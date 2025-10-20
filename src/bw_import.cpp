// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <list>
#include <memory>

extern "C" {
  #include "bigWig.h"
  #include <R_ext/Rdynload.h> // for R_unload_* signature
}

using namespace Rcpp;

// -------------------- one-time libBigWig init (env-aware) --------------------
static std::atomic<bool> bw_ready(false);

inline uint32_t pick_buffer_bytes_from_env() {
  size_t kb = 0;
  if (const char* s = std::getenv("BWIMPORT_BUFSZ_KB")) {
    long v = std::strtol(s, nullptr, 10);
    if (v > 0) kb = static_cast<size_t>(v);
  }
  size_t bytes = kb ? kb * 1024u : (1u << 20); // default 1 MiB
  if (bytes < (1u << 16)) bytes = (1u << 16);  // min 64 KiB
  if (bytes > (1u << 23)) bytes = (1u << 23);  // max 8 MiB
  return static_cast<uint32_t>(bytes);
}

inline void ensure_bw_init() {
  bool expected = false;
  if (bw_ready.compare_exchange_strong(expected, true)) {
    uint32_t buf = pick_buffer_bytes_from_env();
    if (bwInit(buf) != 0) {
      bw_ready.store(false);
      stop("Failed to initialize libBigWig.");
    }
  }
}

// -------------------- small LRU cache of open handles -----------------------
struct BWHandle {
  bigWigFile_t* bw{nullptr};
  // cache a mapping from stripped->actual chrom name for faster matching
  std::unordered_map<std::string, std::string> chrom_map;
};

static std::mutex g_mtx;
static std::unordered_map<std::string, std::pair<BWHandle*, std::list<std::string>::iterator>> g_map;
static std::list<std::string> g_lru; // front = most recently used
static std::vector<std::unique_ptr<BWHandle>> g_pool;

inline size_t max_open_handles() {
  size_t n = 8;
  if (const char* s = std::getenv("BWIMPORT_MAX_OPEN")) {
    long v = std::strtol(s, nullptr, 10);
    if (v > 0) n = static_cast<size_t>(v);
  }
  if (n < 1) n = 1;
  if (n > 64) n = 64;
  return n;
}

inline std::string strip_chr_prefix(const std::string& s) {
  if (s.size() >= 3) {
    char c0 = s[0], c1 = s[1], c2 = s[2];
    if ((c0=='c'||c0=='C') && (c1=='h'||c1=='H') && (c2=='r'||c2=='R'))
      return s.substr(3);
  }
  return s;
}

inline void close_and_free(BWHandle* h) {
  if (!h) return;
  if (h->bw) { bwClose(h->bw); h->bw = nullptr; }
  h->chrom_map.clear();
}

inline BWHandle* open_or_get_cached(const std::string& key) {
  std::lock_guard<std::mutex> lk(g_mtx);

  // hit
  auto it = g_map.find(key);
  if (it != g_map.end()) {
    // move to front
    g_lru.erase(it->second.second);
    g_lru.push_front(key);
    it->second.second = g_lru.begin();
    return it->second.first;
  }

  // miss â†’ open
  auto h = std::make_unique<BWHandle>();
  h->bw = bwOpen(key.c_str(), NULL, "r");
  if (!h->bw) {
    // do not cache failed opens
    return nullptr;
  }

  // build chrom map (stripped -> actual, and actual -> actual)
  if (h->bw->cl && h->bw->cl->nKeys > 0) {
    for (uint32_t i = 0; i < h->bw->cl->nKeys; ++i) {
      std::string actual = h->bw->cl->chrom[i];
      h->chrom_map[strip_chr_prefix(actual)] = actual; // "12" -> "chr12" or "12"
      h->chrom_map[actual] = actual;                   // "chr12" -> "chr12"
    }
  }

  // evict if needed
  if (g_pool.size() >= max_open_handles()) {
    // evict least-recently-used
    const std::string& victim_key = g_lru.back();
    auto mit = g_map.find(victim_key);
    if (mit != g_map.end()) {
      close_and_free(mit->second.first);
      // keep the slot but reuse it for new handle
      // replace the unique_ptr content
      // find which unique_ptr holds this BWHandle*
      for (auto& up : g_pool) {
        if (up.get() == mit->second.first) { up.reset(); break; }
      }
      g_map.erase(mit);
    }
    g_lru.pop_back();
  }

  g_pool.emplace_back(std::move(h));
  BWHandle* ret = g_pool.back().get();
  g_lru.push_front(key);
  g_map.emplace(key, std::make_pair(ret, g_lru.begin()));
  return ret;
}

inline void clear_cache() {
  std::lock_guard<std::mutex> lk(g_mtx);
  for (auto& up : g_pool) {
    close_and_free(up.get());
  }
  g_pool.clear();
  g_map.clear();
  g_lru.clear();
}

// ----------------------------- main import ----------------------------------

// fill a constant span quickly
inline void fill_span(NumericVector& out, uint32_t qStart, uint32_t s, uint32_t e, double v) {
  uint32_t paintStart = std::max(s, qStart);
  uint32_t paintEnd   = e; // caller ensures e <= qEnd
  if (paintEnd <= paintStart) return;
  size_t off  = static_cast<size_t>(paintStart - qStart);
  size_t len  = static_cast<size_t>(paintEnd   - paintStart);
  // zero/NaN handling
  double val = std::isnan(v) ? 0.0 : v;
  std::fill(out.begin() + off, out.begin() + off + len, val);
}

// [[Rcpp::export]]
NumericVector bw_import_impl(std::string bw_file, std::string chrom, int start, int end) {
  if (start < 1 || end < start)
    stop("Invalid coordinates: start must be >= 1 and end >= start.");

  ensure_bw_init();

  BWHandle* H = open_or_get_cached(bw_file);
  if (!H || !H->bw)
    stop("Cannot open BigWig file: %s", bw_file.c_str());

  // resolve chromosome with cache
  std::string want = chrom;
  std::string key  = want;
  auto it = H->chrom_map.find(key);
  if (it == H->chrom_map.end()) {
    key = strip_chr_prefix(want);
    it  = H->chrom_map.find(key);
    if (it == H->chrom_map.end()) {
      // rebuild chrom_map if file was reopened/changed (defensive)
      H->chrom_map.clear();
      if (H->bw->cl && H->bw->cl->nKeys > 0) {
        for (uint32_t i = 0; i < H->bw->cl->nKeys; ++i) {
          std::string actual = H->bw->cl->chrom[i];
          H->chrom_map[strip_chr_prefix(actual)] = actual;
          H->chrom_map[actual] = actual;
        }
      }
      it = H->chrom_map.find(key);
      if (it == H->chrom_map.end()) {
        // build preview list for error
        std::string available;
        if (H->bw->cl && H->bw->cl->nKeys > 0) {
          uint32_t lim = std::min<uint32_t>(H->bw->cl->nKeys, 5);
          for (uint32_t i = 0; i < lim; ++i) {
            if (i) available += ", ";
            available += H->bw->cl->chrom[i];
          }
          if (H->bw->cl->nKeys > lim) available += ", ...";
        }
        stop("Chromosome '%s' not found in BigWig file '%s'. Available examples: [%s]",
             chrom.c_str(), bw_file.c_str(), available.c_str());
      }
    }
  }
  const std::string& chrom_match = it->second;

  // query bounds
  const uint32_t qStart = static_cast<uint32_t>(start - 1); // 0-based inclusive
  const uint32_t qEnd   = static_cast<uint32_t>(end);       // 0-based exclusive
  const int out_len = end - start + 1;
  NumericVector out(out_len, 0.0);

  // fetch intervals
  bwOverlappingIntervals_t* iv = bwGetOverlappingIntervals(H->bw, chrom_match.c_str(), qStart, qEnd);
  if (iv) {
    for (uint32_t i = 0; i < iv->l; ++i) {
      uint32_t s = iv->start[i];
      uint32_t e = iv->end[i];
      if (e <= qStart || s >= qEnd) continue;
      if (e > qEnd) e = qEnd;
      fill_span(out, qStart, s, e, iv->value[i]);
    }
    bwDestroyOverlappingIntervals(iv);
  }
  return out;
}

// [[Rcpp::export]]
void bw_cleanup() {
  clear_cache();
  bwCleanup();
  bw_ready.store(false);
}

// called when the DLL/SO unloads
extern "C" void R_unload_bwimport(DllInfo* /*dll*/) {
  clear_cache();
  bwCleanup();
  bw_ready.store(false);
}
