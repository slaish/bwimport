// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <atomic>

#ifdef _WIN32
#include <windows.h>  // for GetShortPathNameA
#endif

extern "C" {
  #include "bigWig.h"
  #include <R_ext/Rdynload.h>
}

using namespace Rcpp;

// --- one-time libBigWig init guard ------------------------------------------
// FIX: single atomic flag at file scope — no shadowing
static std::atomic<bool> bw_ready(false);

inline void ensure_bw_init() {
  // FIX: removed the local `static std::atomic<bool> bw_ready` that was
  //      shadowing the file-scope variable. bw_cleanup() and R_unload_bwimport()
  //      reset the file-scope bw_ready, so the guard must use the same variable.
  bool expected = false;
  if (bw_ready.compare_exchange_strong(expected, true)) {
    size_t kb = 0;
    if (const char* s = std::getenv("BWIMPORT_BUFSZ_KB")) {
      long v = std::strtol(s, nullptr, 10);
      if (v > 0) kb = static_cast<size_t>(v);
    }
    size_t bytes = kb ? kb * 1024u : (1u << 20); // 1 MiB default
    if (bytes < (1u << 16)) bytes = (1u << 16);  // min 64 KiB
    if (bytes > (1u << 23)) bytes = (1u << 23);  // max 8 MiB

    if (bwInit(bytes) != 0) {
      bw_ready.store(false);
      Rcpp::stop("Failed to initialize libBigWig.");
    }
  }
}


// Helper: strip any leading "chr", case-insensitive
inline std::string strip_chr_prefix(const std::string& s) {
  if (s.size() >= 3 &&
     (s.rfind("chr", 0) == 0 || s.rfind("Chr", 0) == 0 || s.rfind("CHR", 0) == 0)) {
    return s.substr(3);
  }
  return s;
}


// FIX: On Windows, convert paths with spaces (or non-ASCII chars) to the
//      8.3 short-path form so that C fopen() inside libBigWig never chokes.
//      Also normalise forward slashes — harmless but defensive.
inline std::string safe_local_path(const std::string& path) {
#ifdef _WIN32
  // Normalise separators: forward slashes → backslashes
  std::string norm = path;
  std::replace(norm.begin(), norm.end(), '/', '\\');

  // If the path contains a space (or any extended char), try the 8.3 form
  if (norm.find(' ') != std::string::npos ||
      norm.find('~') != std::string::npos) {
    char shortbuf[MAX_PATH];
    DWORD len = GetShortPathNameA(norm.c_str(), shortbuf, MAX_PATH);
    if (len > 0 && len < MAX_PATH) {
      return std::string(shortbuf);
    }
    // Fall through: short-name generation may be disabled on the volume.
    // In that case, just return the normalised path and hope for the best.
  }
  return norm;
#else
  return path;  // no-op on Unix
#endif
}


// [[Rcpp::export]]
NumericVector bw_import_impl(std::string bw_file, std::string chrom, int start, int end) {
  if (start < 1 || end < start)
    stop("Invalid coordinates: start must be >= 1 and end >= start.");

  ensure_bw_init();

  // FIX: sanitise local paths on Windows before handing to libBigWig
  std::string open_path = safe_local_path(bw_file);

  bigWigFile_t* bw = bwOpen(open_path.c_str(), NULL, "r");
  if (!bw)
    stop("Cannot open BigWig file: %s", bw_file.c_str());

  // --- Chromosome name matching (handles both 'chr12' <-> '12') ---
  std::string chrom_match = "";
  std::string chrom_stripped = strip_chr_prefix(chrom);

  if (bw->cl && bw->cl->nKeys > 0) {
    for (uint32_t i = 0; i < bw->cl->nKeys; ++i) {
      std::string bw_chrom = bw->cl->chrom[i];
      std::string bw_chrom_stripped = strip_chr_prefix(bw_chrom);
      if (chrom == bw_chrom || chrom_stripped == bw_chrom_stripped) {
        chrom_match = bw_chrom;
        break;
      }
    }
  }

  if (chrom_match.empty()) {
    std::string available = "";
    if (bw->cl && bw->cl->nKeys > 0) {
      for (uint32_t i = 0; i < std::min<uint32_t>(bw->cl->nKeys, 5); ++i) {
        available += bw->cl->chrom[i];
        if (i < std::min<uint32_t>(bw->cl->nKeys, 5) - 1) available += ", ";
      }
      if (bw->cl->nKeys > 5) available += ", ...";
    }
    bwClose(bw);
    stop("Chromosome '%s' not found in BigWig file '%s'. Available examples: [%s]",
         chrom.c_str(), bw_file.c_str(), available.c_str());
  }

  // --- Query region ---
  const uint32_t qStart = static_cast<uint32_t>(start - 1); // 0-based inclusive
  const uint32_t qEnd   = static_cast<uint32_t>(end);       // 0-based exclusive

  const int out_len = end - start + 1;
  NumericVector out(out_len, 0.0);

  bwOverlappingIntervals_t* iv = bwGetOverlappingIntervals(bw, chrom_match.c_str(), qStart, qEnd);
  if (iv) {
    for (uint32_t i = 0; i < iv->l; ++i) {
      uint32_t s = iv->start[i];
      uint32_t e = iv->end[i];
      if (e <= qStart || s >= qEnd) continue;

      uint32_t paintStart = std::max(s, qStart);
      uint32_t paintEnd   = std::min(e, qEnd);

      for (uint32_t p = paintStart; p < paintEnd; ++p) {
        int idx = static_cast<int>(p - qStart);
        double v = iv->value[i];
        out[idx] = std::isnan(v) ? 0.0 : v;
      }
    }
    bwDestroyOverlappingIntervals(iv);
  }

  bwClose(bw);
  return out;
}

// [[Rcpp::export]]
void bw_cleanup() {
  bwCleanup();
  bw_ready.store(false);  // now correctly resets the SAME flag ensure_bw_init() checks
}

// --- Cleanup hook: called when the DLL/SO unloads ---------------------------
extern "C" void R_unload_bwimport(DllInfo* /*dll*/) {
  bwCleanup();
  bw_ready.store(false);
}
