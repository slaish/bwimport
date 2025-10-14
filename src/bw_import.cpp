// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>

extern "C" {
  #include "bigWig.h"
}

using namespace Rcpp;

// Helper: strip any leading "chr", case-insensitive
inline std::string strip_chr_prefix(const std::string& s) {
  if (s.size() >= 3 && 
     (s.rfind("chr", 0) == 0 || s.rfind("Chr", 0) == 0 || s.rfind("CHR", 0) == 0)) {
    return s.substr(3);
  }
  return s;
}

// [[Rcpp::export]]
NumericVector bw_import_impl(std::string bw_file, std::string chrom, int start, int end) {
  if (start < 1 || end < start)
    stop("Invalid coordinates: start must be >= 1 and end >= start.");

  if (bwInit(1 << 17) != 0)
    stop("Failed to initialize libBigWig.");

  bigWigFile_t* bw = bwOpen(bw_file.c_str(), NULL, "r");
  if (!bw)
    stop("Cannot open BigWig file: %s", bw_file.c_str());

  // --- Chromosome name matching (handles both 'chr12' ↔ '12') ---
  std::string chrom_match = "";
  std::string chrom_stripped = strip_chr_prefix(chrom);

  if (bw->cl && bw->cl->nKeys > 0) {
    for (uint32_t i = 0; i < bw->cl->nKeys; ++i) {
      std::string bw_chrom = bw->cl->chrom[i];
      std::string bw_chrom_stripped = strip_chr_prefix(bw_chrom);

      // Match if they are identical *with or without* the "chr" prefix
      if (chrom == bw_chrom || chrom_stripped == bw_chrom_stripped) {
        chrom_match = bw_chrom;  // use actual name from BigWig
        break;
      }
    }
  }

  if (chrom_match.empty()) {
    // Could not find any match — build a short list of what's in the file
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

