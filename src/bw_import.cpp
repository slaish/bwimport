// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
 
extern "C" {
  #include "bigWig.h"
}
 
using namespace Rcpp;
 
// Fill [start,end] inclusive in 1-based R coordinates
// with values from a BigWig. Uses bwGetOverlappingIntervals
// and paints ranges into the output buffer.
// [[Rcpp::export]]
NumericVector bw_import_impl(std::string bw_file, std::string chrom, int start, int end) {
  if (start < 1 || end < start) {
    stop("Invalid coordinates: start must be >= 1 and end >= start.");
  }
 
  // Init libBigWig once-per-process (harmless if called multiple times)
  if (bwInit(1 << 17) != 0) {
    stop("Failed to initialize libBigWig.");
  }
 
  bigWigFile_t* bw = bwOpen(bw_file.c_str(), NULL, "r");
  if (!bw) {
    stop("Cannot open BigWig file: %s", bw_file.c_str());
  }
 
  const uint32_t qStart = static_cast<uint32_t>(start - 1); // 0-based inclusive
  const uint32_t qEnd   = static_cast<uint32_t>(end);       // 0-based exclusive
 
  // Pre-fill with zeros (or NA_REAL if you prefer)
  const int out_len = end - start + 1;
  NumericVector out(out_len, 0.0);
 
  // Get overlapping intervals, not per-base values
  bwOverlappingIntervals_t* iv = bwGetOverlappingIntervals(bw, chrom.c_str(), qStart, qEnd);
  if (iv) {
    for (uint32_t i = 0; i < iv->l; ++i) {
      // Interval [iv->start[i], iv->end[i]) with value iv->value[i]
      uint32_t s = iv->start[i];
      uint32_t e = iv->end[i];
 
      if (e <= qStart || s >= qEnd) continue; // no overlap
 
      uint32_t paintStart = s < qStart ? qStart : s;
      uint32_t paintEnd   = e > qEnd   ? qEnd   : e;
 
      // Paint inclusive of start, exclusive of end
      for (uint32_t p = paintStart; p < paintEnd; ++p) {
        int idx = static_cast<int>(p - qStart); // 0-based index into out
        double v = iv->value[i];
        out[idx] = std::isnan(v) ? 0.0 : v;     // keep zero for missing
      }
    }
    bwDestroyOverlappingIntervals(iv);
  }
 
  bwClose(bw);
  return out;
}
