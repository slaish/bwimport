bwimport

R package that uses C++ to import data from bigwig files either local or remote

The R package should work out-of-the-box, so the Rccp link should be static.
[Makevars should hard-lock compilation to the same Rcpp used by the current R, right?]

Working directory:
/Users/au103725/Dropbox/Lab_stuff/Bioinformatics/R_packages

Folder structure:

bwimport/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   └── bw_import.R
└── src/
    ├── bw_import.cpp
    ├── Makevars
    └── libBigWig.a


Here are the file contents:

DESCRIPTION
Package: bwimport
Type: Package
Title: Fast BigWig Region Import
Version: 0.1.0
Authors@R: person("Søren", "Lykke-Andersen", email = "sla@mbg.au.dk", role = c("cre", "aut"), comment = c(ORCID = "0000-0001-9357-2910"))
Description: Import numeric base-resolution values from a BigWig region using libBigWig.
License: MIT + file LICENSE
Encoding: UTF-8
LinkingTo: Rcpp
Imports: Rcpp (>= 1.0.10)
SystemRequirements: zlib, libcurl
LazyData: true

NAMESPACE
useDynLib(bwimport, .registration=TRUE)
importFrom(Rcpp, evalCpp)
export(bw_import)

bw_import.R
#' Import BigWig region
#'
#' @param bw_file path or URL
#' @param chrom chromosome (e.g., "chr12")
#' @param start 1-based start
#' @param end 1-based end inclusive
#' @return numeric vector of length end-start+1
#' @export
bw_import <- function(bw_file, chrom, start, end) {
  stopifnot(is.character(bw_file), is.character(chrom), length(bw_file) == 1L, length(chrom) == 1L)
  bw_import_impl(bw_file, chrom, as.integer(start), as.integer(end))
}

bw_import.cpp
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
 
Makevars
PKG_CPPFLAGS += -I. $(shell "${R_HOME}/bin/Rscript" -e "cat(Rcpp:::CxxFlags())")
PKG_LIBS     += libBigWig.a -lz -lcurl $(shell "${R_HOME}/bin/Rscript" -e "cat(Rcpp:::LdFlags())")

The package was installed from inside R forcing a rebuild against the currently loaded R and Rcpp so that symbols match exactly (a problem that took a long time to fix):
system("R CMD INSTALL --preclean --no-multiarch .")


Sanity check:
getLoadedDLLs()[["bwimport"]]

DLL name: bwimport
Filename: /Users/au103725/Library/R/arm64/4.5/library/bwimport/libs/bwimport.so
Dynamic lookup: TRUE


> library(bwimport)
> .bw.file = "/Users/au103725/Dropbox/Lab_stuff/Bioinformatics/R_packages/seqNdisplayR/test/siRRP40_xPAP_in_batch3_plus.bw"
> .chrom = "chr12"
> .chrom.start=6531808
> .chrom.end=6541078
> bw_import(.bw.file, .chrom, .chrom.start, .chrom.end)
Error in bw_import(.bw.file, .chrom, .chrom.start, .chrom.end) : 
  function 'Rcpp_precious_remove' not provided by package 'Rcpp'


It seems to me that the static/dynamic issue is mixed up in this code.
Moreover, based on previous troubleshooting "Rcpp_precious_remove" was present, but with some additional characters "__ZN4Rcpp20Rcpp_precious_removeEP7SEXPREC". I don't know if this is same in this version, but it could be.




cmake .. \
  -DBUILD_SHARED_LIBS=OFF \
  -DWITH_CURL=ON \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DCMAKE_OSX_ARCHITECTURES=arm64 \
  -DCURL_INCLUDE_DIR="$(brew --prefix curl)/include" \
  -DCURL_LIBRARY="$(brew --prefix curl)/lib/libcurl.dylib" \
  -DZLIB_INCLUDE_DIR="$(brew --prefix zlib)/include" \
  -DZLIB_LIBRARY="$(brew --prefix zlib)/lib/libz.dylib"