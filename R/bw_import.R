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
