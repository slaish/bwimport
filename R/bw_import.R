#' Import BigWig region
#'
#' @param bw_file Character scalar: path to a local BigWig or a URL (http/https/ftp)
#' @param chrom   Character scalar: chromosome name (e.g., "chr1")
#' @param start   Integer(1): 1-based start (inclusive)
#' @param end     Integer(1): 1-based end (inclusive)
#' @return Numeric vector of length end - start + 1
#' @export
bw_import <- function(bw_file, chrom, start, end) {
  stopifnot(
    is.character(bw_file), length(bw_file) == 1L,
    is.character(chrom),   length(chrom)   == 1L
  )
  start <- as.integer(start)
  end   <- as.integer(end)
  if (!is.finite(start) || !is.finite(end) || start < 1L || end < start) {
    stop("Invalid coordinates: start must be >= 1 and end >= start.", call. = FALSE)
  }
 
  is_url <- grepl("^(https?|ftp)://", bw_file, ignore.case = TRUE)
 
  # On Windows we build libBigWig without libcurl; download URL to a temp file first.
  if (.Platform$OS.type == "windows" && is_url) {
    tf <- tempfile(fileext = ".bw")
    curl::curl_download(url = bw_file, destfile = tf, mode = "wb")
    on.exit(try(unlink(tf), silent = TRUE), add = TRUE)
    bw_file <- tf
  }
 
  bw_import_impl(bw_file, chrom, start, end)
}
