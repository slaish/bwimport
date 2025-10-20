#' Import BigWig region
#'
#' @param bw_file Character scalar: path to a local BigWig or a URL (http/https/ftp)
#' @param chrom   Character scalar: chromosome name (e.g., "chr1")
#' @param start   Integer(1): 1-based start (inclusive)
#' @param end     Integer(1): 1-based end (inclusive)
#' @return Numeric vector of length end - start + 1
#' @export
#' @examples
#' bw_URL = "http://genome-ftp.mbg.au.dk/public/THJ/seqNdisplayR/examples/tracks/HeLa_3pseq/siGFP_noPAP_in_batch1_plus.bw"
#' chrom = "chr12"
#' chrom_start = 6531808
#' chrom_end = 6541078
#' vals = bw_import(bw_URL, chrom, chrom_start, chrom_end)
#' max(vals) # [1] 9503.29
#' 
#' On Windows - tune import timing through the following
#' 1) Bigger buffer
#' Sys.setenv(BWIMPORT_BUFSZ_KB = "4096")  # 4 MiB
#' 2) Force IPv4
#' Sys.setenv(BWIMPORT_IPRESOLVE = "4")
#' 3) Optional: fresh sockets (only if still slow)
#' # Sys.setenv(BWIMPORT_FRESH = "1")
#' If there is a big drop, BUFSZ can be dialed down (to e.g. 2048) to balance memory/throughput.
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
