# Per-session cache: URL -> path of a downloaded copy on local disk. Keeps a
# remote bigwig alive for the whole R session so subsequent bw_import() calls
# on the same URL reuse the local copy instead of re-downloading. Cleared on
# session exit (files under tempdir()).
.bw_url_cache <- new.env(parent = emptyenv())

#' Import BigWig region
#'
#' @param bw_file Character scalar: path to a local BigWig or a URL (http/https/ftp)
#' @param chrom   Character scalar: chromosome name (e.g., "chr1")
#' @param start   Integer(1): 1-based start (inclusive)
#' @param end     Integer(1): 1-based end (inclusive)
#' @return Numeric vector of length end - start + 1
#'
#' @details
#' On Windows, remote URLs used to be *always* downloaded to a temp file
#' before being passed to libBigWig, which caused very slow first plots
#' (whole-file transfer per bigwig instead of a ~50 kB byte-range fetch).
#' As of bwimport 0.2.2:
#'
#' - The default is direct URL fetching (libBigWig's own libcurl support
#'   is linked in `Makevars.win`). If direct fetching fails at open time,
#'   the code silently falls back to a download-then-open path.
#' - Downloaded files are cached per URL for the R session, so the fallback
#'   pays the transfer cost at most once per bigwig even in the failure
#'   case.
#' - Force the download path (skip direct fetch) with
#'   `Sys.setenv(BWIMPORT_WINDOWS_DOWNLOAD = "1")` -- useful behind
#'   restrictive corporate proxies where libcurl can't reach GitHub /
#'   remote hosts but R's `curl` package can.
#' @export
#' @examples
#' bw_URL <- "http://genome-ftp.mbg.au.dk/public/THJ/seqNdisplayR/examples/tracks/HeLa_3pseq/siGFP_noPAP_in_batch1_plus.bw"
#' chrom <- "chr12"
#' chrom_start <- 6531808
#' chrom_end <- 6541078
#' vals <- bw_import(bw_URL, chrom, chrom_start, chrom_end)
#' max(vals) # [1] 9503.29
#'
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

  is_url  <- grepl("^(https?|ftp)://", bw_file, ignore.case = TRUE)
  is_win  <- .Platform$OS.type == "windows"

  # Local Windows paths: normalise so libBigWig's fopen sees clean slashes.
  if (is_win && !is_url) {
    return(bwimport::bw_import_impl(
      normalizePath(bw_file, winslash = "/", mustWork = TRUE),
      chrom, start, end
    ))
  }

  # Non-Windows, or non-URL: direct pass-through.
  if (!is_win || !is_url) {
    return(bwimport::bw_import_impl(bw_file, chrom, start, end))
  }

  # Windows + URL. Two paths, tried in order:
  #   (1) Direct URL fetch via libBigWig+libcurl -- fast byte-range reads.
  #   (2) Fall back to downloading the whole file (cached per session) --
  #       correct but slow, only pays the transfer cost once per URL.
  force_download <- nzchar(Sys.getenv("BWIMPORT_WINDOWS_DOWNLOAD"))

  if (!force_download) {
    direct <- tryCatch(
      bwimport::bw_import_impl(bw_file, chrom, start, end),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
    if (!is.null(direct)) return(direct)
    # Fall through to download path.
  }

  # Download-and-open fallback. Cache the local copy per URL for the session.
  local_path <- get0(bw_file, envir = .bw_url_cache, inherits = FALSE)
  if (is.null(local_path) || !file.exists(local_path)) {
    td <- tempdir()
    short_td <- tryCatch(utils::shortPathName(td), error = function(e) td)
    local_path <- tempfile(tmpdir = short_td, fileext = ".bw")
    tryCatch(
      curl::curl_download(url = bw_file, destfile = local_path, mode = "wb"),
      error = function(e) {
        stop(sprintf(
          paste0("Failed to download BigWig file from URL:\n  %s\nError: %s\n\n",
                 "Tip: check your internet connection; try Sys.setenv(",
                 "BWIMPORT_IPRESOLVE = \"4\") to force IPv4."),
          bw_file, conditionMessage(e)
        ), call. = FALSE)
      }
    )
    assign(bw_file, local_path, envir = .bw_url_cache)
  }

  bwimport::bw_import_impl(local_path, chrom, start, end)
}

#' Clear the per-session bigwig URL download cache
#'
#' @description
#' Deletes any bigwig files bwimport downloaded to `tempdir()` on Windows
#' as a fallback for URLs that libBigWig couldn't fetch directly, and
#' empties the URL -> local-path map. Safe to call any time; no-op on
#' non-Windows or when nothing has been cached.
#' @return `invisible(NULL)`.
#' @export
bw_clear_url_cache <- function() {
  for (url in ls(.bw_url_cache, all.names = TRUE)) {
    p <- get(url, envir = .bw_url_cache)
    try(unlink(p), silent = TRUE)
    rm(list = url, envir = .bw_url_cache)
  }
  invisible(NULL)
}
