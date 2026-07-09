# bwimport 0.2.3

## Bug fixes

* **Windows: URLs no longer mangled to `http:\\...`.** `safe_local_path()`
  in `src/bw_import.cpp` used to unconditionally convert `/` → `\` on
  Windows to make local paths fopen()-safe. That also mangled `http://`
  into `http:\\...`, and libBigWig then failed to parse it as a URL --
  every remote track failed silently, and some hosts caused a ~50 s
  hang per file (fopen retrying an invalid Windows path with the URL
  host treated as a UNC share). The fix restricts the slash-swap to
  paths that don't start with `http://`, `https://`, or `ftp://`,
  matching how R-level `bw_import()` already routes URLs vs local
  paths. Remote plots on Windows now hit libBigWig's byte-range fast
  path directly, same as macOS/Linux.

# bwimport 0.2.2

## Bug fixes

* Windows: dramatically faster remote BigWig reads. `bw_import()` on
  Windows used to always download the entire remote file to disk
  before opening it, even though `Makevars.win` links libcurl into
  libBigWig and byte-range fetches would work. First plots with many
  tracks paid the transfer cost of every full bigwig, appearing to
  hang for minutes on slow connections.

  Now: on Windows + URL, `bw_import()` tries libBigWig's own
  URL-open first (byte-range reads, typically ~50 kB per region).
  Only on failure does it fall back to the old download path -- and
  when it does, downloaded files are cached per URL for the R session,
  so at most one transfer per bigwig.

  Set `Sys.setenv(BWIMPORT_WINDOWS_DOWNLOAD = "1")` to force the
  download path if you're behind a proxy where libcurl fails but
  R's `curl` package succeeds.

## New features

* New exported helper `bw_clear_url_cache()` releases any bigwigs the
  Windows fallback path cached under `tempdir()` during the session.

# bwimport 0.2.1

## New features

* All libBigWig `stderr` messages are now gated behind the `BWIMPORT_QUIET`
  environment variable. When `BWIMPORT_QUIET` is set (to any value), the
  vendored libBigWig prints (`[bwHdrRead] ...`, `[bwOpen] ...`, etc.) are
  suppressed at the C level. Leave the variable unset to see full diagnostics.
  This is intended for callers that already handle failures at the R level
  (retry with alt chromosome name, all-zero fallback) and don't want the
  cosmetic noise from transient first-attempt failures.

# bwimport 0.2.0

## New features

* `bw_import()` now works on Windows, including importing BigWig files from
  remote URLs (http/https/ftp). Remote files are downloaded to a temporary
  location before reading, since the Windows build of libBigWig is configured
  without libcurl.
* New exported function `bw_cleanup()` releases libBigWig resources.

## Bug fixes

* Fixed a shadow-variable bug in `bw_import()` and improved Windows path
  handling (forward slashes, short-path expansion to avoid spaces).
* Cleaned up the `\examples{}` block in `?bw_import` — removed prose that was
  mixed into example R code.

## Build / internals

* Added `configure.win` and reworked `Makevars.win` so the package builds
  cleanly on Windows.
* Added `curl` to `Imports`.
* `bw_import()` now uses one-shot libBigWig initialisation with proper DLL
  cleanup, and continues to match chromosome names with or without the `chr`
  prefix.
