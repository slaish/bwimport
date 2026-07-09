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
