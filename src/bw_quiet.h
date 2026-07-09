#ifndef BW_QUIET_H
#define BW_QUIET_H

#include <stdio.h>
#include <stdlib.h>

/* Gate all libBigWig stderr messages behind the BWIMPORT_QUIET env var.
 * When BWIMPORT_QUIET is set to a non-empty value, suppress the noisy
 * transient "[bwHdrRead] ...", "[bwOpen] ..." etc. that fire on retryable
 * failures. Unset OR empty means "print" -- matches R's Sys.getenv/nzchar
 * semantics, so `Sys.setenv(BWIMPORT_QUIET = "")` re-enables diagnostics. */
#define BW_STDERR(...) do {                                          \
    const char *_bw_q = getenv("BWIMPORT_QUIET");                    \
    if (!_bw_q || !*_bw_q) {                                         \
        fprintf(stderr, __VA_ARGS__);                                \
    }                                                                \
} while (0)

#endif /* BW_QUIET_H */
