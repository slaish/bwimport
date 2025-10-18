#ifndef NOCURL
#include <curl/curl.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "bigWigIO.h"
#include <inttypes.h>
#include <errno.h>
 
size_t GLOBAL_DEFAULTBUFFERSIZE;
 
#ifndef NOCURL
/* Apply common curl options before each perform() */
static void bw_curl_apply_common_opts(CURL *h) {
    /* Force HTTP/1.1; some Windows stacks/proxies stall on HTTP/2 */
    curl_easy_setopt(h, CURLOPT_HTTP_VERSION, CURL_HTTP_VERSION_1_1);
 
    /* Follow redirects + sane timeouts */
    curl_easy_setopt(h, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(h, CURLOPT_MAXREDIRS,     5L);
    curl_easy_setopt(h, CURLOPT_CONNECTTIMEOUT,10L);  /* seconds */
    curl_easy_setopt(h, CURLOPT_TIMEOUT,       60L);  /* seconds total */
 
    /* Keepalive for long sequences of range requests */
    curl_easy_setopt(h, CURLOPT_TCP_KEEPALIVE, 1L);
    curl_easy_setopt(h, CURLOPT_TCP_KEEPIDLE,  30L);
    curl_easy_setopt(h, CURLOPT_TCP_KEEPINTVL, 15L);
 
    /* Safer in threaded/Windows envs */
    curl_easy_setopt(h, CURLOPT_NOSIGNAL, 1L);
 
    /* Last resort if a proxy/CDN misbehaves:
       curl_easy_setopt(h, CURLOPT_FORBID_REUSE, 1L);
       curl_easy_setopt(h, CURLOPT_FRESH_CONNECT, 1L);
    */
}
 
uint64_t getContentLength(const URL_t *URL) {
    size_t size = 0;
    if (curl_easy_getinfo(URL->x.curl, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &size) != CURLE_OK) {
        return 0;
    }
    if (size == (size_t)-1) return 0;
    return (uint64_t)size;
}
 
/* Fill the buffer; URL may be left unusable on error */
CURLcode urlFetchData(URL_t *URL, unsigned long bufSize) {
    CURLcode rv;
    char range[128];
 
    if (URL->filePos != (size_t)-1) URL->filePos += URL->bufLen;
    else URL->filePos = 0;
 
    URL->bufPos = URL->bufLen = 0; /* reset buffer window */
 
    /* Build Range header using %zu for size_t */
    (void)snprintf(range, sizeof(range), "%zu-%zu",
                   URL->filePos, URL->filePos + (size_t)bufSize - 1U);
    rv = curl_easy_setopt(URL->x.curl, CURLOPT_RANGE, range);
    if (rv != CURLE_OK) {
        fprintf(stderr, "[urlFetchData] Couldn't set the range (%s)\n", range);
        return rv;
    }
 
    /* Apply common opts then perform */
    bw_curl_apply_common_opts(URL->x.curl);
 
    rv = curl_easy_perform(URL->x.curl);
    errno = 0; /* clear remnant errno */
    return rv;
}
 
/* Read data into a buffer, ideally from an in-memory buffer */
size_t url_fread(void *obuf, size_t obufSize, URL_t *URL) {
    size_t remaining = obufSize, fetchSize;
    unsigned char *p = (unsigned char*)obuf;
    CURLcode rv;
 
    while (remaining) {
        if (!URL->bufLen) {
            rv = urlFetchData(URL, URL->bufSize);
            if (rv != CURLE_OK) {
                fprintf(stderr, "[url_fread] urlFetchData (A) returned %s\n", curl_easy_strerror(rv));
                return 0;
            }
        } else if (URL->bufLen < URL->bufPos + remaining) {
            size_t chunk = URL->bufLen - URL->bufPos;
            memcpy(p, URL->memBuf + URL->bufPos, chunk);
            p += chunk;
            remaining -= chunk;
 
            if (remaining) {
                fetchSize = URL->isCompressed
                            ? ((remaining < URL->bufSize) ? remaining : URL->bufSize)
                            : URL->bufSize;
                rv = urlFetchData(URL, fetchSize);
                if (rv != CURLE_OK) {
                    fprintf(stderr, "[url_fread] urlFetchData (B) returned %s\n", curl_easy_strerror(rv));
                    return 0;
                }
            }
        } else {
            memcpy(p, URL->memBuf + URL->bufPos, remaining);
            URL->bufPos += remaining;
            remaining = 0;
        }
    }
    return obufSize;
}
#endif /* !NOCURL */
 
/* Returns bytes requested or fewer on error.
   For remote files, actual bytes read may be less than the return value. */
size_t urlRead(URL_t *URL, void *buf, size_t bufSize) {
#ifndef NOCURL
    if (URL->type == 0) {
        return fread(buf, bufSize, 1, URL->x.fp) * bufSize;
    } else {
        return url_fread(buf, bufSize, URL);
    }
#else
    return fread(buf, bufSize, 1, URL->x.fp) * bufSize;
#endif
}
 
size_t bwFillBuffer(const void *inBuf, size_t l, size_t nmemb, void *pURL) {
    URL_t *URL = (URL_t*)pURL;
    size_t copied = l * nmemb;
    if (!URL->memBuf) return 0;
 
    if (copied > URL->bufSize - URL->bufLen) /* cap to remaining buffer space */
        copied = URL->bufSize - URL->bufLen;
 
    memcpy(URL->memBuf + URL->bufLen, inBuf, copied);
    URL->bufLen += copied;
    return copied;
}
 
/* Seek to an arbitrary location, returning a CURLcode.
   Local file: CURLE_OK on success or CURLE_FAILED_INIT on error. */
CURLcode urlSeek(URL_t *URL, size_t pos) {
#ifndef NOCURL
    char range[128];
    CURLcode rv;
 
    if (URL->type == BWG_FILE) {
#endif
        if (fseek(URL->x.fp, (long)pos, SEEK_SET) == 0) {
            errno = 0;
            return CURLE_OK;
        } else {
            return CURLE_FAILED_INIT; /* arbitrary */
        }
#ifndef NOCURL
    } else {
        /* If the location is covered by the buffer then don't seek */
        if (pos < URL->filePos || pos >= URL->filePos + URL->bufLen) {
            URL->filePos = pos;
            URL->bufLen = 0; /* so next read won’t increment filePos wrongly */
            URL->bufPos = 0;
 
            (void)snprintf(range, sizeof(range), "%zu-%zu", pos, pos + URL->bufSize - 1U);
            rv = curl_easy_setopt(URL->x.curl, CURLOPT_RANGE, range);
            if (rv != CURLE_OK) {
                fprintf(stderr, "[urlSeek] Couldn't set the range (%s)\n", range);
                return rv;
            }
 
            bw_curl_apply_common_opts(URL->x.curl);
 
            rv = curl_easy_perform(URL->x.curl);
            if (rv != CURLE_OK) {
                fprintf(stderr, "[urlSeek] curl_easy_perform received an error!\n");
            }
            errno = 0;  /* clear remnant errno */
            return rv;
        } else {
            URL->bufPos = pos - URL->filePos;
            return CURLE_OK;
        }
    }
#endif
}
 
URL_t *urlOpen(const char *fname, CURLcode (*callBack)(CURL*), const char *mode) {
    URL_t *URL = (URL_t*)calloc(1, sizeof(URL_t));
    if (!URL) return NULL;
    char *url = NULL, *req = NULL;
#ifndef NOCURL
    CURLcode code;
    char range[128];
#endif
 
    URL->fname = fname;
 
    if ((!mode) || (strchr(mode, 'w') == 0)) {
        /* detect protocol */
#ifndef NOCURL
        if      (strncmp(fname, "http://",  7) == 0) URL->type = BWG_HTTP;
        else if (strncmp(fname, "https://", 8) == 0) URL->type = BWG_HTTPS;
        else if (strncmp(fname, "ftp://",   6) == 0) URL->type = BWG_FTP;
        else                                         URL->type = BWG_FILE;
#else
        URL->type = BWG_FILE;
#endif
 
        if (URL->type == BWG_FILE) {
            /* local file */
            URL->filePos = (size_t)-1; /* nothing read yet */
            URL->x.fp = fopen(fname, "rb");
            if (!(URL->x.fp)) {
                free(URL);
                fprintf(stderr, "[urlOpen] Couldn't open %s for reading\n", fname);
                return NULL;
            }
#ifndef NOCURL
        } else {
            /* remote file */
            URL->memBuf = (unsigned char*)malloc(GLOBAL_DEFAULTBUFFERSIZE);
            if (!(URL->memBuf)) {
                free(URL);
                fprintf(stderr, "[urlOpen] Couldn't allocate file buffer!\n");
                return NULL;
            }
            URL->bufSize = GLOBAL_DEFAULTBUFFERSIZE;
 
            URL->x.curl = curl_easy_init();
            if (!(URL->x.curl)) {
                fprintf(stderr, "[urlOpen] curl_easy_init() failed!\n");
                goto error;
            }
 
            if (curl_easy_setopt(URL->x.curl, CURLOPT_HTTPAUTH, CURLAUTH_ANY) != CURLE_OK) {
                fprintf(stderr, "[urlOpen] Failed to set CURLOPT_HTTPAUTH\n");
                goto error;
            }
            if (curl_easy_setopt(URL->x.curl, CURLOPT_FOLLOWLOCATION, 1L) != CURLE_OK) {
                fprintf(stderr, "[urlOpen] Failed to set CURLOPT_FOLLOWLOCATION\n");
                goto error;
            }
            if (curl_easy_setopt(URL->x.curl, CURLOPT_URL, fname) != CURLE_OK) {
                fprintf(stderr, "[urlOpen] Couldn't set CURLOPT_URL!\n");
                goto error;
            }
 
            /* initial warm-up range */
            (void)snprintf(range, sizeof(range), "0-%zu", URL->bufSize - 1U);
            if (curl_easy_setopt(URL->x.curl, CURLOPT_RANGE, range) != CURLE_OK) {
                fprintf(stderr, "[urlOpen] Couldn't set CURLOPT_RANGE (%s)!\n", range);
                goto error;
            }
 
            /* write callback */
            if (curl_easy_setopt(URL->x.curl, CURLOPT_WRITEFUNCTION, bwFillBuffer) != CURLE_OK) {
                fprintf(stderr, "[urlOpen] Couldn't set CURLOPT_WRITEFUNCTION!\n");
                goto error;
            }
            if (curl_easy_setopt(URL->x.curl, CURLOPT_WRITEDATA, (void*)URL) != CURLE_OK) {
                fprintf(stderr, "[urlOpen] Couldn't set CURLOPT_WRITEDATA!\n");
                goto error;
            }
 
            /* HTTPS: ignore cert errors (keeps upstream behavior) */
            if (curl_easy_setopt(URL->x.curl, CURLOPT_SSL_VERIFYPEER, 0L) != CURLE_OK) {
                fprintf(stderr, "[urlOpen] Couldn't set CURLOPT_SSL_VERIFYPEER\n");
                goto error;
            }
            if (curl_easy_setopt(URL->x.curl, CURLOPT_SSL_VERIFYHOST, 0L) != CURLE_OK) {
                fprintf(stderr, "[urlOpen] Couldn't set CURLOPT_SSL_VERIFYHOST\n");
                goto error;
            }
 
            if (callBack) {
                code = callBack(URL->x.curl);
                if (code != CURLE_OK) {
                    fprintf(stderr, "[urlOpen] user callback error: %s\n", curl_easy_strerror(code));
                    goto error;
                }
            }
 
            /* common Windows-friendly opts before first perform */
            bw_curl_apply_common_opts(URL->x.curl);
 
            code = curl_easy_perform(URL->x.curl);
            errno = 0;
            if (code != CURLE_OK) {
                fprintf(stderr, "[urlOpen] curl_easy_perform error: %s\n", curl_easy_strerror(code));
                goto error;
            }
#endif /* NOCURL */
        }
    } else {
        /* writing mode (local only) */
        URL->type = BWG_FILE;
        URL->x.fp = fopen(fname, mode);
        if (!(URL->x.fp)) {
            free(URL);
            fprintf(stderr, "[urlOpen] Couldn't open %s for writing\n", fname);
            return NULL;
        }
    }
 
    if (url) free(url);
    if (req) free(req);
    return URL;
 
#ifndef NOCURL
error:
    if (url) free(url);
    if (req) free(req);
    free(URL->memBuf);
    curl_easy_cleanup(URL->x.curl);
    free(URL);
    return NULL;
#endif
}
 
/* Free resources and cleanup curl */
void urlClose(URL_t *URL) {
    if (URL->type == BWG_FILE) {
        fclose(URL->x.fp);
#ifndef NOCURL
    } else {
        free(URL->memBuf);
        curl_easy_cleanup(URL->x.curl);
#endif
    }
    free(URL);
}

