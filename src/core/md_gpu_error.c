/* md_gpu_error.c -- see md_gpu_error.h. */

#include "md_gpu_error.h"
#include "md_gpu.h"

#include <core/md_format.h>
#include <core/md_log.h>

#include <stdarg.h>

#define MD_GPU_ERROR_CAP 2560u   /* room for an adapter list with every missing feature */

static MD_THREAD_LOCAL char md_gpu_error_buf[MD_GPU_ERROR_CAP];
static MD_THREAD_LOCAL bool md_gpu_has_error;

bool md_gpu_fail(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    md_vsnprintf(md_gpu_error_buf, sizeof(md_gpu_error_buf), fmt, ap);
    va_end(ap);
    md_gpu_has_error = true;
    MD_LOG_ERROR("md_gpu: %s", md_gpu_error_buf);
    return false;
}

void md_gpu_clear_error(void) {
    md_gpu_has_error = false;
}

const char* md_gpu_last_error(void) {
    return md_gpu_has_error ? md_gpu_error_buf : NULL;
}
