/*
md_gpu_error.h -- internal to md_gpu: how its parts report a failure.

md_gpu_fail formats a message, keeps it for md_gpu_last_error (per thread),
logs it through md_log and returns false, so that a failing path can end in
`return md_gpu_fail(...)`. Shared by the backends and md_gpu_heap.
*/

#ifndef MD_GPU_ERROR_H
#define MD_GPU_ERROR_H

#include <core/md_common.h>

#include <stdbool.h>

bool md_gpu_fail(const char* fmt, ...) MD_PRINTF_FORMAT(1, 2);

/* Forget this thread's last error, at the start of an operation whose
   failure is reported through md_gpu_last_error (device creation). */
void md_gpu_clear_error(void);

#endif
