#pragma once

#include <stdint.h>
#include "md_str.h"

struct md_allocator_i;
typedef uint64_t timestamp_t;

#ifdef __cplusplus
extern "C" {
#endif

/*
int64_t     md_os_physical_ram_in_bytes();
*/
str_t       md_os_current_working_directory();

str_t       md_os_path_make_canonical(str_t path, struct md_allocator_i* alloc);
str_t       md_os_path_make_relative(str_t path_from, str_t path_to, struct md_allocator_i* alloc);

timestamp_t md_os_time_current();
double      md_os_time_delta_in_ms(timestamp_t t0, timestamp_t t1);
double      md_os_time_delta_in_s (timestamp_t t0, timestamp_t t1);

void        md_os_sleep(int64_t milliseconds);

#ifdef __cplusplus
}
#endif
