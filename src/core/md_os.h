#pragma once

#include <stdint.h>
#include "md_str.h"

struct  md_allocator_i;
typedef int64_t timestamp_t;

#ifdef __cplusplus
extern "C" {
#endif

/*
int64_t     md_os_physical_ram_in_bytes();
*/
str_t       md_os_current_working_directory();

str_t       md_os_path_make_canonical(str_t path, struct md_allocator_i* alloc);
str_t       md_os_path_make_relative(str_t path_from, str_t path_to, struct md_allocator_i* alloc);

bool		md_os_path_exists(str_t path);
bool		md_os_path_is_directory(str_t path);

timestamp_t md_os_time_current();
timestamp_t md_os_time_from_milliseconds(int64_t milliseconds);
double      md_os_time_as_milliseconds(timestamp_t t);
double      md_os_time_as_seconds(timestamp_t t);

void        md_os_sleep(int64_t milliseconds);

#ifdef __cplusplus
}
#endif
