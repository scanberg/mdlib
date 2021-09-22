#pragma once

#include <stdint.h>
#include "md_str.h"

typedef uint64_t timestamp_t;

int64_t     md_os_physical_ram_in_bytes();
str_t       md_os_current_working_directory();

timestamp_t md_os_time_current();
double      md_os_time_delta_in_ms(timestamp_t t0, timestamp_t t1);
double      md_os_time_delta_in_s (timestamp_t t0, timestamp_t t1);