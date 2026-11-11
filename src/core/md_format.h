#pragma once

#include <core/md_common.h>

#include <stddef.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

// snprintf and vsnprintf as C99 has them, in the "C" locale whatever the process or the calling thread
// has set: the decimal point is always '.'. The process locale is not ours to rely on. A GTK file dialog
// switches it to the user's, and so may a shell extension a Windows file dialog loads, and under a decimal
// comma locale the C library writes 1.5 as "1,5" - which md_parse, as it does not follow the locale either,
// reads back as 1 and 5. Everything mdlib formats goes through these (md_alloc_printf, the writers,
// md_file_printf and the log too), so what it writes is the same text in every locale.
//
// Implemented in md_os.c: each call is handed a "C" locale of its own through the platform's API, rather
// than setlocale, which is process wide and would change what every other thread prints meanwhile.
int md_snprintf (char* buf, size_t cap, const char* format, ...) MD_PRINTF_FORMAT(3, 4);
int md_vsnprintf(char* buf, size_t cap, const char* format, va_list args);

#ifdef __cplusplus
}
#endif
