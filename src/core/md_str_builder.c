#include "md_str_builder.h"

#include "md_allocator.h"
#include "md_log.h"
#include "md_array.h"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

void md_str_builder_init(md_str_builder_t* sb, struct md_allocator_i* alloc) {
	ASSERT(sb);
	ASSERT(alloc);
#if DEBUG
	if (sb->buf != 0) {
		md_print(MD_LOG_TYPE_DEBUG, "Initializing non-zero struct, string builder may be leaking here.");
	}
#endif
	sb->buf = NULL;
	sb->alloc = alloc;
	md_array_ensure(sb->buf, 32, alloc);
}

void md_str_builder_free(md_str_builder_t* sb) {
	ASSERT(sb);
	if (!sb->buf) {
		return;
	}
	ASSERT(sb->alloc);
	md_array_free(sb->buf, sb->alloc);
	sb->buf = 0;
	sb->alloc = NULL;
}

void md_str_builder_append_cstr(md_str_builder_t* sb, const char* cstr) {
	md_str_builder_append_cstr_len(sb, cstr, (int64_t)strlen(cstr));
}

void md_str_builder_append_cstr_len(md_str_builder_t* sb, const char* cstr, int64_t len) {
	ASSERT(sb);
	if (len > 0) {
		if (md_array_size(sb->buf)) md_array_pop(sb->buf); // Remove zero terminator
		md_array_push_array(sb->buf, cstr, len, sb->alloc);
		md_array_push(sb->buf, '\0', sb->alloc);
	}
}

void md_str_builder_printf(md_str_builder_t* sb, const char* format, ...) {
	ASSERT(sb);
	va_list args;
	va_start(args, format);
	int len = vsnprintf(NULL, 0, format, args);
	va_end(args);

	if (len > 0) {
		if (md_array_size(sb->buf)) md_array_pop(sb->buf); // Remove zero terminator
		int64_t offset = md_array_size(sb->buf);
		md_array_resize(sb->buf, md_array_size(sb->buf) + len + 1, sb->alloc);

		va_start(args, format);
		vsnprintf(sb->buf + offset, len + 1, format, args);
		va_end(args);
	}
}

void md_str_builder_append_str(md_str_builder_t* sb, str_t str) {
	ASSERT(sb);
	if (str.len > 0) {
		if (md_array_size(sb->buf)) md_array_pop(sb->buf); // Remove zero terminator
		md_array_push_array(sb->buf, str.ptr, str.len, sb->alloc);
		md_array_push(sb->buf, '\0', sb->alloc);
	}
}

void md_str_builder_pop(md_str_builder_t* sb, int64_t n) {
	ASSERT(sb);
	ASSERT(n > 0);
	md_array_shrink(sb->buf, MAX(0, md_array_size(sb->buf) - n));
	md_array_push(sb->buf, '\0', sb->alloc);
}

void md_str_builder_append_char(md_str_builder_t* sb, char c) {
    ASSERT(sb);
    if (md_array_size(sb->buf)) md_array_pop(sb->buf); // Remove zero terminator
    md_array_push(sb->buf, c, sb->alloc);
    md_array_push(sb->buf, '\0', sb->alloc);
}

void md_str_builder_reset(md_str_builder_t* sb) {
	ASSERT(sb);
	md_array_shrink(sb->buf, 0);
}

const char* md_str_builder_cstr(const md_str_builder_t* sb) {
	return sb->buf;
}

int64_t md_str_builder_len(const md_str_builder_t* sb) {
	return md_array_size(sb->buf);
}

str_t md_str_builder_to_str(const md_str_builder_t* sb) {
	ASSERT(sb);
	return (str_t) { md_array_size(sb->buf) ? sb->buf : 0, md_array_size(sb->buf) };
}

