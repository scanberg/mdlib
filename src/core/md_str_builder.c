#include <core/md_str_builder.h>

#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_array.h>

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

void md_strb_init(md_strb_t* sb, struct md_allocator_i* alloc) {
	ASSERT(sb);
	ASSERT(alloc);
#if DEBUG
	if (sb->buf != 0) {
		md_log(MD_LOG_TYPE_DEBUG, "Initializing non-zero struct, string builder may be leaking here.");
	}
#endif
	sb->buf = NULL;
	sb->alloc = alloc;
	md_array_ensure(sb->buf, 32, alloc);
}

void md_strb_free(md_strb_t* sb) {
	ASSERT(sb);
	if (!sb->buf) {
		return;
	}
	ASSERT(sb->alloc);
	md_array_free(sb->buf, sb->alloc);
	sb->buf = 0;
	sb->alloc = NULL;
}

void md_strb_push_cstr(md_strb_t* sb, const char* cstr) {
	md_strb_push_cstrl(sb, cstr, strlen(cstr));
}

void md_strb_push_cstrl(md_strb_t* sb, const char* cstr, size_t len) {
	ASSERT(sb);
	if (len > 0) {
		if (md_array_size(sb->buf)) md_array_pop(sb->buf); // Remove zero terminator
		md_array_push_array(sb->buf, cstr, len, sb->alloc);
		md_array_push(sb->buf, '\0', sb->alloc);
	}
}

void md_strb_fmt(md_strb_t* sb, const char* format, ...) {
	ASSERT(sb);
	va_list args;
	va_start(args, format);
	int len = vsnprintf(NULL, 0, format, args);
	va_end(args);

	if (len > 0) {
		if (md_array_size(sb->buf)) md_array_pop(sb->buf); // Remove zero terminator
		size_t offset = md_array_size(sb->buf);
		md_array_resize(sb->buf, md_array_size(sb->buf) + len + 1, sb->alloc);

		va_start(args, format);
		vsnprintf(sb->buf + offset, len + 1, format, args);
		va_end(args);
	}
}

void md_strb_push_str(md_strb_t* sb, str_t str) {
	ASSERT(sb);
	if (str.len > 0) {
		if (md_array_size(sb->buf)) md_array_pop(sb->buf); // Remove zero terminator
		md_array_push_array(sb->buf, str.ptr, str.len, sb->alloc);
		md_array_push(sb->buf, '\0', sb->alloc);
	}
}

void md_strb_pop(md_strb_t* sb, size_t n) {
	ASSERT(sb);
	ASSERT(n > 0);
	n = MIN(n, md_array_size(sb->buf));
	md_array_shrink(sb->buf, n);
	md_array_push(sb->buf, '\0', sb->alloc);
}

void md_strb_push_char(md_strb_t* sb, char c) {
    ASSERT(sb);
    if (md_array_size(sb->buf)) md_array_pop(sb->buf); // Remove zero terminator
    md_array_push(sb->buf, c, sb->alloc);
    md_array_push(sb->buf, '\0', sb->alloc);
}

void md_strb_reset(md_strb_t* sb) {
	ASSERT(sb);
	md_array_shrink(sb->buf, 0);
}

const char* md_strb_to_cstr(const md_strb_t* sb) {
	return sb->buf;
}

size_t md_strb_size(const md_strb_t* sb) {
	return md_array_size(sb->buf);
}

size_t md_strb_len(const md_strb_t* sb) {
    const size_t size = md_array_size(sb->buf);
    return size ? size - 1 : 0;
}

str_t md_strb_to_str(const md_strb_t* sb) {
	ASSERT(sb);
	const size_t len = md_strb_len(sb);
	return (str_t) { len ? sb->buf : NULL, len };
}

