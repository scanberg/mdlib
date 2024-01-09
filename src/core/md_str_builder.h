#pragma once

#include <core/md_str.h>
#include <core/md_array.h>

// String Builder API

struct md_allocator_i;

typedef struct md_strb_t {
    md_array(char) buf;
    struct md_allocator_i* alloc;

#ifdef __cplusplus
	md_strb_t& operator += (str_t str);
	md_strb_t& operator += (const char* cstr);
	md_strb_t& operator += (int c);
	operator str_t() const;
#endif

} md_strb_t;

#ifdef __cplusplus
extern "C" {
#endif

void md_strb_init(md_strb_t* sb, struct md_allocator_i* alloc);
void md_strb_free(md_strb_t* sb);

static inline md_strb_t md_strb_create(struct md_allocator_i* alloc) {
	ASSERT(alloc);
	md_strb_t sb = {0};
	md_strb_init(&sb, alloc);
	return sb;
}

void md_strb_fmt(md_strb_t* sb, const char* format, ...);
void md_strb_push_char(md_strb_t* sb, char c);
void md_strb_push_cstr(md_strb_t* sb, const char* cstr);
void md_strb_push_cstrl(md_strb_t* sb, const char* cstr, size_t len);
void md_strb_push_str(md_strb_t* sb, str_t str);

void md_strb_reset(md_strb_t* sb);
void md_strb_pop(md_strb_t* sb, size_t n);

// Pointer to string buffer
char* md_strb_ptr(md_strb_t sb);

// Length of string in bytes (excluding zero terminator)
size_t md_strb_len(md_strb_t sb);

// Capacity of buffer in bytes
size_t md_strb_cap(md_strb_t sb);

// Size of buffer in bytes (including zero terminator)
size_t md_strb_size(md_strb_t sb);

const char* md_strb_to_cstr(md_strb_t sb);
str_t       md_strb_to_str (md_strb_t sb);

#ifdef __cplusplus

inline md_strb_t& md_strb_t::operator += (str_t str) {
	md_strb_push_str(this, str);
	return *this;
}

inline md_strb_t& md_strb_t::operator += (const char* cstr) {
	md_strb_push_cstr(this, cstr);
	return *this;
}

inline md_strb_t& md_strb_t::operator += (int c) {
	md_strb_push_char(this, (char)c);
	return *this;
}

inline md_strb_t::operator str_t() const {
	return md_strb_to_str(*this);
}

}
#endif
