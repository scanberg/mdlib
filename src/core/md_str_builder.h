#pragma once

#include <core/md_str.h>

// String Builder API

struct md_allocator_i;

typedef struct md_strb_t {
    char* buf;
    struct md_allocator_i* alloc;

#ifdef __cplusplus
	md_strb_t& operator += (str_t str);
	md_strb_t& operator += (const char* cstr);
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
void md_strb_char(md_strb_t* sb, char c);
void md_strb_cstr(md_strb_t* sb, const char* cstr);
void md_strb_cstrl(md_strb_t* sb, const char* cstr, int64_t len);
void md_strb_str(md_strb_t* sb, str_t str);
void md_strb_reset(md_strb_t* sb);

void md_strb_pop(md_strb_t* sb, int64_t n);

int64_t     md_strb_len(const md_strb_t* sb);
const char* md_strb_to_cstr(const md_strb_t* sb);
str_t       md_strb_to_str(const md_strb_t* sb);

#ifdef __cplusplus

inline md_strb_t& md_strb_t::operator += (str_t str) {
	md_strb_str(this, str);
	return *this;
}

inline md_strb_t& md_strb_t::operator += (const char* cstr) {
	md_strb_cstr(this, cstr);
	return *this;
}

inline md_strb_t::operator str_t() const {
	return md_strb_to_str(this);
}

}
#endif
