#pragma once

#include <core/md_intrinsics.h>
#include <core/md_str.h>
#include <core/md_os.h>

#include <string.h>

/*
// This is a worse implementation of memchr

static inline uint32_t tzcnt_u32(uint32_t x) {
#ifdef _MSC_VER
    return _tzcnt_u32(x);
#else
    return __builtin_ctz(x);
#endif
}

static inline __m128i _load(const char* p) {
    return _mm_loadu_si128((const __m128i*)p);
}

static inline __m128i _load_partial(const char* p, uint32_t n) {
    __m128i a0 = _mm_set1_epi8((char)n);
    __m128i a1 = _mm_set_epi8(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
    __m128i m  = _mm_cmpgt_epi8(a0,a1);
    return _mm_and_si128(_load(p), m);
}

static inline int64_t find_char(const char* ptr, int64_t len, char character) {
    int64_t off = 0;
    __m128i c = _mm_set1_epi8(character);
    while (len > 16) {
        __m128i d = _load(ptr + off);
        __m128i m = _mm_cmpeq_epi8(d, c);
        int r = _mm_movemask_epi8(m);
        if (r)
            return off + tzcnt_u32(r);
        off += 16;
        len -= 16;
    }
    if (len > 0) {
        __m128i d = _load_partial(ptr + off, (uint32_t)len);
        __m128i m = _mm_cmpeq_epi8(d, c);
        int r = _mm_movemask_epi8(m);
        if (r)
            return off + tzcnt_u32(r);
    }
    return -1;
}
*/

typedef struct md_buffered_reader_t {
    str_t   str;
    char*   buf;
    int64_t cap;
    md_file_o* file;
} md_buffered_reader_t;

static inline md_buffered_reader_t md_buffered_reader_from_file(char* buf, int64_t cap, md_file_o* file) {
    ASSERT(buf);
    ASSERT(cap > 0);
    ASSERT(file);
    
    md_buffered_reader_t lr = {
        .str = {0, 0},
        .buf = buf,
        .cap = cap,
        .file = file,
    };

    return lr;
}

static inline md_buffered_reader_t md_buffered_line_reader_from_str(str_t str) {
    md_buffered_reader_t reader = {
        .str = str,
    };
    return reader;
}

static inline bool md_buffered_line_reader_extract_line(str_t* line, md_buffered_reader_t* lr) {
    ASSERT(lr);
    ASSERT(line);
    if (lr->file && !lr->str.len) {
        ASSERT(lr->buf);
        ASSERT(lr->cap > 0);
        const int64_t bytes_read = md_file_read_lines(lr->file, lr->buf, lr->cap);
        if (bytes_read > 0) {
            lr->str.ptr = lr->buf;
            lr->str.len = bytes_read;
        }
    }
    return str_extract_line(line, &lr->str); 
}

static inline double parse_float_dec(const char* ptr, int64_t len) {
    static const double pow10[22] = {
        1.0 / 1e+0,  1.0 / 1e+1,  1.0 / 1e+2,  1.0 / 1e+3,  1.0 / 1e+4,  1.0 / 1e+5,  1.0 / 1e+6,  1.0 / 1e+7,
        1.0 / 1e+8,  1.0 / 1e+9,  1.0 / 1e+10, 1.0 / 1e+11, 1.0 / 1e+12, 1.0 / 1e+13, 1.0 / 1e+14, 1.0 / 1e+15,
        1.0 / 1e+16, 1.0 / 1e+17, 1.0 / 1e+18, 1.0 / 1e+19, 1.0 / 1e+20, 1.0 / 1e+21};

    const char* c = ptr;
    const char* end = ptr + len;

    bool neg = (*c == '-');
	if (*c == '-' || *c == '+') ++c;

    double val = 0.0;
	while (c < end && '0' <= *c && *c <= '9') {
		val = val * 10 + (*c - '0');
		++c;
	}

    if (c < end && *c == '.') {
        const char* dec = ++c;
        while (c < end && '0' <= *c && *c <= '9') {
            val = val * 10 + ((int)(*c) - '0');
            ++c;
        }
        val *= pow10[c - dec];
    }

    return neg ? -val : val;
}

static inline double parse_float(str_t str) {
    ASSERT(str.ptr);
    static const double pow10[32] = {
        1e+0,  1e+1,  1e+2,  1e+3,  1e+4,  1e+5,  1e+6,  1e+7,
        1e+8,  1e+9,  1e+10, 1e+11, 1e+12, 1e+13, 1e+14, 1e+15,
        1e+16, 1e+17, 1e+18, 1e+19, 1e+20, 1e+21, 1e+22, 1e+23,
        1e+24, 1e+25, 1e+26, 1e+27, 1e+28, 1e+29, 1e+30, 1e+31,
    };

    const char* c = str.ptr;
    const char* end = str.ptr + str.len;

    while (c < end && is_whitespace(*c)) ++c;

    double val = 0;
    double sign = 1;
    if (*c == '-') {
        ++c;
        sign = -1;
    }
    while (c < end && is_digit(*c)) {
        val = val * 10 + ((int)(*c) - '0');
        ++c;
    }

    if (c < end && *c == '.') {
        const char* dec = ++c;
        while (c < end && (c - dec) < (int64_t)ARRAY_SIZE(pow10) - 1 && is_digit(*c)) {
            val = val * 10 + ((int)(*c) - '0');
            ++c;
        }
        int64_t count = c - dec;
        val /= pow10[count];
    }

    if (c < end && (*c == 'e' || *c == 'E')) {
        ++c;
        int exp_sign = 1;
        if (c < end && (*c == '+' || *c == '-')) {
            exp_sign = *c == '+' ? 1 : -1;
            ++c;
        }
        int exp_val = 0;
        while (c < end && is_digit(*c)) {
            exp_val = exp_val * 10 + ((int)(*c) - '0');
            ++c;
        }
        while (exp_val) {
            int ev = MIN(exp_val, (int)ARRAY_SIZE(pow10)-1);
            val = exp_sign > 0 ? val * pow10[ev] : val / pow10[ev];
            exp_val -= ev;
        }
    }

    return sign * val;
}

static inline int64_t parse_int(str_t str) {
    if (str.len <= 0) return 0;

    int64_t val = 0;
    const char* c   = str.ptr;
    const char* end = str.ptr + str.len;
    bool neg = (*c == '-');

    if (neg)
        ++c;

    while (c < end && is_digit(*c)) {
        val = val * 10 + ((int64_t)(*c) - (int64_t)'0');
        ++c;
    }
    return neg ? -val : val;
}

static inline __m128i mask(int64_t n) {
    //ALIGNAS(32) static const char mask[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    //return _mm_loadu_si128((const __m128i*)(mask + n));
    return _mm_cmpgt_epi8(_mm_set1_epi8((char)n), _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
}

static inline int find_first_char(__m128i v, char c) {
    __m128i m = _mm_cmpeq_epi8(v, _mm_set1_epi8(c));
    int r = ctz32(_mm_movemask_epi8(m));
    return MIN(r, 16);
}

static inline uint64_t parse_uint64_si128(__m128i v, __m128i m) {
    v = _mm_subs_epu8(v, _mm_set1_epi8('0'));
    v = _mm_and_si128(v, m);

    v = _mm_maddubs_epi16(v, _mm_set_epi8(1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10));
    v = _mm_madd_epi16(v, _mm_set_epi16(1, 100, 1, 100, 1, 100, 1, 100));
    v = _mm_packus_epi32(v, v);
    v = _mm_madd_epi16(v, _mm_set_epi16(1, 10000, 1, 10000, 1, 10000, 1, 10000));

    return (uint64_t)_mm_cvtsi128_si32(v) * 100000000 + (int64_t)_mm_extract_epi32(v, 1);
}

static inline uint64_t parse_uint64_simd(const char* ptr, int64_t len) {
    __m128i m = mask(len);
    __m128i v = _mm_lddqu_si128((const __m128i*)(ptr+len-16));
    return parse_uint64_si128(v, m);
}

static inline int64_t parse_int_simd(const char* ptr, int64_t len) {
    len = MIN(len, 16);
	bool neg = (*ptr == '-');
    if (neg) {
        ++ptr;
        --len;
    }
    int64_t val = parse_uint64_simd(ptr, len);
    return neg ? -val : val;
}

static inline int64_t str_parse_int_simd(str_t s) {
    return parse_int_simd(s.ptr, s.len);
}

static inline double parse_float_simd(const char* ptr, int64_t len) {
    static const double pow10[16] = { 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16 };
    len = MIN(len, 16);
    bool neg = (*ptr == '-');
    __m128i v = _mm_lddqu_si128((const __m128i*)ptr);
    int dec = find_first_char(v, '.');
	double whole = (double)parse_uint64_si128(_mm_lddqu_si128((const __m128i*)(ptr + dec - 16)), mask(dec));
    double fract = 0.0;
    if (dec < 16) {
        for (int64_t i = dec + 1; i < len; ++i) {
            fract = fract * 10 + ((int)(ptr[i]) - '0');
        }
        fract /= pow10[len-dec];
    }

    //uint64_t fract = parse_uint64_si128(_mm_lddqu_si128((const __m128i*)(ptr + len - 16)), mask(len - (dec + 1)));

    return (double)whole + fract;
}

static inline double str_parse_float_simd(str_t str) {
    return parse_float_simd(str.ptr, str.len);
}

// Extracts token with whitespace as delimiter
static inline bool extract_token(str_t* tok, str_t* str) {
    ASSERT(tok);
    ASSERT(str);
    if (str_empty(*str)) return false;

    const char* end = str->ptr + str->len;
    const char* c = str->ptr;
    while (c < end && is_whitespace(*c)) ++c;
    if (c >= end) return false;

    const char* tok_beg = c;
    while (c < end && !is_whitespace(*c)) ++c;

    tok->ptr = tok_beg;
    tok->len = c - tok_beg;

    str->ptr = c < end ? c + 1 : end;
    str->len = end - str->ptr;

    return true;
}


// Extract multiple tokens with whitespace as delimiter
static inline int64_t extract_tokens(str_t token_arr[], int64_t token_cap, str_t* str) {
    ASSERT(token_arr);
    ASSERT(token_cap >= 0);
    ASSERT(str);

    int64_t num_tokens = 0;
    while (num_tokens < token_cap && extract_token(&token_arr[num_tokens], str)) {
        num_tokens += 1;
    }
    return num_tokens;
}

// Extracts token with specific delimiter
static inline bool extract_token_delim(str_t* tok, str_t* str, char delim) {
    ASSERT(tok);
    ASSERT(str);
    if (!str->ptr || str->len == 0) return false;

    const char* beg = str->ptr;
    const char* end = str->ptr + str->len;
    const char* c = str->ptr;
    while (c != end && *c != delim) {
        ++c;
    }
    tok->ptr = beg;
    tok->len = c - beg;

    str->ptr = c != end ? c + 1 : end;
    str->len = end - str->ptr;

    return true;
}