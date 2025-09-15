#pragma once

#include <core/md_intrinsics.h>
#include <core/md_simd.h>
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

static inline md_128i _load(const char* p) {
    return _mm_loadu_si128((const md_128i*)p);
}

static inline md_128i _load_partial(const char* p, uint32_t n) {
    md_128i a0 = _mm_set1_epi8((char)n);
    md_128i a1 = _mm_set_epi8(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
    md_128i m  = _mm_cmpgt_epi8(a0,a1);
    return _mm_and_si128(_load(p), m);
}

static inline int64_t find_char(const char* ptr, int64_t len, char character) {
    int64_t off = 0;
    md_128i c = _mm_set1_epi8(character);
    while (len > 16) {
        md_128i d = _load(ptr + off);
        md_128i m = _mm_cmpeq_epi8(d, c);
        int r = _mm_movemask_epi8(m);
        if (r)
            return off + tzcnt_u32(r);
        off += 16;
        len -= 16;
    }
    if (len > 0) {
        md_128i d = _load_partial(ptr + off, (uint32_t)len);
        md_128i m = _mm_cmpeq_epi8(d, c);
        int r = _mm_movemask_epi8(m);
        if (r)
            return off + tzcnt_u32(r);
    }
    return -1;
}
*/

// Attempts to read complete lines from a file (As many that can fit into buf with given capacity)
static inline size_t md_parse_read_lines(md_file_o* file, char* buf, size_t cap) {
    if (!file || !buf || cap < 1) return 0;
    size_t len = md_file_read(file, buf, cap);
    if (len == cap) {
        size_t loc;
        const str_t str = {buf, len};
        if (str_rfind_char(&loc, str, '\n')) {
            const int64_t offset = (int64_t)loc + 1 - (int64_t)len;
            // Set file pointer to the beginning of the next line
            md_file_seek(file, offset, MD_FILE_CUR);
            len = loc + 1;
        }
    }
    return len;
}

typedef struct md_buffered_reader_t {
    str_t   str;
    char*   buf;
    size_t  cap;
    md_file_o* file;
} md_buffered_reader_t;

static inline md_buffered_reader_t md_buffered_reader_from_file(char* buf, size_t cap, md_file_o* file) {
    ASSERT(buf);
    ASSERT(file);
    
    md_buffered_reader_t lr = {
        .str = {0, 0},
        .buf = buf,
        .cap = cap,
        .file = file,
    };

    return lr;
}

static inline md_buffered_reader_t md_buffered_reader_from_str(str_t str) {
    md_buffered_reader_t reader = {
        .str = str,
        .cap = str.len,
    };
    return reader;
}

static inline void md_buffered_reader_ensure_lines(md_buffered_reader_t* r) {
    ASSERT(r);
    if (r->file && !r->str.len) {
        ASSERT(r->buf);
        const size_t bytes_read = md_parse_read_lines(r->file, r->buf, r->cap);
        if (bytes_read > 0) {
            r->str.ptr = r->buf;
            r->str.len = bytes_read;
        }
    }
}

static inline bool md_buffered_reader_extract_line(str_t* line, md_buffered_reader_t* r) {
    ASSERT(r);
    ASSERT(line);
    md_buffered_reader_ensure_lines(r);
    return str_extract_line(line, &r->str);
}

static inline bool md_buffered_reader_peek_line(str_t* line, md_buffered_reader_t* r) {
    ASSERT(r);
    ASSERT(line);
    md_buffered_reader_ensure_lines(r);
    return str_peek_line(line, &r->str);
}

static inline bool md_buffered_reader_skip_line(md_buffered_reader_t* r) {
    ASSERT(r);
    md_buffered_reader_ensure_lines(r);
    return str_skip_line(&r->str);
}

static inline void md_buffered_reader_reset(md_buffered_reader_t* r) {
    ASSERT(r);
    if (r->file) {
        md_file_seek(r->file, 0, MD_FILE_BEG);
        r->str.ptr = NULL;
        r->str.len = 0;
    } else {
        ASSERT(r->str.ptr && "Cannot reset uninitialized reader");
        r->str.ptr += r->str.len - r->cap;
        r->str.len = r->cap;
    }
}

/*
static inline bool md_buffered_reader_eof(const md_buffered_reader_t* r) {
    if (r->file) {
        return md_file_eof(r->file) && r->str.len == 0;
    } else {
        return r->str.len == 0;
    }
}
*/

#if 0
static inline void md_buffered_reader_seekg(md_buffered_reader_t* r, int64_t pos) {
    if (r->file) {
        return md_file_seek(r->file, pos, MD_FILE_BEG);
        r->str.ptr = NULL;
        r->str.len = 0;
    } else {
        // This is not correct, simply copied from reset
        ASSERT(r->str.ptr && "Cannot seekg uninitialized reader");
        r->str.ptr += r->str.len - r->cap;
        r->str.len = r->cap;
    }
}
#endif

// Gets the current position within the buffer
static inline int64_t md_buffered_reader_tellg(const md_buffered_reader_t* r) {
    if (r->file) {
        return md_file_tell(r->file) - r->str.len;
    } else {
        return r->cap - r->str.len;
    }
}

// Check if token is a valid float
// Make sure str is trimmed from whitespace
static inline bool is_float(str_t str) {
    const char* c = str.ptr;
    const char* end = str.ptr + str.len;

    if (c >= end) return false;
    if (*c == '-' || *c == '+') ++c;
    while (c < end && is_digit(*c)) ++c;
    if (*c == '.') {
        ++c;
        while (c < end && is_digit(*c)) ++c;
    }
    if (*c == 'e' || *c == 'E') {
        ++c;
        if (*c == '-' || *c == '+') ++c;
        while (c < end && is_digit(*c)) ++c;
    }
    return c == end;
}

// Check if token is a valid decimal float
// Make sure str is trimmed from whitespace
static inline bool is_float_dec(str_t str) {
    const char* c = str.ptr;
    const char* end = str.ptr + str.len;

    if (c >= end) return false;
    if (*c == '-' || *c == '+') ++c;
    while (c < end && is_digit(*c)) ++c;
    if (*c == '.') {
        ++c;
        while (c < end && is_digit(*c)) ++c;
    }
    return c == end;
}

// Check if token is a valid integer
// Make sure str is trimmed from whitespace
static inline bool is_int(str_t str) {
    const char* c = str.ptr;
    const char* end = str.ptr + str.len;

    if (c >= end) return false;
    if (*c == '-' || *c == '+') ++c;
    while (c < end && is_digit(*c)) ++c;
    return c == end;
}

static const double pow10_table[32] = {
    1e+0,  1e+1,  1e+2,  1e+3,  1e+4,  1e+5,  1e+6,  1e+7,
    1e+8,  1e+9,  1e+10, 1e+11, 1e+12, 1e+13, 1e+14, 1e+15,
    1e+16, 1e+17, 1e+18, 1e+19, 1e+20, 1e+21, 1e+22, 1e+23,
    1e+24, 1e+25, 1e+26, 1e+27, 1e+28, 1e+29, 1e+30, 1e+31,
};

// Parsing routine for floating point numbers without scientific notation (exponent)
static inline double parse_float_dec(const char* ptr, int64_t len) {

    const char* c = ptr;
    const char* end = ptr + len;

    bool neg = (*c == '-');
    if (*c == '-' || *c == '+') ++c;

    double val = 0.0;
    while (c < end && ('0' <= *c && *c <= '9')) {
        val = val * 10 + (*c - '0');
        ++c;
    }

    if (c < end && *c == '.') {
        const char* dec = ++c;
        while (c < end && '0' <= *c && *c <= '9') {
            val = val * 10 + ((int)(*c) - '0');
            ++c;
        }
        val *= pow10_table[c - dec];
    }

    return neg ? -val : val;
}

// It is not exact beyond a certain number of decimal digits, but it is fast and will suffice
// for most text based formats we deal with
static inline double parse_float(str_t str) {
    ASSERT(str.ptr);

    const char* c = str.ptr;
    const char* end = str.ptr + str.len;

    while (c < end && is_whitespace(*c)) ++c;

    bool negate = (*c == '-');
    uint64_t whole = 0;
    uint64_t fract = 0;
    uint64_t fract_len = 0;
    
    if (*c == '-')
        ++c;
    
    while (c < end && is_digit(*c)) {
        whole = whole * 10 + ((int)(*c) - '0');
        ++c;
    }

    if (c < end && *c == '.') {
        ++c;
        const char* fbeg = c;
        const char* fend = MIN(end, c + 19); // Will overflow if we go beyond this
        while (c < fend && is_digit(*c)) {
            fract = fract * 10 + ((int)(*c) - '0');
            ++c;
        }
        fract_len = c - fbeg;
    }

    double val = (double)whole + (double)fract / pow10_table[fract_len];

    if (c < end && (*c == 'e' || *c == 'E')) {
        ++c;
        bool exp_neg = *c == '-';
        if (c < end && (*c == '+' || *c == '-')) {
            ++c;
        }
        int exp_val = 0;
        while (c < end && is_digit(*c)) {
            exp_val = exp_val * 10 + ((int)(*c) - '0');
            ++c;
        }
        while (exp_val) {
            int ev = MIN(exp_val, (int)ARRAY_SIZE(pow10_table)-1);
            val = exp_neg ? val / pow10_table[ev] : val * pow10_table[ev];
            exp_val -= ev;
        }
    }

    return negate ? -val : val;
}

static inline int64_t parse_int(str_t str) {
    if (str.len <= 0) return 0;

    int64_t val = 0;
    const char* c   = str.ptr;
    const char* end = str.ptr + str.len;

    // Skip space
    while (c < end && is_whitespace(*c)) ++c;

    bool neg = false;
    if (*c == '-') {
		neg = true;
		++c;
	}

    while (c < end && is_digit(*c)) {
        val = val * 10 + ((int64_t)(*c) - (int64_t)'0');
        ++c;
    }
    return neg ? -val : val;
}

#if 0
static inline md_128i mask(size_t n) {
    return md_mm_cmpgt_epi8(md_mm_set1_epi8((char)n), md_mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
}

static inline uint64_t
#if MD_COMPILER_GCC || MD_COMPILER_CLANG
__attribute__((target("sse4.1")))
#endif
parse_uint64_si128(md_128i v, md_128i m) {
    v = _mm_subs_epu8(v, _mm_set1_epi8('0'));
    v = _mm_and_si128(v, m);

    v = _mm_maddubs_epi16(v, _mm_set_epi8(1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10));
    v = _mm_madd_epi16(v, _mm_set_epi16(1, 100, 1, 100, 1, 100, 1, 100));
    v = _mm_packus_epi32(v, v);
    v = _mm_madd_epi16(v, _mm_set_epi16(1, 10000, 1, 10000, 1, 10000, 1, 10000));

    return (uint64_t)_mm_cvtsi128_si32(v) * 100000000 + (int64_t)_mm_extract_epi32(v, 1);
}

static inline uint64_t
#if MD_COMPILER_GCC || MD_COMPILER_CLANG
__attribute__((target("sse4.1")))
#endif
parse_uint32x2_si128(md_128i v, md_128i m) {
    v = _mm_subs_epu8(v, md_mm_set1_epi8('0'));
    v = md_mm_and_si128(v, m);
    v = _mm_maddubs_epi16(v, md_mm_set_epi8(1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10));
    v = md_mm_madd_epi16(v, md_mm_set_epi16(1, 100, 1, 100, 1, 100, 1, 100));
    v = _mm_packus_epi32(v, v);
    v = md_mm_madd_epi16(v, _mm_set_epi16(1, 10000, 1, 10000, 1, 10000, 1, 10000));
    return md_mm_cvtsi128_si64(v);
}

static inline uint64_t parse_u64(const char* ptr, int64_t len) {
    md_128i m = mask(len);
    md_128i v = _mm_loadu_si128((const md_128i*)(ptr+len-16));
    return parse_uint64_si128(v, m);
}
#endif

static inline uint64_t load_u64(const char* ptr) {
    uint64_t val;
    MEMCPY(&val, ptr, sizeof(val));
    return val;
}

static inline uint32_t parse_u32(const char* ptr, size_t len) {
    const uint64_t mask = 0x000000FF000000FF;
    const uint64_t mul1 = 0x000F424000000064; // 100 + (1000000ULL << 32)
    const uint64_t mul2 = 0x0000271000000001; // 1 + (10000ULL << 32)
    uint64_t shift = (64 - (len << 3));
    uint64_t val = load_u64(ptr);
    val = val << shift;
    val -= 0x3030303030303030 << shift; // '0'
    val = (val * 10) + (val >> 8); // val = (val * 2561) >> 8;
    val = (((val & mask) * mul1) + (((val >> 16) & mask) * mul2)) >> 32;
    return (uint32_t)val;
}

#if 0
static inline uint64_t parse_uint_wide(const char* ptr, size_t len) {
    uint64_t val;
    if (len > 8) {
        return parse_u64(ptr, len);
    } else {
        return parse_u32(ptr, len);
    }
}

static inline int64_t parse_int_wide(const char* ptr, size_t len) {
	bool neg = (*ptr == '-');
    if (*ptr == '-') {
        ++ptr;
        --len;
    }

    int64_t val = parse_uint_wide(ptr, len);
    return neg ? -val : val;
}

static inline md_128i shift_left(md_128i in, int n) {
    static const char shift_shuffle_lookup[32] = {
        -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128,
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
    };

    const char*   lookup_center = shift_shuffle_lookup + 16;
    const md_128i shuffle = md_mm_loadu_si128(lookup_center - n);
    return md_mm_shuffle_epi8(in, shuffle);
}
#endif

// https://github.com/fastfloat/fast_float/blob/main/include/fast_float/ascii_number.h
static inline uint32_t parse_eight_digits_unrolled(uint64_t val) {
    const uint64_t mask = 0x000000FF000000FF;
    const uint64_t mul1 = 0x000F424000000064; // 100 + (1000000ULL << 32)
    const uint64_t mul2 = 0x0000271000000001; // 1 + (10000ULL << 32)
    val -= 0x3030303030303030; // '0'
    val = (val * 10) + (val >> 8); // val = (val * 2561) >> 8;
    val = (((val & mask) * mul1) + (((val >> 16) & mask) * mul2)) >> 32;
    return (uint32_t)val;
}

static inline int find_dot(const char* ptr) {
    uint64_t val = load_u64(ptr);
    uint64_t msk = val & 0x2E2E2E2E2E2E2E2E; // '.'
    return msk ? (int)clz64(msk) >> 1 : 0;
}

struct float_token_t {
    const char* ptr;
    uint8_t whole_beg;
    uint8_t whole_end;
    uint8_t fract_beg;
    uint8_t fract_end;
    uint8_t exp_beg;
    uint8_t exp_end;
    uint8_t flags;
    uint8_t _unused;
};

static inline int find_first_char(md_128i v, char c) {
    md_128i m = md_mm_cmpeq_epi8(v, md_mm_set1_epi8(c));
    return ctz32(md_mm_movemask_epi8(m));
}

// Specialized version where the the integer and fractional part each
// are expected to fit into 8 characters. ptr is expected to have 16 readable characters from its loacation in memory.
static inline double parse_float_wide(const char* ptr, size_t len) {
    int neg = (*ptr == '-');
    if (*ptr == '-') {
        ++ptr;
        --len;
    }
    // Find decimal point by 128-bit intrinsics
    int dec = find_first_char(md_mm_loadu_si128((const md_128i*)ptr), '.');
    dec = MIN(dec, (int)len);
    int frac_len = (int)len - (dec + 1);
    frac_len = MAX(0, frac_len);

	uint32_t whole = parse_u32(ptr, dec);
    uint32_t fract = parse_u32(ptr + dec + 1, frac_len);

    double result = whole + fract / pow10_table[frac_len];
    return neg ? -result : result;
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
static inline size_t extract_tokens(str_t tok_arr[], size_t tok_cap, str_t* str) {
    ASSERT(tok_arr);
    ASSERT(str);

    size_t num_tokens = 0;
    while (num_tokens < tok_cap && extract_token(&tok_arr[num_tokens], str)) {
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

// Extracts token with specific delimiter
static inline size_t extract_tokens_delim(str_t tok_arr[], size_t tok_cap, str_t* str, char delim) {
    ASSERT(tok_arr);
    ASSERT(str);

    size_t num_tokens = 0;
    while (num_tokens < tok_cap && extract_token_delim(&tok_arr[num_tokens], str, delim)) {
        num_tokens += 1;
    }
    return num_tokens;
}
