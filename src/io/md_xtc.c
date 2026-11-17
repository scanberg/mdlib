#include <md_xtc.h>
#include <md_xdr.h>

#include <md_system.h>
#include <md_run.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_file.h>
#include <core/md_path.h>
#include <core/md_vec_math.h>

#include <stdio.h>

#define MD_XTC_CACHE_MAGIC   0x8281237612371
#define MD_XTC_CACHE_VERSION 6   // 6: the shared run cache header; 5: the step and box of every frame

#define XTC_MAGIC 1995
// GROMACS 2023 and later: the same frame, with the byte count of the compressed coordinates as a 64 bit integer (XDR's
// hyper: the high 32 bits first). Written for systems above 298261617 atoms (XTC_1995_MAX_NATOMS), where the
// compressed bytes of a frame no longer fit the 32 bit count, and by GROMACS's own tests at any size.
#define XTC_NEW_MAGIC 2023

static inline bool xtc_magic_ok(int32_t magic) {
    return magic == XTC_MAGIC || magic == XTC_NEW_MAGIC;
}

// The size of the byte count which follows the header of a compressed frame (natoms >= 10)
static inline size_t xtc_count_size(const uint8_t* frame_ptr) {
    return md_xdr_load_i32(frame_ptr) == XTC_NEW_MAGIC ? 8 : 4;
}

// The byte count of a compressed frame, at the end of its header
static inline uint64_t xtc_count(const uint8_t* count_ptr, size_t count_size) {
    return count_size == 8 ? md_xdr_load_u64(count_ptr) : (uint64_t)md_xdr_load_u32(count_ptr);
}

/* XTC small header size (natoms<=9).
*  > int(4) magic
*  > int(4) natoms
*  > int(4) step
*  > float(4) time
*  > 9xfloat(4) box
*  > int(4) natoms (again)
*/
#define XTC_SMALL_HEADER_SIZE 56

/* Size of uncompressed coordinates for one atom.
* 3xfloat(4) x
*/
#define XTC_SMALL_COORDS_SIZE 12

/* XTC header size (natoms>=10).
* Compressed trajectories contain some additional values:
*  > float(4) precision
*  > 3xint(4) minint
*  > 3xint(4) maxint
*  > int(4) smallidx
* See `xdrfile_compress_coord_double()`.
*/
#define XTC_HEADER_SIZE (XTC_SMALL_HEADER_SIZE + 32)

#if defined(_MSC_VER)
#define XTC_NOINLINE __declspec(noinline)
#else
#define XTC_NOINLINE __attribute__((noinline))
#endif

#if defined(__GNUC__) || defined(__clang__)
#define XTC_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define XTC_UNLIKELY(x) (x)
#endif

typedef md_128i v4i_t;
#define v4i_set(x, y, z, w) md_mm_set_epi32(w, z, y, x)
#define v4i_set1(x)         md_mm_set1_epi32(x)
#define v4i_add(a, b)       md_mm_add_epi32(a, b)
#define v4i_sub(a, b)       md_mm_sub_epi32(a, b)
#define v4i_load(addr)      md_mm_loadu_epi32(addr)

// ### Division by a constant ###
//
// The digits of a coordinate take divisions by divisors that are fixed for the frame, done as a multiply
// by a reciprocal and a shift: the method of T. Granlund and P. L. Montgomery, "Division by invariant
// integers using multiplication" (PLDI 1994), here in the round-down form with an incremented dividend
// of A. D. Robison, "N-bit unsigned division via N-bit multiply-add" (ARITH 2005), which takes a divisor
// of 1 as well. The decoder had it from libdivide (https://libdivide.com, by ridiculousfish, zlib
// license) before it was cut down to the one case it needs.
//
// For d >= 1 and 2^sh <= d: m = ceil(2^(64+sh) / d) - 1 < 2^64, m*d = 2^(64+sh) - f with 0 < f <= d, and
//     floor(v / d) == mulhi(v + 1, m) >> sh    for every v with (v + 1) * f <= 2^(64+sh),
// as with v = q*d + r, 0 <= r < d, (v+1)*m / 2^(64+sh) = q + (r+1)/d - (v+1)*f / (d * 2^(64+sh)), whose
// last two terms then sum to [0, 1). The largest shift, sh = max(ceil(log2(d)) - 1, 0), takes every
// v < 2^63 (f <= d <= 2^(sh+1)): the first coordinate of a group. A shift of 0 saves the shift where the
// dividends are small enough: the small coordinates.
typedef struct xtc_div_t {
    uint64_t m;
    uint64_t sh;
} xtc_div_t;

static MD_FORCE_INLINE uint64_t xtc_mulhi(uint64_t a, uint64_t b) {
    uint64_t hi;
    md_umul128(a, b, &hi);
    return hi;
}

static MD_FORCE_INLINE uint64_t xtc_div(uint64_t v, xtc_div_t div) {
    return xtc_mulhi(v + 1, div.m) >> div.sh;
}

// For 1 <= d < 2^63, with the largest shift: ceil(2^(64+sh) / d) by long division. The top sh + 2 bits of
// 2^(64+sh), 2^(sh+1), give a quotient of 1 (d <= 2^(sh+1) < 2d), and the other 63 go as many at a time as
// the remainder, below d, has room for shifted within 64 bits.
static inline xtc_div_t xtc_div_init(uint64_t d) {
    MD_ASSERT(1 <= d && d < (1ull << 63));
    if (d == 1) {
        const xtc_div_t one = {UINT64_MAX, 0};
        return one;
    }
    const uint32_t len = 64 - md_clz64(d - 1);     // sh + 1
    uint64_t q = 1;
    uint64_t r = (1ull << len) - d;
    for (uint32_t left = 63; left > 0;) {
        const uint32_t k = MD_MIN(left, 64 - len);
        r <<= k;
        q = (q << k) | (r / d);
        r %= d;
        left -= k;
    }
    const xtc_div_t div = {r ? q : q - 1, len - 1};
    return div;
}

// ### Bits ###
//
// The stream holds what the writer sent in order, each byte from its most significant bit. Every read
// is an unaligned 8 byte load, which is what the guard bytes past the end of the frame
// (MD_XTC_DECODE_PADDING) are for.

// num_bits (1..57) bits at bit_offset, the first the most significant (receivebits in xdrfile)
static MD_FORCE_INLINE uint64_t xtc_bits(const uint8_t* stream, size_t bit_offset, uint32_t num_bits) {
    return (md_xdr_load_u64(stream + (bit_offset >> 3)) << (bit_offset & 7)) >> (64 - num_bits);
}

// The stream from bit_offset on as if it were byte aligned there, its bytes in the order of an integer:
// byte i of the result holds the bits [bit_offset + 8i, bit_offset + 8i + 8), the first the most
// significant. The first seven bytes are whole and the eighth has its top 8 - (bit_offset & 7) bits: 57.
static MD_FORCE_INLINE uint64_t xtc_bytes(const uint8_t* stream, size_t bit_offset) {
    return md_bswap64(md_xdr_load_u64(stream + (bit_offset >> 3)) << (bit_offset & 7));
}

// As xtc_bytes with all eight bytes whole, from a second load a byte on (and no shift by 64)
static MD_FORCE_INLINE uint64_t xtc_bytes_wide(const uint8_t* stream, size_t bit_offset) {
    const uint8_t* p = stream + (bit_offset >> 3);
    const uint32_t b = bit_offset & 7;
    return md_bswap64((md_xdr_load_u64(p) << b) | (md_xdr_load_u64(p + 1) >> (8 - b)));
}

// ### Integers ###
//
// A coordinate (x, y, z) in [0, size_x) x [0, size_y) x [0, size_z) is sent as the one number
// v = (x * size_y + y) * size_z + z of num_bits bits (sendints in xdrfile): its bytes, least significant
// first, the last holding only the num_bits % 8 bits left over, if any. So v is xtc_bytes with the
// leftover bits moved down next to the whole bytes, from the top of their byte where the stream has them.
typedef struct xtc_bits_t {
    uint64_t lo_mask;      // the whole bytes
    uint64_t hi_mask;      // the leftover bits, moved down
    uint32_t hi_shift;     // how far they are moved: 8 less their number
    uint32_t num_bits;
} xtc_bits_t;

// The layout as a constant expression. Beyond 64 bits it is that of 64, which the readers of such numbers
// (xtc_read_long, xtc_read_slow) do not look at.
#define XTC_BYTES(n)  ((MD_MIN(n, 64) - 1) / 8)
#define XTC_LAST(n)   (MD_MIN(n, 64) - 8 * XTC_BYTES(n))
#define XTC_BITS(n)   {(1ull << (8 * XTC_BYTES(n))) - 1, ((1ull << XTC_LAST(n)) - 1) << (8 * XTC_BYTES(n)), 8 - XTC_LAST(n), (n)}

static inline xtc_bits_t xtc_bits_init(uint32_t num_bits) {
    MD_ASSERT(0 < num_bits && num_bits <= 72);
    const xtc_bits_t bits = XTC_BITS(num_bits);
    return bits;
}

// v of up to 57 bits (xtc_bytes), or up to 64 (wide, xtc_bytes_wide). Everything here by value: in locals
// the parameters stay in registers, through a pointer each is loaded again for every atom, as the
// coordinate stores might change them for all the compiler can tell.
static MD_FORCE_INLINE uint64_t xtc_value(const uint8_t* stream, size_t bit_offset, xtc_bits_t bits, const bool wide) {
    const uint64_t bytes = wide ? xtc_bytes_wide(stream, bit_offset) : xtc_bytes(stream, bit_offset);
    return (bytes & bits.lo_mask) | ((bytes >> bits.hi_shift) & bits.hi_mask);
}

// Small coordinates beyond XTC_SMALL_LIMB_MAX (64 to 72 bits, s > 2^21), which only a precision of 10^7
// or so reaches, with plain divisions: v / s in limbs of 32 bits, the top one first, as xtc_read_long.
// Out of line, as its run loop is inlined into the loop of every file: inlined, its divisions cost the
// ordinary files 1-3% (gcc 11) and gained the files that take it nothing.
static XTC_NOINLINE v4i_t xtc_read_slow(const uint8_t* stream, size_t bit_offset, uint32_t num_bits, uint64_t s) {
    MD_ASSERT(num_bits <= 72);
    uint64_t lo, hi = 0;
    if (num_bits <= 64) {
        lo = xtc_value(stream, bit_offset, xtc_bits_init(num_bits), true);
    } else {
        lo = xtc_bytes_wide(stream, bit_offset);
        hi = xtc_bits(stream, bit_offset + 64, num_bits - 64);
    }
    const uint64_t a = (hi << 32) | (lo >> 32);             // < 2^40
    const uint64_t b = ((a % s) << 32) | (lo & 0xFFFFFFFF);  // < s * 2^32
    const uint64_t t = ((a / s) << 32) | (b / s);            // x * s + y
    return v4i_set((uint32_t)(t / s), (uint32_t)(t % s), (uint32_t)(b % s), 0);
}

// Ranges beyond 24 bits: then the three are sent one after the other (sendbits), num_bits[i] bits each
static MD_FORCE_INLINE v4i_t xtc_read_three(const uint8_t* stream, size_t bit_offset, const uint32_t num_bits[3]) {
    uint32_t xyz[4] = {0};
    for (int i = 0; i < 3; ++i) {
        xyz[i] = (uint32_t)xtc_bits(stream, bit_offset, num_bits[i]);
        bit_offset += num_bits[i];
    }
    return v4i_load(xyz);
}

// The bits sendints takes for a number below size_x * size_y * size_z: the bit length of the product
// itself (sizeofints in xdrfile), one more than needed when it is a power of two.
static uint32_t xtc_sizeofints(const uint32_t sizes[3]) {
    uint64_t hi;
    const uint64_t lo = md_umul128((uint64_t)sizes[0] * sizes[1], sizes[2], &hi);
    return hi ? 128 - md_clz64(hi) : 64 - md_clz64(lo);
}

// The bits sendbits takes for a value below size: its bit length (sizeofint in xdrfile)
static uint32_t xtc_sizeofint(uint32_t size) {
    return 32 - md_clz32(size);
}

// ### The first coordinate of a group ###
//
// Its x and x * size_y + y each in one division of v, rather than one after the other, which leaves the
// next group less to wait for. Below 2^63 (63 bits); more, which takes ranges of 2^21 along all three
// axes, take xtc_read_long.
#define XTC_NARROW_BITS 57
#define XTC_WIDE_BITS   63

typedef struct xtc_big_t {
    xtc_bits_t bits;
    xtc_div_t  div_yz;      // by size_y * size_z: x
    xtc_div_t  div_z;       // by size_z: x * size_y + y
    xtc_div_t  div_y;       // by size_y: only beyond XTC_WIDE_BITS
    uint64_t   size_y;
    uint64_t   size_z;
} xtc_big_t;

static inline xtc_big_t xtc_big_init(uint32_t num_bits, uint32_t size_y, uint32_t size_z) {
    const xtc_big_t big = {
        .bits   = xtc_bits_init(num_bits),
        .div_yz = xtc_div_init((uint64_t)size_y * size_z),
        .div_z  = xtc_div_init(size_z),
        .div_y  = num_bits > XTC_WIDE_BITS ? xtc_div_init(size_y) : xtc_div_init(1),
        .size_y = size_y,
        .size_z = size_z,
    };
    return big;
}

static MD_FORCE_INLINE v4i_t xtc_digits_big(uint64_t v, xtc_big_t big) {
    const uint64_t x = xtc_div(v, big.div_yz);
    const uint64_t t = xtc_div(v, big.div_z);
    return v4i_set((uint32_t)x, (uint32_t)(t - x * big.size_y), (uint32_t)(v - t * big.size_z), 0);
}

// Up to XTC_WIDE_BITS. The width is the same for the whole frame, so the branch is predicted.
static MD_FORCE_INLINE v4i_t xtc_read_big(const uint8_t* stream, size_t bit_offset, xtc_big_t big) {
    if (big.bits.num_bits <= XTC_NARROW_BITS) {
        return xtc_digits_big(xtc_value(stream, bit_offset, big.bits, false), big);
    }
    return xtc_digits_big(xtc_value(stream, bit_offset, big.bits, true), big);
}

// 64 to 72 bits, v = hi * 2^64 + lo: v / size_z in two limbs of 32 bits, the top one first, then the
// quotient, x * size_y + y < 2^48, by size_y. With size_z <= 2^24 every dividend is below 2^56.
static MD_FORCE_INLINE v4i_t xtc_read_long(const uint8_t* stream, size_t bit_offset, xtc_big_t big) {
    uint64_t lo, hi = 0;
    if (big.bits.num_bits == 64) {
        lo = xtc_value(stream, bit_offset, big.bits, true);
    } else {
        lo = xtc_bytes_wide(stream, bit_offset);
        hi = xtc_bits(stream, bit_offset + 64, big.bits.num_bits - 64);
    }
    const uint64_t a  = (hi << 32) | (lo >> 32);                                 // < 2^40
    const uint64_t qa = xtc_div(a, big.div_z);
    const uint64_t b  = ((a - qa * big.size_z) << 32) | (lo & 0xFFFFFFFF);       // < size_z * 2^32
    const uint64_t qb = xtc_div(b, big.div_z);
    const uint64_t t  = (qa << 32) | qb;                                         // x * size_y + y
    const uint64_t x  = xtc_div(t, big.div_y);
    return v4i_set((uint32_t)x, (uint32_t)(t - x * big.size_y), (uint32_t)(b - qb * big.size_z), 0);
}

// ### Small coordinates ###
//
// Three values in [0, s), relative to the coordinate before, as one number of smallidx bits: s is
// magicints[smallidx] of xdrfile, s^3 <= 2^smallidx. Their divisions take a shift of 0, which makes the
// dividers m = ceil(2^64 / d) - 1 = floor((2^64 - 1) / d), and the condition (v + 1) * f <= 2^64 for
// v < s^3. For x = v / s^2 directly, it holds up to XTC_SMALL_PAIR_MAX; past it x comes from v / s, one
// division after the other, which holds up to XTC_SMALL_CHAIN_MAX (s^4 <= 2^64 would do). In one step x
// is ready sooner, which counts on cores whose window is too short to cover the chain (Skylake). Beyond,
// where only files of a precision of 10^6 or so come, up to XTC_SMALL_LIMB_MAX (s <= 2^21, v < 2^63), v / s
// is taken in two limbs of 21 bits, the top one first, and then x from v / s < s^2: every dividend u is
// below 2^42, where the same multiplier holds ((u + 1) * f <= 2^42 * 2^21). Beyond, xtc_read_slow.
typedef struct xtc_small_t {
    xtc_bits_t bits;
    uint64_t   size;
    uint64_t   m;           // for s
    uint64_t   m2;          // for s^2
} xtc_small_t;

#define XTC_SMALL_PAIR_MAX  38
#define XTC_SMALL_CHAIN_MAX 48
#define XTC_SMALL_LIMB_MAX  63
#define XTC_SMALL_LIMB_BITS 21
#define XTC_SMALL(i, s) {XTC_BITS(i), (s), UINT64_MAX / (s), UINT64_MAX / ((uint64_t)(s) * (s))}

#define FIRSTIDX 9
static const xtc_small_t xtc_small[] = {
    XTC_SMALL( 9,        8), XTC_SMALL(10,       10), XTC_SMALL(11,       12), XTC_SMALL(12,       16),
    XTC_SMALL(13,       20), XTC_SMALL(14,       25), XTC_SMALL(15,       32), XTC_SMALL(16,       40),
    XTC_SMALL(17,       50), XTC_SMALL(18,       64), XTC_SMALL(19,       80), XTC_SMALL(20,      101),
    XTC_SMALL(21,      128), XTC_SMALL(22,      161), XTC_SMALL(23,      203), XTC_SMALL(24,      256),
    XTC_SMALL(25,      322), XTC_SMALL(26,      406), XTC_SMALL(27,      512), XTC_SMALL(28,      645),
    XTC_SMALL(29,      812), XTC_SMALL(30,     1024), XTC_SMALL(31,     1290), XTC_SMALL(32,     1625),
    XTC_SMALL(33,     2048), XTC_SMALL(34,     2580), XTC_SMALL(35,     3250), XTC_SMALL(36,     4096),
    XTC_SMALL(37,     5060), XTC_SMALL(38,     6501), XTC_SMALL(39,     8192), XTC_SMALL(40,    10321),
    XTC_SMALL(41,    13003), XTC_SMALL(42,    16384), XTC_SMALL(43,    20642), XTC_SMALL(44,    26007),
    XTC_SMALL(45,    32768), XTC_SMALL(46,    41285), XTC_SMALL(47,    52015), XTC_SMALL(48,    65536),
    XTC_SMALL(49,    82570), XTC_SMALL(50,   104031), XTC_SMALL(51,   131072), XTC_SMALL(52,   165140),
    XTC_SMALL(53,   208063), XTC_SMALL(54,   262144), XTC_SMALL(55,   330280), XTC_SMALL(56,   416127),
    XTC_SMALL(57,   524287), XTC_SMALL(58,   660561), XTC_SMALL(59,   832255), XTC_SMALL(60,  1048576),
    XTC_SMALL(61,  1321122), XTC_SMALL(62,  1664510), XTC_SMALL(63,  2097152), XTC_SMALL(64,  2642245),
    XTC_SMALL(65,  3329021), XTC_SMALL(66,  4194304), XTC_SMALL(67,  5284491), XTC_SMALL(68,  6658042),
    XTC_SMALL(69,  8388607), XTC_SMALL(70, 10568983), XTC_SMALL(71, 13316085), XTC_SMALL(72, 16777216),
};
#define LASTIDX (FIRSTIDX + (int)MD_COUNTOF(xtc_small))

enum {
    XTC_SMALL_PAIR,
    XTC_SMALL_CHAIN,
    XTC_SMALL_LIMB,
    XTC_SMALL_SLOW,
};

static MD_FORCE_INLINE v4i_t xtc_read_small(const uint8_t* stream, size_t bit_offset, xtc_small_t sml, const int kind) {
    if (kind == XTC_SMALL_SLOW) {
        return xtc_read_slow(stream, bit_offset, sml.bits.num_bits, sml.size);
    }
    if (kind == XTC_SMALL_LIMB) {
        const uint64_t v  = xtc_value(stream, bit_offset, sml.bits, true);
        const uint64_t a  = v >> XTC_SMALL_LIMB_BITS;
        const uint64_t qa = xtc_mulhi(a + 1, sml.m);
        const uint64_t b  = ((a - qa * sml.size) << XTC_SMALL_LIMB_BITS) | (v & ((1ull << XTC_SMALL_LIMB_BITS) - 1));
        const uint64_t qb = xtc_mulhi(b + 1, sml.m);
        const uint64_t t  = (qa << XTC_SMALL_LIMB_BITS) | qb;      // x * s + y
        const uint64_t x  = xtc_mulhi(t + 1, sml.m);
        return v4i_set((uint32_t)x, (uint32_t)(t - x * sml.size), (uint32_t)(b - qb * sml.size), 0);
    }
    const uint64_t v = xtc_value(stream, bit_offset, sml.bits, false);
    const uint64_t t = xtc_mulhi(v + 1, sml.m);     // x * s + y
    const uint64_t x = kind == XTC_SMALL_PAIR ? xtc_mulhi(v + 1, sml.m2) : xtc_mulhi(t + 1, sml.m);
    return v4i_set((uint32_t)x, (uint32_t)(t - x * sml.size), (uint32_t)(v - t * sml.size), 0);
}

// ### Frames ###

// The whole vector, fourth lane included: the caller guarantees 16 writable bytes at dst.
static MD_FORCE_INLINE void write_coord16(float* dst, v4i_t coord, md_128 invp) {
    md_mm_storeu_ps(dst, md_mm_mul_ps(md_mm_cvtepi32_ps(coord), invp));
}

// A run of small coordinates after the big one, each relative to the one before, the first two written
// in swapped order: a water's O, the big one, comes between its Hs. Force inlined with 'kind' a constant:
// the loop for the small coordinates of every ordinary file holds no other.
static MD_FORCE_INLINE v4i_t xtc_decode_run(const uint8_t* stream, size_t* bit_offset, v4i_t thiscoord, v4i_t vsmall, xtc_small_t sml,
                                            md_128 invp, float** lfp, int run_count, const int kind) {
    const size_t sml_bits = sml.bits.num_bits;
    size_t offset = *bit_offset;
    float* out = *lfp;
    const v4i_t prevcoord = thiscoord;
    thiscoord = v4i_add(xtc_read_small(stream, offset, sml, kind), v4i_sub(thiscoord, vsmall));
    offset += sml_bits;
    write_coord16(out, thiscoord, invp); out += 3;
    write_coord16(out, prevcoord, invp); out += 3;
    for (int i = 1; i < run_count; ++i) {
        thiscoord = v4i_add(xtc_read_small(stream, offset, sml, kind), v4i_sub(thiscoord, vsmall));
        offset += sml_bits;
        write_coord16(out, thiscoord, invp); out += 3;
    }
    *bit_offset = offset;
    *lfp = out;
    return thiscoord;
}

enum {
    XTC_FORM_ONE,       // one number of up to XTC_WIDE_BITS
    XTC_FORM_LONG,      // one number of more
    XTC_FORM_THREE,     // three numbers
};

// The coordinates of a compressed frame from its stream of bits, in groups: the first coordinate of a
// group is one number (big), or when a range is beyond 24 bits three numbers of bitsizeint[i] bits, and
// a run of small coordinates may follow. Force inlined into its three callers below with 'form' a
// constant in each: the loop that decodes the frames of every ordinary trajectory is then compiled with
// no trace of the forms none of them has.
static MD_FORCE_INLINE bool xtc_decode_coords(const uint8_t* stream, uint64_t num_bytes, int natoms, int smallidx, const int32_t minint[3], float scale,
                                              float precision, xtc_big_t big, const uint32_t bitsizeint[3], float* out_coords, const int form) {
    const size_t max_bits = (size_t)num_bytes * 8;
    const size_t big_bits = form == XTC_FORM_THREE ? (size_t)bitsizeint[0] + bitsizeint[1] + bitsizeint[2] : big.bits.num_bits;
    size_t bit_offset = 0;

    xtc_small_t sml = xtc_small[smallidx - FIRSTIDX];
    v4i_t vsmall = v4i_set1((int)(sml.size / 2));

    float* lfp = out_coords;
    const md_128 invp = md_mm_set1_ps(scale / precision);
    const v4i_t vminint = v4i_set(minint[0], minint[1], minint[2], 0);
    int run_count = 0;
    int atom_idx = 0;

    // Every coordinate is written as a whole 16 byte store whose fourth lane the next one overwrites,
    // a single store instead of the two (or a trip through the stack) a 12 byte copy turns into. A
    // group writes at most 11 atoms, so the body stops 11 atoms short of the end and the last groups
    // decode into tail, from where they are copied out.
    float tail[3 * 11 + 1];
    int tail_beg = natoms;
    int limit = natoms - 11;

    for (int phase = 0; phase < 2; ++phase) {
        if (phase == 1) {
            tail_beg = atom_idx;
            lfp = tail;
            limit = natoms;
        }
        while (atom_idx < limit) {
            v4i_t thiscoord = form == XTC_FORM_THREE ? xtc_read_three(stream, bit_offset, bitsizeint) :
                              form == XTC_FORM_LONG  ? xtc_read_long(stream, bit_offset, big) : xtc_read_big(stream, bit_offset, big);
            thiscoord = v4i_add(thiscoord, vminint);
            bit_offset += big_bits;

            // How far the stream advances depends on the flag. Computed with a select, the next group
            // cannot start before this group's flag has been loaded, so every group waits on the
            // last; as a branch the advance is predicted and the groups overlap.
            const uint32_t flag = (uint32_t)xtc_bits(stream, bit_offset, 6);
            int is_smaller = 0;
            bit_offset += 1;
            if (XTC_UNLIKELY(flag & 32)) {
                bit_offset += 5;
                const int run = flag & 31;
                run_count  = run / 3;
                is_smaller = run % 3 - 1;
            }

            const int batch_size = run_count + 1;
            if (atom_idx + batch_size > natoms) {
                MD_LOG_ERROR("XTC: Buffer overrun during decompression.");
                goto done;
            }
            atom_idx += batch_size;

            if (run_count > 0) {
                if (smallidx <= XTC_SMALL_PAIR_MAX) {
                    thiscoord = xtc_decode_run(stream, &bit_offset, thiscoord, vsmall, sml, invp, &lfp, run_count, XTC_SMALL_PAIR);
                } else if (smallidx <= XTC_SMALL_CHAIN_MAX) {
                    thiscoord = xtc_decode_run(stream, &bit_offset, thiscoord, vsmall, sml, invp, &lfp, run_count, XTC_SMALL_CHAIN);
                } else if (smallidx <= XTC_SMALL_LIMB_MAX) {
                    thiscoord = xtc_decode_run(stream, &bit_offset, thiscoord, vsmall, sml, invp, &lfp, run_count, XTC_SMALL_LIMB);
                } else {
                    thiscoord = xtc_decode_run(stream, &bit_offset, thiscoord, vsmall, sml, invp, &lfp, run_count, XTC_SMALL_SLOW);
                }
            } else {
                write_coord16(lfp, thiscoord, invp); lfp += 3;
            }

            if (is_smaller) {
                smallidx += is_smaller;
                if (smallidx < FIRSTIDX || smallidx >= LASTIDX) {
                    MD_LOG_ERROR("XTC: Invalid size found in 'xdrfile_decompress_coord_float'.");
                    goto done;
                }
                sml = xtc_small[smallidx - FIRSTIDX];
                vsmall = v4i_set1((int)(sml.size / 2));
            }

            // The next group starts before the end of the stream, and with MD_XTC_DECODE_PADDING bytes past
            // it every read of a group stays in the buffer however its flags and sizes are corrupted.
            if (XTC_UNLIKELY(bit_offset > max_bits)) {
                MD_LOG_ERROR("XTC: Decompression ran past the end of the frame.");
                goto done;
            }
        }
    }

done:
    if (atom_idx > tail_beg) {
        md_memcpy(out_coords + (size_t)tail_beg * 3, tail, (size_t)(atom_idx - tail_beg) * 3 * sizeof(float));
    }
    return atom_idx == natoms;
}

// Out of line, each with its copy of the loop: inlined into xtc_decode_frame_data_scaled (as gcc does when
// let), the ordinary files gained 0-1% and frames of precision 1e6 lost some 10% (clang 18).
static XTC_NOINLINE bool xtc_decode_coords_three(const uint8_t* stream, uint64_t num_bytes, int natoms, int smallidx, const int32_t minint[3],
                                                 float scale, float precision, const uint32_t bitsizeint[3], float* out_coords) {
    const xtc_big_t no_big = {0};
    return xtc_decode_coords(stream, num_bytes, natoms, smallidx, minint, scale, precision, no_big, bitsizeint, out_coords, XTC_FORM_THREE);
}

static XTC_NOINLINE bool xtc_decode_coords_long(const uint8_t* stream, uint64_t num_bytes, int natoms, int smallidx, const int32_t minint[3],
                                                float scale, float precision, uint32_t num_bits, uint32_t size_y, uint32_t size_z, float* out_coords) {
    const xtc_big_t big = xtc_big_init(num_bits, size_y, size_z);
    const uint32_t no_bitsizeint[3] = {0};
    return xtc_decode_coords(stream, num_bytes, natoms, smallidx, minint, scale, precision, big, no_bitsizeint, out_coords, XTC_FORM_LONG);
}

static inline bool decode_header(const uint8_t* frame_ptr, md_xtc_header_t* out_header) {
    // Extract header
    int magic;
    magic = md_xdr_load_i32(frame_ptr);
    if (!xtc_magic_ok(magic)) {
        MD_LOG_ERROR("XTC: Magic number did not match");
        return false;
	}

    out_header->natoms = md_xdr_load_i32(frame_ptr + 4);
    out_header->step = md_xdr_load_i32(frame_ptr + 8);
	out_header->time = md_xdr_load_f32(frame_ptr + 12);
	md_xdr_load_f32_array((float*)out_header->box, frame_ptr + 16, 9);

    return true;
}

// The scan reads every frame header on its way through the file, and a header holds the step and
// the box as well as the time. steps and boxes may be NULL; boxes takes nine floats per frame, nm.
static size_t xtc_scan(md_file_t xdr, md_array(int64_t)* frame_offsets, md_array(double)* frame_times,
                       md_array(int64_t)* frame_steps, md_array(float)* frame_boxes, md_allocator_i* alloc);

size_t md_xtc_read_frame_offsets_and_times(md_file_t xdr, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, md_allocator_i* alloc) {
    return xtc_scan(xdr, frame_offsets, frame_times, NULL, NULL, alloc);
}

static void xtc_scan_push_header(md_array(int64_t)* frame_steps, md_array(float)* frame_boxes, const md_xtc_header_t* h, md_allocator_i* alloc) {
    if (frame_steps) {
        md_array_push(*frame_steps, (int64_t)h->step, alloc);
    }
    if (frame_boxes) {
        const float* box = &h->box[0][0];
        for (int i = 0; i < 9; ++i) {
            md_array_push(*frame_boxes, box[i], alloc);
        }
    }
}

static size_t xtc_scan(md_file_t xdr, md_array(int64_t)* frame_offsets, md_array(double)* frame_times,
                       md_array(int64_t)* frame_steps, md_array(float)* frame_boxes, md_allocator_i* alloc) {
    size_t filesize = (size_t)md_file_size(xdr);

    if (filesize == 0) {
        MD_LOG_ERROR("XTC: Failed extract filesize");
        return 0;
    }

    /* Go to file beg */
    if (!md_file_seek(xdr, 0, MD_FILE_BEG)) {
        MD_LOG_ERROR("XTC: Failed to seek to beginning of file");
        return 0;
    }

    uint8_t frame_header_data[XTC_HEADER_SIZE];
    md_xtc_header_t xtc_header;
    size_t read_bytes;

    read_bytes = md_file_read(xdr, frame_header_data, XTC_SMALL_HEADER_SIZE);
    if (read_bytes != XTC_SMALL_HEADER_SIZE || !decode_header(frame_header_data, &xtc_header)) {
        MD_LOG_ERROR("XTC: File does not appear to be a valid xtc trajectory");
        return 0;
    }

    if (xtc_header.natoms <= 0) {
        MD_LOG_ERROR("XTC: Invalid number of atoms in header");
        return 0;
    }

    size_t num_frames = 0;

    // Push first frame
    md_array_push(*frame_offsets, 0, alloc);
    md_array_push(*frame_times, xtc_header.time, alloc);
    xtc_scan_push_header(frame_steps, frame_boxes, &xtc_header, alloc);
    num_frames += 1;

    /* Dont bother with compression for nine atoms or less */
    if (xtc_header.natoms <= 9) {
        const size_t framebytes = XTC_SMALL_HEADER_SIZE + XTC_SMALL_COORDS_SIZE * xtc_header.natoms;
        const size_t est_frames = (filesize / framebytes); /* Should we complain if framesize doesn't divide filesize? */

        md_array_ensure(*frame_offsets, est_frames, alloc);
        md_array_ensure(*frame_times,   est_frames, alloc);

        for (size_t i = 1; i < est_frames; i++) {
            const size_t offset = i * framebytes;
            
            bool success = false;
            if (md_file_seek(xdr, offset, MD_FILE_BEG)) {
                read_bytes = md_file_read(xdr, frame_header_data, XTC_SMALL_HEADER_SIZE);
                success = (read_bytes == XTC_SMALL_HEADER_SIZE) && decode_header(frame_header_data, &xtc_header);
            }

            // Push frame i
            if (success) {
                md_array_push(*frame_offsets, offset, alloc);
                md_array_push(*frame_times, xtc_header.time, alloc);
                xtc_scan_push_header(frame_steps, frame_boxes, &xtc_header, alloc);
                num_frames += 1;
            } else {
               MD_LOG_DEBUG("XTC: encountered corrupted frame header");
               break;
            }
        }
        md_array_push(*frame_offsets, num_frames * framebytes, alloc);
    } else {
        // Every frame is its header, the byte count of its compressed coordinates and those bytes,
        // padded to 4. One positioned read per frame, at an offset kept here: no tell or seek, and
        // nothing that can step past the end of the file.
        uint8_t raw[XTC_HEADER_SIZE + 8];

        const size_t first_count = xtc_count_size(frame_header_data);
        if (md_file_read_at(xdr, XTC_HEADER_SIZE, raw, first_count) != first_count) {
            MD_LOG_ERROR("XTC: Failed to read framebytes");
            return 0;
        }
        const uint64_t first_bytes = xtc_count(raw, first_count);
        if (first_bytes > filesize) {
            MD_LOG_ERROR("XTC: The first frame is cut short");
            return 0;
        }
        int64_t offset = XTC_HEADER_SIZE + (int64_t)first_count + (int64_t)md_xdr_padded_size((size_t)first_bytes);
        if (offset > (int64_t)filesize) {
            MD_LOG_ERROR("XTC: The first frame is cut short");
            return 0;
        }

        size_t est_nframes = filesize / (size_t)offset + 1;
        /* First `framebytes` might be larger than average, so we would underestimate `est_nframes`*/
        est_nframes += est_nframes / 5;

        md_array_ensure(*frame_offsets, est_nframes, alloc);
        md_array_ensure(*frame_times,   est_nframes, alloc);

        while (offset < (int64_t)filesize) {
            // The header and a 32 bit count: what every frame has, and all a frame of the old magic has at the end of a file
            if (md_file_read_at(xdr, offset, raw, XTC_HEADER_SIZE + 4) != XTC_HEADER_SIZE + 4 || !decode_header(raw, &xtc_header)) {
                MD_LOG_DEBUG("XTC: encountered corrupted frame header");
                break;
            }
            const size_t count_size = xtc_count_size(raw);
            if (count_size == 8 && md_file_read_at(xdr, offset + XTC_HEADER_SIZE + 4, raw + XTC_HEADER_SIZE + 4, 4) != 4) {
                MD_LOG_INFO("XTC: The last frame is cut short and is left out");
                break;
            }
            const uint64_t bytes = xtc_count(raw + XTC_HEADER_SIZE, count_size);
            const int64_t next = bytes > filesize ? (int64_t)filesize + 1 :
                offset + XTC_HEADER_SIZE + (int64_t)count_size + (int64_t)md_xdr_padded_size((size_t)bytes);
            if (next > (int64_t)filesize) {
                // A run killed mid write leaves its last frame incomplete; it is not a frame.
                MD_LOG_INFO("XTC: The last frame is cut short and is left out");
                break;
            }
            md_array_push(*frame_offsets, offset, alloc);
            md_array_push(*frame_times, xtc_header.time, alloc);
            xtc_scan_push_header(frame_steps, frame_boxes, &xtc_header, alloc);
            num_frames += 1;
            offset = next;
        }
        // The end of the last whole frame: the file size, or where a damaged tail begins. Always
        // there, as every reader takes the size of frame i to be offsets[i + 1] - offsets[i].
        md_array_push(*frame_offsets, offset, alloc);
    }
    return num_frames;
}

// The coordinates xyz packed, the layout the file holds them in, scaled on the way out: scale 10 turns the
// file's nm into Angstrom with no second pass over the coordinates.
static bool xtc_decode_frame_data_scaled(const uint8_t* frame_ptr, size_t frame_bytes, md_xtc_header_t* out_header, float* out_coords, size_t num_atoms, float scale) {
    if (frame_ptr == NULL || frame_bytes == 0) {
        return false;
    }

    if (frame_bytes < XTC_SMALL_HEADER_SIZE) {
        MD_LOG_ERROR("XTC: Frame size is too small to contain header");
        return false;
	}

    if (out_header) {
        if (!decode_header(frame_ptr, out_header)) {
            return false;
        }
        float* box = &out_header->box[0][0];
        for (int i = 0; i < 9; ++i) {
            box[i] *= scale;
        }
	}

    if (!out_coords) {
        return true;
    }

    int natoms;
	natoms = md_xdr_load_i32(frame_ptr + XTC_SMALL_HEADER_SIZE - 4);

    if (natoms != (int)num_atoms) {
        MD_LOG_ERROR("XTC: Number of atoms in frame header does not match expected number of atoms");
        return false;
    }

    size_t offset = XTC_SMALL_HEADER_SIZE;
    if (natoms <= 9) {
        if (frame_bytes < offset + (size_t)natoms * XTC_SMALL_COORDS_SIZE) {
            MD_LOG_ERROR("XTC: Frame size is too small to contain its coordinates");
            return false;
        }
		// No compression for 9 atoms or less, just read the coordinates directly
		md_xdr_load_f32_array(out_coords, frame_ptr + offset, (size_t)natoms * 3);
        for (size_t i = 0; i < (size_t)natoms * 3; ++i) {
            out_coords[i] *= scale;
        }
        return true;
    }

    const size_t count_size = xtc_count_size(frame_ptr);
    if (frame_bytes < XTC_HEADER_SIZE + count_size) {
        MD_LOG_ERROR("XTC: Frame size is too small to contain header");
        return false;
    }

    float precision;
    int32_t minint[3], maxint[3], smallidx;
    uint64_t num_bytes;
	precision = md_xdr_load_f32(frame_ptr + offset); offset += 4;
	md_xdr_load_i32_array(minint, frame_ptr + offset, 3); offset += 12;
	md_xdr_load_i32_array(maxint, frame_ptr + offset, 3); offset += 12;
	smallidx  = md_xdr_load_i32(frame_ptr + offset); offset += 4;
    num_bytes = xtc_count(frame_ptr + offset, count_size); offset += count_size;

    // A corrupt header fails here rather than indexing past the tables or sending the decoder past the
    // end of the frame.
    if (smallidx < FIRSTIDX || smallidx >= (int)LASTIDX) {
        MD_LOG_ERROR("XTC: Invalid small index in frame header");
        return false;
    }
    if (!(precision > 0.0f) || !isfinite(precision)) {
        MD_LOG_ERROR("XTC: Invalid precision in frame header");
        return false;
    }
    if (num_bytes > frame_bytes - offset) {
        MD_LOG_ERROR("XTC: Frame holds fewer bytes than its header states");
        return false;
    }

    uint32_t sizeint[3];
    for (int i = 0; i < 3; ++i) {
        const int64_t range = (int64_t)maxint[i] - (int64_t)minint[i] + 1;
        if (range < 1 || range > (int64_t)UINT32_MAX) {
            MD_LOG_ERROR("XTC: Invalid coordinate range in frame header");
            return false;
        }
        sizeint[i] = (uint32_t)range;
    }

    const uint8_t* stream = frame_ptr + offset;
    if (XTC_UNLIKELY((sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff)) {
        const uint32_t bitsizeint[3] = {xtc_sizeofint(sizeint[0]), xtc_sizeofint(sizeint[1]), xtc_sizeofint(sizeint[2])};
        return xtc_decode_coords_three(stream, num_bytes, natoms, smallidx, minint, scale, precision, bitsizeint, out_coords);
    }
    // Three ranges of at most 24 bits: one number of at most 72
    const uint32_t num_bits = xtc_sizeofints(sizeint);
    if (XTC_UNLIKELY(num_bits > XTC_WIDE_BITS)) {
        return xtc_decode_coords_long(stream, num_bytes, natoms, smallidx, minint, scale, precision, num_bits, sizeint[1], sizeint[2], out_coords);
    }
    const xtc_big_t big = xtc_big_init(num_bits, sizeint[1], sizeint[2]);
    const uint32_t no_bitsizeint[3] = {0};
    return xtc_decode_coords(stream, num_bytes, natoms, smallidx, minint, scale, precision, big, no_bitsizeint, out_coords, XTC_FORM_ONE);
}

bool md_xtc_decode_frame_data(const uint8_t* frame_ptr, size_t frame_bytes, md_xtc_header_t* out_header, float* out_coords, size_t num_atoms) {
    return xtc_decode_frame_data_scaled(frame_ptr, frame_bytes, out_header, out_coords, num_atoms, 1.0f);
}

// Everything the scan learns about the file, which is everything short of the coordinates: where
// each frame is, and the time, step and box its header states. Kept beside the trajectory as
// '<file>.cache' so the scan is paid once per file rather than once per load.
typedef struct xtc_cache_t {
    md_run_cache_header_t header;
    int64_t* frame_offsets;   // num_frames + 1: the last one is the file size
    double*  frame_times;     // ps
    int64_t* frame_steps;
    float*   frame_boxes;     // 9 per frame, nm, row major
} xtc_cache_t;

static bool cache_read_block(md_file_t file, void** dst, size_t bytes, md_allocator_i* alloc) {
    *dst = md_alloc(alloc, bytes);
    return *dst && md_file_read(file, *dst, bytes) == bytes;
}

// On failure nothing is freed here: alloc is always an arena owned by the caller, which goes as a
// whole, so there is no partial state to unwind.
static bool try_read_cache(xtc_cache_t* cache, md_str_t path, md_allocator_i* alloc) {
    md_file_t file = {0};
    if (!md_run_cache_open(&file, &cache->header, path, MD_XTC_CACHE_MAGIC, MD_XTC_CACHE_VERSION)) {
        return false;
    }
    const size_t n = cache->header.num_frames;
    const bool ok =
        cache_read_block(file, (void**)&cache->frame_offsets, (n + 1) * sizeof(int64_t), alloc) &&
        cache_read_block(file, (void**)&cache->frame_times,   n * sizeof(double),        alloc) &&
        cache_read_block(file, (void**)&cache->frame_steps,   n * sizeof(int64_t),       alloc) &&
        cache_read_block(file, (void**)&cache->frame_boxes,   n * 9 * sizeof(float),     alloc) &&
        md_file_tell(file) == (int64_t)md_file_size(file);
    if (!ok) {
        MD_LOG_ERROR("XTC: the cache of '" MD_STR_FMT "' is incomplete", MD_STR_ARG(path));
    }
    md_file_close(&file);
    return ok;
}

static bool write_cache(const xtc_cache_t* cache, md_str_t path, const md_file_info_t* scanned) {
    md_file_t file = {0};
    const size_t n = cache->header.num_frames;
    if (!md_run_cache_create(&file, path, scanned, MD_XTC_CACHE_MAGIC, MD_XTC_CACHE_VERSION, cache->header.num_atoms, n)) {
        return false;
    }
    const bool ok =
        md_file_write(file, cache->frame_offsets, (n + 1) * sizeof(int64_t)) == (n + 1) * sizeof(int64_t) &&
        md_file_write(file, cache->frame_times,   n * sizeof(double))        == n * sizeof(double) &&
        md_file_write(file, cache->frame_steps,   n * sizeof(int64_t))       == n * sizeof(int64_t) &&
        md_file_write(file, cache->frame_boxes,   n * 9 * sizeof(float))     == n * 9 * sizeof(float);
    if (!ok) {
        MD_LOG_ERROR("XTC: failed to write the cache of '" MD_STR_FMT "'", MD_STR_ARG(path));
    }
    md_file_close(&file);
    return ok;
}

// What the file is, from its cache when that is current and from a scan otherwise (writing the
// cache unless told not to). path must be canonical; every array is allocated from alloc.
static bool xtc_index_load(xtc_cache_t* cache, md_str_t path, uint32_t flags, md_allocator_i* alloc) {
    md_memset(cache, 0, sizeof(*cache));

    md_file_t file = {0};
    if (!md_file_open(&file, path, MD_FILE_READ)) {
        MD_LOG_ERROR("XTC: Failed to open file '" MD_STR_FMT "'", MD_STR_ARG(path));
        return false;
    }

    bool result = false;
    md_file_info_t scanned = {0};
    md_file_info_extract(file, &scanned);

    uint8_t frame_header_data[XTC_SMALL_HEADER_SIZE];
    md_xtc_header_t xtc_header = { 0 };
    if (md_file_read(file, frame_header_data, XTC_SMALL_HEADER_SIZE) != XTC_SMALL_HEADER_SIZE || !decode_header(frame_header_data, &xtc_header)) {
        MD_LOG_ERROR("XTC: Failed to read header of first frame, file may be corrupt or not a valid xtc trajectory");
        goto done;
    }
    if (xtc_header.natoms <= 0) {
        MD_LOG_ERROR("XTC: Number of atoms in trajectory was zero");
        goto done;
    }

    if (try_read_cache(cache, path, alloc)) {
        result = true;
        goto done;
    }

    md_memset(cache, 0, sizeof(*cache));
    cache->header.num_atoms = xtc_header.natoms;

    md_array(int64_t) offsets = 0;
    md_array(double)  times   = 0;
    md_array(int64_t) steps   = 0;
    md_array(float)   boxes   = 0;
    cache->header.num_frames = xtc_scan(file, &offsets, &times, &steps, &boxes, alloc);
    if (!cache->header.num_frames || !offsets || !times) {
        MD_LOG_DEBUG("XTC: frame offsets or frame times was empty");
        goto done;
    }
    cache->frame_offsets = offsets;
    cache->frame_times   = times;
    cache->frame_steps   = steps;
    cache->frame_boxes   = boxes;

    if (!(flags & MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
        // If we fail to write the cache, that's ok, we can inform about it, but do not halt
        if (write_cache(cache, path, &scanned)) {
            MD_LOG_INFO("XTC: Successfully created cache file for '" MD_STR_FMT "'", MD_STR_ARG(path));
        }
    }
    result = true;

done:
    md_file_close(&file);
    return result;
}

// ### RUN ###

// <run>/atom/position. Everything it needs beyond the slice is in the table: the run is the path
// minus "/atom/position", and the file, offset and size of the frame are its source attributes.
static size_t xtc_position_provider(void* dst, size_t cap, const md_attributes_t* table, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    (void)user_data;
    const md_attributes_t* attributes = table;

    // A whole temporal virtual attribute is refused before it gets here, so the frame is fixed.
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) {
        return 0;
    }

    md_run_source_t src;
    if (!md_run_source(&src, attributes, attr, MD_STR_LIT("atom/position"))) {
        return 0;
    }
    const int64_t* offsets = src.offset;
    const int64_t* sizes   = src.size;

    const uint32_t frame = slice->idx[0];
    const size_t num_atoms = attr->format.shape[1];
    const size_t frame_size = (size_t)sizes[frame];

    md_temp_scope_t temp = md_temp_begin();
    size_t written = 0;

    // The decoder reads the bit stream a word at a time and may look past the end of the frame.
    uint8_t* frame_data = md_temp_alloc(temp, MD_ALIGN_TO(frame_size, 16) + MD_XTC_DECODE_PADDING);
    if (frame_data) {
        md_memset(frame_data + frame_size, 0, MD_ALIGN_TO(frame_size, 16) + MD_XTC_DECODE_PADDING - frame_size);
    }
    // Through io: inside an extraction context the file stays open from frame to frame, which is
    // what makes streaming from a cluster's file server bearable. Without one it is opened here.
    const md_str_t file_path = src.path;
    const size_t read = frame_data ? md_attribute_io_read_at(io, file_path, offsets[frame], frame_data, frame_size) : 0;
    if (frame_data && read != frame_size) {
        MD_LOG_ERROR("XTC: Failed to read frame %u from '" MD_STR_FMT "', expected %zu bytes, got %zu", frame, MD_STR_ARG(file_path), frame_size, read);
    }
    if (frame_data && read == frame_size) {
        if (slice->num_idx == 1 && cap == num_atoms * 3) {
            // The whole frame, which is the case that matters: decoded straight into the caller.
            if (xtc_decode_frame_data_scaled(frame_data, frame_size, NULL, (float*)dst, num_atoms, 10.0f)) {
                written = cap;
            }
        } else if (slice->num_idx == 2 && cap == 3 && slice->idx[1] < num_atoms) {
            // One atom: the stream is not seekable, so the frame is decoded and one value kept.
            float* xyz = md_temp_alloc(temp, num_atoms * 3 * sizeof(float));
            if (xyz && xtc_decode_frame_data_scaled(frame_data, frame_size, NULL, xyz, num_atoms, 10.0f)) {
                md_memcpy(dst, xyz + (size_t)slice->idx[1] * 3, 3 * sizeof(float));
                written = cap;
            }
        }
    }

    md_temp_end(temp);
    return written;
}

bool md_xtc_system_publish_run(md_system_t* sys, md_str_t filename, md_str_t run, uint32_t flags) {
    MD_ASSERT(sys);
    if (md_str_empty(run)) {
        MD_LOG_ERROR("XTC: no run to publish into");
        return false;
    }
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    const md_str_t path = {path_buf, path_len};

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    bool result = false;

    xtc_cache_t index;
    if (!xtc_index_load(&index, path, flags, arena)) {
        goto done;
    }

    const size_t F = index.header.num_frames;
    const size_t N = index.header.num_atoms;

    int64_t* sizes = md_alloc(arena, F * sizeof(int64_t));
    float*   boxes = md_alloc(arena, F * 9 * sizeof(float));
    for (size_t i = 0; i < F; ++i) {
        sizes[i] = index.frame_offsets[i + 1] - index.frame_offsets[i];
    }
    for (size_t i = 0; i < F * 9; ++i) {
        boxes[i] = index.frame_boxes[i] * 10.0f;   // nm to Angstrom, as the coordinates
    }

    const md_attribute_virtual_t virt = { .provider = xtc_position_provider };
    const md_run_desc_t desc = {
        .num_frames    = F,
        .num_atoms     = N,
        .time          = index.frame_times,
        .time_unit     = md_unit_picosecond(),
        .step          = index.frame_steps,
        .unitcell      = boxes,
        .source_path   = path,
        .source_offset = index.frame_offsets,
        .source_size   = sizes,
        .position_virt = &virt,
    };
    result = md_run_publish(sys, run, &desc);

done:
    md_arena_allocator_destroy(arena);
    return result;
}

size_t md_xtc_file_atom_count(md_str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return 0;
    uint8_t head[8];
    const bool ok = md_file_read(file, head, sizeof(head)) == sizeof(head);
    md_file_close(&file);
    if (!ok || !xtc_magic_ok(md_xdr_load_i32(head))) return 0;
    const int32_t natoms = md_xdr_load_i32(head + 4);
    return natoms > 0 ? (size_t)natoms : 0;
}
