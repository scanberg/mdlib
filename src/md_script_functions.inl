// This is to mark that the procedure supports a varying length
#define ANY_LENGTH -1

#ifndef MD_DIST_BINS
#define MD_DIST_BINS 1024
#endif

#ifndef MD_VOL_DIM
#define MD_VOL_DIM 128
#endif

#define STATIC_VALIDATION_ERROR -2

#define BAKE_STR(str) {str"", sizeof(str)-1}

// Type info declaration helpers
#define TI_BOOL         {TYPE_BOOL, {1}, 0}

#define TI_FLOAT        {TYPE_FLOAT, {1}, 0}
#define TI_FLOAT_ARR    {TYPE_FLOAT, {ANY_LENGTH}, 0}
#define TI_FLOAT2       {TYPE_FLOAT, {2,1}, 1}
#define TI_FLOAT2_ARR   {TYPE_FLOAT, {2,ANY_LENGTH}, 1}
#define TI_FLOAT3       {TYPE_FLOAT, {3,1}, 1}
#define TI_FLOAT3_ARR   {TYPE_FLOAT, {3,ANY_LENGTH}, 1}
#define TI_FLOAT4       {TYPE_FLOAT, {4,1}, 1}
#define TI_FLOAT4_ARR   {TYPE_FLOAT, {4,ANY_LENGTH}, 1}
#define TI_FLOAT44      {TYPE_FLOAT, {4,4,1}, 2}
#define TI_FLOAT44_ARR  {TYPE_FLOAT, {4,4,ANY_LENGTH}, 2}

// The second dimension in the distribution encodes weights for each bin
#define TI_DISTRIBUTION {TYPE_FLOAT, {MD_DIST_BINS,2}, 1}
#define TI_VOLUME       {TYPE_FLOAT, {MD_VOL_DIM, MD_VOL_DIM, MD_VOL_DIM,1}, 3}

#define TI_INT          {TYPE_INT, {1}, 0}
#define TI_INT_ARR      {TYPE_INT, {ANY_LENGTH}, 0}

#define TI_FRANGE       {TYPE_FRANGE, {1}, 0}
#define TI_FRANGE_ARR   {TYPE_FRANGE, {ANY_LENGTH}, 0}

#define TI_IRANGE       {TYPE_IRANGE, {1}, 0}
#define TI_IRANGE_ARR   {TYPE_IRANGE, {ANY_LENGTH}, 0}

#define TI_STRING       {TYPE_STRING, {1}, 0}
#define TI_STRING_ARR   {TYPE_STRING, {ANY_LENGTH}, 0}

#define TI_BITFIELD     {TYPE_BITFIELD, {1}, 0}
#define TI_BITFIELD_ARR {TYPE_BITFIELD, {ANY_LENGTH}, 0}

#define TI_COORDINATE     {TYPE_COORDINATE, {1}, 0}
#define TI_COORDINATE_ARR {TYPE_COORDINATE, {ANY_LENGTH}, 0}

#define as_float(arg) (*((float*)((arg).ptr)))
#define as_float_arr(arg) ((float*)((arg).ptr))

#define as_int(arg) (*((int*)((arg).ptr)))
#define as_int_arr(arg) ((int*)((arg).ptr))

#define as_frange(arg) (*((frange_t*)((arg).ptr)))
#define as_frange_arr(arg) ((frange_t*)((arg).ptr))

#define as_irange(arg) (*((irange_t*)((arg).ptr)))
#define as_irange_arr(arg) ((irange_t*)((arg).ptr))

#define as_string(arg) (*((str_t*)((arg).ptr)))
#define as_string_arr(arg) ((str_t*)((arg).ptr))

#define as_bitfield(arg) ((md_bitfield_t*)((arg).ptr))

#define as_vec2(arg) (*((vec2_t*)((arg).ptr)))
#define as_vec2_arr(arg) (((vec2_t*)((arg).ptr)))

#define as_vec3(arg) (*((vec3_t*)((arg).ptr)))
#define as_vec3_arr(arg) (((vec3_t*)((arg).ptr)))

#define as_vec4(arg) (*((vec4_t*)((arg).ptr)))
#define as_vec4_arr(arg) (((vec4_t*)((arg).ptr)))

// Macros to bake standard functions into something callable through our interface
#define BAKE_FUNC_F__F(prefix, func) \
    static int prefix##func(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        as_float(*dst) = func(as_float(arg[0])); \
        return 0; \
    }

#define BAKE_FUNC_F__F_F(prefix, func) \
    static int prefix##func(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        as_float(*dst) = func(as_float(arg[0]), as_float(arg[1])); \
        return 0; \
    }

#define BAKE_FUNC_FARR__FARR(prefix, func) \
    static int prefix##func(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        if (dst) { \
            for (size_t i = 0; i < element_count(*dst); ++i) { \
                as_float_arr(*dst)[i] = func(as_float_arr(arg[0])[i]); \
            } \
            return 0; \
        } \
        else { \
            return (int)element_count(arg[0]); \
        } \
    }

// float func(float)
BAKE_FUNC_F__F(_, sqrtf)
BAKE_FUNC_F__F(_, cbrtf)
BAKE_FUNC_F__F(_, fabsf)
BAKE_FUNC_F__F(_, floorf)
BAKE_FUNC_F__F(_, ceilf)
BAKE_FUNC_F__F(_, cosf)
BAKE_FUNC_F__F(_, sinf)
BAKE_FUNC_F__F(_, tanf)
BAKE_FUNC_F__F(_, acosf)
BAKE_FUNC_F__F(_, asinf)
BAKE_FUNC_F__F(_, atanf)
BAKE_FUNC_F__F(_, logf)
BAKE_FUNC_F__F(_, expf)
BAKE_FUNC_F__F(_, log2f)
BAKE_FUNC_F__F(_, exp2f)
BAKE_FUNC_F__F(_, log10f)

// float func(float, float)
BAKE_FUNC_F__F_F(_, atan2f)
BAKE_FUNC_F__F_F(_, powf)
BAKE_FUNC_F__F_F(_, fminf)
BAKE_FUNC_F__F_F(_, fmaxf)

// array versions
// @TODO: Fill in these when needed. It's a bit uncelar up front what should be supported and exposed as array operations?
BAKE_FUNC_FARR__FARR(_arr_, fabsf)
BAKE_FUNC_FARR__FARR(_arr_, floorf)
BAKE_FUNC_FARR__FARR(_arr_, ceilf)

// MACROS TO BAKE OPERATORS INTO CALLABLE FUNCTIONS
#define BAKE_OP_UNARY_S(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        if (!dst) return 0; \
        *((base_type*)dst->ptr) = op *((base_type*)arg[0].ptr); \
        return 0; \
    }

#define BAKE_OP_UNARY_M(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        if (!dst) return 0; \
        for (size_t i = 0; i < element_count(*dst); ++i) { \
            ((base_type*)dst->ptr)[i] = op ((base_type*)arg[0].ptr)[i]; \
        } \
        return 0; \
    }

#define BAKE_OP_S_S(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        if (!dst) return 0; \
        *((base_type*)dst->ptr) = *((base_type*)arg[0].ptr) op *((base_type*)arg[1].ptr); \
        return 0; \
    }

#define BAKE_OP_M_S(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        if (!dst) return 0; \
        for (size_t i = 0; i < element_count(*dst); ++i) { \
            ((base_type*)dst->ptr)[i] = ((base_type*)arg[0].ptr)[i] op *((base_type*)arg[1].ptr); \
        } \
        return 0; \
    }

#define BAKE_OP_M_M(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        if (!dst) return 0; \
        for (size_t i = 0; i < element_count(*dst); ++i) { \
            ((base_type*)dst->ptr)[i] = ((base_type*)arg[0].ptr)[i] op ((base_type*)arg[1].ptr)[i]; \
        } \
        return 0; \
    }

#define BAKE_SIMDF_OP_M_M(name, op) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        ASSERT(dst); \
        ASSERT(dst->type.base_type == TYPE_FLOAT); \
        ASSERT(arg[0].type.base_type == TYPE_FLOAT); \
        ASSERT(arg[1].type.base_type == TYPE_FLOAT); \
        const size_t count = type_info_element_stride_count(arg[0].type); \
        ASSERT(count == type_info_element_stride_count(arg[1].type)); \
        ASSERT(count % 8 == 0); \
        const float* src_a = as_float_arr(arg[0]); \
        const float* src_b = as_float_arr(arg[1]); \
        float* dst_arr = as_float_arr(*dst); \
        for (size_t i = 0; i < count; i += 8) { \
            md_256 a = md_mm256_loadu_ps(src_a + i); \
            md_256 b = md_mm256_loadu_ps(src_b + i); \
            md_mm256_storeu_ps(dst_arr + i, op(a, b)); \
        } \
        return 0; \
    }

#define BAKE_SIMDF_OP_M_S(name, op) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        ASSERT(dst); \
        ASSERT(dst->type.base_type == TYPE_FLOAT); \
        ASSERT(arg[0].type.base_type == TYPE_FLOAT); \
        ASSERT(is_type_equivalent(arg[1].type, (type_info_t)TI_FLOAT)); \
        const size_t count = type_info_element_stride_count(arg[0].type); \
        ASSERT(count % 8 == 0); \
        const float* src_arr = as_float_arr(arg[0]); \
        float* dst_arr = as_float_arr(*dst); \
        md_256 s = md_mm256_set1_ps(as_float(arg[1])); \
        for (size_t i = 0; i < count; i += 8) { \
            md_mm256_storeu_ps(dst_arr + i, op(md_mm256_loadu_ps(src_arr + i), s)); \
        } \
        return 0; \
    }

#define BAKE_SIMDF_OP_M(name, op) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        ASSERT(dst); \
        ASSERT(dst->type.base_type == TYPE_FLOAT); \
        ASSERT(arg[0].type.base_type == TYPE_FLOAT); \
        const size_t count = type_info_element_stride_count(arg[0].type); \
        ASSERT(count % 8 == 0); \
        const float* src_arr = as_float_arr(arg[0]); \
        float* dst_arr = as_float_arr(*dst); \
        for (size_t i = 0; i < count; i += 8) { \
            md_mm256_storeu_ps(dst_arr + i, op(md_mm256_loadu_ps(src_arr + i))); \
        } \
        return 0; \
    }

BAKE_OP_S_S(_op_and_b_b,  &&, bool)
BAKE_OP_S_S(_op_or_b_b,   ||, bool)
BAKE_OP_S_S(_op_xor_b_b,   ^, bool)
BAKE_OP_UNARY_S(_op_not_b, !, bool)

BAKE_OP_S_S(_op_add_i_i, +, int)
BAKE_OP_S_S(_op_sub_i_i, -, int)
BAKE_OP_S_S(_op_mul_i_i, *, int)
BAKE_OP_S_S(_op_div_i_i, /, int)
BAKE_OP_UNARY_S(_op_neg_i, -, int)
BAKE_OP_UNARY_M(_op_neg_iarr, -, int)

BAKE_OP_M_S(_op_add_iarr_i, +, int)
BAKE_OP_M_S(_op_sub_iarr_i, -, int)
BAKE_OP_M_S(_op_mul_iarr_i, *, int)
BAKE_OP_M_S(_op_div_iarr_i, /, int)

BAKE_OP_M_M(_op_add_iarr_iarr, +, int)
BAKE_OP_M_M(_op_sub_iarr_iarr, -, int)
BAKE_OP_M_M(_op_mul_iarr_iarr, *, int)
BAKE_OP_M_M(_op_div_iarr_iarr, /, int)

BAKE_SIMDF_OP_M_M(_op_simd_add_farr_farr, md_mm256_add_ps)
BAKE_SIMDF_OP_M_S(_op_simd_add_farr_f,    md_mm256_add_ps)

BAKE_SIMDF_OP_M_M(_op_simd_sub_farr_farr, md_mm256_sub_ps)
BAKE_SIMDF_OP_M_S(_op_simd_sub_farr_f,    md_mm256_sub_ps)

BAKE_SIMDF_OP_M_M(_op_simd_mul_farr_farr, md_mm256_mul_ps)
BAKE_SIMDF_OP_M_S(_op_simd_mul_farr_f,    md_mm256_mul_ps)

BAKE_SIMDF_OP_M_M(_op_simd_div_farr_farr, md_mm256_div_ps)
BAKE_SIMDF_OP_M_S(_op_simd_div_farr_f,    md_mm256_div_ps)

BAKE_SIMDF_OP_M(_op_simd_abs_farr,        md_mm256_abs_ps)

BAKE_SIMDF_OP_M_M(_op_simd_min_farr_farr, md_mm256_min_ps)
BAKE_SIMDF_OP_M_S(_op_simd_min_farr_f,    md_mm256_min_ps)
BAKE_SIMDF_OP_M_M(_op_simd_max_farr_farr, md_mm256_max_ps)
BAKE_SIMDF_OP_M_S(_op_simd_max_farr_f,    md_mm256_max_ps)

BAKE_OP_S_S(_op_add_f_f,        +, float)
BAKE_OP_S_S(_op_sub_f_f,        -, float)
BAKE_OP_S_S(_op_mul_f_f,        *, float)
BAKE_OP_S_S(_op_div_f_f,        /, float)
BAKE_OP_UNARY_S(_op_neg_f,      -, float)
BAKE_OP_UNARY_M(_op_neg_farr,   -, float)

BAKE_OP_M_S(_op_add_farr_f, +, float)
BAKE_OP_M_S(_op_sub_farr_f, -, float)
BAKE_OP_M_S(_op_mul_farr_f, *, float)
BAKE_OP_M_S(_op_div_farr_f, /, float)

BAKE_OP_M_M(_op_add_farr_farr, +, float)
BAKE_OP_M_M(_op_sub_farr_farr, -, float)
BAKE_OP_M_M(_op_mul_farr_farr, *, float)
BAKE_OP_M_M(_op_div_farr_farr, /, float)

// Forward declarations of functions
// @TODO: Add your declarations here

// Casts (Implicit ones)

static int _cast_int_to_flt             (data_t*, data_t[], eval_context_t*);
static int _cast_int_to_irng            (data_t*, data_t[], eval_context_t*);
static int _cast_irng_to_frng           (data_t*, data_t[], eval_context_t*);
static int _cast_int_arr_to_irng_arr    (data_t*, data_t[], eval_context_t*);
static int _cast_int_arr_to_flt_arr     (data_t*, data_t[], eval_context_t*);
static int _cast_irng_arr_to_frng_arr   (data_t*, data_t[], eval_context_t*);
static int _cast_int_arr_to_bf          (data_t*, data_t[], eval_context_t*);
static int _cast_irng_arr_to_bf         (data_t*, data_t[], eval_context_t*);
static int _join_bf_arr              (data_t*, data_t[], eval_context_t*);

// Basic operations
static int _min_farr  (data_t*, data_t[], eval_context_t*); // (float[]) -> float
static int _max_farr  (data_t*, data_t[], eval_context_t*); // (float[]) -> float

// Logical operators for custom types
static int _not  (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _and  (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _or   (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _xor  (data_t*, data_t[], eval_context_t*); // -> bitfield

// Selectors
// Atom level selectors
static int _all             (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_x        (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_y        (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_z        (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_xyz      (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_expl_flt (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_expl_frng(data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_impl_flt (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_impl_frng(data_t*, data_t[], eval_context_t*); // -> bitfield
static int _name            (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _element_str     (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _element_irng    (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _atom_irng       (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _atom_int        (data_t*, data_t[], eval_context_t*); // -> bitfield

static int _ion             (data_t*, data_t[], eval_context_t*); // -> bitfield

static int _ring            (data_t*, data_t[], eval_context_t*); // -> bitfield[]

// Residue level selectors
static int _water   (data_t*, data_t[], eval_context_t*);   // -> bitfield[]
static int _protein   (data_t*, data_t[], eval_context_t*); // -> bitfield[]
static int _nucleic   (data_t*, data_t[], eval_context_t*); // -> bitfield[]
static int _resname (data_t*, data_t[], eval_context_t*);   // (str[])          -> bitfield[]
static int _resid   (data_t*, data_t[], eval_context_t*);   // (int[]/irange[]) -> bitfield[]
static int _residue (data_t*, data_t[], eval_context_t*);   // (irange[])       -> bitfield[]

static int _fill_residue(data_t*, data_t[], eval_context_t*);   // (bitfield[]) -> bitfield
static int _fill_chain  (data_t*, data_t[], eval_context_t*);   // (bitfield[]) -> bitfield

// Chain level selectors
static int _chain_irng  (data_t*, data_t[], eval_context_t*); // (irange[]) -> bitfield
static int _chain_str   (data_t*, data_t[], eval_context_t*); // (str[]) -> bitfield

// Property Compute
static int _distance        (data_t*, data_t[], eval_context_t*); // (position[], position[]) -> float
static int _distance_min    (data_t*, data_t[], eval_context_t*); // (position[], position[]) -> float
static int _distance_max    (data_t*, data_t[], eval_context_t*); // (position[], position[]) -> float
static int _distance_pair   (data_t*, data_t[], eval_context_t*); // (position[N], position[M]) -> float[N*M]

static int _angle   (data_t*, data_t[], eval_context_t*); // (position[], position[], position[]) -> float
static int _dihedral(data_t*, data_t[], eval_context_t*); // (position[], position[], position[]), position[]) -> float

static int _rmsd    (data_t*, data_t[], eval_context_t*); // (bitfield) -> float

// Radial distribution function: The idea is that we use a fixed (high) amount of bins, then we let the user choose some kernel to smooth it.
static int _rdf_flt (data_t*, data_t[], eval_context_t*); // (position[], position[], float)  -> float[1024] (Histogram).
static int _rdf_frng(data_t*, data_t[], eval_context_t*); // (position[], position[], frange) -> float[1024] (Histogram).

static int _sdf     (data_t*, data_t[], eval_context_t*); // (bitfield, bitfield, float) -> float[128][128][128]. This one cannot be stored explicitly as one copy per frame, but is rather accumulated.

// Geometric operations
static int _com     (data_t*, data_t[], eval_context_t*);  // (position[]) -> float[3]
static int _plane   (data_t*, data_t[], eval_context_t*);  // (position[]) -> float[4]
static int _coordinate(data_t*, data_t[], eval_context_t*);  // (position[]) -> float[3][]
static int _coordinate_x(data_t*, data_t[], eval_context_t*);  // (position[]) -> float[]
static int _coordinate_y(data_t*, data_t[], eval_context_t*);  // (position[]) -> float[]
static int _coordinate_z(data_t*, data_t[], eval_context_t*);  // (position[]) -> float[]

static int _coordinate_xy(data_t*, data_t[], eval_context_t*);  // (position[]) -> float2[]
static int _coordinate_xz(data_t*, data_t[], eval_context_t*);  // (position[]) -> float2[]

static int _coordinate_yz(data_t*, data_t[], eval_context_t*);  // (position[]) -> float2[]

static int _shape_weights(data_t*, data_t[], eval_context_t*); // (position[]) -> float

// Linear algebra
static int _dot           (data_t*, data_t[], eval_context_t*); // (float[], float[]) -> float
static int _cross         (data_t*, data_t[], eval_context_t*); // (float[3], float[3]) -> float[3]
static int _length        (data_t*, data_t[], eval_context_t*); // (float[]) -> float
static int _mat4_mul_mat4 (data_t*, data_t[], eval_context_t*); // (float[4][4], float[4][4]) -> float[4][4]
static int _mat4_mul_vec4 (data_t*, data_t[], eval_context_t*); // (float[4][4], float[4]) -> float[4]

static int _vec2 (data_t*, data_t[], eval_context_t*); // (float, float) -> float[2]
static int _vec3 (data_t*, data_t[], eval_context_t*); // (float, float, float) -> float[4]
static int _vec4 (data_t*, data_t[], eval_context_t*); // (float, float, float, float) -> float[4]

// Misc
// Count the number of bits set (popcount)
static int _count (data_t*, data_t[], eval_context_t*); // (bitfield) -> float

// Specialized count function with supplied argument
static int _count_with_arg(data_t* , data_t[], eval_context_t*); // (bitfield) -> float

static int _op_simd_neg_farr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)ctx;
    ASSERT(dst);
    int64_t total_count = type_info_element_stride_count(arg[0].type);
    ASSERT(total_count % 8 == 0);
    const float* src_arr = as_float_arr(arg[0]);
    float* dst_arr = as_float_arr(*dst);

    for (int64_t i = 0; i < total_count; i += 8) {
        md_256 val = md_mm256_loadu_ps(src_arr + i);
        md_mm256_storeu_ps(dst_arr + i, md_mm256_sub_ps(md_mm256_setzero_ps(), val));
    }
    return 0;
}

// Predefined constants
#define _PI  3.14159265358f
#define _TAU 6.28318530718f
#define _E   2.71828182845f

#define CSTR(cstr) {cstr"", sizeof(cstr)-1}

// @TODO: Add your values here
static constant_t constants[] = {
    {CSTR("PI"),     _PI,  {.base = {.dim.angle = 1 }, .mult = 1.0}},
    {CSTR("TAU"),    _TAU, {.base = {.dim.angle = 1 }, .mult = 1.0}},
    {CSTR("E"),      _E},
};

// IMPLICIT CASTS/CONVERSIONS
static procedure_t casts[] = {
    {CSTR("cast"),    TI_FLOAT,         1,  {TI_INT},           _cast_int_to_flt},
    {CSTR("cast"),    TI_IRANGE,        1,  {TI_INT},           _cast_int_to_irng},
    {CSTR("cast"),    TI_FRANGE,        1,  {TI_IRANGE},        _cast_irng_to_frng},
    {CSTR("cast"),    TI_IRANGE_ARR,    1,  {TI_INT_ARR},       _cast_int_arr_to_irng_arr,  FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("cast"),    TI_FLOAT_ARR,     1,  {TI_INT_ARR},       _cast_int_arr_to_flt_arr,   FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("cast"),    TI_FRANGE_ARR,    1,  {TI_IRANGE_ARR},    _cast_irng_arr_to_frng_arr, FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("cast"),    TI_BITFIELD,      1,  {TI_INT_ARR},       _cast_int_arr_to_bf},
    {CSTR("cast"),    TI_BITFIELD,      1,  {TI_IRANGE_ARR},    _cast_irng_arr_to_bf},
    {CSTR("cast"),    TI_BITFIELD,      1,  {TI_BITFIELD_ARR},  _join_bf_arr},

    //{CSTR("cast"),    TI_FLOAT3_ARR,    1,  {TI_COORDINATE_ARR},  _coordinate, FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_DYNAMIC_LENGTH},
};

static procedure_t operators[] = {
    {CSTR("not"),    TI_BOOL,           1,  {TI_BOOL},                  _op_not_b},
    {CSTR("or"),     TI_BOOL,           2,  {TI_BOOL,   TI_BOOL},       _op_or_b_b},
    {CSTR("xor"),    TI_BOOL,           2,  {TI_BOOL,   TI_BOOL},       _op_xor_b_b},
    {CSTR("and"),    TI_BOOL,           2,  {TI_BOOL,   TI_BOOL},       _op_and_b_b},

    // BITFIELD NOT
    {CSTR("not"),    TI_BITFIELD,       1,  {TI_BITFIELD_ARR},                  _not,   FLAG_FLATTEN},
    {CSTR("and"),    TI_BITFIELD,       2,  {TI_BITFIELD_ARR, TI_BITFIELD_ARR}, _and,   FLAG_FLATTEN},
    {CSTR("or"),     TI_BITFIELD,       2,  {TI_BITFIELD_ARR, TI_BITFIELD_ARR}, _or,    FLAG_FLATTEN},
    {CSTR("xor"),    TI_BITFIELD,       2,  {TI_BITFIELD_ARR, TI_BITFIELD_ARR}, _xor,   FLAG_FLATTEN},

    // Binary add
    {CSTR("+"),      TI_FLOAT,          2,  {TI_FLOAT,      TI_FLOAT},              _op_add_f_f},
    {CSTR("+"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT},              _op_add_farr_f,         FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("+"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},          _op_add_farr_farr,      FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {CSTR("+"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_DISTRIBUTION},    _op_simd_add_farr_farr, FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("+"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_FLOAT},           _op_simd_add_farr_f,    FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("+"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_VOLUME},                _op_simd_add_farr_farr, FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("+"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_FLOAT},                 _op_simd_add_farr_f,    FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("+"),      TI_INT,            2,  {TI_INT,        TI_INT},                _op_add_i_i},
    {CSTR("+"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT},                _op_add_iarr_i,         FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("+"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT_ARR},            _op_add_iarr_iarr,      FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    // Unary negation
    {CSTR("-"),      TI_FLOAT,          1,  {TI_FLOAT},                             _op_neg_f},
    {CSTR("-"),      TI_FLOAT_ARR,      1,  {TI_FLOAT_ARR},                         _op_neg_farr,           FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("-"),      TI_INT,            1,  {TI_INT},                               _op_neg_i},
    {CSTR("-"),      TI_INT_ARR,        1,  {TI_INT_ARR},                           _op_neg_iarr,           FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("-"),      TI_DISTRIBUTION,   1,  {TI_DISTRIBUTION},                      _op_simd_neg_farr,      FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("-"),      TI_VOLUME,         1,  {TI_VOLUME},                            _op_simd_neg_farr,      FLAG_RET_AND_ARG_EQUAL_LENGTH},

    // Binary sub
    {CSTR("-"),      TI_FLOAT,          2,  {TI_FLOAT,      TI_FLOAT},              _op_sub_f_f},
    {CSTR("-"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT},              _op_sub_farr_f,         FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("-"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},          _op_sub_farr_farr,      FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {CSTR("-"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_DISTRIBUTION},    _op_simd_sub_farr_farr, FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("-"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_FLOAT},           _op_simd_sub_farr_f,    FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("-"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_VOLUME},                _op_simd_sub_farr_farr, FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("-"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_FLOAT},                 _op_simd_sub_farr_f,    FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("-"),      TI_INT,            2,  {TI_INT,        TI_INT},                _op_sub_i_i},
    {CSTR("-"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT},                _op_sub_iarr_i,         FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("-"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT_ARR},            _op_sub_iarr_iarr,      FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    // Binary mul
    {CSTR("*"),      TI_FLOAT,          2,  {TI_FLOAT,      TI_FLOAT},              _op_mul_f_f},
    {CSTR("*"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT},              _op_mul_farr_f,         FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("*"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},          _op_mul_farr_farr,      FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {CSTR("*"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_DISTRIBUTION},    _op_simd_mul_farr_farr, FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("*"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_FLOAT},           _op_simd_mul_farr_f,    FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("*"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_VOLUME},                _op_simd_mul_farr_farr, FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("*"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_FLOAT},                 _op_simd_mul_farr_f,    FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("*"),      TI_INT,            2,  {TI_INT,        TI_INT},                _op_mul_i_i},
    {CSTR("*"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT},                _op_mul_iarr_i,         FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("*"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT_ARR},            _op_mul_iarr_iarr,      FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    // Binary div
    {CSTR("/"),      TI_FLOAT,          2,  {TI_FLOAT,      TI_FLOAT},              _op_div_f_f},
    {CSTR("/"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT},              _op_div_farr_f,         FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("/"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},          _op_div_farr_farr,      FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {CSTR("/"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_DISTRIBUTION},    _op_simd_div_farr_farr, FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("/"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_FLOAT},           _op_simd_div_farr_f,    FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    
    {CSTR("/"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_VOLUME},                _op_simd_div_farr_farr, FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("/"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_FLOAT},                 _op_simd_div_farr_f,    FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("/"),      TI_INT,            2,  {TI_INT,        TI_INT},                _op_div_i_i},
    {CSTR("/"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT},                _op_div_iarr_i,         FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("/"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT_ARR},            _op_div_iarr_iarr,      FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
};

static procedure_t procedures[] = {
    // NATIVE FUNCS
    {CSTR("sqrt"),   TI_FLOAT, 1, {TI_FLOAT}, _sqrtf},
    {CSTR("cbrt"),   TI_FLOAT, 1, {TI_FLOAT}, _cbrtf},
    {CSTR("abs"),    TI_FLOAT, 1, {TI_FLOAT}, _fabsf},
    {CSTR("floor"),  TI_FLOAT, 1, {TI_FLOAT}, _floorf},
    {CSTR("ceil"),   TI_FLOAT, 1, {TI_FLOAT}, _ceilf},
    {CSTR("cos"),    TI_FLOAT, 1, {TI_FLOAT}, _cosf},
    {CSTR("sin"),    TI_FLOAT, 1, {TI_FLOAT}, _sinf},
    {CSTR("asin"),   TI_FLOAT, 1, {TI_FLOAT}, _asinf},
    {CSTR("acos"),   TI_FLOAT, 1, {TI_FLOAT}, _acosf},
    {CSTR("atan"),   TI_FLOAT, 1, {TI_FLOAT}, _atanf},
    {CSTR("log"),    TI_FLOAT, 1, {TI_FLOAT}, _logf},
    {CSTR("exp"),    TI_FLOAT, 1, {TI_FLOAT}, _expf},
    {CSTR("log2"),   TI_FLOAT, 1, {TI_FLOAT}, _log2f},
    {CSTR("exp2"),   TI_FLOAT, 1, {TI_FLOAT}, _exp2f},
    {CSTR("log10"),  TI_FLOAT, 1, {TI_FLOAT}, _log10f},

    {CSTR("atan"),   TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _atan2f},
    {CSTR("atan2"),  TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _atan2f},
    {CSTR("pow"),    TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _powf},
    {CSTR("min"),    TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _fminf},
    {CSTR("max"),    TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _fmaxf},

    // ARRAY VERSIONS OF NATIVE FUNCS
    {CSTR("abs"),    TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_fabsf,    FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("floor"),  TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_floorf,   FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("ceil"),   TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_ceilf,    FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("min"),    TI_FLOAT,     1, {TI_FLOAT_ARR}, _min_farr},
    {CSTR("max"),    TI_FLOAT,     1, {TI_FLOAT_ARR}, _max_farr},

    // VECTORIZED VERSIONS FOR ARRAYS OF BIGGER DATA
    {CSTR("abs"),    TI_VOLUME, 1,  {TI_VOLUME}, _op_simd_abs_farr},
    {CSTR("min"),    TI_VOLUME, 2,  {TI_VOLUME, TI_FLOAT},  _op_simd_min_farr_f,    FLAG_SYMMETRIC_ARGS},
    {CSTR("min"),    TI_VOLUME, 2,  {TI_VOLUME, TI_VOLUME}, _op_simd_min_farr_farr},
    {CSTR("max"),    TI_VOLUME, 2,  {TI_VOLUME, TI_FLOAT},  _op_simd_max_farr_f,    FLAG_SYMMETRIC_ARGS},
    {CSTR("max"),    TI_VOLUME, 2,  {TI_VOLUME, TI_VOLUME}, _op_simd_max_farr_farr},

    // LINEAR ALGEBRA
    {CSTR("dot"),    TI_FLOAT,   2, {TI_FLOAT_ARR,   TI_FLOAT_ARR},  _dot},
    {CSTR("cross"),  TI_FLOAT3,  2, {TI_FLOAT3,      TI_FLOAT3},     _cross},
    {CSTR("length"), TI_FLOAT,   1, {TI_FLOAT_ARR},                  _length},
    {CSTR("mul"),    TI_FLOAT44, 2, {TI_FLOAT44,     TI_FLOAT44},    _mat4_mul_mat4},
    {CSTR("mul"),    TI_FLOAT4,  2, {TI_FLOAT44,     TI_FLOAT4},     _mat4_mul_vec4},

    // CONSTRUCTORS
    {CSTR("vec2"),   TI_FLOAT2,  2, {TI_FLOAT, TI_FLOAT},                       _vec2},
    {CSTR("vec3"),   TI_FLOAT3,  3, {TI_FLOAT, TI_FLOAT, TI_FLOAT},             _vec3},
    {CSTR("vec4"),   TI_FLOAT4,  4, {TI_FLOAT, TI_FLOAT, TI_FLOAT, TI_FLOAT},   _vec4},

    // --- SELECTORS ---

    // Atom level
    {CSTR("all"),       TI_BITFIELD, 0, {0},                _all},
    {CSTR("ion"),       TI_BITFIELD, 0, {0},                _ion},
    {CSTR("type"),      TI_BITFIELD, 1, {TI_STRING_ARR},    _name,              FLAG_STATIC_VALIDATION},
    {CSTR("name"),      TI_BITFIELD, 1, {TI_STRING_ARR},    _name,              FLAG_STATIC_VALIDATION},
    {CSTR("label"),     TI_BITFIELD, 1, {TI_STRING_ARR},    _name,              FLAG_STATIC_VALIDATION},
    {CSTR("element"),   TI_BITFIELD, 1, {TI_STRING_ARR},    _element_str,       FLAG_STATIC_VALIDATION},
    {CSTR("element"),   TI_BITFIELD, 1, {TI_IRANGE_ARR},    _element_irng,      FLAG_STATIC_VALIDATION},
    {CSTR("atom"),      TI_BITFIELD, 1, {TI_IRANGE_ARR},    _atom_irng,         FLAG_STATIC_VALIDATION},
    {CSTR("atom"),      TI_BITFIELD, 1, {TI_INT_ARR},       _atom_int,          FLAG_STATIC_VALIDATION},
    
    {CSTR("ring"),      TI_BITFIELD_ARR, 0, {0},            _ring,              FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},

    // Residue level
    {CSTR("protein"),   TI_BITFIELD_ARR, 0, {0},                _protein,       FLAG_QUERYABLE_LENGTH},
    {CSTR("nucleic"),   TI_BITFIELD_ARR, 0, {0},                _nucleic,       FLAG_QUERYABLE_LENGTH},
    {CSTR("water"),     TI_BITFIELD_ARR, 0, {0},                _water,         FLAG_QUERYABLE_LENGTH},
    {CSTR("resname"),   TI_BITFIELD_ARR, 1, {TI_STRING_ARR},    _resname,       FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {CSTR("residue"),   TI_BITFIELD_ARR, 1, {TI_STRING_ARR},    _resname,       FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {CSTR("resid"),     TI_BITFIELD_ARR, 1, {TI_IRANGE_ARR},    _resid,         FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {CSTR("residue"),   TI_BITFIELD_ARR, 1, {TI_IRANGE_ARR},    _residue,       FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},

    // Chain level
    {CSTR("chain"),     TI_BITFIELD_ARR,  1,  {TI_STRING_ARR},    _chain_str,   FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {CSTR("chain"),     TI_BITFIELD_ARR,  1,  {TI_IRANGE_ARR},    _chain_irng,  FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},

    // Dynamic selectors (depend on atomic position, therefore marked as dynamic which means the values cannot be determined at compile-time)
    // Also have variable result (well its a single bitfield, but the number of atoms within is not fixed)
    {CSTR("within_x"),      TI_BITFIELD, 1, {TI_FRANGE},    _within_x,     FLAG_DYNAMIC},
    {CSTR("within_y"),      TI_BITFIELD, 1, {TI_FRANGE},    _within_y,     FLAG_DYNAMIC},
    {CSTR("within_z"),      TI_BITFIELD, 1, {TI_FRANGE},    _within_z,     FLAG_DYNAMIC},
    {CSTR("within_xyz"),    TI_BITFIELD, 3, {TI_FRANGE, TI_FRANGE, TI_FRANGE}, _within_xyz,     FLAG_DYNAMIC},

    {CSTR("within"),    TI_BITFIELD, 2, {TI_BITFIELD_ARR,   TI_FLOAT},  _within_expl_flt,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_FLATTEN | FLAG_SPATIAL_QUERY | FLAG_SYMMETRIC_ARGS },
    {CSTR("within"),    TI_BITFIELD, 2, {TI_FLOAT3_ARR,     TI_FLOAT},  _within_expl_flt,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_FLATTEN | FLAG_SPATIAL_QUERY | FLAG_SYMMETRIC_ARGS },
    {CSTR("within"),    TI_BITFIELD, 2, {TI_BITFIELD_ARR,   TI_FRANGE}, _within_expl_frng,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_FLATTEN | FLAG_SPATIAL_QUERY | FLAG_SYMMETRIC_ARGS },
    {CSTR("within"),    TI_BITFIELD, 2, {TI_FLOAT3_ARR,     TI_FRANGE}, _within_expl_frng,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_FLATTEN | FLAG_SPATIAL_QUERY | FLAG_SYMMETRIC_ARGS },

    {CSTR("within"),    TI_BITFIELD, 1, {TI_FLOAT},                     _within_impl_flt,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_FLATTEN | FLAG_SPATIAL_QUERY},
    {CSTR("within"),    TI_BITFIELD, 1, {TI_FRANGE},                    _within_impl_frng,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_FLATTEN | FLAG_SPATIAL_QUERY},

    // --- PROPERTY COMPUTE ---
    {CSTR("distance"),          TI_FLOAT,       2,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR}, _distance,          FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_NO_FLATTEN},
    {CSTR("distance_min"),      TI_FLOAT,       2,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR}, _distance_min,      FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_NO_FLATTEN},
    {CSTR("distance_max"),      TI_FLOAT,       2,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR}, _distance_max,      FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_NO_FLATTEN},
    {CSTR("distance_pair"),     TI_FLOAT_ARR,   2,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR}, _distance_pair,     FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_QUERYABLE_LENGTH | FLAG_NO_FLATTEN},

    {CSTR("angle"),     TI_FLOAT,   3,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR, TI_COORDINATE_ARR},                    _angle,     FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_NO_FLATTEN},
    {CSTR("dihedral"),  TI_FLOAT,   4,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR, TI_COORDINATE_ARR, TI_COORDINATE_ARR},   _dihedral,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_NO_FLATTEN},

    {CSTR("rmsd"),      TI_FLOAT,   1,  {TI_BITFIELD},    _rmsd,     FLAG_DYNAMIC | FLAG_VISUALIZE},

    {CSTR("rdf"), TI_DISTRIBUTION, 3, {TI_COORDINATE_ARR, TI_COORDINATE_ARR, TI_FLOAT},  _rdf_flt,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_NO_FLATTEN},
    {CSTR("rdf"), TI_DISTRIBUTION, 3, {TI_COORDINATE_ARR, TI_COORDINATE_ARR, TI_FRANGE}, _rdf_frng,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_NO_FLATTEN},

    {CSTR("sdf"),       TI_VOLUME, 3, {TI_BITFIELD_ARR, TI_BITFIELD, TI_FLOAT},      _sdf,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_SDF | FLAG_VISUALIZE | FLAG_NO_FLATTEN},

    // --- GEOMETRICAL OPERATIONS ---
    {CSTR("com"),           TI_FLOAT3,   1,  {TI_COORDINATE_ARR},  _com,           FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_NO_FLATTEN},
    {CSTR("plane"),         TI_FLOAT4,   1,  {TI_COORDINATE_ARR},  _plane,         FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_NO_FLATTEN},
    {CSTR("shape_weights"), TI_FLOAT3,   1,  {TI_COORDINATE_ARR},  _shape_weights, FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_FLATTEN},

    {CSTR("coord"),      TI_FLOAT3_ARR,  1,  {TI_COORDINATE_ARR},  _coordinate,      FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE | FLAG_NO_FLATTEN},
    {CSTR("coord_x"),    TI_FLOAT_ARR,   1,  {TI_COORDINATE_ARR},  _coordinate_x,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE | FLAG_NO_FLATTEN},
    {CSTR("coord_y"),    TI_FLOAT_ARR,   1,  {TI_COORDINATE_ARR},  _coordinate_y,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE | FLAG_NO_FLATTEN},
    {CSTR("coord_z"),    TI_FLOAT_ARR,   1,  {TI_COORDINATE_ARR},  _coordinate_z,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE | FLAG_NO_FLATTEN},

    {CSTR("coord_xy"),   TI_FLOAT2_ARR,  1,  {TI_COORDINATE_ARR},  _coordinate_xy,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE | FLAG_NO_FLATTEN},
    {CSTR("coord_xz"),   TI_FLOAT2_ARR,  1,  {TI_COORDINATE_ARR},  _coordinate_xz,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE | FLAG_NO_FLATTEN},

    {CSTR("coord_yz"),   TI_FLOAT2_ARR,  1,  {TI_COORDINATE_ARR},  _coordinate_yz,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE | FLAG_NO_FLATTEN},


    // --- MISC ---
    {CSTR("count"),     TI_FLOAT,           1,  {TI_BITFIELD},              _count},
    {CSTR("count"),     TI_FLOAT,           2,  {TI_BITFIELD, TI_STRING},   _count_with_arg, FLAG_STATIC_VALIDATION},

    {CSTR("join"),      TI_BITFIELD,        1,  {TI_BITFIELD_ARR},  _join_bf_arr,   FLAG_FLATTEN},
    {CSTR("flatten"),   TI_BITFIELD,        1,  {TI_BITFIELD_ARR},  _join_bf_arr,   FLAG_FLATTEN},

    {CSTR("residue"),   TI_BITFIELD_ARR,    1,  {TI_BITFIELD_ARR},  _fill_residue,  FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_NO_FLATTEN},
    {CSTR("chain"),     TI_BITFIELD_ARR,    1,  {TI_BITFIELD_ARR},  _fill_chain,    FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_NO_FLATTEN},

};

#undef CSTR

static inline md_spatial_hash_t* get_spatial_hash(eval_context_t* ctx) {
	ASSERT(ctx);
    if (!ctx->spatial_hash) {
        // Lazily generate the spatial hash if needed
        ctx->spatial_hash = md_spatial_hash_create_soa(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, 0, ctx->mol->atom.count, &ctx->mol->unit_cell, ctx->temp_alloc);
    }
	return ctx->spatial_hash;
}

static inline void visualize_atom_mask(const md_bitfield_t* mask, eval_context_t* ctx) {
    ASSERT(ctx->vis);
    if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        md_bitfield_or_inplace(ctx->vis_structure, mask);
    }
}

static inline void visualize_atom_range(irange_t range, eval_context_t* ctx) {
    ASSERT(ctx->vis);
    if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        md_bitfield_set_range(ctx->vis_structure, range.beg, range.end);
    }
}

static inline void visualize_atom_index(int64_t index, eval_context_t* ctx) {
    ASSERT(ctx->vis);
    if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        md_bitfield_set_bit(ctx->vis_structure, index);
    }
}

/*
static inline void visualize_atom_indices32(const int32_t* indices, int64_t num_indices, eval_context_t* ctx) {
    ASSERT(ctx->vis);
    if (ctx->vis_flag & MD_SCRIPT_VISUALIZE_ATOMS) {
        for (int64_t i = 0; i < num_indices; ++i) {
            md_bitfield_set_bit(vis->atom_mask, indices[i]);
        }
    }
}

static inline void visualize_atom_indices64(const int64_t* indices, int64_t num_indices, eval_context_t* ctx) {
    ASSERT(ctx->vis);
    if (ctx->vis_flag & MD_SCRIPT_VISUALIZE_ATOMS) {
        for (int64_t i = 0; i < num_indices; ++i) {
            md_bitfield_set_bit(ctx->vis->atom_mask, indices[i]);
        }
    }
}
*/

#define COLOR_WHITE 0xFFFFFFFFU

static inline md_script_vis_vertex_t vertex(vec3_t p, uint32_t color) {
    return (md_script_vis_vertex_t){p, color};
}

static inline void push_point(md_script_vis_vertex_t v, md_script_vis_t* vis) {
    ASSERT(vis);
    ASSERT(vis->alloc);
    md_array_push(vis->points, v, vis->alloc);
}

static inline void push_line(md_script_vis_vertex_t v0, md_script_vis_vertex_t v1, md_script_vis_t* vis) {
    ASSERT(vis);
    ASSERT(vis->alloc);
    md_array_push(vis->lines, v0, vis->alloc);
    md_array_push(vis->lines, v1, vis->alloc);
}

static inline void push_triangle(md_script_vis_vertex_t v0, md_script_vis_vertex_t v1, md_script_vis_vertex_t v2, md_script_vis_t* vis) {
    ASSERT(vis);
    ASSERT(vis->alloc);
    md_array_push(vis->triangles, v0, vis->alloc);
    md_array_push(vis->triangles, v1, vis->alloc);
    md_array_push(vis->triangles, v2, vis->alloc);
}

static inline void push_sphere(vec3_t xyz, float r, uint32_t color, md_script_vis_t* vis) {
    ASSERT(vis);
    ASSERT(vis->alloc);
    md_script_vis_sphere_t sphere = {xyz, r, color};
    md_array_push(vis->spheres, sphere, vis->alloc);
}

static inline md_bitfield_t* push_structure(md_script_vis_t* vis) {
    ASSERT(vis);
    ASSERT(vis->alloc);
    md_bitfield_t* bf = md_array_push(vis->structures, (md_bitfield_t){0}, vis->alloc);
    md_bitfield_init(bf, vis->alloc);
    return bf;
}

static inline bool idx_in_range(int idx, irange_t range) {
    // I think a range should be inclusive in this context... Since we should be 1 based and not 0 based on indices
    return range.beg <= idx && idx <= range.end;
}

static inline bool range_in_range(irange_t small_range, irange_t big_range) {
    return big_range.beg <= small_range.beg && small_range.end <= big_range.end;
}

static inline vec3_t extract_com(const float* x, const float* y, const float* z, const float* w, const md_bitfield_t* bitfield) {
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    ASSERT(bitfield);

    vec4_t sum = {0};
    md_bitfield_iter_t it = md_bitfield_iter_create(bitfield);
    while (md_bitfield_iter_next(&it)) {
        const int64_t idx = md_bitfield_iter_idx(&it);
        const float weight = w ? w[idx] : 1.0f;
        vec4_t v = vec4_set(x[idx], y[idx], z[idx], 1.0f);
        v = vec4_mul_f(v, weight);
        sum = vec4_add(sum, v);
    }
    if (sum.w == 0) sum.w = 1.0f;
    return vec3_div_f(vec3_from_vec4(sum), sum.w);
}

static inline vec3_t* extract_vec3(const float* x, const float* y, const float* z, const md_bitfield_t* bf, md_allocator_i* alloc) {
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    ASSERT(bf);
    ASSERT(alloc);

    md_array(vec3_t) result = md_array_create(vec3_t, md_bitfield_popcount(bf), alloc);
    md_bitfield_iter_t it = md_bitfield_iter_create(bf);
    size_t i = 0;
    while (md_bitfield_iter_next(&it)) {
        const uint64_t idx = md_bitfield_iter_idx(&it);
        result[i++] = vec3_set(x[idx], y[idx], z[idx]);
    }
    return result;
}

static inline int64_t extract_xyz(float* dst_x, float* dst_y, float* dst_z, const float* src_x, const float* src_y, const float* src_z, const md_bitfield_t* bf) {
    ASSERT(dst_x);
    ASSERT(dst_y);
    ASSERT(dst_z);
    ASSERT(src_x);
    ASSERT(src_y);
    ASSERT(src_z);
    ASSERT(bf);

    int64_t count = 0;

    md_bitfield_iter_t it = md_bitfield_iter_create(bf);
    while (md_bitfield_iter_next(&it)) {
        const int64_t idx = md_bitfield_iter_idx(&it);
        dst_x[count] = src_x[idx];
        dst_y[count] = src_y[idx];
        dst_z[count] = src_z[idx];
        count += 1;
    }

    return count;
}

static inline int64_t extract_xyz_vec3(vec3_t* dst_xyz, const float* src_x, const float* src_y, const float* src_z, const md_bitfield_t* bf) {
    ASSERT(dst_xyz);
    ASSERT(src_x);
    ASSERT(src_y);
    ASSERT(src_z);
    ASSERT(bf);

    int64_t count = 0;

    md_bitfield_iter_t it = md_bitfield_iter_create(bf);
    while (md_bitfield_iter_next(&it)) {
        const int64_t idx = md_bitfield_iter_idx(&it);
        dst_xyz[count] = vec3_set(src_x[idx], src_y[idx], src_z[idx]);
        count += 1;
    }

    return count;
}

static inline int64_t extract_xyzw(float* dst_x, float* dst_y, float* dst_z, float* dst_w, const float* src_x, const float* src_y, const float* src_z, const float* src_w, const md_bitfield_t* bf) {
    ASSERT(dst_x);
    ASSERT(dst_y);
    ASSERT(dst_z);
    ASSERT(dst_w);
    ASSERT(src_x);
    ASSERT(src_y);
    ASSERT(src_z);
    ASSERT(src_w);
    ASSERT(bf);

    int64_t count = 0;

    int64_t beg_bit = bf->beg_bit;
    int64_t end_bit = bf->end_bit;
    while ((beg_bit = md_bitfield_scan(bf, beg_bit, end_bit)) != 0) {
        const int64_t idx = beg_bit - 1;
        dst_x[count] = src_x[idx];
        dst_y[count] = src_y[idx];
        dst_z[count] = src_z[idx];
        dst_w[count] = src_w[idx];
        count += 1;
    }

    return count;
}

static inline int64_t extract_xyzw_vec4(vec4_t* dst_xyzw, const float* src_x, const float* src_y, const float* src_z, const float* src_w, const md_bitfield_t* bf) {
    ASSERT(dst_xyzw);
    ASSERT(src_x);
    ASSERT(src_y);
    ASSERT(src_z);
    ASSERT(src_w);
    ASSERT(bf);

    int64_t count = 0;
    
    md_bitfield_iter_t it = md_bitfield_iter_create(bf);
    while (md_bitfield_iter_next(&it)) {
        const int64_t idx = md_bitfield_iter_idx(&it);
        dst_xyzw[count] = vec4_set(src_x[idx], src_y[idx], src_z[idx], src_w[idx]);
        count += 1;
    }

    return count;
}


// @TODO: All these have poor names

static inline irange_t remap_range_to_context(irange_t range, irange_t context) {
    if (range.beg == INT32_MIN)
        range.beg = 0;
    else
        range.beg = range.beg - 1;

    if (range.end == INT32_MAX)
        range.end = context.end - context.beg;
    else
        range.end = range.end;

    range.beg = context.beg + range.beg;
    range.end = context.beg + range.end;

    return range;
}

static inline int32_t remap_index_to_context(int32_t index, irange_t context) {
    return context.beg + index - 1;
}

static inline irange_t clamp_range(irange_t range, irange_t context) {
    range.beg = CLAMP(range.beg, context.beg, context.end);
    range.end = CLAMP(range.end, context.beg, context.end);
    return range;
}

static int32_t find_val(const int32_t* arr, int64_t count, int32_t val) {
    for (int32_t i = 0; i < (int32_t)count; ++i) {
        if (arr[i] == val) {
			return i;
		}
    }
    return -1;
}

static int32_t find_label(const md_label_t* arr, int64_t count, str_t lbl) {
    for (int32_t i = 0; i < (int32_t)count; ++i) {
		if (str_eq(LBL_TO_STR(arr[i]), lbl)) {
			return i;
		}
	}
	return -1;
}

static md_array(int32_t) get_residue_indices_in_context(const md_molecule_t* mol, const md_bitfield_t* bitfield, md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);

    md_array(int32_t) arr = 0;

    if (mol->residue.count) {
        if (bitfield) {
            for (size_t i = 0; i < mol->residue.count; ++i) {
                md_range_t res = md_residue_atom_range(mol->residue, i);
                if (md_bitfield_popcount_range(bitfield, res.beg, res.end)) {
                    md_array_push(arr, (int32_t)i, alloc);
                }
            }
        }
        else {
            for (int32_t i = 0; i < (int32_t)mol->residue.count; ++i) {
                md_array_push(arr, i, alloc);
            }
        }
    }

    return arr;
}

static md_array(int32_t) get_chain_indices_in_context(const md_molecule_t* mol, const md_bitfield_t* bitfield, md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);

    md_array(int32_t) arr = 0;
    if (mol->chain.count) {
        if (bitfield) {
            for (size_t i = 0; i < mol->chain.count; ++i) {
                md_range_t range = md_chain_atom_range(mol->chain, i);
                if (md_bitfield_popcount_range(bitfield, range.beg, range.end)) {
                    md_array_push(arr, (int32_t)i, alloc);
                }
            }
        }
        else {
            for (size_t i = 0; i < mol->chain.count; ++i) {
                md_array_push(arr, (int32_t)i, alloc);
            }
        }
    }

    return arr;
}

static inline irange_t get_atom_range_in_context(const md_molecule_t* mol, const md_bitfield_t* mol_ctx) {
    ASSERT(mol);
    irange_t range = { 0, (int32_t)mol->atom.count };
    if (mol_ctx) {
        range.beg = mol_ctx->beg_bit;
        range.end = mol_ctx->end_bit;
    }
    return range;
}

static inline irange_t get_residue_range_in_context(const md_molecule_t* mol, const md_bitfield_t* mol_ctx) {
    ASSERT(mol);
    irange_t range = { 0, (int32_t)mol->residue.count };
    if (mol_ctx) {
        range.beg = INT_MAX;
        range.end = 0;
        for (size_t i = 0; i < mol->residue.count; ++i) {
            md_range_t res = md_residue_atom_range(mol->residue, i);
			if (md_bitfield_popcount_range(mol_ctx, res.beg, res.end)) {
				range.beg = MIN(range.beg, (int32_t)i);
                range.end = MAX(range.end, (int32_t)i + 1);
				break;
			}
        }
    }
    return range;
}

static inline irange_t get_chain_range_in_context(const md_molecule_t* mol, const md_bitfield_t* mol_ctx) {
    ASSERT(mol);
    irange_t range = { 0, (int32_t)mol->chain.count };
    if (mol_ctx) {
        range.beg = INT_MAX;
        range.end = 0;
        for (size_t i = 0; i < mol->chain.count; ++i) {
            md_range_t atom_range = md_chain_atom_range(mol->chain, i);
            if (md_bitfield_popcount_range(mol_ctx, atom_range.beg, atom_range.end)) {
                range.beg = MIN(range.beg, (int32_t)i);
                range.end = MAX(range.end, (int32_t)i + 1);
                break;
            }
        }
    }
    return range;
}

static int coordinate_validate(data_t arg, int arg_idx, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_COORDINATE_ARR));
    if (element_count(arg) == 0) return 0;
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);
    const int64_t ctx_size = ctx_range.end - ctx_range.beg;

    switch (arg.type.base_type) {
    case TYPE_FLOAT:
        ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_FLOAT3_ARR));
        return (int)type_info_array_len(arg.type);
    case TYPE_INT: {
        int* indices = as_int_arr(arg);
        int count = 0;
        for (size_t i = 0; i < element_count(arg); ++i) {
            const int64_t idx = (int64_t)ctx_range.beg + (int64_t)indices[i] - 1;

            if (!(ctx_range.beg <= idx && idx < ctx_range.end)) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[arg_idx], "supplied index (%i) is not within expected range (%i:%i)",
                    (int)idx, 1, ctx_size);
                return STATIC_VALIDATION_ERROR;
            }
            if (ctx->mol_ctx) {
                if (!md_bitfield_test_bit(ctx->mol_ctx, idx)) {
                    LOG_ERROR(ctx->ir, ctx->arg_tokens[arg_idx], "supplied index (%i) is not represented within the supplied context", (int)idx);
                    return STATIC_VALIDATION_ERROR;
                }
            }
            count += 1;
        }
        return count;
    }
    case TYPE_BITFIELD:
    {
        if (ctx->backchannel && ctx->arg_flags[arg_idx] & FLAG_DYNAMIC) {
            ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
        }
        if (element_count(arg) == 1) {
            return (int)md_bitfield_popcount(as_bitfield(arg));
        }
        else {
            return (int)type_info_array_len(arg.type);
        }
    }
    case TYPE_IRANGE:
    {
        if (ctx->backchannel && ctx->arg_flags[arg_idx] & FLAG_DYNAMIC) {
            ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
        }
        irange_t* ranges = as_irange_arr(arg);
        int count = 0;
        for (size_t i = 0; i < element_count(arg); ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            if (range_in_range(range, ctx_range)) {
                if (ctx->mol_ctx) {
                    count += (int)md_bitfield_popcount_range(ctx->mol_ctx, range.beg, range.end);
                }
                else {
                    count += range.end - range.beg;
                }
            }
            else {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "supplied range (%i:%i) is not within expected range (%i:%i)",
                    ranges[i].beg, ranges[i].end, 1, ctx_size);
                return STATIC_VALIDATION_ERROR;
            }
        }
        return count;
    }
    default:
        ASSERT(false);
        return STATIC_VALIDATION_ERROR;
    }
}

static void coordinate_visualize(data_t arg, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(ctx->vis);

    if (element_count(arg) == 0) return;
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    switch (arg.type.base_type) {
    case TYPE_FLOAT:
        ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_FLOAT3_ARR));
        if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
            const vec3_t* pos = as_vec3_arr(arg);
            const int64_t len = element_count(arg);
            for (int64_t i = 0; i < len; ++i) {
                md_script_vis_vertex_t v = vertex(pos[i], COLOR_WHITE);
                push_point(v, ctx->vis);
            }
        }
        break;
    case TYPE_INT: {
        if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
            const int* indices = as_int_arr(arg);
            const size_t len = element_count(arg);
            for (size_t i = 0; i < len; ++i) {
                const int64_t idx = (int64_t)ctx_range.beg + (int64_t)indices[i] - 1;
                visualize_atom_index(idx, ctx);
            }
        }
        break;
    }
    case TYPE_IRANGE:
        if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
            const irange_t* ranges = as_irange_arr(arg);
            for (size_t i = 0; i < element_count(arg); ++i) {
                irange_t range = remap_range_to_context(ranges[i], ctx_range);
                if (range_in_range(range, ctx_range)) {
                    if (ctx->mol_ctx) {
                        md_bitfield_iter_t it = md_bitfield_iter_range_create(ctx->mol_ctx, range.beg, range.end);
                        while (md_bitfield_iter_next(&it)) {
                            const int64_t idx = md_bitfield_iter_idx(&it);
						    visualize_atom_index(idx, ctx);
                        }
                    }
                    else {
                        visualize_atom_range(range, ctx);
                    }
                }
            }
        }
        break;
    case TYPE_BITFIELD:
        if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
            md_bitfield_t tmp_bf = {0};
            const md_bitfield_t* in_bf = as_bitfield(arg);

            if (ctx->mol_ctx) {
                md_bitfield_init(&tmp_bf, ctx->temp_alloc);
            }

            for (size_t i = 0; i < element_count(arg); ++i) {
                const md_bitfield_t* bf = &in_bf[i];
                if (ctx->mol_ctx) {
                    md_bitfield_and(&tmp_bf, bf, ctx->mol_ctx);
                    bf = &tmp_bf;
                }
                visualize_atom_mask(bf, ctx);
            }

            if (ctx->mol_ctx) {
                md_bitfield_free(&tmp_bf);
            }
        }
        break;
    default:
        ASSERT(false);
    }
}

// Atom indices is optional and if supplied, will be filled with the atom indices corresponding to the positions
static vec3_t* coordinate_extract(data_t arg, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_COORDINATE_ARR));
    vec3_t* positions = 0;
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    switch (arg.type.base_type) {
    case TYPE_FLOAT:
        ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_FLOAT3_ARR));
        md_array_ensure(positions, element_count(arg), ctx->temp_alloc);
        md_array_push_array(positions, as_vec3_arr(arg), element_count(arg), ctx->temp_alloc);
        break;
    case TYPE_INT: {
        md_array_ensure(positions, element_count(arg), ctx->temp_alloc);
        int* indices = as_int_arr(arg);
        size_t num_idx = element_count(arg);
        for (size_t i = 0; i < num_idx; ++i) {
            // Shift here since we use 1 based indices for atoms
            const int idx = ctx_range.beg + indices[i] - 1;
            vec3_t pos = { ctx->mol->atom.x[idx], ctx->mol->atom.y[idx], ctx->mol->atom.z[idx] };
            md_array_push(positions, pos, ctx->temp_alloc);
        }
        break;
    }
    case TYPE_IRANGE: {
        irange_t* ranges = as_irange_arr(arg);
        size_t num_ranges = element_count(arg);
        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            range = clamp_range(range, ctx_range);

            if (ctx->mol_ctx) {
                md_array_ensure(positions, md_bitfield_popcount_range(ctx->mol_ctx, (uint64_t)range.beg, (uint64_t)range.end), ctx->temp_alloc);
                md_bitfield_iter_t it = md_bitfield_iter_range_create(ctx->mol_ctx, range.beg, range.end);
                while (md_bitfield_iter_next(&it)) {
                    const int64_t idx = md_bitfield_iter_idx(&it);
                    vec3_t pos = { ctx->mol->atom.x[idx], ctx->mol->atom.y[idx], ctx->mol->atom.z[idx] };
                    md_array_push(positions, pos, ctx->temp_alloc);
                }
            }
            else {
                md_array_ensure(positions, (size_t)MAX(0, range.end - range.beg), ctx->temp_alloc);
                for (int j = range.beg; j < range.end; ++j) {
                    vec3_t pos = { ctx->mol->atom.x[j], ctx->mol->atom.y[j], ctx->mol->atom.z[j] };
                    md_array_push(positions, pos, ctx->temp_alloc);
                }
            }
        }
        break;
    }
    case TYPE_BITFIELD: {
        md_bitfield_t tmp_bf = { 0 };
        if (ctx->mol_ctx) {
            md_bitfield_init(&tmp_bf, ctx->temp_alloc);
        }

        md_bitfield_t* bf_arr = as_bitfield(arg);
        size_t num_bf = element_count(arg);

        if (num_bf == 1) {
            md_bitfield_t* bf = bf_arr;
            if (ctx->mol_ctx) {
                md_bitfield_init(&tmp_bf, ctx->temp_alloc);
                md_bitfield_and(&tmp_bf, bf_arr, ctx->mol_ctx);
                bf = &tmp_bf;
            }
            positions = extract_vec3(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, bf, ctx->temp_alloc);
        }
        else {
            for (size_t i = 0; i < num_bf; ++i) {
                md_bitfield_t* bf = &bf_arr[i];
                if (ctx->mol_ctx) {
                    md_bitfield_and(&tmp_bf, bf, ctx->mol_ctx);
                    bf = &tmp_bf;
                }
                vec3_t com = extract_com(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol->atom.mass, bf);
                md_array_push(positions, com, ctx->temp_alloc);
            }
            ASSERT(num_bf == md_array_size(positions));
        }

        if (ctx->mol_ctx) {
            md_bitfield_free(&tmp_bf);
        }
        break;
    }
    default:
        break;
    }

    return positions;
}

static md_array(int) coordinate_extract_indices(data_t arg, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_COORDINATE_ARR));
    md_array(int) out_indices = 0;
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    switch (arg.type.base_type) {
    case TYPE_FLOAT:
        break;
    case TYPE_INT: {
        md_array_resize(out_indices, element_count(arg), ctx->temp_alloc);
        const int* indices = as_int_arr(arg);
        const size_t num_idx = element_count(arg);
        for (size_t i = 0; i < num_idx; ++i) {
            // Shift here since we use 1 based indices for atoms
            out_indices[i] = ctx_range.beg + indices[i] - 1;
        }
        break;
    }
    case TYPE_IRANGE: {
        irange_t* ranges = as_irange_arr(arg);
        size_t num_ranges = element_count(arg);
        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            range = clamp_range(range, ctx_range);
            if (ctx->mol_ctx) {
                md_bitfield_iter_t it = md_bitfield_iter_range_create(ctx->mol_ctx, range.beg, range.end);
                while (md_bitfield_iter_next(&it)) {
                    const int idx = (int)md_bitfield_iter_idx(&it);
                    md_array_push(out_indices, (int)idx, ctx->temp_alloc);
                }
            }
            else {
                md_array_ensure(out_indices, md_array_size(out_indices) + (range.end - range.beg), ctx->temp_alloc);
                for (int j = range.beg; j < range.end; ++j) {
                    md_array_push(out_indices, j, ctx->temp_alloc);
                }
            }
        }
        break;
    }
    case TYPE_BITFIELD: {
        const md_bitfield_t* bf = as_bitfield(arg);
        const size_t num_bf = element_count(arg);

        if (num_bf == 1) {
            md_bitfield_t tmp_bf = { 0 };
            if (ctx->mol_ctx) {
                md_bitfield_init(&tmp_bf, ctx->temp_alloc);
                md_bitfield_and(&tmp_bf, bf, ctx->mol_ctx);
                bf = &tmp_bf;
            }

            const size_t count = md_bitfield_popcount(bf);
            md_array_resize(out_indices, count, ctx->temp_alloc);
            md_bitfield_iter_extract_indices(out_indices, count, md_bitfield_iter_create(bf));

            if (ctx->mol_ctx) {
                md_bitfield_free(&tmp_bf);
            }
        }
        else {
            // @NOTE(Robin): We want to be consistent with the coordinate_extract function which in the case of multiple bitfields returns the com of each bitfield
            // Therefore we do not want in this case to return some index as we cannot map the single position to a single atom
        }

        break;
    }
    default:
        break;
    }

    return out_indices;
}

static vec3_t position_extract_com(data_t arg, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_COORDINATE_ARR));
    vec3_t com = {0};

    md_vm_arena_temp_t tmp_pos = md_vm_arena_temp_begin(ctx->temp_arena);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    switch (arg.type.base_type) {
    case TYPE_FLOAT: {
        // This corresponds to user supplied coordinates
        ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_FLOAT3_ARR));
        vec3_t* in_pos = as_vec3_arr(arg);
        size_t len = element_count(arg);
        vec4_t* xyzw = md_vm_arena_push(ctx->temp_arena, len * sizeof(vec4_t));
        for (size_t i = 0; i < len; ++i) {
            xyzw[i] = vec4_from_vec3(in_pos[i], 1.0f);
        }
        com = md_util_com_compute_vec4(xyzw, len, &ctx->mol->unit_cell);
        break;
    }
    case TYPE_INT: {
        const int32_t* in_idx = as_int_arr(arg);
        size_t num_idx = element_count(arg);
        if (ctx->mol_ctx) {
            int32_t* indices = md_vm_arena_push(ctx->temp_arena, num_idx * sizeof(int32_t));
            for (size_t i = 0; i < num_idx; ++i) {
                indices[i] = remap_index_to_context(in_idx[i], ctx_range);
            }
            in_idx = indices;
        }
        
        com = md_util_com_compute(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol->atom.mass, in_idx, num_idx, &ctx->mol->unit_cell);
        break;
    }
    case TYPE_IRANGE: {
        const irange_t* ranges = as_irange_arr(arg);
        size_t num_ranges = element_count(arg);
        md_array(int32_t) indices = 0;
        // We have multiple ranges which may be overlapping
        // We therefore compute the com for each range and then compute the com from the sub-com
        vec4_t* xyzw_arr = md_vm_arena_push(ctx->temp_arena, num_ranges * sizeof(vec4_t));
        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            range = clamp_range(range, ctx_range);
            if (ctx->mol_ctx) {
                size_t len = md_bitfield_popcount_range(ctx->mol_ctx, range.beg, range.end);
                md_bitfield_iter_extract_indices(indices, len, md_bitfield_iter_range_create(ctx->mol_ctx, range.beg, range.end));
                vec4_t xyzw = vec4_from_vec3(md_util_com_compute(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol->atom.mass, indices, len, &ctx->mol->unit_cell), 1.0f);
                md_array_push(xyzw_arr, xyzw, ctx->temp_alloc);
            }
            else {
                size_t len = range.end - range.beg;
                md_array_ensure(indices, len, ctx->temp_alloc);
                md_array_shrink(indices, len);
                for (int32_t j = range.beg; j < range.end; ++j) {
					md_array_push(indices, j, ctx->temp_alloc);
				}
                xyzw_arr[i] = vec4_from_vec3(md_util_com_compute(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol->atom.mass, indices, len, &ctx->mol->unit_cell), 1.0f);
            }
        }

        com = md_util_com_compute_vec4(xyzw_arr, num_ranges, &ctx->mol->unit_cell);
        break;
    }
    case TYPE_BITFIELD: {
        
        md_bitfield_t tmp_bf = { 0 };
        if (ctx->mol_ctx) {
            md_bitfield_init(&tmp_bf, ctx->temp_alloc);
        }

        md_bitfield_t* in_bf = as_bitfield(arg);
        const size_t  num_bf = element_count(arg);

        if (num_bf == 1) {
            md_bitfield_t* bf = in_bf;
            if (ctx->mol_ctx) {
                md_bitfield_init(&tmp_bf, ctx->temp_alloc);
                md_bitfield_and(&tmp_bf, in_bf, ctx->mol_ctx);
                bf = &tmp_bf;
            }
            size_t len = md_bitfield_popcount(bf);
            int32_t* indices = md_vm_arena_push(ctx->temp_arena, md_bitfield_popcount(bf) * sizeof(int32_t));
            md_bitfield_iter_extract_indices(indices, len, md_bitfield_iter_create(bf));
            com = md_util_com_compute(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol->atom.mass, indices, len, &ctx->mol->unit_cell);
        }
        else {
            // If we have multiple bitfields we compute the center of mass for each bitfield before computing a single com from the sub-coms
            md_array(int32_t) indices = 0;
            vec4_t* xyzw_arr = md_vm_arena_push(ctx->temp_arena, num_bf * sizeof(vec4_t));
            for (size_t i = 0; i < num_bf; ++i) {
                md_bitfield_t* bf = &in_bf[i];
                if (ctx->mol_ctx) {
                    md_bitfield_and(&tmp_bf, bf, ctx->mol_ctx);
                    bf = &tmp_bf;
                }
                size_t len = md_bitfield_popcount(bf);
                md_array_ensure(indices, len, ctx->temp_alloc);
                md_array_shrink(indices, len);
                md_bitfield_iter_extract_indices(indices, len, md_bitfield_iter_create(bf));
                xyzw_arr[i] = vec4_from_vec3(md_util_com_compute(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol->atom.mass, indices, len, &ctx->mol->unit_cell), 1.0f);
            }
            com = md_util_com_compute_vec4(xyzw_arr, num_bf, &ctx->mol->unit_cell);
        }
        break;
    }
    default:
        break;
    }

    md_vm_arena_temp_end(tmp_pos);

    return com;
}

static inline float dihedral_angle(vec3_t p0, vec3_t p1, vec3_t p2, vec3_t p3) {
    const vec3_t b1 = vec3_normalize(vec3_sub(p1, p0));
    const vec3_t b2 = vec3_normalize(vec3_sub(p2, p1));
    const vec3_t b3 = vec3_normalize(vec3_sub(p3, p2));
    const vec3_t c1 = vec3_cross(b1, b2);
    const vec3_t c2 = vec3_cross(b2, b3);
    return atan2f(vec3_dot(vec3_cross(c1, c2), b2), vec3_dot(c1, c2));
}

// This should ideally be some light weight regular expression matching.
// For now we just use joker characters '*'
// "*AD" matches anything that ends with AD
// "AD*" matches anything that begins with AD
// @TODO: "*AD*" matches anything that contains AD
// @TODO: "AD*DA" matches anything that begins with AD and ends with DA

static bool match_query(str_t query, str_t str) {
    if (query.len == 0 || str.len == 0) return false;

    int joker_count = 0;
    for (size_t i = 0; i < query.len; ++i) {
        if (query.ptr[i] == '*') joker_count += 1;
    }

    if (joker_count == 0) return str_eq(query, str);

    size_t loc;
    str_t beg_match =  str_find_char(&loc, query, '*') ? str_substr(query, 0, loc) : (str_t){0};
    str_t end_match = str_rfind_char(&loc, query, '*') ? str_substr(query, loc + 1, SIZE_MAX) : (str_t){0};

    if (!str_empty(beg_match)) {
        if (!str_eq_n(beg_match, str, beg_match.len)) return false;
    }

    if (!str_empty(end_match)) {
        str = str_substr(str, str.len - end_match.len, end_match.len);
        if (!str_eq(end_match, str)) return false;
    }

    return true;
}

static bool validate_query(str_t query) {
    if (query.len == 0) return false;
    int joker_count = 0;
    for (size_t i = 0; i < query.len; ++i) {
        if (query.ptr[i] == '*') {
            joker_count += 1;
        }
    }
    if (joker_count > 2) return false;

    if (joker_count == 2) {
        return query.ptr[0] == '*' && query.ptr[query.len-1] == '*';
    }

    return true;
}

// IMPLEMENTATIONS
// @TODO: Add more here

static int _min_farr  (data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT_ARR));
    (void)ctx;

    if (dst) {
        const float* src = as_float_arr(arg[0]);
        float min_val = 0;

        if (element_count(arg[0]) > 0) {
            min_val = FLT_MAX;
            for (size_t i = 0; i < element_count(arg[0]); ++i) {
                min_val = MIN(min_val, src[i]);
            }
        }

        as_float(*dst) = min_val;
    }
    return 0;
}

static int _max_farr  (data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT_ARR));
    (void)ctx;

    if (dst) {
        const float* src = as_float_arr(arg[0]);
        float max_val = 0;

        if (element_count(arg[0]) > 0) {
            max_val = -FLT_MAX;
            for (size_t i = 0; i < element_count(arg[0]); ++i) {
                max_val = MAX(max_val, src[i]);
            }
        }

        as_float(*dst) = max_val;
    }
    return 0;
}

static inline md_bitfield_t _internal_flatten_bf(const md_bitfield_t* bf_arr, size_t count, md_allocator_i* alloc) {
    if (count == 1) {
        return bf_arr[0];
    }
    
    md_bitfield_t bf = md_bitfield_create(alloc);
    for (size_t i = 0; i < count; ++i) {
        md_bitfield_or_inplace(&bf, &bf_arr[i]);
    }
    return bf;
}

static int _not(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD_ARR));

    if (dst) {
        md_bitfield_t* bf_dst = dst->ptr;
        const md_bitfield_t* bf_src = arg[0].ptr;

        md_bitfield_clear_range(bf_dst, 0, ctx->mol->atom.count);
        const size_t count = element_count(arg[0]);
        for (size_t i = 0; i < count; ++i) {
            md_bitfield_or_inplace(bf_dst, &bf_src[i]);
        }
        md_bitfield_not_inplace(bf_dst, 0, ctx->mol->atom.count);

        if (ctx->mol_ctx) {
            // This is a bit conceptually strange,
            // But if you perform a bitwise negation within a context, only the bits within the context are flipped
            md_bitfield_and_inplace(bf_dst, ctx->mol_ctx); // Store intermediate result in dst
        }
    }

    return 0;
}

static int _and(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_BITFIELD_ARR));

    if (dst) {
        md_bitfield_t* bf_dst = dst->ptr;
        const md_bitfield_t* bf_src[2] = {arg[0].ptr, arg[1].ptr};
        const size_t bf_src_len[2] = {element_count(arg[0]), element_count(arg[1])};

        if (bf_src_len[0] == 1 && bf_src_len[1] == 1) {
            md_bitfield_and(bf_dst, bf_src[0], bf_src[1]);
        }
        else {
            for (size_t i = 0; i < bf_src_len[0]; ++i) {
                md_bitfield_or_inplace(bf_dst, bf_src[0] + i);
            }

            md_bitfield_t bf = _internal_flatten_bf(bf_src[1], bf_src_len[1], ctx->temp_alloc);
            md_bitfield_and_inplace(bf_dst, &bf);
        }

        if (ctx->mol_ctx) {
            md_bitfield_and_inplace(bf_dst, ctx->mol_ctx);
        }
    }

    return 0;
}

static int _or(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_BITFIELD_ARR));

    if (dst) {
        md_bitfield_t* bf_dst = dst->ptr;
        const md_bitfield_t* bf_src[2] = {arg[0].ptr, arg[1].ptr};
        const size_t bf_src_len[2] = {element_count(arg[0]), element_count(arg[1])};

        for (size_t i = 0; i < bf_src_len[0]; ++i) {
            md_bitfield_or_inplace(bf_dst, &bf_src[0][i]);
        }

        for (size_t i = 0; i < bf_src_len[1]; ++i) {
            md_bitfield_or_inplace(bf_dst, &bf_src[1][i]);
        }
        
        if (ctx->mol_ctx) {
            md_bitfield_and_inplace(bf_dst, ctx->mol_ctx);
        }
    }

    return 0;
}

static int _xor(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_BITFIELD_ARR));

    if (dst) {
        md_bitfield_t* bf_dst = dst->ptr;
        const md_bitfield_t* bf_src[2] = { arg[0].ptr, arg[1].ptr };
        const size_t bf_src_len[2] = {element_count(arg[0]), element_count(arg[1])};

        if (bf_src_len[0] == 1 && bf_src_len[1] == 1) {
            md_bitfield_xor(bf_dst, bf_src[0], bf_src[1]);
        } else {
            for (size_t i = 0; i < bf_src_len[0]; ++i) {
                md_bitfield_or_inplace(bf_dst, &bf_src[0][i]);
            }

            md_bitfield_t bf = _internal_flatten_bf(bf_src[1], bf_src_len[1], ctx->temp_alloc);
            md_bitfield_xor_inplace(bf_dst, &bf);
        }

        if (ctx->mol_ctx) {
            md_bitfield_and_inplace(bf_dst, ctx->mol_ctx);
        }
    }

    return 0;
}

static int _dot(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT_ARR));
    (void)ctx;

    if (dst) {
        float* a = (float*)arg[0].ptr;
        float* b = (float*)arg[1].ptr;
        double res = 0; // Accumulate in double, then cast to float
        for (size_t i = 0; i < element_count(arg[0]); ++i) {
            res += (double)a[i] * (double)b[i]; 
        }
        as_float(*dst) = (float)res;
    }
    return 0;
}

static int _cross(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT3));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT3));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT3));
    (void)ctx;

    if (dst) {
        float* a = (float*)arg[0].ptr;
        float* b = (float*)arg[1].ptr;
        float* c = (float*)dst->ptr;
        c[0] = a[1]*b[2] - b[1]*a[2];
        c[1] = a[2]*b[0] - b[2]*a[0];
        c[2] = a[0]*b[1] - b[0]*a[1];
    }

    return 0;
}

static int _length(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT_ARR));
    (void)ctx;

    if (dst) {
        float* a = (float*)arg[0].ptr;
        double res = 0; // Accumulate in double, then cast to float
        for (size_t i = 0; i < element_count(arg[0]); ++i) {
            res += (double)a[i] * (double)a[i];
        }
        as_float(*dst) = (float)sqrt(res);
    }

    return 0;
}


static int _mat4_mul_mat4(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT44));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT44));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT44));
    (void)ctx;

    if (dst) {
        mat4_t* A = (mat4_t*) arg[0].ptr;
        mat4_t* B = (mat4_t*) arg[1].ptr;
        mat4_t* C = (mat4_t*) dst->ptr;
        *C = mat4_mul(*A, *B);
    }
    return 0;
}

static int _mat4_mul_vec4(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type,   (type_info_t)TI_FLOAT4));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT44));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT4));
    (void)ctx;

    if (dst) {
        mat4_t* M = (mat4_t*) arg[0].ptr;
        vec4_t* v = (vec4_t*) arg[1].ptr;
        vec4_t* r = (vec4_t*) dst->ptr;
        *r = mat4_mul_vec4(*M, *v);
    }
    return 0;
}

static int _vec2(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT2));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT));

    (void)ctx;
    if (dst) {
        float (*res) = (float*)dst->ptr;
        res[0] = as_float(arg[0]);
        res[1] = as_float(arg[1]);
    }
    return 0;
}

static int _vec3(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT3));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FLOAT));

    (void)ctx;
    if (dst) {
        float (*res) = (float*)dst->ptr;
        res[0] = as_float(arg[0]);
        res[1] = as_float(arg[1]);
        res[2] = as_float(arg[2]);
    }
    return 0;
}

static int _vec4(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT4));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[3].type, (type_info_t)TI_FLOAT));

    (void)ctx;
    if (dst) {
        float (*res) = (float*)dst->ptr;
        res[0] = as_float(arg[0]);
        res[1] = as_float(arg[1]);
        res[2] = as_float(arg[2]);
        res[3] = as_float(arg[3]);
    }
    return 0;
}

static int _all(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (dst) {
        ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
        (void)arg;
        md_bitfield_t* bf = (md_bitfield_t*)dst->ptr;
        md_bitfield_set_range(bf, 0, ctx->mol->atom.count);
        if (ctx->mol_ctx) {
            md_bitfield_and_inplace(bf, ctx->mol_ctx);
        }
    }
    return 0;
}

static int _name(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));
    ASSERT(ctx && ctx->mol && ctx->mol->atom.type);

    const str_t* str = as_string_arr(arg[0]);
    const size_t num_str = element_count(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = (md_bitfield_t*)dst->ptr;

        for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
            if (ctx->mol_ctx && !md_bitfield_test_bit(ctx->mol_ctx, i)) continue;
            
            str_t atom_str = LBL_TO_STR(ctx->mol->atom.type[i]);
            for (size_t j = 0; j < num_str; ++j) {
                if (match_query(str[j], atom_str)) {
                    md_bitfield_set_bit(bf, i);
                }
            }
        }
    } else {
        // We are only validating the arguments here, making sure that they are represented within the potential context
        for (size_t j = 0; j < num_str; ++j) {
            if (!validate_query(str[j])) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' is not a valid query", str[j].len, str[j].ptr);
                return -1;
            }
            bool match = false;
            for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
                if (ctx->mol_ctx && !md_bitfield_test_bit(ctx->mol_ctx, i)) continue;
                str_t atom_str = LBL_TO_STR(ctx->mol->atom.type[i]);
                if (match_query(str[j], atom_str)) {
                    match = true;
                    break;
                }
            }
            if (!match) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' did not match any atom label within the structure", str[j].len, str[j].ptr);
                return -1;
            }
        }
    }

    return 0;
}

static int _element_str(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));
    ASSERT(ctx && ctx->mol);

    if (!ctx->mol->atom.element) {
        if (!dst) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "Dataset does not contain elements, perhaps it is coarse grained?");
            return -1;
        }
        return 0;
    }

    uint8_t* elem_idx = 0;

    const size_t num_str = element_count(arg[0]);
    const str_t* str = as_string_arr(arg[0]);

    for (size_t i = 0; i < num_str; ++i) {
        md_element_t elem = md_util_element_lookup(str[i]);
        if (elem)
            md_array_push(elem_idx, elem, ctx->temp_alloc);
        else if (!dst) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "Failed to map '%.*s' into any Element.", str[i].len, str[i].ptr);
            return -1;
        }
    }
    const size_t num_elem = md_array_size(elem_idx);
    if (!dst && num_elem == 0) {
        LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "No valid arguments in Element");
        return -1;
    }

    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);

        md_bitfield_reserve_range(bf, 0, ctx->mol->atom.count);

        for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
            if (ctx->mol_ctx && !md_bitfield_test_bit(ctx->mol_ctx, i)) continue;
            for (size_t j = 0; j < num_elem; ++j) {
                if (elem_idx[j] == ctx->mol->atom.element[i]) {
                    md_bitfield_set_bit(bf, i);
                    //bit_set_idx(result.bits, i);
                    break;
                }
            }
        }
    } else {
        for (size_t j = 0; j < num_elem; ++j) {
            bool found = false;
            for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
                if (elem_idx[j] == ctx->mol->atom.element[i]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "Element '%.*s' was not found within structure.", str[j].len, str[j].ptr);
                return -1;
            }
        }
    }

    return 0;
}

static int _element_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    ASSERT(ctx && ctx->mol && ctx->mol->atom.element);

    irange_t* ranges = as_irange_arr(arg[0]);
    const size_t num_ranges = element_count(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    if (dst) {
        ASSERT(dst->ptr && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);
        for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
            if (ctx->mol_ctx && !md_bitfield_test_bit(ctx->mol_ctx, i)) continue;
            for (size_t j = 0; j < num_ranges; ++j) {
                if (idx_in_range(ctx->mol->atom.element[i], ranges[j])) {
                    md_bitfield_set_bit(bf, i);
                    break;
                }
            }
        }
    }
    else {
        for (size_t j = 0; j < num_ranges; ++j) {
            bool match = false;
            for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
                if (ctx->mol_ctx && !md_bitfield_test_bit(ctx->mol_ctx, i)) continue;
                if (idx_in_range(ctx->mol->atom.element[i], ranges[j])) {
                    match = true;
                    break;
                }
            }
            if (!match) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "No element within range (%i:%i) was found within structure.", ranges[j].beg, ranges[j].end);
                return -1;
            }
        }
    }

    return 0;
}

static int coordinate_range(data_t* dst, eval_context_t* ctx, frange_t range_x, frange_t range_y, frange_t range_z) {
    ASSERT(ctx && ctx->mol);
    const md_bitfield_t* src_bf = ctx->mol_ctx;

    if (dst) {
        md_bitfield_t* dst_bf = as_bitfield(*dst);
        if (src_bf) {
            int64_t beg_bit = src_bf->beg_bit;
            int64_t end_bit = src_bf->end_bit;
            while ((beg_bit = md_bitfield_scan(src_bf, beg_bit, end_bit)) != 0) {
                const int64_t idx = beg_bit - 1;

                if (range_x.beg <= ctx->mol->atom.x[idx] && ctx->mol->atom.x[idx] <= range_x.end &&
                    range_y.beg <= ctx->mol->atom.y[idx] && ctx->mol->atom.y[idx] <= range_y.end &&
                    range_z.beg <= ctx->mol->atom.z[idx] && ctx->mol->atom.z[idx] <= range_z.end)
                {
                    md_bitfield_set_bit(dst_bf, idx);
                }
            }
        }
        else {
            md_bitfield_reserve_range(dst_bf, 0, ctx->mol->atom.count);
            for (size_t idx = 0; idx < ctx->mol->atom.count; ++idx) {
                if (range_x.beg <= ctx->mol->atom.x[idx] && ctx->mol->atom.x[idx] <= range_x.end &&
                    range_y.beg <= ctx->mol->atom.y[idx] && ctx->mol->atom.y[idx] <= range_y.end &&
                    range_z.beg <= ctx->mol->atom.z[idx] && ctx->mol->atom.z[idx] <= range_z.end)
                {
                    md_bitfield_set_bit(dst_bf, idx);
                }
            }
        }
        return 0;
    }
    else {
        int count = 0;
        if (src_bf) {
            int64_t beg_bit = src_bf->beg_bit;
            int64_t end_bit = src_bf->end_bit;
            while ((beg_bit = md_bitfield_scan(src_bf, beg_bit, end_bit)) != 0) {
                const int64_t idx = beg_bit - 1;
                if (range_x.beg <= ctx->mol->atom.x[idx] && ctx->mol->atom.x[idx] <= range_x.end &&
                    range_y.beg <= ctx->mol->atom.y[idx] && ctx->mol->atom.y[idx] <= range_y.end &&
                    range_z.beg <= ctx->mol->atom.z[idx] && ctx->mol->atom.z[idx] <= range_z.end)
                {
                    count += 1;
                }
            }
        }
        else {
            for (size_t idx = 0; idx < ctx->mol->atom.count; ++idx) {
                if (range_x.beg <= ctx->mol->atom.x[idx] && ctx->mol->atom.x[idx] <= range_x.end &&
                    range_y.beg <= ctx->mol->atom.y[idx] && ctx->mol->atom.y[idx] <= range_y.end &&
                    range_z.beg <= ctx->mol->atom.z[idx] && ctx->mol->atom.z[idx] <= range_z.end)
                {
                    count += 1;
                }
            }
        }
        return count;
    }
}

static int _within_x(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FRANGE));
    return coordinate_range(dst, ctx, as_frange(arg[0]), (frange_t){-FLT_MAX, FLT_MAX}, (frange_t){-FLT_MAX, FLT_MAX});
}

static int _within_y(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FRANGE));
    return coordinate_range(dst, ctx, (frange_t){-FLT_MAX, FLT_MAX}, as_frange(arg[0]), (frange_t){-FLT_MAX, FLT_MAX});
}

static int _within_z(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FRANGE));
    return coordinate_range(dst, ctx, (frange_t){-FLT_MAX, FLT_MAX}, (frange_t){-FLT_MAX, FLT_MAX}, as_frange(arg[0]));
}

static int _within_xyz(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FRANGE));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FRANGE));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FRANGE));
    return coordinate_range(dst, ctx, as_frange(arg[0]), as_frange(arg[1]), as_frange(arg[2]));
}

static bool within_float_iter(const md_spatial_hash_elem_t* elem, int mask, void* user_data) {
    md_bitfield_t* bf = (md_bitfield_t*)user_data;
    while (mask) {
        const int idx = ctz32(mask);
        md_bitfield_set_bit(bf, elem[idx].idx);
        mask &= ~(1 << idx);
    }
    return true;
}

static int _within_expl_flt(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x && ctx->mol->atom.y && ctx->mol->atom.z);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT));

    const float radius = as_float(arg[1]);

    if (dst || ctx->vis) {
        const vec3_t* in_pos = coordinate_extract(arg[0], ctx);
        const size_t num_pos = md_array_size(in_pos);
        const md_spatial_hash_t* sh = get_spatial_hash(ctx);
        const md_bitfield_t* bf_mask = 0;
        md_bitfield_t* bf_dst = 0;
        
        if (is_type_equivalent(arg[0].type, (type_info_t)TI_BITFIELD)) {           
            bf_mask = as_bitfield(arg[0]);
        }

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
            bf_dst = as_bitfield(*dst);
            for (size_t i = 0; i < num_pos; ++i) {
                md_spatial_hash_query_batch(sh, in_pos[i], radius, within_float_iter, bf_dst);
            }
        } else {
            // Visualize
            ASSERT(ctx->vis);
            bf_dst = ctx->vis_structure;
            coordinate_visualize(arg[0], ctx);
            
            for (size_t i = 0; i < num_pos; ++i) {
                push_sphere(in_pos[i], radius, COLOR_WHITE, ctx->vis);
                md_spatial_hash_query_batch(sh, in_pos[i], radius, within_float_iter, bf_dst);
            }
        }
            
        if (bf_mask) {
            ASSERT(bf_dst);
            // Exclude the input atoms from the result
            md_bitfield_andnot_inplace(bf_dst, bf_mask);
        }
    }
    else {
        if (coordinate_validate(arg[0], 1, ctx) < 0) return -1;
        if (radius <= 0) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[1], "The supplied radius is negative or zero, please supply a positive value");
            return -1;
        }
    }

    return 0;
}

static int _within_impl_flt(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x && ctx->mol->atom.y && ctx->mol->atom.z);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT));

    const float radius = as_float(arg[0]);

    if (dst || ctx->vis) {
        ASSERT(ctx->mol_ctx);
        const md_spatial_hash_t* sh = get_spatial_hash(ctx);
        const vec3_t* in_pos = extract_vec3(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol_ctx, ctx->temp_alloc);
        const size_t num_pos = md_array_size(in_pos);
            
        md_bitfield_t* bf_dst = 0;

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
            bf_dst = as_bitfield(*dst);
            for (size_t i = 0; i < num_pos; ++i) {
                md_spatial_hash_query_batch(sh, in_pos[i], radius, within_float_iter, bf_dst);
            }
        }
        else {
            // Visualize
            ASSERT(ctx->vis);
            bf_dst = ctx->vis_structure;
            visualize_atom_mask(ctx->mol_ctx, ctx);
            for (size_t i = 0; i < num_pos; ++i) {
                push_sphere(in_pos[i], radius, COLOR_WHITE, ctx->vis);
                md_spatial_hash_query_batch(sh, in_pos[i], radius, within_float_iter, bf_dst);
            }
        }
        // Exclude the input atoms from the result
        md_bitfield_andnot_inplace(bf_dst, ctx->mol_ctx);
    }
    else {
        if (radius <= 0) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The supplied radius is negative or zero, please supply a positive value");
            return -1;
        }
        if (!ctx->mol_ctx) {
            LOG_ERROR(ctx->ir, ctx->op_token, "The operation is missing its context");
            return -1;
        }
    }

    return 0;
}

typedef struct within_frng_payload_t {
    md_bitfield_t* bf;
    const float min_r2;
    const float max_r2;
    vec4_t ref_pos;
    vec4_t pbc_ext;
} within_frng_payload_t;

static bool within_frng_iter(const md_spatial_hash_elem_t* elem, void* user_data) {
    within_frng_payload_t* data = (within_frng_payload_t*)user_data;
    
    const float d2 = vec4_periodic_distance_squared(data->ref_pos, vec4_from_vec3(elem->xyz, 0), data->pbc_ext);
    if (data->min_r2 < d2 && d2 < data->max_r2) {
        md_bitfield_set_bit(data->bf, elem->idx);
    }
    return true;
}

static int _within_expl_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.z);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FRANGE));

    const frange_t rad_range = as_frange(arg[1]);

    if (dst || ctx->vis) {
        const md_spatial_hash_t* sh = get_spatial_hash(ctx);
        const vec3_t* in_pos = coordinate_extract(arg[0], ctx);
        const size_t num_pos = md_array_size(in_pos);
        const md_bitfield_t* bf_mask = 0;
        md_bitfield_t* bf_dst = 0;
            
        if (is_type_equivalent(arg[0].type, (type_info_t)TI_BITFIELD)) {           
            bf_mask = as_bitfield(arg[0]);
        }
        
        const vec4_t pbc_ext = vec4_from_vec3(mat3_mul_vec3(ctx->mol->unit_cell.basis, vec3_set1(1.0f)), 0);
        within_frng_payload_t payload = {
            .min_r2 = rad_range.beg * rad_range.beg,
            .max_r2 = rad_range.end * rad_range.end,
            .pbc_ext = pbc_ext,
        };
        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
            bf_dst = as_bitfield(*dst);
            payload.bf = bf_dst;
            for (size_t i = 0; i < num_pos; ++i) {
                const float rad = rad_range.end;
                payload.ref_pos = vec4_from_vec3(in_pos[i], 0);
                md_spatial_hash_query(sh, in_pos[i], rad, within_frng_iter, &payload);
            }
        }
        else {
            // Visualize
            bf_dst = ctx->vis_structure;
            payload.bf = bf_dst;
            for (size_t i = 0; i < num_pos; ++i) {
                const float rad = rad_range.end;
                payload.ref_pos = vec4_from_vec3(in_pos[i], 0);
                md_spatial_hash_query(sh, in_pos[i], rad, within_frng_iter, &payload);
                push_sphere(in_pos[i], rad, COLOR_WHITE, ctx->vis);
            }
        }
        if (bf_mask) {
            ASSERT(bf_dst);
            // Exclude the input atoms from the result
            md_bitfield_andnot_inplace(bf_dst, bf_mask);
        }
    }
    else {
        if (coordinate_validate(arg[0], 1, ctx) < 0) return -1;
        if (rad_range.beg < 0 || rad_range.end < rad_range.beg) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[1], "The supplied radius range is invalid, ");
            return -1;
        }
    }

    return 0;
}

static int _within_impl_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.z);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FRANGE));

    const frange_t rad_range = as_frange(arg[0]);

    if (dst || ctx->vis) {
        ASSERT(ctx->mol_ctx);
        const md_spatial_hash_t* sh = get_spatial_hash(ctx);
        const vec3_t* in_pos = extract_vec3(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol_ctx, ctx->temp_alloc);
        const size_t num_pos = md_array_size(in_pos);
        md_bitfield_t* bf_dst = 0;
           
        const vec4_t pbc_ext = vec4_from_vec3(mat3_mul_vec3(ctx->mol->unit_cell.basis, vec3_set1(1.0f)), 0);
        within_frng_payload_t payload = {
            .min_r2 = rad_range.beg * rad_range.beg,
            .max_r2 = rad_range.end * rad_range.end,
            .pbc_ext = pbc_ext,
        };
        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
            bf_dst = as_bitfield(*dst);
            payload.bf = bf_dst;
            for (size_t i = 0; i < num_pos; ++i) {
                const float rad = rad_range.end;
                payload.ref_pos = vec4_from_vec3(in_pos[i], 0);
                md_spatial_hash_query(sh, in_pos[i], rad, within_frng_iter, &payload);
            }
        }
        else {
            // Visualize
            bf_dst = ctx->vis_structure;
            payload.bf = bf_dst;
            for (size_t i = 0; i < num_pos; ++i) {
                const float rad = rad_range.end;
                payload.ref_pos = vec4_from_vec3(in_pos[i], 0);
                md_spatial_hash_query(sh, in_pos[i], rad, within_frng_iter, &payload);
                push_sphere(in_pos[i], rad, COLOR_WHITE, ctx->vis);
            }
        }
        md_bitfield_andnot_inplace(bf_dst, ctx->mol_ctx);
    } else {
        if (rad_range.beg < 0 || rad_range.end < rad_range.beg) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The supplied radius range is invalid, ");
            return -1;
        }
        if (!ctx->mol_ctx) {
            LOG_ERROR(ctx->ir, ctx->op_token, "The operation is missing its context");
            return -1;
        }
    }

    return 0;
}

static int _atom_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    ASSERT(ctx && ctx->mol);

    const size_t num_ranges = element_count(arg[0]);
    const irange_t* ranges = as_irange_arr(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    if (dst) {
        ASSERT(dst->ptr && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);

        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            range = clamp_range(range, ctx_range);
            if (range.beg < range.end) {
                md_bitfield_set_range(bf, range.beg, range.end);
            }
        }
        // Apply context if supplied
        if (ctx->mol_ctx) {
            md_bitfield_and_inplace(bf, ctx->mol_ctx);
        }
    }
    else {
        ASSERT(ctx->arg_tokens);
        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = ranges[i];
            if (!range_in_range(remap_range_to_context(range, ctx_range), ctx_range)) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "supplied range (%i:%i) is not within expected range (%i:%i)",
                    range.beg, range.end, 1, ctx_range.end - ctx_range.beg);
                return -1;
            }
        }
    }

    return 0;
}

static int _atom_int(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_INT_ARR));
    ASSERT(ctx && ctx->mol);

    const size_t num_indices = element_count(arg[0]);
    const int* indices = as_int_arr(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    if (dst) {
        ASSERT(dst->ptr && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);
        for (size_t i = 0; i < num_indices; ++i) {
            int32_t idx = remap_index_to_context(indices[i], ctx_range);
            ASSERT(idx_in_range(idx, ctx_range));
            md_bitfield_set_bit(bf, idx);
        }
    }
    else {
        ASSERT(ctx->arg_tokens);
        for (size_t i = 0; i < num_indices; ++i) {
            if (!idx_in_range(remap_index_to_context(indices[i], ctx_range), ctx_range)) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "supplied index (%i) is not within expected range (%i:%i)",
                    indices[i], 1, ctx_range.end - ctx_range.beg);
                return -1;
            }
        }
    }

    return 0;
}

static int _ring(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    (void)arg;

    int result = 0;

    if (dst || ctx->vis) {
        if (dst) {
            if (!dst->ptr) return 0;
            
            ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD_ARR));
            md_bitfield_t* bf_arr = as_bitfield(*dst);
        
            const size_t num_rings = md_index_data_count(ctx->mol->rings);
            int64_t dst_idx = 0;
            for (size_t i = 0; i < num_rings; ++i) {
                md_bitfield_t* bf = &bf_arr[dst_idx];
            
                // Only accept the ring if it is fully within the given context
                const md_atom_idx_t* ring_beg = md_index_range_beg(ctx->mol->rings, i);
                const md_atom_idx_t* ring_end = md_index_range_end(ctx->mol->rings, i);

                if (ctx->mol_ctx) {
                    bool discard = false;
                    for (const md_atom_idx_t* it = ring_beg; it != ring_end; ++it) {
                        if (!md_bitfield_test_bit(ctx->mol_ctx, *it)) {
                            discard = true;
                            break;
                        }
                    }
                    if (discard) continue;
                }

                for (const md_atom_idx_t* it = ring_beg; it != ring_end; ++it) {
                    md_bitfield_set_bit(bf, *it);
                }
                dst_idx += 1;
            }
        } else {
            ASSERT(ctx->vis);
            const size_t num_rings = md_index_data_count(ctx->mol->rings);
            for (size_t i = 0; i < num_rings; ++i) {
                // Only visualize the ring if it is fully within the given context
                const md_atom_idx_t* ring_beg = md_index_range_beg(ctx->mol->rings, i);
                const md_atom_idx_t* ring_end = md_index_range_end(ctx->mol->rings, i);

                if (ctx->mol_ctx) {
                    bool discard = false;
                    for (const md_atom_idx_t* it = ring_beg; it != ring_end; ++it) {
                        if (!md_bitfield_test_bit(ctx->mol_ctx, *it)) {
                            discard = true;
                            break;
                        }
                    }
                    if (discard) continue;
                }

                md_script_vis_vertex_t vbeg = vertex(vec3_set(ctx->mol->atom.x[*ring_beg], ctx->mol->atom.y[*ring_beg], ctx->mol->atom.z[*ring_beg]), COLOR_WHITE);
                md_script_vis_vertex_t v0 = vbeg;
                for (const md_atom_idx_t* it = ring_beg+1; it != ring_end; ++it) {
                    md_script_vis_vertex_t v1 = vertex(vec3_set(ctx->mol->atom.x[*it], ctx->mol->atom.y[*it], ctx->mol->atom.z[*it]), COLOR_WHITE);
                    push_line(v0, v1, ctx->vis);
                    v0 = v1;
                }
                push_line(v0, vbeg, ctx->vis);
            }
        }
    } else {
        // We need to check if the ring is within the given context
        if (ctx->mol_ctx) {
            int count = 0;
            const size_t num_rings = md_index_data_count(ctx->mol->rings);
            for (size_t i = 0; i < num_rings; ++i) {
                bool discard = false;
                const md_atom_idx_t* ring_beg = md_index_range_beg(ctx->mol->rings, i);
                const md_atom_idx_t* ring_end = md_index_range_end(ctx->mol->rings, i);
                for (const md_atom_idx_t* it = ring_beg; it != ring_end; ++it) {
                    if (!md_bitfield_test_bit(ctx->mol_ctx, *it)) {
                        discard = true;
                        break;
                    }
                }
                if (!discard) {
                    count += 1;
                }
            }
            result = count;
        } else {
            return (int)md_index_data_count(ctx->mol->rings);
        }
    }

    return result;
}

static int _water(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    (void)arg;

    if (!ctx->mol->residue.atom_offset || !ctx->mol->residue.flags) {
        return 0;
    }

    int result = 0;
    int32_t* res_indices = get_residue_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);

    if (dst) {
        ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD_ARR));
        const int64_t cap = type_info_array_len(dst->type);
        if (cap == 0) return 0;

        md_bitfield_t* bf = (md_bitfield_t*)dst->ptr;
        int64_t dst_idx = 0;
        for (size_t i = 0; i < md_array_size(res_indices); ++i) {
            int32_t res_idx = res_indices[i];
            const md_flags_t flags = ctx->mol->residue.flags[res_idx];
            const md_range_t range = md_residue_atom_range(ctx->mol->residue, i);
            if (flags & MD_FLAG_WATER) {
                md_bitfield_set_range(&bf[dst_idx], range.beg, range.end);
                // Do not progress this if we are evaluating in a filter context (we want a single bitfield then)
                dst_idx = (cap == 1) ? dst_idx : dst_idx + 1;
            }
        }
    } else {
        // Length query
        int count = 0;
        for (size_t i = 0; i < md_array_size(res_indices); ++i) {
            int32_t res_idx = res_indices[i];
            const md_flags_t flags = ctx->mol->residue.flags[res_idx];
            count += (flags & MD_FLAG_WATER) ? 1 : 0;
        }
        if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
            count = MIN(1, count);
        }
        result = count;
    }

    md_array_free(res_indices, ctx->temp_alloc);
    
    return result;
}

static int _protein(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    (void)arg;

    int result = 0;

    int32_t* res_indices = get_residue_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);

    if (dst) {
        ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD_ARR));
        md_bitfield_t* bf = (md_bitfield_t*)dst->ptr;
        const int64_t cap = type_info_array_len(dst->type);
        if (cap == 0) return 0;

        // Its a pretty good guess that the bitfield should be in the ballpark of the size of the entire thing.
        md_bitfield_reserve_range(bf, 0, ctx->mol->atom.count);
        
        int64_t dst_idx = 0;
        for (size_t i = 0; i < md_array_size(res_indices); ++i) {
            int32_t res_idx = res_indices[i];
            if (ctx->mol->residue.flags[res_idx] & MD_FLAG_AMINO_ACID) {
                const md_range_t range = md_residue_atom_range(ctx->mol->residue, res_idx);
                ASSERT(dst_idx < cap);
                md_bitfield_set_range(&bf[dst_idx], range.beg, range.end);
                dst_idx = (cap == 1) ? dst_idx : dst_idx + 1;
            }
        }
    }
    else {
        if (!ctx->mol->residue.count) {
            LOG_ERROR(ctx->ir, ctx->op_token, "The system does not contain any residues");
            return -1;
        }

        int count = 0;
        for (size_t i = 0; i < md_array_size(res_indices); ++i) {
            int32_t res_idx = res_indices[i];
            if (ctx->mol->residue.flags[res_idx] & MD_FLAG_AMINO_ACID) {
                count += 1;
            }
        }
        if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
            count = MIN(1, count);
        }
        result = count;
    }

    md_array_free(res_indices, ctx->temp_alloc);

    return result;
}

static int _nucleic(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    (void)arg;

    int result = 0;

    int32_t* res_indices = get_residue_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);

    if (dst) {
        ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD_ARR));
        md_bitfield_t* bf = (md_bitfield_t*)dst->ptr;
        const int64_t cap = type_info_array_len(dst->type);
        if (cap == 0) return 0;

        int64_t dst_idx = 0;
        for (size_t i = 0; i < md_array_size(res_indices); ++i) {
            int32_t res_idx = res_indices[i];
            if (ctx->mol->residue.flags[res_idx] & MD_FLAG_NUCLEOTIDE) {
                const md_range_t range = md_residue_atom_range(ctx->mol->residue, res_idx);
                ASSERT(dst_idx < cap);
                md_bitfield_set_range(&bf[dst_idx], range.beg, range.end);
                dst_idx = (cap == 1) ? dst_idx : dst_idx + 1;
            }
        }
    }
    else {
        if (!ctx->mol->residue.count) {
            LOG_ERROR(ctx->ir, ctx->op_token, "The system does not contain any residues");
            return -1;
        }

        int count = 0;
        for (size_t i = 0; i < md_array_size(res_indices); ++i) {
            int32_t res_idx = res_indices[i];
            if (ctx->mol->residue.flags[res_idx] & MD_FLAG_NUCLEOTIDE) {
                count += 1;
            }
        }
        if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
            count = MIN(1, count);
        }
        result = count;
    }

    md_array_free(res_indices, ctx->temp_alloc);

    return result;
}

static int _ion(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    (void)arg;

    int result = 0;

    if (dst) {
        ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);

        if (ctx->mol_ctx) {
            md_bitfield_iter_t it = md_bitfield_iter_create(ctx->mol_ctx);
            while (md_bitfield_iter_next(&it)) {
				if (ctx->mol->atom.flags[it.idx] & MD_FLAG_ION) {
					md_bitfield_set_bit(bf, it.idx);
				}
			}
        } else {
            for (size_t i = 0; i < ctx->mol->atom.count; ++i) {
                if (ctx->mol->atom.flags[i] & MD_FLAG_ION) {
                    md_bitfield_set_bit(bf, i);
                }
            }
        }
    }

    return result;
}



static int _residue(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));

    int result = 0;

    const size_t num_ranges = element_count(arg[0]);
    const irange_t* ranges = as_irange_arr(arg[0]);
    // Here we use the implicit range given by the context and use that to select substructures within it
    // The supplied iranges will be used as relative indices into the context
    const irange_t ctx_range = get_residue_range_in_context(ctx->mol, ctx->mol_ctx);
    const int32_t ctx_size = ctx_range.end - ctx_range.beg;

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD_ARR));
        if (dst->ptr) {
            md_bitfield_t* bf_arr = (md_bitfield_t*)dst->ptr;
            const int64_t cap = type_info_array_len(dst->type);
            if (cap == 0) return 0;

            int64_t dst_idx = 0;
            for (size_t i = 0; i < num_ranges; ++i) {
                irange_t range = remap_range_to_context(ranges[i], ctx_range);
                range = clamp_range(range, ctx_range);
                for (int64_t j = range.beg; j < range.end; ++j) {
                    //const uint64_t offset = ctx->mol->residue.atom_range[j].beg;
                    //const uint64_t length = ctx->mol->residue.atom_range[j].end - ctx->mol->residue.atom_range[j].beg;
                    //bit_set(result.bits, offset, length);
                    const md_range_t atom_range = md_residue_atom_range(ctx->mol->residue, j);
                    ASSERT(dst_idx < cap);
                    md_bitfield_t* bf = &bf_arr[dst_idx];
                    dst_idx = (cap == 1) ? dst_idx : dst_idx + 1;
                    md_bitfield_set_range(bf, atom_range.beg, atom_range.end);
                    if (ctx->mol_ctx) {
                        md_bitfield_and_inplace(bf, ctx->mol_ctx);
                    }
                }
            }
        }
    }
    else {
        int count = 0;
        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            if (!range_in_range(range, ctx_range)) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "supplied range (%i:%i) is not within expected range (%i:%i)",
                    ranges[i].beg, ranges[i].end, 1, ctx_size);
                return -1;
            }
            count += range.end - range.beg;
        }
        if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
            count = MIN(1, count);
        }
        result = count;
    }

    return result;
}

static int _fill_residue(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD_ARR));
    const size_t num_bf = element_count(arg[0]);
    const md_bitfield_t* src_bf = as_bitfield(arg[0]);

    md_bitfield_t tmp_bf = {0};
    if (num_bf > 1) {
        md_bitfield_init(&tmp_bf, ctx->temp_alloc);
        for (size_t i = 0; i < num_bf; ++i) {
            md_bitfield_or_inplace(&tmp_bf, &src_bf[i]);
        }
        src_bf = &tmp_bf;
    }

    int result = STATIC_VALIDATION_ERROR;
    if (ctx->mol && ctx->mol->residue.atom_offset) {
        if (dst) {
            md_bitfield_t* dst_bf = as_bitfield(*dst);
            const int cap = type_info_array_len(dst->type);
            if (cap == 0) return 0;
            ASSERT(cap > 0);

            const size_t capacity = cap;
            const size_t inc = capacity > 1 ? 1 : 0;

            size_t dst_idx = 0;
            for (size_t i = 0; i < ctx->mol->residue.count; ++i) {
                const md_range_t range = md_residue_atom_range(ctx->mol->residue, i);
                ASSERT(dst_idx <= capacity);
                size_t popcount = md_bitfield_popcount_range(src_bf, range.beg, range.end);
                if (popcount) {
                    md_bitfield_set_range(&dst_bf[dst_idx], range.beg, range.end);
                    dst_idx += inc;
                }
            }

            result = 0;
        } else {
            int count = 0;
            
            if (ctx->backchannel) {
                if (ctx->arg_flags[0] & FLAG_DYNAMIC) {
                    ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
                }
            }

            for (size_t i = 0; i < ctx->mol->residue.count; ++i) {
                const md_range_t range = md_residue_atom_range(ctx->mol->residue, i);
				size_t popcount = md_bitfield_popcount_range(src_bf, range.beg, range.end);
				if (popcount) {
					count += 1;
				}
			}

            if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
                count = MIN(1, count);
            }
            result = count;
        }
    }

    if (num_bf > 1) {
        md_bitfield_free(&tmp_bf);
    }

    return result;
}

static int _fill_chain(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD_ARR));

    const int64_t num_bf = element_count(arg[0]);
    const md_bitfield_t* src_bf = as_bitfield(arg[0]);
    
    md_bitfield_t tmp_bf = {0};
    if (num_bf > 1) {
        md_bitfield_init(&tmp_bf, ctx->temp_alloc);
        for (int64_t i = 0; i < num_bf; ++i) {
            md_bitfield_or_inplace(&tmp_bf, &src_bf[i]);
        }
        src_bf = &tmp_bf;
    }

    int result = 0;
    if (ctx->mol && ctx->mol->chain.atom_offset) {
        if (dst) {
            md_bitfield_t* dst_bf = as_bitfield(*dst);
            const int64_t cap = type_info_array_len(dst->type);
            if (cap == 0) return 0;
            const int64_t inc = cap > 1 ? 1 : 0;

            int64_t dst_idx = 0;
            for (size_t i = 0; i < ctx->mol->chain.count; ++i) {
                ASSERT(dst_idx <= cap);
                md_range_t range = md_chain_atom_range(ctx->mol->chain, i);
                if (md_bitfield_popcount_range(src_bf, range.beg, range.end)) {
                    md_bitfield_set_range(&dst_bf[dst_idx], range.beg, range.end);
                    dst_idx += inc;
                }
            }

            result = 0;
        } else {
            int count = 0;

            if (ctx->backchannel) {
                if (ctx->arg_flags[0] & FLAG_DYNAMIC) {
                    ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
                }
            }

            for (size_t i = 0; i < ctx->mol->chain.count; ++i) {
                md_range_t range = md_chain_atom_range(ctx->mol->chain, i);
                size_t popcount = md_bitfield_popcount_range(src_bf, range.beg, range.end);
                if (popcount) {
                    count += 1;
                }
            }

            if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
                count = MIN(1, count);
            }
            result = count;
        }
    }

    if (num_bf > 1) {
        md_bitfield_free(&tmp_bf);
    }

    return result;
}

static int _resname(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));

    int result = 0;

    if (!ctx->mol->residue.name) {
        LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The molecule does not contain any residue names");
        return -1;
    }

    const size_t num_queries = element_count(arg[0]);
    const str_t* queries = as_string_arr(arg[0]);
    
    // Here we only pick the residues which are represented within the context (by having bits set)
    int32_t* res_indices = get_residue_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);
    const size_t res_count = md_array_size(res_indices);

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD_ARR));
        if (dst->ptr) {
            md_bitfield_t* bf_arr = (md_bitfield_t*)dst->ptr;
            const int64_t cap = type_info_array_len(dst->type);
            if (cap == 0) return 0;
            int64_t dst_idx = 0;

            for (size_t i = 0; i < res_count; ++i) {
                const int64_t res_idx = res_indices[i];
                for (size_t j = 0; j < num_queries; ++j) {
                    if (match_query(queries[j], LBL_TO_STR(ctx->mol->residue.name[res_idx]))) {
                        const md_range_t atom_range = md_residue_atom_range(ctx->mol->residue, res_idx);
                        ASSERT(dst_idx < cap);
                        md_bitfield_t* bf = &bf_arr[dst_idx];
                        dst_idx = (cap == 1) ? dst_idx : dst_idx + 1;
                        // @TODO: These two operations could be replaced with a single copy_range routine if implemented in the future
                        md_bitfield_set_range(bf, atom_range.beg, atom_range.end);
                        if (ctx->mol_ctx) {
                            md_bitfield_and_inplace(bf, ctx->mol_ctx);
                        }

                        break;
                    }
                }
            }
        }
    }
    else {
        int count = 0;
        for (size_t j = 0; j < num_queries; ++j) {
            if (!validate_query(queries[j])) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' is not a valid query", queries[j].len, queries[j].ptr);
                return -1;
            }
            bool match = false;
            for (size_t i = 0; i < res_count; ++i) {
                const int64_t res_idx = res_indices[i];
                if (match_query(queries[j], LBL_TO_STR(ctx->mol->residue.name[res_idx]))) {
                    count += 1;
                    match = true;
                }
            }
            if (!match) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' did not match any residue within the context", queries[j].len, queries[j].ptr);
                return -1;
            }
        }
        if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
            count = MIN(1, count);
        }
        result = count;
    }

    md_array_free(res_indices, ctx->temp_alloc);

    return result;
}

static int _resid(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));

    int result = 0;

    if (!ctx->mol->residue.id) {
        LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The molecule does not contain any residue ids");
        return -1;
    }

    const size_t num_rid = element_count(arg[0]);
    const irange_t*  rid = arg[0].ptr;

    // Here we only pick the residues which are represented within the context (by having bits set)
    int32_t* res_indices = get_residue_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD_ARR));
        if (dst->ptr) {
            md_bitfield_t* bf_arr = (md_bitfield_t*)dst->ptr;
            const int64_t cap = type_info_array_len(dst->type);
            if (cap == 0) return 0;

            int64_t dst_idx = 0;
            for (size_t i = 0; i < md_array_size(res_indices); ++i) {
                const int64_t res_idx = res_indices[i];
                for (size_t j = 0; j < num_rid; ++j) {
                    if (idx_in_range((int)ctx->mol->residue.id[res_idx], rid[j])) {
                        const md_range_t atom_range = md_residue_atom_range(ctx->mol->residue, res_idx);
                        md_bitfield_t* bf = &bf_arr[dst_idx];
                        dst_idx = (cap == 1) ? dst_idx : dst_idx + 1;
                        md_bitfield_set_range(bf, atom_range.beg, atom_range.end);
                        if (ctx->mol_ctx) {
                            md_bitfield_and_inplace(bf, ctx->mol_ctx);
                        }

                        ASSERT(dst_idx <= cap);
                        break;
                    }
                }
            }
        }
    }
    else {
        int count = 0;
        for (size_t j = 0; j < num_rid; ++j) {
            bool match = false;
            for (size_t i = 0; i < md_array_size(res_indices); ++i) {
                const int64_t res_idx = res_indices[i];
                if (idx_in_range((int)ctx->mol->residue.id[res_idx], rid[j])) {
                    count += 1;
                    match = true;
                }
            }
            if (!match) {
                // @TODO: Should this just be a soft warning instead?
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "No matching residue id was found within the range (%i:%i)", rid[j].beg, rid[j].end);
                return -1;
            }
        }
        if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
            count = MIN(1, count);
        }
        result = count;
    }

    return result;
}

static int _chain_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));

    if (ctx->mol->chain.count == 0) {
        LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The molecule does not contain any chains");
        return -1;
    }

    int result = 0;

    const size_t num_ranges = element_count(arg[0]);
    const irange_t* ranges = as_irange_arr(arg[0]);
    const irange_t ctx_range = get_chain_range_in_context(ctx->mol, ctx->mol_ctx);
    const int32_t ctx_size = ctx_range.end - ctx_range.beg;

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD_ARR));
        if (dst->ptr) {
            md_bitfield_t* bf_arr = (md_bitfield_t*)dst->ptr;
            const int64_t cap = type_info_array_len(dst->type);
            if (cap == 0) return 0;

            int64_t dst_idx = 0;
            for (size_t i = 0; i < num_ranges; ++i) {
                irange_t range = remap_range_to_context(ranges[i], ctx_range);
                for (int64_t j = range.beg; j < range.end; ++j) {
                    const md_range_t atom_range = md_chain_atom_range(ctx->mol->chain, j);
                    ASSERT(dst_idx < cap);
                    md_bitfield_t* bf = &bf_arr[dst_idx];
                    dst_idx = (cap == 1) ? dst_idx : dst_idx + 1;
                    md_bitfield_set_range(bf, atom_range.beg, atom_range.end);
                    if (ctx->mol_ctx) {
                        md_bitfield_and_inplace(bf, ctx->mol_ctx);
                    }

                }
            }
        }
    }
    else {
        int count = 0;
        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            if (!range_in_range(range, ctx_range)) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "supplied range (%i:%i) is not within expected range (%i:%i)",
                    ranges[i].beg, ranges[i].end, 1, ctx_size);
                return -1;
            }

            count += range.end - range.beg;
        }
        if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
            count = MIN(1, count);
        }
        result = count;
    }

    return result;
}

static int _chain_str(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));

    if (ctx->mol->chain.count == 0 || !ctx->mol->chain.id) {
        LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The molecule does not contain any chains");
        return -1;
    }

    int result = 0;

    const size_t num_str = element_count(arg[0]);
    const str_t* str     = arg[0].ptr;
    int32_t* chain_indices = get_chain_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD_ARR));
        if (dst->ptr) {
            md_bitfield_t* bf_arr = (md_bitfield_t*)dst->ptr;
            const int64_t cap = type_info_array_len(dst->type);
            if (cap == 0) return 0;

            int64_t dst_idx = 0;
            for (size_t i = 0; i < md_array_size(chain_indices); ++i) {
                for (size_t j = 0; j < num_str; ++j) {
                    if (match_query(str[j], LBL_TO_STR(ctx->mol->chain.id[i]))) {
                        ASSERT(dst_idx < cap);
                        md_bitfield_t* bf = &bf_arr[dst_idx];
                        dst_idx = (cap == 1) ? dst_idx : dst_idx + 1;
                        const md_range_t atom_range = md_chain_atom_range(ctx->mol->chain, i);
                        md_bitfield_set_range(bf, atom_range.beg, atom_range.end);
                        if (ctx->mol_ctx) {
                            md_bitfield_and_inplace(bf, ctx->mol_ctx);
                        }

                        break;
                    }
                }
            }
        }
    }
    else {
        int count = 0;
        for (size_t j = 0; j < num_str; ++j) {
            if (!validate_query(str[j])) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' is not a valid query", str[j].len, str[j].ptr);
                return -1;
            }
            int pre_count = count;
            for (size_t i = 0; i < md_array_size(chain_indices); ++i) {
                if (match_query(str[j], LBL_TO_STR(ctx->mol->chain.id[i]))) {
                    count += 1;
                }
            }
            if (pre_count == count) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' did not match any chain within the structure", str[j].len, str[j].ptr);
                return -1;
            }
        }
        if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
            count = MIN(1, count);
        }
        result = count;
    }

    md_array_free(chain_indices, ctx->temp_alloc);

    return result;
}

// Property Compute

static int _distance(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        vec3_t a = position_extract_com(arg[0], ctx);
        vec3_t b = position_extract_com(arg[1], ctx);

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT));
            float* dst_arr = as_float_arr(*dst);
            md_util_unit_cell_distance_array(dst_arr, &a, 1, &b, 1, &ctx->mol->unit_cell);
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
            coordinate_visualize(arg[1], ctx);

            if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                md_script_vis_vertex_t va = vertex(a, COLOR_WHITE);
                md_script_vis_vertex_t vb = vertex(b, COLOR_WHITE);
                push_point(va, ctx->vis);
                push_point(vb, ctx->vis);
                push_line(va, vb, ctx->vis);
            }
        }
    }
    else {
        int res_a = coordinate_validate(arg[0], 0, ctx);
        int res_b = coordinate_validate(arg[1], 1, ctx);
        if (res_a < 0) return res_a;
        if (res_b < 0) return res_b;
        if (ctx->backchannel) {
            ctx->backchannel->unit = md_unit_angstrom();
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
        return 1;
    }

    return 0;
}

static int _distance_min(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        const vec3_t* a_pos = coordinate_extract(arg[0], ctx);
        const vec3_t* b_pos = coordinate_extract(arg[1], ctx);
        const size_t  a_len = md_array_size(a_pos);
        const size_t  b_len = md_array_size(b_pos);

        int64_t min_i, min_j;
        float min_dist = md_util_unit_cell_min_distance(&min_i, &min_j, a_pos, a_len, b_pos, b_len, &ctx->mol->unit_cell);

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT));
            as_float(*dst) = min_dist;
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
            coordinate_visualize(arg[1], ctx);

            if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                if (a_pos && b_pos) {
                    md_script_vis_vertex_t va = vertex(a_pos[min_i], COLOR_WHITE);
                    md_script_vis_vertex_t vb = vertex(b_pos[min_j], COLOR_WHITE);
                    push_point(va, ctx->vis);
                    push_point(vb, ctx->vis);
                    push_line(va, vb, ctx->vis);
                }
            }
        }
    }
    else {
        int res_a = coordinate_validate(arg[0], 0, ctx);
        int res_b = coordinate_validate(arg[1], 1, ctx);
        if (res_a < 0) return res_a;
        if (res_b < 0) return res_b;
        if (ctx->backchannel) {
            ctx->backchannel->unit = md_unit_angstrom();
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
        return 1;
    }

    return 0;
}

static int _distance_max(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        const vec3_t* a_pos = coordinate_extract(arg[0], ctx);
        const vec3_t* b_pos = coordinate_extract(arg[1], ctx);
        const size_t  a_len = md_array_size(a_pos);
        const size_t  b_len = md_array_size(b_pos);

        int64_t max_i, max_j;
        float max_dist = md_util_unit_cell_min_distance(&max_i, &max_j, a_pos, a_len, b_pos, b_len, &ctx->mol->unit_cell);

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT));
            as_float(*dst) = max_dist;
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
            coordinate_visualize(arg[1], ctx);

            if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                if (a_pos && b_pos) {
                    md_script_vis_vertex_t va = vertex(a_pos[max_i], COLOR_WHITE);
                    md_script_vis_vertex_t vb = vertex(a_pos[max_j], COLOR_WHITE);
                    push_point(va, ctx->vis);
                    push_point(vb, ctx->vis);
                    push_line(va, vb, ctx->vis);
                }
            }
        }
    } else {
        int res_a = coordinate_validate(arg[0], 0, ctx);
        int res_b = coordinate_validate(arg[1], 1, ctx);
        if (res_a < 0) return res_a;
        if (res_b < 0) return res_b;
        if (ctx->backchannel) {
            ctx->backchannel->unit = md_unit_angstrom();
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
        return 1;
    }

    return 0;
}

static int _distance_pair(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));
    
    int result = 0;

    if (dst || ctx->vis) {
        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(ctx->temp_arena);

        vec3_t* a_pos = coordinate_extract(arg[0], ctx);
        vec3_t* b_pos = coordinate_extract(arg[1], ctx);
        const int64_t a_len = md_array_size(a_pos);
        const int64_t b_len = md_array_size(b_pos);

        if (dst) {
            ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
            ASSERT(element_count(*dst) == md_array_size(a_pos) * md_array_size(b_pos));
            float* dst_arr = as_float_arr(*dst);

            md_util_unit_cell_distance_array(dst_arr, a_pos, a_len, b_pos, b_len, &ctx->mol->unit_cell);
        }
        if (ctx->vis) {
            if (ctx->subscript_ranges) {
                // @NOTE(Robin): This is not trivial...
                // We want to support visualization of subranges, which can be imposed by for example array subscripts or hovering a line plot corresponding to one subelement.
                // In such case we have this encoded in ctx->subscript_ranges, which gives us the range of elements to visualize.
                // This means we have to do a bit of extra work to figure out which elements to visualize, since they need to map 1:1 to the elements in the result array (dst)
                int32_t* a_idx = coordinate_extract_indices(arg[0], ctx);
                int32_t* b_idx = coordinate_extract_indices(arg[1], ctx);

                for (size_t ri = 0; ri < md_array_size(ctx->subscript_ranges); ++ri) {
                    irange_t range = ctx->subscript_ranges[ri];
                    range.beg -= 1; // @NOTE: We subtract 1 here because the ranges are 1-indexed

                    int prev_a = -1;
                    md_script_vis_vertex_t va = {0};

                    for (int i = range.beg; i < range.end; ++i) {
                        int a = i % (int)a_len;
                        int b = i / (int)a_len;

                        if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
                            if (a_idx) {
                                visualize_atom_index(a_idx[a], ctx);
                            } else if (arg[0].type.base_type == TYPE_BITFIELD && element_count(arg[0]) > 1) {
                                if (i < (int)element_count(arg[0])) {
                                	const md_bitfield_t* bf_arr = as_bitfield(arg[0]);
                                	visualize_atom_mask(&bf_arr[i], ctx);
                                }
                            }

                            if (b_idx) {
                                visualize_atom_index(b_idx[b], ctx);
                            } else if (arg[1].type.base_type == TYPE_BITFIELD && element_count(arg[1]) > 1) {
                                if (i < (int)element_count(arg[1])) {
                                    const md_bitfield_t* bf_arr = as_bitfield(arg[1]);
                                    visualize_atom_mask(&bf_arr[i], ctx);
                                }
                            }
                        }

                        if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                            if (a != prev_a) {
                                prev_a = a;
                                va = vertex(a_pos[a], COLOR_WHITE);
                                push_point(va, ctx->vis);
                            }
                            md_script_vis_vertex_t vb = vertex(b_pos[b], COLOR_WHITE);
                            push_point(vb, ctx->vis);
                            push_line(va, vb, ctx->vis);
                        }
                    }
                }
            } else {
                coordinate_visualize(arg[0], ctx);
                coordinate_visualize(arg[1], ctx);

                if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                    for (int64_t i = 0; i < a_len; ++i) {
                        md_script_vis_vertex_t va = vertex(a_pos[i], COLOR_WHITE);
                        push_point(va, ctx->vis);
                        for (int64_t j = 0; j < b_len; ++j) {
                            md_script_vis_vertex_t vb = vertex(b_pos[j], COLOR_WHITE);
                            push_point(vb, ctx->vis);
                            push_line(va, vb, ctx->vis);
                        }
                    }
                }
            }
        }
        md_vm_arena_temp_end(temp);
    } else {
        int res0 = coordinate_validate(arg[0], 0, ctx);
        int res1 = coordinate_validate(arg[1], 1, ctx);
        if (res0 < 0) return res0;
        if (res1 < 0) return res1;
        if (ctx->backchannel) {
            ctx->backchannel->unit = md_unit_angstrom();
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
        int64_t count = (int64_t)res0 * (int64_t)res1;
        if (count > 1000000) {
            LOG_ERROR(ctx->ir, ctx->op_token, "The size produced by the operation is %"PRId64", which exceeds the upper limit of 1'000'000", count);
            return STATIC_VALIDATION_ERROR;
        }
        result = (int)count;
    }

    return result;
}

static int _angle(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_COORDINATE_ARR));
    
    if (dst || ctx->vis) {
        const vec3_t a = position_extract_com(arg[0], ctx);
        const vec3_t b = position_extract_com(arg[1], ctx);
        const vec3_t c = position_extract_com(arg[2], ctx);
        const vec3_t v0 = vec3_normalize(vec3_sub(a, b));
        const vec3_t v1 = vec3_normalize(vec3_sub(c, b));

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT));
            as_float(*dst) = acosf(vec3_dot(v0, v1));
        }

        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
            coordinate_visualize(arg[1], ctx);
            coordinate_visualize(arg[2], ctx);

            if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                md_script_vis_vertex_t va = vertex(a, COLOR_WHITE);
                md_script_vis_vertex_t vb = vertex(b, COLOR_WHITE);
                md_script_vis_vertex_t vc = vertex(c, COLOR_WHITE);

                push_point(va, ctx->vis);
                push_point(vb, ctx->vis);
                push_point(vc, ctx->vis);
                push_line(va, vb, ctx->vis);
                push_line(vb, vc, ctx->vis);

                // This is the angle arc
                // @TODO: Insert more vertices here and make it a bit smoother
                md_script_vis_vertex_t vd = vertex(vec3_add(b, vec3_mul_f(v0, 0.5f)), COLOR_WHITE);
                md_script_vis_vertex_t ve = vertex(vec3_add(b, vec3_mul_f(v1, 0.5f)), COLOR_WHITE);

                push_triangle(vb, vd, ve, ctx->vis);
            }
        }
    } else {
        // Validate input
        if (coordinate_validate(arg[0], 0, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (coordinate_validate(arg[1], 1, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (coordinate_validate(arg[2], 2, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (ctx->backchannel) {
            ctx->backchannel->unit = md_unit_radian();
            //ctx->backchannel->value_range = (frange_t){0, PI};
        }
    }

    return 0;
}

static int _dihedral(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[3].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        const vec3_t a = position_extract_com(arg[0], ctx);
        const vec3_t b = position_extract_com(arg[1], ctx);
        const vec3_t c = position_extract_com(arg[2], ctx);
        const vec3_t d = position_extract_com(arg[3], ctx);

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT));
            as_float(*dst) = (float)dihedral_angle(a,b,c,d);
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
            coordinate_visualize(arg[1], ctx);
            coordinate_visualize(arg[2], ctx);
            coordinate_visualize(arg[3], ctx);

            if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {

                md_script_vis_vertex_t va = vertex(a, COLOR_WHITE);
                md_script_vis_vertex_t vb = vertex(b, COLOR_WHITE);
                md_script_vis_vertex_t vc = vertex(c, COLOR_WHITE);
                md_script_vis_vertex_t vd = vertex(d, COLOR_WHITE);

                push_point(va, ctx->vis);
                push_point(vb, ctx->vis);
                push_point(vc, ctx->vis);
                push_point(vd, ctx->vis);

                push_line(va, vb, ctx->vis);
                push_line(vb, vc, ctx->vis);
                push_line(vc, vd, ctx->vis);

                // @TODO: Draw planes and the angle between them
            }
        }
    } else {
        // Validate input
        if (coordinate_validate(arg[0], 0, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (coordinate_validate(arg[1], 1, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (coordinate_validate(arg[2], 2, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (coordinate_validate(arg[3], 3, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (ctx->backchannel) {
            ctx->backchannel->unit = md_unit_radian();
            //ctx->backchannel->value_range = (frange_t){-PI, PI};
        }
    }

    return 0;
}

static int _rmsd(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD_ARR));
    ASSERT(ctx && ctx->mol && ctx->mol->atom.mass);

    bool result = 0;

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT));
        ASSERT(ctx->initial_configuration.x && ctx->initial_configuration.y && ctx->initial_configuration.z);

        if (dst->ptr) {
            const md_bitfield_t* src_bf = as_bitfield(arg[0]);
            md_bitfield_t bf = _internal_flatten_bf(src_bf, element_count(arg[0]), ctx->temp_alloc);
            md_bitfield_t tmp_bf = {0};

            if (ctx->mol_ctx) {
                md_bitfield_init(&tmp_bf, ctx->temp_alloc);
                md_bitfield_and (&tmp_bf, &bf, ctx->mol_ctx);
                bf = tmp_bf;
            }
            const size_t count = md_bitfield_popcount(&bf);
            if (count > 0) {
                int32_t* indices = md_alloc(ctx->temp_alloc, sizeof(int32_t) * count);
                md_bitfield_iter_extract_indices(indices, count, md_bitfield_iter_create(&bf));

                const float* const in_x[2] = {ctx->initial_configuration.x, (const float*)ctx->mol->atom.x};
                const float* const in_y[2] = {ctx->initial_configuration.y, (const float*)ctx->mol->atom.y};
                const float* const in_z[2] = {ctx->initial_configuration.z, (const float*)ctx->mol->atom.z};
                const float* in_w = ctx->mol->atom.mass;

                const vec3_t com[2] = {
                    md_util_com_compute(in_x[0], in_y[0], in_z[0], in_w, NULL, count, NULL),
                    md_util_com_compute(in_x[1], in_y[1], in_z[1], in_w, NULL, count, NULL),
                };

                as_float(*dst) = (float)md_util_rmsd_compute(in_x, in_y, in_z, in_w, indices, count, com);
                md_free(ctx->temp_alloc, indices, sizeof(int32_t) * count);
            }
        }
    }
    else {
        if (ctx->vis) {
            // Visualize
            // I don't think we need another visualization than the bitfield highlighting the atoms involved...
            coordinate_visualize(arg[0], ctx);
        } else {
            // Validate args
            // Nothing really to validate, arguments are of type bitfields and if the bitfield is empty, that would be ok, since that would yield a valid rmsd -> 0.
			// And the empty bitfield must be valid in the case of dynamic selection.
        }
    }

    return result;
}


// #################
// ###   CASTS   ###
// #################

static int _cast_int_to_flt(data_t* dst, data_t arg[], eval_context_t* ctx){
    if (!dst) return 0;
    
    ASSERT(type_info_equal(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(type_info_equal(arg[0].type, (type_info_t)TI_INT));
    (void)ctx;
    as_float(*dst) = (float)as_int(arg[0]);
    return 0;
}

static int _cast_int_to_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    
    ASSERT(type_info_equal(dst->type, (type_info_t)TI_IRANGE));
    ASSERT(type_info_equal(arg[0].type, (type_info_t)TI_INT));
    (void)ctx;
    as_irange(*dst) = (irange_t){as_int(arg[0]), as_int(arg[0])};
    return 0;
}

static int _cast_irng_to_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    
    ASSERT(type_info_equal(dst->type, (type_info_t)TI_FRANGE));
    ASSERT(type_info_equal(arg[0].type, (type_info_t)TI_IRANGE));
    (void)ctx;
    as_frange(*dst) = (frange_t){(float)as_irange(arg[0]).beg, (float)as_irange(arg[0]).end};
    return 0;
}

static int _cast_int_arr_to_irng_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;

    ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_IRANGE_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_INT_ARR));
    (void)ctx;
    const size_t dst_len = element_count(*dst);
    ASSERT(dst_len == element_count(arg[0]));

    irange_t* d = dst->ptr;
    int*      s = arg[0].ptr;

    for (size_t i = 0; i < dst_len; ++i) {
        d[i].beg = s[i];
        d[i].end = s[i];
    }
    return 0;
}

static int _cast_int_arr_to_flt_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    
    ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_INT_ARR));
    (void)ctx;
    const size_t dst_len = element_count(*dst);
    ASSERT(dst_len == element_count(arg[0]));

    float* d = dst->ptr;
    int*   s = arg[0].ptr;

    for (size_t i = 0; i < dst_len; ++i) {
        d[i] = (float)s[i]; 
    }
    return 0;
}

static int _cast_irng_arr_to_frng_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    
    ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FRANGE_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    (void)ctx;
    const size_t dst_len = element_count(*dst);
    ASSERT(dst_len == element_count(arg[0]));

    frange_t* d = dst->ptr;
    irange_t* s = arg[0].ptr;

    for (size_t i = 0; i < dst_len; ++i) {
        d[i].beg = (float)s[i].beg;
        d[i].end = (float)s[i].end; 
    }
    return 0;
}

static int _cast_int_arr_to_bf(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_INT_ARR));
    ASSERT(arg[0].ptr);

    int result = 0;

    const int32_t* indices = as_int_arr(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);
    const int64_t ctx_size = ctx_range.end - ctx_range.beg;
    const size_t num_idx = element_count(arg[0]);

    // The idea here is that if we have a context, we only make sure to add the indices which are represented within the context.

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);
       
        for (size_t i = 0; i < num_idx; ++i) {
            // Shift here since we use 1 based indices for atoms
            const int64_t idx = (int64_t)ctx_range.beg + (int64_t)indices[i] - 1;
            if (ctx->mol_ctx) {
                ASSERT(md_bitfield_test_bit(ctx->mol_ctx, idx));
                // This is pre-checked in the static check bellow in the else statement
            }
            md_bitfield_set_bit(bf, idx);
        }
    } else {
        for (size_t i = 0; i < num_idx; ++i) {
            const int64_t idx = (int64_t)ctx_range.beg + (int64_t)indices[i] - 1;

            if (!(ctx_range.beg <= idx && idx < ctx_range.end)) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "supplied index (%i) is not within expected range (%i:%i)", (int)idx, 1, ctx_size);
                return -1;
            }
            if (ctx->mol_ctx) {
                if (!md_bitfield_test_bit(ctx->mol_ctx, idx)) {
                    LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "supplied index (%i) is not represented within the supplied context", (int)idx);
                    return -1;
                }
            }
        }
    }

    return result;
}

static int _cast_irng_arr_to_bf(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    ASSERT(arg[0].ptr);

    int result = 0;

    const size_t num_ranges = element_count(arg[0]);
    const irange_t* ranges = as_irange_arr(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);
    const int32_t ctx_size = ctx_range.end - ctx_range.beg;

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);
        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = clamp_range(remap_range_to_context(ranges[i], ctx_range), ctx_range);

            if (ctx->mol_ctx) {
                int64_t beg_bit = range.beg;
                int64_t end_bit = range.end;
                while ((beg_bit = md_bitfield_scan(ctx->mol_ctx, beg_bit, end_bit)) != 0) {
                    const int64_t idx = beg_bit - 1;
                    md_bitfield_set_bit(bf, idx);
                }
            }
            else {
                if (range.beg >= 0 && range.end > range.beg) {
                    md_bitfield_set_range(bf, range.beg, range.end);
                }
            }
        }
    } else {
        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            if (!range_in_range(range, ctx_range)) {
                LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "supplied range (%i:%i) is not within expected range (%i:%i)",
                    ranges[i].beg, ranges[i].end, 1, ctx_size);
                return -1;
            }
        }
    }

    return result;
}

static int _join_bf_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)ctx;
    if (!dst) return 0;

    ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD_ARR));

    md_bitfield_t* dst_bf = as_bitfield(*dst);
    const md_bitfield_t* src_bf_arr = (md_bitfield_t*)arg[0].ptr;
    const int64_t count = type_info_array_len(arg[0].type);

    for (int64_t i = 0; i < count; ++i) {
        md_bitfield_or_inplace(dst_bf, &src_bf_arr[i]);
    }

    return 0;
}

static int _com(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));

    vec3_t com;
    if (dst || ctx->vis) {
        com = position_extract_com(arg[0], ctx);
        if (dst) {
            ASSERT(type_info_equal(dst->type, (type_info_t)TI_FLOAT3));
            as_vec3(*dst) = com;
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
            if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                vec3_t* pos = coordinate_extract(arg[0], ctx);
                md_script_vis_vertex_t c = vertex(com, COLOR_WHITE);
                push_point(c, ctx->vis);
                for (size_t i = 0; i < md_array_size(pos); ++i) {
                    md_script_vis_vertex_t v = vertex(pos[i], COLOR_WHITE);
                    push_line(c, v, ctx->vis);
                }
            }
        }
    }
    else {
        return coordinate_validate(arg[0], 0, ctx);
    }
    return 0;
}

static int _plane(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    (void)ctx;

    if (dst || ctx->vis) {
        // Compute eigen vectors from covariance matrix
        // normal should be the third axis (i.e. the smallest eigen value)
        md_array(int) indices = coordinate_extract_indices(arg[0], ctx);
        size_t count = md_array_size(indices);

        vec4_t* xyzw = md_alloc(ctx->temp_alloc, sizeof(vec4_t) * count);
        for (size_t i = 0; i < count; ++i) {
            int32_t idx = indices[i];
            float w = ctx->mol->atom.mass ? ctx->mol->atom.mass[idx] : 1.0f;
			xyzw[i] = (vec4_t){ctx->mol->atom.x[idx], ctx->mol->atom.y[idx], ctx->mol->atom.z[idx], w};
		}

        md_util_unwrap_vec4(xyzw, count, &ctx->mol->unit_cell);
        vec3_t com = md_util_com_compute_vec4(xyzw, count, 0); // @NOTE: No need to supply the unit cell here since we already unwrapped the structure
        mat3_t M = mat3_covariance_matrix_vec4(xyzw, 0, count, com);
        mat3_eigen_t eigen = mat3_eigen(M);

        vec3_t normal = vec3_normalize(eigen.vectors.col[2]);
        float d = vec3_dot(normal, com);
        vec4_t plane = {normal.x, normal.y, normal.z, d};

        if (dst) {
            ASSERT(type_info_equal(dst->type, (type_info_t)TI_FLOAT4));
            as_vec4(*dst) = plane;
        }

        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
            if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                vec3_t v = vec3_mul_f(eigen.vectors.col[0], eigen.values.elem[0]);
                vec3_t u = vec3_mul_f(eigen.vectors.col[1], eigen.values.elem[1]);

                // Draw something plane-ish
                md_script_vis_vertex_t ci = vertex(com, COLOR_WHITE);
                md_script_vis_vertex_t ni = vertex(vec3_add(com, normal), COLOR_WHITE);

                md_script_vis_vertex_t corner[4] = {
                    vertex(vec3_sub(vec3_sub(com, v), u), COLOR_WHITE),
                    vertex(vec3_sub(vec3_add(com, v), u), COLOR_WHITE),
                    vertex(vec3_add(vec3_add(com, v), u), COLOR_WHITE),
                    vertex(vec3_add(vec3_sub(com, v), u), COLOR_WHITE),
                };

                push_line(corner[0], corner[1], ctx->vis);
                push_line(corner[1], corner[2], ctx->vis);
                push_line(corner[2], corner[3], ctx->vis);
                push_line(corner[3], corner[0], ctx->vis);

                push_line(ci, ni, ctx->vis);
            }
        }

        md_free(ctx->temp_alloc, xyzw, sizeof(vec4_t) * count);
    } else {
        int count = coordinate_validate(arg[0], 0, ctx);
        if (count < 0) return -1;
        if (count < 3) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "Invalid number of positions, need at least 3 to compute a plane");
            return -1;
        }
    }

    return 0;
}

static int _coordinate(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        if (dst) {
            ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT3_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            ASSERT(dst->size == md_array_bytes(pos));
            MEMCPY(dst->ptr, pos, md_array_size(pos) * sizeof(vec3_t));
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
        }
    } else {
        return coordinate_validate(arg[0], 0, ctx);
    }
    return 0;
}

static int _coordinate_x(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        if (dst) {
            ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            ASSERT(element_count(*dst) == md_array_size(pos));
            float* x = as_float_arr(*dst);
            for (size_t i = 0; i < element_count(*dst); ++i) {
                x[i] = pos[i].x;
            }
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
        }
    } else {
        return coordinate_validate(arg[0], 0, ctx);
    }
    return 0;
}

static int _coordinate_y(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        if (dst) {
            ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            ASSERT(element_count(*dst) == md_array_size(pos));
            float* y = as_float_arr(*dst);
            for (size_t i = 0; i < element_count(*dst); ++i) {
                y[i] = pos[i].y;
            }
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
        }
    } else {
        return coordinate_validate(arg[0], 0, ctx);
    }
    return 0;
}

static int _coordinate_z(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        if (dst) {
            ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            ASSERT(element_count(*dst) == md_array_size(pos));
            float* z = as_float_arr(*dst);
            for (size_t i = 0; i < element_count(*dst); ++i) {
                z[i] = pos[i].z;
            }
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
        }
    } else {
        return coordinate_validate(arg[0], 0, ctx);
    }
    return 0;
}

static int _coordinate_xy(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        if (dst) {
            ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT2_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            ASSERT(element_count(*dst) == md_array_size(pos));
        
            vec2_t* xy = as_vec2_arr(*dst);
            for (size_t i = 0; i < element_count(*dst); ++i) {
                xy[i] = (vec2_t){pos[i].x, pos[i].y};
            }
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
        }
    } else {
        return coordinate_validate(arg[0], 0, ctx);
    }
    return 0;
}

static int _coordinate_xz(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        if (dst) {
            ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT2_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            ASSERT(element_count(*dst) == md_array_size(pos));

            vec2_t* xz = as_vec2_arr(*dst);
            for (size_t i = 0; i < element_count(*dst); ++i) {
                xz[i] = (vec2_t){pos[i].x, pos[i].z};
            }
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
        }
    } else {
        return coordinate_validate(arg[0], 0, ctx);
    }
    return 0;
}

static int _coordinate_yz(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        if (dst) {
            ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT2_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            ASSERT(element_count(*dst) == md_array_size(pos));

            vec2_t* yz = as_vec2_arr(*dst);
            for (size_t i = 0; i < element_count(*dst); ++i) {
                yz[i] = (vec2_t){pos[i].y, pos[i].z};
            }
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
        }
    } else {
        return coordinate_validate(arg[0], 0, ctx);
    }
    return 0;
}

typedef struct {
    vec4_t pos;
    vec4_t pbc_ext;
    const md_bitfield_t* exclusion_mask;
    float min_cutoff;
    float min_cutoff2;
    float max_cutoff;
    float max_cutoff2;
    float inv_cutoff_range;
    float* bins;
    int32_t num_bins;
    uint32_t idx;
    uint64_t total_count;
} rdf_payload_t;

bool rdf_iter(const md_spatial_hash_elem_t* elem_arr, int mask, void* user_param) {
    rdf_payload_t* data = user_param;
    
    while (mask) {
        const int idx = ctz32(mask);
        const float d2 = vec4_periodic_distance_squared(vec4_from_vec3(elem_arr[idx].xyz, 0), data->pos, data->pbc_ext);
        if (data->min_cutoff2 < d2 && d2 < data->max_cutoff2) {
            const float d = sqrtf(d2);
            int32_t bin_idx = (int32_t)(((d - data->min_cutoff) * data->inv_cutoff_range) * data->num_bins);
            bin_idx = CLAMP(bin_idx, 0, data->num_bins - 1);
            data->bins[bin_idx] += 1.0f;
            data->total_count += 1;
        }
        mask &= ~(1 << idx);
    }
    return true;
}

// iterator version which excludes the supplied index
bool rdf_iter_excl_idx(const md_spatial_hash_elem_t* elem_arr, int mask, void* user_param) {
    rdf_payload_t* data = user_param;

    while (mask) {
        const int idx = ctz32(mask);
        if (elem_arr[idx].idx != data->idx) {
            const float d2 = vec4_periodic_distance_squared(vec4_from_vec3(elem_arr[idx].xyz, 0), data->pos, data->pbc_ext);
            if (data->min_cutoff2 < d2 && d2 < data->max_cutoff2) {
                const float d = sqrtf(d2);
                int32_t bin_idx = (int32_t)(((d - data->min_cutoff) * data->inv_cutoff_range) * data->num_bins);
                bin_idx = CLAMP(bin_idx, 0, data->num_bins - 1);
                data->bins[bin_idx] += 1.0f;
                data->total_count += 1;
            }
        }
        mask &= ~(1 << idx);
    }
    return true;
}

// iterator version which excludes indices that are present in the supplied exclusion_mask
bool rdf_iter_excl_mask(const md_spatial_hash_elem_t* elem_arr, int mask, void* user_param) {
    rdf_payload_t* data = user_param;

    while (mask) {
        const int idx = ctz32(mask);
        if (!md_bitfield_test_bit(data->exclusion_mask, elem_arr[idx].idx)) {
            const float d2 = vec4_periodic_distance_squared(vec4_from_vec3(elem_arr[idx].xyz, 0), data->pos, data->pbc_ext);
            if (data->min_cutoff2 < d2 && d2 < data->max_cutoff2) {
                const float d = sqrtf(d2);
                int32_t bin_idx = (int32_t)(((d - data->min_cutoff) * data->inv_cutoff_range) * data->num_bins);
                bin_idx = CLAMP(bin_idx, 0, data->num_bins - 1);
                data->bins[bin_idx] += 1.0f;
                data->total_count += 1;
            }
        }
        mask &= ~(1 << idx);
    }
    return true;
}

#define RDF_BRUTE_FORCE_LIMIT (3000)

static inline double sphere_volume(double r) {
    return (4.0 / 3.0) * PI * (r * r * r);
}

static void compute_rdf(float* bins, float* weights, int num_bins, const data_t arg[2], float min_cutoff, float max_cutoff, const md_unit_cell_t* unit_cell, md_allocator_i* alloc, eval_context_t* ctx) {
    const float inv_cutoff_range = 1.0f / (max_cutoff - min_cutoff);

    uint64_t total_count = 0;

    // Add a small bias to avoid adding 'zero' distances
    if (min_cutoff == 0) {
        min_cutoff = 1.0e-3f;
    }

    const md_array(int32_t) ref_idx = coordinate_extract_indices(arg[0], ctx);
    md_array(vec3_t) ref_pos = 0;

    const bool use_exclusion_masks = (arg[0].type.base_type == TYPE_BITFIELD) && (element_count(arg[0]) > 1);
    const bool use_exclusion_indices = (ref_idx != 0);

    if (ref_idx) {
        md_array_resize(ref_pos, md_array_size(ref_idx), ctx->temp_alloc);
        for (size_t i = 0; i < md_array_size(ref_idx); ++i) {
            const int64_t idx = ref_idx[i];
			ref_pos[i] = (vec3_t){ctx->mol->atom.x[idx], ctx->mol->atom.y[idx], ctx->mol->atom.z[idx]};
		}
    } else {
		ref_pos = coordinate_extract(arg[0], ctx);
	}
    const size_t ref_len = md_array_size(ref_pos);

    md_array(vec3_t) trg_pos = coordinate_extract(arg[1], ctx);
    const size_t trg_len = md_array_size(trg_pos);

    const vec4_t ext = unit_cell ? vec4_from_vec3(mat3_mul_vec3(unit_cell->basis, vec3_set1(1.0f)), 0) : vec4_zero();
    if ((ref_len * trg_len) < RDF_BRUTE_FORCE_LIMIT) {
        for (size_t i = 0; i < ref_len; ++i) {
            for (size_t j = 0; j < trg_len; ++j) {
                const float d2 = vec4_periodic_distance_squared(vec4_from_vec3(ref_pos[i], 0), vec4_from_vec3(trg_pos[j], 0), ext);
                if (min_cutoff * min_cutoff < d2 && d2 < max_cutoff * max_cutoff) {
                    const float d = sqrtf(d2);
                    int32_t bin_idx = (int32_t)(((d - min_cutoff) * inv_cutoff_range) * num_bins);
                    bin_idx = CLAMP(bin_idx, 0, num_bins - 1);
                    bins[bin_idx] += 1.0f;
                    total_count += 1;
                }
            }
        }
    } else {
        md_spatial_hash_t* hash = md_spatial_hash_create_vec3(trg_pos, NULL, trg_len, unit_cell, alloc);

        rdf_payload_t payload = {
            .pbc_ext = ext,
            .min_cutoff = min_cutoff,
            .min_cutoff2 = min_cutoff * min_cutoff,
            .max_cutoff = max_cutoff,
            .max_cutoff2 = max_cutoff * max_cutoff,
            .inv_cutoff_range = inv_cutoff_range,
            .bins = bins,
            .num_bins = num_bins,
            .total_count = 0,
        };

        if (use_exclusion_indices) {
            for (size_t i = 0; i < ref_len; ++i) {
                payload.pos = vec4_from_vec3(ref_pos[i], 0);
                payload.idx = (uint32_t)i;
                md_spatial_hash_query_batch(hash, ref_pos[i], max_cutoff, rdf_iter_excl_idx, &payload);
            }
        } else if (use_exclusion_masks) {
            for (size_t i = 0; i < ref_len; ++i) {
                payload.pos = vec4_from_vec3(ref_pos[i], 0);
                payload.exclusion_mask = as_bitfield(arg[0]) + i;
                md_spatial_hash_query_batch(hash, ref_pos[i], max_cutoff, rdf_iter_excl_mask, &payload);
            }
        } else {
            for (size_t i = 0; i < ref_len; ++i) {
                payload.pos = vec4_from_vec3(ref_pos[i], 0);
                md_spatial_hash_query_batch(hash, ref_pos[i], max_cutoff, rdf_iter, &payload);
            }
        }

        total_count = payload.total_count;
    }

    // Reference rho
    const double total_vol = sphere_volume(max_cutoff) - sphere_volume(min_cutoff);
    const double ref_rho   = total_count / total_vol;

    // Compute normalization factor of bins
    // With respect to the nominal distribution rho
    // Each bin is normalized with respect to the nominal distribution
    const double dr = (max_cutoff - min_cutoff) / (float)num_bins;
    double prev_sphere_vol = 0;
    for (int64_t i = 0; i < num_bins; ++i) {
        const double sphere_vol = sphere_volume(min_cutoff + (i + 0.5) * dr);
        const double bin_vol = sphere_vol - prev_sphere_vol;
        prev_sphere_vol = sphere_vol;
        weights[i] = (float)(ref_rho * bin_vol);
    }
}

#if 0
static int visualize_rdf(const data_t arg[2], float min_cutoff, float max_cutoff, eval_context_t* ctx) {
    ASSERT(ctx->vis);

    //coordinate_visualize(arg[0], ctx);
    //coordinate_visualize(arg[1], ctx);
    // Visualize
    if (ctx->vis && ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
        const float min_cutoff2 = min_cutoff * min_cutoff;
        const float max_cutoff2 = max_cutoff * max_cutoff;
        // This is stupid: Will do an n^2 visualization of all references and targets of each reference.
        // Clutters the visual space and tells you nothing in the end.
        const int64_t max_vertices = 10000;
        int64_t vertex_count = 0;
        for (int64_t i = 0; i < ref_len; ++i) {
            md_script_vis_vertex_t r_v = vertex(ref_pos[i], COLOR_WHITE);
            vertex_count += 1;
            push_point(r_v, ctx->vis);
            for (int64_t j = 0; j < trg_len; ++j) {
                const vec3_t d = vec3_sub(ref_pos[i], trg_pos[j]);
                const float d2 = vec3_dot(d, d);
                if (min_cutoff2 < d2 && d2 < max_cutoff2) {
                    md_script_vis_vertex_t t_v = vertex(trg_pos[j], COLOR_WHITE);
                    vertex_count += 1;
                    push_line(r_v, t_v, ctx->vis);
                    if (vertex_count > max_vertices) return 0;
                }
            }
        }
    }

    return 0;
}
#endif

static int internal_rdf(data_t* dst, data_t arg[], float min_cutoff, float max_cutoff, eval_context_t* ctx) {
    if (dst || ctx->vis) {
        const int num_bins = MD_DIST_BINS;
       
        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_DISTRIBUTION));
            ASSERT(dst->ptr);
            float* bins    = as_float_arr(*dst);
			float* weights = as_float_arr(*dst) + num_bins;
            compute_rdf(bins, weights, num_bins, arg, min_cutoff, max_cutoff, &ctx->mol->unit_cell, ctx->temp_alloc, ctx);
        }
        if (ctx->vis) {
            //return visualize_rdf(arg, min_cutoff, max_cutoff, ctx);
        }
        return 0;
    }
    else {
        int ref_len = coordinate_validate(arg[0], 0, ctx);
        int trg_len = coordinate_validate(arg[1], 1, ctx);

        // Validate input
        if (ref_len <= 0) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "empty reference positions");
            return STATIC_VALIDATION_ERROR;
        }
        if (trg_len <= 0) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[1], "empty target positions");
            return STATIC_VALIDATION_ERROR;
        }
        if (min_cutoff < 0.0f || max_cutoff <= min_cutoff) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[2], "Invalid cutoff");
            return STATIC_VALIDATION_ERROR;
        }
        if (ctx->arg_flags[2] & FLAG_DYNAMIC) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[2], "Cutoff needs to have a static value in order to determine the upper bound if the distribution");
            return STATIC_VALIDATION_ERROR;
        }
        if (ctx->backchannel) {
            ctx->backchannel->unit = (md_unit_t){ 0 };
            ctx->backchannel->value_range = (frange_t){min_cutoff, max_cutoff};
        }
    }
    return 0;
}

static int _rdf_flt(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FLOAT));

    const float cutoff = as_float(arg[2]);
    return internal_rdf(dst, arg, 0, cutoff, ctx);
}

static int _rdf_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_BITFIELD_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FRANGE));

    const frange_t cutoff = as_frange(arg[2]);
    return internal_rdf(dst, arg, cutoff.beg, cutoff.end, ctx);
}

typedef enum {
    COUNT_TYPE_UNKNOWN = 0,
    COUNT_TYPE_ATOM = 1,
    COUNT_TYPE_RESIDUE = 2,
    COUNT_TYPE_CHAIN = 3,
    COUNT_TYPE_STRUCTURE = 4,
    COUNT_TYPE_COUNT
} count_type_t;

static const str_t count_type_str[] = {
	BAKE_STR(""),
	BAKE_STR("atom"),
	BAKE_STR("residue"),
	BAKE_STR("chain"),
	BAKE_STR("structure"),
};

count_type_t count_type_from_str(str_t str) {
	for (size_t i = 0; i < COUNT_TYPE_COUNT; ++i) {
		if (str_eq(str, count_type_str[i])) {
			return (count_type_t)i;
		}
	}
	return COUNT_TYPE_UNKNOWN;
}

size_t internal_count(const md_bitfield_t* bf, count_type_t count_type, eval_context_t* ctx) {
    ASSERT(bf);
    ASSERT(ctx);

    md_bitfield_t tmp = {0};
    if (ctx->mol_ctx) {
        md_bitfield_and(&tmp, bf, ctx->mol_ctx);
		bf = &tmp;
    }
    
    size_t count = 0;

    switch (count_type) {
    	case COUNT_TYPE_ATOM: {
			return md_bitfield_popcount(bf);
		}
		case COUNT_TYPE_RESIDUE:
			for (size_t i = 0; i < ctx->mol->residue.count; ++i) {
                md_range_t range = md_residue_atom_range(ctx->mol->residue, i);
                if (md_bitfield_popcount_range(bf, range.beg, range.end) > 0) {
					count += 1;
				}
            }
            break;
		case COUNT_TYPE_CHAIN: {
            for (size_t i = 0; i < ctx->mol->chain.count; ++i) {
                md_range_t range = md_chain_atom_range(ctx->mol->chain, i);
                if (md_bitfield_popcount_range(bf, range.beg, range.end) > 0) {
                    count += 1;
                }
            }
            break;
		}
		case COUNT_TYPE_STRUCTURE: {
            size_t num_structures = md_index_data_count(ctx->mol->structures);
            md_bitfield_t tmp2 = md_bitfield_create(ctx->temp_alloc);
			for (size_t i = 0; i < num_structures; ++i) {
                const int32_t* idx = md_index_range_beg(ctx->mol->structures, i);
                const size_t num_idx = md_index_range_size(ctx->mol->structures, i);
            
                md_bitfield_clear(&tmp2);
                md_bitfield_set_indices_u32(&tmp2, (const uint32_t*)idx, num_idx);
                md_bitfield_and_inplace(&tmp2, bf);
                if (md_bitfield_popcount(&tmp2) > 0) {
                	count += 1;
                }
            }
            break;
		}
		default: {
			ASSERT(false);
			break;
		}
    }
    return count;
}

static int _count(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD));

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
        const md_bitfield_t* bf = as_bitfield(arg[0]);

        as_float(*dst) = (float)internal_count(bf, COUNT_TYPE_ATOM, ctx);
    }

    return 0;
}

static int _count_with_arg(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_STRING));

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));

        const md_bitfield_t* bf = as_bitfield(arg[0]);
        count_type_t type = count_type_from_str(as_string(arg[1]));

        as_float(*dst) = (float)internal_count(bf, type, ctx);
    } else {
        str_t str = as_string(arg[1]);
        count_type_t type = count_type_from_str(str);
        if (type == COUNT_TYPE_UNKNOWN) {
            md_strb_t sb = md_strb_create(ctx->temp_alloc);
            for (int i = 1; i < COUNT_TYPE_COUNT; ++i) {
                md_strb_fmt(&sb, "\""STR_FMT"\"\n", STR_ARG(count_type_str[i]));
            }
            LOG_ERROR(ctx->ir, ctx->arg_tokens[1], "Unknown argument: '"STR_FMT"', valid arguments are:\n"STR_FMT, STR_ARG(str), STR_ARG(md_strb_to_str(sb)));
			return -1;
		}
    }

    return 0;
}

// This is quirky and wierd.
// It should essentially check if all the supplied bitfields are 'equivalent'
// Equivalent in this case means that all of them represent the same type of structures within the system
// We do this now by checking against the atom elements if available (which is not the case in CoarseGrained simulations)
// In such case we revert to matching Atom labels
// To resolve this PROPERLY, one needs to determine the topology of each supplied substructure and match it as a graph.

static inline bool are_bitfields_equivalent(const md_bitfield_t bitfields[], int64_t num_bitfields, const md_molecule_t* mol) {
    // Number of bits should match.
    // The atomic element of each set bit should match.

    const md_bitfield_t* ref_bf = &bitfields[0];
    const int64_t ref_count = md_bitfield_popcount(ref_bf);

    // We compare against the first one, which is considered to be the reference
    for (int64_t i = 1; i < num_bitfields; ++i) {
        const md_bitfield_t* bf = &bitfields[i];
        const int64_t count = md_bitfield_popcount(bf);
        if (count != ref_count) {
            return false;
        }

        int64_t beg_bit = bf->beg_bit;
        int64_t end_bit = bf->end_bit;
        while ((beg_bit = md_bitfield_scan(bf, beg_bit, end_bit)) != 0) {
            const int64_t idx = beg_bit - 1;
            const int64_t ref_idx = ref_bf->beg_bit + (idx - bf->beg_bit);
            if (!md_bitfield_test_bit(ref_bf, ref_idx)) {
                return false;
            }
            if (mol->atom.element) {
                if (mol->atom.element[idx] != mol->atom.element[ref_idx]) {
                    return false;
                }
            } else {
                str_t a = LBL_TO_STR(mol->atom.type[idx]);
                str_t b = LBL_TO_STR(mol->atom.type[ref_idx]);
                if (!str_eq(a, b)) {
                    return false;
                }

            }
        }
    }

    return true;
}

static inline void populate_volume(float* vol, const vec3_t* xyz, size_t count, mat4_t M) {
    const vec4_t dim = vec4_set(MD_VOL_DIM, MD_VOL_DIM, MD_VOL_DIM, 0);
    for (size_t i = 0; i < count; ++i) {
        const vec4_t coord = mat4_mul_vec4(M, vec4_from_vec3(xyz[i], 1.0f));
        
        const md_128 a = md_mm_cmplt_ps(md_mm_setzero_ps(), coord.m128);
        const md_128 b = md_mm_cmplt_ps(coord.m128, dim.m128);
        const md_128 c = md_mm_and_ps(a, b);

        // Count the number of lanes which is not zero
        // 0x7 = 1 + 2 + 4 means that first three lanes (x,y,z) are within min and max
        if (_mm_movemask_ps(c) == 0x7) {
            vol[(uint32_t)coord.z * (MD_VOL_DIM * MD_VOL_DIM) + (uint32_t)coord.y * MD_VOL_DIM + (uint32_t)coord.x] += 1.0f;
        }
    }
}

static inline mat4_t compute_volume_matrix(float radius) {
    // We have the cutoff as a radius, meaning our volume radius has the length 'r' for each axis.
    // Thus this means the diameter is 2*r.
    // We have the resolution MD_VOL_DIM for each axis, meaning each voxel has the extent of 2*r / VOL_RES units
    const float voxel_ext = (2*radius) / MD_VOL_DIM;
    const float s = 1.0f / voxel_ext;
    mat4_t S = mat4_scale(s, s, s);

    const float t = (MD_VOL_DIM / 2);
    mat4_t T = mat4_translate(t, t, t);

    return mat4_mul(T, S);
}

typedef struct sdf_payload_t {
    mat4_t M;
    float* vol;
} sdf_payload_t;

bool sdf_iter(const md_spatial_hash_elem_t* elem_arr, int mask, void* user_param) {
    sdf_payload_t* data = user_param;
    
    while (mask) {
        const int idx = ctz32(mask);
        const vec4_t coord = mat4_mul_vec4(data->M, vec4_from_vec3(elem_arr[idx].xyz, 1.0f));

        // Dwelling into intrinsic land here, is a bit unsure if we should expose these as vec4 operations
        const md_128 a = md_mm_cmplt_ps(md_mm_setzero_ps(), coord.m128);
        const md_128 b = md_mm_cmplt_ps(coord.m128, md_mm_set_ps(0, MD_VOL_DIM, MD_VOL_DIM, MD_VOL_DIM));
        const md_128 c = md_mm_and_ps(a, b);

        // Count the number of lanes which are not zero
        // 0x7 = 1 + 2 + 4 means that first three lanes (x,y,z) are within min and max
        if (_mm_movemask_ps(c) == 0x7) {
            data->vol[(uint32_t)coord.z * (MD_VOL_DIM * MD_VOL_DIM) + (uint32_t)coord.y * MD_VOL_DIM + (uint32_t)coord.x] += 1.0f;
        }
        mask &= ~(1 << idx);
    }

    return true;
}

static int _sdf(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FLOAT));
    ASSERT(arg[0].ptr);
    ASSERT(arg[1].ptr);
    ASSERT(arg[2].ptr);

    const int64_t num_ref_bitfields = element_count(arg[0]);
    const md_bitfield_t* ref_bf_arr = as_bitfield(arg[0]);
    const md_bitfield_t* trg_bf     = as_bitfield(arg[1]);
    float cutoff = as_float(arg[2]);

    const md_bitfield_t* ref_bf = &ref_bf_arr[0];

    // We currently only support bitfields which represent equivalent structures. 
    // his each bitfield should contain the same number of atoms, and the n'th atom within the bitfield should match the n'th atom in all other bitfields.
    // This may be loosened in the future to allow aligning structures which are not strictly equivalent.
    int result = 0;

    if (dst || ctx->vis) {
        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(ctx->temp_arena);

        ASSERT(ctx->initial_configuration.x);
        ASSERT(ctx->initial_configuration.y);
        ASSERT(ctx->initial_configuration.z);

        // This could happen if we have dynamic length as input
        if (num_ref_bitfields == 0) return 0;

        const size_t ref_size = md_bitfield_popcount(ref_bf);
        const size_t trg_size = md_bitfield_popcount(trg_bf);

        if (ref_size == 0) return 0;
        if (trg_size == 0) return 0; 

        float* vol = 0;
        if (dst) {
            ASSERT(dst->ptr);
            ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_VOLUME));
            vol = as_float_arr(*dst);
        }

        // Allocate temporary memory
        int32_t* ref_idx[2] = {
            md_alloc(ctx->temp_alloc, sizeof(int32_t) * ref_size),
            md_alloc(ctx->temp_alloc, sizeof(int32_t) * ref_size)
        };

        // Extract indices
        md_bitfield_iter_extract_indices(ref_idx[0], ref_size, md_bitfield_iter_create(ref_bf));

        /*
        const int64_t mem_size = sizeof(float) * (ref_size * 7);
        float* mem = md_alloc(ctx->temp_alloc, mem_size);

        float* ref_x[2] = { mem + ref_size * 0, mem + ref_size * 1 };
        float* ref_y[2] = { mem + ref_size * 2, mem + ref_size * 3 };
        float* ref_z[2] = { mem + ref_size * 4, mem + ref_size * 5 };
        float* ref_w = mem + ref_size * 6;
        */
        const float* const ref_x[2] = {ctx->initial_configuration.x, (const float*)ctx->mol->atom.x};
        const float* const ref_y[2] = {ctx->initial_configuration.y, (const float*)ctx->mol->atom.y};
        const float* const ref_z[2] = {ctx->initial_configuration.z, (const float*)ctx->mol->atom.z};
        const float* const ref_w[2] = {ctx->mol->atom.mass, ctx->mol->atom.mass};
        vec3_t ref_com[2] = {0};

        // Fetch initial reference positions
        //extract_xyzw(ref_x[0], ref_y[0], ref_z[0], ref_w, ctx->initial_configuration.x, ctx->initial_configuration.y, ctx->initial_configuration.z, ctx->mol->atom.mass, ref_bf);

        // Fetch target positions
        vec3_t* trg_xyz = extract_vec3(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, trg_bf, ctx->temp_alloc);
        ref_com[0] = md_util_com_compute(ref_x[0], ref_y[0], ref_z[0], ref_w[0], ref_idx[0], ref_size, &ctx->mol->unit_cell);

        // @TODO(Robin): This should be measured
        const bool brute_force = trg_size < 5000;

        md_spatial_hash_t* spatial_hash = NULL;
        const float spatial_hash_radius = sqrtf(cutoff * cutoff * 3); // distance from center of volume to corners of the volume.

        if (!brute_force) {
            spatial_hash = md_spatial_hash_create_vec3(trg_xyz, 0, trg_size, &ctx->mol->unit_cell, ctx->temp_alloc);
        }

        // A for alignment matrix, Align eigen vectors with axis x,y,z etc.
        mat3_eigen_t eigen = mat3_eigen(mat3_covariance_matrix(ref_x[0], ref_y[0], ref_z[0], ref_w[0], ref_idx[0], ref_size, ref_com[0]));
        mat4_t A = mat4_from_mat3(mat3_transpose(eigen.vectors));

        // V for volume matrix scale and align with the volume which we aim to populate with density
        mat4_t V = compute_volume_matrix(cutoff);
        mat4_t VA = mat4_mul(V, A);

        if (ctx->vis && ctx->vis_flags & MD_SCRIPT_VISUALIZE_SDF) {
            ctx->vis->sdf.extent = cutoff;
        }

        for (int64_t i = 0; i < num_ref_bitfields; ++i) {
            const md_bitfield_t* bf = &ref_bf_arr[i];
            md_bitfield_iter_extract_indices(ref_idx[1], ref_size, md_bitfield_iter_create(bf));

            ref_com[1] = md_util_com_compute(ref_x[1], ref_y[1], ref_z[1], ref_w[1], ref_idx[1], ref_size, &ctx->mol->unit_cell);
            mat3_t R  = mat3_optimal_rotation(ref_x, ref_y, ref_z, ref_w, (const int32_t* const*)ref_idx, ref_size, ref_com);
            mat4_t RT = mat4_mul(mat4_from_mat3(R), mat4_translate(-ref_com[1].x, -ref_com[1].y, -ref_com[1].z));

            if (vol) {
                if (brute_force) {
                    populate_volume(vol, trg_xyz, trg_size, mat4_mul(VA, RT));
                } else {
                    sdf_payload_t payload = {
                        .M = mat4_mul(VA, RT),
                        .vol = vol,
                    };
                    md_spatial_hash_query_batch(spatial_hash, ref_com[1], spatial_hash_radius, sdf_iter, &payload);
                }
            }
            if (ctx->vis && ctx->vis_flags & MD_SCRIPT_VISUALIZE_SDF) {
                md_allocator_i* alloc = ctx->vis->alloc;
                md_bitfield_t bf_cpy = {0};
                md_bitfield_init(&bf_cpy, alloc);
                md_bitfield_copy(&bf_cpy, bf);

                mat4_t ART = mat4_mul(A, RT);
                md_array_push(ctx->vis->sdf.matrices, ART, alloc);
                md_array_push(ctx->vis->sdf.structures, bf_cpy, alloc);
            }
        }

        md_vm_arena_temp_end(temp);
    } else {
        // Validation
        if (num_ref_bitfields < 1) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "Number of bitfields which serve as reference frame must be 1 or more");
            return -1;
        }

        const int64_t ref_bit_count = md_bitfield_popcount(ref_bf);
        if (ref_bit_count <= 0) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The supplied reference bitfield(s) are empty");
            return -1;
        }

        // Test for equivalence
        if (!are_bitfields_equivalent(ref_bf_arr, num_ref_bitfields, ctx->mol)) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[0], "The supplied reference bitfields are not identical: the number of atoms and their corresponding elements do not match between all supplied bitfields");
            return -1;
        }

        const int64_t target_bit_count = md_bitfield_popcount(trg_bf);
        if (target_bit_count <= 0) {
            LOG_ERROR(ctx->ir, ctx->arg_tokens[1], "The supplied target bitfield is empty");
            return -1;
        }

        if (ctx->backchannel) {
            ctx->backchannel->unit = (md_unit_t){0};
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
    }

    return result;
}

static vec3_t compute_shape_weights(data_t arg[], eval_context_t* ctx) {
    vec3_t weights = {0};
    md_vm_arena_temp_t pos = md_vm_arena_temp_begin(ctx->temp_arena);

    md_array(int) indices = 0;

    if (arg[0].type.base_type == TYPE_BITFIELD && element_count(arg[0]) > 0) {
        md_bitfield_t tmp_bf = {0};
        md_bitfield_init(&tmp_bf, ctx->temp_alloc);

        const size_t count = element_count(arg[0]);
        const md_bitfield_t* bf_arr = as_bitfield(arg[0]);
        for (size_t i = 0; i < count; ++i) {
            const md_bitfield_t* bf = &bf_arr[i];
            md_bitfield_or_inplace(&tmp_bf, bf);
        }

        if (ctx->mol_ctx) {
            md_bitfield_and_inplace(&tmp_bf, ctx->mol_ctx);
        }

        const size_t size = md_bitfield_popcount(&tmp_bf);
        if (size > 0) {
            md_array_resize(indices, size, ctx->temp_alloc);
            md_bitfield_iter_extract_indices(indices, size, md_bitfield_iter_create(&tmp_bf));
        }
    } else {
        indices = coordinate_extract_indices(arg[0], ctx);
    }

    if (indices) {
        // @TODO: This needs to be deperiodized on the com
        const size_t count = md_array_size(indices);
        const vec3_t com = md_util_com_compute(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol->atom.mass, indices, count, &ctx->mol->unit_cell);
        const mat3_t M = mat3_covariance_matrix(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol->atom.mass, indices, count, com);
        weights = md_util_shape_weights(&M);
    }

    md_vm_arena_temp_end(pos);
    return weights;
}

static int _shape_weights(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));    
    ASSERT(arg[0].ptr);

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_FLOAT3));
        vec3_t* weights = (vec3_t*)dst->ptr;
        *weights = compute_shape_weights(arg, ctx);
    } else {
        if (ctx->backchannel) {
			ctx->backchannel->unit = (md_unit_t){0};
			ctx->backchannel->value_range = (frange_t){0, 1.0};
		}
    }

    return 0;
}
/*
// This is some experimental future work, for matching structures using maximum common subgraph

static int32_t greedy_align(const md_bitfield_t* a, const md_bitfield_t* b, const md_molecule_t* mol) {
    const int64_t a_size = md_bitfield_popcount(a);
    const int64_t b_size = md_bitfield_popcount(b);
    const md_bitfield_t* small;
    const md_bitfield_t* big;

    if (a_size < b_size) {
        small = a;
        big = b;
    } else {
        small = b;
        big = a;
    }
}

static int _align(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_BITFIELD));

    if (dst) {
        ASSERT(dst->ptr);
        ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_INT_ARR));
    }
}
*/
