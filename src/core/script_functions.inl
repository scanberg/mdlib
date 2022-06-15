// This is to mark that the procedure supports a varying length
#define ANY_LENGTH -1
//#define ANY_LEVEL -1

#define DIST_BINS 1024
#define VOL_DIM 128

#define STATIC_VALIDATION_ERROR -2

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

#define TI_DISTRIBUTION {TYPE_FLOAT, {DIST_BINS,1}, 1}
#define TI_VOLUME       {TYPE_FLOAT, {VOL_DIM,VOL_DIM,VOL_DIM,1}, 3}

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

#define TI_POSITION {TYPE_POSITION, {1}, 0}
#define TI_POSITION_ARR {TYPE_POSITION, {ANY_LENGTH}, 0}

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
            for (int64_t i = 0; i < element_count(*dst); ++i) { \
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
        for (int64_t i = 0; i < element_count(*dst); ++i) { \
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
        for (int64_t i = 0; i < element_count(*dst); ++i) { \
            ((base_type*)dst->ptr)[i] = ((base_type*)arg[0].ptr)[i] op *((base_type*)arg[1].ptr); \
        } \
        return 0; \
    }

#define BAKE_OP_M_M(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        if (!dst) return 0; \
        for (int64_t i = 0; i < element_count(*dst); ++i) { \
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
        int64_t count_a = type_info_element_stride_count(arg[0].type); \
        int64_t count_b = type_info_element_stride_count(arg[1].type); \
        ASSERT(count_a == count_b); \
        ASSERT(count_a % md_simd_widthf == 0); \
        const float* src_a = as_float_arr(arg[0]); \
        const float* src_b = as_float_arr(arg[1]); \
        float* dst_arr = as_float_arr(*dst); \
        for (int64_t i = 0; i < count_a; i += md_simd_widthf) { \
            md_simd_typef a = md_simd_loadf(src_a + i); \
            md_simd_typef b = md_simd_loadf(src_b + i); \
            md_simd_storef(dst_arr + i, op(a, b)); \
        } \
        return 0; \
    }

#define BAKE_SIMDF_OP_M_S(name, op) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        ASSERT(dst); \
        ASSERT(dst->type.base_type == TYPE_FLOAT); \
        ASSERT(arg[0].type.base_type == TYPE_FLOAT); \
        ASSERT(is_type_equivalent(arg[1].type, (md_type_info_t)TI_FLOAT)); \
        int64_t count = type_info_element_stride_count(arg[0].type); \
        ASSERT(count % md_simd_widthf == 0); \
        const float* src_arr = as_float_arr(arg[0]); \
        float* dst_arr = as_float_arr(*dst); \
        md_simd_typef s = md_simd_set1f(as_float(arg[1])); \
        for (int64_t i = 0; i < count; i += md_simd_widthf) { \
            md_simd_storef(dst_arr + i, op(md_simd_loadf(src_arr + i), s)); \
        } \
        return 0; \
    }

#define BAKE_SIMDF_OP_M(name, op) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        ASSERT(dst); \
        ASSERT(dst->type.base_type == TYPE_FLOAT); \
        ASSERT(arg[0].type.base_type == TYPE_FLOAT); \
        int64_t count = type_info_element_stride_count(arg[0].type); \
        ASSERT(count % md_simd_widthf == 0); \
        const float* src_arr = as_float_arr(arg[0]); \
        float* dst_arr = as_float_arr(*dst); \
        for (int64_t i = 0; i < count; i += md_simd_widthf) { \
            md_simd_storef(dst_arr + i, op(md_simd_loadf(src_arr + i))); \
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

BAKE_SIMDF_OP_M_M(_op_simd_add_farr_farr, md_simd_addf)
BAKE_SIMDF_OP_M_S(_op_simd_add_farr_f,    md_simd_addf)

BAKE_SIMDF_OP_M_M(_op_simd_sub_farr_farr, md_simd_subf)
BAKE_SIMDF_OP_M_S(_op_simd_sub_farr_f,    md_simd_subf)

BAKE_SIMDF_OP_M_M(_op_simd_mul_farr_farr, md_simd_mulf)
BAKE_SIMDF_OP_M_S(_op_simd_mul_farr_f,    md_simd_mulf)

BAKE_SIMDF_OP_M_M(_op_simd_div_farr_farr, md_simd_divf)
BAKE_SIMDF_OP_M_S(_op_simd_div_farr_f,    md_simd_divf)

BAKE_SIMDF_OP_M(_op_simd_abs_farr, md_simd_absf)

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
static int _flatten_bf_arr              (data_t*, data_t[], eval_context_t*);

// Basic operations
static int _min_farr  (data_t*, data_t[], eval_context_t*); // (float[]) -> float
static int _max_farr  (data_t*, data_t[], eval_context_t*); // (float[]) -> float

// Logical operators for custom types
static int _not  (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _and  (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _or   (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _xor  (data_t*, data_t[], eval_context_t*); // -> bitfield

// Selectors
// Atomic level selectors
static int _all     (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _x       (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _y       (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _z       (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_expl_flt  (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_expl_frng (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_impl_flt  (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_impl_frng (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _name    (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _element_str (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _element_irng(data_t*, data_t[], eval_context_t*); // -> bitfield
static int _atom_irng    (data_t*, data_t[], eval_context_t*);   // -> bitfield
static int _atom_int     (data_t*, data_t[], eval_context_t*);   // -> bitfield

// Residue level selectors
static int _water   (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _protein   (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _ion     (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _resname (data_t*, data_t[], eval_context_t*); // (str[]) -> bitfield
static int _resid   (data_t*, data_t[], eval_context_t*); // (irange[]) -> bitfield   // irange also covers integers, since integers are implicitly convertible to irange
static int _residue (data_t*, data_t[], eval_context_t*); // (irange[]) -> bitfield

static int _fill_residue(data_t*, data_t[], eval_context_t*);   // (bitfield[]) -> bitfield
static int _fill_chain  (data_t*, data_t[], eval_context_t*);    // (bitfield[]) -> bitfield

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
static int _position(data_t*, data_t[], eval_context_t*);  // (position[]) -> float[3][]
static int _position_x(data_t*, data_t[], eval_context_t*);  // (position[]) -> float[]
static int _position_y(data_t*, data_t[], eval_context_t*);  // (position[]) -> float[]
static int _position_z(data_t*, data_t[], eval_context_t*);  // (position[]) -> float[]

// Linear algebra
static int _dot           (data_t*, data_t[], eval_context_t*); // (float[], float[]) -> float
static int _cross         (data_t*, data_t[], eval_context_t*); // (float[3], float[3]) -> float[3]
static int _length        (data_t*, data_t[], eval_context_t*); // (float[]) -> float
static int _mat4_mul_mat4 (data_t*, data_t[], eval_context_t*); // (float[4][4], float[4][4]) -> float[4][4]
static int _mat4_mul_vec4 (data_t*, data_t[], eval_context_t*); // (float[4][4], float[4]) -> float[4]

static int _vec2 (data_t*, data_t[], eval_context_t*); // (float, float) -> float[2]
static int _vec3 (data_t*, data_t[], eval_context_t*); // (float, float, float) -> float[4]
static int _vec4 (data_t*, data_t[], eval_context_t*); // (float, float, float, float) -> float[4]

static int _op_simd_neg_farr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)ctx;
    ASSERT(dst);
    int64_t total_count = type_info_element_stride_count(arg[0].type);
    ASSERT(total_count % md_simd_widthf == 0);
    const float* src_arr = as_float_arr(arg[0]);
    float* dst_arr = as_float_arr(*dst);

    for (int64_t i = 0; i < total_count; i += md_simd_widthf) {
        md_simd_typef val = md_simd_loadf(src_arr + i);
        md_simd_storef(dst_arr + i, md_simd_subf(md_simd_zerof(), val));
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
    {CSTR("PI"),     _PI,  {.dim.angle = 1 }},
    {CSTR("TAU"),    _TAU, {.dim.angle = 1 }},
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
    {CSTR("cast"),    TI_BITFIELD,      1,  {TI_BITFIELD_ARR},  _flatten_bf_arr},

    //{CSTR("cast"),    TI_FLOAT3_ARR,    1,  {TI_POSITION_ARR},  _position, FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_DYNAMIC_LENGTH},
};

static procedure_t operators[] = {
    {CSTR("not"),    TI_BOOL,            1,  {TI_BOOL},              _op_not_b},
    {CSTR("or"),     TI_BOOL,            2,  {TI_BOOL,   TI_BOOL},   _op_or_b_b},
    {CSTR("xor"),    TI_BOOL,            2,  {TI_BOOL,   TI_BOOL},   _op_xor_b_b},
    {CSTR("and"),    TI_BOOL,            2,  {TI_BOOL,   TI_BOOL},   _op_and_b_b},

    // BITFIELD NOT
    {CSTR("not"),    TI_BITFIELD,   1,  {TI_BITFIELD},  _not},
    {CSTR("and"),    TI_BITFIELD,   2,  {TI_BITFIELD, TI_BITFIELD}, _and},
    {CSTR("or"),     TI_BITFIELD,   2,  {TI_BITFIELD, TI_BITFIELD}, _or},
    {CSTR("xor"),    TI_BITFIELD,   2,  {TI_BITFIELD, TI_BITFIELD}, _xor},

    // Binary add
    {CSTR("+"),      TI_FLOAT,       2,  {TI_FLOAT,      TI_FLOAT},      _op_add_f_f},
    {CSTR("+"),      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT},      _op_add_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("+"),      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},  _op_add_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {CSTR("+"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_DISTRIBUTION},    _op_simd_add_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("+"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_FLOAT},           _op_simd_add_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("+"),      TI_VOLUME, 2,  {TI_VOLUME,  TI_VOLUME},    _op_simd_add_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("+"),      TI_VOLUME, 2,  {TI_VOLUME,  TI_FLOAT},     _op_simd_add_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("+"),      TI_INT,         2,  {TI_INT,        TI_INT},        _op_add_i_i},
    {CSTR("+"),      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT},        _op_add_iarr_i,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("+"),      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT_ARR},    _op_add_iarr_iarr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    // Unary negation
    {CSTR("-"),      TI_FLOAT,       1,  {TI_FLOAT},                     _op_neg_f},
    {CSTR("-"),      TI_FLOAT_ARR,   1,  {TI_FLOAT_ARR},                 _op_neg_farr,       FLAG_RET_AND_ARG_EQUAL_LENGTH},

    {CSTR("-"),      TI_INT,         1,  {TI_INT},                       _op_neg_i},
    {CSTR("-"),      TI_INT_ARR,     1,  {TI_INT_ARR},                   _op_neg_iarr,       FLAG_RET_AND_ARG_EQUAL_LENGTH},

    {CSTR("-"),      TI_DISTRIBUTION,   1,  {TI_DISTRIBUTION},           _op_simd_neg_farr,     FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("-"),      TI_VOLUME,         1,  {TI_DISTRIBUTION},           _op_simd_neg_farr,     FLAG_RET_AND_ARG_EQUAL_LENGTH},


    // Binary sub
    {CSTR("-"),      TI_FLOAT,       2,  {TI_FLOAT,      TI_FLOAT},      _op_sub_f_f},
    {CSTR("-"),      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT},      _op_sub_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("-"),      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},  _op_sub_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {CSTR("-"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_DISTRIBUTION},    _op_simd_sub_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("-"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_FLOAT},           _op_simd_sub_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("-"),      TI_VOLUME, 2,  {TI_VOLUME,  TI_VOLUME},    _op_simd_sub_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("-"),      TI_VOLUME, 2,  {TI_VOLUME,  TI_FLOAT},     _op_simd_sub_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("-"),      TI_INT,         2,  {TI_INT,        TI_INT},        _op_sub_i_i},
    {CSTR("-"),      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT},        _op_sub_iarr_i,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("-"),      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT_ARR},    _op_sub_iarr_iarr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    // Binary mul
    {CSTR("*"),      TI_FLOAT,       2,  {TI_FLOAT,      TI_FLOAT},      _op_mul_f_f},
    {CSTR("*"),      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT},      _op_mul_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("*"),      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},  _op_mul_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {CSTR("*"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_DISTRIBUTION},    _op_simd_mul_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("*"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_FLOAT},           _op_simd_mul_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("*"),      TI_VOLUME,     2,  {TI_VOLUME,  TI_VOLUME},    _op_simd_mul_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("*"),      TI_VOLUME,     2,  {TI_VOLUME,  TI_FLOAT},     _op_simd_mul_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("*"),      TI_INT,         2,  {TI_INT,        TI_INT},        _op_mul_i_i},
    {CSTR("*"),      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT},        _op_mul_iarr_i,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("*"),      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT_ARR},    _op_mul_iarr_iarr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    // Binary div
    {CSTR("/"),      TI_FLOAT,       2,  {TI_FLOAT,      TI_FLOAT},      _op_div_f_f},
    {CSTR("/"),      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT},      _op_div_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("/"),      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},  _op_div_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {CSTR("/"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_DISTRIBUTION},    _op_simd_div_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("/"),      TI_DISTRIBUTION,   2,  {TI_DISTRIBUTION,  TI_FLOAT},           _op_simd_div_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("/"),      TI_VOLUME, 2,  {TI_VOLUME,  TI_VOLUME},    _op_simd_div_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
    {CSTR("/"),      TI_VOLUME, 2,  {TI_VOLUME,  TI_FLOAT},     _op_simd_div_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},

    {CSTR("/"),      TI_INT,         2,  {TI_INT,        TI_INT},        _op_div_i_i},
    {CSTR("/"),      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT},        _op_div_iarr_i,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {CSTR("/"),      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT_ARR},    _op_div_iarr_iarr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
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

    // ARRAY VERSIONS OF NATIVE FUNCS
    {CSTR("abs"),    TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_fabsf,    FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("floor"),  TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_floorf,   FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {CSTR("ceil"),   TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_ceilf,    FLAG_RET_AND_ARG_EQUAL_LENGTH},

    // VECTORIZED VERSIONS FOR ARRAYS OF BIGGER DATA
    {CSTR("abs"),    TI_VOLUME, 1,  {TI_VOLUME}, _op_simd_abs_farr},

    {CSTR("min"),    TI_FLOAT,  1, {TI_FLOAT_ARR}, _min_farr},
    {CSTR("max"),    TI_FLOAT,  1, {TI_FLOAT_ARR}, _max_farr},    

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
    {CSTR("type"),      TI_BITFIELD, 1, {TI_STRING_ARR},    _name,              FLAG_STATIC_VALIDATION},
    {CSTR("name"),      TI_BITFIELD, 1, {TI_STRING_ARR},    _name,              FLAG_STATIC_VALIDATION},
    {CSTR("label"),     TI_BITFIELD, 1, {TI_STRING_ARR},    _name,              FLAG_STATIC_VALIDATION},
    {CSTR("element"),   TI_BITFIELD, 1, {TI_STRING_ARR},    _element_str,       FLAG_STATIC_VALIDATION},
    {CSTR("element"),   TI_BITFIELD, 1, {TI_IRANGE_ARR},    _element_irng,      FLAG_STATIC_VALIDATION},
    {CSTR("atom"),      TI_BITFIELD, 1, {TI_IRANGE_ARR},    _atom_irng,         FLAG_STATIC_VALIDATION},
    {CSTR("atom"),      TI_BITFIELD, 1, {TI_INT_ARR},       _atom_int,          FLAG_STATIC_VALIDATION},

    // Residue level
    {CSTR("protein"),   TI_BITFIELD_ARR, 0, {0},                _protein,       FLAG_QUERYABLE_LENGTH},
    {CSTR("water"),     TI_BITFIELD_ARR, 0, {0},                _water,         FLAG_QUERYABLE_LENGTH},
    {CSTR("ion"),       TI_BITFIELD_ARR, 0, {0},                _ion,           FLAG_QUERYABLE_LENGTH},
    {CSTR("resname"),   TI_BITFIELD_ARR, 1, {TI_STRING_ARR},    _resname,       FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {CSTR("residue"),   TI_BITFIELD_ARR, 1, {TI_STRING_ARR},    _resname,       FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {CSTR("resid"),     TI_BITFIELD_ARR, 1, {TI_IRANGE_ARR},    _resid,         FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {CSTR("residue"),   TI_BITFIELD_ARR, 1, {TI_IRANGE_ARR},    _residue,       FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},

    // Chain level
    {CSTR("chain"),     TI_BITFIELD_ARR,  1,  {TI_STRING_ARR},    _chain_str,   FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {CSTR("chain"),     TI_BITFIELD_ARR,  1,  {TI_IRANGE_ARR},    _chain_irng,  FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},

    // Dynamic selectors (depend on atomic position, therefore marked as dynamic which means the values cannot be determined at compile-time)
    // Also have variable result (well its a single bitfield, but the number of atoms within is not fixed)
    {CSTR("x"),         TI_BITFIELD, 1, {TI_FRANGE},                    _x,                 FLAG_DYNAMIC},
    {CSTR("y"),         TI_BITFIELD, 1, {TI_FRANGE},                    _y,                 FLAG_DYNAMIC},
    {CSTR("z"),         TI_BITFIELD, 1, {TI_FRANGE},                    _z,                 FLAG_DYNAMIC},
    {CSTR("within"),    TI_BITFIELD, 2, {TI_FLOAT,  TI_POSITION_ARR},   _within_expl_flt,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},
    {CSTR("within"),    TI_BITFIELD, 2, {TI_FRANGE, TI_POSITION_ARR},   _within_expl_frng,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},
    {CSTR("within"),    TI_BITFIELD, 1, {TI_FLOAT},                     _within_impl_flt,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},
    {CSTR("within"),    TI_BITFIELD, 1, {TI_FRANGE},                    _within_impl_frng,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},

    // --- PROPERTY COMPUTE ---
    {CSTR("distance"),          TI_FLOAT,       2,  {TI_POSITION_ARR, TI_POSITION_ARR}, _distance,          FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},
    {CSTR("distance_min"),      TI_FLOAT,       2,  {TI_POSITION_ARR, TI_POSITION_ARR}, _distance_min,      FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},
    {CSTR("distance_max"),      TI_FLOAT,       2,  {TI_POSITION_ARR, TI_POSITION_ARR}, _distance_max,      FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},
    {CSTR("distance_pair"),     TI_FLOAT_ARR,   2,  {TI_POSITION_ARR, TI_POSITION_ARR}, _distance_pair,     FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_QUERYABLE_LENGTH},

    {CSTR("angle"),     TI_FLOAT,   3,  {TI_POSITION_ARR, TI_POSITION_ARR, TI_POSITION_ARR},                    _angle,     FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},
    {CSTR("dihedral"),  TI_FLOAT,   4,  {TI_POSITION_ARR, TI_POSITION_ARR, TI_POSITION_ARR, TI_POSITION_ARR},   _dihedral,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},

    {CSTR("rmsd"),      TI_FLOAT,   1,  {TI_BITFIELD},    _rmsd,     FLAG_DYNAMIC},

    {CSTR("rdf"), TI_DISTRIBUTION, 3, {TI_POSITION_ARR, TI_POSITION_ARR, TI_FLOAT},  _rdf_flt,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},
    {CSTR("rdf"), TI_DISTRIBUTION, 3, {TI_POSITION_ARR, TI_POSITION_ARR, TI_FRANGE}, _rdf_frng,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},

    {CSTR("sdf"),       TI_VOLUME, 3, {TI_BITFIELD_ARR, TI_BITFIELD, TI_FLOAT},      _sdf,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_SDF | FLAG_VISUALIZE},

    // --- GEOMETRICAL OPERATIONS ---
    {CSTR("com"),       TI_FLOAT3,      1,  {TI_POSITION_ARR},  _com,       FLAG_DYNAMIC | FLAG_VISUALIZE | FLAG_STATIC_VALIDATION},
    {CSTR("plane"),     TI_FLOAT4,      1,  {TI_POSITION_ARR},  _plane,     FLAG_DYNAMIC | FLAG_VISUALIZE | FLAG_STATIC_VALIDATION},

    {CSTR("position"),      TI_FLOAT3_ARR,  1,  {TI_POSITION_ARR},  _position,      FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH},
    {CSTR("position_x"),    TI_FLOAT_ARR,   1,  {TI_POSITION_ARR},  _position_x,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH},
    {CSTR("position_y"),    TI_FLOAT_ARR,   1,  {TI_POSITION_ARR},  _position_y,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH},
    {CSTR("position_z"),    TI_FLOAT_ARR,   1,  {TI_POSITION_ARR},  _position_z,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH},

    // --- MISC ---
    {CSTR("flatten"),   TI_BITFIELD,        1,  {TI_BITFIELD_ARR},  _flatten_bf_arr},
    {CSTR("residue"),   TI_BITFIELD_ARR,    1,  {TI_BITFIELD_ARR},  _fill_residue,  FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH},
    {CSTR("chain"),     TI_BITFIELD_ARR,    1,  {TI_BITFIELD_ARR},  _fill_chain,    FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH},
};

#undef CSTR

static inline void visualize_atom_mask(const md_bitfield_t* mask, md_script_visualization_t* vis) {
    ASSERT(vis);
    if (vis->o->flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        md_bitfield_or_inplace(vis->atom_mask, mask);
    }
}

static inline void visualize_atom_range(irange_t range, md_script_visualization_t* vis) {
    ASSERT(vis);
    if (vis->o->flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        md_bitfield_set_range(vis->atom_mask, range.beg, range.end);
    }
}

static inline void visualize_atom_index(int64_t index, md_script_visualization_t* vis) {
    ASSERT(vis);
    if (vis->o->flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        md_bitfield_set_bit(vis->atom_mask, index);
    }
}

static inline void visualize_atom_indices32(const int32_t* indices, int64_t num_indices, md_script_visualization_t* vis) {
    ASSERT(vis);
    if (vis->o->flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        for (int64_t i = 0; i < num_indices; ++i) {
            md_bitfield_set_bit(vis->atom_mask, indices[i]);
        }
    }
}

static inline void visualize_atom_indices64(const int64_t* indices, int64_t num_indices, md_script_visualization_t* vis) {
    ASSERT(vis);
    if (vis->o->flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        for (int64_t i = 0; i < num_indices; ++i) {
            md_bitfield_set_bit(vis->atom_mask, indices[i]);
        }
    }
}

static inline uint16_t push_vertex(vec3_t pos, md_script_visualization_t* vis) {
    ASSERT(vis);
    ASSERT(vis->o->alloc);
    const int64_t idx = vis->vertex.count;
    md_array_push(vis->vertex.pos, pos, vis->o->alloc);
    vis->vertex.count += 1;
    ASSERT(idx < UINT16_MAX);
    return (uint16_t)idx;
}

static inline void push_point(uint16_t v, md_script_visualization_t* vis) {
    ASSERT(vis);
    ASSERT(vis->o->alloc);
    md_array_push(vis->point.idx, v, vis->o->alloc);
    vis->point.count += 1;
}

static inline void push_line(uint16_t v0, uint16_t v1, md_script_visualization_t* vis) {
    ASSERT(vis);
    ASSERT(vis->o->alloc);
    md_array_push(vis->line.idx, v0, vis->o->alloc);
    md_array_push(vis->line.idx, v1, vis->o->alloc);
    vis->line.count += 1;
}

static inline void push_triangle(uint16_t v0, uint16_t v1, uint16_t v2, md_script_visualization_t* vis) {
    ASSERT(vis);
    ASSERT(vis->o->alloc);
    md_array_push(vis->triangle.idx, v0, vis->o->alloc);
    md_array_push(vis->triangle.idx, v1, vis->o->alloc);
    md_array_push(vis->triangle.idx, v2, vis->o->alloc);
    vis->triangle.count += 1;

}

static inline void push_sphere(vec4_t pos_rad, md_script_visualization_t* vis) {
    ASSERT(vis);
    ASSERT(vis->o->alloc);
    md_array_push(vis->sphere.pos_rad, pos_rad, vis->o->alloc);
    vis->sphere.count += 1;
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
    int64_t beg_bit = bitfield->beg_bit;
    int64_t end_bit = bitfield->end_bit;
    while ((beg_bit = md_bitfield_scan(bitfield, beg_bit, end_bit)) != 0) {
        const int64_t idx = beg_bit - 1;
        const float weight = w[idx];
        vec4_t v = {x[idx], y[idx], z[idx], 1.0f};
        v = vec4_mul_f(v, weight);
        sum = vec4_add(sum, v);
    }
    if (sum.w == 0) sum.w = 1.0f;
    return vec3_div_f(vec3_from_vec4(sum), sum.w);
}

static inline vec3_t* extract_vec3(const float* x, const float* y, const float* z, const md_bitfield_t* bitfield, md_allocator_i* alloc) {
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    ASSERT(bitfield);
    ASSERT(alloc);

    vec3_t* pos = 0;

    int64_t beg_bit = bitfield->beg_bit;
    int64_t end_bit = bitfield->end_bit;
    while ((beg_bit = md_bitfield_scan(bitfield, beg_bit, end_bit)) != 0) {
        const int64_t idx = beg_bit - 1;
        vec3_t p = { x[idx], y[idx], z[idx] };
        md_array_push(pos, p, alloc);
    }
    return pos;
}

static inline int64_t extract_xyz(float* dst_x, float* dst_y, float* dst_z, const float* src_x, const float* src_y, const float* src_z, const md_bitfield_t* bitfield) {
    ASSERT(dst_x);
    ASSERT(dst_y);
    ASSERT(dst_z);
    ASSERT(src_x);
    ASSERT(src_y);
    ASSERT(src_z);
    ASSERT(bitfield);

    int64_t count = 0;

    int64_t beg_bit = bitfield->beg_bit;
    int64_t end_bit = bitfield->end_bit;
    while ((beg_bit = md_bitfield_scan(bitfield, beg_bit, end_bit)) != 0) {
        const int64_t idx = beg_bit - 1;
        dst_x[count] = src_x[idx];
        dst_y[count] = src_y[idx];
        dst_z[count] = src_z[idx];
        count += 1;
    }

    return count;
}

static inline int64_t extract_xyzw(float* dst_x, float* dst_y, float* dst_z, float* dst_w, const float* src_x, const float* src_y, const float* src_z, const float* src_w, const md_bitfield_t* bitfield) {
    ASSERT(dst_x);
    ASSERT(dst_y);
    ASSERT(dst_z);
    ASSERT(dst_w);
    ASSERT(src_x);
    ASSERT(src_y);
    ASSERT(src_z);
    ASSERT(src_w);
    ASSERT(bitfield);

    int64_t count = 0;

    int64_t beg_bit = bitfield->beg_bit;
    int64_t end_bit = bitfield->end_bit;
    while ((beg_bit = md_bitfield_scan(bitfield, beg_bit, end_bit)) != 0) {
        const int64_t idx = beg_bit - 1;
        dst_x[count] = src_x[idx];
        dst_y[count] = src_y[idx];
        dst_z[count] = src_z[idx];
        dst_w[count] = src_w[idx];
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

static int64_t* get_residue_indices_in_context(const md_molecule_t* mol, const md_bitfield_t* bitfield, md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(mol->atom.residue_idx);
    ASSERT(mol->residue.atom_range);
    ASSERT(alloc);

    int64_t* arr = 0;

    if (bitfield) {
        int64_t beg = bitfield->beg_bit;
        int64_t end = bitfield->end_bit;
        while ((beg = md_bitfield_scan(bitfield, beg, end)) != 0) {
            const int64_t i = beg - 1;
            const int64_t res_idx = mol->atom.residue_idx[i];
            md_array_push(arr, res_idx, alloc);
            beg = mol->residue.atom_range[res_idx].end;
        }
    }
    else {
        for (int64_t i = 0; i < mol->residue.count; ++i) {
            md_array_push(arr, i, alloc);
        }
    }

    return arr;
}

static int64_t* get_chain_indices_in_context(const md_molecule_t* mol, const md_bitfield_t* bitfield, md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(mol->atom.chain_idx);
    ASSERT(mol->chain.atom_range);
    ASSERT(alloc);
    int64_t* arr = 0;

    if (bitfield) {
        int64_t beg = bitfield->beg_bit;
        int64_t end = bitfield->end_bit;
        while ((beg = md_bitfield_scan(bitfield, beg, end)) != 0) {
            const int64_t i = beg - 1;
            const int64_t chain_idx = mol->atom.chain_idx[i];
            md_array_push(arr, chain_idx, alloc);
            beg = mol->chain.atom_range[chain_idx].end;
        }
    }
    else {
        for (int64_t i = 0; i < mol->chain.count; ++i) {
            md_array_push(arr, i, alloc);
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
        int64_t first = CLAMP((int64_t)mol_ctx->beg_bit, 0, mol->atom.count - 1);
        int64_t last = CLAMP((int64_t)mol_ctx->end_bit - 1, first, mol->atom.count - 1);
        range.beg = mol->atom.residue_idx[first];
        range.end = mol->atom.residue_idx[last] + 1;
    }
    return range;
}

static inline irange_t get_chain_range_in_context(const md_molecule_t* mol, const md_bitfield_t* mol_ctx) {
    ASSERT(mol);
    irange_t range = { 0, (int32_t)mol->chain.count };
    if (mol_ctx) {
        int64_t first = CLAMP((int64_t)mol_ctx->beg_bit, 0, mol->atom.count - 1);
        int64_t last = CLAMP((int64_t)mol_ctx->end_bit - 1, first, mol->atom.count - 1);
        range.beg = mol->atom.chain_idx[first];
        range.end = mol->atom.chain_idx[last] + 1;
    }
    return range;
}

static int position_validate(data_t arg, int arg_idx, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (md_type_info_t)TI_POSITION_ARR));
    if (element_count(arg) == 0) return 0;
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);
    const int64_t ctx_size = ctx_range.end - ctx_range.beg;

    switch (arg.type.base_type) {
    case TYPE_FLOAT:
        ASSERT(is_type_directly_compatible(arg.type, (md_type_info_t)TI_FLOAT3_ARR));
        return (int)type_info_array_len(arg.type);
    case TYPE_INT: {
        int* indices = as_int_arr(arg);
        int count = 0;
        for (int64_t i = 0; i < element_count(arg); ++i) {
            const int64_t idx = (int64_t)ctx_range.beg + (int64_t)indices[i] - 1;

            if (!(ctx_range.beg <= idx && idx < ctx_range.end)) {
                create_error(ctx->ir, ctx->arg_tokens[arg_idx], "supplied index (%i) is not within expected range (%i:%i)",
                    (int)idx, 1, ctx_size);
                return STATIC_VALIDATION_ERROR;
            }
            if (ctx->mol_ctx) {
                if (!md_bitfield_test_bit(ctx->mol_ctx, idx)) {
                    create_error(ctx->ir, ctx->arg_tokens[arg_idx], "supplied index (%i) is not represented within the supplied context", (int)idx);
                    return STATIC_VALIDATION_ERROR;
                }
            }
            count += 1;
        }
        return count;
    }
    case TYPE_BITFIELD:
    {
        if (ctx->arg_flags[arg_idx] & FLAG_DYNAMIC) {
            ASSERT(ctx->backchannel);
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
        if (ctx->arg_flags[arg_idx] & FLAG_DYNAMIC) {
            ASSERT(ctx->backchannel);
            ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
        }
        irange_t* ranges = as_irange_arr(arg);
        int count = 0;
        for (int64_t i = 0; i < element_count(arg); ++i) {
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
                create_error(ctx->ir, ctx->arg_tokens[0], "supplied range (%i:%i) is not within expected range (%i:%i)",
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

static void position_visualize(data_t arg, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(ctx->vis);

    if (element_count(arg) == 0) return;
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    switch (arg.type.base_type) {
    case TYPE_FLOAT:
        ASSERT(is_type_directly_compatible(arg.type, (md_type_info_t)TI_FLOAT3_ARR));
        if (ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
            const vec3_t* pos = as_vec3_arr(arg);
            const int64_t len = element_count(arg);
            for (int64_t i = 0; i < len; ++i) {
                uint16_t v = push_vertex(pos[i], ctx->vis);
                push_point(v, ctx->vis);
            }
        }
        break;
    case TYPE_INT: {
        int* indices = as_int_arr(arg);
        const int64_t len = element_count(arg);
        for (int64_t i = 0; i < len; ++i) {
            const int64_t idx = (int64_t)ctx_range.beg + (int64_t)indices[i] - 1;
            visualize_atom_index(idx, ctx->vis);
        }
        break;
    }
    case TYPE_BITFIELD:
    {
        md_bitfield_t tmp_bf = {0};
        const md_bitfield_t* in_bf = as_bitfield(arg);

        if (ctx->mol_ctx) {
            md_bitfield_init(&tmp_bf, ctx->temp_alloc);
        }

        for (int64_t i = 0; i < element_count(arg); ++i) {
            const md_bitfield_t* bf = &in_bf[i];
            if (ctx->mol_ctx) {
                md_bitfield_and(&tmp_bf, bf, ctx->mol_ctx);
                bf = &tmp_bf;
            }
            visualize_atom_mask(bf, ctx->vis);
        }

        if (ctx->mol_ctx) {
            md_bitfield_free(&tmp_bf);
        }
        break;
    }
    case TYPE_IRANGE:
    {
        irange_t* ranges = as_irange_arr(arg);
        for (int64_t i = 0; i < element_count(arg); ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            if (range_in_range(range, ctx_range)) {
                if (ctx->mol_ctx) {
                    int64_t beg_bit = range.beg;
                    int64_t end_bit = range.end;
                    while ((beg_bit = md_bitfield_scan(ctx->mol_ctx, beg_bit, end_bit)) != 0) {
                        const int64_t idx = beg_bit - 1;
                        visualize_atom_index(idx, ctx->vis);
                    }
                }
                else {
                    visualize_atom_range(range, ctx->vis);
                }
            }
        }
        break;
    }
    default:
        ASSERT(false);
    }
}

static vec3_t* position_extract(data_t arg, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (md_type_info_t)TI_POSITION_ARR));
    vec3_t* positions = 0;
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    switch (arg.type.base_type) {
    case TYPE_FLOAT:
        ASSERT(is_type_directly_compatible(arg.type, (md_type_info_t)TI_FLOAT3_ARR));
        md_array_push_array(positions, as_vec3_arr(arg), element_count(arg), ctx->temp_alloc);
        break;
    case TYPE_INT: {
        int* indices = as_int_arr(arg);
        int64_t num_idx = element_count(arg);
        for (int64_t i = 0; i < num_idx; ++i) {
            // Shift here since we use 1 based indices for atoms
            const int64_t idx = (int64_t)ctx_range.beg + (int64_t)indices[i] - 1;
            vec3_t pos = { ctx->mol->atom.x[idx], ctx->mol->atom.y[idx], ctx->mol->atom.z[idx] };
            md_array_push(positions, pos, ctx->temp_alloc);
        }
        break;
    }
    case TYPE_IRANGE: {
        irange_t* ranges = as_irange_arr(arg);
        int64_t num_ranges = element_count(arg);
        for (int64_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            range = clamp_range(range, ctx_range);

            if (ctx->mol_ctx) {
                int64_t beg_bit = range.beg;
                int64_t end_bit = range.end;
                while ((beg_bit = md_bitfield_scan(ctx->mol_ctx, beg_bit, end_bit)) != 0) {
                    const int64_t idx = beg_bit - 1;
                    vec3_t pos = { ctx->mol->atom.x[idx], ctx->mol->atom.y[idx], ctx->mol->atom.z[idx] };
                    md_array_push(positions, pos, ctx->temp_alloc);
                }
            }
            else {
                for (int64_t j = range.beg; j < range.end; ++j) {
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
        int64_t num_bf = element_count(arg);

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
            for (int64_t i = 0; i < num_bf; ++i) {
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

static vec3_t position_extract_com(data_t arg, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (md_type_info_t)TI_POSITION_ARR));
    vec4_t sum = { 0 };
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    switch (arg.type.base_type) {
    case TYPE_FLOAT: {
        ASSERT(is_type_directly_compatible(arg.type, (md_type_info_t)TI_FLOAT3_ARR));
        vec3_t* positions = as_vec3_arr(arg);
        for (int64_t i = 0; i < element_count(arg); ++i) {
            vec4_add(sum, vec4_from_vec3(positions[i], 1.0f));
        }
        break;
    }
    case TYPE_INT: {
        int* indices = as_int_arr(arg);
        int64_t num_idx = element_count(arg);
        for (int64_t i = 0; i < num_idx; ++i) {
            // Shift here since we use 1 based indices for atoms
            const int64_t idx = (int64_t)ctx_range.beg + (int64_t)indices[i] - 1;
            vec4_t xyzw = { ctx->mol->atom.x[idx], ctx->mol->atom.y[idx], ctx->mol->atom.z[idx], 1.0f };
            xyzw = vec4_mul_f(xyzw, ctx->mol->atom.mass[idx]);
            sum = vec4_add(sum, xyzw);
        }
        break;
    }
    case TYPE_IRANGE: {
        irange_t* ranges = as_irange_arr(arg);
        int64_t num_ranges = element_count(arg);
        for (int64_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            range = clamp_range(range, ctx_range);

            if (ctx->mol_ctx) {
                int64_t beg_bit = range.beg;
                int64_t end_bit = range.end;
                while ((beg_bit = md_bitfield_scan(ctx->mol_ctx, beg_bit, end_bit)) != 0) {
                    const int64_t idx = beg_bit - 1;
                    vec4_t xyzw = { ctx->mol->atom.x[idx], ctx->mol->atom.y[idx], ctx->mol->atom.z[idx], 1.0f };
                    xyzw = vec4_mul_f(xyzw, ctx->mol->atom.mass[idx]);
                    sum = vec4_add(sum, xyzw);
                }
            }
            else {
                for (int64_t j = range.beg; j < range.end; ++j) {
                    vec4_t xyzw = { ctx->mol->atom.x[j], ctx->mol->atom.y[j], ctx->mol->atom.z[j], 1.0f };
                    xyzw = vec4_mul_f(xyzw, ctx->mol->atom.mass[j]);
                    sum = vec4_add(sum, xyzw);
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
        int64_t num_bf = element_count(arg);

        if (num_bf == 1) {
            md_bitfield_t* bf = bf_arr;
            if (ctx->mol_ctx) {
                md_bitfield_init(&tmp_bf, ctx->temp_alloc);
                md_bitfield_and(&tmp_bf, bf_arr, ctx->mol_ctx);
                bf = &tmp_bf;
            }
            vec3_t com = extract_com(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol->atom.mass, bf);
            sum = vec4_from_vec3(com, 1.0f);
        }
        else {
            for (int64_t i = 0; i < num_bf; ++i) {
                md_bitfield_t* bf = &bf_arr[i];
                if (ctx->mol_ctx) {
                    md_bitfield_and(&tmp_bf, bf, ctx->mol_ctx);
                    bf = &tmp_bf;
                }
                vec3_t com = extract_com(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol->atom.mass, bf);
                sum = vec4_add(sum, vec4_from_vec3(com, 1.0f));
            }
        }

        if (ctx->mol_ctx) {
            md_bitfield_free(&tmp_bf);
        }
        break;
    }
    default:
        break;
    }

    return vec3_div_f(vec3_from_vec4(sum), sum.w ? sum.w : 1.0f);
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
    for (int64_t i = 0; i < query.len; ++i) {
        if (query.ptr[i] == '*') joker_count += 1;
    }

    if (joker_count == 0) return compare_str(query, str);

    int64_t first_joker = find_char(query, '*');
    int64_t last_joker = rfind_char(query, '*');

    str_t beg_match = first_joker != 0 ? substr(query, 0, first_joker) : (str_t){0};
    str_t end_match = last_joker  != query.len - 1 ? substr(query, last_joker + 1, -1) : (str_t){0};

    if (!str_empty(beg_match)) {
        if (!compare_str_n(beg_match, str, beg_match.len)) return false;
    }

    if (!str_empty(end_match)) {
        str = substr(str, str.len - end_match.len, end_match.len);
        if (!compare_str(end_match, str)) return false;
    }

    return true;
}

static bool validate_query(str_t query) {
    if (query.len == 0) return false;
    int joker_count = 0;
    for (int64_t i = 0; i < query.len; ++i) {
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
    ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT_ARR));
    (void)ctx;

    const float* src = as_float_arr(arg[0]);
    float min_val = 0;

    if (element_count(arg[0]) > 0) {
        min_val = FLT_MAX;
        for (int64_t i = 0; i < element_count(arg[0]); ++i) {
            min_val = MIN(min_val, src[i]);
        }
    }

    as_float(*dst) = min_val;
    return 0;
}

static int _max_farr  (data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT_ARR));
    (void)ctx;

    const float* src = as_float_arr(arg[0]);
    float max_val = 0;

    if (element_count(arg[0]) > 0) {
        max_val = -FLT_MAX;
        for (int64_t i = 0; i < element_count(arg[0]); ++i) {
            max_val = MAX(max_val, src[i]);
        }
    }

    as_float(*dst) = max_val;
    return 0;
}

static int _not  (data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_BITFIELD));
    (void)ctx;

    md_bitfield_t* bf_dst = dst->ptr;
    const md_bitfield_t* bf_src = arg[0].ptr;

    //ASSERT(bf_dst->num_bits == bf_src->num_bits);

    md_bitfield_not(bf_dst, bf_src, 0, ctx->mol->atom.count);

    if (ctx->mol_ctx) {
        // This is a bit conceptually strange,
        // But if you perform a bitwise negation within a context, only the bits within the context are flipped
        md_bitfield_and_inplace(bf_dst, ctx->mol_ctx); // Store intermediate result in dst
    }

    return 0;
}

static int _and  (data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_BITFIELD));
    (void)ctx;

    md_bitfield_t* bf_dst = dst->ptr;
    const md_bitfield_t* bf_src[2] = {arg[0].ptr, arg[1].ptr};

    md_bitfield_and(bf_dst, bf_src[0], bf_src[1]);
    if (ctx->mol_ctx) {
        md_bitfield_and_inplace(bf_dst, ctx->mol_ctx);
    }

    return 0;
}

static int _or   (data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_BITFIELD));
    (void)ctx;

    md_bitfield_t* bf_dst = dst->ptr;
    md_bitfield_t* bf_src[2] = {arg[0].ptr, arg[1].ptr};

    md_bitfield_or(bf_dst, bf_src[0], bf_src[1]);
    if (ctx->mol_ctx) {
        md_bitfield_and_inplace(bf_dst, ctx->mol_ctx);
    }

    return 0;
}

static int _xor(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_BITFIELD));
    (void)ctx;

    md_bitfield_t* bf_dst = dst->ptr;
    md_bitfield_t* bf_src[2] = { arg[0].ptr, arg[1].ptr };

    md_bitfield_xor(bf_dst, bf_src[0], bf_src[1]);
    if (ctx->mol_ctx) {
        md_bitfield_and_inplace(bf_dst, ctx->mol_ctx);
    }

    return 0;
}

static int _dot(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_FLOAT_ARR));
    (void)ctx;
    float* a = (float*)arg[0].ptr;
    float* b = (float*)arg[1].ptr;
    double res = 0; // Accumulate in double, then cast to float
    for (int64_t i = 0; i < element_count(arg[0]); ++i) {
        res += (double)a[i] * (double)b[i]; 
    }
    as_float(*dst) = (float)res;
    return 0;
}

static int _cross(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT3));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT3));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_FLOAT3));
    (void)ctx;

    float* a = (float*)arg[0].ptr;
    float* b = (float*)arg[1].ptr;
    float* c = (float*)dst->ptr;
    c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = a[2]*b[0] - b[2]*a[0];
    c[2] = a[0]*b[1] - b[0]*a[1];

    return 0;
}

static int _length(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT_ARR));
    (void)ctx;

    float* a = (float*)arg[0].ptr;
    double res = 0; // Accumulate in double, then cast to float
    for (int64_t i = 0; i < element_count(arg[0]); ++i) {
        res += (double)a[i] * (double)a[i];
    }
    as_float(*dst) = (float)sqrt(res);

    return 0;
}


static int _mat4_mul_mat4(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT44));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT44));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_FLOAT44));
    (void)ctx;

    // The type system should already have covered this, we are only reading data we know exists.
    mat4_t* A = (mat4_t*) arg[0].ptr;
    mat4_t* B = (mat4_t*) arg[1].ptr;
    mat4_t* C = (mat4_t*) dst->ptr;

    *C = mat4_mul(*A, *B);
    return 0;
}

static int _mat4_mul_vec4(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type,   (md_type_info_t)TI_FLOAT4));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT44));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_FLOAT4));
    (void)ctx;

    // The type system should already have covered this, we are only reading data we know exists.
    mat4_t* M = (mat4_t*) arg[0].ptr;
    vec4_t* v = (vec4_t*) arg[1].ptr;
    vec4_t* r = (vec4_t*) dst->ptr;

    *r = mat4_mul_vec4(*M, *v);
    return 0;
}

static int _vec2(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT2));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_FLOAT));

    (void)ctx;
    float (*res) = (float*)dst->ptr;
    res[0] = as_float(arg[0]);
    res[1] = as_float(arg[1]);
    return 0;
}

static int _vec3(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT3));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[2].type, (md_type_info_t)TI_FLOAT));

    (void)ctx;
    float (*res) = (float*)dst->ptr;
    res[0] = as_float(arg[0]);
    res[1] = as_float(arg[1]);
    res[2] = as_float(arg[2]);
    return 0;
}

static int _vec4(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT4));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[2].type, (md_type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[3].type, (md_type_info_t)TI_FLOAT));

    (void)ctx;
    float (*res) = (float*)dst->ptr;
    res[0] = as_float(arg[0]);
    res[1] = as_float(arg[1]);
    res[2] = as_float(arg[2]);
    res[3] = as_float(arg[3]);
    return 0;
}

static int _all(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (dst) {
        ASSERT(dst && is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
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
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_STRING_ARR));
    ASSERT(ctx && ctx->mol && ctx->mol->atom.name);

    const str_t* str = as_string_arr(arg[0]);
    const int64_t num_str = element_count(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = (md_bitfield_t*)dst->ptr;

        for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
            if (ctx->mol_ctx && !md_bitfield_test_bit(ctx->mol_ctx, i)) continue;
            
            str_t atom_str = label_to_str(&ctx->mol->atom.name[i]);
            for (int64_t j = 0; j < num_str; ++j) {
                if (match_query(str[j], atom_str)) {
                    md_bitfield_set_bit(bf, i);
                }
            }
        }
    } else {
        // We are only validating the arguments here, making sure that they are represented within the potential context
        for (int64_t j = 0; j < num_str; ++j) {
            if (!validate_query(str[j])) {
                create_error(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' is not a valid query", str[j].len, str[j].ptr);
                return -1;
            }
            bool match = false;
            for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
                if (ctx->mol_ctx && !md_bitfield_test_bit(ctx->mol_ctx, i)) continue;
                str_t atom_str = label_to_str(&ctx->mol->atom.name[i]);
                if (match_query(str[j], atom_str)) {
                    match = true;
                    break;
                }
            }
            if (!match) {
                create_error(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' did not match any atom label within the structure", str[j].len, str[j].ptr);
                return -1;
            }
        }
    }

    return 0;
}

static int _element_str(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_STRING_ARR));
    ASSERT(ctx && ctx->mol && ctx->mol->atom.element);

    uint8_t* elem_idx = 0;

    const int64_t num_str = element_count(arg[0]);
    const str_t* str = as_string_arr(arg[0]);

    for (int64_t i = 0; i < num_str; ++i) {
        md_element_t elem = md_util_element_lookup(str[i]);
        if (elem)
            md_array_push(elem_idx, elem, ctx->temp_alloc);
        else {
            create_error(ctx->ir, ctx->arg_tokens[0], "Failed to map '%.*s' into any Element.", str[i].len, str[i].ptr);
            return -1;
        }
    }
    const int64_t num_elem = md_array_size(elem_idx);
    if (num_elem == 0) {
        create_error(ctx->ir, ctx->arg_tokens[0], "No valid arguments in Element");
        return -1;
    }

    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);

        for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
            if (ctx->mol_ctx && !md_bitfield_test_bit(ctx->mol_ctx, i)) continue;
            for (int64_t j = 0; j < num_elem; ++j) {
                if (elem_idx[j] == ctx->mol->atom.element[i]) {
                    md_bitfield_set_bit(bf, i);
                    //bit_set_idx(result.bits, i);
                    break;
                }
            }
        }
    } else {
        for (int64_t j = 0; j < num_elem; ++j) {
            bool found = false;
            for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
                if (elem_idx[j] == ctx->mol->atom.element[i]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                create_error(ctx->ir, ctx->arg_tokens[0], "Element '%.*s' was not found within structure.", str[j].len, str[j].ptr);
                return -1;
            }
        }
    }

    return 0;
}

static int _element_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_IRANGE_ARR));
    ASSERT(ctx && ctx->mol && ctx->mol->atom.element);

    md_bitfield_t* bf = as_bitfield(*dst);
    irange_t* ranges = as_irange_arr(arg[0]);
    const int64_t num_ranges = element_count(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    if (dst) {
        ASSERT(dst->ptr && is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
        for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
            if (ctx->mol_ctx && !md_bitfield_test_bit(ctx->mol_ctx, i)) continue;
            for (int64_t j = 0; j < num_ranges; ++j) {
                if (idx_in_range(ctx->mol->atom.element[i], ranges[j])) {
                    md_bitfield_set_bit(bf, i);
                    break;
                }
            }
        }
    }
    else {
        for (int64_t j = 0; j < num_ranges; ++j) {
            bool match = false;
            for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
                if (ctx->mol_ctx && !md_bitfield_test_bit(ctx->mol_ctx, i)) continue;
                if (idx_in_range(ctx->mol->atom.element[i], ranges[j])) {
                    match = true;
                    break;
                }
            }
            if (!match) {
                create_error(ctx->ir, ctx->arg_tokens[0], "No element within range (%i:%i) was found within structure.", ranges[j].beg, ranges[j].end);
                return -1;
            }
        }
    }

    return 0;
}

static int coordinate_range(data_t* dst, data_t arg[], eval_context_t* ctx, const float* src_arr) {
    ASSERT(ctx && ctx->mol);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FRANGE));
    const frange_t range = as_frange(arg[0]);
    const md_bitfield_t* src_bf = ctx->mol_ctx;

    if (dst) {
        md_bitfield_t* dst_bf = as_bitfield(*dst);
        if (src_bf) {
            int64_t beg_bit = src_bf->beg_bit;
            int64_t end_bit = src_bf->end_bit;
            while ((beg_bit = md_bitfield_scan(src_bf, beg_bit, end_bit)) != 0) {
                const int64_t i = beg_bit - 1;
                if (range.beg <= src_arr[i] && src_arr[i] <= range.end) {
                    md_bitfield_set_bit(dst_bf, i);
                }
            }
        }
        else {
            for (int64_t i = 0; i < ctx->mol->atom.count; ++i) {
                if (range.beg <= src_arr[i] && src_arr[i] <= range.end) {
                    md_bitfield_set_bit(dst_bf, i);
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
                const int64_t i = beg_bit - 1;
                if (range.beg <= src_arr[i] && src_arr[i] <= range.end) {
                    count += 1;
                }
            }
        }
        else {
            for (int64_t i = 0; i < ctx->mol->atom.count; ++i) {
                if (range.beg <= src_arr[i] && src_arr[i] <= range.end) {
                    count += 1;
                }
            }
        }
        return count;
    }
}

static int _x(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return coordinate_range(dst, arg, ctx, ctx->mol->atom.x);
}

static int _y(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return coordinate_range(dst, arg, ctx, ctx->mol->atom.y);
}

static int _z(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return coordinate_range(dst, arg, ctx, ctx->mol->atom.z);
}

static bool within_float_iter(uint32_t idx, vec3_t pos, void* user_data) {
    (void)pos;
    md_bitfield_t* bf = (md_bitfield_t*)user_data;
    md_bitfield_set_bit(bf, (int64_t)idx);
    return true;
}

static int _within_expl_flt(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x && ctx->mol->atom.y && ctx->mol->atom.z);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_POSITION_ARR));

    const float radius = as_float(arg[0]);

    if (dst || ctx->vis) {
        const vec3_t* in_pos = position_extract(arg[1], ctx);
        const int64_t num_pos = md_array_size(in_pos);

        md_spatial_hash_t spatial_hash = {0};
        md_spatial_hash_args_t args = {
            .coords = {
                .count = ctx->mol->atom.count,
                .x = ctx->mol->atom.x,
                .y = ctx->mol->atom.y,
                .z = ctx->mol->atom.z,
            },
            .cell_ext = radius / 1.5f,
            .alloc = ctx->temp_alloc,
            .temp_alloc = ctx->temp_alloc
        };
        md_spatial_hash_init(&spatial_hash, &args);

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
            md_bitfield_t* dst_bf = as_bitfield(*dst);
            for (int64_t i = 0; i < num_pos; ++i) {
                md_spatial_hash_query(&spatial_hash, in_pos[i], radius, within_float_iter, dst_bf);
            }
        } else {
            ASSERT(ctx->vis);
            position_visualize(arg[1], ctx);
            // Visualize
            md_bitfield_t* dst_bf = ctx->vis->atom_mask;
            for (int64_t i = 0; i < num_pos; ++i) {
                push_sphere(vec4_from_vec3(in_pos[i], radius), ctx->vis);
                md_spatial_hash_query(&spatial_hash, in_pos[i], radius, within_float_iter, dst_bf);
            }
        }
    }
    else {
        if (position_validate(arg[1], 1, ctx) < 0) return -1;
        if (radius <= 0) {
            create_error(ctx->ir, ctx->arg_tokens[0], "The supplied radius is negative or zero, please supply a positive value");
            return -1;
        }
    }

    return 0;
}

static int _within_impl_flt(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x && ctx->mol->atom.y && ctx->mol->atom.z);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FLOAT));

    const float radius = as_float(arg[0]);

    if (dst || ctx->vis) {
        ASSERT(ctx->mol_ctx);

        vec3_t* in_pos = extract_vec3(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol_ctx, ctx->temp_alloc);
        const int64_t num_pos = md_array_size(in_pos);

        md_spatial_hash_t spatial_hash = { 0 };
        md_spatial_hash_args_t args = {
            .coords = {
                .count = ctx->mol->atom.count,
                .x = ctx->mol->atom.x,
                .y = ctx->mol->atom.y,
                .z = ctx->mol->atom.z,
            },
            .cell_ext = radius / 1.5f,
            .alloc = ctx->temp_alloc,
            .temp_alloc = ctx->temp_alloc
        };
        md_spatial_hash_init(&spatial_hash, &args);

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
            md_bitfield_t* dst_bf = as_bitfield(*dst);
            for (int64_t i = 0; i < num_pos; ++i) {
                md_spatial_hash_query(&spatial_hash, in_pos[i], radius, within_float_iter, dst_bf);
            }
        }
        else {
            // Visualize
            ASSERT(ctx->vis);
            visualize_atom_mask(ctx->mol_ctx, ctx->vis);
            md_bitfield_t* dst_bf = ctx->vis->atom_mask;
            for (int64_t i = 0; i < num_pos; ++i) {
                push_sphere(vec4_from_vec3(in_pos[i], radius), ctx->vis);
                md_spatial_hash_query(&spatial_hash, in_pos[i], radius, within_float_iter, dst_bf);
            }
        }
    }
    else {
        if (radius <= 0) {
            create_error(ctx->ir, ctx->arg_tokens[0], "The supplied radius is negative or zero, please supply a positive value");
            return -1;
        }
        if (!ctx->mol_ctx) {
            create_error(ctx->ir, ctx->op_token, "The operation is missing its context");
            return -1;
        }
    }

    return 0;
}

typedef struct within_frng_payload_t {
    md_bitfield_t* bf;
    frange_t rad_range;
    vec3_t ref_pos;
} within_frng_payload_t;

static bool within_frng_iter(uint32_t idx, vec3_t pos, void* user_data) {
    (void)pos;
    within_frng_payload_t* data = (within_frng_payload_t*)user_data;
    float d2 = vec3_distance_squared(data->ref_pos, pos);
    if (data->rad_range.beg * data->rad_range.beg < d2 && d2 < data->rad_range.end * data->rad_range.end) {
        md_bitfield_set_bit(data->bf, (int64_t)idx);
    }
    return true;
}

static int _within_expl_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.z);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FRANGE));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_POSITION_ARR));

    const frange_t rad_range = as_frange(arg[0]);

    if (dst || ctx->vis) {
        const vec3_t* in_pos = as_vec3_arr(arg[1]);
        const int64_t num_pos = element_count(arg[1]);

        md_spatial_hash_t spatial_hash = {0};
        md_spatial_hash_args_t args = {
            .coords = {
                .count = ctx->mol->atom.count,
                .x = ctx->mol->atom.x,
                .y = ctx->mol->atom.y,
                .z = ctx->mol->atom.z,
            },
            .cell_ext = rad_range.end / 1.5f,
            .alloc = ctx->temp_alloc,
            .temp_alloc = ctx->temp_alloc
        };
        md_spatial_hash_init(&spatial_hash, &args);

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
            md_bitfield_t* dst_bf = as_bitfield(*dst);
            within_frng_payload_t payload = {
                .bf = dst_bf,
                .rad_range = rad_range,
            };
            for (int64_t i = 0; i < num_pos; ++i) {
                float  rad = rad_range.end;
                payload.ref_pos = in_pos[i];
                md_spatial_hash_query(&spatial_hash, in_pos[i], rad, within_frng_iter, &payload);
            }
        } else {
            // Visualize
            md_bitfield_t* dst_bf = ctx->vis->atom_mask;
            within_frng_payload_t payload = {
                .bf = dst_bf,
                .rad_range = rad_range,
            };
            for (int64_t i = 0; i < num_pos; ++i) {
                float  rad = rad_range.end;
                payload.ref_pos = in_pos[i];
                md_spatial_hash_query(&spatial_hash, in_pos[i], rad, within_frng_iter, &payload);
                push_sphere(vec4_from_vec3(in_pos[i], rad), ctx->vis);
            }
        }
    }
    else {
        if (position_validate(arg[0], 0, ctx) < 0) return -1;
        if (rad_range.beg < 0 || rad_range.end < rad_range.beg) {
            create_error(ctx->ir, ctx->arg_tokens[0], "The supplied radius range is invalid, ");
            return -1;
        }
    }

    return 0;
}

static int _within_impl_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.z);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_FRANGE));

    const frange_t rad_range = as_frange(arg[0]);

    if (dst || ctx->vis) {
        ASSERT(ctx->mol_ctx);
        vec3_t* in_pos = extract_vec3(ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol_ctx, ctx->temp_alloc);
        const int64_t num_pos = md_array_size(in_pos);

        md_spatial_hash_t spatial_hash = { 0 };
        md_spatial_hash_args_t args = {
            .coords = {
                .count = ctx->mol->atom.count,
                .x = ctx->mol->atom.x,
                .y = ctx->mol->atom.y,
                .z = ctx->mol->atom.z,
            },
            .cell_ext = rad_range.end / 1.5f,
            .alloc = ctx->temp_alloc,
            .temp_alloc = ctx->temp_alloc
        };
        md_spatial_hash_init(&spatial_hash, &args);

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
            md_bitfield_t* dst_bf = as_bitfield(*dst);
            within_frng_payload_t payload = {
                .bf = dst_bf,
                .rad_range = rad_range,
            };
            for (int64_t i = 0; i < num_pos; ++i) {
                float  rad = rad_range.end;
                payload.ref_pos = in_pos[i];
                md_spatial_hash_query(&spatial_hash, in_pos[i], rad, within_frng_iter, &payload);
            }
        }
        else {
            // Visualize
            md_bitfield_t* dst_bf = ctx->vis->atom_mask;
            within_frng_payload_t payload = {
                .bf = dst_bf,
                .rad_range = rad_range,
            };
            for (int64_t i = 0; i < num_pos; ++i) {
                float  rad = rad_range.end;
                payload.ref_pos = in_pos[i];
                md_spatial_hash_query(&spatial_hash, in_pos[i], rad, within_frng_iter, &payload);
                push_sphere(vec4_from_vec3(in_pos[i], rad), ctx->vis);
            }
        }
    }
    else {
        if (rad_range.beg < 0 || rad_range.end < rad_range.beg) {
            create_error(ctx->ir, ctx->arg_tokens[0], "The supplied radius range is invalid, ");
            return -1;
        }
        if (!ctx->mol_ctx) {
            create_error(ctx->ir, ctx->op_token, "The operation is missing its context");
            return -1;
        }
    }

    return 0;
}

static int _atom_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_IRANGE_ARR));
    ASSERT(ctx && ctx->mol);

    const int64_t num_ranges = element_count(arg[0]);
    const irange_t* ranges = as_irange_arr(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    if (dst) {
        ASSERT(dst->ptr && is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);

        for (int64_t i = 0; i < num_ranges; ++i) {
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
        for (int64_t i = 0; i < num_ranges; ++i) {
            irange_t range = ranges[i];
            if (!range_in_range(remap_range_to_context(range, ctx_range), ctx_range)) {
                create_error(ctx->ir, ctx->arg_tokens[0], "supplied range (%i:%i) is not within expected range (%i:%i)",
                    range.beg, range.end, 1, ctx_range.end - ctx_range.beg);
                return -1;
            }
        }
    }

    return 0;
}

static int _atom_int(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_INT_ARR));
    ASSERT(ctx && ctx->mol);

    const int64_t num_indices = element_count(arg[0]);
    const int* indices = as_int_arr(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);

    if (dst) {
        ASSERT(dst->ptr && is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);
        for (int64_t i = 0; i < num_indices; ++i) {
            int32_t idx = remap_index_to_context(indices[i], ctx_range);
            ASSERT(idx_in_range(idx, ctx_range));
            md_bitfield_set_bit(bf, idx);
        }
    }
    else {
        ASSERT(ctx->arg_tokens);
        for (int64_t i = 0; i < num_indices; ++i) {
            if (!idx_in_range(remap_index_to_context(indices[i], ctx_range), ctx_range)) {
                create_error(ctx->ir, ctx->arg_tokens[0], "supplied index (%i) is not within expected range (%i:%i)",
                    indices[i], 1, ctx_range.end - ctx_range.beg);
                return -1;
            }
        }
    }

    return 0;
}

static int _water(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->residue.atom_range && ctx->mol->atom.element);
    (void)arg;

    int result = 0;
    int64_t* res_indices = get_residue_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);
    const md_element_t* elem = ctx->mol->atom.element;

    if (dst) {
        ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD_ARR));
        const int64_t capacity = type_info_array_len(dst->type);
        (void)capacity;

        md_bitfield_t* bf = (md_bitfield_t*)dst->ptr;
        int64_t dst_idx = 0;
        for (int64_t i = 0; i < md_array_size(res_indices); ++i) {
            const md_range_t range = ctx->mol->residue.atom_range[res_indices[i]];
            if (range.end - range.beg == 3) {
                int32_t j = range.beg;
                // Square each element index in the composition and sum up.
                // Water should be: 1*1 + 1*1 + 8*8 = 66
                int magic = elem[j] * elem[j] + elem[j+1] * elem[j+1] + elem[j+2] * elem[j+2];
                if (magic == 66) {
                    //bit_set(result.bits, (uint64_t)range.beg, (uint64_t)range.end - (uint64_t)range.beg);
                    md_bitfield_set_range(&bf[dst_idx++], range.beg, range.end);
                    ASSERT(dst_idx <= capacity);
                }
            }
        }
    } else {
        // Length query
        int count = 0;
        for (int64_t i = 0; i < md_array_size(res_indices); ++i) {
            const md_range_t range = ctx->mol->residue.atom_range[res_indices[i]];
            if (range.end - range.beg == 3) {
                int32_t j = range.beg;
                // Square each element index in the composition and sum up.
                // Water should be: 1*1 + 1*1 + 8*8 = 66
                int magic = elem[j] * elem[j] + elem[j+1] * elem[j+1] + elem[j+2] * elem[j+2];
                if (magic == 66) {
                    ++count;
                }
            }
        }
        result = count;
    }

    md_array_free(res_indices, ctx->temp_alloc);
    
    return result;
}

static int _protein(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->residue.name && ctx->mol->residue.atom_range);
    (void)arg;

    int result = 0;
    int64_t* res_indices = get_residue_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);

    if (dst) {
        ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD_ARR));
        md_bitfield_t* bf = (md_bitfield_t*)dst->ptr;
        const int64_t capacity = type_info_array_len(dst->type);
        int64_t dst_idx = 0;

        for (int64_t i = 0; i < md_array_size(res_indices); ++i) {
            int64_t ri = res_indices[i];
            str_t str = label_to_str(&ctx->mol->residue.name[ri]);
            if (md_util_resname_amino_acid(str)) {
                md_range_t range = ctx->mol->residue.atom_range[ri];
                //bit_set(result.bits, (uint64_t)range.beg, (uint64_t)range.end - (uint64_t)range.beg);
                
                md_bitfield_set_range(&bf[dst_idx++], range.beg, range.end);
                ASSERT(dst_idx <= capacity);
            }
        }
    }
    else {
        int count = 0;
        for (int64_t i = 0; i < md_array_size(res_indices); ++i) {
            int64_t ri = res_indices[i];
            str_t str = label_to_str(&ctx->mol->residue.name[ri]);
            if (md_util_resname_amino_acid(str)) {
                count += 1;
            }
        }
        result = count;
    }

    md_array_free(res_indices, ctx->temp_alloc);

    return result;
}

static int _ion(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->residue.atom_range && ctx->mol->atom.element);
    (void)arg;

    int result = 0;
    int64_t* res_indices = get_residue_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);

    if (dst) {
        ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD_ARR));
        md_bitfield_t* bf = as_bitfield(*dst);
        for (int64_t i = 0; i < md_array_size(res_indices); ++i) {
            const md_range_t range = ctx->mol->residue.atom_range[res_indices[i]];
            // Currently only check for monatomic ions
            // We don't know the charges.
            // We make the assumption that if it is a residue of size 1 and of the correct element then its an ion
            if (range.end - range.beg == 1) {
                md_element_t elem = ctx->mol->atom.element[range.beg];
                switch (elem) {
                case 1:  goto set_bit;  // H+ / H-
                case 3:  goto set_bit;  // Li+
                case 4:  goto set_bit;  // Be2+
                case 7:  goto set_bit;  // N3-
                case 8:  goto set_bit;  // O2-
                case 9:  goto set_bit;  // F-
                case 11: goto set_bit;  // Na+
                case 12: goto set_bit;  // Mg2+
                case 13: goto set_bit;  // Al3+
                case 15: goto set_bit;  // P3-
                case 16: goto set_bit;  // S2-
                case 17: goto set_bit;  // Cl-
                case 19: goto set_bit;  // K+
                case 20: goto set_bit;  // Ca2+
                case 24: goto set_bit;  // Cr2+ / Cr3+
                case 25: goto set_bit;  // Mn2+ / Mn3+
                case 26: goto set_bit;  // Fe2+ / Fe3+
                case 27: goto set_bit;  // Co2+ / Co3+
                case 29: goto set_bit;  // Cu1+ / Cu2+
                case 30: goto set_bit;  // Zn2+
                case 34: goto set_bit;  // Se2-
                case 35: goto set_bit;  // Br-
                case 37: goto set_bit;  // Rb+
                case 38: goto set_bit;  // Sr2+
                case 47: goto set_bit;  // Ag+
                case 48: goto set_bit;  // Cd2+
                case 50: goto set_bit;  // Sn2+ / Sn4+
                case 53: goto set_bit;  // I-
                case 55: goto set_bit;  // Cs+
                case 56: goto set_bit;  // Ba2+
                case 66: goto set_bit;  // Hg2+
                case 68: goto set_bit;  // Pb2+ / Pb4+
                default: continue;
                }
            set_bit:
                md_bitfield_set_bit(bf, range.beg);
            }
        }
    }
    else {
        int count = 0;
        for (int64_t i = 0; i < md_array_size(res_indices); ++i) {
            const md_range_t range = ctx->mol->residue.atom_range[res_indices[i]];
            // Currently only check for monatomic ions
            // We don't know the charges.
            // We make the assumption that if it is a residue of size 1 and of the correct element then its an ion
            if (range.end - range.beg == 1) {
                md_element_t elem = ctx->mol->atom.element[range.beg];
                switch (elem) {
                case 1:  goto count_res;  // H+ / H-
                case 3:  goto count_res;  // Li+
                case 4:  goto count_res;  // Be2+
                case 8:  goto count_res;  // O2-
                case 9:  goto count_res;  // F-
                case 11: goto count_res;  // Na+
                case 12: goto count_res;  // Mg2+
                case 13: goto count_res;  // Al3+
                case 16: goto count_res;  // S2-
                case 17: goto count_res;  // Cl-
                case 19: goto count_res;  // K+
                case 20: goto count_res;  // Ca2+
                case 35: goto count_res;  // Br-
                case 47: goto count_res;  // Ag+
                case 53: goto count_res;  // I-
                case 55: goto count_res;  // Cs+
                case 56: goto count_res;  // Ba2+
                default: continue;
                }
            count_res:
                count += 1;
            }
        }
        result = count;
    }

    md_array_free(res_indices, ctx->temp_alloc);

    return result;
}



static int _residue(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->residue.atom_range);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_IRANGE_ARR));

    int result = 0;

    const int64_t num_ranges = element_count(arg[0]);
    const irange_t* ranges = as_irange_arr(arg[0]);
    // Here we use the implicit range given by the context and use that to select substructures within it
    // The supplied iranges will be used as relative indices into the context
    const irange_t ctx_range = get_residue_range_in_context(ctx->mol, ctx->mol_ctx);
    const int32_t ctx_size = ctx_range.end - ctx_range.beg;

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD_ARR));
        if (dst->ptr) {
            md_bitfield_t* bf_arr = (md_bitfield_t*)dst->ptr;
            const int64_t capacity = type_info_array_len(dst->type);

            int64_t dst_idx = 0;
            for (int64_t i = 0; i < num_ranges; ++i) {
                irange_t range = remap_range_to_context(ranges[i], ctx_range);
                range = clamp_range(range, ctx_range);
                for (int64_t j = range.beg; j < range.end; ++j) {
                    //const uint64_t offset = ctx->mol->residue.atom_range[j].beg;
                    //const uint64_t length = ctx->mol->residue.atom_range[j].end - ctx->mol->residue.atom_range[j].beg;
                    //bit_set(result.bits, offset, length);
                    const md_range_t atom_range = ctx->mol->residue.atom_range[j];
                    ASSERT(dst_idx < capacity);
                    md_bitfield_t* bf = &bf_arr[dst_idx++];
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
        for (int64_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            if (!range_in_range(range, ctx_range)) {
                create_error(ctx->ir, ctx->arg_tokens[0], "supplied range (%i:%i) is not within expected range (%i:%i)",
                    ranges[i].beg, ranges[i].end, 1, ctx_size);
                return -1;
            }
            count += range.end - range.beg;
        }
        result = count;
    }

    return result;
}

static int _fill_residue(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_BITFIELD_ARR));
    const int64_t num_bf = element_count(arg[0]);
    const md_bitfield_t* src_bf = as_bitfield(arg[0]);

    md_bitfield_t tmp_bf = {0};
    md_bitfield_init(&tmp_bf, ctx->temp_alloc);
    for (int64_t i = 0; i < num_bf; ++i) {
        md_bitfield_or_inplace(&tmp_bf, &src_bf[i]);
    }

    int result = 0;
    if (ctx->mol && ctx->mol->atom.residue_idx && ctx->mol->residue.atom_range) {
        if (dst) {
            md_bitfield_t* dst_bf = as_bitfield(*dst);
            const int64_t capacity = element_count(*dst);

            int i = 0;
            int64_t beg = tmp_bf.beg_bit;
            int64_t end = tmp_bf.end_bit;
            while ((beg = md_bitfield_scan(&tmp_bf, beg, end)) != 0) {
                const int64_t res_idx = ctx->mol->atom.residue_idx[beg - 1];
                md_range_t range = ctx->mol->residue.atom_range[res_idx];
                ASSERT(i < capacity);
                md_bitfield_set_range(&dst_bf[i], range.beg, range.end);
                i += 1;
                beg = range.end;
            }
        } else {
            int count = 0;

            int64_t beg = tmp_bf.beg_bit;
            int64_t end = tmp_bf.end_bit;
            while ((beg = md_bitfield_scan(&tmp_bf, beg, end)) != 0) {
                const int64_t res_idx = ctx->mol->atom.residue_idx[beg - 1];
                count += 1;
                // Skip rest of atoms in residue
                beg = ctx->mol->residue.atom_range[res_idx].end;
            }
            result = count;
        }
    }
    return result;
}

static int _fill_chain(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_BITFIELD_ARR));

    const int64_t num_bf = element_count(arg[0]);
    const md_bitfield_t* src_bf = as_bitfield(arg[0]);

    md_bitfield_t tmp_bf = {0};
    md_bitfield_init(&tmp_bf, ctx->temp_alloc);
    for (int64_t i = 0; i < num_bf; ++i) {
        md_bitfield_or_inplace(&tmp_bf, &src_bf[i]);
    }

    int result = 0;
    if (ctx->mol && ctx->mol->atom.chain_idx && ctx->mol->chain.atom_range) {
        if (dst) {
            md_bitfield_t* dst_bf = as_bitfield(*dst);
            const int64_t capacity = element_count(*dst);

            int i = 0;
            int64_t beg = tmp_bf.beg_bit;
            int64_t end = tmp_bf.end_bit;
            while ((beg = md_bitfield_scan(&tmp_bf, beg, end)) != 0) {
                const int64_t chain_idx = ctx->mol->atom.chain_idx[beg - 1];
                md_range_t range = ctx->mol->chain.atom_range[chain_idx];
                ASSERT(i < capacity);
                md_bitfield_set_range(&dst_bf[i], range.beg, range.end);
                i += 1;
                beg = range.end;
            }
            ASSERT(i == element_count(*dst));
        } else {
            int count = 0;

            int64_t beg = tmp_bf.beg_bit;
            int64_t end = tmp_bf.end_bit;
            while ((beg = md_bitfield_scan(&tmp_bf, beg, end)) != 0) {
                const int64_t chain_idx = ctx->mol->atom.chain_idx[beg - 1];
                count += 1;
                // Skip rest of atoms in chain
                beg = ctx->mol->chain.atom_range[chain_idx].end;
            }
            result = count;
        }
    }
    return result;
}

static int _resname(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_STRING_ARR));

    int result = 0;

    if (!ctx->mol->residue.name) {
        create_error(ctx->ir, ctx->arg_tokens[0], "The molecule does not contain any residue names");
        return -1;
    }

    const int64_t num_queries = element_count(arg[0]);
    const str_t* queries = as_string_arr(arg[0]);
    
    // Here we only pick the residues which are represented within the context (by having bits set)
    int64_t* res_indices = get_residue_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);
    const int64_t res_count = md_array_size(res_indices);

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD_ARR));
        if (dst->ptr) {
            md_bitfield_t* bf_arr = (md_bitfield_t*)dst->ptr;
            const int64_t capacity = type_info_array_len(dst->type);

            int64_t dst_idx = 0;

            for (int64_t i = 0; i < res_count; ++i) {
                const int64_t res_idx = res_indices[i];
                for (int64_t j = 0; j < num_queries; ++j) {
                    if (match_query(queries[j], label_to_str(&ctx->mol->residue.name[res_idx]))) {
                        //uint64_t offset = ctx->mol->residue.atom_range[i].beg;
                        //uint64_t length = ctx->mol->residue.atom_range[i].end - ctx->mol->residue.atom_range[i].beg;
                        //bit_set(result.bits, offset, length);
                        const md_range_t atom_range = ctx->mol->residue.atom_range[res_idx];
                        ASSERT(dst_idx < capacity);
                        md_bitfield_t* bf = &bf_arr[dst_idx++];

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
        for (int64_t j = 0; j < num_queries; ++j) {
            if (!validate_query(queries[j])) {
                create_error(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' is not a valid query", queries[j].len, queries[j].ptr);
                return -1;
            }
            bool match = false;
            for (int64_t i = 0; i < res_count; ++i) {
                const int64_t res_idx = res_indices[i];
                if (match_query(queries[j], label_to_str(&ctx->mol->residue.name[res_idx]))) {
                    count += 1;
                    match = true;
                }
            }
            if (!match) {
                create_error(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' did not match any residue within the context", queries[j].len, queries[j].ptr);
                return -1;
            }
        }
        result = count;
    }

    md_array_free(res_indices, ctx->temp_alloc);

    return result;
}

static int _resid(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_IRANGE_ARR));

    int result = 0;

    if (!ctx->mol->residue.id) {
        create_error(ctx->ir, ctx->arg_tokens[0], "The molecule does not contain any residue ids");
        return -1;
    }

    const int64_t   num_rid = element_count(arg[0]);
    const irange_t*     rid = arg[0].ptr;

    // Here we only pick the residues which are represented within the context (by having bits set)
    int64_t* res_indices = get_residue_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD_ARR));
        if (dst->ptr) {
            md_bitfield_t* bf_arr = (md_bitfield_t*)dst->ptr;
            const int64_t capacity = type_info_array_len(dst->type);

            int64_t dst_idx = 0;
            for (int64_t i = 0; i < md_array_size(res_indices); ++i) {
                const int64_t res_idx = res_indices[i];
                for (int64_t j = 0; j < num_rid; ++j) {
                    if (idx_in_range((int)ctx->mol->residue.id[res_idx], rid[j])) {
                        //uint64_t offset = ctx->mol->residue.atom_range[i].beg;
                        //uint64_t length = ctx->mol->residue.atom_range[i].end - ctx->mol->residue.atom_range[i].beg;
                        //bit_set(bf->bits, offset, length);
                        const md_range_t atom_range = ctx->mol->residue.atom_range[res_idx];
                        md_bitfield_t* bf = &bf_arr[dst_idx++];
                        md_bitfield_set_range(bf, atom_range.beg, atom_range.end);
                        if (ctx->mol_ctx) {
                            md_bitfield_and_inplace(bf, ctx->mol_ctx);
                        }

                        ASSERT(dst_idx <= capacity);
                        break;
                    }
                }
            }
        }
    }
    else {
        int count = 0;
        for (int64_t j = 0; j < num_rid; ++j) {
            bool match = false;
            for (int64_t i = 0; i < md_array_size(res_indices); ++i) {
                const int64_t res_idx = res_indices[i];
                if (idx_in_range((int)ctx->mol->residue.id[res_idx], rid[j])) {
                    count += 1;
                    match = true;
                }
            }
            if (!match) {
                // @TODO: Should this just be a soft warning instead?
                create_error(ctx->ir, ctx->arg_tokens[0], "No matching residue id was found within the range (%i:%i)", rid[j].beg, rid[j].end);
                return -1;
            }
        }
        result = count;
    }

    return result;
}

static int _chain_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_IRANGE_ARR));

    if (ctx->mol->chain.count == 0 || !ctx->mol->chain.atom_range) {
        create_error(ctx->ir, ctx->arg_tokens[0], "The molecule does not contain any chains");
        return -1;
    }

    int result = 0;

    const int64_t num_ranges = element_count(arg[0]);
    const irange_t* ranges = as_irange_arr(arg[0]);
    const irange_t ctx_range = get_chain_range_in_context(ctx->mol, ctx->mol_ctx);
    const int32_t ctx_size = ctx_range.end - ctx_range.beg;

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD_ARR));
        if (dst->ptr) {
            md_bitfield_t* bf_arr = (md_bitfield_t*)dst->ptr;
            const int64_t capacity = type_info_array_len(dst->type);

            int64_t dst_idx = 0;
            for (int64_t i = 0; i < num_ranges; ++i) {
                irange_t range = remap_range_to_context(ranges[i], ctx_range);
                for (int64_t j = range.beg; j < range.end; ++j) {
                    //const uint64_t offset = ctx->mol->chain.atom_range[j].beg;
                    //const uint64_t length = ctx->mol->chain.atom_range[j].end - ctx->mol->chain.atom_range[j].beg;
                    //bit_set(result.bits, offset, length);

                    const md_range_t atom_range = ctx->mol->chain.atom_range[j];
                    ASSERT(dst_idx < capacity);
                    md_bitfield_t* bf = &bf_arr[dst_idx++];

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
        for (int64_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            if (!range_in_range(range, ctx_range)) {
                create_error(ctx->ir, ctx->arg_tokens[0], "supplied range (%i:%i) is not within expected range (%i:%i)",
                    ranges[i].beg, ranges[i].end, 1, ctx_size);
                return -1;
            }

            count += range.end - range.beg;
        }
        result = count;
    }

    return result;
}

static int _chain_str(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_STRING_ARR));

    if (ctx->mol->chain.count == 0 || !ctx->mol->chain.id) {
        create_error(ctx->ir, ctx->arg_tokens[0], "The molecule does not contain any chains");
        return -1;
    }

    int result = 0;

    const int64_t num_str = element_count(arg[0]);
    const str_t*  str     = arg[0].ptr;
    int64_t* chain_indices = get_chain_indices_in_context(ctx->mol, ctx->mol_ctx, ctx->temp_alloc);

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (md_type_info_t)TI_BITFIELD_ARR));
        if (dst->ptr) {
            md_bitfield_t* bf_arr = (md_bitfield_t*)dst->ptr;
            const int64_t capacity = type_info_array_len(dst->type);

            int64_t dst_idx = 0;
            for (int64_t i = 0; i < md_array_size(chain_indices); ++i) {
                for (int64_t j = 0; j < num_str; ++j) {
                    if (match_query(str[j], label_to_str(&ctx->mol->chain.id[i]))) {
                        //uint64_t offset = ctx->mol->chain.atom_range[i].beg;
                        //uint64_t length = ctx->mol->chain.atom_range[i].end - ctx->mol->chain.atom_range[i].beg;
                        //bit_set(bf->bits, offset, length);
                        const md_range_t atom_range = ctx->mol->chain.atom_range[i];
                        md_bitfield_t* bf = &bf_arr[dst_idx++];
                        md_bitfield_set_range(bf, atom_range.beg, atom_range.end);
                        if (ctx->mol_ctx) {
                            md_bitfield_and_inplace(bf, ctx->mol_ctx);
                        }

                        ASSERT(dst_idx <= capacity);
                        break;
                    }
                }
            }
        }
    }
    else {
        int count = 0;
        for (int64_t j = 0; j < num_str; ++j) {
            if (!validate_query(str[j])) {
                create_error(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' is not a valid query", str[j].len, str[j].ptr);
                return -1;
            }
            bool match = false;
            for (int64_t i = 0; i < md_array_size(chain_indices); ++i) {
                if (match_query(str[j], label_to_str(&ctx->mol->chain.id[i]))) {
                    match = true;
                    break;
                }
            }
            if (!match) {
                create_error(ctx->ir, ctx->arg_tokens[0], "The string '%.*s' did not match any chain within the structure", str[j].len, str[j].ptr);
                return -1;
            }
            count += 1;
        }
        result = count;
    }

    md_array_free(chain_indices, ctx->temp_alloc);

    return result;
}

// Property Compute

static int _distance(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_POSITION_ARR));

    if (dst || ctx->vis) {
        vec3_t a = position_extract_com(arg[0], ctx);
        vec3_t b = position_extract_com(arg[1], ctx);

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT));            
            if (ctx->frame_header) {
                vec4_t vp = {ctx->frame_header->box[0][0], ctx->frame_header->box[1][1], ctx->frame_header->box[2][2], 1.0f};
                vec4_t va = vec4_from_vec3(a, 0.0f);
                vec4_t vb = vec4_from_vec3(b, 0.0f);
                as_float(*dst) = vec4_periodic_distance(va, vb, vp);
            } else {
                as_float(*dst) = vec3_distance(a, b);
            }
        } else {
            ASSERT(ctx->vis);
            position_visualize(arg[0], ctx);
            position_visualize(arg[1], ctx);

            if (ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                uint16_t va = push_vertex(a, ctx->vis);
                uint16_t vb = push_vertex(b, ctx->vis);
                push_point(va, ctx->vis);
                push_point(vb, ctx->vis);
                push_line(va, vb, ctx->vis);
            }
        }
    }
    else {
        int res_a = position_validate(arg[0], 0, ctx);
        int res_b = position_validate(arg[1], 1, ctx);
        if (res_a < 0) return res_a;
        if (res_b < 0) return res_b;
        if (ctx->backchannel) {
            ctx->backchannel->unit = (md_unit_t){ .dim.length = 1};
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
        return 1;
    }

    return 0;
}

static int _distance_min(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_POSITION_ARR));

    if (dst || ctx->vis) {
        const vec3_t* a_pos = position_extract(arg[0], ctx);
        const vec3_t* b_pos = position_extract(arg[1], ctx);
        const int64_t a_len = md_array_size(a_pos);
        const int64_t b_len = md_array_size(b_pos);

        float min_dist = FLT_MAX;
        int64_t min_i = 0;
        int64_t min_j = 0;
        for (int64_t i = 0; i < a_len; ++i) {
            for (int64_t j = 0; j < b_len; ++j) {
                const float dist = vec3_distance(a_pos[i], b_pos[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT));
            as_float(*dst) = min_dist;
        }
        else {
            ASSERT(ctx->vis);
            position_visualize(arg[0], ctx);
            position_visualize(arg[1], ctx);

            if (ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                if (a_pos && b_pos) {
                    uint16_t va = push_vertex(a_pos[min_i], ctx->vis);
                    uint16_t vb = push_vertex(b_pos[min_j], ctx->vis);
                    push_point(va, ctx->vis);
                    push_point(vb, ctx->vis);
                    push_line(va, vb, ctx->vis);
                }
            }
        }
    }
    else {
        int res_a = position_validate(arg[0], 0, ctx);
        int res_b = position_validate(arg[1], 1, ctx);
        if (res_a < 0) return res_a;
        if (res_b < 0) return res_b;
        if (ctx->backchannel) {
            ctx->backchannel->unit = (md_unit_t){ .dim.length = 1 };
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
        return 1;
    }

    return 0;
}

static int _distance_max(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_POSITION_ARR));

    if (dst || ctx->vis) {
        const vec3_t* a_pos = position_extract(arg[0], ctx);
        const vec3_t* b_pos = position_extract(arg[1], ctx);
        const int64_t a_len = md_array_size(a_pos);
        const int64_t b_len = md_array_size(b_pos);

        float max_dist = 0;
        int64_t max_i = 0;
        int64_t max_j = 0;
        for (int64_t i = 0; i < a_len; ++i) {
            for (int64_t j = 0; j < b_len; ++j) {
                const float dist = vec3_distance(a_pos[i], b_pos[j]);
                if (dist > max_dist) {
                    max_dist = dist;
                    max_i = i;
                    max_j = j;
                }
            }
        }

        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT));
            as_float(*dst) = max_dist;
        } else {
            ASSERT(ctx->vis);
            position_visualize(arg[0], ctx);
            position_visualize(arg[1], ctx);

            if (ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                if (a_pos && b_pos) {
                    uint16_t va = push_vertex(a_pos[max_i], ctx->vis);
                    uint16_t vb = push_vertex(b_pos[max_j], ctx->vis);
                    push_point(va, ctx->vis);
                    push_point(vb, ctx->vis);
                    push_line(va, vb, ctx->vis);
                }
            }
        }
    } else {
        int res_a = position_validate(arg[0], 0, ctx);
        int res_b = position_validate(arg[1], 1, ctx);
        if (res_a < 0) return res_a;
        if (res_b < 0) return res_b;
        if (ctx->backchannel) {
            ctx->backchannel->unit = (md_unit_t){ .dim.length = 1};
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
        return 1;
    }

    return 0;
}

static int _distance_pair(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_POSITION_ARR));
    
    int result = 0;

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (md_type_info_t)TI_FLOAT_ARR));

        vec3_t* a_pos = position_extract(arg[0], ctx);
        vec3_t* b_pos = position_extract(arg[1], ctx);
        const int64_t a_len = md_array_size(a_pos);
        const int64_t b_len = md_array_size(b_pos);

        float* dst_arr = as_float_arr(*dst);
        const int64_t dst_len = element_count(*dst);
        ASSERT(dst_len == a_len * b_len);

        for (int64_t i = 0; i < md_array_size(a_pos); ++i) {
            for (int64_t j = 0; j < md_array_size(b_pos); ++j) {
                const float dist = vec3_distance(a_pos[i], b_pos[j]);
                dst_arr[i * md_array_size(b_pos) + j] = dist;
            }
        }
    } else if (ctx->vis) {
        position_visualize(arg[0], ctx);
        position_visualize(arg[1], ctx);

        if (ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
            vec3_t* a_pos = position_extract(arg[0], ctx);
            vec3_t* b_pos = position_extract(arg[1], ctx);
            const int64_t a_len = md_array_size(a_pos);
            const int64_t b_len = md_array_size(b_pos);

            for (int64_t i = 0; i < a_len; ++i) {
                uint16_t va = push_vertex(a_pos[i], ctx->vis);
                push_point(va, ctx->vis);
                for (int64_t j = 0; j < b_len; ++j) {
                    uint16_t vb = push_vertex(b_pos[j], ctx->vis);
                    push_point(vb, ctx->vis);
                    push_line(va, vb, ctx->vis);
                }
            }
        }
    } else {
        int res0 = position_validate(arg[0], 0, ctx);
        int res1 = position_validate(arg[1], 1, ctx);
        if (res0 < 0) return res0;
        if (res1 < 0) return res1;
        if (ctx->backchannel) {
            ctx->backchannel->unit = (md_unit_t){ .dim.length = 1};
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
        int64_t count = (int64_t)res0 * (int64_t)res1;
        if (count > 1000000) {
            create_error(ctx->ir, ctx->op_token, "The size produced by the operation is %"PRId64", which exceeds the upper limit of 1000000", count);
            return STATIC_VALIDATION_ERROR;
        }
        result = (int)count;
    }

    return result;
}

static int _angle(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (md_type_info_t)TI_POSITION_ARR));
    
    if (dst) {
        const vec3_t a = position_extract_com(arg[0], ctx);
        const vec3_t b = position_extract_com(arg[1], ctx);
        const vec3_t c = position_extract_com(arg[2], ctx);
        const vec3_t v0 = vec3_normalize(vec3_sub(a, b));
        const vec3_t v1 = vec3_normalize(vec3_sub(c, b));

        ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT));
        as_float(*dst) = acosf(vec3_dot(v0, v1));
    } else if (ctx->vis) {
        position_visualize(arg[0], ctx);
        position_visualize(arg[1], ctx);
        position_visualize(arg[2], ctx);

        if (ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
            const vec3_t a = position_extract_com(arg[0], ctx);
            const vec3_t b = position_extract_com(arg[1], ctx);
            const vec3_t c = position_extract_com(arg[2], ctx);

            uint16_t va = push_vertex(a, ctx->vis);
            uint16_t vb = push_vertex(b, ctx->vis);
            uint16_t vc = push_vertex(c, ctx->vis);

            push_point(va, ctx->vis);
            push_point(vb, ctx->vis);
            push_point(vc, ctx->vis);
            push_line(va, vb, ctx->vis);
            push_line(vb, vc, ctx->vis);

            // This is the angle arc
            // @TODO: Insert more vertices here and make it a bit smoother
            uint16_t vd = push_vertex(vec3_lerp(b, a, 0.3f), ctx->vis);
            uint16_t ve = push_vertex(vec3_lerp(b, c, 0.3f), ctx->vis);

            push_triangle(vb, vd, ve, ctx->vis);
        }
    } else {
        // Validate input
        if (position_validate(arg[0], 0, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (position_validate(arg[1], 1, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (position_validate(arg[2], 2, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (ctx->backchannel) {
            ctx->backchannel->unit = (md_unit_t){ .dim.angle = 1};
            ctx->backchannel->value_range = (frange_t){0, PI};
        }
    }

    return 0;
}

static int _dihedral(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(is_type_directly_compatible(arg[3].type, (md_type_info_t)TI_POSITION_ARR));

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT));
        const vec3_t a = position_extract_com(arg[0], ctx);
        const vec3_t b = position_extract_com(arg[1], ctx);
        const vec3_t c = position_extract_com(arg[2], ctx);
        const vec3_t d = position_extract_com(arg[3], ctx);

        as_float(*dst) = (float)dihedral_angle(a,b,c,d);
    }
    else if (ctx->vis) {
        position_visualize(arg[0], ctx);
        position_visualize(arg[1], ctx);
        position_visualize(arg[2], ctx);
        position_visualize(arg[3], ctx);

        if (ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
            const vec3_t a = position_extract_com(arg[0], ctx);
            const vec3_t b = position_extract_com(arg[1], ctx);
            const vec3_t c = position_extract_com(arg[2], ctx);
            const vec3_t d = position_extract_com(arg[3], ctx);

            uint16_t va = push_vertex(a, ctx->vis);
            uint16_t vb = push_vertex(b, ctx->vis);
            uint16_t vc = push_vertex(c, ctx->vis);
            uint16_t vd = push_vertex(d, ctx->vis);

            push_point(va, ctx->vis);
            push_point(vb, ctx->vis);
            push_point(vc, ctx->vis);
            push_point(vd, ctx->vis);

            push_line(va, vb, ctx->vis);
            push_line(vb, vc, ctx->vis);
            push_line(vc, vd, ctx->vis);

            // @TODO: Draw planes and the angle between them
        }
    } else {
        // Validate input
        if (position_validate(arg[0], 0, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (position_validate(arg[1], 1, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (position_validate(arg[2], 2, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (position_validate(arg[3], 3, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (ctx->backchannel) {
            ctx->backchannel->unit = (md_unit_t){.dim.angle = 1};
            ctx->backchannel->value_range = (frange_t){-PI, PI};
        }
    }

    return 0;
}

static int _rmsd(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_equivalent(arg[0].type, (md_type_info_t)TI_BITFIELD));
    ASSERT(ctx && ctx->mol && ctx->mol->atom.mass);

    bool result = 0;

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_FLOAT));
        ASSERT(ctx->initial_configuration.x && ctx->initial_configuration.y && ctx->initial_configuration.z);

        if (dst->ptr) {
            const md_bitfield_t* bf = as_bitfield(arg[0]);
            md_bitfield_t tmp_bf = {0};

            if (ctx->mol_ctx) {
                md_bitfield_init(&tmp_bf, ctx->temp_alloc);
                md_bitfield_and(&tmp_bf, bf, ctx->mol_ctx);
                bf = &tmp_bf;
            }
            const int64_t count = md_bitfield_popcount(bf);

            if (count > 0) {
                const int64_t stride = ROUND_UP(count, md_simd_widthf);
                const int64_t coord_bytes = stride * 7 * sizeof(float);
                float* coord_data = md_alloc(ctx->temp_alloc, coord_bytes);
                float* x0 = coord_data + stride * 0;
                float* y0 = coord_data + stride * 1;
                float* z0 = coord_data + stride * 2;
                float* x  = coord_data + stride * 3;
                float* y  = coord_data + stride * 4;
                float* z  = coord_data + stride * 5;
                float* w  = coord_data + stride * 6;

                extract_xyz (x0, y0, z0,    ctx->initial_configuration.x, ctx->initial_configuration.y, ctx->initial_configuration.z, bf);
                extract_xyzw(x,  y,  z,  w, ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, ctx->mol->atom.mass, bf);

                as_float(*dst) = (float)md_util_compute_rmsd(x0, y0, z0, x, y, z, w, count);
            }
        }

        // No need to free anything, this will be taken care of in the procedure call
    }
    else {
        // Validate args
        // Nothing really to validate, arguments are of type bitfields and if the bitfield is empty, that would be ok, since that would yield a valid rmsd -> 0.
        // @TODO: Possibly warn if a bitfield is empty, but currently we only have hard validation errors...

        // Visualize
        // I don't think we need another visualization than the bitfield highlighting the atoms involved...
    }

    return result;
}


// #################
// ###   CASTS   ###
// #################

static int _cast_int_to_flt(data_t* dst, data_t arg[], eval_context_t* ctx){
    ASSERT(dst && compare_type_info(dst->type, (md_type_info_t)TI_FLOAT));
    ASSERT(compare_type_info(arg[0].type, (md_type_info_t)TI_INT));
    (void)ctx;
    as_float(*dst) = (float)as_int(arg[0]);
    return 0;
}

static int _cast_int_to_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && compare_type_info(dst->type, (md_type_info_t)TI_IRANGE));
    ASSERT(compare_type_info(arg[0].type, (md_type_info_t)TI_INT));
    (void)ctx;
    as_irange(*dst) = (irange_t){as_int(arg[0]), as_int(arg[0])};
    return 0;
}

static int _cast_irng_to_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && compare_type_info(dst->type, (md_type_info_t)TI_FRANGE));
    ASSERT(compare_type_info(arg[0].type, (md_type_info_t)TI_IRANGE));
    (void)ctx;
    as_frange(*dst) = (frange_t){(float)as_irange(arg[0]).beg, (float)as_irange(arg[0]).end};
    return 0;
}

static int _cast_int_arr_to_irng_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_IRANGE_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_INT_ARR));
    (void)ctx;
    const int64_t dst_len = element_count(*dst);
    const int64_t arg_len = element_count(arg[0]);
    ASSERT(dst_len == arg_len);

    irange_t* d = dst->ptr;
    int*      s = arg[0].ptr;

    for (int64_t i = 0; i < dst_len; ++i) {
        d[i].beg = s[i];
        d[i].end = s[i];
    }
    return 0;
}

static int _cast_int_arr_to_flt_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_FLOAT_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_INT_ARR));
    (void)ctx;
    const int64_t dst_len = element_count(*dst);
    const int64_t arg_len = element_count(arg[0]);
    ASSERT(dst_len == arg_len);

    float* d = dst->ptr;
    int*   s = arg[0].ptr;

    for (int64_t i = 0; i < dst_len; ++i) {
        d[i] = (float)s[i]; 
    }
    return 0;
}

static int _cast_irng_arr_to_frng_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (md_type_info_t)TI_FRANGE_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_IRANGE_ARR));
    (void)ctx;
    const int64_t dst_len = element_count(*dst);
    const int64_t arg_len = element_count(arg[0]);
    ASSERT(dst_len == arg_len);

    frange_t* d = dst->ptr;
    irange_t* s = arg[0].ptr;

    for (int64_t i = 0; i < dst_len; ++i) {
        d[i].beg = (float)s[i].beg;
        d[i].end = (float)s[i].end; 
    }
    return 0;
}

static int _cast_int_arr_to_bf(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_INT_ARR));
    ASSERT(arg[0].ptr);

    int result = 0;

    const int32_t* indices = as_int_arr(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);
    const int64_t ctx_size = ctx_range.end - ctx_range.beg;
    const int64_t num_idx = element_count(arg[0]);

    // The idea here is that if we have a context, we only make sure to add the indices which are represented within the context.

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);
       
        for (int64_t i = 0; i < num_idx; ++i) {
            // Shift here since we use 1 based indices for atoms
            const int64_t idx = (int64_t)ctx_range.beg + (int64_t)indices[i] - 1;
            if (ctx->mol_ctx) {
                ASSERT(md_bitfield_test_bit(ctx->mol_ctx, idx));
                // This is pre-checked in the static check bellow in the else statement
            }
            md_bitfield_set_bit(bf, idx);
        }
    } else {
        for (int64_t i = 0; i < num_idx; ++i) {
            const int64_t idx = (int64_t)ctx_range.beg + (int64_t)indices[i] - 1;

            if (!(ctx_range.beg <= idx && idx < ctx_range.end)) {
                create_error(ctx->ir, ctx->arg_tokens[0], "supplied index (%i) is not within expected range (%i:%i)", (int)idx, 1, ctx_size);
                return -1;
            }
            if (ctx->mol_ctx) {
                if (!md_bitfield_test_bit(ctx->mol_ctx, idx)) {
                    create_error(ctx->ir, ctx->arg_tokens[0], "supplied index (%i) is not represented within the supplied context", (int)idx);
                    return -1;
                }
            }
        }
    }

    return result;
}

static int _cast_irng_arr_to_bf(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_IRANGE_ARR));
    ASSERT(arg[0].ptr);

    int result = 0;

    const int64_t num_ranges = element_count(arg[0]);
    const irange_t* ranges = as_irange_arr(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->mol, ctx->mol_ctx);
    const int32_t ctx_size = ctx_range.end - ctx_range.beg;

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
        md_bitfield_t* bf = as_bitfield(*dst);
        for (int64_t i = 0; i < num_ranges; ++i) {
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
        for (int64_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            if (!range_in_range(range, ctx_range)) {
                create_error(ctx->ir, ctx->arg_tokens[0], "supplied range (%i:%i) is not within expected range (%i:%i)",
                    ranges[i].beg, ranges[i].end, 1, ctx_size);
                return -1;
            }
        }
    }

    return result;
}

static int _flatten_bf_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)ctx;
    ASSERT(dst && is_type_equivalent(dst->type, (md_type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_BITFIELD_ARR));

    md_bitfield_t* dst_bf = as_bitfield(*dst);
    const md_bitfield_t* src_bf_arr = (md_bitfield_t*)arg[0].ptr;
    const int64_t count = type_info_array_len(arg[0].type);

    for (int64_t i = 0; i < count; ++i) {
        md_bitfield_or_inplace(dst_bf, &src_bf_arr[i]);
    }

    return 0;
}

static int _com(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));

    if (dst) {
        ASSERT(compare_type_info(dst->type, (md_type_info_t)TI_FLOAT3));
        as_vec3(*dst) = position_extract_com(arg[0], ctx);
    } else if (ctx->vis) {
        position_visualize(arg[0], ctx);
        if (ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
            vec3_t com = position_extract_com(arg[0], ctx);
            vec3_t* pos = position_extract(arg[0], ctx);
            uint16_t c = push_vertex(com, ctx->vis);
            push_point(c, ctx->vis);
            for (int64_t i = 0; i < md_array_size(pos); ++i) {
                uint16_t v = push_vertex(pos[i], ctx->vis);
                push_line(c, v, ctx->vis);
            }
        }
    } else {
        return position_validate(arg[0], 0, ctx);
    }
    return 0;
}

static int64_t extract_xyz_from_arg(vec3_t** xyz, data_t arg, eval_context_t* ctx) {    
    if (is_type_directly_compatible(arg.type, (md_type_info_t)TI_FLOAT3_ARR)) {
        // Nothing to do!
        if (*xyz) *xyz = arg.ptr;
        return element_count(arg);
    }
    else if (is_type_directly_compatible(arg.type, (md_type_info_t)TI_BITFIELD_ARR)) {
        const md_bitfield_t* bf_arr = as_bitfield(arg);
        int64_t count = 0;
        md_bitfield_t tmp_bf = {0};
        md_bitfield_init(&tmp_bf, ctx->temp_alloc);
        for (int64_t i = 0; i < element_count(arg); ++i) {
            const md_bitfield_t* bf = &bf_arr[0];
            if (ctx->mol_ctx) {
                md_bitfield_and(&tmp_bf, bf, ctx->mol_ctx);
                bf = &tmp_bf;
            }
            int64_t popcount = md_bitfield_popcount(bf);
            if (popcount) {
                md_array_ensure(*xyz, md_array_size(*xyz) + popcount, ctx->temp_alloc);
                count += popcount;
                int64_t beg_bit = bf->beg_bit;
                int64_t end_bit = bf->end_bit;
                while ((beg_bit = md_bitfield_scan(bf, beg_bit, end_bit)) != 0) {
                    int64_t idx = beg_bit - 1;
                    vec3_t p = {ctx->mol->atom.x[idx], ctx->mol->atom.y[idx], ctx->mol->atom.z[idx]};
                    md_array_push(*xyz, p, ctx->temp_alloc);
                }
            }
        }
    }
    else if (is_type_directly_compatible(arg.type, (md_type_info_t)TI_INT_ARR)) {
        int32_t* indices = as_int_arr(arg);
        for (int64_t i = 0; i < element_count(arg); ++i) {
            int32_t idx = indices[i];
            vec3_t p = {ctx->mol->atom.x[idx], ctx->mol->atom.y[idx], ctx->mol->atom.z[idx]};
            md_array_push(*xyz, p, ctx->temp_alloc);
        }
    }
    else if (is_type_directly_compatible(arg.type, (md_type_info_t)TI_IRANGE_ARR)) {
        for (int64_t i = 0; i < element_count(arg); ++i) {
            
        }
    }

    // Bad.
    return -1;
}

static int _plane(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));
    (void)ctx;

    const vec3_t* in_pos = position_extract(arg[0], ctx);
    const int64_t num_pos = md_array_size(in_pos);

    if (dst || ctx->vis) {
        ASSERT(num_pos >= 3);
        // Compute eigen vectors from covariance matrix
        // normal should be the third axis (i.e. the smallest eigen value)

        vec3_t com = {0,0,0};
        for (int64_t i = 0; i < num_pos; ++i) {
            com = vec3_add(com, in_pos[i]);
        }
        com = vec3_div_f(com, (float)num_pos);

        vec3_t eigen_vec[3];
        float  eigen_val[3];
        mat3_eigen(mat3_covariance_matrix_vec3(in_pos, com, num_pos), eigen_vec, eigen_val);

        vec3_t normal = vec3_normalize(eigen_vec[2]);
        float d = vec3_dot(normal, com);
        vec4_t plane = {normal.x, normal.y, normal.z, d};

        if (dst) {
            ASSERT(compare_type_info(dst->type, (md_type_info_t)TI_FLOAT4));
            as_vec4(*dst) = plane;
        }

        if (ctx->vis) {
            position_visualize(arg[0], ctx);
            if (ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                vec3_t v = vec3_mul_f(eigen_vec[0], eigen_val[0]);
                vec3_t u = vec3_mul_f(eigen_vec[1], eigen_val[1]);

                // Draw something plane-ish
                uint16_t ci = push_vertex(com, ctx->vis);
                uint16_t ni = push_vertex(vec3_add(com, normal), ctx->vis);

                uint16_t corner[4] = {
                    push_vertex(vec3_sub(vec3_sub(com, v), u), ctx->vis),
                    push_vertex(vec3_sub(vec3_add(com, v), u), ctx->vis),
                    push_vertex(vec3_add(vec3_add(com, v), u), ctx->vis),
                    push_vertex(vec3_add(vec3_sub(com, v), u), ctx->vis),
                };

                push_line(corner[0], corner[1], ctx->vis);
                push_line(corner[1], corner[2], ctx->vis);
                push_line(corner[2], corner[3], ctx->vis);
                push_line(corner[3], corner[0], ctx->vis);

                push_line(ci, ni, ctx->vis);
            }
        }
    } else {
        if (position_validate(arg[0], 0, ctx) < 0) return -1;
        if (num_pos < 3) {
            create_error(ctx->ir, ctx->arg_tokens[0], "Invalid number of positions, need at least 3 to compute a plane");
            return -1;
        }
    }

    return 0;
}

static int _position(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (md_type_info_t)TI_FLOAT3_ARR));
        const vec3_t* pos = position_extract(arg[0], ctx);
        ASSERT(element_count(*dst) == md_array_size(pos));
        memcpy(dst->ptr, pos, md_array_size(pos) * sizeof(vec3_t));
    }
    else if (ctx->vis) {
        position_visualize(arg[0], ctx);
    } else {
        return position_validate(arg[0], 0, ctx);
    }
    return 0;
}

static int _position_x(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (md_type_info_t)TI_FLOAT_ARR));
        const vec3_t* pos = position_extract(arg[0], ctx);
        ASSERT(element_count(*dst) == md_array_size(pos));
        float* x = as_float_arr(*dst);
        for (int64_t i = 0; i < element_count(*dst); ++i) {
            x[i] = pos[i].x;
        }
    }
    else if (ctx->vis) {
        position_visualize(arg[0], ctx);
    } else {
        return position_validate(arg[0], 0, ctx);
    }
    return 0;
}

static int _position_y(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (md_type_info_t)TI_FLOAT_ARR));
        const vec3_t* pos = position_extract(arg[0], ctx);
        ASSERT(element_count(*dst) == md_array_size(pos));
        float* y = as_float_arr(*dst);
        for (int64_t i = 0; i < element_count(*dst); ++i) {
            y[i] = pos[i].y;
        }
    }
    else if (ctx->vis) {
        position_visualize(arg[0], ctx);
    } else {
        return position_validate(arg[0], 0, ctx);
    }
    return 0;
}

static int _position_z(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (md_type_info_t)TI_FLOAT_ARR));
        const vec3_t* pos = position_extract(arg[0], ctx);
        ASSERT(element_count(*dst) == md_array_size(pos));
        float* z = as_float_arr(*dst);
        for (int64_t i = 0; i < element_count(*dst); ++i) {
            z[i] = pos[i].z;
        }
    }
    else if (ctx->vis) {
        position_visualize(arg[0], ctx);
    } else {
        return position_validate(arg[0], 0, ctx);
    }
    return 0;
}

typedef struct {
    vec3_t pos;
    float min_cutoff;
    float max_cutoff;
    float inv_cutoff_range;
    float* bins;
    int32_t num_bins;
    uint32_t idx;
    int32_t total_count;
} rdf_payload_t;

bool rdf_iter(uint32_t idx, vec3_t coord, void* user_param) {
    (void)idx;
    rdf_payload_t* data = user_param;
    const float d = vec3_distance(coord, data->pos);
    if (d == 0) {
        while(0);
    }
    if (data->min_cutoff < d && d < data->max_cutoff) {
        int32_t bin_idx = (int32_t)(((d - data->min_cutoff) * data->inv_cutoff_range) * data->num_bins);
        bin_idx = CLAMP(bin_idx, 0, data->num_bins - 1);
        data->bins[bin_idx] += 1.0f;
        data->total_count += 1;
    }
    return true;
}

#define RDF_BRUTE_FORCE_LIMIT 1000

static inline double spherical_shell_volume(double r_inner, double r_outer) {
    return (4.0f / 3.0f) * PI * (r_outer * r_outer * r_outer - r_inner * r_inner * r_inner);
}

static void compute_rdf(float* bins, int num_bins, const vec3_t* ref_pos, int64_t ref_size, const vec3_t* target_pos, int64_t target_size, float min_cutoff, float max_cutoff, vec3_t pbc_ext, md_allocator_i* alloc) {
    const int64_t sum = ref_size * target_size;
    const float inv_cutoff_range = 1.0f / (max_cutoff - min_cutoff);

    int total_count = 0;

    if (sum < RDF_BRUTE_FORCE_LIMIT) {
        for (int64_t i = 0; i < ref_size; ++i) {
            for (int64_t j = 0; j < target_size; ++j) {

                const float d = vec3_distance(ref_pos[i], target_pos[j]);
                if (d == 0) {
                    while (0);
                }
                if (min_cutoff < d && d < max_cutoff) {
                    int32_t bin_idx = (int32_t)(((d - min_cutoff) * inv_cutoff_range) * num_bins);
                    bin_idx = CLAMP(bin_idx, 0, num_bins - 1);
                    bins[bin_idx] += 1.0f;
                    total_count += 1;
                }
            }
        }
    }
    else {
        md_spatial_hash_args_t args = {
            .alloc = alloc,
            .temp_alloc = alloc,
            .cell_ext = CLAMP(max_cutoff / 1.5f, 1.0f, 16.0f),
            .coords = {
                .count = target_size,
                .stride = sizeof(vec3_t),
                .x = &target_pos->x,
                .y = &target_pos->y,
                .z = &target_pos->z,
            },
        };

        md_spatial_hash_t hash = { 0 };
        md_spatial_hash_init(&hash, &args);

        rdf_payload_t payload = {
            .min_cutoff = min_cutoff,
            .max_cutoff = max_cutoff,
            .inv_cutoff_range = inv_cutoff_range,
            .bins = bins,
            .num_bins = num_bins,
            .total_count = 0,
        };

        if (vec3_dot(pbc_ext, pbc_ext) > 0) {
            vec3_t min_box = { 0, 0, 0 };
            vec3_t max_box = pbc_ext;

            for (int64_t i = 0; i < ref_size; ++i) {
                payload.pos = ref_pos[i];
                payload.idx = (uint32_t)i;
                md_spatial_hash_query_periodic(&hash, ref_pos[i], max_cutoff, min_box, max_box, rdf_iter, &payload);
            }
            total_count = payload.total_count;
        }
        else {
            for (int64_t i = 0; i < ref_size; ++i) {
                payload.pos = ref_pos[i];
                payload.idx = (uint32_t)i;
                md_spatial_hash_query(&hash, ref_pos[i], max_cutoff, rdf_iter, &payload);
            }
            total_count = payload.total_count;
        }
    }

    // Reference rho
    const double ref_rho = total_count / spherical_shell_volume(min_cutoff, max_cutoff);

    // Normalize bins
    // With respect to the nominal distribution rho
    // Each bin is normalized with respect to the nominal distribution
    const double dr = (max_cutoff - min_cutoff) / (float)num_bins;
    for (int64_t i = 0; i < num_bins; ++i) {
        const double r_inner = min_cutoff + (i + 0) * dr;
        const double r_outer = min_cutoff + (i + 1) * dr;
        const double rho = bins[i] / spherical_shell_volume(r_inner, r_outer);
        bins[i] = (float)(rho / ref_rho);
    }
}

static int rdf_visualize(const vec3_t* ref_pos, int64_t ref_size, const vec3_t* target_pos, int64_t target_size, float min_cutoff, float max_cutoff, eval_context_t* ctx) {
    const float min_cutoff2 = min_cutoff * min_cutoff;
    const float max_cutoff2 = max_cutoff * max_cutoff;

    // Visualize
    if (ctx->vis && ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
        const int64_t max_vertices = 10000;
        int64_t vertex_count = 0;

        for (int64_t i = 0; i < ref_size; ++i) {
            uint16_t r_v = push_vertex(ref_pos[i], ctx->vis);
            vertex_count += 1;
            push_point(r_v, ctx->vis);
            for (int64_t j = 0; j < target_size; ++j) {
                const vec3_t d = vec3_sub(ref_pos[i], target_pos[j]);
                const float d2 = vec3_dot(d, d);
                if (min_cutoff2 < d2 && d2 < max_cutoff2) {
                    uint16_t t_v = push_vertex(target_pos[j], ctx->vis);
                    vertex_count += 1;
                    push_line(r_v, t_v, ctx->vis);
                    if (vertex_count > max_vertices) return 0;
                }
            }
        }
    }

    return 0;
}

static int internal_rdf(data_t* dst, data_t arg[], float min_cutoff, float max_cutoff, eval_context_t* ctx) {
    if (dst || ctx->vis) {
        const int num_bins = DIST_BINS;

        const vec3_t* ref_pos = position_extract(arg[0], ctx);
        const int64_t ref_len = md_array_size(ref_pos);

        const vec3_t* trg_pos = position_extract(arg[1], ctx);
        const int64_t trg_len = md_array_size(trg_pos);
        
        if (dst) {
            ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_DISTRIBUTION));
            ASSERT(dst->ptr);
            float* bins = as_float_arr(*dst);
            if (ref_len > 0 && trg_len > 0) {
                vec3_t pbc_ext = { 0,0,0 };
                if (ctx->frame_header) {
                    pbc_ext = (vec3_t){ ctx->frame_header->box[0][0], ctx->frame_header->box[1][1], ctx->frame_header->box[2][2] };
                }
                compute_rdf(bins, num_bins, ref_pos, ref_len, trg_pos, trg_len, min_cutoff, max_cutoff, pbc_ext, ctx->temp_alloc);
            }
        } else {
            return rdf_visualize(ref_pos, ref_len, trg_pos, trg_len, min_cutoff, max_cutoff, ctx);
        }
        return 0;
    }
    else {
        int ref_len = position_validate(arg[0], 0, ctx);
        int trg_len = position_validate(arg[1], 1, ctx);

        // Validate input
        if (ref_len <= 0) {
            create_error(ctx->ir, ctx->arg_tokens[0], "empty reference positions");
            return STATIC_VALIDATION_ERROR;
        }
        if (trg_len <= 0) {
            create_error(ctx->ir, ctx->arg_tokens[1], "empty target positions");
            return STATIC_VALIDATION_ERROR;
        }
        if (min_cutoff < 0.0f || max_cutoff <= min_cutoff) {
            create_error(ctx->ir, ctx->arg_tokens[2], "Invalid cutoff");
            return STATIC_VALIDATION_ERROR;
        }
        if (ctx->arg_flags[2] & FLAG_DYNAMIC) {
            create_error(ctx->ir, ctx->arg_tokens[2], "Cutoff needs to have a static value in order to determine the upper bound if the distribution");
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
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_POSITION_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (md_type_info_t)TI_FLOAT));

    const float cutoff = as_float(arg[2]);
    return internal_rdf(dst, arg, 0, cutoff, ctx);
}

static int _rdf_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_BITFIELD_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_BITFIELD_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (md_type_info_t)TI_FRANGE));

    const frange_t cutoff = as_frange(arg[2]);
    return internal_rdf(dst, arg, cutoff.beg, cutoff.end, ctx);
}

static inline bool are_bitfields_equivalent(const md_bitfield_t bitfields[], int64_t num_bitfields, const md_element_t atom_elements[]) {
    // Number of bits should match.
    // The atomic element of each set bit should match.

    const md_bitfield_t* ref_bf = &bitfields[0];
    const int64_t ref_count = md_bitfield_popcount(ref_bf);

    // We compare against the first one, which is considered to be the reference
    for (int64_t i = 1; i < num_bitfields; ++i) {
        const md_bitfield_t* bf = &bitfields[i];
        const int64_t count = md_bitfield_popcount(bf);
        if (count != ref_count) return false;

        int64_t beg_bit = bf->beg_bit;
        int64_t end_bit = bf->end_bit;
        while ((beg_bit = md_bitfield_scan(bf, beg_bit, end_bit)) != 0) {
            const int64_t idx = beg_bit - 1;
            const int64_t ref_idx = ref_bf->beg_bit + (idx - bf->beg_bit);
            if (!md_bitfield_test_bit(ref_bf, ref_idx)) return false;
            if (atom_elements[idx] != atom_elements[ref_idx]) return false;
        }
    }

    return true;
}

static inline void populate_volume(float* vol, const float* x, const float* y, const float* z, int64_t num_pos, vec4_t max) {
    for (int64_t i = 0; i < num_pos; ++i) {
        const md_simd_f128_t min = md_simd_zero_f128();
        const md_simd_f128_t val = md_simd_set_f128(x[i], y[i], z[i], 0.0f);

        md_simd_f128_t a = md_simd_cmp_lt_f128(min, val);
        md_simd_f128_t b = md_simd_cmp_lt_f128(val, max.f128);
        md_simd_f128_t c = md_simd_and_f128(a, b);

        int res = _mm_movemask_ps(c);
        if (res == 0x7) {
            vol[(int)z[i] * VOL_DIM * VOL_DIM + (int)y[i] * VOL_DIM + (int)x[i]] += 1.0f;
        }

        /*
        if ((0 < x[i] && x[i] < max.x) &&
            (0 < y[i] && y[i] < max.y) &&
            (0 < z[i] && z[i] < max.z)) {
            
        }
        */
    }
}

static inline mat4_t compute_volume_matrix(float radius) {
    // We have the cutoff as a radius, meaning our volume radius has the length 'r' for each axis.
    // Thus this means the diameter is 2*r.
    // We have the resolution VOL_DIM for each axis, meaning each voxel has the extent of 2*r / VOL_RES units
    const float voxel_ext = (2*radius) / VOL_DIM;
    const float s = 1.0f / voxel_ext;
    mat4_t S = mat4_scale(s, s, s);

    const float t = (VOL_DIM / 2);
    mat4_t T = mat4_translate(t, t, t);

    return mat4_mul(T, S);
}

static int _sdf(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_BITFIELD_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[2].type, (md_type_info_t)TI_FLOAT));
    ASSERT(arg[0].ptr);
    ASSERT(arg[1].ptr);
    ASSERT(arg[2].ptr);

    const int64_t num_ref_bitfields = element_count(arg[0]);
    const md_bitfield_t* ref_bitfields = as_bitfield(arg[0]);
    const md_bitfield_t* target_bitfield = as_bitfield(arg[1]);
    float cutoff = as_float(arg[2]);

    const md_bitfield_t* ref_bf = &ref_bitfields[0];

    // We currently only support bitfields which represent equivalent structures. 
    // his each bitfield should contain the same number of atoms, and the n'th atom within the bitfield should match the n'th atom in all other bitfields.
    // This may be loosened in the future to allow aligning structures which are not strictly equivalent.
    int result = 0;

    if (dst || ctx->vis) {
        ASSERT(num_ref_bitfields > 0);
        ASSERT(ctx->initial_configuration.x);
        ASSERT(ctx->initial_configuration.y);
        ASSERT(ctx->initial_configuration.z);

        float* vol = 0;
        if (dst) {
            ASSERT(dst->ptr);
            ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_VOLUME));
            vol = as_float_arr(*dst);
        }

        const int64_t target_size = md_bitfield_popcount(target_bitfield);
        ASSERT(target_size > 0);

        const int64_t ref_size = md_bitfield_popcount(ref_bf);
        ASSERT(ref_size > 0);

        const int64_t mem_size = sizeof(float) * (ref_size * 7 + target_size * 6);
        float* mem = md_alloc(ctx->temp_alloc, mem_size);

        float* ref_x  = mem + ref_size * 0;
        float* ref_y  = mem + ref_size * 1;
        float* ref_z  = mem + ref_size * 2;
        float* ref_x0 = mem + ref_size * 3;
        float* ref_y0 = mem + ref_size * 4;
        float* ref_z0 = mem + ref_size * 5;
        float* ref_w  = mem + ref_size * 6;
        float* target_x  = mem + ref_size * 7 + target_size * 0;
        float* target_y  = mem + ref_size * 7 + target_size * 1;
        float* target_z  = mem + ref_size * 7 + target_size * 2;
        float* target_x0 = mem + ref_size * 7 + target_size * 3;
        float* target_y0 = mem + ref_size * 7 + target_size * 4;
        float* target_z0 = mem + ref_size * 7 + target_size * 5;

        // Fetch initial reference positions
        extract_xyzw(ref_x0, ref_y0, ref_z0, ref_w, ctx->initial_configuration.x, ctx->initial_configuration.y, ctx->initial_configuration.z, ctx->mol->atom.mass, ref_bf);

        // Fetch target positions
        extract_xyz(target_x0, target_y0, target_z0, ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, target_bitfield);

        const vec3_t ref_com0 = md_util_compute_com(ref_x0, ref_y0, ref_z0, ref_w, ref_size);

        // A for alignment matrix, Align eigen vectors with axis x,y,z etc.
        mat3_t eigen_vecs;
        vec3_t eigen_vals;
        mat3_eigen(mat3_covariance_matrix(ref_x0, ref_y0, ref_z0, ref_com0, ref_size), eigen_vecs.col, eigen_vals.elem);
        mat4_t A = mat4_from_mat3(mat3_transpose(eigen_vecs));

        // V for volume matrix scale and align with the volume which we aim to populate with density
        mat4_t V = compute_volume_matrix(cutoff);
        mat4_t VA = mat4_mul(V, A);

        if (ctx->vis && ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_SDF) {
            ctx->vis->sdf.count = num_ref_bitfields;
            ctx->vis->sdf.extent = cutoff;
        }

        const vec4_t max_ext = {VOL_DIM, VOL_DIM, VOL_DIM, VOL_DIM};

        for (int64_t i = 0; i < num_ref_bitfields; ++i) {
            const md_bitfield_t* bf = &ref_bitfields[i];

            extract_xyz(ref_x, ref_y, ref_z, ctx->mol->atom.x, ctx->mol->atom.y, ctx->mol->atom.z, bf);
            vec3_t ref_com = md_util_compute_com(ref_x, ref_y, ref_z, ref_w, ref_size);
            mat3_t R = md_util_compute_optimal_rotation(ref_x0, ref_y0, ref_z0, ref_com0, ref_x, ref_y, ref_z, ref_com, ref_w, ref_size);
            mat4_t RT = mat4_mul(mat4_from_mat3(R), mat4_translate(-ref_com.x, -ref_com.y, -ref_com.z));

            if (vol) {
                mat4_t VART = mat4_mul(VA, RT);
                mat4_batch_transform2(target_x, target_y, target_z, target_x0, target_y0, target_z0, 1.0f, target_size, VART);
                populate_volume(vol, target_x, target_y, target_z, target_size, max_ext);
            }
            if (ctx->vis && ctx->vis->o->flags & MD_SCRIPT_VISUALIZE_SDF) {
                md_allocator_i* alloc = ctx->vis->o->alloc;
                md_bitfield_t bf_cpy = {0};
                md_bitfield_init(&bf_cpy, alloc);
                md_bitfield_copy(&bf_cpy, bf);

                mat4_t ART = mat4_mul(A, RT);
                md_array_push(ctx->vis->sdf.matrices, ART, alloc);
                md_array_push(ctx->vis->sdf.structures, bf_cpy, alloc);
            }
        }

        md_free(ctx->temp_alloc, mem, mem_size);

    } else {
        // Validation
        if (num_ref_bitfields < 1) {
            create_error(ctx->ir, ctx->arg_tokens[0], "Number of bitfields which serve as reference frame must be 1 or more");
            return -1;
        }

        const int64_t ref_bit_count = md_bitfield_popcount(ref_bf);
        if (ref_bit_count <= 0) {
            create_error(ctx->ir, ctx->arg_tokens[0], "The supplied reference bitfield(s) are empty");
            return -1;
        }

        // Test for equivalence
        if (!are_bitfields_equivalent(ref_bitfields, num_ref_bitfields, ctx->mol->atom.element)) {
            create_error(ctx->ir, ctx->arg_tokens[0], "The supplied reference bitfields are not identical: the number of atoms and their corresponding elements do not match between all supplied bitfields");
            return -1;
        }

        const int64_t target_bit_count = md_bitfield_popcount(target_bitfield);
        if (target_bit_count <= 0) {
            create_error(ctx->ir, ctx->arg_tokens[1], "The supplied target bitfield is empty");
            return -1;
        }

        if (ctx->backchannel) {
            ctx->backchannel->unit = (md_unit_t){ .dim = {.count = 1, .length = -3} };
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
    }

    return result;
}


/*
// This is some experimental future work, for matching structures using maximum overlapping subgraph
typedef struct node_t {
    md_element_t elem;
} node_t;

typedef struct edge_t {
    int32_t idx[2];
} edge_t;

typedef struct graph_t {
    node_t* nodes;
    edge_t* edges;
} graph_t;

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
    ASSERT(is_type_directly_compatible(arg[0].type, (md_type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[1].type, (md_type_info_t)TI_BITFIELD));

    if (dst) {
        ASSERT(dst->ptr);
        ASSERT(is_type_equivalent(dst->type, (md_type_info_t)TI_INT_ARR));
    }
}
*/

#undef ANY_LENGTH
//#undef ANY_LEVEL
