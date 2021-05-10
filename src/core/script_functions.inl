#ifndef __MD_SCRIPT_FUNCTIONS_INL__
#define __MD_SCRIPT_FUNCTIONS_INL__


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

#define as_bitrange(arg) (*((bitrange_t*)((arg).ptr)))
#define as_bitrange_arr(arg) ((bitrange_t*)((arg).ptr))

#define as_bitfield(arg) (*((bitfield_t*)((arg).ptr)))

#define as_vec3(arg) (*((vec3*)((arg).ptr)))
#define as_vec3_arr(arg) (((vec3*)((arg).ptr)))

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
        *((base_type*)dst->ptr) = op *((base_type*)arg[0].ptr); \
        return 0; \
    }

#define BAKE_OP_UNARY_M(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        for (int64_t i = 0; i < element_count(*dst); ++i) { \
            ((base_type*)dst->ptr)[i] = op ((base_type*)arg[0].ptr)[i]; \
        } \
        return 0; \
    }

#define BAKE_OP_S_S(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        *((base_type*)dst->ptr) = *((base_type*)arg[0].ptr) op *((base_type*)arg[1].ptr); \
        return 0; \
    }

#define BAKE_OP_M_S(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        for (int64_t i = 0; i < element_count(*dst); ++i) { \
            ((base_type*)dst->ptr)[i] = ((base_type*)arg[0].ptr)[i] op *((base_type*)arg[1].ptr); \
        } \
        return 0; \
    }

#define BAKE_OP_M_M(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        for (int64_t i = 0; i < element_count(*dst); ++i) { \
            ((base_type*)dst->ptr)[i] = ((base_type*)arg[0].ptr)[i] op ((base_type*)arg[1].ptr)[i]; \
        } \
        return 0; \
    }

BAKE_OP_S_S(_op_and_b_b, &&,  bool)
BAKE_OP_S_S(_op_or_b_b,  ||,  bool)
BAKE_OP_UNARY_S(_op_not_b, !, bool)

BAKE_OP_S_S(_op_add_f_f, +, float)
BAKE_OP_S_S(_op_sub_f_f, -, float)
BAKE_OP_S_S(_op_mul_f_f, *, float)
BAKE_OP_S_S(_op_div_f_f, /, float)
BAKE_OP_UNARY_S(_op_neg_f, -, float)
BAKE_OP_UNARY_M(_op_neg_farr, -, float)

BAKE_OP_M_S(_op_add_farr_f, +, float)
BAKE_OP_M_S(_op_sub_farr_f, -, float)
BAKE_OP_M_S(_op_mul_farr_f, *, float)
BAKE_OP_M_S(_op_div_farr_f, /, float)

BAKE_OP_M_M(_op_add_farr_farr, +, float)
BAKE_OP_M_M(_op_sub_farr_farr, -, float)
BAKE_OP_M_M(_op_mul_farr_farr, *, float)
BAKE_OP_M_M(_op_div_farr_farr, /, float)

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

// Forward declarations of functions
// @TODO: Add your declarations here

// Casts (Implicit ones)

static int _cast_int_to_flt             (data_t*, data_t[], eval_context_t*);
static int _cast_int_to_irng            (data_t*, data_t[], eval_context_t*);
static int _cast_irng_to_frng           (data_t*, data_t[], eval_context_t*);
static int _cast_int_arr_to_flt_arr     (data_t*, data_t[], eval_context_t*);
static int _cast_irng_arr_to_frng_arr   (data_t*, data_t[], eval_context_t*);

static int _cast_bf_to_atom     (data_t*, data_t[], eval_context_t*);
static int _cast_bf_to_residue  (data_t*, data_t[], eval_context_t*);
static int _cast_bf_to_chain    (data_t*, data_t[], eval_context_t*);

// Logical operators for custom types
static int _not  (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _and  (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _or   (data_t*, data_t[], eval_context_t*); // -> bitfield

// Selectors
// Atomic level selectors
static int _all     (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _x       (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _y       (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _z       (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_flt  (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _within_frng (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _name    (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _element_str (data_t*, data_t[], eval_context_t*); // -> bitfield
static int _element_irng(data_t*, data_t[], eval_context_t*); // -> bitfield
static int _atom    (data_t*, data_t[], eval_context_t*);   // -> bitfield
//static int _water   (data_t*, data_t[], eval_context_t*); // -> bitfield

// Residue level selectors
static int _resname (data_t*, data_t[], eval_context_t*); // (str[]) -> bitfield
static int _resid   (data_t*, data_t[], eval_context_t*); // (irange[]) -> bitfield   // irange also covers integers, since integers are implicitly convertible to irange
static int _residue (data_t*, data_t[], eval_context_t*); // (irange[]) -> bitfield

// Chain level selectors
static int _chain_irng  (data_t*, data_t[], eval_context_t*); // (irange[]) -> bitfield
static int _chain_str   (data_t*, data_t[], eval_context_t*); // (str[]) -> bitfield

// Property Compute
static int _distance        (data_t*, data_t[], eval_context_t*); // (float[3], float[3])   -> float
static int _distance_min    (data_t*, data_t[], eval_context_t*); // (float[3][], float[3][]) -> float
static int _distance_max    (data_t*, data_t[], eval_context_t*); // (float[3][], float[3][]) -> float
static int _distance_pair   (data_t*, data_t[], eval_context_t*); // (float[3][], float[3][]) -> float

static int _angle   (data_t*, data_t[], eval_context_t*); // (float[3], float[3], float[3]) -> float
static int _dihedral(data_t*, data_t[], eval_context_t*); // (float[3], float[3], float[3], float[3]) -> float

static int _rmsd    (data_t*, data_t[], eval_context_t*); // (float[3][]) -> float

static int _rdf     (data_t*, data_t[], eval_context_t*); // (bitfield, bitfield, float) -> float[128] (Histogram). The idea is that we use a fixed amount of bins, then we let the user choose some kernel to smooth it.
static int _sdf     (data_t*, data_t[], eval_context_t*); // (bitfield, bitfield, float) -> float[128][128][128]. This one cannot be stored explicitly as one copy per frame, but is rather accumulated.

// Geometric operations
static int _com     (data_t*, data_t[], eval_context_t*); // (float[3][]) -> float
static int _plane   (data_t*, data_t[], eval_context_t*);  // (float[3][]) -> float[4]

static int _position_int    (data_t*, data_t[], eval_context_t*);   // (int)      -> float[3]
static int _position_irng   (data_t*, data_t[], eval_context_t*);   // (irange)   -> float[3][]
static int _position_bf     (data_t*, data_t[], eval_context_t*);   // (bitfield) -> float[3][]

static int _ref_frame_bf    (data_t*, data_t[], eval_context_t*);   // (bitfield) -> float[4][4][]

// Linear algebra
static int _dot           (data_t*, data_t[], eval_context_t*); // (float[], float[]) -> float
static int _cross         (data_t*, data_t[], eval_context_t*); // (float[3], float[3]) -> float[3]
static int _mat4_mul_mat4 (data_t*, data_t[], eval_context_t*); // (float[4][4], float[4][4]) -> float[4][4]
static int _mat4_mul_vec4 (data_t*, data_t[], eval_context_t*); // (float[4][4], float[4]) -> float[4]

static int _vec2 (data_t*, data_t[], eval_context_t*); // (float, float) -> float[2]
static int _vec3 (data_t*, data_t[], eval_context_t*); // (float, float, float) -> float[4]
static int _vec4 (data_t*, data_t[], eval_context_t*); // (float, float, float, float) -> float[4]


// This is to mark that the procedure supports a varying length
#define ANY_LENGTH -1
#define ANY_LEVEL -1

#define DIST_BINS 128
#define VOL_DIM 128

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

#define TI_BITRANGE     {TYPE_BITRANGE, {1}, 0}
#define TI_BITRANGE_ARR {TYPE_BITRANGE, {ANY_LENGTH}, 0}

#define TI_BITFIELD         {TYPE_BITFIELD, {1}, 0, ANY_LEVEL}
#define TI_BITFIELD_ATOM    {TYPE_BITFIELD, {1}, 0, LEVEL_ATOM}
#define TI_BITFIELD_RESIDUE {TYPE_BITFIELD, {1}, 0, LEVEL_RESIDUE}
#define TI_BITFIELD_CHAIN   {TYPE_BITFIELD, {1}, 0, LEVEL_CHAIN}

// Predefined constants
static const float _PI  = 3.14159265358f;
static const float _TAU = 6.28318530718f;
static const float _E   = 2.71828182845f;

// @TODO: Add your values here
static identifier_t constants[] = {
    {{"PI", 2},     {TI_FLOAT, (void*)(&_PI),    sizeof(float)}},
    {{"TAU", 3},    {TI_FLOAT, (void*)(&_TAU),   sizeof(float)}},
    {{"E", 1},      {TI_FLOAT, (void*)(&_E),     sizeof(float)}},
};

// IMPLICIT CASTS/CONVERSIONS
static procedure_t casts[] = {
    {{"cast", 4},    TI_FLOAT,              1,  {TI_INT},           _cast_int_to_flt},
    {{"cast", 4},    TI_IRANGE,             1,  {TI_INT},           _cast_int_to_irng},
    {{"cast", 4},    TI_FRANGE,             1,  {TI_IRANGE},        _cast_irng_to_frng},
    {{"cast", 4},    TI_FLOAT_ARR,          1,  {TI_INT_ARR},       _cast_int_arr_to_flt_arr,   FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {{"cast", 4},    TI_FRANGE_ARR,         1,  {TI_IRANGE_ARR},    _cast_irng_arr_to_frng_arr, FLAG_RET_AND_ARG_EQUAL_LENGTH},

    {{"cast", 4},    TI_BITFIELD_ATOM,      1,  {TI_BITFIELD},          _cast_bf_to_atom},
    {{"cast", 4},    TI_BITFIELD_RESIDUE,   1,  {TI_BITFIELD_CHAIN},    _cast_bf_to_residue},

    // This does the heavy lifting for implicitly converting every compatible argument into a position (vec3) if the procedure is marked with FLAG_POSITION
    {{"extract pos", 11},   TI_FLOAT3_ARR,  1,  {TI_INT_ARR},       _position_int,  FLAG_DYNAMIC | FLAG_POSITION | FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {{"extract pos", 11},   TI_FLOAT3_ARR,  1,  {TI_IRANGE_ARR},    _position_irng, FLAG_DYNAMIC | FLAG_POSITION | FLAG_QUERYABLE},
    {{"extract pos", 11},   TI_FLOAT3_ARR,  1,  {TI_BITFIELD},      _position_bf,   FLAG_DYNAMIC | FLAG_POSITION | FLAG_QUERYABLE},
};

static procedure_t operators[] = {
    {{"not", 3},    TI_BOOL,            1,  {TI_BOOL},              _op_not_b},
    {{"or", 2},     TI_BOOL,            2,  {TI_BOOL,   TI_BOOL},   _op_or_b_b},
    {{"and", 3},    TI_BOOL,            2,  {TI_BOOL,   TI_BOOL},   _op_and_b_b},

    // BITFIELD NOT
    {{"not", 3},    TI_BITFIELD_ATOM,   1,  {TI_BITFIELD_ATOM},     _not},
    {{"not", 3},    TI_BITFIELD_RESIDUE,1,  {TI_BITFIELD_RESIDUE},  _not},
    {{"not", 3},    TI_BITFIELD_CHAIN,  1,  {TI_BITFIELD_CHAIN},    _not},

    // BITFIELD AND -> MAINTAIN THE HIGHEST LEVEL OF CONTEXT AMONG OPERANDS
    {{"and", 3},    TI_BITFIELD_ATOM,   2,  {TI_BITFIELD_ATOM,      TI_BITFIELD_ATOM},      _and},
    {{"and", 3},    TI_BITFIELD_RESIDUE,2,  {TI_BITFIELD_RESIDUE,   TI_BITFIELD_RESIDUE},   _and},
    {{"and", 3},    TI_BITFIELD_CHAIN,  2,  {TI_BITFIELD_CHAIN,     TI_BITFIELD_CHAIN},     _and},

    {{"and", 3},    TI_BITFIELD_RESIDUE,2,  {TI_BITFIELD_ATOM,      TI_BITFIELD_RESIDUE},   _and,   FLAG_SYMMETRIC_ARGS},
    {{"and", 3},    TI_BITFIELD_CHAIN,  2,  {TI_BITFIELD_ATOM,      TI_BITFIELD_CHAIN},     _and,   FLAG_SYMMETRIC_ARGS},
    {{"and", 3},    TI_BITFIELD_CHAIN,  2,  {TI_BITFIELD_RESIDUE,   TI_BITFIELD_CHAIN},     _and,   FLAG_SYMMETRIC_ARGS},

    {{"or", 2},     TI_BITFIELD_ATOM,   2,  {TI_BITFIELD_ATOM,      TI_BITFIELD_ATOM},      _or},
    {{"or", 2},     TI_BITFIELD_RESIDUE,2,  {TI_BITFIELD_RESIDUE,   TI_BITFIELD_RESIDUE},   _or},
    {{"or", 2},     TI_BITFIELD_CHAIN,  2,  {TI_BITFIELD_CHAIN,     TI_BITFIELD_CHAIN},     _or},

    {{"or", 2},     TI_BITFIELD_ATOM,   2,  {TI_BITFIELD_ATOM,      TI_BITFIELD_RESIDUE},   _or,    FLAG_SYMMETRIC_ARGS},
    {{"or", 2},     TI_BITFIELD_ATOM,   2,  {TI_BITFIELD_ATOM,      TI_BITFIELD_CHAIN},     _or,    FLAG_SYMMETRIC_ARGS},
    {{"or", 2},     TI_BITFIELD_RESIDUE,2,  {TI_BITFIELD_RESIDUE,   TI_BITFIELD_CHAIN},     _or,    FLAG_SYMMETRIC_ARGS},

    // Binary add
    {{"+", 1},      TI_FLOAT,       2,  {TI_FLOAT,      TI_FLOAT},      _op_add_f_f},
    {{"+", 1},      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT},      _op_add_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {{"+", 1},      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},  _op_add_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {{"+", 1},      TI_INT,         2,  {TI_INT,        TI_INT},        _op_add_i_i},
    {{"+", 1},      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT},        _op_add_iarr_i,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {{"+", 1},      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT_ARR},    _op_add_iarr_iarr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    // Unary negation
    {{"-", 1},      TI_FLOAT,       1,  {TI_FLOAT},                     _op_neg_f},
    {{"-", 1},      TI_FLOAT_ARR,   1,  {TI_FLOAT_ARR},                 _op_neg_farr,       FLAG_RET_AND_ARG_EQUAL_LENGTH},

    {{"-", 1},      TI_INT,         1,  {TI_INT},                       _op_neg_i},
    {{"-", 1},      TI_INT_ARR,     1,  {TI_INT_ARR},                   _op_neg_iarr,       FLAG_RET_AND_ARG_EQUAL_LENGTH},

    // Binary sub
    {{"-", 1},      TI_FLOAT,       2,  {TI_FLOAT,      TI_FLOAT},      _op_sub_f_f},
    {{"-", 1},      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT},      _op_sub_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {{"-", 1},      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},  _op_sub_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {{"-", 1},      TI_INT,         2,  {TI_INT,        TI_INT},        _op_sub_i_i},
    {{"-", 1},      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT},        _op_sub_iarr_i,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {{"-", 1},      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT_ARR},    _op_sub_iarr_iarr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    // Binary mul
    {{"*", 1},      TI_FLOAT,       2,  {TI_FLOAT,      TI_FLOAT},      _op_mul_f_f},
    {{"*", 1},      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT},      _op_mul_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {{"*", 1},      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},  _op_mul_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {{"*", 1},      TI_INT,         2,  {TI_INT,        TI_INT},        _op_mul_i_i},
    {{"*", 1},      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT},        _op_mul_iarr_i,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {{"*", 1},      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT_ARR},    _op_mul_iarr_iarr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    // Binary div
    {{"/", 1},      TI_FLOAT,       2,  {TI_FLOAT,      TI_FLOAT},      _op_div_f_f},
    {{"/", 1},      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT},      _op_div_farr_f,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {{"/", 1},      TI_FLOAT_ARR,   2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},  _op_div_farr_farr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {{"/", 1},      TI_INT,         2,  {TI_INT,        TI_INT},        _op_div_i_i},
    {{"/", 1},      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT},        _op_div_iarr_i,     FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_SYMMETRIC_ARGS},
    {{"/", 1},      TI_INT_ARR,     2,  {TI_INT_ARR,    TI_INT_ARR},    _op_div_iarr_iarr,  FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_ARGS_EQUAL_LENGTH},
};

static procedure_t procedures[] = {
    // NATIVE FUNCS
    {{"sqrt", 4},   TI_FLOAT, 1, {TI_FLOAT}, _sqrtf},
    {{"cbrt", 4},   TI_FLOAT, 1, {TI_FLOAT}, _cbrtf},
    {{"abs", 3},    TI_FLOAT, 1, {TI_FLOAT}, _fabsf},
    {{"floor", 5},  TI_FLOAT, 1, {TI_FLOAT}, _floorf},
    {{"ceil", 4},   TI_FLOAT, 1, {TI_FLOAT}, _ceilf},
    {{"cos", 3},    TI_FLOAT, 1, {TI_FLOAT}, _cosf},
    {{"sin", 3},    TI_FLOAT, 1, {TI_FLOAT}, _sinf},
    {{"asin", 4},   TI_FLOAT, 1, {TI_FLOAT}, _asinf},
    {{"acos", 4},   TI_FLOAT, 1, {TI_FLOAT}, _acosf},
    {{"atan", 4},   TI_FLOAT, 1, {TI_FLOAT}, _atanf},
    {{"log", 3},    TI_FLOAT, 1, {TI_FLOAT}, _logf},
    {{"exp", 3},    TI_FLOAT, 1, {TI_FLOAT}, _expf},
    {{"log2", 4},   TI_FLOAT, 1, {TI_FLOAT}, _log2f},
    {{"exp2", 4},   TI_FLOAT, 1, {TI_FLOAT}, _exp2f},
    {{"log10", 5},  TI_FLOAT, 1, {TI_FLOAT}, _log10f},

    {{"atan", 4},   TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _atan2f},
    {{"atan2", 5},  TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _atan2f},
    {{"pow", 3},    TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _powf},

    // ARRAY VERSIONS OF NATIVE FUNCS
    {{"abs", 3},    TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_fabsf,    FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {{"floor", 5},  TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_floorf,   FLAG_RET_AND_ARG_EQUAL_LENGTH},
    {{"ceil", 4},   TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_ceilf,    FLAG_RET_AND_ARG_EQUAL_LENGTH},

    // LINEAR ALGEBRA
    {{"dot", 3},    TI_FLOAT,   2, {TI_FLOAT_ARR,   TI_FLOAT_ARR},  _dot},
    {{"cross", 5},  TI_FLOAT3,  2, {TI_FLOAT3,      TI_FLOAT3},     _cross},
    {{"mul", 3},    TI_FLOAT44, 2, {TI_FLOAT44,     TI_FLOAT44},    _mat4_mul_mat4},
    {{"mul", 3},    TI_FLOAT4,  2, {TI_FLOAT44,     TI_FLOAT4},     _mat4_mul_vec4},

    // CONSTRUCTORS
    {{"vec2", 4},   TI_FLOAT2,  2, {TI_FLOAT, TI_FLOAT},                       _vec2},
    {{"vec3", 4},   TI_FLOAT3,  3, {TI_FLOAT, TI_FLOAT, TI_FLOAT},             _vec3},
    {{"vec4", 4},   TI_FLOAT4,  4, {TI_FLOAT, TI_FLOAT, TI_FLOAT, TI_FLOAT},   _vec4},

    // BITFIELD CONVERSION (Explicit casts of bitfield types)
    {{"atoms", 5},      TI_BITFIELD_ATOM, 1, {TI_BITFIELD}, _cast_bf_to_atom},
    {{"residues", 8},   TI_BITFIELD_ATOM, 1, {TI_BITFIELD}, _cast_bf_to_residue},
    {{"chains", 6},     TI_BITFIELD_ATOM, 1, {TI_BITFIELD}, _cast_bf_to_chain},

    // --- SELECTION ---

    // Atom level
    {{"all", 3},    TI_BITFIELD_ATOM, 0, {0},               _all},
    {{"type", 4},   TI_BITFIELD_ATOM, 1, {TI_STRING_ARR},   _name,              FLAG_QUERYABLE},
    {{"name", 4},   TI_BITFIELD_ATOM, 1, {TI_STRING_ARR},   _name,              FLAG_QUERYABLE},
    {{"label", 5},  TI_BITFIELD_ATOM, 1, {TI_STRING_ARR},   _name,              FLAG_QUERYABLE},
    {{"element", 7},TI_BITFIELD_ATOM, 1, {TI_STRING_ARR},   _element_str,       FLAG_QUERYABLE},
    {{"element", 7},TI_BITFIELD_ATOM, 1, {TI_IRANGE_ARR},   _element_irng,      FLAG_QUERYABLE},
    {{"atom", 4},   TI_BITFIELD_ATOM, 1, {TI_IRANGE_ARR},   _atom,              FLAG_QUERYABLE},

    // Residue level
    {{"resname", 7},TI_BITFIELD_RESIDUE, 1, {TI_STRING_ARR},    _resname,       FLAG_QUERYABLE},
    {{"residue", 7},TI_BITFIELD_RESIDUE, 1, {TI_STRING_ARR},    _resname,       FLAG_QUERYABLE},
    {{"resid", 5},  TI_BITFIELD_RESIDUE, 1, {TI_IRANGE_ARR},    _resid,         FLAG_QUERYABLE},
    {{"residue", 7},TI_BITFIELD_RESIDUE, 1, {TI_IRANGE_ARR},    _residue,       FLAG_QUERYABLE},

    // Chain level
    {{"chain", 5},  TI_BITFIELD_CHAIN,  1,  {TI_STRING_ARR},    _chain_str,     FLAG_QUERYABLE},
    {{"chain", 5},  TI_BITFIELD_CHAIN,  1,  {TI_IRANGE_ARR},    _chain_irng,    FLAG_QUERYABLE},

    // Dynamic selectors (depend on atomic position, therefore marked as dynamic which means the values cannot be determined at compile-time)
    {{"x", 1},      TI_BITFIELD_ATOM, 1, {TI_FRANGE},  _x, FLAG_DYNAMIC},
    {{"y", 1},      TI_BITFIELD_ATOM, 1, {TI_FRANGE},  _y, FLAG_DYNAMIC},
    {{"z", 1},      TI_BITFIELD_ATOM, 1, {TI_FRANGE},  _z, FLAG_DYNAMIC},
    {{"within", 6}, TI_BITFIELD_ATOM, 2, {TI_FLOAT,  TI_FLOAT3_ARR},      _within_flt,    FLAG_DYNAMIC | FLAG_POSITION},
    {{"within", 6}, TI_BITFIELD_ATOM, 2, {TI_FRANGE, TI_FLOAT3_ARR},      _within_frng,   FLAG_DYNAMIC | FLAG_POSITION},

    // --- PROPERTY COMPUTE ---
    {{"distance", 8},       TI_FLOAT,   2,  {TI_FLOAT3, TI_FLOAT3},         _distance,      FLAG_DYNAMIC | FLAG_POSITION},
    {{"distance_min", 12},  TI_FLOAT,   2,  {TI_FLOAT3_ARR, TI_FLOAT3_ARR}, _distance_min,  FLAG_DYNAMIC | FLAG_POSITION},
    {{"distance_max", 12},  TI_FLOAT,   2,  {TI_FLOAT3_ARR, TI_FLOAT3_ARR}, _distance_max,  FLAG_DYNAMIC | FLAG_POSITION},
    {{"distance_pair", 13}, TI_FLOAT,   2,  {TI_FLOAT3_ARR, TI_FLOAT3_ARR}, _distance_pair, FLAG_RET_AND_ARG_EQUAL_LENGTH | FLAG_DYNAMIC | FLAG_POSITION},

    {{"angle", 5},      TI_FLOAT,   3,  {TI_FLOAT3, TI_FLOAT3, TI_FLOAT3},              _angle,     FLAG_DYNAMIC | FLAG_POSITION},
    {{"dihedral", 8},   TI_FLOAT,   4,  {TI_FLOAT3, TI_FLOAT3, TI_FLOAT3, TI_FLOAT3},   _dihedral,  FLAG_DYNAMIC | FLAG_POSITION},

    {{"rmsd", 4},       TI_FLOAT,   1,  {TI_FLOAT3_ARR},    _rmsd,     FLAG_DYNAMIC | FLAG_POSITION},

    {{"rdf", 3},        TI_DISTRIBUTION,    3,  {TI_FLOAT3_ARR, TI_FLOAT3_ARR, TI_FLOAT},   _rdf,  FLAG_DYNAMIC | FLAG_POSITION},
    {{"sdf", 3},        TI_VOLUME,          3,  {TI_BITFIELD_RESIDUE, TI_FLOAT3_ARR, TI_FLOAT},  _sdf,  FLAG_DYNAMIC | FLAG_POSITION},

    // --- GEOMETRICAL OPERATIONS ---
    {{"com", 3},        TI_FLOAT3,  1,  {TI_FLOAT3_ARR},    _com,           FLAG_DYNAMIC | FLAG_POSITION},
    {{"plane", 5},      TI_FLOAT4,  1,  {TI_FLOAT3_ARR},    _plane,         FLAG_DYNAMIC | FLAG_POSITION},
    //{{"referenceframe", 14}, TI_FLOAT44, 1, {TI_FLOAT3_ARR}, _ref_frame_bf, FLAG_DYNAMIC | FLAG_POSITION},
};

// Tables
enum {
    Unknown = 0,
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V,
    Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru,
    Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb,
    Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh,
    Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og
};

static const char* element_symbols[] = {
    "Xx", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
    "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og" };

static const char* amino_acids[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU",
    "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASC", "GLX", "XLE",
};

static const char* dna[] = {"DA", "DA3", "DA5", "DC", "DC3", "DC5", "DG", "DG3", "DG5", "DT", "DT3", "DT5"};
static const char* acidic[] = { "ASP", "GLU" };
static const char* basic[] = { "ARG", "HIS", "LYS" };

static const char* neutral[] = { "VAL", "PHE", "GLN", "TYR", "HIS", "CYS", "MET", "TRP", "ASX", "GLX", "PCA", "HYP" };
static const char* water[] = { "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP", "TIP2", "TIP3", "TIP4" };
static const char* hydrophobic[] = { "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "CYX" };

static inline uint8_t element_symbol_to_index(str_t str) {
    if (str.len == 1 || str.len == 2) {
        for (uint8_t i = 0; i < ARRAY_SIZE(element_symbols); ++i) {
            const char* sym = element_symbols[i];
            if ( (str.ptr[0] != sym[0]) || (str.len == 2 && str.ptr[1] != sym[1])) continue;
            return i;
        }
    }
    return 0;
}

static inline bool compare_str_cstr(str_t str, const char* cstr) {
    if (!str.ptr || !str.len || !cstr) return false;
    for (uint64_t i = 0; i < str.len; ++i) {
        if (cstr[i] == '\0' || str.ptr[i] != cstr[i]) return false;
    }
    return cstr[str.len] == '\0';
}

static inline bool name_exists_in_array(str_t name, const char* arr[], uint64_t arr_len) {
    for (uint64_t i = 0; i < arr_len; ++i) {
        if (compare_str_cstr(name, arr[i])) return true;
    }
    return false;
}

static inline bool in_range(irange_t range, int idx) {
    // I think a range should be inclusive in this context... Since we should be 1 based and not 0 based on indices
    return range.beg <= idx && idx <= range.end;
}

typedef struct vec3 {
    float x, y, z;
} vec3;

static inline vec3 vec3_add(vec3 a, vec3 b) {
    return (vec3){a.x+b.x, a.y+b.y, a.z+b.z};
}

static inline vec3 vec3_sub(vec3 a, vec3 b) {
    return (vec3){a.x-b.x, a.y-b.y, a.z-b.z};
}

static inline vec3 vec3_mul(vec3 a, vec3 b) {
    return (vec3){a.x*b.x, a.y*b.y, a.z*b.z};
}

static inline vec3 vec3_div(vec3 a, vec3 b) {
    return (vec3){a.x/b.x, a.y/b.y, a.z/b.z};
}

static inline float vec3_dot(vec3 a, vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline float vec3_len(vec3 v) {
    return sqrtf(vec3_dot(v, v));
}

static inline float vec3_dist(vec3 a, vec3 b) {
    return vec3_len(vec3_sub(a, b));
}

static inline vec3 vec3_normalize(vec3 v) {
    const float len = vec3_len(v);
    return (vec3) {
        .x = v.x / len,
        .y = v.y / len,
        .z = v.z / len
    };
}


// IMPLEMENTATIONS
// @TODO: Add more here

static int _not  (data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD));
    (void)ctx;

    bitfield_t* bf_dst = dst->ptr;
    bitfield_t* bf_src = arg[0].ptr;

    ASSERT(bf_dst->num_bits == bf_src->num_bits);

    bit_not(bf_dst->bits, bf_src->bits, 0, bf_dst->num_bits);
    return 0;
}

static int _and  (data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_BITFIELD));
    (void)ctx;

    bitfield_t* bf_dst = dst->ptr;
    bitfield_t* bf_src[2] = {arg[0].ptr, arg[1].ptr};

    ASSERT(bf_dst->num_bits == bf_src[0]->num_bits);
    ASSERT(bf_dst->num_bits == bf_src[1]->num_bits);

    bit_and(bf_dst->bits, bf_src[0]->bits, bf_src[1]->bits, 0, bf_dst->num_bits);
    return 0;
}

static int _or   (data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_BITFIELD));
    (void)ctx;

    bitfield_t* bf_dst = dst->ptr;
    bitfield_t* bf_src[2] = {arg[0].ptr, arg[1].ptr};

    ASSERT(bf_dst->num_bits == bf_src[0]->num_bits);
    ASSERT(bf_dst->num_bits == bf_src[1]->num_bits);

    bit_or(bf_dst->bits, bf_src[0]->bits, bf_src[1]->bits, 0, bf_dst->num_bits);
    return 0;
}

static int _dot(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT_ARR));
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
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_FLOAT3));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT3));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT3));
    (void)ctx;

    float* a = (float*)arg[0].ptr;
    float* b = (float*)arg[1].ptr;
    float* c = (float*)dst->ptr;
    c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = a[2]*b[0] - b[2]*a[0];
    c[2] = a[0]*b[1] - b[0]*a[1];

    return 0;
}

static int _mat4_mul_mat4(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_FLOAT44));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT44));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT44));
    (void)ctx;

    // The type system should already have covered this, we are only reading data we know exists.
    md_mat4* A = (md_mat4*) arg[0].ptr;
    md_mat4* B = (md_mat4*) arg[1].ptr;
    md_mat4* C = (md_mat4*) dst->ptr;

    *C = md_mat4_mul(*A, *B);
    return 0;
}

static int _mat4_mul_vec4(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type,   (type_info_t)TI_FLOAT4));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT44));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT4));
    (void)ctx;

    // The type system should already have covered this, we are only reading data we know exists.
    md_mat4* M = (md_mat4*) arg[0].ptr;
    md_vec4* v = (md_vec4*) arg[1].ptr;
    md_vec4* r = (md_vec4*) dst->ptr;

    *r = md_mat4_mul_vec4(*M, *v);
    return 0;
}

static int _vec2(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_FLOAT2));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT));

    (void)ctx;
    float (*res) = (float*)dst->ptr;
    res[0] = as_float(arg[0]);
    res[1] = as_float(arg[1]);
    return 0;
}

static int _vec3(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_FLOAT3));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FLOAT));

    (void)ctx;
    float (*res) = (float*)dst->ptr;
    res[0] = as_float(arg[0]);
    res[1] = as_float(arg[1]);
    res[2] = as_float(arg[2]);
    return 0;
}

static int _vec4(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_FLOAT4));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[3].type, (type_info_t)TI_FLOAT));

    (void)ctx;
    float (*res) = (float*)dst->ptr;
    res[0] = as_float(arg[0]);
    res[1] = as_float(arg[1]);
    res[2] = as_float(arg[2]);
    res[3] = as_float(arg[3]);
    return 0;
}

static int _all(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_ATOM));
    (void)arg;
    (void)ctx;
    bitfield_t result = *((bitfield_t*)dst->ptr);
    bit_clear(result.bits, 0, result.num_bits);
    bit_set(result.bits, ctx->mol_ctx.atom.beg, ctx->mol_ctx.atom.end - ctx->mol_ctx.atom.beg);
    return 0;
}

static int _name(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));
    ASSERT(ctx && ctx->mol && ctx->mol->atom.name);

    const str_t* str = as_string_arr(arg[0]);
    const uint64_t num_str = (uint64_t)element_count(arg[0]);

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_ATOM));
        bitfield_t result = *((bitfield_t*)dst->ptr);

        for (uint64_t i = ctx->mol_ctx.atom.beg; i < ctx->mol_ctx.atom.end; ++i) {
            const char* atom_str = ctx->mol->atom.name[i];
            for (uint64_t j = 0; j < num_str; ++j) {
                if (compare_str_cstr(str[j], atom_str)) {
                    bit_set_idx(result.bits, i);
                }
            }
        }
    } else {
        int count = 0;
        for (uint64_t i = ctx->mol_ctx.atom.beg; i < ctx->mol_ctx.atom.end; ++i) {
            const char* atom_str = ctx->mol->atom.name[i];
            for (uint64_t j = 0; j < num_str; ++j) {
                if (compare_str_cstr(str[j], atom_str)) {
                    ++count;
                } else {
                    return -1;
                }
            }
        }
        return count;
    }

    return 0;
}

static int _element_str(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));
    ASSERT(ctx && ctx->mol && ctx->mol->atom.element);

    uint8_t* elem_idx = 0;

    const uint64_t num_str = (uint64_t)element_count(arg[0]);
    const str_t* str = as_string_arr(arg[0]);

    for (uint64_t i = 0; i < num_str; ++i) {
        uint8_t idx = element_symbol_to_index(str[i]);
        if (idx)
            md_array_push(elem_idx, idx, ctx->temp_alloc);
        else
            return -1;
    }

    if (dst) {
        ASSERT(is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_ATOM));
        bitfield_t result = as_bitfield(*dst);

        const uint64_t num_elem = md_array_size(elem_idx);
        for (uint64_t i = ctx->mol_ctx.atom.beg; i < ctx->mol_ctx.atom.end; ++i) {
            for (uint64_t j = 0; j < num_elem; ++j) {
                if (elem_idx[j] == ctx->mol->atom.element[i]) {
                    bit_set_idx(result.bits, i);
                    break;
                }
            }
        }
    
    } else {
        int count = 0;
        const uint64_t num_elem = md_array_size(elem_idx);
        for (uint64_t i = ctx->mol_ctx.atom.beg; i < ctx->mol_ctx.atom.end; ++i) {
            bool found = false;
            for (uint64_t j = 0; j < num_elem; ++j) {
                if (elem_idx[j] == ctx->mol->atom.element[i]) {
                    ++count;
                    found = true;
                    break;
                }
            }
            if (!found) {
                return -1;
            }
        }
        return count;
    }

    return 0;
}

static int _element_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_ATOM));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    ASSERT(ctx && ctx->mol && ctx->mol->atom.element);

    bitfield_t result = as_bitfield(*dst);
    irange_t* ranges = as_irange_arr(arg[0]);
    const uint64_t num_ranges = element_count(arg[0]);

    for (uint64_t i = ctx->mol_ctx.atom.beg; i < ctx->mol_ctx.atom.end; ++i) {
        for (uint64_t j = 0; j < num_ranges; ++j) {
            if (in_range(ranges[j], ctx->mol->atom.element[i])) {
                bit_set_idx(result.bits, i);
                break;
            }
        }
    }

    return 0;
}

static int _atom(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_ATOM));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    ASSERT(ctx && ctx->mol);

    bitfield_t result = as_bitfield(*dst);
    irange_t* ranges = as_irange_arr(arg[0]);
    const uint64_t num_ranges = element_count(arg[0]);

    for (uint64_t j = 0; j < num_ranges; ++j) {
        irange_t range = ranges[j];
        range.beg = CLAMP(range.beg + ctx->mol_ctx.atom.beg, (int32_t)ctx->mol_ctx.atom.beg, (int32_t)ctx->mol_ctx.atom.end);
        range.end = CLAMP(range.end + ctx->mol_ctx.atom.beg, (int32_t)ctx->mol_ctx.atom.beg, (int32_t)ctx->mol_ctx.atom.end);
        bit_set(result.bits, range.beg, range.end-range.beg);
    }

    return 0;
}

// @PERFORMANCE: Improve these by using spatial hashing
static int _x(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x);
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_ATOM));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FRANGE));
    bitfield_t result = as_bitfield(*dst);
    const frange_t range = as_frange(arg[0]);

    for (uint64_t i = ctx->mol_ctx.atom.beg; i < ctx->mol_ctx.atom.end; ++i) {
        if (range.beg <= ctx->mol->atom.x[i] && ctx->mol->atom.x[i] <= range.end) {
            bit_set_idx(result.bits, i);
        }
    }

    return 0;
}

static int _y(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.y);
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_ATOM));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FRANGE));
    bitfield_t result = as_bitfield(*dst);
    const frange_t range = as_frange(arg[0]);

    for (uint64_t i = ctx->mol_ctx.atom.beg; i < ctx->mol_ctx.atom.end; ++i) {
        if (range.beg <= ctx->mol->atom.y[i] && ctx->mol->atom.y[i] <= range.end) {
            bit_set_idx(result.bits, i);
        }
    }
    return 0;
}

static int _z(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.z);
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_ATOM));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FRANGE));
    bitfield_t result = as_bitfield(*dst);
    const frange_t range = as_frange(arg[0]);

    for (uint64_t i = ctx->mol_ctx.atom.beg; i < ctx->mol_ctx.atom.end; ++i) {
        if (range.beg <= ctx->mol->atom.z[i] && ctx->mol->atom.z[i] <= range.end) {
            bit_set_idx(result.bits, i);
        }
    }
    return 0;
}

static int _within_flt(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.z);
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_ATOM));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_BITFIELD));

    ASSERT(false);

    return 0;
}

static int _within_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.z);
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_ATOM));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FRANGE));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_BITFIELD));

    ASSERT(false);

    return 0;
}

static int _resname(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->residue.name);
    ASSERT(dst && dst->ptr && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_RESIDUE));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));
    bitfield_t* bf = dst->ptr;

    const uint64_t num_str = element_count(arg[0]);
    const str_t*   str     = arg[0].ptr;

    for (uint64_t i = ctx->mol_ctx.residue.beg; i < ctx->mol_ctx.residue.end; ++i) {
        for (uint64_t j = 0; j < num_str; ++j) {
            if (compare_str_cstr(str[j], ctx->mol->residue.name[i])) {
                uint64_t offset = ctx->mol->residue.atom_range[i].beg;
                uint64_t length = ctx->mol->residue.atom_range[i].end - ctx->mol->residue.atom_range[i].beg;
                bit_set(bf->bits, offset, length);
                break;
            }
        }
    }

    return 0;
}

static int _resid(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->residue.name);
    ASSERT(dst && dst->ptr && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_RESIDUE));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    bitfield_t* bf = dst->ptr;

    const uint64_t  num_rid = element_count(arg[0]);
    const irange_t*     rid = arg[0].ptr;

    for (uint64_t i = ctx->mol_ctx.residue.beg; i < ctx->mol_ctx.residue.end; ++i) {
        for (uint64_t j = 0; j < num_rid; ++j) {
            if (in_range(rid[j], (int)ctx->mol->residue.id[i])) {
                uint64_t offset = ctx->mol->residue.atom_range[i].beg;
                uint64_t length = ctx->mol->residue.atom_range[i].end - ctx->mol->residue.atom_range[i].beg;
                bit_set(bf->bits, offset, length);
                break;
            }
        }
    }

    return 0;
}

static int _residue(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->residue.name);
    ASSERT(dst && dst->ptr && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_RESIDUE));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    bitfield_t* bf = dst->ptr;

    const uint64_t  num_ranges = element_count(arg[0]);
    const irange_t*     ranges = arg[0].ptr;

    for (uint64_t i = 0; i < num_ranges; ++i) {
        irange_t range = ranges[i];
        range.beg = CLAMP(range.beg, (int32_t)ctx->mol_ctx.residue.beg, (int32_t)ctx->mol_ctx.residue.end - 1);
        range.end = CLAMP(range.end, (int32_t)ctx->mol_ctx.residue.beg, (int32_t)ctx->mol_ctx.residue.end - 1);
        for (uint64_t j = range.beg; j <= range.end; ++j) {
            const uint64_t offset = ctx->mol->residue.atom_range[j].beg;
            const uint64_t length = ctx->mol->residue.atom_range[j].end - ctx->mol->residue.atom_range[j].beg;
            bit_set(bf->bits, offset, length);
        }
    }

    return 0;
}

static int _chain_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->chain.atom_range);
    ASSERT(dst && dst->ptr && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_CHAIN));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    bitfield_t* bf = dst->ptr;

    const uint64_t  num_ranges = element_count(arg[0]);
    const irange_t*     ranges = arg[0].ptr;

    for (uint64_t i = 0; i < num_ranges; ++i) {
        irange_t range = ranges[i];
        range.beg = CLAMP(range.beg, (int32_t)ctx->mol_ctx.chain.beg, (int32_t)ctx->mol_ctx.chain.end - 1);
        range.end = CLAMP(range.end, (int32_t)ctx->mol_ctx.chain.beg, (int32_t)ctx->mol_ctx.chain.end - 1);
        for (uint64_t j = range.beg; j <= range.end; ++j) {
            const uint64_t offset = ctx->mol->chain.atom_range[j].beg;
            const uint64_t length = ctx->mol->chain.atom_range[j].end - ctx->mol->chain.atom_range[j].beg;
            bit_set(bf->bits, offset, length);
        }
    }

    return 0;
}

static int _chain_str(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->chain.id && ctx->mol->chain.atom_range);
    ASSERT(dst && dst->ptr && is_type_equivalent(dst->type, (type_info_t)TI_BITFIELD_RESIDUE));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));
    bitfield_t* bf = dst->ptr;

    const uint64_t num_str = element_count(arg[0]);
    const str_t*   str     = arg[0].ptr;

    for (uint64_t i = ctx->mol_ctx.residue.beg; i < ctx->mol_ctx.residue.end; ++i) {
        for (uint64_t j = 0; j < num_str; ++j) {
            if (compare_str_cstr(str[j], ctx->mol->chain.id[i])) {
                uint64_t offset = ctx->mol->chain.atom_range[i].beg;
                uint64_t length = ctx->mol->chain.atom_range[i].end - ctx->mol->chain.atom_range[i].beg;
                bit_set(bf->bits, offset, length);
                break;
            }
        }
    }

    return 0;
}

// Property Compute

static int _distance(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x && ctx->mol->atom.y && ctx->mol->atom.z);
    ASSERT(dst);
    ASSERT(compare_type_info(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_FLOAT3));
    ASSERT(compare_type_info(arg[1].type, (type_info_t)TI_FLOAT3));

    const vec3 a = as_vec3(arg[0]);
    const vec3 b = as_vec3(arg[1]);
    const float d = vec3_dist(a, b);
    as_float(*dst) = d;
    dst->unit = MD_SCRIPT_UNIT_ANGSTROM;

    return 0;
}

static int _distance_min(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x && ctx->mol->atom.y && ctx->mol->atom.z);
    ASSERT(dst);
    ASSERT(compare_type_info(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_FLOAT3_ARR));
    ASSERT(compare_type_info(arg[1].type, (type_info_t)TI_FLOAT3_ARR));

    const vec3* a = as_vec3_arr(arg[0]);
    const vec3* b = as_vec3_arr(arg[1]);
    int64_t a_len = element_count(arg[0]);
    int64_t b_len = element_count(arg[1]);

    float min_dist = 0;
    for (int64_t i = 0; i < a_len; ++i) {
        for (int64_t j = 0; j < b_len; ++j) {
            const float dist = vec3_dist(a[i], b[j]);
            min_dist = MIN(min_dist, dist);
        }
    }

    as_float(*dst) = min_dist;
    dst->unit = MD_SCRIPT_UNIT_ANGSTROM;

    return 0;
}

static int _distance_max(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x && ctx->mol->atom.y && ctx->mol->atom.z);
    ASSERT(dst);
    ASSERT(compare_type_info(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_FLOAT3_ARR));
    ASSERT(compare_type_info(arg[1].type, (type_info_t)TI_FLOAT3_ARR));

    const vec3* a = as_vec3_arr(arg[0]);
    const vec3* b = as_vec3_arr(arg[1]);
    int64_t a_len = element_count(arg[0]);
    int64_t b_len = element_count(arg[1]);

    float max_dist = 0;
    for (int64_t i = 0; i < a_len; ++i) {
        for (int64_t j = 0; j < b_len; ++j) {
            const float dist = vec3_dist(a[i], b[j]);
            max_dist = MAX(max_dist, dist);
        }
    }

    as_float(*dst) = max_dist;
    dst->unit = MD_SCRIPT_UNIT_ANGSTROM;
    return 0;
}

static int _distance_pair(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x && ctx->mol->atom.y && ctx->mol->atom.z);
    ASSERT(dst);
    ASSERT(compare_type_info(dst->type, (type_info_t)TI_FLOAT_ARR));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_FLOAT3_ARR));
    ASSERT(compare_type_info(arg[1].type, (type_info_t)TI_FLOAT3_ARR));

    const int64_t a_len = element_count(arg[0]);
    const int64_t b_len = element_count(arg[1]);
    if (dst) {
        const vec3* a = as_vec3_arr(arg[0]);
        const vec3* b = as_vec3_arr(arg[1]);
        float* dst_arr = as_float_arr(*dst);
        const int64_t dst_len = element_count(*dst);
        ASSERT(dst_len == a_len * b_len);

        for (int64_t i = 0; i < a_len; ++i) {
            for (int64_t j = 0; j < b_len; ++j) {
                const float dist = vec3_dist(a[i], b[j]);
                dst_arr[i * b_len + j] = dist;
            }
        }

        dst->unit = MD_SCRIPT_UNIT_ANGSTROM;
        return 0;
    } else {
        return (int)(a_len * b_len);
    }
}

static int _angle(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x && ctx->mol->atom.y && ctx->mol->atom.z);
    ASSERT(dst);
    ASSERT(compare_type_info(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_FLOAT3));
    ASSERT(compare_type_info(arg[1].type, (type_info_t)TI_FLOAT3));
    ASSERT(compare_type_info(arg[2].type, (type_info_t)TI_FLOAT3));
    
    const vec3 a = as_vec3(arg[0]);
    const vec3 b = as_vec3(arg[1]);
    const vec3 c = as_vec3(arg[2]);
    const vec3 v0 = vec3_normalize(vec3_sub(a, b));
    const vec3 v1 = vec3_normalize(vec3_sub(c, b));
    
    as_float(*dst) = (float)RAD_TO_DEG(acosf(vec3_dot(v0, v1)));
    dst->unit = MD_SCRIPT_UNIT_DEGREES;
    return 0;
}

static int _dihedral(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x && ctx->mol->atom.y && ctx->mol->atom.z);
    ASSERT(dst);
    ASSERT(compare_type_info(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_FLOAT3));
    ASSERT(compare_type_info(arg[1].type, (type_info_t)TI_FLOAT3));
    ASSERT(compare_type_info(arg[2].type, (type_info_t)TI_FLOAT3));
    ASSERT(compare_type_info(arg[3].type, (type_info_t)TI_FLOAT3));

    ASSERT(false);

    return 0;
}

static int _rmsd(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->mol && ctx->mol->atom.x && ctx->mol->atom.y && ctx->mol->atom.z);
    ASSERT(dst);
    ASSERT(compare_type_info(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));

    ASSERT(false);

    return 0;
}


// #################
// ###   CASTS   ###
// #################

static int _cast_int_to_flt(data_t* dst, data_t arg[], eval_context_t* ctx){
    ASSERT(dst && compare_type_info(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_INT));
    (void)ctx;
    as_int(*dst) = (int)as_float(arg[0]);
    return 0;
}

static int _cast_int_to_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && compare_type_info(dst->type, (type_info_t)TI_IRANGE));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_INT));
    (void)ctx;
    as_irange(*dst) = (irange_t){as_int(arg[0]), as_int(arg[0])};
    return 0;
}

static int _cast_irng_to_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && compare_type_info(dst->type, (type_info_t)TI_FRANGE));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_IRANGE));
    (void)ctx;
    as_frange(*dst) = (frange_t){(float)as_irange(arg[0]).beg, (float)as_irange(arg[0]).end};
    return 0;
}

static int _cast_int_arr_to_flt_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_INT_ARR));
    (void)ctx;
    const uint64_t dst_len = element_count(*dst);
    const uint64_t arg_len = element_count(arg[0]);
    ASSERT(dst_len == arg_len);

    float* d = dst->ptr;
    int*   s = arg[0].ptr;

    for (uint64_t i = 0; i < dst_len; ++i) {
        d[i] = (float)s[i]; 
    }
    return 0;
}
static int _cast_irng_arr_to_frng_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_FRANGE_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    (void)ctx;
    const uint64_t dst_len = element_count(*dst);
    const uint64_t arg_len = element_count(arg[0]);
    ASSERT(dst_len == arg_len);

    frange_t* d = dst->ptr;
    irange_t* s = arg[0].ptr;

    for (uint64_t i = 0; i < dst_len; ++i) {
        d[i].beg = (float)s[i].beg;
        d[i].end = (float)s[i].end; 
    }
    return 0;
}



static int _cast_bf_to_atom(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && compare_type_info(dst->type, (type_info_t)TI_BITFIELD_ATOM));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_BITFIELD));
    (void)dst;
    (void)arg;
    (void)ctx;
    return 0;
}

static int _cast_bf_to_residue(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && compare_type_info(dst->type, (type_info_t)TI_BITFIELD_RESIDUE));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_BITFIELD));
    (void)dst;
    (void)arg;
    (void)ctx;
    return 0;
}

static int _cast_bf_to_chain(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && compare_type_info(dst->type, (type_info_t)TI_BITFIELD_CHAIN));
    ASSERT(compare_type_info(arg[0].type, (type_info_t)TI_BITFIELD));
    (void)dst;
    (void)arg;
    (void)ctx;
    return 0;
}

static int _com(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && compare_type_info(dst->type, (type_info_t)TI_FLOAT3));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT3_ARR));
    (void)dst;
    (void)arg;
    (void)ctx;

    const vec3* pos = as_vec3_arr(arg[0]);
    const int64_t pos_size = element_count(arg[0]);

    vec3 com = {0};
    if (pos_size > 0) {
        for (int64_t i = 0; i < pos_size; ++i) {
            com = vec3_add(com, pos[i]);
        }
        const double div = 1.0 / (double)pos_size;
        com.x = (float)(com.x / div);
        com.y = (float)(com.y / div);
        com.z = (float)(com.z / div);
    }

    as_vec3(*dst) = com;

    return 0;
}

static int _plane(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(dst && compare_type_info(dst->type, (type_info_t)TI_FLOAT4));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT3_ARR));
    (void)dst;
    (void)arg;
    (void)ctx;

    ASSERT(false);
    return 0;
}

static int _position_int(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_INT_ARR));
    ASSERT(arg[0].ptr);
    const int* indices = as_int_arr(arg[0]);
    const uint64_t num_idx = element_count(arg[0]);

    if (dst) {
        ASSERT(element_count(*dst) == element_count(arg[0]));
        ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT3_ARR));
        vec3* position = as_vec3_arr(*dst);

        for (uint64_t i = 0; i < num_idx; ++i) {
            // Shift here since we use 1 based indices for atoms
            const uint64_t idx = ctx->mol_ctx.atom.beg + indices[i] - 1;
            ASSERT(0 <= idx && idx < ctx->mol_ctx.atom.end);
            position[i] = (vec3){ctx->mol->atom.x[idx], ctx->mol->atom.y[idx], ctx->mol->atom.z[idx]};
        }
    } else {
        for (uint64_t i = 0; i < num_idx; ++i) {
            if (indices[i] < 1 || (int)ctx->mol_ctx.atom.end < indices[i]) {
                return -1;
            }
        }
        return (int)num_idx;
    }

    return 0;
}

static int _position_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    ASSERT(arg[0].ptr);
    const irange_t* ranges = as_irange_arr(arg[0]);
    const uint64_t num_ranges = element_count(arg[0]);

    if (dst) {
        ASSERT(element_count(*dst) == element_count(arg[0]));
        ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT3_ARR));
        vec3* position = as_vec3_arr(*dst);

        uint64_t dst_idx = 0;
        for (uint64_t i = 0; i < num_ranges; ++i) {
            for (uint64_t j = ranges[i].beg; j < ranges[j].end; ++j) {
                const uint64_t src_idx = ctx->mol_ctx.atom.beg + j;
                position[dst_idx++] = (vec3){ctx->mol->atom.x[src_idx], ctx->mol->atom.y[src_idx], ctx->mol->atom.z[src_idx]};
            }
        }
    } else {
        int count = 0;
        for (uint64_t i = 0; i < num_ranges; ++i) {
            if ((int32_t)ctx->mol_ctx.atom.beg <= ranges[i].beg && ranges[i].end < (int32_t)ctx->mol_ctx.atom.end) {
                // We good!
                count += (ranges[i].end - ranges[i].beg);
            }
            else {
                return -1;
            }
        }
        return count;
    }

    return 0;
}

static int _position_bf(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD));
    ASSERT(arg[0].ptr);
    const bitfield_t bf = as_bitfield(arg[0]);
    const uint64_t offset = ctx->mol_ctx.atom.beg;
    const uint64_t length = ctx->mol_ctx.atom.end - ctx->mol_ctx.atom.beg;

    if (dst) {
        ASSERT(element_count(*dst) == (int64_t)bit_count(bf.bits, offset, length));
        ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT3_ARR));
        vec3* position = as_vec3_arr(*dst);

        int64_t dst_idx = 0;
        uint64_t bit_idx = 0;
        while ((bit_idx = bit_scan(bf.bits, offset + bit_idx, length - bit_idx)) != 0) {
            uint64_t src_idx = bit_idx - 1;
            position[dst_idx++] = (vec3){ctx->mol->atom.x[src_idx], ctx->mol->atom.y[src_idx], ctx->mol->atom.z[src_idx]};
        }
        ASSERT(dst_idx == element_count(*dst));
    } else {
        return (int)bit_count(bf.bits, offset, length);
    }

    return 0;
}

static int _rdf(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)ctx;
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT3_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT3_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FLOAT));
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_DISTRIBUTION));
    ASSERT(arg[0].ptr);
    ASSERT(arg[1].ptr);
    ASSERT(arg[2].ptr);
    ASSERT(dst->ptr);

    const vec3* ref_pos = as_vec3_arr(arg[0]);
    const int64_t ref_size = element_count(arg[0]);

    const vec3* target_pos = as_vec3_arr(arg[1]);
    const int64_t target_size = element_count(arg[1]);

    const float cutoff = as_float(arg[2]);
    const int32_t num_bins = DIST_BINS;
    float* bins = as_float_arr(*dst);

    memset(bins, 0, num_bins * sizeof(float));

    for (int64_t i = 0; i < ref_size; ++i) {
        for (int64_t j = 0; j < target_size; ++j) {
            const float d = vec3_dist(ref_pos[i], target_pos[j]);
            if (d < cutoff) {
                const int32_t bin_idx = (int32_t)((d / cutoff) * DIST_BINS);
                bins[bin_idx] += 1.0f;
            }
        }
    }

    // Normalize the distribution
    const float scl = 1.0f / num_bins;
    for (int64_t i = 0; i < num_bins; ++i) {
        bins[i] *= scl;
    }

    return 0;
}

static int _sdf(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_BITFIELD_RESIDUE));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT3_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FLOAT));
    ASSERT(dst && is_type_equivalent(dst->type, (type_info_t)TI_VOLUME));
    ASSERT(arg[0].ptr);
    ASSERT(arg[1].ptr);
    ASSERT(arg[2].ptr);
    ASSERT(dst->ptr);

    if (dst) {
        
    } else {
        
    }

    return 0;
}

#undef ANY_LENGTH
#undef ANY_LEVEL

#endif