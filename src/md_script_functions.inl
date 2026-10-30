// This is to mark that the procedure supports a varying length
#define ANY_LENGTH -1

#ifndef MD_DIST_BINS
#define MD_DIST_BINS 1024
#endif

#ifndef MD_VOL_DIM
#define MD_VOL_DIM 128
#endif

#define STATIC_VALIDATION_ERROR -2


// Type info declaration helpers
#define TI_BOOL         {TYPE_BOOL,  {1}}
#define TI_BOOL_ARR     {TYPE_BOOL,  {ANY_LENGTH}}
#define TI_FLOAT        {TYPE_FLOAT, {1}}
#define TI_FLOAT_ARR    {TYPE_FLOAT, {ANY_LENGTH}}
#define TI_FLOAT2       {TYPE_FLOAT, {2}}
#define TI_FLOAT2_ARR   {TYPE_FLOAT, {ANY_LENGTH,2}}
#define TI_FLOAT3       {TYPE_FLOAT, {3}}
#define TI_FLOAT3_ARR   {TYPE_FLOAT, {ANY_LENGTH,3}}
#define TI_FLOAT4       {TYPE_FLOAT, {4}}
#define TI_FLOAT4_ARR   {TYPE_FLOAT, {ANY_LENGTH,4}}
#define TI_FLOAT44      {TYPE_FLOAT, {4,4}}
#define TI_FLOAT44_ARR  {TYPE_FLOAT, {ANY_LENGTH,4,4}}

// One distribution, and any number of them, over MD_DIST_BINS bins. The weights are not in the type (data_t.weight).
#define TI_DISTRIBUTION     {TYPE_DISTRIBUTION, {1,MD_DIST_BINS}}
#define TI_DISTRIBUTION_ARR {TYPE_DISTRIBUTION, {ANY_LENGTH,MD_DIST_BINS}}
#define TI_VOLUME       {TYPE_FLOAT, {MD_VOL_DIM, MD_VOL_DIM, MD_VOL_DIM}}

#define TI_INT          {TYPE_INT, {1}}
#define TI_INT_ARR      {TYPE_INT, {ANY_LENGTH}}

#define TI_FRANGE       {TYPE_FRANGE, {1}}
#define TI_FRANGE_ARR   {TYPE_FRANGE, {ANY_LENGTH}}

#define TI_IRANGE       {TYPE_IRANGE, {1}}
#define TI_IRANGE_ARR   {TYPE_IRANGE, {ANY_LENGTH}}

#define TI_STRING       {TYPE_STRING, {1}}
#define TI_STRING_ARR   {TYPE_STRING, {ANY_LENGTH}}

#define TI_ATOMSET     {TYPE_ATOMSET, {1}}
#define TI_ATOMSET_ARR {TYPE_ATOMSET, {ANY_LENGTH}}

#define TI_CONTACT        {TYPE_CONTACT, {1}}
#define TI_PLANE          {TYPE_PLANE, {1}}
#define TI_PLANE_ARR      {TYPE_PLANE, {ANY_LENGTH}}
#define TI_PUCKER         {TYPE_PUCKER, {1}}
#define TI_PUCKER_ARR     {TYPE_PUCKER, {ANY_LENGTH}}

#define TI_BONDSET        {TYPE_BONDSET, {1}}
#define TI_BONDSET_ARR    {TYPE_BONDSET, {ANY_LENGTH}}
#define TI_BOND_DELTA     {TYPE_BOND_DELTA, {1}}

#define TI_COORDINATE     {TYPE_COORDINATE, {1}}
#define TI_COORDINATE_ARR {TYPE_COORDINATE, {ANY_LENGTH}}

#define as_bool(arg) (*((bool*)((arg).ptr)))
#define as_bool_arr(arg) ((bool*)((arg).ptr))

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

#define as_bitset(arg) ((md_bitset_t*)((arg).ptr))
#define as_bondset(arg) ((bondset_t*)((arg).ptr))

// Sets written as results are built in a mask in the temp arena (over the atoms of the system: it allocates only
// the words it writes), then extracted into the allocator of the destination: their words live exactly as long as
// the data holding them (see data_bitset_alloc). Function-local scratch scopes are therefore opened with
// md_temp_begin_avoid(ctx->temp_alloc): rewinding one must never release a result, or a destination, which lives
// in ctx->temp_alloc.
static inline md_bitmask_t result_mask(eval_context_t* ctx) {
    return md_bitmask_create(0, atom_domain_end(ctx), ctx->temp_alloc);
}

static inline void finish_result(data_t* dst, size_t idx, md_bitmask_t* b, eval_context_t* ctx) {
    as_bitset(*dst)[idx] = md_bitmask_extract(b, data_bitset_alloc(dst, ctx));
    md_bitmask_free(b);
}

// Restricts a mask to the evaluation context, if there is one
static inline void restrict_to_context(md_bitmask_t* b, const eval_context_t* ctx) {
    if (ctx->mol_ctx) md_bitmask_and(b, *ctx->mol_ctx);
}

// Writes one set per element of dst, in order. A destination which holds a single set is flattened: every
// set pushed is added to it. Separate ranges are RUNs and cost nothing; a union is built, so it stays linear.
typedef struct set_writer_t {
    data_t* dst;
    eval_context_t* ctx;
    md_bitmask_t union_b;
    size_t cap;
    size_t idx;
    bool restrict_ctx;  // intersect each set with the evaluation context
} set_writer_t;

static inline set_writer_t set_writer_begin(data_t* dst, eval_context_t* ctx, bool restrict_ctx) {
    set_writer_t w = {0};
    w.dst = dst;
    w.ctx = ctx;
    w.cap = (size_t)MAX(0, type_info_array_len(dst->type));
    w.restrict_ctx = restrict_ctx && ctx->mol_ctx;
    if (w.cap == 1) w.union_b = result_mask(ctx);
    return w;
}

static inline void set_writer_push(set_writer_t* w, md_bitset_t s) {
    if (w->cap == 1) {
        md_bitmask_or(&w->union_b, s);
        return;
    }
    ASSERT(w->idx < w->cap);
    if (w->idx >= w->cap) return;
    md_allocator_i* alloc = data_bitset_alloc(w->dst, w->ctx);
    as_bitset(*w->dst)[w->idx++] = w->restrict_ctx ? md_bitset_and(s, *w->ctx->mol_ctx, alloc) : md_bitset_copy(s, alloc);
}

static inline void set_writer_push_range(set_writer_t* w, uint32_t beg, uint32_t end) {
    set_writer_push(w, md_bitset_range(beg, end));
}

static inline void set_writer_end(set_writer_t* w) {
    if (w->cap == 1) {
        if (w->restrict_ctx) restrict_to_context(&w->union_b, w->ctx);
        finish_result(w->dst, 0, &w->union_b, w->ctx);
    }
}

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

// The element wise operations below count the elements by the bytes there are, of the destination and of the
// arguments alike: the dimensions of an empty array (none, all 0) make for a size of 1 (see dim_size). For the same
// reason the type of an empty destination is not that of an array, and the asserts on it let one through
// (dst->size == 0).
#define ELEMENTWISE_COUNT(d, a, type) (MIN((d).size, (a).size) / sizeof(type))

#define BAKE_FUNC_FARR__FARR(prefix, func) \
    static int prefix##func(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        if (dst) { \
            const size_t n = ELEMENTWISE_COUNT(*dst, arg[0], float); \
            for (size_t i = 0; i < n; ++i) { \
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
        const size_t n = ELEMENTWISE_COUNT(*dst, arg[0], base_type); \
        for (size_t i = 0; i < n; ++i) { \
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
        const size_t n = ELEMENTWISE_COUNT(*dst, arg[0], base_type); \
        for (size_t i = 0; i < n; ++i) { \
            ((base_type*)dst->ptr)[i] = ((base_type*)arg[0].ptr)[i] op *((base_type*)arg[1].ptr); \
        } \
        return 0; \
    }

// The scalar on the left: 's - a', 's / a'. The operators which are not commutative have one of these for an array on
// the right, rather than FLAG_SYMMETRIC_ARGS, which swaps the operands into 'a - s'.
#define BAKE_OP_S_M(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        if (!dst) return 0; \
        const base_type s = *((base_type*)arg[0].ptr); \
        const size_t n = ELEMENTWISE_COUNT(*dst, arg[1], base_type); \
        for (size_t i = 0; i < n; ++i) { \
            ((base_type*)dst->ptr)[i] = s op ((base_type*)arg[1].ptr)[i]; \
        } \
        return 0; \
    }

#define BAKE_OP_M_M(name, op, base_type) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        if (!dst) return 0; \
        const size_t n = MIN(ELEMENTWISE_COUNT(*dst, arg[0], base_type), arg[1].size / sizeof(base_type)); \
        for (size_t i = 0; i < n; ++i) { \
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
        const size_t count = (size_t)type_info_dim_size(dst->type); \
        ASSERT(count == (size_t)type_info_dim_size(arg[0].type)); \
        ASSERT(count == (size_t)type_info_dim_size(arg[1].type)); \
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
        ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT)); \
        const size_t count = (size_t)type_info_dim_size(dst->type); \
        ASSERT(count == (size_t)type_info_dim_size(arg[0].type)); \
        ASSERT(count % 8 == 0); \
        const float* src_arr = as_float_arr(arg[0]); \
        float* dst_arr = as_float_arr(*dst); \
        md_256 s = md_mm256_set1_ps(as_float(arg[1])); \
        for (size_t i = 0; i < count; i += 8) { \
            md_mm256_storeu_ps(dst_arr + i, op(md_mm256_loadu_ps(src_arr + i), s)); \
        } \
        return 0; \
    }

#define BAKE_SIMDF_OP_S_M(name, op) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        ASSERT(dst); \
        ASSERT(dst->type.base_type == TYPE_FLOAT); \
        ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT)); \
        ASSERT(arg[1].type.base_type == TYPE_FLOAT); \
        const size_t count = (size_t)type_info_dim_size(dst->type); \
        ASSERT(count == (size_t)type_info_dim_size(arg[1].type)); \
        ASSERT(count % 8 == 0); \
        const float* src_arr = as_float_arr(arg[1]); \
        float* dst_arr = as_float_arr(*dst); \
        md_256 s = md_mm256_set1_ps(as_float(arg[0])); \
        for (size_t i = 0; i < count; i += 8) { \
            md_mm256_storeu_ps(dst_arr + i, op(s, md_mm256_loadu_ps(src_arr + i))); \
        } \
        return 0; \
    }

#define BAKE_SIMDF_OP_M(name, op) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        ASSERT(dst); \
        ASSERT(dst->type.base_type == TYPE_FLOAT); \
        ASSERT(arg[0].type.base_type == TYPE_FLOAT); \
        const size_t count = (size_t)type_info_dim_size(dst->type); \
        ASSERT(count == (size_t)type_info_dim_size(arg[0].type)); \
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

// Element-wise logic on arrays of bool (for example the result of a comparison)
BAKE_OP_UNARY_M(_op_not_barr, !, bool)
BAKE_OP_M_S(_op_and_barr_b, &&, bool)
BAKE_OP_M_S(_op_or_barr_b,  ||, bool)
BAKE_OP_M_S(_op_xor_barr_b,  ^, bool)
BAKE_OP_M_M(_op_and_barr_barr, &&, bool)
BAKE_OP_M_M(_op_or_barr_barr,  ||, bool)
BAKE_OP_M_M(_op_xor_barr_barr,  ^, bool)

BAKE_OP_S_S(_op_add_i_i, +, int)
BAKE_OP_S_S(_op_sub_i_i, -, int)
BAKE_OP_S_S(_op_mul_i_i, *, int)
// Integer division traps (SIGFPE) on a zero divisor, and on INT32_MIN // -1, whose result does not fit: an error
static bool int_div_valid(int a, int b, eval_context_t* ctx) {
    if (b == 0) {
        EVAL_LOG_ERROR(ctx, ctx->op_token, "Integer division by zero");
        return false;
    }
    if (a == INT32_MIN && b == -1) {
        EVAL_LOG_ERROR(ctx, ctx->op_token, "Integer division overflows: %i // %i", a, b);
        return false;
    }
    return true;
}

static int _op_div_i_i(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    const int a = as_int(arg[0]);
    const int b = as_int(arg[1]);
    if (!int_div_valid(a, b, ctx)) return -1;
    as_int(*dst) = a / b;
    return 0;
}
BAKE_OP_UNARY_S(_op_neg_i, -, int)
BAKE_OP_UNARY_M(_op_neg_iarr, -, int)

BAKE_OP_M_S(_op_add_iarr_i, +, int)
BAKE_OP_M_S(_op_sub_iarr_i, -, int)
BAKE_OP_S_M(_op_sub_i_iarr, -, int)
BAKE_OP_M_S(_op_mul_iarr_i, *, int)
static int _op_div_iarr_i(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    const int* a = (const int*)arg[0].ptr;
    const int b = as_int(arg[1]);
    int* d = (int*)dst->ptr;
    const size_t n = ELEMENTWISE_COUNT(*dst, arg[0], int);
    for (size_t i = 0; i < n; ++i) {
        if (!int_div_valid(a[i], b, ctx)) return -1;
        d[i] = a[i] / b;
    }
    return 0;
}

static int _op_div_i_iarr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    const int a = as_int(arg[0]);
    const int* b = (const int*)arg[1].ptr;
    int* d = (int*)dst->ptr;
    const size_t n = ELEMENTWISE_COUNT(*dst, arg[1], int);
    for (size_t i = 0; i < n; ++i) {
        if (!int_div_valid(a, b[i], ctx)) return -1;
        d[i] = a / b[i];
    }
    return 0;
}

BAKE_OP_M_M(_op_add_iarr_iarr, +, int)
BAKE_OP_M_M(_op_sub_iarr_iarr, -, int)
BAKE_OP_M_M(_op_mul_iarr_iarr, *, int)
static int _op_div_iarr_iarr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    const int* a = (const int*)arg[0].ptr;
    const int* b = (const int*)arg[1].ptr;
    int* d = (int*)dst->ptr;
    const size_t n = MIN(ELEMENTWISE_COUNT(*dst, arg[0], int), arg[1].size / sizeof(int));
    for (size_t i = 0; i < n; ++i) {
        if (!int_div_valid(a[i], b[i], ctx)) return -1;
        d[i] = a[i] / b[i];
    }
    return 0;
}

BAKE_SIMDF_OP_M_M(_op_simd_add_farr_farr, md_mm256_add_ps)
BAKE_SIMDF_OP_M_S(_op_simd_add_farr_f,    md_mm256_add_ps)

BAKE_SIMDF_OP_M_M(_op_simd_sub_farr_farr, md_mm256_sub_ps)
BAKE_SIMDF_OP_M_S(_op_simd_sub_farr_f,    md_mm256_sub_ps)
BAKE_SIMDF_OP_S_M(_op_simd_sub_f_farr,    md_mm256_sub_ps)

BAKE_SIMDF_OP_M_M(_op_simd_mul1arr_farr, md_mm256_mul_ps)
BAKE_SIMDF_OP_M_S(_op_simd_mul1arr_f,    md_mm256_mul_ps)

BAKE_SIMDF_OP_M_M(_op_simd_div_farr_farr, md_mm256_div_ps)
BAKE_SIMDF_OP_M_S(_op_simd_div_farr_f,    md_mm256_div_ps)
BAKE_SIMDF_OP_S_M(_op_simd_div_f_farr,    md_mm256_div_ps)

BAKE_SIMDF_OP_M(_op_simd_abs_farr,        md_mm256_abs_ps)

BAKE_SIMDF_OP_M_M(_op_simd_min_farr_farr, md_mm256_min_ps)
BAKE_SIMDF_OP_M_S(_op_simd_min_farr_f,    md_mm256_min_ps)
BAKE_SIMDF_OP_M_M(_op_simd_max_farr_farr, md_mm256_max_ps)
BAKE_SIMDF_OP_M_S(_op_simd_max_farr_f,    md_mm256_max_ps)

BAKE_OP_S_S(_op_add_f_f,        +, float)
BAKE_OP_S_S(_op_sub_f_f,        -, float)
BAKE_OP_S_S(_op_mul1_f,        *, float)
BAKE_OP_S_S(_op_div_f_f,        /, float)
BAKE_OP_UNARY_S(_op_neg_f,      -, float)
BAKE_OP_UNARY_M(_op_neg_farr,   -, float)

BAKE_OP_M_S(_op_add_farr_f, +, float)
BAKE_OP_M_S(_op_sub_farr_f, -, float)
BAKE_OP_S_M(_op_sub_f_farr, -, float)
BAKE_OP_M_S(_op_mul1arr_f, *, float)
BAKE_OP_M_S(_op_div_farr_f, /, float)
BAKE_OP_S_M(_op_div_f_farr, /, float)

BAKE_OP_M_M(_op_add_farr_farr, +, float)
BAKE_OP_M_M(_op_sub_farr_farr, -, float)
BAKE_OP_M_M(_op_mul1arr_farr, *, float)
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
static int _join_bf_arr                 (data_t*, data_t[], eval_context_t*);
static int _join_deprecated            (data_t*, data_t[], eval_context_t*);
static int _split_bf_int                (data_t*, data_t[], eval_context_t*);

// Basic operations
static int _min_farr  (data_t*, data_t[], eval_context_t*); // (float[]) -> float
static int _max_farr  (data_t*, data_t[], eval_context_t*); // (float[]) -> float
static int _sum_farr  (data_t*, data_t[], eval_context_t*); // (float[]) -> float

// Logical operators for custom types
static int _not  (data_t*, data_t[], eval_context_t*); // -> atomset
static int _and  (data_t*, data_t[], eval_context_t*); // -> atomset
static int _or   (data_t*, data_t[], eval_context_t*); // -> atomset
static int _xor  (data_t*, data_t[], eval_context_t*); // -> atomset

// Boolean operators
static int _cmp_lt  (data_t*, data_t[], eval_context_t*); // -> bool
static int _cmp_gt  (data_t*, data_t[], eval_context_t*); // -> bool
static int _cmp_le  (data_t*, data_t[], eval_context_t*); // -> bool
static int _cmp_ge  (data_t*, data_t[], eval_context_t*); // -> bool
static int _cmp_eq  (data_t*, data_t[], eval_context_t*); // -> bool
static int _cmp_lt_f_f(data_t*, data_t[], eval_context_t*); // -> bool
static int _cmp_gt_f_f(data_t*, data_t[], eval_context_t*); // -> bool
static int _cmp_le_f_f(data_t*, data_t[], eval_context_t*); // -> bool
static int _cmp_ge_f_f(data_t*, data_t[], eval_context_t*); // -> bool
static int _cmp_eq_f_f(data_t*, data_t[], eval_context_t*); // -> bool

// Selectors
// Atom level selectors
static int _all             (data_t*, data_t[], eval_context_t*); // -> atomset
static int _within_x        (data_t*, data_t[], eval_context_t*); // -> atomset
static int _within_y        (data_t*, data_t[], eval_context_t*); // -> atomset
static int _within_z        (data_t*, data_t[], eval_context_t*); // -> atomset
static int _within_xyz      (data_t*, data_t[], eval_context_t*); // -> atomset
static int _within_expl_flt (data_t*, data_t[], eval_context_t*); // -> atomset
static int _within_expl_frng(data_t*, data_t[], eval_context_t*); // -> atomset
static int _within_impl_flt (data_t*, data_t[], eval_context_t*); // -> atomset
static int _within_impl_frng(data_t*, data_t[], eval_context_t*); // -> atomset
static int _name            (data_t*, data_t[], eval_context_t*); // -> atomset
static int _element_str     (data_t*, data_t[], eval_context_t*); // -> atomset
static int _element_irng    (data_t*, data_t[], eval_context_t*); // -> atomset
static int _atom_irng       (data_t*, data_t[], eval_context_t*); // -> atomset
static int _atom_int        (data_t*, data_t[], eval_context_t*); // -> atomset

static int _ion             (data_t*, data_t[], eval_context_t*); // -> atomset
static int _backbone        (data_t*, data_t[], eval_context_t*); // -> atomset
static int _side            (data_t*, data_t[], eval_context_t*); // -> atomset
static int _sidechain       (data_t*, data_t[], eval_context_t*); // -> atomset
static int _nucleoside      (data_t*, data_t[], eval_context_t*); // -> atomset
static int _nucleobase      (data_t*, data_t[], eval_context_t*); // -> atomset
static int _head            (data_t*, data_t[], eval_context_t*); // -> atomset
static int _tail            (data_t*, data_t[], eval_context_t*); // -> atomset[]

static int _ring            (data_t*, data_t[], eval_context_t*); // -> atomset[]
static int _ring_irng       (data_t*, data_t[], eval_context_t*); // (irange[], ...) -> atomset[]

// Structure matching, see md_match.h
static int _smiles          (data_t*, data_t[], eval_context_t*); // (str, level: str, mode: str) -> atomset[]
static int _match           (data_t*, data_t[], eval_context_t*); // (atomset, by: str, level: str, mode: str) -> atomset[]

// Residue level selectors
static int _water       (data_t*, data_t[], eval_context_t*);   // -> atomset[]
static int _lipid       (data_t*, data_t[], eval_context_t*);   // -> atomset[]
static int _protein     (data_t*, data_t[], eval_context_t*);   // -> atomset[]
static int _nucleic     (data_t*, data_t[], eval_context_t*);   // -> atomset[]
static int _polymer     (data_t*, data_t[], eval_context_t*);   // -> atomset[]
static int _protein_irng(data_t*, data_t[], eval_context_t*);   // (irange[]) -> atomset[]
static int _nucleic_irng(data_t*, data_t[], eval_context_t*);   // (irange[]) -> atomset[]
static int _polymer_irng(data_t*, data_t[], eval_context_t*);   // (irange[]) -> atomset[]
static int _nucleotide  (data_t*, data_t[], eval_context_t*);   // -> atomset[]

// Regions of a QM calculation (MD_ATOM_FLAG_QM), one selection per residue which has atoms in it
static int _qm          (data_t*, data_t[], eval_context_t*);   // -> atomset[]
static int _environment (data_t*, data_t[], eval_context_t*);   // -> atomset[]
static int _res_name    (data_t*, data_t[], eval_context_t*);   // (str[])          -> atomset[]
static int _res_seq_id (data_t*, data_t[], eval_context_t*);    // (int[]/irange[]) -> atomset[]
static int _res         (data_t*, data_t[], eval_context_t*);   // (irange[])       -> atomset[]

static int _fill_res(data_t*, data_t[], eval_context_t*);    // (atomset[]) -> atomset
static int _fill_chain(data_t*, data_t[], eval_context_t*);  // (atomset[]) -> atomset[]

// Instance level selectors
static int _chain_irng(data_t*, data_t[], eval_context_t*);     // (irange[]) -> atomset[]
static int _chain_id(data_t*, data_t[], eval_context_t*);       // (str[]) -> atomset
static int _molecule_irng(data_t*, data_t[], eval_context_t*);  // (irange[]) -> atomset[]
static int _fill_molecule(data_t*, data_t[], eval_context_t*);  // (atomset[]) -> atomset[]
static int _part_irng(data_t*, data_t[], eval_context_t*);      // (irange[]) -> atomset[]
static int _fill_part(data_t*, data_t[], eval_context_t*);      // (atomset[]) -> atomset[]
static int _glycan(data_t*, data_t[], eval_context_t*);         // -> atomset[]

// Property Compute
static int _distance        (data_t*, data_t[], eval_context_t*); // (position[], position[]) -> float
static int _distance_min    (data_t*, data_t[], eval_context_t*); // (position[N], position[]) -> float[N]
static int _distance_max    (data_t*, data_t[], eval_context_t*); // (position[N], position[]) -> float[N]
static int _distance_pair   (data_t*, data_t[], eval_context_t*); // (position[N], position[M]) -> float[N*M]

static int _angle   (data_t*, data_t[], eval_context_t*); // (position[], position[], position[]) -> float
static int _dihedral(data_t*, data_t[], eval_context_t*); // (position[], position[], position[]), position[]) -> float

static int _rmsd    (data_t*, data_t[], eval_context_t*); // (atomset[N]) -> float[N]

// Radial distribution function: The idea is that we use a fixed (high) amount of bins, then we let the user choose some kernel to smooth it.
static int _rdf_flt (data_t*, data_t[], eval_context_t*); // (position[], position[], float)  -> float[1024] (Histogram).
static int _rdf_frng(data_t*, data_t[], eval_context_t*); // (position[], position[], frange) -> float[1024] (Histogram).

// Density along the axis of the unit cell
static int _density  (data_t*, data_t[], eval_context_t*);
static int _dist_scale(data_t*, data_t[], eval_context_t*); // (distribution[], float) -> distribution[]
static int _dist_div  (data_t*, data_t[], eval_context_t*); // (distribution[], float) -> distribution[]
static int _dist_pool (data_t*, data_t[], eval_context_t*); // (distribution[], distribution[]) -> distribution[]
static int _density_x(data_t*, data_t[], eval_context_t*);
static int _density_y(data_t*, data_t[], eval_context_t*);
static int _density_z(data_t*, data_t[], eval_context_t*);

static int _sdf     (data_t*, data_t[], eval_context_t*); // (atomset, atomset, float) -> float[128][128][128]. This one cannot be stored explicitly as one copy per frame, but is rather accumulated.

// Misc
static int _contact_count(data_t*, data_t[], eval_context_t*); // (atomset[], atomset[], float) -> int
static int _contacts(data_t*, data_t[], eval_context_t*);        // (a, b, cutoff, exclude_bonds, min_separation, parent, exclude_within) -> contact
static int _count_contact(data_t*, data_t[], eval_context_t*);   // (contact) -> float
static int _count_contact_unit(data_t*, data_t[], eval_context_t*); // (contact, str) -> float
static int _degree_contact(data_t*, data_t[], eval_context_t*);  // (contact) -> float[num_a]
static int _chunks(data_t*, data_t[], eval_context_t*);          // (atomset[], int) -> atomset[]
static int _porosity   (data_t*, data_t[], eval_context_t*); // (atomset[]) -> float

// Bonds, see ### BONDS ###
static int _bond        (data_t*, data_t[], eval_context_t*); // (a, b, kind, order, resonance, origin) -> bondset[]
static int _bond_atoms  (data_t*, data_t[], eval_context_t*); // atom(bondset[]) -> atomset[]
static int _count_bonds (data_t*, data_t[], eval_context_t*); // (bondset[]) -> float
static int _bond_length (data_t*, data_t[], eval_context_t*); // (bondset[]) -> float[]
static int _bond_not    (data_t*, data_t[], eval_context_t*); // (bondset[]) -> bondset
static int _bond_and    (data_t*, data_t[], eval_context_t*); // (bondset[], bondset[]) -> bondset
static int _bond_or     (data_t*, data_t[], eval_context_t*); // (bondset[], bondset[]) -> bondset
static int _bond_xor    (data_t*, data_t[], eval_context_t*); // (bondset[], bondset[]) -> bondset
static int _bond_delta  (data_t*, data_t[], eval_context_t*); // (scope, tolerance, breakable) -> bond_delta
static int _bond_delta_formed(data_t*, data_t[], eval_context_t*); // (bond_delta) -> bondset
static int _bond_delta_broken(data_t*, data_t[], eval_context_t*); // (bond_delta) -> bondset
static int _bond_of_delta    (data_t*, data_t[], eval_context_t*); // bond(bond_delta, b, kind, order, resonance, origin) -> bondset

// Geometric operations
static int _com     (data_t*, data_t[], eval_context_t*);  // (position[]) -> float[3]
static int _plane_groups    (data_t*, data_t[], eval_context_t*);  // (atomset[]) -> plane[], one per group
static int _plane_points    (data_t*, data_t[], eval_context_t*);  // (position[]) -> plane, one through the positions
static int _plane_centroid  (data_t*, data_t[], eval_context_t*);  // (plane[]) -> float[][3]
static int _plane_normal    (data_t*, data_t[], eval_context_t*);  // (plane[]) -> float[][3]
static int _plane_radius    (data_t*, data_t[], eval_context_t*);  // (plane[]) -> float[]
static int _plane_deviation (data_t*, data_t[], eval_context_t*);  // (plane[]) -> float[]
static int _plane_angle     (data_t*, data_t[], eval_context_t*);  // (plane[], plane[]) -> float[]
static int _pucker_rings    (data_t*, data_t[], eval_context_t*);  // (atomset[]) -> pucker[], one per ring
static int _pucker_atoms    (data_t*, data_t[], eval_context_t*);  // (int[]) -> pucker[1], the atoms in order
static int _pucker_amplitude(data_t*, data_t[], eval_context_t*);  // (pucker[]) -> float[]
static int _pucker_phase    (data_t*, data_t[], eval_context_t*);  // (pucker[]) -> float[]
static int _pucker_theta    (data_t*, data_t[], eval_context_t*);  // (pucker[]) -> float[]
static int _pucker_pseudorotation(data_t*, data_t[], eval_context_t*);  // (pucker[]) -> float[]
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
static int _normalize     (data_t*, data_t[], eval_context_t*); // (float[]) -> float[]
static int _cross         (data_t*, data_t[], eval_context_t*); // (float[3], float[3]) -> float[3]
static int _length        (data_t*, data_t[], eval_context_t*); // (float[]) -> float
static int _mat4_mul_mat4 (data_t*, data_t[], eval_context_t*); // (float[4][4], float[4][4]) -> float[4][4]
static int _mat4_mul_vec4 (data_t*, data_t[], eval_context_t*); // (float[4][4], float[4]) -> float[4]

static int _vec2 (data_t*, data_t[], eval_context_t*); // (float, float) -> float[2]
static int _vec3 (data_t*, data_t[], eval_context_t*); // (float, float, float) -> float[4]
static int _vec4 (data_t*, data_t[], eval_context_t*); // (float, float, float, float) -> float[4]

// Misc
// Count the number of bits set (popcount)
static int _count (data_t*, data_t[], eval_context_t*); // (atomset) -> float

// Specialized count function with supplied argument
static int _count_with_arg(data_t* , data_t[], eval_context_t*); // (atomset) -> float

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


// @TODO: Add your values here
static constant_t constants[] = {
    {STR_INIT("PI"),     _PI,  {.base = {.dim.angle = 1 }, .mult = 1.0}},
    {STR_INIT("TAU"),    _TAU, {.base = {.dim.angle = 1 }, .mult = 1.0}},
    {STR_INIT("E"),      _E,   {.mult = 1.0}},
};

// IMPLICIT CASTS/CONVERSIONS
static procedure_t casts[] = {
    {STR_INIT("cast"),    TI_FLOAT,         1,  {TI_INT},           _cast_int_to_flt},
    {STR_INIT("cast"),    TI_IRANGE,        1,  {TI_INT},           _cast_int_to_irng},
    {STR_INIT("cast"),    TI_FRANGE,        1,  {TI_IRANGE},        _cast_irng_to_frng},
    {STR_INIT("cast"),    TI_IRANGE_ARR,    1,  {TI_INT_ARR},       _cast_int_arr_to_irng_arr,  FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("cast"),    TI_FLOAT_ARR,     1,  {TI_INT_ARR},       _cast_int_arr_to_flt_arr,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("cast"),    TI_FRANGE_ARR,    1,  {TI_IRANGE_ARR},    _cast_irng_arr_to_frng_arr, FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("cast"),    TI_ATOMSET,      1,  {TI_INT_ARR},       _cast_int_arr_to_bf},
    {STR_INIT("cast"),    TI_ATOMSET,      1,  {TI_IRANGE_ARR},    _cast_irng_arr_to_bf},
    {STR_INIT("cast"),    TI_ATOMSET,      1,  {TI_ATOMSET_ARR},  _join_bf_arr},

    //{STR_INIT("cast"),    TI_FLOAT3_ARR,    1,  {TI_COORDINATE_ARR},  _coordinate, FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_DYNAMIC_LENGTH},
};

static procedure_t operators[] = {
    {STR_INIT("not"),    TI_BOOL,           1,  {TI_BOOL},                  _op_not_b},
    {STR_INIT("or"),     TI_BOOL,           2,  {TI_BOOL,   TI_BOOL},       _op_or_b_b},
    {STR_INIT("xor"),    TI_BOOL,           2,  {TI_BOOL,   TI_BOOL},       _op_xor_b_b},
    {STR_INIT("and"),    TI_BOOL,           2,  {TI_BOOL,   TI_BOOL},       _op_and_b_b},

    // Element-wise logic on bool arrays
    {STR_INIT("not"),    TI_BOOL_ARR,       1,  {TI_BOOL_ARR},                  _op_not_barr,       FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("or"),     TI_BOOL_ARR,       2,  {TI_BOOL_ARR, TI_BOOL},         _op_or_barr_b,      FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_SYMMETRIC_ARGS},
    {STR_INIT("xor"),    TI_BOOL_ARR,       2,  {TI_BOOL_ARR, TI_BOOL},         _op_xor_barr_b,     FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_SYMMETRIC_ARGS},
    {STR_INIT("and"),    TI_BOOL_ARR,       2,  {TI_BOOL_ARR, TI_BOOL},         _op_and_barr_b,     FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_SYMMETRIC_ARGS},
    {STR_INIT("or"),     TI_BOOL_ARR,       2,  {TI_BOOL_ARR, TI_BOOL_ARR},     _op_or_barr_barr,   FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},
    {STR_INIT("xor"),    TI_BOOL_ARR,       2,  {TI_BOOL_ARR, TI_BOOL_ARR},     _op_xor_barr_barr,  FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},
    {STR_INIT("and"),    TI_BOOL_ARR,       2,  {TI_BOOL_ARR, TI_BOOL_ARR},     _op_and_barr_barr,  FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},

    // BITFIELD NOT
    {STR_INIT("not"),    TI_ATOMSET,       1,  {TI_ATOMSET_ARR},                  _not,   FLAG_FLATTEN},
    {STR_INIT("and"),    TI_ATOMSET,       2,  {TI_ATOMSET_ARR, TI_ATOMSET_ARR}, _and,   FLAG_FLATTEN},
    {STR_INIT("or"),     TI_ATOMSET,       2,  {TI_ATOMSET_ARR, TI_ATOMSET_ARR}, _or,    FLAG_FLATTEN},
    {STR_INIT("xor"),    TI_ATOMSET,       2,  {TI_ATOMSET_ARR, TI_ATOMSET_ARR}, _xor,   FLAG_FLATTEN},

    // BONDSET LOGIC, on the union of each side
    {STR_INIT("not"),    TI_BONDSET,       1,  {TI_BONDSET_ARR},                  _bond_not},
    {STR_INIT("and"),    TI_BONDSET,       2,  {TI_BONDSET_ARR, TI_BONDSET_ARR},  _bond_and},
    {STR_INIT("or"),     TI_BONDSET,       2,  {TI_BONDSET_ARR, TI_BONDSET_ARR},  _bond_or},
    {STR_INIT("xor"),    TI_BONDSET,       2,  {TI_BONDSET_ARR, TI_BONDSET_ARR},  _bond_xor},

    // Comparison. The scalar overload comes first and the array/array overload last, so that a scalar operand is
    // matched against the broadcasting overloads before it is treated as a length-one array.
    {STR_INIT("<"),     TI_BOOL,           2,  {TI_FLOAT,     TI_FLOAT},     _cmp_lt_f_f},
    {STR_INIT("<"),     TI_BOOL_ARR,       2,  {TI_FLOAT_ARR, TI_FLOAT},     _cmp_lt,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("<"),     TI_BOOL_ARR,       2,  {TI_FLOAT,     TI_FLOAT_ARR}, _cmp_lt,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("<"),     TI_BOOL_ARR,       2,  {TI_FLOAT_ARR, TI_FLOAT_ARR}, _cmp_lt,   FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},

    {STR_INIT(">"),     TI_BOOL,           2,  {TI_FLOAT,     TI_FLOAT},     _cmp_gt_f_f},
    {STR_INIT(">"),     TI_BOOL_ARR,       2,  {TI_FLOAT_ARR, TI_FLOAT},     _cmp_gt,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT(">"),     TI_BOOL_ARR,       2,  {TI_FLOAT,     TI_FLOAT_ARR}, _cmp_gt,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT(">"),     TI_BOOL_ARR,       2,  {TI_FLOAT_ARR, TI_FLOAT_ARR}, _cmp_gt,   FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},

    {STR_INIT("<="),    TI_BOOL,           2,  {TI_FLOAT,     TI_FLOAT},     _cmp_le_f_f},
    {STR_INIT("<="),    TI_BOOL_ARR,       2,  {TI_FLOAT_ARR, TI_FLOAT},     _cmp_le,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("<="),    TI_BOOL_ARR,       2,  {TI_FLOAT,     TI_FLOAT_ARR}, _cmp_le,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("<="),    TI_BOOL_ARR,       2,  {TI_FLOAT_ARR, TI_FLOAT_ARR}, _cmp_le,   FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},

    {STR_INIT(">="),    TI_BOOL,           2,  {TI_FLOAT,     TI_FLOAT},     _cmp_ge_f_f},
    {STR_INIT(">="),    TI_BOOL_ARR,       2,  {TI_FLOAT_ARR, TI_FLOAT},     _cmp_ge,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT(">="),    TI_BOOL_ARR,       2,  {TI_FLOAT,     TI_FLOAT_ARR}, _cmp_ge,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT(">="),    TI_BOOL_ARR,       2,  {TI_FLOAT_ARR, TI_FLOAT_ARR}, _cmp_ge,   FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},

    {STR_INIT("=="),    TI_BOOL,           2,  {TI_FLOAT,     TI_FLOAT},     _cmp_eq_f_f},
    {STR_INIT("=="),    TI_BOOL_ARR,       2,  {TI_FLOAT_ARR, TI_FLOAT},     _cmp_eq,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("=="),    TI_BOOL_ARR,       2,  {TI_FLOAT,     TI_FLOAT_ARR}, _cmp_eq,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("=="),    TI_BOOL_ARR,       2,  {TI_FLOAT_ARR, TI_FLOAT_ARR}, _cmp_eq,   FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},

    // Binary add
    {STR_INIT("+"),      TI_FLOAT,          2,  {TI_FLOAT,      TI_FLOAT},              _op_add_f_f},
    {STR_INIT("+"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT},              _op_add_farr_f,         FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_SYMMETRIC_ARGS},
    {STR_INIT("+"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},          _op_add_farr_farr,      FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},


    {STR_INIT("+"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_VOLUME},                _op_simd_add_farr_farr, FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},
    {STR_INIT("+"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_FLOAT},                 _op_simd_add_farr_f,    FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_SYMMETRIC_ARGS},

    {STR_INIT("+"),      TI_INT,            2,  {TI_INT,        TI_INT},                _op_add_i_i},
    {STR_INIT("+"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT},                _op_add_iarr_i,         FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_SYMMETRIC_ARGS},
    {STR_INIT("+"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT_ARR},            _op_add_iarr_iarr,      FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},

    // Unary negation
    {STR_INIT("-"),      TI_FLOAT,          1,  {TI_FLOAT},                             _op_neg_f},
    {STR_INIT("-"),      TI_FLOAT_ARR,      1,  {TI_FLOAT_ARR},                         _op_neg_farr,           FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("-"),      TI_INT,            1,  {TI_INT},                               _op_neg_i},
    {STR_INIT("-"),      TI_INT_ARR,        1,  {TI_INT_ARR},                           _op_neg_iarr,           FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("-"),      TI_VOLUME,         1,  {TI_VOLUME},                            _op_simd_neg_farr,      FLAG_DEDUCE_LENGTH_FROM_ARG},

    // Binary sub
    {STR_INIT("-"),      TI_FLOAT,          2,  {TI_FLOAT,      TI_FLOAT},              _op_sub_f_f},
    // Not commutative: the scalar on either side has an overload of its own (FLAG_SYMMETRIC_ARGS would swap the operands)
    {STR_INIT("-"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT},              _op_sub_farr_f,         FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("-"),      TI_FLOAT_ARR,      2,  {TI_FLOAT,      TI_FLOAT_ARR},          _op_sub_f_farr,         FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("-"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},          _op_sub_farr_farr,      FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},


    {STR_INIT("-"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_VOLUME},                _op_simd_sub_farr_farr, FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},
    {STR_INIT("-"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_FLOAT},                 _op_simd_sub_farr_f,    FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("-"),      TI_VOLUME,         2,  {TI_FLOAT,   TI_VOLUME},                _op_simd_sub_f_farr,    FLAG_DEDUCE_LENGTH_FROM_ARG},

    {STR_INIT("-"),      TI_INT,            2,  {TI_INT,        TI_INT},                _op_sub_i_i},
    {STR_INIT("-"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT},                _op_sub_iarr_i,         FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("-"),      TI_INT_ARR,        2,  {TI_INT,        TI_INT_ARR},            _op_sub_i_iarr,         FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("-"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT_ARR},            _op_sub_iarr_iarr,      FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},

    // Binary mul
    {STR_INIT("*"),      TI_FLOAT,          2,  {TI_FLOAT,      TI_FLOAT},              _op_mul1_f},
    {STR_INIT("*"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT},              _op_mul1arr_f,         FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_SYMMETRIC_ARGS},
    {STR_INIT("*"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},          _op_mul1arr_farr,      FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},

    // A distribution only scales (its values, not its weights): other arithmetic on it is not defined (see pool)
    {STR_INIT("*"),      TI_DISTRIBUTION_ARR, 2, {TI_DISTRIBUTION_ARR, TI_FLOAT},       _dist_scale,           FLAG_QUERYABLE_LENGTH | FLAG_SYMMETRIC_ARGS},

    {STR_INIT("*"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_VOLUME},                _op_simd_mul1arr_farr, FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},
    {STR_INIT("*"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_FLOAT},                 _op_simd_mul1arr_f,    FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_SYMMETRIC_ARGS},

    {STR_INIT("*"),      TI_INT,            2,  {TI_INT,        TI_INT},                _op_mul_i_i},
    {STR_INIT("*"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT},                _op_mul_iarr_i,         FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_SYMMETRIC_ARGS},
    {STR_INIT("*"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT_ARR},            _op_mul_iarr_iarr,      FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},

    // Binary div
    {STR_INIT("/"),      TI_FLOAT,          2,  {TI_FLOAT,      TI_FLOAT},              _op_div_f_f},
    {STR_INIT("/"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT},              _op_div_farr_f,         FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("/"),      TI_FLOAT_ARR,      2,  {TI_FLOAT,      TI_FLOAT_ARR},          _op_div_f_farr,         FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("/"),      TI_FLOAT_ARR,      2,  {TI_FLOAT_ARR,  TI_FLOAT_ARR},          _op_div_farr_farr,      FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},

    {STR_INIT("/"),      TI_DISTRIBUTION_ARR, 2, {TI_DISTRIBUTION_ARR, TI_FLOAT},       _dist_div,             FLAG_QUERYABLE_LENGTH},
    
    {STR_INIT("/"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_VOLUME},                _op_simd_div_farr_farr, FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},
    {STR_INIT("/"),      TI_VOLUME,         2,  {TI_VOLUME,  TI_FLOAT},                 _op_simd_div_farr_f,    FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("/"),      TI_VOLUME,         2,  {TI_FLOAT,   TI_VOLUME},                _op_simd_div_f_farr,    FLAG_DEDUCE_LENGTH_FROM_ARG},

    {STR_INIT("//"),      TI_INT,            2,  {TI_INT,        TI_INT},                _op_div_i_i},
    {STR_INIT("//"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT},                _op_div_iarr_i,         FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("//"),      TI_INT_ARR,        2,  {TI_INT,        TI_INT_ARR},            _op_div_i_iarr,         FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("//"),      TI_INT_ARR,        2,  {TI_INT_ARR,    TI_INT_ARR},            _op_div_iarr_iarr,      FLAG_DEDUCE_LENGTH_FROM_ARG | FLAG_ARGS_EQUAL_LENGTH},
};

static procedure_t procedures[] = {
    // NATIVE FUNCS
    {STR_INIT("sqrt"),   TI_FLOAT, 1, {TI_FLOAT}, _sqrtf},
    {STR_INIT("cbrt"),   TI_FLOAT, 1, {TI_FLOAT}, _cbrtf},
    {STR_INIT("abs"),    TI_FLOAT, 1, {TI_FLOAT}, _fabsf},
    {STR_INIT("floor"),  TI_FLOAT, 1, {TI_FLOAT}, _floorf},
    {STR_INIT("ceil"),   TI_FLOAT, 1, {TI_FLOAT}, _ceilf},
    {STR_INIT("cos"),    TI_FLOAT, 1, {TI_FLOAT}, _cosf},
    {STR_INIT("sin"),    TI_FLOAT, 1, {TI_FLOAT}, _sinf},
    {STR_INIT("asin"),   TI_FLOAT, 1, {TI_FLOAT}, _asinf},
    {STR_INIT("acos"),   TI_FLOAT, 1, {TI_FLOAT}, _acosf},
    {STR_INIT("atan"),   TI_FLOAT, 1, {TI_FLOAT}, _atanf},
    {STR_INIT("log"),    TI_FLOAT, 1, {TI_FLOAT}, _logf},
    {STR_INIT("exp"),    TI_FLOAT, 1, {TI_FLOAT}, _expf},
    {STR_INIT("log2"),   TI_FLOAT, 1, {TI_FLOAT}, _log2f},
    {STR_INIT("exp2"),   TI_FLOAT, 1, {TI_FLOAT}, _exp2f},
    {STR_INIT("log10"),  TI_FLOAT, 1, {TI_FLOAT}, _log10f},

    {STR_INIT("atan"),   TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _atan2f},
    {STR_INIT("atan2"),  TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _atan2f},
    {STR_INIT("pow"),    TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _powf},
    {STR_INIT("min"),    TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _fminf},
    {STR_INIT("max"),    TI_FLOAT, 2, {TI_FLOAT, TI_FLOAT}, _fmaxf},

    // ARRAY VERSIONS OF NATIVE FUNCS
    {STR_INIT("abs"),    TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_fabsf,    FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("floor"),  TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_floorf,   FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("ceil"),   TI_FLOAT_ARR, 1, {TI_FLOAT_ARR}, _arr_ceilf,    FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("min"),    TI_FLOAT,     1, {TI_FLOAT_ARR}, _min_farr, FLAG_UNIT_OF_ARG},
    {STR_INIT("max"),    TI_FLOAT,     1, {TI_FLOAT_ARR}, _max_farr, FLAG_UNIT_OF_ARG},
    {STR_INIT("sum"),    TI_FLOAT,     1, {TI_FLOAT_ARR}, _sum_farr, FLAG_UNIT_OF_ARG},

    // VECTORIZED VERSIONS FOR ARRAYS OF BIGGER DATA
    {STR_INIT("abs"),    TI_VOLUME, 1,  {TI_VOLUME}, _op_simd_abs_farr},
    {STR_INIT("min"),    TI_VOLUME, 2,  {TI_VOLUME, TI_FLOAT},  _op_simd_min_farr_f,    FLAG_SYMMETRIC_ARGS},
    {STR_INIT("min"),    TI_VOLUME, 2,  {TI_VOLUME, TI_VOLUME}, _op_simd_min_farr_farr},
    {STR_INIT("max"),    TI_VOLUME, 2,  {TI_VOLUME, TI_FLOAT},  _op_simd_max_farr_f,    FLAG_SYMMETRIC_ARGS},
    {STR_INIT("max"),    TI_VOLUME, 2,  {TI_VOLUME, TI_VOLUME}, _op_simd_max_farr_farr},

    // LINEAR ALGEBRA
    {STR_INIT("dot"),       TI_FLOAT,       2, {TI_FLOAT_ARR,   TI_FLOAT_ARR},  _dot},
    {STR_INIT("cross"),     TI_FLOAT3,      2, {TI_FLOAT3,      TI_FLOAT3},     _cross},
    {STR_INIT("length"),    TI_FLOAT,       1, {TI_FLOAT_ARR},                  _length},
    {STR_INIT("normalize"), TI_FLOAT_ARR,   1, {TI_FLOAT_ARR},                  _normalize, FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("mul"),       TI_FLOAT44,     2, {TI_FLOAT44,     TI_FLOAT44},    _mat4_mul_mat4},
    {STR_INIT("mul"),       TI_FLOAT4,      2, {TI_FLOAT44,     TI_FLOAT4},     _mat4_mul_vec4},

    // CONSTRUCTORS
    {STR_INIT("vec2"),   TI_FLOAT2,  2, {TI_FLOAT, TI_FLOAT},                       _vec2},
    {STR_INIT("vec3"),   TI_FLOAT3,  3, {TI_FLOAT, TI_FLOAT, TI_FLOAT},             _vec3},
    {STR_INIT("vec4"),   TI_FLOAT4,  4, {TI_FLOAT, TI_FLOAT, TI_FLOAT, TI_FLOAT},   _vec4},

    // --- SELECTORS ---

    // Atom level
    {STR_INIT("all"),           TI_ATOMSET, 0, {0},                _all},
    {STR_INIT("ion"),           TI_ATOMSET, 0, {0},                _ion},
    {STR_INIT("backbone"),      TI_ATOMSET, 0, {0},                _backbone},
    {STR_INIT("side"),          TI_ATOMSET, 0, {0},                _side},
    {STR_INIT("sidechain"),     TI_ATOMSET, 0, {0},                _sidechain},
    {STR_INIT("nucleoside"),    TI_ATOMSET, 0, {0},                _nucleoside},
    {STR_INIT("nucleobase"),    TI_ATOMSET, 0, {0},                _nucleobase},
    // The parts of lipids (md_lipid_data_t): their heads, and one selection per tail
    {STR_INIT("head"),          TI_ATOMSET, 0, {0},                _head},
    {STR_INIT("tail"),          TI_ATOMSET_ARR, 0, {0},            _tail,          FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE},

    {STR_INIT("type"),      TI_ATOMSET, 1, {TI_STRING_ARR},        _name,          FLAG_STATIC_VALIDATION},
    {STR_INIT("name"),      TI_ATOMSET, 1, {TI_STRING_ARR},        _name,          FLAG_STATIC_VALIDATION},
    {STR_INIT("label"),     TI_ATOMSET, 1, {TI_STRING_ARR},        _name,          FLAG_STATIC_VALIDATION},
    {STR_INIT("element"),   TI_ATOMSET, 1, {TI_STRING_ARR},        _element_str,   FLAG_STATIC_VALIDATION},
    {STR_INIT("element"),   TI_ATOMSET, 1, {TI_IRANGE_ARR},        _element_irng,  FLAG_STATIC_VALIDATION},
    {STR_INIT("atom"),      TI_ATOMSET, 1, {TI_IRANGE_ARR},        _atom_irng,     FLAG_STATIC_VALIDATION},
    {STR_INIT("atom"),      TI_ATOMSET, 1, {TI_INT_ARR},           _atom_int,      FLAG_STATIC_VALIDATION},
    
    {STR_INIT("ring"),      TI_ATOMSET_ARR, 3, {TI_STRING_ARR, TI_STRING_ARR, TI_IRANGE}, _ring, FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},
    // By position: listed after the pattern, so that ring() (no pattern, an absent argument) stays the one above
    {STR_INIT("ring"),      TI_ATOMSET_ARR, 3, {TI_IRANGE_ARR, TI_STRING_ARR, TI_IRANGE}, _ring_irng, FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},

    // Structure matching: one selection per occurrence of a SMILES pattern, or of a reference structure
    {STR_INIT("smiles"),    TI_ATOMSET_ARR, 3, {TI_STRING, TI_STRING, TI_STRING},              _smiles, FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("match"),     TI_ATOMSET_ARR, 4, {TI_ATOMSET, TI_STRING, TI_STRING, TI_STRING}, _match,  FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},

    // Parts by their kind (md_part_kind_t), one set per part: the strands of a polymer, the waters
    {STR_INIT("polymer"),   TI_ATOMSET_ARR, 0, {0},                _polymer,       FLAG_QUERYABLE_LENGTH},
    {STR_INIT("polymer"),   TI_ATOMSET_ARR, 1, {TI_IRANGE_ARR},    _polymer_irng,  FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("protein"),   TI_ATOMSET_ARR, 0, {0},                _protein,       FLAG_QUERYABLE_LENGTH},
    {STR_INIT("protein"),   TI_ATOMSET_ARR, 1, {TI_IRANGE_ARR},    _protein_irng,  FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("nucleic"),   TI_ATOMSET_ARR, 0, {0},                _nucleic,       FLAG_QUERYABLE_LENGTH},
    {STR_INIT("nucleic"),   TI_ATOMSET_ARR, 1, {TI_IRANGE_ARR},    _nucleic_irng,  FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("water"),     TI_ATOMSET_ARR, 0, {0},                _water,         FLAG_QUERYABLE_LENGTH},

    // Residue level
    {STR_INIT("nucleotide"),TI_ATOMSET_ARR, 0, {0},                _nucleotide,    FLAG_QUERYABLE_LENGTH},
    {STR_INIT("lipid"),     TI_ATOMSET_ARR, 0, {0},                _lipid,         FLAG_QUERYABLE_LENGTH},
    {STR_INIT("glycan"),    TI_ATOMSET_ARR, 0, {0},                _glycan,        FLAG_QUERYABLE_LENGTH},

    // A compile error where the system has no such region: a QM region, and atoms which are not in it
    {STR_INIT("qm"),          TI_ATOMSET_ARR, 0, {0},              _qm,            FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("environment"), TI_ATOMSET_ARR, 0, {0},              _environment,   FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("resname"),   TI_ATOMSET_ARR, 1, {TI_STRING_ARR},    _res_name,      FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("residue"),   TI_ATOMSET_ARR, 1, {TI_STRING_ARR},    _res_name,      FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("component"), TI_ATOMSET_ARR, 1, {TI_STRING_ARR},    _res_name,      FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("resid"),     TI_ATOMSET_ARR, 1, {TI_IRANGE_ARR},    _res_seq_id,    FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("residue"),   TI_ATOMSET_ARR, 1, {TI_IRANGE_ARR},    _res,           FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("component"), TI_ATOMSET_ARR, 1, {TI_IRANGE_ARR},    _res,           FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},

    // Instance level
    // A chain id is a label: every atom which carries it, as one set (see md_chain_data_t). chain_id() and auth_id()
    // are other names of it. The position of a chain is that of a structure, one set per chain, as residue(i) is.
    {STR_INIT("chain"),        TI_ATOMSET,      1,  {TI_STRING_ARR},    _chain_id,         FLAG_STATIC_VALIDATION},
    {STR_INIT("chain_id"),     TI_ATOMSET,      1,  {TI_STRING_ARR},    _chain_id,         FLAG_STATIC_VALIDATION},
    {STR_INIT("auth_id"),      TI_ATOMSET,      1,  {TI_STRING_ARR},    _chain_id,         FLAG_STATIC_VALIDATION},
    {STR_INIT("chain"),        TI_ATOMSET_ARR,  1,  {TI_IRANGE_ARR},    _chain_irng,       FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    // Molecules (md_molecule_data_t): by position, one set each
    {STR_INIT("molecule"),     TI_ATOMSET_ARR,  1,  {TI_IRANGE_ARR},    _molecule_irng,    FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    // Parts (md_part_data_t): a strand, a glycan tree, a ligand, an ion, a water; by position, one set each
    {STR_INIT("part"),         TI_ATOMSET_ARR,  1,  {TI_IRANGE_ARR},    _part_irng,        FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},

    // Dynamic selectors (depend on atomic position, therefore marked as dynamic which means the values cannot be determined at compile-time)
    // Also have variable result (well its a single atomset, but the number of atoms within is not fixed)
    {STR_INIT("within_x"),      TI_ATOMSET, 1, {TI_FRANGE},    _within_x,     FLAG_DYNAMIC},
    {STR_INIT("within_y"),      TI_ATOMSET, 1, {TI_FRANGE},    _within_y,     FLAG_DYNAMIC},
    {STR_INIT("within_z"),      TI_ATOMSET, 1, {TI_FRANGE},    _within_z,     FLAG_DYNAMIC},
    {STR_INIT("within_xyz"),    TI_ATOMSET, 3, {TI_FRANGE, TI_FRANGE, TI_FRANGE}, _within_xyz,     FLAG_DYNAMIC},

    {STR_INIT("within"),    TI_ATOMSET, 2, {TI_FLOAT,  TI_COORDINATE_ARR}, _within_expl_flt,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_FLATTEN | FLAG_SPATIAL_QUERY },
    {STR_INIT("within"),    TI_ATOMSET, 2, {TI_FRANGE, TI_COORDINATE_ARR}, _within_expl_frng,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_FLATTEN | FLAG_SPATIAL_QUERY },

    {STR_INIT("within"),    TI_ATOMSET, 1, {TI_FLOAT},                     _within_impl_flt,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_FLATTEN | FLAG_SPATIAL_QUERY },
    {STR_INIT("within"),    TI_ATOMSET, 1, {TI_FRANGE},                    _within_impl_frng,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_FLATTEN | FLAG_SPATIAL_QUERY },

    // --- PROPERTY COMPUTE ---
    {STR_INIT("distance"),          TI_FLOAT,       2,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR}, _distance,          FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_FLATTEN },
    {STR_INIT("distance_min"),      TI_FLOAT_ARR,   2,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR}, _distance_min,      FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("distance_max"),      TI_FLOAT_ARR,   2,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR}, _distance_max,      FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("distance_pair"),     TI_FLOAT_ARR,   2,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR}, _distance_pair,     FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE | FLAG_QUERYABLE_LENGTH },

    {STR_INIT("angle"),     TI_FLOAT,   3,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR, TI_COORDINATE_ARR},                      _angle,     FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE },
    {STR_INIT("dihedral"),  TI_FLOAT,   4,  {TI_COORDINATE_ARR, TI_COORDINATE_ARR, TI_COORDINATE_ARR, TI_COORDINATE_ARR},   _dihedral,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE },

    {STR_INIT("rmsd"),      TI_FLOAT_ARR, 1, {TI_ATOMSET_ARR}, _rmsd,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE},

    {STR_INIT("rdf"),       TI_DISTRIBUTION, 3, {TI_COORDINATE_ARR, TI_COORDINATE_ARR, TI_FLOAT},  _rdf_flt,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE },
    {STR_INIT("rdf"),       TI_DISTRIBUTION, 3, {TI_COORDINATE_ARR, TI_COORDINATE_ARR, TI_FRANGE}, _rdf_frng,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE },

    {STR_INIT("density"),   {TYPE_DISTRIBUTION, {3,MD_DIST_BINS}}, 1, {TI_ATOMSET_ARR}, _density,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_FLATTEN},
    {STR_INIT("density_x"), TI_DISTRIBUTION, 1, {TI_ATOMSET_ARR}, _density_x,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_FLATTEN},
    {STR_INIT("density_y"), TI_DISTRIBUTION, 1, {TI_ATOMSET_ARR}, _density_y,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_FLATTEN},
    {STR_INIT("density_z"), TI_DISTRIBUTION, 1, {TI_ATOMSET_ARR}, _density_z,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_FLATTEN},
    {STR_INIT("pool"),      TI_DISTRIBUTION_ARR, 2, {TI_DISTRIBUTION_ARR, TI_DISTRIBUTION_ARR}, _dist_pool, FLAG_QUERYABLE_LENGTH | FLAG_ARGS_EQUAL_LENGTH},

    {STR_INIT("sdf"),       TI_VOLUME, 3,   {TI_ATOMSET_ARR, TI_ATOMSET, TI_FLOAT},      _sdf,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_SDF | FLAG_VISUALIZE },

    {STR_INIT("count"),     TI_FLOAT,  1,   {TI_ATOMSET_ARR},              _count},
    {STR_INIT("count"),     TI_FLOAT,  2,   {TI_ATOMSET_ARR, TI_STRING},   _count_with_arg, FLAG_STATIC_VALIDATION},

    {STR_INIT("porosity"), TI_FLOAT,  1,   {TI_ATOMSET},   _porosity, FLAG_DYNAMIC | FLAG_STATIC_VALIDATION},

    {STR_INIT("contact_count"), TI_FLOAT_ARR,  3,   {TI_ATOMSET_ARR, TI_ATOMSET, TI_FLOAT}, _contact_count, FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE},

    // Contacts between groups, see md_contact.h
    {STR_INIT("contacts"),  TI_CONTACT,    7,   {TI_ATOMSET_ARR, TI_ATOMSET_ARR, TI_FLOAT, TI_INT, TI_INT, TI_ATOMSET_ARR, TI_ATOMSET_ARR}, _contacts, FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE},
    {STR_INIT("count"),     TI_FLOAT,      1,   {TI_CONTACT},               _count_contact},
    {STR_INIT("count"),     TI_FLOAT,      2,   {TI_CONTACT, TI_STRING},    _count_contact_unit, FLAG_STATIC_VALIDATION},
    {STR_INIT("degree"),    TI_FLOAT_ARR,  1,   {TI_CONTACT},               _degree_contact, FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("chunks"),    TI_ATOMSET_ARR, 2, {TI_ATOMSET_ARR, TI_INT},  _chunks, FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},

    // --- GEOMETRICAL OPERATIONS ---
    {STR_INIT("com"),           TI_FLOAT3,      1,  {TI_COORDINATE_ARR},  _com,           FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_VISUALIZE },
    // One plane per group, or one through positions. The groups come first: an atomset is a position as well, and the
    // first overload which takes the arguments as they are wins. One plane through all the groups is plane(union(...)).
    {STR_INIT("plane"),         TI_PLANE_ARR,   1,  {TI_ATOMSET_ARR},    _plane_groups,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE },
    {STR_INIT("plane"),         TI_PLANE_ARR,   1,  {TI_COORDINATE_ARR},  _plane_points,  FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE },
    {STR_INIT("centroid"),      TI_FLOAT3_ARR,  1,  {TI_PLANE_ARR},       _plane_centroid,  FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("normal"),        TI_FLOAT3_ARR,  1,  {TI_PLANE_ARR},       _plane_normal,    FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("radius"),        TI_FLOAT_ARR,   1,  {TI_PLANE_ARR},       _plane_radius,    FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("deviation"),     TI_FLOAT_ARR,   1,  {TI_PLANE_ARR},       _plane_deviation, FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("angle"),         TI_FLOAT_ARR,   2,  {TI_PLANE_ARR, TI_PLANE_ARR}, _plane_angle, FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("pucker"),        TI_PUCKER_ARR,  1,  {TI_ATOMSET_ARR},    _pucker_rings,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE },
    {STR_INIT("pucker"),        TI_PUCKER_ARR,  1,  {TI_INT_ARR},         _pucker_atoms,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE },
    {STR_INIT("amplitude"),     TI_FLOAT_ARR,   1,  {TI_PUCKER_ARR},      _pucker_amplitude, FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("phase"),         TI_FLOAT_ARR,   1,  {TI_PUCKER_ARR},      _pucker_phase,     FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("theta"),         TI_FLOAT_ARR,   1,  {TI_PUCKER_ARR},      _pucker_theta,     FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("pseudorotation"), TI_FLOAT_ARR,  1,  {TI_PUCKER_ARR},      _pucker_pseudorotation, FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("shape_weights"), TI_FLOAT3_ARR,  1,  {TI_COORDINATE_ARR},  _shape_weights, FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE},

    {STR_INIT("coord"),      TI_FLOAT3_ARR,  1,  {TI_COORDINATE_ARR},  _coordinate,      FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE },
    {STR_INIT("coord_x"),    TI_FLOAT_ARR,   1,  {TI_COORDINATE_ARR},  _coordinate_x,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE },
    {STR_INIT("coord_y"),    TI_FLOAT_ARR,   1,  {TI_COORDINATE_ARR},  _coordinate_y,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE },
    {STR_INIT("coord_z"),    TI_FLOAT_ARR,   1,  {TI_COORDINATE_ARR},  _coordinate_z,    FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE },

    {STR_INIT("coord_xy"),   TI_FLOAT2_ARR,  1,  {TI_COORDINATE_ARR},  _coordinate_xy,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE },
    {STR_INIT("coord_xz"),   TI_FLOAT2_ARR,  1,  {TI_COORDINATE_ARR},  _coordinate_xz,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE },
    {STR_INIT("coord_yz"),   TI_FLOAT2_ARR,  1,  {TI_COORDINATE_ARR},  _coordinate_yz,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH | FLAG_VISUALIZE },

    // --- MISC ---
    {STR_INIT("union"),     TI_ATOMSET,     1,  {TI_ATOMSET_ARR},      _join_bf_arr,  FLAG_FLATTEN},
    // The old name of union: it compiles, with a warning (as flatten of selections does, see static_check_flatten)
    {STR_INIT("join"),      TI_ATOMSET,     1,  {TI_ATOMSET_ARR},      _join_deprecated, FLAG_FLATTEN | FLAG_STATIC_VALIDATION},
    {STR_INIT("flatten"),   TI_ATOMSET,     1,  {TI_ATOMSET_ARR},      _join_bf_arr,  FLAG_FLATTEN},

    {STR_INIT("split"),     TI_ATOMSET_ARR, 2,  {TI_ATOMSET, TI_INT},  _split_bf_int, FLAG_FLATTEN | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH},

    {STR_INIT("residue"),   TI_ATOMSET_ARR, 1,  {TI_ATOMSET_ARR},  _fill_res,   FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("chain"),     TI_ATOMSET_ARR, 1,  {TI_ATOMSET_ARR},  _fill_chain, FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("molecule"),  TI_ATOMSET_ARR, 1,  {TI_ATOMSET_ARR},  _fill_molecule, FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },
    {STR_INIT("part"),      TI_ATOMSET_ARR, 1,  {TI_ATOMSET_ARR},  _fill_part,     FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH },

    // --- BONDS ---
    // The atoms overload first: bond() without arguments is that one
    {STR_INIT("bond"),       TI_BONDSET_ARR, 6, {TI_ATOMSET_ARR, TI_ATOMSET_ARR, TI_STRING_ARR, TI_IRANGE, TI_STRING_ARR, TI_STRING_ARR}, _bond, FLAG_QUERYABLE_LENGTH | FLAG_STATIC_VALIDATION},
    {STR_INIT("bond"),       TI_BONDSET,     6, {TI_BOND_DELTA,  TI_ATOMSET_ARR, TI_STRING_ARR, TI_IRANGE, TI_STRING_ARR, TI_STRING_ARR}, _bond_of_delta, FLAG_STATIC_VALIDATION},
    {STR_INIT("atom"),       TI_ATOMSET_ARR, 1, {TI_BONDSET_ARR},   _bond_atoms,    FLAG_DEDUCE_LENGTH_FROM_ARG},
    {STR_INIT("count"),      TI_FLOAT,       1, {TI_BONDSET_ARR},   _count_bonds},
    {STR_INIT("length"),     TI_FLOAT_ARR,   1, {TI_BONDSET_ARR},   _bond_length,   FLAG_DYNAMIC | FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH},
    {STR_INIT("bond_delta"), TI_BOND_DELTA,  3, {TI_ATOMSET_ARR, TI_FLOAT, TI_STRING_ARR}, _bond_delta, FLAG_DYNAMIC | FLAG_STATIC_VALIDATION},
    {STR_INIT("formed"),     TI_BONDSET,     1, {TI_BOND_DELTA},    _bond_delta_formed},
    {STR_INIT("broken"),     TI_BONDSET,     1, {TI_BOND_DELTA},    _bond_delta_broken},
};

// Parameter names, declared per procedure NAME and shared by all of its overloads (see NAMED ARGUMENTS in
// md_script.c). A procedure without an entry here accepts positional arguments only. Declaring names does not
// change positional calls, so a procedure can be given names without breaking existing scripts, provided the
// names are listed in the order its arguments are already taken.
//
// Rules, checked by the unittest 'script.named_args_signature_table':
//   - A name means the same position in every overload of the procedure.
//   - Every parameter which some overload does not take is PARAM_OPTIONAL, and PARAM_OPTIONAL parameters are trailing.
//   - Every overload takes a PARAM_NULLABLE parameter (it receives the absent argument).
//   - As in Python, a parameter which can be given by position and has to be given does not follow one which can be
//     left out: required parameters after a PARAM_DEFAULT or PARAM_NULLABLE one are PARAM_KW_ONLY.
//   - Procedures with FLAG_SYMMETRIC_ARGS get no entry: swapping arguments given by name has no meaning.
//   - Names cannot be keywords of the language (in, of, out, and, or, xor, not).
#define REQ(name) {STR_INIT(name), PARAM_REQUIRED}
#define OPT(name) {STR_INIT(name), PARAM_OPTIONAL}
#define NUL(name) {STR_INIT(name), PARAM_NULLABLE}
#define KW_REQ(name) {STR_INIT(name), PARAM_REQUIRED | PARAM_KW_ONLY}
#define KW_NUL(name) {STR_INIT(name), PARAM_NULLABLE | PARAM_KW_ONLY}
#define KW_INT(name, value) {STR_INIT(name), PARAM_DEFAULT | PARAM_KW_ONLY, TI_INT, {._int = (value)}}
#define KW_STR(name, value) {STR_INIT(name), PARAM_DEFAULT | PARAM_KW_ONLY, TI_STRING, {._string = STR_INIT(value)}}

static const proc_sig_t signatures[] = {
    {STR_INIT("distance"),      2, {REQ("a"), REQ("b")}},
    {STR_INIT("distance_min"),  2, {REQ("a"), REQ("b")}},
    {STR_INIT("distance_max"),  2, {REQ("a"), REQ("b")}},
    {STR_INIT("distance_pair"), 2, {REQ("a"), REQ("b")}},
    // angle(a, b, c) of three positions, angle(a, b) between two sets of planes
    {STR_INIT("angle"),         3, {REQ("a"), REQ("b"), OPT("c")}},
    {STR_INIT("dihedral"),      4, {REQ("a"), REQ("b"), REQ("c"), REQ("d")}},

    // within(radius) takes its centre from the context ('within(3) in x'), within(radius, around) explicitly
    {STR_INIT("within"),        2, {REQ("radius"), OPT("around")}},
    {STR_INIT("within_xyz"),    3, {REQ("x"), REQ("y"), REQ("z")}},

    {STR_INIT("rdf"),           3, {REQ("a"), REQ("b"), REQ("cutoff")}},
    // The structures are superimposed onto the first, the density of target is accumulated around them within extent
    {STR_INIT("sdf"),           3, {REQ("structures"), REQ("target"), REQ("extent")}},

    {STR_INIT("count"),         2, {REQ("sel"), OPT("unit")}},
    {STR_INIT("split"),         2, {REQ("sel"), REQ("parts")}},
    {STR_INIT("contact_count"), 3, {REQ("a"), REQ("b"), REQ("cutoff")}},

    // contacts(a, b = none, *, cutoff, exclude_bonds = 3, min_separation = 0, parent = none, exclude_within = none)
    {STR_INIT("contacts"),      7, {REQ("a"), NUL("b"), KW_REQ("cutoff"), KW_INT("exclude_bonds", 3), KW_INT("min_separation", 0), KW_NUL("parent"), KW_NUL("exclude_within")}},
    {STR_INIT("degree"),        1, {REQ("c")}},
    {STR_INIT("chunks"),        2, {REQ("sel"), REQ("size")}},

    // smiles(pattern, *, level = "molecule", mode = "unique"), match(ref, *, by = "element", level = "molecule", mode = "unique")
    {STR_INIT("smiles"),        3, {REQ("pattern"), KW_STR("level", "molecule"), KW_STR("mode", "unique")}},
    {STR_INIT("match"),         4, {REQ("ref"), KW_STR("by", "element"), KW_STR("level", "molecule"), KW_STR("mode", "unique")}},

    // ring(pattern = none, *, kind = none, size = none), see _ring. The first is a pattern or, in the other overload, an
    // index (_ring_irng): 'pattern' names it in both.
    {STR_INIT("ring"),          3, {NUL("pattern"), KW_NUL("kind"), KW_NUL("size")}},

    // bond(a = none, b = none, *, kind = none, order = none, resonance = none, origin = none), see _bond, and bond(d, ...),
    // whose 'a' is a bond_delta (_bond_of_delta)
    {STR_INIT("bond"),          6, {NUL("a"), NUL("b"), KW_NUL("kind"), KW_NUL("order"), KW_NUL("resonance"), KW_NUL("origin")}},
    // bond_delta(scope = none, *, tolerance = none, breakable = none), see _bond_delta
    {STR_INIT("bond_delta"),    3, {NUL("scope"), KW_NUL("tolerance"), KW_NUL("breakable")}},
};

#undef REQ
#undef OPT
#undef NUL
#undef KW_REQ
#undef KW_NUL
#undef KW_INT
#undef KW_STR

// ### VALUE DOMAINS ###
// The values that string parameters take, for completion (md_script_complete). The domains are read from the system
// with the same accessors as the procedures use to match their arguments (see value_domain_t and
// completion_domain_values in md_script_complete.inl). Like the signatures, an entry is per procedure NAME and
// parameter position, and holds for every overload taking a string there.
// Checked by the unittest 'script.completion_domain_table': every procedure with a string parameter has an entry.
// attr() is parsed rather than looked up in the procedure table, but its path is declared here all the same.
static const param_domain_t param_domains[] = {
    {STR_INIT("name"),      0, DOMAIN_ATOM_NAME},
    {STR_INIT("type"),      0, DOMAIN_ATOM_NAME},
    {STR_INIT("label"),     0, DOMAIN_ATOM_NAME},
    {STR_INIT("element"),   0, DOMAIN_ELEMENT},
    {STR_INIT("resname"),   0, DOMAIN_RES_NAME},
    {STR_INIT("residue"),   0, DOMAIN_RES_NAME},
    {STR_INIT("component"), 0, DOMAIN_RES_NAME},
    {STR_INIT("chain"),     0, DOMAIN_CHAIN_ID},
    {STR_INIT("chain_id"),  0, DOMAIN_CHAIN_ID},
    {STR_INIT("auth_id"),   0, DOMAIN_CHAIN_ID},
    // count(contacts, unit) takes 'atom' or 'group' instead, which cannot be told apart without knowing the type of
    // the first argument. Counting selections is the common case.
    {STR_INIT("count"),     1, DOMAIN_COUNT_UNIT},
    {STR_INIT("attr"),      0, DOMAIN_ATTR_PATH},
    {STR_INIT("smiles"),    0, DOMAIN_SMILES},
    {STR_INIT("smiles"),    1, DOMAIN_MATCH_LEVEL},
    {STR_INIT("smiles"),    2, DOMAIN_MATCH_MODE},
    {STR_INIT("match"),     1, DOMAIN_MATCH_LABEL},
    {STR_INIT("match"),     2, DOMAIN_MATCH_LEVEL},
    {STR_INIT("match"),     3, DOMAIN_MATCH_MODE},
    {STR_INIT("ring"),      0, DOMAIN_SMILES},
    {STR_INIT("ring"),      1, DOMAIN_RING_KIND},
    {STR_INIT("bond"),      2, DOMAIN_BOND_KIND},
    {STR_INIT("bond"),      4, DOMAIN_BOND_RESONANCE},
    {STR_INIT("bond"),      5, DOMAIN_BOND_ORIGIN},
    {STR_INIT("bond_delta"), 2, DOMAIN_BOND_ORIGIN},
};

// The values of the string parameters of smiles() and match()
typedef struct match_name_t {
    str_t name;
    int   value;
} match_name_t;

// "chain": within one chain (md_chain_data_t), the chains given to md_match as its units
#define MATCH_LEVEL_CHAIN -1
static const match_name_t match_levels[] = {
    {STR_INIT("molecule"),  MD_MATCH_LEVEL_MOLECULE},
    {STR_INIT("residue"),   MD_MATCH_LEVEL_RESIDUE},
    {STR_INIT("component"), MD_MATCH_LEVEL_RESIDUE},
    {STR_INIT("part"),      MD_MATCH_LEVEL_PART},
    {STR_INIT("chain"),     MATCH_LEVEL_CHAIN},
};

// 'whole' is UNIQUE with MD_MATCH_FLAG_WHOLE: a match which is all of its unit is the only one of the unit, so the
// other modes come to the same, and ALL would only add the symmetric mappings of the same atoms to a selection.
#define MATCH_MODE_WHOLE -1
static const match_name_t match_modes[] = {
    {STR_INIT("unique"),       MD_MATCH_MODE_UNIQUE},
    {STR_INIT("all"),          MD_MATCH_MODE_ALL},
    {STR_INIT("one_per_unit"), MD_MATCH_MODE_ONE_PER_UNIT},
    {STR_INIT("disjoint"),     MD_MATCH_MODE_DISJOINT},
    {STR_INIT("whole"),        MATCH_MODE_WHOLE},
};

static const match_name_t match_labels[] = {
    {STR_INIT("element"), MD_MATCH_LABEL_ELEMENT},
    {STR_INIT("name"),    MD_MATCH_LABEL_NAME},
};

static inline md_spatial_acc_t* get_spatial_acc(eval_context_t* ctx, double max_cutoff) {
    ASSERT(ctx);
    
    // Built for a cutoff at least as large as asked for: reused, otherwise rebuilt
    if (ctx->spatial_acc_cutoff < max_cutoff) {
        if (ctx->spatial_acc.alloc) {
			md_spatial_acc_free(&ctx->spatial_acc);
        }
        ctx->spatial_acc.alloc = ctx->temp_alloc;

        // Rounded up to a multiple of 6 A, so a few different cutoffs within one evaluation share a structure
		const double cutoff = ceil(max_cutoff / 6.0) * 6.0;
        md_coord_stream_t coords = md_coord_stream_from_aos((const float*)ctx->cur_state->xyz, sizeof(vec3_t), NULL, ctx->sys->atom.count);
        md_spatial_acc_init(&ctx->spatial_acc, &(md_spatial_acc_desc_t){ .coords = &coords, .cutoff = cutoff, .unitcell = &ctx->cur_state->unitcell });
		ctx->spatial_acc_cutoff = cutoff;
    }

    return &ctx->spatial_acc;
}

static inline void visualize_atom_mask(const md_bitset_t* mask, eval_context_t* ctx) {
    ASSERT(ctx->vis);
    if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        md_bitmask_or(ctx->vis_structure, *mask);
    }
}

static inline void visualize_atom_range(irange_t range, eval_context_t* ctx) {
    ASSERT(ctx->vis);
    if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        md_bitmask_set_range(ctx->vis_structure, (uint32_t)range.beg, (uint32_t)range.end);
    }
}

static inline void visualize_atom_index(int64_t index, eval_context_t* ctx) {
    ASSERT(ctx->vis);
    if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        md_bitmask_set_bit(ctx->vis_structure, (uint32_t)index);
    }
}

#define COLOR_WHITE 0xFFFFFFFFU
#define COLOR_GRAY  0x88888888U

#define COLOR_U32(R,G,B,A) (((uint32_t)A << 24) | ((uint32_t)B << 16) | ((uint32_t)G << 8) | ((uint32_t)R))

static inline md_script_vis_vertex_t vertex(vec3_t p, uint32_t color) {
    return (md_script_vis_vertex_t){p, color, vec3_set(0,0,1), 0xFFFFFFFFu};
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

static inline void push_text(vec3_t xyz, const char* cstr, md_script_vis_t* vis) {
    ASSERT(vis);
    ASSERT(vis->alloc);
    md_script_vis_text_t text = {xyz, str_copy_cstr(cstr, vis->alloc) };
    md_array_push(vis->text, text, vis->alloc);
} 

// A label which is a quantity. 'str' is still filled in with the default spelling so the label can
// be drawn without looking at anything else, but the value and its unit are carried alongside for a
// viewer which would rather show it in a unit of its own choosing.
static inline void push_quantity_text(vec3_t xyz, double value, md_unit_t unit, md_script_vis_t* vis) {
    ASSERT(vis);
    ASSERT(vis->alloc);

    char buf[64];
    size_t len = (size_t)md_snprintf(buf, sizeof(buf), "%.2f", value);
    if (len < sizeof(buf) && !md_unit_is_none(unit)) {
        char unit_buf[32];
        const size_t unit_len = md_unit_print(unit_buf, sizeof(unit_buf), unit);
        if (unit_len > 0) {
            // The degree sign hugs the number the way the convention has it, everything else takes a space.
            const char* sep = strcmp(unit_buf, "\xC2\xB0") == 0 ? "" : " ";
            len += (size_t)md_snprintf(buf + len, sizeof(buf) - len, "%s%s", sep, unit_buf);
        }
    }
    len = MIN(len, sizeof(buf) - 1);

    md_script_vis_text_t text = {
        .pos   = xyz,
        .str   = str_copy((str_t){buf, len}, vis->alloc),
        .value = value,
        .unit  = unit,
    };
    md_array_push(vis->text, text, vis->alloc);
}

// Adds an empty structure to the visualization and returns its index. Filled in when it is complete, as an
// md_array may move its elements when it grows.
static inline size_t push_structure(md_script_vis_t* vis) {
    ASSERT(vis);
    ASSERT(vis->alloc);
    md_array_push(vis->structures, md_bitset_empty_set(), vis->alloc);
    return md_array_size(vis->structures) - 1;
}

static inline bool idx_in_range(int idx, irange_t range) {
    // I think a range should be inclusive in this context... Since we should be 1 based and not 0 based on indices
    return range.beg <= idx && idx <= range.end;
}

static inline bool range_in_range(irange_t small_range, irange_t big_range) {
    return big_range.beg <= small_range.beg && small_range.end <= big_range.end;
}

static inline vec3_t extract_com(const vec3_t* xyz, const float* w, const md_bitset_t* bitfield) {
    ASSERT(xyz);
    ASSERT(bitfield);

    vec4_t sum = {0};
    md_bitset_iter_t it = md_bitset_iter(*bitfield);
    while (md_bitset_iter_next(&it)) {
        const int64_t idx = it.idx;
        const float weight = w ? w[idx] : 1.0f;
        vec4_t v = vec4_from_vec3(xyz[idx], 1.0f);
        v = vec4_mul1(v, weight);
        sum = vec4_add(sum, v);
    }
    if (sum.w == 0) sum.w = 1.0f;
    return vec3_div1(vec3_from_vec4(sum), sum.w);
}

static inline vec3_t* extract_vec3(const vec3_t* xyz, const md_bitset_t* bf, md_allocator_i* alloc) {
    ASSERT(xyz);
    ASSERT(bf);
    ASSERT(alloc);

    size_t count = md_bitset_count(*bf);
    if (count == 0) {
        return NULL;
    }

    md_array(vec3_t) result = md_array_create(vec3_t, count, alloc);
    md_bitset_iter_t it = md_bitset_iter(*bf);
    size_t i = 0;
    while (md_bitset_iter_next(&it)) {
        const uint64_t idx = it.idx;
        result[i++] = xyz[idx];
    }
    return result;
}

static inline size_t extract_xyzw_vec4(vec4_t* dst_xyzw, const vec3_t* src_xyz, const float* src_w, const md_bitset_t* bf) {
    ASSERT(dst_xyzw);
    ASSERT(src_xyz);
    ASSERT(bf);

    size_t count = 0;
    
    md_bitset_iter_t it = md_bitset_iter(*bf);
    while (md_bitset_iter_next(&it)) {
        const uint64_t idx = it.idx;
        dst_xyzw[count] = vec4_from_vec3(src_xyz[idx], src_w ? src_w[idx] : 1.0f);
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

// The atoms of a chain: the atom ranges of its runs of residues (see md_chain_data_t), into b
static inline void chain_set_bits(md_bitmask_t* b, const md_system_t* sys, size_t chain_idx) {
    const md_urange_t* end = md_system_chain_run_end(sys, chain_idx);
    for (const md_urange_t* it = md_system_chain_run_beg(sys, chain_idx); it && it < end; ++it) {
        const md_urange_t r = md_system_chain_run_atom_range(sys, *it);
        if (r.beg < r.end) md_bitmask_set_range(b, r.beg, r.end);
    }
}

static inline bool chain_in_context(const md_system_t* sys, size_t chain_idx, const md_bitset_t* bitfield) {
    if (!bitfield) return true;
    const md_urange_t* end = md_system_chain_run_end(sys, chain_idx);
    for (const md_urange_t* it = md_system_chain_run_beg(sys, chain_idx); it && it < end; ++it) {
        const md_urange_t r = md_system_chain_run_atom_range(sys, *it);
        if (md_bitset_count_range(*bitfield, r.beg, r.end)) return true;
    }
    return false;
}

// The residues with an atom in the context (all of them without one), in their order: what residue(index) counts
// within it. Only those within the span of the context [beg, end) can have one, so they are found from a binary search
// rather than a walk over all of the residues, which residue(1) in residue(:) would otherwise make once per residue.
static md_array(md_residue_idx_t) get_res_indices_in_context(const md_system_t* sys, const md_bitset_t* bitfield, md_allocator_i* alloc) {
    ASSERT(sys);
    ASSERT(alloc);

    md_array(int32_t) arr = 0;
    const md_residue_data_t* res = &sys->residue;
    if (!res->count) return arr;

    if (bitfield) {
        if (!res->atom_offset || bitfield->end <= bitfield->beg) return arr;
        // The first residue which ends past the beginning of the context
        size_t lo = 0, hi = res->count;
        while (lo < hi) {
            const size_t mid = lo + (hi - lo) / 2;
            if (res->atom_offset[mid + 1] <= bitfield->beg) lo = mid + 1;
            else hi = mid;
        }
        for (size_t i = lo; i < res->count && res->atom_offset[i] < bitfield->end; ++i) {
            if (md_bitset_count_range(*bitfield, res->atom_offset[i], res->atom_offset[i + 1])) {
                md_array_push(arr, (int32_t)i, alloc);
            }
        }
    }
    else {
        for (int32_t i = 0; i < (int32_t)res->count; ++i) {
            md_array_push(arr, i, alloc);
        }
    }

    return arr;
}

// The chains with an atom in the context (all of them without one), in their order: what chain(index) counts
static md_array(int32_t) get_chain_indices_in_context(const md_system_t* sys, const md_bitset_t* bitfield, md_allocator_i* alloc) {
    ASSERT(sys);
    ASSERT(alloc);

    md_array(int32_t) arr = 0;
    for (size_t c = 0; c < md_system_chain_count(sys); ++c) {
        if (chain_in_context(sys, c, bitfield)) md_array_push(arr, (int32_t)c, alloc);
    }
    return arr;
}

static inline irange_t get_atom_range_in_context(const md_system_t* sys, const md_bitset_t* mol_ctx) {
    ASSERT(sys);
    irange_t range = { 0, (int32_t)sys->atom.count };
    if (mol_ctx) {
        range.beg = mol_ctx->beg;
        range.end = mol_ctx->end;
    }
    return range;
}

// The positions given to an index selector (residue(1:5), ring(3), chain(2), molecule(1:2:7)): 1-based and inclusive,
// an open end runs to the first or the last, a step takes every step-th. out (optional) receives them 0-based, in the
// order given, into the num things the selector counts; returns how many, or STATIC_VALIDATION_ERROR (reported when
// report is set). A position before the first is an error. One past the last is an error without a context, and
// within one selects nothing: the selector counts what lies in the context, and a chain need not have a tenth residue,
// nor a residue a second ring.
static int select_positions(md_array(int32_t)* out, const irange_t* ranges, size_t num_ranges, int64_t num, const char* unit, bool report, eval_context_t* ctx) {
    int count = 0;
    for (size_t i = 0; i < num_ranges; ++i) {
        const int64_t beg  = ranges[i].beg == INT32_MIN ? 1 : ranges[i].beg;
        const int64_t end  = ranges[i].end == INT32_MAX ? num : ranges[i].end;
        const int64_t step = ranges[i].step > 0 ? ranges[i].step : 1;
        if (beg < 1) {
            if (report) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "%s: index %" PRId64 " is before the first %s, they are counted from 1", unit, beg, unit);
            return STATIC_VALIDATION_ERROR;
        }
        if (end > num && !ctx->mol_ctx) {
            if (report) {
                if (num == 0) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "%s: index %" PRId64 " is past the last %s, there are none", unit, end, unit);
                else EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "%s: index %" PRId64 " is past the last %s (%" PRId64 ")", unit, end, unit, num);
            }
            return STATIC_VALIDATION_ERROR;
        }
        for (int64_t k = beg; k <= MIN(end, num); k += step) {
            if (out) md_array_push(*out, (int32_t)(k - 1), ctx->temp_alloc);
            count += 1;
        }
    }
    return count;
}

static bool validate_atom_index_in_context(int32_t in_idx, token_t token, const eval_context_t* ctx) {
    ASSERT(ctx);

    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);
    const int ctx_size = ctx_range.end - ctx_range.beg;
    const int32_t ctx_idx = remap_index_to_context(in_idx, ctx_range);
    if (ctx_idx < ctx_range.beg || ctx_idx >= ctx_range.end) {
        EVAL_LOG_ERROR(ctx, token, "supplied index (%i) is not within range of its context (%i:%i)", in_idx, 1, ctx_size);
        return false;
    }
    if (ctx->mol_ctx) {
        if (!md_bitset_test(*ctx->mol_ctx, ctx_idx)) {
            EVAL_LOG_ERROR(ctx, token, "supplied index (%i) is not represented within the supplied context", in_idx);
            return false;
        }
    }
    return true;
}

static bool validate_atom_range_in_context(irange_t in_range, token_t token, const eval_context_t* ctx) {
    ASSERT(ctx);
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);
    const int ctx_size = ctx_range.end - ctx_range.beg;
    const irange_t range = remap_range_to_context(in_range, ctx_range);
    if (!range_in_range(range, ctx_range)) {
        EVAL_LOG_ERROR(ctx, token, "supplied range (%i:%i) is not within range of its context (%i:%i)", in_range.beg, in_range.end, 1, ctx_size);
        return false;
    }
    if (ctx->mol_ctx) {
        if (!range_in_range(range, ctx_range)) {
            EVAL_LOG_ERROR(ctx, token, "supplied range (%i:%i) is not contained within its context (%i:%i)", in_range.beg, in_range.end, 1, ctx_size);
            return false;
        }
    }
    return true;
}

static int coordinate_validate(data_t arg, int arg_idx, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_COORDINATE_ARR));
    if (element_count(arg) == 0) return 0;
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

    if (ctx->backchannel) {
        ctx->backchannel->unit[0] = md_unit_none();
        ctx->backchannel->unit[1] = md_unit_angstrom();
    }

    switch (arg.type.base_type) {
    case TYPE_FLOAT:
        ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_FLOAT3_ARR));
        return (int)type_info_array_len(arg.type);
    case TYPE_INT: {
        int* indices = as_int_arr(arg);
        int count = 0;
        for (size_t i = 0; i < element_count(arg); ++i) {
            if (!validate_atom_index_in_context(indices[i], ctx->arg_tokens[arg_idx], ctx)) {
                return STATIC_VALIDATION_ERROR;
            }                
            count += 1;
        }
        return count;
    }
    case TYPE_IRANGE: {
        if (ctx->backchannel && ctx->arg_flags[arg_idx] & FLAG_DYNAMIC) {
            ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
        }
        const irange_t* ranges = as_irange_arr(arg);
        size_t num_ranges = element_count(arg);
        int count = 0;
        for (size_t i = 0; i < num_ranges; ++i) {
            if (!validate_atom_range_in_context(ranges[i], ctx->arg_tokens[arg_idx], ctx)) {
                return STATIC_VALIDATION_ERROR;
            }
            const irange_t range = remap_range_to_context(ranges[i], ctx_range);
            const int64_t step = ranges[i].step > 0 ? ranges[i].step : 1;
            if (ctx->mol_ctx) {
                int64_t pop = (int64_t)md_bitset_count_range(*ctx->mol_ctx, range.beg, range.end);
                count += (int)((pop + step - 1) / step);
            } else {
                int64_t len = MAX(0, range.end - range.beg);
                count += (int)((len + step - 1) / step);
            }
        }
        return count;
    }
    case TYPE_ATOMSET:
    {
        if (ctx->backchannel && ctx->arg_flags[arg_idx] & FLAG_DYNAMIC) {
            ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
        }
        if (element_count(arg) == 1) {
            return (int)md_bitset_count(*as_bitset(arg));
        }
        else {
            return (int)type_info_array_len(arg.type);
        }
    }
    default:
        ASSERT(false);
        return STATIC_VALIDATION_ERROR;
    }
}

// The union of the sets. A single set is returned as is (an alias), more are combined into alloc.
static inline md_bitset_t _internal_flatten_bf(const md_bitset_t* bf_arr, size_t count, md_allocator_i* alloc) {
    if (count == 1) {
        return bf_arr[0];
    }

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_bitmask_t b = md_bitmask_create(0, MD_BITSET_MAX_INDEX, md_temp_allocator(temp));
    for (size_t i = 0; i < count; ++i) {
        md_bitmask_or(&b, bf_arr[i]);
    }
    md_bitset_t bf = md_bitmask_extract(&b, alloc);
    md_temp_end(temp);
    return bf;
}

static void coordinate_visualize(data_t arg, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(ctx->vis);

    if (element_count(arg) == 0) return;
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

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
                    const int64_t step = ranges[i].step > 0 ? ranges[i].step : 1;
                    if (ctx->mol_ctx) {
                        md_bitset_iter_t it = md_bitset_iter_range(*ctx->mol_ctx, range.beg, range.end);
                        int64_t k = 0;
                        while (md_bitset_iter_next(&it)) {
                            const int64_t idx = it.idx;
                            if ((k % step) == 0) {
                                visualize_atom_index(idx, ctx);
                            }
                            ++k;
                        }
                    }
                    else {
                        for (int64_t j = range.beg; j < range.end; j += step) {
                            visualize_atom_index(j, ctx);
                        }
                    }
                }
            }
        }
        break;
    case TYPE_ATOMSET:
        if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
            md_bitset_t tmp_bf = {0};
            const md_bitset_t* in_bf = as_bitset(arg);
            size_t bf_len = element_count(arg);


            if (ctx->subscript_ranges) {
                // For atomsets, we only care about the first subscript dimension [0]
                irange_t range = ctx->subscript_ranges[0];
                for (int i = range.beg; i < range.end; ++i) {
                    if (i < (int)bf_len) {
                        const md_bitset_t* bf = &in_bf[i];
                        if (ctx->mol_ctx) {
                            tmp_bf = md_bitset_and(*bf, *ctx->mol_ctx, ctx->temp_alloc);
                            bf = &tmp_bf;
                        }
                        visualize_atom_mask(bf, ctx);
                    }
                }
            } else {
                for (size_t i = 0; i < bf_len; ++i) {
                    const md_bitset_t* bf = &in_bf[i];
                    if (ctx->mol_ctx) {
                        tmp_bf = md_bitset_and(*bf, *ctx->mol_ctx, ctx->temp_alloc);
                        bf = &tmp_bf;
                    }
                    visualize_atom_mask(bf, ctx);
                }
            }

        }
        break;
    default:
        ASSERT(false);
    }
}

static md_array(vec3_t) coordinate_extract(data_t arg, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_COORDINATE_ARR));
    md_array(vec3_t) positions = 0;
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

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
            vec3_t pos = { ctx->cur_state->xyz[idx].x, ctx->cur_state->xyz[idx].y, ctx->cur_state->xyz[idx].z };
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
            const int64_t step = ranges[i].step > 0 ? ranges[i].step : 1;

            if (ctx->mol_ctx) {
                int64_t pop = (int64_t)md_bitset_count_range(*ctx->mol_ctx, (uint64_t)range.beg, (uint64_t)range.end);
                md_array_ensure(positions, (size_t)((pop + step - 1) / step), ctx->temp_alloc);
                md_bitset_iter_t it = md_bitset_iter_range(*ctx->mol_ctx, range.beg, range.end);
                int64_t k = 0;
                while (md_bitset_iter_next(&it)) {
                    const int64_t idx = it.idx;
                    if ((k % step) == 0) {
                        vec3_t pos = { ctx->cur_state->xyz[idx].x, ctx->cur_state->xyz[idx].y, ctx->cur_state->xyz[idx].z };
                        md_array_push(positions, pos, ctx->temp_alloc);
                    }
                    ++k;
                }
            }
            else {
                int64_t len = MAX(0, range.end - range.beg);
                md_array_ensure(positions, (size_t)((len + step - 1) / step), ctx->temp_alloc);
                for (int64_t j = range.beg; j < range.end; j += step) {
                    vec3_t pos = { ctx->cur_state->xyz[j].x, ctx->cur_state->xyz[j].y, ctx->cur_state->xyz[j].z };
                    md_array_push(positions, pos, ctx->temp_alloc);
                }
            }
        }
        break;
    }
    case TYPE_ATOMSET: {
        md_bitset_t tmp_bf = { 0 };

        md_bitset_t* bf_arr = as_bitset(arg);
        size_t num_bf = element_count(arg);

        if (num_bf == 1) {
            md_bitset_t* bf = bf_arr;
            if (ctx->mol_ctx) {
                tmp_bf = md_bitset_and(*bf_arr, *ctx->mol_ctx, ctx->temp_alloc);
                bf = &tmp_bf;
            }
            positions = extract_vec3(ctx->cur_state->xyz, bf, ctx->temp_alloc);
        }
        else {
            for (size_t i = 0; i < num_bf; ++i) {
                md_bitset_t* bf = &bf_arr[i];
                if (ctx->mol_ctx) {
                    tmp_bf = md_bitset_and(*bf, *ctx->mol_ctx, ctx->temp_alloc);
                    bf = &tmp_bf;
                }
                
                if (ctx->proc_flags & FLAG_FLATTEN) {
                    md_bitset_iter_t it = md_bitset_iter(*bf);
                    size_t count = md_bitset_count(*bf);
                    md_array_ensure(positions, md_array_size(positions) + count, ctx->temp_alloc);
                    while (md_bitset_iter_next(&it)) {
                        const uint64_t idx = it.idx;
                        vec3_t pos = vec3_set(ctx->cur_state->xyz[idx].x, ctx->cur_state->xyz[idx].y, ctx->cur_state->xyz[idx].z);
                        md_array_push(positions, pos, ctx->temp_alloc);
                    }
                } else {
                    vec3_t com = extract_com(ctx->cur_state->xyz, ctx->atom_mass, bf);
                    md_array_push(positions, com, ctx->temp_alloc);
                }
            }
            ASSERT(num_bf == md_array_size(positions));
        }

        break;
    }
    default:
        break;
    }

    return positions;
}

static md_array(vec4_t) coordinate_extract_xyzw(data_t arg, float default_weight, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_COORDINATE_ARR));
    md_array(vec4_t) out_xyzw = 0;
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

    switch (arg.type.base_type) {
    case TYPE_FLOAT:
        ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_FLOAT3_ARR));
        size_t num_xyz = element_count(arg);
        md_array_resize(out_xyzw, num_xyz, ctx->temp_alloc);
        const vec3_t* src_xyz = as_vec3_arr(arg);
        for (size_t i = 0; i < num_xyz; ++i) {
            out_xyzw[i] = vec4_from_vec3(src_xyz[i], default_weight); 
        }
        md_array_push_array(out_xyzw, as_vec3_arr(arg), element_count(arg), ctx->temp_alloc);
        break;
    case TYPE_INT: {
        ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_INT_ARR));
        size_t num_idx = element_count(arg);
        md_array_resize(out_xyzw, num_idx, ctx->temp_alloc);
        const int* in_idx = as_int_arr(arg);
        for (size_t i = 0; i < num_idx; ++i) {
            // Shift here since we use 1 based indices for atoms
            const int idx = ctx_range.beg + in_idx[i] - 1;
            out_xyzw[i] = vec4_set(ctx->cur_state->xyz[idx].x, ctx->cur_state->xyz[idx].y, ctx->cur_state->xyz[idx].z, ctx->atom_mass ? ctx->atom_mass[idx] : default_weight);
        }
        break;
    }
    case TYPE_IRANGE: {
        ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_IRANGE_ARR));
        const irange_t* ranges = as_irange_arr(arg);
        size_t num_ranges = element_count(arg);
        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            range = clamp_range(range, ctx_range);
            const int64_t step = ranges[i].step > 0 ? ranges[i].step : 1;

            if (ctx->mol_ctx) {
                md_bitset_iter_t it = md_bitset_iter_range(*ctx->mol_ctx, range.beg, range.end);
                int64_t k = 0;
                while (md_bitset_iter_next(&it)) {
                    const uint64_t idx = it.idx;
                    if ((k % step) == 0) {
                        vec4_t xyzw = vec4_set(ctx->cur_state->xyz[idx].x, ctx->cur_state->xyz[idx].y, ctx->cur_state->xyz[idx].z, ctx->atom_mass ? ctx->atom_mass[idx] : default_weight);
                        md_array_push(out_xyzw, xyzw, ctx->temp_alloc);
                    }
                    ++k;
                }
            }
            else {
                for (int64_t j = range.beg; j < range.end; j += step) {
                    vec4_t xyzw = vec4_set(ctx->cur_state->xyz[j].x, ctx->cur_state->xyz[j].y, ctx->cur_state->xyz[j].z, ctx->atom_mass ? ctx->atom_mass[j] : default_weight);
                    md_array_push(out_xyzw, xyzw, ctx->temp_alloc);
                }
            }
        }
        break;
    }
    case TYPE_ATOMSET: {
        md_bitset_t tmp_bf = { 0 };

        md_bitset_t* bf_arr = as_bitset(arg);
        size_t num_bf = element_count(arg);

        if (num_bf == 1) {
            md_bitset_t* bf = bf_arr;
            if (ctx->mol_ctx) {
                tmp_bf = md_bitset_and(*bf_arr, *ctx->mol_ctx, ctx->temp_alloc);
                bf = &tmp_bf;
            }
            md_array_resize(out_xyzw, md_bitset_count(*bf), ctx->temp_alloc);
            extract_xyzw_vec4(out_xyzw, ctx->cur_state->xyz, ctx->atom_mass, bf);
        }
        else {
            md_array(int32_t) indices = 0;
            for (size_t i = 0; i < num_bf; ++i) {
                md_bitset_t* bf = &bf_arr[i];
                if (ctx->mol_ctx) {
                    tmp_bf = md_bitset_and(*bf, *ctx->mol_ctx, ctx->temp_alloc);
                    bf = &tmp_bf;
                }
                md_array_resize(indices, md_bitset_count(*bf), ctx->temp_alloc);
                md_bitset_iter_t it = md_bitset_iter(*bf);
                md_bitset_iter_extract_indices(indices, md_array_size(indices), it);
                vec3_t xyz = md_util_com_compute(ctx->cur_state->xyz, ctx->atom_mass, indices, md_array_size(indices), &ctx->cur_state->unitcell);
                double w = 0;
                if (ctx->atom_mass) {
                for (size_t j = 0; j < md_array_size(indices); ++j) {
                    w += ctx->atom_mass[indices[j]];
                }
                } else {
                    w = (double)md_array_size(indices);
                }
                md_array_push(out_xyzw, vec4_from_vec3(xyz, (float)w), ctx->temp_alloc);
            }
        }

        break;
    }
    default:
        break;
    }

    return out_xyzw;
}

static md_array(int) coordinate_extract_indices(data_t arg, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_COORDINATE_ARR));
    md_array(int) out_indices = NULL;
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

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
            const int step = ranges[i].step > 0 ? ranges[i].step : 1;
            if (ctx->mol_ctx) {
                md_bitset_iter_t it = md_bitset_iter_range(*ctx->mol_ctx, range.beg, range.end);
                int k = 0;
                while (md_bitset_iter_next(&it)) {
                    const int idx = (int)it.idx;
                    if ((k % step) == 0) {
                        md_array_push(out_indices, (int)idx, ctx->temp_alloc);
                    }
                    ++k;
                }
            }
            else {
                int len = MAX(0, range.end - range.beg);
                md_array_ensure(out_indices, md_array_size(out_indices) + (size_t)((len + step - 1) / step), ctx->temp_alloc);
                for (int j = range.beg; j < range.end; j += step) {
                    md_array_push(out_indices, j, ctx->temp_alloc);
                }
            }
        }
        break;
    }
    case TYPE_ATOMSET: {
        const md_bitset_t* bf = as_bitset(arg);
        const size_t num_bf = element_count(arg);

        if (num_bf == 1) {
            md_bitset_t tmp_bf = { 0 };
            if (ctx->mol_ctx) {
                tmp_bf = md_bitset_and(*bf, *ctx->mol_ctx, ctx->temp_alloc);
                bf = &tmp_bf;
            }

            const size_t count = md_bitset_count(*bf);
            md_array_resize(out_indices, count, ctx->temp_alloc);
            md_bitset_iter_extract_indices(out_indices, count, md_bitset_iter(*bf));

        }
        else {
            // @NOTE(Robin): We want to be consistent with the coordinate_extract function which in the case of multiple atomsets returns the com of each atomset
            // Therefore we do not want in this case to return some index as we cannot map the single position to a single atom
        }

        break;
    }
    default:
        break;
    }

    return out_indices;
}

static vec3_t coordinate_extract_com(data_t arg, eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_COORDINATE_ARR));
    vec3_t com = {0};

    md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

    switch (arg.type.base_type) {
    case TYPE_FLOAT: {
        // This corresponds to user supplied coordinates
        ASSERT(is_type_directly_compatible(arg.type, (type_info_t)TI_FLOAT3_ARR));
        vec3_t* in_pos = as_vec3_arr(arg);
        size_t len = element_count(arg);

        if (len == 1) {
            com = in_pos[0];
            goto done;
        }

        vec4_t* xyzw = md_temp_alloc_array(temp, vec4_t, len);
        for (size_t i = 0; i < len; ++i) {
            xyzw[i] = vec4_from_vec3(in_pos[i], 1.0f);
        }
        com = md_util_com_compute_vec4(xyzw, 0, len, &ctx->cur_state->unitcell);
        break;
    }
    case TYPE_INT: {
        const int32_t* in_idx = as_int_arr(arg);
        size_t num_idx = element_count(arg);
        int32_t* idx = md_temp_alloc_array(temp, int32_t, num_idx);

        for (size_t i = 0; i < num_idx; ++i) {
            idx[i] = remap_index_to_context(in_idx[i], ctx_range);
        }

        if (num_idx == 1) {
            com = ctx->cur_state->xyz[idx[0]];
            goto done;
        }
        
        com = md_util_com_compute(ctx->cur_state->xyz, ctx->atom_mass, idx, num_idx, &ctx->cur_state->unitcell);
        break;
    }
    case TYPE_IRANGE: {
        const irange_t* ranges = as_irange_arr(arg);
        size_t num_ranges = element_count(arg);
        md_array(int32_t) indices = 0;
        // We have multiple ranges which may be overlapping
        // We therefore compute the com for each range and then compute the com from the sub-com
        vec4_t* xyzw_arr = md_temp_alloc(temp, num_ranges * sizeof(vec4_t));
        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            range = clamp_range(range, ctx_range);
            const int64_t step = ranges[i].step > 0 ? ranges[i].step : 1;
            if (ctx->mol_ctx) {
                size_t pop = md_bitset_count_range(*ctx->mol_ctx, range.beg, range.end);
                size_t len = (size_t)((pop + step - 1) / step);
                md_array_ensure(indices, len, ctx->temp_alloc);
                md_array_shrink(indices, 0);
                md_bitset_iter_t it = md_bitset_iter_range(*ctx->mol_ctx, range.beg, range.end);
                int64_t k = 0;
                while (md_bitset_iter_next(&it)) {
                    int32_t idx = (int32_t)it.idx;
                    if ((k % step) == 0) {
                        md_array_push(indices, idx, ctx->temp_alloc);
                    }
                    ++k;
                }
                xyzw_arr[i] = vec4_from_vec3(md_util_com_compute(ctx->cur_state->xyz, ctx->atom_mass, indices, len, &ctx->cur_state->unitcell), 1.0f);
            }
            else {
                size_t len = (size_t)((MAX(0, range.end - range.beg) + step - 1) / step);
                md_array_ensure(indices, len, ctx->temp_alloc);
                md_array_shrink(indices, 0);
                for (int32_t j = range.beg; j < range.end; j += (int32_t)step) {
                    // 1 based indexing to 0 based indexing
                    md_array_push(indices, j, ctx->temp_alloc);
                }
                xyzw_arr[i] = vec4_from_vec3(md_util_com_compute(ctx->cur_state->xyz, ctx->atom_mass, indices, len, &ctx->cur_state->unitcell), 1.0f);
            }
        }

        com = md_util_com_compute_vec4(xyzw_arr, 0, num_ranges, &ctx->cur_state->unitcell);
        break;
    }
    case TYPE_ATOMSET: {
        
        md_bitset_t tmp_bf = { 0 };

        md_bitset_t* in_bf = as_bitset(arg);
        const size_t  num_bf = element_count(arg);

        if (num_bf == 1) {
            md_bitset_t* bf = in_bf;
            if (ctx->mol_ctx) {
                tmp_bf = md_bitset_and(*in_bf, *ctx->mol_ctx, ctx->temp_alloc);
                bf = &tmp_bf;
            }
            size_t len = md_bitset_count(*bf);
            int32_t* indices = md_temp_alloc_array(temp, int32_t, len);
            md_bitset_iter_extract_indices(indices, len, md_bitset_iter(*bf));
            com = md_util_com_compute(ctx->cur_state->xyz, ctx->atom_mass, indices, len, &ctx->cur_state->unitcell);
        }
        else {
            // If we have multiple atomsets we compute the center of mass for each atomset before computing a single com from the sub-coms
            md_array(int32_t) indices = 0;
            vec4_t* xyzw_arr = md_temp_alloc_array(temp, vec4_t, num_bf);
            for (size_t i = 0; i < num_bf; ++i) {
                md_bitset_t* bf = &in_bf[i];
                if (ctx->mol_ctx) {
                    tmp_bf = md_bitset_and(*bf, *ctx->mol_ctx, ctx->temp_alloc);
                    bf = &tmp_bf;
                }
                size_t len = md_bitset_count(*bf);
                md_array_ensure(indices, len, ctx->temp_alloc);
                md_array_shrink(indices, len);
                md_bitset_iter_extract_indices(indices, len, md_bitset_iter(*bf));
                xyzw_arr[i] = vec4_from_vec3(md_util_com_compute(ctx->cur_state->xyz, ctx->atom_mass, indices, len, &ctx->cur_state->unitcell), 1.0f);
            }
            com = md_util_com_compute_vec4(xyzw_arr, 0, num_bf, &ctx->cur_state->unitcell);
        }
        break;
    }
    default:
        break;
    }

done:
    md_temp_end(temp);

    return com;
}

// The names that name(), resname() and chain() are given match what they select by glob_match: '*' stands for any
// run of characters, "C*" (begins with C), "*A" (ends with A), "*A*" (contains A), "C*A". Nothing matches an empty name.
static bool match_query(str_t query, str_t str) {
    if (query.len == 0 || str.len == 0) return false;
    return glob_match(query, str);
}

// Any pattern is one, except the empty one, which matches nothing
static bool validate_query(str_t query) {
    return query.len > 0;
}

// IMPLEMENTATIONS
// @TODO: Add more here

static int _min_farr  (data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
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
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
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

// The sum of the elements, 0 for none. Accumulated in double: a float sum of many terms of mixed sign (energies)
// loses digits that the result keeps.
static int _sum_farr  (data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT_ARR));
    (void)ctx;

    if (dst) {
        const float* src = as_float_arr(arg[0]);
        double sum = 0;
        for (size_t i = 0; i < element_count(arg[0]); ++i) {
            sum += src[i];
        }
        as_float(*dst) = (float)sum;
    }
    return 0;
}

static int _not(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));

    if (dst) {
        const md_bitset_t* bf_src = arg[0].ptr;

        md_bitmask_t b = result_mask(ctx);
        const size_t count = element_count(arg[0]);
        for (size_t i = 0; i < count; ++i) {
            md_bitmask_or(&b, bf_src[i]);
        }
        md_bitmask_invert_range(&b, 0, (uint32_t)ctx->sys->atom.count);

        // This is a bit conceptually strange,
        // But if you perform a bitwise negation within a context, only the bits within the context are flipped
        restrict_to_context(&b, ctx);
        finish_result(dst, 0, &b, ctx);
    }

    return 0;
}

static int _and(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_ATOMSET_ARR));

    if (dst) {
        const md_bitset_t* bf_src[2] = {arg[0].ptr, arg[1].ptr};
        const size_t bf_src_len[2] = {element_count(arg[0]), element_count(arg[1])};

        md_bitmask_t b = result_mask(ctx);
        for (size_t i = 0; i < bf_src_len[0]; ++i) {
            md_bitmask_or(&b, bf_src[0][i]);
        }
        md_bitset_t bf = _internal_flatten_bf(bf_src[1], bf_src_len[1], ctx->temp_alloc);
        md_bitmask_and(&b, bf);
        restrict_to_context(&b, ctx);
        finish_result(dst, 0, &b, ctx);
    }

    return 0;
}

static int _or(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_ATOMSET_ARR));

    if (dst) {
        const md_bitset_t* bf_src[2] = {arg[0].ptr, arg[1].ptr};
        const size_t bf_src_len[2] = {element_count(arg[0]), element_count(arg[1])};

        md_bitmask_t b = result_mask(ctx);
        for (size_t i = 0; i < bf_src_len[0]; ++i) {
            md_bitmask_or(&b, bf_src[0][i]);
        }
        for (size_t i = 0; i < bf_src_len[1]; ++i) {
            md_bitmask_or(&b, bf_src[1][i]);
        }
        restrict_to_context(&b, ctx);
        finish_result(dst, 0, &b, ctx);
    }

    return 0;
}

static int _xor(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_ATOMSET_ARR));

    if (dst) {
        const md_bitset_t* bf_src[2] = { arg[0].ptr, arg[1].ptr };
        const size_t bf_src_len[2] = {element_count(arg[0]), element_count(arg[1])};

        md_bitmask_t b = result_mask(ctx);
        for (size_t i = 0; i < bf_src_len[0]; ++i) {
            md_bitmask_or(&b, bf_src[0][i]);
        }
        md_bitset_t bf = _internal_flatten_bf(bf_src[1], bf_src_len[1], ctx->temp_alloc);
        md_bitmask_xor(&b, bf);
        restrict_to_context(&b, ctx);
        finish_result(dst, 0, &b, ctx);
    }

    return 0;
}

static int _dot(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
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
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT3));
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
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT_ARR));
    (void)ctx;

    if (dst) {
        float* a = (float*)arg[0].ptr;
        double len2 = 0; // Accumulate in double, then cast to float
        for (size_t i = 0; i < arg[0].size / sizeof(float); ++i) {
            len2 += (double)a[i] * (double)a[i];
        }
        as_float(*dst) = (float)sqrt(len2);
    }

    return 0;
}

static int _normalize(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT_ARR));
    (void)ctx;

    if (dst) {
        float* val = (float*)arg[0].ptr;
        double len2 = 0;
        const size_t count = ELEMENTWISE_COUNT(*dst, arg[0], float);
        for (size_t i = 0; i < count; ++i) {
            len2 += (double)val[i] * (double)val[i];
        }
        if (len2 < DBL_EPSILON) {
            MEMSET(dst->ptr, 0, dst->size);
            return 0;
        }
        double scl = 1.0 / sqrt(len2);
        for (size_t i = 0; i < count; ++i) {
            as_float_arr(*dst)[i] = (float)(val[i] * scl);
        }
    }

    return 0;
}

static int _mat4_mul_mat4(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT44));
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
    ASSERT(!dst || is_type_directly_compatible(dst->type,   (type_info_t)TI_FLOAT4));
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
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT2));
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
    ASSERT(!dst || is_type_directly_compatible(dst->type,  (type_info_t)TI_FLOAT3));
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
    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT4));
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
        ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
        (void)arg;
        const md_bitset_t all = md_bitset_range(0, (uint32_t)ctx->sys->atom.count);
        as_bitset(*dst)[0] = ctx->mol_ctx ? md_bitset_and(all, *ctx->mol_ctx, data_bitset_alloc(dst, ctx)) : all;
    }
    return 0;
}

static int _name(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));
    ASSERT(ctx && ctx->sys);

    const str_t* str = as_string_arr(arg[0]);
    const size_t num_str = element_count(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

    if (dst) {
        ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
        md_bitmask_t b = result_mask(ctx);

        for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
            if (ctx->mol_ctx && !md_bitset_test(*ctx->mol_ctx, i)) continue;
            
            str_t atom_id = md_atom_name(&ctx->sys->atom, i);
            for (size_t j = 0; j < num_str; ++j) {
                if (match_query(str[j], atom_id)) {
                    md_bitmask_set_bit(&b, (uint32_t)i);
                }
            }
        }
        finish_result(dst, 0, &b, ctx);
    } else {
        // We are only validating the arguments here, making sure that they are represented within the potential context
        for (size_t j = 0; j < num_str; ++j) {
            if (!validate_query(str[j])) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The string '%.*s' is not a valid query", str[j].len, str[j].ptr);
                return -1;
            }
            bool match = false;
            for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
                if (ctx->mol_ctx && !md_bitset_test(*ctx->mol_ctx, i)) continue;
                str_t atom_id = md_atom_name(&ctx->sys->atom, i);
                if (match_query(str[j], atom_id)) {
                    match = true;
                    break;
                }
            }
            if (!match) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The string '%.*s' did not match any atom label within the structure", str[j].len, str[j].ptr);
                return -1;
            }
        }
    }

    return 0;
}

static int _element_str(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));
    ASSERT(ctx && ctx->sys);

    md_atomic_number_t* query_z = 0;

    const size_t num_str = element_count(arg[0]);
    const str_t* str = as_string_arr(arg[0]);

    for (size_t i = 0; i < num_str; ++i) {
        md_atomic_number_t z = md_atomic_number_from_symbol(str[i], false);
        if (z)
            md_array_push(query_z, z, ctx->temp_alloc);
        else if (!dst) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "Failed to map '%.*s' into any Element.", str[i].len, str[i].ptr);
            return -1;
        }
    }
    const size_t num_query_z = md_array_size(query_z);
    if (!dst && num_query_z == 0) {
        EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "No valid arguments in Element");
        return -1;
    }

    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

    if (dst) {
        ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
        md_bitmask_t b = result_mask(ctx);
        if (ctx_range.beg < ctx_range.end) md_bitmask_reserve(&b, (uint32_t)ctx_range.beg, (uint32_t)ctx_range.end);

        for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
            if (ctx->mol_ctx && !md_bitset_test(*ctx->mol_ctx, i)) continue;
            md_atomic_number_t z = md_atom_atomic_number(&ctx->sys->atom, i);

            for (size_t j = 0; j < num_query_z; ++j) {
                if (query_z[j] == z) {
                    md_bitmask_set_bit(&b, (uint32_t)i);
                    break;
                }
            }
        }
        finish_result(dst, 0, &b, ctx);
    } else {
        for (size_t j = 0; j < num_query_z; ++j) {
            md_atomic_number_t z_j = query_z[j];
            bool found = false;
            for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
                md_atomic_number_t z_i = md_atom_atomic_number(&ctx->sys->atom, i);
                if (z_i == z_j) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "Element '" STR_FMT "' was not found within structure.", STR_ARG(str[j]));
                return -1;
            }
        }
    }

    return 0;
}

static int _element_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    ASSERT(ctx && ctx->sys);

    irange_t* ranges = as_irange_arr(arg[0]);
    const size_t num_ranges = element_count(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

    if (dst) {
        ASSERT(dst->ptr);
        ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
        md_bitmask_t b = result_mask(ctx);
        for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
            if (ctx->mol_ctx && !md_bitset_test(*ctx->mol_ctx, i)) continue;
            md_atomic_number_t z_i = md_atom_atomic_number(&ctx->sys->atom, i);
            for (size_t j = 0; j < num_ranges; ++j) {
                const int64_t step = ranges[j].step > 0 ? ranges[j].step : 1;
                if (idx_in_range(z_i, ranges[j])) {
                    int64_t diff = (int64_t)z_i - (int64_t)ranges[j].beg;
                    if ((diff % step) == 0) {
                        md_bitmask_set_bit(&b, (uint32_t)i);
                        break;
                    }
                }
            }
        }
        finish_result(dst, 0, &b, ctx);
    }
    else {
        for (size_t j = 0; j < num_ranges; ++j) {
            bool match = false;
            const int64_t step = ranges[j].step > 0 ? ranges[j].step : 1;
            for (int64_t i = ctx_range.beg; i < ctx_range.end; ++i) {
                if (ctx->mol_ctx && !md_bitset_test(*ctx->mol_ctx, i)) continue;
                md_atomic_number_t z_i = md_atom_atomic_number(&ctx->sys->atom, i);
                if (idx_in_range(z_i, ranges[j])) {
                    int64_t diff = (int64_t)z_i - (int64_t)ranges[j].beg;
                    if ((diff % step) == 0) {
                        match = true;
                        break;
                    }
                }
            }
            if (!match) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "No element within range (%i:%i) was found within structure.", ranges[j].beg, ranges[j].end);
                return -1;
            }
        }
    }

    return 0;
}

static int coordinate_range(data_t* dst, eval_context_t* ctx, frange_t range_x, frange_t range_y, frange_t range_z) {
    ASSERT(ctx && ctx->sys);
    const md_bitset_t* src_bf = ctx->mol_ctx;

    if (dst) {
        md_bitmask_t b = result_mask(ctx);
        if (src_bf) {
            md_bitset_iter_t it = md_bitset_iter(*src_bf);
            while (md_bitset_iter_next(&it)) {
                const int64_t idx = it.idx;

                if (range_x.beg <= ctx->cur_state->xyz[idx].x && ctx->cur_state->xyz[idx].x <= range_x.end &&
                    range_y.beg <= ctx->cur_state->xyz[idx].y && ctx->cur_state->xyz[idx].y <= range_y.end &&
                    range_z.beg <= ctx->cur_state->xyz[idx].z && ctx->cur_state->xyz[idx].z <= range_z.end)
                {
                    md_bitmask_set_bit(&b, (uint32_t)idx);
                }
            }
        }
        else {
            if (ctx->sys->atom.count) md_bitmask_reserve(&b, 0, (uint32_t)ctx->sys->atom.count);
            for (size_t idx = 0; idx < ctx->sys->atom.count; ++idx) {
                if (range_x.beg <= ctx->cur_state->xyz[idx].x && ctx->cur_state->xyz[idx].x <= range_x.end &&
                    range_y.beg <= ctx->cur_state->xyz[idx].y && ctx->cur_state->xyz[idx].y <= range_y.end &&
                    range_z.beg <= ctx->cur_state->xyz[idx].z && ctx->cur_state->xyz[idx].z <= range_z.end)
                {
                    md_bitmask_set_bit(&b, (uint32_t)idx);
                }
            }
        }
        finish_result(dst, 0, &b, ctx);
        return 0;
    }
    else {
        int count = 0;
        if (src_bf) {
            md_bitset_iter_t it = md_bitset_iter(*src_bf);
            while (md_bitset_iter_next(&it)) {
                const int64_t idx = it.idx;
                if (range_x.beg <= ctx->cur_state->xyz[idx].x && ctx->cur_state->xyz[idx].x <= range_x.end &&
                    range_y.beg <= ctx->cur_state->xyz[idx].y && ctx->cur_state->xyz[idx].y <= range_y.end &&
                    range_z.beg <= ctx->cur_state->xyz[idx].z && ctx->cur_state->xyz[idx].z <= range_z.end)
                {
                    count += 1;
                }
            }
        }
        else {
            for (size_t idx = 0; idx < ctx->sys->atom.count; ++idx) {
                if (range_x.beg <= ctx->cur_state->xyz[idx].x && ctx->cur_state->xyz[idx].x <= range_x.end &&
                    range_y.beg <= ctx->cur_state->xyz[idx].y && ctx->cur_state->xyz[idx].y <= range_y.end &&
                    range_z.beg <= ctx->cur_state->xyz[idx].z && ctx->cur_state->xyz[idx].z <= range_z.end)
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

static void within_float_cb(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    (void)i_idx;    // i idx is external point indcies
    (void)ij_dist2;
    md_bitmask_t* b = (md_bitmask_t*)user_param;
    md_bitmask_set_indices(b, j_idx, num_pairs);
}

static int _within_expl_flt(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys && ctx->cur_state);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));

    const float radius = as_float(arg[0]);

    if (dst || ctx->vis) {
        const vec3_t* in_pos = coordinate_extract(arg[1], ctx);
        const size_t num_pos = md_array_size(in_pos);
        const md_spatial_acc_t* sa = get_spatial_acc(ctx, radius);
        const md_bitset_t* bf_mask = 0;
        md_bitmask_t result = {0};
        md_bitmask_t* bf_dst = 0;
        
        if (is_type_directly_compatible(arg[1].type, (type_info_t)TI_ATOMSET)) {           
            bf_mask = as_bitset(arg[1]);
        }

        md_coord_stream_t stream = md_coord_stream_from_aos((const float*)in_pos, sizeof(vec3_t), NULL, num_pos);
        if (dst) {
            ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
            result = result_mask(ctx);
            bf_dst = &result;
        } else {
            // Visualize
            ASSERT(ctx->vis);
            bf_dst = ctx->vis_structure;
            coordinate_visualize(arg[1], ctx);
            
            for (size_t i = 0; i < num_pos; ++i) {
                push_sphere(in_pos[i], radius, COLOR_WHITE, ctx->vis);
            }
        }

        // Fill atomset with indices
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(sa, &stream, radius, within_float_cb, bf_dst, 0);
            
        if (bf_mask) {
            ASSERT(bf_dst);
            // Exclude the input atoms from the result
            md_bitmask_andnot(bf_dst, *bf_mask);
        }
        if (dst) finish_result(dst, 0, &result, ctx);
    }
    else {
        if (coordinate_validate(arg[1], 1, ctx) < 0) return -1;
        if (radius <= 0) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[1], "The supplied radius is negative or zero, please supply a positive value");
            return -1;
        }
    }

    return 0;
}

static int _within_impl_flt(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys && ctx->cur_state);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT));

    const float radius = as_float(arg[0]);

    if (dst || ctx->vis) {
        if (!ctx->mol_ctx) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "The operation is missing its context");
            return -1;
        }
        md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);
        const md_spatial_acc_t* acc = get_spatial_acc(ctx, radius);

        size_t num_idx = md_bitset_count(*ctx->mol_ctx);
        int32_t* idx = md_temp_alloc_array(temp, int32_t, num_idx);
        size_t num_written = md_bitset_iter_extract_indices(idx, num_idx, md_bitset_iter(*ctx->mol_ctx));
        (void)num_written;
        ASSERT(num_written == num_idx);
            
        md_bitmask_t result = {0};
        md_bitmask_t* bf_dst = 0;

        md_coord_stream_t stream = md_coord_stream_from_aos((const float*)ctx->cur_state->xyz, sizeof(vec3_t), idx, num_idx);
        if (dst) {
            ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
            result = result_mask(ctx);
            bf_dst = &result;
        }
        else {
            // Visualize
            ASSERT(ctx->vis);
            bf_dst = ctx->vis_structure;
            visualize_atom_mask(ctx->mol_ctx, ctx);
            for (size_t i = 0; i < num_idx; ++i) {
                vec3_t pos = ctx->cur_state->xyz[idx[i]];
                push_sphere(pos, radius, COLOR_WHITE, ctx->vis);
            }
        }
        
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(acc, &stream, radius, within_float_cb, bf_dst, 0);

        // Exclude the input atoms from the result
        md_bitmask_andnot(bf_dst, *ctx->mol_ctx);
        if (dst) finish_result(dst, 0, &result, ctx);

        md_temp_end(temp);
    }
    else {
        if (radius <= 0) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The supplied radius is negative or zero, please supply a positive value");
            return -1;
        }
        if (!ctx->mol_ctx) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "The operation is missing its context");
            return -1;
        }
    }

    return 0;
}

typedef struct within_frng_payload_t {
    md_bitmask_t* bf;
    float min_r2;
} within_frng_payload_t;

static void within_frng_cb(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    (void)i_idx;    // i idx is external point indcies
    within_frng_payload_t* data = (within_frng_payload_t*)user_param;
    for (size_t k = 0; k < num_pairs; ++k) {
        if (ij_dist2[k] >= data->min_r2) {
            md_bitmask_set_bit(data->bf, j_idx[k]);
        }
    }
}

static int _within_expl_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys && ctx->cur_state);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FRANGE));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));

    const frange_t rad_range = as_frange(arg[0]);

    if (dst || ctx->vis) {
        const md_spatial_acc_t* sa = get_spatial_acc(ctx, rad_range.end);
        const vec3_t* in_pos = coordinate_extract(arg[1], ctx);
        const size_t num_pos = md_array_size(in_pos);
        const md_bitset_t* bf_mask = 0;
        md_bitmask_t result = {0};
        md_bitmask_t* bf_dst = 0;
            
        if (is_type_directly_compatible(arg[1].type, (type_info_t)TI_ATOMSET)) {           
            bf_mask = as_bitset(arg[1]);
        }
        
        if (dst) {
            ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
            result = result_mask(ctx);
            bf_dst = &result;
        }
        else {
            // Visualize
            bf_dst = ctx->vis_structure;
            const float rad = rad_range.end;
            for (size_t i = 0; i < num_pos; ++i) {
                push_sphere(in_pos[i], rad, COLOR_WHITE, ctx->vis);
            }
        }
        within_frng_payload_t payload = {
            .bf = bf_dst,
            .min_r2 = rad_range.beg * rad_range.beg,
        };
        md_coord_stream_t stream = md_coord_stream_from_aos((const float*)in_pos, sizeof(vec3_t), NULL, num_pos);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(sa, &stream, rad_range.end, within_frng_cb, &payload, 0);

        if (bf_mask) {
            ASSERT(bf_dst);
            // Exclude the input atoms from the result
            md_bitmask_andnot(bf_dst, *bf_mask);
        }
        if (dst) finish_result(dst, 0, &result, ctx);
    }
    else {
        if (coordinate_validate(arg[1], 1, ctx) < 0) return -1;
        if (rad_range.beg < 0 || rad_range.end < rad_range.beg) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[1], "The supplied radius range is invalid, ");
            return -1;
        }
    }

    return 0;
}

static int _within_impl_frng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys && ctx->cur_state);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FRANGE));

    const frange_t rad_range = as_frange(arg[0]);

    if (dst || ctx->vis) {
        if (!ctx->mol_ctx) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "The operation is missing its context");
            return -1;
        }
        md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);
        const md_spatial_acc_t* sa = get_spatial_acc(ctx, rad_range.end);

        size_t num_idx = md_bitset_count(*ctx->mol_ctx);
        int32_t* idx = md_temp_alloc_array(temp, int32_t, num_idx);
        size_t num_written = md_bitset_iter_extract_indices(idx, num_idx, md_bitset_iter(*ctx->mol_ctx));
        (void)num_written;
        ASSERT(num_written == num_idx);

        md_coord_stream_t stream = md_coord_stream_from_aos((const float*)ctx->cur_state->xyz, sizeof(vec3_t), idx, num_idx);
        md_bitmask_t result = {0};
        md_bitmask_t* bf_dst = 0;
        
        if (dst) {
            ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
            result = result_mask(ctx);
            bf_dst = &result;
        }
        else {
            // Visualize
            bf_dst = ctx->vis_structure;
            const float rad = rad_range.end;
            for (size_t i = 0; i < num_idx; ++i) {
                vec3_t pos = ctx->cur_state->xyz[idx[i]];
                push_sphere(pos, rad, COLOR_WHITE, ctx->vis);
            }
        }

        within_frng_payload_t payload = {
            .bf = bf_dst,
            .min_r2 = rad_range.beg * rad_range.beg,
        };

        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(sa, &stream, rad_range.end, within_frng_cb, &payload, 0);
        
        // Exclude the input atoms from the result
        md_bitmask_andnot(bf_dst, *ctx->mol_ctx);
        if (dst) finish_result(dst, 0, &result, ctx);

        md_temp_end(temp);
    } else {
        if (rad_range.beg < 0 || rad_range.end < rad_range.beg) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The supplied radius range is invalid, ");
            return -1;
        }
        if (!ctx->mol_ctx) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "The operation is missing its context");
            return -1;
        }
    }

    return 0;
}

typedef struct contact_count_callback_data_t {
    const md_bitmask_t* exclusion_bf;
    size_t count;
    eval_context_t* ctx;
} contact_count_callback_data_t;

void contact_count_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    (void)i_idx;
    (void)ij_dist2;
	contact_count_callback_data_t* data = (contact_count_callback_data_t*)user_param;
    for (size_t k = 0; k < num_pairs; ++k) {
        if (data->exclusion_bf && md_bitmask_test_bit(data->exclusion_bf, j_idx[k])) continue;
        data->count += 1;
	}
}

void contact_count_vis_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    (void)ij_dist2;
	contact_count_callback_data_t* data = (contact_count_callback_data_t*)user_param;
    eval_context_t* ctx = data->ctx;
    for (size_t k = 0; k < num_pairs; ++k) {
        uint32_t i = i_idx[k];
        uint32_t j = j_idx[k];
        if (data->exclusion_bf && md_bitmask_test_bit(data->exclusion_bf, j)) continue;
        const vec3_t pos_i = md_state_coord(ctx->cur_state, i);
        const vec3_t pos_j = md_state_coord(ctx->cur_state, j);
        md_script_vis_vertex_t v1 = vertex(pos_i, COLOR_U32(0, 0, 0, 255));
        md_script_vis_vertex_t v2 = vertex(pos_j, COLOR_U32(0, 0, 0, 0));
        push_line(v1, v2, ctx->vis);
        visualize_atom_index(i, ctx);
        visualize_atom_index(j, ctx);
    }
}

static int _contact_count(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_ATOMSET_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FLOAT));
    const float cutoff = as_float(arg[2]);
	int path_length = 4; // Default minimum path length

    if (is_type_directly_compatible(arg[3].type, (type_info_t)TI_INT)) {
		path_length = as_int(arg[3]);
        if (path_length < 0) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[3], "The minimum path length cannot be negative");
            return STATIC_VALIDATION_ERROR;
		}
    }

    const md_bitset_t* bf_a = as_bitset(arg[0]);
    const md_bitset_t* bf_b = as_bitset(arg[1]);
	int dim_a = (int)element_count(arg[0]);
	int dim_b = (int)element_count(arg[1]);

    if (dst || ctx->vis) {
        float* out_counts = NULL;
        if (dst) {
            ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
            int dst_len = type_info_array_len(dst->type);
            if (dst_len != dim_a) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "Mismatching size in dst array");
                return -1;
            } 
            out_counts = (float*)dst->ptr;
        }

		md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);
        md_bitset_t tmp_bf = {0};

        if (dim_b > 1) {
            tmp_bf = _internal_flatten_bf(bf_b, dim_b, ctx->temp_alloc);
			bf_b = &tmp_bf;
        }

        size_t num_indices = md_bitset_count(*bf_b);
		int32_t* indices = md_temp_alloc_array(temp, int32_t, num_indices);

		md_bitset_iter_t it = md_bitset_iter(*bf_b);
		md_bitset_iter_extract_indices(indices, num_indices, it);

        md_coord_stream_t stream = md_coord_stream_from_aos((const float*)ctx->cur_state->xyz, sizeof(vec3_t), indices, num_indices);
		md_spatial_acc_t acc = { .alloc = ctx->temp_alloc };
		md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .coords = &stream, .cutoff = cutoff, .unitcell = &ctx->cur_state->unitcell, .flags = MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX });

		md_array(int32_t) a_indices = 0;

        // We must maintain an exclusion atomset to avoid counting atoms within the same set
        // Also we want to avoid counting contacts between bonded atoms
        // How to achieve this efficiently:
        // 1. For each atom in set A, we create an exclusion atomset that includes all atoms in set A
        // 2. We then add all bonded atoms to the exclusion atomset
        // 3. When counting contacts from atom in set A to atoms in set B, we skip those in the exclusion atomset

		// What if set A and set B overlap?
        // 

        md_bitmask_t exclusion_bf = md_bitmask_create(0, atom_domain_end(ctx), ctx->temp_alloc);

        contact_count_callback_data_t data = {
            .exclusion_bf = &exclusion_bf,
            .count = 0,
            .ctx = ctx,
        };

        int beg = 0;
        int end = dim_a;
        if (ctx->vis && ctx->subscript_ranges) {
            // We only care about the first dimension here
            beg = ctx->subscript_ranges[0].beg;
            end = ctx->subscript_ranges[0].end;
        }

		md_spatial_acc_pair_callback_t cb = ctx->vis ? contact_count_vis_callback : contact_count_callback;

        for (int i = beg; i < end; ++i) {
			const md_bitset_t* bf = &bf_a[i];
			const size_t a_length = md_bitset_count(*bf);
			md_array_ensure(a_indices, a_length, ctx->temp_alloc);
			md_bitset_iter_extract_indices(a_indices, a_length, md_bitset_iter(*bf));

            data.count = 0; // Each element of a gets its own count
            md_bitmask_assign(&exclusion_bf, *bf);
            md_bitmask_and(&exclusion_bf, *bf_b); // Start with overlap
			md_util_mask_grow_by_bonds(&exclusion_bf, ctx->sys, path_length, NULL);

            // Iterate over atoms in set A and exclude those potential contact points
            md_coord_stream_t a_stream = md_coord_stream_from_aos((const float*)ctx->cur_state->xyz, sizeof(vec3_t), a_indices, a_length);
            md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &a_stream, cutoff, cb, &data, MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX);
            if (out_counts) {
                out_counts[i] = (float)data.count;
            }
		}

		md_temp_end(temp);
    } else {
        // Validation
        if (cutoff <= 0) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The cutoff distance must be positive.");
            return -1;
        }
        return dim_a;
    }
    return 0;
}

static int _atom_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    ASSERT(ctx && ctx->sys);

    const size_t num_ranges = element_count(arg[0]);
    const irange_t* ranges = as_irange_arr(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

    if (dst) {
        ASSERT(dst->ptr && is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
        // The common case, a single contiguous range, is a RUN: no payload at all
        if (num_ranges == 1 && ranges[0].step <= 1) {
            irange_t range = clamp_range(remap_range_to_context(ranges[0], ctx_range), ctx_range);
            const md_bitset_t set = md_bitset_range((uint32_t)MAX(range.beg, 0), (uint32_t)MAX(range.end, 0));
            as_bitset(*dst)[0] = ctx->mol_ctx ? md_bitset_and(set, *ctx->mol_ctx, data_bitset_alloc(dst, ctx)) : set;
            return 0;
        }

        md_bitmask_t b = result_mask(ctx);
        for (size_t i = 0; i < num_ranges; ++i) {
            irange_t range = remap_range_to_context(ranges[i], ctx_range);
            range = clamp_range(range, ctx_range);
            const int64_t step = ranges[i].step > 0 ? ranges[i].step : 1;
            if (range.beg < range.end) {
                if (step == 1) {
                    md_bitmask_set_range(&b, (uint32_t)range.beg, (uint32_t)range.end);
                } else {
                    for (int64_t j = range.beg; j < range.end; j += step) {
                        md_bitmask_set_bit(&b, (uint32_t)j);
                    }
                }
            }
        }
        // Apply context if supplied
        restrict_to_context(&b, ctx);
        finish_result(dst, 0, &b, ctx);
    }
    else {
        ASSERT(ctx->arg_tokens);
        for (size_t i = 0; i < num_ranges; ++i) {
            if (!validate_atom_range_in_context(ranges[i], ctx->arg_tokens[0], ctx)) {
                return STATIC_VALIDATION_ERROR;
            }
        }
    }

    return 0;
}

static int _atom_int(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_INT_ARR));
    ASSERT(ctx && ctx->sys);

    const size_t num_indices = element_count(arg[0]);
    const int* indices = as_int_arr(arg[0]);
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

    if (dst) {
        ASSERT(dst->ptr && is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
        md_bitmask_t b = result_mask(ctx);
        for (size_t i = 0; i < num_indices; ++i) {
            const int32_t idx = remap_index_to_context(indices[i], ctx_range);
            // Checked by the static validation, against the context known then
            if (!idx_in_range(idx, ctx_range)) {
                EVAL_LOG_ERROR(ctx, ctx->op_token, "supplied index (%i) is not within range of its context", indices[i]);
                return -1;
            }
            md_bitmask_set_bit(&b, (uint32_t)idx);
        }
        finish_result(dst, 0, &b, ctx);
    }
    else {
        ASSERT(ctx->arg_tokens);
        for (size_t i = 0; i < num_indices; ++i) {
            if (!validate_atom_index_in_context(indices[i], ctx->arg_tokens[0], ctx)) {
                return STATIC_VALIDATION_ERROR;
            }
        }
    }

    return 0;
}

// ### STRUCTURE MATCHING ###
// smiles() and match(), see md_match.h. The matches depend on the topology only, so both are static; within an 'in'
// context only the atoms of the context are searched.

static bool match_lookup(int* out, const match_name_t* table, size_t count, str_t str) {
    for (size_t i = 0; i < count; ++i) {
        if (str_eq(table[i].name, str)) {
            *out = table[i].value;
            return true;
        }
    }
    return false;
}

static void match_report_unknown(eval_context_t* ctx, int arg, const char* what, const match_name_t* table, size_t count, str_t str) {
    md_str_writer_t sb = md_str_writer_create(ctx->temp_alloc);
    for (size_t i = 0; i < count; ++i) {
        md_str_writer_fmt(&sb, "%s'" STR_FMT "'", i ? ", " : "", STR_ARG(table[i].name));
    }
    EVAL_LOG_ERROR(ctx, ctx->arg_tokens[arg], "Unknown %s '" STR_FMT "', valid ones are: " STR_FMT, what, STR_ARG(str), STR_ARG(md_str_writer_to_str(&sb)));
}

// Runs the query in the context and fills dst, or returns the number of matches when validating
// The steps a match of a script may take (md_match_desc_t.max_steps). A connected reference or SMILES stays far below
// it (a few 10^5 steps in a system of 10^5 atoms); what goes past it is a query in several parts, each placed anywhere
// in its unit, which fit in more ways than anyone wants listed (see HOW THE SEARCH RUNS in md_match.h). Some 10^7 steps
// are a fraction of a second, which a compilation can wait for and an evaluation can take per frame.
#define MATCH_MAX_STEPS(num_atoms) (10000000ull + 64ull * (uint64_t)(num_atoms))

// The parts of a query which are not bonded to each other
static size_t match_query_parts(const md_match_query_t* q, md_allocator_i* alloc) {
    if (q->num_atoms == 0) return 0;
    int32_t* parent = md_alloc(alloc, sizeof(int32_t) * q->num_atoms);
    for (size_t i = 0; i < q->num_atoms; ++i) parent[i] = (int32_t)i;
    for (size_t i = 0; i < q->num_bonds; ++i) {
        int32_t a = (int32_t)q->bonds[i].a, b = (int32_t)q->bonds[i].b;
        while (parent[a] != a) a = parent[a] = parent[parent[a]];
        while (parent[b] != b) b = parent[b] = parent[parent[b]];
        if (a != b) parent[MAX(a, b)] = MIN(a, b);
    }
    size_t num = 0;
    for (size_t i = 0; i < q->num_atoms; ++i) num += parent[i] == (int32_t)i;
    md_free(alloc, parent, sizeof(int32_t) * q->num_atoms);
    return num;
}

// proc: the procedure, what: what the query was made from, as the error names them ("match", "reference")
static int match_eval(data_t* dst, const md_match_query_t* query, data_t arg[], int level_arg, int mode_arg, uint32_t select_flags, const char* proc, const char* what, eval_context_t* ctx) {
    int level = 0, mode = 0;
    if (!match_lookup(&level, match_levels, ARRAY_SIZE(match_levels), as_string(arg[level_arg]))) {
        if (!dst) match_report_unknown(ctx, level_arg, "level", match_levels, ARRAY_SIZE(match_levels), as_string(arg[level_arg]));
        return STATIC_VALIDATION_ERROR;
    }
    if (!match_lookup(&mode, match_modes, ARRAY_SIZE(match_modes), as_string(arg[mode_arg]))) {
        if (!dst) match_report_unknown(ctx, mode_arg, "mode", match_modes, ARRAY_SIZE(match_modes), as_string(arg[mode_arg]));
        return STATIC_VALIDATION_ERROR;
    }
    if (level == MD_MATCH_LEVEL_RESIDUE && ctx->sys->residue.count == 0) {
        if (!dst) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[level_arg], "The system has no residues to match within");
        return STATIC_VALIDATION_ERROR;
    }
    if (level == MD_MATCH_LEVEL_PART && md_system_part_count(ctx->sys) == 0) {
        if (!dst) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[level_arg], "The system has no parts to match within");
        return STATIC_VALIDATION_ERROR;
    }
    if (level == MATCH_LEVEL_CHAIN && md_system_chain_count(ctx->sys) == 0) {
        if (!dst) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[level_arg], "The data gives no chain ids, so the system has no chains to match within: level='part' matches within each strand");
        return STATIC_VALIDATION_ERROR;
    }
    // Within a chain: its atoms are the unit, the bonds leaving it (a disulfide to another chain) not searched. The
    // atoms of no chain (a molecule the data gives no chain id) are in no unit.
    md_bitset_t* chain_units = NULL;
    size_t num_chain_units = 0;
    if (level == MATCH_LEVEL_CHAIN) {
        num_chain_units = md_system_chain_count(ctx->sys);
        if (num_chain_units == 0) {
            if (!dst) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[level_arg], "The system has no chains to match within");
            return STATIC_VALIDATION_ERROR;
        }
        chain_units = md_alloc(ctx->temp_alloc, sizeof(md_bitset_t) * num_chain_units);
        for (size_t c = 0; c < num_chain_units; ++c) {
            md_bitmask_t b = md_bitmask_create(0, (uint32_t)ctx->sys->atom.count, ctx->temp_alloc);
            chain_set_bits(&b, ctx->sys, c);
            chain_units[c] = md_bitmask_view(&b);
        }
    }

    const md_match_desc_t desc = {
        .query = query,
        .level = level == MATCH_LEVEL_CHAIN ? MD_MATCH_LEVEL_MOLECULE : (md_match_level_t)level,
        .units = chain_units,
        .num_units = num_chain_units,
        .mode  = mode == MATCH_MODE_WHOLE ? MD_MATCH_MODE_UNIQUE : (md_match_mode_t)mode,
        .flags = mode == MATCH_MODE_WHOLE ? MD_MATCH_FLAG_WHOLE : MD_MATCH_FLAG_NONE,
        .mask  = ctx->mol_ctx,
        .max_steps = MATCH_MAX_STEPS(ctx->sys->atom.count),
    };
    md_match_result_t res = {0};
    if (!md_match_find(&res, &desc, ctx->sys, ctx->temp_alloc)) {
        if (!dst) EVAL_LOG_ERROR(ctx, ctx->op_token, "Structure matching failed, see the log");
        return STATIC_VALIDATION_ERROR;
    }
    if (res.incomplete) {
        // Reported by the evaluation as well: the search depends on the frame through the context
        const size_t num_parts = match_query_parts(query, ctx->temp_alloc);
        const char* unit = level == MD_MATCH_LEVEL_RESIDUE ? "residue" : level == MD_MATCH_LEVEL_PART ? "part" : level == MATCH_LEVEL_CHAIN ? "chain" : "molecule";
        if (num_parts > 1) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "%s: the search was given up, the %s is in %zu parts which are not bonded to each other. "
                "Each is placed anywhere in the %s, and together they fit in too many ways: give a %s whose atoms are bonded "
                "(a chain with a gap is two parts)", proc, what, num_parts, unit, what);
        } else {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "%s: the search was given up, the %s fits in too many ways", proc, what);
        }
        md_match_result_free(&res);
        return STATIC_VALIDATION_ERROR;
    }

    int result = 0;
    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        if (dst->ptr && res.count > 0) {
            md_bitset_t* bf_arr = as_bitset(*dst);
            md_allocator_i* alloc = data_bitset_alloc(dst, ctx);
            // Restricted to the context, the selections are made in temp and intersected into alloc
            md_allocator_i* sel_alloc = ctx->mol_ctx ? ctx->temp_alloc : alloc;
            const int64_t cap = type_info_array_len(dst->type);
            size_t n = 0;
            if (cap == 1) {
                // Flattened: all of the matches in one
                bf_arr[0] = md_match_result_select_all(&res, ctx->sys, select_flags, sel_alloc);
                n = 1;
            } else if (cap > 1) {
                n = md_match_result_select(bf_arr, (size_t)cap, &res, ctx->sys, select_flags, sel_alloc);
            }
            if (ctx->mol_ctx) {
                for (size_t i = 0; i < n; ++i) {
                    bf_arr[i] = md_bitset_and(bf_arr[i], *ctx->mol_ctx, alloc);
                }
            }
        }
    } else {
        result = (int)res.count;
        if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
            result = MIN(1, result);
        }
    }

    md_match_result_free(&res);
    return result;
}

static int _smiles(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING));

    const str_t pattern = as_string(arg[0]);
    md_match_query_t query = {0};
    md_smiles_error_t err = {0};
    if (!md_match_query_init_smiles(&query, pattern, ctx->temp_alloc, &err)) {
        if (!dst) {
            // Point at the offending character when the pattern is written out as a string literal
            token_t tok = ctx->arg_tokens[0];
            if (tok.type == TOKEN_STRING && tok.str.len == pattern.len + 2) {
                tok.beg = tok.beg + 1 + (int)err.offset;
                tok.end = MIN(tok.beg + 1, tok.end - 1);
                tok.end = MAX(tok.end, tok.beg);
            }
            EVAL_LOG_ERROR(ctx, tok, "Invalid SMILES: %s", err.message);
        }
        return STATIC_VALIDATION_ERROR;
    }

    // A SMILES describes its atoms with their hydrogens: the selections have them
    const int result = match_eval(dst, &query, arg, 1, 2, MD_MATCH_SELECT_HYDROGENS, "smiles", "SMILES", ctx);
    md_match_query_free(&query);
    return result;
}

static int _match(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET));

    if (!dst && !ctx->vis && (ctx->arg_flags[0] & FLAG_DYNAMIC)) {
        EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "match: the reference cannot depend on the frame");
        return STATIC_VALIDATION_ERROR;
    }

    int label = 0;
    if (!match_lookup(&label, match_labels, ARRAY_SIZE(match_labels), as_string(arg[1]))) {
        if (!dst) match_report_unknown(ctx, 1, "label", match_labels, ARRAY_SIZE(match_labels), as_string(arg[1]));
        return STATIC_VALIDATION_ERROR;
    }

    const md_bitset_t ref = *as_bitset(arg[0]);
    const size_t count = md_bitset_count(ref);
    if (count == 0) {
        if (!dst) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "match: the reference is empty");
        return STATIC_VALIDATION_ERROR;
    }

    // The reference in ascending order of its atoms, which is the order of the atoms of each match
    int32_t* atom_idx = md_alloc(ctx->temp_alloc, sizeof(int32_t) * count);
    const size_t num_idx = md_bitset_iter_extract_indices(atom_idx, count, md_bitset_iter(ref));

    md_match_query_t query = {0};
    int result = STATIC_VALIDATION_ERROR;
    if (md_match_query_init_atoms(&query, atom_idx, num_idx, (md_match_label_t)label, ctx->sys, ctx->temp_alloc)) {
        // The atoms the reference maps, its hydrogens if it has them, and nothing else
        result = match_eval(dst, &query, arg, 2, 3, MD_MATCH_SELECT_NONE, "match", "reference", ctx);
        md_match_query_free(&query);
    } else if (!dst) {
        EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "match: the reference cannot be made into a query, see the log");
    }
    md_free(ctx->temp_alloc, atom_idx, sizeof(int32_t) * count);
    return result;
}

// ### RINGS ###
// ring(pattern = none, *, kind = none, size = none): the rings of the system (md_system_t::ring) which all of the
// given filters let through, each fully within the context. ring(index, *, kind = none, size = none) takes rings by
// their position among those instead (_ring_irng).
//   pattern  SMILES, one or several: a ring which has any of them along it (md_match with the rings as units). A chain
//            is consecutive ring atoms in either direction (CCNCC), a chain closed onto its first atom the whole ring
//            (c1ccncc1), parts separated by '.' are anywhere in the ring (N.O). A branch is not along a ring: an error.
//   kind     One or several of ring_kinds: a ring of any of them.
//   size     The number of atoms, an int or a range (5:6).
typedef enum ring_kind_t {
    RING_KIND_AROMATIC,
    RING_KIND_ALIPHATIC,
    RING_KIND_HETERO,
    RING_KIND_HETEROAROMATIC,
    RING_KIND_SUGAR,
} ring_kind_t;

static const match_name_t ring_kinds[] = {
    {STR_INIT("aromatic"),       RING_KIND_AROMATIC},          // Every bond of the ring aromatic
    {STR_INIT("aliphatic"),      RING_KIND_ALIPHATIC},         // Not aromatic
    {STR_INIT("hetero"),         RING_KIND_HETERO},            // An atom other than C in the ring
    {STR_INIT("heteroaromatic"), RING_KIND_HETEROAROMATIC},    // Both
    {STR_INIT("sugar"),          RING_KIND_SUGAR},             // A furanose or pyranose, see MD_RING_FLAG_SUGAR
};

static bool ring_is_kind(md_ring_flags_t f, int kind) {
    switch (kind) {
    case RING_KIND_AROMATIC:       return (f & MD_RING_FLAG_AROMATIC) != 0;
    case RING_KIND_ALIPHATIC:      return (f & MD_RING_FLAG_AROMATIC) == 0;
    case RING_KIND_HETERO:         return (f & MD_RING_FLAG_HETERO) != 0;
    case RING_KIND_HETEROAROMATIC: return (f & MD_RING_FLAG_AROMATIC) && (f & MD_RING_FLAG_HETERO);
    case RING_KIND_SUGAR:          return (f & MD_RING_FLAG_SUGAR) != 0;
    default:                       return false;
    }
}

// Reports a SMILES syntax error, at the offending character when the pattern is a string literal of its own
static void report_smiles_error(eval_context_t* ctx, int arg, str_t pattern, const md_smiles_error_t* err) {
    token_t tok = ctx->arg_tokens[arg];
    if (tok.type == TOKEN_STRING && tok.str.len == pattern.len + 2) {
        tok.beg = tok.beg + 1 + (int)err->offset;
        tok.end = MIN(tok.beg + 1, tok.end - 1);
        tok.end = MAX(tok.end, tok.beg);
    }
    EVAL_LOG_ERROR(ctx, tok, "Invalid SMILES: %s", err->message);
}

// The query of a ring pattern, or false (reported when validating) if it cannot be along a ring
static bool ring_pattern_query(md_match_query_t* query, str_t pattern, bool report, eval_context_t* ctx) {
    md_smiles_error_t err = {0};
    if (!md_match_query_init_smiles(query, pattern, ctx->temp_alloc, &err)) {
        if (report) report_smiles_error(ctx, 0, pattern, &err);
        return false;
    }
    // Hydrogens written as atoms of their own ([H]C) are hydrogen counts, not ring atoms
    const size_t n = query->num_atoms;
    uint32_t* degree = md_alloc(ctx->temp_alloc, sizeof(uint32_t) * MAX(n, 1));
    MEMSET(degree, 0, sizeof(uint32_t) * MAX(n, 1));
    size_t num_heavy = 0;
    for (size_t i = 0; i < n; ++i) {
        const md_match_atom_t* a = &query->atoms[i];
        num_heavy += !((a->flags & MD_MATCH_ATOM_ELEMENT) && a->z == MD_Z_H);
    }
    for (size_t i = 0; i < query->num_bonds; ++i) {
        const md_match_atom_t* a = &query->atoms[query->bonds[i].a];
        const md_match_atom_t* b = &query->atoms[query->bonds[i].b];
        if (((a->flags & MD_MATCH_ATOM_ELEMENT) && a->z == MD_Z_H) || ((b->flags & MD_MATCH_ATOM_ELEMENT) && b->z == MD_Z_H)) continue;
        degree[query->bonds[i].a]++;
        degree[query->bonds[i].b]++;
    }
    for (size_t i = 0; i < n; ++i) {
        if (degree[i] > 2) {
            if (report) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "ring: '" STR_FMT "' branches, and a ring atom has two neighbors along the ring: use smiles() for substituents", STR_ARG(pattern));
            md_match_query_free(query);
            return false;
        }
    }
    if (num_heavy > MD_RING_MAX_SIZE) {
        if (report) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "ring: '" STR_FMT "' has %zu atoms, more than the largest ring (%d)", STR_ARG(pattern), num_heavy, MD_RING_MAX_SIZE);
        md_match_query_free(query);
        return false;
    }
    return true;
}

// The rings which pass the filters of arg (see above), within the context. Returns the count, or
// STATIC_VALIDATION_ERROR (reported when validating, dst == NULL). out_ring (optional) receives their indices.
static int ring_select(md_array(uint32_t)* out_ring, data_t arg[], bool report, eval_context_t* ctx) {
    const md_system_t* sys = ctx->sys;
    const size_t num_rings = md_ring_data_count(&sys->ring);

    // Kinds
    uint32_t kind_mask = 0;
    if (!is_absent_type(arg[1].type)) {
        ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_STRING_ARR));
        const str_t* kinds = as_string_arr(arg[1]);
        const size_t num_kinds = element_count(arg[1]);
        for (size_t i = 0; i < num_kinds; ++i) {
            int kind = 0;
            if (!match_lookup(&kind, ring_kinds, ARRAY_SIZE(ring_kinds), kinds[i])) {
                if (report) match_report_unknown(ctx, 1, "ring kind", ring_kinds, ARRAY_SIZE(ring_kinds), kinds[i]);
                return STATIC_VALIDATION_ERROR;
            }
            kind_mask |= 1u << kind;
        }
    }

    // Size
    int size_lo = INT32_MIN, size_hi = INT32_MAX;
    if (!is_absent_type(arg[2].type)) {
        ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_IRANGE));
        const irange_t r = as_irange(arg[2]);
        size_lo = r.beg;
        size_hi = r.end;
        if (size_hi < size_lo) {
            if (report) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[2], "ring: the size range is empty");
            return STATIC_VALIDATION_ERROR;
        }
    }

    // Patterns: the rings which have any of them
    uint8_t* has_pattern = NULL;
    if (!is_absent_type(arg[0].type)) {
        ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));
        const str_t* patterns = as_string_arr(arg[0]);
        const size_t num_patterns = element_count(arg[0]);
        has_pattern = md_alloc(ctx->temp_alloc, MAX(num_rings, 1));
        MEMSET(has_pattern, 0, MAX(num_rings, 1));
        // The rings as the units to match within
        md_bitset_t* units = md_alloc(ctx->temp_alloc, sizeof(md_bitset_t) * MAX(num_rings, 1));
        for (size_t r = 0; r < num_rings; ++r) {
            units[r] = md_bitset_from_indices((const uint32_t*)md_ring_atoms(&sys->ring, r), md_ring_size(&sys->ring, r), ctx->temp_alloc);
        }
        for (size_t p = 0; p < num_patterns; ++p) {
            md_match_query_t query = {0};
            if (!ring_pattern_query(&query, patterns[p], report, ctx)) return STATIC_VALIDATION_ERROR;
            if (num_rings > 0) {
                const md_match_desc_t desc = {
                    .query     = &query,
                    .mode      = MD_MATCH_MODE_ONE_PER_UNIT,
                    .units     = units,
                    .num_units = num_rings,
                };
                md_match_result_t res = {0};
                const bool ok = md_match_find(&res, &desc, sys, ctx->temp_alloc);
                md_match_query_free(&query);
                if (!ok) {
                    if (report) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "Ring matching failed, see the log");
                    return STATIC_VALIDATION_ERROR;
                }
                for (size_t i = 0; i < res.count; ++i) {
                    if (res.unit[i] < num_rings) has_pattern[res.unit[i]] = 1;
                }
                md_match_result_free(&res);
            } else {
                md_match_query_free(&query);
            }
        }
    }

    int count = 0;
    for (size_t r = 0; r < num_rings; ++r) {
        const int n = (int)md_ring_size(&sys->ring, r);
        if (n < size_lo || n > size_hi) continue;
        if (has_pattern && !has_pattern[r]) continue;
        if (kind_mask) {
            const md_ring_flags_t f = md_system_ring_flags(sys, r);
            bool any = false;
            for (int k = 0; k < (int)ARRAY_SIZE(ring_kinds) && !any; ++k) {
                any = (kind_mask & (1u << k)) && ring_is_kind(f, k);
            }
            if (!any) continue;
        }
        // Only a ring fully within the context
        if (ctx->mol_ctx) {
            const int32_t* ring = md_ring_atoms(&sys->ring, r);
            bool inside = true;
            for (int k = 0; k < n && inside; ++k) inside = md_bitset_test(*ctx->mol_ctx, (uint32_t)ring[k]);
            if (!inside) continue;
        }
        if (out_ring) md_array_push(*out_ring, (uint32_t)r, ctx->temp_alloc);
        count += 1;
    }
    return count;
}

// The result of ring(): one selection per ring of rings (indices into md_system_t::ring), in their order, or all of them
// in one when flattened; their outlines when visualizing; their count when validating.
static int ring_result(data_t* dst, const uint32_t* rings, int count, eval_context_t* ctx) {
    const md_system_t* sys = ctx->sys;

    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        if (!dst->ptr) return 0;
        md_bitset_t* bf_arr = as_bitset(*dst);
        const int64_t cap = type_info_array_len(dst->type);
        if (cap <= 0) return 0;
        // Flattened (cap 1): all of the rings in one. Rings beyond the capacity go into the last set.
        md_allocator_i* alloc = data_bitset_alloc(dst, ctx);
        md_bitmask_t b = result_mask(ctx);
        int64_t slot = 0;
        for (int i = 0; i < count; ++i) {
            const int64_t s = cap == 1 ? 0 : MIN(i, cap - 1);
            if (s != slot) {
                bf_arr[slot] = md_bitmask_extract(&b, alloc);
                md_bitmask_clear(&b);
                slot = s;
            }
            const int32_t* ring = md_ring_atoms(&sys->ring, rings[i]);
            const size_t n = md_ring_size(&sys->ring, rings[i]);
            for (size_t k = 0; k < n; ++k) md_bitmask_set_bit(&b, (uint32_t)ring[k]);
        }
        if (count > 0) bf_arr[slot] = md_bitmask_extract(&b, alloc);
        md_bitmask_free(&b);
        return 0;
    }

    if (ctx->vis) {
        for (int i = 0; i < count; ++i) {
            const int32_t* ring = md_ring_atoms(&sys->ring, rings[i]);
            const size_t n = md_ring_size(&sys->ring, rings[i]);
            if (n == 0) continue;
            const vec3_t* xyz = ctx->cur_state->xyz;
            md_script_vis_vertex_t vbeg = vertex(vec3_set(xyz[ring[0]].x, xyz[ring[0]].y, xyz[ring[0]].z), COLOR_WHITE);
            md_script_vis_vertex_t v0 = vbeg;
            for (size_t k = 1; k < n; ++k) {
                md_script_vis_vertex_t v1 = vertex(vec3_set(xyz[ring[k]].x, xyz[ring[k]].y, xyz[ring[k]].z), COLOR_WHITE);
                push_line(v0, v1, ctx->vis);
                v0 = v1;
            }
            push_line(v0, vbeg, ctx->vis);
        }
        return 0;
    }

    if (ctx->eval_flags & EVAL_FLAG_FLATTEN) return MIN(1, count);
    return count;
}

static int _ring(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);

    md_array(uint32_t) rings = 0;
    const int count = ring_select(&rings, arg, !dst && !ctx->vis, ctx);
    if (count < 0) return count;
    return ring_result(dst, rings, count, ctx);
}

// ring(index, *, kind = none, size = none): rings by their 1-based position among those ring(kind = kind, size = size)
// gives, in the order of md_system_t::ring: within a context the rings fully within it, as residue(index) counts the
// residues of its context. So without a context or filters, ring(k) is the ring k - 1 of the system, which is what the
// circle of an aromatic ring picks as (md_gl_picking_offset_t). An index before the first is an error; one past the
// last is an error without a context, and selects nothing within one (a residue need not have a second ring).
static int _ring_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    const bool report = !dst && !ctx->vis;

    // The rings to index into: those the filters let through. The index takes the place of the pattern.
    data_t filter_arg[3] = { {0}, arg[1], arg[2] };
    md_array(uint32_t) rings = 0;
    const int num_rings = ring_select(&rings, filter_arg, report, ctx);
    if (num_rings < 0) return num_rings;

    md_array(int32_t) pos = 0;
    const int count = select_positions(&pos, as_irange_arr(arg[0]), element_count(arg[0]), num_rings, "ring", report, ctx);
    if (count < 0) return count;
    uint32_t* picked = md_alloc(ctx->temp_alloc, sizeof(uint32_t) * MAX(count, 1));
    for (int i = 0; i < count; ++i) picked[i] = rings[pos[i]];
    return ring_result(dst, picked, count, ctx);
}

static int _select_atoms_with_flags(data_t* dst, data_t arg[], eval_context_t* ctx, md_atom_flags_t flags) {
    ASSERT(ctx && ctx->sys);
    (void)arg;

    int result = 0;

    if (dst) {
        ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
        md_bitmask_t b = result_mask(ctx);

        if (ctx->mol_ctx) {
            md_bitset_iter_t it = md_bitset_iter(*ctx->mol_ctx);
            while (md_bitset_iter_next(&it)) {
                uint64_t idx = it.idx;
                if (md_atom_flags(&ctx->sys->atom, idx) & flags) {
                    md_bitmask_set_bit(&b, (uint32_t)idx);
                }
            }
        } else {
            for (size_t i = 0; i < ctx->sys->atom.count; ++i) {
                if (md_atom_flags(&ctx->sys->atom, i) & flags) {
                    md_bitmask_set_bit(&b, (uint32_t)i);
                }
            }
        }
        finish_result(dst, 0, &b, ctx);
    }

    return result;
}

static int _select_residues_of_kind(data_t* dst, data_t arg[], eval_context_t* ctx, md_residue_kind_t kind) {
    ASSERT(ctx && ctx->sys);
    (void)arg;

    int result = 0;

    int* res_indices = get_res_indices_in_context(ctx->sys, ctx->mol_ctx, ctx->temp_alloc);
	size_t num_res = md_array_size(res_indices);

    if (dst) {
        ASSERT(dst && is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        const int cap = type_info_array_len(dst->type);
        if (cap == 0) return 0;

        set_writer_t w = set_writer_begin(dst, ctx, false);
        for (size_t i = 0; i < num_res; ++i) {
            int res_idx = res_indices[i];
            if (md_system_residue_kind(ctx->sys, res_idx) == kind) {
                const md_urange_t range = md_system_residue_atom_range(ctx->sys, res_idx);
                set_writer_push_range(&w, range.beg, range.end);
            }
        }
        set_writer_end(&w);
    }
    else {
        if (!ctx->sys->residue.count) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "The system does not contain any components");
            return -1;
        }

        int count = 0;
        for (size_t i = 0; i < num_res; ++i) {
            int32_t res_idx = res_indices[i];
            if (md_system_residue_kind(ctx->sys, res_idx) == kind) {
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

// ### PARTS BY KIND ###
// polymer(), protein(), nucleic(), water(): one set per part of the kind (md_part_kind_t), a strand of a polymer with its
// caps and modified residues, a water. Within a context the parts with an atom in it, whole, as glycan() takes them; by
// position (protein(2)) the parts of the kind in the context, each cut to it, as part(2) counts them. A system whose parts
// were not inferred has a set per residue of the kind instead, as lipid() has without its lipids.

typedef bool (*part_kind_pred_t)(md_part_kind_t kind);
static bool part_kind_polymer(md_part_kind_t kind) { return md_part_kind_is_polymer(kind); }
static bool part_kind_peptide(md_part_kind_t kind) { return kind == MD_PART_KIND_PEPTIDE; }
static bool part_kind_nucleic(md_part_kind_t kind) { return md_part_kind_is_nucleic_acid(kind); }
static bool part_kind_water  (md_part_kind_t kind) { return kind == MD_PART_KIND_WATER; }

// The parts of the kind with an atom in the context (all of them without one), in their order
static md_array(int32_t) get_parts_of_kind(const md_system_t* sys, part_kind_pred_t pred, const md_bitset_t* bitfield, md_allocator_i* alloc) {
    md_array(int32_t) arr = 0;
    for (size_t p = 0; p < md_system_part_count(sys); ++p) {
        if (!pred(md_system_part_kind(sys, p))) continue;
        const md_urange_t r = md_system_part_atom_range(sys, p);
        if (!bitfield || md_bitset_count_range(*bitfield, r.beg, r.end)) md_array_push(arr, (int32_t)p, alloc);
    }
    return arr;
}

static int select_parts_of_kind(data_t* dst, data_t arg[], eval_context_t* ctx, part_kind_pred_t pred, md_residue_kind_t res_kind) {
    ASSERT(ctx && ctx->sys);
    const md_system_t* sys = ctx->sys;
    if (md_system_part_count(sys) == 0) {
        if (res_kind == MD_RESIDUE_KIND_OTHER) {
            // polymer(): no residue kind stands for it
            if (!dst && sys->residue.count == 0) {
                EVAL_LOG_ERROR(ctx, ctx->op_token, "The system does not contain any components");
                return -1;
            }
            return 0;
        }
        return _select_residues_of_kind(dst, arg, ctx, res_kind);
    }

    md_array(int32_t) parts = get_parts_of_kind(sys, pred, ctx->mol_ctx, ctx->temp_alloc);
    const int count = (int)md_array_size(parts);
    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        if (!dst->ptr || type_info_array_len(dst->type) == 0) return 0;
        set_writer_t w = set_writer_begin(dst, ctx, false);
        for (int i = 0; i < count; ++i) {
            const md_urange_t r = md_system_part_atom_range(sys, (size_t)parts[i]);
            set_writer_push_range(&w, r.beg, r.end);
        }
        set_writer_end(&w);
        return 0;
    }
    return (ctx->eval_flags & EVAL_FLAG_FLATTEN) ? MIN(1, count) : count;
}

static int select_parts_of_kind_irng(data_t* dst, data_t arg[], eval_context_t* ctx, part_kind_pred_t pred, const char* unit) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    const md_system_t* sys = ctx->sys;
    if (md_system_part_count(sys) == 0) {
        if (!dst) EVAL_LOG_ERROR(ctx, ctx->op_token, "The system has no parts to count the %ss of", unit);
        return -1;
    }

    md_array(int32_t) parts = get_parts_of_kind(sys, pred, ctx->mol_ctx, ctx->temp_alloc);
    md_array(int32_t) pos = 0;
    const int count = select_positions(dst ? &pos : NULL, as_irange_arr(arg[0]), element_count(arg[0]), (int64_t)md_array_size(parts), unit, !dst, ctx);
    if (count < 0) return count;

    if (dst) {
        if (!dst->ptr || type_info_array_len(dst->type) == 0) return 0;
        set_writer_t w = set_writer_begin(dst, ctx, true);
        for (int i = 0; i < count; ++i) {
            const md_urange_t r = md_system_part_atom_range(sys, (size_t)parts[pos[i]]);
            set_writer_push_range(&w, r.beg, r.end);
        }
        set_writer_end(&w);
        return 0;
    }
    return (ctx->eval_flags & EVAL_FLAG_FLATTEN) ? MIN(1, count) : count;
}

static int _polymer(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return select_parts_of_kind(dst, arg, ctx, part_kind_polymer, MD_RESIDUE_KIND_OTHER);
}

static int _water(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return select_parts_of_kind(dst, arg, ctx, part_kind_water, MD_RESIDUE_KIND_WATER);
}

static int _polymer_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return select_parts_of_kind_irng(dst, arg, ctx, part_kind_polymer, "polymer");
}

static int _protein_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return select_parts_of_kind_irng(dst, arg, ctx, part_kind_peptide, "protein");
}

static int _nucleic_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return select_parts_of_kind_irng(dst, arg, ctx, part_kind_nucleic, "nucleic acid");
}

// lipid(): one selection per lipid molecule (md_lipid_data_t), all of its residues: a lipid split into several residues
// (the AMBER lipids) is one, as a part is. What is attached to a lipid (the further sugars of a glycolipid, the core
// and O antigen of a lipopolysaccharide) is not of it. Within a context the lipids with an atom in it, whole, as
// protein() takes its residues. A system whose lipids were not inferred has a lipid per residue of the kind.
static int _lipid(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    const md_system_t* sys = ctx->sys;
    const md_lipid_data_t* lipid = &sys->lipid;
    if (lipid->count == 0) {
        return _select_residues_of_kind(dst, arg, ctx, MD_RESIDUE_KIND_LIPID);
    }

    // The atoms of each lipid: a run when its residues are consecutive, as they mostly are
    md_array(md_bitset_t) sets = 0;
    for (size_t l = 0; l < lipid->count; ++l) {
        const int32_t* residues = md_lipid_residues(lipid, l);
        const size_t num_residues = md_lipid_residue_count(lipid, l);
        if (num_residues == 0) continue;
        md_bitset_t set = md_bitset_empty_set();
        bool in_ctx = !ctx->mol_ctx;
        md_urange_t run = md_system_residue_atom_range(sys, (size_t)residues[0]);
        for (size_t i = 1; i <= num_residues; ++i) {
            const md_urange_t r = i < num_residues ? md_system_residue_atom_range(sys, (size_t)residues[i]) : (md_urange_t){0, 0};
            if (i < num_residues && r.beg == run.end) {
                run.end = r.end;
                continue;
            }
            in_ctx |= ctx->mol_ctx && md_bitset_count_range(*ctx->mol_ctx, run.beg, run.end) > 0;
            set = md_bitset_empty(set) ? md_bitset_range(run.beg, run.end) : md_bitset_or(set, md_bitset_range(run.beg, run.end), ctx->temp_alloc);
            run = r;
        }
        if (!in_ctx) continue;
        md_array_push(sets, set, ctx->temp_alloc);
    }
    const int count = (int)md_array_size(sets);

    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        if (!dst->ptr || type_info_array_len(dst->type) == 0) return 0;
        set_writer_t w = set_writer_begin(dst, ctx, false);
        for (int i = 0; i < count; ++i) set_writer_push(&w, sets[i]);
        set_writer_end(&w);
        return 0;
    }
    if (ctx->eval_flags & EVAL_FLAG_FLATTEN) return MIN(1, count);
    return count;
}

static int _protein(data_t* dst, data_t arg[], eval_context_t* ctx) {
	return select_parts_of_kind(dst, arg, ctx, part_kind_peptide, MD_RESIDUE_KIND_AMINO_ACID);
}

static int _nucleic(data_t* dst, data_t arg[], eval_context_t* ctx) {
	return select_parts_of_kind(dst, arg, ctx, part_kind_nucleic, MD_RESIDUE_KIND_NUCLEOTIDE);
}

static int _nucleotide(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _select_residues_of_kind(dst, arg, ctx, MD_RESIDUE_KIND_NUCLEOTIDE);
}

static int _ion(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _select_residues_of_kind(dst, arg, ctx, MD_RESIDUE_KIND_ION);
}

// The atoms of the QM region (MD_ATOM_FLAG_QM), or of its environment: one selection per residue with
// atoms in the region, holding those atoms. A residue is not assumed to lie on one side: a QM/MM
// boundary can cut through a residue. Without residues the region is a single selection.
//
// Whether the region exists is a property of the SYSTEM, not of the context: 'qm' needs a QM region,
// 'environment' a QM region and atoms outside it, and without them either is a compile error rather
// than an empty selection - an empty 'environment' would read as "nothing around the QM region", where
// the truth is that the system does not say. Within a context ('qm in water') it is the part of the
// region that lies there, which can be nothing.
static int _select_region(data_t* dst, data_t arg[], eval_context_t* ctx, bool qm) {
    ASSERT(ctx && ctx->sys);
    (void)arg;

    const md_system_t* sys = ctx->sys;
    const size_t num_atoms = sys->atom.count;

    if (!dst) {
        size_t num_qm = 0;
        for (size_t i = 0; i < num_atoms; ++i) {
            num_qm += (md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_QM) ? 1 : 0;
        }
        if (num_qm == 0) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, qm ?
                "The system has no QM region: it was not loaded from a quantum chemistry calculation" :
                "The system has no QM region, and so no environment of one");
            return -1;
        }
        if (!qm && num_qm == num_atoms) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "The system is QM throughout: it has no environment");
            return -1;
        }
    }

    const md_atom_flags_t want = qm ? MD_ATOM_FLAG_QM : MD_ATOM_FLAG_NONE;

    // The groups: the residues in the context, or the whole system as one when it has none
    md_array(md_urange_t) groups = 0;
    if (sys->residue.count) {
        int* res_indices = get_res_indices_in_context(sys, ctx->mol_ctx, ctx->temp_alloc);
        for (size_t i = 0; i < md_array_size(res_indices); ++i) {
            md_array_push(groups, md_system_residue_atom_range(sys, res_indices[i]), ctx->temp_alloc);
        }
        md_array_free(res_indices, ctx->temp_alloc);
    } else if (!ctx->mol_ctx || md_bitset_test_any_range(*ctx->mol_ctx, 0, (uint32_t)num_atoms)) {
        md_urange_t all = {0, (uint32_t)num_atoms};
        md_array_push(groups, all, ctx->temp_alloc);
    }

    int result = 0;
    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        md_bitset_t* bf = as_bitset(*dst);
        const int cap = type_info_array_len(dst->type);
        if (cap > 0) {
            // One set per group with any such atom, or all of them in one (cap 1)
            md_allocator_i* alloc = data_bitset_alloc(dst, ctx);
            md_bitmask_t b = result_mask(ctx);
            int dst_idx = 0;
            for (size_t g = 0; g < md_array_size(groups); ++g) {
                bool any = false;
                for (uint32_t i = groups[g].beg; i < groups[g].end; ++i) {
                    if ((md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_QM) == want) {
                        md_bitmask_set_bit(&b, i);
                        any = true;
                    }
                }
                if (any && cap > 1) {
                    ASSERT(dst_idx < cap);
                    if (dst_idx < cap) bf[dst_idx++] = md_bitmask_extract(&b, alloc);
                    md_bitmask_clear(&b);
                }
            }
            if (cap == 1) bf[0] = md_bitmask_extract(&b, alloc);
            md_bitmask_free(&b);
        }
    } else {
        int count = 0;
        for (size_t g = 0; g < md_array_size(groups); ++g) {
            for (uint32_t i = groups[g].beg; i < groups[g].end; ++i) {
                if ((md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_QM) == want) {
                    count += 1;
                    break;
                }
            }
        }
        if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
            count = MIN(1, count);
        }
        result = count;
    }

    md_array_free(groups, ctx->temp_alloc);
    return result;
}

static int _qm(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _select_region(dst, arg, ctx, true);
}

static int _environment(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _select_region(dst, arg, ctx, false);
}

static int _backbone(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _select_atoms_with_flags(dst, arg, ctx, MD_ATOM_FLAG_BACKBONE);
}

static int _side(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _select_atoms_with_flags(dst, arg, ctx, MD_ATOM_FLAG_SIDE_CHAIN | MD_ATOM_FLAG_NUCLEOSIDE);
}

static int _sidechain(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _select_atoms_with_flags(dst, arg, ctx, MD_ATOM_FLAG_SIDE_CHAIN);
}

static int _nucleoside(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _select_atoms_with_flags(dst, arg, ctx, MD_ATOM_FLAG_NUCLEOSIDE);
}

static int _nucleobase(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _select_atoms_with_flags(dst, arg, ctx, MD_ATOM_FLAG_NUCLEOBASE);
}

// head(): the atoms of the heads of the lipids (md_lipid_data_t), hydrogens included: all of a lipid but its tails, the
// head group, the phosphate, the glycerol and the linkages of a phospholipid, the hydroxyl of a sterol, both heads of a
// bolalipid. A lipid whose heads and tails are not known has none.
static int _head(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    (void)arg;
    if (dst) {
        ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
        const md_system_t* sys = ctx->sys;
        md_bitmask_t b = result_mask(ctx);
        for (size_t l = 0; l < sys->lipid.count; ++l) {
            if (sys->lipid.head_offset[l + 1] == sys->lipid.head_offset[l]) continue;
            const int32_t* residues = md_lipid_residues(&sys->lipid, l);
            const size_t num_residues = md_lipid_residue_count(&sys->lipid, l);
            for (size_t c = 0; c < num_residues; ++c) {
                const md_urange_t r = md_system_residue_atom_range(sys, (size_t)residues[c]);
                for (uint32_t i = r.beg; i < r.end; ++i) {
                    if (!(md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_LIPID_TAIL)) md_bitmask_set_bit(&b, i);
                }
            }
        }
        restrict_to_context(&b, ctx);
        finish_result(dst, 0, &b, ctx);
    }
    return 0;
}

// tail(): one selection per tail of the lipids (md_lipid_data_t), in their order, with their branches and hydrogens: the
// acyl, alkyl and sphingoid chains, and the ring skeleton of a sterol. Within a context the tails fully within it, as
// ring() takes the rings. Visualized as the spines of the tails.
static int _tail(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    (void)arg;
    const md_system_t* sys = ctx->sys;
    const md_lipid_data_t* lipid = &sys->lipid;

    // The tails taken
    md_array(uint32_t) tails = 0;
    for (size_t t = 0; t < lipid->tail.count; ++t) {
        if (ctx->mol_ctx) {
            const int32_t* atoms = md_lipid_tail_atoms(lipid, t);
            const size_t n = md_lipid_tail_atom_count(lipid, t);
            bool inside = true;
            for (size_t k = 0; k < n && inside; ++k) inside = md_bitset_test(*ctx->mol_ctx, (uint32_t)atoms[k]);
            if (!inside) continue;
        }
        md_array_push(tails, (uint32_t)t, ctx->temp_alloc);
    }
    const int count = (int)md_array_size(tails);

    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        if (!dst->ptr) return 0;
        md_bitset_t* bf_arr = as_bitset(*dst);
        const int64_t cap = type_info_array_len(dst->type);
        if (cap <= 0) return 0;
        // Flattened (cap 1): all of the tails in one. Tails beyond the capacity go into the last set.
        md_allocator_i* alloc = data_bitset_alloc(dst, ctx);
        md_bitmask_t b = result_mask(ctx);
        int64_t slot = 0;
        for (int i = 0; i < count; ++i) {
            const int64_t s = cap == 1 ? 0 : MIN(i, cap - 1);
            if (s != slot) {
                bf_arr[slot] = md_bitmask_extract(&b, alloc);
                md_bitmask_clear(&b);
                slot = s;
            }
            const int32_t* atoms = md_lipid_tail_atoms(lipid, tails[i]);
            const size_t n = md_lipid_tail_atom_count(lipid, tails[i]);
            for (size_t k = 0; k < n; ++k) md_bitmask_set_bit(&b, (uint32_t)atoms[k]);
        }
        if (count > 0) bf_arr[slot] = md_bitmask_extract(&b, alloc);
        md_bitmask_free(&b);
        return 0;
    }

    if (ctx->vis) {
        const vec3_t* xyz = ctx->cur_state->xyz;
        for (int i = 0; i < count; ++i) {
            const int32_t* spine = md_lipid_tail_spine(lipid, tails[i]);
            const size_t n = md_lipid_tail_spine_size(lipid, tails[i]);
            for (size_t k = 1; k < n; ++k) {
                const vec3_t a = xyz[spine[k - 1]];
                const vec3_t c = xyz[spine[k]];
                push_line(vertex(vec3_set(a.x, a.y, a.z), COLOR_WHITE), vertex(vec3_set(c.x, c.y, c.z), COLOR_WHITE), ctx->vis);
            }
        }
        return 0;
    }

    if (ctx->eval_flags & EVAL_FLAG_FLATTEN) return MIN(1, count);
    return count;
}

// residue(index), component(index): residues by their 1-based position, in a context among the residues with an atom
// in it (get_res_indices_in_context), each cut to the context: residue(2) in chain(1) is the second residue of the
// chain. A position before the first is an error; one past the last is an error without a context, and selects
// nothing within one (a chain need not have a tenth residue), as ring(index) does.
static int _res(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    const md_system_t* sys = ctx->sys;

    // The residues counted: those of the context, or all of the system (no list)
    md_array(int32_t) residues = ctx->mol_ctx ? get_res_indices_in_context(sys, ctx->mol_ctx, ctx->temp_alloc) : NULL;
    const int64_t num_residues = ctx->mol_ctx ? (int64_t)md_array_size(residues) : (int64_t)sys->residue.count;

    md_array(int32_t) pos = 0;
    const int count = select_positions(dst ? &pos : NULL, as_irange_arr(arg[0]), element_count(arg[0]), num_residues, "residue", !dst, ctx);
    if (count < 0) return count;

    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        if (!dst->ptr || type_info_array_len(dst->type) == 0) return 0;
        set_writer_t w = set_writer_begin(dst, ctx, true);
        for (int i = 0; i < count; ++i) {
            const md_urange_t atom_range = md_system_residue_atom_range(sys, residues ? (size_t)residues[pos[i]] : (size_t)pos[i]);
            set_writer_push_range(&w, atom_range.beg, atom_range.end);
        }
        set_writer_end(&w);
        return 0;
    }
    return (ctx->eval_flags & EVAL_FLAG_FLATTEN) ? MIN(1, count) : count;
}

static int _fill_res(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    const size_t num_bf = element_count(arg[0]);
    const md_bitset_t* src_bf = as_bitset(arg[0]);

    md_bitset_t tmp_bf = {0};
    // None (an empty array, 'residue(x)' of an empty x) is flattened as well: src_bf would point at nothing
    if (num_bf != 1) {
        tmp_bf = _internal_flatten_bf(src_bf, num_bf, ctx->temp_alloc);
        src_bf = &tmp_bf;
    }

    int result = STATIC_VALIDATION_ERROR;
    if (ctx->sys && ctx->sys->residue.atom_offset) {
        if (dst) {
            const int cap = type_info_array_len(dst->type);
            if (cap == 0) return 0;
            ASSERT(cap > 0);

            set_writer_t w = set_writer_begin(dst, ctx, false);
            for (size_t i = 0; i < ctx->sys->residue.count; ++i) {
                const md_urange_t range = md_system_residue_atom_range(ctx->sys, i);
                if (md_bitset_test_any_range(*src_bf, range.beg, range.end)) {
                    set_writer_push_range(&w, range.beg, range.end);
                }
            }
            set_writer_end(&w);

            result = 0;
        } else {
            int count = 0;
            
            if (ctx->backchannel) {
                if (ctx->arg_flags[0] & FLAG_DYNAMIC) {
                    ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
                }
            }

            for (size_t i = 0; i < ctx->sys->residue.count; ++i) {
                const md_urange_t range = md_system_residue_atom_range(ctx->sys, i);
                size_t popcount = md_bitset_count_range(*src_bf, range.beg, range.end);
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


    return result;
}

// chain(sel): the chains the selection has an atom of, whole, one set each. An atom of no chain (a molecule the data
// gives no chain id, which is not a polymer) adds none.
static int _fill_chain(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    ASSERT(ctx && ctx->sys);
    const md_system_t* sys = ctx->sys;

    if (md_system_chain_count(sys) == 0) {
        if (!dst) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The data gives no chain ids, so the system has no chains: a strand of a polymer is polymer(i), or part(sel)");
        return -1;
    }

    const int64_t num_bf = element_count(arg[0]);
    const md_bitset_t* src_bf = as_bitset(arg[0]);

    md_bitset_t tmp_bf = {0};
    // None (an empty array, 'chain(x)' of an empty x) is flattened as well: src_bf would point at nothing
    if (num_bf != 1) {
        tmp_bf = _internal_flatten_bf(src_bf, (size_t)num_bf, ctx->temp_alloc);
        src_bf = &tmp_bf;
    }

    md_array(int32_t) chains = 0;
    for (size_t c = 0; c < md_system_chain_count(sys); ++c) {
        if (chain_in_context(sys, c, src_bf)) md_array_push(chains, (int32_t)c, ctx->temp_alloc);
    }
    const int count = (int)md_array_size(chains);

    if (dst) {
        if (!dst->ptr || type_info_array_len(dst->type) == 0) return 0;
        set_writer_t w = set_writer_begin(dst, ctx, false);
        for (int i = 0; i < count; ++i) {
            md_bitmask_t b = result_mask(ctx);
            chain_set_bits(&b, sys, (size_t)chains[i]);
            set_writer_push(&w, md_bitmask_view(&b));
            md_bitmask_free(&b);
        }
        set_writer_end(&w);
        return 0;
    }

    if (ctx->backchannel && (ctx->arg_flags[0] & FLAG_DYNAMIC)) {
        ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
    }
    return (ctx->eval_flags & EVAL_FLAG_FLATTEN) ? MIN(1, count) : count;
}

static int _res_name(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));

    int result = 0;

    if (!ctx->sys->residue.name) {
        EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The molecule does not contain any residue names");
        return -1;
    }

    const size_t num_queries = element_count(arg[0]);
    const str_t* queries = as_string_arr(arg[0]);
    
    // Here we only pick the residues which are represented within the context (by having bits set)
    int32_t* res_indices = get_res_indices_in_context(ctx->sys, ctx->mol_ctx, ctx->temp_alloc);
    const size_t res_count = md_array_size(res_indices);

    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        if (dst->ptr) {
            const int64_t cap = type_info_array_len(dst->type);
            if (cap == 0) return 0;
            set_writer_t w = set_writer_begin(dst, ctx, true);

            for (size_t i = 0; i < res_count; ++i) {
                const int64_t res_idx = res_indices[i];
                for (size_t j = 0; j < num_queries; ++j) {
                    if (match_query(queries[j], LBL_TO_STR(ctx->sys->residue.name[res_idx]))) {
                        const md_urange_t atom_range = md_system_residue_atom_range(ctx->sys, res_idx);
                        set_writer_push_range(&w, atom_range.beg, atom_range.end);

                        break;
                    }
                }
            }
            set_writer_end(&w);
        }
    }
    else {
        int count = 0;
        for (size_t j = 0; j < num_queries; ++j) {
            if (!validate_query(queries[j])) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The string '%.*s' is not a valid query", queries[j].len, queries[j].ptr);
                return -1;
            }
            bool match = false;
            for (size_t i = 0; i < res_count; ++i) {
                const int64_t res_idx = res_indices[i];
                if (match_query(queries[j], LBL_TO_STR(ctx->sys->residue.name[res_idx]))) {
                    count += 1;
                    match = true;
                }
            }
            if (!match) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The string '%.*s' did not match any residue within the context", queries[j].len, queries[j].ptr);
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

static int _res_seq_id(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));

    int result = 0;

    if (!ctx->sys->residue.seq_id) {
        EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The molecule does not contain any residue ids");
        return -1;
    }

    const size_t num_rid = element_count(arg[0]);
    const irange_t*  rid = arg[0].ptr;

    // Here we only pick the residues which are represented within the context (by having bits set)
    int32_t* res_indices = get_res_indices_in_context(ctx->sys, ctx->mol_ctx, ctx->temp_alloc);

    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        if (dst->ptr) {
            const int64_t cap = type_info_array_len(dst->type);
            if (cap == 0) return 0;
            set_writer_t w = set_writer_begin(dst, ctx, true);

            for (size_t i = 0; i < md_array_size(res_indices); ++i) {
                const int64_t res_idx = res_indices[i];
                for (size_t j = 0; j < num_rid; ++j) {
                    const int64_t step = rid[j].step > 0 ? rid[j].step : 1;
                    int seq_id = (int)ctx->sys->residue.seq_id[res_idx];
                    if (idx_in_range(seq_id, rid[j])) {
                        int64_t diff = (int64_t)seq_id - (int64_t)rid[j].beg;
                        if ((diff % step) != 0) continue;
                        const md_urange_t atom_range = md_system_residue_atom_range(ctx->sys, res_idx);
                        set_writer_push_range(&w, atom_range.beg, atom_range.end);

                        break;
                    }
                }
            }
            set_writer_end(&w);
        }
    }
    else {
        int count = 0;
        for (size_t j = 0; j < num_rid; ++j) {
            bool match = false;
            const int64_t step = rid[j].step > 0 ? rid[j].step : 1;
            for (size_t i = 0; i < md_array_size(res_indices); ++i) {
                const int64_t res_idx = res_indices[i];
                int seq = (int)ctx->sys->residue.seq_id[res_idx];
                if (idx_in_range(seq, rid[j])) {
                    int64_t diff = (int64_t)seq - (int64_t)rid[j].beg;
                    if ((diff % step) == 0) {
                        count += 1;
                        match = true;
                    }
                }
            }
            if (!match) {
                // @TODO: Should this just be a soft warning instead?
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "No matching seq id was found within the range (%i:%i)", rid[j].beg, rid[j].end);
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
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));

    if (md_system_chain_count(ctx->sys) == 0) {
        if (!dst) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The data gives no chain ids, so the system has no chains: a strand of a polymer is polymer(i)");
        return -1;
    }

    // Those with an atom in the context, or all of the system
    md_array(int32_t) indices = get_chain_indices_in_context(ctx->sys, ctx->mol_ctx, ctx->temp_alloc);
    md_array(int32_t) pos = 0;
    const int count = select_positions(dst ? &pos : NULL, as_irange_arr(arg[0]), element_count(arg[0]), (int64_t)md_array_size(indices), "chain", !dst, ctx);
    if (count < 0) return count;

    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        if (!dst->ptr || type_info_array_len(dst->type) == 0) return 0;
        set_writer_t w = set_writer_begin(dst, ctx, true);
        for (int i = 0; i < count; ++i) {
            md_bitmask_t b = result_mask(ctx);
            chain_set_bits(&b, ctx->sys, (size_t)indices[pos[i]]);
            set_writer_push(&w, md_bitmask_view(&b));
            md_bitmask_free(&b);
        }
        set_writer_end(&w);
        return 0;
    }
    return (ctx->eval_flags & EVAL_FLAG_FLATTEN) ? MIN(1, count) : count;
}

// chain("A"): every atom which carries the chain id, as one set (see md_chain_data_t): a chain id is a label of the
// data, as an atom name is, rather than a structure counted one by one. Several ids, or a pattern matching several, are
// all of their atoms together. Within a context, its atoms which carry the id. An id which no chain in the context has
// is a compile error, and so is any id where the data gives none: the system then has no chains, and a strand of a
// polymer is found as one (polymer(i)).
static int _chain_id(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_STRING_ARR));
    const md_system_t* sys = ctx->sys;

    const size_t num_str = element_count(arg[0]);
    const str_t* str = as_string_arr(arg[0]);
    md_array(int32_t) chains = get_chain_indices_in_context(sys, ctx->mol_ctx, ctx->temp_alloc);

    if (dst) {
        ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
        md_bitmask_t b = result_mask(ctx);
        for (size_t i = 0; i < md_array_size(chains); ++i) {
            const str_t id = md_system_chain_id(sys, (size_t)chains[i]);
            if (str_empty(id)) continue;
            for (size_t j = 0; j < num_str; ++j) {
                if (match_query(str[j], id)) {
                    chain_set_bits(&b, sys, (size_t)chains[i]);
                    break;
                }
            }
        }
        restrict_to_context(&b, ctx);
        finish_result(dst, 0, &b, ctx);
        return 0;
    }

    if (md_system_chain_count(sys) == 0) {
        EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The data gives no chain ids, so the system has no chains: a strand of a polymer is polymer(i)");
        return -1;
    }
    for (size_t j = 0; j < num_str; ++j) {
        if (!validate_query(str[j])) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The string '" STR_FMT "' is not a valid query", STR_ARG(str[j]));
            return -1;
        }
        bool match = false;
        for (size_t i = 0; i < md_array_size(chains) && !match; ++i) {
            const str_t id = md_system_chain_id(sys, (size_t)chains[i]);
            match = !str_empty(id) && match_query(str[j], id);
        }
        if (!match) {
            // The ids there are, when they are few enough to list
            md_str_writer_t sb = md_str_writer_create(ctx->temp_alloc);
            size_t listed = 0;
            for (size_t i = 0; i < md_array_size(chains) && listed <= 16; ++i) {
                const str_t id = md_system_chain_id(sys, (size_t)chains[i]);
                if (str_empty(id)) continue;
                md_str_writer_fmt(&sb, "%s'" STR_FMT "'", listed ? ", " : "", STR_ARG(id));
                listed += 1;
            }
            const char* where = ctx->mol_ctx ? " within the context" : "";
            if (listed == 0 || listed > 16) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The string '" STR_FMT "' did not match any chain id%s", STR_ARG(str[j]), where);
            } else {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The string '" STR_FMT "' did not match any chain id%s, they are " STR_FMT, STR_ARG(str[j]), where, STR_ARG(md_str_writer_to_str(&sb)));
            }
            return -1;
        }
    }
    return 0;
}

// ### MOLECULES ###
// The molecules of the system (md_molecule_data_t): its atoms joined by bonds, within residues and along polymer
// chains. In the order of their first atoms; their atoms need not be one range.

static inline void molecule_set_bits(md_bitmask_t* b, const md_system_t* sys, size_t mol_idx) {
    const md_molecule_t mol = md_system_molecule(sys, mol_idx);
    for (size_t k = 0; k < mol.count; ++k) md_bitmask_set_bit(b, (uint32_t)mol.atom_idx[k]);
}

static inline bool molecule_in_context(const md_system_t* sys, size_t mol_idx, const md_bitset_t* bitfield) {
    if (!bitfield) return true;
    const md_molecule_t mol = md_system_molecule(sys, mol_idx);
    for (size_t k = 0; k < mol.count; ++k) {
        if (md_bitset_test(*bitfield, (uint32_t)mol.atom_idx[k])) return true;
    }
    return false;
}

static md_array(int32_t) get_molecule_indices_in_context(const md_system_t* sys, const md_bitset_t* bitfield, md_allocator_i* alloc) {
    md_array(int32_t) arr = 0;
    const size_t num = md_system_molecule_count(sys);
    for (size_t m = 0; m < num; ++m) {
        if (molecule_in_context(sys, m, bitfield)) md_array_push(arr, (int32_t)m, alloc);
    }
    return arr;
}

// molecule(i): the molecules by their 1-based position, inside a context among those with an atom in it, each cut to it
static int _molecule_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));

    md_array(int32_t) indices = get_molecule_indices_in_context(ctx->sys, ctx->mol_ctx, ctx->temp_alloc);
    md_array(int32_t) pos = 0;
    const int count = select_positions(dst ? &pos : NULL, as_irange_arr(arg[0]), element_count(arg[0]), (int64_t)md_array_size(indices), "molecule", !dst, ctx);
    if (count < 0) return count;

    if (dst) {
        if (!dst->ptr || type_info_array_len(dst->type) == 0) return 0;
        set_writer_t w = set_writer_begin(dst, ctx, true);
        for (int i = 0; i < count; ++i) {
            md_bitmask_t b = result_mask(ctx);
            molecule_set_bits(&b, ctx->sys, (size_t)indices[pos[i]]);
            set_writer_push(&w, md_bitmask_view(&b));
            md_bitmask_free(&b);
        }
        set_writer_end(&w);
        return 0;
    }
    return (ctx->eval_flags & EVAL_FLAG_FLATTEN) ? MIN(1, count) : count;
}

// molecule(sel): the molecules the selection has an atom of, whole, one set each
static int _fill_molecule(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    ASSERT(ctx && ctx->sys);
    const md_system_t* sys = ctx->sys;

    const int64_t num_bf = element_count(arg[0]);
    const md_bitset_t* src_bf = as_bitset(arg[0]);
    md_bitset_t tmp_bf = {0};
    if (num_bf != 1) {
        tmp_bf = _internal_flatten_bf(src_bf, (size_t)num_bf, ctx->temp_alloc);
        src_bf = &tmp_bf;
    }

    // The molecule of each selected atom, by its slot: the molecules are in the order of their slots
    md_array(int32_t) mols = 0;
    const size_t num_mols = md_system_molecule_count(sys);
    if (num_mols && sys->molecule.atom_slot) {
        uint8_t* taken = (uint8_t*)md_alloc(ctx->temp_alloc, num_mols);
        MEMSET(taken, 0, num_mols);
        md_bitset_iter_t it = md_bitset_iter(*src_bf);
        while (md_bitset_iter_next(&it)) {
            if (it.idx >= sys->atom.count) break;
            const md_molecule_idx_t mol = md_system_molecule_find_by_atom_idx(sys, it.idx);
            if (mol >= 0) taken[mol] = 1;
        }
        for (size_t m = 0; m < num_mols; ++m) {
            if (taken[m]) md_array_push(mols, (int32_t)m, ctx->temp_alloc);
        }
    }
    const int count = (int)md_array_size(mols);

    if (dst) {
        if (!dst->ptr || type_info_array_len(dst->type) == 0) return 0;
        set_writer_t w = set_writer_begin(dst, ctx, false);
        for (int i = 0; i < count; ++i) {
            md_bitmask_t b = result_mask(ctx);
            molecule_set_bits(&b, sys, (size_t)mols[i]);
            set_writer_push(&w, md_bitmask_view(&b));
            md_bitmask_free(&b);
        }
        set_writer_end(&w);
        return 0;
    }

    if (ctx->backchannel && (ctx->arg_flags[0] & FLAG_DYNAMIC)) {
        ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
    }
    return (ctx->eval_flags & EVAL_FLAG_FLATTEN) ? MIN(1, count) : count;
}

// ### PARTS ###
// The parts of the system (md_part_data_t): a polymer strand, a glycan tree, a ligand, an ion, a water. Each is a run of
// consecutive residues, so one range of atoms.

// The parts with an atom in the context (all of them without one), in their order
static md_array(int32_t) get_part_indices_in_context(const md_system_t* sys, const md_bitset_t* bitfield, md_allocator_i* alloc) {
    md_array(int32_t) arr = 0;
    for (size_t p = 0; p < md_system_part_count(sys); ++p) {
        const md_urange_t r = md_system_part_atom_range(sys, p);
        if (!bitfield || md_bitset_count_range(*bitfield, r.beg, r.end)) md_array_push(arr, (int32_t)p, alloc);
    }
    return arr;
}

// part(i): the parts by their 1-based position, inside a context among those with an atom in it, each cut to it
static int _part_irng(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));

    md_array(int32_t) indices = get_part_indices_in_context(ctx->sys, ctx->mol_ctx, ctx->temp_alloc);
    md_array(int32_t) pos = 0;
    const int count = select_positions(dst ? &pos : NULL, as_irange_arr(arg[0]), element_count(arg[0]), (int64_t)md_array_size(indices), "part", !dst, ctx);
    if (count < 0) return count;

    if (dst) {
        if (!dst->ptr || type_info_array_len(dst->type) == 0) return 0;
        set_writer_t w = set_writer_begin(dst, ctx, true);
        for (int i = 0; i < count; ++i) {
            const md_urange_t r = md_system_part_atom_range(ctx->sys, (size_t)indices[pos[i]]);
            set_writer_push_range(&w, r.beg, r.end);
        }
        set_writer_end(&w);
        return 0;
    }
    return (ctx->eval_flags & EVAL_FLAG_FLATTEN) ? MIN(1, count) : count;
}

// part(sel): the parts the selection has an atom of, whole, one set each
static int _fill_part(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    ASSERT(ctx && ctx->sys);
    const md_system_t* sys = ctx->sys;

    const int64_t num_bf = element_count(arg[0]);
    const md_bitset_t* src_bf = as_bitset(arg[0]);
    md_bitset_t tmp_bf = {0};
    if (num_bf != 1) {
        tmp_bf = _internal_flatten_bf(src_bf, (size_t)num_bf, ctx->temp_alloc);
        src_bf = &tmp_bf;
    }
    md_array(int32_t) parts = get_part_indices_in_context(sys, src_bf, ctx->temp_alloc);
    const int count = (int)md_array_size(parts);

    if (dst) {
        if (!dst->ptr || type_info_array_len(dst->type) == 0) return 0;
        set_writer_t w = set_writer_begin(dst, ctx, false);
        for (int i = 0; i < count; ++i) {
            const md_urange_t r = md_system_part_atom_range(sys, (size_t)parts[i]);
            set_writer_push_range(&w, r.beg, r.end);
        }
        set_writer_end(&w);
        return 0;
    }

    if (ctx->backchannel && (ctx->arg_flags[0] & FLAG_DYNAMIC)) {
        ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
    }
    return (ctx->eval_flags & EVAL_FLAG_FLATTEN) ? MIN(1, count) : count;
}

// glycan(): one selection per glycan, a part of saccharides only: each of its residues has the ring of a monosaccharide
// (MD_RING_FLAG_SUGAR) and is of no polymer. A glycan tree attached to a protein, an oligosaccharide on its own, a
// single sugar. Within a context the glycans with an atom in it, whole, as lipid() takes them.
static int _glycan(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)arg;
    ASSERT(ctx && ctx->sys);
    const md_system_t* sys = ctx->sys;
    const size_t num_res = sys->residue.count;

    md_array(int32_t) glycans = 0;
    if (num_res && md_system_part_count(sys)) {
        // The residues with a sugar ring
        uint8_t* sugar = (uint8_t*)md_alloc(ctx->temp_alloc, num_res);
        MEMSET(sugar, 0, num_res);
        const size_t num_rings = md_ring_data_count(&sys->ring);
        for (size_t r = 0; r < num_rings; ++r) {
            if (!(md_system_ring_flags(sys, r) & MD_RING_FLAG_SUGAR)) continue;
            const int32_t* atoms = md_ring_atoms(&sys->ring, r);
            const md_residue_idx_t c = md_system_residue_find_by_atom_idx(sys, (size_t)atoms[0]);
            if (c >= 0) sugar[c] = 1;
        }
        for (size_t p = 0; p < md_system_part_count(sys); ++p) {
            const md_urange_t residues = md_system_part_residue_range(sys, p);
            bool all = residues.beg < residues.end;
            for (uint32_t c = residues.beg; c < residues.end && all; ++c) {
                const md_residue_kind_t kind = md_system_residue_kind(sys, c);
                all = sugar[c] && kind != MD_RESIDUE_KIND_AMINO_ACID && kind != MD_RESIDUE_KIND_NUCLEOTIDE;
            }
            if (!all) continue;
            const md_urange_t r = md_system_part_atom_range(sys, p);
            if (ctx->mol_ctx && !md_bitset_count_range(*ctx->mol_ctx, r.beg, r.end)) continue;
            md_array_push(glycans, (int32_t)p, ctx->temp_alloc);
        }
    }
    const int count = (int)md_array_size(glycans);

    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET_ARR));
        if (!dst->ptr || type_info_array_len(dst->type) == 0) return 0;
        set_writer_t w = set_writer_begin(dst, ctx, false);
        for (int i = 0; i < count; ++i) {
            const md_urange_t r = md_system_part_atom_range(sys, (size_t)glycans[i]);
            set_writer_push_range(&w, r.beg, r.end);
        }
        set_writer_end(&w);
        return 0;
    }
    return (ctx->eval_flags & EVAL_FLAG_FLATTEN) ? MIN(1, count) : count;
}

// Property Compute

static inline void draw_distance_text(const vec3_t c_a, const vec3_t c_b, const float value, md_script_vis_t* vis) {
    vec3_t c = vec3_mul1(vec3_add(c_a, c_b), 0.5f);
    // Coordinates are Angstrom, so that is what the separation between two of them is.
    push_quantity_text(c, value, md_unit_angstrom(), vis);
}

// Draw single distance between two points
static inline void draw_distance(const vec3_t c_a, const vec3_t c_b, const float value, md_script_vis_t* vis, uint32_t vis_flags) {
    if (vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
        md_script_vis_vertex_t va = { c_a, COLOR_WHITE };
        md_script_vis_vertex_t vb = { c_b, COLOR_WHITE };
        push_point(va, vis);
        push_point(vb, vis);
        push_line(va, vb, vis);
    }
    if (vis_flags & MD_SCRIPT_VISUALIZE_TEXT) {
        draw_distance_text(c_a, c_b, value, vis);
    }
}

static void draw_distances(const vec3_t* a, size_t num_a, const vec3_t* b, size_t num_b, const float* values, md_script_vis_t* vis, uint32_t vis_flags) {
    for (size_t i = 0; i < num_a; ++i) {
        md_script_vis_vertex_t va = {a[i], COLOR_WHITE};
        if (vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
            push_point(va, vis);
        }
        for (size_t j = 0; j < num_b; ++j) {
            md_script_vis_vertex_t vb = {b[j], COLOR_WHITE};
            if (vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                push_point(vb, vis);
                push_line(va, vb, vis);
            }
            if (vis_flags & MD_SCRIPT_VISUALIZE_TEXT) {
                size_t idx = i * num_b + j;
                draw_distance_text(a[i], b[j], values[idx], vis);
            }
        }
    }
}

static int _distance(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        vec3_t a = coordinate_extract_com(arg[0], ctx);
        vec3_t b = coordinate_extract_com(arg[1], ctx);

        vec4_t b4 = vec4_from_vec3(b, 0);
        md_util_deperiodize_vec4(&b4, 1, a, &ctx->cur_state->unitcell);
        b = vec3_from_vec4(b4);
        float dist = vec3_distance(a, b);

        if (dst) {
            ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
            float* ptr = dst->ptr;
            *ptr = dist;
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
            coordinate_visualize(arg[1], ctx);
            draw_distance(a, b, dist, ctx->vis, ctx->vis_flags);
        }
    }
    else {
        int res_a = coordinate_validate(arg[0], 0, ctx);
        int res_b = coordinate_validate(arg[1], 1, ctx);
        if (res_a < 0) return res_a;
        if (res_b < 0) return res_b;
        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_angstrom();
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
        return 1;
    }

    return 0;
}

// distance_min and distance_max: one value per element of a, the shortest (longest) distance from the positions of that
// element to any position of b. An element of a is a selection, an atom index, a range of atoms or a point, and is
// measured by its atoms, not by its centre of mass. b is taken as a whole: an array of selections is merged into one.
static int distance_extent(data_t* dst, data_t arg[], eval_context_t* ctx, bool largest) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));

    const size_t num_elem = element_count(arg[0]);

    if (dst || ctx->vis) {
        md_temp_scope_t temp = md_temp_begin_in(ctx->temp_alloc);

        data_t b = arg[1];
        md_bitset_t b_bf = {0};
        if (b.type.base_type == TYPE_ATOMSET && element_count(b) > 1) {
            b_bf = _internal_flatten_bf(as_bitset(b), element_count(b), ctx->temp_alloc);
            b.ptr  = &b_bf;
            b.size = sizeof(md_bitset_t);
            b.type = type_info_element_type(b.type);
        }
        const vec3_t* b_pos = coordinate_extract(b, ctx);
        const size_t  b_len = md_array_size(b_pos);

        // Total count, not element_count: evaluated in a context the result is laid out as [context][N] (see rmsd)
        float* out = NULL;
        size_t out_len = 0;
        if (dst) {
            ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
            out = as_float_arr(*dst);
            out_len = type_info_total_element_count(dst->type);
            MEMSET(out, 0, out_len * sizeof(float));
        }

        size_t beg = 0;
        size_t end = num_elem;
        md_array(irange_t) subscript_ranges = ctx->subscript_ranges;
        if (ctx->vis && subscript_ranges) {
            // Only the elements of the result which are shown, e.g. a hovered element of its plot
            beg = (size_t)CLAMP(subscript_ranges[0].beg, 0, (int)num_elem);
            end = (size_t)CLAMP(subscript_ranges[0].end, (int)beg, (int)num_elem);
            // Which is a subscript of the result, not of the single elements and of b which are visualized below
            ctx->subscript_ranges = NULL;
        }

        const type_info_t elem_type   = type_info_element_type(arg[0].type);
        const size_t      elem_stride = type_info_element_byte_stride(arg[0].type);

        // The positions of the elements one after another, all measured against b in one pass: each element starts
        // from the pair of the one before, which is what makes it cheap (see md_util_min_distance_groups)
        const size_t num = end - beg;
        md_array(vec3_t) a_pos = 0;
        size_t* a_off = md_temp_alloc_array(temp, size_t, num + 1);
        for (size_t i = beg; i < end; ++i) {
            const data_t elem = {
                .type = elem_type,
                .ptr  = (char*)arg[0].ptr + i * elem_stride,
                .size = elem_stride,
            };
            const vec3_t* pos = coordinate_extract(elem, ctx);
            a_off[i - beg] = md_array_size(a_pos);
            md_array_push_array(a_pos, pos, md_array_size(pos), ctx->temp_alloc);
            if (ctx->vis) {
                coordinate_visualize(elem, ctx);
            }
        }
        a_off[num] = md_array_size(a_pos);

        float*   dist = md_temp_alloc_array(temp, float,   num);
        int64_t* ia   = md_temp_alloc_array(temp, int64_t, num);
        int64_t* ib   = md_temp_alloc_array(temp, int64_t, num);
        if (largest) {
            md_util_max_distance_groups(dist, ia, ib, a_pos, a_off, num, b_pos, b_len, &ctx->cur_state->unitcell);
        } else {
            md_util_min_distance_groups(dist, ia, ib, a_pos, a_off, num, b_pos, b_len, &ctx->cur_state->unitcell);
        }

        for (size_t k = 0; k < num; ++k) {
            // An empty selection - a dynamic one can be empty in some frames - has no distance to anything: 0
            if (ib[k] < 0) continue;

            const size_t i = beg + k;
            if (out && i < out_len) {
                out[i] = dist[k];
            }
            if (ctx->vis) {
                // The end in b drawn in the image it was measured in
                const vec3_t a = a_pos[ia[k]];
                vec4_t b4 = vec4_from_vec3(b_pos[ib[k]], 0);
                md_util_deperiodize_vec4(&b4, 1, a, &ctx->cur_state->unitcell);
                draw_distance(a, vec3_from_vec4(b4), dist[k], ctx->vis, ctx->vis_flags);
            }
        }

        if (ctx->vis) {
            coordinate_visualize(b, ctx);
            ctx->subscript_ranges = subscript_ranges;
        }

        md_temp_end(temp);
    } else {
        // The positions are validated as for distance, but without its length: the length of the result is the number
        // of elements of a, which a selection changing size from frame to frame does not change
        static_backchannel_t* backchannel = ctx->backchannel;
        ctx->backchannel = NULL;
        const int res_a = coordinate_validate(arg[0], 0, ctx);
        const int res_b = coordinate_validate(arg[1], 1, ctx);
        ctx->backchannel = backchannel;
        if (res_a < 0) return res_a;
        if (res_b < 0) return res_b;
        if (ctx->backchannel) {
            // Unless the number of elements of a is what changes, e.g. residue(within(...))
            if (ctx->arg_flags && (ctx->arg_flags[0] & FLAG_DYNAMIC_LENGTH)) {
                ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
            }
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_angstrom();
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
        return (int)num_elem;
    }

    return 0;
}

static int _distance_min(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return distance_extent(dst, arg, ctx, false);
}

static int _distance_max(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return distance_extent(dst, arg, ctx, true);
}

static int _distance_pair(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));
    
    int result = 0;

    if (dst || ctx->vis) {
        md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);

        vec3_t* a_pos = coordinate_extract(arg[0], ctx);
        vec3_t* b_pos = coordinate_extract(arg[1], ctx);
        const size_t a_len = md_array_size(a_pos);
        const size_t b_len = md_array_size(b_pos);

        float* dst_arr = 0;

        if (dst) {
            ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
            // As for coord_x, dst can be of another length than there are pairs: then they are computed aside, and
            // what fits is written with the rest zeroed
            dst_arr = (dst->size == a_len * b_len * sizeof(float)) ? as_float_arr(*dst) : md_temp_alloc_array(temp, float, a_len * b_len);
        } else {
            dst_arr = md_temp_alloc_array(temp, float, a_len * b_len);
        }
        md_util_distance_array(dst_arr, a_pos, a_len, b_pos, b_len, &ctx->cur_state->unitcell);
        if (dst && dst_arr != as_float_arr(*dst)) {
            const size_t n = MIN(dst->size / sizeof(float), a_len * b_len);
            MEMCPY(dst->ptr, dst_arr, n * sizeof(float));
            MEMSET((char*)dst->ptr + n * sizeof(float), 0, dst->size - n * sizeof(float));
        }
        if (ctx->vis) {
            if (ctx->subscript_ranges) {
                // @NOTE(Robin): This is not trivial...
                // We want to support visualization of subranges, which can be imposed by for example array subscripts or hovering a line plot corresponding to one subelement.
                // In such case we have this encoded in ctx->subscript_ranges, which gives us the range of elements to visualize.
                // This means we have to do a bit of extra work to figure out what elements to visualize, since they need to map 1:1 to the elements in the result array (dst)
                int32_t* a_idx = coordinate_extract_indices(arg[0], ctx);
                int32_t* b_idx = coordinate_extract_indices(arg[1], ctx);

                for (size_t ri = 0; ri < md_array_size(ctx->subscript_ranges); ++ri) {
                    irange_t range = ctx->subscript_ranges[ri];

                    for (int i = range.beg; i < range.end && i < (int)(a_len * b_len); ++i) {
                        // i represents the linear index which we would like to visualize
                        // But we need to map this back to the corresponding indices in a and b
                        // As the distances are essentially calculated and stored like dist[a][b]
                        int a = i / (int)b_len;
                        int b = i % (int)b_len;

                        if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
                            if (a_idx) {
                                visualize_atom_index(a_idx[a], ctx);
                            } else if (arg[0].type.base_type == TYPE_ATOMSET && element_count(arg[0]) > 1) {
                                if (a < (int)element_count(arg[0])) {
                                    const md_bitset_t* bf_arr = as_bitset(arg[0]);
                                    visualize_atom_mask(&bf_arr[a], ctx);
                                }
                            }

                            if (b_idx) {
                                visualize_atom_index(b_idx[b], ctx);
                            } else if (arg[1].type.base_type == TYPE_ATOMSET && element_count(arg[1]) > 1) {
                                if (b < (int)element_count(arg[1])) {
                                    const md_bitset_t* bf_arr = as_bitset(arg[1]);
                                    visualize_atom_mask(&bf_arr[b], ctx);
                                }
                            }
                        }

                        draw_distance(a_pos[a], b_pos[b], dst_arr[i], ctx->vis, ctx->vis_flags);
                    }
                }
            } else {
                coordinate_visualize(arg[0], ctx);
                coordinate_visualize(arg[1], ctx);
                draw_distances(a_pos, a_len, b_pos, b_len, dst_arr, ctx->vis, ctx->vis_flags);
            }
        }
        md_temp_end(temp);
    } else {
        int res0 = coordinate_validate(arg[0], 0, ctx);
        int res1 = coordinate_validate(arg[1], 1, ctx);
        if (res0 < 0) return res0;
        if (res1 < 0) return res1;
        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_angstrom();
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
        size_t count = (size_t)res0 * (size_t)res1;
        if (count > 1000000) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "The size produced by the operation is %zu, which exceeds the upper limit of 1'000'000", count);
            return STATIC_VALIDATION_ERROR;
        }
        result = (int)count;
    }

    return result;
}

static void draw_angle_arc(vec3_t c, vec3_t v, vec3_t axis, float angle, uint32_t color, md_script_vis_t* vis) {
    ASSERT(vis);
    const float len = 0.5f;
    /*if (angle < 0.0f) {
        axis = vec3_mul1(axis, -1.0f);
        angle = -angle;
    }*/
    const size_t num_segments = MAX(1, (size_t)(fabs(angle) / DEG_TO_RAD(10)));
    const float scl = 1.0f / (float)num_segments;
    md_script_vis_vertex_t vb = vertex(c, color);
    md_script_vis_vertex_t vp = vertex(vec3_add(c, vec3_mul1(v, len)), color);
    push_line(vb, vp, vis);
    for (size_t i = 1; i <= num_segments; ++i) {
        const float t = i * scl;
        quat_t q = quat_axis_angle(axis, angle * t);
        vec3_t n = quat_mul_vec3(q, v);
        md_script_vis_vertex_t vn = vertex(vec3_add(c, vec3_mul1(n, len)), color);
        push_triangle(vb, vp, vn, vis);
        push_line(vp, vn, vis);
        vp = vn;
    }
    push_line(vb, vp, vis);

    // Add text. The angle is carried in radians with its unit rather than converted here, so a
    // viewer showing angles in something else is not fighting a conversion baked into the string.
    quat_t q = quat_axis_angle(axis, angle * 0.5f);
    vec3_t n = quat_mul_vec3(q, v);
    c = vec3_add(c, vec3_mul1(n, len + 0.1f));
    push_quantity_text(c, angle, md_unit_radian(), vis);
}

static int _angle(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_COORDINATE_ARR));
    
    if (dst || ctx->vis) {
        const vec3_t b = coordinate_extract_com(arg[1], ctx);

        // Both arms under the minimum image convention, as distance and dihedral do. Taken from the
        // raw coordinates, an angle whose atoms sit on either side of a cell boundary came out as the
        // angle to an arm a whole cell long.
        vec3_t dx[2] = {
            vec3_sub(coordinate_extract_com(arg[0], ctx), b),
            vec3_sub(coordinate_extract_com(arg[2], ctx), b),
        };
        md_util_min_image_vec3(dx, ARRAY_SIZE(dx), &ctx->cur_state->unitcell);

        // The ends placed in the images the arms were taken in, so the drawn angle is the measured one
        const vec3_t a = vec3_add(b, dx[0]);
        const vec3_t c = vec3_add(b, dx[1]);
        const vec3_t v0 = vec3_normalize(dx[0]);
        const vec3_t v1 = vec3_normalize(dx[1]);
        const float  cos_angle = CLAMP(vec3_dot(v0, v1), -1.0f, 1.0f);

        if (dst) {
            ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
            as_float(*dst) = acosf(cos_angle);
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
                const vec3_t axis = vec3_normalize(vec3_cross(v0, v1));
                const float angle = acosf(cos_angle);

                draw_angle_arc(b, v0, axis, angle, COLOR_WHITE, ctx->vis);

                /*
                const size_t num_segments = 10;
                const float scl = 1.0 / num_segments;
                md_script_vis_vertex_t vp = vertex(vec3_add(b, vec3_mul1(v0, 0.5f)), COLOR_WHITE);
                push_line(vb, vp, ctx->vis);
                for (size_t i = 1; i <= num_segments; ++i) {
                    const float t = i * scl;
                    quat_t q = quat_axis_angle(axis, angle * t);
                    vec3_t n = quat_mul_vec3(q, v0);
                    md_script_vis_vertex_t v = vertex(vec3_add(b, vec3_mul1(n, 0.5f)), COLOR_WHITE);
                    push_triangle(vb, vp, v, ctx->vis);
                    push_line(vp, v, ctx->vis);
                    vp = v;
                }
                push_line(vb, vp, ctx->vis);
                */
            }
        }
    } else {
        // Validate input
        if (coordinate_validate(arg[0], 0, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (coordinate_validate(arg[1], 1, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (coordinate_validate(arg[2], 2, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_radian();
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
        vec3_t x[4] = {
            coordinate_extract_com(arg[0], ctx),
            coordinate_extract_com(arg[1], ctx),
            coordinate_extract_com(arg[2], ctx),
            coordinate_extract_com(arg[3], ctx),
        };

        vec3_t dx[3] = {
            vec3_sub(x[1], x[0]),
            vec3_sub(x[2], x[1]),
            vec3_sub(x[3], x[2]),
        };
        md_util_min_image_vec3(dx, ARRAY_SIZE(dx), &ctx->cur_state->unitcell);
        float angle = (float)vec3_dihedral_angle(dx[0], dx[1], dx[2]);

        if (dst) {
            ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
            as_float(*dst) = angle;
        }
        if (ctx->vis) {
            coordinate_visualize(arg[0], ctx);
            coordinate_visualize(arg[1], ctx);
            coordinate_visualize(arg[2], ctx);
            coordinate_visualize(arg[3], ctx);

            if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) {
                md_script_vis_vertex_t va = vertex(x[0], COLOR_WHITE);
                md_script_vis_vertex_t vb = vertex(x[1], COLOR_WHITE);
                md_script_vis_vertex_t vc = vertex(x[2], COLOR_WHITE);
                md_script_vis_vertex_t vd = vertex(x[3], COLOR_WHITE);

                push_point(va, ctx->vis);
                push_point(vb, ctx->vis);
                push_point(vc, ctx->vis);
                push_point(vd, ctx->vis);

                push_line(va, vb, ctx->vis);
                push_line(vb, vc, ctx->vis);
                push_line(vc, vd, ctx->vis);

#if 1
                vec3_t b = x[1];
                vec3_t c = x[2];

                vec3_t ba = vec3_sub(x[0], x[1]);
                vec3_t bc = vec3_sub(x[2], x[1]);
                vec3_t cd = vec3_sub(x[3], x[2]);

                vec3_t u  = vec3_normalize(bc); // Axis
                vec3_t vp = vec3_mul1(u, vec3_dot(ba, u));
                vec3_t vo = vec3_sub(ba, vp);
                vec3_t wp = vec3_mul1(u, vec3_dot(cd, u));
                vec3_t wo = vec3_sub(cd, wp);

                // Plane 1
                {
                    vec3_t d = vec3_mul1(vec3_normalize(vo), 1.0f);
                    md_script_vis_vertex_t v[4] = {
                        vb,
                        vc,
                        vertex(vec3_add(c, d), COLOR_WHITE),
                        vertex(vec3_add(b, d), COLOR_WHITE),
                    };

                    push_line(v[0], v[1], ctx->vis);
                    push_line(v[1], v[2], ctx->vis);
                    push_line(v[2], v[3], ctx->vis);
                    push_line(v[3], v[0], ctx->vis);
                }

                // Plane 2
                {
                    vec3_t d = vec3_mul1(vec3_normalize(wo), 1.0f);
                    md_script_vis_vertex_t v[4] = {
                        vc,
                        vb,
                        vertex(vec3_add(b, d), COLOR_WHITE),
                        vertex(vec3_add(c, d), COLOR_WHITE),
                    };

                    push_line(v[0], v[1], ctx->vis);
                    push_line(v[1], v[2], ctx->vis);
                    push_line(v[2], v[3], ctx->vis);
                    push_line(v[3], v[0], ctx->vis);
                }
                // Angle
                {
                    vec3_t p0 = vec3_lerp(b, c, 0.5f);
                    draw_angle_arc(p0, vo, u, angle, COLOR_WHITE, ctx->vis);
                }
#endif
            }
        }
    } else {
        // Validate input
        if (coordinate_validate(arg[0], 0, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (coordinate_validate(arg[1], 1, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (coordinate_validate(arg[2], 2, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (coordinate_validate(arg[3], 3, ctx) < 0) return STATIC_VALIDATION_ERROR;
        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_radian();
            //ctx->backchannel->value_range = (frange_t){-PI, PI};
            // In [-180, 180[, cut at 0, where torsions are least likely to be: at 0 the outer bonds are eclipsed,
            // the highest barrier of an sp3-sp3 torsion (chi, the nucleic acid backbone), the rare cis of a peptide
            // omega, and a phi which the backbone does not reach. At 180 they would be cut through the trans of
            // all of these, and through the beta strands of psi.
            ctx->backchannel->period = period_make(-PI, PI, 0.0);
        }
    }

    return 0;
}

static int _rmsd(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    ASSERT(ctx);

    // One value per atomset of the argument: rmsd(residue(:)) gives the deviation of every residue on its
    // own, each fitted separately. rmsd(union(...)) pools them into a single fit and a single value.
    const size_t num_bf = element_count(arg[0]);

    if (dst) {
        ASSERT(ctx->sys);
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
        // Without masses (outside the evaluation of frames, when it is visualized) the atoms weigh the same, as they
        // do for every other procedure
        if (!ctx->ref_state || !ctx->ref_state->xyz) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "rmsd: the system has no reference coordinates to compare with");
            return -1;
        }

        // Total count, not element_count: evaluated in a context the result is laid out as [context][N],
        // so this call sees float[1][N]
        // The two agree except for an argument whose length varies over frames evaluated inside a context:
        // the context fixes the length of each of its parts at compile time. Write what fits, zero the rest.
        float* out = as_float_arr(*dst);
        const size_t out_len = type_info_total_element_count(dst->type);
        for (size_t i = num_bf; i < out_len; ++i) {
            out[i] = 0.0f;
        }

        const md_bitset_t* bf_arr = as_bitset(arg[0]);

        md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);

        md_bitset_t tmp_bf = {0};

        md_array(vec4_t) xyzw[2] = {0};

        for (size_t i = 0; i < MIN(num_bf, out_len); ++i) {
            const md_bitset_t* bf = &bf_arr[i];
            if (ctx->mol_ctx) {
                tmp_bf = md_bitset_and(*bf, *ctx->mol_ctx, ctx->temp_alloc);
                bf = &tmp_bf;
            }

            // An empty selection is valid - a dynamic one can be empty in some frames - and deviates by nothing
            out[i] = 0.0f;
            const size_t count = md_bitset_count(*bf);
            if (count == 0) continue;

            md_array_resize(xyzw[0], count, ctx->temp_alloc);
            md_array_resize(xyzw[1], count, ctx->temp_alloc);

            extract_xyzw_vec4(xyzw[0], ctx->ref_state->xyz, ctx->atom_mass, bf);
            extract_xyzw_vec4(xyzw[1], ctx->cur_state->xyz, ctx->atom_mass, bf);

            // Each set is placed into mutually consistent images against its own cell, and its centre is
            // the plain weighted mean of the placed points - which is what the fit below subtracts.
            // Coordinates stay absolute: md_util_rmsd_compute_vec4 subtracts the centres itself.
            // Correct while the selection spans less than half a cell, as any per point image choice.
            vec3_t com[2] = {0};
            md_util_deperiodize_self_vec4(xyzw[0], count, &ctx->ref_state->unitcell, &com[0]);
            md_util_deperiodize_self_vec4(xyzw[1], count, &ctx->cur_state->unitcell, &com[1]);

            out[i] = (float)md_util_rmsd_compute_vec4((const vec4_t* const*)xyzw, 0, count, com);
        }

        md_temp_end(temp);
        return 0;
    }

    if (ctx->vis) {
        // The atomset highlighting the atoms involved is all there is to show
        coordinate_visualize(arg[0], ctx);
        return 0;
    }

    // Query / validation. Nothing to validate: an empty atomset is a valid (zero) deviation, and has to be
    // for dynamic selections. The length follows the argument, including when that varies over frames.
    if (ctx->backchannel) {
        if (ctx->arg_flags && (ctx->arg_flags[0] & FLAG_DYNAMIC_LENGTH)) {
            ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
        }
        ctx->backchannel->unit[0] = md_unit_none();
        ctx->backchannel->unit[1] = md_unit_angstrom();
        ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
    }
    return (int)num_bf;
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

    ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_IRANGE_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_INT_ARR));
    (void)ctx;
    // The length is deduced from the argument, so the two agree: the smaller is what can be read and written
    const size_t dst_len = MIN(element_count(*dst), element_count(arg[0]));

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
    
    ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_INT_ARR));
    (void)ctx;
    // The length is deduced from the argument, so the two agree: the smaller is what can be read and written
    const size_t dst_len = MIN(element_count(*dst), element_count(arg[0]));

    float* d = dst->ptr;
    int*   s = arg[0].ptr;

    for (size_t i = 0; i < dst_len; ++i) {
        d[i] = (float)s[i]; 
    }
    return 0;
}

static int _cast_irng_arr_to_frng_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    
    ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FRANGE_ARR));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_IRANGE_ARR));
    (void)ctx;
    // The length is deduced from the argument, so the two agree: the smaller is what can be read and written
    const size_t dst_len = MIN(element_count(*dst), element_count(arg[0]));

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
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);
    const size_t num_idx = element_count(arg[0]);

    // The idea here is that if we have a context, we only make sure to add the indices which are represented within the context.

    if (dst) {
        ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
        md_bitmask_t b = result_mask(ctx);
        for (size_t i = 0; i < num_idx; ++i) {
            // Shift here since we use 1 based indices for atoms
            const int64_t idx = remap_index_to_context(indices[i], ctx_range);
            // This is pre-checked in the static check, against the context known then
            if (!idx_in_range((int32_t)idx, ctx_range) || (ctx->mol_ctx && !md_bitset_test(*ctx->mol_ctx, (uint32_t)idx))) {
                EVAL_LOG_ERROR(ctx, ctx->op_token, "supplied index (%i) is not represented within the supplied context", indices[i]);
                return -1;
            }
            md_bitmask_set_bit(&b, (uint32_t)idx);
        }
        finish_result(dst, 0, &b, ctx);
    } else {
        for (size_t i = 0; i < num_idx; ++i) {
            if (!validate_atom_index_in_context(indices[i], ctx->arg_tokens[0], ctx)) {
                return STATIC_VALIDATION_ERROR;
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
    const irange_t ctx_range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);

    if (dst) {
        ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
        md_bitmask_t b = result_mask(ctx);
        for (size_t i = 0; i < num_ranges; ++i) {
            const irange_t range = clamp_range(remap_range_to_context(ranges[i], ctx_range), ctx_range);
            const int64_t step = ranges[i].step > 0 ? ranges[i].step : 1;
            if (ctx->mol_ctx) {
                md_bitset_iter_t it = md_bitset_iter_range(*ctx->mol_ctx, (uint32_t)MAX(range.beg, 0), (uint32_t)MAX(range.end, 0));
                int64_t k = 0;
                while (md_bitset_iter_next(&it)) {
                    if ((k % step) == 0) {
                        md_bitmask_set_bit(&b, it.idx);
                    }
                    ++k;
                }
            }
            else {
                if (range.beg >= 0 && range.end > range.beg) {
                    if (step == 1) {
                        md_bitmask_set_range(&b, (uint32_t)range.beg, (uint32_t)range.end);
                    } else {
                        for (int64_t j = range.beg; j < range.end; j += step) {
                            md_bitmask_set_bit(&b, (uint32_t)j);
                        }
                    }
                }
            }
        }
        finish_result(dst, 0, &b, ctx);
    } else {
        for (size_t i = 0; i < num_ranges; ++i) {
            if (!validate_atom_range_in_context(ranges[i], ctx->arg_tokens[0], ctx)) {
                return STATIC_VALIDATION_ERROR;
            }
        }
    }

    return result;
}

static int _join_bf_arr(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)ctx;
    if (!dst) return 0;

    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_ATOMSET));
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));

    const md_bitset_t* src_bf_arr = (md_bitset_t*)arg[0].ptr;
    const int64_t count = type_info_array_len(arg[0].type);

    md_bitmask_t b = result_mask(ctx);
    for (int64_t i = 0; i < count; ++i) {
        md_bitmask_or(&b, src_bf_arr[i]);
    }
    finish_result(dst, 0, &b, ctx);

    return 0;
}

// join(): the old name of union(), which says what becomes of an atom in more than one of the selections
static int _join_deprecated(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst && !ctx->vis) {
        EVAL_LOG_WARNING(ctx, ctx->op_token, "join() is deprecated, use union(): every atom of the selections, once");
        return 0;
    }
    return _join_bf_arr(dst, arg, ctx);
}

// Split an atomset into an array of atomsets given a criterion (residue, chain, structure or integer > 0)
static int _split_bf_int(data_t* dst, data_t arg[], eval_context_t* ctx) {
    // Implement integer split for now
	(void)ctx;

    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_INT));

	int in_N = as_int(arg[1]);

    if (dst) {
        if (in_N <= 0) {
            return 0;
        }

		md_bitset_t* dst_bf_arr = as_bitset(*dst);
		const md_bitset_t* src_bf = (md_bitset_t*)arg[0].ptr;
		md_allocator_i* dst_alloc = data_bitset_alloc(dst, ctx);

		// Modify atomset if we have a context
		md_bitset_t tmp_bf = { 0 };
		if (ctx->mol_ctx) {
			tmp_bf = md_bitset_and(*src_bf, *ctx->mol_ctx, ctx->temp_alloc);
			src_bf = &tmp_bf;
        }

        uint32_t N = (uint32_t)in_N;

        uint32_t src_beg_bit = (uint32_t)src_bf->beg;
		uint32_t src_end_bit = (uint32_t)src_bf->end;
		uint32_t src_size = src_end_bit - src_beg_bit;

		// Split atomset into N atomsets
        for (uint32_t i = 0; i < N; ++i) {
			uint32_t range_beg = (src_size * i) / N;
			uint32_t range_end = (i + 1 == N) ? src_size : (src_size * (i + 1)) / N;

			dst_bf_arr[i] = md_bitset_and(md_bitset_range(src_beg_bit + range_beg, src_beg_bit + range_end), *src_bf, dst_alloc);
		}
        return 0;
    } else {
        // Static check, return number of atomsets returned
        if (in_N <= 0) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[1], "Split count must be greater than zero");
            return STATIC_VALIDATION_ERROR;
		}
		return in_N;
    }
}

static int _com(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));

    vec3_t com;
    if (dst || ctx->vis) {
        com = coordinate_extract_com(arg[0], ctx);
        if (dst) {
            ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT3));
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

// ### PLANE ###
// plane(groups): the least squares plane of each group of atoms, a plane record (plane_t) per group, which centroid(),
// normal(), radius(), deviation() and angle() read. plane(positions): one plane through positions, atom indices in the
// order given or float[3] (from coord(), ...). One plane through all of the atoms of the groups is
// plane(union(groups)).
//
// The normal of a plane is oriented, so that it does not flip from frame to frame:
//   - A group which is a ring of the system (md_system_t::ring): by the order of the ring, the normal on the side from
//     which the atoms go around anticlockwise. A sugar ring goes from the ring O to the anomeric C (see MD_RING_FLAG_SUGAR).
//   - Any other group: by two of its atoms chosen in the reference structure, the one farthest from its centroid and
//     the one which spans the largest area with it. The normal is on the side from which the first turns towards the
//     second anticlockwise.
//   - Positions: by their order, as a ring is, the normal on the side from which they go around anticlockwise.

// The ring of the system whose atoms are those of bf, -1 if none. first_ring[k] .. holds the rings of atom k.
typedef struct plane_rings_t {
    uint32_t* off;      // [num_atoms + 1]
    uint32_t* ring;
} plane_rings_t;

static plane_rings_t plane_rings_build(const md_system_t* sys, md_allocator_i* alloc) {
    plane_rings_t pr = {0};
    const size_t num_rings = md_ring_data_count(&sys->ring);
    const size_t N = sys->atom.count;
    if (num_rings == 0 || N == 0) return pr;
    pr.off = md_alloc(alloc, sizeof(uint32_t) * (N + 1));
    MEMSET(pr.off, 0, sizeof(uint32_t) * (N + 1));
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(&sys->ring, r);
        const size_t n = md_ring_size(&sys->ring, r);
        for (size_t k = 0; k < n; ++k) if ((size_t)ring[k] < N) pr.off[ring[k] + 1] += 1;
    }
    for (size_t i = 0; i < N; ++i) pr.off[i + 1] += pr.off[i];
    pr.ring = md_alloc(alloc, sizeof(uint32_t) * MAX(pr.off[N], 1));
    uint32_t* fill = md_alloc(alloc, sizeof(uint32_t) * N);
    MEMCPY(fill, pr.off, sizeof(uint32_t) * N);
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(&sys->ring, r);
        const size_t n = md_ring_size(&sys->ring, r);
        for (size_t k = 0; k < n; ++k) if ((size_t)ring[k] < N) pr.ring[fill[ring[k]]++] = (uint32_t)r;
    }
    return pr;
}

static int64_t plane_rings_find(const plane_rings_t* pr, const md_system_t* sys, const int32_t* idx, size_t count) {
    if (!pr->off || count < MD_RING_MIN_SIZE || count > MD_RING_MAX_SIZE) return -1;
    const int32_t first = idx[0];
    for (uint32_t t = pr->off[first]; t < pr->off[first + 1]; ++t) {
        const uint32_t r = pr->ring[t];
        if (md_ring_size(&sys->ring, r) != count) continue;
        const int32_t* ring = md_ring_atoms(&sys->ring, r);
        bool same = true;
        for (size_t k = 0; k < count && same; ++k) {
            bool found = false;
            for (size_t j = 0; j < count && !found; ++j) found = ring[k] == idx[j];
            same = found;
        }
        if (same) return r;
    }
    return -1;
}

// The points xyzw made mutually consistent across the cell, and their centroid
static vec3_t plane_place(vec4_t* xyzw, size_t count, const md_unitcell_t* cell) {
    vec3_t c = {0};
    md_util_deperiodize_self_vec4(xyzw, count, cell, &c);
    // The centroid of the points as placed
    vec3_t sum = {0};
    for (size_t k = 0; k < count; ++k) sum = vec3_add(sum, vec3_from_vec4(xyzw[k]));
    return vec3_mul1(sum, 1.0f / (float)count);
}

// Positions of the atoms idx in state, made mutually consistent across the cell, and their centroid
static vec3_t plane_positions(vec4_t* xyzw, const int32_t* idx, size_t count, const md_system_state_t* state) {
    for (size_t k = 0; k < count; ++k) {
        const vec3_t x = state->xyz[idx[k]];
        xyzw[k] = vec4_set(x.x, x.y, x.z, 1.0f);
    }
    return plane_place(xyzw, count, &state->unitcell);
}

// The side from which the points go around anticlockwise, in their order (twice the area of a planar polygon)
static vec3_t plane_turn(const vec4_t* xyzw, size_t count, vec3_t c) {
    vec3_t ref = {0};
    for (size_t k = 0; k < count; ++k) {
        const vec3_t a = vec3_sub(vec3_from_vec4(xyzw[k]), c);
        const vec3_t b = vec3_sub(vec3_from_vec4(xyzw[(k + 1) % count]), c);
        ref = vec3_add(ref, vec3_cross(a, b));
    }
    return ref;
}

// The least squares plane of count >= 3 points, placed, with centroid c. The normal is on the side of ref.
static void plane_from_points(plane_t* out, const vec4_t* xyzw, size_t count, vec3_t c, vec3_t ref) {
    const mat3_t M = mat3_covariance_matrix_vec4(xyzw, 0, count, c);
    const mat3_eigen_t E = mat3_eigen(M);

    // Eigenvalues in descending order: the two in the plane, then the normal
    int order[3] = {0, 1, 2};
    for (int i = 0; i < 3; ++i) for (int j = i + 1; j < 3; ++j) {
        if (E.values.elem[order[j]] > E.values.elem[order[i]]) { int t = order[i]; order[i] = order[j]; order[j] = t; }
    }
    vec3_t n    = vec3_normalize(E.vectors.col[order[2]]);
    vec3_t axis = vec3_normalize(E.vectors.col[order[0]]);

    if (vec3_dot(n, ref) < 0.0f) n = vec3_mul1(n, -1.0f);
    // A right handed frame: axis, normal x axis, normal
    if (vec3_dot(vec3_cross(axis, vec3_cross(n, axis)), n) < 0.0f) axis = vec3_mul1(axis, -1.0f);

    float radius = 0.0f, dev = 0.0f;
    for (size_t k = 0; k < count; ++k) {
        const vec3_t d = vec3_sub(vec3_from_vec4(xyzw[k]), c);
        radius += vec3_length(d);
        const float h = vec3_dot(d, n);
        dev += h * h;
    }
    out->count     = (int)count;
    out->centroid  = c;
    out->normal    = n;
    out->axis      = axis;
    out->extent[0] = sqrtf(MAX(E.values.elem[order[0]], 0.0f));
    out->extent[1] = sqrtf(MAX(E.values.elem[order[1]], 0.0f));
    out->radius    = radius / (float)count;
    out->deviation = sqrtf(dev / (float)count);
}

static void plane_fit(plane_t* out, const int32_t* idx, size_t count, int64_t ring, const eval_context_t* ctx, md_allocator_i* temp) {
    MEMSET(out, 0, sizeof(plane_t));
    out->count = (int)count;
    if (count < 3) return;
    const md_system_t* sys = ctx->sys;

    vec4_t* xyzw = md_alloc(temp, sizeof(vec4_t) * count);
    const vec3_t c = plane_positions(xyzw, idx, count, ctx->cur_state);

    // Orientation
    vec3_t ref = {0};
    if (ring >= 0) {
        // Around the ring, in its order
        const int32_t* ring_idx = md_ring_atoms(&sys->ring, (size_t)ring);
        const size_t n_ring = md_ring_size(&sys->ring, (size_t)ring);
        for (size_t k = 0; k < n_ring; ++k) {
            size_t a = 0, b = 0;
            for (size_t j = 0; j < count; ++j) {
                if (idx[j] == ring_idx[k]) a = j;
                if (idx[j] == ring_idx[(k + 1) % n_ring]) b = j;
            }
            ref = vec3_add(ref, vec3_cross(vec3_sub(vec3_from_vec4(xyzw[a]), c), vec3_sub(vec3_from_vec4(xyzw[b]), c)));
        }
    } else {
        // Two anchors chosen in the reference structure, the same in every frame
        const md_system_state_t* rs = (ctx->ref_state && ctx->ref_state->xyz) ? ctx->ref_state : ctx->cur_state;
        vec4_t* r = md_alloc(temp, sizeof(vec4_t) * count);
        const vec3_t rc = plane_positions(r, idx, count, rs);
        size_t a = 0;
        float best = -1.0f;
        for (size_t k = 0; k < count; ++k) {
            const float d = vec3_length_squared(vec3_sub(vec3_from_vec4(r[k]), rc));
            if (d > best + 1.0e-6f) { best = d; a = k; }
        }
        size_t b = a;
        best = 1.0e-8f;
        const vec3_t ra = vec3_sub(vec3_from_vec4(r[a]), rc);
        for (size_t k = 0; k < count; ++k) {
            const float area = vec3_length(vec3_cross(ra, vec3_sub(vec3_from_vec4(r[k]), rc)));
            if (area > best * (1.0f + 1.0e-4f)) { best = area; b = k; }
        }
        if (b != a) {
            ref = vec3_cross(vec3_sub(vec3_from_vec4(xyzw[a]), c), vec3_sub(vec3_from_vec4(xyzw[b]), c));
        }
    }
    plane_from_points(out, xyzw, count, c, ref);
}

static void plane_visualize(const plane_t* p, eval_context_t* ctx) {
    if (!(ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY) || p->count < 3) return;
    const vec3_t u = vec3_mul1(p->axis, p->radius);
    const vec3_t v = vec3_mul1(vec3_cross(p->normal, p->axis), p->radius);
    const int segments = 32;
    md_script_vis_vertex_t prev = vertex(vec3_add(p->centroid, u), COLOR_WHITE);
    for (int i = 1; i <= segments; ++i) {
        const float t = (float)i / (float)segments * 2.0f * PI;
        const vec3_t x = vec3_add(p->centroid, vec3_add(vec3_mul1(u, cosf(t)), vec3_mul1(v, sinf(t))));
        md_script_vis_vertex_t cur = vertex(x, COLOR_WHITE);
        push_line(prev, cur, ctx->vis);
        prev = cur;
    }
    // The normal, as long as the radius, from the centroid
    push_point(vertex(p->centroid, COLOR_WHITE), ctx->vis);
    push_line(vertex(p->centroid, COLOR_WHITE), vertex(vec3_add(p->centroid, vec3_mul1(p->normal, p->radius)), COLOR_WHITE), ctx->vis);
}

static int _plane_groups(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    const md_bitset_t* groups = as_bitset(arg[0]);
    const size_t num = element_count(arg[0]);

    if (!dst && !ctx->vis) {
        // The groups have to be fixed, and have the atoms a plane needs
        if (ctx->arg_flags[0] & FLAG_DYNAMIC) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "plane: the groups cannot depend on the frame (plane(coord(...)) fits the positions of each frame)");
            return STATIC_VALIDATION_ERROR;
        }
        for (size_t g = 0; g < num; ++g) {
            const size_t n = ctx->mol_ctx ? md_bitset_count(md_bitset_and(groups[g], *ctx->mol_ctx, ctx->temp_alloc)) : md_bitset_count(groups[g]);
            if (n < 3) {
                if (num == 1) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "plane: the group has %zu atom%s, a plane needs at least 3", n, n == 1 ? "" : "s");
                else          EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "plane: group %zu of %zu has %zu atom%s, a plane needs at least 3", g + 1, num, n, n == 1 ? "" : "s");
                return STATIC_VALIDATION_ERROR;
            }
        }
        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_none();
        }
        return (int)num;
    }

    if (dst && ctx->backchannel) {
        // Evaluated within the static check, for a procedure asking about its argument (centroid() its length):
        // there are no coordinates there, only the shape
        MEMSET(dst->ptr, 0, element_count(*dst) * sizeof(plane_t));
        return 0;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(ctx->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    const plane_rings_t pr = plane_rings_build(ctx->sys, temp_alloc);

    plane_t* out = dst ? (plane_t*)dst->ptr : NULL;
    const size_t cap = dst ? element_count(*dst) : num;
    size_t beg = 0, end = MIN(num, cap);
    if (!dst && ctx->vis && ctx->subscript_ranges) {
        beg = (size_t)MAX(ctx->subscript_ranges[0].beg, 0);
        end = MIN((size_t)MAX(ctx->subscript_ranges[0].end, 0), end);
    }

    for (size_t g = beg; g < end; ++g) {
        md_bitset_t group = groups[g];
        if (ctx->mol_ctx) {
            group = md_bitset_and(group, *ctx->mol_ctx, temp_alloc);
        }
        const md_bitset_t* bf = &group;
        const size_t count = md_bitset_count(group);
        int32_t* idx = md_alloc(temp_alloc, sizeof(int32_t) * MAX(count, 1));
        md_bitset_iter_extract_indices(idx, count, md_bitset_iter(group));
        plane_t p;
        plane_fit(&p, idx, count, plane_rings_find(&pr, ctx->sys, idx, count), ctx, temp_alloc);
        if (out) out[g] = p;
        if (ctx->vis) {
            visualize_atom_mask(bf, ctx);
            plane_visualize(&p, ctx);
        }
    }

    md_temp_end(temp);
    return 0;
}

static int _plane_points(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));

    if (!dst && !ctx->vis) {
        // Positions only counted when they are evaluated (coord(within(...))) are checked then
        if (!(arg[0].type.base_type == TYPE_FLOAT && is_variable_length(arg[0].type))) {
            const int count = coordinate_validate(arg[0], 0, ctx);
            if (count < 0) return STATIC_VALIDATION_ERROR;
            if (count < 3) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "plane: %d position%s, a plane needs at least 3", count, count == 1 ? "" : "s");
                return STATIC_VALIDATION_ERROR;
            }
        }
        if (ctx->backchannel) {
            // Always one plane, however many positions
            ctx->backchannel->flags &= ~FLAG_DYNAMIC_LENGTH;
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_none();
        }
        return 1;
    }

    if (dst && ctx->backchannel) {
        // Evaluated within the static check, for a procedure asking about its argument: the shape only
        MEMSET(dst->ptr, 0, element_count(*dst) * sizeof(plane_t));
        return 0;
    }

    plane_t p = {0};
    md_array(vec3_t) xyz = coordinate_extract(arg[0], ctx);
    const size_t count = md_array_size(xyz);
    p.count = (int)count;
    if (count >= 3) {
        vec4_t* xyzw = md_alloc(ctx->temp_alloc, sizeof(vec4_t) * count);
        for (size_t k = 0; k < count; ++k) {
            xyzw[k] = vec4_from_vec3(xyz[k], 1.0f);
        }
        const vec3_t c = plane_place(xyzw, count, &ctx->cur_state->unitcell);
        plane_from_points(&p, xyzw, count, c, plane_turn(xyzw, count, c));
        md_free(ctx->temp_alloc, xyzw, sizeof(vec4_t) * count);
    }

    if (dst) {
        *(plane_t*)dst->ptr = p;
    }
    if (ctx->vis) {
        coordinate_visualize(arg[0], ctx);
        plane_visualize(&p, ctx);
    }
    return 0;
}

static inline const plane_t* as_planes(data_t arg) {
    return (const plane_t*)arg.ptr;
}

// Static: the length, one value per plane, and the unit of what is read
static int plane_accessor_static(data_t arg[], md_unit_t unit, eval_context_t* ctx) {
    if (ctx->backchannel) {
        ctx->backchannel->unit[0] = md_unit_none();
        ctx->backchannel->unit[1] = unit;
    }
    return (int)element_count(arg[0]);
}

static int _plane_centroid(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return plane_accessor_static(arg, md_unit_angstrom(), ctx);
    const plane_t* p = as_planes(arg[0]);
    vec3_t* out = as_vec3_arr(*dst);
    const size_t n = MIN(element_count(arg[0]), type_info_total_element_count(dst->type) / 3);
    for (size_t i = 0; i < n; ++i) out[i] = p ? p[i].centroid : (vec3_t){0};
    return 0;
}

static int _plane_normal(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return plane_accessor_static(arg, md_unit_none(), ctx);
    const plane_t* p = as_planes(arg[0]);
    vec3_t* out = as_vec3_arr(*dst);
    const size_t n = MIN(element_count(arg[0]), type_info_total_element_count(dst->type) / 3);
    for (size_t i = 0; i < n; ++i) out[i] = p ? p[i].normal : (vec3_t){0};
    return 0;
}

static int _plane_radius(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return plane_accessor_static(arg, md_unit_angstrom(), ctx);
    const plane_t* p = as_planes(arg[0]);
    float* out = as_float_arr(*dst);
    const size_t n = MIN(element_count(arg[0]), element_count(*dst));
    for (size_t i = 0; i < n; ++i) out[i] = p ? p[i].radius : 0.0f;
    return 0;
}

static int _plane_deviation(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return plane_accessor_static(arg, md_unit_angstrom(), ctx);
    const plane_t* p = as_planes(arg[0]);
    float* out = as_float_arr(*dst);
    const size_t n = MIN(element_count(arg[0]), element_count(*dst));
    for (size_t i = 0; i < n; ++i) out[i] = p ? p[i].deviation : 0.0f;
    return 0;
}

// angle(p, q): the angle between planes, 0 to 90 degrees, element by element (a single plane against each of the other)
static int _plane_angle(data_t* dst, data_t arg[], eval_context_t* ctx) {
    const size_t na = element_count(arg[0]);
    const size_t nb = element_count(arg[1]);
    if (!dst) {
        if (na != nb && na != 1 && nb != 1) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "angle: %zu planes against %zu, they have to be as many or one of them a single plane", na, nb);
            return STATIC_VALIDATION_ERROR;
        }
        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_radian();
            ctx->backchannel->value_range = (frange_t){0, PI / 2};
        }
        return (int)MAX(na, nb);
    }
    const plane_t* a = as_planes(arg[0]);
    const plane_t* b = as_planes(arg[1]);
    float* out = as_float_arr(*dst);
    const size_t n = MIN(MAX(na, nb), element_count(*dst));
    for (size_t i = 0; i < n; ++i) {
        const plane_t* pa = &a[na == 1 ? 0 : i];
        const plane_t* pb = &b[nb == 1 ? 0 : i];
        const float c = fabsf(vec3_dot(pa->normal, pb->normal));
        out[i] = acosf(CLAMP(c, 0.0f, 1.0f));
    }
    return 0;
}

// ### PUCKER ###
// pucker(rings): how each ring of 5 to 7 atoms is puckered (md_pucker.h), a pucker record per ring, which amplitude(),
// phase(), theta() and pseudorotation() read. The rings are those of the system, in the order of md_pucker_ring_order
// (a sugar from its ring O and anomeric C, proline from N and CA). pucker({i, j, ...}) takes the atoms of one ring in
// the order given instead.

static void pucker_visualize(const int32_t* order, size_t n, eval_context_t* ctx) {
    if (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
        for (size_t k = 0; k < n; ++k) md_bitmask_set_bit(ctx->vis_structure, (uint32_t)order[k]);
    }
    if (!(ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY)) return;
    vec4_t xyzw[MD_RING_MAX_SIZE];
    plane_positions(xyzw, order, n, ctx->cur_state);
    for (size_t k = 0; k < n; ++k) {
        push_line(vertex(vec3_from_vec4(xyzw[k]), COLOR_WHITE), vertex(vec3_from_vec4(xyzw[(k + 1) % n]), COLOR_WHITE), ctx->vis);
    }
}

// The ring of each group, -1 where it is not one of 5 to 7 atoms of the system
static int64_t pucker_group_ring(const plane_rings_t* pr, const md_bitset_t* bf, const eval_context_t* ctx, md_allocator_i* temp) {
    const size_t count = md_bitset_count(*bf);
    if (count < 5 || count > 7) return -1;
    int32_t* idx = md_alloc(temp, sizeof(int32_t) * count);
    md_bitset_iter_extract_indices(idx, count, md_bitset_iter(*bf));
    return plane_rings_find(pr, ctx->sys, idx, count);
}

static int _pucker_rings(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    const md_bitset_t* groups = as_bitset(arg[0]);
    const size_t num = element_count(arg[0]);
    const md_system_t* sys = ctx->sys;

    md_temp_scope_t temp = md_temp_begin_avoid(ctx->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    const plane_rings_t pr = plane_rings_build(sys, temp_alloc);
    int result = 0;

    if (!dst && !ctx->vis) {
        if (ctx->arg_flags[0] & FLAG_DYNAMIC) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "pucker: the rings cannot depend on the frame");
            result = STATIC_VALIDATION_ERROR;
            goto done;
        }
        for (size_t g = 0; g < num; ++g) {
            if (pucker_group_ring(&pr, &groups[g], ctx, temp_alloc) < 0) {
                const size_t n = md_bitset_count(groups[g]);
                char what[32] = "";
                if (num > 1) md_snprintf(what, sizeof(what), " %zu of %zu", g + 1, num);
                if (n < 5 || n > 7) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "pucker: group%s has %zu atoms, a ring of 5 to 7 is needed (ring(size=5:7) gives those)", what, n);
                else                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "pucker: group%s is not a ring of the system; pucker({i, j, ...}) takes atoms in the order of a ring", what);
                result = STATIC_VALIDATION_ERROR;
                goto done;
            }
        }
        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_none();
        }
        result = (int)num;
        goto done;
    }

    pucker_rec_t* out = dst ? (pucker_rec_t*)dst->ptr : NULL;
    const size_t end = MIN(num, dst ? element_count(*dst) : num);
    size_t beg = 0, vend = end;
    if (!dst && ctx->vis && ctx->subscript_ranges) {
        beg = (size_t)MAX(ctx->subscript_ranges[0].beg, 0);
        vend = MIN((size_t)MAX(ctx->subscript_ranges[0].end, 0), end);
    }
    for (size_t g = beg; g < vend; ++g) {
        const int64_t r = pucker_group_ring(&pr, &groups[g], ctx, temp_alloc);
        if (r < 0) {
            if (out) MEMSET(&out[g], 0, sizeof(pucker_rec_t));
            continue;
        }
        int32_t order[MD_RING_MAX_SIZE];
        bool anchored = false;
        const size_t n = md_pucker_ring_order(order, &anchored, sys, (size_t)r);
        const bool sugar = (md_system_ring_flags(sys, (size_t)r) & MD_RING_FLAG_SUGAR) != 0;
        if (out) {
            if (ctx->backchannel) {
                // The static check asking about the argument of an accessor: the shape only
                MEMSET(&out[g], 0, sizeof(pucker_rec_t));
                out[g].cp.size = (uint32_t)n;
                out[g].size = (uint8_t)n;
                out[g].anchored = anchored;
                out[g].sugar = sugar;
                out[g].has_pseudo = sugar && n == 5;
            } else {
                md_pucker_ring(&out[g], sys, ctx->cur_state, (size_t)r);
            }
        }
        if (ctx->vis) pucker_visualize(order, n, ctx);
    }

done:
    md_temp_end(temp);
    return result;
}

static int _pucker_atoms(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    const int* in = as_int_arr(arg[0]);
    const size_t n = element_count(arg[0]);
    const irange_t range = get_atom_range_in_context(ctx->sys, ctx->mol_ctx);
    int32_t order[MD_RING_MAX_SIZE];

    if (!dst && !ctx->vis) {
        if (ctx->arg_flags[0] & FLAG_DYNAMIC) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "pucker: the atoms cannot depend on the frame");
            return STATIC_VALIDATION_ERROR;
        }
        if (n < 5 || n > 7) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "pucker: %zu atoms, a ring of 5 to 7 atoms in order is needed", n);
            return STATIC_VALIDATION_ERROR;
        }
    }
    if (n < 5 || n > 7) return dst ? 0 : 1;
    for (size_t k = 0; k < n; ++k) {
        const int idx = range.beg + in[k] - 1;
        if (idx < range.beg || idx >= range.end) {
            if (!dst && !ctx->vis) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "pucker: atom %d is out of range", in[k]);
            return STATIC_VALIDATION_ERROR;
        }
        order[k] = idx;
    }
    if (!dst && !ctx->vis) {
        // Around a ring: each atom bonded to the next, and the last to the first
        for (size_t k = 0; k < n; ++k) {
            if (md_bond_find(&ctx->sys->bond, order[k], order[(k + 1) % n]) < 0) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "pucker: atoms %d and %d are not bonded, the atoms have to go around a ring", in[k], in[(k + 1) % n]);
                return STATIC_VALIDATION_ERROR;
            }
        }
        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_none();
        }
        return 1;
    }
    if (dst) {
        pucker_rec_t* out = (pucker_rec_t*)dst->ptr;
        if (ctx->backchannel) {
            MEMSET(out, 0, sizeof(pucker_rec_t));
            out->cp.size = (uint32_t)n;
            out->size = (uint8_t)n;
            out->anchored = true;
            out->has_pseudo = n == 5;
        } else {
            md_pucker_ring_atoms(out, ctx->sys, ctx->cur_state, order, n);
        }
    }
    if (ctx->vis) pucker_visualize(order, n, ctx);
    return 0;
}

static inline const pucker_rec_t* as_puckers(data_t arg) {
    return (const pucker_rec_t*)arg.ptr;
}

// What an accessor needs of the rings, checked on their shape (the static check fills it in): 6 membered for theta,
// a furanose for pseudorotation. A phase of a ring without a first atom only warns.
typedef enum pucker_need_t { PUCKER_NEED_ANY, PUCKER_NEED_PHASE, PUCKER_NEED_THETA, PUCKER_NEED_PSEUDO } pucker_need_t;

static int pucker_accessor_static(data_t arg[], pucker_need_t need, md_unit_t unit, eval_context_t* ctx) {
    const pucker_rec_t* p = as_puckers(arg[0]);
    const size_t n = element_count(arg[0]);
    if (p) {
        bool warned = false;
        for (size_t i = 0; i < n; ++i) {
            if (p[i].cp.size == 0) continue;
            if (need == PUCKER_NEED_THETA && p[i].cp.size != 6) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "theta: ring %zu has %u atoms, theta is of six membered rings (ring(size=6) gives those)", i + 1, p[i].cp.size);
                return STATIC_VALIDATION_ERROR;
            }
            if (need == PUCKER_NEED_PSEUDO && !p[i].has_pseudo) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "pseudorotation: ring %zu is not a furanose (ring(kind='sugar', size=5) gives those)", i + 1);
                return STATIC_VALIDATION_ERROR;
            }
            if ((need == PUCKER_NEED_PHASE || need == PUCKER_NEED_THETA) && !p[i].anchored && !p[i].has_pseudo && !warned) {
                EVAL_LOG_WARNING(ctx, ctx->arg_tokens[0], "ring %zu has no atom to start from (a ring of C alike all around): its %s depends on which atom is taken as the first", i + 1, need == PUCKER_NEED_THETA ? "theta" : "phase");
                warned = true;
            }
        }
    }
    if (ctx->backchannel) {
        ctx->backchannel->unit[0] = md_unit_none();
        ctx->backchannel->unit[1] = unit;
        if (need == PUCKER_NEED_PHASE) {
            // Cremer-Pople phi2 in [0, 360[, cut at 0: it moves around the whole circle (the pseudorotation of a five
            // membered ring, the boats and skews of a six), and where it is least likely to be depends on the ring
            ctx->backchannel->period = period_make(0.0, TWO_PI, 0.0);
        } else if (need == PUCKER_NEED_PSEUDO) {
            // P in [0, 360[, cut at 270: the west of the pseudorotation wheel (O4'-exo), the highest barrier of a
            // furanose. North (C3'-endo, P near 18) and south (C2'-endo, near 162) meet through the east. At 0 they
            // would be cut through the north.
            ctx->backchannel->period = period_make(0.0, TWO_PI, 1.5 * PI);
        }
    }
    return (int)n;
}

static int _pucker_amplitude(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return pucker_accessor_static(arg, PUCKER_NEED_ANY, md_unit_angstrom(), ctx);
    const pucker_rec_t* p = as_puckers(arg[0]);
    float* out = as_float_arr(*dst);
    const size_t n = MIN(element_count(arg[0]), element_count(*dst));
    for (size_t i = 0; i < n; ++i) out[i] = p ? p[i].cp.amplitude : 0.0f;
    return 0;
}

static int _pucker_phase(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return pucker_accessor_static(arg, PUCKER_NEED_PHASE, md_unit_radian(), ctx);
    const pucker_rec_t* p = as_puckers(arg[0]);
    float* out = as_float_arr(*dst);
    const size_t n = MIN(element_count(arg[0]), element_count(*dst));
    for (size_t i = 0; i < n; ++i) out[i] = p ? DEG_TO_RAD(p[i].cp.phi2) : 0.0f;
    return 0;
}

static int _pucker_theta(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return pucker_accessor_static(arg, PUCKER_NEED_THETA, md_unit_radian(), ctx);
    const pucker_rec_t* p = as_puckers(arg[0]);
    float* out = as_float_arr(*dst);
    const size_t n = MIN(element_count(arg[0]), element_count(*dst));
    for (size_t i = 0; i < n; ++i) out[i] = p ? DEG_TO_RAD(p[i].cp.theta) : 0.0f;
    return 0;
}

static int _pucker_pseudorotation(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return pucker_accessor_static(arg, PUCKER_NEED_PSEUDO, md_unit_radian(), ctx);
    const pucker_rec_t* p = as_puckers(arg[0]);
    float* out = as_float_arr(*dst);
    const size_t n = MIN(element_count(arg[0]), element_count(*dst));
    for (size_t i = 0; i < n; ++i) out[i] = p ? DEG_TO_RAD(p[i].pseudo_phase) : 0.0f;
    return 0;
}

static int _internal_density(data_t* dst, data_t arg[], eval_context_t* ctx, int axis) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));

    if (dst) {
        float* density[3];
        float* weight[3];

        // Set the distribution weight vector for each dimension
        // And zero the contributuion of each bin
        ASSERT(dst->weight);
        if (axis == -1) {
            // Three distributions: x, y, z
            float* ptr = dst->ptr;
            density[0] = ptr + MD_DIST_BINS * 0;
            weight[0]  = dst->weight + MD_DIST_BINS * 0;
            density[1] = ptr + MD_DIST_BINS * 1;
            weight[1]  = dst->weight + MD_DIST_BINS * 1;
            density[2] = ptr + MD_DIST_BINS * 2;
            weight[2]  = dst->weight + MD_DIST_BINS * 2;

            for (size_t i = 0; i < MD_DIST_BINS; ++i) {
                weight [0][i] = 1.0f;
                density[0][i] = 0.0f;
                weight [1][i] = 1.0f;
                density[1][i] = 0.0f;
                weight [2][i] = 1.0f;
                density[2][i] = 0.0f;
            }
        } else {
            density[axis] = (float*)dst->ptr;
            weight[axis]  = dst->weight;
            for (size_t i = 0; i < MD_DIST_BINS; ++i) {
                weight [axis][i] = 1.0f;
                density[axis][i] = 0.0f;
            }
        }

        const md_bitset_t* src_bf = as_bitset(arg[0]);
        md_bitset_t bf = _internal_flatten_bf(src_bf, element_count(arg[0]), ctx->temp_alloc);
        if (ctx->mol_ctx) {
            bf = md_bitset_and(bf, *ctx->mol_ctx, ctx->temp_alloc);
        }

        mat3_t A = { 0 };
        md_unitcell_A_extract(A.elem, &ctx->ref_state->unitcell);

        // Reference point
        vec3_t rc = mat3_mul_vec3(A, vec3_set1(0.5f));
        // Reference extent
        vec3_t re = mat3_diag(A);

        vec4_t ext = vec4_from_vec3(re, 0);
        if (axis == -1) {
            // We use the max axis as we cannot output three individual value_ranges
            ext = vec4_set1(MAX(re.x, MAX(re.y, re.z)));
            ext.w = 0;
        }

        md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);

        size_t count = md_bitset_count(bf);
        vec4_t* xyzw = md_temp_alloc_array(temp, vec4_t, count);

        extract_xyzw_vec4(xyzw, ctx->cur_state->xyz, ctx->atom_mass, &bf);
        //md_util_pbc_vec4(xyzw, count, &ctx->sys->unit_cell);

        vec4_t mid_point = vec4_from_vec3(rc, 0);
        vec4_t period = vec4_from_vec3(re, 0);
        vec4_t inv_ext = {
            ext.x > 0.0f ? 1.0f / ext.x : 0.0f,
            ext.y > 0.0f ? 1.0f / ext.y : 0.0f,
            ext.z > 0.0f ? 1.0f / ext.z : 0.0f,
            0.0f
        };

        vec4_t min_point = vec4_sub(mid_point, vec4_mul1(ext, 0.5f));

        // @TODO: This is conceptually not correct in the periferi of the system as it does not count multiple occurences of atoms
        // This should be resolved with spatial hashing and counting occurrences within the cells
        // But the current inplementation of the spatial hash does not facilitate it properly...


        for (size_t i = 0; i < count; ++i) {
            vec4_t xyzw_dp = vec4_deperiodize_ortho(xyzw[i], mid_point, period);
            vec4_t fc = vec4_mul(vec4_sub(xyzw_dp, min_point), inv_ext);
            float w = xyzw[i].w;

            int bin_idx[3] = {
                CLAMP((int)(fc.x * MD_DIST_BINS), 0, MD_DIST_BINS - 1),
                CLAMP((int)(fc.y * MD_DIST_BINS), 0, MD_DIST_BINS - 1),
                CLAMP((int)(fc.z * MD_DIST_BINS), 0, MD_DIST_BINS - 1),
            };

            if (axis == -1) {
                density[0][bin_idx[0]] += w;
                density[1][bin_idx[1]] += w;
                density[2][bin_idx[2]] += w;
            } else {
                density[axis][bin_idx[axis]] += w;
            }
        }

        // The slab of a bin along an axis: the bin's thickness times the cross section of the cell. With all three
        // axes (axis == -1) every axis is binned over the longest extent, so the thickness is that over the bin
        // count, while the cross section is still the cell's: the cell's volume over the bins would not be.
        const double slab_vol[3] = {
            (double)ext.x / MD_DIST_BINS * re.y * re.z,
            (double)ext.y / MD_DIST_BINS * re.x * re.z,
            (double)ext.z / MD_DIST_BINS * re.x * re.y,
        };

        // Convert atomic mass / Ångström^3 to SI units KG/M^3
        const double factor[3] = {
            slab_vol[0] > 0 ? 1660.5390666 / slab_vol[0] : 0.0,
            slab_vol[1] > 0 ? 1660.5390666 / slab_vol[1] : 0.0,
            slab_vol[2] > 0 ? 1660.5390666 / slab_vol[2] : 0.0,
        };

        if (axis == -1) {
            for (size_t i = 0; i < MD_DIST_BINS; ++i) {
                density[0][i] = (float)(density[0][i] * factor[0]);
                density[1][i] = (float)(density[1][i] * factor[1]);
                density[2][i] = (float)(density[2][i] * factor[2]);
            }
        } else {
            for (size_t i = 0; i < MD_DIST_BINS; ++i) {
                density[axis][i] = (float)(density[axis][i] * factor[axis]);
            }
        }

        md_temp_end(temp);
    } else {
        // STATIC CHECK
        const md_unitcell_t* unit_cell = &ctx->ref_state->unitcell;

        if (axis == -1) {
            if (unit_cell->flags == 0) {
                EVAL_LOG_ERROR(ctx, ctx->op_token, "The unitcell is not defined");
                return STATIC_VALIDATION_ERROR;
            }
        } else {
            uint32_t required_flags = 0;
            if (axis == 0)      required_flags = MD_UNITCELL_PBC_X;
            else if (axis == 1) required_flags = MD_UNITCELL_PBC_Y;
            else if (axis == 2) required_flags = MD_UNITCELL_PBC_Z;

            if (!(ctx->cur_state->unitcell.flags & required_flags)) {
                EVAL_LOG_ERROR(ctx, ctx->op_token, "The unitcell is not defined along the given axis");
                return STATIC_VALIDATION_ERROR;
            }
        }
        if (ctx->backchannel) {
            // Reference half extent
            mat3_t A = { 0 };
            md_unitcell_A_extract_float(A.elem, unit_cell);
            vec3_t re = vec3_mul1(mat3_diag(A), 0.5f);

            // We use the max axis half ext as we cannot output three individual value_ranges
            float rad = 0;
            if (axis == -1) {
                rad  = MAX(re.x, MAX(re.y, re.z));
            } else {
                rad = re.elem[axis];
            }

            ctx->backchannel->unit[0] = md_unit_angstrom();
            ctx->backchannel->unit[1] = md_unit_div(md_unit_kilogram(), md_unit_pow(md_unit_meter(), 3));
            ctx->backchannel->value_range = (frange_t){-rad, rad};
        }
    }
    return 0;
}

static int _density(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _internal_density(dst, arg, ctx, -1);
}

static int _density_x(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _internal_density(dst, arg, ctx, 0);
}

static int _density_y(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _internal_density(dst, arg, ctx, 1);
}

static int _density_z(data_t* dst, data_t arg[], eval_context_t* ctx) {
    return _internal_density(dst, arg, ctx, 2);
}

static int _coordinate(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));

    if (dst || ctx->vis) {
        if (dst) {
            ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT3_ARR));
            md_array(vec3_t) pos = coordinate_extract(arg[0], ctx);
            // As for coord_x: what fits is written and the rest zeroed
            const size_t bytes = MIN(dst->size, md_array_bytes(pos));
            MEMCPY(dst->ptr, pos, bytes);
            MEMSET((char*)dst->ptr + bytes, 0, dst->size - bytes);
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
            ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            float* x = as_float_arr(*dst);
            // Laid out [context][N] when evaluated in a context, where an argument whose length varies over frames can
            // have another length than the one fixed at compile time: what fits is written and the rest zeroed
            const size_t n = MIN(dst->size / sizeof(float), md_array_size(pos));
            MEMSET((char*)dst->ptr + n * sizeof(float), 0, dst->size - n * sizeof(float));
            for (size_t i = 0; i < n; ++i) {
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
            ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            float* y = as_float_arr(*dst);
            // Laid out [context][N] when evaluated in a context, where an argument whose length varies over frames can
            // have another length than the one fixed at compile time: what fits is written and the rest zeroed
            const size_t n = MIN(dst->size / sizeof(float), md_array_size(pos));
            MEMSET((char*)dst->ptr + n * sizeof(float), 0, dst->size - n * sizeof(float));
            for (size_t i = 0; i < n; ++i) {
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
            ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            float* z = as_float_arr(*dst);
            // Laid out [context][N] when evaluated in a context, where an argument whose length varies over frames can
            // have another length than the one fixed at compile time: what fits is written and the rest zeroed
            const size_t n = MIN(dst->size / sizeof(float), md_array_size(pos));
            MEMSET((char*)dst->ptr + n * sizeof(float), 0, dst->size - n * sizeof(float));
            for (size_t i = 0; i < n; ++i) {
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
            ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT2_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            vec2_t* xy = as_vec2_arr(*dst);
            // Laid out [context][N] when evaluated in a context, where an argument whose length varies over frames can
            // have another length than the one fixed at compile time: what fits is written and the rest zeroed
            const size_t n = MIN(dst->size / sizeof(vec2_t), md_array_size(pos));
            MEMSET((char*)dst->ptr + n * sizeof(vec2_t), 0, dst->size - n * sizeof(vec2_t));
            for (size_t i = 0; i < n; ++i) {
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
            ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT2_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            vec2_t* xz = as_vec2_arr(*dst);
            // Laid out [context][N] when evaluated in a context, where an argument whose length varies over frames can
            // have another length than the one fixed at compile time: what fits is written and the rest zeroed
            const size_t n = MIN(dst->size / sizeof(vec2_t), md_array_size(pos));
            MEMSET((char*)dst->ptr + n * sizeof(vec2_t), 0, dst->size - n * sizeof(vec2_t));
            for (size_t i = 0; i < n; ++i) {
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
            ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT2_ARR));
            const vec3_t* pos = coordinate_extract(arg[0], ctx);
            vec2_t* yz = as_vec2_arr(*dst);
            // Laid out [context][N] when evaluated in a context, where an argument whose length varies over frames can
            // have another length than the one fixed at compile time: what fits is written and the rest zeroed
            const size_t n = MIN(dst->size / sizeof(vec2_t), md_array_size(pos));
            MEMSET((char*)dst->ptr + n * sizeof(vec2_t), 0, dst->size - n * sizeof(vec2_t));
            for (size_t i = 0; i < n; ++i) {
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
    const md_bitset_t* exclusion_masks;
    float min_cutoff;
    float inv_cutoff_range;
    float* bins;
    int32_t num_bins;
    uint64_t total_count;
} rdf_payload_t;

#if 0
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

// iterator version which excludes indices that are present in the supplied exclusion_mask
bool rdf_iter_excl_mask(const md_spatial_hash_elem_t* elem_arr, int mask, void* user_param) {
    rdf_payload_t* data = user_param;

    while (mask) {
        const int idx = ctz32(mask);
        if (!md_bitset_test(*data->exclusion_mask, elem_arr[idx].idx)) {
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
#endif

static inline void rdf_increment_bin(rdf_payload_t* data, float d) {
    int32_t bin_idx = (int32_t)(((d - data->min_cutoff) * data->inv_cutoff_range) * data->num_bins);
    bin_idx = CLAMP(bin_idx, 0, data->num_bins - 1);
    data->bins[bin_idx] += 1.0f;
    data->total_count += 1;
}

// @NOTE: This can be vectorized
void rdf_cb(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    (void)i_idx;
    (void)j_idx;
    rdf_payload_t* data = user_param;
    const float min_r2 = data->min_cutoff * data->min_cutoff;
    for (size_t k = 0; k < num_pairs; ++k) {
        if (ij_dist2[k] < min_r2) {
            continue;
        }
        const float d = sqrtf(ij_dist2[k]);
        rdf_increment_bin(data, d);
    }
}

void rdf_cb_excl_mask(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    rdf_payload_t* data = user_param;
    const float min_r2 = data->min_cutoff * data->min_cutoff;
    for (size_t k = 0; k < num_pairs; ++k) {
        if (ij_dist2[k] < min_r2) {
            continue;
        }
        uint32_t i = i_idx[k];
        uint32_t j = j_idx[k];
        if (!md_bitset_test(data->exclusion_masks[i], j)) {
            const float d = sqrtf(ij_dist2[k]);
            rdf_increment_bin(data, d);
        }
    }
}

static inline double sphere_volume(double r) {
    return (4.0 / 3.0) * PI * (r * r * r);
}

static void compute_rdf(float bins[], float weights[], int num_bins, const data_t arg[2], float min_cutoff, float max_cutoff, eval_context_t* ctx) {
    const float inv_cutoff_range = 1.0f / (max_cutoff - min_cutoff);

    size_t total_count = 0;

    // Add a small bias to avoid adding 'zero' distances, also to remove *self* when counting atomic RDFs
    min_cutoff = MAX(min_cutoff, 1e-3f);

    md_array(int32_t) ref_idx = coordinate_extract_indices(arg[0], ctx);
    size_t ref_len  = 0;

    // IF we have multiple atomsets supplied as reference, we calculate a RDF based on the COM of each atomset and excluded *self* when counting.
    const bool use_exclusion_masks = (arg[0].type.base_type == TYPE_ATOMSET) && (element_count(arg[0]) > 1);

    md_coord_stream_t ref_stream = { 0 };
    md_spatial_acc_flags_t ref_flags = 0;
    if (ref_idx) {
        ref_len = md_array_size(ref_idx);
        ref_stream = md_coord_stream_from_aos((const float*)ctx->cur_state->xyz, sizeof(vec3_t), ref_idx, ref_len);
        ref_flags = MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX;
    } else {
        md_array(vec3_t) ref_pos = coordinate_extract(arg[0], ctx);
        ref_len = md_array_size(ref_pos);
        ref_stream = md_coord_stream_from_aos((const float*)ref_pos, sizeof(vec3_t), NULL, ref_len);
    }

    md_array(int32_t) trg_idx = coordinate_extract_indices(arg[1], ctx);
    size_t trg_len = 0;

    md_coord_stream_t trg_stream = { 0 };
    md_spatial_acc_flags_t trg_flags = 0;
    if (trg_idx) {
        trg_len = md_array_size(trg_idx);
        trg_stream = md_coord_stream_from_aos((const float*)ctx->cur_state->xyz, sizeof(vec3_t), trg_idx, trg_len);
        trg_flags = MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX;
    } else {
        md_array(vec3_t) trg_pos = coordinate_extract(arg[1], ctx);
        trg_len = md_array_size(trg_pos);
        trg_stream = md_coord_stream_from_aos((const float*)trg_pos, sizeof(vec3_t), NULL, trg_len);
    }

    {
        md_spatial_acc_t sa = { .alloc = ctx->temp_alloc };
        md_spatial_acc_init(&sa, &(md_spatial_acc_desc_t){ .coords = &trg_stream, .cutoff = max_cutoff, .unitcell = &ctx->cur_state->unitcell, .flags = trg_flags });

        rdf_payload_t payload = {
            .exclusion_masks = as_bitset(arg[0]),
            .min_cutoff = min_cutoff,
            .inv_cutoff_range = inv_cutoff_range,
            .bins = bins,
            .num_bins = num_bins,
            .total_count = 0,
        };

        md_spatial_acc_pair_callback_t cb = use_exclusion_masks ? rdf_cb_excl_mask : rdf_cb;
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&sa, &ref_stream, max_cutoff, cb, &payload, ref_flags);

        total_count = payload.total_count;
    }

    // Reference rho
    const double total_vol = sphere_volume(max_cutoff) - sphere_volume(min_cutoff);
    const double ref_rho   = total_count / total_vol;

    // Compute normalization factor of bins
    // With respect to the nominal distribution rho
    // Each bin is normalized with respect to the nominal distribution
    // The shell of bin i is [min + i dr, min + (i + 1) dr), as rdf_increment_bin bins it
    const double dr = (max_cutoff - min_cutoff) / (float)num_bins;
    double prev_sphere_vol = sphere_volume(min_cutoff);
    for (int64_t i = 0; i < num_bins; ++i) {
        const double sphere_vol = sphere_volume(min_cutoff + (i + 1) * dr);
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

// ### DISTRIBUTIONS ###
// A bin is shown as value / weight, and the evaluation sums both over frames. Only operations which give the same
// result whether they are applied per frame or to those sums are defined: scaling, which scales the values and
// leaves the weights, and pooling, which adds both.

static size_t dist_bin_count(const data_t* d) {
    return type_info_total_element_count(d->type);
}

static float dist_weight(const data_t* d, size_t i) {
    return d->weight ? d->weight[i] : 1.0f;
}

// The number of distributions in the result: as many as in the first argument. Asked for with dst NULL (the
// length is not deduced from the arguments' dimensions, which would take [1][bins] for a leading one to replicate).
static int dist_count(const data_t arg[]) {
    return arg[0].type.dim[0];
}

static int dist_scaled(data_t* dst, data_t arg[], float scl) {
    if (!dst) {
        return dist_count(arg);
    }
    {
        ASSERT(dst->weight);
        const size_t n = dist_bin_count(&arg[0]);
        ASSERT(n == dist_bin_count(dst));
        const float* v = (const float*)arg[0].ptr;
        float* out = (float*)dst->ptr;
        for (size_t i = 0; i < n; ++i) {
            out[i] = v[i] * scl;
            dst->weight[i] = dist_weight(&arg[0], i);
        }
    }
    return 0;
}

static int _dist_scale(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)ctx;
    return dist_scaled(dst, arg, as_float(arg[1]));
}

static int _dist_div(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)ctx;
    return dist_scaled(dst, arg, dst ? 1.0f / as_float(arg[1]) : 0.0f);
}

static int _dist_pool(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)ctx;
    if (!dst) {
        return dist_count(arg);
    }
    {
        ASSERT(dst->weight);
        const size_t n = dist_bin_count(&arg[0]);
        ASSERT(n == dist_bin_count(&arg[1]) && n == dist_bin_count(dst));
        const float* a = (const float*)arg[0].ptr;
        const float* b = (const float*)arg[1].ptr;
        float* out = (float*)dst->ptr;
        for (size_t i = 0; i < n; ++i) {
            out[i] = a[i] + b[i];
            dst->weight[i] = dist_weight(&arg[0], i) + dist_weight(&arg[1], i);
        }
    }
    return 0;
}

static int internal_rdf(data_t* dst, data_t arg[], float min_cutoff, float max_cutoff, eval_context_t* ctx) {
    if (dst || ctx->vis) {
        const int num_bins = MD_DIST_BINS;
       
        if (dst) {
            ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_DISTRIBUTION));
            ASSERT(dst->ptr && dst->weight);
            float* bins    = (float*)dst->ptr;
            float* weights = dst->weight;
            compute_rdf(bins, weights, num_bins, arg, min_cutoff, max_cutoff, ctx);
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
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "empty reference positions");
            return STATIC_VALIDATION_ERROR;
        }
        if (trg_len <= 0) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[1], "empty target positions");
            return STATIC_VALIDATION_ERROR;
        }
        if (min_cutoff < 0.0f || max_cutoff <= min_cutoff) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[2], "Invalid cutoff");
            return STATIC_VALIDATION_ERROR;
        }
        if (ctx->arg_flags[2] & FLAG_DYNAMIC) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[2], "Cutoff needs to have a static value in order to determine the upper bound if the distribution");
            return STATIC_VALIDATION_ERROR;
        }
        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_angstrom();
            ctx->backchannel->unit[1] = md_unit_none();
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
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_COORDINATE_ARR));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FRANGE));

    const frange_t cutoff = as_frange(arg[2]);
    return internal_rdf(dst, arg, cutoff.beg, cutoff.end, ctx);
}

typedef enum {
    COUNT_TYPE_UNKNOWN = 0,
    COUNT_TYPE_ATOM = 1,
    COUNT_TYPE_RESIDUE = 2,
    COUNT_TYPE_PART = 3,
    COUNT_TYPE_CHAIN = 4,
    COUNT_TYPE_MOLECULE = 5,
    COUNT_TYPE_COUNT
} count_type_t;

static const str_t count_type_str[] = {
    STR_INIT(""),
    STR_INIT("atom"),
    STR_INIT("residue"),
    STR_INIT("part"),
    STR_INIT("chain"),
    STR_INIT("molecule"),
};

count_type_t count_type_from_str(str_t str) {
    for (size_t i = 0; i < COUNT_TYPE_COUNT; ++i) {
        if (str_eq(str, count_type_str[i])) {
            return (count_type_t)i;
        }
    }
    return COUNT_TYPE_UNKNOWN;
}

size_t internal_count(const md_bitset_t* bf_arr, size_t bf_len, count_type_t count_type, eval_context_t* ctx) {
    ASSERT(ctx);
    if (!bf_arr || bf_len == 0) {
        return 0;
    }

    md_bitset_t bf = _internal_flatten_bf(bf_arr, bf_len, ctx->temp_alloc);

    if (ctx->mol_ctx) {
        bf = md_bitset_and(bf, *ctx->mol_ctx, ctx->temp_alloc);
    }
    
    size_t count = 0;

    switch (count_type) {
        case COUNT_TYPE_ATOM: {
            return md_bitset_count(bf);
        }
        case COUNT_TYPE_RESIDUE:
            for (size_t i = 0; i < ctx->sys->residue.count; ++i) {
                md_urange_t range = md_system_residue_atom_range(ctx->sys, i);
                if (md_bitset_count_range(bf, range.beg, range.end) > 0) {
                    count += 1;
                }
            }
            break;
        case COUNT_TYPE_PART: {
            for (size_t p = 0; p < md_system_part_count(ctx->sys); ++p) {
                const md_urange_t range = md_system_part_atom_range(ctx->sys, p);
                count += md_bitset_count_range(bf, range.beg, range.end) > 0;
            }
            break;
        }
        case COUNT_TYPE_CHAIN: {
            for (size_t c = 0; c < md_system_chain_count(ctx->sys); ++c) {
                count += chain_in_context(ctx->sys, c, &bf) ? 1 : 0;
            }
            break;
        }
        case COUNT_TYPE_MOLECULE: {
            const size_t num_molecules = md_system_molecule_count(ctx->sys);
            for (size_t i = 0; i < num_molecules; ++i) {
                const md_molecule_t mol = md_system_molecule(ctx->sys, i);

                // Counted when any of its atoms is in the set
                for (size_t k = 0; k < mol.count; ++k) {
                    if (md_bitset_test(bf, (uint32_t)mol.atom_idx[k])) {
                        count += 1;
                        break;
                    }
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
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));

    if (dst) {
        ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
        const md_bitset_t* bf = as_bitset(arg[0]);
        size_t len = element_count(arg[0]);
        size_t count = internal_count(bf, len, COUNT_TYPE_ATOM, ctx);
        as_float(*dst) = (float)count;
    }

    return 0;
}

static int _count_with_arg(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_STRING));

    if (dst) {
        ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));
        const md_bitset_t* bf = as_bitset(arg[0]);
        size_t len = element_count(arg[0]);
        count_type_t type = count_type_from_str(as_string(arg[1]));
        as_float(*dst) = (float)internal_count(bf, len, type, ctx);
    } else {
        str_t str = as_string(arg[1]);
        count_type_t type = count_type_from_str(str);
        if (type == COUNT_TYPE_UNKNOWN) {
            md_str_writer_t sb = md_str_writer_create(ctx->temp_alloc);
            for (int i = 1; i < COUNT_TYPE_COUNT; ++i) {
                md_str_writer_fmt(&sb, "\""STR_FMT"\"\n", STR_ARG(count_type_str[i]));
            }
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[1], "Unknown argument: '"STR_FMT"', valid arguments are:\n"STR_FMT, STR_ARG(str), STR_ARG(md_str_writer_to_str(&sb)));
            return -1;
        }
    }

    return 0;
}

// This is quirky and wierd.
// It should essentially check if all the supplied atomsets are 'equivalent'
// Equivalent in this case means that all of them represent the same type of structures within the system
// We do this now by checking against the atom elements if available (which is not the case in CoarseGrained simulations)
// In such case we revert to matching Atom labels
// To resolve this PROPERLY, one needs to determine the topology of each supplied substructure and match it as a graph.

static inline bool are_bitfields_equivalent(const md_bitset_t bitfields[], size_t num_bitfields, const md_system_t* sys) {
    // Number of bits should match.
    // The atomic element of each set bit should match.

    const md_bitset_t* ref_bf = &bitfields[0];
    const int64_t ref_count = md_bitset_count(*ref_bf);

    // We compare against the first one, which is considered to be the reference
    for (size_t i = 1; i < num_bitfields; ++i) {
        const md_bitset_t* bf = &bitfields[i];
        const int64_t count = md_bitset_count(*bf);
        if (count != ref_count) {
            return false;
        }

        md_bitset_iter_t scan_it = md_bitset_iter(*bf);
        while (md_bitset_iter_next(&scan_it)) {
            const int64_t idx = scan_it.idx;
            const int64_t ref_idx = ref_bf->beg + (idx - bf->beg);
            if (!md_bitset_test(*ref_bf, ref_idx)) {
                return false;
            }

            if (sys->atom.type_idx[idx] == sys->atom.type_idx[ref_idx]) {
                continue;
            }

            md_atomic_number_t z_idx = md_atom_atomic_number(&sys->atom, idx);
            md_atomic_number_t z_ref = md_atom_atomic_number(&sys->atom, ref_idx);
            if (z_idx != 0 && z_ref != 0 && z_idx != z_ref) {
                return false;
            } else {
                // Test the type labels
                str_t lbl_idx = md_atom_name(&sys->atom, idx);
                str_t lbl_ref = md_atom_name(&sys->atom, ref_idx);
                if (!str_eq(lbl_idx, lbl_ref)) {
                    return false;
                }
            }
        }
    }

    return true;
}

static inline void populate_volume(float* vol, const vec3_t* xyz, const int* idx, size_t count, const md_bitset_t* exclusion_mask, mat4_t M) {
    const vec4_t dim = vec4_set(MD_VOL_DIM, MD_VOL_DIM, MD_VOL_DIM, 0);
    for (size_t i = 0; i < count; ++i) {
        if (exclusion_mask) {
            if (md_bitset_test(*exclusion_mask, idx[i])) {
                continue;
            }
        }

        const vec4_t coord = mat4_mul_vec4(M, vec4_from_vec3(xyz[i], 1.0f));
        
        const md_128 a = md_mm_cmplt_ps(md_mm_setzero_ps(), coord.m128);
        const md_128 b = md_mm_cmplt_ps(coord.m128, dim.m128);
        const md_128 c = md_mm_and_ps(a, b);

        // Count the number of lanes which is not zero
        // 0x7 = 1 + 2 + 4 means that first three lanes (x,y,z) are within min and max
        if (md_mm_movemask_ps(c) == 0x7) {
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
    const md_bitset_t* exclusion_mask;
} sdf_payload_t;

// Spatial acceleration callback for point query
void sdf_cb(const uint32_t in_idx[], const float in_x[], const float in_y[], const float in_z[], size_t num_points, void* user_param) {
    // The points come in the periodic images around the query centre the spatial acc works in, and M is
    // built against exactly that centre (see md_spatial_acc_aabb_query_center), so M maps them straight
    // onto the grid. The query box encloses the rotated volume, so its corners hold points outside the
    // volume: those are dropped. Clamping them piled their counts onto the faces of the volume.

    sdf_payload_t* data = user_param;

    for (size_t i = 0; i < num_points; ++i) {
        if (data->exclusion_mask && md_bitset_test(*data->exclusion_mask, in_idx[i])) {
            continue;
        }

        const vec4_t c = mat4_mul_vec4(data->M, vec4_set(in_x[i], in_y[i], in_z[i], 1.0f));
        if (!(c.x >= 0.0f && c.x < (float)MD_VOL_DIM &&
              c.y >= 0.0f && c.y < (float)MD_VOL_DIM &&
              c.z >= 0.0f && c.z < (float)MD_VOL_DIM)) {
            continue;
        }

        const uint32_t ix = MIN((uint32_t)c.x, MD_VOL_DIM - 1);
        const uint32_t iy = MIN((uint32_t)c.y, MD_VOL_DIM - 1);
        const uint32_t iz = MIN((uint32_t)c.z, MD_VOL_DIM - 1);
        const size_t vol_idx = iz * (MD_VOL_DIM * MD_VOL_DIM) + iy * MD_VOL_DIM + ix;
        ASSERT(vol_idx < MD_VOL_DIM * MD_VOL_DIM * MD_VOL_DIM);
        data->vol[vol_idx] += 1.0f;
    }
}

static int _sdf(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET_ARR));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_ATOMSET));
    ASSERT(is_type_directly_compatible(arg[2].type, (type_info_t)TI_FLOAT));

    const size_t num_ref_bitfields = element_count(arg[0]);
    const md_bitset_t* ref_bf_arr = as_bitset(arg[0]);
    const md_bitset_t* trg_bf     = as_bitset(arg[1]);
    float cutoff = as_float(arg[2]);

    const md_bitset_t* ref_bf = &ref_bf_arr[0];

    // We currently only support atomsets which represent equivalent structures. 
    // his each atomset should contain the same number of atoms, and the n'th atom within the atomset
    // should match the n'th atom in all other atomsets.
    // This may be loosened in the future to allow aligning structures which are not strictly equivalent.
    int result = 0;

    if (dst || ctx->vis) {
        if (!ctx->ref_state || !ctx->ref_state->xyz) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "sdf: the system has no reference coordinates to align with");
            return -1;
        }
        md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);

        // This could happen if we have dynamic length as input
        if (num_ref_bitfields == 0) goto done;

        const size_t ref_size = md_bitset_count(*ref_bf);
        const size_t trg_size = md_bitset_count(*trg_bf);

        if (ref_size == 0) goto done;
        if (trg_size == 0) goto done;

        float* vol = 0;
        if (dst) {
            ASSERT(dst->ptr);
            ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_VOLUME));
            vol = as_float_arr(*dst);
        }

        // Allocate temporary memory.
        // [0] is the reference template, taken from the configuration the topology was inferred from.
        // [1] is refilled below from the current frame for each of the equivalent structures.
        vec4_t* ref_xyzw[2] = {
            md_temp_alloc_array(temp, vec4_t, ref_size),
            md_temp_alloc_array(temp, vec4_t, ref_size)
        };

        const float* ref_w = (const float*)ctx->atom_mass;

        vec3_t ref_com[2] = {0};

        // The reference template. Placed into mutually consistent images against its own circular
        // mean, which needs no image assignment to compute and so cannot be skewed by whichever
        // images the stored configuration happens to use.
        // @NOTE: against ref_state's OWN cell. The reference configuration and the current frame need
        // not share one, which matters under a varying cell.
        extract_xyzw_vec4(ref_xyzw[0], ctx->ref_state->xyz, ref_w, ref_bf);
        md_util_deperiodize_self_vec4(ref_xyzw[0], ref_size, &ctx->ref_state->unitcell, &ref_com[0]);

        // Fetch target positions
        int*    trg_idx = md_temp_alloc_array(temp, int, trg_size);
        md_bitset_iter_extract_indices(trg_idx, trg_size, md_bitset_iter(*trg_bf));

        md_coord_stream_t stream = md_coord_stream_from_aos((const float*)ctx->cur_state->xyz, sizeof(vec3_t), trg_idx, trg_size);
		md_spatial_acc_t spatial_acc = { .alloc = temp.arena };
        md_spatial_acc_init(&spatial_acc, &(md_spatial_acc_desc_t){ .coords = &stream, .cutoff = cutoff, .unitcell = &ctx->cur_state->unitcell, .flags = MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX });

        // A for alignment matrix, Align eigen vectors with axis x,y,z etc.
        mat3_eigen_t eigen = mat3_eigen(mat3_covariance_matrix_vec4(ref_xyzw[0], 0, ref_size, ref_com[0]));
        const mat3_t A3 = mat3_transpose(eigen.vectors);
        mat4_t A = mat4_from_mat3(A3);

        // V for volume matrix scale and align with the volume which we aim to populate with density
        mat4_t V = compute_volume_matrix(cutoff);
        mat4_t VA = mat4_mul(V, A);

        if (ctx->vis && ctx->vis_flags & MD_SCRIPT_VISUALIZE_SDF) {
            ctx->vis->sdf.extent = cutoff;
        }

        // Scale against which an image ambiguity is judged. Left at FLT_MAX when there is no cell,
        // which disables the guard below - without periodicity there is no ambiguity to guard against.
        float min_cell_extent = FLT_MAX;
        if (md_unitcell_flags(&ctx->cur_state->unitcell) & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) {
            vec3_t ext = {0};
            md_unitcell_diag_extract_float(ext.elem, &ctx->cur_state->unitcell);
            for (int k = 0; k < 3; ++k) {
                if (ext.elem[k] > 0.0f) min_cell_extent = MIN(min_cell_extent, ext.elem[k]);
            }
        }

        for (size_t i = 0; i < num_ref_bitfields; ++i) {
            const md_bitset_t* bf = &ref_bf_arr[i];

            // The atomsets are validated as equivalent during the static check, but a dynamic
            // expression can still hand us a differently sized one at evaluation time. The buffers
            // below are sized for ref_size, so a mismatch would be an overrun rather than a bad result.
            if (md_bitset_count(*bf) != ref_size) continue;

            extract_xyzw_vec4(ref_xyzw[1], ctx->cur_state->xyz, ref_w, bf);

            // Rotation, centre and each point's periodic image are solved for together, so the images
            // are chosen to minimise the very residual this alignment is about to be judged by. No
            // topology is consulted: the copies need not be whole structures, or even bonded.
            mat3_t R = {0};
            const float residual = md_util_optimal_rotation_pbc_vec4_iter(&R, &ref_com[1], ref_xyzw[1], ref_xyzw[0], ref_com[0], ref_xyzw[1], ref_size, &ctx->cur_state->unitcell, 8, 1.0e-6f);

            // The residual is the margin to an ambiguous image choice. Approaching half a cell means
            // no per point image assignment can represent this copy and the alignment is meaningless,
            // so drop the copy rather than smear a wrong orientation into the volume.
            if (residual > 0.25f * min_cell_extent) continue;

            // ref_com[1] is in the image the copy arrived in (see md_util_optimal_rotation_pbc_vec4_iter),
            // which is where it is drawn, so the visualized frame sits on the atoms.
            mat4_t RT = mat4_mul(mat4_from_mat3(R), mat4_translate(-ref_com[1].x, -ref_com[1].y, -ref_com[1].z));

            if (vol) {
                // The neighbours come back around the query centre as the spatial acc folds it into the
                // cell, which is a whole lattice vector off ref_com[1] whenever the copy sits outside the
                // cell - trajectories with whole molecules, or unwrapped ones. Built against ref_com[1],
                // every one of those frames was mapped a cell away from the grid and clamped onto its faces.
                const double aabb_cen[3] = { ref_com[1].x, ref_com[1].y, ref_com[1].z };
                double query_cen[3];
                md_spatial_acc_aabb_query_center(query_cen, &spatial_acc, aabb_cen);
                const mat4_t RT_query = mat4_mul(mat4_from_mat3(R), mat4_translate(-(float)query_cen[0], -(float)query_cen[1], -(float)query_cen[2]));

                sdf_payload_t payload = {
                    .M = mat4_mul(VA, RT_query),
                    .vol = vol,
                    .exclusion_mask = bf,
                };

                // The volume is a cube of half extent 'cutoff' in the aligned frame. Its extent along each
                // world axis is the query box, which must enclose all of it: a box of half extent 'cutoff'
                // in world space misses the corners of any rotated volume.
                const mat3_t Q = mat3_mul(A3, R);   // world -> aligned
                double aabb_rad[3];
                for (int k = 0; k < 3; ++k) {
                    aabb_rad[k] = cutoff * (fabs(Q.elem[k][0]) + fabs(Q.elem[k][1]) + fabs(Q.elem[k][2]));
                }
				md_spatial_acc_for_each_point_in_aabb(&spatial_acc, aabb_cen, aabb_rad, sdf_cb, &payload);
            }
            if (ctx->vis && ctx->vis_flags & MD_SCRIPT_VISUALIZE_SDF) {
                md_allocator_i* alloc = ctx->vis->alloc;
                ASSERT(alloc);
                md_bitset_t bf_cpy = md_bitset_copy(*bf, alloc);

                mat4_t ART = mat4_mul(A, RT);
                md_array_push(ctx->vis->sdf.matrices, ART, alloc);
                md_array_push(ctx->vis->sdf.structures, bf_cpy, alloc);
            }
        }

done:
        md_temp_end(temp);
    } else {
        // Validation
        if (num_ref_bitfields < 1) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "Number of atomsets which serve as reference frame must be 1 or more");
            return -1;
        }

        const size_t ref_bit_count = md_bitset_count(*ref_bf);
        if (ref_bit_count == 0) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The supplied reference atomset(s) are empty");
            return -1;
        }

        // Test for equivalence
        if (!are_bitfields_equivalent(ref_bf_arr, num_ref_bitfields, ctx->sys)) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "The supplied reference atomsets are not identical: the number of atoms and their corresponding elements do not match between all supplied atomsets");
            return -1;
        }

        const size_t target_bit_count = md_bitset_count(*trg_bf);
        if (target_bit_count == 0) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[1], "The supplied target atomset is empty");
            return -1;
        }

        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_angstrom();
            ctx->backchannel->unit[1] = md_unit_div(md_unit_count(), md_unit_pow(md_unit_angstrom(), 3));
            ctx->backchannel->value_range = (frange_t){0, FLT_MAX};
        }
    }

    return result;
}

// Construct a ratio of void space over total space within the bounding box of the selected atoms
static int _porosity(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET));

    // Static/validation path
    if (!dst) {
        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_none();
            ctx->backchannel->value_range = (frange_t){0, 1}; // void ratio e in [0, 1]
        }
        return 1; // Produces a single scalar
    }

    ASSERT(!dst || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT));

    // Only support orthogonal cells for now
    if (ctx->cur_state->unitcell.flags & MD_UNITCELL_TRICLINIC) {
        EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "void_ratio: triclinic unit cells not supported yet");
        *(float*)dst->ptr = 0.0f;
        return 0;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);

    const md_bitset_t* bf = as_bitset(arg[0]);
    size_t bf_dim = element_count(arg[0]);

    md_bitset_t tmp = {0};
    if (bf_dim > 1) {
        tmp = _internal_flatten_bf(bf, bf_dim, ctx->temp_alloc);
        bf = &tmp;
    }

    const size_t count = md_bitset_count(*bf);
    float* out = (float*)dst->ptr;
    *out = 0.0f;

    if (count == 0) {
        EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "void_ratio: empty selection");
        goto done;
    }

    // Extract indices of selected atoms
    int* idx = md_temp_alloc_array(temp, int, count);
    md_bitset_iter_extract_indices(idx, count, md_bitset_iter(*bf));

    // Per atom radii, indexed by atom index as extract_xyzw_vec4 reads them. This used to be the radii of
    // atoms [0, count) - read by atom index, so past its end for most selections - and those values
    // set the extent of the grid through the aabb below.
    const float* rad = ctx->atom_radius;
    if (!rad) {
        float* all = md_temp_alloc_array(temp, float, ctx->sys->atom.count);
        md_atom_extract_radii(all, 0, ctx->sys->atom.count, &ctx->sys->atom);
        rad = all;
    }

    // Build a working copy of coordinates and deperiodize to make selection contiguous
    md_array(vec4_t) xyzr = 0;
    md_array_resize(xyzr, count, ctx->temp_alloc);
    extract_xyzw_vec4(xyzr, ctx->cur_state->xyz, rad, bf);
    md_util_deperiodize_self_vec4(xyzr, count, &ctx->cur_state->unitcell, NULL);

    vec3_t bmin, bmax;
    md_util_aabb_compute_vec4(bmin.elem, bmax.elem, xyzr, 0, count);

    // Compute voxel sizes from ROI extents; guard zero extents
    vec3_t ext = vec3_max(vec3_sub(bmax, bmin), vec3_set1(1.0));
    float max_ext = MAX(ext.x, MAX(ext.y, ext.z));

    const int max_dim = 512; // Use this for largest extent then try to construct a uniform voxel size
    const float target_voxel_ext = max_ext / max_dim;

    int dim[3] = {
        MAX(1, (int)(ext.x / target_voxel_ext)),
        MAX(1, (int)(ext.y / target_voxel_ext)),
        MAX(1, (int)(ext.z / target_voxel_ext)),
    };

    const int dim2 = dim[0] * dim[1];
    const size_t num_bits = (size_t)dim[0] * (size_t)dim[1] * (size_t)dim[2];
    const size_t num_u64 = (num_bits + 63) / 64;
    uint64_t* vol_bits = md_temp_alloc_zero_array(temp, uint64_t, num_u64);

    const float dx = ext.x / (float)dim[0];
    const float dy = ext.y / (float)dim[1];
    const float dz = ext.z / (float)dim[2];

    // Voxelize occupancy
    for (size_t i = 0; i < count; ++i) {
        const vec4_t p = xyzr[i];
        const float x = p.x, y = p.y, z = p.z;
        const float r  = md_atom_radius(&ctx->sys->atom, idx[i]);

        // Bounding voxel range for this sphere (in ROI grid coordinates)
        int ix_min = (int)floorf((x - r - bmin.x) / dx);
        int iy_min = (int)floorf((y - r - bmin.y) / dy);
        int iz_min = (int)floorf((z - r - bmin.z) / dz);
        int ix_max = (int)floorf((x + r - bmin.x) / dx);
        int iy_max = (int)floorf((y + r - bmin.y) / dy);
        int iz_max = (int)floorf((z + r - bmin.z) / dz);

        // Clamp to grid
        ix_min = MAX(0, MIN(dim[0] - 1, ix_min));
        iy_min = MAX(0, MIN(dim[1] - 1, iy_min));
        iz_min = MAX(0, MIN(dim[2] - 1, iz_min));
        ix_max = MAX(0, MIN(dim[0] - 1, ix_max));
        iy_max = MAX(0, MIN(dim[1] - 1, iy_max));
        iz_max = MAX(0, MIN(dim[2] - 1, iz_max));

        const float r2 = r * r;

        for (int iz = iz_min; iz <= iz_max; ++iz) {
            const float vz = bmin.z + ((float)iz + 0.5f) * dz;
            const float dzv = vz - z;
            for (int iy = iy_min; iy <= iy_max; ++iy) {
                const float vy = bmin.y + ((float)iy + 0.5f) * dy;
                const float dyv = vy - y;
                for (int ix = ix_min; ix <= ix_max; ++ix) {
                    const float vx = bmin.x + ((float)ix + 0.5f) * dx;
                    const float dxv = vx - x;
                    const float d2  = fmaf(dxv, dxv, fmaf(dyv, dyv, dzv * dzv));
                    if (d2 <= r2) {
                        const size_t bit_idx = (size_t)iz * dim2 + (size_t)iy * dim[0] + (size_t)ix;
                        vol_bits[bit_idx >> 6] |= ((uint64_t)1u << (bit_idx & 63));
                    }
                }
            }
        }
    }

    // Count occupied voxels
    size_t set_bits = 0;
    for (size_t i = 0; i < num_u64; ++i) set_bits += popcnt64(vol_bits[i]);

    // Compute void ratio e = V_void / V_solid
    if (set_bits == 0) {
        // No solid found in ROI — undefined e; report 0 and warn
        EVAL_LOG_ERROR(ctx, ctx->arg_tokens[0], "void_ratio: zero occupied volume in ROI");
        *out = 0.0f;
    } else {
        const double v_void  = (double)num_bits - set_bits;
        const double v_total = (double)num_bits;
        *out = (float)(v_void / v_total);
    }

done:
    md_temp_end(temp);
    return 0;
}

static int _shape_weights(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_COORDINATE_ARR));    
    ASSERT(arg[0].ptr);

    if (dst) {
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_FLOAT3_ARR));
        md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);

        vec3_t* out_weights = as_vec3_arr(*dst);
        // Bounded by the destination: evaluated inside a context, an argument whose length varies over frames
        // can hold more atomsets than the length the context fixed at compile time
        const size_t out_len = type_info_total_element_count(dst->type) / 3;

        md_array(vec4_t)  xyzw = 0;
        if (arg[0].type.base_type == TYPE_ATOMSET) {
            const size_t bf_len = MIN(element_count(arg[0]), out_len);
            const md_bitset_t* bf_arr = as_bitset(arg[0]);
            md_bitset_t tmp_bf = {0};

            for (size_t i = 0; i < bf_len; ++i) {
                const md_bitset_t* bf = &bf_arr[i];
                if (ctx->mol_ctx) {
                    tmp_bf = md_bitset_and(*bf, *ctx->mol_ctx, ctx->temp_alloc);
                    bf = &tmp_bf;
                }
                const size_t count = md_bitset_count(*bf);
                if (count > 0) {
                    md_array_resize(xyzw, count, ctx->temp_alloc);

                    extract_xyzw_vec4(xyzw, ctx->cur_state->xyz, ctx->atom_mass, bf);
                    // The covariance is about the mean of the placed points. The circular mean the
                    // placement is seeded from is not that mean, and the offset between the two enters
                    // the covariance squared.
                    vec3_t com = {0};
                    md_util_deperiodize_self_vec4(xyzw, count, &ctx->cur_state->unitcell, &com);

                    const mat3_t M = mat3_covariance_matrix_vec4(xyzw, 0, count, com);
                    out_weights[i] = md_util_shape_weights(&M);
                }
            }
        } else {
            xyzw = coordinate_extract_xyzw(arg[0], 1.0f, ctx);
            const size_t count = md_array_size(xyzw);
            vec3_t com = {0};
            md_util_deperiodize_self_vec4(xyzw, count, &ctx->cur_state->unitcell, &com);
            const mat3_t M = mat3_covariance_matrix_vec4(xyzw, 0, count, com);
            out_weights[0] = md_util_shape_weights(&M);
        }

        md_temp_end(temp);
    } else if (ctx->vis) {
        coordinate_visualize(arg[0], ctx);
    } else {
        int count = 0;
        if (arg[0].type.base_type == TYPE_ATOMSET) {
            count = (int)element_count(arg[0]);
        } else {
            count = 1;
        }

        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_none();
            ctx->backchannel->value_range = (frange_t){0, 1.0};
        }
        return count;
    }

    return 0;
}

// Comparison of floats. A scalar is broadcast against an array, and two arrays of equal length are compared
// element by element. The same implementation serves the array overloads, and when 'dst' is NULL it reports the
// length of the result (the length of the longest argument) so that dynamic arguments can be validated.
#define BAKE_CMP_ARR(name, op) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_FLOAT_ARR)); \
        ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_FLOAT_ARR)); \
        const float* a = as_float_arr(arg[0]); \
        const float* b = as_float_arr(arg[1]); \
        int a_len = type_info_array_len(arg[0].type); \
        int b_len = type_info_array_len(arg[1].type); \
        if (dst) { \
            ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_BOOL_ARR)); \
            int out_len = type_info_array_len(dst->type); \
            bool* out = as_bool_arr(*dst); \
            if (a_len < b_len) { \
                ASSERT(a_len == 1); \
                for (int i = 0; i < out_len; ++i) out[i] = a[0] op b[i]; \
            } else if (a_len > b_len) { \
                ASSERT(b_len == 1); \
                ASSERT(out_len == a_len); \
                for (int i = 0; i < out_len; ++i) out[i] = a[i] op b[0]; \
            } else { \
                ASSERT(a_len == out_len && b_len == out_len); \
                for (int i = 0; i < a_len; ++i) out[i] = a[i] op b[i]; \
            } \
            return 0; \
        } \
        return MAX(a_len, b_len); \
    }

#define BAKE_CMP_S_S(name, op) \
    static int name(data_t* dst, data_t arg[], eval_context_t* ctx) { \
        (void)ctx; \
        if (!dst) return 0; \
        *((bool*)dst->ptr) = *((float*)arg[0].ptr) op *((float*)arg[1].ptr); \
        return 0; \
    }

BAKE_CMP_ARR(_cmp_lt, <)
BAKE_CMP_ARR(_cmp_gt, >)
BAKE_CMP_ARR(_cmp_le, <=)
BAKE_CMP_ARR(_cmp_ge, >=)
BAKE_CMP_ARR(_cmp_eq, ==)

BAKE_CMP_S_S(_cmp_lt_f_f, <)
BAKE_CMP_S_S(_cmp_gt_f_f, >)
BAKE_CMP_S_S(_cmp_le_f_f, <=)
BAKE_CMP_S_S(_cmp_ge_f_f, >=)
BAKE_CMP_S_S(_cmp_eq_f_f, ==)

/*
// This is some experimental future work, for matching structures using maximum common subgraph

static int32_t greedy_align(const md_bitset_t* a, const md_bitset_t* b, const md_system_t* sys) {
    const int64_t a_size = md_bitset_count(*a);
    const int64_t b_size = md_bitset_count(*b);
    const md_bitset_t* small;
    const md_bitset_t* big;

    if (a_size < b_size) {
        small = a;
        big = b;
    } else {
        small = b;
        big = a;
    }
}

static int _align(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(is_type_directly_compatible(arg[0].type, (type_info_t)TI_ATOMSET));
    ASSERT(is_type_directly_compatible(arg[1].type, (type_info_t)TI_ATOMSET));

    if (dst) {
        ASSERT(dst->ptr);
        ASSERT(!dst || dst->size == 0 || is_type_directly_compatible(dst->type, (type_info_t)TI_INT_ARR));
    }
}
*/

// ### CONTACTS ###
// contacts() is the script face of md_contact: a set of group pairs in contact per frame, an opaque value of
// TYPE_CONTACT which the procedures below reduce. Its groups cannot depend on the frame, so everything which
// depends on the topology is prepared once per call site, during the static check (contact_query_entry_t).

// The element of arr which holds each particle. Particles in no element get a label of their own (num + index),
// so they never share one. The first element holding a particle wins.
static uint32_t* contacts_label_particles(const md_bitset_t* arr, size_t num, uint32_t num_atoms, md_allocator_i* alloc) {
    uint32_t* label = md_alloc(alloc, sizeof(uint32_t) * MAX(num_atoms, 1));
    for (uint32_t k = 0; k < num_atoms; ++k) label[k] = (uint32_t)num + k;
    for (size_t e = num; e-- > 0;) {
        md_bitset_iter_t it = md_bitset_iter(arr[e]);
        while (md_bitset_iter_next(&it)) {
            const uint64_t k = it.idx;
            if (k < num_atoms) label[k] = (uint32_t)e;
        }
    }
    return label;
}

// Argument order: a, b, cutoff, exclude_bonds, min_separation, parent, exclude_within
static bool contacts_prepare(md_contact_query_t* q, const data_t arg[], const eval_context_t* ctx, md_allocator_i* alloc) {
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    const uint32_t N = (uint32_t)ctx->sys->atom.count;

    md_contact_desc_t desc = {
        .group_a = as_bitset(arg[0]),
        .num_a   = element_count(arg[0]),
        .cutoff  = as_float(arg[2]),
        .exclude_bonds  = (uint32_t)as_int(arg[3]),
        .min_separation = (uint32_t)as_int(arg[4]),
    };
    if (!is_absent_type(arg[1].type)) {
        // An empty b is an empty set, which is not the same as no b: that would be the contacts within a
        static const md_bitset_t no_groups = {0};
        desc.num_b   = element_count(arg[1]);
        desc.group_b = desc.num_b ? as_bitset(arg[1]) : &no_groups;
    }
    if (!is_absent_type(arg[5].type) && desc.num_a) {
        // The parent of a group is the parent of its first particle
        const uint32_t* label = contacts_label_particles(as_bitset(arg[5]), element_count(arg[5]), N, temp_alloc);
        uint32_t* parent = md_alloc(temp_alloc, sizeof(uint32_t) * desc.num_a);
        for (size_t g = 0; g < desc.num_a; ++g) {
            md_bitset_iter_t it = md_bitset_iter(desc.group_a[g]);
            parent[g] = (md_bitset_iter_next(&it) && it.idx < N) ? label[it.idx] : UINT32_MAX - (uint32_t)g;
        }
        desc.group_parent = parent;
    }
    if (!is_absent_type(arg[6].type)) {
        desc.particle_label = contacts_label_particles(as_bitset(arg[6]), element_count(arg[6]), N, temp_alloc);
    }

    const bool result = md_contact_query_init(q, &desc, ctx->sys, alloc);
    md_temp_end(temp);
    return result;
}

static const md_contact_query_t* contacts_find_query(const eval_context_t* ctx) {
    if (!ctx->op_node || !ctx->ir) return NULL;
    for (size_t i = 0; i < md_array_size(ctx->ir->contact_queries); ++i) {
        const contact_query_entry_t* e = &ctx->ir->contact_queries[i];
        if (e->node == ctx->op_node && e->sys == ctx->sys) {
            return &e->query;
        }
    }
    return NULL;
}

static int _contacts(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);

    if (!dst && !ctx->vis) {
        // Static validation: everything a query is prepared from must be the same for every frame
        static const char* names[] = { "a", "b", "cutoff", "exclude_bonds", "min_separation", "parent", "exclude_within" };
        for (int i = 0; i < 7; ++i) {
            if (!is_absent_type(arg[i].type) && (ctx->arg_flags[i] & FLAG_DYNAMIC)) {
                EVAL_LOG_ERROR(ctx, ctx->arg_tokens[i], "contacts: '%s' cannot depend on the frame, the groups and their parameters have to be fixed", names[i]);
                return STATIC_VALIDATION_ERROR;
            }
        }
        if (!(as_float(arg[2]) > 0.0f)) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[2], "contacts: the cutoff has to be positive");
            return STATIC_VALIDATION_ERROR;
        }
        if (as_int(arg[3]) < 0 || as_int(arg[4]) < 0) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "contacts: exclude_bonds and min_separation cannot be negative");
            return STATIC_VALIDATION_ERROR;
        }
        // A plain selection is one group. Legal, but rarely what is meant, and the result looks broken
        const size_t num_a = element_count(arg[0]);
        const size_t num_b = is_absent_type(arg[1].type) ? 0 : element_count(arg[1]);
        if (is_absent_type(arg[1].type) && num_a <= 1) {
            EVAL_LOG_WARNING(ctx, ctx->arg_tokens[0], "contacts: 'a' is a single group, and a group is never in contact with itself, so there are no contacts. "
                "Give the groups as an array, e.g. residue(...), or a second set as 'b'");
        } else if (!is_absent_type(arg[1].type) && num_a == 1 && num_b == 1) {
            EVAL_LOG_WARNING(ctx, ctx->op_token, "contacts: 'a' and 'b' are single groups, so at most one pair of groups is in contact and count(c) is 0 or 1. "
                "count(c, 'atom') counts the particle pairs; groups are given as arrays, e.g. residue(...)");
        }
        if (ctx->backchannel) {
            ctx->backchannel->unit[1] = md_unit_none();
        }
        // Prepare the call site, once, when compiling. An evaluation does not write to the IR: one which finds
        // no prepared query makes one for its call only (below).
        if (ctx->compile_ir && ctx->op_node && !contacts_find_query(ctx)) {
            contact_query_entry_t entry = { .node = ctx->op_node, .sys = ctx->sys };
            if (!contacts_prepare(&entry.query, arg, ctx, ctx->compile_ir->arena)) {
                EVAL_LOG_ERROR(ctx, ctx->op_token, "contacts: failed to prepare the query");
                return STATIC_VALIDATION_ERROR;
            }
            md_array_push(ctx->compile_ir->contact_queries, entry, ctx->compile_ir->arena);
        }
        return 0;
    }

    const md_bitset_t* A = as_bitset(arg[0]);
    const md_bitset_t* B = is_absent_type(arg[1].type) ? A : as_bitset(arg[1]);

    if (!dst) {
        // Visualization only. Finding the pairs means a neighbour search over every particle of the groups, which in a
        // large system stalls whatever asked for a picture, so it is never done for one: the groups are shown instead.
        for (size_t g = 0; g < element_count(arg[0]); ++g) visualize_atom_mask(&A[g], ctx);
        if (B != A) {
            for (size_t g = 0; g < element_count(arg[1]); ++g) visualize_atom_mask(&B[g], ctx);
        }
        return 0;
    }

    if (ctx->backchannel) {
        // Evaluated within the static check, for a procedure which is asking about its argument (degree asks for the
        // number of groups). The shape is known without searching, and a search here would stall the compilation.
        md_contact_set_t* set = (md_contact_set_t*)dst->ptr;
        MEMSET(set, 0, sizeof(md_contact_set_t));
        set->num_a = (uint32_t)element_count(arg[0]);
        set->num_b = is_absent_type(arg[1].type) ? set->num_a : (uint32_t)element_count(arg[1]);
        set->flags = is_absent_type(arg[1].type) ? MD_CONTACT_FLAG_SELF : MD_CONTACT_FLAG_NONE;
        return 0;
    }

    int result = 0;
    md_temp_scope_t temp = md_temp_begin_avoid(ctx->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    const md_contact_query_t* q = contacts_find_query(ctx);
    md_contact_query_t local = {0};
    if (!q) {
        // Not prepared by a static check of this call site, so for this call only
#ifdef MD_SCRIPT_TEST_HOOKS
        test_hook_contact_local_prepares += 1;
#endif
        if (!contacts_prepare(&local, arg, ctx, temp_alloc)) {
            result = -1;
            goto done;
        }
        q = &local;
    }

    // The pairs live for the evaluation
    md_contact_set_t* set = (md_contact_set_t*)dst->ptr;
#ifdef MD_SCRIPT_TEST_HOOKS
    test_hook_contact_searches += 1;
#endif
    if (!md_contact_query_eval(set, q, ctx->cur_state, ctx->alloc)) {
        result = -1;
        goto done;
    }

    if (ctx->vis) {
        // Computed anyway: show the groups actually in contact
        for (size_t k = 0; k < set->count; ++k) {
            visualize_atom_mask(&A[set->i[k]], ctx);
            visualize_atom_mask(&B[set->j[k]], ctx);
        }
    }

done:
    md_temp_end(temp);
    return result;
}

static inline const md_contact_set_t* as_contact_set(data_t arg) {
    return (const md_contact_set_t*)arg.ptr;
}

// count(c): the number of group pairs in contact
static int _count_contact(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)ctx;
    if (dst) {
        const md_contact_set_t* set = as_contact_set(arg[0]);
        as_float(*dst) = set ? (float)set->count : 0.0f;
    }
    return 0;
}

// count(c, unit): 'group' counts group pairs, 'atom' the particle pairs behind them
static int _count_contact_unit(data_t* dst, data_t arg[], eval_context_t* ctx) {
    const str_t unit = as_string(arg[1]);
    const bool atom  = str_eq(unit, STR_LIT("atom"));
    const bool group = str_eq(unit, STR_LIT("group"));
    if (!dst) {
        if (!atom && !group) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[1], "count: the unit of a contact set is 'group' (pairs of groups) or 'atom' (pairs of particles)");
            return STATIC_VALIDATION_ERROR;
        }
        return 0;
    }
    const md_contact_set_t* set = as_contact_set(arg[0]);
    double count = 0.0;
    if (set) {
        if (atom) {
            for (size_t k = 0; k < set->count; ++k) count += set->atom_pairs[k];
        } else {
            count = (double)set->count;
        }
    }
    as_float(*dst) = (float)count;
    return 0;
}

// degree(c): per group of A, the number of groups it is in contact with
static int _degree_contact(data_t* dst, data_t arg[], eval_context_t* ctx) {
    (void)ctx;
    const md_contact_set_t* set = as_contact_set(arg[0]);
    if (!dst) {
        // The length is the number of groups of A, which is fixed
        return set ? (int)set->num_a : 0;
    }
    float* out = as_float_arr(*dst);
    const size_t len = element_count(*dst);
    MEMSET(out, 0, len * sizeof(float));
    if (set) {
        const bool self = set->flags & MD_CONTACT_FLAG_SELF;
        for (size_t k = 0; k < set->count; ++k) {
            if (set->i[k] < len) out[set->i[k]] += 1.0f;
            if (self && set->j[k] < len) out[set->j[k]] += 1.0f;
        }
    }
    return 0;
}

// chunks(sel, size): each element of sel cut into consecutive runs of size particles, in index order. A last run
// shorter than size is kept. Groups which the topology does not name: the slices of a fibril which is one residue.
static int _chunks(data_t* dst, data_t arg[], eval_context_t* ctx) {
    const md_bitset_t* src = as_bitset(arg[0]);
    const size_t num_src = element_count(arg[0]);
    const int size = as_int(arg[1]);
    if (size <= 0) {
        EVAL_LOG_ERROR(ctx, ctx->arg_tokens[1], "chunks: the size has to be positive");
        return STATIC_VALIDATION_ERROR;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    md_bitset_t* dst_bf = dst ? as_bitset(*dst) : NULL;
    const size_t cap = dst ? element_count(*dst) : 0;
    md_allocator_i* dst_alloc = dst ? data_bitset_alloc(dst, ctx) : NULL;
    // The last chunk also takes everything beyond the capacity (a flattened evaluation), so it is built
    md_bitmask_t last = md_bitmask_create(0, atom_domain_end(ctx), temp_alloc);
    size_t count = 0;
    md_array(int32_t) idx = 0;
    for (size_t e = 0; e < num_src; ++e) {
        md_bitset_t tmp = {0};
        const md_bitset_t* bf = &src[e];
        if (ctx->mol_ctx) {
            tmp = md_bitset_and(*bf, *ctx->mol_ctx, temp_alloc);
            bf = &tmp;
        }
        const size_t pop = md_bitset_count(*bf);
        if (!pop) continue;
        if (!dst) {
            count += DIV_UP(pop, (size_t)size);
            continue;
        }
        md_array_resize(idx, pop, temp_alloc);
        md_bitset_iter_extract_indices(idx, pop, md_bitset_iter(*bf));
        for (size_t k0 = 0; k0 < pop; k0 += (size_t)size) {
            const size_t k1 = MIN(k0 + (size_t)size, pop);
            const size_t c = count + k0 / (size_t)size;
            if (c + 1 < cap) {
                dst_bf[c] = md_bitset_from_indices((const uint32_t*)idx + k0, k1 - k0, dst_alloc);
            } else {
                // Beyond the capacity (a flattened evaluation), everything goes into the last one
                md_bitmask_set_indices(&last, (const uint32_t*)idx + k0, k1 - k0);
            }
        }
        count += DIV_UP(pop, (size_t)size);
    }
    if (dst && cap > 0 && !md_bitmask_empty(&last)) {
        dst_bf[cap - 1] = md_bitmask_extract(&last, dst_alloc);
    }
    md_temp_end(temp);

    if (!dst) {
        if (ctx->backchannel && (ctx->arg_flags[0] & FLAG_DYNAMIC)) {
            ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
        }
        if (ctx->eval_flags & EVAL_FLAG_FLATTEN) {
            count = MIN(1, count);
        }
        return (int)count;
    }
    return 0;
}

// ### BONDS ###
// The bonds as values of their own (bondset_t): a set of bonds is not a set of atoms, which brings in every bond among
// its atoms. bond() selects persistent bonds by their atoms and flags, bond_delta() the bonds formed and broken in a
// frame, and the procedures below turn them into atoms and numbers.

static const match_name_t bond_resonances[] = {
    {STR_INIT("aromatic"),    MD_BOND_FLAG_AROMATIC},
    {STR_INIT("delocalized"), MD_BOND_FLAG_DELOCALIZED},
    {STR_INIT("localized"),   MD_BOND_FLAG_NONE},
};

// The persistent bonds of the system: the graph, then the non-covalent bonds (md_bond_data_t)
static inline uint32_t bond_total(const md_system_t* sys) {
    return (uint32_t)(sys->bond.count + sys->bond.noncovalent_count);
}

static inline uint64_t bond_pair_key(md_atom_pair_t p) {
    return ((uint64_t)(uint32_t)MIN(p.idx[0], p.idx[1]) << 32) | (uint32_t)MAX(p.idx[0], p.idx[1]);
}

static inline md_atom_pair_t bond_key_pair(uint64_t key) {
    const md_atom_pair_t p = {{ (md_atom_idx_t)(key >> 32), (md_atom_idx_t)(key & 0xFFFFFFFFu) }};
    return p;
}

// Both atoms of a pair within the context, if there is one
static inline bool bond_in_context(md_atom_pair_t p, const eval_context_t* ctx) {
    return !ctx->mol_ctx || (md_bitset_test(*ctx->mol_ctx, (uint32_t)p.idx[0]) && md_bitset_test(*ctx->mol_ctx, (uint32_t)p.idx[1]));
}

// The extra pairs of a set being built, kept as keys with their flags
typedef struct bond_extras_t {
    md_array(uint64_t)        key;
    md_array(md_bond_flags_t) flags;
} bond_extras_t;

static int compare_bond_key(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

typedef struct bond_keyed_t {
    uint64_t key;       // First: compare_bond_key reads it
    uint64_t idx;
} bond_keyed_t;

// The extras of a set, sorted and without duplicates (the first of equal keys keeps its flags)
static void bond_extras_tidy(bond_extras_t* e, md_allocator_i* temp) {
    const size_t n = md_array_size(e->key);
    if (n < 2) return;
    bond_keyed_t* k = md_alloc(temp, sizeof(bond_keyed_t) * n);
    for (size_t i = 0; i < n; ++i) {
        k[i].key = e->key[i];
        k[i].idx = i;
    }
    qsort(k, n, sizeof(bond_keyed_t), compare_bond_key);
    md_bond_flags_t* flags = md_alloc(temp, sizeof(md_bond_flags_t) * n);
    MEMCPY(flags, e->flags, sizeof(md_bond_flags_t) * n);
    size_t w = 0;
    for (size_t i = 0; i < n; ++i) {
        if (w > 0 && e->key[w - 1] == k[i].key) continue;
        e->key[w]   = k[i].key;
        e->flags[w] = flags[k[i].idx];
        w += 1;
    }
    md_array_shrink(e->key, w);
    md_array_shrink(e->flags, w);
}

static bondset_t bondset_make(const md_bitmask_t* bonds, const bond_extras_t* extras, md_allocator_i* alloc) {
    bondset_t r = { .bonds = md_bitmask_extract(bonds, alloc) };
    const size_t n = extras ? md_array_size(extras->key) : 0;
    if (n) {
        md_atom_pair_t*  pair  = md_alloc(alloc, sizeof(md_atom_pair_t)  * n);
        md_bond_flags_t* flags = md_alloc(alloc, sizeof(md_bond_flags_t) * n);
        for (size_t i = 0; i < n; ++i) {
            pair[i]  = bond_key_pair(extras->key[i]);
            flags[i] = extras->flags[i];
        }
        r.num_extra   = (uint32_t)n;
        r.extra_pair  = pair;
        r.extra_flags = flags;
    }
    return r;
}

static void bond_extras_push(bond_extras_t* e, md_atom_pair_t pair, md_bond_flags_t flags, md_allocator_i* temp) {
    md_array_push(e->key, bond_pair_key(pair), temp);
    md_array_push(e->flags, flags, temp);
}

// The union of the sets of an argument: its bonds into b, its extras into e
static void bondset_union_into(md_bitmask_t* b, bond_extras_t* e, const bondset_t* sets, size_t count, md_allocator_i* temp) {
    for (size_t i = 0; i < count; ++i) {
        md_bitmask_or(b, sets[i].bonds);
        for (uint32_t k = 0; k < sets[i].num_extra; ++k) {
            bond_extras_push(e, sets[i].extra_pair[k], sets[i].extra_flags[k], temp);
        }
    }
    bond_extras_tidy(e, temp);
}

// The bonds of a set within the context (both atoms in it)
static void bondset_restrict_to_context(md_bitmask_t* b, bond_extras_t* e, const eval_context_t* ctx) {
    if (!ctx->mol_ctx) return;
    md_bitset_t set = md_bitmask_extract(b, ctx->temp_alloc);
    md_bitset_iter_t it = md_bitset_iter(set);
    while (md_bitset_iter_next(&it)) {
        const md_atom_pair_t p = ctx->sys->bond.pairs[it.idx];
        if (!bond_in_context(p, ctx)) md_bitmask_clear_bit(b, (uint32_t)it.idx);
    }
    size_t w = 0;
    for (size_t k = 0; k < md_array_size(e->key); ++k) {
        if (!bond_in_context(bond_key_pair(e->key[k]), ctx)) continue;
        e->key[w] = e->key[k];
        e->flags[w] = e->flags[k];
        w += 1;
    }
    md_array_shrink(e->key, w);
    md_array_shrink(e->flags, w);
}

static void visualize_bond(md_atom_pair_t p, bool atoms, bool lines, eval_context_t* ctx) {
    if (atoms) {
        md_bitmask_set_bit(ctx->vis_atom_mask, (uint32_t)p.idx[0]);
        md_bitmask_set_bit(ctx->vis_atom_mask, (uint32_t)p.idx[1]);
    }
    if (lines) {
        const vec3_t* xyz = ctx->cur_state->xyz;
        push_line(vertex(xyz[p.idx[0]], COLOR_WHITE), vertex(xyz[p.idx[1]], COLOR_WHITE), ctx->vis);
    }
}

// The atoms of the bonds of a set, and a line along each
static void visualize_bondset(const bondset_t* b, eval_context_t* ctx) {
    const md_system_t* sys = ctx->sys;
    const bool atoms = ctx->vis_atom_mask && (ctx->vis_flags == MD_SCRIPT_VISUALIZE_DEFAULT || (ctx->vis_flags & MD_SCRIPT_VISUALIZE_ATOMS));
    const bool lines = ctx->vis && ctx->cur_state && ctx->cur_state->xyz && (ctx->vis_flags == MD_SCRIPT_VISUALIZE_DEFAULT || (ctx->vis_flags & MD_SCRIPT_VISUALIZE_GEOMETRY));
    const uint32_t total = bond_total(sys);
    md_bitset_iter_t it = md_bitset_iter(b->bonds);
    while (md_bitset_iter_next(&it)) {
        if (it.idx >= total) break;
        visualize_bond(sys->bond.pairs[it.idx], atoms, lines, ctx);
    }
    for (uint32_t k = 0; k < b->num_extra; ++k) {
        visualize_bond(b->extra_pair[k], atoms, lines, ctx);
    }
}

// The filters of bond() (kind, order, resonance, origin): which flags pass
typedef struct bond_filter_t {
    md_bond_kind_name_t kinds[16];
    size_t   num_kinds;         // 0: any kind
    bool     has_order;         // Bonds of the graph whose order is in [order_lo, order_hi]
    int      order_lo, order_hi;
    uint32_t resonance;         // Bits of bond_resonances, 0: any. Bonds of the graph only
    uint32_t origin;            // Bits of 1 << md_bond_origin_t, 0: any
} bond_filter_t;

static bool bond_origin_lookup(int* out, str_t name) {
    for (int o = MD_BOND_ORIGIN_FILE; o <= MD_BOND_ORIGIN_INFERRED; ++o) {
        if (str_eq_cstr(name, md_bond_origin_name((md_bond_origin_t)o))) {
            *out = o;
            return true;
        }
    }
    return false;
}

static void bond_report_unknown_origin(eval_context_t* ctx, int arg, str_t name) {
    match_name_t names[4];
    for (int o = 0; o < 4; ++o) names[o] = (match_name_t){ str_from_cstr(md_bond_origin_name((md_bond_origin_t)o)), o };
    match_report_unknown(ctx, arg, "bond origin", names, ARRAY_SIZE(names), name);
}

// The origins of a string argument as a mask, or false (reported when validating)
static bool bond_origin_mask(uint32_t* out, data_t arg, int arg_idx, bool report, eval_context_t* ctx) {
    *out = 0;
    const str_t* names = as_string_arr(arg);
    for (size_t i = 0; i < element_count(arg); ++i) {
        int o = 0;
        if (!bond_origin_lookup(&o, names[i])) {
            if (report) bond_report_unknown_origin(ctx, arg_idx, names[i]);
            return false;
        }
        *out |= 1u << o;
    }
    return true;
}

static bool bond_filter_parse(bond_filter_t* f, const data_t arg[], bool report, eval_context_t* ctx) {
    MEMSET(f, 0, sizeof(bond_filter_t));
    if (!is_absent_type(arg[2].type)) {
        size_t count = 0;
        const md_bond_kind_name_t* table = md_bond_kind_names(&count);
        const str_t* names = as_string_arr(arg[2]);
        for (size_t i = 0; i < element_count(arg[2]); ++i) {
            size_t k = 0;
            while (k < count && !str_eq_cstr(names[i], table[k].name)) ++k;
            if (k == count) {
                if (report) {
                    match_name_t valid[16];
                    for (size_t j = 0; j < count && j < ARRAY_SIZE(valid); ++j) valid[j] = (match_name_t){ str_from_cstr(table[j].name), (int)j };
                    match_report_unknown(ctx, 2, "bond kind", valid, MIN(count, ARRAY_SIZE(valid)), names[i]);
                }
                return false;
            }
            if (f->num_kinds < ARRAY_SIZE(f->kinds)) f->kinds[f->num_kinds++] = table[k];
        }
    }
    if (!is_absent_type(arg[3].type)) {
        const irange_t r = as_irange(arg[3]);
        if (r.end < r.beg || r.beg < 0) {
            if (report) EVAL_LOG_ERROR(ctx, ctx->arg_tokens[3], "bond: the order is a number from 0 (unknown) to %d, or a range of them", MD_BOND_ORDER_QUADRUPLE);
            return false;
        }
        f->has_order = true;
        f->order_lo = r.beg;
        f->order_hi = r.end;
    }
    if (!is_absent_type(arg[4].type)) {
        const str_t* names = as_string_arr(arg[4]);
        for (size_t i = 0; i < element_count(arg[4]); ++i) {
            size_t k = 0;
            while (k < ARRAY_SIZE(bond_resonances) && !str_eq(names[i], bond_resonances[k].name)) ++k;
            if (k == ARRAY_SIZE(bond_resonances)) {
                if (report) match_report_unknown(ctx, 4, "bond resonance", bond_resonances, ARRAY_SIZE(bond_resonances), names[i]);
                return false;
            }
            f->resonance |= 1u << k;
        }
    }
    if (!is_absent_type(arg[5].type) && !bond_origin_mask(&f->origin, arg[5], 5, report, ctx)) return false;
    return true;
}

static bool bond_filter_pass(const bond_filter_t* f, md_bond_flags_t fl) {
    if (f->num_kinds) {
        bool any = false;
        for (size_t k = 0; k < f->num_kinds && !any; ++k) any = (fl & f->kinds[k].mask) == f->kinds[k].value;
        if (!any) return false;
    }
    // Order and resonance are of the bonds of the graph: neither a non-covalent bond, a spring (elastic) nor the bond of
    // a metal (metallic) has them
    const bool graph = md_bond_is_graph(fl);
    if (f->has_order) {
        const int o = md_bond_order(fl);
        if (!graph || o < f->order_lo || o > f->order_hi) return false;
    }
    if (f->resonance) {
        if (!graph) return false;
        const uint32_t r = fl & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED);
        bool any = false;
        for (size_t k = 0; k < ARRAY_SIZE(bond_resonances) && !any; ++k) any = (f->resonance & (1u << k)) && r == (uint32_t)bond_resonances[k].value;
        if (!any) return false;
    }
    if (f->origin && !(f->origin & (1u << md_bond_origin(fl)))) return false;
    return true;
}

// Adds the bonds of an atom which pass: the other atom in other (NULL: any), within the context
static void bond_collect_atom(md_bitmask_t* b, uint32_t atom, const md_bitset_t* group, const md_bitset_t* other, const bond_filter_t* f, const eval_context_t* ctx) {
    const md_system_t* sys = ctx->sys;
    md_bond_iter_t it = md_bond_iter(&sys->bond, atom);
    while (md_bond_iter_has_next(&it)) {
        const uint32_t o = (uint32_t)md_bond_iter_atom_index(&it);
        const uint32_t bi = (uint32_t)md_bond_iter_bond_index(&it);
        const md_bond_flags_t fl = (md_bond_flags_t)md_bond_iter_bond_flags(&it);
        md_bond_iter_next(&it);
        if (other) {
            if (!md_bitset_test(*other, o)) continue;
        } else if (group && !md_bitset_test(*group, o)) {
            continue;
        }
        const md_atom_pair_t p = {{ (md_atom_idx_t)atom, (md_atom_idx_t)o }};
        if (!bond_in_context(p, ctx) || !bond_filter_pass(f, fl)) continue;
        md_bitmask_set_bit(b, bi);
    }
}

// bond(a = none, b = none, *, kind = none, order = none, resonance = none, origin = none): persistent bonds of the
// system (both of its ranges, md_bond_data_t), within the context (both atoms in it), which all the filters let through.
//   a, b       Without either: all bonds. With a: the bonds with both atoms in a, one set per element of a (bond(residue(:))
//              is the bonds within each residue). With b as well: the bonds with one atom in the element of a and the other
//              in b. With b alone: the bonds with an atom in b.
//   kind       Any of md_bond_kind_names: 'covalent', 'coordinate', 'elastic', 'metallic', 'noncovalent', 'interaction',
//              'hydrogen'.
//   order      The order (an int or a range, 0 for unknown) of a bond of the graph: non-covalent bonds have none.
//   resonance  Any of 'aromatic', 'delocalized', 'localized' (neither), of a bond of the graph.
//   origin     Any of 'file', 'topology', 'user', 'inferred'.
// bond(d, ...) with a bond_delta in place of a is _bond_of_delta.
static int _bond(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    const bool report = !dst && !ctx->vis;
    bond_filter_t filter;
    if (!bond_filter_parse(&filter, arg, report, ctx)) return STATIC_VALIDATION_ERROR;

    const bool has_a = !is_absent_type(arg[0].type);
    const bool has_b = !is_absent_type(arg[1].type);
    const int num_groups = has_a ? (int)element_count(arg[0]) : 1;
    if (!dst) {
        if (ctx->vis) return 0;
        return (ctx->eval_flags & EVAL_FLAG_FLATTEN) ? MIN(1, num_groups) : num_groups;
    }
    ASSERT(is_type_directly_compatible(dst->type, (type_info_t)TI_BONDSET_ARR));
    const int64_t cap = type_info_array_len(dst->type);
    if (cap <= 0 || !dst->ptr) return 0;

    const md_system_t* sys = ctx->sys;
    const uint32_t total = bond_total(sys);
    md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_allocator_i* alloc = data_bitset_alloc(dst, ctx);

    // b as one set
    md_bitset_t B = {0};
    if (has_b) B = _internal_flatten_bf(as_bitset(arg[1]), element_count(arg[1]), temp_alloc);

    for (int64_t slot = 0; slot < cap; ++slot) {
        md_bitmask_t b = md_bitmask_create(0, total, temp_alloc);
        // The groups of this slot: one each, or all of them when there is one slot for many (flattened)
        const int g_beg = cap == 1 ? 0 : (int)slot;
        const int g_end = cap == 1 ? num_groups : MIN((int)slot + 1, num_groups);
        for (int g = g_beg; g < g_end; ++g) {
            const md_bitset_t* A = has_a ? &as_bitset(arg[0])[g] : NULL;
            if (A) {
                md_bitset_iter_t it = md_bitset_iter(*A);
                while (md_bitset_iter_next(&it)) {
                    if (it.idx >= sys->atom.count) break;
                    bond_collect_atom(&b, (uint32_t)it.idx, A, has_b ? &B : NULL, &filter, ctx);
                }
            } else if (has_b) {
                md_bitset_iter_t it = md_bitset_iter(B);
                while (md_bitset_iter_next(&it)) {
                    if (it.idx >= sys->atom.count) break;
                    bond_collect_atom(&b, (uint32_t)it.idx, NULL, NULL, &filter, ctx);
                }
            } else {
                for (uint32_t i = 0; i < sys->bond.count; ++i) {
                    if (bond_in_context(sys->bond.pairs[i], ctx) && bond_filter_pass(&filter, sys->bond.flags[i])) md_bitmask_set_bit(&b, i);
                }
            }
            // The non-covalent bonds have no connectivity: they are few, and looked at one by one
            for (uint32_t i = (uint32_t)sys->bond.count; i < total; ++i) {
                const md_atom_pair_t p = sys->bond.pairs[i];
                const uint32_t x = (uint32_t)p.idx[0], y = (uint32_t)p.idx[1];
                bool in = true;
                if (A && has_b)   in = (md_bitset_test(*A, x) && md_bitset_test(B, y)) || (md_bitset_test(*A, y) && md_bitset_test(B, x));
                else if (A)       in = md_bitset_test(*A, x) && md_bitset_test(*A, y);
                else if (has_b)   in = md_bitset_test(B, x) || md_bitset_test(B, y);
                if (in && bond_in_context(p, ctx) && bond_filter_pass(&filter, sys->bond.flags[i])) md_bitmask_set_bit(&b, i);
            }
        }
        as_bondset(*dst)[slot] = bondset_make(&b, NULL, alloc);
    }
    md_temp_end(temp);
    return 0;
}

// The atoms of the bonds of a set, both of each, into b
static void bondset_atoms_into(md_bitmask_t* b, const bondset_t* set, const md_system_t* sys) {
    const uint32_t total = bond_total(sys);
    md_bitset_iter_t it = md_bitset_iter(set->bonds);
    while (md_bitset_iter_next(&it)) {
        if (it.idx >= total) break;
        md_bitmask_set_bit(b, (uint32_t)sys->bond.pairs[it.idx].idx[0]);
        md_bitmask_set_bit(b, (uint32_t)sys->bond.pairs[it.idx].idx[1]);
    }
    for (uint32_t k = 0; k < set->num_extra; ++k) {
        md_bitmask_set_bit(b, (uint32_t)set->extra_pair[k].idx[0]);
        md_bitmask_set_bit(b, (uint32_t)set->extra_pair[k].idx[1]);
    }
}

// atom(b): the atoms of the bonds of each set, as residue(sel) grows a selection to its residues
static int _bond_atoms(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return (int)element_count(arg[0]);
    const bondset_t* src = as_bondset(arg[0]);
    const int64_t len = MIN((int64_t)element_count(arg[0]), type_info_array_len(dst->type));
    for (int64_t i = 0; i < len; ++i) {
        md_bitmask_t b = result_mask(ctx);
        bondset_atoms_into(&b, &src[i], ctx->sys);
        finish_result(dst, (size_t)i, &b, ctx);
    }
    return 0;
}

// The bonds of the union of the sets of an argument, in their order: persistent by index, then the extras
static size_t bondset_union_pairs(md_atom_pair_t** out, const bondset_t* sets, size_t count, const eval_context_t* ctx, md_allocator_i* temp) {
    const md_system_t* sys = ctx->sys;
    md_bitmask_t b = md_bitmask_create(0, bond_total(sys), temp);
    bond_extras_t e = {0};
    bondset_union_into(&b, &e, sets, count, temp);
    const size_t num_bonds = (size_t)md_bitmask_count(&b);
    const size_t n = num_bonds + md_array_size(e.key);
    md_atom_pair_t* pairs = md_alloc(temp, sizeof(md_atom_pair_t) * MAX(n, 1));
    size_t w = 0;
    md_bitset_iter_t it = md_bitmask_iter(&b);
    while (md_bitset_iter_next(&it) && w < num_bonds) pairs[w++] = sys->bond.pairs[it.idx];
    for (size_t k = 0; k < md_array_size(e.key); ++k) pairs[w++] = bond_key_pair(e.key[k]);
    *out = pairs;
    return w;
}

// count(b): the number of bonds in the union of the sets
static int _count_bonds(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);
    md_atom_pair_t* pairs = NULL;
    as_float(*dst) = (float)bondset_union_pairs(&pairs, as_bondset(arg[0]), element_count(arg[0]), ctx, md_temp_allocator(temp));
    md_temp_end(temp);
    return 0;
}

// length(b): the length of each bond of the union of the sets in the frame (minimum image), in the order of the set
static int _bond_length(data_t* dst, data_t arg[], eval_context_t* ctx) {
    md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);
    md_atom_pair_t* pairs = NULL;
    const size_t n = arg[0].ptr ? bondset_union_pairs(&pairs, as_bondset(arg[0]), element_count(arg[0]), ctx, md_temp_allocator(temp)) : 0;
    int result = 0;
    if (!dst) {
        if (ctx->backchannel) {
            ctx->backchannel->unit[0] = md_unit_none();
            ctx->backchannel->unit[1] = md_unit_angstrom();
            if (ctx->arg_flags[0] & FLAG_DYNAMIC) ctx->backchannel->flags |= FLAG_DYNAMIC_LENGTH;
        }
        result = ctx->vis ? 0 : (int)n;
    } else if (ctx->cur_state && ctx->cur_state->xyz) {
        float* out = as_float_arr(*dst);
        const size_t len = MIN(n, (size_t)element_count(*dst));
        for (size_t k = 0; k < len; ++k) {
            vec3_t dx = vec3_sub(ctx->cur_state->xyz[pairs[k].idx[1]], ctx->cur_state->xyz[pairs[k].idx[0]]);
            md_util_min_image_vec3(&dx, 1, &ctx->cur_state->unitcell);
            out[k] = vec3_length(dx);
        }
    }
    md_temp_end(temp);
    return result;
}

// The logic of bond sets, on the union of each side and within the context. 'not' is the bonds of the system which are
// not in the set: of the persistent bonds, as the bonds computed for a frame have no system to be the complement in.
static int bond_logic(data_t* dst, data_t arg[], int num_args, char op, eval_context_t* ctx) {
    if (!dst) return 0;
    const md_system_t* sys = ctx->sys;
    const uint32_t total = bond_total(sys);
    md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    md_bitmask_t x = md_bitmask_create(0, total, temp_alloc);
    bond_extras_t ex = {0};
    bondset_union_into(&x, &ex, as_bondset(arg[0]), element_count(arg[0]), temp_alloc);

    bond_extras_t out = {0};
    if (op == '!') {
        md_bitset_t set = md_bitmask_extract(&x, temp_alloc);
        md_bitmask_clear(&x);
        md_bitmask_invert_range(&x, 0, total);
        md_bitmask_andnot(&x, set);
    } else {
        ASSERT(num_args == 2);
        md_bitmask_t y = md_bitmask_create(0, total, temp_alloc);
        bond_extras_t ey = {0};
        bondset_union_into(&y, &ey, as_bondset(arg[1]), element_count(arg[1]), temp_alloc);
        const md_bitset_t ys = md_bitmask_extract(&y, temp_alloc);
        switch (op) {
        case '&': md_bitmask_and(&x, ys); break;
        case '|': md_bitmask_or(&x, ys); break;
        case '^': {
            const md_bitset_t xs = md_bitmask_extract(&x, temp_alloc);
            md_bitmask_or(&x, ys);
            md_bitmask_andnot(&x, md_bitset_and(xs, ys, temp_alloc));
            break;
        }
        default: ASSERT(false); break;
        }
        // The extras: a merge of the two sorted lists
        size_t i = 0, j = 0;
        const size_t nx = md_array_size(ex.key), ny = md_array_size(ey.key);
        while (i < nx || j < ny) {
            const bool take_x = j >= ny || (i < nx && ex.key[i] < ey.key[j]);
            const bool take_y = i >= nx || (j < ny && ey.key[j] < ex.key[i]);
            if (take_x) {
                if (op != '&') { md_array_push(out.key, ex.key[i], temp_alloc); md_array_push(out.flags, ex.flags[i], temp_alloc); }
                i += 1;
            } else if (take_y) {
                if (op != '&') { md_array_push(out.key, ey.key[j], temp_alloc); md_array_push(out.flags, ey.flags[j], temp_alloc); }
                j += 1;
            } else {
                // In both
                if (op != '^') { md_array_push(out.key, ex.key[i], temp_alloc); md_array_push(out.flags, ex.flags[i], temp_alloc); }
                i += 1;
                j += 1;
            }
        }
    }
    bondset_restrict_to_context(&x, &out, ctx);
    as_bondset(*dst)[0] = bondset_make(&x, &out, data_bitset_alloc(dst, ctx));
    md_temp_end(temp);
    return 0;
}

static int _bond_not(data_t* dst, data_t arg[], eval_context_t* ctx) { return bond_logic(dst, arg, 1, '!', ctx); }
static int _bond_and(data_t* dst, data_t arg[], eval_context_t* ctx) { return bond_logic(dst, arg, 2, '&', ctx); }
static int _bond_or (data_t* dst, data_t arg[], eval_context_t* ctx) { return bond_logic(dst, arg, 2, '|', ctx); }
static int _bond_xor(data_t* dst, data_t arg[], eval_context_t* ctx) { return bond_logic(dst, arg, 2, '^', ctx); }

// bond_delta(scope = none, *, tolerance = none, breakable = none): the bonds formed and broken in the frame against the
// bond graph of the system, for the bonds with an atom in the scope (md_util_bond_delta). A value of its own, which
// formed(), broken() and bond(d) read.
//   scope      The atoms whose bonds are looked at, all of them without it.
//   tolerance  Half the width of the band around the cutoffs, a fraction of them in [0, 0.5): 0.1 without it.
//   breakable  The origins whose bonds may break: 'file' and 'inferred' without it.
static int _bond_delta(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    const bool has_scope = !is_absent_type(arg[0].type);
    const float tolerance = is_absent_type(arg[1].type) ? MD_BOND_DELTA_TOLERANCE : as_float(arg[1]);
    uint32_t breakable = MD_BOND_DELTA_BREAKABLE;

    if (!dst && !ctx->vis) {
        if (!(tolerance >= 0.0f && tolerance < 0.5f)) {
            EVAL_LOG_ERROR(ctx, ctx->arg_tokens[1], "bond_delta: the tolerance is a fraction of the cutoffs, from 0 up to (not including) 0.5");
            return STATIC_VALIDATION_ERROR;
        }
        if (!is_absent_type(arg[2].type) && !bond_origin_mask(&breakable, arg[2], 2, true, ctx)) return STATIC_VALIDATION_ERROR;
        if (md_system_is_coarse_grained(ctx->sys)) {
            EVAL_LOG_ERROR(ctx, ctx->op_token, "bond_delta: the bonds of a coarse grained system come from its templates, not from distances, so none form or break");
            return STATIC_VALIDATION_ERROR;
        }
        if (ctx->backchannel) ctx->backchannel->unit[1] = md_unit_none();
        return 0;
    }
    if (!is_absent_type(arg[2].type)) bond_origin_mask(&breakable, arg[2], 2, false, ctx);
    if (!dst) return 0;

    const md_system_t* sys = ctx->sys;
    const uint32_t total = bond_total(sys);
    bond_delta_value_t* val = (bond_delta_value_t*)dst->ptr;
    MEMSET(val, 0, sizeof(bond_delta_value_t));
    // Asked about within the static check (the shape only): no search
    if (ctx->backchannel) return 0;

    md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_allocator_i* alloc = data_bitset_alloc(dst, ctx);

    md_bitset_t scope = {0};
    if (has_scope) scope = _internal_flatten_bf(as_bitset(arg[0]), element_count(arg[0]), temp_alloc);
    // Within a context, the scope is what of it is in the context
    if (ctx->mol_ctx) scope = has_scope ? md_bitset_and(scope, *ctx->mol_ctx, temp_alloc) : *ctx->mol_ctx;
    const md_bitset_t* scope_ptr = (has_scope || ctx->mol_ctx) ? &scope : NULL;

    // The persistent bonds touching the scope, but for the elastic and the metallic bonds: springs and the shells of a
    // metal, not bonds of a molecule, which neither break nor form
    md_bitmask_t touched = md_bitmask_create(0, total, temp_alloc);
    for (uint32_t i = 0; i < total; ++i) {
        const md_atom_pair_t p = sys->bond.pairs[i];
        if (md_bond_kind(sys->bond.flags[i]) == MD_BOND_KIND_ELASTIC || md_bond_kind(sys->bond.flags[i]) == MD_BOND_KIND_METALLIC) continue;
        if (!scope_ptr || md_bitset_test(scope, (uint32_t)p.idx[0]) || md_bitset_test(scope, (uint32_t)p.idx[1])) md_bitmask_set_bit(&touched, i);
    }
    val->touched = md_bitmask_extract(&touched, alloc);

    md_bond_delta_t d = {0};
    const md_bond_delta_desc_t desc = { .scope = scope_ptr, .tolerance = tolerance, .breakable = breakable };
    if (ctx->cur_state && md_util_bond_delta(&d, sys, ctx->cur_state, &desc, temp_alloc)) {
        val->broken = md_bitset_copy(d.broken, alloc);
        bond_extras_t e = {0};
        for (size_t k = 0; k < d.num_formed; ++k) bond_extras_push(&e, d.formed_pair[k], d.formed_flags[k], temp_alloc);
        md_bitmask_t none = md_bitmask_create(0, total, temp_alloc);
        val->formed = bondset_make(&none, &e, alloc);
    }
    md_temp_end(temp);
    return 0;
}

static inline const bond_delta_value_t* as_bond_delta(data_t arg) {
    return (const bond_delta_value_t*)arg.ptr;
}

// formed(d): the bonds formed in the frame, which the system does not have
static int _bond_delta_formed(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    const bond_delta_value_t* d = as_bond_delta(arg[0]);
    as_bondset(*dst)[0] = d ? bondset_copy(d->formed, data_bitset_alloc(dst, ctx)) : (bondset_t){0};
    return 0;
}

// broken(d): the bonds of the graph broken in the frame
static int _bond_delta_broken(data_t* dst, data_t arg[], eval_context_t* ctx) {
    if (!dst) return 0;
    const bond_delta_value_t* d = as_bond_delta(arg[0]);
    as_bondset(*dst)[0] = (bondset_t){ .bonds = d ? md_bitset_copy(d->broken, data_bitset_alloc(dst, ctx)) : md_bitset_empty_set() };
    return 0;
}

// bond(d, b = none, *, kind = none, order = none, resonance = none, origin = none): the bonds of the frame within the
// scope of d: the persistent bonds which touch it and are not broken, and those formed, within the context (both atoms
// in it) as bond() is, which b (an atom in it) and the filters let through. The filters read the flags of the system
// for a persistent bond, and those md_util_bond_delta gives a bond formed (inferred; covalent or coordinate, of no
// known order).
static int _bond_of_delta(data_t* dst, data_t arg[], eval_context_t* ctx) {
    ASSERT(ctx && ctx->sys);
    const bool report = !dst && !ctx->vis;
    bond_filter_t filter;
    if (!bond_filter_parse(&filter, arg, report, ctx)) return STATIC_VALIDATION_ERROR;
    if (!dst) return 0;

    const bond_delta_value_t* d = as_bond_delta(arg[0]);
    if (!d) {
        as_bondset(*dst)[0] = (bondset_t){0};
        return 0;
    }
    const md_system_t* sys = ctx->sys;
    md_temp_scope_t temp = md_temp_begin_avoid(ctx->temp_alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    const bool has_b = !is_absent_type(arg[1].type);
    md_bitset_t B = {0};
    if (has_b) B = _internal_flatten_bf(as_bitset(arg[1]), element_count(arg[1]), temp_alloc);

    md_bitmask_t b = md_bitmask_create(0, bond_total(sys), temp_alloc);
    md_bitset_iter_t it = md_bitset_iter(md_bitset_andnot(d->touched, d->broken, temp_alloc));
    while (md_bitset_iter_next(&it)) {
        if (it.idx >= bond_total(sys)) break;
        const md_atom_pair_t p = sys->bond.pairs[it.idx];
        if (has_b && !md_bitset_test(B, (uint32_t)p.idx[0]) && !md_bitset_test(B, (uint32_t)p.idx[1])) continue;
        if (bond_in_context(p, ctx) && bond_filter_pass(&filter, sys->bond.flags[it.idx])) md_bitmask_set_bit(&b, (uint32_t)it.idx);
    }
    bond_extras_t e = {0};
    for (uint32_t k = 0; k < d->formed.num_extra; ++k) {
        const md_atom_pair_t p = d->formed.extra_pair[k];
        if (has_b && !md_bitset_test(B, (uint32_t)p.idx[0]) && !md_bitset_test(B, (uint32_t)p.idx[1])) continue;
        if (bond_in_context(p, ctx) && bond_filter_pass(&filter, d->formed.extra_flags[k])) bond_extras_push(&e, p, d->formed.extra_flags[k], temp_alloc);
    }
    as_bondset(*dst)[0] = bondset_make(&b, &e, data_bitset_alloc(dst, ctx));
    md_temp_end(temp);
    return 0;
}
