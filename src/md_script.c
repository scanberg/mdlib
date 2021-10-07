#include "md_script.h"
#include "md_molecule.h"
#include "md_trajectory.h"
#include "md_frame_cache.h"
#include "md_filter.h"
#include "md_util.h"

#include "core/md_log.h"
#include "core/md_allocator.h"
#include "core/md_common.h"
#include "core/md_str.h"
#include "core/md_array.inl"
#include "core/md_file.h"
#include "core/md_arena_allocator.h"
#include "core/md_stack_allocator.h"
#include "core/md_compiler.h"
#include "core/md_bitop.h"
#include "core/md_vec_math.h"
#include "core/md_simd.h"
#include "core/md_platform.h"
#include "core/md_compiler.h"
#include "core/md_os.h"

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if MD_COMPILER_MSVC
#pragma warning(disable:4063) // Single character tokens not being valid tokens
#endif

#if MD_COMPILER_GCC || MD_COMPILER_CLANG
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wswitch"   // Single character tokens not part of enumeration
#endif

#define SCRIPT_IR_MAGIC 0x371bacfe8274910a
#define SCRIPT_EVAL_MAGIC 0x89bca715287bcaff

#define MAX_SUPPORTED_PROC_ARGS 8
#define MAX_SUPPORTED_TYPE_DIMS 4

// #############################
// ###   TYPE DECLARATIONS   ###
// #############################

typedef enum token_type_t {
    TOKEN_UNDEF = 0,
    // Reserve the first indices for character literals
    TOKEN_IDENT = 128, // idenfitier
    TOKEN_LE, // '<='
    TOKEN_GE, // '>='
    TOKEN_NE, // '!='
    TOKEN_EQ, // '=='
    TOKEN_AND,
    TOKEN_OR,
    TOKEN_NOT,
    TOKEN_OF,     // Reserved
    TOKEN_IN,     // Operator for setting the evaluation context.
    TOKEN_FLOAT,  // Floating point number
    TOKEN_INT,    // integer
    TOKEN_STRING, // String literal, defined between two quotation marks "LIKE THIS" or 'LIKE THIS'
    TOKEN_END     // End of tokenstream
} token_type_t;

typedef enum ast_type_t {
    AST_UNDEFINED = 0,
    AST_EXPRESSION,         // Parenthesis wrapped expression, We need to have it up here to ensure its precedence over things
    AST_PROC_CALL,          // Procedure call, Operators are directly translated into procedure calls as well.
    AST_CONSTANT_VALUE,
    AST_IDENTIFIER,
    AST_ARRAY,              // This is also just to declare an array in tree form, the node tells us the types, the children contains the values... inefficient, yes.
    AST_ARRAY_SUBSCRIPT,    // [..]
    AST_CAST,
    AST_UNARY_NEG,
    AST_NOT,
    AST_MUL,
    AST_DIV,
    AST_ADD,
    AST_SUB,
    AST_LT,
    AST_GT,
    AST_LE,
    AST_GE,
    AST_NE,
    AST_EQ,
    AST_AND,
    AST_OR,
    AST_CONTEXT,
    AST_ASSIGNMENT,
} ast_type_t;

typedef enum base_type_t {
    TYPE_UNDEFINED = 0,
    TYPE_FLOAT,
    TYPE_INT,
    TYPE_BOOL,
    TYPE_FRANGE,
    TYPE_IRANGE,
    TYPE_BITFIELD,      // Bitfield used to represent a selection of atoms, also contains a level which indicates the context (atom, residue, chain)
    TYPE_STRING,
} base_type_t;

typedef enum flags_t {
    // Function Flags
    FLAG_SYMMETRIC_ARGS             = 0x002, // Indicates that the arguments are symmetric, meaning the arguments can be swapped
    FLAG_ARGS_EQUAL_LENGTH          = 0x004, // Arguments should have the same array length
    FLAG_RET_AND_ARG_EQUAL_LENGTH   = 0x008, // Return type array length matches argument arrays' length
    FLAG_DYNAMIC_LENGTH             = 0x010, // Return type has a varying length which is not constant over frames
    FLAG_QUERYABLE_LENGTH           = 0x020, // Marks a procedure as queryable, meaning it can be called with NULL as dst to query the length of the resulting array
    FLAG_STATIC_VALIDATION          = 0x040, // Marks a procedure as validatable, meaning it can be called with NULL as dst to validate it in a static context during compilation
    FLAG_POSITION                   = 0x080, // Hints that float3 argument(s) is atomic position(s), this means we can implicitly compute this from alot of input types using the position procedure
    FLAG_CONSTANT                   = 0x100, // Hints that the identifier is constant and should not be modified.
    // Flags from 0x1000 and upwards are automatically propagated upwards the AST_TREE
    FLAG_DYNAMIC                    = 0x1000, // Indicates that it needs to be reevaluated for every frame of the trajectory (it has a dependency on atomic positions)
    FLAG_SDF                        = 0x2000, // Indicates that the expression involves an sdf computation
    FLAG_VISUALIZE                  = 0x4000, // Hints that a procedure can produce visualizations
} flags_t;

static const uint32_t FLAG_PROPAGATION_MASK = ~0xFFFU;

typedef struct token_t {
    token_type_t type;
    int32_t line;
    int32_t col_beg;
    int32_t col_end;
    str_t str;
} token_t;

typedef struct frange_t {
    float beg;
    float end;
} frange_t;

typedef struct irange_t {
    int32_t beg;
    int32_t end;
} irange_t;

typedef struct md_type_info_t {
    base_type_t  base_type;

    // This is how we encode the types dimensionality in dim
    // The example is given with base_type == float
    // 
    // dim = [0][0][0][0], len_dim = 0,  Uninitialized
    // dim = [1][0][0][0], len_dim = 0,  float[1]       (scalar)
    // 
    // dim = [2][0][0][0], len_dim = 0   float[2]       (float of length 2)
    // dim = [2][0][0][0], len_dim = 1   float[2][0]    (float[2] of length 0)
    // dim = [2][2][0][0], len_dim = 1   float[2][2]    (float[2] of length 2)
    // 
    // dim = [4][0][0][0], len_dim = 0   float[4]       (float of length 4)
    // dim = [4][1][0][0], len_dim = 1   float[4][1]    (float[4] of length 1)
    // 
    // dim = [4][4][0][0], len_dim = 1   float[4][4]    (float[4] of length 4)
    // dim = [4][4][1][0], len_dim = 2   float[4][4][1] (mat4 of length 1)
    // 
    // dim = [4][4][3][0], len_dim = 2   float[4][4][3] (mat4 of length 3)
    // 
    // Note that there is no differentiation between a true scalar and an array of length [1]. This is intentional to implicitly support scalars
    // as input arguments to procedures which operate on arrays of varying length.
    // 
    // Procedures encode their support of varying length in arguments and return values by encoding that particular dimension with -1
    // dim = [-1][0][0][0] Would mean that this argument accepts arrays of base_type of any length
    // dim = [4][-1][0][0] Would mean that this argument accepts arrays of base_type[4] of any length

    int32_t     dim[MAX_SUPPORTED_TYPE_DIMS]; // The dimensionality of a single element, each entry in dim represents a multidimensional length
    int32_t     len_dim;                      // Tells us which dimension in dim that encodes the length
//    context_level_t level;                    // Is only used for bitfields to signify which context they operate on. Should always be 0 otherwise
} md_type_info_t;

// There is a conceptual twist here if len should belong to the type_info_t or if it should belong to data_t
// We need to differentiate between different lengths of arrays with the same type in our system, so for now it should probably reside within the type_info
// I.e. the array length should be part of the fundamental type. This only strenghtens the type-system

// This is the struct which is passed as arguments into the procedures
typedef struct data_t {
    md_type_info_t      type;
    void*               ptr;    // Pointer to the data
    int64_t             size;   // Size in bytes of data (This we want to determine during the static check so we can allocate the data when evaluating)
    md_script_unit_t    unit;

    // Only used for distributions
    float               min_range;
    float               max_range;
} data_t;

// This is data stored directly in the nodes to hold scalar values
typedef union value_t {
    str_t       _string;
    float       _float;
    int32_t     _int;
    frange_t    _frange;
    irange_t    _irange;
    bool        _bool;
    md_exp_bitfield_t  _bitfield;
} value_t;

typedef struct identifier_t {
    str_t       name;
    data_t      data;
    flags_t     flags;
} identifier_t;

typedef struct md_script_visualization_o {
    uint64_t magic;
    md_allocator_i* alloc;
} md_script_visualization_o;

typedef struct eval_context {
    struct md_script_ir_o* ir;
    const md_molecule_t* mol;
    const md_exp_bitfield_t* mol_ctx;   // The atomic bit context in which we perform the operation, this can be null, in that case we are not limited to a smaller context and the full molecule is our context

    md_stack_allocator_t* stack_alloc;  // This is the same allocator as temp alloc, just with the raw interface
    md_allocator_i* temp_alloc;         // For allocating transient data
    md_allocator_i* alloc;              // For allocating persistent data (for the duration of the evaluation)

    token_t* arg_tokens;                // Tokens to arguments for contextual information when reporting errors
    identifier_t* identifiers;          // Evaluated identifiers for references
                                      
    md_script_visualization_t* vis;     // These are used when calling a procedure flagged with the VISUALIZE flag so the procedure can fill in the geometry

    struct {
        md_trajectory_frame_header_t* header;
        const float* x;
        const float* y;
        const float* z;
    } initial_configuration;
} eval_context_t;

typedef struct procedure {
    str_t name;
    md_type_info_t return_type;
    int64_t num_args;
    md_type_info_t arg_type[MAX_SUPPORTED_PROC_ARGS];
    int  (*proc_ptr)(data_t*, data_t[], eval_context_t*); // Function pointer to the "C" procedure
    flags_t flags;
} procedure_t;

typedef struct ast_node {
    // @OPTIMIZE: Make specific types for each type of Ast_node, where the first member is the type which we can use and cast after.
    ast_type_t      type;    
    token_t         token;          // Corresponding token from which the node was created (used to tie errors back into src)
    flags_t         flags;

    // PAYLOAD
    struct ast_node**   children;       // For AND, OR, NOT, procedure arguments etc.
    value_t             value;          // Scalar values for storing data directly inside the node
    data_t              data;           // Structure for passing as argument into procedure. (Holds pointer, length and type)

    procedure_t*    proc;           // Procedure reference
    str_t           ident;          // Identifier reference

    // CONTEXT
    md_type_info_t* lhs_context_types; // For static type checking
    const md_exp_bitfield_t* context;  // Since contexts have to be statically known at compile time, we store a reference to it. This enables independent evaluation of a node without missing the context.
} ast_node_t;

typedef struct tokenizer {
    str_t str;
    int64_t cur;
    int64_t line;
    int64_t line_offset;
} tokenizer_t;

typedef struct expression {
    ast_node_t*   node;     
    identifier_t* ident;    // Every expression should end up being assigned to an identifier. This is that identifier
    str_t str;
} expression_t;

typedef struct md_script_ir_o {
    uint64_t magic;

    // We could use a direct raw interface here to save some function pointer indirections
    struct md_allocator_i *arena;

    str_t str;  // Original string containing the 'source'
    
    // These are resizable arrays
    ast_node_t          **nodes;
    expression_t        **expressions;
    expression_t        **type_checked_expressions;     // List of expressions which passed type checking
    expression_t        **eval_targets;                 // List of dynamic expressions which needs to be evaluated per frame
    expression_t        **prop_expressions;             // List of expressions which are meant for exposure as properties
    
    identifier_t        *identifiers;                   // List of identifiers, notice that the data in a const context should only be used if it is flagged as

    // These are the final products which can be read through the public part of the structure
    md_script_error_t   *errors;
    md_script_token_t   *tokens;

    bool record_errors; // This is to toggle if new errors should be registered... We don't want to flood the errors
    const char* stage;  // This is just to supply a context for the errors i.e. which stage the error occured
} md_script_ir_o;

typedef struct md_script_eval_o {
    uint64_t magic;

    struct md_allocator_i *arena;
    int num_frames_completed;
    volatile bool interrupt;

    md_script_property_t    *properties;
} md_script_eval_o;

typedef struct parse_context_t {
    md_script_ir_o*    ir;
    tokenizer_t*    tokenizer;
    ast_node_t*     node;   // This represents the current root node in the tree
    md_allocator_i* temp_alloc;
} parse_context_t;


// ##########################
// ###   CORE FUNCTIONS   ###
// ##########################

static inline uint64_t generate_fingerprint() {
    return md_os_time_current();
}

static uint32_t operator_precedence(ast_type_t type) {
    switch(type) {
    case AST_EXPRESSION:
    case AST_PROC_CALL:
    case AST_CONSTANT_VALUE:
    case AST_IDENTIFIER:
    case AST_ARRAY:
    case AST_ARRAY_SUBSCRIPT:
        return 1;
    case AST_CAST:
    case AST_UNARY_NEG:
    case AST_NOT:
        return 2;
    case AST_MUL:
    case AST_DIV:
        return 3;
    case AST_ADD:
    case AST_SUB:
        return 4;
    case AST_LT:
    case AST_GT:
    case AST_LE:
    case AST_GE:
        return 5;
    case AST_NE:
    case AST_EQ:
        return 6;
    case AST_AND:
        return 7;
    case AST_OR:
        return 8;
    case AST_CONTEXT:
        return 9;
    case AST_ASSIGNMENT:
        return 10;
    default:
        ASSERT(false);
    }
    return 0;
}

static inline bool compare_type_info_dim(md_type_info_t a, md_type_info_t b) {
    return (memcmp(&a.dim, &b.dim, sizeof(a.dim)) == 0) && (a.len_dim == b.len_dim);
}

static inline bool compare_type_info(md_type_info_t a, md_type_info_t b) {
    /*if (a.base_type == TYPE_BITFIELD && b.base_type == TYPE_BITFIELD) {
        if ((a.level == -1 || b.level == -1) || a.level == b.level) {
            return compare_type_info_dim(a, b);
        }
        return false;
    }*/
    return memcmp(&a, &b, sizeof(md_type_info_t)) == 0;
}

static inline bool is_undefined_type(md_type_info_t ti) {
    return memcmp(&ti, &(md_type_info_t){0}, sizeof(md_type_info_t)) == 0;
}

static inline bool is_array(md_type_info_t ti) {
    return ti.dim[0] > 1;
}

static inline bool is_scalar(md_type_info_t ti) {
    return ti.dim[0] == 1 && ti.dim[1] == 0 && ti.len_dim == 0;
}

static inline bool is_variable_length(md_type_info_t ti) {
    return ti.dim[ti.len_dim] == -1;
}

/*
static inline bool is_incomplete_bitfield(md_type_info_t ti) {
    return ti.base_type == TYPE_BITFIELD && ti.level == -1;
}
*/

static inline int64_t type_info_array_len(md_type_info_t ti) {
    ASSERT(ti.len_dim < MAX_SUPPORTED_TYPE_DIMS);
    return ti.dim[ti.len_dim];
}

static inline int64_t element_count(data_t arg) {
    ASSERT(arg.type.dim[arg.type.len_dim] >= 0);
    return arg.type.dim[arg.type.len_dim];
}

// If the type is a scalar, then it is its own element type.
// If the type is an array, then its element type is the type of one element within the array
static inline md_type_info_t type_info_element_type(md_type_info_t ti) {
    if (is_scalar(ti)) return ti;
    md_type_info_t elem_ti = ti;
    elem_ti.dim[elem_ti.len_dim] = 1;
    return elem_ti;
}

// Size of the base type structure
static inline uint64_t base_type_element_byte_size(base_type_t type) {
    switch (type) {
    case TYPE_UNDEFINED: return 0;
    case TYPE_FLOAT:    return sizeof(float);
    case TYPE_INT:      return sizeof(int);
    case TYPE_BOOL:     return sizeof(bool);        
    case TYPE_FRANGE:   return sizeof(frange_t);
    case TYPE_IRANGE:   return sizeof(irange_t);
    case TYPE_BITFIELD: return sizeof(md_exp_bitfield_t);
    case TYPE_STRING:   return sizeof(str_t);
    default:            return 0;
    }
}

static inline bool is_type_dim_compatible(md_type_info_t from, md_type_info_t to) {

    // Some examples of matching, ^ shows len_dim
    // float[4][0][0][0] should match     float[-1][0][0][0]
    //       ^                                   ^
    // float[4][4][0][0] should match     float[4][-1][0][0]
    //          ^                                   ^
    // float[4][4][0][0] should not match float[4][-1][0][0]
    //             ^                                ^
    // 
    // This case should be nice to have, but is perhaps more confusing and too loose to support
    // float[4][4][1][0] should     match float[4][-1][0][0]
    //             ^                                ^


    // Is this always the case??? (Not if we want to support the last special case
    if (from.len_dim != to.len_dim) return false;

    for (int i = 0; i < MAX_SUPPORTED_TYPE_DIMS; ++i) {
        if (from.dim[i] == to.dim[i] && from.dim[i] > 0) continue;
        else if (from.dim[i] == -1 || to.dim[i] == -1) return true; // We consider a zero length array to be a valid type as well
    }
    return false;
}

static inline bool is_type_equivalent(md_type_info_t a, md_type_info_t b) {
    return memcmp(&a, &b, sizeof(md_type_info_t)) == 0;
}

static inline bool is_type_directly_compatible(md_type_info_t from, md_type_info_t to) {
    if (compare_type_info(from, to)) return true;

    if (from.base_type == to.base_type) {
        // This is essentially a logical XOR, we only want to support this operation if we have one unspecified array dimension.
        if (is_variable_length(from) != is_variable_length(to)) {
            return is_type_dim_compatible(from, to);
        }
    }

    return false;
}

static inline bool compare_type_info_array(const md_type_info_t a[], const md_type_info_t b[], int64_t num_arg_types) {
    for (int64_t i = 0; i < num_arg_types; ++i) {
        if (!is_type_directly_compatible(a[i], b[i])) {
            return false;
        }
    }
    return true;
}

// Returns the count of elements for this type, i.e an array of [4][4][1] would return 16
// NOTE: It does not include the last dimension which encodes the length of the array
// 
static inline int64_t type_info_element_stride_count(md_type_info_t ti) {
    if (ti.len_dim == 0) return 1;
    int64_t stride = ti.dim[0];
    for (int32_t i = 1; i < ti.len_dim; ++i) {
        stride *= ti.dim[i];
    }
    return stride;
}

static inline int64_t type_info_byte_stride(md_type_info_t ti) {
    return type_info_element_stride_count(ti) * base_type_element_byte_size(ti.base_type);
}

static inline int64_t type_info_total_byte_size(md_type_info_t ti) {
    return type_info_byte_stride(ti) * type_info_array_len(ti);
}

static inline int64_t type_info_total_element_count(md_type_info_t ti) {
    return type_info_element_stride_count(ti) * type_info_array_len(ti);
}

static inline int64_t bitfield_byte_size(int64_t num_bits) {
    return DIV_UP(num_bits, 64) * sizeof(int64_t);
}

static void* allocate_type(md_type_info_t type, md_allocator_i* alloc) {
    ASSERT(!is_undefined_type(type));
    ASSERT(alloc);

    void* data = 0;

    // The convention should be that we only allocate data for known lengths
    const int64_t array_len = type_info_array_len(type);
    ASSERT(array_len > -1);

    const int64_t bytes = type_info_byte_stride(type) * array_len;
    if (bytes > 0) {
        data = md_alloc(alloc, bytes);
        if (!data) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to allocate data in script!");
            return NULL;
        }
        memset(data, 0, bytes);

        if (type.base_type == TYPE_BITFIELD) {
            md_exp_bitfield_t* bf = data;
            for (int64_t i = 0; i < array_len; ++i) {
                md_bitfield_init(&bf[i], alloc);
            }
        }
    }

    return data;
}

static bool allocate_data(data_t* data, md_allocator_i* alloc) {
    ASSERT(data && !is_undefined_type(data->type));
    ASSERT(alloc);
    
    // The convention should be that we only allocate data for known lengths
    const int64_t array_len = type_info_array_len(data->type);
    ASSERT(array_len > -1);

    // Do the base type allocation (array)
    const int64_t bytes = type_info_byte_stride(data->type) * array_len;

    if (bytes > 0) {
        data->ptr = md_alloc(alloc, bytes);
        if (!data->ptr) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to allocate data in script!");
            return false;
        }
        memset(data->ptr, 0, bytes);
    }
    data->size = bytes;

    if (data->type.base_type == TYPE_BITFIELD) {
        md_exp_bitfield_t* bf = data->ptr;
        for (int64_t i = 0; i < array_len; ++i) {
            md_bitfield_init(&bf[i], alloc);
        }
    }

    return true;
}

static void free_data(data_t* data, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);    
    if (data->ptr && data->size) {
        if (data->type.base_type == TYPE_BITFIELD) {
            md_exp_bitfield_t* bf = data->ptr;
            const uint64_t num_elem = element_count(*data);
            for (uint64_t i = 0; i < num_elem; ++i) {
                md_bitfield_free(&bf[i]);
            }
        }
        md_free(alloc, data->ptr, data->size);
    }
    data->ptr = 0;
    data->size = 0;
}

static void copy_data(data_t* dst, const data_t* src) {
    ASSERT(dst);
    ASSERT(src);
    ASSERT(src->size);
    ASSERT(src->ptr);
    ASSERT(dst->ptr);
    ASSERT(dst->size == src->size);

    dst->type = src->type;
    dst->unit = src->unit;

    if (dst->type.base_type == TYPE_BITFIELD) {
        const uint64_t num_elem = element_count(*dst);

        md_exp_bitfield_t* dst_bf = dst->ptr;
        const md_exp_bitfield_t* src_bf = src->ptr;
        // Set bit pointers into destination memory
        for (uint64_t i = 0; i < num_elem; ++i) {
            md_bitfield_copy(&dst_bf[i], &src_bf[i]);
            //dst_bf->bits = (uint64_t*)((char*)dst->ptr + ((uint64_t)src_bf->bits - (uint64_t)src->ptr));
        }
    } else {
        memcpy(dst->ptr, src->ptr, src->size);
    }
}

static inline bool is_number(token_type_t type) {
    return type == TOKEN_INT || type == TOKEN_FLOAT;
}

static inline bool is_operator(ast_type_t type) {
    return (type == AST_ADD || type == AST_SUB || type == AST_MUL || type == AST_DIV ||
            type == AST_AND || type == AST_OR || type == AST_NOT ||
            type == AST_EQ || type == AST_NE || type == AST_LE || type == AST_GE || type == AST_LT || type == AST_GT);
}

static const char* get_token_type_str(token_type_t type) {
    switch (type) {
    case TOKEN_UNDEF: return "undefined";
    case '(': return "( (opening parentheses)";
    case ')': return ") (closing parentheses)";
    case '[': return "[ (opening bracket)";
    case ']': return "] (closing bracket)";
    case '<': return "< (less than)";
    case '>': return "> (greater than)";
    case '+': return "+ (addition)";
    case '-': return "- (subtraction or unary minus)";
    case '*': return "* (multiplication)";
    case '/': return "/ (division)";
    case '=': return "= (assignment)";
    case ':': return ": (range delimiter)";
    case ';': return "; (end of expression)";
    case ',': return ", (argument separator)";
    case TOKEN_IDENT: return "identifier";
    case TOKEN_LE: return "<= (less than or equal)";
    case TOKEN_GE: return ">= (greater than or equal)";
    case TOKEN_EQ: return "== (equals)";
    case TOKEN_AND: return "and";
    case TOKEN_OR: return "or";
    case TOKEN_NOT: return "not";
    case TOKEN_OF: return "of";
    case TOKEN_IN: return "in";
    case TOKEN_FLOAT: return "float";
    case TOKEN_INT: return "integer";
    case TOKEN_STRING: return "string";
    case TOKEN_END: return "end of token stream";
    default: return "symbol";
    }
}

static const char* get_value_type_str(base_type_t type) {
    switch(type) {
    case TYPE_UNDEFINED: return "undefined";
    case TYPE_FLOAT: return "float";
    case TYPE_INT: return "int";
    case TYPE_FRANGE: return "frange";
    case TYPE_IRANGE: return "irange";
    case TYPE_BOOL: return "boolean";
    case TYPE_STRING: return "string";
    case TYPE_BITFIELD: return "bitfield";
    default: return "type out of range";
    }
}

/*
static int64_t count_residues_in_mask(md_bitfield_t mask, const md_molecule* mol) {
    ASSERT(mol);
    ASSERT(mask.num_bits == mol->atom.count);
    ASSERT(mol->atom.residue_idx);
    ASSERT(mol->residue.atom_range);

    int64_t count = 0;
    uint64_t idx = 0;
    while (idx = bit_scan(mask.bits, idx, mask.num_bits - idx)) {
        ++count;
        const int64_t i = idx - 1;
        idx = mol->residue.atom_range[mol->atom.residue_idx[i]].end;
    }

    return count;
}

static int64_t count_chains_in_mask(md_bitfield_t mask, const md_molecule* mol) {
    ASSERT(mol);
    ASSERT(mask.num_bits == mol->atom.count);
    ASSERT(mol->atom.chain_idx);
    ASSERT(mol->chain.atom_range);

    int64_t count = 0;
    uint64_t idx = 0;
    while (idx = bit_scan(mask.bits, idx, mask.num_bits - idx)) {
        ++count;
        const int64_t i = idx - 1;
        idx = mol->chain.atom_range[mol->atom.chain_idx[i]].end;
    }

    return count;
}
*/

/*
static int64_t* get_atom_indices_in_mask(const md_exp_bitfield_t* mask, const md_molecule* mol, md_allocator_i* alloc) {
    ASSERT(mask);
    ASSERT(mol);
    ASSERT(alloc);
    ASSERT(mask->end_bit <= mol->atom.count);

    int64_t* result = NULL;

    int64_t beg_bit = mask->beg_bit;
    int64_t end_bit = mask->end_bit;
    while ((beg_bit = md_bitfield_scan(mask, beg_bit, end_bit)) != 0) {
        const int64_t idx = beg_bit - 1;
        md_array_push(result, idx, alloc);
    }
    return result;
}

static int64_t* get_residue_indices_in_mask(const md_exp_bitfield_t* mask, const md_molecule* mol, md_allocator_i* alloc) {
    ASSERT(mask);
    ASSERT(mol);
    ASSERT(alloc);
    ASSERT(mask->end_bit <= mol->atom.count);
    ASSERT(mol->atom.residue_idx);
    ASSERT(mol->residue.atom_range);

    int64_t* result = NULL;

    int64_t beg_bit = mask->beg_bit;
    int64_t end_bit = mask->end_bit;
    while ((beg_bit = md_bitfield_scan(mask, beg_bit, end_bit)) != 0) {
        const int64_t atom_idx = beg_bit - 1;
        int64_t res_idx = mol->atom.residue_idx[atom_idx];
        md_array_push(result, res_idx, alloc);
        // Skip to end of residue
        beg_bit = mol->residue.atom_range[res_idx].end;
    }
    return result;
}

static int64_t* get_chain_indices_in_mask(const md_exp_bitfield_t* mask, const md_molecule* mol, md_allocator_i* alloc) {
    ASSERT(mask);
    ASSERT(mol);
    ASSERT(alloc);
    ASSERT(mask->end_bit <= mol->atom.count);
    ASSERT(mol->atom.chain_idx);
    ASSERT(mol->chain.atom_range);

    int64_t* result = NULL;

    int64_t beg_bit = mask->beg_bit;
    int64_t end_bit = mask->end_bit;
    while ((beg_bit = md_bitfield_scan(mask, beg_bit, end_bit)) != 0) {
        const int64_t atom_idx = beg_bit - 1;
        int64_t chain_idx = mol->atom.chain_idx[atom_idx];
        md_array_push(result, chain_idx, alloc);
        // skip to end of chain
        beg_bit = mol->chain.atom_range[chain_idx].end;
    }
    return result;
}
*/

/*
static mol_context_t* get_residue_contexts_in_mask(const md_exp_bitfield_t* mask, const md_molecule* mol, md_allocator_i* alloc) {
    ASSERT(mask);
    ASSERT(mol);
    ASSERT(alloc);
    ASSERT(mol->residue.atom_range);    

    int64_t* res_idx = get_residue_indices_in_mask(mask, mol, alloc);
    mol_context_t* res_ctx = NULL;

    for (int64_t i = 0; i < md_array_size(res_idx); ++i) {
        const int64_t idx = res_idx[i];
        mol_context_t ctx = {
            .atom = {mol->residue.atom_range[idx].beg, mol->residue.atom_range[idx].end},
            .residue = {(uint32_t)idx, (uint32_t)idx+1},
            .chain = {0},
        };
        md_array_push(res_ctx, ctx, alloc);
    }

    return res_ctx;
}

static mol_context_t* get_chain_contexts_in_mask(const md_exp_bitfield_t* mask, const md_molecule* mol, md_allocator_i* alloc) {
    ASSERT(mask);
    ASSERT(mol);
    ASSERT(alloc);
    ASSERT(mol->chain.atom_range);
    ASSERT(mol->chain.residue_range);

    int64_t* chain_idx = get_chain_indices_in_mask(mask, mol, alloc);
    mol_context_t* chain_ctx = NULL;

    for (int64_t i = 0; i < md_array_size(chain_idx); ++i) {
        const int64_t idx = chain_idx[i];
        mol_context_t ctx = {
            .atom = {mol->chain.atom_range[idx].beg, mol->chain.atom_range[idx].end},
            .residue = {mol->chain.residue_range[idx].beg, mol->chain.residue_range[idx].end},
            .chain = {(uint32_t)idx, (uint32_t)idx+1},
        };
        md_array_push(chain_ctx, ctx, alloc);
    }

    return chain_ctx;
}
*/

static void create_error(md_script_ir_o* ir, token_t token, const char* format, ...) {
    if (!ir->record_errors) return;
    ir->record_errors = false;

    char buffer[512] = {0};
    va_list args;
    va_start(args, format);
    int len = vsnprintf(buffer, ARRAY_SIZE(buffer), format, args);
    va_end(args);

    ASSERT(len > 0);
    ASSERT(len < (int)ARRAY_SIZE(buffer));

    char* err_str = md_alloc(ir->arena, (uint64_t)len + 1);
    memcpy(err_str, buffer, len);
    err_str[len] = '\0';

    md_script_error_t error = {
        .line = (uint32_t)token.line,
        .col_beg = (uint32_t)(token.col_beg),
        .length = (uint32_t)token.str.len,
        .error = {.ptr = err_str, .len = (uint64_t)len},
    };
    md_array_push(ir->errors, error, ir->arena);

    // @TODO: Remove at some point
    md_printf(MD_LOG_TYPE_DEBUG, "%s line %llu: %s", ir->stage, token.line, buffer);

    if (token.str.ptr && token.str.len) {
        memset(buffer, 0, ARRAY_SIZE(buffer));
        const char* src_beg = ir->str.ptr;
        const char* src_end = ir->str.ptr + ir->str.len;
        const char* beg = token.str.ptr;
        const char* end = token.str.ptr + token.str.len;

        for (int i = 0; i < 50; ++i) {
            if (beg - 1 <= src_beg || beg[-1] == '\n') break;
            --beg;
        }
        for (int i = 0; i < 50; ++i) {
            if (end + 1 >= src_end || end[1] == '\n') break;
            ++end;
        }
        static const char* long_ass_carret = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
        md_printf(MD_LOG_TYPE_DEBUG, "%.*s", (end-beg), beg);
        md_printf(MD_LOG_TYPE_DEBUG, "%*s^%.*s", (token.str.ptr - beg), "", token.str.len-1, long_ass_carret);
    }
}

// Include all procedures, operators and defines
#include "core/script_functions.inl"




// ############################
// ###   HELPER FUNCTIONS   ###
// ############################

// Returns if the value type can be operated on using logical operators
static inline bool is_value_type_logical_operator_compatible(base_type_t type) {
    return type == TYPE_BOOL || type == TYPE_BITFIELD;
}

static inline bool is_identifier_procedure(str_t ident) {
    for (uint64_t i = 0; i < ARRAY_SIZE(procedures); ++i) {
        if (compare_str(ident, procedures[i].name)) {
            return true;
        }
    }
    return false;
}

static ast_node_t* create_node(md_script_ir_o* ir, ast_type_t type, token_t token) {
    ast_node_t* node = md_alloc(ir->arena, sizeof(ast_node_t));
    memset(node, 0, sizeof(ast_node_t));
    node->type = type;
    node->token = token;
    md_array_push(ir->nodes, node, ir->arena);
    return node;
}

static identifier_t* get_identifier(md_script_ir_o* ir, str_t name) {
    for (int64_t i = 0; i < (int64_t)ARRAY_SIZE(constants); ++i) {
        if (compare_str(constants[i].name, name)) {
            return &constants[i];
        }
    }

    for (int64_t i = 0; i < md_array_size(ir->identifiers); ++i) {
        if (compare_str(ir->identifiers[i].name, name)) {
            return &ir->identifiers[i];
        }
    }
    return NULL;
}

static identifier_t* create_identifier(md_script_ir_o* ir, str_t name) {
    ASSERT(get_identifier(ir, name) == NULL);

    identifier_t ident = {
        .name = copy_str(name, ir->arena),
        .data = {0},
        .flags = 0,
    };

    return md_array_push(ir->identifiers, ident, ir->arena);
}

static void fix_precedence(ast_node_t** node) {
    // We parse from left to right thus unbalancing the tree in its right child node (child[1])
    // Therefore we only need to fix potential precedence issues down this path
    // This is unfortunately not the case with keyword in, since that can be unbalanced down in the left subtree...
    // Therefore we check both

    if (!node) return;
    ast_node_t* parent = *node;
    if (!parent) return;

    const uint64_t num_children = md_array_size(parent->children);
    for (uint64_t i = 0; i < num_children; ++i) {
        ast_node_t* child = parent->children[i];
        if (!child || md_array_size(child->children) == 0) continue;
        ast_node_t* grand_child = child->children[0];

        if (operator_precedence(child->type) > operator_precedence(parent->type)) {
            child->children[0] = parent;
            parent->children[i] = grand_child;
            *node = child;
        }
    }
}

static bool is_type_implicitly_convertible(md_type_info_t from, md_type_info_t to, flags_t flags) {
    flags &= FLAG_POSITION; // This is the only flag we care about from the perspective of matching casting procedures
    for (uint64_t i = 0; i < ARRAY_SIZE(casts); ++i) {
        // if (flags && !(casts[i].flags & flags)) continue;
        if ((casts[i].flags & FLAG_POSITION) && !(flags & FLAG_POSITION)) continue; // Make sure we only consider casts which are marked with POSITION when the procedure is marked with POSITION
        if (is_type_directly_compatible(from, casts[i].arg_type[0]) &&
            is_type_directly_compatible(casts[i].return_type, to)) {
            return true;
        }
    }
    return false;
}

typedef struct procedure_match_result_t {
    bool success;
    procedure_t* procedure;
    uint32_t flags;
} procedure_match_result_t;

static procedure_match_result_t find_cast_procedure(md_type_info_t from, md_type_info_t to) {

    procedure_match_result_t res = {0};
    uint32_t lowest_cost = ~0U;

    // In the future, we might need to do two casts to actually get to the proper type.

    for (uint64_t i = 0; i < ARRAY_SIZE(casts); ++i) {
        if (is_type_directly_compatible(from, casts[i].arg_type[0]) &&
            is_type_directly_compatible(casts[i].return_type, to)) {

            uint32_t cost = 1;
            if (type_info_array_len(to) != -1 && type_info_array_len(casts[i].return_type) == -1) cost += 1; // Penalty for not being a perfect fit, there might be a better match which directly matches the array length
            if (type_info_array_len(to) == -1 && type_info_array_len(casts[i].return_type) != -1) cost += 1; // Same as above

            if (cost < lowest_cost) {
                res.success = true;
                res.procedure = &casts[i];
                lowest_cost = cost;
            }
        }
    }

    return res;
}

static procedure_match_result_t find_procedure_supporting_arg_types_in_candidates(str_t name, const md_type_info_t arg_types[], int64_t num_arg_types, procedure_t* candidates, int64_t num_cantidates, bool allow_implicit_conversions) {
    procedure_match_result_t res = {0};

    for (int64_t i = 0; i < num_cantidates; ++i) {
        procedure_t* proc = &candidates[i];
        if (compare_str(proc->name, name)) {
            if (num_arg_types == proc->num_args) {
                if (compare_type_info_array(arg_types, proc->arg_type, num_arg_types)) {
                    res.success = true;
                    res.procedure = proc;
                    return res;
                }
                else if (proc->flags & FLAG_SYMMETRIC_ARGS) {
                    // @TODO: Does this make sense for anything else than two arguments???
                    // I cannot come up with such a scenario.
                    ASSERT(proc->num_args == 2);
                    md_type_info_t swapped_args[2] = {arg_types[1], arg_types[0]};
                    if (compare_type_info_array(swapped_args, proc->arg_type, 2)) {
                        res.success = true;
                        res.procedure = proc;
                        res.flags |= FLAG_SYMMETRIC_ARGS;   // Flag this match as symmetric
                        return res;
                    }
                }
            }
        }
    }

    if (allow_implicit_conversions) {
        // Try using implicit casts on types where this is allowed.
        // int to float and irange to frange.
        // Store the best candidate with the lowest cost heuristic (number of conversions needed)
        procedure_t*    best_proc  = NULL;
        uint32_t        best_flags = 0;
        uint32_t        best_cost  = 0xFFFFFFFFU;

        for (int64_t i = 0; i < num_cantidates; ++i) {
            procedure_t* proc = &candidates[i];
            if (compare_str(proc->name, name)) {
                if (num_arg_types == proc->num_args) {
                    // Same name and same number of args...

                    uint32_t cost = 0;
                    uint32_t flags = 0;
                    for (int64_t j = 0; j < proc->num_args; ++j) {
                        if (is_type_directly_compatible(arg_types[j], proc->arg_type[j])) {
                            // No conversion needed for this argument (0 cost)
                        }
                        else if (is_type_implicitly_convertible(arg_types[j], proc->arg_type[j], proc->flags)) {
                            ++cost;
                        }
                        else {
                            // We are smoked.. This particular matching procedure cannot be resolved using implicit conversions...
                            cost = 0xFFFFFFFFU;
                            break;
                        }
                    }

                    if (cost < best_cost) {
                        best_cost = cost;
                        best_proc = proc;
                        best_flags = flags;
                    }

                    if (proc->flags & FLAG_SYMMETRIC_ARGS) {
                        ASSERT(proc->num_args == 2);
                        // Test if we can get a (better) result by swapping the arguments
                        md_type_info_t swapped_args[2] = {arg_types[1], arg_types[0]};

                        cost = 0;
                        flags = FLAG_SYMMETRIC_ARGS;
                        for (int64_t j = 0; j < proc->num_args; ++j) {
                            if (is_type_directly_compatible(swapped_args[j], proc->arg_type[j])) {
                                // No conversion needed for this argument (0 cost)
                            }
                            else if (is_type_implicitly_convertible(swapped_args[j], proc->arg_type[j], proc->flags)) {
                                ++cost;
                            }
                            else {
                                // We are smoked.. This particular matching procedure cannot be resolved using implicit conversions...
                                cost = 0xFFFFFFFFU;
                                break;
                            }
                        }

                        if (cost < best_cost) {
                            best_cost = cost;
                            best_proc = proc;
                            best_flags = flags;
                        }
                    }
                }
            }
        }
        if (best_cost < 0xFFFFFFFFU) {
            res.success = true;
            res.procedure = best_proc;
            res.flags = best_flags;
        }
    }

    return res;
}

static procedure_match_result_t find_procedure_supporting_arg_types(str_t name, const md_type_info_t arg_types[], uint64_t num_arg_types, bool allow_implicit_conversions) {
    return find_procedure_supporting_arg_types_in_candidates(name, arg_types, num_arg_types, procedures, ARRAY_SIZE(procedures), allow_implicit_conversions);
}

static procedure_match_result_t find_operator_supporting_arg_types(ast_type_t op, const md_type_info_t arg_types[], uint64_t num_arg_types, bool allow_implicit_conversions) {
    // Map token type to string which we use to identify operator procedures
    str_t name = {0};
    switch(op) {
    case AST_ADD: name = make_cstr("+");   break;
    case AST_SUB: name = make_cstr("-");   break;
    case AST_MUL: name = make_cstr("*");   break;
    case AST_DIV: name = make_cstr("/");   break;
    case AST_AND: name = make_cstr("and"); break;
    case AST_OR:  name = make_cstr("or");  break;
    case AST_NOT: name = make_cstr("not"); break;
    case AST_EQ:  name = make_cstr("==");  break;
    case AST_NE:  name = make_cstr("!=");  break;
    case AST_LT:  name = make_cstr("<");   break;
    case AST_GT:  name = make_cstr(">");   break;
    case AST_LE:  name = make_cstr("<=");  break;
    case AST_GE:  name = make_cstr(">=");  break;

    default:
        ASSERT(false);
    }

    return find_procedure_supporting_arg_types_in_candidates(name, arg_types, num_arg_types, operators, ARRAY_SIZE(operators), allow_implicit_conversions);
}

static identifier_t* find_constant(str_t name) {
    for (uint64_t i = 0; i < ARRAY_SIZE(constants); ++i) {
        if (compare_str(constants[i].name, name)) {
            return &constants[i];
        }
    }
    return NULL;
}

static inline bool is_token_type_comparison(token_type_t type) {
    return type == '<' || type == TOKEN_LE || type == '>' || type == TOKEN_GE || type == TOKEN_EQ;
}

static bool expect_token_type(md_script_ir_o* ir, token_t token, token_type_t type) {
    if (token.type != type) {
        create_error(ir, token, "Unexpected token '%.*s', expected token '%s'", token.str.len, token.str.ptr, get_token_type_str(type));
        return false;
    }
    return true;
}

// #####################
// ###   TOKENIZER   ###
// #####################

static token_t tokenizer_get_next_from_buffer(tokenizer_t* tokenizer) {
    token_t token = {0};

    if (tokenizer->cur >= tokenizer->str.len) {
        tokenizer->cur = tokenizer->str.len;
        token.type = TOKEN_END;
        token.line = (int32_t)tokenizer->line;
        token.col_beg = (int32_t)(tokenizer->cur - tokenizer->line_offset);
        return token;
    }
    
    int64_t line = tokenizer->line;
    const char* buf = tokenizer->str.ptr;
    const int64_t len = tokenizer->str.len;
    for (int64_t i = tokenizer->cur; i != len; ++i) {
        int64_t j = 0;
        
        if (buf[i] == '#') {
            // Skip to end of line
            while (++i < len) {
                if (buf[i] == '\n') {
                    break;
                }
            }
        }

        if (buf[i] == '\n') {
            ++line;
            tokenizer->line_offset = i + 1;
        }
        
        // Alpha numeric block
        else if (is_alpha(buf[i]) || buf[i] == '_') {
            for (j = i+1; j != len; ++j) {
                if (!is_alpha(buf[j]) && buf[j] != '_' && !is_digit(buf[j])) {
                    break;
                }
            }
            
            // Keywords
            const int64_t n = j - i;
            if (n == 2) {
                str_t str = {buf+i, 2};
                if (compare_str(str, (str_t){"or", 2})) {
                    token.type = TOKEN_OR;
                }
                else if (compare_str(str, (str_t){"in", 2})) {
                    token.type = TOKEN_IN;
                }
                else if (compare_str(str, (str_t){"of", 2})) {
                    token.type = TOKEN_OF;
                }
            }
            else if (n == 3) {
                str_t str = {buf+i, 3};
                if (compare_str(str, (str_t){"and", 3})) {
                    token.type = TOKEN_AND;
                }
                else if (compare_str(str, (str_t){"not", 3})) {
                    token.type = TOKEN_NOT;
                }
            }
            
            if (token.type == TOKEN_UNDEF) {
                token.type = TOKEN_IDENT;
            }
        }
        
        // Numeric literal
        else if (is_digit(buf[i]) || (i+1 < len && buf[i] == '-' && is_digit(buf[i+1]))) {
            bool is_float = false;
            for (j = i+1; j != len; ++j) {
                if (buf[j] == '.') {
                    is_float = true;
                    continue;
                } else if (!is_digit(buf[j])) {
                    break;
                }
            }
            if (is_float)
                token.type = TOKEN_FLOAT;
            else
                token.type = TOKEN_INT;
        }
        
        // String literal
        else if (buf[i] == '\'' || buf[i] == '\"') {
            char terminator = buf[i]; // Store the termination character to search for
            for (j = i+1; j != len; ++j) {
                if (buf[j] == terminator) {
                    token.type = TOKEN_STRING;
                    ++j;
                    break;
                }
                else if (buf[j] == ';' || (buf[j] == '\r' && buf[j+1] == '\n')) { // We do not want to leak return carry '\r' into the text since that will malform any output
                    break;
                }
            }
        }
        
        else if (is_symbol(buf[i])) {
            typedef struct Symbol {
                const char* str;
                token_type_t type;
            } Symbol;
            
            static Symbol symbols_2[] = {
                {"<=", TOKEN_LE},
                {">=", TOKEN_GE},
                {"==", TOKEN_EQ}
            };
            
            static char symbols_1[] = {
                '(',
                ')',
                '[',
                ']',
                '{',
                '}',
                '<',
                '>',
                '*',
                '+',
                '-',
                '/',
                '=',
                ':',
                ';',
                ',',
            };
            
            // Potential symbol count in buf to match (only care about 1 or 2)
            const int64_t n = (i + 1 < len && is_symbol(buf[i + 1])) ? 2 : 1;
            
            // Match buf + i against the longest accepted symbol that we can find.
            if (n == 2) {
                str_t str = {buf + i, 2};
                for (int k = 0; k < (int)ARRAY_SIZE(symbols_2); ++k) {
                    if (compare_str(str, (str_t){symbols_2[k].str, 2})) {
                        token.type = symbols_2[k].type;
                        j = i + 2;
                        break;
                    }
                }
            }
            if (!token.type) {
                for (int k = 0; k < (int)ARRAY_SIZE(symbols_1); ++k) {
                    if (symbols_1[k] == buf[i]) {
                        // The character is the type
                        token.type = symbols_1[k];
                        j = i + 1;
                        break;
                    }
                }
            }
        }
        
        if (j != 0) {
            token.str.ptr = tokenizer->str.ptr + i;
            token.str.len = j - i;
            token.line = (int32_t)line;
            token.col_beg = (int32_t)(i - tokenizer->line_offset);
            token.col_end = (int32_t)(j - tokenizer->line_offset);
            break;
        }

        if (!token.type && i >= len - 1) {
            // We are going to hit the end in next iteration and have nothing to show for.
            tokenizer->cur = tokenizer->str.len;
            token.type = TOKEN_END;
            token.line = (int32_t)line;
            token.col_beg = (int32_t)(i - tokenizer->line_offset);
            token.col_end = token.col_beg;
            break;
        }
    }

    return token;
}

static token_t tokenizer_consume_next(tokenizer_t* tokenizer) {
    token_t token = tokenizer_get_next_from_buffer(tokenizer);
    if (token.str.ptr && token.str.len) {
        // Advance current of tokenizer
        tokenizer->cur = token.str.ptr + token.str.len - tokenizer->str.ptr;
        tokenizer->line = token.line;
    }
    return token;
}

static token_t tokenizer_peek_next(tokenizer_t* tokenizer) {
    token_t token = tokenizer_get_next_from_buffer(tokenizer);
    return token;
}

static token_t tokenizer_consume_until_type(tokenizer_t* tokenizer, token_type_t* types, uint64_t num_types) {
    token_t token = {0};
    while (token = tokenizer_get_next_from_buffer(tokenizer), token.type != TOKEN_END) {
        for (uint64_t i = 0; i < num_types; ++i) {
            if (token.type == types[i]) goto end;
        }
        if (token.str.ptr && token.str.len) {
            // Advance current of tokenizer
            tokenizer->cur = token.str.ptr + token.str.len - tokenizer->str.ptr;
            tokenizer->line = token.line;
        }
    }
end:
    return token;
}

static tokenizer_t tokenizer_init(str_t str) {
    tokenizer_t tokenizer = {0};
    tokenizer.str = str;
    tokenizer.cur = 0;
    tokenizer.line = 1; // We start on line 1 not 0
    return tokenizer;
}


// ####################
// ###   PRINTING   ###
// ####################

static int print_type_info(char* buf, int buf_size, md_type_info_t info) {
    int len = snprintf(buf, buf_size, "%s", get_value_type_str(info.base_type));

    for (int i = 0; i < MAX_SUPPORTED_TYPE_DIMS; ++i) {
        if (info.dim[i] == -1) {
            len += snprintf(buf + len, MAX(0, buf_size - len), "[..]");
        }
        else if (info.dim[i] >= 1) {
            len += snprintf(buf + len, MAX(0, buf_size - len), "[%i]", (int)info.dim[i]);
        }
        else if (info.dim[i] == 0) {
            break;
        }
        else {
            ASSERT(false);
        }
    }

    return len;
}

#define PRINT(fmt, ...) len += snprintf(buf + len, MAX(0, buf_size - len), fmt, ##__VA_ARGS__)

static int print_bitfield(char* buf, int buf_size, const md_exp_bitfield_t* bitfield) {
    ASSERT(bitfield);
    if (buf_size <= 0) return 0;

    const int64_t max_bits = 64;
    const int64_t num_bits = MIN((int64_t)bitfield->end_bit - (int64_t)bitfield->beg_bit, max_bits);
    int len = 0;

    PRINT("<%i,%i>[", bitfield->beg_bit, bitfield->end_bit);
    for (int64_t i = 0; i < num_bits; ++i) {
        PRINT("%i", md_bitfield_test_bit(bitfield, i) ? 1 : 0);
    }
    PRINT("%s]", num_bits == max_bits ? "..." : "");

    return len;
}

static int print_value(char* buf, int buf_size, data_t data) {
    int len = 0;
    if (is_variable_length(data.type)) return 0;
    if (buf_size <= 0) return 0;

    if (data.ptr) {
        if (is_array(data.type) && data.type.base_type != TYPE_BITFIELD) {
            int64_t arr_len = data.type.dim[data.type.len_dim];
            if (arr_len == 1) {
                data.type.dim[data.type.len_dim] = 0;
                data.type.len_dim = MAX(0, data.type.len_dim-1);
                arr_len = data.type.dim[data.type.len_dim];
            }
            md_type_info_t type = type_info_element_type(data.type);
            int64_t stride = type_info_byte_stride(data.type);

            PRINT("{");
            for (int64_t i = 0; i < arr_len; ++i) {
                data_t elem_data = {
                    .ptr = (char*)data.ptr + stride * i,
                    .size = data.size,
                    .type = type,
                    .unit = data.unit,
                };
                len += print_value(buf + len, buf_size - len, elem_data);
                if (i < arr_len - 1) PRINT(",");
            }
            PRINT("}");
        } else {
            switch(data.type.base_type) {
            case TYPE_BITFIELD:
                //print_bitfield(file, (md_exp_bitfield_t*)data.ptr);
                break;
            case TYPE_BOOL:
                PRINT("%s", *(bool*)data.ptr ? "true" : "false");
                break;
            case TYPE_INT:
                PRINT("%i", *(int*)data.ptr);
                break;
            case TYPE_FLOAT:
                PRINT("%.2f", *(float*)data.ptr);
                break;
            case TYPE_IRANGE:
            {
                irange_t rng = *(irange_t*)data.ptr;
                if (rng.beg == INT32_MIN)
                    PRINT("int_min");
                else
                    PRINT("%i", rng.beg);
                PRINT(":");
                if (rng.end == INT32_MAX)
                    PRINT("int_max");
                else
                    PRINT("%i", rng.end);
                break;
            }
            case TYPE_FRANGE:
            {
                frange_t rng = *(frange_t*)data.ptr;
                if (rng.beg == -FLT_MAX)
                    PRINT("-flt_max");
                else
                    PRINT("%.2f", rng.beg);
                PRINT(":");
                if (rng.end == FLT_MAX)
                    PRINT("flt_max");
                else
                    PRINT("%.2f", rng.end);
                break;
            }
            case TYPE_STRING:
            {
                str_t str = *(str_t*)data.ptr;
                PRINT("%.*s", (int)str.len, str.ptr);
                break;
            }
            default:
                ASSERT(false);
            }
        }
    }
    else {
        PRINT("NULL");
    }
    return len;
}

#undef PRINT

#if MD_DEBUG

/*
static void print_value(FILE* file, data_t data) {
    if (is_variable_length(data.type)) return;

    if (is_array(data.type)) {
        int64_t len = data.type.dim[data.type.len_dim];
        if (len == 1) {
            data.type.dim[data.type.len_dim] = 0;
            data.type.len_dim = MAX(0, data.type.len_dim-1);
            len = data.type.dim[data.type.len_dim];
        }
        md_type_info_t type = type_info_element_type(data.type);
        int64_t stride = type_info_byte_stride(data.type);

        fprintf(file, "[");
        for (int64_t i = 0; i < len; ++i) {
            data_t new_data = {
                .ptr = (char*)data.ptr + stride * i,
                .size = data.size,
                .type = type,
                .unit = data.unit,
            };
            print_value(file, new_data);
            if (i < len - 1) fprintf(file, ", ");
        }
        fprintf(file, "]");
    } else {
        if (data.ptr) {
            switch(data.type.base_type) {
            case TYPE_BITFIELD:
                print_bitfield(file, (md_exp_bitfield_t*)data.ptr);
                break;
            case TYPE_BOOL:
                fprintf(file, "%s", *(bool*)data.ptr ? "true" : "false");
                break;
            case TYPE_INT:
                fprintf(file, "%i", *(int*)data.ptr);
                break;
            case TYPE_FLOAT:
                fprintf(file, "%.1f", *(float*)data.ptr);
                break;
            case TYPE_IRANGE:
            {
                irange_t rng = *(irange_t*)data.ptr;
                if (rng.beg == INT32_MIN)
                    fprintf(file, "int_min");
                else
                    fprintf(file, "%i", rng.beg);
                fprintf(file, ":");
                if (rng.end == INT32_MAX)
                    fprintf(file, "int_max");
                else
                    fprintf(file, "%i", rng.end);
                break;
            }
            case TYPE_FRANGE:
            {
                frange_t rng = *(frange_t*)data.ptr;
                if (rng.beg == -FLT_MAX)
                    fprintf(file, "-flt_max");
                else
                    fprintf(file, "%f", rng.beg);
                fprintf(file, ":");
                if (rng.end == FLT_MAX)
                    fprintf(file, "flt_max");
                else
                    fprintf(file, "%f", rng.end);
                break;
            }
            case TYPE_STRING:
            {
                str_t str = *(str_t*)data.ptr;
                fprintf(file, "%.*s", (int)str.len, str.ptr);
                break;
            }
            default:
                ASSERT(false);
            }
        }
        else {
            fprintf(file, "NULL");
        }
    }
}
*/

static void print_label(FILE* file, const ast_node_t* node) {
    char buf[1024];
    switch(node->type) {
    case AST_CONSTANT_VALUE:
        print_value(buf, sizeof(buf), node->data);
        fprintf(file, "%s", buf);
        break;
    case AST_PROC_CALL:
        if (node->proc) {
            fprintf(file, "%.*s", (int)node->proc->name.len, node->proc->name.ptr);
        }
        else {
            fprintf(file, "NULL");
        }
        fprintf(file, "(proc call)");
        break;
    case AST_ASSIGNMENT:
        fprintf(file, "assignment =");
        break;
    case AST_CAST:
        ASSERT(node->children);
        print_type_info(buf, ARRAY_SIZE(buf), node->children[0]->data.type);
        fprintf(file, "cast (%s) ->", buf);
        break;
    case AST_IDENTIFIER:
        fprintf(file, "%.*s", (int)node->ident.len, node->ident.ptr);
        fprintf(file, " (identifier)");
        break;
    case AST_ARRAY:
        fprintf(file, "array");
        break;
    case AST_ARRAY_SUBSCRIPT:
        fprintf(file, "array subscript");
        break;
    case AST_AND:
    case AST_OR:
    case AST_NOT:
        fprintf(file, "logical op %s", get_token_type_str(node->token.type));
        break;
    case AST_ADD:
    case AST_SUB:
    case AST_MUL:
    case AST_DIV:
        fprintf(file, "binary op %s", get_token_type_str(node->token.type));
        break;
    default:
        fprintf(file, "%s", get_token_type_str(node->token.type));
        break;
    }

    print_type_info(buf, sizeof(buf), node->data.type);
    fprintf(file, " -> %s", buf);

    if (node->flags) {
        fprintf(file, " {");
        if (node->flags & FLAG_DYNAMIC) fprintf(file, "D");
        fprintf(file, "}");
    }

    if (node->type != AST_ASSIGNMENT) {
        if (node->data.ptr) {
            fprintf(file, " = ");
            print_value(buf, sizeof(buf), node->data);
            fprintf(file, "%s", buf);
        }
    }
}


static inline void indent(FILE* file, int amount) {
    for (int i = 0; i < amount; ++i) {
        fprintf(file, "\t");
    }
}

static void print_node(FILE* file, const ast_node_t* node, int depth) {
    if (!node || !file) return;

    indent(file, depth);
    fprintf(file, "{\"name\": \"");
    print_label(file, node);
    fprintf(file, "\",\n");

    const uint64_t num_children = md_array_size(node->children);
    if (num_children) {
        indent(file, depth);
        fprintf(file, "\"children\": [\n");
        for (uint64_t i = 0; i < num_children; ++i) {
            print_node(file, node->children[i], depth + 1);
        }
        indent(file, depth);
        fprintf(file, "]\n");
    }

    indent(file, depth);
    fprintf(file, "},\n");
}

static void print_expr(FILE* file, str_t str) {

    // We need to add escape character to quotation marks, since they are not representable within a json string...
    char buf[1024] = {0};
    uint64_t at = 0;
    for (const char* c = str.ptr; c != str.ptr + str.len; ++c) {
        if (*c == '\"') {
            buf[at++] = '\\';
            buf[at++] = *c;
        }
        else if (*c == '\r') {

        }
        else if (*c == '\n') {
            buf[at++] = '\\';
            buf[at++] = 'n';
        }
        else {
            buf[at++] = *c;
        }
    }

    fprintf(file, "\"%.*s\"\n", (int)at, buf);
}

static void save_expressions_to_json(expression_t** expr, uint64_t num_expr, str_t filename) {
    FILE* file = (FILE*)md_file_open(filename, MD_FILE_WRITE);

    if (file) {
        fprintf(file, "var treeData = [\n");
        for (uint64_t i = 0; i < num_expr; ++i) {
            fprintf(file, "{\n\"node\" :\n");
            print_node(file, expr[i]->node, 1);
            fprintf(file, "\"flags\" : %u,\n", expr[i]->node->flags);
            fprintf(file, "\"expr\" : ");
            print_expr(file, expr[i]->str);
            fprintf(file, "},\n");

        }
        fprintf(file, "];");
        md_file_close((md_file_o*)file);
    }
}

#endif

// ###################
// ###   PARSING   ###
// ###################

static void expand_node_token_range_with_children(ast_node_t* node) {
    ASSERT(node);
    for (int64_t i = 0; i < md_array_size(node->children); ++i) {
        if (!node->children[i]) continue;
        expand_node_token_range_with_children(node->children[i]);
        node->token.col_beg = MIN(node->token.col_beg, node->children[i]->token.col_beg);
        node->token.col_end = MAX(node->token.col_end, node->children[i]->token.col_end);
    }
}

static ast_node_t* parse_expression(parse_context_t* ctx);

static ast_node_t** parse_comma_separated_arguments_until_token(parse_context_t* ctx, token_type_t token_type) {
    ast_node_t** args = 0;
    token_t next = {0};
    while (next = tokenizer_peek_next(ctx->tokenizer), next.type != TOKEN_END) {
        if (next.type == token_type) goto done;
        if (next.type == ',') {
            create_error(ctx->ir, next, "Empty argument in argument list");
            return NULL;
        } 
        ast_node_t* arg = parse_expression(ctx);
        next = tokenizer_peek_next(ctx->tokenizer);
        if (arg && (next.type == ',' || next.type == token_type)) {
            md_array_push(args, arg, ctx->temp_alloc);
            if (next.type == ',')
                tokenizer_consume_next(ctx->tokenizer);
            else // (next.type == token_type)
                goto done;
        } else {
            create_error(ctx->ir, next, "Unexpected token in argument list");
            return NULL;
        }
    }
done:
    return args;
}

static ast_node_t* parse_procedure_call(parse_context_t* ctx, token_t token) {
    // Here we know atleast that the token refers to something which is listed as a procedure!
    token_t next = tokenizer_peek_next(ctx->tokenizer);
    ast_node_t* node = 0;
    
    if (next.type == '(') {
        tokenizer_consume_next(ctx->tokenizer); // '('

        ast_node_t **args = parse_comma_separated_arguments_until_token(ctx, ')');
        next = tokenizer_consume_next(ctx->tokenizer);
        if (expect_token_type(ctx->ir, next, ')')) {
            node = create_node(ctx->ir, AST_PROC_CALL, token);
            const int64_t num_args = md_array_size(args);
            if (num_args) {
                md_array_push_array(node->children, args, num_args, ctx->ir->arena);
            }
            // Expand proc call to contain entire argument list ')'
            node->token.col_end = next.col_end;
        } else {
            create_error(ctx->ir, token, "Unexpected end of argument list");
        }
    } else {
        node = create_node(ctx->ir, AST_PROC_CALL, token);
    }

    return node;
}

ast_node_t* parse_identifier(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == TOKEN_IDENT);
    
    const str_t ident = token.str;
    ast_node_t* node = 0;
    identifier_t* c = 0;
    
    if (is_identifier_procedure(ident)) {
        // The identifier has matched some procedure name
        node = parse_procedure_call(ctx, token);
    } else if ((c = find_constant(ident)) != 0) {
        node = create_node(ctx->ir, AST_IDENTIFIER, token);
        //node->ident =;
        node->data = c->data;
    } else {
        // Identifier
        token_t next = tokenizer_peek_next(ctx->tokenizer);
        if (next.type == '(') {
            // This is an intended procedure call, but has no matching procedure by that name
            create_error(ctx->ir, token, "Undefined procedure '%.*s'", ident.len, ident.ptr);

            token_type_t token_types[] = {'(', ')', ';'};
            int paren_bal = 0;
            while (next = tokenizer_consume_until_type(ctx->tokenizer, token_types, ARRAY_SIZE(token_types)), next.type != TOKEN_END) {
                switch (next.type) {
                case '(':
                    tokenizer_consume_next(ctx->tokenizer);
                    ++paren_bal;
                    break;
                case ')':
                    tokenizer_consume_next(ctx->tokenizer);
                    if (--paren_bal == 0) {
                        goto done;
                    }
                    break;
                default: goto done;
                }
            }
        done:;
        }
        else if (next.type == '=') {
            // Assignment, therefore a new identifier
            if (!get_identifier(ctx->ir, ident)) {
                node = create_node(ctx->ir, AST_IDENTIFIER, token);
                create_identifier(ctx->ir, ident);
                node->ident = copy_str(ident, ctx->ir->arena);
            }
            else {
                create_error(ctx->ir, token, "The identifier is already taken. Variables cannot be reassigned.");
            }
        } else {
            // Identifier reference, resolve this later in the static check
            node = create_node(ctx->ir, AST_IDENTIFIER, token);
        }
    }
    return node;
}

ast_node_t* parse_logical(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == TOKEN_AND || token.type == TOKEN_OR || token.type == TOKEN_NOT);
    
    ast_node_t* arg[2] = {ctx->node, parse_expression(ctx)};   // lhs and rhs
    static const char* label[2] = {"Left hand side", "Right hand side"};
    const int num_operands = (token.type == TOKEN_NOT) ? 1 : 2;
    const int col_beg = (2 - num_operands);
    ast_node_t* node = 0;
    
    for (int i = col_beg; i < 2; ++i) {
        if (!arg[i]) {
            create_error(ctx->ir, token,
                         "%s of token '%s' did not evaluate to a valid expression.", label[i], get_token_type_str(token.type));
            return NULL;
        }
    }

    ast_type_t type = 0;
    switch (token.type) {
    case TOKEN_AND: type = AST_AND; break;
    case TOKEN_OR:  type = AST_OR;  break;
    case TOKEN_NOT: type = AST_NOT; break;
    default:
        ASSERT(false);
    }

    node = create_node(ctx->ir, type, token);
    md_array_push_array(node->children, arg + col_beg, num_operands, ctx->ir->arena);
    
    return node;
}

ast_node_t* parse_assignment(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == '=');
    
    ast_node_t* lhs = ctx->node;
    ast_node_t* rhs = parse_expression(ctx);
    
    if (!lhs || lhs->type != AST_IDENTIFIER) {
        create_error(ctx->ir, token,
                     "Left hand side of '%s' is not a valid identifier.", get_token_type_str(token.type));
        return NULL;
    }
    if (!rhs) {
        create_error(ctx->ir, token,
                     "Right hand side of token '%s' did not evaluate to a valid expression.", get_token_type_str(token.type));
        return NULL;
    }

    if (rhs->type == AST_ASSIGNMENT) {
        create_error(ctx->ir, token, "Syntax error: Invalid assignment");
        return NULL;
    }
    
    ast_node_t* node = create_node(ctx->ir, AST_ASSIGNMENT, token);
    md_array_push(node->children, lhs, ctx->ir->arena);
    md_array_push(node->children, rhs, ctx->ir->arena);
    node->data.type = rhs->data.type;
    return node;
}

ast_node_t* parse_comparison(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == '<' || token.type == '>' || token.type == TOKEN_EQ || token.type == TOKEN_LE || token.type == TOKEN_GE);
    
    ast_node_t* lhs = ctx->node;
    ast_node_t* rhs = parse_expression(ctx);
    ast_node_t* node = 0;

    if (lhs && rhs) {
        ast_type_t type = 0;
        switch(token.type) {
            case TOKEN_EQ: type = AST_EQ; break;
            case TOKEN_NE: type = AST_NE; break;
            case TOKEN_LE: type = AST_LE; break;
            case TOKEN_GE: type = AST_GE; break;
            case '<': type = AST_LT; break;
            case '>': type = AST_GT; break;
            default: ASSERT(false);
        }
        node = create_node(ctx->ir, type, token);
        ast_node_t* args[2] = {lhs, rhs};
        md_array_push_array(node->children, args, 2, ctx->ir->arena);
    }
   
    return node;
}

ast_node_t* parse_arithmetic(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == '-' || token.type == '+' || token.type == '*' || token.type == '/');

    // We don't support unary negation of arbitrary types yet...
    // If it is for a constant number, it is already baked into the number.
    
    ast_node_t* lhs = ctx->node;
    ast_node_t* rhs = parse_expression(ctx);
    ast_node_t* node = 0;
   
    if (lhs && rhs) {
        ast_type_t type = 0;
        switch(token.type) {
            case '+': type = AST_ADD; break;
            case '-': type = AST_SUB; break;
            case '*': type = AST_MUL; break;
            case '/': type = AST_DIV; break;
            default: ASSERT(false);
        }
        node = create_node(ctx->ir, type, token);
        ast_node_t* args[2] = {lhs, rhs};
        md_array_push_array(node->children, args, 2, ctx->ir->arena);
    }
    
    return node;
}

ast_node_t* parse_value(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == TOKEN_FLOAT || token.type == TOKEN_INT || token.type == TOKEN_STRING || token.type == ':');
    token_t next = tokenizer_peek_next(ctx->tokenizer);
    
    ast_node_t* node = create_node(ctx->ir, AST_CONSTANT_VALUE, token);
    
    if (token.type == ':' || (next.type == ':' && (is_number(token.type)))) {
        // Expand the nodes token to contain the range
        node->token.col_end = next.col_end;

        // RANGE
        // ... ugh! alot to parse!
        // Parse everything as double for now and convert to integer in the end if this bool is not set.
        bool is_float = (token.type == TOKEN_FLOAT);
        double beg = 0;
        double end = 0;

        if (is_number(token.type)) {
            beg = parse_float(token.str);
            expect_token_type(ctx->ir, tokenizer_consume_next(ctx->tokenizer), ':');
        } else { // == ':'
            beg = -FLT_MAX;
        }
        
        next = tokenizer_peek_next(ctx->tokenizer);
        if (is_number(next.type)) {
            token = tokenizer_consume_next(ctx->tokenizer);
            end = parse_float(token.str);
            is_float |= (token.type == TOKEN_FLOAT);
            // Expand the nodes token to contain the number
            node->token.col_end = next.col_end;
        } else {
            end = FLT_MAX;
        }

        if (is_float) {
            node->data.type = (md_type_info_t)TI_FRANGE;
            node->value._frange.beg = (float)beg;
            node->value._frange.end = (float)end;
        } else {
            node->data.type = (md_type_info_t)TI_IRANGE;
            node->value._irange.beg = (beg == -FLT_MAX) ?  INT32_MIN : (int32_t)beg;
            node->value._irange.end = (end ==  FLT_MAX) ?  INT32_MAX : (int32_t)end;
        }
    }
    else if (token.type == TOKEN_FLOAT) {
        node->data.type = (md_type_info_t)TI_FLOAT;
        node->value._float = (float)parse_float(token.str);
    }
    else if (token.type == TOKEN_INT) {
        node->data.type = (md_type_info_t)TI_INT;
        node->value._int = (int32_t)parse_int(token.str);
    }
    else if (token.type == TOKEN_STRING) {
        node->data.type = (md_type_info_t)TI_STRING;
        node->value._string = copy_str(substr(token.str, 1, MAX(0, ((int)token.str.len - 2))), ctx->ir->arena); // remove quotation marks
    }
    else {
        ASSERT(false);
    }
    return node;
}

ast_node_t* parse_array_subscript(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == '[');

    if (ctx->node) {
        ast_node_t* node = ctx->node;
        // We need to consume the node and turn it into a subscript
        // Put the original node as the first child, and insert whatever arguments we have as the following children.

        ast_node_t** elements = parse_comma_separated_arguments_until_token(ctx, ']');
        token_t next = tokenizer_consume_next(ctx->tokenizer);
        if (expect_token_type(ctx->ir, next, ']')) {
            const int64_t num_elements = md_array_size(elements);
            if (num_elements) {

                ast_node_t* node_copy = create_node(ctx->ir, node->type, token);
                memcpy(node_copy, node, sizeof(ast_node_t));

                memset(node, 0, sizeof(ast_node_t));
                node->type = AST_ARRAY_SUBSCRIPT;
                node->token = token;

                md_array_push(node->children, node_copy, ctx->ir->arena);
                md_array_push_array(node->children, elements, num_elements, ctx->ir->arena);

                // Expand token range to contain ']'
                node->token.col_end = next.col_end;

                return node;
            } else {
                create_error(ctx->ir, token, "Empty array subscripts are not allowed.");
            }
        } else {
            create_error(ctx->ir, token, "Unexpected end of argument list");
        }
    } else {
        create_error(ctx->ir, token, "Missing left hand side expression to apply subscript operator");
    }
    return NULL;
}

ast_node_t* parse_array(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == '{');

    ast_node_t* node = 0;

    ast_node_t** elements = parse_comma_separated_arguments_until_token(ctx, '}');
    token_t next = tokenizer_consume_next(ctx->tokenizer);
    if (expect_token_type(ctx->ir, next, '}')) {
        const int64_t num_elements = md_array_size(elements);
        if (num_elements) {
            // We only commit the results if everything at least parsed ok.
            node = create_node(ctx->ir, AST_ARRAY, token);
            md_array_push_array(node->children, elements, num_elements, ctx->ir->arena);
            // Expand token range with '}'
            node->token.col_end = next.col_end;
        } else {
            create_error(ctx->ir, token, "Empty arrays are not allowed.");
        }
    }
    else {
        create_error(ctx->ir, token, "Unexpected end of argument list");
    }
    return node;
}

ast_node_t* parse_in(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == TOKEN_IN);
    
    ast_node_t* lhs = ctx->node;
    ast_node_t* rhs = parse_expression(ctx);
    ast_node_t* node = 0;

    if (!lhs) {
        create_error(ctx->ir, token, "Left hand side of 'in' did not evaluate into a valid expression.");
    }
    else if(!rhs) {
        create_error(ctx->ir, token, "Right hand side of 'in' did not evaluate into a valid expression.");
    }
    else {
        node = create_node(ctx->ir, AST_CONTEXT, token);
        md_array_push(node->children, lhs, ctx->ir->arena);
        md_array_push(node->children, rhs, ctx->ir->arena);
    }

    return node;
}

ast_node_t* parse_parentheses(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == '(');
    ast_node_t* expr = parse_expression(ctx);
    ast_node_t* node = create_node(ctx->ir, AST_EXPRESSION, token);
    md_array_push(node->children, expr, ctx->ir->arena);
    token_t next = tokenizer_consume_next(ctx->tokenizer);
    if (expect_token_type(ctx->ir, next, ')')) {
        // Expand token with ')'
        node->token.col_end = next.col_end;
        return node;
    }
    return NULL;
}

static ast_node_t* parse_argument(parse_context_t* ctx) {
    token_t next = tokenizer_peek_next(ctx->tokenizer);
    switch (next.type) {
    case TOKEN_NOT:
        ctx->node = parse_logical(ctx);
        fix_precedence(&ctx->node);
        break;
    case TOKEN_IDENT:
        ctx->node = parse_identifier(ctx);
        break;
    case TOKEN_FLOAT:
    case TOKEN_INT:
    case TOKEN_STRING:
    case ':':
        ctx->node = parse_value(ctx);
        break;
    case '(':
        ctx->node = parse_parentheses(ctx);
        break;
    case '{':
        ctx->node = parse_array(ctx);
        break;
    case '-':
        ctx->node = parse_arithmetic(ctx);
        fix_precedence(&ctx->node);
        break;
    case TOKEN_UNDEF:
        create_error(ctx->ir, next, "Undefined token!");
        tokenizer_consume_next(ctx->tokenizer);
        return NULL;
    default:
        create_error(ctx->ir, next, "Unexpected token value! '%.*s'",
            next.str.len, next.str.ptr);
        return NULL;
    }
    return ctx->node;
}

ast_node_t* parse_expression(parse_context_t* ctx) {
    token_t token = {0};
    while (token = tokenizer_peek_next(ctx->tokenizer), token.type != TOKEN_END) {
        switch (token.type) {
            case TOKEN_AND:
            case TOKEN_OR:
            case TOKEN_NOT:
            ctx->node = parse_logical(ctx);
            fix_precedence(&ctx->node);
            break;
            case '=':
            ctx->node = parse_assignment(ctx);
            fix_precedence(&ctx->node);
            break;
            case TOKEN_GE:
            case TOKEN_LE:
            case TOKEN_EQ:
            ctx->node = parse_comparison(ctx);
            break;
            case TOKEN_IDENT:
            ctx->node = parse_identifier(ctx);
            break;
            case TOKEN_IN:
            ctx->node = parse_in(ctx);
            fix_precedence(&ctx->node);
            break;
            case TOKEN_FLOAT:
            case TOKEN_INT:
            case TOKEN_STRING:
            case ':':
            ctx->node = parse_value(ctx);
            break;
            case '(':
            ctx->node = parse_parentheses(ctx);
            break;
            case '[':
            ctx->node = parse_array_subscript(ctx);
            break;
            case '{':
            ctx->node = parse_array(ctx);
            break;
            case '-':
            case '+':
            case '*':
            case '/':
            ctx->node = parse_arithmetic(ctx);
            fix_precedence(&ctx->node);
            break;
            case ')':
            case ']':
            case '}':
            case ',':
            case ';':
            goto done;
            case TOKEN_UNDEF:
            create_error(ctx->ir, token, "Undefined token!");
            tokenizer_consume_next(ctx->tokenizer);
            return NULL;
            default:
            create_error(ctx->ir, token, "Unexpected token value! '%.*s'", token.str.len, token.str.ptr);
            return NULL;
        }
        if (!ctx->node) goto done;
    }
    
done:
    if (ctx->node) {
        expand_node_token_range_with_children(ctx->node);
    }
    return ctx->node;
}

// ####################
// ###   EVALUATE   ###
// ####################

static bool finalize_type(md_type_info_t* type, const ast_node_t* node, eval_context_t* ctx);

static bool evaluate_node(data_t*, const ast_node_t*, eval_context_t*);

static bool evaluate_proc_call(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(node && node->type == AST_PROC_CALL);
    ASSERT(node->proc);

    // This procedure can be called within the static check to 'query' the result of the procedure.
    // In such case, dst will be NULL.

    if (dst == NULL && ctx->vis) {
        // We are called within a visualization context, to visualize things
        // If our function is not marked as FLAG_VISUALIZE, we just stop here.
        if (!(node->proc->flags & FLAG_VISUALIZE)) {
            return true;
        }
    }

    const int64_t num_args = md_array_size(node->children);
    ASSERT(num_args < MAX_SUPPORTED_PROC_ARGS);

    data_t arg_data[MAX_SUPPORTED_PROC_ARGS] = {0};
    token_t arg_tokens[MAX_SUPPORTED_PROC_ARGS] = {0};

    for (int64_t i = 0; i < num_args; ++i) {
        // We need to evaluate the argument nodes first before we make the proc call.
        // In this context we are not interested in storing any data (since it is only used to be passed as arguments)
        // so we can allocate the data required for the node with the temp alloc
        arg_data[i] = node->children[i]->data;
        arg_tokens[i] = node->children[i]->token;

        arg_data[i].ptr = 0;

        if (is_variable_length(arg_data[i].type)) {
            if (!finalize_type(&arg_data[i].type, node->children[i], ctx)) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to finalize dynamic type in procedure call");
                return false;
            }
        }

        allocate_data(&arg_data[i], ctx->temp_alloc);
        if (!evaluate_node(&arg_data[i], node->children[i], ctx)) {
            return false;
        }
    }

    token_t* old_arg_tokens = ctx->arg_tokens;
    ctx->arg_tokens = arg_tokens;
    int result = node->proc->proc_ptr(dst, arg_data, ctx);
    ctx->arg_tokens = old_arg_tokens;

    for (int64_t i = num_args - 1; i >= 0; --i) {
        free_data(&arg_data[i], ctx->temp_alloc);
    }

    return result >= 0;
}

static bool evaluate_constant_value(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_CONSTANT_VALUE);
    ASSERT(node->data.ptr && node->data.size > 0);  // Assume that the static check set the data already
    (void)ctx;

    if (dst) {
        ASSERT(dst->ptr && dst->size >= node->data.size); // Make sure we have the expected size
        copy_data(dst, &node->data);
    }

    return true;
}

// This should only called in an evaluation context
static identifier_t* eval_find_identifier(str_t name, eval_context_t* ctx) {
    ASSERT(ctx);
    
    // Constants
    for (int64_t i = 0; i < (int64_t)ARRAY_SIZE(constants); ++i) {
        if (compare_str(name, constants[i].name)) {
            return &constants[i];
        }
    }

    // Static Identifiers from Compilation
    if (ctx->ir) {
        for (int64_t i = 0; i < md_array_size(ctx->ir->identifiers); ++i) {
            // Only return IR identifiers if they are constant
            if ((ctx->ir->identifiers[i].flags & FLAG_CONSTANT) && compare_str(name, ctx->ir->identifiers[i].name)) {
                return &ctx->ir->identifiers[i];
            }
        }
    }

    // Dynamic Identifiers from Evaluation
    for (int64_t i = 0; i < md_array_size(ctx->identifiers); ++i) {
        if (compare_str(name, ctx->identifiers[i].name)) {
            return &ctx->identifiers[i];
        }
    }

    return NULL;
}

static bool evaluate_identifier_reference(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_IDENTIFIER);
    (void)ctx;
    
    // We assume in this context that we have a reference to the identifier, we copy the data
    if (dst) {
        identifier_t* ident = eval_find_identifier(node->ident, ctx);
        if (ident) {
            ASSERT(dst->ptr && dst->size >= ident->data.size); // Make sure we have the expected size
            copy_data(dst, &ident->data);
        }
    }

    return true;
}

static bool evaluate_assignment(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ASSIGNMENT);
    ASSERT(node->children);

    ASSERT(md_array_size(node->children) == 2);
    ast_node_t* lhs = node->children[0];
    ast_node_t* rhs = node->children[1];

    ASSERT(lhs && lhs->type == AST_IDENTIFIER && lhs->ident.ptr);
    ASSERT(rhs);

    // We know from the static check that the lhs (child[0]) is an identifier and the rhs is a valid expression with a known data size
    // the data size could be zero if it is an array of length 0, in that case we don't need to allocate anything.
    // This should also be the conceptual entry point for all evaluations in this declarative language

    identifier_t* ident = eval_find_identifier(lhs->ident, ctx);

    if (!ident && ctx->alloc) {
        identifier_t id = {
            .name = lhs->ident,
            .data = lhs->data,
            .flags = lhs->flags
        };
        ident = md_array_push(ctx->identifiers, id, ctx->alloc);
        allocate_data(&ident->data, ctx->alloc);
    }

    bool result = false;
    if (dst || ctx->vis) {
        result = evaluate_node(dst, rhs, ctx);
        if (ident && (ident->flags & FLAG_CONSTANT) == 0 && dst) {
            ident->data = *dst;
        }
    }
    else if ((ident->flags & FLAG_CONSTANT) == 0) {
        result = evaluate_node(&ident->data, rhs, ctx);
    }

    return result;
}

static bool evaluate_array(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY);

    // Even if we do not have a dst parameter we still want to propagate the evaluation to the children
    const int64_t num_args = md_array_size(node->children);
    ast_node_t** const args = node->children;

    if (dst) {
        ASSERT(compare_type_info(dst->type, node->data.type));

        int64_t byte_offset = 0;
        for (int64_t i = 0; i < num_args; ++i) {
            data_t data = {
                .type = args[i]->data.type,
                .ptr = (char*)dst->ptr + byte_offset,
                .size = type_info_total_byte_size(args[i]->data.type)
            };
            if (!evaluate_node(&data, args[i], ctx)) {
                return false;
            }
            byte_offset += data.size;
        }
        ASSERT(byte_offset == dst->size);
    }
    else {
        for (int64_t i = 0; i < num_args; ++i) {
            if (!evaluate_node(NULL, args[i], ctx)) {
                return false;
            }
        }
    }

    return true;
}

static bool evaluate_array_subscript(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY_SUBSCRIPT);
    ASSERT(md_array_size(node->children) == 2);

    const ast_node_t* arr = node->children[0];
    const ast_node_t* idx = node->children[1];

    data_t arr_data = {arr->data.type};
    data_t idx_data = {idx->data.type};

    if (!allocate_data(&arr_data, ctx->temp_alloc)) return false;
    if (!allocate_data(&idx_data, ctx->temp_alloc)) return false;

    if (!evaluate_node(&arr_data, arr, ctx)) return false;
    if (!evaluate_node(&idx_data, idx, ctx)) return false;

    ASSERT(arr_data.ptr);
    ASSERT(idx_data.ptr);

    int offset = 0;
    int length = 1;
    if (compare_type_info(idx_data.type, (md_type_info_t)TI_INT)) {
        offset = *((int*)(idx_data.ptr));
    } else if (compare_type_info(idx_data.type, (md_type_info_t)TI_IRANGE)) {
        irange_t range = *((irange_t*)(idx_data.ptr));
        offset = range.beg;
        length = range.end - range.beg + 1;
    } else {
        ASSERT(false);
    }

    // We have front end notion of 1 based indexing to be coherent with the molecule conventions
    offset -= 1;

    if (dst) {
        ASSERT(compare_type_info(dst->type, node->data.type));
        ASSERT(dst->ptr);
        const int64_t elem_size = type_info_byte_stride(arr_data.type);
        const data_t src = {dst->type, (char*)arr_data.ptr + elem_size * offset, elem_size * length};
        copy_data(dst, &src);
    }

    return true;
}

static bool evaluate_context(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(node && node->type == AST_CONTEXT);
    ASSERT(md_array_size(node->children) == 2);

    const ast_node_t* expr_node = node->children[0];
    const ast_node_t* ctx_node  = node->children[1];

    // No need to reevaluate CTX_NODE here since it has already been done during static type checking and its results should be stored in the data ptr

    ASSERT(ctx_node->data.type.base_type == TYPE_BITFIELD);
    ASSERT(ctx_node->data.ptr);

    const int64_t num_ctx = type_info_array_len(ctx_node->data.type);
    const md_exp_bitfield_t* ctx_bf = (const md_exp_bitfield_t*)ctx_node->data.ptr;

    ASSERT(node->lhs_context_types);
    const md_type_info_t* lhs_types = node->lhs_context_types;
    ASSERT(md_array_size(lhs_types) == num_ctx);

    // Evaluate the expression within each context.
    int64_t dst_idx = 0;
    for (int64_t i = 0; i < num_ctx; ++i) {
        eval_context_t sub_ctx = *ctx;
        sub_ctx.mol_ctx = &ctx_bf[i];

        if (dst) {
            const uint64_t elem_size = type_info_byte_stride(dst->type);

            data_t data = {
                .type = lhs_types[i],
                .ptr = (char*)dst->ptr + elem_size * dst_idx,
                .size = elem_size
            };
            int64_t offset = type_info_array_len(lhs_types[i]);
            ASSERT(offset > 0);
            dst_idx += offset;

            if (!evaluate_node(&data, expr_node, &sub_ctx)) return false;

            if (i == 0) {
                dst->min_range = data.min_range;
                dst->max_range = data.max_range;
            }
        }
        else {
            if (!evaluate_node(NULL, expr_node, &sub_ctx)) return false;
        }
    }

    return true;
}

static bool evaluate_node(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);

    switch (node->type) {
        case AST_PROC_CALL:
            return evaluate_proc_call(dst, node, ctx);
        case AST_ASSIGNMENT:
            return evaluate_assignment(dst, node, ctx);
        case AST_ARRAY:
            return evaluate_array(dst, node, ctx);
        case AST_ARRAY_SUBSCRIPT:
            return evaluate_array_subscript(dst, node, ctx);
        case AST_CONSTANT_VALUE:
            return evaluate_constant_value(dst, node, ctx);
        case AST_IDENTIFIER:
            return evaluate_identifier_reference(dst, node, ctx);
        case AST_CONTEXT:
            return evaluate_context(dst, node, ctx);

        case AST_NOT:       // All operators should have been resolved during static check and converted into proc-calls
        case AST_UNARY_NEG:
        case AST_MUL:
        case AST_DIV:
        case AST_ADD:
        case AST_SUB:
        case AST_AND:
        case AST_OR:
        case AST_EQ:
        case AST_NE:
        case AST_LE:
        case AST_GE:
        case AST_LT:
        case AST_GT:
        case AST_CAST:      // Casts should have been resolved during static check and converted into a proc-call
        case AST_UNDEFINED:
        default:
            ASSERT(false);
    }
    return false;
}

// #############################
// ###   STATIC TYPE CHECK   ###
// #############################
// Determine the type of all nodes

static bool static_check_node(ast_node_t*, eval_context_t*);
static bool finalize_proc_call(ast_node_t*, procedure_match_result_t*, eval_context_t*);

static bool convert_node(ast_node_t* node, md_type_info_t new_type, eval_context_t* ctx) {
    const md_type_info_t from = node->data.type;
    md_type_info_t to   = new_type;

    procedure_match_result_t res = find_cast_procedure(from, to);
    if (res.success) {
        // We need to update the return type here
        to = res.procedure->return_type;

        // We need to check for the position flag here since that is not a true constant value, but a derived value.
        if (node->type == AST_CONSTANT_VALUE && !(res.procedure->flags & FLAG_POSITION)) {
            value_t val = node->value;
            data_t old = {
                .type = node->data.type,
                .size = node->data.size,
                .ptr = &val,
            };

            ASSERT(type_info_array_len(to) > -1);

            // Convert the type
            node->data.type = to;
            node->data.size = type_info_total_byte_size(to);

            if (to.base_type == TYPE_BITFIELD) {
                // Initialize that sucker (only type that needs to be initialized)
                md_bitfield_init(&node->value._bitfield, ctx->alloc);
            }

            // Perform the data conversion
            if (res.procedure->proc_ptr(&node->data, &old, ctx) == 0) {
                return true;
            }
        } else {
            // We need to convert this node into a cast node and add the original node data as a child
            ast_node_t* node_copy = create_node(ctx->ir, node->type, node->token);
            memcpy(node_copy, node, sizeof(ast_node_t));
            node->type = AST_PROC_CALL;
            node->children = 0; // node_copy have taken over the children, we need to zero this to trigger a proper allocation in next step
            md_array_push(node->children, node_copy, ctx->ir->arena);

            return finalize_proc_call(node, &res, ctx);
        }
    }

    ASSERT(false);
    return false;
}

static bool finalize_type(md_type_info_t* type, const ast_node_t* node, eval_context_t* ctx) {
    //ASSERT(is_variable_length(node->data.type));
    ASSERT(node->type == AST_PROC_CALL);
    ASSERT(type);

    *type = node->data.type;

    ast_node_t** const args = node->children;
    int64_t num_args = md_array_size(node->children);

    if (node->proc->flags & FLAG_RET_AND_ARG_EQUAL_LENGTH) {
        ASSERT(num_args > 0);
        // We can deduce the length of the array from the input type (take the array which is the longest???)
        int64_t max_len = 0;
        for (int64_t i = 0; i < num_args; ++i) {
            int64_t len = (int32_t)type_info_array_len(args[i]->data.type);
            ASSERT(len > -1);
            max_len = MAX(len, max_len);
        }
        type->dim[type->len_dim] = (int32_t)max_len;
    } else {
        ASSERT(node->proc->flags & FLAG_QUERYABLE_LENGTH);

        data_t arg_data[MAX_SUPPORTED_PROC_ARGS] = {0};
        token_t arg_tokens[MAX_SUPPORTED_PROC_ARGS] = {0};

        for (int64_t i = 0; i < num_args; ++i) {
            // We need to evaluate the argument nodes first before we make the proc call.
            // In this context we are not interested in storing any data (since it is only used to be passed as arguments)
            // so we can allocate the data required for the node with the temp alloc

            arg_data[i] = node->children[i]->data;
            arg_tokens[i] = node->children[i]->token;

            if (!arg_data[i].ptr) {
                if (is_variable_length(arg_data[i].type)) {
                    if (!finalize_type(&arg_data[i].type, node->children[i], ctx)) {
                        md_print(MD_LOG_TYPE_ERROR, "Failed to finalize dynamic type in procedure call");
                        return false;
                    }
                }
                allocate_data(&arg_data[i], ctx->temp_alloc);
                if (!evaluate_node(&arg_data[i], node->children[i], ctx)) {
                    return false;
                }
            }
            /*
            arg_data[i].type = node->children[i]->data.type;
            arg_tokens[i] = node->children[i]->token;
            allocate_data(&arg_data[i], ctx->temp_alloc);
            if (!evaluate_node(&arg_data[i], node->children[i], ctx)) {
                return false;
            }
            */
        }

        eval_context_t eval_ctx = *ctx;
        eval_ctx.arg_tokens = arg_tokens;

        // Perform the call
        int query_result = node->proc->proc_ptr(NULL, arg_data, &eval_ctx);

        // Free
        for (int64_t i = num_args - 1 ; i >= 0; --i) {
            free_data(&arg_data[i], ctx->temp_alloc);
        }

        if (query_result >= 0) { // Zero length is valid
            type->dim[type->len_dim] = query_result;
        } else if (query_result == 0) {
            create_error(ctx->ir, node->token, "The supplied argument is of zero length.", query_result);
            return false;
        } else {
            create_error(ctx->ir, node->token, "Unexpected return value (%i) when querying procedure for array length.", query_result);
            return false;
        }
    }

    return true;
}

static bool finalize_proc_call(ast_node_t* node, procedure_match_result_t* res, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(res);
    ASSERT(ctx);

    node->proc = res->procedure;
    node->data.type = node->proc->return_type;

    const int64_t num_args = md_array_size(node->children);
    ast_node_t** arg = node->children;

    if (num_args > 0) {
        if (res->flags & FLAG_SYMMETRIC_ARGS) {
            // If the symmetric flag is set, we need to swap the arguments
            ASSERT(num_args == 2);
            ast_node_t* tmp_child = node->children[0];
            node->children[0] = node->children[1];
            node->children[1] = tmp_child;
        }

        // Make sure we cast the arguments into the expected types of the procedure.
        for (int64_t i = 0; i < num_args; ++i) {
            if (!is_type_directly_compatible(arg[i]->data.type, node->proc->arg_type[i])) {
                // Types are not directly compatible, but should be implicitly convertible (otherwise the match should never have been found)
                ASSERT(is_type_implicitly_convertible(arg[i]->data.type, node->proc->arg_type[i], node->proc->flags));
                if (!convert_node(arg[i], node->proc->arg_type[i], ctx)) {
                    return false;
                }
                node->flags |= arg[i]->flags & FLAG_PROPAGATION_MASK;
            }
        }

        if (node->proc->flags & FLAG_ARGS_EQUAL_LENGTH) {
            ASSERT(num_args > 1); // This flag should only be set for procedures with multiple arguments
                                  // Make sure that the input arrays are of equal length. Otherwise create an error.
            const int64_t expected_len = type_info_array_len(arg[0]->data.type);
            for (int64_t i = 1; i < num_args; ++i) {
                int64_t len = type_info_array_len(arg[i]->data.type);
                if (len != expected_len) {
                    create_error(ctx->ir, node->token, "Expected array-length of arguments to match. arg 0 has length %i, arg %i has length %i.", (int)expected_len, (int)i, (int)len);
                    return false;
                }
            }
        }
    }

    // Propagate procedure flags
    node->flags |= node->proc->flags & FLAG_PROPAGATION_MASK;

    // Need to flag DYNAMIC_LENGTH if the node has QUERYABLE_LENGTH and depends on DYNAMIC arguments
    // Or if the node is evaluated within a context and has QUERYABLE_LENGTH
    if (node->proc->flags & FLAG_QUERYABLE_LENGTH) {
        if (ctx->mol_ctx)
            node->flags |= FLAG_DYNAMIC_LENGTH;
        else {
            for (int64_t i = 0; i < num_args; ++i) {
                if (arg[i]->flags & FLAG_DYNAMIC) {
                    node->flags |= FLAG_DYNAMIC_LENGTH;
                    break;
                }
            }
        }
    }

    if (is_variable_length(node->data.type)) {
        // The nodes type has a variable length. We would like to resolve it statically in order to determine the type
        if (node->flags & FLAG_DYNAMIC_LENGTH) {
            // This node has been marked as returning a variable length, which means it has a dependency on something which cannot be resolved statically during compile time.
            // Therefore we cannot deduce the length of the type statically and we need to yield this operation until evaluation.
            // This is conceptually fine for intermediate values, but should never be stored into actual properties.
            // Example: distance(com(x(1:10)), vec3(1,0,0));
            return true;
        }
        else {
            return finalize_type(&node->data.type, node, ctx);
        }
    }

    return true;
}

static bool static_check_children(ast_node_t* node, eval_context_t* ctx) {
    const int64_t num_children = md_array_size(node->children);
    for (int64_t i = 0; i < num_children; ++i) {
        ASSERT(node != node->children[i]);
        if (static_check_node(node->children[i], ctx)) {
            node->flags |= node->children[i]->flags & FLAG_PROPAGATION_MASK;
        } else {
            return false;
        }
    }
    return true;
}

static bool static_check_operator(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(is_operator(node->type));
    ASSERT(node->children);

    if (static_check_children(node, ctx)) {
        const int64_t num_args = md_array_size(node->children);
        ast_node_t** arg = node->children;

        md_type_info_t arg_type[2] = {arg[0]->data.type, num_args > 1 ? arg[1]->data.type : (md_type_info_t){0}};
        procedure_match_result_t res = find_operator_supporting_arg_types(node->type, arg_type, num_args, true);
        if (res.success) {
            node->type = AST_PROC_CALL;
            return finalize_proc_call(node, &res, ctx);
        } else {
            char arg_type_str[2][64] = {"empty", "empty"};
            for (int i = 0; i < num_args; ++i) {
                if (arg[i]) print_type_info(arg_type_str[i], ARRAY_SIZE(arg_type_str[i]), arg[i]->data.type);                
            }
            create_error(ctx->ir, node->token,
                "Could not find supporting operator '%s' with left hand side type (%s) and right hand side type (%s)",
                get_token_type_str(node->token.type), arg_type_str[0], arg_type_str[1]);
        }
    }

    return false;
}

static bool static_check_proc_call(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node->type == AST_PROC_CALL);

    bool result = true;

    if (!node->proc) {
        result = false;
        if (static_check_children(node, ctx)) {
            const uint64_t num_args = md_array_size(node->children);
            const str_t proc_name = node->token.str;

            if (num_args == 0) {
                // Zero arguments
                procedure_match_result_t res = find_procedure_supporting_arg_types(proc_name, 0, 0, false);
                if (res.success) {
                    result = finalize_proc_call(node, &res, ctx);
                }
                else {
                    create_error(ctx->ir, node->token,
                        "Could not find matching procedure '%.*s' which takes no arguments", (int)proc_name.len, proc_name.ptr);
                }
            }
            else {
                // One or more arguments
                ASSERT(num_args < MAX_SUPPORTED_PROC_ARGS);
                md_type_info_t arg_type[MAX_SUPPORTED_PROC_ARGS] = {0};

                for (uint64_t i = 0; i < num_args; ++i) {
                    arg_type[i] = node->children[i]->data.type;
                }

                procedure_match_result_t res = find_procedure_supporting_arg_types(proc_name, arg_type, num_args, true);
                if (res.success) {
                    result = finalize_proc_call(node, &res, ctx);
                }
                else {
                    char buf[512] = {0};
                    int  cur = 0;
                    for (uint64_t i = 0; i < num_args; ++i) {
                        cur += print_type_info(buf + cur, ARRAY_SIZE(buf) - cur, arg_type[i]);
                        if (i + 1 != num_args) {
                            cur += snprintf(buf + cur, ARRAY_SIZE(buf) - cur, ",");
                        }
                    }

                    create_error(ctx->ir, node->token,
                        "Could not find matching procedure '%.*s' which takes the following argument(s) (%.*s)",
                        (int)proc_name.len, proc_name.ptr, cur, buf);
                }
            }
        }
    }

    if (result) {
        // Perform new child check here since children may have changed due to conversions etc.
        result &= static_check_children(node, ctx);

        ASSERT(node->proc);
        // Perform static validation by evaluating with NULL ptr
        if (node->proc->flags & FLAG_STATIC_VALIDATION) {
            result &= evaluate_proc_call(NULL, node, ctx);
        }
    }

    return result;
}

static bool static_check_constant_value(ast_node_t* node, eval_context_t* ctx) {
    (void)ctx;
    ASSERT(node && node->type == AST_CONSTANT_VALUE);
    ASSERT(node->data.type.base_type != TYPE_UNDEFINED);
    ASSERT(is_scalar(node->data.type));

    node->data.ptr = &node->value;
    node->data.size = base_type_element_byte_size(node->data.type.base_type);

    if (node->data.type.base_type == TYPE_IRANGE) {
        // Make sure this range is not negative
        irange_t rng = node->value._irange;
        if (rng.beg > rng.end) {
            create_error(ctx->ir, node->token, "The range is invalid, a range must have an ascending format, did you mean '%i:%i'?", rng.end, rng.beg);
            return false;
        }
    }
    
    return true;
}

static bool static_check_array(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY && node->children);

    // I think the rules should be as following:
    // Allow integers, floats, ranges within the array.
    // If any float is present, then the array's base type gets promoted to float
    // if any intrange is present, all integers are promoted to intranges
    // We do not promote floats into float ranges for now, while possible, it may introduce errors where floats are incorrectly passed as arguments
    // via implicit conversion to floatrange. This would then be a zero extent range (since float).

    if (static_check_children(node, ctx)) {

        const int64_t  num_elem = md_array_size(node->children);
        ast_node_t**        elem = node->children;

        md_type_info_t array_type = {0};
        int64_t array_len = 0;
        
        // Pass 1: find a suitable base_type for the array
        for (int64_t i = 0; i < num_elem; ++i) {
            md_type_info_t elem_type = type_info_element_type(elem[i]->data.type);
            array_len += type_info_array_len(elem[i]->data.type);

            if (array_type.base_type == TYPE_UNDEFINED) {
                // Assign type from first element
                //array_type = elem[i]->data.type;
                array_type = elem_type;
            }
            else {
                if (!is_type_directly_compatible(elem_type, array_type)) {
                    if (is_type_implicitly_convertible(elem[i]->data.type, array_type, 0)) {
                        // We can convert the element to the array type.   
                    }
                    else if (is_type_implicitly_convertible(array_type, elem_type, 0)) {
                        // Option 2: We can convert the array type into the elements type (e.g int to float)
                        // Promote the type of the array
                        array_type = elem_type;

                        // Just to make sure, we recheck all elements up until this one
                        for (int64_t j = 0; j < i; ++j) {
                            if (!is_type_directly_compatible(elem_type, array_type) &&
                                !is_type_implicitly_convertible(elem_type, array_type, 0)) {
                                create_error(ctx->ir, elem[i]->token, "Incompatible types wihin array construct");
                                return false;
                            }
                        }
                    }
                    else {
                        // Incompatible types...
                        create_error(ctx->ir, elem[i]->token, "Incompatible types wihin array construct");
                        return false;
                    }
                }
            }
        }

        // Pass 2: Perform implicit conversions of nodes if required
        for (int64_t i = 0; i < num_elem; ++i) {
            md_type_info_t converted_type = array_type;
            converted_type.dim[converted_type.len_dim] = (int32_t)type_info_array_len(elem[i]->data.type);

            if (!is_type_directly_compatible(elem[i]->data.type, converted_type) &&
                is_type_implicitly_convertible(elem[i]->data.type, converted_type, 0)) {
                if (!convert_node(elem[i], converted_type, ctx)) {
                    return false;
                }
            }
        }

        // Finalize the type for the array
        array_type.dim[array_type.len_dim] = (int32_t)array_len;
        node->data.type = array_type;

        return true;
    }
    return false;
}

static bool static_check_array_subscript(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY_SUBSCRIPT && node->children);

    // In array subscripts, only integers and ranges should be allowed.
    // For now, we only support SINGLE integers or ranges
    // @TODO: In the future, se should expand this to support mixed integers and ranges and create a proper subset of the original data
    // @TODO: evaluate the children in order to determine the length and values (to make sure they are in range)

    if (static_check_children(node, ctx)) {
        const uint64_t  num_elem = md_array_size(node->children);
        ast_node_t**        elem = node->children;

        if (num_elem < 2) {
            create_error(ctx->ir, node->token, "Missing arguments in array subscript");
            return false;
        }
        else if (num_elem > 2) {
            create_error(ctx->ir, elem[2]->token, "Only single entries are allowed inside array subscript");
            return false;
        }

        ast_node_t* lhs = elem[0];
        ast_node_t* arg = elem[1];

        if (is_scalar(arg->data.type)) {
            if (arg->data.type.base_type == TYPE_INT || arg->data.type.base_type == TYPE_IRANGE) {
                irange_t range = {0};
                if (arg->data.type.base_type == TYPE_INT) {
                    range.beg = arg->value._int;
                    range.end = arg->value._int;
                } else {
                    range = arg->value._irange;
                }
                if (range.beg <= range.end && 1 <= range.beg && range.end <= element_count(lhs->data)) {
                    ASSERT(lhs->data.type.dim[0] > 0);
                    // SUCCESS!
                    node->data.type = lhs->data.type;
                    node->data.type.dim[node->data.type.len_dim] = range.end - range.beg + 1;
                    return true;
                } else {
                    create_error(ctx->ir, arg->token, "Invalid Array subscript range");
                    return false;
                }
            } else {
                create_error(ctx->ir, arg->token, "Only int and int-ranges are allowed inside array subscript");
                return false;
            }
        } else {
            create_error(ctx->ir, arg->token, "No arrays are allowed inside array subscript");
            return false;
        }
    }

    return false;
}

static bool static_check_identifier_reference(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_IDENTIFIER);
    
    if (!node->ident.ptr) {
        node->ident = node->token.str;
    }

    identifier_t* ident = get_identifier(ctx->ir, node->ident);

    if (ident) {
        node->data = ident->data;
        node->flags |= ident->flags & FLAG_PROPAGATION_MASK;    // This means it has a dependency on an identifier, and that identifier might be dynamic... we need to propagate flags

        if (ident->data.type.base_type != TYPE_UNDEFINED) {
            return true;
        } else {
            create_error(ctx->ir, node->token, "Identifier (%.*s) has an unresolved type", ident->name.len, ident->name.ptr);
        }
    } else {
        create_error(ctx->ir, node->token, "Unresolved reference to identifier (%.*s)", node->ident.len, node->ident.ptr);
    }
    return false;
}

static bool static_check_assignment(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ASSIGNMENT && node->children && md_array_size(node->children) == 2);
    ast_node_t* lhs = node->children[0];
    ast_node_t* rhs = node->children[1];

    ASSERT(lhs && lhs->ident.ptr);

    // We don't statically check lhs here since it is an assignment and the identifiers type is deduced from rhs.
    if (lhs->type != AST_IDENTIFIER) {
        create_error(ctx->ir, node->token, "Left hand side of assignment is not an identifier");
        return false;
    }

    if (static_check_node(rhs, ctx)) {
        ASSERT(lhs->data.type.base_type == TYPE_UNDEFINED);  // Identifiers type should always be undefined until explicitly assigned.
        ASSERT(rhs->data.type.base_type != TYPE_UNDEFINED);  // Right hand side type must be known

        identifier_t* ident = get_identifier(ctx->ir, lhs->ident);
        ASSERT(ident);

        if (!(rhs->flags & FLAG_DYNAMIC) && !rhs->data.ptr) {
            // If it is not a dynamic node, we evaluate it directly and store the data.
            allocate_data(&rhs->data, ctx->alloc);
            evaluate_node(&rhs->data, rhs, ctx);
            rhs->flags |= FLAG_CONSTANT;
        }

        node->data   = rhs->data;
        node->flags  = rhs->flags;

        lhs->data    = rhs->data;
        lhs->flags   = rhs->flags;

        ident->data  = rhs->data;
        ident->flags = rhs->flags;
        return true;
    }

    return false;
}

static void propagate_context(ast_node_t* node, const md_exp_bitfield_t* context) {
    ASSERT(node);
    node->context = context;
    for (int64_t i = 0; i < md_array_size(node->children); ++i) {
        propagate_context(node->children[i], context);
    }
}

static bool static_check_context(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_CONTEXT && node->children && md_array_size(node->children) == 2);

    // The rules are as follows:
    // A context RHS must be a bitfield type.
    //     
    // If the the basetype of LHS is (float, int, float[4], bool or whatever), then the result is of the same type but with an extra dimension with the same length as the RHS.
    // I.e.
    // If RHS is bitfield with 8 matching structures:
    // An LHS of float[4] would then result in a final type of float[4][8], one for each evaluated bitrange.
    
    // WE MUST BE ABLE TO EVALUATE RHS EXPRESSION DURING COMPILE TIME IN ORDER TO KNOW THE LENGTH
    // If the RHS is a dynamic type, we need to push the 'length' type deduction to 'run-time'
    ast_node_t* lhs = node->children[0];
    ast_node_t* rhs = node->children[1];

    if (static_check_node(rhs, ctx)) {
        if (rhs->data.type.base_type == TYPE_BITFIELD) {
            if (!(rhs->flags & FLAG_DYNAMIC)) {
                const int64_t num_contexts = type_info_array_len(rhs->data.type);
                if (num_contexts > 0) {
                    // Store this persistently and set this as a context for all child nodes
                    // Store as md_array so we can query its length later
                    md_exp_bitfield_t* contexts = 0;
                    md_array_resize(contexts, num_contexts, ctx->ir->arena);
                    for (int64_t i = 0; i < num_contexts; ++i) {
                        md_bitfield_init(&contexts[i], ctx->ir->arena);
                    }
                    rhs->data.ptr = contexts;
                    rhs->data.size = num_contexts * sizeof(md_exp_bitfield_t);

                    ASSERT(md_array_size(contexts) == num_contexts);
                    //allocate_data(&rhs->data, ctx->ir->arena);
                    if (evaluate_node(&rhs->data, rhs, ctx)) {

                        // We differentiate here if the LHS is a bitfield or not
                        // if LHS is a bitfield of length 1, the resulting bitfield length is the same as RHS (M)
                        // if LHS is a bitfield of length N, the resulting bitfield length is (N*M), N may vary for each context.
                        
                        node->lhs_context_types = 0;
                        int64_t arr_len = 0;
                        for (int64_t i = 0; i < num_contexts; ++i) {
                            eval_context_t local_ctx = {
                                .ir = ctx->ir,
                                .mol = ctx->mol,
                                .mol_ctx = &contexts[i],
                                .temp_alloc = ctx->temp_alloc
                            };
                            if (!static_check_node(lhs, &local_ctx)) {
                                return false;
                            }

                            md_type_info_t type = lhs->data.type;

                            if (lhs->flags & FLAG_DYNAMIC_LENGTH) {
                                if (!finalize_type(&type, lhs, &local_ctx)) {
                                    create_error(ctx->ir, lhs->token, "Failed to deduce length of type which is required for determining type and size of context");
                                    return false;
                                }
                            }

                            int64_t len = type_info_array_len(type);
                            ASSERT(len > 0);

                            md_array_push(node->lhs_context_types, type, ctx->ir->arena);
                            arr_len += len;
                        }

                        node->flags |= lhs->flags & FLAG_PROPAGATION_MASK;
                        node->data.type = lhs->data.type;
                        node->data.type.dim[node->data.type.len_dim] = (int32_t)arr_len;

                        propagate_context(lhs, contexts);

                        return true;

                    } else {
                        create_error(ctx->ir, node->token, "Right hand side of 'in' failed to evaluate at compile time.");
                    }
                } else {
                    create_error(ctx->ir, node->token, "The context is empty.");
                }
            } else {
                create_error(ctx->ir, node->token, "Right hand side of 'in' cannot be evaluated at compile time.");
            }
        } else {
            char buf[128] = {0};
            print_type_info(buf, ARRAY_SIZE(buf), rhs->data.type);
            create_error(ctx->ir, node->token, "Right hand side of keyword 'in' has an incompatible type: Expected bitfield, got '%s'.", buf);
        }
    }
    return false;
}

static bool static_check_node(ast_node_t* node, eval_context_t* ctx) {
    // This is probably more like a static check which encompass more than just checking the types...
    // The idea here is that we have an syntax tree given by node and we want to map the function calls into concrete functions and check for the type.
    // The type will propagate back throughout the tree and we will remap nodes into concrete function calls.
    // This is not done in the parsing step as we want the full context before attempting this operation.

    ASSERT(node);

    switch (node->type) {
    case AST_ASSIGNMENT:
        return static_check_assignment(node, ctx);
    case AST_CONSTANT_VALUE:
        return static_check_constant_value(node, ctx);
    case AST_ARRAY:
        return static_check_array(node, ctx);
    case AST_ARRAY_SUBSCRIPT:
        return static_check_array_subscript(node, ctx);
    case AST_IDENTIFIER:
        return static_check_identifier_reference(node, ctx);
    case AST_PROC_CALL:
        return static_check_proc_call(node, ctx);
    case AST_CONTEXT:
        return static_check_context(node, ctx);
    case AST_ADD:
    case AST_SUB:
    case AST_MUL:
    case AST_DIV:
    case AST_AND:
    case AST_OR:
    case AST_NOT:
    case AST_EQ:
    case AST_NE:
    case AST_LE:
    case AST_GE:
    case AST_LT:
    case AST_GT:
        return static_check_operator(node, ctx);
    case AST_CAST: // Should never happen since we insert casts here in type checking phase.
    default:
        ASSERT(false);
    }

    return false;
}

// Prunes the tree from AST_EXPRESSION, since that is just an unnecessary indirection for sorting the tree correctly on precedence
static ast_node_t* prune_expressions(ast_node_t* node) {
    if (node) {
        for (int64_t i = 0; i < md_array_size(node->children); ++i) {
            if (node->children[i]) {
                ASSERT(node->children[i] != node);
                node->children[i] = prune_expressions(node->children[i]);
            }
        }
        if (node->type == AST_EXPRESSION) {
            ASSERT(node->children && md_array_size(node->children) == 1);
            return node->children[0];
        }
    }
    return node;
}

static bool parse_script(md_script_ir_o* ir) {
    ASSERT(ir);

    bool result = true;

    tokenizer_t tokenizer = tokenizer_init(ir->str);

    parse_context_t ctx = {
        .ir = ir,
        .tokenizer = &tokenizer,
        .node = 0,
        .temp_alloc = default_temp_allocator,
    };

    ir->stage = "Parsing";

    // Parse statements until we have exhausted the tokenizer

    token_t tok = {0};
    while (tok = tokenizer_peek_next(ctx.tokenizer), tok.type != TOKEN_END) {
        const char* beg = tok.str.ptr;

        ir->record_errors = true;   // We reset the error recording flag before each statement
        ast_node_t* node = parse_expression(&ctx);
        node = prune_expressions(node);
        if (node) {
            tok = tokenizer_consume_next(ctx.tokenizer);
            if (tok.type != ';') {
                create_error(ir, tok, "Missing ';' to mark end of statement.");
                continue;
            }

            const char* end = tok.str.ptr + tok.str.len;
            str_t expr_str = {beg, (uint64_t)(end - beg)};

            if (node->type == AST_ASSIGNMENT &&
                md_array_size(node->children) == 2 && node->children[0]->type == AST_IDENTIFIER && node->children[0]->ident.ptr && node->children[1]) {
                identifier_t* ident = get_identifier(ir, node->children[0]->ident);
                ASSERT(ident);
                expression_t* expr = md_alloc(ir->arena, sizeof(expression_t));
                memset(expr, 0, sizeof(expression_t));
                expr->ident = ident;
                expr->node = node;
                expr->str = expr_str;
                md_array_push(ir->expressions, expr, ir->arena);
            }
            else {
                create_error(ir, node->token, "Every parsed expression should end up being assigned to a unique identifier.");
                result = false;
            }
        } else {
            token_type_t types[1] = {';'};
            tokenizer_consume_until_type(ctx.tokenizer, types, ARRAY_SIZE(types)); // Goto next expression
            tokenizer_consume_next(ctx.tokenizer); 
            result = false;
        }
    }

    return result;
}

static bool static_type_check(md_script_ir_o* ir, const md_molecule_t* mol) {
    ASSERT(ir);
    ASSERT(mol);

    bool result = true;

    const int64_t mem_size = MEGABYTES(64);
    void* mem = md_alloc(default_allocator, mem_size);

    md_stack_allocator_t stack_alloc = {0};
    md_stack_allocator_init(&stack_alloc, mem, mem_size);

    md_allocator_i temp_alloc = md_stack_allocator_create_interface(&stack_alloc);

    eval_context_t static_check_ctx = {
        .ir = ir,
        .temp_alloc = &temp_alloc,
        .alloc = ir->arena,
        .mol = mol,
    };
    ir->stage = "Static Type Checking";

    for (int64_t i = 0; i < md_array_size(ir->expressions); ++i) {
        ir->record_errors = true;   // We reset the error recording flag before each statement
        if (static_check_node(ir->expressions[i]->node, &static_check_ctx)) {
            md_array_push(ir->type_checked_expressions, ir->expressions[i], ir->arena);
        } else {
            result = false;
        }
    }

    md_free(default_allocator, mem, mem_size);

    return result;
}

static bool extract_dynamic_evaluation_targets(md_script_ir_o* ir) {
    ASSERT(ir);

    // Check all type checked expressions for dynamic flag
    for (int64_t i = 0; i < md_array_size(ir->type_checked_expressions); ++i) {
        expression_t* expr = ir->type_checked_expressions[i];
        ASSERT(expr);
        ASSERT(expr->node);
        if (expr->node->flags & FLAG_DYNAMIC) {
            md_array_push(ir->eval_targets, ir->type_checked_expressions[i], ir->arena);
        }
    }
    return true;
}

static inline bool is_temporal_type(md_type_info_t ti) {
    return (ti.base_type == TYPE_FLOAT && ti.dim[0] > 0 && ti.dim[1] == 0);
}

static inline bool is_distribution_type(md_type_info_t ti) {
    return (ti.base_type == TYPE_FLOAT && ti.dim[0] == DIST_BINS && ti.dim[1] > 0 && ti.dim[2] == 0);
}

static inline bool is_volume_type(md_type_info_t ti) {
    return (ti.base_type == TYPE_FLOAT && ti.dim[0] == VOL_DIM && ti.dim[1] == VOL_DIM && ti.dim[2] == VOL_DIM && ti.dim[3] > 0);
}

static inline bool is_property_type(md_type_info_t ti) {
    return is_temporal_type(ti) || is_distribution_type(ti) || is_volume_type(ti);
}

static bool extract_property_expressions(md_script_ir_o* ir) {
    ASSERT(ir);

    for (int64_t i = 0; i < md_array_size(ir->eval_targets); ++i) {
        expression_t* expr = ir->eval_targets[i];
        ASSERT(expr);
        ASSERT(expr->node);
        if (is_property_type(expr->node->data.type)) {
            md_array_push(ir->prop_expressions, expr, ir->arena);
        }
    }

    return true;
}

//static bool static_eval_node(ast_node_t* node, eval_context_t* ctx);

/*
static bool static_eval_context_node(ast_node_t* node, eval_context_t* ctx) {
    // We need to make a special case here because if we want to statically evaluate a node which is a context node.
    // Then we properly need to supply a context i.e. mol_context and all that before evaluating those nodes...

    ASSERT(ctx);
    ASSERT(node && node->type == AST_CONTEXT);
    ASSERT(md_array_size(node->children) == 2);

    const ast_node_t* expr_node = node->children[0];
    const ast_node_t* ctx_node  = node->children[1];

    if (!static_eval_node(ctx_node, ctx)) return false;

    // These should be available by now...
    ASSERT(ctx_node->data.type.base_type == TYPE_BITFIELD);
    ASSERT(ctx_node->data.ptr);

    const int64_t num_ctx = type_info_array_len(ctx_node->data.type);
    const md_exp_bitfield_t* ctx_bf = (const md_exp_bitfield_t*)ctx_node->data.ptr;

    ASSERT(node->lhs_context_types);
    const md_type_info_t* lhs_types = node->lhs_context_types;
    ASSERT(md_array_size(lhs_types) == num_ctx);

    // Evaluate the expression within each context.
    for (int64_t i = 0; i < num_ctx; ++i) {
        eval_context_t sub_ctx = *ctx;
        sub_ctx.mol_ctx = &ctx_bf[i];
        if (!static_eval_node(expr_node, &sub_ctx)) return false;
    }

    return true;
}
*/

static bool static_eval_node(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(ctx);
    ASSERT(ctx->ir);

    const int64_t num_children = md_array_size(node->children);
    bool result = true;
    
    // Propagate the evaluation all the way down to the children.
    if (node->children) {
        int64_t offset = 0;
        if (node->type == AST_CONTEXT) offset = 1;  // We don't want to evaluate the lhs in a context on its own since then it will loose its context
        for (int64_t i = offset; i < num_children; ++i) {
            result &= static_eval_node(node->children[i], ctx);
        }
    }

    // Only evaluate the node if it is not flagged as dynamic
    // Only evaluate if data.ptr is not already set (which can happen during static check)
    if (!(node->flags & FLAG_DYNAMIC) && !node->data.ptr) {
        if (allocate_data(&node->data, ctx->ir->arena)) {
            result &= evaluate_node(&node->data, node, ctx);
            if (node->type == AST_CONTEXT) {
                ASSERT(node->children);
                // If its a context node, we copy the data to the child as well
                node->children[0]->data = node->data;
            }
        } else {
            create_error(ctx->ir, node->token, "Could not allocate data for node during static evaluation");
        }
    }

    return result;
}

static bool static_evaluation(md_script_ir_o* ir, const md_molecule_t* mol) {
    ASSERT(mol);

    const int64_t mem_size = MEGABYTES(64);
    void* mem = md_alloc(default_allocator, mem_size);

    md_stack_allocator_t stack_alloc = {0};
    md_stack_allocator_init(&stack_alloc, mem, mem_size);

    md_allocator_i temp_alloc = md_stack_allocator_create_interface(&stack_alloc);

    eval_context_t ctx = {
        .ir = ir,
        .mol = mol,
        .temp_alloc = &temp_alloc,
    };

    ir->stage = "Static Evaluation";
    bool result = true;

    // Evaluate every node which is not flagged with FLAG_DYNAMIC and store its value.
    for (int64_t i = 0; i < md_array_size(ir->type_checked_expressions); ++i) {
        result &= static_eval_node(ir->type_checked_expressions[i]->node, &ctx);
        md_stack_allocator_reset(&stack_alloc);
    }

    md_free(default_allocator, mem, mem_size);

    return result;
}

static bool allocate_property_data(md_script_property_t* prop, md_type_info_t type, int64_t num_frames, md_allocator_i* alloc) {
    int64_t num_values = 0;
    memcpy(prop->data.dim, type.dim, sizeof(type.dim));

    int64_t aggregate_size = 0;

    switch (prop->type) {
    case MD_SCRIPT_PROPERTY_TYPE_TEMPORAL:
        // For temporal data, we store and expose all values, this enables filtering to be performed afterwards to create distributions
        prop->data.dim[1] = (int32_t)num_frames;
        num_values = prop->data.dim[0] * prop->data.dim[1];
        aggregate_size = (int32_t)num_frames;
        break;
    case MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION:
        // For distributions we only store and expose the aggregate
        num_values = DIST_BINS;
        prop->data.dim[0] = DIST_BINS;
        aggregate_size = num_values;
        break;
    case MD_SCRIPT_PROPERTY_TYPE_VOLUME:
        // For volumes we only store and expose the aggregate
        num_values = VOL_DIM * VOL_DIM * VOL_DIM;
        prop->data.dim[0] = VOL_DIM;
        prop->data.dim[1] = VOL_DIM;
        prop->data.dim[2] = VOL_DIM;
        aggregate_size = num_values;
        break;
    default:
        ASSERT(false);
    }

    prop->data.values = md_alloc(alloc, num_values * sizeof(float));
    memset(prop->data.values, 0, num_values * sizeof(float));
    prop->data.num_values = num_values;

    bool aggregate = type_info_array_len(type) > 1;
    if (aggregate) {
        prop->data.aggregate = md_alloc(alloc, sizeof(md_script_property_data_aggregate_t));

        memset(prop->data.aggregate, 0, sizeof(md_script_property_data_aggregate_t));
        prop->data.aggregate->num_values = aggregate_size;

        if (aggregate_size != num_values) {
            prop->data.aggregate->mean = md_alloc(alloc, aggregate_size * sizeof(float));
            memset(prop->data.aggregate->mean, 0, aggregate_size * sizeof(float));
        } else {
            prop->data.aggregate->mean = prop->data.values;
        }

        prop->data.aggregate->variance = md_alloc(alloc, aggregate_size * sizeof(float));
        memset(prop->data.aggregate->variance, 0, aggregate_size * sizeof(float));
    }

    return true;
}

static void create_properties(md_script_property_t** properties, int64_t num_frames, expression_t** expressions, int64_t num_expressions, md_allocator_i* alloc) {
    ASSERT(properties);

    md_script_property_t* props = NULL;

    for (int64_t i = 0; i < num_expressions; ++i) {
        const expression_t* expr = expressions[i];
        ASSERT(expr);
        ASSERT(expr->node);

        md_script_property_t prop = {
            .ident = copy_str(expr->ident->name, alloc),
            .type = 0,
            .unit = expr->node->data.unit,
            .data = {0},
            .vis_token = (struct md_script_vis_token_t*)expr->node,
        };

        if (is_temporal_type(expr->node->data.type)) {
            prop.type = MD_SCRIPT_PROPERTY_TYPE_TEMPORAL;
        }
        else if (is_distribution_type(expr->node->data.type)) {
            prop.type = MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION;
        }
        else if (is_volume_type(expr->node->data.type)) {
            prop.type = MD_SCRIPT_PROPERTY_TYPE_VOLUME;
        }
        else {
            ASSERT(false);
        }

        allocate_property_data(&prop, expr->node->data.type, num_frames, alloc);
        md_array_push(props, prop, alloc);
    }

    *properties = props;
}

static void compute_min_max(float* min, float* max, const float* values, int64_t num_values) {
    ASSERT(min);
    ASSERT(max);
    ASSERT(values);

    *min = FLT_MAX;
    *max = -FLT_MAX;
    for (int64_t i = 0; i < num_values; ++i) {
        *min = MIN(*min, values[i]);
        *max = MAX(*max, values[i]);
    }
}

static void compute_distribution(float* bins, int64_t num_bins, float min_bin_range, float max_bin_range, const float* values, int64_t num_values) {
    ASSERT(bins);
    ASSERT(num_bins > 0);
    ASSERT(max_bin_range >= min_bin_range);
    ASSERT(values);

    memset(bins, 0, num_bins * sizeof(float));

    const float bin_range = max_bin_range - min_bin_range;
    ASSERT(bin_range > 0.0f);
    for (int64_t i = 0; i < num_values; ++i) {
        const float val = values[i];
        const int32_t idx = (int32_t)CLAMP(((val - min_bin_range) / bin_range) * num_bins, 0, num_bins-1);
        ASSERT(0 <= idx && idx < num_bins);
        bins[idx] += 1;
    }

    const float scl = 1.0f / (float)num_values;
    for (int64_t i = 0; i < num_bins; ++i) {
        bins[i] *= scl;
    }
}

static bool eval_properties(md_script_property_t* props, int64_t num_props, const md_molecule_t* mol, const md_trajectory_i* traj, const md_exp_bitfield_t* mask, md_script_ir_o* ir, md_script_eval_o* eval_o) {
    ASSERT(props);
    ASSERT(mol);
    ASSERT(traj);
    ASSERT(ir);
    ASSERT(eval_o);
    
    const int64_t num_expr = md_array_size(ir->eval_targets);
    expression_t** const expr = ir->eval_targets;
    
    ASSERT(md_array_size(ir->prop_expressions) == num_props);

    const int64_t mem_size = MEGABYTES(128);
    void* mem = md_alloc(default_allocator, mem_size);

    md_stack_allocator_t stack_alloc = {0};
    md_stack_allocator_init(&stack_alloc, mem, mem_size);
    md_allocator_i temp_alloc = md_stack_allocator_create_interface(&stack_alloc);

    // coordinate data for reading trajectory frames into
    const int64_t stride = ROUND_UP(mol->atom.count, md_simd_width);    // Round up allocation size to simd width to allow for vectorized operations
    const int64_t coord_bytes = stride * 3 * sizeof(float);
    float* init_coords = md_alloc(&temp_alloc, coord_bytes);
    float* curr_coords = md_alloc(&temp_alloc, coord_bytes);

    md_exp_bitfield_t tmp_bf = {0};
    if (!mask) {
        md_bitfield_init(&tmp_bf, &temp_alloc);
        md_bitfield_set_range(&tmp_bf, 0, md_trajectory_num_frames(traj));
        mask = &tmp_bf;
    }
    
    // This data is meant to hold the evaluated expressions
    data_t* data = md_stack_allocator_alloc(&stack_alloc, num_expr * sizeof(data_t));

    const uint64_t STACK_RESET_POINT = md_stack_allocator_get_offset(&stack_alloc);

    md_trajectory_frame_header_t init_header = {0};

    float* init_x = init_coords + 0 * stride;
    float* init_y = init_coords + 1 * stride;
    float* init_z = init_coords + 2 * stride;

    float* curr_x = curr_coords + 0 * stride;
    float* curr_y = curr_coords + 1 * stride;
    float* curr_z = curr_coords + 2 * stride;

    // We want a mutable molecule which we can modify the atom coordinate section of
    md_molecule_t mutable_mol = *mol;
    mutable_mol.atom.x = curr_x;
    mutable_mol.atom.y = curr_y;
    mutable_mol.atom.z = curr_z;

    eval_context_t ctx = {
        .ir = ir,
        .mol = &mutable_mol,
        .temp_alloc = &temp_alloc,
        .alloc = &temp_alloc,
        .initial_configuration = {
            .header = &init_header,
            .x = init_x,
            .y = init_y,
            .z = init_z,
        },
    };

    // Fill the data for the initial configuration (needed by rmsd and SDF as a 'reference')
    md_trajectory_load_frame(traj, 0, &init_header, init_x, init_y, init_z);

    bool result = true;

    const int64_t num_evaluated_frames = md_bitfield_popcount(mask);

    int64_t beg_bit = mask->beg_bit;
    int64_t end_bit = mask->end_bit;

    int64_t dst_idx = 0;

    // Preprocess the property data
    for (int64_t p_idx = 0; p_idx < num_props; ++p_idx) {
        md_script_property_t* prop = &props[p_idx];
        switch (prop->type) {
        case MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION:
        case MD_SCRIPT_PROPERTY_TYPE_VOLUME:
            memset(prop->data.values, 0, prop->data.num_values * sizeof(float));
            if (prop->data.aggregate) {
                memset(prop->data.aggregate->mean, 0, prop->data.aggregate->num_values * sizeof(float));
                memset(prop->data.aggregate->variance, 0, prop->data.aggregate->num_values * sizeof(float));
            }
        case MD_SCRIPT_PROPERTY_TYPE_TEMPORAL:
            break;
        default:
            ASSERT(false);
        }
    }

    // We evaluate each frame, one at a time
    while ((beg_bit = md_bitfield_scan(mask, beg_bit, end_bit)) != 0) {
        if (eval_o->interrupt) {
            goto done;
        }

        int64_t f_idx = beg_bit - 1;
        ASSERT(f_idx < md_trajectory_num_frames(traj));
        
        result = md_trajectory_load_frame(traj, f_idx, NULL, curr_x, curr_y, curr_z);

        if (!result) {
            md_printf(MD_LOG_TYPE_ERROR, "Something went wrong when loading the frames during evaluation");
            goto done;
        }
        
        md_stack_allocator_set_offset(&stack_alloc, STACK_RESET_POINT);
        ctx.identifiers = 0;

        for (int64_t i = 0; i < num_expr; ++i) {
            data[i] = expr[i]->node->data;
            data[i].ptr = 0;
            allocate_data(&data[i], &temp_alloc);
            data[i].min_range = -FLT_MAX;
            data[i].max_range = FLT_MAX;
            if (!evaluate_node(&data[i], expr[i]->node, &ctx)) {
                // @TODO: Report error
                md_printf(MD_LOG_TYPE_ERROR, "Evaluation error when evaluating property '%.*s' at frame %lli", expr[i]->ident->name.len, expr[i]->ident->name.ptr, f_idx);
                result = false;
                goto done;
            }
        }

        for (int64_t p_idx = 0; p_idx < num_props; ++p_idx) {
            md_script_property_t* prop = &props[p_idx];
            const int32_t size = (int32_t)type_info_array_len(data[p_idx].type);
            float* values = (float*)data[p_idx].ptr;

            if (prop->type == MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) {
                ASSERT(prop->data.values);
                
                if (prop->data.aggregate) {
                    // Compute mean + variance values for this frame
                    ASSERT(prop->data.aggregate->mean);
                    ASSERT(prop->data.aggregate->variance);

                    float mean = 0;
                    for (int32_t i = 0 ; i < size; ++i) {
                        mean += values[i];
                    }
                    mean /= (float)size;

                    float variance = 0;
                    for (int32_t i = 0 ; i < size; ++i) {
                        const float x = values[i] - mean;
                        variance += x * x;
                    }
                    variance /= (float)size;
                    prop->data.aggregate->mean[dst_idx] = mean;
                    prop->data.aggregate->variance[dst_idx] = variance;
                }
                memcpy(prop->data.values + dst_idx * size, values, size * sizeof(float));
            }
            else if (prop->type == MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION) {
                // Accumulate values
                ASSERT(prop->data.values);

                if (prop->data.aggregate) {
                    // Compute mean + variance values for this frame
                    ASSERT(prop->data.aggregate->mean);
                    ASSERT(prop->data.aggregate->variance);

                    float mean[DIST_BINS] = {0};
                    for (int32_t i = 0 ; i < size; ++i) {
                        for (int32_t j = 0; j < DIST_BINS; ++j) {
                            mean[j] += values[i * DIST_BINS + j];                    
                        }
                    }
                    for (int32_t i = 0; i < DIST_BINS; ++i) {
                        mean[i] /= (float)size;
                    }

                    float variance[DIST_BINS] = {0};
                    for (int32_t i = 0 ; i < size; ++i) {
                        for (int32_t j = 0; j < DIST_BINS; ++j) {
                            const float x = values[i * DIST_BINS + j] - mean[j];
                            variance[j] += x * x;
                        }
                    }
                    for (int32_t i = 0; i < DIST_BINS; ++i) {
                        variance[i] /= (float)size;
                    }

                    for (int32_t i = 0; i < prop->data.num_values; ++i) {
                        prop->data.aggregate->mean[i] = mean[i];
                        prop->data.aggregate->variance[i] = variance[i];
                    }
                }
                else {
                    const md_simd_typef scl = md_simd_set1f(1.0f / (float)num_evaluated_frames);
                    for (int64_t i = 0; i < ROUND_UP(prop->data.num_values, md_simd_width); i += md_simd_width) {
                        md_simd_typef old_val = md_simd_loadf(prop->data.values + i);
                        md_simd_typef new_val = md_simd_mulf(md_simd_loadf(values + i), scl);
                        md_simd_storef(prop->data.values + i, md_simd_addf(old_val, new_val));
                    }
                }
            }
            else if (prop->type == MD_SCRIPT_PROPERTY_TYPE_VOLUME) {
                // Accumulate values
                ASSERT(prop->data.values);
                const md_simd_typef scl = md_simd_set1f(1.0f / (float)num_evaluated_frames);
                for (int64_t i = 0; i < ROUND_UP(prop->data.num_values, md_simd_width); i += md_simd_width) {
                    md_simd_typef old_val = md_simd_loadf(prop->data.values + i);
                    md_simd_typef new_val = md_simd_mulf(md_simd_loadf(values + i), scl);
                    md_simd_storef(prop->data.values + i, md_simd_addf(old_val, new_val));
                }
                /*for (int64_t i = 0; i < prop->data.num_values; ++i) {
                    prop->data.values[i] += values[i];
                }
                */
                // @TODO: Compute variance here as well for volumes
            }
            else {
                ASSERT(false);
            }
        }

        dst_idx += 1;
    }

    ASSERT(dst_idx == num_evaluated_frames);

    // Postprocess the data
    for (int64_t p_idx = 0; p_idx < num_props; ++p_idx) {
        md_script_property_t* prop = &props[p_idx];

        if (prop->type == MD_SCRIPT_PROPERTY_TYPE_TEMPORAL) {
            // This needs to be set to a lower number if we are evaluating a subset
            prop->data.num_values = num_evaluated_frames * prop->data.dim[0];

            // This is for the mean if the property is an aggregate
            compute_min_max(&prop->data.min_value, &prop->data.max_value, prop->data.values, prop->data.num_values);
            prop->data.min_range[0] = data[p_idx].min_range == -FLT_MAX ? prop->data.min_value : data[p_idx].min_range;
            prop->data.max_range[0] = data[p_idx].max_range ==  FLT_MAX ? prop->data.max_value : data[p_idx].max_range;
        }
        else if (prop->type == MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION) {
        }
        else if (prop->type == MD_SCRIPT_PROPERTY_TYPE_VOLUME) {
        }
        else {
            ASSERT(false);
        }
    }

done:
    md_free(default_allocator, mem, mem_size);

    return result;
}

static inline bool validate_ir_o(const md_script_ir_o* o) {
    if (!o) {
        md_print(MD_LOG_TYPE_ERROR, "Script object is NULL.");
        return false;
    }

    if (o->magic != SCRIPT_IR_MAGIC) {
        md_print(MD_LOG_TYPE_ERROR, "Script object is corrupt or invalid.");
        return false;
    }

    return true;
}

static inline bool validate_eval_o(const md_script_eval_o* o) {
    if (!o) {
        md_print(MD_LOG_TYPE_ERROR, "Eval object is NULL.");
        return false;
    }

    if (o->magic != SCRIPT_EVAL_MAGIC) {
        md_print(MD_LOG_TYPE_ERROR, "Eval object is corrupt or invalid.");
        return false;
    }

    return true;
}

static md_script_ir_o* create_ir(md_allocator_i* alloc) {
    md_script_ir_o* ir = md_alloc(alloc, sizeof(md_script_ir_o));
    memset(ir, 0, sizeof(md_script_ir_o));
    ir->magic = SCRIPT_IR_MAGIC;
    ir->arena = md_arena_allocator_create(alloc, MEGABYTES(1));
    
    return ir;
}

static bool init_ir(md_script_ir_o* o, str_t str) {
    if (validate_ir_o(o)) {
        md_arena_allocator_reset(o->arena);
        md_allocator_i* arena = o->arena;
        memset(o, 0, sizeof(md_script_ir_o));
        o->magic = SCRIPT_IR_MAGIC;
        o->arena = arena;
        o->str = copy_str(str, o->arena);
        return true;
    }
    return false;
}

static bool free_ir(md_script_ir_o* o) {
    if (validate_ir_o(o)) {
        md_arena_allocator_destroy(o->arena);
        memset(o, 0, sizeof(md_script_ir_o));
        return true;
    }
    return false;
}

static void create_tokens(md_script_ir_o* ir, const ast_node_t* node, const ast_node_t* node_override, int32_t depth) {
    md_script_token_t token = {0};

    ASSERT(node->type != AST_UNDEFINED);

    int len = 0;
    char buf[256] = {0};

    if (!(node->flags & FLAG_DYNAMIC)) {
        if (node->data.type.base_type != TYPE_BITFIELD) {
            char val_buf[128] = {0};
            int val_len = print_value(val_buf, sizeof(val_buf), node->data);
            len += snprintf(buf + len, sizeof(buf) - len, "%.*s ", val_len, val_buf);
        }
    } else {
        len += snprintf(buf + len, sizeof(buf) - len, "[d] ");
    }

    char type_buf[128] = {0};
    int type_len = print_type_info(type_buf, (int)sizeof(type_buf), node->data.type);
    len += snprintf(buf + len, sizeof(buf) - len, "(%.*s)", type_len, type_buf);

    token.line = node->token.line;
    token.col_beg = node->token.col_beg;
    token.col_end = node->token.col_end;
    token.depth = depth;
    token.vis_token = (const struct md_script_vis_token_t*)(node_override ? node_override : node);
    token.text = copy_str((str_t){.ptr = buf, .len = len}, ir->arena);

    // Parent token should be last token added (unless empty)
    md_script_token_t* last_token = md_array_last(ir->tokens);
    if (last_token && token.line == last_token->line && token.col_beg == last_token->col_beg && token.col_end == last_token->col_end) {
        // If the current token has the exact same line, col_beg and col_end, we only want to update the text, not the payload (the node)
        // since we will evaluate that top down anyways. Otherwise, visualizations for position extraction casts will not be triggered.
        last_token->text = token.text;
    } else {
        md_array_push(ir->tokens, token, ir->arena);
    }

    // Recurse and add children
    const int64_t num_children = md_array_size(node->children);
    if (node->type == AST_CONTEXT) {
        ASSERT(num_children == 2);
        create_tokens(ir, node->children[0], node, depth + 1);
        create_tokens(ir, node->children[1], NULL, depth + 1);
    } else {
        for (int64_t i = 0; i < num_children; ++i) {
            create_tokens(ir, node->children[i], NULL, depth + 1);
        }
    }
}

bool extract_tokens(md_script_ir_o* ir) {
    const int64_t num_expr = md_array_size(ir->type_checked_expressions);
    for (int64_t i = 0; i < num_expr; ++i) {
        expression_t* expr = ir->type_checked_expressions[i];
        create_tokens(ir, expr->node, NULL, 0);
    }

    return true;
}

bool md_script_ir_compile(md_script_ir_t* ir, md_script_ir_compile_args_t args) {
    ASSERT(ir);
    if (!args.src.ptr) {
        md_print(MD_LOG_TYPE_ERROR, "Script Compile: Source string was null");
        return false;
    }

    if (!args.mol) {
        md_print(MD_LOG_TYPE_ERROR, "Script Compile: Molecule was null");
        return false;
    }

    if (!args.alloc) {
        md_print(MD_LOG_TYPE_ERROR, "Script Compile: Allocator was null");
        return false;
    }

    if (!ir->o) {
        ir->o = create_ir(args.alloc);
    }
    init_ir(ir->o, args.src);

    bool result = parse_script(ir->o) &&
        static_type_check(ir->o, args.mol) &&
        static_evaluation(ir->o, args.mol) &&
        extract_dynamic_evaluation_targets(ir->o) &&
        extract_property_expressions(ir->o) &&
        extract_tokens(ir->o);

    ir->errors = ir->o->errors;
    ir->num_errors = md_array_size(ir->o->errors);

    ir->tokens = ir->o->tokens;
    ir->num_tokens = md_array_size(ir->o->tokens);

    ir->fingerprint = generate_fingerprint();

#if MD_DEBUG
    //save_expressions_to_json(ir->o->expressions, md_array_size(ir->o->expressions), make_cstr("tree.json"));
#endif
    return result;
}

bool md_script_ir_free(md_script_ir_t* ir) {
    bool result = false;
    if (ir) {
        result = ir->o ? free_ir(ir->o) : true;
        memset(ir, 0, sizeof(md_script_ir_t));
    }
    return result;
}

static struct md_script_eval_o* create_eval(md_allocator_i* alloc) {
    md_script_eval_o* o = md_alloc(alloc, sizeof(md_script_eval_o));
    o->magic = SCRIPT_EVAL_MAGIC;
    o->arena = md_arena_allocator_create(alloc, MEGABYTES(1));
    o->properties = 0;
    return o;
}

static bool free_eval(md_script_eval_o* o) {
    if (validate_eval_o(o)) {
        md_arena_allocator_destroy(o->arena);
        memset(o, 0, sizeof(md_script_eval_o));
        return true;
    }
    return false;
}

static bool init_eval(md_script_eval_o* o) {
    if (validate_eval_o(o)) {
        md_arena_allocator_reset(o->arena);
        o->properties = NULL;
        return true;
    }
    return false;
}

bool md_script_eval_init(md_script_eval_t* eval, int64_t num_frames, const md_script_ir_t* ir, md_allocator_i* alloc) {
    ASSERT(eval);
    ASSERT(ir);
    ASSERT(alloc);

    // Set these to zero until we have initialized the memory for the new properties
    eval->properties = NULL;
    eval->num_properties = 0;

    if (!eval->o) {
        eval->o = create_eval(alloc);
    }
    init_eval(eval->o);

    const int64_t num_prop_expr = md_array_size(ir->o->prop_expressions);
    if (num_prop_expr > 0) {
        create_properties(&eval->properties, num_frames, ir->o->prop_expressions, num_prop_expr, eval->o->arena);
        eval->num_properties = num_prop_expr;
    }

    return true;
}

bool md_script_eval_compute(md_script_eval_t* eval, md_script_eval_args_t args) {
    ASSERT(eval);

    bool result = true;
    if (!args.ir) {
        md_print(MD_LOG_TYPE_ERROR, "Script eval: Immediate representation was null");
        result = false;
    }
    if (!args.mol) {
        md_print(MD_LOG_TYPE_ERROR, "Script eval: Molecule was null");
        result = false;
    }
    if (!args.traj) {
        md_print(MD_LOG_TYPE_ERROR, "Script eval: Trajectory was null");
        result = false;
    }
    if (md_trajectory_num_frames(args.traj) == 0) {
        md_print(MD_LOG_TYPE_ERROR, "Script eval: Trajectory was empty");
        result = false;
    }

    if (result) {    
        if (eval->num_properties > 0) {
            ASSERT(eval->num_properties == md_array_size(args.ir->o->prop_expressions));
            eval->o->interrupt = false;
            eval->fingerprint = generate_fingerprint();
            result = eval_properties(eval->properties, eval->num_properties, args.mol, args.traj, args.filter_mask, args.ir->o, eval->o);
            eval->fingerprint = generate_fingerprint();
        }
    }

    return result;
}

bool md_script_eval_free(md_script_eval_t* eval) {
    bool result = false;
    if (eval) {
        result = eval->o ? free_eval(eval->o) : true;
        memset(eval, 0, sizeof(md_script_eval_t));
    }
    return result;
}

int md_script_eval_completed_frame_count(const md_script_eval_t* eval) {
    if (eval && eval->o)
        return eval->o->num_frames_completed;
    return 0;
}

void md_script_eval_interrupt(const md_script_eval_t* eval) {
    if (eval && eval->o) {
        eval->o->interrupt = true;
    }
}

static bool eval_expression(data_t* dst, str_t expr, md_molecule_t* mol, md_allocator_i* alloc) {
    md_script_ir_o* ir = create_ir(default_temp_allocator);
    init_ir(ir, expr);
    tokenizer_t tokenizer = tokenizer_init(ir->str);
    ast_node_t* node = parse_expression(&(parse_context_t){ .ir = ir, .tokenizer = &tokenizer, .temp_alloc = default_temp_allocator});
    if (node) {
        eval_context_t ctx = {
            .ir = ir,
            .mol = mol,
            .temp_alloc = default_temp_allocator,
        };

        if (static_check_node(node, &ctx)) {
            dst->type = node->data.type;
            allocate_data(dst, alloc);
            if (evaluate_node(dst, node, &ctx)) {
                return true;
            }
        }
    }

    if (ir->errors) {
        for (int64_t i = 0; i < md_array_size(ir->errors); ++i) {
            md_printf(MD_LOG_TYPE_ERROR, "%.*s", ir->errors[i].error.len, ir->errors[i].error.ptr);
        }
    }

    return false;
}

bool md_filter_evaluate(str_t expr, md_exp_bitfield_t* target, const md_filter_context_t* filter_ctx, md_filter_additional_info_t* info) {
    ASSERT(target);
    ASSERT(filter_ctx);
    ASSERT(filter_ctx->mol);
    ASSERT(filter_ctx->alloc);

    bool result = false;

    if (info && info->error_buf) {
        memset(info->error_buf, 0, info->error_cap);
    }

    md_allocator_i* alloc = md_arena_allocator_create(filter_ctx->alloc, MEGABYTES(1));
    md_allocator_i* temp_alloc = alloc;

    md_script_ir_o* ir = create_ir(alloc);
    init_ir(ir, expr);

    tokenizer_t tokenizer = tokenizer_init(ir->str);

    parse_context_t parse_ctx = {
        .ir = ir,
        .tokenizer = &tokenizer,
        .node = 0,
        .temp_alloc = alloc,
    };

    // @TODO: Find a way to copy the AST tree from the ir within filter_ctx and use that
    if (filter_ctx->ir && filter_ctx->ir->o) {
        if (filter_ctx->ir->o->identifiers) {
            md_array_push_array(ir->identifiers, filter_ctx->ir->o->identifiers, md_array_size(filter_ctx->ir->o->identifiers), alloc);
        }
        if (filter_ctx->ir->o->eval_targets) {
            md_array_push_array(ir->eval_targets, filter_ctx->ir->o->eval_targets, md_array_size(filter_ctx->ir->o->eval_targets), alloc);
        }
    }

    ir->stage = "Filter evaluate";
    ir->record_errors = true;

    ast_node_t* node = parse_expression(&parse_ctx);
    node = prune_expressions(node);
    if (node) {
        eval_context_t static_ctx = {
            .ir = ir,
            .mol = filter_ctx->mol,
            .temp_alloc = temp_alloc,
        };

        for (int64_t i = 0; i < filter_ctx->selection.count; ++i) {
            const md_filter_stored_selection_t sel = filter_ctx->selection.ptr[i];
            identifier_t ident = {
                .name = sel.ident,
                .data = {
                    .ptr = (void*)(&sel.bitfield),      // We are casting away const here, but internally this data is never modified because of the FLAG_CONSTANT
                    .size = sizeof(md_exp_bitfield_t),
                    .type = {.base_type = TYPE_BITFIELD, .dim = {1}},
                },
                .flags = FLAG_CONSTANT
            };
            md_array_push(static_ctx.identifiers, ident, alloc);
        }
        if (static_check_node(node, &static_ctx)) {
            if (node->data.type.base_type == TYPE_BITFIELD) {
                md_bitfield_clear(target);

                eval_context_t eval_ctx = {
                    .ir = ir,
                    .mol = filter_ctx->mol,
                    .temp_alloc = temp_alloc,
                    .initial_configuration = {
                        .x = filter_ctx->mol->atom.x,
                        .y = filter_ctx->mol->atom.y,
                        .z = filter_ctx->mol->atom.z,
                    }
                };

                // Evaluate dynamic targets if available
                for (int64_t i = 0; i < md_array_size(ir->eval_targets); ++i) {
                    data_t data = { .type = ir->eval_targets[i]->node->data.type };
                    allocate_data(&data, alloc);
                    evaluate_node(&data, ir->eval_targets[i]->node, &eval_ctx);
                }

                data_t data = { .type = node->data.type };
                allocate_data(&data, alloc);

                result = evaluate_node(&data, node, &eval_ctx);
                if (result) {
                    const int64_t len = type_info_array_len(data.type);
                    const md_exp_bitfield_t* bf_arr = data.ptr;
                    if (bf_arr) {
                        for (int64_t i = 0; i < len; ++i) {
                            md_bitfield_or(target, target, &bf_arr[i]);
                        }
                    }
                    if (info) {
                        if (info->is_dynamic) *info->is_dynamic = (bool)(node->flags & FLAG_DYNAMIC);
                    }
                }
            }
        }
    }
    if (!result) {
        if (info && info->error_buf) {
            snprintf(info->error_buf, info->error_cap, "Filter did not evaluate to a bitfield\n");
        } else {
            md_print(MD_LOG_TYPE_ERROR, "Filter did not evaluate to a bitfield");
        }
    }

    int64_t len = 0;
    for (int64_t i = 0; i < md_array_size(ir->errors); ++i) {
        if (info && info->error_buf) {
            int64_t space_left = MAX(0, info->error_cap - len);
            if (space_left) len += snprintf(info->error_buf + len, (size_t)space_left, "%.*s", (int)ir->errors[i].error.len, ir->errors[i].error.ptr);
        } else {
            md_printf(MD_LOG_TYPE_ERROR, "%.*s", ir->errors[i].error.len, ir->errors[i].error.ptr);
        }
    }

    md_arena_allocator_destroy(alloc);
    
    return result;
}

#define VIS_MAGIC 0xbc6169abd9628947

static md_range_t node_range(const ast_node_t* node) {
    md_range_t range = {node->token.col_beg, node->token.col_end};
    const int64_t num_children = md_array_size(node->children);
    ast_node_t** const children = node->children;
    for (int64_t i = 0; i < num_children; ++i) {
        md_range_t child_range = node_range(children[i]);
        range.beg = MIN(range.beg, child_range.beg);
        range.end = MAX(range.end, child_range.end);
    }

    return range;
}

static const ast_node_t* get_node(const ast_node_t* node, int32_t col) {
    ASSERT(node);

    const ast_node_t* result = 0;

    const int64_t num_children = md_array_size(node->children);
    if (num_children > 0) {
        for (int64_t i = 0; i < num_children; ++i) {
            md_range_t child_range = node_range(node->children[i]);
            if (child_range.beg <= col && col < child_range.end) {
                result = get_node(node->children[i], col);
            }
        }
    }

    if (!result) {
        md_range_t range = node_range(node);
        if (range.beg <= col && col < range.end) {
            result = node;
        }
    }

    return result;
}

static void do_vis_eval(const ast_node_t* node, eval_context_t* ctx) {
    if (node->data.type.base_type == TYPE_BITFIELD) {
        data_t data = { .type = node->data.type };
        allocate_data(&data, ctx->temp_alloc);
        evaluate_node(&data, node, ctx);

        ASSERT(data.ptr);
        const md_exp_bitfield_t* bf_arr = data.ptr;
        for (int64_t i = 0; i < element_count(data); ++i) {
            md_bitfield_or(&ctx->vis->atom_mask, &ctx->vis->atom_mask, &bf_arr[i]);
        }

        free_data(&data, ctx->temp_alloc);
    } else {
        evaluate_node(NULL, node, ctx);
    }
}

static void visualize_node(const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(ctx->vis);
    if (node->type != AST_CONTEXT && node->context) {
        for (int64_t i = 0; i < md_array_size(node->context); ++i) {
            ctx->mol_ctx = &node->context[i];
            do_vis_eval(node, ctx);
        }
    }
    else {
        do_vis_eval(node, ctx);
    }

    /*
    // Recurse and add children
    const int64_t num_children = md_array_size(node->children);
    for (int64_t i = 0; i < num_children; ++i) {
        if (node->type == AST_ASSIGNMENT && i == 0) continue;   // Skip LHS in assignment
        visualize_node(node->children[i], ctx);
    }
    */
}

bool md_script_visualization_init(md_script_visualization_t* vis, md_script_visualization_args_t args) {
    ASSERT(vis);
    if (!args.token) return false;
    if (!args.mol) return false;
    if (!args.ir) return false;
    if (!args.ir->o) return false;
    if (!args.alloc) return false;

    if (vis->o) {
        ASSERT(vis->o->magic == VIS_MAGIC);
        ASSERT(vis->o->alloc);
        md_arena_allocator_reset(vis->o->alloc);
    } else {
        md_allocator_i* arena = md_arena_allocator_create(args.alloc, MEGABYTES(1));
        vis->o = md_alloc(arena, sizeof(md_script_visualization_o));
        memset(vis->o, 0, sizeof(md_script_visualization_o));
        vis->o->alloc = arena;
        vis->o->magic = VIS_MAGIC;
    }

    const int64_t stack_size = MEGABYTES(32);
    void* stack_ptr = md_alloc(default_allocator, stack_size);
    md_stack_allocator_t stack = {0};
    md_stack_allocator_init(&stack, stack_ptr, stack_size);
    md_allocator_i stack_alloc = md_stack_allocator_create_interface(&stack);

    int64_t num_frames = md_trajectory_num_atoms(args.traj);
    float* x = md_stack_allocator_alloc(&stack, num_frames * sizeof(float));
    float* y = md_stack_allocator_alloc(&stack, num_frames * sizeof(float));
    float* z = md_stack_allocator_alloc(&stack, num_frames * sizeof(float));
    
    md_trajectory_frame_header_t header = { 0 };
    if (args.traj->inst) {
        md_trajectory_load_frame(args.traj, 0, &header, x, y, z);
    } else {
        x = args.mol->atom.x;
        y = args.mol->atom.y;
        z = args.mol->atom.z;
    }

    eval_context_t ctx = {
        .ir = args.ir->o,
        .mol = args.mol,
        .temp_alloc = &stack_alloc,
        .vis = vis,
        .initial_configuration = {
            .header = &header,
            .x = x,
            .y = y,
            .z = z
        }
    };
    
    md_bitfield_init(&vis->atom_mask, vis->o->alloc);

    visualize_node((ast_node_t*)args.token, &ctx);

    md_free(default_allocator, stack_ptr, stack_size);

    return true;
}

bool md_script_visualization_free(md_script_visualization_t* vis) {
    ASSERT(vis);
    if (vis->o) {
        ASSERT(vis->o->magic == VIS_MAGIC);
        md_arena_allocator_destroy(vis->o->alloc);
    }
    memset(vis, 0, sizeof(md_script_visualization_t));

    return true;
}
