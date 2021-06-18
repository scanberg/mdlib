#include "md_script.h"
#include "md_molecule.h"
#include "md_trajectory.h"
#include "md_filter.h"
#include "md_util.h"

#include "core/md_log.h"
#include "core/md_allocator.h"
#include "core/md_common.h"
#include "core/md_str.h"
#include "core/md_array.inl"
#include "core/md_file.h"
#include "core/md_arena_allocator.h"
#include "core/md_compiler.h"
#include "core/md_bitop.h"
#include "core/md_vec_math.h"

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#pragma warning(disable:4063) // Single character tokens not being valid tokens

#define SCRIPT_IR_MAGIC 0x371bacfe8274910a
#define SCRIPT_EVAL_MAGIC 0x89bca715287bcaff

#define MAX_SUPPORTED_PROC_ARGS 8
#define MAX_SUPPORTED_TYPE_DIMS 4

#define DIV_UP(x, y) ((x + (y-1)) / y) // Round up division

// #############################
// ###   TYPE DECLARATIONS   ###
// #############################

typedef enum token_type {
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

typedef enum ast_type {
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

typedef enum base_type {
    TYPE_UNDEFINED = 0,
    TYPE_FLOAT,
    TYPE_INT,
    TYPE_BOOL,
    TYPE_FRANGE,
    TYPE_IRANGE,
    TYPE_BITFIELD,      // Bitfield used to represent a selection of atoms, also contains a level which indicates the context (atom, residue, chain)
    TYPE_STRING,
} base_type_t;

typedef enum context_level {
    LEVEL_ATOM,
    LEVEL_RESIDUE,
    LEVEL_CHAIN,
} context_level_t;

typedef enum flags {
    FLAG_DYNAMIC                    = 0x001, // Indicates that it needs to be reevaluated for every frame of the trajectory (it has a dependence on atomic positions)
    FLAG_SYMMETRIC_ARGS             = 0x002, // Indicates that the arguments are symmetric
    FLAG_ARGS_EQUAL_LENGTH          = 0x004, // Arguments should have the same array length
    FLAG_RET_AND_ARG_EQUAL_LENGTH   = 0x008, // Return type array length matches argument arrays' length
    FLAG_RET_VARYING_LENGTH         = 0x010, // Return type has a varying length which is not constant over frames
    FLAG_QUERYABLE_LENGTH           = 0x020, // Marks a procedure as queryable, meaning it can be called with NULL as dst to query the length of the resulting array
    FLAG_VALIDATABLE_ARGS           = 0x040, // Marks a procedure as validatable, meaning it can be called with NULL as dst to validate the input arguments
    FLAG_POSITION                   = 0x080, // Hints that float3 argument(s) is atomic position(s), this means we can implicitly compute this from alot of input types using the position procedure
    FLAG_REF_FRAME                  = 0x100,
} flags_t;

typedef struct token {
    token_type_t type;
    str_t str;
    int64_t line;
    int64_t offset;
} token_t;

typedef struct frange {
    float beg;
    float end;
} frange_t;

typedef struct irange {
    int32_t beg;
    int32_t end;
} irange_t;

typedef struct mol_context {
    irange_t atom;
    irange_t residue;
    irange_t chain;
} mol_context_t;

typedef struct eval_context {
    struct md_script_ir_o* ir;
    token_t* arg_tokens; // Tokens to arguments for contextual information when reporting errors
    md_allocator_i* temp_alloc;
    const md_molecule* mol;
    mol_context_t mol_ctx;  // This gives us the subcontext ranges in which we operate
} eval_context_t;

typedef struct md_type_info {
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
    context_level_t level;                    // Is only used for bitfields to signify which context they operate on. Should always be 0 otherwise
} md_type_info_t;

// There is a conceptual twist here if len should belong to the type_info_t or if it should belong to data_t
// We need to differentiate between different lengths of arrays with the same type in our system, so for now it should probably reside within the type_info
// I.e. the array length should be part of the fundamental type. This only strenghtens the type-system

// This is the struct which is passed as arguments into the procedures
typedef struct data {
    md_type_info_t      type;
    void*               ptr;    // Pointer to the data
    int64_t             size;   // Size in bytes of data (This we want to determine during the static check so we can allocate the data when evaluating)
    md_script_unit_t    unit;

    // Only used for distributions
    float               min_range;
    float               max_range;
} data_t;

// This is data stored directly in the nodes to hold scalar values
typedef union value {
    str_t       _string;
    float       _float;
    int32_t     _int;
    frange_t    _frange;
    irange_t    _irange;
    bool        _bool;
    md_bitfield_t  _bitfield;
} value_t;

typedef struct identifier {
    str_t       name;
    data_t      data;
    flags_t     flags;
} identifier_t;

typedef struct procedure {
    str_t name;
    md_type_info_t return_type;
    uint32_t num_args;
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
    struct ast_node**    children;       // For AND, OR, NOT, procedure arguments etc.
    value_t         value;          // Scalar values for storing data directly inside the node
    data_t          data;           // Structure for passing as argument into procedure. (Holds pointer, length and type)

    procedure_t*    proc;           // Procedure reference
    identifier_t*   ident;          // Identifier reference
} ast_node_t;

typedef struct tokenizer {
    str_t str;
    int64_t cur;
    int64_t line;
} tokenizer_t;

typedef struct expression {
    ast_node_t*   node;     
    identifier_t* ident;    // Every expression should end up being assigned to an identifier. This is that identifier
    str_t str;
} expression_t;

typedef struct md_script_ir_o {
    uint64_t magic;

    // We could use a direct raw interface here to save some function pointer indirections
    struct md_allocator *arena;

    str_t str;  // Original string containing the 'source'
    
    // These are resizable arrays
    expression_t        **expressions;
    identifier_t        **identifiers;
    ast_node_t          **nodes;

    expression_t        **type_checked_expressions;     // List of expressions which passed type checking
    expression_t        **eval_targets;                 // List of dynamic expressions which needs to be evaluated per frame
    expression_t        **prop_expressions;             // List of expressions which are meant for exposure as properties

    // These are the final products which can be read through the public part of the structure
    md_script_error_t   *errors;

    bool record_errors; // This is to toggle if new errors should be registered... We don't want to flood the errors
    const char* stage;  // This is just to supply a context for the errors i.e. which stage the error occured
} md_script_ir_o;

typedef struct md_script_eval_o {
    uint64_t magic;

    struct md_allocator *arena;

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
    if (a.base_type == TYPE_BITFIELD && b.base_type == TYPE_BITFIELD) {
        if ((a.level == -1 || b.level == -1) || a.level == b.level) {
            return compare_type_info_dim(a, b);
        }
        return false;
    }
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

static inline bool is_incomplete_bitfield(md_type_info_t ti) {
    return ti.base_type == TYPE_BITFIELD && ti.level == -1;
}

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
    case TYPE_BITFIELD: return sizeof(md_bitfield_t);
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

static bool allocate_data(data_t* data, int64_t atom_count, md_allocator_i* alloc) {
    ASSERT(data && !is_undefined_type(data->type));
    ASSERT(alloc);
    
    // The convention should be that we only allocate data for known lengths
    const int64_t array_len = type_info_array_len(data->type);
    ASSERT(array_len > -1);

    // Do the base type allocation (array)
    const int64_t extra_bytes = (data->type.base_type == TYPE_BITFIELD) ? bitfield_byte_size(atom_count) * array_len : 0;
    const int64_t base_bytes = type_info_byte_stride(data->type) * array_len;
    const int64_t bytes = base_bytes + extra_bytes;

    data->ptr = md_alloc(alloc, bytes);
    if (!data->ptr) {
        return false;
    }
    data->size = bytes;
    memset(data->ptr, 0, bytes);

    if (data->type.base_type == TYPE_BITFIELD) {
        // Configure bitfields to point in memory
        md_bitfield_t* bf = data->ptr;
        const int64_t bf_data_stride = bitfield_byte_size(atom_count);
        const int64_t num_bits = atom_count;
        for (int64_t i = 0; i < array_len; ++i) {
            bf[i].bits = (uint64_t*)((char*)data->ptr + base_bytes + bf_data_stride * i);
            bf[i].num_bits = num_bits;
        }
    }

    return true;
}

static void free_data(data_t* data, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);
    if (data->ptr && data->size) md_free(alloc, data->ptr, data->size);
    data->ptr = 0;
    data->size = 0;
}

static inline void copy_data(data_t* dst, const data_t* src) {
    ASSERT(dst);
    ASSERT(src);
    ASSERT(src->size);
    ASSERT(src->ptr);
    ASSERT(dst->ptr);
    ASSERT(dst->size == src->size);

    dst->type = src->type;
    dst->unit = src->unit;
    memcpy(dst->ptr, src->ptr, src->size);

    if (dst->type.base_type == TYPE_BITFIELD) {
        const uint64_t num_elem = element_count(*dst);
        md_bitfield_t* dst_bf = dst->ptr;
        const md_bitfield_t* src_bf = src->ptr;
        // Set bit pointers into destination memory
        for (uint64_t i = 0; i < num_elem; ++i) {
            dst_bf->bits = (uint64_t*)((char*)dst->ptr + ((uint64_t)src_bf->bits - (uint64_t)src->ptr));
        }
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

static int64_t count_atoms_in_mask(md_bitfield_t mask, const md_molecule* mol) {
    ASSERT(mol);
    ASSERT(mask.num_bits == mol->atom.count);

    return (int64_t)bit_count(mask.bits, 0, mask.num_bits);
}

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

static uint64_t* get_atom_indices_in_mask(md_bitfield_t mask, const md_molecule* mol, md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);
    ASSERT(mask.num_bits == mol->atom.count);

    uint64_t* result = NULL;

    uint64_t idx = 0;
    while (idx = bit_scan(mask.bits, idx, mask.num_bits - idx)) {
        md_array_push(result, idx-1, alloc);
    }
    return result;
}

static uint64_t* get_residue_indices_in_mask(md_bitfield_t mask, const md_molecule* mol, md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);
    ASSERT(mask.num_bits == mol->atom.count);
    ASSERT(mol->atom.residue_idx);
    ASSERT(mol->residue.atom_range);

    uint64_t* result = NULL;

    uint64_t idx = 0;
    while ((idx = bit_scan(mask.bits, idx, mask.num_bits - idx)) != 0) {
        uint64_t res_idx = mol->atom.residue_idx[idx-1];
        md_array_push(result, res_idx, alloc);
        idx = mol->residue.atom_range[res_idx].end;
    }
    return result;
}

static uint64_t* get_chain_indices_in_mask(md_bitfield_t mask, const md_molecule* mol, md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);
    ASSERT(mask.num_bits == mol->atom.count);
    ASSERT(mol->atom.chain_idx);
    ASSERT(mol->chain.atom_range);

    uint64_t* result = NULL;

    uint64_t idx = 0;
    while ((idx = bit_scan(mask.bits, idx, mask.num_bits - idx)) != 0) {
        uint64_t chain_idx = mol->atom.chain_idx[idx-1];
        md_array_push(result, chain_idx, alloc);
        idx = mol->chain.atom_range[chain_idx].end;
    }
    return result;
}

static mol_context_t* get_residue_contexts_in_mask(md_bitfield_t mask, const md_molecule* mol, md_allocator_i* alloc) {
    ASSERT(alloc);
    ASSERT(mol);
    ASSERT(mol->residue.atom_range);    

    uint64_t* res_idx = get_residue_indices_in_mask(mask, mol, alloc);
    mol_context_t* res_ctx = NULL;

    for (int64_t i = 0; i < md_array_size(res_idx); ++i) {
        const uint64_t idx = res_idx[i];
        mol_context_t ctx = {
            .atom = {mol->residue.atom_range[idx].beg, mol->residue.atom_range[idx].end},
            .residue = {(uint32_t)idx, (uint32_t)idx+1},
            .chain = {0},
        };
        md_array_push(res_ctx, ctx, alloc);
    }

    return res_ctx;
}

static mol_context_t* get_chain_contexts_in_mask(md_bitfield_t mask, const md_molecule* mol, md_allocator_i* alloc) {
    ASSERT(alloc);
    ASSERT(mol);
    ASSERT(mol->chain.atom_range);
    ASSERT(mol->chain.residue_range);

    uint64_t* chain_idx = get_chain_indices_in_mask(mask, mol, alloc);
    mol_context_t* chain_ctx = NULL;

    for (int64_t i = 0; i < md_array_size(chain_idx); ++i) {
        const uint64_t idx = chain_idx[i];
        mol_context_t ctx = {
            .atom = {mol->chain.atom_range[idx].beg, mol->chain.atom_range[idx].end},
            .residue = {mol->chain.residue_range[idx].beg, mol->chain.residue_range[idx].end},
            .chain = {(uint32_t)idx, (uint32_t)idx+1},
        };
        md_array_push(chain_ctx, ctx, alloc);
    }

    return chain_ctx;
}



static void create_error(md_script_ir_o* ir, token_t token, const char* format, ...) {
    if (!ir->record_errors) return;
    ir->record_errors = false;

    char buffer[512] = {0};
    va_list args;
    va_start(args, format);
    int len = vsnprintf(buffer, ARRAY_SIZE(buffer), format, args);
    va_end(args);

    ASSERT(len > 0);
    ASSERT(len < ARRAY_SIZE(buffer));

    char* err_str = md_alloc(ir->arena, (uint64_t)len + 1);
    memcpy(err_str, buffer, len);
    err_str[len] = '\0';

    md_script_error_t error = {
        .line = (uint32_t)token.line,
        .offset = (uint32_t)(token.offset),
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
    const uint64_t num_identifiers = md_array_size(ir->identifiers);
    for (uint64_t i = 0; i < num_identifiers; ++i) {
        if (compare_str(ir->identifiers[i]->name, name)) {
            return ir->identifiers[i];
        }
    }
    return NULL;
}

static identifier_t* create_identifier(md_script_ir_o* ir, str_t name) {
    ASSERT(get_identifier(ir, name) == NULL);

    identifier_t* ident = md_alloc(ir->arena, sizeof(identifier_t));
    memset(ident, 0, sizeof(identifier_t));
    ident->name = copy_str(name, ir->arena);
    md_array_push(ir->identifiers, ident, ir->arena);
    return ident;
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
        if (!child || md_array_size(child->children) == 0) return;
        ast_node_t* grand_child = child->children[0];

        if (operator_precedence(child->type) > operator_precedence(parent->type)) {
            child->children[0] = parent;
            parent->children[i] = grand_child;
            *node = child;
        }
    }
}

static bool is_type_implicitly_convertible(md_type_info_t from, md_type_info_t to, flags_t flags) {
    flags &= FLAG_POSITION; // This is the only flag we care about from the perspective of matching casing procedures
    for (uint64_t i = 0; i < ARRAY_SIZE(casts); ++i) {
        if (flags && !(casts[i].flags & flags)) continue;
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

    for (uint64_t i = 0; i < ARRAY_SIZE(casts); ++i) {
        if (is_type_directly_compatible(casts[i].arg_type[0], from) &&
            is_type_directly_compatible(casts[i].return_type, to)) {
            res.success = true;
            res.procedure = &casts[i];
            break;
        }
    }

    return res;
}

static procedure_match_result_t find_procedure_supporting_arg_types_in_candidates(str_t name, const md_type_info_t arg_types[], uint64_t num_arg_types, procedure_t* candidates, uint64_t num_cantidates, bool allow_implicit_conversions) {
    procedure_match_result_t res = {0};

    for (uint64_t i = 0; i < num_cantidates; ++i) {
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

        for (uint64_t i = 0; i < num_cantidates; ++i) {
            procedure_t* proc = &candidates[i];
            if (compare_str(proc->name, name)) {
                if (num_arg_types == proc->num_args) {
                    // Same name and same number of args...

                    uint32_t cost = 0;
                    uint32_t flags = 0;
                    for (uint64_t j = 0; j < proc->num_args; ++j) {
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
                        for (uint64_t j = 0; j < proc->num_args; ++j) {
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
        token.line = tokenizer->line;
        token.offset = tokenizer->cur;
        return token;
    }
    
    int64_t line = tokenizer->line;
    const char* buf = tokenizer->str.ptr;
    const int64_t len = tokenizer->str.len;
    for (int64_t i = tokenizer->cur; i != len; ++i) {
        int64_t j = 0;
        
        if (buf[i] == '#') {
            // Skip to end of line
            while (++i != len) {
                if (buf[i] == '\n') {
                    ++line;
                    break;
                }
            }
        }
        
        else if (buf[i] == '\n') {
            ++line;
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
                for (int k = 0; k < ARRAY_SIZE(symbols_2); ++k) {
                    if (compare_str(str, (str_t){symbols_2[k].str, 2})) {
                        token.type = symbols_2[k].type;
                        j = i + 2;
                        break;
                    }
                }
            }
            if (!token.type) {
                for (int k = 0; k < ARRAY_SIZE(symbols_1); ++k) {
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
            token.line = line;
            token.offset = i;
            break;
        }

        if (!token.type && i == len - 1) {
            // We are going to hit the end in next iteration and have nothing to show for.
            tokenizer->cur = tokenizer->str.len;
            token.type = TOKEN_END;
            token.line = line;
            token.offset = i;
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

static token_t tokenizer_expect(tokenizer_t* tokenizer, token_type_t type) {
    token_t token = tokenizer_get_next_from_buffer(tokenizer);
    if (token.type == type) {
        if (token.str.ptr && token.str.len) {
            // Advance current of tokenizer
            tokenizer->cur = token.str.ptr + token.str.len - tokenizer->str.ptr;
            tokenizer->line = token.line;
            return token;
        }
        ASSERT(false);
    }
    ASSERT(false);
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
    int result = snprintf(buf, buf_size, "%s", get_value_type_str(info.base_type));
    int i = 0;
    while (info.dim[i] != 0 && i < MAX_SUPPORTED_TYPE_DIMS) {
        if (info.dim[i] == -1)
            result += snprintf(buf + result, buf_size, "[..]");
        else
            result += snprintf(buf + result, buf_size, "[%i]", (int)info.dim[i]);
        ++i;
    }
    return result;
}

#if MD_DEBUG

static void print_bitfield(FILE* file, md_bitfield_t bitfield) {
    const int64_t max_bits = 64;
    const int64_t num_bits = MIN(bitfield.num_bits, max_bits);
    fprintf(file, "[");
    for (int64_t i = 0; i < num_bits; ++i) {
        fprintf(file, "%i", bit_test(bitfield.bits, i) ? 1 : 0);
    }
    if (num_bits == max_bits) {
        fprintf(file, "...");
    }
    fprintf(file, "]");
}

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
                print_bitfield(file, *(md_bitfield_t*)data.ptr);
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

static void print_label(FILE* file, const ast_node_t* node) {
    switch(node->type) {
    case AST_CONSTANT_VALUE:
        print_value(file, node->data);
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
    {
        ASSERT(node->children);
        char buf[64] = {0};
        print_type_info(buf, ARRAY_SIZE(buf), node->children[0]->data.type);
        fprintf(file, "cast (%s) ->", buf);
        break;
    }
    case AST_IDENTIFIER:
        if (node->ident) {
            fprintf(file, "%.*s", (int)node->ident->name.len, node->ident->name.ptr);
        }
        else {
            fprintf(file, "NULL");
        }
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

    char buf[64] = {0};
    print_type_info(buf, ARRAY_SIZE(buf), node->data.type);
    fprintf(file, " -> %s", buf);

    if (node->flags) {
        fprintf(file, " {");
        if (node->flags & FLAG_DYNAMIC) fprintf(file, "D");
        fprintf(file, "}");
    }

    if (node->type != AST_ASSIGNMENT) {
        if (node->data.ptr) {
            fprintf(file, " = ");
            print_value(file, node->data);
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

static ast_node_t* parse_expression(parse_context_t* ctx);

static ast_node_t** parse_comma_separated_arguments_until_token(parse_context_t* ctx, token_type_t token_type) {
    ast_node_t** args = 0;
    token_t next = {0};
    while (next = tokenizer_peek_next(ctx->tokenizer), next.type != TOKEN_END) {
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
        if (expect_token_type(ctx->ir, tokenizer_consume_next(ctx->tokenizer), ')')) {
            node = create_node(ctx->ir, AST_PROC_CALL, token);
            const int64_t num_args = md_array_size(args);
            if (num_args) {
                md_array_push_array(node->children, args, num_args, ctx->ir->arena);
            }
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
        node->ident = c;
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
                node->ident = create_identifier(ctx->ir, ident);
            }
            else {
                create_error(ctx->ir, token, "The identifier is already taken. Variables cannot be reassigned.");
            }
        } else {
            // Identifier reference, resolve this later
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
    const int offset = (2 - num_operands);
    ast_node_t* node = 0;
    
    for (int i = offset; i < 2; ++i) {
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
    md_array_push_array(node->children, arg + offset, num_operands, ctx->ir->arena);
    
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
    
    if (token.type == ':' || next.type == ':' && (is_number(token.type))) {
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
        if (expect_token_type(ctx->ir, tokenizer_consume_next(ctx->tokenizer), ']')) {
            const int64_t num_elements = md_array_size(elements);
            if (num_elements) {

                ast_node_t* node_copy = create_node(ctx->ir, node->type, token);
                memcpy(node_copy, node, sizeof(ast_node_t));

                memset(node, 0, sizeof(ast_node_t));
                node->type = AST_ARRAY_SUBSCRIPT;
                node->token = token;

                md_array_push(node->children, node_copy, ctx->ir->arena);
                md_array_push_array(node->children, elements, num_elements, ctx->ir->arena);

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
    if (expect_token_type(ctx->ir, tokenizer_consume_next(ctx->tokenizer), '}')) {
        const int64_t num_elements = md_array_size(elements);
        if (num_elements) {
            // We only commit the results if everything at least parsed ok.
            node = create_node(ctx->ir, AST_ARRAY, token);
            md_array_push_array(node->children, elements, num_elements, ctx->ir->arena);
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
    token = tokenizer_consume_next(ctx->tokenizer);
    expect_token_type(ctx->ir, token, ')');
    return node;
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
    return ctx->node;
}

// ####################
// ###   EVALUATE   ###
// ####################

static bool finalize_type(md_type_info_t* type, const ast_node_t* node, eval_context_t* ctx);

static bool evaluate_node(data_t*, const ast_node_t*, eval_context_t*);

static bool evaluate_proc_call(data_t* dst, const ast_node_t* node, eval_context_t* eval_ctx) {
    ASSERT(node && node->type == AST_PROC_CALL);
    ASSERT(node->proc);

    // This procedure can be called within the static check to 'query' the result of the procedure.
    // In such case, dst will be NULL.

    const int64_t num_args = md_array_size(node->children);
    ASSERT(num_args < MAX_SUPPORTED_PROC_ARGS);

    data_t arg_data[MAX_SUPPORTED_PROC_ARGS] = {0};
    token_t arg_tokens[MAX_SUPPORTED_PROC_ARGS] = {0};

    for (int64_t i = 0; i < num_args; ++i) {
        // We need to evaluate the argument nodes first before we make the proc call.
        // In this context we are not interested in storing any data (since it is only used to be passed as arguments)
        // so we can allocate the data required for the node with the temp alloc
        arg_data[i].type = node->children[i]->data.type;
        arg_tokens[i] = node->children[i]->token;

        if (is_variable_length(arg_data[i].type)) {
            if (!finalize_type(&arg_data[i].type, node->children[i], eval_ctx)) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to finalize dynamic type in procedure call");
                return false;
            }
        }

        allocate_data(&arg_data[i], eval_ctx->mol->atom.count, eval_ctx->temp_alloc);
        if (!evaluate_node(&arg_data[i], node->children[i], eval_ctx))
            return false;
    }
    
    eval_context_t ctx = {
        .ir = eval_ctx->ir,
        .arg_tokens = arg_tokens,
        .mol = eval_ctx->mol,
        .mol_ctx = eval_ctx->mol_ctx,
        .temp_alloc = eval_ctx->temp_alloc
    };

    int result = node->proc->proc_ptr(dst, arg_data, &ctx);

    for (int64_t i = num_args - 1; i >= 0; --i) {
        free_data(&arg_data[i], eval_ctx->temp_alloc);
    }

    return result >= 0;
}

static bool evaluate_constant_value(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_CONSTANT_VALUE);
    ASSERT(node->data.ptr && node->data.size > 0);  // Assume that the static check set the data already
    ASSERT(dst && dst->ptr && dst->size == node->data.size); // Make sure we have room for the data
    (void)ctx;

    copy_data(dst, &node->data);

    return true;
}

static bool evaluate_identifier(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_IDENTIFIER);
    ASSERT(dst && dst->ptr && dst->size >= node->ident->data.size); // Make sure we have the expected size
    (void)ctx;
    
    // We assume in this context that we have a reference to the identifier, we copy the contents
    copy_data(dst, &node->ident->data);

    return true;
}

static bool evaluate_assignment(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ASSIGNMENT);
    ASSERT(node->children);

    ASSERT(md_array_size(node->children) == 2);
    ast_node_t* lhs = node->children[0];
    ast_node_t* rhs = node->children[1];

    ASSERT(lhs && lhs->type == AST_IDENTIFIER && lhs->ident);
    ASSERT(rhs);

    if (strncmp(lhs->ident->name.ptr, "p1", 2) == 0) {
        while(0);
    }

    // We know from the static check that the lhs (child[0]) is an identifier and the rhs is a valid expression with a known data size
    // the data size could be zero if it is an array of length 0, in that case we don't need to allocate anything.
    // This should also be the conceptual entry point for all evaluations in this declarative language

    if (evaluate_node(dst, rhs, ctx)) {
        // Good stuff
        lhs->data = *dst;
        lhs->ident->data = *dst;
        return true;
    }

    return false;
}

static bool evaluate_array(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY);
    ASSERT(dst && dst->ptr);
    ASSERT(compare_type_info(dst->type, node->data.type));

    const uint64_t elem_size = type_info_byte_stride(dst->type);
    const md_type_info_t elem_type = type_info_element_type(dst->type);

    for (int64_t i = 0; i < md_array_size(node->children); ++i) {
        data_t data = {elem_type, (char*)dst->ptr + elem_size * i, elem_size};
        if (!evaluate_node(&data, node->children[i], ctx)) {
            return false;
        }
    }

    return true;
}

static bool evaluate_array_subscript(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY_SUBSCRIPT);
    ASSERT(dst && dst->ptr);
    ASSERT(compare_type_info(dst->type, node->data.type));
    ASSERT(md_array_size(node->children) == 2);

    const ast_node_t* arr = node->children[0];
    const ast_node_t* idx = node->children[1];

    data_t arr_data = {arr->data.type};
    data_t idx_data = {idx->data.type};

    if (!allocate_data(&arr_data, ctx->mol->atom.count, ctx->temp_alloc)) return false;
    if (!allocate_data(&idx_data, ctx->mol->atom.count, ctx->temp_alloc)) return false;

    if (!evaluate_node(&arr_data, arr, ctx)) return false;
    if (!evaluate_node(&idx_data, idx, ctx)) return false;

    ASSERT(arr_data.ptr);
    ASSERT(idx_data.ptr);

    const int64_t elem_size = type_info_byte_stride(arr_data.type);
    int offset = 0;
    int length = 1;
    if (compare_type_info(idx_data.type, (md_type_info_t)TI_INT)) {
        offset = *((int*)(idx_data.ptr));
    } else if (compare_type_info(idx_data.type, (md_type_info_t)TI_IRANGE)) {
        irange_t rng = *((irange_t*)(idx_data.ptr));
        offset = rng.beg;
        length = rng.end - rng.beg + 1;
    } else {
        ASSERT(false);
    }

    // @TODO: Implement support for bitfield here.
    // For bitfields we operate on the level of the bitfield (atoms, residues, chains)
    // the index or range corresponds to the n'th atom / residue / chain
    // mask these with an bit_and operation

    ASSERT(elem_size == dst->size);
    const data_t src = {dst->type, (char*)arr_data.ptr + elem_size * offset, elem_size * length};
    copy_data(dst, &src);

    return true;
}

static bool evaluate_context(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_CONTEXT);
    ASSERT(md_array_size(node->children) == 2);

    const ast_node_t* expr_node = node->children[0];
    const ast_node_t* ctx_node  = node->children[1];

    data_t ctx_data  = {ctx_node->data.type};

    if (!allocate_data(&ctx_data,  ctx->mol->atom.count, ctx->temp_alloc)) return false;

    if (!evaluate_node(&ctx_data, ctx_node, ctx)) return false;
    ASSERT(ctx_data.type.base_type == TYPE_BITFIELD);
    ASSERT(ctx_data.ptr);

    md_bitfield_t ctx_bf = *((md_bitfield_t*)ctx_data.ptr);

    mol_context_t* mol_contexts = 0;

    switch (ctx_data.type.level) {
    case LEVEL_RESIDUE:
        mol_contexts = get_residue_contexts_in_mask(ctx_bf, ctx->mol, ctx->temp_alloc);
        break;
    case LEVEL_CHAIN:
        mol_contexts = get_chain_contexts_in_mask(ctx_bf, ctx->mol, ctx->temp_alloc);
        break;
    case LEVEL_ATOM:
    default:
        ASSERT(false);  // This should have been handeled by the static check
    }

    const uint64_t num_ctx = md_array_size(mol_contexts);
    const uint64_t dst_len = type_info_array_len(dst->type);

    if (dst->type.base_type == TYPE_BITFIELD) {
        ASSERT(dst_len == 1);
        // If destination type is bitfield, we know that the lhs in the context is also a bitfield
        // we evaluate each context and apply bitwise OR between all the results.

        data_t tmp_data = {
            .type = expr_node->data.type,
        };
        if (!allocate_data(&tmp_data,  ctx->mol->atom.count, ctx->temp_alloc)) return false;

        md_bitfield_t dst_bf = *((md_bitfield_t*)dst->ptr);
        md_bitfield_t tmp_bf = *((md_bitfield_t*)tmp_data.ptr);

        ASSERT(dst_bf.bits && dst_bf.num_bits);
        ASSERT(tmp_bf.bits && tmp_bf.num_bits == dst_bf.num_bits);

        // Evaluate the expression within each context.
        for (uint64_t i = 0; i < num_ctx; ++i) {
            eval_context_t sub_ctx = {
                .mol = ctx->mol,
                .mol_ctx = mol_contexts[i],
                .temp_alloc = ctx->temp_alloc,
            };
            if (!evaluate_node(&tmp_data, expr_node, &sub_ctx)) return false;
            bit_or(dst_bf.bits, dst_bf.bits, tmp_bf.bits, 0, dst_bf.num_bits);
        }

        free_data(&tmp_data, ctx->temp_alloc);
    }
    else {
        ASSERT(num_ctx == dst_len);

        const uint64_t elem_size = type_info_byte_stride(dst->type);
        const md_type_info_t elem_type = type_info_element_type(dst->type);

        // Evaluate the expression within each context.
        for (uint64_t i = 0; i < num_ctx; ++i) {
            eval_context_t sub_ctx = {
                .mol = ctx->mol,
                .mol_ctx = mol_contexts[i],
                .temp_alloc = ctx->temp_alloc,
            };
            data_t data = {
                .type = elem_type,
                .ptr = (char*)dst->ptr + elem_size * i,
                .size = elem_size
            };
            if (!evaluate_node(&data, expr_node, &sub_ctx)) return false;
            if (i == 0) {
                dst->min_range = data.min_range;
                dst->max_range = data.max_range;
            }
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
            return evaluate_identifier(dst, node, ctx);
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
            ASSERT(is_scalar(from));
            ASSERT(is_scalar(to));
            // Convert node directly
            switch (to.base_type) {
            case TYPE_FLOAT:
                ASSERT(from.base_type == TYPE_INT);
                node->value._float = (float)node->value._int;
                node->data.size = sizeof(node->value._float);
                break;
            case TYPE_IRANGE:
                ASSERT(from.base_type == TYPE_INT);
                node->value._irange = (irange_t){node->value._int, node->value._int};
                node->data.size = sizeof(node->value._irange);
                break;
            case TYPE_FRANGE:
                ASSERT(from.base_type == TYPE_IRANGE);
                node->value._frange = (frange_t){(float)node->value._irange.beg, (float)node->value._irange.end};
                node->data.size = sizeof(node->value._frange);
                break;
            default:
                ASSERT(false);
            }
            node->data.type.base_type = to.base_type;
            return true;
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
    ASSERT(is_variable_length(node->data.type));
    ASSERT(node->type == AST_PROC_CALL);
    ASSERT(type);

    *type = node->data.type;

    const ast_node_t** args = node->children;
    int64_t num_args = md_array_size(node->children);

    if (node->proc->flags & FLAG_RET_AND_ARG_EQUAL_LENGTH) {
        ASSERT(num_args > 0);
        // We can deduce the length of the array from the input type (take the array which is the longest???)
        int64_t max_len = 0;
        for (uint64_t i = 0; i < num_args; ++i) {
            int64_t len = (int32_t)type_info_array_len(args[i]->data.type);
            ASSERT(len > -1);
            max_len = MAX(len, max_len);
        }
        type->dim[type->len_dim] = (int32_t)max_len;
    } else {
        ASSERT(node->proc->flags & FLAG_QUERYABLE_LENGTH);

        // Find which dimension we should fill in
        uint64_t dim = 0;
        for (dim = 0; dim < MAX_SUPPORTED_TYPE_DIMS; ++dim) {
            if (node->data.type.dim[dim] == -1) break;
        }
        ASSERT(dim < MAX_SUPPORTED_TYPE_DIMS);

        eval_context_t eval_ctx = {
            .ir = ctx->ir,
            .mol = ctx->mol,
            .mol_ctx = ctx->mol_ctx,
            .temp_alloc = ctx->temp_alloc
        };

        data_t args[MAX_SUPPORTED_PROC_ARGS] = {0};

        for (uint64_t i = 0; i < num_args; ++i) {
            // We need to evaluate the argument nodes first before we make the proc call.
            // In this context we are not interested in storing any data (since it is only used to be passed as arguments)
            // so we can allocate the data required for the node with the temp alloc
            args[i].type = node->children[i]->data.type;
            allocate_data(&args[i], eval_ctx.mol->atom.count, eval_ctx.temp_alloc);
            if (!evaluate_node(&args[i], node->children[i], &eval_ctx)) {
                return false;
            }
        }

        // Perform the call
        int query_result = node->proc->proc_ptr(NULL, args, &eval_ctx);

        // Free
        for (int64_t i = num_args - 1 ; i >= 0; --i) {
            free_data(&args[i], eval_ctx.temp_alloc);
        }

        if (query_result >= 0) { // Zero length is valid
            type->dim[dim] = query_result;
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

    const uint64_t num_args = md_array_size(node->children);
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
        for (uint64_t i = 0; i < num_args; ++i) {
            if (!is_type_directly_compatible(arg[i]->data.type, node->proc->arg_type[i])) {
                // Types are not directly compatible, but should be implicitly convertible (otherwise the match should never have been found)
                ASSERT(is_type_implicitly_convertible(arg[i]->data.type, node->proc->arg_type[i], node->proc->flags));
                if (!convert_node(arg[i], node->proc->arg_type[i], ctx)) {
                    return false;
                }
                if (arg[i]->flags & FLAG_DYNAMIC) node->flags |= FLAG_DYNAMIC;
            }
        }

        if (node->proc->flags & FLAG_ARGS_EQUAL_LENGTH) {
            ASSERT(num_args > 1); // This flag should only be set for procedures with multiple arguments
                                  // Make sure that the input arrays are of equal length. Otherwise create an error.
            const uint64_t expected_len = type_info_array_len(arg[0]->data.type);
            for (uint64_t i = 1; i < num_args; ++i) {
                uint64_t len = type_info_array_len(arg[i]->data.type);
                if (len != expected_len) {
                    create_error(ctx->ir, node->token, "Expected array-length of arguments to match. arg 0 has length %i, arg %i has length %i.", (int)expected_len, (int)i, (int)len);
                    return false;
                }
            }
        }
    }

    // Make sure this flag is set in the tree so it can be propagated upwards later
    if (node->proc->flags & FLAG_DYNAMIC) node->flags |= FLAG_DYNAMIC;
    if (node->proc->flags & FLAG_RET_VARYING_LENGTH) node->flags |= FLAG_RET_VARYING_LENGTH;

    if (is_variable_length(node->data.type)) {
        // The nodes type has a variable length. We would like to resolve it statically in order to determine the type
        if (node->flags & FLAG_RET_VARYING_LENGTH) {
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
    const uint64_t num_children = md_array_size(node->children);
    for (uint64_t i = 0; i < num_children; ++i) {
        ASSERT(node != node->children[i]);
        if (static_check_node(node->children[i], ctx)) {
            node->flags |= node->children[i]->flags;
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
        const uint64_t num_args = md_array_size(node->children);
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

    if (node->token.line == 92) {
        while(0);
    }

    if (!node->proc) {
        if (static_check_children(node, ctx)) {
            const uint64_t num_args = md_array_size(node->children);
            const str_t proc_name = node->token.str;

            if (num_args == 0) {
                // Zero arguments
                procedure_match_result_t res = find_procedure_supporting_arg_types(proc_name, 0, 0, false);
                if (res.success) {
                    return finalize_proc_call(node, &res, ctx);
                }
                else {
                    create_error(ctx->ir, node->token,
                        "Could not find matching procedure '.*s' which takes no arguments", (int)proc_name.len, proc_name.ptr);
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
                    return finalize_proc_call(node, &res, ctx);
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
    } else {
        // We already have gotten the procedure now we just possibly want to statically check the function by querying
        if (node->proc->flags & FLAG_QUERYABLE_LENGTH) {
            return evaluate_proc_call(NULL, node, ctx);
        } else {
            return static_check_children(node, ctx);
        }
    }

    return false;
}

static bool static_check_constant_value(ast_node_t* node, eval_context_t* ctx) {
    (void)ctx;
    ASSERT(node && node->type == AST_CONSTANT_VALUE);
    ASSERT(node->data.type.base_type != TYPE_UNDEFINED);
    ASSERT(is_scalar(node->data.type));

    node->data.ptr = &node->value;
    node->data.size = base_type_element_byte_size(node->data.type.base_type);
    
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

        const uint64_t  num_elem = md_array_size(node->children);
        ast_node_t**        elem = node->children;

        md_type_info_t array_type = {0};
        
        // Pass 1: find a suitable base_type for the array
        for (uint64_t i = 0; i < num_elem; ++i) {
            if (array_type.base_type == TYPE_UNDEFINED) {
                // Assign type from first element
                array_type = elem[i]->data.type;
            }
            else {
                if (!is_type_directly_compatible(elem[i]->data.type, array_type)) {
                    if (is_type_implicitly_convertible(elem[i]->data.type, array_type, 0)) {
                        // We can convert the element to the array type.   
                    }
                    else if (is_type_implicitly_convertible(array_type, elem[i]->data.type, 0)) {
                        // Option 2: We can convert the array type into the elements type (e.g int to float)
                        // Promote the type of the array
                        array_type = elem[i]->data.type;

                        // Just to make sure, we recheck all elements up until this one
                        for (uint64_t j = 0; j < i; ++j) {
                            if (!is_type_directly_compatible(elem[i]->data.type, array_type) &&
                                !is_type_implicitly_convertible(elem[i]->data.type, array_type, 0)) {
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
        for (uint64_t i = 0; i < num_elem; ++i) {
            if (!is_type_directly_compatible(elem[i]->data.type, array_type) &&
                is_type_implicitly_convertible(elem[i]->data.type, array_type, 0)) {
                if (!convert_node(elem[i], array_type, ctx)) {
                    return false;
                }
            }
        }

        // Finalize the type for the array
        array_type.dim[array_type.len_dim] = (int32_t)num_elem;
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
                if (0 <= range.beg && range.end < element_count(lhs->data)) {
                    ASSERT(lhs->data.type.dim[0] > 0);
                    // SUCCESS!
                    node->data.type = lhs->data.type;
                    node->data.type.dim[node->data.type.len_dim] = range.end - range.beg + 1;
                    return true;
                } else {
                    create_error(ctx->ir, arg->token, "Array subscript out of range");
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
    
    if (!node->ident) {   
        node->ident = get_identifier(ctx->ir, node->token.str);
        if (node->ident) {
            node->data = node->ident->data;
            node->flags |= node->ident->flags;    // This means it has a dependency on an identifier, and that identifier might be dynamic... we need to propagate flags
        }
    }

    if (node->ident) {
        if (node->ident->data.type.base_type != TYPE_UNDEFINED) {
            return true;
        } else {
            create_error(ctx->ir, node->token, "Identifier (%.*s) has an unresolved type", node->ident->name.len, node->ident->name.ptr);
        }
    } else {
        create_error(ctx->ir, node->token, "Unresolved reference to identifier (%.*s)", node->token.str.len, node->token.str.ptr);
    }
    return false;
}

static bool static_check_assignment(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ASSIGNMENT && node->children && md_array_size(node->children) == 2);
    ast_node_t* lhs = node->children[0];
    ast_node_t* rhs = node->children[1];

    ASSERT(lhs && lhs->ident);

    // We don't statically check lhs here since it is an assignment and the identifiers type is deduced from rhs.
    if (lhs->type != AST_IDENTIFIER) {
        create_error(ctx->ir, node->token, "Left hand side of assignment is not an identifier");
        return false;
    }

    if (static_check_node(rhs, ctx)) {
        ASSERT(lhs->data.type.base_type == TYPE_UNDEFINED);  // Identifiers type should always be undefined until explicitly assigned.
        if (rhs->data.type.base_type != TYPE_UNDEFINED) {
            node->data.type         = rhs->data.type;
            lhs->data.type          = rhs->data.type;
            lhs->ident->data.type   = rhs->data.type;
            node->flags             = rhs->flags;
            lhs->flags              = rhs->flags;
            lhs->ident->flags       = rhs->flags;
            return true;
        } else {
            create_error(ctx->ir, node->token, "Right hand side of assignment has an unresolved type");
        }
    }

    return false;
}

static bool static_check_context(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_CONTEXT && node->children && md_array_size(node->children) == 2);

    // The rules are as follows:
    // A context RHS must be a bitfield type.
    // 
    // If the LHS is a bitfield, then we just modify the bits with an AND operation.
    // If the the basetype of LHS is (float, int, float[4], bool or whatever), then the result is of the same type but with an extra dimension with the same length as the RHS.
    // I.e.
    // If RHS is bitfield with 8 matching structures:
    // An LHS of float[4] would then result in a final type of float[4][8], one for each evaluated bitrange.
    
    // WE MUST BE ABLE TO EVALUATE RHS EXPRESSION DURING COMPILE TIME IN ORDER TO KNOW THE LENGTH

    if (node->token.line == 10) {
        while(0);
    }

    if (static_check_children(node, ctx)) {
        ast_node_t* lhs = node->children[0];
        ast_node_t* rhs = node->children[1];

        if (rhs->data.type.base_type == TYPE_BITFIELD && rhs->data.type.level > 0) {
            if (!(rhs->flags & FLAG_DYNAMIC)) {
                eval_context_t eval_ctx = {
                    .mol = ctx->mol,
                    .mol_ctx = {
                        .atom = {0, (uint32_t)ctx->mol->atom.count},
                        .residue = {0, (uint32_t)ctx->mol->residue.count},
                        .chain = {0, (uint32_t)ctx->mol->chain.count}
                    },
                    .temp_alloc = ctx->temp_alloc,
                };
                allocate_data(&rhs->data, ctx->mol->atom.count, ctx->temp_alloc);
                if (evaluate_node(&rhs->data, rhs, &eval_ctx)) {
                    md_bitfield_t bf = *((md_bitfield_t*)rhs->data.ptr);
                    mol_context_t* mol_contexts = NULL;
                    switch(rhs->data.type.level) {
                        case LEVEL_RESIDUE: mol_contexts = get_residue_contexts_in_mask(bf, ctx->mol, ctx->temp_alloc); break;
                        case LEVEL_CHAIN:   mol_contexts = get_chain_contexts_in_mask(bf, ctx->mol, ctx->temp_alloc); break;
                        default: ASSERT(false);
                    }

                    if (mol_contexts) {
                        // We differentiate here if the LHS is a bitfield or not
                        // If it is a bitfield, we don't create an array of bitfields, we just aim to OR them all into the same single bitfield
                        // This is why we don't set the length in that case

                        node->data.type = lhs->data.type;
                        if (lhs->data.type.base_type != TYPE_BITFIELD) {
                            node->data.type.dim[node->data.type.len_dim] = (int32_t)md_array_size(mol_contexts);
                        }

                        for (int64_t i = 0; i < md_array_size(mol_contexts); ++i) {
                            eval_context_t local_ctx = {
                                .ir = ctx->ir,
                                .mol = ctx->mol,
                                .mol_ctx = mol_contexts[i],
                                .temp_alloc = ctx->temp_alloc
                            };
                            if (!static_check_node(lhs, &local_ctx)) {
                                return false;
                            }
                        }

                        return true;
                    } else {
                        create_error(ctx->ir, node->token, "The context is empty.");
                    }
                } else {
                    create_error(ctx->ir, node->token, "Right hand side of 'in' failed to evaluate at compile time.");
                }
            } else {
                create_error(ctx->ir, node->token, "Right hand side of 'in' cannot be evaluated at compile time.");
            }
        } else {
            char buf[128] = {0};
            print_type_info(buf, ARRAY_SIZE(buf), rhs->data.type);
            create_error(ctx->ir, node->token, "Right hand side of in has an incompatible type, expected residue or chain level structure got (%s)", buf);
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
static void prune_ast_expressions(ast_node_t* node) {
    ASSERT(node);
    
    for (int64_t i = 0; i < md_array_size(node->children); ++i) {
        ast_node_t* child = node->children[i];
        prune_ast_expressions(child);

        if (child && child->type == AST_EXPRESSION) {
            ASSERT(child->children && md_array_size(child->children) == 1);
            node->children[i] = child->children[0];
        }
    }
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
        if (node) {
            tok = tokenizer_consume_next(ctx.tokenizer);
            if (tok.type != ';') {
                create_error(ir, tok, "Missing ';' to mark end of statement.");
                continue;
            }

            const char* end = tok.str.ptr + tok.str.len;
            str_t expr_str = {beg, (uint64_t)(end - beg)};

            if (node->type == AST_ASSIGNMENT &&
                md_array_size(node->children) == 2 && node->children[0]->type == AST_IDENTIFIER && node->children[0]->ident && node->children[1]) {
                prune_ast_expressions(node);
                identifier_t* ident = node->children[0]->ident;
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

static bool static_type_check(md_script_ir_o* ir, const md_molecule* mol) {
    ASSERT(ir);
    ASSERT(mol);

    bool result = true;

    eval_context_t static_check_ctx = {
        .ir = ir,
        .temp_alloc = default_temp_allocator,
        .mol = mol,
        .mol_ctx = {
            .atom = {0, (uint32_t)mol->atom.count},
            .residue = {0, (uint32_t)mol->residue.count},
            .chain = {0, (uint32_t)mol->chain.count},
        },
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

static bool static_eval_node(ast_node_t* node, md_script_ir_o* ir, eval_context_t* ctx) {
    ASSERT(node);

    bool result = true;

    if (node->flags & FLAG_DYNAMIC) {
        // This node cannot be resolved at compile time.
        // However, its children might be candidates for compile time evaluation
        // @TODO: Allocate the required data for the top level dynamic nodes here. (size * num_frames)
        if (node->children) {
            for (int64_t i = 0; i < md_array_size(node->children); ++i) {
                result &= static_eval_node(node->children[i], ir, ctx);
            }
        }
    } else {
        // Allocate persistent data for the 'top' node within the tree/subtree which is saved for later evaluations.
        // Use temp malloc for everything within the subtree
        // Constant value expressions will already have their ptr and size set
        if (!node->data.ptr) {
            if (allocate_data(&node->data, ctx->mol->atom.count, ir->arena)) {
                result &= evaluate_node(&node->data, node, ctx);
            } else {
                create_error(ir, node->token, "Could not allocate data for node during static evaluation");
            }
        }
    }

    return result;
}

static bool static_evaluation(md_script_ir_o* ir, const md_molecule* mol) {
    ASSERT(mol);
    eval_context_t ctx = {
        .mol = mol,
        .mol_ctx = {
            .atom = {0, (uint32_t)mol->atom.count},
            .residue = {0, (uint32_t)mol->residue.count},
            .chain = {0, (uint32_t)mol->chain.count}
        },
        .temp_alloc = default_temp_allocator,   // We might need a linear allocator here which can handle big data. Something which does not wrap around and we can reset after each evaluation.
    };

    ir->stage = "Static Evaluation";
    bool result = true;

    // Evaluate every node which is not flagged with FLAG_DYNAMIC and store its value.
    for (int64_t i = 0; i < md_array_size(ir->type_checked_expressions); ++i) {
        result &= static_eval_node(ir->type_checked_expressions[i]->node, ir, &ctx);
    }

    return result;
}



static bool allocate_property_data(md_script_property_t* prop, int64_t array_len, int64_t num_frames, md_allocator_i* alloc) {
    switch (prop->type) {
    case MD_SCRIPT_TYPE_TEMPORAL:
    {
        prop->temporal = md_alloc(alloc, sizeof(md_script_scalar_t));
        memset(prop->temporal, 0, sizeof(md_script_scalar_t));

        const int64_t num_values = num_frames;
        prop->temporal->values = md_alloc(alloc, num_values * sizeof(float));
        memset(prop->temporal->values, 0, num_values * sizeof(float));
        prop->temporal->num_values = num_values;

        if (array_len > 1) {
            prop->temporal->variance = md_alloc(alloc, num_values * sizeof(float));
            memset(prop->temporal->variance, 0, num_values * sizeof(float));
        }
        // FALLTHROUGH
    }
    case MD_SCRIPT_TYPE_DISTRIBUTION:
    {
        prop->distribution = md_alloc(alloc, sizeof(md_script_distribution_t));
        memset(prop->distribution, 0, sizeof(md_script_distribution_t));

        const int64_t num_values = DIST_BINS;
        prop->distribution->values = md_alloc(alloc, num_values * sizeof(float));
        memset(prop->distribution->values, 0, num_values * sizeof(float));
        prop->distribution->num_values = num_values;

        if (array_len > 1) {
            prop->distribution->variance = md_alloc(alloc, num_values * sizeof(float));
            memset(prop->distribution->variance, 0, num_values * sizeof(float));
        }
        break; // NO FALLTHROUGH
    }
    case MD_SCRIPT_TYPE_VOLUME:
    {
        prop->volume = md_alloc(alloc, sizeof(md_script_volume_t));
        memset(prop->volume, 0, sizeof(md_script_volume_t));

        const int64_t num_values = VOL_DIM * VOL_DIM * VOL_DIM;
        prop->volume->values = md_alloc(alloc, num_values * sizeof(float));
        memset(prop->volume->values, 0, num_values * sizeof(float));
        prop->volume->num_values = num_values;

        if (array_len > 1) {
            prop->volume->variance = md_alloc(alloc, num_values * sizeof(float));
            memset(prop->volume->variance, 0, num_values * sizeof(float));
        }
        break;
    }
    default:
        ASSERT(false);
    }

    return true;
}



static bool init_properties(md_script_property_t properties[], int64_t num_frames, const expression_t** expressions, int64_t num_expressions, md_allocator_i* alloc) {
    ASSERT(properties);
    for (int64_t i = 0; i < num_expressions; ++i) {
        expression_t* expr = expressions[i];
        ASSERT(expr);
        ASSERT(expr->node);

        md_script_property_t prop = {
            .ident = copy_str(expr->ident->name, alloc),
            .type = 0,
            .unit = expr->node->data.unit,
            .dim = {0},
        };

        if (is_temporal_type(expr->node->data.type)) {
            prop.type = MD_SCRIPT_TYPE_TEMPORAL;
        }
        else if (is_distribution_type(expr->node->data.type)) {
            prop.type = MD_SCRIPT_TYPE_DISTRIBUTION;
        }
        else if (is_volume_type(expr->node->data.type)) {
            prop.type = MD_SCRIPT_TYPE_VOLUME;
        }
        else {
            ASSERT(false);
        }

        ASSERT(ARRAY_SIZE(prop.dim) == ARRAY_SIZE(expr->node->data.type.dim));
        ASSERT(sizeof(prop.dim) == sizeof(expr->node->data.type.dim));
        memcpy(prop.dim, expr->node->data.type.dim, sizeof(prop.dim));

        const int64_t array_len = type_info_array_len(expr->node->data.type);
        allocate_property_data(&prop, array_len, num_frames, alloc),

        properties[i] = prop;
    }
    return true;
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

static void compute_distribution(md_script_distribution_t* dist, const float* values, int64_t num_values) {
    ASSERT(dist);
    ASSERT(values);

    const float bin_range = dist->max_range - dist->min_range;
    const int32_t num_bins = dist->num_values;
    ASSERT(bin_range > 0.0f);
    for (int64_t i = 0; i < num_values; ++i) {
        const float val = values[i];
        const int32_t bin = (val - dist->min_range) / bin_range;
        if (0 <= bin && bin < num_bins) {
            dist->values[bin] += 1;
        }
    }

    for (int64_t i = 0; i < num_bins; ++i) {
        dist->values[i] /= (float)num_values;
    }
}

static bool eval_properties(md_script_property_t* props, int64_t num_props, const md_molecule* mol, const md_trajectory_i* traj, md_script_ir_o* ir, md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(traj);
    ASSERT(ir);
    ASSERT(alloc);

    // Modifiable position xyz data for reading trajectory frames into
    const int64_t xyz_size = mol->atom.count * 3 * sizeof(float);
    float* xyz_ptr = md_alloc(alloc, xyz_size);

    // We want a version which we can modify the atom coordinate section of
    md_molecule mutable_mol = *mol;
    mutable_mol.atom.x = xyz_ptr;
    mutable_mol.atom.y = xyz_ptr + 1 * mol->atom.count;
    mutable_mol.atom.z = xyz_ptr + 2 * mol->atom.count;

    md_allocator_i* eval_arena = md_arena_allocator_create(alloc, KILOBYTES(16));

    eval_context_t ctx = {
        .mol = &mutable_mol,
        .mol_ctx = {
            .atom    = {0, (uint32_t)mutable_mol.atom.count},
            .residue = {0, (uint32_t)mutable_mol.residue.count},
            .chain   = {0, (uint32_t)mutable_mol.chain.count}
        },
        .temp_alloc = eval_arena
    };

    bool result = true;

    const int64_t num_expr = md_array_size(ir->eval_targets);
    const expression_t** expr = ir->eval_targets;

    ASSERT(md_array_size(ir->prop_expressions) == num_props);

    // We evaluate each frame, one at a time
    for (int64_t f_idx = 0; f_idx < traj->num_frames; ++f_idx) {
        md_trajectory_load_frame(traj, f_idx, MD_TRAJ_FIELD_XYZ, &(md_trajectory_data_t){
            .num_atoms = mutable_mol.atom.count,
            .x = mutable_mol.atom.x,
            .y = mutable_mol.atom.y,
            .z = mutable_mol.atom.z,
        });

        md_arena_allocator_reset(eval_arena);

        for (int64_t i = 0; i < num_expr; ++i) {
            allocate_data(&expr[i]->node->data, ctx.mol->atom.count, eval_arena);
            expr[i]->node->data.min_range = -FLT_MAX;
            expr[i]->node->data.max_range = FLT_MAX;
            if (!evaluate_node(&expr[i]->node->data, expr[i]->node, &ctx)) {
                // @TODO: Report error
                md_printf(MD_LOG_TYPE_ERROR, "Evaluation error when evaluating property '%.*s' at frame %lli", expr[i]->ident->name.len, expr[i]->ident->name.ptr, f_idx);
                result = false;
                goto done;
            }
        }

        for (int64_t p_idx = 0; p_idx < num_props; ++p_idx) {
            ast_node_t* node = ir->prop_expressions[p_idx]->node;
            md_script_property_t* prop = &props[p_idx];
            const int32_t size = type_info_array_len(node->data.type);
            float* values = (float*)node->data.ptr;

            if (prop->type == MD_SCRIPT_TYPE_TEMPORAL) {
                ASSERT(prop->temporal);
                ASSERT(prop->temporal->values);
                
                if (size > 1) {
                    // Compute mean + variance!
                    ASSERT(prop->temporal->variance);
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
                    prop->temporal->values[f_idx] = mean;
                    prop->temporal->variance[f_idx] = variance;
                }
                else {
                    prop->temporal->values[f_idx] = values[0];
                }
            }
            else if (prop->type == MD_SCRIPT_TYPE_DISTRIBUTION) {
                // Accumulate values
                ASSERT(prop->distribution);
                ASSERT(prop->distribution->values);

                if (size > 1) {
                    // Compute mean + variance!
                    ASSERT(prop->distribution->variance);
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

                    for (int32_t i = 0; i < prop->distribution->num_values; ++i) {
                        prop->distribution->values[i] = mean[i];
                        prop->distribution->variance[i] = variance[i];
                    }
                }
                else {
                    for (int64_t i = 0; i < prop->distribution->num_values; ++i) {
                        prop->distribution->values[i] += values[i];
                    }
                }
            }
            else if (prop->type == MD_SCRIPT_TYPE_VOLUME) {
                // Accumulate values
                ASSERT(prop->volume);
                ASSERT(prop->volume->values);
                for (int64_t i = 0; i < prop->volume->num_values; ++i) {
                    prop->volume->values[i] += values[i];
                }
            }
            else {
                ASSERT(false);
            }
        }
    }

    // Postprocess the data
    for (int64_t p_idx = 0; p_idx < num_props; ++p_idx) {
        md_script_property_t* prop = &props[p_idx];
        ast_node_t* node = ir->prop_expressions[p_idx]->node;

        if (prop->type == MD_SCRIPT_TYPE_TEMPORAL) {
            compute_min_max(&prop->temporal->min_value, &prop->temporal->max_value, prop->temporal->values, prop->temporal->num_values);
            ASSERT(prop->distribution);
            prop->distribution->min_range = node->data.min_range == -FLT_MAX ? prop->temporal->min_value : node->data.min_range;
            prop->distribution->max_range = node->data.max_range ==  FLT_MAX ? prop->temporal->max_value : node->data.max_range;

            compute_distribution(prop->distribution, prop->temporal->values, prop->temporal->num_values);
        }
        else if (prop->type == MD_SCRIPT_TYPE_DISTRIBUTION) {
            const double scl = 1.0 / (double)traj->num_frames;
            for (int64_t i = 0; i < prop->distribution->num_values; ++i) {
                prop->distribution->values[i] *= scl;
            }
        }
        else if (prop->type == MD_SCRIPT_TYPE_VOLUME) {
            const double scl = 1.0 / (double)traj->num_frames;
            for (int64_t i = 0; i < prop->volume->num_values; ++i) {
                prop->volume->values[i] *= scl;
            }
        }
        else {
            ASSERT(false);
        }
    }

done:
    md_arena_allocator_destroy(eval_arena);
    md_free(alloc, xyz_ptr, xyz_size);

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
    ir->arena = md_arena_allocator_create(alloc, KILOBYTES(16));
    
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

bool md_script_ir_compile(struct md_script_ir* ir, struct md_script_ir_compile_args args) {
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
        extract_property_expressions(ir->o);

    ir->errors = ir->o->errors;
    ir->num_errors = md_array_size(ir->o->errors);

#if MD_DEBUG
    //save_expressions_to_json(ir->expressions, md_array_size(ir->expressions), make_cstr("tree.json"));
#endif
    return result;
}

bool md_script_ir_free(struct md_script_ir* ir) {
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
    o->arena = md_arena_allocator_create(alloc, KILOBYTES(16));
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
        o->properties = 0;
        return true;
    }
    return false;
}

bool md_script_eval(struct md_script_eval_result* eval, struct md_script_eval_args args) {
    ASSERT(eval);
    if (!args.ir) {
        md_print(MD_LOG_TYPE_ERROR, "Script eval: Immediate representation was null");
        return false;
    }
    if (!args.mol) {
        md_print(MD_LOG_TYPE_ERROR, "Script eval: Molecule was null");
        return false;
    }
    if (!args.traj) {
        md_print(MD_LOG_TYPE_ERROR, "Script eval: Trajectory was null");
        return false;
    }
    if (!args.alloc) {
        md_print(MD_LOG_TYPE_ERROR, "Script eval: Allocator was null");
        return false;
    }

    if (!eval->o) {
        eval->o = create_eval(args.alloc);
    }
    init_eval(eval->o);

    bool result = false;
    
    const int64_t num_props = md_array_size(args.ir->o->prop_expressions);
    if (num_props) {
        eval->properties = md_alloc(args.alloc, num_props * sizeof(md_script_property_t));
        eval->num_properties = num_props;
        result = init_properties(eval->properties, args.traj->num_frames, args.ir->o->prop_expressions, num_props, eval->o->arena) &&
                 eval_properties(eval->properties, num_props, args.mol, args.traj, args.ir->o, eval->o->arena);
    }

    if (!result) {
        eval->properties = 0;
        eval->num_properties = 0;
    }

    return result;
}

bool md_script_eval_free(struct md_script_eval_result* eval) {
    bool result = false;
    if (eval) {
        result = eval->o ? free_eval(eval->o) : true;
        memset(eval, 0, sizeof(md_script_eval_result_t));
    }
    return result;
}

static bool eval_expression(data_t* dst, str_t expr, md_molecule* mol, md_allocator_i* alloc) {
    md_script_ir_o* ir = create_ir(default_temp_allocator);
    init_ir(ir, expr);
    tokenizer_t tokenizer = tokenizer_init(ir->str);
    ast_node_t* node = parse_expression(&(parse_context_t){ .ir = ir, .tokenizer = &tokenizer, .temp_alloc = default_temp_allocator});
    if (node) {
        eval_context_t ctx = {
            .ir = ir,
            .mol = mol,
            .mol_ctx = {
                .atom = {0, (uint32_t)mol->atom.count},
                .residue = {0, (uint32_t)mol->residue.count},
                .chain = {0, (uint32_t)mol->chain.count},
            },
            .temp_alloc = default_temp_allocator,
        };

        if (static_check_node(node, &ctx)) {
            dst->type = node->data.type;
            allocate_data(dst, mol->atom.count, alloc);
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

bool md_filter_evaluate(str_t expr, md_bitfield_t target, md_filter_context_t filter_ctx) {
    ASSERT(target.bits);
    ASSERT(target.num_bits);
    ASSERT(filter_ctx.mol);

    if (target.num_bits < filter_ctx.mol->atom.count) {
        md_print(MD_LOG_TYPE_ERROR, "Not enough bits to store result of filter evaluation");
        return false;
    }

    bool result = true;

    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, KILOBYTES(16));
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

    ir->stage = "Filter evaluate";
    ir->record_errors = true;

    ast_node_t* node = parse_expression(&parse_ctx);
    if (node) {
        prune_ast_expressions(node);
        identifier_t* stored_selections = 0;
        for (int64_t i = 0; i < filter_ctx.selection.count; ++i) {
            const md_filter_stored_selection_t sel = filter_ctx.selection.ptr[i];
            identifier_t ident = {
                .name = sel.ident,
                .data = {
                    .ptr = &sel.bitfield,
                    .size = DIV_UP(sel.bitfield.num_bits, 64) * 8,
                    .type = {.base_type = TYPE_BITFIELD, .dim = {1}},
            },
            };
            md_array_push(stored_selections, ident, alloc);
        }

        for (int64_t i = 0; i < md_array_size(stored_selections); ++i) {
            md_array_push(ir->identifiers, &stored_selections[i], alloc);
        }

        mol_context_t mol_ctx = {
            .atom = { 0, (uint32_t)filter_ctx.mol->atom.count },
            .residue = { 0, (uint32_t)filter_ctx.mol->residue.count },
            .chain = { 0, (uint32_t)filter_ctx.mol->chain.count },
        };

        eval_context_t static_ctx = {
            .ir = ir,
            .mol = filter_ctx.mol,
            .mol_ctx = mol_ctx,
            .temp_alloc = temp_alloc,
        };

        if (static_check_node(node, &static_ctx)) {
            if (node->data.type.base_type == TYPE_BITFIELD && is_scalar(node->data.type)) {
                bit_clear(target.bits, 0, target.num_bits);
                data_t data = {
                    .ptr = &target,
                    .size = DIV_UP(target.num_bits, 64) * 8,
                    .type = node->data.type,
                };
                eval_context_t eval_ctx = {
                    .mol = filter_ctx.mol,
                    .mol_ctx = mol_ctx,
                    .temp_alloc = temp_alloc,
                };

                result = evaluate_node(&data, node, &eval_ctx);
            } else {
                md_print(MD_LOG_TYPE_ERROR, "Filter did not evaluate to a bitfield");
            }
        }
    }

    for (int64_t i = 0; i < md_array_size(ir->errors); ++i) {
        md_printf(MD_LOG_TYPE_ERROR, "%.*s", ir->errors[i].error.len, ir->errors[i].error.ptr);
    }

    md_arena_allocator_destroy(alloc);
    return result;
}