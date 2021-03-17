#include "md_script.h"
#include "md_allocator.h"
#include "md_log.h"
#include "md_molecule.h"
#include "core/common.h"
#include "core/str_util.h"
#include "core/array.inl"
#include "core/file.inl"
#include "core/arena_allocator.h"
#include "core/compiler.h"
#include "core/bitop.h"

#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include <float.h>

#pragma warning(disable:4063) // Single character tokens not being valid tokens

#define MAX_SUPPORTED_PROC_ARGS 8
#define MAX_SUPPORTED_TYPE_DIMS 4

#define DIV_UP(x, y) ((x + (y-1)) / y) // Round up division

// #############################
// ###   TYPE DECLARATIONS   ###
// #############################

typedef struct eval_context_t {
    md_allocator_i* temp_alloc;
    const md_molecule* mol;
} eval_context_t;

typedef struct val_context_t {
    md_script_ir* ir;
    md_allocator_i* temp_alloc;
    const md_molecule* mol;
} val_context_t;

typedef enum token_type_t {
    TOKEN_UNDEF = 0,
    // Reserve the first indices for character literals
    TOKEN_IDENT = 128, // idenfitier
    TOKEN_LE, // '<='
    TOKEN_GE, // '>='
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

typedef enum context_level_t {
    LEVEL_ATOM,
    LEVEL_RESIDUE,
    LEVEL_CHAIN,
} context_level_t;

typedef enum flags_t {
    FLAG_DYNAMIC = 0x1,     // Indicates that it needs to be reevaluated for every frame of the trajectory (it has a dependence on atomic positions)
} flags_t;

typedef struct token_t {
    token_type_t type;
    str_t str;
    uint64_t line;
} token_t;

typedef struct bitfield_t {
    uint64_t* bits;
    uint32_t  num_bits;
} bitfield_t;

typedef struct frange_t {
    float beg;
    float end;
} frange_t;

typedef struct irange_t {
    int32_t beg;
    int32_t end;
} irange_t;

typedef struct type_info_t {
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
} type_info_t;

// There is a conceptual twist here if len should belong to the type_info_t or if it should belong to data_t
// We need to differentiate between different lengths of arrays with the same type in our system, so for now it should probably reside within the type_info
// I.e. the array length should be part of the fundamental type. This only strenghtens the type-system

// This is the struct which is passed as arguments into the procedures
typedef struct data_t {
    type_info_t type;
    void*       ptr;    // Pointer to the data
    uint64_t    size;   // Size in bytes of data (This we want to determine during the static check so we can allocate the data when evaluating)
    str_t       unit;   // The unit of the particular data, e.g rad as in radians, or ° as in degrees, or Å as in Ångström (This is added as a suffix when displaying the data)
} data_t;

// This is data stored directly in the nodes to hold scalar values
typedef union value_t {
    str_t       _string;
    float       _float;
    int32_t     _int;
    frange_t    _frange;
    irange_t    _irange;
    bool        _bool;
    bitfield_t  _bitfield;
} value_t;

typedef struct identifier_t {
    str_t       name;
    data_t      data;
    flags_t     flags;
} identifier_t;

typedef struct procedure_t {
    str_t name;
    type_info_t return_type;
    uint32_t num_args;
    type_info_t arg_type[MAX_SUPPORTED_PROC_ARGS];
    int  (*proc_ptr)(data_t*, data_t[], eval_context_t*); // Function pointer to the "C" procedure
    flags_t flags;
} procedure_t;

typedef struct ast_node_t ast_node_t;
struct ast_node_t {
    // @OPTIMIZE: Make specific types for each type of Ast_node, where the first member is the type which we can use and cast after.
    ast_type_t      type;    
    token_t         token;          // Corresponding token from which the node was created (used to tie errors back into src)
    flags_t         flags;

    // PAYLOAD
    ast_node_t**    children;       // For AND, OR, NOT, procedure arguments etc.
    value_t         value;          // Scalar values for storing data directly inside the node
    data_t          data;           // Structure for passing as argument into procedure. (Holds pointer, length and type)

    procedure_t*    proc;           // Procedure reference
    identifier_t*   ident;          // Identifier reference
};

typedef struct tokenizer_t {
    str_t str;
    uint64_t cur;
    uint64_t line;
} tokenizer_t;

typedef struct expression_t {
    ast_node_t*   node;     
    identifier_t* ident;    // Every expression should end up being assigned to an identifier. This is that identifier
    str_t str;
} expression_t;

typedef struct md_script_ir {
    struct md_allocator_i *arena; // We could use a direct raw interface here to save some pointer indirections
    str_t str;  // Original string containing the 'source'
    
    // These are resizable arrays
    md_script_error *errors;
    expression_t    **expressions;
    identifier_t    **identifiers;
    ast_node_t      **nodes;

    expression_t** eval_targets; // A list of all expressions which need to be evaluated every frame

    bool record_errors; // This is to toggle if new errors should be registered... We don't want to flood the errors
    const char* stage;  // This is just to supply a context for the errors i.e. which stage the error occured
} md_script_ir;

typedef struct parse_context_t {
    md_script_ir* ir;
    ast_node_t* node;   // This represents the current root node in the tree
    tokenizer_t tokenizer;
    md_allocator_i* temp_alloc;
} parse_context_t;


// ##########################
// ###   CORE FUNCTIONS   ###
// ##########################

static inline bool compare_type_info_dim(type_info_t a, type_info_t b) {
    return (memcmp(&a.dim, &b.dim, sizeof(a.dim)) == 0) && (a.len_dim == b.len_dim);
}

static inline bool compare_type_info(type_info_t a, type_info_t b) {
    if (a.base_type == TYPE_BITFIELD && b.base_type == TYPE_BITFIELD) {
        if ((a.level == -1 || b.level == -1) || a.level == b.level) {
            return compare_type_info_dim(a, b);
        }
        return false;
    }
    return memcmp(&a, &b, sizeof(type_info_t)) == 0;
}

static inline bool is_undefined_type(type_info_t ti) {
    return memcmp(&ti, &(type_info_t){0}, sizeof(type_info_t)) == 0;
}

static inline bool is_array(type_info_t ti) {
    return ti.dim[0] > 1;
}

static inline bool is_scalar(type_info_t ti) {
    return ti.dim[0] == 1 && ti.dim[1] == 0 && ti.len_dim == 0;
}

static inline bool is_variable_length(type_info_t ti) {
    return ti.dim[ti.len_dim] == -1;
}

static inline uint64_t type_info_array_len(type_info_t ti) {
    return ti.dim[ti.len_dim];
}

static inline uint64_t element_count(data_t arg) {
    ASSERT(arg.type.dim[arg.type.len_dim] >= 0);
    return arg.type.dim[arg.type.len_dim];
}

// If the type is a scalar, then it is its own element type.
// If the type is an array, then its element type is the type of one element within the array
static inline type_info_t type_info_element_type(type_info_t ti) {
    if (is_scalar(ti)) return ti;
    type_info_t elem_ti = ti;
    elem_ti.dim[elem_ti.len_dim] = 1;
    return elem_ti;
}

// Size of the base type structure
static inline uint64_t base_type_byte_size(base_type_t type) {
    switch (type) {
    case TYPE_UNDEFINED: return 0;
    case TYPE_FLOAT:    return sizeof(float);
    case TYPE_INT:      return sizeof(int);
    case TYPE_BOOL:     return sizeof(bool);        
    case TYPE_FRANGE:   return sizeof(frange_t);
    case TYPE_IRANGE:   return sizeof(irange_t);
    case TYPE_BITFIELD: return sizeof(bitfield_t);
    case TYPE_STRING:   return sizeof(str_t);
    default:            return 0;
    }
}

static inline bool is_type_dim_compatible(type_info_t from, type_info_t to) {

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
        else if (to.dim[i] == -1) return true; // We consider a zero length array to be a valid type as well
    }
    return false;
}

static inline bool is_type_equivalent(type_info_t a, type_info_t b) {
    return memcmp(&a, &b, sizeof(type_info_t)) == 0;
}

static inline bool is_type_directly_compatible(type_info_t from, type_info_t to) {
    if (compare_type_info(from, to)) return true;

    if (is_scalar(from)) {
        if (from.base_type == to.base_type) {
            // To support varying length convention used by procedure args
            if (to.dim[0] == -1) return true;
        }
    }
    else {
        if (from.base_type == to.base_type) {
            if (is_variable_length(to)) {
                return is_type_dim_compatible(from, to);
            }
        }
    }
    return false;
}

static inline bool compare_type_info_array(const type_info_t a[], const type_info_t b[], uint64_t num_arg_types) {
    for (uint64_t i = 0; i < num_arg_types; ++i) {
        if (!is_type_directly_compatible(a[i], b[i])) {
            return false;
        }
    }
    return true;
}

// Returns the count of elements for this type, i.e an array of [4][4][1] would return 16
// NOTE: It does not include the last dimension which encodes the length of the array
// 
static inline uint64_t type_info_element_stride(type_info_t ti) {
    if (ti.len_dim == 0) return 1;
    uint64_t stride = (uint64_t)ti.dim[0];
    for (int32_t i = 1; i < ti.len_dim; ++i) {
        stride *= ti.dim[i];
    }
    return stride;
}

static inline uint64_t type_info_byte_stride(type_info_t ti) {
    return type_info_element_stride(ti) * base_type_byte_size(ti.base_type);
}

static bool allocate_data(data_t* data, const md_molecule* mol, md_allocator_i* alloc) {
    ASSERT(data && data->ptr == NULL && !is_undefined_type(data->type));
    ASSERT(mol);
    ASSERT(alloc);
    
    const uint64_t array_len = type_info_array_len(data->type);
        // The convention should be that this function should have never been called in the first place since the arrays length was 0.
    ASSERT(array_len > 0);

    // Do the array allocation
    const uint64_t array_bytes = type_info_byte_stride(data->type) * array_len;
    data->ptr = md_alloc(alloc, array_bytes);
    data->size = array_bytes;

    // Perform extra allocations if necessary
    if (data->type.base_type == TYPE_BITFIELD) {
        bitfield_t* bf = data->ptr;
        const uint64_t num_bits = mol->atom.count;
        const uint64_t bit_size = DIV_UP(mol->atom.count, 64); // Compute how many uint64_t we need to allocate the bitfield, 64-bits per uint64_t.
        for (uint64_t i = 0; i < array_len; ++i) {
            bf[i].bits = md_alloc(alloc, bit_size);
            bf[i].num_bits = num_bits;
        }
    }

    return true;
}

static inline bool is_number(token_type_t type) {
    return type == TOKEN_INT || type == TOKEN_FLOAT;
}

static inline bool is_operator(ast_type_t type) {
    return (type == AST_ADD || type == AST_SUB || type == AST_MUL || type == AST_DIV ||
        type == AST_AND || type == AST_OR || type == AST_NOT);
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
    case TYPE_FRANGE: return "float range";
    case TYPE_IRANGE: return "int range";
    case TYPE_BOOL: return "boolean";
    case TYPE_STRING: return "string";
    case TYPE_BITFIELD: return "bitfield";
    default: return "type out of range";
    }
}

static uint64_t count_atoms(bitfield_t mask, md_molecule* mol) {
    ASSERT(mol);
    ASSERT(mask.num_bits == mol->atom.count);

    return bit_count(mask.bits, 0, mask.num_bits);
}

static uint32_t count_residues(bitfield_t mask, md_molecule* mol) {
    ASSERT(mol);
    ASSERT(mask.num_bits == mol->atom.count);
    ASSERT(mol->atom.residue_idx);
    ASSERT(mol->residue.atom_range);

    uint64_t count = 0;
    // @PERFORMANCE: Improve this by bitscanning
    for (uint64_t i = 0; i < mask.num_bits; ++i) {
        if (bit_test(mask.bits, i)) {
            ++count;
            i = mol->residue.atom_range[mol->atom.residue_idx[i]].end; // Skip the remaining residue
        }
    }

    return count;
}

static uint32_t count_chains(bitfield_t mask, md_molecule* mol) {
    ASSERT(mol);
    ASSERT(mask.num_bits == mol->atom.count);
    ASSERT(mol->atom.chain_idx);
    ASSERT(mol->chain.atom_range);

    uint64_t count = 0;
    // @PERFORMANCE: Improve this by bitscanning
    for (uint64_t i = 0; i < mask.num_bits; ++i) {
        if (bit_test(mask.bits, i)) {
            ++count;
            i = mol->chain.atom_range[mol->atom.chain_idx[i]].end; // Skip the remaining chain
        }
    }

    return count;
}

// Include all procedures, operators and defines
#include "core/script_functions.inl"




// ############################
// ###   HELPER FUNCTIONS   ###
// ############################

static inline bool is_type_implicitly_convertible(type_info_t from, type_info_t to) {
    for (uint64_t i = 0; i < ARRAY_SIZE(casts); ++i) {
        if (is_type_directly_compatible(casts[i].arg_type[0], from) &&
            is_type_directly_compatible(casts[i].return_type, to)) {
            return true;
        }
    }

    return false;
}

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

static ast_node_t* create_node(md_script_ir* ir, ast_type_t type, token_t token) {
    ast_node_t* node = md_alloc(ir->arena, sizeof(ast_node_t));
    memset(node, 0, sizeof(ast_node_t));
    node->type = type;
    node->token = token;
    md_array_push(ir->nodes, node, ir->arena);
    return node;
}

static identifier_t* get_identifier(md_script_ir* ir, str_t name) {
    const uint64_t num_identifiers = md_array_size(ir->identifiers);
    for (uint64_t i = 0; i < num_identifiers; ++i) {
        if (compare_str(ir->identifiers[i]->name, name)) {
            return ir->identifiers[i];
        }
    }
    return NULL;
}

static identifier_t* create_identifier(md_script_ir* ir, str_t name) {
    ASSERT(get_identifier(ir, name) == NULL);

    identifier_t* ident = md_alloc(ir->arena, sizeof(identifier_t));
    memset(ident, 0, sizeof(identifier_t));
    ident->name = copy_str(name, ir->arena);
    md_array_push(ir->identifiers, ident, ir->arena);
    return ident;
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

        if (child->type > parent->type) { // type is sorted on precedence
            child->children[0] = parent;
            parent->children[i] = grand_child;
            *node = child;
        }
    }
}

static bool convert_node(ast_node_t* node, type_info_t new_type, md_script_ir* ir) {
    // We need to support the following conversions
    // Direct value conversions of Constant Value nodes, which are always scalar
    // Cast conversions which are performed on results of expressions / operators / proc_calls where we cannot directly convert the nodes

    const type_info_t from = node->data.type;
    const type_info_t to   = new_type;

    if (is_scalar(from)) {
        if (is_scalar(to)) {
            if (node->type == AST_CONSTANT_VALUE) {
                // If its a constant value, we modify the node directly
                ASSERT(is_scalar(node->data.type)); // Constant value nodes should always be scalars
                switch (to.base_type) {
                case TYPE_FLOAT:
                    ASSERT(from.base_type == TYPE_INT);
                    node->value._float = (float)node->value._int;
                    node->data.type.base_type = TYPE_FLOAT;
                    return true;
                case TYPE_FRANGE:
                    ASSERT(from.base_type == TYPE_IRANGE);
                    node->value._frange = (frange_t){(float)node->value._irange.beg, (float)node->value._irange.end};
                    node->data.type.base_type = TYPE_FRANGE;
                    return true;
                case TYPE_IRANGE:
                    ASSERT(from.base_type == TYPE_INT);
                    node->value._irange = (irange_t){node->value._int, node->value._int};
                    node->data.type.base_type = TYPE_IRANGE;
                    return true;
                default:
                    ASSERT(false);
                }
            } else {
                // We need to convert this node into a cast node and add the original node data as a child
                ast_node_t* node_copy = create_node(ir, node->type, node->token);
                memcpy(node_copy, node, sizeof(ast_node_t));
                node->type = AST_CAST;
                node->data.type = to;
                node->children = 0; // node_copy have taken over the children, we need to zero this to trigger a proper allocation in next step
                md_array_push(node->children, node_copy, ir->arena);
                return true;
            }
        }
    } else { // is_array(from)
        if (is_array(to)) {
            ASSERT(false);
        }
    }

    // Unsupported conversion
    ASSERT(false); // Should not reach this since the implicit conversion check should have not allowed this
    return false;
}

static void create_error(md_script_ir* ir, token_t token, const char* format, ...) {
    if (!ir->record_errors) return;
    ir->record_errors = false;

    char buffer[512] = {0};
    va_list args;
    va_start(args, format);
    int len = vsnprintf(buffer, ARRAY_SIZE(buffer), format, args);
    va_end(args);

    ASSERT(len > 0);
    ASSERT(len < ARRAY_SIZE(buffer));

    char* err_str = md_alloc(ir->arena, len + 1);
    memcpy(err_str, buffer, len);
    err_str[len] = '\0';

    md_script_error error = {0};
    error.line = (uint32_t)token.line;
    error.offset = (uint32_t)(token.str.ptr - ir->str.ptr);
    error.length = (uint32_t)token.str.len;
    error.str_ptr = err_str;
    error.str_len = (uint64_t)len;
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
        const char* long_ass_carret = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
        md_printf(MD_LOG_TYPE_DEBUG, "%.*s", (end-beg), beg);
        md_printf(MD_LOG_TYPE_DEBUG, "%*s^%.*s", (token.str.ptr - beg), "", token.str.len-1, long_ass_carret);
    }
}

static procedure_t* find_procedure_supporting_arg_types_in_candidates(str_t name, const type_info_t arg_types[], uint64_t num_arg_types, procedure_t* candidates, uint64_t num_cantidates, bool allow_implicit_conversions) {
    for (uint64_t i = 0; i < num_cantidates; ++i) {
        procedure_t* proc = &candidates[i];
        if (compare_str(proc->name, name)) {
            if (num_arg_types == proc->num_args) {
                if (compare_type_info_array(arg_types, proc->arg_type, num_arg_types)) {
                    return proc;
                }
            }
        }
    }

    if (allow_implicit_conversions) {
        // Try using implicit casts on types where this is allowed.
        // int to float and irange to frange.
        // Store the best candidate with the lowest cost heuristic (number of conversions needed)
        procedure_t* best_proc = NULL;
        uint32_t     best_cost = 0xFFFFFFFFU;
        for (uint64_t i = 0; i < num_cantidates; ++i) {
            procedure_t* proc = &candidates[i];
            if (compare_str(proc->name, name)) {
                if (num_arg_types == proc->num_args) {
                    // Same name and same number of args...

                    uint32_t cost = 0;
                    for (uint64_t j = 0; j < proc->num_args; ++j) {
                        if (is_type_directly_compatible(arg_types[j], proc->arg_type[j])) {
                            // No conversion needed for this argument (0 cost)
                        }
                        else if (is_type_implicitly_convertible(arg_types[j], proc->arg_type[j])) {
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
                    }
                }
            }
        }
        return best_proc;
    }

    return NULL;
}

static procedure_t* find_procedure_supporting_arg_types(str_t name, const type_info_t arg_types[], uint64_t num_arg_types, bool allow_implicit_conversions) {
    return find_procedure_supporting_arg_types_in_candidates(name, arg_types, num_arg_types, procedures, ARRAY_SIZE(procedures), allow_implicit_conversions);
}

static procedure_t* find_operator_supporting_arg_types(token_type_t op, const type_info_t arg_types[], uint64_t num_arg_types, bool allow_implicit_conversions) {
    // Map token type to string which we use to identify operator procedures
    str_t name = {0};
    switch(op) {
    case '+': name.ptr = "+"; name.len = 1; break;
    case '-': name.ptr = "-"; name.len = 1; break;
    case '*': name.ptr = "*"; name.len = 1; break;
    case '/': name.ptr = "/"; name.len = 1; break;
    case TOKEN_AND: name.ptr = "and"; name.len = 3; break;
    case TOKEN_OR:  name.ptr = "or";  name.len = 2; break;
    case TOKEN_NOT: name.ptr = "not"; name.len = 3; break;
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

static bool expect_token_type(md_script_ir* ir, token_t token, token_type_t type) {
    if (token.type != type) {
        create_error(ir, token, "Unexpected token '%.*s', expected token '%s'", token.str.len, token.str.ptr, get_token_type_str(type));
        return false;
    }
    return true;
}



/*
static void allocate_node_data(ast_node_t* node, md_allocator_i* alloc) {
ASSERT(node);
ASSERT(alloc);
ASSERT(node->data.type.base_type != TYPE_UNDEFINED);

if (is_scalar(node->data.type)) {
// We use the value field within the node.
} else {
uint64_t element_payload = 0; // If each element has an additional payload, e.g. bitfields (need to allocate all the bits)
uint64_t element_size = type_info_stride(node->data.type);
}

if (node->type == AST_CONSTANT_VALUE) {
// we don't need to allocate any data for simple types (float, int, bool, frange, irange, string*, bitrange)
// *string has its data allocated when parsed.

switch (node->data.type.base_type) {
case TYPE_BITFIELD:
{
ASSERT(node->value._bitfield.bits = 0);
node->value._bitfield.bits = md_alloc(alloc, DIV_UP(node->value._bitfield.num_bits, 64));
break;
}
default:
break;
}
node->data.ptr = &node->value;
} else {
uint64_t required_size = base_type_byte_size();
}
}
*/



// #####################
// ###   TOKENIZER   ###
// #####################

static token_t tokenizer_get_next_from_buffer(tokenizer_t* tokenizer) {
    token_t token = {0};

    if (tokenizer->cur >= tokenizer->str.len) {
        tokenizer->cur = tokenizer->str.len;
        token.type = TOKEN_END;
        return token;
    }
    
    uint64_t line = tokenizer->line;
    const char* buf = tokenizer->str.ptr;
    const uint64_t len = tokenizer->str.len;
    for (uint64_t i = tokenizer->cur; i != len; ++i) {
        uint64_t j = 0;
        
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
            const uint64_t n = j - i;
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
            const uint64_t n = (i + 1 < len && is_symbol(buf[i + 1])) ? 2 : 1;
            
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
            break;
        }

        if (!token.type && i == len - 1) {
            // We are going to hit the end in next iteration and have nothing to show for.
            tokenizer->cur = tokenizer->str.len;
            token.type = TOKEN_END;
            token.line = line;
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

static int print_type_info(char* buf, int buf_size, type_info_t info) {
    int result = snprintf(buf, buf_size, "%s", get_value_type_str(info.base_type));
    int i = 0;
    while (info.dim[i] > 0 && i < MAX_SUPPORTED_TYPE_DIMS) {
        result += snprintf(buf + result, buf_size, "[%i]", (int)info.dim[i]);
        ++i;
    }
    return result;
}

#if MD_DEBUG

static void print_bitfield(FILE* file, bitfield_t bitfield) {
    const uint64_t max_bits = 64;
    const uint64_t num_bits = MIN(bitfield.num_bits, max_bits); // Don't want to show more than the first 128 bits ...
    fprintf(file, "[");
    for (uint64_t i = 0; i < num_bits; ++i) {
        fprintf(file, "%i", bit_test(bitfield.bits, i) ? 1 : 0);
    }
    if (num_bits == max_bits) {
        fprintf(file, "...");
    }
    fprintf(file, "]");
}

static void print_array(FILE* file, data_t data) {
    fprintf(file, "Array");
}

static void print_value(FILE* file, data_t data) {
    if (is_array(data.type)) {
        print_array(file, data);
    } else {
        if (data.ptr) {
            switch(data.type.base_type) {
            case TYPE_BITFIELD:
                print_bitfield(file, *(bitfield_t*)data.ptr);
                break;
            case TYPE_BOOL:
                fprintf(file, "%s", *(bool*)data.ptr ? "true" : "false");
                break;
            case TYPE_INT:
                fprintf(file, "%i", *(int*)data.ptr);
                break;
            case TYPE_FLOAT:
                fprintf(file, "%f", *(float*)data.ptr);
                break;
            case TYPE_IRANGE:
            {
                irange_t rng = *(irange_t*)data.ptr;
                fprintf(file, "%i:%i", rng.beg, rng.end);
                break;
            }
            case TYPE_FRANGE:
            {
                frange_t rng = *(frange_t*)data.ptr;
                fprintf(file, "%f:%f", rng.beg, rng.end);
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

static void save_expressions_to_json(const expression_t** expr, uint64_t num_expr, str_t filename) {
    FILE* file = (FILE*)md_file_open(filename.ptr, filename.len, "w");

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
        md_file_close(file);
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
    while (next = tokenizer_peek_next(&ctx->tokenizer), next.type != TOKEN_END) {
        ast_node_t* arg = parse_expression(ctx);
        next = tokenizer_peek_next(&ctx->tokenizer);
        if (arg && (next.type == ',' || next.type == token_type)) {
            md_array_push(args, arg, ctx->temp_alloc);
            if (next.type == ',')
                tokenizer_consume_next(&ctx->tokenizer);
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
    token_t next = tokenizer_peek_next(&ctx->tokenizer);
    ast_node_t* node = 0;
    
    if (next.type == '(') {
        tokenizer_consume_next(&ctx->tokenizer); // '('

        ast_node_t **args = parse_comma_separated_arguments_until_token(ctx, ')');
        if (expect_token_type(ctx->ir, tokenizer_consume_next(&ctx->tokenizer), ')')) {
            const uint64_t num_args = md_array_size(args);
            if (num_args) {
                node = create_node(ctx->ir, AST_PROC_CALL, token);
                md_array_push_array(node->children, args, num_args, ctx->ir->arena);
            }
        } else {
            create_error(ctx->ir, token, "Unexpected end of argument list");
        }
    }

    return node;
}

ast_node_t* parse_identifier(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(&ctx->tokenizer);
    ASSERT(token.type == TOKEN_IDENT);
    
    const str_t ident = token.str;
    ast_node_t* node = 0;
    identifier_t* c = 0;
    
    if (is_identifier_procedure(ident)) {
        // The identifier has matched some procedure name
        node = parse_procedure_call(ctx, token);
    } else if (c = find_constant(ident)) {
        node = create_node(ctx->ir, AST_IDENTIFIER, token);
        node->ident = c;
        node->data = c->data;
    } else {
        // Identifier
        token_t next = tokenizer_peek_next(&ctx->tokenizer);
        if (next.type == '(') {
            // This is an intended procedure call, but has no matching procedure by that name
            create_error(ctx->ir, token, "Undefined procedure '%.*s'", ident.len, ident.ptr);

            token_type_t token_types[] = {'(', ')', ';'};
            int paren_bal = 0;
            while (next = tokenizer_consume_until_type(&ctx->tokenizer, token_types, ARRAY_SIZE(token_types)), next.type != TOKEN_END) {
                switch (next.type) {
                case '(':
                    tokenizer_consume_next(&ctx->tokenizer);
                    ++paren_bal;
                    break;
                case ')':
                    tokenizer_consume_next(&ctx->tokenizer);
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
    token_t token = tokenizer_consume_next(&ctx->tokenizer);
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
    token_t token = tokenizer_consume_next(&ctx->tokenizer);
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
    
    ast_node_t* node = create_node(ctx->ir, AST_ASSIGNMENT, token);
    md_array_push(node->children, lhs, ctx->ir->arena);
    md_array_push(node->children, rhs, ctx->ir->arena);
    node->data.type = rhs->data.type;
    return node;
}

ast_node_t* parse_comparison(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(&ctx->tokenizer);
    ASSERT(token.type == '<' || token.type == '>' || token.type == TOKEN_EQ || token.type == TOKEN_LE || token.type == TOKEN_GE);
    
    ast_node_t* lhs = ctx->node;
    ast_node_t* rhs = parse_expression(ctx);
    ast_node_t* node = 0;

    if (!lhs) {
        create_error(ctx->ir, token,
            "Left hand side of operator '%s' did not evaluate to a proper expression.", get_token_type_str(token.type));
        return NULL;
    }

    if (!rhs) {
        create_error(ctx->ir, token,
            "Right hand side of operator '%s' did not evaluate to a proper expression.", get_token_type_str(token.type));
        return NULL;
    }

    const type_info_t arg_type[2] = {lhs->data.type, rhs->data.type};
    procedure_t* op = find_operator_supporting_arg_types(token.type, arg_type, 2, false);

    if (op) {
        node = create_node(ctx->ir, AST_PROC_CALL, token);
        md_array_push(node->children, lhs, ctx->ir->arena);
        md_array_push(node->children, rhs, ctx->ir->arena);
        node->data.type = op->return_type;
    } else {
        char arg_type_str[2][64] = {0};
        print_type_info(arg_type_str[0], ARRAY_SIZE(arg_type_str[0]), lhs->data.type);
        print_type_info(arg_type_str[1], ARRAY_SIZE(arg_type_str[1]), rhs->data.type);
        create_error(ctx->ir, token,
            "Could not find supporting operator '%s' with left hand side type (%s) and right hand side type (%s)",
            get_token_type_str(token.type), arg_type_str[0], arg_type_str[1]);
    }

    return node;
}

ast_node_t* parse_arithmetic(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(&ctx->tokenizer);
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
    token_t token = tokenizer_consume_next(&ctx->tokenizer);
    ASSERT(token.type == TOKEN_FLOAT || token.type == TOKEN_INT || token.type == TOKEN_STRING || token.type == ':');
    token_t next = tokenizer_peek_next(&ctx->tokenizer);
    
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
            expect_token_type(ctx->ir, tokenizer_consume_next(&ctx->tokenizer), ':');
        } else { // == ':'
            beg = -FLT_MAX;
        }
        
        next = tokenizer_peek_next(&ctx->tokenizer);
        if (is_number(next.type)) {
            token = tokenizer_consume_next(&ctx->tokenizer);
            end = parse_float(token.str);
            is_float |= (token.type == TOKEN_FLOAT);
        } else {
            end = FLT_MAX;
        }

        if (is_float) {
            node->data.type = (type_info_t)TI_FRANGE;
            node->value._frange.beg = (float)beg;
            node->value._frange.end = (float)end;
            node->data.ptr = &node->value;
            node->data.size = base_type_byte_size(node->data.type.base_type);
        } else {
            node->data.type = (type_info_t)TI_IRANGE;
            node->value._irange.beg = (beg == -FLT_MAX) ?  INT32_MIN : (int32_t)beg;
            node->value._irange.end = (end ==  FLT_MAX) ?  INT32_MAX : (int32_t)end;
            node->data.ptr = &node->value;
            node->data.size = base_type_byte_size(node->data.type.base_type);
        }
    }
    else if (token.type == TOKEN_FLOAT) {
        node->data.type = (type_info_t)TI_FLOAT;
        node->value._float = (float)parse_float(token.str);
        node->data.ptr = &node->value;
        node->data.size = base_type_byte_size(node->data.type.base_type);
    }
    else if (token.type == TOKEN_INT) {
        node->data.type = (type_info_t)TI_INT;
        node->value._int = (int32_t)parse_int(token.str);
        node->data.ptr = &node->value;
        node->data.size = base_type_byte_size(node->data.type.base_type);
    }
    else if (token.type == TOKEN_STRING) {
        node->data.type = (type_info_t)TI_STRING;
        node->value._string = copy_str(substr(token.str, 1, MAX(0, ((int)token.str.len - 2))), ctx->ir->arena); // remove quotation marks
        node->data.ptr = &node->value;
        node->data.size = base_type_byte_size(node->data.type.base_type);
    }
    else {
        ASSERT(false);
    }
    return node;
}

ast_node_t* parse_array_subscript(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(&ctx->tokenizer);
    ASSERT(token.type == '[');

    if (ctx->node) {
        ast_node_t* node = ctx->node;
        // We need to consume the node and turn it into a subscript
        // Put the original node as the first child, and insert whatever arguments we have as the following children.

        ast_node_t** elements = parse_comma_separated_arguments_until_token(ctx, ']');
        if (expect_token_type(ctx->ir, tokenizer_consume_next(&ctx->tokenizer), ']')) {
            const uint64_t num_elements = md_array_size(elements);
            if (num_elements) {

                ast_node_t* node_copy = create_node(ctx->ir, node->type, token);
                memcpy(node_copy, node, sizeof(ast_node_t));

                memset(node, 0, sizeof(ast_node_t));
                node->type = AST_ARRAY_SUBSCRIPT;
                node->token = token;

                md_array_push(node->children, node_copy, ctx->ir->arena);
                md_array_push_array(node->children, elements, num_elements, ctx->ir->arena);

                return ctx->node;
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
    token_t token = tokenizer_consume_next(&ctx->tokenizer);
    ASSERT(token.type == '{');

    ast_node_t* node = 0;

    ast_node_t** elements = parse_comma_separated_arguments_until_token(ctx, '}');
    if (expect_token_type(ctx->ir, tokenizer_consume_next(&ctx->tokenizer), '}')) {
        const uint64_t num_elements = md_array_size(elements);
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
    token_t token = tokenizer_consume_next(&ctx->tokenizer);
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
    token_t token = tokenizer_consume_next(&ctx->tokenizer);
    ASSERT(token.type == '(');
    ast_node_t* expr = parse_expression(ctx);
    ast_node_t* node = create_node(ctx->ir, AST_EXPRESSION, token);
    md_array_push(node->children, expr, ctx->ir->arena);
    token = tokenizer_consume_next(&ctx->tokenizer);
    expect_token_type(ctx->ir, token, ')');
    return node;
}

static ast_node_t* parse_argument(parse_context_t* ctx) {
    token_t next = tokenizer_peek_next(&ctx->tokenizer);
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
        tokenizer_consume_next(&ctx->tokenizer);
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
    while (token = tokenizer_peek_next(&ctx->tokenizer), token.type != TOKEN_END) {
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
            tokenizer_consume_next(&ctx->tokenizer);
            return NULL;
            default:
            create_error(ctx->ir, token, "Unexpected token value! '%.*s'", token.str.len, token.str.ptr);
            return NULL;
        }
    }
    
    done:
    return ctx->node;
}

// ####################
// ###   EVALUATE   ###
// ####################

static int evaluate_node(data_t*, const ast_node_t*, eval_context_t*);

static int evaluate_proc_call(data_t* dst, const ast_node_t* node, eval_context_t* eval_ctx) {
    ASSERT(node && node->type == AST_PROC_CALL);
    ASSERT(node->proc);

    // This procedure can be called within the static check to 'query' the result of the procedure.
    // In such case, dst will be NULL.

    const uint64_t num_args = md_array_size(node->children);
    data_t args[MAX_SUPPORTED_PROC_ARGS] = {0};

    for (uint64_t i = 0; i < num_args; ++i) {
        // We need to evaluate the argument nodes first before we make the proc call.
        // In this context we are not interested in storing any data (since it is only used to be passed as arguments)
        // so we can allocate the data required for the node with the temp alloc
        args[i].type = node->children[i]->data.type;
        allocate_data(&args[i], eval_ctx->mol, eval_ctx->temp_alloc);
        evaluate_node(&args[i], node->children[i], eval_ctx);
    }

    int result = node->proc->proc_ptr(dst, args, eval_ctx);
    return result;
}

static int evaluate_constant_value(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_CONSTANT_VALUE);
    ASSERT(node->data.ptr && node->data.size > 0);  // Assume that the static check set the data already
    ASSERT(dst && dst->ptr && dst->size == node->data.size); // Make sure we have room for the data
    (void)ctx;

    copy_data(dst, &node->data);

    return 0;
}

static int evaluate_identifier(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_IDENTIFIER);
    ASSERT(node->data.ptr && node->data.size > 0); // Assume that this has been evaluated already!
    ASSERT(dst && dst->ptr && dst->size == node->data.size); // Make sure we have room for the data
    (void)ctx;
    
    // We assume in this context that we have a reference to the identifier, we copy the contents
    copy_data(dst, &node->data);

    return 0;
}

static int evaluate_assignment(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ASSIGNMENT);
    ASSERT(node->children);

    const uint64_t num_children = md_array_size(node->children);
    ASSERT(num_children == 2);
    ast_node_t* lhs = node->children[0];
    ast_node_t* rhs = node->children[1];

    ASSERT(lhs && lhs->type == AST_IDENTIFIER && lhs->ident);
    ASSERT(rhs);

    // We know from the static check that the lhs (child[0]) is an identifier and the rhs is a valid expression with a known data size
    // the data size could be zero if it is an array of length 0, in that case we don't need to allocate anything.
    // This should also be the conceptual entry point for all evaluations in this declarative language

    return -1;
}

static int evaluate_array(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY);
    ASSERT(dst && dst->ptr && dst->size == node->data.size);
    ASSERT(compare_type_info(dst->type, node->data.type));

    const uint64_t elem_size = type_info_byte_stride(dst->type);
    const type_info_t elem_type = type_info_element_type(dst->type);

    for (uint64_t i = 0; i < md_array_size(node->children); ++i) {
        data_t data = {elem_type, (char*)dst->ptr + elem_size * i, elem_size};
        evaluate_node(&data, node->children[i], ctx);
    }

    return 0;
}

static int evaluate_array_subscript(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY_SUBSCRIPT);
    ASSERT(dst && dst->ptr && dst->size == node->data.size);
    ASSERT(compare_type_info(dst->type, node->data.type));
    ASSERT(md_array_size(node->children) == 2);

    const ast_node_t* arr = node->children[0];
    const ast_node_t* idx = node->children[1];

    data_t arr_data = {arr->data.type};
    data_t idx_data = {idx->data.type};

    if (!allocate_data(&arr_data, ctx->mol, ctx->temp_alloc)) return -1;
    if (!allocate_data(&idx_data, ctx->mol, ctx->temp_alloc)) return -1;

    int result = 0;
    if ((result = evaluate_node(&arr_data, arr, ctx)) != 0) return result;
    if ((result = evaluate_node(&idx_data, idx, ctx)) != 0) return result;

    ASSERT(arr_data.ptr);
    ASSERT(idx_data.ptr);

    const uint64_t elem_size = type_info_byte_stride(arr_data.type);
    int offset = 0;
    int length = 1;
    if (compare_type_info(idx_data.type, (type_info_t)TI_INT)) {
        offset = *((int*)(idx_data.ptr));
    } else if (compare_type_info(idx_data.type, (type_info_t)TI_IRANGE)) {
        irange_t rng = *((irange_t*)(idx_data.ptr));
        offset = rng.beg;
        length = rng.end - rng.beg + 1;
    } else {
        ASSERT(false);
    }

    const data_t src = {dst->type, (char*)arr_data.ptr + elem_size * offset, elem_size * length};
    copy_data(dst, &src);

    return 0;
}

static int evaluate_cast(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_CAST);
    ASSERT(dst && dst->ptr && dst->size == node->data.size);
    ASSERT(compare_type_info(dst->type, node->data.type));
    ASSERT(md_array_size(node->children) == 1);

    

    return -1;
}

static int evaluate_context(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_CONTEXT);

    return -1;
}

static int evaluate_node(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);

    switch (node->type) {
        case AST_PROC_CALL:
        case AST_NOT:
        case AST_UNARY_NEG:
        case AST_MUL:
        case AST_DIV:
        case AST_ADD:
        case AST_SUB:
        case AST_AND:
        case AST_OR:
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
        case AST_CAST:
            return evaluate_cast(dst, node, ctx);         // Should probably be baked into function calls
        case AST_CONTEXT:
            return evaluate_context(dst, node, ctx);
        case AST_UNDEFINED:
        default:
            ASSERT(false);
    }
    return false;
}

// ########################
// ###   STATIC CHECK   ###
// ########################

typedef struct static_check_ctx_t {
    md_allocator_i* temp_alloc;
    md_script_ir* ir;
    md_molecule* mol;
} static_check_ctx_t;

// This function performs an actual procedure call with a NULL ptr as the destination value, thus querying the length of the result without anything being written.
// This result is then used to deduce the full array type of our node.
static type_info_t static_check_query_proc_result_type(ast_node_t* node, static_check_ctx_t* ctx) {
    ASSERT(node && node->proc);
    eval_context_t eval_ctx;
    eval_ctx.mol = ctx->mol;
    eval_ctx.temp_alloc = ctx->temp_alloc;

    type_info_t result_type = {0};

    int proc_result = evaluate_proc_call(NULL, node, &eval_ctx);
    if (proc_result >= 0) { // Zero length is valid
        type_info_t proc_type = node->proc->return_type;
        ASSERT(is_variable_length(proc_type));
        // Find which dimension we should fill in
        uint64_t i = 0;
        for (i = 0; i < MAX_SUPPORTED_TYPE_DIMS; ++i) {
            if (proc_type.dim[i] == -1) break;
        }
        if (i < MAX_SUPPORTED_TYPE_DIMS) {
            result_type = proc_type; 
            result_type.dim[i] = proc_result;
        } else {
            create_error(ctx->ir, node->token, "Unexpected array type. Is not varying in length?");
        }
    } else {
        create_error(ctx->ir, node->token, "Unexpected return value when calling procedure.");
    }

    return result_type;
}

static bool static_check(ast_node_t*, static_check_ctx_t*);

static bool static_check_children(ast_node_t* node, static_check_ctx_t* ctx) {
    bool result = true;
    if (node->children) {
        const uint64_t size = md_array_size(node->children);
        for (uint64_t i = 0; i < size; ++i) {
            ASSERT(node != node->children[i]);
            if (static_check(node->children[i], ctx)) {
                node->flags |= node->children[i]->flags;
            } else {
                result = false;
                break;
            }
        }
    }
    return result;
}

static bool static_check_operator(ast_node_t* node, static_check_ctx_t* ctx) {
    ASSERT(node);
    ASSERT(is_operator(node->type));
    ASSERT(node->children);

    if (static_check_children(node, ctx)) {
        const uint64_t num_args = md_array_size(node->children);
        ast_node_t* arg[2] = {0};
        if (node->type == AST_NOT) {
            ASSERT(num_args == 1);
            arg[0] = node->children[0];
        }
        else {
            ASSERT(num_args == 2);
            arg[0] = node->children[0];
            arg[1] = node->children[1];
        }

        const type_info_t arg_type[2] = {arg[0]->data.type, num_args > 1 ? arg[1]->data.type : (type_info_t){0}};
        procedure_t* op = find_operator_supporting_arg_types(node->token.type, arg_type, num_args, true);

        if (op) {
            // Convert node into a proc call
            node->type = AST_PROC_CALL;
            node->proc = op;

            // Make sure we cast the arguments into the expected types of the procedure.
            for (uint64_t i = 0; i < num_args; ++i) {
                if (!is_type_directly_compatible(arg_type[i], op->arg_type[i])) {
                    // Types are not identical, but should be implicitly convertible
                    ASSERT(is_type_implicitly_convertible(arg_type[i], op->arg_type[i]));
                    convert_node(arg[i], op->arg_type[i], ctx->ir);
                }
            }

            if (is_variable_length(op->return_type)) {
                // Any variable length Operator result means that all array arguments must match in length or be scalars
                // The length of the result is strictly determined by the length of the array input
                if (num_args == 1) {
                    ASSERT(is_array(arg_type[0]) && !is_variable_length(arg_type[0]));
                    node->data.type = arg_type[0];
                } else {
                    if (is_array(arg_type[0]) && !is_array(arg_type[1])) {
                        node->data.type = arg_type[0];
                    } else if (!is_array(arg_type[0]) && is_array(arg_type[1])) {
                        node->data.type = arg_type[1];
                    } else {
                        ASSERT(is_array(arg_type[0]) && is_array(arg_type[1]));
                        if (compare_type_info(arg_type[0], arg_type[1])) {
                            node->data.type = arg_type[0];
                        }
                        else {
                            create_error(ctx->ir, node->token, "The array dimensions of the operands did not match.");
                        }
                    }
                }
            } else {
                node->data.type = op->return_type;
            }
            return true;
        } else {
            char arg_type_str[2][64] = {0};
            for (int i = 0; i < num_args; ++i) {
                if (arg[i])
                    print_type_info(arg_type_str[i], ARRAY_SIZE(arg_type_str[i]), arg[i]->data.type);
                else
                    snprintf(arg_type_str[i], ARRAY_SIZE(arg_type_str[i]), "empty");
            }

            create_error(ctx->ir, node->token,
                "Could not find supporting operator '%s' with left hand side type (%s) and right hand side type (%s)",
                get_token_type_str(node->token.type), arg_type_str[0], arg_type_str[1]);
        }
    }

    return false;
}

static bool static_check_proc_call(ast_node_t* node, static_check_ctx_t* ctx) {
    ASSERT(node->type == AST_PROC_CALL);

    if (static_check_children(node, ctx)) {
        const uint64_t num_args = md_array_size(node->children);
        const str_t proc_name = node->token.str;

        if (num_args == 0) {
            // Zero arguments
            procedure_t* proc = find_procedure_supporting_arg_types(proc_name, 0, 0, false);
            if (proc) {
                node->proc = proc;
                if (is_variable_length(proc->return_type)) {
                    node->data.type = static_check_query_proc_result_type(node, ctx);
                } else {
                    node->data.type = proc->return_type;
                }
                return true;
            }
            else {
                create_error(ctx->ir, node->token,
                    "Could not find matching procedure '.*s' which takes no arguments", (int)proc_name.len, proc_name.ptr);
            }
        }
        else {
            ASSERT(num_args < MAX_SUPPORTED_PROC_ARGS);
            ast_node_t** arg = node->children;
            // One or more arguments
            type_info_t arg_type[MAX_SUPPORTED_PROC_ARGS] = {0};

            for (uint64_t i = 0; i < num_args; ++i) {
                arg_type[i] = arg[i]->data.type;
            }

            procedure_t* proc = find_procedure_supporting_arg_types(proc_name, arg_type, num_args, true);
            if (proc) {
                node->proc = proc;
                // Make sure this flag is set in the tree so it can be propagated upwards later
                if (proc->flags & FLAG_DYNAMIC) node->flags |= FLAG_DYNAMIC;

                // Make sure we cast the arguments into the expected types of the procedure.
                for (uint64_t i = 0; i < num_args; ++i) {
                    if (!is_type_directly_compatible(arg_type[i], proc->arg_type[i])) {
                        // Types are not directly compatible, but should be implicitly convertible
                        ASSERT(is_type_implicitly_convertible(arg_type[i], proc->arg_type[i]));
                        convert_node(arg[i], proc->arg_type[i], ctx->ir);
                    }
                }

                if (is_variable_length(proc->return_type)) {
                    node->data.type = static_check_query_proc_result_type(node, ctx);
                } else {
                    node->data.type = proc->return_type;
                }

                return true;
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
    return false;
}

static bool static_check_constant_value(ast_node_t* node, static_check_ctx_t* ctx) {
    (void)ctx;
    ASSERT(node && node->type == AST_CONSTANT_VALUE);
    ASSERT(node->data.type.base_type != TYPE_UNDEFINED);
    ASSERT(is_scalar(node->data.type));
    
    return true;
}

static bool static_check_array(ast_node_t* node, static_check_ctx_t* ctx) {
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

        type_info_t array_type = {0};
        
        // Pass 1: find a suitable base_type for the array
        for (uint64_t i = 0; i < num_elem; ++i) {
            if (is_scalar(elem[i]->data.type)) {
                if (array_type.base_type == TYPE_UNDEFINED) {
                    // Assign type from first element
                    array_type = elem[i]->data.type;
                }
                else {
                    if (!is_type_directly_compatible(elem[i]->data.type, array_type)) {
                        if (is_type_implicitly_convertible(elem[i]->data.type, array_type)) {
                            // We can convert the element to the array type.   
                        }
                        else if (is_type_implicitly_convertible(array_type, elem[i]->data.type)) {
                            // Option 2: We can convert the array type into the elements type (e.g int to float)
                            // Promote the type of the array
                            array_type = elem[i]->data.type;

                            // Just to make sure, we recheck all elements up until this one
                            for (uint64_t j = 0; j < i; ++j) {
                                if (!is_type_directly_compatible(elem[i]->data.type, array_type) &&
                                    !is_type_implicitly_convertible(elem[i]->data.type, array_type)) {
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
            else {
                create_error(ctx->ir, elem[i]->token, "Only scalar types are allowed within array construct");
                return false;
            }
        }

        // Pass 2: Perform implicit conversions of nodes if required
        for (uint64_t i = 0; i < num_elem; ++i) {
            if (!is_type_directly_compatible(elem[i]->data.type, array_type) &&
                is_type_implicitly_convertible(elem[i]->data.type, array_type)) {
                convert_node(elem[i], array_type, ctx->ir);
            }
        }

        // Finalize the type for the array
        array_type.dim[array_type.len_dim] = (int32_t)num_elem;
        node->data.type = array_type;

        return true;
    }
    return false;
}

static bool static_check_array_subscript(ast_node_t* node, static_check_ctx_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY_SUBSCRIPT && node->children);

    // In array subscripts, only integers and ranges should be allowed.
    // For now, we only support SINGLE integers or ranges
    // @TODO: In the future, se should expand this to support mixed integers and ranges and create a proper subset of the original data

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

static bool static_check_identifier_reference(ast_node_t* node, static_check_ctx_t* ctx) {
    ASSERT(node && node->type == AST_IDENTIFIER);
    
    if (!node->ident) {
        node->ident = get_identifier(ctx->ir, node->token.str);
        node->data = node->ident->data;
        node->flags |= node->ident->flags;    // This means it has a dependency on an identifier, and that identifier might be dynamic... we need to propagate flags
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

static bool static_check_assignment(ast_node_t* node, static_check_ctx_t* ctx) {
    ASSERT(node && node->type == AST_ASSIGNMENT && node->children && md_array_size(node->children) == 2);
    ast_node_t* lhs = node->children[0];
    ast_node_t* rhs = node->children[1];

    ASSERT(lhs && lhs->ident);

    // We don't statically check lhs here since it is an assignment and the identifiers type is deduced from rhs.

    if (static_check(rhs, ctx)) {
        ASSERT(lhs->data.type.base_type == TYPE_UNDEFINED);  // Identifiers type should always be undefined until explicitly assigned.
        if (rhs->data.type.base_type != TYPE_UNDEFINED) {
            node->data          = rhs->data;
            lhs->data           = rhs->data;
            lhs->ident->data    = rhs->data;
            node->flags         = rhs->flags;
            lhs->flags          = rhs->flags;
            lhs->ident->flags   = rhs->flags;
            return true;
        } else {
            create_error(ctx->ir, node->token, "Right hand side of assignment has an unresolved type");
        }
    }

    return false;
}

static bool static_check_context(ast_node_t* node, static_check_ctx_t* ctx) {
    ASSERT(node && node->type == AST_CONTEXT && node->children && md_array_size(node->children) == 2);

    // The rules are as follows:
    // A context RHS must be a bitfield type.
    // 
    // If the LHS is a bitfield, then we just modify the bits with an AND operation.
    // If the the basetype of LHS is (float, int, float[4], bool or whatever), then the result is of the same type but with an extra dimension with the same length as the RHS.
    // I.e.
    // If RHS is bitrange[8]:
    // An LHS of float[4] would then result in a final type of float[4][8], one for each evaluated bitrange.
    
    // WE MUST BE ABLE TO EVALUATE RHS EXPRESSION DURING COMPILE TIME IN ORDER TO KNOW THE LENGTH


    if (static_check_children(node, ctx)) {
        ast_node_t* lhs = node->children[0];
        ast_node_t* rhs = node->children[1];

        if (rhs->data.type.base_type == TYPE_BITFIELD && rhs->data.type.level > 0) {
            if (!rhs->flags & FLAG_DYNAMIC) {
                eval_context_t eval_ctx = {0};
                eval_ctx.mol = ctx->mol;
                eval_ctx.temp_alloc = ctx->temp_alloc;
                allocate_data(&rhs->data, ctx->mol, ctx->temp_alloc);
                if (evaluate_node(&rhs->data, rhs, &eval_ctx) == 0) {
                    bitfield_t* bf = rhs->data.ptr;
                    uint64_t count = 0;
                    switch(rhs->data.type.level) {
                        case LEVEL_RESIDUE: count = count_residues(*bf, ctx->mol); break;
                        case LEVEL_CHAIN:   count = count_chains(*bf, ctx->mol); break;
                        default: ASSERT(false);
                    }

                    type_info_t new_type = lhs->data.type;
                    if (lhs->data.type.base_type == TYPE_BITFIELD) {
                        new_type.level = rhs->data.type.level;
                    } else {
                        new_type = lhs->data.type;
                        new_type.dim[new_type.len_dim] = (int32_t)count;
                    }

                    node->data.type = new_type;
                    return true;
                } else {
                    create_error(ctx->ir, node->token, "Right hand side of 'in' failed to evaluate at compile time.");
                }
            } else {
                create_error(ctx->ir, node->token, "Right hand side of 'in' cannot be evaluated at compile time.");
            }
        } else {
            create_error(ctx->ir, node->token, "Right hand side of in has an incompatible type (expected residue or chain level structure)");
        }
    }
    return false;
}

static bool static_check(ast_node_t* node, static_check_ctx_t* ctx) {
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
    
    for (uint64_t i = 0; i < md_array_size(node->children); ++i) {
        ast_node_t* child = node->children[i];
        prune_ast_expressions(child);

        if (child && child->type == AST_EXPRESSION) {
            ASSERT(child->children && md_array_size(child->children) == 1);
            node->children[i] = child->children[0];
        }
    }
}

static bool parse_script(md_script_ir* ir) {
    ASSERT(ir);

    bool result = true;

    parse_context_t ctx = {
        .temp_alloc = default_temp_allocator,
        .ir = ir,
        .node = 0,
        .tokenizer = tokenizer_init(ir->str),
    };

    ir->stage = "Parsing";

    // Parse statements until we have exhausted the tokenizer

    token_t tok = {0};
    while (tok = tokenizer_peek_next(&ctx.tokenizer), tok.type != TOKEN_END) {
        const char* beg = tok.str.ptr;

        ir->record_errors = true;   // We reset the error recording flag before each statement
        ast_node_t* node = parse_expression(&ctx);
        tok = tokenizer_consume_next(&ctx.tokenizer);
        expect_token_type(ir, tok, ';');

        const char* end = tok.str.ptr + tok.str.len;
        str_t expr_str = {beg, (uint64_t)(end - beg)};

        if (node) {
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
            result = false;
        }
    }

    return result;
}

static bool static_type_check(md_script_ir* ir, md_molecule* mol) {
    ASSERT(ir);
    ASSERT(mol);

    bool result = true;

    static_check_ctx_t static_check_ctx = {
        .temp_alloc = default_temp_allocator,
        .ir = ir,
        .mol = mol
    };
    ir->stage = "Static Type Checking";

    for (uint64_t i = 0; i < md_array_size(ir->expressions); ++i) {
        ir->record_errors = true;   // We reset the error recording flag before each statement
        if (!static_check(ir->expressions[i]->node, &static_check_ctx)) {
            result = false;
        }
    }

    return result;
}

static void extract_dynamic_evaluation_targets(md_script_ir* ir) {
    ASSERT(ir);

    // Check all expressions for dynamic flag
    for (uint64_t i = 0; i < md_array_size(ir->expressions); ++i) {
        if (ir->expressions[i]->node->flags & FLAG_DYNAMIC) {          
            md_array_push(ir->eval_targets, ir->expressions[i], ir->arena);
        }
    }
}

static bool static_eval_node(ast_node_t* node, md_script_ir* ir, eval_context_t* ctx) {
    ASSERT(node);

    bool result = true;

    if (node->flags & FLAG_DYNAMIC) {
        // This node cannot be resolved at compile time.
        // However, its children might be candidates for compile time evaluation
        // @TODO: Allocate the required data for the top level dynamic nodes here. (size * num_frames)
        if (node->children) {
            for (uint64_t i = 0; i < md_array_size(node->children); ++i) {
                result &= static_eval_node(node->children[i], ir, ctx);
            }
        }
    } else {
        // Allocate persistent data for the 'top' node within the tree/subtree which is saved for later evaluations.
        // Use temp malloc for everything within the subtree
        if (allocate_data(&node->data, ctx->mol, ir->arena)) {
            result &= (evaluate_node(&node->data, node, ctx) > 0);
        } else {
            create_error(ir, node->token, "Could not allocate data for node during static evaluation");
        }
    }

    return result;
}

static bool static_evaluation(md_script_ir* ir, md_molecule* mol) {
    eval_context_t ctx = {0};
    ctx.mol = mol;
    ctx.temp_alloc = default_temp_allocator; // We might need a linear allocator here which can handle big data. Something which does not wrap around and we can reset after each evaluation.
    ir->stage = "Static Evaluation";

    bool result = true;

    // Evaluate every node which is not flagged with FLAG_DYNAMIC and store its value.
    for (uint64_t i = 0; i < md_array_size(ir->expressions); ++i) {
        result &= static_eval_node(ir->expressions[i]->node, ir, &ctx);
    }

    return result;
}

struct md_script_ir* md_script_compile(const char* str, uint64_t len, struct md_molecule* mol, md_allocator_i* alloc) {
    md_script_ir* ir = md_alloc(alloc, sizeof(md_script_ir));
    memset(ir, 0, sizeof(md_script_ir));
    ir->arena = md_arena_allocator_create(alloc, 4096);
    ir->str = copy_str((str_t){str, len}, ir->arena);

    parse_script(ir);
    static_type_check(ir, mol);
    extract_dynamic_evaluation_targets(ir);

#if MD_DEBUG
    save_expressions_to_json(ir->expressions, md_array_size(ir->expressions), (str_t){"tree.json", 9});
#endif
    
    return ir;
}

void md_script_free(md_script_ir* ir) {
    md_arena_allocator_destroy(ir->arena);
}

bool md_script_compile_success(const md_script_ir* ir) {
    // What do we return here?
    return true;
}

uint64_t md_script_get_num_errors(const md_script_ir* ir) {
    if (ir) {
        return md_array_size(ir->errors);
    }
    return 0;
}

const md_script_error* md_script_get_errors(const md_script_ir* ir) {
    if (ir) {
        return ir->errors;
    }
    return NULL;
}