/**
*   The script system allows for evaluation of expressions in the context of molecules.
*   It is a statically typed language with a small set of types and operations.
*   It is used both for 'filtering' i.e. selecting subsets of atoms through textual queries,
*   used in representations among other things.
*   
*   The supported types are
*   A subset of types 
* 
* 
*   There are alot of things to improve here.
*   
*   PERIODS
*   - We should support periods. or at least propagate periods within our system and expose it.
*   - The question is where it belongs, is it part of the type or part of the value?
*   
*   UNITS
*   Units are half arse supported atm. It only works for the trivial cases and the trivial operators (+,-,/,*)
*   - There is no way to specify units for user supplied values.
*   - The mathematical functions abs, min, max, sin, cos, log, pow etc. needs to be overloaded to properly modify units.
*   - The units are resolved during the static check phase, which is concepcually correct. But it is currently very cluncky and not very robust to deal with units in the underlying functions.
*   
*   PROPERTY DATA
*   The property data is a bit of a mess in the API, this should be simplified towards the user. It is fine that we only expose a small fixed quantity of property Types and make their types more concrete.
*    
**/

#include <md_script.h>

#include <md_molecule.h>
#include <md_trajectory.h>
#include <md_filter.h>
#include <md_util.h>
#include <md_xvg.h>
#include <md_edr.h>
#include <md_csv.h>

#include <core/md_common.h>
#include <core/md_compiler.h>
#include <core/md_log.h>
#include <core/md_bitfield.h>
#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_array.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_unit.h>
#include <core/md_spatial_hash.h>
#include <core/md_vec_math.h>
#include <core/md_str_builder.h>
#include <core/md_hash.h>

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <inttypes.h>

#if MD_COMPILER_MSVC
#pragma warning(disable:4063) // Single character tokens not being valid tokens
#endif

#if MD_COMPILER_GCC || MD_COMPILER_CLANG
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wswitch"   // Single character tokens not part of enumeration
#endif
#if MD_COMPILER_CLANG
#pragma GCC diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#endif

#define SCRIPT_IR_MAGIC 0x371bacfe8274910a
#define SCRIPT_EVAL_MAGIC 0x89bca715287bcaff

#define MAX_SUPPORTED_PROC_ARGS 8
#define MAX_NUM_DIMS 4
#define TIMESTAMP_ERROR_MARGIN 1.0e-4

#define SETUP_TEMP_ALLOC(reserve_size) \
    md_vm_arena_t vm_arena; \
    md_vm_arena_init(&vm_arena, reserve_size); \
    md_allocator_i temp_alloc = md_vm_arena_create_interface(&vm_arena);

#define FREE_TEMP_ALLOC() md_vm_arena_free(&vm_arena);

// ################################
// ###   FORWARD DECLARATIONS   ###
// ################################

typedef struct tokenizer_t tokenizer_t;
typedef struct token_t token_t;
typedef struct ast_node_t ast_node_t;
typedef struct procedure_t procedure_t;
typedef struct identifier_t identifier_t;
typedef struct frange_t frange_t;
typedef struct irange_t irange_t;
typedef struct type_info_t type_info_t;
typedef struct data_t data_t;
typedef struct constant_t constant_t;
typedef struct static_backchannel_t static_backchannel_t;
typedef struct parse_context_t parse_context_t;
typedef struct eval_context_t eval_context_t;

typedef union value_t value_t;

// #############################
// ###   TYPE DECLARATIONS   ###
// #############################

typedef enum log_level_t {
    LOG_LEVEL_WARNING,
    LOG_LEVEL_ERROR
} log_level_t;

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
    TOKEN_XOR,
    TOKEN_NOT,
    TOKEN_OF,     // Reserved
    TOKEN_IN,     // Operator for setting the evaluation context.
    TOKEN_OUT,
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
    AST_TABLE,
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
    AST_XOR,
    AST_OR,
    AST_OUT,
    AST_CONTEXT,
    AST_ASSIGNMENT,
    AST_FLATTEN,
    AST_TRANSPOSE,
} ast_type_t;

typedef enum base_type_t {
    TYPE_UNDEFINED = 0,
    TYPE_FLOAT,
    TYPE_INT,
    TYPE_BOOL,
    TYPE_FRANGE,
    TYPE_IRANGE,
    TYPE_BITFIELD,      // Bitfield used to represent a selection of atoms
    TYPE_STRING,
    TYPE_COORDINATE,    // This is a pseudo type which signifies that the underlying type is something that can be interpereted as a coordinate (INT, IRANGE, BITFIELD or float[3])
} base_type_t;

typedef enum flags_t {
    // Procedure Flags
    FLAG_SYMMETRIC_ARGS             = 0x00001, // Indicates that the first two arguments are symmetric, meaning they can be swapped
    FLAG_ARGS_EQUAL_LENGTH          = 0x00002, // Arguments should have the same array length
    FLAG_RET_AND_FIRST_ARG_EQUAL_LENGTH   = 0x00004, // Return type array length matches argument arrays' length
    FLAG_DYNAMIC_LENGTH             = 0x00008, // Return type has a varying length which is not constant over frames
    FLAG_QUERYABLE_LENGTH           = 0x00010, // Marks a procedure as queryable, meaning it can be called with NULL as dst to query the length of the resulting array
    FLAG_STATIC_VALIDATION          = 0x00020, // Marks a procedure as validatable, meaning it can be called with NULL as dst to validate it in a static context during compilation
    FLAG_FLATTEN                    = 0x00040, // Hints that a procedure want flattened bitfields as input, e.g. within, this can be propagated to the arguments during static type checking
    //FLAG_NO_FLATTEN                 = 0x00080, // Hints that a procedure do not want flattened bitfields as input, e.g. com
    
    // Node flags
    FLAG_CONSTANT                   = 0x00100, // Marks the node as constant and should not be modified. It should also have data in its ptr
    FLAG_CONTEXTS_EQUIVALENT        = 0x00200, // Marks the contexts in the node as equivalent, meaning it
    FLAG_COORDINATE                 = 0x00400, // Marks the node as a coordinate, meaning it can be visualized as a coordinate

    // Flags from 0x10000 and upwards are automatically propagated upwards the Nodes of the AST_TREE
    FLAG_DYNAMIC                    = 0x10000, // Indicates that it needs to be reevaluated for every frame of the trajectory (it has a dependency on atomic positions)
    FLAG_SDF                        = 0x20000, // Indicates that the expression involves an sdf computation
    FLAG_VISUALIZE                  = 0x40000, // Hints that a procedure can produce visualizations

    // Flags from 0x10000000 and upwards are automatically propagated to the IR flags
    FLAG_SPATIAL_QUERY              = 0x10000000, // This means that the expression involves some form of spatial query. If we have this, we can pre-compute a spatial hash acceleration structure and provide it through the context
} flags_t;

typedef enum eval_flags_t {
    //EVAL_FLAG_STATIC_TYPE_CHECK     = 0x1,
    //EVAL_FLAG_EVALUATE              = 0x2,
    EVAL_FLAG_FLATTEN               = 0x100,
    EVAL_FLAG_NO_STATIC_EVAL        = 0x200,
    EVAL_FLAG_NO_LENGTH_CHECK       = 0x400,
} eval_flags_t;

// Propagate upwards to parent nodes within the AST tree
#define FLAG_AST_PROPAGATION_MASK 0xFFFF0000U
#define FLAG_IR_PROPAGATION_MASK  0xF0000000U

struct token_t {
    token_type_t type;
    int line_beg;
    int line_end;
    int col_beg;
    int col_end;
    int beg;
    int end;
    str_t str;
};

struct frange_t {
    float beg;
    float end;
};

struct irange_t {
    int beg;
    int end;
};

struct type_info_t {
    base_type_t  base_type;
    int dim[MAX_NUM_DIMS]; // The dimensionality of a single element, each entry in dim represents a multidimensional length
};

// There is a conceptual twist here if len should belong to the type_info_t or if it should belong to data_t
// We need to differentiate between different lengths of arrays with the same type in our system, so for now it should probably reside within the type_info
// I.e. the array length should be part of the fundamental type. This only strenghtens the type-system

// This is the struct which is passed as arguments into the procedures
struct data_t {
    type_info_t type;
    void*       ptr;    // Pointer to the data
    size_t      size;   // Size in bytes of data (This we want to determine during the static check so we can allocate the data when evaluating)

    frange_t    value_range;
    md_unit_t   unit;
};

// This is data stored directly in the nodes of the AST tree
union value_t {
    str_t           _string;
    float           _float;
    int             _int;
    frange_t        _frange;
    irange_t        _irange;
    bool            _bool;
    md_bitfield_t   _bitfield;
};

struct identifier_t {
    str_t       name;
    ast_node_t* node;    // This is the node to evaluate in order to generate the data for the identifier
    data_t*     data;    // This is the data to fill in...
};

struct constant_t {
    str_t     name;
    float     value;
    md_unit_t unit;
};

// The backchannel is used to provide additional information of the data during static analysis.
// Ideally this should perhaps reside together with the data_t struct, but for now we keep it separate.
// The reason is because we do not pass in data_t during static analysis. which may be a mistake.
typedef struct static_backchannel_t {
    flags_t   flags;
    md_unit_t unit;
    frange_t  value_range;
} static_backchannel_t;

typedef struct eval_context_t {
    md_script_ir_t* ir;
    const md_molecule_t* mol;
    const md_bitfield_t* mol_ctx;       // The atomic bit context in which we perform the operation, this can be null, in that case we are not limited to a smaller context and the full molecule is our context
    const md_trajectory_i* traj;

    md_vm_arena_t* temp_arena;          // For allocating transient data  Raw interface, prefer this
    md_allocator_i* temp_alloc;         // For allocating transient data  (Generic interface to temp_arena)
    md_allocator_i* alloc;              // For allocating persistent data (for the duration of the evaluation)

    // Contextual information for static checking 
    token_t  op_token;                  // Token for the operation which is evaluated
    token_t* arg_tokens;                // Tokens to arguments for contextual information when reporting errors
    flags_t* arg_flags;                 // Flags of arguments
    identifier_t* identifiers;          // Evaluated identifiers for references
                                      
    md_script_vis_t* vis;               // These are used when calling a procedure flagged with the VISUALIZE flag so the procedure can fill in the geometry
    uint32_t vis_flags;                 // Visualization flags
    md_bitfield_t* vis_structure;       // Current visualization structure to mark atoms in

    md_trajectory_frame_header_t* frame_header;

    struct {
        md_trajectory_frame_header_t* header;
        const float* x;
        const float* y;
        const float* z;
    } initial_configuration;

    md_spatial_hash_t* spatial_hash;

    // During evaluations, whenever we hit a array subscript operator, we store the indices
    // In order to propagate what will be visible and or used by the array subscript operator
    // This allows us to only visualize parts of an expression
    md_array(irange_t) subscript_ranges;

    // This is information which is only present in the static check and serves as a backchannel for the procedure in order to deduce value range, units, flags and such
    static_backchannel_t* backchannel;

    uint32_t eval_flags; // Flags containing information in what context the evaluation is occurring
} eval_context_t;

struct procedure_t {
    str_t name;
    type_info_t return_type;
    size_t num_args;
    type_info_t arg_type[MAX_SUPPORTED_PROC_ARGS];
    int  (*proc_ptr)(data_t*, data_t[], eval_context_t*); // Function pointer to the "C" procedure    
    flags_t flags;
};

typedef struct table_t {
    str_t name;
    size_t num_fields;
    size_t num_values;
    
    md_array(str_t) field_names;
    md_array(md_unit_t) field_units;
    md_array(md_array(double)) field_values;
} table_t;

void table_push_field_d(table_t* table, str_t name, md_unit_t unit, const double* data, size_t count, md_allocator_i* alloc) {
    (void)count;
    ASSERT(count == table->num_values);
    md_array_push(table->field_names, str_copy(name, alloc), alloc);
    md_array_push(table->field_units, unit, alloc);
    md_array(double) field_data = md_array_create(double, table->num_values, alloc);
    MEMCPY(field_data, data, sizeof(double) * table->num_values);
    md_array_push(table->field_values, field_data, alloc);
    table->num_fields = md_array_size(table->field_names);
}

void table_push_field_f(table_t* table, str_t name, md_unit_t unit, const float* data, size_t count, md_allocator_i* alloc) {
    (void)count;
    ASSERT(count == table->num_values);
    md_array_push(table->field_names, str_copy(name, alloc), alloc);
    md_array_push(table->field_units, unit, alloc);
    md_array(double) field_data = md_array_create(double, table->num_values, alloc);
    for (size_t i = 0; i < table->num_values; ++i) {
        field_data[i] = data[i];
    }
    md_array_push(table->field_values, field_data, alloc);
    table->num_fields = md_array_size(table->field_names);
}

struct ast_node_t {
    // @TODO, @PERF: Make specific types for each type of Ast_node, where the first member is the type which we can use and cast by.

    ast_type_t          type;    
    flags_t             flags;      // Flags for node (set during static type checking)
    token_t             token;      // Corresponding token from which the node was created (used to tie errors back into src)

    // PAYLOAD
    ast_node_t**        children;   // For AND, OR, NOT, procedure arguments etc.
    value_t             value;      // Scalar values for storing data directly inside the node
    data_t              data;       // Structure for passing as argument into procedure. (Holds pointer, length and type)

    // PROCEDURE
    const procedure_t*  proc;       // Procedure reference
    uint32_t            proc_flags; // Additional flags used for e.g. swapping arguments which may be required for the procedure call

    // IDENTIFIER
    str_t               ident;      // Identifier reference (also used for procedure name lookup)

    // ARRAY SUBSCRIPT
    irange_t            subscript_ranges[MAX_NUM_DIMS]; // The ranges for the array subscript operator
    size_t              subscript_dim;

    // TABLE
    const table_t*      table;
    md_array(int)       table_field_indices;  // Indices into the table fields which are used for this node

    // CONTEXT
    size_t              num_contexts;   // Number of arguments for procedure calls
    type_info_t*        lhs_context_types; // For static type checking
    const md_bitfield_t* contexts;  // Since contexts have to be statically known at compile time, we store a reference to it.
};

typedef struct tokenizer_t {
    str_t   str;
    int     cur;
    int     line;
    int     line_offset;  // offset to the current line
} tokenizer_t;

struct md_script_ir_t {
    uint64_t magic;
    uint64_t fingerprint;

    // We could use a direct raw interface here to save some function pointer indirections
    struct md_allocator_i *arena;

    str_t str;      // Original string containing the 'source'
    uint32_t flags; // special flags which have been assigned during compile time
    
    md_array(ast_node_t*) nodes;                        // All nodes in the AST tree

    md_array(ast_node_t*) parsed_expressions;           // List of parsed expressions in the script
    md_array(ast_node_t*) type_checked_expressions;     // List of expressions which passed static type checking
    md_array(ast_node_t*) eval_targets;                 // List of dynamic expressions which needs to be evaluated per frame

    // These form the basis for 'static' common subexpression elimination
    md_array(data_t)    static_expression_data;         // Static expressions which are constant and evaluated during static analysis.
    md_array(uint64_t)  static_expression_hash;         // hash of expression
    md_array(str_t)     static_expression_str;          // string for debugging
    
    md_array(identifier_t)  identifiers;                // List of identifiers, notice that the data in a const context should only be used if it is flagged as
    md_array(table_t) tables;                           // List of tables which are used in the script

    md_array(md_script_property_flags_t) property_flags;    // List of property infos
    md_array(const ast_node_t*)          property_nodes;    // List of property nodes;
    md_array(str_t)                      property_names;    // List of property names

    // These are the final products which can be read through the public part of the structure
    md_array(md_log_token_t)        warnings;
    md_array(md_log_token_t)        errors;
    md_array(md_script_vis_token_t) vis_tokens;
    md_array(str_t)                 identifier_names;

    const char* stage;      // This is just to supply a context for the errors i.e. in which compilation stage the error occured
    bool record_log;        // This is to toggle if new errors should be registered... We don't want to flood the errors
    bool compile_success;
};

struct md_script_eval_t {
    uint64_t magic;
    uint64_t ir_fingerprint;

    struct md_allocator_i *arena;
    volatile bool interrupt;

    size_t        frame_count;
    md_bitfield_t frame_mask;
    md_mutex_t    frame_lock;

    md_array(str_t)                     property_names;
    md_array(md_script_property_data_t) property_data;
    md_array(volatile uint32_t)         property_dist_count;   // Counters property distributions
    md_array(md_mutex_t)                property_dist_mutex;   // Protect the data when writing to it in a threaded context (Distributions)
};

struct parse_context_t {
    md_script_ir_t* ir;
    tokenizer_t*    tokenizer;
    ast_node_t*     node;   // This represents the current root node in the tree
    md_allocator_i* temp_alloc;
};


// ##########################
// ###   CORE FUNCTIONS   ###
// ##########################

static uint64_t generate_fingerprint() {
    return md_time_current();
}

static uint32_t operator_precedence(ast_type_t type) {
    switch(type) {
    case AST_EXPRESSION:
    case AST_PROC_CALL:
    case AST_CONSTANT_VALUE:
    case AST_TABLE:
    case AST_IDENTIFIER:
    case AST_ARRAY:
    case AST_ARRAY_SUBSCRIPT:
    case AST_FLATTEN:
    case AST_TRANSPOSE:
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
    case AST_XOR:
        return 8;
    case AST_OR:
        return 9;
    case AST_OUT:
        return 10;
    case AST_CONTEXT:
        return 11;
    case AST_ASSIGNMENT:
        return 12;
    default:
        ASSERT(false);
    }
    return 0;
}

// This is used to get the dim of an element in an array for example
static void dim_shift_left(int dim[MAX_NUM_DIMS]) {
    for (int i = 0; i < MAX_NUM_DIMS - 1; ++i) {
        dim[i] = dim[i+1];
    }
    dim[MAX_NUM_DIMS-1] = 0;
}

// The semantics are a bit weird, but the idea is that if we shift right,
// that means we are creating an array of length 1, which is the same as a scalar
static void dim_shift_right(int dim[MAX_NUM_DIMS]) {
    for (int i = MAX_NUM_DIMS - 1; i > 0; --i) {
        dim[i] = dim[i-1];
    }
    dim[0] = 1;
}

// returns the number of non zero dimensions
static int dim_ndims(const int dim[MAX_NUM_DIMS]) {
    for (int i = 0; i < 4; ++i) {
        if (dim[i] == 0) return i;
    }
    return 4;
}

// Return the product of the dimensions (excluding zero dimensions)
static int dim_size(const int dim[MAX_NUM_DIMS]) {
	int size = 1;
	for (int i = 0; i < MAX_NUM_DIMS; ++i) {
    	if (dim[i] == 0) break;
    	size *= dim[i];
    }
	return size;
}

static void dim_flatten(int dim[MAX_NUM_DIMS]) {
    int size = dim_size(dim);
    MEMSET(dim, 0, sizeof(int) * MAX_NUM_DIMS);
    dim[0] = size;
}

static void dim_prune_leading_ones(int dim[MAX_NUM_DIMS]) {
    int num_dim = dim_ndims(dim);
    if (num_dim < 2) return;
    for (int i = 0; i < num_dim - 1; ++i) {
        if (dim[i] == 1) {
            dim_shift_left(dim);
            num_dim--;
        }
    }
}

static bool type_info_equal(type_info_t a, type_info_t b) {
    return MEMCMP(&a, &b, sizeof(type_info_t)) == 0;
}

static bool is_undefined_type(type_info_t ti) {
    return MEMCMP(&ti, &(type_info_t){0}, sizeof(type_info_t)) == 0;
}

static bool is_array(type_info_t ti) {
    return ti.dim[0] > 1;
}

static bool is_scalar(type_info_t ti) {
    return ti.dim[0] == 1 && (ti.dim[1] == 0 || ti.dim[1] == 1);
}

static bool is_variable_length(type_info_t ti) {
    return ti.dim[0] == -1;
}

static int type_info_array_len(type_info_t ti) {
    return ti.dim[0];
}

static size_t element_count(data_t arg) {
    return arg.type.dim[0];
}

static type_info_t type_info_element_type(type_info_t ti) {
    ti.dim[0] = 1;
    return ti;
}

// Size of the base type structure
static size_t base_type_element_byte_size(base_type_t type) {
    switch (type) {
    case TYPE_FLOAT:    return sizeof(float);
    case TYPE_INT:      return sizeof(int);
    case TYPE_BOOL:     return sizeof(bool);        
    case TYPE_FRANGE:   return sizeof(frange_t);
    case TYPE_IRANGE:   return sizeof(irange_t);
    case TYPE_BITFIELD: return sizeof(md_bitfield_t);
    case TYPE_STRING:   return sizeof(str_t);
    case TYPE_UNDEFINED:
    default:            return 0;
    }
}

static bool is_type_dim_compatible(type_info_t from, type_info_t to) {

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
    //             ^                        
    // 
    // float[3][1] should match float[3][0]


    // Is this always the case??? (Not if we want to support the last special case
    //if (from.len_dim != to.len_dim) return false;

    for (int i = 0; i < MAX_NUM_DIMS; ++i) {
        if (i == 0) {
            if (from.dim[i] == to.dim[i]) continue;
            if (from.dim[i] == -1 && to.dim[i] > 0) continue;
            if (from.dim[i] > 0 && to.dim[i] == -1) continue;
            return false;
        }
        else {
            if (from.dim[i] != to.dim[i]) return false;
        }
    }
    return true;
}

static bool is_type_equivalent(type_info_t a, type_info_t b) {
    if (a.base_type != b.base_type) return false;
    for (int i = 0; i < MAX_NUM_DIMS; ++i) {
        if (a.dim[i] == 0) a.dim[i] = 1;
        if (b.dim[i] == 0) b.dim[i] = 1;
    }
    return MEMCMP(&a, &b, sizeof(type_info_t)) == 0;
}

static bool is_type_position_compatible(type_info_t type) {
    switch(type.base_type) {
        case TYPE_INT: return true;
        case TYPE_IRANGE: return true;
        case TYPE_BITFIELD: return true;
        case TYPE_FLOAT: 
            return (type.dim[0] == 3 && type.dim[1] == 0) || (type.dim[1] == 3 && type.dim[2] == 0);
        default: return false;
    }
}

static bool is_type_directly_compatible(type_info_t from, type_info_t to) {
    if (to.base_type == TYPE_COORDINATE && is_type_position_compatible(from)) {
        to.base_type = from.base_type;
        if (from.base_type == TYPE_FLOAT) {
            if (from.dim[1] != 0) {
                to.dim[0] = from.dim[0];
                to.dim[1] = 3;
            } else {
                to.dim[0] = 3;
            }
        }
    }

    int from_ndim = dim_ndims(from.dim);
    int to_ndim   = dim_ndims(to.dim);
    int max_ndim  = MAX(from_ndim, to_ndim);
    if (from.dim[0] != -1) {
        for (int i = 0; i < max_ndim - from_ndim; ++i) {
            dim_shift_right(from.dim);
        }
    }
    if (to.dim[0] != -1) {
        for (int i = 0; i < max_ndim - to_ndim; ++i) {
            dim_shift_right(to.dim);
        }
    }
    for (int i = from_ndim; i < MAX_NUM_DIMS; ++i) {
        if (from.dim[i] == 0) from.dim[i] = 1;
    }
    for (int i = to_ndim; i < MAX_NUM_DIMS; ++i) {
        if (to.dim[i] == 0) to.dim[i] = 1;
    }

    if (type_info_equal(from, to)) return true;

    if (from.base_type == to.base_type) {
        // This is essentially a logical XOR, we only want to support this operation if we have one unspecified array dimension.
        if (!is_variable_length(from) && is_variable_length(to)) {
            return is_type_dim_compatible(from, to);
        }
    }

    return false;
}

static bool compare_type_info_array(const type_info_t a[], const type_info_t b[], size_t num_arg_types) {
    for (size_t i = 0; i < num_arg_types; ++i) {
        if (!is_type_directly_compatible(a[i], b[i])) {
            return false;
        }
    }
    return true;
}

// Returns the count of elements for this type, i.e an array of [1][4][4] would return 16
// NOTE: It does not include the last dimension which encodes the length of the array
static size_t type_info_element_stride_count(type_info_t ti) {
    if (ti.dim[1] == 0) return 1;
    size_t stride = ti.dim[1];
    for (int32_t i = 2; i < MAX_NUM_DIMS; ++i) {
		if (ti.dim[i] == 0) break;
		stride *= ti.dim[i];
	}
    /*
    if (ti.len_dim == 0) return 1;
    ASSERT(ti.dim[0] >= 0);
    size_t stride = ti.dim[0];
    for (int32_t i = 1; i < ti.len_dim; ++i) {
        ASSERT(ti.dim[i] >= 0);
        stride *= ti.dim[i];
    }
    */
    return stride;
}

static size_t type_info_element_byte_stride(type_info_t ti) {
    return type_info_element_stride_count(ti) * base_type_element_byte_size(ti.base_type);
}

static size_t type_info_total_byte_size(type_info_t ti) {
    return type_info_element_byte_stride(ti) * type_info_array_len(ti);
}

static size_t type_info_total_element_count(type_info_t ti) {
    return type_info_element_stride_count(ti) * type_info_array_len(ti);
}

static size_t bitfield_byte_size(size_t num_bits) {
    return DIV_UP(num_bits, 64) * sizeof(int64_t);
}

static bool timestamps_approx_equal(double t0, double t1) {
    return (t0 - t1) < TIMESTAMP_ERROR_MARGIN;
}

static bool allocate_data(data_t* data, type_info_t type, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(!is_undefined_type(type));
    ASSERT(alloc);
    
    // The convention should be that we only allocate data for known lengths
    const int64_t array_len = type_info_array_len(type);
    ASSERT(array_len > -1);

    // Do the base type allocation (array)
    const int64_t bytes = type_info_element_byte_stride(type) * array_len;

    if (bytes > 0) {
        data->ptr = md_alloc(alloc, bytes);
        if (!data->ptr) {
            MD_LOG_ERROR("Failed to allocate data in script!");
            return false;
        }
        // @NOTE (Robin): This could be a very wasteful operation most times as the memory will be overwritten anyways.
        MEMSET(data->ptr, 0, bytes);
    }
    data->size = bytes;
    data->type = type;
    data->value_range = (frange_t){0, 0};

    if (data->unit.base.raw_bits == 0 && data->unit.mult == 0) {
        data->unit.base.raw_bits = 0;
        data->unit.mult = 1.0;
    }

    if (type.base_type == TYPE_BITFIELD) {
        md_bitfield_t* bf = data->ptr;
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
            md_bitfield_t* bf = data->ptr;
            const int64_t num_elem = element_count(*data);
            for (int64_t i = 0; i < num_elem; ++i) {
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

    if (dst->type.base_type == TYPE_BITFIELD) {
        const uint64_t num_elem = element_count(*dst);

        md_bitfield_t* dst_bf = dst->ptr;
        const md_bitfield_t* src_bf = src->ptr;
        // Set bit pointers into destination memory
        for (uint64_t i = 0; i < num_elem; ++i) {
            md_bitfield_copy(&dst_bf[i], &src_bf[i]);
        }
    } else {
        MEMCPY(dst->ptr, src->ptr, src->size);
    }
}



static bool is_operator(ast_type_t type) {
    return (type == AST_ADD || type == AST_SUB || type == AST_MUL || type == AST_DIV || type == AST_UNARY_NEG ||
            type == AST_AND || type == AST_OR || type == AST_XOR || type == AST_NOT ||
            type == AST_EQ || type == AST_NE || type == AST_LE || type == AST_GE || type == AST_LT || type == AST_GT);
}

static bool is_bitwise_operator(ast_type_t type) {
    return (type == AST_AND || type == AST_OR || type == AST_XOR || type == AST_NOT);
}

static bool is_number(token_type_t type) {
    return type == TOKEN_INT || type == TOKEN_FLOAT;
}

// Returns if a token should be considered an operand
// Mainly used to disambiguate unary operator from binary operator
// Then we need to look at the previous token
static bool is_operand(token_type_t type) {
    switch (type) {
    case ')':
    case ']':
    case ':':
    case TOKEN_IDENT:
    case TOKEN_OF:
    case TOKEN_IN:
    case TOKEN_FLOAT:
    case TOKEN_INT:
    case TOKEN_STRING:
        return true;
    case TOKEN_EQ:
    case TOKEN_AND:
    case TOKEN_OR:
    case TOKEN_XOR:
    case TOKEN_GE:
    case TOKEN_NOT:
    case TOKEN_LE:
    case '+':
    case '>':
    case '*':
    case '-':
    case '=':
    case '/':
    case '<':
    case '[':
    case '(':
    case ',':
    case ';':
    case TOKEN_UNDEF:
    case TOKEN_END:
        return false;
    default:
        ASSERT(false);
    }
    return false;
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
    case TOKEN_XOR: return "xor";
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

static token_t concat_tokens(token_t a, token_t b) {
    token_t t;
    if (a.line_beg < b.line_beg) {
        t.col_beg  = a.col_beg;
        t.col_end  = b.col_end;
        t.line_beg = a.line_beg;
        t.line_end = b.line_end;
    } else if (b.line_beg < a.line_beg) {
        t.col_beg  = b.col_beg;
        t.col_end  = a.col_end;
        t.line_beg = b.line_beg;
        t.line_end = a.line_end;
    } else {
        t.col_beg  = MIN(a.col_beg, b.col_beg);
        t.col_end  = MAX(a.col_end, b.col_end);
        t.line_beg = a.line_beg;
        t.line_end = b.line_end;
    }
    t.type = a.beg < b.beg ? a.type : b.type;
    t.beg = MIN(a.beg, b.beg);
    t.end = MAX(a.end, b.end);

    t.str.ptr = MIN(str_beg(a.str), str_beg(b.str));
    t.str.len = MAX(str_end(a.str), str_end(b.str)) - str_beg(t.str);
    return t;
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

#define LOG_WARNING(ir, tok, fmt, ...)  create_log_token(LOG_LEVEL_WARNING, ir, tok, fmt, ##__VA_ARGS__)
#define LOG_ERROR(ir, tok, fmt, ...)    create_log_token(LOG_LEVEL_ERROR, ir, tok, fmt, ##__VA_ARGS__)

static void create_log_token(log_level_t level, md_script_ir_t* ir, token_t token, const char* format, ...) {
    ASSERT(ir);
    if (!ir->record_log) return;

    if (level >= LOG_LEVEL_ERROR) {
        // Stop after first error
        ir->record_log = false;
    }
    
    char buf[1024];
    va_list args;
    va_start(args, format);
    int len = vsnprintf(buf, ARRAY_SIZE(buf), format, args);
    va_end(args);

    md_log_token_t log_tok = {
        .range = {
            .beg = token.beg,
            .end = token.end,
        },
        .text = str_copy_cstrn(buf, len, ir->arena),
    };

    switch (level) {
    case LOG_LEVEL_WARNING:
        md_array_push(ir->warnings, log_tok, ir->arena);
        break;
    case LOG_LEVEL_ERROR:
        md_array_push(ir->errors, log_tok, ir->arena);
        break;
    default:
        break;
    }

    /*
    md_logf(MD_LOG_TYPE_DEBUG, "%s line %i: %s", ir->stage, token.line_beg, buffer);

    if (token.str.ptr && token.str.len) {
        MEMSET(buffer, 0, ARRAY_SIZE(buffer));
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
        static const char long_ass_carret[] = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
        md_logf(MD_LOG_TYPE_DEBUG, ""STR_FMT"", (end-beg), beg);
        md_logf(MD_LOG_TYPE_DEBUG, "%*s^"STR_FMT"", (token.str.ptr - beg), "", token.str.len-1, long_ass_carret);
    }
    */
}

// Include all procedures, operators and defines
#include "md_script_functions.inl"

// ############################
// ###   HELPER FUNCTIONS   ###
// ############################

// Returns if the value type can be operated on using logical operators
static inline bool is_value_type_logical_operator_compatible(base_type_t type) {
    return type == TYPE_BOOL || type == TYPE_BITFIELD;
}

static inline bool is_identifier_static_procedure(str_t ident) {
    if (str_eq(ident, STR_LIT("import"))) {
        return true;   
    }
    else if (str_eq(ident, STR_LIT("flatten"))) {
        return true;
    }
    return false;
}

static inline bool is_identifier_procedure(str_t ident) {
    for (size_t i = 0; i < ARRAY_SIZE(procedures); ++i) {
        if (str_eq(ident, procedures[i].name)) {
            return true;
        }
    }
    return false;
}

static ast_node_t* create_node(md_script_ir_t* ir, ast_type_t type, token_t token) {
    ast_node_t* node = md_alloc(ir->arena, sizeof(ast_node_t));
    MEMSET(node, 0, sizeof(ast_node_t));
    node->type = type;
    node->token = token;
    node->data.unit = md_unit_none();
    md_array_push(ir->nodes, node, ir->arena);
    return node;
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

static bool is_type_implicitly_convertible(type_info_t from, type_info_t to) {
    for (size_t i = 0; i < ARRAY_SIZE(casts); ++i) {
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

static procedure_match_result_t find_cast_procedure(type_info_t from, type_info_t to) {

    procedure_match_result_t res = {0};
    uint32_t lowest_cost = ~0U;

    // In the future, we might need to do two casts to actually get to the proper type.

    for (size_t i = 0; i < ARRAY_SIZE(casts); ++i) {
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

static uint32_t compute_cost(const procedure_t* proc, const type_info_t arg_types[]) {
    uint32_t cost = 0;
    for (size_t j = 0; j < proc->num_args; ++j) {
        if (type_info_equal(arg_types[j], proc->arg_type[j])) {
            // No conversion needed for this argument (0 cost)
        }
        else if (is_type_directly_compatible(arg_types[j], proc->arg_type[j])) {
            // It is not a perfect match, but we can do a direct conversion (probably variable array length)
            cost += 1;
        }
        else if (is_type_implicitly_convertible(arg_types[j], proc->arg_type[j])) {
            // Implicit conversions are more expensive than direct conversions
            cost += 4;
        }
        else {
            // We are smoked.. This particular matching procedure cannot be resolved using implicit conversions...
            cost = 0xFFFFFFFFU;
            break;
        }
    }
    return cost;
}

static procedure_match_result_t find_procedure_supporting_arg_types_in_candidates(str_t name, const type_info_t arg_types[], size_t num_arg_types, procedure_t* candidates, size_t num_cantidates, bool allow_implicit_conversions) {
    procedure_match_result_t res = {0};

    for (size_t i = 0; i < num_cantidates; ++i) {
        procedure_t* proc = &candidates[i];
        if (str_eq(proc->name, name)) {
            if (num_arg_types == proc->num_args) {
                if (compare_type_info_array(arg_types, proc->arg_type, num_arg_types)) {
                    res.success = true;
                    res.procedure = proc;
                    return res;
                }
                else if (proc->flags & FLAG_SYMMETRIC_ARGS) {
                    // @TODO: Does this make sense for anything else than two arguments???
                    // The only exception here is the RDF function which should be symmetric on the first two, and then have the third cutoff parameter.
                    ASSERT(proc->num_args >= 2);
                    type_info_t swapped_args[2] = {arg_types[1], arg_types[0]};
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

        for (size_t i = 0; i < num_cantidates; ++i) {
            procedure_t* proc = &candidates[i];
            if (str_eq(proc->name, name)) {
                if (num_arg_types == proc->num_args) {
                    // Same name and same number of args...

                    uint32_t cost = compute_cost(proc, arg_types);
                    if (cost < best_cost) {
                        best_cost = cost;
                        best_proc = proc;
                        best_flags = 0;
                    }

                    if (proc->flags & FLAG_SYMMETRIC_ARGS) {
                        ASSERT(proc->num_args >= 2);
                        // Test if we can get a (better) result by swapping the arguments
                        type_info_t swapped_args[2] = {arg_types[1], arg_types[0]};

                        cost = compute_cost(proc, swapped_args);
                        if (cost < best_cost) {
                            best_cost = cost;
                            best_proc = proc;
                            best_flags = FLAG_SYMMETRIC_ARGS;
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

static procedure_match_result_t find_procedure_supporting_arg_types(str_t name, const type_info_t arg_types[], size_t num_arg_types, bool allow_implicit_conversions) {
    return find_procedure_supporting_arg_types_in_candidates(name, arg_types, num_arg_types, procedures, ARRAY_SIZE(procedures), allow_implicit_conversions);
}

static procedure_match_result_t find_operator_supporting_arg_types(ast_type_t op, const type_info_t arg_types[], size_t num_arg_types, bool allow_implicit_conversions) {
    // Map marker type to string which we use to identify operator procedures
    str_t name = {0};
    switch(op) {
    case AST_ADD: name = STR_LIT("+");   break;
    case AST_UNARY_NEG:
    case AST_SUB: name = STR_LIT("-");   break;
    case AST_MUL: name = STR_LIT("*");   break;
    case AST_DIV: name = STR_LIT("/");   break;
    case AST_AND: name = STR_LIT("and"); break;
    case AST_OR:  name = STR_LIT("or");  break;
    case AST_XOR: name = STR_LIT("xor");  break;
    case AST_NOT: name = STR_LIT("not"); break;
    case AST_EQ:  name = STR_LIT("==");  break;
    case AST_NE:  name = STR_LIT("!=");  break;
    case AST_LT:  name = STR_LIT("<");   break;
    case AST_GT:  name = STR_LIT(">");   break;
    case AST_LE:  name = STR_LIT("<=");  break;
    case AST_GE:  name = STR_LIT(">=");  break;

    default:
        ASSERT(false);
    }

    return find_procedure_supporting_arg_types_in_candidates(name, arg_types, num_arg_types, operators, ARRAY_SIZE(operators), allow_implicit_conversions);
}

static constant_t* find_constant(str_t name) {
    for (size_t i = 0; i < ARRAY_SIZE(constants); ++i) {
        if (str_eq(constants[i].name, name)) {
            return &constants[i];
        }
    }
    return NULL;
}

static identifier_t* find_identifier(str_t name, identifier_t* identifiers, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        if (str_eq(name, identifiers[i].name)) {
            return &identifiers[i];
        }
    }
    return NULL;
}

static identifier_t* get_identifier(md_script_ir_t* ir, str_t name) {
    return find_identifier(name, ir->identifiers, md_array_size(ir->identifiers));
}

static identifier_t* create_identifier(md_script_ir_t* ir, str_t name) {
    if (find_constant(name)) {
        return NULL;
    }

    if (get_identifier(ir, name)) {
        return NULL;
    }

    identifier_t ident = {
        .name = str_copy(name, ir->arena),
        .node = 0,
        .data = 0,
    };

    return md_array_push(ir->identifiers, ident, ir->arena);
}

static inline bool is_token_type_comparison(token_type_t type) {
    return type == '<' || type == TOKEN_LE || type == '>' || type == TOKEN_GE || type == TOKEN_EQ;
}

static bool expect_token_type(md_script_ir_t* ir, token_t token, token_type_t type) {
    if (token.type != type) {
        LOG_ERROR(ir, token, "Unexpected token '"STR_FMT"', expected token '%s'", token.str.len, token.str.ptr, get_token_type_str(type));
        return false;
    }
    return true;
}

// #####################
// ###   TOKENIZER   ###
// #####################

static token_t tokenizer_get_next_from_buffer(tokenizer_t* tokenizer) {
    token_t token = {0};

    if (tokenizer->cur >= (int)tokenizer->str.len) {
        tokenizer->cur = (int)tokenizer->str.len;
        token.type = TOKEN_END;
        token.line_beg = tokenizer->line;
        token.line_end = tokenizer->line;
        token.col_beg  = tokenizer->cur - tokenizer->line_offset;
        token.col_end  = tokenizer->cur - tokenizer->line_offset;
        token.beg      = tokenizer->cur;
        token.end      = tokenizer->cur;
        return token;
    }
    
    int        line = tokenizer->line;
    const char* buf = tokenizer->str.ptr;
    const int   len = (int)tokenizer->str.len;
    for (int i = (int)tokenizer->cur; i != len; ++i) {
        int j = 0;
        
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
            const int n = j - i;
            if (n == 2) {
                str_t str = {buf+i, 2};
                if (str_eq(str, STR_LIT("or"))) {
                    token.type = TOKEN_OR;
                }
                else if (str_eq(str, STR_LIT("in"))) {
                    token.type = TOKEN_IN;
                }
                else if (str_eq(str, STR_LIT("of"))) {
                    token.type = TOKEN_OF;
                }
            }
            else if (n == 3) {
                str_t str = {buf+i, 3};
                if (str_eq(str, STR_LIT("and"))) {
                    token.type = TOKEN_AND;
                }
                else if (str_eq(str, STR_LIT("not"))) {
                    token.type = TOKEN_NOT;
                }
                else if (str_eq(str, STR_LIT("xor"))) {
                    token.type = TOKEN_XOR;
                }
                else if (str_eq(str, STR_LIT("out"))) {
                    token.type = TOKEN_OUT;
                }
            }
            
            if (token.type == TOKEN_UNDEF) {
                token.type = TOKEN_IDENT;
            }
        }
        
        // Numeric literal
        else if (is_digit(buf[i])) {
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
            typedef struct symbol_t {
                const char* str;
                token_type_t type;
            } symbol_t;
            
            static const symbol_t symbols_2[] = {
                {"<=", TOKEN_LE},
                {">=", TOKEN_GE},
                {"==", TOKEN_EQ}
            };
            
            static const char symbols_1[] = {
                '(', ')', '[', ']', '{', '}', '<', '>', '*', '+', '-', '/', '=', ':', ';', ',', '.'
            };
            
            // Potential symbol count in buf to match (only care about 1 or 2)
            const int n = (i + 1 < len && is_symbol(buf[i + 1])) ? 2 : 1;
            
            // Match buf + i against the longest accepted symbol that we can find.
            if (n == 2) {
                str_t str = {buf + i, 2};
                for (size_t k = 0; k < ARRAY_SIZE(symbols_2); ++k) {
                    if (str_eq(str, (str_t){symbols_2[k].str, 2})) {
                        token.type = symbols_2[k].type;
                        j = i + 2;
                        break;
                    }
                }
            }
            if (!token.type) {
                for (size_t k = 0; k < ARRAY_SIZE(symbols_1); ++k) {
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
            token.line_beg = line;
            token.line_end = line;
            token.col_beg = i - tokenizer->line_offset;
            token.col_end = j - tokenizer->line_offset;
            token.beg = i;
            token.end = j;
            break;
        }

        if (!token.type && i >= len - 1) {
            // We are going to hit the end in next iteration and have nothing to show for.
            tokenizer->cur = (int)tokenizer->str.len;
            token.type = TOKEN_END;
            token.line_beg = line;
            token.line_end = line;
            token.col_beg  = i - tokenizer->line_offset;
            token.col_end  = i - tokenizer->line_offset;
            token.beg = i;
            token.end = i;
            break;
        }
    }

    return token;
}

static token_t tokenizer_consume_next(tokenizer_t* tokenizer) {
    token_t token = tokenizer_get_next_from_buffer(tokenizer);
    if (token.str.ptr && token.str.len) {
        // Advance current of tokenizer
        tokenizer->cur = (int)(token.str.ptr + token.str.len - tokenizer->str.ptr);
        tokenizer->line = token.line_end;
    }
    return token;
}

static token_t tokenizer_peek_next(tokenizer_t* tokenizer) {
    token_t token = tokenizer_get_next_from_buffer(tokenizer);
    return token;
}

static token_t tokenizer_consume_until_type(tokenizer_t* tokenizer, token_type_t* types, size_t num_types) {
    token_t token = {0};
    while (token = tokenizer_get_next_from_buffer(tokenizer), token.type != TOKEN_END) {
        for (size_t i = 0; i < num_types; ++i) {
            if (token.type == types[i]) goto end;
        }
        if (token.str.ptr && token.str.len) {
            // Advance current of tokenizer
            tokenizer->cur = (int)(token.str.ptr + token.str.len - tokenizer->str.ptr);
            tokenizer->line = token.line_end;
        }
    }
end:
    return token;
}

static tokenizer_t tokenizer_init(str_t str) {
    tokenizer_t tokenizer = {
        .str = str,
        .line = 1
    };
    return tokenizer;
}


// ####################
// ###   PRINTING   ###
// ####################

static size_t print_type_info(char* buf, size_t cap, type_info_t info) {
    size_t len = snprintf(buf, cap, "%s", get_value_type_str(info.base_type));

    for (size_t i = 0; i < MAX_NUM_DIMS; ++i) {
        if (info.dim[i] == -1) {
            len += snprintf(buf + len, cap - MIN(len, cap), "[..]");
        }
        else if (info.dim[i] >= 1) {
            len += snprintf(buf + len, cap - MIN(len, cap), "[%i]", (int)info.dim[i]);
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

#define PRINT(...) len += snprintf(buf + len, cap - MIN(len, cap), ##__VA_ARGS__)

static size_t print_bitfield(char* buf, size_t cap, const md_bitfield_t* bitfield) {
    ASSERT(bitfield);
    if (cap <= 0) return 0;

    const size_t max_bits = 64;
    const size_t num_bits = MIN(bitfield->end_bit - bitfield->beg_bit, max_bits);
    size_t len = 0;

    PRINT("<%i,%i>[", bitfield->beg_bit, bitfield->end_bit);
    for (size_t i = 0; i < num_bits; ++i) {
        PRINT("%i", md_bitfield_test_bit(bitfield, i) ? 1 : 0);
    }
    PRINT("%s]", num_bits == max_bits ? "..." : "");

    return len;
}

static size_t print_data_value(char* buf, size_t cap, data_t data) {
    size_t len = 0;
    if (is_variable_length(data.type)) return 0;
    if (cap <= 0) return 0;

    if (data.ptr) {
        if (is_array(data.type) && data.type.base_type != TYPE_BITFIELD) {
            int arr_len = data.type.dim[0];
            if (arr_len == 1) {
                dim_shift_left(data.type.dim);
                arr_len = data.type.dim[0];
            }
            type_info_t type = type_info_element_type(data.type);
            int64_t stride = (int64_t)type_info_element_byte_stride(data.type);

            PRINT("{");
            for (int64_t i = 0; i < arr_len; ++i) {
                data_t elem_data = {
                    .ptr = (char*)data.ptr + stride * i,
                    .size = data.size,
                    .type = type,
                };
                len += print_data_value(buf + len, cap - MIN(len, cap), elem_data);
                if (i < arr_len - 1) PRINT(",");
            }
            PRINT("}");
        } else {
            switch(data.type.base_type) {
            case TYPE_BITFIELD:
                //print_bitfield(file, (md_bitfield_t*)data.ptr);
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
                PRINT(STR_FMT, STR_ARG(str));
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
    if (!md_unit_empty(data.unit)) {
        len += md_unit_print(buf + len, cap - MIN(len, cap), data.unit);
    }
    return len;
}

#undef PRINT

#if MD_DEBUG

/*
static void print_data_value(FILE* file, data_t data) {
    if (is_variable_length(data.type)) return;

    if (is_array(data.type)) {
        int64_t len = data.type.dim[data.type.len_dim];
        if (len == 1) {
            data.type.dim[data.type.len_dim] = 0;
            data.type.len_dim = MAX(0, data.type.len_dim-1);
            len = data.type.dim[data.type.len_dim];
        }
        type_info_t type = type_info_element_type(data.type);
        int64_t stride = type_info_element_byte_stride(data.type);

        fprintf(file, "[");
        for (int64_t i = 0; i < len; ++i) {
            data_t new_data = {
                .ptr = (char*)data.ptr + stride * i,
                .size = data.size,
                .type = type,
                .unit = data.unit,
            };
            print_data_value(file, new_data);
            if (i < len - 1) fprintf(file, ", ");
        }
        fprintf(file, "]");
    } else {
        if (data.ptr) {
            switch(data.type.base_type) {
            case TYPE_BITFIELD:
                print_bitfield(file, (md_bitfield_t*)data.ptr);
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
                fprintf(file, ""STR_FMT"", (int)str.len, str.ptr);
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
        print_data_value(buf, sizeof(buf), node->data);
        fprintf(file, "%s", buf);
        break;
    case AST_PROC_CALL:
        if (node->proc) {
            fprintf(file, ""STR_FMT"", (int)node->proc->name.len, node->proc->name.ptr);
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
        fprintf(file, ""STR_FMT"", (int)node->ident.len, node->ident.ptr);
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
        fprintf(file, "logical op %s", get_token_type_str(node->marker.type));
        break;
    case AST_ADD:
    case AST_SUB:
    case AST_MUL:
    case AST_DIV:
        fprintf(file, "binary op %s", get_token_type_str(node->marker.type));
        break;
    default:
        fprintf(file, "%s", get_token_type_str(node->marker.type));
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
            print_data_value(buf, sizeof(buf), node->data);
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

    fprintf(file, "\""STR_FMT"\"\n", (int)at, buf);
}

static void save_expressions_to_json(expression_t** expr, uint64_t num_expr, str_t filename) {
    md_file_o* file = md_file_open(filename, MD_FILE_WRITE);

    if (file) {
        md_file_printf(file, "var treeData = [\n");
        for (uint64_t i = 0; i < num_expr; ++i) {
            md_file_printf(file, "{\n\"node\" :\n");
            print_node(file, expr[i]->node, 1);
            md_file_printf(file, "\"flags\" : %u,\n", expr[i]->node->flags);
            md_file_printf(file, "\"expr\" : ");
            print_expr(file, expr[i]->str);
            md_file_printf(file, "},\n");

        }
        md_file_printf(file, "];");
        md_file_close((md_file_o*)file);
    }
}

#endif

// ###################
// ###   PARSING   ###
// ###################

static void expand_node_token_range_with_children(ast_node_t* node) {
    ASSERT(node);
    for (size_t i = 0; i < md_array_size(node->children); ++i) {
        ASSERT(node->children[i] != node);
        if (!node->children[i]) continue;
        expand_node_token_range_with_children(node->children[i]);
        node->token = concat_tokens(node->token, node->children[i]->token);
        //node->token.col_beg = MIN(node->token.col_beg, node->children[i]->token.col_beg);
        //node->token.col_end = MAX(node->token.col_end, node->children[i]->token.col_end);
    }
}

static ast_node_t* parse_expression(parse_context_t* ctx);

static ast_node_t** parse_comma_separated_arguments_until_token(parse_context_t* ctx, token_type_t token_type) {
    ast_node_t** args = 0;
    token_t next = {0};
    while (next = tokenizer_peek_next(ctx->tokenizer), next.type != TOKEN_END) {
        if (next.type == token_type) goto done;
        if (next.type == ',') {
            LOG_ERROR(ctx->ir, next, "Empty argument in argument list");
            return NULL;
        }
        ctx->node = 0;
        ast_node_t* arg = parse_expression(ctx);
        next = tokenizer_peek_next(ctx->tokenizer);
        if (arg && (next.type == ',' || next.type == token_type)) {
            md_array_push(args, arg, ctx->temp_alloc);
            if (next.type == ',')
                tokenizer_consume_next(ctx->tokenizer);
            else // (next.type == token_type)
                goto done;
        } else {
            LOG_ERROR(ctx->ir, next, "Unexpected token in argument list");
            return NULL;
        }
    }
done:
    return args;
}

static ast_node_t* parse_procedure_call(parse_context_t* ctx, token_t token) {
    // Here we know atleast that the marker refers to something which is listed as a procedure!
    token_t next = tokenizer_peek_next(ctx->tokenizer);
    ast_node_t* node = 0;
    
    if (next.type == '(') {
        tokenizer_consume_next(ctx->tokenizer); // '('

        ast_node_t **args = parse_comma_separated_arguments_until_token(ctx, ')');
        next = tokenizer_consume_next(ctx->tokenizer);
        if (expect_token_type(ctx->ir, next, ')')) {
            node = create_node(ctx->ir, AST_PROC_CALL, token);
            node->ident = str_copy(token.str, ctx->ir->arena);
            const size_t num_args = md_array_size(args);
            if (num_args) {
                md_array_push_array(node->children, args, num_args, ctx->ir->arena);
            }
            // Expand proc call to contain entire argument list ')'
            node->token = concat_tokens(node->token, next);
        } else {
            LOG_ERROR(ctx->ir, token, "Unexpected end of argument list");
        }
    } else {
        node = create_node(ctx->ir, AST_PROC_CALL, token);
        node->ident = str_copy(token.str, ctx->ir->arena);
    }

    return node;
}

ast_node_t* parse_identifier(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == TOKEN_IDENT);

    if (ctx->node && (ctx->node->token.type == TOKEN_IDENT)) {
        LOG_ERROR(ctx->ir, concat_tokens(ctx->node->token, token), "Invalid syntax, identifiers must be separated by operation");
        return NULL;
    }
    
    const str_t ident = token.str;
    ast_node_t* node = 0;
    constant_t* c = 0;
    
    if (str_eq(ident, STR_LIT("import"))) {
        node = parse_procedure_call(ctx, token);
        if (node) {
            node->type = AST_TABLE;
        }
    } else if (str_eq(ident, STR_LIT("flatten"))) {
        node = parse_procedure_call(ctx, token);
        if (node) {
            node->type = AST_FLATTEN;
        }
    } else if (str_eq(ident, STR_LIT("transpose"))) {
        node = parse_procedure_call(ctx, token);
        if (node) {
            node->type = AST_TRANSPOSE;
        }
    } else if (is_identifier_procedure(ident)) {
        // The identifier has matched some procedure name
        node = parse_procedure_call(ctx, token);
    } else if ((c = find_constant(ident)) != 0) {
        node = create_node(ctx->ir, AST_CONSTANT_VALUE, token);
        node->value._float = c->value;
        node->data.unit = c->unit;
        node->data.ptr = &node->value._float;
        node->data.size = sizeof(float);
        node->data.type = (type_info_t)TI_FLOAT;
        node->flags |= FLAG_CONSTANT;
        node->ident = str_copy(ident, ctx->ir->arena);
    } else {
        // Identifier
        token_t next = tokenizer_peek_next(ctx->tokenizer);
        if (next.type == '(') {
            // This is an intended procedure call, but has no matching procedure by that name
            LOG_ERROR(ctx->ir, token, "Undefined procedure '"STR_FMT"'", ident.len, ident.ptr);

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
        } else {
            // Identifier reference, resolve this later in the static check
            node = create_node(ctx->ir, AST_IDENTIFIER, token);
            node->ident = str_copy(ident, ctx->ir->arena);
        }
    }
    return node;
}

ast_node_t* parse_logical(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == TOKEN_AND || token.type == TOKEN_OR || token.type == TOKEN_NOT || token.type == TOKEN_XOR);

    ast_type_t type = 0;

    switch (token.type) {
    case TOKEN_AND: type = AST_AND; break;
    case TOKEN_OR:  type = AST_OR;  break;
    case TOKEN_NOT: type = AST_NOT; break;
    case TOKEN_XOR: type = AST_XOR; break;
    default:
        ASSERT(false);
    }

    ast_node_t* lhs = ctx->node;
    ctx->node = 0;
    ast_node_t* rhs = parse_expression(ctx);

    if (type != AST_NOT && !lhs) {
        LOG_ERROR(ctx->ir, token, "Left hand side of '%s' did not evaluate to a valid expression.", get_token_type_str(token.type));
        return NULL;
    }

    if (!rhs) {
        LOG_ERROR(ctx->ir, token, "Right hand side of '%s' did not evaluate to a valid expression.", get_token_type_str(token.type));
        return NULL;
    }

    ast_node_t* node = create_node(ctx->ir, type, token);
    if (type != AST_NOT) md_array_push(node->children, lhs, ctx->ir->arena);
    md_array_push(node->children, rhs, ctx->ir->arena);
    
    return node;
}

ast_node_t* parse_assignment(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == '=');
    
    ast_node_t* lhs = ctx->node;
    if (!lhs) {
		LOG_ERROR(ctx->ir, token, "Left hand side of assignment (=) did not evaluate to a valid expression.");
		return NULL;
	}

    ctx->node = 0;
    ast_node_t* rhs = parse_expression(ctx);
    
    if (lhs->type != AST_IDENTIFIER) {
        if (lhs->type == AST_ARRAY) {
            for (size_t i = 0; i < md_array_size(lhs->children); ++i) {
				if (lhs->children[i]->type != AST_IDENTIFIER) {
					LOG_ERROR(ctx->ir, lhs->children[i]->token, "'%s' is not a valid identifier.", lhs->children[i]->token.str.ptr);
					return NULL;
				}
			}
        } else {
            LOG_ERROR(ctx->ir, token, "Left hand side of assignment (=) is not a valid identifier.");
            return NULL;
        }
    }
    if (!rhs) {
        LOG_ERROR(ctx->ir, token, "Right hand side of token '%s' did not evaluate to a valid expression.", get_token_type_str(token.type));
        return NULL;
    }

    if (rhs->type == AST_ASSIGNMENT) {
        LOG_ERROR(ctx->ir, token, "Syntax error: Invalid assignment");
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

    ast_node_t* lhs = ctx->node;
    ctx->node = 0;
    ast_node_t* rhs = parse_expression(ctx);
   
    if (!rhs) {
        LOG_ERROR(ctx->ir, token, "Invalid artihmetic expression, right operand is undefined.");
    } else {
        ast_type_t type = 0;
        switch(token.type) {
            case '+': type = AST_ADD; break;
            case '-': type = lhs ? AST_SUB : AST_UNARY_NEG; break;
            case '*': type = AST_MUL; break;
            case '/': type = AST_DIV; break;
            default: ASSERT(false);
        }
        ast_node_t* node = create_node(ctx->ir, type, token);
        if (type != AST_UNARY_NEG) {
            md_array_push(node->children, lhs, ctx->ir->arena);
        }
        md_array_push(node->children, rhs, ctx->ir->arena);
        return node;
    }
    
    return NULL;
}

ast_node_t* parse_value(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == TOKEN_FLOAT || token.type == TOKEN_INT || token.type == TOKEN_STRING || token.type == ':');
    token_t next = tokenizer_peek_next(ctx->tokenizer);
    
    ast_node_t* node = create_node(ctx->ir, AST_CONSTANT_VALUE, token);
    
    if (token.type == ':' || (next.type == ':' && (is_number(token.type)))) {
        // Expand the nodes marker to contain the range
        //node->token.col_end = next.col_end;
        node->token = concat_tokens(node->token, next);

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
            // Expand the nodes marker to contain the number
            //node->token.col_end = next.col_end;
            node->token = concat_tokens(node->token, next);

        } else {
            end = FLT_MAX;
        }

        if (is_float) {
            node->data.type = (type_info_t)TI_FRANGE;
            node->value._frange.beg = (float)beg;
            node->value._frange.end = (float)end;
        } else {
            node->data.type = (type_info_t)TI_IRANGE;
            node->value._irange.beg = (beg == -FLT_MAX) ?  INT32_MIN : (int32_t)beg;
            node->value._irange.end = (end ==  FLT_MAX) ?  INT32_MAX : (int32_t)end;
        }
    }
    else if (token.type == TOKEN_FLOAT) {
        node->data.type = (type_info_t)TI_FLOAT;
        node->value._float = (float)parse_float(token.str);
    }
    else if (token.type == TOKEN_INT) {
        node->data.type = (type_info_t)TI_INT;
        node->value._int = (int32_t)parse_int(token.str);
    }
    else if (token.type == TOKEN_STRING) {
        node->data.type = (type_info_t)TI_STRING;
        node->value._string = str_copy(str_substr(token.str, 1, MAX(0, ((int)token.str.len - 2))), ctx->ir->arena); // remove quotation marks
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
            const size_t num_elements = md_array_size(elements);
            if (num_elements) {

                ast_node_t* node_copy = create_node(ctx->ir, node->type, token);
                MEMCPY(node_copy, node, sizeof(ast_node_t));

                MEMSET(node, 0, sizeof(ast_node_t));
                node->type = AST_ARRAY_SUBSCRIPT;
                node->token = token;

                md_array_push(node->children, node_copy, ctx->ir->arena);
                md_array_push_array(node->children, elements, num_elements, ctx->ir->arena);

                // Expand marker range to contain ']'
                //node->token.col_end = next.col_end;
                node->token = concat_tokens(node->token, next);

                return node;
            } else {
                LOG_ERROR(ctx->ir, token, "Empty array subscripts are not allowed.");
            }
        } else {
            LOG_ERROR(ctx->ir, token, "Unexpected end of argument list");
        }
    } else {
        LOG_ERROR(ctx->ir, token, "Missing left hand side expression to apply subscript operator");
    }
    return NULL;
}

ast_node_t* parse_array(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == '{');

    ast_node_t** elements = parse_comma_separated_arguments_until_token(ctx, '}');
    token_t next = tokenizer_consume_next(ctx->tokenizer);
    if (expect_token_type(ctx->ir, next, '}')) {
        const size_t num_elements = md_array_size(elements);
        if (num_elements) {
            // We only commit the results if everything at least parsed ok.
            ast_node_t* node = create_node(ctx->ir, AST_ARRAY, token);
            md_array_push_array(node->children, elements, num_elements, ctx->ir->arena);
            // Expand marker range with '}'
            //node->token.col_end = next.col_end;
            node->token = concat_tokens(node->token, next);
            return node;
        } else {
            LOG_ERROR(ctx->ir, token, "Empty arrays are not allowed.");
        }
    }
    else {
        LOG_ERROR(ctx->ir, token, "Unexpected end of argument list");
    }
    return NULL;
}

ast_node_t* parse_in(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == TOKEN_IN);
    
    ast_node_t* lhs = ctx->node;
    ctx->node = 0;
    ast_node_t* rhs = parse_expression(ctx);

    if (!lhs) {
        LOG_ERROR(ctx->ir, token, "Left hand side of 'in' did not evaluate into a valid expression.");
    }
    else if(!rhs) {
        LOG_ERROR(ctx->ir, token, "Right hand side of 'in' did not evaluate into a valid expression.");
    }
    else {
        ast_node_t* node = create_node(ctx->ir, AST_CONTEXT, token);
        md_array_push(node->children, lhs, ctx->ir->arena);
        md_array_push(node->children, rhs, ctx->ir->arena);
        return node;
    }
    return NULL;
}

ast_node_t* parse_out(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == TOKEN_OUT);

    ctx->node = 0;
    ast_node_t* rhs = parse_expression(ctx);

    if(!rhs) {
        LOG_ERROR(ctx->ir, token, "Right hand side of 'out' did not evaluate into a valid expression.");
    }
    else {
        ast_node_t* node = create_node(ctx->ir, AST_OUT, token);
        md_array_push(node->children, rhs, ctx->ir->arena);
        return node;
    }
    return NULL;
}

ast_node_t* parse_parentheses(parse_context_t* ctx) {
    token_t token = tokenizer_consume_next(ctx->tokenizer);
    ASSERT(token.type == '(');
    ast_node_t* expr = parse_expression(ctx);
    if (!expr) {
        LOG_ERROR(ctx->ir, token, "Expression inside parentheses did not evaluate into a valid expression.");
        return NULL;
    }
    token_t next = tokenizer_consume_next(ctx->tokenizer);
    if (expect_token_type(ctx->ir, next, ')')) {
        ast_node_t* node = create_node(ctx->ir, AST_EXPRESSION, token);
        md_array_push(node->children, expr, ctx->ir->arena);
        // Expand token marker
        node->token = concat_tokens(node->token, next);
        return node;
    }
    return NULL;
}

/*
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
        LOG_ERROR(ctx->ir, next, "Undefined marker!");
        tokenizer_consume_next(ctx->tokenizer);
        return NULL;
    default:
        LOG_ERROR(ctx->ir, next, "Unexpected marker value! '"STR_FMT"'",
            next.str.len, next.str.ptr);
        return NULL;
    }
    return ctx->node;
}
*/

ast_node_t* parse_expression(parse_context_t* ctx) {
    token_t token = {0};
    while (token = tokenizer_peek_next(ctx->tokenizer), token.type != TOKEN_END) {
        switch (token.type) {
            case TOKEN_AND:
            case TOKEN_OR:
            case TOKEN_XOR:
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
            case TOKEN_OUT:
            ctx->node = parse_out(ctx);
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
            LOG_ERROR(ctx->ir, token, "Invalid token!");
            tokenizer_consume_next(ctx->tokenizer);
            return NULL;
            default:
            LOG_ERROR(ctx->ir, token, "Unexpected token: '"STR_FMT"'", token.str.len, token.str.ptr);
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

static bool finalize_type(type_info_t* type, const ast_node_t* node, eval_context_t* ctx);
static bool evaluate_node(data_t*, const ast_node_t*, eval_context_t*);

static int do_proc_call(data_t* dst, const procedure_t* proc,  ast_node_t** const args, size_t num_args, eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(proc);
    ASSERT(num_args < MAX_SUPPORTED_PROC_ARGS);

    int result = 0;

    const bool visualize_only = ctx->vis && !dst;
    if (visualize_only && !(proc->flags & FLAG_VISUALIZE)) {
        // If the node has not been flagged with visualize, it should not be called during visualization
        // Unless there is some data (dst) which should be filled in.
        // Do however propagate the evaluation to the children
        for (size_t i = 0; i < num_args; ++i) {
            if (!evaluate_node(NULL, args[i], ctx)) {
                return -1;
            }
        }
        return result;
    }

    data_t  arg_data  [MAX_SUPPORTED_PROC_ARGS] = {0};
    token_t arg_tokens[MAX_SUPPORTED_PROC_ARGS] = {0};
    flags_t arg_flags [MAX_SUPPORTED_PROC_ARGS] = {0};

    //uint64_t stack_reset_point = md_stack_allocator_get_pos(ctx->stack_alloc);

    md_script_vis_t* old_vis = ctx->vis;
    // If we are in visualization mode and the procedure is flagged with visualization, we do not want to propagate the visualization to the arguments
    if (ctx->vis && proc->flags & FLAG_VISUALIZE) {
        ctx->vis = 0;
    }

    for (size_t i = 0; i < num_args; ++i) {
        // We need to evaluate the argument nodes first before we make the proc call.
        // In this context we are not interested in storing any data (since it is only used to be passed as arguments)
        // so we can allocate the data required for the node with the temp alloc

        arg_tokens[i] = args[i]->token;
        arg_flags[i] = args[i]->flags;

        if (args[i]->flags & FLAG_CONSTANT) {
            ASSERT(args[i]->data.ptr);
            MEMCPY(&arg_data[i], &args[i]->data, sizeof(data_t));
        }
        else {
            type_info_t arg_type = args[i]->data.type;

            if (is_variable_length(arg_type)) {
                if (!finalize_type(&arg_type, args[i], ctx)) {
                    MD_LOG_ERROR("Failed to finalize dynamic type in procedure call");
                    result = -1;
                    goto done;
                }
            }
            allocate_data(&arg_data[i], arg_type, ctx->temp_alloc);
            if (!evaluate_node(&arg_data[i], args[i], ctx)) {
                result = -1;
                goto done;
            }
        }
    }

    ctx->vis = old_vis;

    flags_t* old_arg_flags  = ctx->arg_flags;
    token_t* old_arg_tokens = ctx->arg_tokens;
    // Set
    ctx->arg_tokens = arg_tokens;
    ctx->arg_flags  = arg_flags;

    result = proc->proc_ptr(dst, arg_data, ctx);

    // Reset
    ctx->arg_flags  = old_arg_flags;
    ctx->arg_tokens = old_arg_tokens;

done:
    for (int64_t i = (int64_t)num_args - 1; i >= 0; --i) {
        if (!(args[i]->flags & FLAG_CONSTANT)) {
            free_data(&arg_data[i], ctx->temp_alloc);
        }
    }

    // @NOTE(Robin): We cannot simply reset the stack since bitfields are lazily allocated, meaning they allocate data on demand.
    // This means that they will most likely allocate data deeply nested within the calling scope of another procedure call with another reset point for its stack
    // Thus, the data for the bitfields will be reset with that reset point and we will overwrite it later.
    // To solve it I think we need another 'specific' allocator for just the bitfields which may outlive the other stack variables and can be reset separately.
    // Conceptually we could use the top end of the stack for this and have two different stack pointers.
    // 
    //md_stack_allocator_set_pos(ctx->stack_alloc, stack_reset_point);
    return result;
}

static bool evaluate_proc_call(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    //ASSERT(node->type == AST_PROC_CALL);
    ASSERT(node->proc);

    token_t old_token = ctx->op_token;
    ctx->op_token = node->token;
    int result = do_proc_call(dst, node->proc, node->children, md_array_size(node->children), ctx);
    ctx->op_token = old_token;

    return result >= 0;
}

static bool evaluate_constant_value(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_CONSTANT_VALUE);
    (void)ctx;

    if (dst) {
        ASSERT(node->data.ptr && node->data.size > 0);  // Assume that the static check set the data already
        ASSERT(dst->ptr && dst->size >= node->data.size); // Make sure we have the expected size
        copy_data(dst, &node->data);
    }

    return true;
}

static bool evaluate_table_lookup(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_TABLE);

    if (dst) {
        ASSERT(dst->ptr && dst->size >= type_info_total_byte_size(node->data.type));

        const size_t num_frames = md_trajectory_num_frames(ctx->traj);

        if (!num_frames || !ctx->frame_header) {
            MD_LOG_ERROR("Missing trajectory header data, cannot evaluate table");
            return false;
        }
        if (!is_type_equivalent(dst->type, node->data.type)) {
            MD_LOG_ERROR("Type mismatch when evaluating table value");
            return false;
        }
        if (!node->table || !node->table_field_indices) {
            MD_LOG_ERROR("Missing table and or field mappings when evaluating table value");
            return false;
        }

        const int64_t num_rows = (int64_t)node->table->num_values;
        int64_t row_index = -1;
        if (str_eq(node->table->field_names[0], STR_LIT("Time"))) {
            // Complex case, find the matching time
            const double* time = node->table->field_values[0];

            if (!time) {
                MD_LOG_DEBUG("Missing time field when evaluating table value");
                return false;
            }

            const double ref_time = ctx->frame_header->timestamp;

            // Find the row_index in the table which corresponds to the current time
            double t0 = 0.0;
            double t1 = 1.0;
            const double* traj_times = md_trajectory_frame_times(ctx->traj);
            if (traj_times) {
                t0 = traj_times[0];
                t1 = traj_times[num_frames - 1];
            }

            // Make a guess based on the current time
            const double t = (ref_time - t0) / (t1 - t0);

            // Guess row index from t
            row_index = CLAMP((int64_t)(t * num_rows), 0, num_rows - 1);

            if (time[row_index] < ref_time) {
                for (int64_t i = row_index; i < num_rows; ++i) {
                    if (timestamps_approx_equal(time[i], ref_time)) {
                        row_index = i;
                        break;
                    }
                    else if (time[i] - TIMESTAMP_ERROR_MARGIN > ref_time) {
                        row_index = -1;
                    }
                }
            }
            else if (time[row_index] > ref_time) {
                for (int64_t i = row_index; i >= 0; --i) {
                    if (timestamps_approx_equal(time[i], ref_time)) {
                        row_index = i;
                        break;
                    }
                    else if (time[i] + TIMESTAMP_ERROR_MARGIN < ref_time) {
                        row_index = -1;
                    }
                }
            }
        } else {
            // Simple case, find the matching frame
			row_index = ctx->frame_header->index;
        }

        if (row_index == -1) {
            MD_LOG_DEBUG("Unable to lookup corresponding row when evaluating table");
            return false;
        }

        if (row_index >= num_rows) {
			MD_LOG_DEBUG("Row index out of bounds when evaluating table");
			return false;
		}
         
        for (size_t i = 0; i < md_array_size(node->table_field_indices); ++i) {
            int idx = node->table_field_indices[i];
            ((float*)dst->ptr)[i] = (float)node->table->field_values[idx][row_index];
        }
    }

    return true;
}

static bool evaluate_flatten(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_FLATTEN);
    ASSERT(node->children && md_array_size(node->children) == 1);
    const ast_node_t* arg = node->children[0];

    if (dst) {
        if (node->data.type.base_type == TYPE_BITFIELD) {
            data_t data = {0};
            const md_bitfield_t* src_arr = 0;
            size_t src_len = 0;

            if (arg->flags & FLAG_CONSTANT && arg->data.ptr) {
                // Assume that it has data. No need to copy here...
                src_arr = arg->data.ptr;
                src_len = type_info_array_len(arg->data.type);
            } else {
                if (!allocate_data(&data, arg->data.type, ctx->temp_alloc)) {
                    LOG_ERROR(ctx->ir, node->token, "Failed to allocate data to evaluate flatten operation");
                    return false;
                }
                if (!evaluate_node(&data, arg, ctx)) {
                    return false;
                }
                src_arr = data.ptr;
                src_len = type_info_array_len(data.type);
            }

            md_bitfield_t* dst_bf = dst->ptr;
            for (size_t i = 0; i < src_len; ++i) {
                md_bitfield_or_inplace(dst_bf, src_arr + i);
            }
            return true;
        } else {
            data_t new_dst = *dst;
            new_dst.type = arg->data.type;
            return evaluate_node(&new_dst, arg, ctx);
        }
    } else {
        return evaluate_node(dst, arg, ctx);
    }
}

static bool evaluate_transpose(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_TRANSPOSE);
    ASSERT(node->children && md_array_size(node->children) == 1);
    const ast_node_t* arg = node->children[0];
    int ndim = dim_ndims(arg->data.type.dim);

    if (dst && ndim == 2) {
        if (ndim != dim_ndims(dst->type.dim)) {
            LOG_ERROR(ctx->ir, node->token, "Incorrect dimensions in transpose operation");
            return false;
        }
        ASSERT(dst->type.base_type == arg->data.type.base_type);
        ASSERT(dst->type.dim[0] == arg->data.type.dim[1]);
        ASSERT(dst->type.dim[1] == arg->data.type.dim[0]);

        data_t data = arg->data;
        data.ptr = 0;
        data.size = 0;
        /*
        // @NOTE: For now, we just support static length arrays as the transposition throws a bone into the machinery when it comes to determining length
        // That is, it shifts the dimension of the undetermined length
        if (arg->data.type.dim[0] == -1) {
            if (!finalize_type(&data.type, arg, ctx)) {
                LOG_ERROR(ctx->ir, node->token, "Failed to finalize type of subexpression in transpose operation");
                return false;
            }
            dst->type.dim[1] = arg->data.type.dim[0];
        }
        */
        allocate_data(&data, data.type, ctx->temp_alloc);
        if (!evaluate_node(&data, arg, ctx)) {
            LOG_ERROR(ctx->ir, node->token, "Failed to evaluate subexpression in transpose operation");
            return false;
        }
        const size_t num_col = (size_t)dst->type.dim[0];
        const size_t num_row = (size_t)dst->type.dim[1];
        const size_t num_bytes = base_type_element_byte_size(dst->type.base_type);
        char* dst_base = dst->ptr;
        const char* src_base = data.ptr;
        for (int col = 0; col < num_col; ++col) {
            for (int row = 0; row < num_row; ++row) {
                MEMCPY(dst_base + (col * num_row + row) * num_bytes, src_base + (row * num_col + col) * num_bytes, num_bytes);
            }
        }
        return true;
    } else {
        return evaluate_node(dst, arg, ctx);
    }
}

static identifier_t* find_static_identifier(str_t name, eval_context_t* ctx) {
    ASSERT(ctx);

    // Static Identifiers from Compilation
    if (ctx->ir) {
        for (size_t i = 0; i < md_array_size(ctx->ir->identifiers); ++i) {
            // Only return IR identifiers if they are constant
            if (ctx->ir->identifiers[i].node && (ctx->ir->identifiers[i].node->flags & FLAG_CONSTANT) && str_eq(name, ctx->ir->identifiers[i].name)) {
                return &ctx->ir->identifiers[i];
            }
        }
    }

    return NULL;
}

// This should only called in an evaluation context
static identifier_t* find_dynamic_identifier(str_t name, eval_context_t* ctx) {
    ASSERT(ctx);

    // Dynamic Identifiers from Evaluation
    for (size_t i = 0; i < md_array_size(ctx->identifiers); ++i) {
        if (str_eq(name, ctx->identifiers[i].name)) {
            return &ctx->identifiers[i];
        }
    }

    return NULL;
}

static bool evaluate_identifier_reference(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_IDENTIFIER);
    ASSERT(md_array_size(node->children) == 1 && node->children[0]);
    (void)ctx;
    
    // @TODO(Robin): Check if we are properly handling the case where we have already evaluated the identifier
    // And in such case, that should be copied as well.
    
    // Only copy data if the node is constant
    if (dst && node->flags & FLAG_CONSTANT) {
        ASSERT(dst->ptr && node->data.ptr && dst->size >= node->data.size);
        copy_data(dst, &node->data);
    } else if (dst || ctx->vis) {
        evaluate_node(dst, node->children[0], ctx);
    }

    return true;
}

static bool evaluate_assignment(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ASSIGNMENT);
    ASSERT(node->children);

    ASSERT(md_array_size(node->children) == 2);
    ast_node_t* lhs = node->children[0];
    ast_node_t* rhs = node->children[1];

    ASSERT(lhs && (lhs->type == AST_IDENTIFIER || lhs->type == AST_ARRAY));
    ASSERT(rhs);

    // We know from the static check that the lhs (child[0]) is one or more identifiers and the rhs is a valid expression with a known data size
    // the data size could be zero if it is an array of length 0, in that case we don't need to allocate anything.
    // This should also be the conceptual entry point for all evaluations in this declarative language
    if (lhs->type == AST_IDENTIFIER) {
        if (!dst && ctx->vis) {
            return evaluate_node(NULL, rhs, ctx);
        }

        identifier_t* ident = find_static_identifier(lhs->ident, ctx);
        if (ident) {
            ASSERT(ident->node->flags & FLAG_CONSTANT);
            if (dst) {
                copy_data(dst, ident->data);
            }
            if (ctx->vis) {
                return evaluate_node(NULL, rhs, ctx);
            }
            return true;
        }

        ident = find_dynamic_identifier(lhs->ident, ctx);
        if (!ident && ctx->alloc) {
            identifier_t id = {
                .name = lhs->ident,
                .node = rhs,
                .data = dst,
            };
            ident = md_array_push(ctx->identifiers, id, ctx->alloc);
            return evaluate_node(dst, rhs, ctx);
        }
    } else if (lhs->type == AST_ARRAY) {
        // assignment of LHS which is an array of Identifiers
        // e.g {x,y,z} = coord(1);
        int count = (int)md_array_size(lhs->children);
        if (rhs->flags & FLAG_CONSTANT) {
            if (dst) {
                copy_data(dst, &rhs->data);
            }
            if (ctx->vis) {
                return evaluate_node(NULL, rhs, ctx);
            }
            return true;
        }

        if (!evaluate_node(dst, rhs, ctx)) {
            return false;
        }

        for (int i = 0; i < count; ++i) {
            ast_node_t* child = lhs->children[i];
            int stride = (int)type_info_element_byte_stride(child->data.type);
            ASSERT(child->type == AST_IDENTIFIER);
            identifier_t* ident = find_dynamic_identifier(child->ident, ctx);

            if (!ident && ctx->alloc) {
                identifier_t id = {
                    .name = child->ident,
                    .node = rhs,
                    .data = md_alloc(ctx->alloc, sizeof(data_t))
                };
                if (dst && dst->ptr) {
                    *id.data = child->data;
                    id.data->ptr = (char*)dst->ptr + stride * i;
                    id.data->size = stride;
                }
                ident = md_array_push(ctx->identifiers, id, ctx->alloc);
            }
        }
        return true;
    }

    return false;
}

static bool index_in_range(int idx, irange_t range) {
    return range.beg <= idx && idx < range.end;
}

static md_array(irange_t) extract_subscript_ranges(const data_t data, md_allocator_i* alloc) {
    md_array(irange_t) result = 0;

    if (type_info_equal(data.type, (type_info_t)TI_INT)) {
        int len = type_info_array_len(data.type);
        for (int i = 0; i < len; ++i) {
            int value = *(int*)data.ptr;
            irange_t range = {value, value};
            md_array_push(result, range, alloc);
        }
    } else if (type_info_equal(data.type, (type_info_t)TI_IRANGE)) {
        int len = type_info_array_len(data.type);
        irange_t* ranges = (irange_t*)(data.ptr);
        md_array_push_array(result, ranges, len, alloc);
    } else {
        ASSERT(false);
    }

    return result;
}

static bool evaluate_array(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY);

    // Even if we do not have a dst parameter we still want to propagate the evaluation to the children
    const size_t num_args = md_array_size(node->children);
    ast_node_t** const args = node->children;
    bool result = false;

    md_array(irange_t) ranges = ctx->subscript_ranges;
    ctx->subscript_ranges = 0;

    if (dst) {
        ASSERT(type_info_equal(type_info_element_type(dst->type), type_info_element_type(node->data.type)));

        size_t byte_offset = 0;
        for (size_t i = 0; i < num_args; ++i) {
            type_info_t type = args[i]->data.type;
            if (is_variable_length(type)) {
                if (!finalize_type(&type, args[i], ctx)) {
                    md_logf(MD_LOG_TYPE_DEBUG, "evaluate_array: Failed to finalize type for variable length argument...");
                    goto done;
                }
            }
            data_t data = {
                .type = type,
                .ptr = (char*)dst->ptr + byte_offset,
                .size = type_info_total_byte_size(type)
            };
            md_script_vis_t* vis = ctx->vis;
            if (ranges && !index_in_range((int)i, ranges[0])) {
                ctx->vis = 0;
            }
            if (!evaluate_node(&data, args[i], ctx)) {
                goto done;
            }
            ctx->vis = vis;
            byte_offset += data.size;
        }
        ASSERT(byte_offset == dst->size);
    }
    else {

        for (size_t i = 0; i < num_args; ++i) {
            md_script_vis_t* vis = ctx->vis;
            if (ranges && !index_in_range((int)i, ranges[0])) {
                ctx->vis = 0;
            }
            if (!evaluate_node(NULL, args[i], ctx)) {
                goto done;
            }
            ctx->vis = vis;
        }
    }

    result = true;
done:
    ctx->subscript_ranges = ranges;
    return result;
}

static bool evaluate_array_subscript(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY_SUBSCRIPT);
    ASSERT(node->children);

    const ast_node_t* arr = node->children[0];
    data_t arr_data = {0};

    md_array(irange_t) idx_ranges = md_array_create(irange_t, node->subscript_dim, ctx->temp_alloc);
    md_array(irange_t) ext_ranges = ctx->subscript_ranges;

    // @NOTE: This is a bit of a hack, we should not need to allocate the subscript ranges here
    // It is just to serve the poor design of the current propagation of array subscripts in the context
    MEMCPY(idx_ranges, node->subscript_ranges, sizeof(irange_t) * node->subscript_dim);

    if (ext_ranges) {
        for (int i = 0; i < node->subscript_dim; ++i) {
            if (ext_ranges[i].beg == 0 && ext_ranges[i].end == 0) {
                idx_ranges[i].beg = node->subscript_ranges[i].beg;
                idx_ranges[i].end = node->subscript_ranges[i].end;
            } else {
                idx_ranges[i].beg = node->subscript_ranges[i].beg + ext_ranges[i].beg;
                idx_ranges[i].end = node->subscript_ranges[i].beg + ext_ranges[i].end;
            }

        }
    } else {
        MEMCPY(idx_ranges, node->subscript_ranges, sizeof(irange_t) * node->subscript_dim);
    }
    ctx->subscript_ranges = idx_ranges;

    if (dst) {
        if (!allocate_data(&arr_data, arr->data.type, ctx->temp_alloc)) return false;
        if (!evaluate_node(&arr_data, arr, ctx)) return false;
    } else {
        if (!evaluate_node(0, arr, ctx)) return false;
    }
    ctx->subscript_ranges = ext_ranges;

    if (dst) {
        ASSERT(arr_data.ptr);
        // Extract data from the underlying type given the subscript ranges
        ASSERT(type_info_equal(dst->type, node->data.type));
        ASSERT(dst->ptr);

        const size_t elem_size = base_type_element_byte_size(arr_data.type.base_type);
        size_t write_offset = 0;
        const irange_t* src_rng = node->subscript_ranges;
        const int* src_dim = arr_data.type.dim;

        switch (node->subscript_dim) {
            case 1: {
                const size_t slice_size = (src_rng[0].end - src_rng[0].beg) * elem_size;
                const size_t read_offset = src_rng[0].beg * elem_size;
                data_t dst_slice = {dst->type, (char*)dst->ptr + write_offset, slice_size};
                const data_t src_slice = {dst->type, (char*)arr_data.ptr + read_offset, slice_size};
                copy_data(&dst_slice, &src_slice);
            } break;
            case 2:
                for (int i = src_rng[0].beg; i < src_rng[0].end; ++i) {
                    const size_t slice_size = (src_rng[1].end - src_rng[1].beg) * elem_size;
                    const size_t read_offset = (i * src_dim[1] + src_rng[1].beg) * elem_size;
                    data_t dst_slice = {dst->type, (char*)dst->ptr + write_offset, slice_size};
                    const data_t src_slice = {dst->type, (char*)arr_data.ptr + read_offset, slice_size};
                    copy_data(&dst_slice, &src_slice);
                    write_offset += slice_size;
                }
                break;
            case 3:
                // @NOTE: TO BE TESTED
                for (int j = src_rng[0].beg; j < src_rng[0].end; ++j) {
                    for (int i = src_rng[1].beg; i < src_rng[1].end; ++i) {
                        const size_t slice_size = (src_rng[2].end - src_rng[2].beg) * elem_size;
                        const size_t read_offset = (j * src_dim[1] * src_dim[2] + i * src_dim[2] + src_rng[2].beg) * elem_size;
                        data_t dst_slice = {dst->type, (char*)dst->ptr + write_offset, slice_size};
                        const data_t src_slice = {dst->type, (char*)arr_data.ptr + read_offset, slice_size};
                        copy_data(&dst_slice, &src_slice);
                        write_offset += slice_size;
                    }
                }
                break;
            case 4:
                // @NOTE: TO BE TESTED
                for (int k = src_rng[0].beg; k < src_rng[0].end; ++k) {
                    for (int j = src_rng[1].beg; j < src_rng[1].end; ++j) {
                        for (int i = src_rng[2].beg; i < src_rng[2].end; ++i) {
                            const size_t slice_size = (src_rng[3].end - src_rng[3].beg) * elem_size;
                            const size_t read_offset = (k * src_dim[1] * src_dim[2] * src_dim[3] + j * src_dim[2] * src_dim[3] + i * src_dim[3] + src_rng[3].beg) * elem_size;
                            data_t dst_slice = {dst->type, (char*)dst->ptr + write_offset, slice_size};
                            const data_t src_slice = {dst->type, (char*)arr_data.ptr + read_offset, slice_size};
                            copy_data(&dst_slice, &src_slice);
                            write_offset += slice_size;
                        }
                    }
                }
                break;
            default:
                ASSERT(false);
        }
    }

    return true;
}

static bool evaluate_out(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(node && node->type == AST_OUT);
    ASSERT(md_array_size(node->children) == 1);

    const ast_node_t* expr = node->children[0];
    const md_bitfield_t* old_ctx = ctx->mol_ctx;
    ctx->mol_ctx = 0;
    bool result = evaluate_node(dst, expr, ctx);
    ctx->mol_ctx = old_ctx;

    return result;
}

static bool evaluate_context(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(node && node->type == AST_CONTEXT);

    if (dst && node->flags & FLAG_CONSTANT) {
        copy_data(dst, &node->data);
        return true;
    }

    ASSERT(md_array_size(node->children) == 2);

    const ast_node_t* expr_node = node->children[0];
    const ast_node_t* ctx_node  = node->children[1];

    // No need to reevaluate CTX_NODE here since it has already been done during static type checking and its results should be stored in the data ptr

    ASSERT(ctx_node->data.type.base_type == TYPE_BITFIELD);
    ASSERT(ctx_node->data.ptr);

    const int num_ctx = type_info_array_len(ctx_node->data.type);
    if (num_ctx < 0) {
        LOG_ERROR(ctx->ir, node->token, "Invalid number of contexts (%i) in context expression", num_ctx);
		return false;
    }
    const md_bitfield_t* ctx_bf = (const md_bitfield_t*)ctx_node->data.ptr;

    ASSERT(node->lhs_context_types);
    const type_info_t* lhs_types = node->lhs_context_types;
    ASSERT((int)md_array_size(lhs_types) == num_ctx);

    data_t data = {0};
    data.unit = expr_node->data.unit;
    data.value_range = expr_node->data.value_range;
   
    data_t* sub_dst = dst ? &data : NULL;

    bool result = true;

    // Evaluate the expression within each context.
    int64_t dst_idx = 0;
    for (int64_t i = 0; i < num_ctx; ++i) {
        eval_context_t sub_ctx = *ctx;
        sub_ctx.mol_ctx = &ctx_bf[i];
        sub_ctx.subscript_ranges = 0;

        if (dst) {
            const size_t elem_size = type_info_element_byte_stride(dst->type);
            data.type = lhs_types[i];
            data.ptr = (char*)dst->ptr + elem_size * dst_idx;
            data.size = type_info_total_byte_size(lhs_types[i]);;

            int64_t offset = type_info_array_len(lhs_types[i]);
            dst_idx += offset;
        }

        if (ctx->subscript_ranges) {
            if (!index_in_range((int)i, ctx->subscript_ranges[0])) {
                // Clear the visualization context if this index is not present in subscript ranges to prevent visualization
                sub_ctx.vis = 0;
            }
        }

        md_bitfield_t* prev_vis_structure = 0;
        if (sub_ctx.vis) {
            prev_vis_structure = sub_ctx.vis_structure;
            if (sub_ctx.vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
                sub_ctx.vis_structure = push_structure(sub_ctx.vis);
            }
        }

        if (!evaluate_node(sub_dst, expr_node, &sub_ctx)) {
            result = false;
        }

        if (sub_ctx.vis) {
            // Reset atom mask
            if (sub_ctx.vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
                sub_ctx.vis_structure = prev_vis_structure;
            }
        }

        sub_ctx.vis = ctx->vis;

        if (!result) break;
    }

    return result;
}

static bool evaluate_node(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);

    if (dst && node->flags & FLAG_CONSTANT) {
        copy_data(dst, &node->data);
        return true;
    }

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
        case AST_TABLE:
            return evaluate_table_lookup(dst, node, ctx);
        case AST_FLATTEN:
            return evaluate_flatten(dst, node, ctx);
        case AST_TRANSPOSE:
            return evaluate_transpose(dst, node, ctx);
        case AST_IDENTIFIER:
            return evaluate_identifier_reference(dst, node, ctx);
        case AST_CONTEXT:
            return evaluate_context(dst, node, ctx);
        case AST_OUT:
            return evaluate_out(dst, node, ctx);
        case AST_NOT:       // All operators should have been resolved during static check and converted into proc-calls
        case AST_UNARY_NEG:
        case AST_MUL:
        case AST_DIV:
        case AST_ADD:
        case AST_SUB:
        case AST_AND:
        case AST_XOR:
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
static bool finalize_proc_call(ast_node_t*, eval_context_t*);
static bool static_eval_node(ast_node_t*, eval_context_t*);

static bool convert_node(ast_node_t* node, type_info_t new_type, eval_context_t* ctx) {
    const type_info_t from = node->data.type;
    type_info_t to   = new_type;

    procedure_match_result_t res = find_cast_procedure(from, to);
    if (res.success) {
        // We need to update the return type here
        to = res.procedure->return_type;

        if (node->flags & FLAG_CONSTANT) {
            if (to.dim[0] == -1) {
                if (res.procedure->flags & FLAG_RET_AND_FIRST_ARG_EQUAL_LENGTH) {
                    to.dim[0] = from.dim[0];
                } else {
                    ASSERT(res.procedure->flags & FLAG_QUERYABLE_LENGTH);

                    token_t  old_op_token = ctx->op_token;
                    ctx->op_token = node->token;
                    // Perform the call to get length
                    int query_result = do_proc_call(NULL, res.procedure, &node, 1, ctx);
                    ctx->op_token = old_op_token;

                    if (query_result >= 0) { // Zero length is valid
                        to.dim[0] = query_result;
                    }
                    else {
                        LOG_ERROR(ctx->ir, node->token, "Unexpected return value (%i) when querying procedure for array length.", query_result);
                        return false;
                    }
                }
            }

            ASSERT(to.dim[0] > -1);

            data_t new_data = {0};
            if (node->type == AST_CONSTANT_VALUE) {
                new_data = node->data;
                ASSERT(is_scalar(to));
                new_data.type = to;
                new_data.size = type_info_total_byte_size(to);
                new_data.ptr  = &node->value;
            } else {
                allocate_data(&new_data, to, ctx->alloc);        
            }

            // Perform the data conversion
            if (do_proc_call(&new_data, res.procedure, &node, 1, ctx) < 0) {
                LOG_ERROR(ctx->ir, node->token, "Failed to perform data conversion of static data");
                return false;
            }
            // Finalize the node with the new converted data
            node->data = new_data;
            return true;
        } else {
            // We need to convert this node into a cast node and add the original node data as a child
            ast_node_t* node_copy = create_node(ctx->ir, node->type, node->token);
            MEMCPY(node_copy, node, sizeof(ast_node_t));
            node->data = (data_t){0};
            node->data.unit = node_copy->data.unit;
            node->type = AST_PROC_CALL;
            node->children = 0; // node_copy have taken over the children, we need to zero this to trigger a proper allocation in next step
            md_array_push(node->children, node_copy, ctx->ir->arena);
            node->proc = res.procedure;
            node->proc_flags = res.flags;

            return finalize_proc_call(node, ctx);
        }
    }

    ASSERT(false);
    return false;
}

static bool finalize_type_proc(type_info_t* type, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(type);
    ASSERT(node);
    ASSERT(ctx);

    *type = node->data.type;

    ast_node_t** const args = node->children;
    const size_t num_args = md_array_size(node->children);

    if (node->proc->flags & FLAG_RET_AND_FIRST_ARG_EQUAL_LENGTH) {
        if (!args || num_args == 0) {
            LOG_ERROR(ctx->ir, node->token, "Unexpected number of arguments (0), but procedure is marked to extract length from first argument");
            return false;
        }
        // We can deduce the length of the array from the input type (first argument)
        type->dim[0] = type_info_array_len(args[0]->data.type);
        if (type->dim[0] == -1) {
            LOG_ERROR(ctx->ir, node->token, "Procedure is marked to extract length from first argument, but the provided length was not set!");
            return false;
        }
    } else {
        ASSERT(node->proc->flags & FLAG_QUERYABLE_LENGTH);

        // Perform the call
        md_script_vis_t* old_vis = ctx->vis;
        token_t old_token = ctx->op_token;
        ctx->op_token = node->token;
        ctx->vis = NULL;
        int query_result = do_proc_call(NULL, node->proc, node->children, md_array_size(node->children), ctx);
        ctx->op_token = old_token;
        ctx->vis = old_vis;

        if (query_result >= 0) { // Zero length is valid
            type->dim[0] = query_result;
        } else {
            LOG_ERROR(ctx->ir, node->token, "Unexpected return value (%i) when querying procedure for array length.", query_result);
            return false;
        }
    }
    return true;
}

static bool finalize_type_array(type_info_t* type, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(type);
    ASSERT(node);
    ASSERT(ctx);
    ASSERT(ctx->temp_alloc);

    ast_node_t** const children = node->children;
    const size_t num_children = md_array_size(node->children);
    int length = 0;
    for (size_t i = 0; i < num_children; ++i) {
        type_info_t c_type = children[i]->data.type;
        if (is_variable_length(c_type)) {
            if (!finalize_type(&c_type, children[i], ctx)) {
                return false;
            }
        }
        if (!type_info_equal(type_info_element_type(c_type), type_info_element_type(*type))) {
            md_logf(MD_LOG_TYPE_DEBUG, "finalize_type_array: Base element mismatch");
            return false;
        }
        length += (int)type_info_array_len(c_type);
    }

    type->dim[0] = length;
    return true;
}

// This should be invoked whenever we have an variable length type and we need to determine what the node evaluates into within the current context / state
static bool finalize_type(type_info_t* type, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(type);
    ASSERT(node);
    ASSERT(ctx);

    switch (node->type) {
    case AST_ASSIGNMENT:
        ASSERT(node->children);
        ASSERT(md_array_size(node->children) == 2);
        return finalize_type(type, node->children[1], ctx);
    case AST_ARRAY:
        return finalize_type_array(type, node, ctx);
    case AST_IDENTIFIER:
        ASSERT(node->children);
        ASSERT(md_array_size(node->children) == 1);
        return finalize_type(type, node->children[0], ctx);
    case AST_ADD:
    case AST_SUB:
    case AST_MUL:
    case AST_DIV:
    case AST_AND:
    case AST_XOR:
    case AST_OR:
    case AST_NOT:
    case AST_EQ:
    case AST_NE:
    case AST_LE:
    case AST_GE:
    case AST_LT:
    case AST_GT:
    case AST_CAST:
    case AST_PROC_CALL:
        return finalize_type_proc(type, node, ctx);

    case AST_ARRAY_SUBSCRIPT: // Array subscripts should not be variable length, since the input should be statically evaluated so no out of bounds accesses occur...
    case AST_CONSTANT_VALUE:
    case AST_OUT:
    case AST_CONTEXT:   // We don't allow dynamic length expressions for Contexts for now... Perhaps in the future...
    default:
        ASSERT(false);
    }
    return false;
}

static bool finalize_proc_call(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(ctx);
    ASSERT(node->proc);
    ASSERT(node->type == AST_PROC_CALL);

    node->data.type = node->proc->return_type;

    const size_t num_args = md_array_size(node->children);
    ast_node_t** args = node->children;

    if (num_args > 0) {
        if (node->proc_flags & FLAG_SYMMETRIC_ARGS) {
            // If the symmetric flag is set, we swap the first two arguments
            ASSERT(num_args >= 2);
            ast_node_t* tmp_child = node->children[0];
            node->children[0] = node->children[1];
            node->children[1] = tmp_child;
            node->proc_flags &= ~FLAG_SYMMETRIC_ARGS; // Clear flags so we don't swap again
        }

        // Make sure we cast the arguments into the expected types of the procedure.
        for (size_t i = 0; i < num_args; ++i) {
            if (!is_type_directly_compatible(args[i]->data.type, node->proc->arg_type[i])) {
                // Types are not directly compatible, but should be implicitly convertible (otherwise the match should never have been found)
                ASSERT(is_type_implicitly_convertible(args[i]->data.type, node->proc->arg_type[i]));
                if (!convert_node(args[i], node->proc->arg_type[i], ctx)) {
                    return false;
                }
                node->flags    |= args[i]->flags & FLAG_AST_PROPAGATION_MASK;
                ctx->ir->flags |= args[i]->flags & FLAG_IR_PROPAGATION_MASK;
            }
        }

        if (node->proc->flags & FLAG_ARGS_EQUAL_LENGTH) {
            ASSERT(num_args > 1); // This flag should only be set for procedures with multiple arguments
                                  // Make sure that the input arrays are of equal length. Otherwise create an error.
            const int expected_len = type_info_array_len(args[0]->data.type);
            for (size_t i = 1; i < num_args; ++i) {
                int len = type_info_array_len(args[i]->data.type);
                if (len != expected_len) {
                    LOG_ERROR(ctx->ir, node->token, "Expected array-length of arguments to match. arg 0 has length %i, arg %i has length %i.", expected_len, (int)i, len);
                    return false;
                }
            }
        }
    }

    // Propagate procedure flags
    node->flags    |= node->proc->flags & FLAG_AST_PROPAGATION_MASK;
    ctx->ir->flags |= node->proc->flags & FLAG_IR_PROPAGATION_MASK;

    // @TODO: Test if all contexts are equivalent
    // In such case, we can make an exception from the DYNAMIC_LENGTH flag
    if (node->num_contexts > 0 && !(node->flags & FLAG_CONTEXTS_EQUIVALENT) && (node->proc->flags & FLAG_QUERYABLE_LENGTH)) {
        node->flags |= FLAG_DYNAMIC_LENGTH;
        return true;
    }
    
    /*
    // Need to flag with DYNAMIC_LENGTH if we evaluate within a context since the return type length may depend on its context.
    if (ctx->mol_ctx && (node->proc->flags & FLAG_QUERYABLE_LENGTH)) {
        node->flags |= FLAG_DYNAMIC_LENGTH;
        return true;
    }
    */

    if (!(ctx->eval_flags & EVAL_FLAG_NO_LENGTH_CHECK) && is_variable_length(node->data.type)) {
        if (node->proc->flags & FLAG_QUERYABLE_LENGTH) {
            // Try to resolve length by calling the procedure
            static_backchannel_t channel = {0};
            static_backchannel_t* prev_channel = ctx->backchannel;
            token_t prev_op_token = ctx->op_token;
            ctx->backchannel = &channel;
            int query_result = do_proc_call(NULL, node->proc, node->children, md_array_size(node->children), ctx);
            ctx->backchannel = prev_channel;
            ctx->op_token = prev_op_token;

            if (channel.flags & FLAG_DYNAMIC_LENGTH) {
                // If we fail to deduce the length it is still valid in this scenario
                // We just postpone the length-evaluation to run-time
                node->flags |= FLAG_DYNAMIC_LENGTH;
                return true;
            }

            if (query_result >= 0) { // Zero length is valid
                node->data.type.dim[0] = query_result;
                return true;
            } else {
                //node-> flags |= FLAG_DYNAMIC_LENGTH;
                //return true;
                LOG_ERROR(ctx->ir, node->token, "Failed to determine length of procedure return type!");
                return false;
            }
        } else if (node->proc->flags & FLAG_RET_AND_FIRST_ARG_EQUAL_LENGTH) {
            // We can deduce the length of the array by using the type of the first argument
            ASSERT(num_args > 0);
            ASSERT(node->proc->arg_type[0].base_type == args[0]->data.type.base_type);
            node->data.type.dim[0] = type_info_array_len(args[0]->data.type);
            return true;
        } else {
            LOG_ERROR(ctx->ir, node->token, "Procedure returns variable length, but its length cannot be determined.");
            return false;
        }
    }

    return true;
}

static bool static_check_children(ast_node_t* node, eval_context_t* ctx) {
    const size_t num_children = md_array_size(node->children);
    for (size_t i = 0; i < num_children; ++i) {
        ASSERT(node != node->children[i]);
        if (static_check_node(node->children[i], ctx)) {
            node->flags |= node->children[i]->flags & FLAG_AST_PROPAGATION_MASK;
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
    ASSERT(ctx);

    // @NOTE: if we have a bitwise operator in the tree, we flag this
    // so that we can shortcut some ambigious operations such as 'protein' 'residue' etc.
    // That they can be concatenated into a single bitfield instead of an array of bitfields.
    uint32_t backup_flags = ctx->eval_flags;
    if (is_bitwise_operator(node->type)) {
        ctx->eval_flags |= EVAL_FLAG_FLATTEN;
    }

    bool result = false;

    if (static_check_children(node, ctx)) {
        const size_t num_args = md_array_size(node->children);
        ast_node_t** arg = node->children;

        type_info_t arg_type[2] = {arg[0]->data.type, num_args > 1 ? arg[1]->data.type : (type_info_t){0}};
        procedure_match_result_t res = find_operator_supporting_arg_types(node->type, arg_type, num_args, true);
        if (res.success) {
            if (num_args == 2) {
                switch (node->type) {
                // These are the only operators which potentially modify or conserve units.
                case AST_ADD: node->data.unit = md_unit_add(arg[0]->data.unit, arg[1]->data.unit); break;
                case AST_SUB: node->data.unit = md_unit_sub(arg[0]->data.unit, arg[1]->data.unit); break;
                case AST_MUL: node->data.unit = md_unit_mul(arg[0]->data.unit, arg[1]->data.unit); break;
                case AST_DIV: node->data.unit = md_unit_div(arg[0]->data.unit, arg[1]->data.unit); break;
                default: break;
                }
            }

            node->type = AST_PROC_CALL;
            node->proc = res.procedure;
            node->proc_flags = res.flags;
            result = finalize_proc_call(node, ctx);
        } else {
            char arg_type_str[2][64] = {"empty", "empty"};
            for (size_t i = 0; i < num_args; ++i) {
                if (arg[i]) print_type_info(arg_type_str[i], ARRAY_SIZE(arg_type_str[i]), arg[i]->data.type);                
            }
            LOG_ERROR(ctx->ir, node->token,
                "Could not find supporting operator '%s' with left hand side type (%s) and right hand side type (%s)",
                get_token_type_str(node->token.type), arg_type_str[0], arg_type_str[1]);
        }
    }

    // Restore flags
    ctx->eval_flags = backup_flags;
    return result;
}

static size_t extract_argument_types(type_info_t arg_type[], size_t cap, const ast_node_t* node) {
    ASSERT(arg_type);
    ASSERT(node);
    const size_t num_args = MIN(cap, md_array_size(node->children));
    for (size_t i = 0; i < num_args; ++i) {
        arg_type[i] = node->children[i]->data.type;
    }
    return num_args;
}

static size_t print_argument_list(char* buf, size_t cap, const type_info_t arg_type[], size_t num_args) {
    size_t len = 0;
    for (size_t i = 0; i < num_args; ++i) {
        len += print_type_info(buf + len, cap - len, arg_type[i]);
        if (i + 1 != num_args) {
            len += snprintf(buf + len, cap - len, ", ");
        }
    }
    return len;
}

static table_t* find_table(md_script_ir_t* ir, str_t name) {
    for (size_t i = 0; i < md_array_size(ir->tables); ++i) {
        if (str_eq(ir->tables[i].name, name)) {
            return &ir->tables[i];
        }
    }
    return NULL;
}

// There are some heuristics/mumbo jumbo here,
// we try to match the frame times of the trajectory and the table.
// If the table has more frames than the trajectory, we assume that the table is a subset of the trajectory.
static bool frame_times_compatible(const double* traj_times, size_t num_traj_frames, const double* table_times, size_t num_table_frames) {
    if (!traj_times) return false;
    if (!table_times) return false;

    if (num_table_frames < num_traj_frames) {
        return false;
    }

    // We want to ensure that every frame time in the trajectory is present in the table.
    size_t table_idx = 0;
    for (size_t i = 0; i < num_traj_frames; ++i) {
        while (table_idx < num_table_frames && table_times[table_idx] + TIMESTAMP_ERROR_MARGIN < traj_times[i]) {
            ++table_idx;
        }
        if (!timestamps_approx_equal(table_times[table_idx], traj_times[i])) {
            return false;
        }
    }

    return true;
}

static bool frame_times_compatible_f(const double* traj_times, size_t num_traj_frames, const float* table_times, size_t num_table_frames) {
    if (!traj_times) return false;
    if (!table_times) return false;

    if (num_table_frames < num_traj_frames) {
        return false;
    }

    // We want to ensure that every frame time in the trajectory is present in the table.
    size_t table_idx = 0;
    for (size_t i = 0; i < num_traj_frames; ++i) {
        while (table_idx < num_table_frames && (double)table_times[table_idx] + TIMESTAMP_ERROR_MARGIN < traj_times[i]) {
            ++table_idx;
        }
        if (!timestamps_approx_equal((double)table_times[table_idx], traj_times[i])) {
            return false;
        }
    }

    return true;
}

// Find the possible occurrence of parenthesis and try to match its contents to a unit.
static md_unit_t extract_unit_from_label(str_t label) {
    md_unit_t unit = md_unit_none();
    size_t beg_loc, end_loc;
    if (str_find_char(&beg_loc, label, '(') && str_find_char(&end_loc, label, ')')) {
        str_t unit_str = str_substr(label, beg_loc + 1, end_loc - beg_loc - 1);
        md_unit_t xvg_unit = md_unit_from_string(unit_str);
        if (!md_unit_empty(xvg_unit)) {
            unit = xvg_unit;
        }
    }
    return unit;
}

static table_t* import_table(md_script_ir_t* ir, token_t tok, str_t path_to_file, const md_trajectory_i* traj) {
    str_t ext;
    if (!extract_ext(&ext, path_to_file)) {
        LOG_ERROR(ir, tok, "Could not extract extension from file path '"STR_FMT"'", STR_ARG(path_to_file));
        return NULL;
    }
    table_t* table = NULL;
    md_unit_t traj_time_unit = md_trajectory_time_unit(traj);
    const size_t num_frames = md_trajectory_num_frames(traj);
    bool traj_has_time = !md_unit_empty(traj_time_unit);
   
    if (str_eq_cstr_ignore_case(ext, "edr")) {
        md_edr_energies_t edr = {0};
        if (md_edr_energies_parse_file(&edr, path_to_file, md_heap_allocator)) {
            bool success = true;
            if (traj_has_time) {
                if (!frame_times_compatible(md_trajectory_frame_times(traj), num_frames, edr.frame_time, edr.num_frames)) {
                    LOG_ERROR(ir, tok, "EDR file is not compatible with loaded trajectory, could not match timestamps");
                    success = false;
                }
            } else {
                if (num_frames != edr.num_frames) {
					LOG_ERROR(ir, tok, "EDR file is not compatible with loaded trajectory, number of frames did not match");
					success = false;
				}
            }
            if (success) {
                table = md_array_push(ir->tables, (table_t){.num_values = edr.num_frames}, ir->arena);
                table->name = str_copy(path_to_file, ir->arena);

                table_push_field_d(table, STR_LIT("Time"), md_unit_pikosecond(), edr.frame_time, edr.num_frames, ir->arena);
                for (size_t i = 0; i < edr.num_energies; ++i) {
                    table_push_field_f(table, edr.energy[i].name, edr.energy[i].unit, edr.energy[i].values, edr.num_frames, ir->arena);
                }
            }
        }
        md_edr_energies_free(&edr);
    } else if (str_eq_cstr_ignore_case(ext, "xvg")) {
        md_xvg_t xvg = {0};
        if (md_xvg_parse_file(&xvg, path_to_file, md_heap_allocator)) {
            bool success = true;
            bool xvg_has_time = str_eq_cstr_n_ignore_case(xvg.header_info.xaxis_label, "time", 4);
            if (traj_has_time && xvg_has_time) {
                if (!frame_times_compatible_f(md_trajectory_frame_times(traj), num_frames, xvg.fields[0], xvg.num_values)) {
                    LOG_ERROR(ir, tok, "XVG file is not compatible with loaded trajectory, could not match timestamps");
                    success = false;
                }
            } else {
            	if (num_frames != xvg.num_values) {
                	LOG_ERROR(ir, tok, "XVG file is not compatible with loaded trajectory, number of frames did not match");
                	success = false;
                }
            }
            if (success) {
                table = md_array_push(ir->tables, (table_t){.num_values = xvg.num_values}, ir->arena);
                table->name = str_copy(path_to_file, ir->arena);
                
                size_t i = 0;
                if (xvg_has_time) {
                    md_unit_t time_unit = extract_unit_from_label(xvg.header_info.xaxis_label);
                    table_push_field_f(table, STR_LIT("Time"), time_unit, xvg.fields[0], xvg.num_values, ir->arena);
                    i = 1;
                }

                md_unit_t unit = extract_unit_from_label(xvg.header_info.yaxis_label);
                for (; i < xvg.num_fields; ++i) {
                    str_t name = STR_LIT("");
                    if (0 < i && i-1 < md_array_size(xvg.header_info.legends)) {
                        name = xvg.header_info.legends[i-1];
                    }
                    table_push_field_f(table, name, unit, xvg.fields[i], xvg.num_values, ir->arena);
                }
            }
            md_xvg_free(&xvg, md_heap_allocator);
        }
    } else if (str_eq_cstr_ignore_case(ext, "csv")) {
        md_csv_t csv = {0};
        if (md_csv_parse_file(&csv, path_to_file, md_heap_allocator)) {
            bool success = true;
            bool csv_has_time = csv.field_names && str_eq_cstr_n_ignore_case(csv.field_names[0], "time", 4);

            if (traj_has_time && csv_has_time) {
                if (!frame_times_compatible_f(md_trajectory_frame_times(traj), num_frames, csv.field_values[0], csv.num_values)) {
                    LOG_ERROR(ir, tok, "CSV file is not compatible with loaded trajectory, could not match timestamps");
                    success = false;
                }
            } else {
                if (num_frames != csv.num_values) {
					LOG_ERROR(ir, tok, "CSV file is not compatible with loaded trajectory, number of frames did not match");
					success = false;
				}
            }
            if (success) {
                table = md_array_push(ir->tables, (table_t){.num_values = csv.num_values}, ir->arena);
                table->name = str_copy(path_to_file, ir->arena);

                size_t i = 0;
                if (csv_has_time) {
                    md_unit_t time_unit = extract_unit_from_label(csv.field_names[0]);
                    table_push_field_f(table, STR_LIT("Time"), time_unit, csv.field_values[0], csv.num_values, ir->arena);
                    i = 1;
                }

                for (; i < csv.num_fields; ++i) {
                    md_unit_t unit = md_unit_none();
                    str_t name = STR_LIT("");
                    if (csv.field_names) {
                        name = csv.field_names[i];
                        unit = extract_unit_from_label(csv.field_names[i]);
                    }
                    table_push_field_f(table, name, unit, csv.field_values[i], csv.num_values, ir->arena);
                }
            }

            md_csv_free(&csv, md_heap_allocator);
        }
    } else {
        LOG_ERROR(ir, (token_t){0}, "import: unsupported file extension '"STR_FMT"'", (int)ext.len, ext.ptr);
        return NULL;
    }

    return table;
}

static void swap_int(int* a, int* b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

// The number of matches may be more than number of candidates, therefore the true number of matches computed are returned
static size_t str_find_n_best_matches(int match_idx[], size_t num_idx, str_t str, str_t candidates[], size_t num_candidates) {
    num_idx = MIN(num_idx, num_candidates);
    int* distances = md_temp_push(sizeof(int) * num_idx);
    for (size_t i = 0; i < num_idx; ++i) {
        distances[i] = INT32_MAX;
    }
    for (size_t i = 0; i < num_candidates; ++i) {
        str_t can = candidates[i];
        int dist = str_edit_distance(str, can);
        
        // If the first characters differ, then add a penalty
        const size_t min_len = MIN(str.len, can.len);
        for (size_t j = 0; j < MIN(3, min_len); ++j) {
            if (to_lower(str.ptr[j]) != to_lower(can.ptr[j])) {
                dist += 10;
            }
        }
        
        if (dist < distances[num_idx - 1]) {
            // Insertion sort
            distances[num_idx - 1] = dist;
            match_idx[num_idx - 1] = (int)i;
            size_t j = num_idx - 1;
            while (j > 0 && distances[j] < distances[j - 1]) {
                swap_int(&distances[j], &distances[j-1]);
                swap_int(&match_idx[j], &match_idx[j-1]);
                --j;
            }
        }
    }
    return num_idx;
}

static bool static_check_import(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(ctx);

    if (!static_check_children(node, ctx)) {
        return false;
    }

    const size_t num_args = md_array_size(node->children);
    ast_node_t** const args = node->children;

    if (num_args == 0) {
        LOG_ERROR(ctx->ir, node->token, "import: expected at least one argument");
        return false;
    }

    if (num_args > 2) {
        LOG_ERROR(ctx->ir, node->token, "import: to many arguments");
        return false;
    }

    if (!is_type_equivalent(args[0]->data.type, (type_info_t)TI_STRING) || (args[0]->type != AST_CONSTANT_VALUE)) {
        LOG_ERROR(ctx->ir, node->token, "import: requires as minimum a constant string containing the path to the file to import");
        return false;
    }

    str_t full_path = md_path_make_canonical(args[0]->value._string, ctx->temp_alloc);
    if (str_empty(full_path)) {
        LOG_ERROR(ctx->ir, node->token, "import: failed to resolve path '"STR_FMT"'", (int)full_path.len, full_path.ptr);
        return false;
    }

    table_t* table = find_table(ctx->ir, full_path);
    if (!table) {
        table = import_table(ctx->ir, node->token, full_path, ctx->traj);
        if (!table) {
            LOG_ERROR(ctx->ir, node->token, "import: failed to import file '"STR_FMT"'", (int)full_path.len, full_path.ptr);
            return false;
        }
    }

    md_array(int) field_indices = 0;

    // Resolve mappings if supplied (what fields to load)
    if (num_args == 2) {
        const ast_node_t* arg = args[1];
        type_info_t type = arg->data.type;

        if (arg->flags & FLAG_DYNAMIC) {
            LOG_ERROR(ctx->ir, node->token, "import: second argument must be a constant");
            return false;
        }
       
        if (is_type_directly_compatible(type, (type_info_t)TI_INT_ARR)) {
            const size_t count = element_count(arg->data);
            const int* ints = as_int_arr(arg->data);
            for (size_t i = 0; i < count; ++i) {
                int idx = ints[i];
                if (idx < 1 || (int)table->num_fields < idx) {
                    LOG_ERROR(ctx->ir, node->token, "import: second argument must be a valid field index within range [1, %d]",
                        (int)table->num_fields);
                    return false;
                }
                md_array_push(field_indices, idx - 1, ctx->ir->arena);
            }
        } else if (is_type_directly_compatible(type, (type_info_t)TI_IRANGE_ARR)) {
            const size_t count = element_count(arg->data);
            const irange_t* ranges = as_irange_arr(arg->data);

            for (size_t i = 0; i < count; ++i) {
                irange_t range = ranges[i];
                if (range.end < range.beg) {
                    LOG_ERROR(ctx->ir, node->token, "import: invalid range");
                    return false;
                }
                if (range.beg == INT32_MIN) range.beg = 1;
                if (range.end == INT32_MAX) range.end = (int)table->num_fields;
                for (int idx = range.beg; idx <= range.end; ++idx) {
                    if (idx < 1 || (int)table->num_fields < idx) {
                        LOG_ERROR(ctx->ir, node->token, "import: second argument must be a valid field index within range [1, %d]",
                            (int)table->num_fields);
                        return false;
                    }
                    md_array_push(field_indices, idx - 1, ctx->ir->arena);
                }
            }
        } else if (is_type_directly_compatible(type, (type_info_t)TI_STRING_ARR)) {
            const size_t count = element_count(arg->data);
            const str_t* strings = as_string_arr(arg->data);
            for (size_t i = 0; i < count; ++i) {
                str_t str = strings[i];
                if (str.len == 0) {
                    LOG_ERROR(ctx->ir, node->token, "import: empty string not allowed");
                    return false;
                }
                int match_idx = -1;
                for (int j = 0; j < (int)table->num_fields; ++j) {
                    if (str_eq(str, table->field_names[j])) {
                        match_idx = j;
                        break;
                    }
                }
                if (match_idx == -1) {
                    char buf[512];
                    int len = 0;
                    int best_idx[3];
                    const size_t num_best_matches = str_find_n_best_matches(best_idx, ARRAY_SIZE(best_idx), str, table->field_names, table->num_fields);
                    for (size_t j = 0; j < num_best_matches; ++j) {
                        int idx = best_idx[j];
                        str_t name = table->field_names[idx];
                        len += snprintf(buf + len, sizeof(buf) - len, "'"STR_FMT"'\n", (int)name.len, name.ptr);
                    }
                    
                    LOG_ERROR(ctx->ir, node->token, "import: field '"STR_FMT"' not found, closest field names are:\n%s", (int)str.len, str.ptr, buf);
                    return false;
                }
                md_array_push(field_indices, match_idx, ctx->ir->arena);
            }
        } else {
            LOG_ERROR(ctx->ir, node->token, "import: unexpected type of second argument");
            return false;
        }
    } else {
        ASSERT(num_args == 1);
        // Setup node table fields to hold all available fields within the imported table
        for (size_t i = 0; i < table->num_fields; ++i) {
            md_array_push(field_indices, (int)i, ctx->ir->arena);
        }
    }

    if (md_array_size(field_indices) == 0) {
        LOG_ERROR(ctx->ir, node->token, "import: no fields matched");
        return false;
    }

    md_unit_t unit = table->field_units[field_indices[0]];
    for (size_t i = 0; i < md_array_size(field_indices); ++i) {
        if (!md_unit_equal(unit, table->field_units[field_indices[i]])) {
            // If we have conflicting units, we make it unitless and let the user know.
            //LOG_WARNING(ctx->ir, node->token, "import: conflicting units, perhaps separate the import into two separate");
            unit = md_unit_none();
            break;
        }
    }

    node->type = AST_TABLE;
    node->data.type = (type_info_t) {TYPE_FLOAT, {(int)md_array_size(field_indices)}};
    node->data.unit = unit;
    node->table = table;
    node->table_field_indices = field_indices;
    node->flags |= FLAG_DYNAMIC; // It is not a constant value, it will change depending on what time we evaluate it.

    return true;
}

static bool static_check_flatten(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(node->type == AST_FLATTEN);
    ASSERT(ctx);

    size_t num_args = md_array_size(node->children);
    if (num_args != 1) {
        LOG_ERROR(ctx->ir, node->token, "Expected 1 argument to built in procedure flatten, got %zu", num_args);
        return false;
    }

    uint32_t flags = ctx->eval_flags;
    ctx->eval_flags |= EVAL_FLAG_FLATTEN;
    if (!static_check_children(node, ctx)) {
        return false;
    }
    ctx->eval_flags = flags;

    const ast_node_t* arg = node->children[0];
    if (arg->data.type.base_type == TYPE_UNDEFINED) {
        LOG_ERROR(ctx->ir, arg->token, "Argument to built in procedure 'flatten' has undefined base type");
        return false;
    }

    // @NOTE: Flatten conceptually flattens the dimensions into a one dimensional array
    node->data = arg->data;
    MEMSET(node->data.type.dim, 0, sizeof(int) * MAX_NUM_DIMS);
    if (arg->data.type.dim[0] == -1) {
        node->data.type.dim[0] = -1;
    } else {
        node->data.type.dim[0] = dim_size(arg->data.type.dim);
    }

    if (node->data.type.base_type == TYPE_BITFIELD) {
        // If we are dealing with bitfields, we flatten even further by collapsing it into a single bitfield
        node->data.type.dim[0] = 1;
    }

    return true;
}

static bool static_check_transpose(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(node->type == AST_TRANSPOSE);
    ASSERT(ctx);

    size_t num_args = md_array_size(node->children);
    if (num_args != 1) {
        LOG_ERROR(ctx->ir, node->token, "Expected 1 argument to built in procedure 'transpose', got %zu", num_args);
        return false;
    }

    if (!static_check_children(node, ctx)) {
        return false;
    }

    const ast_node_t* arg = node->children[0];
    if (arg->data.type.base_type == TYPE_UNDEFINED) {
        LOG_ERROR(ctx->ir, arg->token, "Argument to built in procedure 'transpose' has undefined base type");
        return false;
    }

    int ndim = dim_ndims(arg->data.type.dim);
    if (ndim == 0) {
        LOG_ERROR(ctx->ir, arg->token, "Undefined number of dimensions of argument supplied to transpose");
        return false;
    }

    if (arg->data.type.dim[0] == -1) {
        LOG_ERROR(ctx->ir, arg->token, "Built in procedure 'transpose' does not support dynamic length types");
        return false;
    }

    else if (ndim == 1) {
        // NOOP
        // @TODO: Remove this node from the tree and directly link the parent to arg
        node->data = arg->data;
        node->data.ptr = 0;
        node->data.size = 0;
        return true;
    }
    else if (ndim == 2) {
        // Do Implicit Transpose (swap row and col)
        node->data = arg->data;
        node->data.ptr = 0;
        node->data.size = 0;
        node->data.type.dim[0] = arg->data.type.dim[1];
        node->data.type.dim[1] = arg->data.type.dim[0];
        return true;
    }
    else {
        // We can conceptually support this if additional axis-permutation arguments are given.
        LOG_ERROR(ctx->ir, arg->token, "objects of dim > 2 are currently not supported in transpose operation");
        return false;
    }
}

static bool static_check_proc_call(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node->type == AST_PROC_CALL);
    uint32_t backup_flags = ctx->eval_flags;

    bool result = true;
    // @NOTE: We do not want to overwrite the procedure if already assigned, as it might have been assigned as an operator (which is a proc call)
    if (!node->proc) {
        // No point in letting expressions statically evaluate and store its data within the tree at this point
        // Conversions and other stuff may occur later
        uint32_t flags = ctx->eval_flags;
        ctx->eval_flags |= EVAL_FLAG_NO_STATIC_EVAL | EVAL_FLAG_NO_LENGTH_CHECK;
        result = static_check_children(node, ctx);
        ctx->eval_flags = flags;
        if (result) {
            const size_t num_args = md_array_size(node->children);
            const str_t proc_name = node->ident;

            // One or more arguments
            ASSERT(num_args < MAX_SUPPORTED_PROC_ARGS);
            type_info_t arg_type[MAX_SUPPORTED_PROC_ARGS];
            
            size_t arg_len = extract_argument_types(arg_type, ARRAY_SIZE(arg_type), node);
            ASSERT(arg_len == num_args);
            (void)arg_len;

            procedure_match_result_t res = find_procedure_supporting_arg_types(proc_name, arg_type, num_args, true);
            result = result && res.success;

            if (res.success) {
                node->proc = res.procedure;
                node->proc_flags = res.flags;
            } else {
                if (num_args == 0) {
                    LOG_ERROR(ctx->ir, node->token,
                        "Could not find matching procedure '"STR_FMT"' which takes no arguments", STR_ARG(proc_name));
                } else {
                    char buf[512];
                    print_argument_list(buf, ARRAY_SIZE(buf), arg_type, num_args);

                    LOG_ERROR(ctx->ir, node->token,
                        "Could not find matching procedure '"STR_FMT"' which takes the following argument(s): %s",
                        STR_ARG(proc_name), buf);
                }
            }
        }
    }

    // Resolve type of returned data
    if (result && node->proc) {
        // Set eval flags from procedures before performing the final type check on arguments
        if (node->proc->flags & FLAG_FLATTEN) {
            ctx->eval_flags |= EVAL_FLAG_FLATTEN;
        }
        /*
        else if (node->proc->flags & FLAG_NO_FLATTEN) {
            ctx->eval_flags = ctx->eval_flags & (~EVAL_FLAG_FLATTEN);
        }
        */
        
        // Perform new child check here since children may have changed due to conversions etc.
        // Also allow for any static evaluation to occur and for length to be determined etc.
        result = result && static_check_children(node, ctx);
            
        const size_t num_args = md_array_size(node->children);
        ast_node_t** args = node->children;
        for (size_t i = 0; i < num_args; ++i) {
            if (node->proc->arg_type[i].base_type == TYPE_COORDINATE) {
                if (args[i]->data.type.base_type == TYPE_FLOAT && dim_ndims(args[i]->data.type.dim) < 2) {
                    // I.e. if the procedures expects a coordinate array, but the input is a single coordinate (float[3]), we need to convert it to an array (float[1][3])
                    // To properly propagate the number of supplied coordinates (1)
                    dim_shift_right(args[i]->data.type.dim);
                }
                // Set coordinate flag, so we can visualize it properly later, since the coordinate type is a pseudo type
                args[i]->flags |= FLAG_COORDINATE;
            }
        }

        // Finalize the procedure call and its type
        result = result && finalize_proc_call(node, ctx);

        // Perform static validation by evaluating with NULL ptr
        if (result && node->proc->flags & FLAG_STATIC_VALIDATION) {
            static_backchannel_t* prev_channel = ctx->backchannel;
            static_backchannel_t channel = {0};
            ctx->backchannel = &channel;
            result = result && evaluate_proc_call(NULL, node, ctx);
            node->data.unit = channel.unit;
            node->data.value_range = channel.value_range;
            ctx->backchannel = prev_channel;
        }
    }

    ctx->eval_flags = backup_flags;
    return result;
}

static bool static_check_constant_value(ast_node_t* node, eval_context_t* ctx) {
    (void)ctx; 
    ASSERT(node && node->type == AST_CONSTANT_VALUE);
    ASSERT(node->data.type.base_type != TYPE_UNDEFINED);
    ASSERT(is_scalar(node->data.type));

    if (node->data.type.base_type == TYPE_IRANGE) {
        // Make sure the range is given in an ascending format, lo:hi
        irange_t rng = node->value._irange;
        if (rng.beg > rng.end) {
            LOG_ERROR(ctx->ir, node->token, "The range is invalid, a range must have an ascending format, did you mean '%i:%i'?", rng.end, rng.beg);
            return false;
        }
    }

    node->data.ptr  = &node->value;
    node->data.size = base_type_element_byte_size(node->data.type.base_type);
    node->flags |= FLAG_CONSTANT;
    
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

        const size_t num_elem = md_array_size(node->children);
        ast_node_t**     elem = node->children;

        int array_len = 0;
        type_info_t array_type = {0};
        md_unit_t array_unit = md_unit_none();
        
        // Pass 1: find a suitable base_type for the array
        for (size_t i = 0; i < num_elem; ++i) {
            int elem_len = type_info_array_len(elem[i]->data.type);
            type_info_t elem_type = type_info_element_type(elem[i]->data.type);
            if (is_variable_length(elem_type)) {
                array_len = -1;
            }

            /*
            // Normalize the type and dimensions here for the length to be computed correctly
            if (elem_type.base_type != TYPE_BITFIELD) {
                int ndim = dim_ndims(elem_type.dim);
                if (ndim == 1 || (ndim > 1 && elem_type.dim[0] != 1)) {
                    dim_shift_right(elem_type.dim);
                }
                // Turn the normalized form into a single entity
                elem_type.dim[0] = 1;
            }
            */
            if (elem_type.base_type == TYPE_BITFIELD) {
                // Bitfield
                dim_flatten(elem_type.dim);
            }

            if (array_len >= 0) array_len += elem_len;

            if (array_type.base_type == TYPE_UNDEFINED) {
                // Assign type from first element
                array_type = elem_type;
            }
            else {
                if (array_type.base_type == TYPE_BITFIELD) {
                    // Bitfields are handled separately and concatenated
                    if (elem_type.base_type != TYPE_BITFIELD) {
                        LOG_ERROR(ctx->ir, elem[i]->token, "Incompatible types wihin array construct");
                        return false;
                    }
                } else {
                    // Other standard types are packed into an array, but require type checking to ensure that the element types are the same
                    // Or can be converted into some common base element
                    if (!is_type_directly_compatible(elem_type, array_type)) {
                        if (is_type_implicitly_convertible(elem[i]->data.type, array_type)) {
                            // We can convert the element to the array type.   
                        }
                        else if (is_type_implicitly_convertible(array_type, elem_type)) {
                            // Option 2: We can convert the array type into the elements type (e.g int to float)
                            // Promote the type of the array
                            array_type = elem_type;

                            // Just to make sure, we recheck all elements up until this one
                            for (size_t j = 0; j < i; ++j) {
                                if (!is_type_directly_compatible(elem_type, array_type) &&
                                    !is_type_implicitly_convertible(elem_type, array_type)) {
                                    LOG_ERROR(ctx->ir, elem[i]->token, "Incompatible types wihin array construct");
                                    return false;
                                }
                            }
                        }
                        else {
                            // Incompatible types...
                            LOG_ERROR(ctx->ir, elem[i]->token, "Incompatible types wihin array construct");
                            return false;
                        }
                    }
                }
            }
        }

        // Pass 2: Perform implicit conversions of nodes if required
        for (size_t i = 0; i < num_elem; ++i) {
            type_info_t converted_type = array_type;
            converted_type.dim[0] = type_info_array_len(elem[i]->data.type);

            if (!is_type_directly_compatible(elem[i]->data.type, converted_type) &&
                is_type_implicitly_convertible(elem[i]->data.type, converted_type)) {
                if (!convert_node(elem[i], converted_type, ctx)) {
                    return false;
                }
            }
        }

        // Pass 3: Deduce the unit, if defined and the same.
        if (num_elem > 0) {
            array_unit = elem[0]->data.unit;
            for (size_t i = 1; i < num_elem; ++i) {
                if (!md_unit_equal(elem[i]->data.unit, array_unit)) {
                    array_unit = md_unit_none();
                    break;
                }
            }
        }

        // Finalize the type of the array (dimensionality)
        array_type.dim[0] = (int32_t)array_len;
        node->data.type = array_type;
        node->data.unit = array_unit;

        bool children_constant = true;
        for (size_t i = 0; i < num_elem; ++i) {
            if (elem[i]->type != AST_CONSTANT_VALUE) {
                children_constant = false;
                break;
            }
        }

        if (children_constant && !(ctx->eval_flags & EVAL_FLAG_NO_STATIC_EVAL)) {
            return static_eval_node(node, ctx);
        }

        return true;
    }
    return false;
}

static bool static_check_array_subscript(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY_SUBSCRIPT && node->children);

    // In array subscripts, only integers and ranges should be allowed.
    // For now, we only support SINGLE integers or ranges
    // If the node has a dynamic length, then we produce an error.
    // @TODO: In the future, se should expand this to support mixed integers and ranges and create a proper subset of the original data
    // @TODO: evaluate the children in order to determine the length and values (to make sure they are in range)

    const size_t num_elem = md_array_size(node->children);
    ast_node_t**     elem = node->children;

    if (num_elem < 2) {
        LOG_ERROR(ctx->ir, node->token, "Missing arguments in array subscript");
        return false;
    }
    //else if (num_elem > 2) {
    //    LOG_ERROR(ctx->ir, elem[2]->token, "Only single entries are allowed inside array subscript");
    //    return false;
    //}

    uint32_t eval_flags = ctx->eval_flags;
    ctx->eval_flags = ctx->eval_flags & ~(EVAL_FLAG_FLATTEN);
    if (!static_check_node(elem[0], ctx)) return false;
    ctx->eval_flags = eval_flags;

    for (size_t i = 1; i < num_elem; ++i) {
        if (!static_check_node(elem[i], ctx)) return false;
    }

    ast_node_t*  lhs  = elem[0];
    ast_node_t** args = elem + 1;
    const size_t num_args = num_elem - 1;

    if (is_variable_length(lhs->data.type)) {
        LOG_ERROR(ctx->ir, lhs->token, "Array subscript operator can only be applied to expressions which have a static length");
        return false;
    }

    int num_dim = dim_ndims(lhs->data.type.dim);
    if (num_args > 1 && num_args != (size_t)num_dim) {
    	LOG_ERROR(ctx->ir, elem[1]->token, "Invalid number of arguments (%i) in array subscript, expected number of arguments to match number of dimensions of lhs (%i)", (int)num_args, num_dim);
		return false;
    }

    type_info_t result_type = lhs->data.type;
    for (size_t i = 0; i < num_args; ++i) {
        int lhs_dim = lhs->data.type.dim[i];
        if (lhs_dim <= 0) {
            LOG_ERROR(ctx->ir, lhs->token, "Unexpected length (%i) in dimension (%i) of lhs", lhs_dim, (int)i);
            return false;
        }
        ast_node_t* arg = args[i];
        if (arg->flags & FLAG_DYNAMIC) {
            LOG_ERROR(ctx->ir, args[i]->token, "Only static expressions are allowed within array subscript");
            return false;
        }

        if (is_scalar(arg->data.type)) {
            if (arg->data.type.base_type == TYPE_INT || arg->data.type.base_type == TYPE_IRANGE) {
                irange_t range = {0};
                if (arg->data.type.base_type == TYPE_INT) {
                    range.beg = arg->value._int;
                    range.end = arg->value._int;
                } else {
                    range = arg->value._irange;
                    if (range.beg == INT32_MIN) {
                        range.beg = 1;
                        arg->value._irange.beg = 1;
                    }
                    if (range.end == INT32_MAX) {
                        range.end = lhs_dim;
                        arg->value._irange.end = range.end;
                    }
                }
                if (range.beg <= range.end && 1 <= range.beg && range.end <= lhs_dim) {
                    node->subscript_ranges[i] = (irange_t){range.beg - 1, range.end};
                    result_type.dim[i] = range.end - range.beg + 1;
                } else {
                    LOG_ERROR(ctx->ir, arg->token, "Invalid Array subscript range");
                    return false;
                }
            } else {
                LOG_ERROR(ctx->ir, arg->token, "Only int and int-ranges are allowed inside array subscript");
                return false;
            }
        } else {
            LOG_ERROR(ctx->ir, arg->token, "No arrays are allowed inside array subscript");
            return false;
        }
    }

    // Fill in missing dimensions with complete ranges
    for (size_t i = num_args; i < num_dim; ++i) {
        node->subscript_ranges[i].beg = 0;
        node->subscript_ranges[i].end = lhs->data.type.dim[i];
    }

    // Remove leading ones in the dimensions
    for (size_t i = 0; i < MAX_NUM_DIMS - 1; ++i) {
        if (result_type.dim[i] == 0) break;
        if (result_type.dim[i] == 1 && result_type.dim[i+1] > 0) {
            dim_shift_left(result_type.dim);
		}
    }

    // SUCCESS!
    // @NOTE Subscript dim is not just the number of args supplied,
    // It is the number of dimensions
    node->subscript_dim = num_dim;
    node->flags = (lhs->flags & FLAG_AST_PROPAGATION_MASK);
    node->data.type = result_type;
    node->data.unit = lhs->data.unit;
    return true;
}

static bool static_check_identifier_reference(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_IDENTIFIER);
    
    if (str_empty(node->ident)) {
        node->ident = node->token.str;
    }

    identifier_t* ident = get_identifier(ctx->ir, node->ident);
    if (ident && ident->node) {
        if (ident->node->data.type.base_type == TYPE_UNDEFINED) {
            LOG_ERROR(ctx->ir, node->token, "Identifier ("STR_FMT") has an unresolved type", ident->name.len, ident->name.ptr);
        } else {
            node->flags = ident->node->flags;
            if (ident->data) {
                node->data = *ident->data;
            }
            if (!node->children) {
                md_array_push(node->children, ident->node, ctx->alloc);
            }
            return true;
        }
    } else {
        LOG_ERROR(ctx->ir, node->token, "Unresolved reference to identifier ("STR_FMT")", node->ident.len, node->ident.ptr);
    }
    return false;
}

static bool static_check_assignment(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ASSIGNMENT && node->children && md_array_size(node->children) == 2);
    ast_node_t* lhs = node->children[0];
    ast_node_t* rhs = node->children[1];

    ASSERT(lhs);
    ASSERT(rhs);

    int num_idents = 1;
    ast_node_t** idents = &lhs;

    if (lhs->type == AST_ARRAY) {
        for (size_t i = 0; i < md_array_size(lhs->children); ++i) {
			if (lhs->children[i]->type != AST_IDENTIFIER) {
				LOG_ERROR(ctx->ir, node->token, "An element on the left hand side of assignment is not an identifier");
				return false;
			}
		}
        num_idents = (int)md_array_size(lhs->children);
        idents = lhs->children;
    } else if (lhs->type != AST_IDENTIFIER) {
        LOG_ERROR(ctx->ir, node->token, "Left hand side of assignment is not an identifier");
		return false;
    }

    for (int i = 0; i < num_idents; ++i) {
        str_t ident = idents[i]->ident;
        // Assignment, therefore a new identifier
        if (find_constant(ident)) {
            LOG_ERROR(ctx->ir, idents[i]->token, "The identifier is occupied by a constant and cannot be assigned.");
        }
        else if (get_identifier(ctx->ir, ident)) {
            LOG_ERROR(ctx->ir, idents[i]->token, "The identifier is already taken. Variables cannot be reassigned.");
        }
    }

    if (static_check_node(rhs, ctx)) {
        ASSERT(rhs->data.type.base_type != TYPE_UNDEFINED);  // Right hand side type must be known
        node->data   = rhs->data;
        node->flags  = rhs->flags;

        if (num_idents > 1) {
            int rhs_len = 0;
            if (rhs->type == AST_ARRAY) {
                rhs_len = (int)md_array_size(rhs->children);
            } else {
                rhs_len = (int)type_info_array_len(rhs->data.type);
            }
            if (rhs_len == -1) {
                LOG_ERROR(ctx->ir, node->token, "Assignment mismatch between left and right hand side: The right hand side has a variable length");
                return false;
            }

            // Sometimes the dimensions of the right had side will be given in normalized form, i.e. [1][3][0][0]
            // In such cases, we still want to match assign this if the left hand side is provided as three identifiers.
            // Therefore we shift the dimensions of the rhs to fit.
            if (rhs->data.type.dim[0] == 1) {
                dim_shift_left(rhs->data.type.dim);
                rhs_len = (int)type_info_array_len(rhs->data.type);
            }

            // This is valid if the array length of rhs is equal to the number of identifiers on the left hand side
            if (rhs_len != num_idents) {
                LOG_ERROR(ctx->ir, node->token, "Assignment mismatch between left and right hand side: The number of identifiers on the left hand side must match the length of the right hand side");
                return false;
            }

            for (int i = 0; i < num_idents; ++i) {
                ASSERT(idents[i]->data.type.base_type == TYPE_UNDEFINED);  // Identifiers type should always be undefined until explicitly assigned.
                identifier_t* ident = create_identifier(ctx->ir, idents[i]->ident);
                if (!ident) {
                    LOG_ERROR(ctx->ir, node->token, "Failed to create identifier. Is the identifier already taken?");
                    return false;
                }

                ident->data = 0;
                ident->node = rhs;

                if (rhs->type == AST_ARRAY) {
                    // Match identifier with the corresponding element in the array
                    ident->data = &rhs->children[i]->data;
                    ident->node = rhs->children[i];
                }
                else {
                    idents[i]->data.type = type_info_element_type(rhs->data.type);
                    const int stride = (int)type_info_element_byte_stride(rhs->data.type);
                    ident->data = md_alloc(ctx->ir->arena, sizeof(data_t));
                    *ident->data = rhs->data;
                    ident->data->type = type_info_element_type(rhs->data.type);
                    dim_prune_leading_ones(ident->data->type.dim);
                    ident->data->size = stride;
                    ident->data->ptr = 0;
                }
                if (!static_check_node(idents[i], ctx)) {
                    return false;
                }
            }
        } else {
            ASSERT(lhs->data.type.base_type == TYPE_UNDEFINED);  // Identifiers type should always be undefined until explicitly assigned.

            identifier_t* ident = create_identifier(ctx->ir, lhs->ident);
            if (!ident) {
                LOG_ERROR(ctx->ir, node->token, "Failed to create identifier. Is the identifier already taken?");
                return false;
            }

            ident->data  = &rhs->data;
            ident->node  = rhs;

            if (!static_check_node(lhs, ctx)) {
                return false;
            }
        }
        lhs->data    = rhs->data;
        lhs->flags   = rhs->flags;
        return true;
    }

    return false;
}

static void propagate_contexts(ast_node_t* node, const md_bitfield_t* contexts, int64_t num_contexts) {
    ASSERT(node);
    node->num_contexts = num_contexts;
    node->contexts = contexts;
    for (size_t i = 0; i < md_array_size(node->children); ++i) {
        propagate_contexts(node->children[i], contexts, num_contexts);
    }
}

static bool static_check_out(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(node->type == AST_OUT);
    ASSERT(md_array_size(node->children) == 1);

    ast_node_t* expr = node->children[0];
    ASSERT(expr);

    const md_bitfield_t* old_ctx = ctx->mol_ctx;
    ctx->mol_ctx = 0;
    bool result = false;
    if (static_check_node(expr, ctx)) {
        node->flags |= expr->flags & FLAG_AST_PROPAGATION_MASK;
        node->data.type = expr->data.type;
        node->data.unit = expr->data.unit;
        node->data.value_range = expr->data.value_range;
        result = true;
    }
    ctx->mol_ctx = old_ctx;

    return result;
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

    // @NOTE: We explicitly dissallow any flattening for type checking of the Right hand side
    // Also we need to allow the length check to occur in order to determine the size of the context (rhs)
    uint32_t backup_flags = ctx->eval_flags;
    ctx->eval_flags &= ~(EVAL_FLAG_FLATTEN | EVAL_FLAG_NO_LENGTH_CHECK);
    bool result = static_check_node(rhs, ctx);
    ctx->eval_flags = backup_flags;
    
    if (result && rhs->data.type.base_type == TYPE_BITFIELD) {
        if (!(rhs->flags & FLAG_DYNAMIC)) {
            const int ctx_len = type_info_array_len(rhs->data.type);
            if (ctx_len > 0) {
                const size_t num_contexts = (size_t)ctx_len;
                // Store this persistently and set this as a context for all child nodes
                md_bitfield_t* contexts = 0;

                if (rhs->flags & FLAG_CONSTANT) {
                    ASSERT(rhs->data.ptr);
                    ASSERT(rhs->data.size == num_contexts * sizeof(md_bitfield_t));
                    contexts = (md_bitfield_t*)rhs->data.ptr;
                } else {
                    md_array_resize(contexts, num_contexts, ctx->ir->arena);
                    ASSERT(md_array_size(contexts) == num_contexts);
                    for (size_t i = 0; i < num_contexts; ++i) {
                        md_bitfield_init(&contexts[i], ctx->ir->arena);
                    }
                    rhs->data.ptr = contexts;
                    rhs->data.size = num_contexts * sizeof(md_bitfield_t);
                    result = evaluate_node(&rhs->data, rhs, ctx);
                }

                //allocate_data(&rhs->data, ctx->ir->arena);
                if (result) {

                    // We differentiate here if the LHS is a bitfield or not
                    // if LHS is a bitfield of length 1, the resulting bitfield length is the same as RHS (M)
                    // if LHS is a bitfield of length N, the resulting bitfield length is (N*M), N may vary for each context.

                    const bool contexts_equivalent = are_bitfields_equivalent(contexts, num_contexts, ctx->mol); 

                    node->lhs_context_types = 0;
                    size_t arr_len = 0;
                    for (size_t i = 0; i < num_contexts; ++i) {
                        eval_context_t local_ctx = *ctx;
                        local_ctx.mol_ctx = &contexts[i];
                        local_ctx.eval_flags |= EVAL_FLAG_NO_STATIC_EVAL;

                        if (!static_check_node(lhs, &local_ctx)) {
                            // @NOTE: This is a bit wild wild west,
                            // We hijack the last error message created and expand it with information about the context
                            // This is not a good solution, but it is a solution for now.
                            md_log_token_t* last_error = md_array_last(ctx->ir->errors);
                            if (last_error) {
                                last_error->context = &contexts[i];
                            }
                            return false;
                        }

                        type_info_t local_type = lhs->data.type;

                        if (lhs->flags & FLAG_DYNAMIC_LENGTH) {
                            if (!finalize_type(&local_type, lhs, &local_ctx)) {
                                LOG_ERROR(ctx->ir, lhs->token, "Failed to deduce length of type which is required for determining type and size of context");
                                return false;
                            }
                        }

                        // Normalize the type and dimensions here for the length to be computed correctly
                        if (local_type.base_type != TYPE_BITFIELD) {
                            int ndim = dim_ndims(local_type.dim);
                            if (ndim == 1 || (ndim > 1 && local_type.dim[0] != 1)) {
                                dim_shift_right(local_type.dim);
                            }
                        }

                        int len = type_info_array_len(local_type);
                        // Some arrays will be zero and that is accepted
                        // We still have to keep it even though it does not contribute
                        // To keep arrays in sync.
                        if (len < 0) {
                            LOG_ERROR(ctx->ir, lhs->token, "The type of the left hand side of 'in' must have a static length");
                            return false;
                        }

                        md_array_push(node->lhs_context_types, local_type, ctx->ir->arena);
                        arr_len += (size_t)len;
                    }

                    node->flags |= lhs->flags & FLAG_AST_PROPAGATION_MASK;
                    node->data.type = node->lhs_context_types[0];
                    if (node->data.type.base_type == TYPE_BITFIELD) {
                        node->data.type.dim[0] *= (int)arr_len;
                        // This is a flattening of the dimensions
                        int len = dim_size(node->data.type.dim);
                        MEMSET(node->data.type.dim, 0, sizeof(node->data.type.dim));
                        node->data.type.dim[0] = len;
                    } else {
                        node->data.type.dim[0] = (int)arr_len;
                    }
                    node->data.unit = lhs->data.unit;
                    node->data.value_range = lhs->data.value_range;

                    propagate_contexts(lhs, contexts, num_contexts);
                    if (contexts_equivalent) {
                        lhs->flags |= FLAG_CONTEXTS_EQUIVALENT;
                    }

                    result = true;

                } else {
                    LOG_ERROR(ctx->ir, node->token, "Right hand side of 'in' failed to evaluate at compile time.");
                }
            } else {
                LOG_ERROR(ctx->ir, node->token, "The context is empty.");
            }
        } else {
            LOG_ERROR(ctx->ir, node->token, "Right hand side of 'in' must be known at compile time.");
        }
    } else {
        char buf[128];
        print_type_info(buf, ARRAY_SIZE(buf), rhs->data.type);
        LOG_ERROR(ctx->ir, node->token, "Right hand side of keyword 'in' has an incompatible type: Expected bitfield, got '%s'.", buf);
    }

    return result;
}

static bool static_check_node(ast_node_t* node, eval_context_t* ctx) {
    // This is probably more like a static check which encompass more than just checking the types...
    // The idea here is that we have an syntax tree given by node and we want to map the function calls into concrete functions and check for the type.
    // The type will propagate back throughout the tree and we will remap nodes into concrete function calls.
    // This is not done in the parsing step as we want the full context before attempting this operation.

    ASSERT(node);

    bool result = false;

    switch (node->type) {
    case AST_ASSIGNMENT:
        result = static_check_assignment(node, ctx);
        break;
    case AST_CONSTANT_VALUE:
        result = static_check_constant_value(node, ctx);
        break;
    case AST_ARRAY:
        result = static_check_array(node, ctx);
        break;
    case AST_ARRAY_SUBSCRIPT:
        result = static_check_array_subscript(node, ctx);
        break;
    case AST_IDENTIFIER:
        result = static_check_identifier_reference(node, ctx);
        break;
    case AST_TABLE:
        result = static_check_import(node, ctx);
        break;
    case AST_FLATTEN:
        result = static_check_flatten(node, ctx);
        break;
    case AST_TRANSPOSE:
        result = static_check_transpose(node, ctx);
        break;
    case AST_PROC_CALL:
        result = static_check_proc_call(node, ctx);
        break;
    case AST_CONTEXT:
        result = static_check_context(node, ctx);
        break;
    case AST_OUT:
        result = static_check_out(node, ctx);
        break;
    case AST_ADD:
    case AST_SUB:
    case AST_MUL:
    case AST_DIV:
    case AST_AND:
    case AST_XOR:
    case AST_OR:
    case AST_NOT:
    case AST_EQ:
    case AST_NE:
    case AST_LE:
    case AST_GE:
    case AST_LT:
    case AST_GT:
    case AST_UNARY_NEG:
        result = static_check_operator(node, ctx);
        break;
    case AST_CAST: // Should never happen since we insert casts here in type checking phase.
    default:
        ASSERT(false);
    }

    // Prune trailing ones in dimensions
    for (int i = MAX_NUM_DIMS - 1; i > 0; --i) {
        if (node->data.type.dim[i] == 1) {
            node->data.type.dim[i] = 0;
        }
    }

    if (result && !(node->flags & FLAG_DYNAMIC) && !(ctx->eval_flags & EVAL_FLAG_NO_STATIC_EVAL)) {
        static_eval_node(node, ctx);
    }

    return result;
}

static uint64_t hash_node(const ast_node_t*, uint64_t);

static uint64_t hash_children(const ast_node_t* node, uint64_t seed) {
    uint64_t hash = seed;
    for (size_t i = 0; i < md_array_size(node->children); ++i) {
        hash = hash_node(node->children[i], seed);
    }
    return hash;
}

static uint64_t hash_table(const table_t* table, uint64_t seed) {
    // Tables are resolved at compile time, thus we can check their contents directly
    uint64_t hash = seed;
    for (size_t i = 0; i < table->num_fields; ++i) {
        hash = md_hash64(table->field_names[i].ptr, table->field_names[i].len, hash);
        hash = md_hash64(table->field_values[i], sizeof(double) * table->num_values, hash);
    }
    return hash;
}

static uint64_t hash_node(const ast_node_t* node, uint64_t seed) {
    ASSERT(node);
    switch (node->type) {
    case AST_CONSTANT_VALUE:
        switch (node->data.type.base_type) {
        case TYPE_BITFIELD:
        {
            const md_bitfield_t* bf = (const md_bitfield_t*)node->data.ptr;
            return md_bitfield_hash64(bf, seed);
        }
        case TYPE_STRING:
        {
            const str_t* str = (const str_t*)node->data.ptr;
            return md_hash64(str->ptr, str->len, seed);
        }
        default:
            return md_hash64(node->data.ptr, node->data.size, seed);
        }
    case AST_IDENTIFIER:
        return md_hash64(node->ident.ptr, node->ident.len, md_hash64(&node->data.type, sizeof(node->data.type), seed));
    case AST_PROC_CALL:
        return md_hash64(node->ident.ptr, node->ident.len, md_hash64(&node->data.type, sizeof(node->data.type), hash_children(node, seed)));
    case AST_TABLE:
        return hash_table(node->table, md_hash64(node->table_field_indices, md_array_bytes(node->table_field_indices), seed));
    case AST_ARRAY:
    case AST_ARRAY_SUBSCRIPT:
    default:
        return md_hash64(&node->data.type, sizeof(node->data.type), hash_children(node, seed));
    }
}

// Prunes the tree from AST_EXPRESSION, since that is just an unnecessary indirection for sorting the tree correctly on precedence
static ast_node_t* prune_expressions(ast_node_t* node) {
    if (node) {
        for (size_t i = 0; i < md_array_size(node->children); ++i) {
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

static bool parse_script(md_script_ir_t* ir) {
    ASSERT(ir);

    bool result = true;

    tokenizer_t tokenizer = tokenizer_init(ir->str);

    parse_context_t ctx = {
        .ir = ir,
        .tokenizer = &tokenizer,
        .node = 0,
        .temp_alloc = md_temp_allocator,
    };

    ir->stage = "Parsing";

    // Parse statements until we have exhausted the tokenizer

    token_t tok = {0};
    while (tok = tokenizer_peek_next(ctx.tokenizer), tok.type != TOKEN_END) {
        ctx.node = 0;
        ir->record_log = true;   // We reset the error recording flag before each statement
        ast_node_t* raw_node = parse_expression(&ctx);
        ast_node_t* pruned_node = prune_expressions(raw_node);
        if (pruned_node) {
            tok = tokenizer_consume_next(ctx.tokenizer);
            if (tok.type != ';') {
                LOG_ERROR(ir, pruned_node->token, "Missing ';' to mark end of statement.");
                result = false;
                continue;
            }
            md_array_push(ir->parsed_expressions, pruned_node, ir->arena);
        } else {
            token_type_t types[1] = {';'};
            tokenizer_consume_until_type(ctx.tokenizer, types, ARRAY_SIZE(types)); // Goto next statement
            tokenizer_consume_next(ctx.tokenizer); 
            result = false;
        }
    }

    return result;
}

static bool static_type_check(md_script_ir_t* ir, const md_molecule_t* mol, const md_trajectory_i* traj) {
    ASSERT(ir);
    ASSERT(mol);

    bool result = true;

    SETUP_TEMP_ALLOC(GIGABYTES(4));

    eval_context_t ctx = {
        .ir = ir,
        .temp_arena = &vm_arena,
        .temp_alloc = &temp_alloc,
        .alloc = ir->arena,
        .mol = mol,
        .traj = traj,
    };

    ir->stage = "Static Type Checking";

    for (size_t i = 0; i < md_array_size(ir->parsed_expressions); ++i) {
        ir->record_log = true;   // We reset the error recording flag before each statement
        if (static_check_node(ir->parsed_expressions[i], &ctx)) {
            md_array_push(ir->type_checked_expressions, ir->parsed_expressions[i], ir->arena);
            ir->flags |= (ir->parsed_expressions[i]->flags & FLAG_IR_PROPAGATION_MASK);
        } else {
            MD_LOG_DEBUG("Static type checking failed for expression: '"STR_FMT"'", STR_ARG(ir->parsed_expressions[i]->token.str));
            result = false;
        }
    }

    FREE_TEMP_ALLOC();

    return result;
}

static bool extract_dynamic_evaluation_targets(md_script_ir_t* ir) {
    ASSERT(ir);

    // Check all type checked expressions for dynamic flag
    for (size_t i = 0; i < md_array_size(ir->type_checked_expressions); ++i) {
        ast_node_t* expr = ir->type_checked_expressions[i];
        ASSERT(expr);
        if (expr->flags & FLAG_DYNAMIC && expr->type == AST_ASSIGNMENT) {
            // If it does not have an identifier, it cannot be referenced and therefore we don't need to evaluate it dynamcally
            md_array_push(ir->eval_targets, ir->type_checked_expressions[i], ir->arena);
        }
    }
    return true;
}

static inline bool is_temporal_type(type_info_t ti) {
    const int ndim = dim_ndims(ti.dim);
    return ((ti.dim[0] != -1) && (ti.base_type == TYPE_FLOAT) && (ndim == 1)) || (ndim == 2 && ti.dim[1] == 1);
}

static inline bool is_distribution_type(type_info_t ti) {
    return ti.dim[0] != -1 && is_type_directly_compatible(ti, (type_info_t)TI_DISTRIBUTION);
}

static inline bool is_volume_type(type_info_t ti) {
    return ti.dim[0] != -1 && is_type_directly_compatible(ti, (type_info_t)TI_VOLUME);
}

static inline bool is_property_type(type_info_t ti) {
    return is_temporal_type(ti) || is_distribution_type(ti) || is_volume_type(ti);
}

static bool static_eval_node(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(ctx);
    ASSERT(ctx->ir);

    const size_t num_children = md_array_size(node->children);
    
    // Propagate the evaluation all the way down to the children.
    if (node->children) {
        size_t offset = 0;
        if (node->type == AST_CONTEXT) offset = 1;  // We don't want to evaluate the lhs in a context on its own since then it will loose its context
        for (size_t i = offset; i < num_children; ++i) {
            if (!static_eval_node(node->children[i], ctx)) {
                return false;
            }
        }
    }

    // Only evaluate the node if it is not flagged as dynamic
    // Only evaluate if data.ptr is not already set (which can happen during static check)
    if (!(node->flags & FLAG_DYNAMIC) && !node->data.ptr) {
        uint64_t hash = hash_node(node, 0);

        // Try to find in existing expressions
        const size_t num_expressions = md_array_size(ctx->ir->static_expression_hash);
        for (size_t i = 0; i < num_expressions; ++i) {
            if (hash == ctx->ir->static_expression_hash[i]) {
                // str_t str = ctx->ir->static_expression_str[i];
                node->data = ctx->ir->static_expression_data[i];
                node->flags |= FLAG_CONSTANT;
                return true;
            }
        }

        data_t data = {0};
        if (allocate_data(&data, node->data.type, ctx->ir->arena)) {
            if (!evaluate_node(&data, node, ctx)) {
                LOG_ERROR(ctx->ir, node->token, "Failed to evaluate node during static evaluation");
                return false;
            }
            if (node->type == AST_CONTEXT) {
                ASSERT(node->children);
                // If its a context node, we copy the data to the child as well
                node->children[0]->data = data;
            }
            node->data = data;
            node->flags |= FLAG_CONSTANT;

            md_array_push(ctx->ir->static_expression_hash, hash, ctx->alloc);
            md_array_push(ctx->ir->static_expression_data, data, ctx->alloc);
            md_array_push(ctx->ir->static_expression_str, node->token.str, ctx->alloc);

            return true;
        } else {
            LOG_ERROR(ctx->ir, node->token, "Could not allocate data for node during static evaluation");
            return false;
        }
    }

    return true;
}

static bool static_evaluation(md_script_ir_t* ir, const md_molecule_t* mol) {
    ASSERT(mol);

    SETUP_TEMP_ALLOC(GIGABYTES(4));

    eval_context_t ctx = {
        .ir = ir,
        .mol = mol,
        .temp_arena = &vm_arena,
        .temp_alloc = &temp_alloc,
    };

    ir->stage = "Static Evaluation";
    bool result = true;

    // Evaluate every node which is not flagged with FLAG_DYNAMIC and store its value.
    for (size_t i = 0; i < md_array_size(ir->type_checked_expressions); ++i) {
        result &= static_eval_node(ir->type_checked_expressions[i], &ctx);
        md_vm_arena_reset(&vm_arena);
    }

    FREE_TEMP_ALLOC();

    return result;
}

 static void allocate_property_data(md_script_property_data_t* data, md_script_property_flags_t flags, type_info_t type, size_t num_frames, md_allocator_i* alloc) {
    ASSERT(data);    
    ASSERT(alloc);

    // @NOTE: We need to 'normalize' the dimensionality of the types in a consistent way, since the properties will be exposed
    // Therefore we make sure all types encode the array length in the first dimension, even if its a scalar, i.e. [1]
    // This simplifies extraction later on.

    // Step 1: Prune trailing ones
    for (int i = MAX_NUM_DIMS - 1; i > 0; --i) {
        if (type.dim[i] == 1) {
            type.dim[i] = 0;
        }
    }

    // @NOTE: At this stage, the flags should only contain the property type
    switch (flags) {
        case MD_SCRIPT_PROPERTY_FLAG_TEMPORAL:
            // For temporal data, we store and expose all values, this enables filtering to be performed afterwards to create distributions
            if (dim_ndims(type.dim) < 2) {
				dim_shift_right(type.dim);
			}
            data->dim[0] = (int32_t)num_frames;
            data->dim[1] = type.dim[1];
            break;
        case MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION:
            if (dim_ndims(type.dim) < 3) {
                dim_shift_right(type.dim);
            }
            MEMCPY(data->dim, type.dim, sizeof(type.dim));
            // @NOTE: Multidimensional distributions are not supported yet
            ASSERT(data->dim[0] == 1);
            break;
        case MD_SCRIPT_PROPERTY_FLAG_VOLUME:
            if (dim_ndims(type.dim) < 4) {
            	dim_shift_right(type.dim);
            }
            MEMCPY(data->dim, type.dim, sizeof(type.dim));
            // @NOTE: Multidimensional volumes are not supported yet
            ASSERT(data->dim[0] == 1);
            break;
        default:
            ASSERT(false);
            break;
    }

    const size_t num_values = (size_t)dim_size(data->dim);
    const size_t num_bytes  = num_values * sizeof(float);
    data->values = md_alloc(alloc, num_bytes);
    MEMSET(data->values, 0, num_bytes);
    data->num_values = num_values;

    if (flags == MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) {
        data->weights = data->values + data->dim[2];
        for (size_t i = 0; i < num_values; ++i) {
            data->weights[i] = 1.0f;
        }
    }

    if (flags == MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
        int ndim = dim_ndims(data->dim);
        if (data->dim[ndim - 1] > 1) {
            // Allocate data for aggregate
            const size_t aggregate_size = num_frames;
            data->aggregate = md_alloc(alloc, sizeof(md_script_aggregate_t));

            MEMSET(data->aggregate, 0, sizeof(md_script_aggregate_t));
            data->aggregate->num_values = aggregate_size;

            if (aggregate_size != num_values) {
                data->aggregate->population_mean = md_alloc(alloc, aggregate_size * sizeof(float));
                MEMSET(data->aggregate->population_mean, 0, aggregate_size * sizeof(float));
            } else {
                data->aggregate->population_mean = data->values;
            }

            data->aggregate->population_var = md_alloc(alloc, aggregate_size * sizeof(float));
            MEMSET(data->aggregate->population_var, 0, aggregate_size * sizeof(float));

            data->aggregate->population_ext = md_alloc(alloc, aggregate_size * sizeof(vec2_t));
            MEMSET(data->aggregate->population_ext, 0, aggregate_size * sizeof(vec2_t));
        }
    }
}

static void compute_min_max_mean_variance(float* out_min, float* out_max, float* out_mean, float* out_var, const float* data, size_t count) {
    ASSERT(out_min);
    ASSERT(out_max);
    ASSERT(out_mean);
    ASSERT(out_var);
    ASSERT(data);

    const float N = (float)count;
    float min = FLT_MAX;
    float max = -FLT_MAX;
    float s1 = 0;
    float s2 = 0;
    size_t i = 0;

    for (i = 0; i < count; ++i) {
        s1 += data[i];
        min = MIN(min, data[i]);
        max = MAX(max, data[i]);
    }

    s1 = s1 / N;

    for (i = 0; i < count; ++i) {
        s2 += (data[i] - s1) * (data[i] - s1);
    }
    s2 = s2 / N;

    *out_min = min;
    *out_max = max;
    *out_mean = s1;
    *out_var = s2;

    /*
    if (count > md_simd_width) {
        md_256 vmin = md_mm256_set1_ps(FLT_MAX);
        md_256 vmax = md_mm256_set1_ps(-FLT_MAX);
        md_256 v1 = md_simd_zero_f32();
        md_256 v2 = md_simd_zero_f32();

        const int simd_count = (count / md_simd_width) * md_simd_width;
        for (; i < simd_count; i += md_simd_width) {
            md_256 val = md_mm256_loadu_ps(data + i);
            vmin = md_simd_min_f32(vmin, val);
            vmax = md_simd_min_f32(vmax, val);
            v1 = md_mm256_add_ps_f32(v1, val);
            v2 = md_mm256_add_ps_f32(v2, md_mm256_mul_ps_f32(val, val));
        }

        s1 = md_simd_hadd_f32(v1);
        s2 = md_simd_hadd_f32(v2);
        min = md_simd_reduce_min_f32(vmin);
        max = md_simd_reduce_max_f32(vmax);
    }

    for (; i < count; ++i) {
        float val = data[i];
        min = MIN(min, val);
        max = MAX(max, val);
        s1 += val;
        s2 += val * val;
    }

    *out_min  = min;
    *out_max  = max;
    *out_mean = s1 / N;
    *out_var  = fabsf((N * s2) - (s1 * s1)) / (N*N);
    * */
}

static void clear_property_data(md_script_property_data_t* data) {
    MEMSET(data->values, 0, data->num_values * sizeof(float));
    if (data->aggregate) {
        MEMSET(data->aggregate->population_mean, 0, data->aggregate->num_values * sizeof(float));
        MEMSET(data->aggregate->population_var,  0, data->aggregate->num_values * sizeof(float));
        MEMSET(data->aggregate->population_ext,  0, data->aggregate->num_values * sizeof(vec2_t));
        /*
        MEMSET(data->aggregate->population_min,  0, data->aggregate->num_values * sizeof(float));
        MEMSET(data->aggregate->population_max,  0, data->aggregate->num_values * sizeof(float));
        */
    }
    data->min_value = +FLT_MAX;
    data->max_value = -FLT_MAX;
}

static bool eval_properties(md_script_eval_t* eval, const md_molecule_t* mol, const md_trajectory_i* traj, const md_script_ir_t* ir, uint32_t frame_beg, uint32_t frame_end) {
    ASSERT(eval);
    ASSERT(mol);
    ASSERT(traj);
    ASSERT(ir);

    // No properties to evaluate!
    if (!eval->property_data) return true;
    
    const size_t num_expr = md_array_size(ir->eval_targets);
    ast_node_t** const expr = ir->eval_targets;
    
    SETUP_TEMP_ALLOC(GIGABYTES(4));

    // coordinate data for reading trajectory frames into
    const size_t stride = ALIGN_TO(mol->atom.count, 8);    // Round up allocation size to simd width to allow for vectorized operations
    const size_t coord_bytes = stride * 3 * sizeof(float);
    float* init_coords = md_vm_arena_push(&vm_arena, coord_bytes);
    float* curr_coords = md_vm_arena_push(&vm_arena, coord_bytes);
    
    // This data is meant to hold the evaluated expressions
    data_t* data = md_vm_arena_push(&vm_arena, num_expr * sizeof(data_t));

    const size_t STACK_RESET_POINT = md_vm_arena_get_pos(&vm_arena);

    //int thread_id = (md_thread_id() & 0xFFFF);
    //md_logf(MD_LOG_TYPE_DEBUG, "Starting evaluation on thread %i, range (%i,%i) arena size: %.2f MB", thread_id, (int)frame_beg, (int)frame_end, (double)vm_arena.commit_pos / (double)MEGABYTES(1));
    //uint64_t max_arena_pos = 0;

    md_trajectory_frame_header_t init_header = { 0 };
    md_trajectory_frame_header_t curr_header = { 0 };

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
        .ir = (md_script_ir_t*)ir,  // We cast away the const here. The evaluation will not modify ir.
        .mol = &mutable_mol,
        .traj = traj,
        .temp_arena = &vm_arena,
        .temp_alloc = &temp_alloc,
        .alloc = &temp_alloc,
        .frame_header = &curr_header,
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

    // We evaluate each frame, one at a time
    for (uint32_t f_idx = frame_beg; f_idx < frame_end; ++f_idx) {
        if (eval->interrupt) {
            goto done;
        }
        
        result = md_trajectory_load_frame(traj, f_idx, &curr_header, curr_x, curr_y, curr_z);

        if (!result) {
            MD_LOG_ERROR("Failed to load frame during evaluation");
            goto done;
        }

        MEMCPY(&mutable_mol.unit_cell, &curr_header.unit_cell, sizeof(md_unit_cell_t));
        
        md_vm_arena_set_pos(&vm_arena, STACK_RESET_POINT);
        ctx.identifiers = 0;
        ctx.spatial_hash = 0;

        for (size_t i = 0; i < num_expr; ++i) {
            type_info_t type = expr[i]->data.type;
            if (is_variable_length(type)) {
                if (!finalize_type(&type, expr[i], &ctx)) {
                    str_t str = expr[i]->token.str;
                    MD_LOG_ERROR("Evaluation error when evaluating the following expression '"STR_FMT"', failed to finalize its type", STR_ARG(str));
                    result = false;
                    goto done;
                }
            }
            allocate_data(&data[i], type, &temp_alloc);
            data[i].unit = expr[i]->data.unit;
            data[i].value_range = expr[i]->data.value_range;
            if (data[i].value_range.beg == 0 && data[i].value_range.end == 0) {
                data[i].value_range.beg = -FLT_MAX;
                data[i].value_range.end = +FLT_MAX;
            }
            if (!evaluate_node(&data[i], expr[i], &ctx)) {
                str_t str = expr[i]->token.str;
                MD_LOG_ERROR("Evaluation error when evaluating the following expression '"STR_FMT"' at frame %i", STR_ARG(str), (int)f_idx);
                result = false;
                goto done;
            }
        }

        size_t num_props = md_array_size(eval->property_data);
        for (size_t p_idx = 0; p_idx < num_props; ++p_idx) {
            ASSERT(str_eq(eval->property_names[p_idx], ir->property_names[p_idx]));

            str_t p_name = eval->property_names[p_idx];
            md_script_property_flags_t p_flags = ir->property_flags[p_idx];
            md_script_property_data_t* p_data = &eval->property_data[p_idx];
            // Find data matching property identifier
            const identifier_t* ident = find_dynamic_identifier(p_name, &ctx);
            if (!ident) {
                MD_LOG_ERROR("Not good!");
                result = false;
                goto done;
            }

            const frange_t value_range = ident->data->value_range;
            const float* values = (const float*)ident->data->ptr;

            if (p_flags & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
                ASSERT(ident->data);
                ASSERT(p_data->values);

                const int32_t size = ident->data->type.dim[0];
                MEMCPY(p_data->values + f_idx * size, values, size * sizeof(float));

                // Determine min max values
                float min, max, mean, var;
                compute_min_max_mean_variance(&min, &max, &mean, &var, values, size);
                p_data->min_value = MIN(p_data->min_value, min);
                p_data->max_value = MAX(p_data->max_value, max);

                if (p_data->aggregate) {
                    p_data->aggregate->population_mean[f_idx] = mean;
                    p_data->aggregate->population_var [f_idx] = var;
                    p_data->aggregate->population_ext [f_idx] = (vec2_t){min, max};
                    /*
                    p_data->aggregate->population_min [f_idx] = min;
                    p_data->aggregate->population_max [f_idx] = max;
                    */
                }

                // Update range if not explicitly set
                p_data->min_range[0] = (value_range.beg == -FLT_MAX) ? p_data->min_value : value_range.beg;
                p_data->max_range[0] = (value_range.end == +FLT_MAX) ? p_data->max_value : value_range.end;
            }
            else if (p_flags & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) {
                // Accumulate values
                ASSERT(p_data->values);

                size_t num_bins = p_data->dim[2];
                float min, max, mean, var;
                compute_min_max_mean_variance(&min, &max, &mean, &var, values, num_bins);

                // ATOMIC WRITE
                md_mutex_lock(&eval->property_dist_mutex[p_idx]);
                {
                    // Cumulative moving average
                    const uint32_t count = eval->property_dist_count[p_idx]++;
                    const md_256 N   = md_mm256_set1_ps((float)(count));
                    const md_256 scl = md_mm256_set1_ps(1.0f / (float)(count + 1));

                    const int64_t length = ALIGN_TO(num_bins, 8);
                    for (int64_t i = 0; i < length; i += 8) {
                        md_256 old_val = md_mm256_mul_ps(md_mm256_loadu_ps(p_data->values + i), N);
                        md_256 new_val = md_mm256_loadu_ps(values + i);
                        md_mm256_storeu_ps(p_data->values + i, md_mm256_mul_ps(md_mm256_add_ps(new_val, old_val), scl));
                    }

                    // Copy weights
                    MEMCPY(p_data->weights, values + num_bins, num_bins * sizeof(float));
                    
                    // Determine min max values
                    p_data->min_value = MIN(p_data->min_value, min);
                    p_data->max_value = MAX(p_data->max_value, max);

                    // Update range if not explicitly set
                    p_data->min_range[0] = (value_range.beg == -FLT_MAX) ? p_data->min_value : value_range.beg;
                    p_data->max_range[0] = (value_range.end == +FLT_MAX) ? p_data->max_value : value_range.end;
                }
                md_mutex_unlock(&eval->property_dist_mutex[p_idx]);
            }
            else if (p_flags & MD_SCRIPT_PROPERTY_FLAG_VOLUME) {
                // Accumulate values
                ASSERT(p_data->values);
                ASSERT(p_data->num_values % 8 == 0); // This should always be the case if we nice powers of 2 for our volumes
                
                // ATOMIC WRITE
                md_mutex_lock(&eval->property_dist_mutex[p_idx]);
                {
                    // Cumulative moving average
                    const uint32_t count = eval->property_dist_count[p_idx]++;
                    const md_256 N = md_mm256_set1_ps((float)(count));
                    const md_256 scl = md_mm256_set1_ps(1.0f / (float)(count + 1));

                    for (size_t i = 0; i < p_data->num_values; i += 8) {
                        md_256 old_val = md_mm256_mul_ps(md_mm256_loadu_ps(p_data->values + i), N);
                        md_256 new_val = md_mm256_loadu_ps(values + i);
                        md_mm256_storeu_ps(p_data->values + i, md_mm256_mul_ps(md_mm256_add_ps(new_val, old_val), scl));
                    }
                }
                md_mutex_unlock(&eval->property_dist_mutex[p_idx]);
            }
            else {
                ASSERT(false);
            }
        }

        md_mutex_lock(&eval->frame_lock);
        md_bitfield_set_bit(&eval->frame_mask, f_idx);
        md_mutex_unlock(&eval->frame_lock);
        
        //max_arena_pos = MAX(max_arena_pos, vm_arena.commit_pos);
    }
done:
    //md_logf(MD_LOG_TYPE_DEBUG, "Finished evaluation on thread %i, max arena size: %.2f MB", thread_id, (double)max_arena_pos / (double)MEGABYTES(1));
    FREE_TEMP_ALLOC();
    return result;
}

static inline bool validate_ir(const md_script_ir_t* ir) {
    if (!ir) {
        return false;
    }

    if (ir->magic != SCRIPT_IR_MAGIC) {
        return false;
    }

    return true;
}

static inline bool validate_eval(const md_script_eval_t* eval) {
    if (!eval) {
        return false;
    }

    if (eval->magic != SCRIPT_EVAL_MAGIC) {
        return false;
    }

    return true;
}

static void reset_ir(md_script_ir_t* ir) {
    md_allocator_i* arena = ir->arena;
    md_arena_allocator_reset(ir->arena);
    MEMSET(ir, 0, sizeof(md_script_ir_t));
    ir->magic = SCRIPT_IR_MAGIC;
    ir->arena = arena;
}

static md_script_ir_t* create_ir(md_allocator_i* alloc) {
    md_script_ir_t* ir = md_alloc(alloc, sizeof(md_script_ir_t));
    md_allocator_i* arena = md_arena_allocator_create(alloc, MEGABYTES(1));
    MEMSET(ir, 0, sizeof(md_script_ir_t));
    ir->magic = SCRIPT_IR_MAGIC;
    ir->arena = arena;
    return ir;
}

static bool add_ir_ctx(md_script_ir_t* ir, const md_script_ir_t* ctx_ir) {
    ASSERT(ir);
    ASSERT(ctx_ir);

    // we want to extract identifiers and subexpressions to be used so we can reference it.
    // @TODO: Check for collisions of identifiers here
    md_array_push_array(ir->identifiers,  ctx_ir->identifiers,  md_array_size(ctx_ir->identifiers),  ir->arena);
    //md_array_push_array(ir->eval_targets, ctx_ir->eval_targets, md_array_size(ctx_ir->eval_targets), ir->arena);
    ir->flags |= ctx_ir->flags;

    return true;
}

static void free_ir(md_script_ir_t* ir) {
    if (validate_ir(ir)) {
        md_arena_allocator_destroy(ir->arena);
        memset(ir, 0, sizeof(md_script_ir_t));
    }
}

static void create_vis_tokens(md_script_ir_t* ir, const ast_node_t* node, const ast_node_t* node_override, int32_t depth) {
    md_script_vis_token_t vis = {0};

    ASSERT(node->type != AST_UNDEFINED);

    md_strb_t sb = {0};
    md_strb_init(&sb, md_temp_allocator);

    char type_buf[128];
    size_t type_len = print_type_info(type_buf, (int)sizeof(type_buf), node->data.type);
    md_strb_push_cstrl(&sb, type_buf, type_len);

    char unit_buf[128];
    int unit_len = (int)md_unit_print(unit_buf, (int)sizeof(unit_buf), node->data.unit);
    if (unit_len) {
        md_strb_fmt(&sb, " ("STR_FMT")", unit_len, unit_buf);
    }

    if (node->type == AST_TABLE) {
        // Write header contents of table to provide an overview of what is there
        md_strb_push_char(&sb, '\n');
        bool matching_units = true;
        for (size_t i = 0; i < md_array_size(node->table_field_indices); ++i) {
            int idx = node->table_field_indices[i];
            if (idx < 0 || (int)node->table->num_fields <= idx) {
                MD_LOG_DEBUG("Attempting to read out of bounds in table_field_indices");
                continue;
            }
            str_t name = node->table->field_names[idx];
            md_strb_fmt(&sb, "[%i]: \""STR_FMT"\"", (int)(i + 1), STR_ARG(name));
            md_unit_t unit = node->table->field_units[idx];
            if (!md_unit_empty(unit) && !md_unit_unitless(unit)) {
                str_t unit_str = md_unit_to_string(unit, md_temp_allocator);
                if (!str_empty(unit_str)) {
                    md_strb_fmt(&sb, " ("STR_FMT")", STR_ARG(unit_str));
                }
            }
            if (!md_unit_equal(node->data.unit, unit)) {
                matching_units = false;
            }
            md_strb_push_char(&sb, '\n');
        }
        if (!matching_units) {
            md_strb_fmt(&sb, "Cannot assign unit to data due to conflicting units in fields.");
        }
    } else if (!(node->flags & FLAG_DYNAMIC)) {
        if (node->data.type.base_type != TYPE_BITFIELD) {
            md_strb_push_char(&sb, '\n');
            char val_buf[128] = {0};
            size_t val_len = print_data_value(val_buf, sizeof(val_buf), node->data);
            md_strb_push_cstrl(&sb, val_buf, val_len);
        }
    } else {
        md_strb_push_str(&sb, STR_LIT(" [d]"));
    }

#if 0
    if (node->data.size) {
        if (node->data.size / MEGABYTES(1)) {
            md_strb_fmt(&sb, " [%.2fMB]", (double)node->data.size / (double)MEGABYTES(1));
        } else if (node->data.size / KILOBYTES(1)) {
            md_strb_fmt(&sb, " [%.2fKB]", (double)node->data.size / (double)KILOBYTES(1));
        } else {
            md_strb_fmt(&sb, " [%iB]", (int)node->data.size);
        }
    }
#endif

    vis.range.beg = node->token.beg;
    vis.range.end = node->token.end;
    vis.depth = depth;
    vis.text = str_copy(md_strb_to_str(sb), ir->arena);
    vis.payload = (const struct md_script_vis_payload_o*)(node_override ? node_override : node);

    // Parent marker should be last marker added (unless empty)
    md_script_vis_token_t* last = md_array_last(ir->vis_tokens);
    if (last &&
        vis.range.beg == last->range.beg &&
        vis.range.end == last->range.end)
    {
        // If the current marker has the exact same line, col_beg and col_end, we only want to update the text, not the payload (the node)
        // since we will evaluate that top down anyways.
        // This is to deal with 'casts' which have the same marker size, but shadow the actual operation we want to visualize
        last->text = vis.text;
    } else {
        md_array_push(ir->vis_tokens, vis, ir->arena);
    }

    // Recurse and add children
    const size_t num_children = md_array_size(node->children);
    if (node->type == AST_CONTEXT) {
        ASSERT(num_children == 2);
        create_vis_tokens(ir, node->children[0], node, depth + 1);
        create_vis_tokens(ir, node->children[1], NULL, depth + 1);
    } else {
        for (size_t i = 0; i < num_children; ++i) {
            create_vis_tokens(ir, node->children[i], NULL, depth + 1);
        }
    }
}

static bool extract_vis_tokens(md_script_ir_t* ir) {
    const size_t num_expr = md_array_size(ir->type_checked_expressions);
    for (size_t i = 0; i < num_expr; ++i) {
        ast_node_t* expr = ir->type_checked_expressions[i];
        create_vis_tokens(ir, expr, NULL, 0);
    }

    return true;
}

static bool extract_identifiers(md_script_ir_t* ir) {
    const size_t num_ident = md_array_size(ir->identifiers);
    for (size_t i = 0; i < num_ident; ++i) {
        md_array_push(ir->identifier_names, ir->identifiers->name, ir->arena);
    }
    return true;
}

static bool create_property(md_script_ir_t* ir, str_t ident, const ast_node_t* node) {
    md_script_property_flags_t flags = 0;

    if (is_temporal_type(node->data.type)) {
        flags |= MD_SCRIPT_PROPERTY_FLAG_TEMPORAL;
    }
    else if (is_distribution_type(node->data.type)) {
        flags |= MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION;
    }
    else if (is_volume_type(node->data.type)) {
        flags |= MD_SCRIPT_PROPERTY_FLAG_VOLUME;
    }

    md_array_push(ir->property_names, ident, ir->arena);
    md_array_push(ir->property_flags, flags, ir->arena);
    md_array_push(ir->property_nodes, node,  ir->arena);

    return true;
}

static bool extract_properties(md_script_ir_t* ir) {
    for (size_t i = 0; i < md_array_size(ir->eval_targets); ++i) {
        ast_node_t* expr = ir->eval_targets[i];
        ASSERT(expr);
        if (expr->type == AST_ASSIGNMENT) {
            ASSERT(expr->children);
            ast_node_t* lhs = expr->children[0];
            if (expr->children[0]->type == AST_IDENTIFIER && is_property_type(expr->data.type)) {
                create_property(ir, lhs->ident, expr);
            } else if (expr->children[0]->type == AST_ARRAY) {
                for (size_t j = 0; j < md_array_size(lhs->children); ++j) {
                    ast_node_t* child = lhs->children[j];
                    if (child->type == AST_IDENTIFIER && is_property_type(child->data.type)) {
                        create_property(ir, child->ident, child);
                    }
                }
            }
        }
    }
    return true;
}

md_script_ir_t* md_script_ir_create(md_allocator_i* alloc) {
    if (!alloc) {
        MD_LOG_ERROR("Script Create: Allocator was null");
        return NULL;
    }
    return create_ir(alloc);
}

bool md_script_ir_compile_from_source(md_script_ir_t* ir, str_t src, const md_molecule_t* mol, const md_trajectory_i* traj, const md_script_ir_t* ctx_ir) {
    if (!validate_ir(ir)) {
        MD_LOG_ERROR("Script Compile: IR is not valid");
        return false;
    }

    if (str_empty(src)) {
        MD_LOG_ERROR("Script Compile: Source string was empty");
        return false;
    }

    if (!mol) {
        MD_LOG_ERROR("Script Compile: Molecule was not supplied");
        return false;
    }
    
    ir->str = str_copy(src, ir->arena);

    if (ctx_ir) {
        add_ir_ctx(ir, ctx_ir);
    }

    ir->compile_success =
        parse_script(ir) &&
        static_type_check(ir, mol, traj) &&
        // optimize ast?
        //static_evaluation(ir, mol) &&
        extract_dynamic_evaluation_targets(ir) &&
        extract_identifiers(ir) &&
        extract_properties(ir);

    extract_vis_tokens(ir);

#if MD_DEBUG
    //save_expressions_to_json(ir->expressions, md_array_size(ir->expressions), make_cstr("tree.json"));
#endif

    uint64_t hash = 0;
    const size_t num_expr = md_array_size(ir->type_checked_expressions);
    for (size_t i = 0; i < num_expr; ++i) {
        hash = hash_node(ir->type_checked_expressions[i], hash);
#if DEBUG
        //md_logf(MD_LOG_TYPE_DEBUG, "%llu", hash);
#endif
    }
    ir->fingerprint = hash;
#if DEBUG
    //md_log(MD_LOG_TYPE_DEBUG, "---");
#endif

    return ir->compile_success;
}

bool md_script_ir_add_identifier_bitfield(md_script_ir_t* ir, str_t name, const md_bitfield_t* bf) {
    if (!validate_ir(ir)) {
        MD_LOG_ERROR("IR is not valid");
        return false;
    }

    if (!md_script_identifier_name_valid(name)) {
        MD_LOG_ERROR("'"STR_FMT"' is not a valid identifier", STR_ARG(name));
        return false;
    }

    if (!md_bitfield_validate(bf)) {
        MD_LOG_ERROR("The supplied bitfield is not a valid bitfield");
        return false;
    }

    identifier_t* ident = create_identifier(ir, name);
    if (!ident) {
        MD_LOG_ERROR("Failed to create identifier with name '"STR_FMT"', it is perhaps already taken.", STR_ARG(name));
        return false;
    }

    ast_node_t* node = create_node(ir, AST_CONSTANT_VALUE, (token_t){0});
    node->data.ptr = &node->value._bitfield;
    node->data.size = sizeof(md_bitfield_t);
    node->data.type = (type_info_t)TI_BITFIELD;
    md_bitfield_init(&node->value._bitfield, ir->arena);
    md_bitfield_copy(&node->value._bitfield, bf);

    ident->data = &node->data;
    ident->node = node;

    return true;
}

static bool node_contains_identifier_reference(const ast_node_t* node, str_t name) {
    if (node->type == AST_IDENTIFIER) {
        return str_eq(node->ident, name);
    } else {
        for (size_t i = 0; i < md_array_size(node->children); ++i) {
            if (node_contains_identifier_reference(node->children[i], name)) {
                return true;
            }
        }
    }
    return false;
}

bool md_script_ir_contains_identifier_reference(const md_script_ir_t* ir, str_t name) {
    if (!validate_ir(ir)) {
        MD_LOG_ERROR("IR is not valid");
        return false;
    }

    if (!md_script_identifier_name_valid(name)) {
        MD_LOG_ERROR("'"STR_FMT"' is not a valid identifier", STR_ARG(name));
        return false;
    }

    if (!find_identifier(name, ir->identifiers, md_array_size(ir->identifiers))) {
        // name is not even registered as an identifier within the script
        return false;
    }

    for (size_t i = 0; i < md_array_size(ir->parsed_expressions); ++i) {
        const ast_node_t* node = ir->parsed_expressions[i];
        if (node_contains_identifier_reference(node, name)) {
            return true;
        }
    }
    return false;
}

void md_script_ir_free(md_script_ir_t* ir) {
    if (validate_ir(ir)) {
        free_ir(ir);
    }
}

void md_script_ir_clear(md_script_ir_t* ir) {
    if (validate_ir(ir)) {
        reset_ir(ir);
    }
}

size_t md_script_ir_num_warnings(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return 0;
    }
    return md_array_size(ir->warnings);
}

const md_log_token_t* md_script_ir_warnings(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return NULL;
    }
    return ir->warnings;
}

size_t md_script_ir_num_errors(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return 0;
    }
    return md_array_size(ir->errors);
}

const md_log_token_t* md_script_ir_errors(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return NULL;
    }
    return ir->errors;
}

size_t md_script_ir_num_vis_tokens(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return 0;
    }
    return md_array_size(ir->vis_tokens);
}

const md_script_vis_token_t* md_script_ir_vis_tokens(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return NULL;
    }
    return ir->vis_tokens;
}

bool md_script_ir_valid(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return false;
    }
    return ir->compile_success;
}

uint64_t md_script_ir_fingerprint(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return 0;
    }
    return ir->fingerprint;
}

size_t md_script_ir_num_identifiers(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return 0;
    }
    return md_array_size(ir->identifier_names);
}

const str_t* md_script_ir_identifiers(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return NULL;
    }
    return ir->identifier_names;
}

size_t md_script_ir_property_count(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return 0;
    }
    return md_array_size(ir->property_names);
}

const str_t* md_script_ir_property_names(const md_script_ir_t* ir) {
    if (!validate_ir(ir)) {
        return 0;
    }
    return ir->property_names;
}

/*
size_t md_script_ir_property_id_filter_on_flags(md_script_property_id_t* out_ids, size_t out_cap, const md_script_ir_t* ir, md_script_property_flags_t flags) {
    if (!validate_ir(ir)) {
        return 0;
    }

    size_t count = 0;
    ASSERT(md_array_size(ir->property_ids) == md_array_size(ir->property_infos));
    for (size_t i = 0; i < md_array_size(ir->property_ids); ++i) {
        if (count == out_cap) return count;
        if (ir->property_infos[i].flags & flags) {
            out_ids[count++] = ir->property_ids[i];
        }
    }
    return count;
}
*/

md_script_property_flags_t md_script_ir_property_flags(const md_script_ir_t* ir, str_t name) {
    if (!validate_ir(ir)) {
        return 0;
    }

    ASSERT(md_array_size(ir->property_names) == md_array_size(ir->property_flags));
    for (size_t i = 0; i < md_array_size(ir->property_names); ++i) {
        if (str_eq(ir->property_names[i], name)) {
            return ir->property_flags[i];
        }
    }

    return MD_SCRIPT_PROPERTY_FLAG_NONE;
}

const md_script_vis_payload_o* md_script_ir_property_vis_payload(const md_script_ir_t* ir, str_t name) {
    if (!validate_ir(ir)) {
        return 0;
    }

    ASSERT(md_array_size(ir->property_names) == md_array_size(ir->property_nodes));
    for (size_t i = 0; i < md_array_size(ir->property_names); ++i) {
        if (str_eq(ir->property_names[i], name)) {
            return (const md_script_vis_payload_o*)ir->property_nodes[i];
        }
    }

    return NULL;
}

static md_script_eval_t* create_eval(md_allocator_i* alloc) {
    md_allocator_i* arena = md_arena_allocator_create(alloc, MEGABYTES(1));
    md_script_eval_t* eval = md_alloc(arena, sizeof(md_script_eval_t));
    MEMSET(eval, 0, sizeof(md_script_eval_t));
    eval->magic = SCRIPT_EVAL_MAGIC;
    eval->arena = arena;
    return eval;
}

md_script_eval_t* md_script_eval_create(size_t num_frames, const md_script_ir_t* ir, md_allocator_i* alloc) {
    if (num_frames == 0) {
        MD_LOG_ERROR("Script eval: Number of frames was 0");
        return NULL;
    }
    if (!ir) {
        MD_LOG_ERROR("Script eval: IR was null");
        return NULL;
    }
    if (!md_script_ir_valid(ir)) {
        MD_LOG_ERROR("Script eval: IR was not valid");
        return NULL;
    }
    if (!alloc) {
        MD_LOG_ERROR("Script eval: Allocator was null");
        return NULL;
    }

    if (md_array_size(ir->property_names) == 0) {
        MD_LOG_INFO("No properties present in ir");
        return NULL;
    }

    md_script_eval_t* eval = create_eval(alloc);

    eval->ir_fingerprint = ir->fingerprint;

    eval->frame_count = num_frames;
    md_bitfield_init(&eval->frame_mask, eval->arena);
    md_bitfield_reserve_range(&eval->frame_mask, 0, num_frames);

    size_t num_props = md_array_size(ir->property_names);
    for (size_t i = 0; i < num_props; ++i) {
        md_array_push(eval->property_names, ir->property_names[i], eval->arena);
        md_script_property_data_t* data = md_array_push(eval->property_data, (md_script_property_data_t){0}, eval->arena);
        const ast_node_t* node = ir->property_nodes[i];
        allocate_property_data(data, ir->property_flags[i], node->data.type, num_frames, eval->arena);
        data->unit = node->data.unit;
    }
    
    md_array_resize(eval->property_dist_count, num_props, eval->arena);
    md_array_resize(eval->property_dist_mutex, num_props, eval->arena);
        
    for (size_t i = 0; i < num_props; ++i) {
        clear_property_data(&eval->property_data[i]);
        eval->property_dist_count[i] = 0;
        md_mutex_init(&eval->property_dist_mutex[i]);
    }

    md_mutex_init(&eval->frame_lock);

    return eval;
}

void md_script_eval_clear_data(md_script_eval_t* eval) {
    ASSERT(eval);
    ASSERT(eval->magic == SCRIPT_EVAL_MAGIC);
    md_bitfield_clear(&eval->frame_mask);
    for (size_t i = 0; i < md_array_size(eval->property_data); ++i) {
        clear_property_data(&eval->property_data[i]);
        eval->property_dist_count[i] = 0;
    }
    eval->interrupt = false;
}

bool md_script_eval_frame_range(md_script_eval_t* eval, const struct md_script_ir_t* ir, const struct md_molecule_t* mol, const struct md_trajectory_i* traj, uint32_t frame_beg, uint32_t frame_end) {
    ASSERT(eval);

    if (!ir) {
        MD_LOG_ERROR("Script eval: Immediate representation was null");
        return false;
    }
    if (!mol) {
        MD_LOG_ERROR("Script eval: Molecule was null");
        return false;
    }
    if (!traj) {
        MD_LOG_ERROR("Script eval: Trajectory was null");
        return false;
    }

    const uint32_t num_frames = (uint32_t)md_trajectory_num_frames(traj);
    if (num_frames == 0) {
        MD_LOG_ERROR("Script eval: Trajectory was empty");
        return false;
    }
    if (frame_beg > frame_end || frame_end > num_frames) {
        MD_LOG_ERROR("Script eval: Invalid frame range");
        return false;
    }

    if (md_array_size(eval->property_data) == 0) {
        MD_LOG_INFO("Script eval: No properties present, nothing to evaluate");
        return false;
    }
    
    bool result = eval_properties(eval, mol, traj, ir, frame_beg, frame_end);

    const uint64_t fingerprint = generate_fingerprint();
    for (size_t i = 0; i < md_array_size(eval->property_data); ++i) {
        eval->property_data[i].fingerprint = fingerprint;
    }

    return result;
}

uint64_t md_script_eval_ir_fingerprint(const md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        return eval->ir_fingerprint;
    }
    return 0;
}

void md_script_eval_free(md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        for (size_t i = 0; i < md_array_size(eval->property_dist_mutex); ++i) {
            md_mutex_destroy(&eval->property_dist_mutex[i]);
        }
        md_mutex_destroy(&eval->frame_lock);
        md_arena_allocator_destroy(eval->arena);
    }
}

size_t md_script_eval_property_count(const md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        return md_array_size(eval->property_names);
    }
    return 0;
}

const md_script_property_data_t* md_script_eval_property_data(const md_script_eval_t* eval, str_t name) {
    if (validate_eval(eval)) {
        for (size_t i = 0; i < md_array_size(eval->property_names); ++i) {
            if (str_eq(eval->property_names[i], name)) {
                return &eval->property_data[i];
            }
        }
    }
    return NULL;
}

size_t md_script_eval_frame_count(const md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        return eval->frame_count;
    }
    return 0;
}

const md_bitfield_t* md_script_eval_frame_mask(const md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        return &eval->frame_mask;
    }
    return NULL;
}

void md_script_eval_interrupt(md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        eval->interrupt = true;
    }
}

static bool eval_expression(data_t* dst, str_t expr, md_molecule_t* mol, md_allocator_i* alloc) {
    SETUP_TEMP_ALLOC(GIGABYTES(4));

    md_script_ir_t* ir = create_ir(&temp_alloc);
    ir->str = str_copy(expr, ir->arena);

    tokenizer_t tokenizer = tokenizer_init(ir->str);
    bool result = false;

    ast_node_t* node = parse_expression(&(parse_context_t){ .ir = ir, .tokenizer = &tokenizer, .temp_alloc = &temp_alloc});
    if (node) {
        eval_context_t ctx = {
            .ir = ir,
            .mol = mol,
            .temp_arena = &vm_arena,
            .temp_alloc = &temp_alloc,
            .alloc = &temp_alloc,
        };

        if (static_check_node(node, &ctx)) {
            allocate_data(dst, node->data.type, alloc);
            if (evaluate_node(dst, node, &ctx)) {
                if (dst->type.base_type == TYPE_STRING) {
                    // @HACK: We reallocate strings here here: If the data type is a str_t, then it gets a shallow copy
                    // Which means that the actual string data is contained within the ir->arena => temp_alloc
                    const uint64_t num_elem = element_count(*dst);
                    str_t* dst_str = dst->ptr;
                    for (uint64_t i = 0; i < num_elem; ++i) {
                        dst_str[i] = str_copy(dst_str[i], alloc);
                    }
                }
                result = true;
            }
        }
    }

    if (ir->errors) {
        for (size_t i = 0; i < md_array_size(ir->errors); ++i) {
            MD_LOG_ERROR(""STR_FMT"", ir->errors[i].text.len, ir->errors[i].text.ptr);
        }
    }

    FREE_TEMP_ALLOC();

    return result;
}

bool md_filter_evaluate(md_array(md_bitfield_t)* bitfields, str_t expr, const md_molecule_t* mol, const md_script_ir_t* ctx_ir, bool* is_dynamic, char* err_buf, size_t err_cap, md_allocator_i* alloc) {
    ASSERT(bitfields);
    ASSERT(mol);
    ASSERT(alloc);

    bool success = false;

    SETUP_TEMP_ALLOC(GIGABYTES(4));

    md_script_ir_t* ir = create_ir(&temp_alloc);
    ir->str = str_copy(expr, ir->arena);

    if (ctx_ir) {
        add_ir_ctx(ir, ctx_ir);
    }

    tokenizer_t tokenizer = tokenizer_init(ir->str);

    parse_context_t parse_ctx = {
        .ir = ir,
        .tokenizer = &tokenizer,
        .node = 0,
        .temp_alloc = &temp_alloc,
    };

    ir->stage = "Filter evaluate";
    ir->record_log = true;

    ast_node_t* node = parse_expression(&parse_ctx);
    node = prune_expressions(node);
    if (node) {
        eval_context_t ctx = {
            .ir = ir,
            .mol = mol,
            .temp_arena = &vm_arena,
            .temp_alloc = &temp_alloc,
            .alloc = &temp_alloc,
            .initial_configuration = {
                .x = mol->atom.x,
                .y = mol->atom.y,
                .z = mol->atom.z,
            },
            .eval_flags = EVAL_FLAG_NO_STATIC_EVAL,
        };

        if (static_check_node(node, &ctx)) {
            if (node->data.type.base_type == TYPE_BITFIELD) {               
                if (type_info_array_len(node->data.type) == -1 && (node->flags & FLAG_DYNAMIC_LENGTH)) {
                    finalize_type(&node->data.type, node, &ctx);
                }

                data_t data = {0};
                allocate_data(&data, node->data.type, &temp_alloc);

                if (evaluate_node(&data, node, &ctx)) {
                    const int64_t len = type_info_array_len(data.type);
                    const md_bitfield_t* bf_arr = data.ptr;
                    if (bf_arr) {
                        for (int64_t i = 0; i < len; ++i) {
                            md_bitfield_t bf = {0};
                            md_bitfield_init(&bf, alloc);
                            md_bitfield_copy(&bf, &bf_arr[i]);
                            md_array_push(*bitfields, bf, alloc);
                        }
                    }
                    if (is_dynamic) *is_dynamic = (bool)(node->flags & FLAG_DYNAMIC);
                    success = true;
                }
            }
            if (!success) {
                if (err_buf) {
                    snprintf(err_buf, err_cap, "Expression did not evaluate to a bitfield\n");
                } else {
                    MD_LOG_ERROR("md_filter: Expression did not evaluate to a valid bitfield\n");
                }
            }
        }
    }

    if (err_buf) {
        size_t len = 0;
        for (size_t i = 0; i < md_array_size(ir->errors); ++i) {
            size_t space_left = err_cap - MIN(len, err_cap);
            if (!space_left) break;
            len += snprintf(err_buf + len, space_left, ""STR_FMT"", (int)ir->errors[i].text.len, ir->errors[i].text.ptr);
        }
    }

    FREE_TEMP_ALLOC();
    return success;
}

bool md_filter(md_bitfield_t* dst_bf, str_t expr, const struct md_molecule_t* mol, const struct md_script_ir_t* ctx_ir, bool* is_dynamic, char* err_buf, size_t err_cap) {
    ASSERT(mol);

    if (!dst_bf || !md_bitfield_validate(dst_bf)) {
        MD_LOG_ERROR("md_filter: Passed in bitfield was NULL or not valid.");
        return false;
    }

    if (!mol) {
        MD_LOG_ERROR("md_filter: Passed in molecule was NULL");
        return false;
    }

    if (mol->atom.count == 0) {
        MD_LOG_ERROR("md_filter: Passed in molecule was empty");
        return false;
    }

    md_bitfield_clear(dst_bf);

    bool success = false;

    SETUP_TEMP_ALLOC(GIGABYTES(4));

    md_script_ir_t* ir = create_ir(&temp_alloc);
    ir->str = str_copy(expr, ir->arena);

    if (ctx_ir) {
        add_ir_ctx(ir, ctx_ir);
    }

    tokenizer_t tokenizer = tokenizer_init(ir->str);

    parse_context_t parse_ctx = {
        .ir = ir,
        .tokenizer = &tokenizer,
        .node = 0,
        .temp_alloc = &temp_alloc,
    };

    ir->stage = "Filter evaluate";
    ir->record_log = true;

    ast_node_t* node = parse_expression(&parse_ctx);
    node = prune_expressions(node);
    if (node) {
        eval_context_t ctx = {
            .ir = ir,
            .mol = mol,
            .temp_arena = &vm_arena,
            .temp_alloc = &temp_alloc,
            .alloc = &temp_alloc,
            .initial_configuration = {
                .x = mol->atom.x,
                .y = mol->atom.y,
                .z = mol->atom.z,
            },
            .eval_flags = EVAL_FLAG_FLATTEN | EVAL_FLAG_NO_STATIC_EVAL,
        };

        if (static_check_node(node, &ctx)) {
            if (node->data.type.base_type == TYPE_BITFIELD) {
                int len = (int)type_info_array_len(node->data.type);
                if (len == -1 && (node->flags & FLAG_DYNAMIC_LENGTH)) {
                    finalize_type(&node->data.type, node, &ctx);
                    len = (int)type_info_array_len(node->data.type);
                }
                
                if (len == 1) {
                    data_t data = {0};
                    data.type = node->data.type;
                    data.ptr = dst_bf;
                    data.size = sizeof(md_bitfield_t);

                    if (evaluate_node(&data, node, &ctx)) {
                        success = true;
                        if (is_dynamic) {
                            *is_dynamic = (bool)(node->flags & FLAG_DYNAMIC);
                        }
                    }
                } else if (len > 1) {
                    data_t data = {0};
                    allocate_data(&data, node->data.type, &temp_alloc);
                    
                    if (evaluate_node(&data, node, &ctx)) {
                        success = true;
                        if (is_dynamic) {
                            *is_dynamic = (bool)(node->flags & FLAG_DYNAMIC);
                        }
                    }

                    if (success) {
                        const md_bitfield_t* src_bf = (const md_bitfield_t*)data.ptr;
                        for (int i = 0; i < len; ++i) {
                            md_bitfield_or_inplace(dst_bf, &src_bf[i]);
                        }
                    }
                }
                else {
                    MD_LOG_DEBUG("md_filter: Expression did not evaluate to a valid bitfield\n");
                    snprintf(err_buf, err_cap, "Expression did not evaluate to a valid bitfield\n");
                }
            } else {
                MD_LOG_DEBUG("md_filter: Expression did not evaluate to a valid bitfield\n");
                snprintf(err_buf, err_cap, "Expression did not evaluate to a bitfield\n");
            }
        }
    }

    if (err_buf) {
        size_t len = 0;
        for (size_t i = 0; i < md_array_size(ir->errors); ++i) {
            int space_left = MAX(0, (int)(err_cap - len));
            if (space_left) len += snprintf(err_buf + len, (size_t)space_left, ""STR_FMT"\n", (int)ir->errors[i].text.len, ir->errors[i].text.ptr);
        }
    }

    FREE_TEMP_ALLOC();
    return success;
}

#define VIS_MAGIC 0xbc6169abd9628947

static void do_vis_eval(const ast_node_t* node, eval_context_t* ctx) {

    if (!node) {
        MD_LOG_DEBUG("Vis Eval: Failed to visualize node: node was NULL\n");
        return;
    }
    
    // Take some shortcuts in specific cases
    switch (node->type) {
    case AST_IDENTIFIER:
        ASSERT(node->children);
        ASSERT(md_array_size(node->children) == 1);
        node = node->children[0];
        break;
    case AST_ASSIGNMENT:
        ASSERT(node->children);
        ASSERT(md_array_size(node->children) == 2);
        node = node->children[1];
        break;
    default:
        break;
    }

    if (node->data.type.base_type == TYPE_BITFIELD) {
        data_t data = {0};

        type_info_t type = node->data.type;
        if (is_variable_length(type)) {
            if (!finalize_type(&type, node, ctx)) {
                md_log(MD_LOG_TYPE_DEBUG, "Vis Eval: Failed to finalize type for variable length expression");
                return;
            }
        }

        allocate_data(&data, type, ctx->temp_alloc);
        if (data.ptr) {
            evaluate_node(&data, node, ctx);
            const md_bitfield_t* bf_arr = data.ptr;
            for (size_t i = 0; i < element_count(data); ++i) {
                md_bitfield_or_inplace(&ctx->vis->atom_mask, &bf_arr[i]);
            }
            free_data(&data, ctx->temp_alloc);
        } else {
            md_log(MD_LOG_TYPE_DEBUG,"Failed to allocate data for bitfield visualization");
        }
    } else if (node->flags & FLAG_COORDINATE) {
        data_t tmp_data = {0};
        const data_t* data = &node->data;
        if (!node->data.ptr) {
            // Evaluate the non constant value expression here into a data_t item and visualize it.
            allocate_data(&tmp_data, node->data.type, ctx->temp_alloc);
            if (!evaluate_node(&tmp_data, node, ctx)) {
                MD_LOG_ERROR("Vis Eval: Failed to evaluate coordinate expression");
                return;
			}
		    data = &tmp_data;
        }
        coordinate_visualize(*data, ctx);
    } else {
        evaluate_node(NULL, node, ctx);
    }
}

static void visualize_node(const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(ctx->vis);
    if (node->type != AST_CONTEXT && node->num_contexts) {
        ASSERT(node->contexts);
        for (size_t i = 0; i < node->num_contexts; ++i) {
            ctx->mol_ctx = &node->contexts[i];
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

bool md_script_vis_eval_payload(md_script_vis_t* vis, const md_script_vis_payload_o* payload, int subidx, const md_script_vis_ctx_t* vis_ctx, md_script_vis_flags_t flags) {
    ASSERT(vis);
    ASSERT(payload);
    ASSERT(vis_ctx);
    
    if (vis->magic != VIS_MAGIC) {
        MD_LOG_ERROR("Visualize: vis object not initialized.");
        return false;
    }

    if (subidx < -1) {
        MD_LOG_ERROR("Visualize: invalid subidx");
        return false;
    }

    if (flags == 0) flags = 0xFFFFFFFFU;
    
    md_bitfield_clear(&vis->atom_mask);

    SETUP_TEMP_ALLOC(GIGABYTES(4));

    int64_t num_atoms = md_trajectory_num_atoms(vis_ctx->traj);
    float* init_x = md_vm_arena_push(&vm_arena, num_atoms * sizeof(float));
    float* init_y = md_vm_arena_push(&vm_arena, num_atoms * sizeof(float));
    float* init_z = md_vm_arena_push(&vm_arena, num_atoms * sizeof(float));

    md_trajectory_frame_header_t header = { 0 };
    if (vis_ctx->traj) {
        md_trajectory_load_frame(vis_ctx->traj, 0, &header, init_x, init_y, init_z);
    } else {
        init_x = vis_ctx->mol->atom.x;
        init_y = vis_ctx->mol->atom.y;
        init_z = vis_ctx->mol->atom.z;
    }

    md_array(irange_t) ranges = 0;
    if (subidx > -1) {
        irange_t range = {subidx, subidx+1};
        md_array_push(ranges, range, &temp_alloc);
    }

    eval_context_t ctx = {
        .ir = (md_script_ir_t*)vis_ctx->ir,
        .mol = vis_ctx->mol,
        .traj = vis_ctx->traj,
        .temp_arena = &vm_arena,
        .temp_alloc = &temp_alloc,
        .vis = vis,
        .vis_flags = flags,
        .vis_structure = &vis->atom_mask,
        .frame_header = &header,
        .initial_configuration = {
            .header = &header,
            .x = init_x,
            .y = init_y,
            .z = init_z
        },
        .subscript_ranges = ranges,
    };

    visualize_node((ast_node_t*)payload, &ctx);

    // This just to see that we conceptually did not make any errors when evaluating sub_contexts
    ASSERT(ctx.vis_structure == &vis->atom_mask);

    // Append all structures into the 'global' atom mask
    for (size_t i = 0; i < md_array_size(vis->structures); ++i) {
        md_bitfield_or_inplace(&vis->atom_mask, &vis->structures[i]);
    }

    FREE_TEMP_ALLOC();

    return true;
}

void md_script_vis_init(md_script_vis_t* vis, md_allocator_i* alloc) {
    ASSERT(vis);
    ASSERT(alloc);
    MEMSET(vis, 0, sizeof(md_script_vis_t));
    vis->magic = VIS_MAGIC;
    vis->alloc = md_arena_allocator_create(alloc, MEGABYTES(1));

    md_bitfield_init(&vis->atom_mask, vis->alloc);
}

bool md_script_vis_free(md_script_vis_t* vis) {
    ASSERT(vis);
    if (vis->magic == VIS_MAGIC) {
        ASSERT(vis->alloc);
        md_arena_allocator_destroy(vis->alloc);
    }
    MEMSET(vis, 0, sizeof(md_script_vis_t));
    return true;
}

bool md_script_vis_clear(md_script_vis_t* vis) {
    if (vis->magic != VIS_MAGIC) {
        MD_LOG_ERROR("Vis: Failed to clear object: not initialized");
        return false;
    }

    md_arena_allocator_reset(vis->alloc);
    md_allocator_i* alloc = vis->alloc;
    MEMSET(vis, 0, sizeof(md_script_vis_t));
    vis->magic = VIS_MAGIC;
    vis->alloc = alloc;

    md_bitfield_init(&vis->atom_mask, vis->alloc);
    return true;
}

bool md_script_identifier_name_valid(str_t ident) {
    if (!ident.ptr) return false;
    if (!ident.len) return false;

    const char* beg = str_beg(ident);
    const char* end = str_end(ident);

    if (!is_alpha(*beg) && *beg != '_') return false;
    for (const char* c = beg + 1; c < end; ++c) {
        if (!is_alpha(*c) && (*c != '_') && !is_digit(*c)) return false;
    }

    return true;
}
