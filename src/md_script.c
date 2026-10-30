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
*   PROPERTIES
*   Evaluated properties are published as attributes under 'script/' (see md_script_eval_attributes).
*    
**/

#include <md_script.h>
#include <md_contact.h>
#include <md_pucker.h>
#include <md_match.h>
#include <md_smiles.h>

#include <md_system.h>
#include <md_filter.h>
#include <md_util.h>

#include <core/md_coord_stream.h>
#include <core/md_common.h>
#include <core/md_compiler.h>
#include <core/md_log.h>
#include <core/md_bitset.h>
#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_array.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_unit.h>
#include <core/md_spatial_acc.h>
#include <core/md_vec_math.h>
#include <core/md_writer.h>
#include <core/md_writer.h>
#include <core/md_hash.h>
#include <core/md_intrinsics.h>

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
#include <intrin.h>
#pragma warning(disable:4063) // Single character tokens not being valid tokens
#endif

// The one atomic operation the IR needs: its id (see next_ir_id)
static inline int64_t atomic_add_i64(volatile int64_t* x, int64_t v) {
#if MD_COMPILER_MSVC
    return _InterlockedExchangeAdd64(x, v) + v;
#else
    return __atomic_add_fetch(x, v, __ATOMIC_ACQ_REL);
#endif
}

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
    TOKEN_INT_DIV, // '//'
    TOKEN_OF,     // Reserved keyword
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
    AST_ATTRIBUTE,          // attr("edr/potential"): a temporal attribute of the system, read at the evaluated frame
    AST_IDENTIFIER,
    AST_ARRAY,              // This is also just to declare an array in tree form, the node tells us the types, the children contains the values... inefficient, yes.
    AST_ARRAY_SUBSCRIPT,    // [..]
    AST_CAST,
    AST_UNARY_NEG,
    AST_NOT,
    AST_MUL,
    AST_DIV,
    AST_INT_DIV,
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

typedef enum associativity_t {
    ASSOCIATIVITY_LEFT_TO_RIGHT,
    ASSOCIATIVITY_RIGHT_TO_LEFT,
} associativity_t;

typedef enum base_type_t {
    TYPE_UNDEFINED = 0,
    TYPE_FLOAT,
    TYPE_INT,
    TYPE_BOOL,
    TYPE_FRANGE,
    TYPE_IRANGE,
    TYPE_ATOMSET,       // A set of atoms (md_bitset_t over the atom indices): what every selector gives
    TYPE_STRING,
    TYPE_COORDINATE,    // This is a pseudo type which signifies that the underlying type is something that can be interpereted as a coordinate (INT, IRANGE, BITFIELD or float[3])
    TYPE_CONTACT,       // Opaque: a set of group pairs in contact (md_contact_set_t). Its arrays live in the evaluation's allocator.
    TYPE_DISTRIBUTION,  // A bin of a histogram: a float. Shaped [n][bins], one distribution per row. The weight of each bin,
                        // what normalizes it, is not part of the value: it sits beside it in data_t.weight.
    TYPE_PLANE,         // The plane fitted to a group of atoms (plane_t), read through centroid(), normal(), radius(), ...
    TYPE_PUCKER,        // The puckering of a ring (pucker_rec_t), read through amplitude(), phase(), theta(), pseudorotation()
    TYPE_BONDSET,       // A set of bonds (bondset_t): what bond() gives
    TYPE_BOND_DELTA,    // Opaque: the bonds formed and broken in a frame (bond_delta_value_t), read through formed(), broken(), bond(d)
} base_type_t;

// Appended to an error where a bond set is given for something which takes atoms: within an expression the two do not
// convert into each other (a filter of atoms alone takes the atoms of the bonds it gives, md_filter)
#define BONDS_DO_NOT_MIX_HINT ". Bonds are not atoms: atom(b) gives the atoms of the bonds b"

typedef enum flags_t {
    // Procedure Flags
    FLAG_SYMMETRIC_ARGS             = 0x00001, // The first two arguments can be swapped to match: the procedure is called with them swapped, so only for commutative ones (a - s is not s - a)
    FLAG_ARGS_EQUAL_LENGTH          = 0x00002, // Arguments should have the same array length
    FLAG_DEDUCE_LENGTH_FROM_ARG     = 0x00004, // Return type array length matches argument arrays' length
    FLAG_DYNAMIC_LENGTH             = 0x00008, // Return type has a varying length which is not constant over frames
    FLAG_QUERYABLE_LENGTH           = 0x00010, // Marks a procedure as queryable, meaning it can be called with NULL as dst to query the length of the resulting array
    FLAG_STATIC_VALIDATION          = 0x00020, // Marks a procedure as validatable, meaning it can be called with NULL as dst to validate it in a static context during compilation
    FLAG_FLATTEN                    = 0x00040, // Hints that a procedure want flattened atomsets as input, e.g. within, this can be propagated to the arguments during static type checking
    FLAG_UNIT_OF_ARG                = 0x00080, // The result is in the unit of the first argument: a reduction over it (sum, min, max of an array)
    
    // Node flags
    FLAG_CONSTANT                   = 0x00100, // Marks the node as constant and should not be modified. It should also have data in its ptr
    FLAG_CONTEXTS_EQUIVALENT        = 0x00200, // Marks the contexts in the node as equivalent, meaning it
    FLAG_COORDINATE                 = 0x00400, // Marks the node as a coordinate, meaning it can be visualized as a coordinate
    FLAG_IMPLICIT                   = 0x00800, // Not written in the source: a default or absent argument the binding of a call supplied

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

// The plane through a group of atoms or through positions (plane()): the least squares plane of the positions
typedef struct plane_t {
    vec3_t centroid;        // Mean position of the atoms (Å)
    vec3_t normal;          // Unit normal, oriented (see plane())
    vec3_t axis;            // Unit direction in the plane along which the atoms spread the most
    float  extent[2];       // RMS spread of the atoms along axis and along normal x axis (Å)
    float  radius;          // Mean distance of the atoms from the centroid (Å)
    float  deviation;       // RMS distance of the atoms from the plane (Å)
    int    count;           // Positions fitted
} plane_t;

// The puckering of a ring (pucker()), as md_pucker_ring makes it
typedef md_pucker_ring_t pucker_rec_t;

// A set of bonds (TYPE_BONDSET), in two parts: the persistent bonds of the system by index (both ranges of
// md_bond_data_t: the graph and the non-covalent bonds), and the bonds which are none of them (computed for a frame: the
// bonds formed in it) by atom pair. A pair with a persistent bond is always in the first part, never in the second, so
// sets compare and combine by their bonds whatever produced them. The indices are those of the topology the set was
// made from. The payloads live in the allocator of the data holding the set (data_bitset_alloc), as those of an atomset.
typedef struct bondset_t {
    md_bitset_t            bonds;
    uint32_t               num_extra;
    const md_atom_pair_t*  extra_pair;      // [num_extra] the lower atom first, sorted
    const md_bond_flags_t* extra_flags;     // [num_extra]
} bondset_t;

// The bonds formed and broken in a frame within a scope (TYPE_BOND_DELTA, bond_delta()), see md_util_bond_delta
typedef struct bond_delta_value_t {
    md_bitset_t touched;    // The persistent bonds with an atom in the scope
    md_bitset_t broken;     // The bonds of the graph broken, of them
    bondset_t   formed;     // The bonds formed, which the system does not have: extra pairs only
} bond_delta_value_t;

struct irange_t {
    int beg;
    int end;
    int step; // Optional stride (>=1). Defaults to 1. Interpreted as [beg, end) with indices beg, beg+step, ... < end after 0-based transform
};

struct type_info_t {
    base_type_t  base_type;
    int dim[MAX_NUM_DIMS]; // The dimensionality of a single element, each entry in dim represents a multidimensional length
};

// There is a conceptual twist here if len should belong to the type_info_t or if it should belong to data_t
// We need to differentiate between different lengths of arrays with the same type in our system, so for now it should probably reside within the type_info
// I.e. the array length should be part of the fundamental type. This only strenghtens the type-system

// PERIODIC VALUES. A value which is the same as itself plus a whole number of periods: an angle about an axis (a
// dihedral is the same at -180 and at 180), the phase of a pucker. [beg, end[ is a period, the one its values are given
// in, and end - beg its length. seam is where, in [beg, end[, a plot of the values best cuts the circle they lie on to lay
// it out flat: where they are least likely to be, so that they cross the cut as seldom as possible. Every value which
// crosses it jumps the length of the period in the plot, and is drawn across the whole of it.
//
// A procedure whose result is periodic states the period and the seam (static_backchannel_t), and the static check
// carries them through the arithmetic done on it (see period_of_operator): dihedral(a,b,c,d) - dihedral(e,f,g,h) is an
// angle as well, and so is dihedral(a,b,c,d) * 180 / PI, in a period of 360 rather than 2 PI. What no longer lies on a
// circle (abs, cos, a comparison) is not periodic. Not periodic: end <= beg, as a zeroed one is.
typedef struct period_t {
    float beg;
    float end;
    float seam;     // in [beg, end[
} period_t;

// This is the struct which is passed as arguments into the procedures
struct data_t {
    type_info_t type;
    void*       ptr;    // Pointer to the data
    size_t      size;   // Size in bytes of data (This we want to determine during the static check so we can allocate the data when evaluating)

    // The allocator ptr (and weight) came from. The payloads of the bitsets held in ptr are allocated with it
    // too, so they live exactly as long as the values.
    md_allocator_i* alloc;

    frange_t    value_range;
    md_unit_t   unit[2];
    period_t    period;     // Of a periodic value (see PERIODIC VALUES), in the unit of the values (unit[1])

    // TYPE_DISTRIBUTION only: the weight of every bin, element for element beside ptr, so at the same byte offsets.
    // A bin is shown as value / weight. Both are additive (over frames, over merged bins, over pooled distributions),
    // which is why the weight is kept rather than divided in: the ratio of sums is the estimate, not the sum of ratios.
    float*      weight;
};

// This is data stored directly in the nodes of the AST tree
union value_t {
    str_t           _string;
    float           _float;
    int             _int;
    frange_t        _frange;
    irange_t        _irange;
    bool            _bool;
    md_bitset_t     _bitset;
};

struct identifier_t {
    str_t       name;
    ast_node_t* node;    // This is the node to evaluate in order to generate the data for the identifier
    data_t*     data;    // This is the data to fill in...
    const md_bitset_t* mol_ctx; // Dynamic identifiers: the context the value was evaluated in
    bool        is_element; // Declared by destructuring a non-array expression, '{x,y,z} = com(...)': the value
    uint32_t    elem_idx;   // is element 'elem_idx' of evaluating 'node', not all of it
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
    md_unit_t unit[2];
    frange_t  value_range;
    period_t  period;       // Set by a procedure whose result is periodic (see PERIODIC VALUES)
} static_backchannel_t;

// The diagnostics of one evaluation. They belong to whoever runs the evaluation, never to the IR,
// which an evaluation only reads.
typedef struct eval_log_t {
    md_array(md_log_token_t) errors;
    md_array(md_log_token_t) warnings;
    md_allocator_i* alloc;
    bool full;      // Set by the first error: as in a compilation, only the first error is kept
} eval_log_t;

typedef struct eval_context_t {
    // The IR being evaluated. Read-only: no evaluation or visualization writes to it, so one IR can be
    // shared by any number of them at once (frame ranges on worker threads, visualizations).
    const md_script_ir_t* ir;

    // The IR being compiled, the same object as ir. Set only by the compilation stages (static type
    // check, static evaluation), which store what they determine in the AST. NULL in every evaluation
    // and visualization: what they produce goes into the context, its allocators or its log.
    md_script_ir_t* compile_ir;

    // Where the diagnostics of an evaluation go (compile_ir is NULL). NULL sends them to md_log, at debug level.
    eval_log_t* log;

    const md_system_t* sys;
    const md_bitset_t* mol_ctx;         // The atomic bit context in which we perform the operation, this can be null, in that case we are not limited to a smaller context and the full molecule is our context

    const float* atom_mass;
    const float* atom_radius;

    md_allocator_i* temp_alloc;         // For allocating transient data  (interface to temp vm arena)
    md_allocator_i* alloc;              // For allocating persistent data (for the duration of the evaluation)

    // Contextual information for static checking 
    token_t  op_token;                  // Token for the operation which is evaluated
    const ast_node_t* op_node;          // Node of the procedure call which is evaluated: identifies its call site
    token_t* arg_tokens;                // Tokens to arguments for contextual information when reporting errors
    flags_t* arg_flags;                 // Flags of arguments
    flags_t proc_flags;                 // Procedure flags of called procedure
    identifier_t* identifiers;          // Evaluated identifiers for references
                                      
    md_script_vis_t* vis;               // These are used when calling a procedure flagged with the VISUALIZE flag so the procedure can fill in the geometry
    uint32_t vis_flags;                 // Visualization flags
    md_bitmask_t* vis_structure;        // Current visualization structure to mark atoms in
    md_bitmask_t* vis_atom_mask;        // The visualization's all encompassing atom mask, while it is built

    // The state being evaluated. Every coordinate, cell and frame read goes through this, so they
    // can never come from different frames.
    const md_system_state_t* cur_state;

    // The frame axis of the run being evaluated, which cur_state->frame indexes. 0 outside an
    // evaluation over a run, where an attr() follows the axis of its own run.
    md_attribute_id_t run_axis_id;

    // The configuration the system's topology was inferred from. Replaces the old
    // initial_configuration, which meant three different things depending on entry point:
    // trajectory frame 0, the system's own coordinates, or whatever the caller passed.
    const md_system_state_t* ref_state;

    // A common spatial acceleration structure which is built lazily.
    // It contains all atoms within the dataset and the cells are at a default resolution ~6Å
    // Will be rebuilt and scaled up if required
    md_spatial_acc_t spatial_acc;
	double spatial_acc_cutoff;     // The cutoff spatial_acc was built for, 0 if it is not built

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

// ### NAMED ARGUMENTS ###
// A procedure may declare names for its parameters, which lets a call supply arguments by name as well as by
// position, as in Python: 'rdf(a, b, cutoff=12)'. Names are declared once per procedure NAME (the 'signatures'
// table in md_script_functions.inl) and mean the same position in every overload of that name. That is what
// allows the arguments of a call to be bound into positional order BEFORE overload resolution, so that
// everything after binding (matching, implicit conversions, evaluation and the C procedures themselves) only
// ever sees positional arguments. A procedure without a signature accepts positional arguments only.

typedef enum param_flags_t {
    PARAM_REQUIRED  = 0,
    PARAM_OPTIONAL  = 1,    // May be omitted as long as no later parameter is given; the call then binds to an overload taking fewer arguments
    PARAM_DEFAULT   = 2,    // If omitted, a constant holding 'def' (of type 'def_type') is supplied in its place
    PARAM_KW_ONLY   = 4,    // Can only be supplied by name
    PARAM_NULLABLE  = 8,    // May be omitted anywhere: the procedure then receives an absent argument in its place, of
                            // type TYPE_UNDEFINED and without data (Python's '= None'). Absent matches any parameter type.
} param_flags_t;

typedef struct param_sig_t {
    str_t       name;
    uint32_t    flags;
    type_info_t def_type;   // PARAM_DEFAULT only, must be a scalar type
    value_t     def;        // PARAM_DEFAULT only
} param_sig_t;

typedef struct proc_sig_t {
    str_t       proc;       // Procedure name, shared by all of its overloads
    size_t      num_params;
    param_sig_t param[MAX_SUPPORTED_PROC_ARGS];
} proc_sig_t;

// ### VALUE DOMAINS ###
// What the strings given to a parameter name, so that completion (md_script_complete) can offer them. Declared per
// procedure NAME and parameter position in the 'param_domains' table in md_script_functions.inl.
typedef enum value_domain_t {
    DOMAIN_NONE = 0,
    DOMAIN_ATOM_NAME,       // Atom (type) names: name(), type(), label()
    DOMAIN_ELEMENT,         // Element symbols of the atoms: element()
    DOMAIN_RES_NAME,        // Component (residue) names: resname(), residue(), component()
    DOMAIN_CHAIN_ID,        // The chain ids the data gives (md_chain_data_t): chain(), chain_id(), auth_id()
    DOMAIN_COUNT_UNIT,      // What count() counts in: 'atom', 'residue', ...
    DOMAIN_ATTR_PATH,       // Paths of the attributes that attr() can read
    DOMAIN_SMILES,          // A SMILES pattern: smiles(). Free text, nothing to offer
    DOMAIN_MATCH_LEVEL,     // Where a match lies: 'structure', 'residue', 'chain', ...
    DOMAIN_MATCH_MODE,      // Which matches are reported: 'unique', 'all', ...
    DOMAIN_MATCH_LABEL,     // What the atoms of a reference are compared by: 'element', 'name'
    DOMAIN_RING_KIND,       // What a ring is: 'aromatic', 'sugar', ...
    DOMAIN_BOND_KIND,       // What a bond is: 'covalent', 'hydrogen', ... (md_bond_kind_names)
    DOMAIN_BOND_RESONANCE,  // 'aromatic', 'delocalized', 'localized'
    DOMAIN_BOND_ORIGIN,     // Where a bond came from: 'file', 'user', ... (md_bond_origin_name)
} value_domain_t;

typedef struct param_domain_t {
    str_t           proc;   // Procedure name, shared by all of its overloads
    uint32_t        param;  // Position of the parameter
    value_domain_t  domain;
} param_domain_t;

// One attribute that attr() reads. Ids rather than pointers: an attribute pointer dies with the next
// create or remove on the table, an id is stable for as long as the path exists.
typedef struct attr_member_t {
    str_t             path;     // absolute
    str_t             label;    // what it is called as a member of a pattern's array (attr_member_labels), empty otherwise
    md_attribute_id_t id;       // what is read
    md_attribute_id_t axis_id;  // its frame axis
} attr_member_t;

// An argument of a procedure call as written in the source
typedef struct named_arg_t {
    str_t   name;   // Empty for a positional argument
    token_t token;  // Spans 'name=value' for a named argument
} named_arg_t;

struct ast_node_t {
    // @TODO, @PERF: Make specific types for each type of Ast_node, where the first member is the type which we can use and cast by.

    ast_type_t          type;    
    uint32_t            index;      // Position in ir->nodes, what a md_script_vis_ref_t refers to it by. Belongs to the node, not to its contents: a copy of a node keeps its own
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

    // ATTRIBUTE
    // Resolved by the static check. A path reads one attribute, a pattern ('*' in it) every one it
    // matches, in path order, as an array over them (also when that is one): its members.
    str_t               attr_path;      // the resolved, absolute path, or the pattern as it was given
    attr_member_t*      attr_members;   // md_array: what is read
    bool                attr_pattern;
    md_attribute_id_t   run_axis_id;    // the frame axis of the run being evaluated
    uint64_t            attr_version;   // the versions of all of it combined, so a reload changes the fingerprint

    // CONTEXT
    size_t              num_contexts;   // Number of arguments for procedure calls
    type_info_t*        lhs_context_types; // For static type checking
    const md_bitset_t* contexts;    // Since contexts have to be statically known at compile time, we store a reference to it.

    // NAMED ARGUMENTS (procedure calls only)
    // Set by the parser when at least one argument was given by name, parallel to 'children'. Consumed by the
    // static check, which binds the arguments into positional order and clears it.
    md_array(named_arg_t) named_args;
};

typedef struct tokenizer_t {
    str_t   str;
    int     cur;
    int     line;
    int     line_offset;  // offset to the current line
} tokenizer_t;

// The topology dependent part of a contacts() call site (md_contact_query_t), prepared once during the static check.
// The groups of a call site cannot depend on the frame, so the evaluation of every frame only reads it.
typedef struct contact_query_entry_t {
    const ast_node_t*  node;
    const md_system_t* sys;
    md_contact_query_t query;
} contact_query_entry_t;

struct md_script_ir_t {
    uint64_t magic;
    uint64_t fingerprint;
    uint64_t id;    // Unique to each compilation (never reused), what a md_script_vis_ref_t is made by. 0 until compiled

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
    md_array(contact_query_entry_t) contact_queries;    // Prepared contacts() call sites, allocated in the arena

    md_array(md_script_property_flags_t) property_flags;    // List of property infos
    md_array(const ast_node_t*)          property_nodes;    // List of property nodes;
    md_array(str_t)                      property_names;    // List of property names
    md_array(uint64_t)                   property_fingerprints; // What each property is computed from (see md_script_ir_property_fingerprint)

    // These are the final products which can be read through the public part of the structure
    md_array(md_log_token_t)        warnings;
    md_array(md_log_token_t)        errors;
    md_array(md_script_vis_token_t) vis_tokens;
    md_array(str_t)                 identifier_names;

    const char* stage;      // This is just to supply a context for the errors i.e. in which compilation stage the error occured
    bool record_log;        // This is to toggle if new errors should be registered... We don't want to flood the errors
    bool compile_success;
};

// One property of an evaluation.
//
// The attribute table OWNS every buffer a property writes, so this is not a second copy of
// anything. It is the part of a property that is not data: the writable views the per frame loop
// needs (resolved once here rather than looked up by path every frame), and the running state an
// evaluation keeps while it accumulates - the observed extremes and the running mean counter.
// Everything a consumer reads is in the table under script/<ident>; see md_script_eval_attributes.
typedef struct eval_property_t {
    str_t                      ident;
    md_script_property_flags_t kind;        // exactly one of TEMPORAL, DISTRIBUTION or VOLUME

    uint32_t count;         // values per frame (temporal), bins (distribution), voxels (volume)
    size_t   num_values;    // floats behind values

    float*  values;         // script/<ident>
    float*  weights;        // script/<ident>/weight     distribution only

    // Distribution only: the sums over the evaluated frames of the values and of the weights. What is published
    // in values and weights are their means, so values / weights is sum(values) / sum(weights) at all times, also
    // while the evaluation runs, and merging bins (sum over values / sum over weights) stays exact.
    double* sum_values;
    double* sum_weights;
    float*  mean;           // script/<ident>/mean       temporal with a population only
    float*  var;            // script/<ident>/variance
    vec2_t* ext;            // script/<ident>/extent
    vec2_t* range;          // script/<ident>/range      temporal and distribution

    // Of a periodic temporal property (MD_SCRIPT_PROPERTY_FLAG_PERIODIC), as script/<ident>/period and /seam publish
    // it: what makes its summary over the population (mean, var, ext) a circular one
    period_t period;

    float   min_value;      // observed over everything evaluated so far
    float   max_value;

    // Distributions and volumes are running means over the evaluated frames. Frame ranges are
    // evaluated concurrently, so the accumulation is serialized per property.
    uint32_t   accum_count;
    md_mutex_t accum_mutex;

    uint64_t     fingerprint;   // What it is computed from (md_script_ir_property_fingerprint)
    // The frames it holds: a temporal property those it has evaluated, an aggregate those it has added. Each
    // property its own, so that one taken over from a previous evaluation (md_script_eval_create_from) keeps its
    // frames while a new one starts without any. Guarded by the eval's frame_lock.
    md_bitmask_t done;

    // A distribution kept per frame (md_script_eval_desc_t::frame_budget): the values and the weights of frame f at
    // f * num_values, for the frames in done. Written before the frame is marked done, and not again.
    float* frame_values;
    float* frame_weights;
} eval_property_t;

// A distribution or a volume over frames of the caller's choosing (md_script_eval_slot_request). Everything in it is
// guarded by its lock, which lives as long as the evaluation, in use or not. Taken after the eval's frame_lock when
// both are held. It is taken into use and let go (active, key, prop, and whether it has sums) under both, so either
// will do to find it, and to tell whether it is one added up as frames are evaluated.
typedef struct eval_slot_t {
    md_mutex_t   lock;
    bool         active;
    uint64_t     key;
    uint32_t     prop;          // the property, an index into the eval's
    uint64_t     version;       // of what it holds (next_version)

    // Reserved over every frame when the slot is first used, so that marking a frame never allocates
    md_bitmask_t wanted;        // the frames asked for
    md_bitmask_t included;      // the frames it holds: those added to the sums, or for a distribution kept per frame
                                // those of the wanted which are evaluated
    uint32_t     count;         // of included

    // Added up over the included frames, for anything not kept per frame. In the heap: let go with the slot.
    double*      sum_values;
    double*      sum_weights;   // distribution only
    size_t       num_values;
} eval_slot_t;

struct md_script_eval_t {
    uint64_t magic;
    uint64_t ir_fingerprint;

    struct md_allocator_i *arena;
    volatile bool interrupt;

    size_t        frame_count;
    md_mutex_t    frame_lock;   // guards the frames each property holds (eval_property_t::done)

    // One per property of the IR, in its order. Sized once at creation and never grown: a published attribute may
    // hold the address of the property it describes.
    md_array(eval_property_t) props;

    // Where the properties above live. Owns every buffer they point at.
    md_attributes_t attributes;

    eval_slot_t* slots;             // MD_SCRIPT_EVAL_MAX_SLOTS of them, in use or not
    uint32_t     num_slots_used;    // slots [0, num_slots_used) have been in use. Guarded by frame_lock.
    volatile uint64_t version_counter;  // of the slots (next_version)

    // Which of the IR's evaluation targets each property needs evaluated for a frame: the one assigning it, and
    // those of every variable it uses, however indirectly. Row p of num_targets, for property p. A frame only
    // some properties still need evaluates what those need, not the whole script.
    uint8_t* needs;
    size_t   num_targets;
};

struct parse_context_t {
    md_script_ir_t* ir;
    tokenizer_t*    tokenizer;
    ast_node_t*     node;   // This represents the current root node in the tree
    md_allocator_i* temp_alloc;
};


// What the parser and the static check guarantee the rest of the compilation and the evaluation: the shape of
// the tree, the types settled on, the sizes of the data. No script should be able to break one of these, and one
// which does is a bug. It is not one to bring the application down with, though: an ASSERT is compiled out of a
// release build, which then carries on into whatever the broken assumption leads to (a NULL dereference, a write
// out of bounds). INTERNAL_CHECK reports the broken invariant as an error at the token instead, and is false, so
// that the caller fails what it was doing:
//     if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 2)) return false;
// INTERNAL_ERROR is the same for a place which should never be reached, such as the default of a switch.
#define INTERNAL_CHECK(ctx, tok, cond) ((cond) ? true : (internal_error((ctx), (tok), "'" #cond "' does not hold", __FILE__, __LINE__), false))
#define INTERNAL_ERROR(ctx, tok, what) internal_error((ctx), (tok), (what), __FILE__, __LINE__)

static void internal_error(const eval_context_t* ctx, token_t tok, const char* what, const char* file, int line);

// ##########################
// ###   CORE FUNCTIONS   ###
// ##########################

static int operator_precedence(ast_type_t type) {
    switch(type) {
    case AST_EXPRESSION:
    case AST_PROC_CALL:
    case AST_CONSTANT_VALUE:
    case AST_ATTRIBUTE:
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
    case AST_INT_DIV:
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

static bool operator_binary(ast_type_t type) {
    switch(type) {
    case AST_MUL:
    case AST_DIV:
    case AST_INT_DIV:
    case AST_ADD:
    case AST_SUB:
    case AST_LT:
    case AST_GT:
    case AST_LE:
    case AST_GE:
    case AST_NE:
    case AST_EQ:
    case AST_AND:
    case AST_XOR:
    case AST_OR:
        return true;
    case AST_EXPRESSION:
    case AST_PROC_CALL:
    case AST_CONSTANT_VALUE:
    case AST_ATTRIBUTE:
    case AST_IDENTIFIER:
    case AST_ARRAY:
    case AST_ARRAY_SUBSCRIPT:
    case AST_FLATTEN:
    case AST_CAST:
    case AST_UNARY_NEG:
    case AST_NOT:
    case AST_OUT:
    case AST_CONTEXT:
    case AST_ASSIGNMENT:
    case AST_TRANSPOSE:
        return false;
    default:
        ASSERT(false);
    }
    return false;
}

static associativity_t operator_associativity(ast_type_t type) {
    switch(type) {
    case AST_EXPRESSION:
    case AST_PROC_CALL:
    case AST_CONSTANT_VALUE:
    case AST_ATTRIBUTE:
    case AST_IDENTIFIER:
    case AST_ARRAY:
    case AST_ARRAY_SUBSCRIPT:
    case AST_FLATTEN:
    case AST_MUL:
    case AST_DIV:
    case AST_INT_DIV:
    case AST_ADD:
    case AST_SUB:
    case AST_LT:
    case AST_GT:
    case AST_LE:
    case AST_GE:
    case AST_NE:
    case AST_EQ:
    case AST_AND:
    case AST_XOR:
    case AST_OR:
        return ASSOCIATIVITY_LEFT_TO_RIGHT;
    case AST_CAST:
    case AST_UNARY_NEG:
    case AST_NOT:
    case AST_OUT:
    case AST_CONTEXT:
    case AST_ASSIGNMENT:
    case AST_TRANSPOSE:
        return ASSOCIATIVITY_RIGHT_TO_LEFT;
    default:
        ASSERT(false);
    }
    return 0;
}

// This is used to get the dim of an element in an array for example
// shifting left by shifting in zeroes
// before [1,3,0,0]
// after  [3,0,0,0]
static void dim_shift_left(int dim[MAX_NUM_DIMS]) {
    for (int i = 0; i < MAX_NUM_DIMS - 1; ++i) {
        dim[i] = dim[i+1];
    }
    dim[MAX_NUM_DIMS-1] = 0;
}

// The semantics are a bit weird, but the idea is that if we shift right, i.e. shifting in dimensions of 1
// that means we are creating an array of length 1, which is the same as a scalar
// before [3,0,0,0]
// after  [1,3,0,0]
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

static int type_info_dim_size(type_info_t ti) {
    return dim_size(ti.dim);
}

static size_t element_count(data_t arg) {
    return arg.type.dim[0];
}

// ### PERIODIC VALUES ###

static inline bool period_valid(period_t p) {
    return p.end > p.beg;
}

static inline double period_length(period_t p) {
    return (double)p.end - (double)p.beg;
}

// x moved by whole periods of length len into [beg, beg + len[
static inline double period_wrap(double x, double beg, double len) {
    const double t = x - beg;
    return beg + (t - len * floor(t / len));
}

// The period [beg, end[ with its seam at seam, moved into it. Not periodic unless end > beg.
static period_t period_make(double beg, double end, double seam) {
    if (!(end > beg)) return (period_t){0};
    return (period_t){ (float)beg, (float)end, (float)period_wrap(seam, beg, end - beg) };
}

// Of the same length, to within what float holds: the two can be told apart from the other's values
static bool period_same_length(period_t a, period_t b) {
    const double la = period_length(a), lb = period_length(b);
    return fabs(la - lb) <= 1.0e-5 * MAX(la, lb);
}

// The period of x + c, for a constant c
static period_t period_shift(period_t p, double c) {
    return period_make(p.beg + c, p.end + c, p.seam + c);
}

// The period of x * c, for a constant c. Its length is |c| times as long: 2x is the same at x and x + p only if it is
// the same at 2x and 2x + 2p. None for c = 0, which is no longer anything but 0.
static period_t period_scale(period_t p, double c) {
    if (c > 0) return period_make(p.beg * c, p.end * c, p.seam * c);
    if (c < 0) return period_make(p.end * c, p.beg * c, p.seam * c);
    return (period_t){0};
}

// The period of a + b (sign 1) or a - b (sign -1), both periodic in periods of the same length: an angle as well, of
// the same length, about the sum (or difference) of the middles of theirs. The seam half a period from that middle,
// where the two are as far apart (or together) as they can be. The difference of two angles is then in [-180, 180[.
static period_t period_combine(period_t a, period_t b, double sign) {
    const double len = period_length(a);
    const double mid = 0.5 * ((double)a.beg + a.end) + sign * 0.5 * ((double)b.beg + b.end);
    return period_make(mid - 0.5 * len, mid + 0.5 * len, mid - 0.5 * len);
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
    case TYPE_ATOMSET: return sizeof(md_bitset_t);
    case TYPE_STRING:   return sizeof(str_t);
    case TYPE_CONTACT:  return sizeof(md_contact_set_t);
    case TYPE_DISTRIBUTION: return sizeof(float);
    case TYPE_PLANE:    return sizeof(plane_t);
    case TYPE_PUCKER:   return sizeof(pucker_rec_t);
    case TYPE_BONDSET:  return sizeof(bondset_t);
    case TYPE_BOND_DELTA: return sizeof(bond_delta_value_t);
    case TYPE_UNDEFINED:
    default:            return 0;
    }
}

static bool is_type_dim_compatible(type_info_t from, type_info_t to) {

    // Some examples of matching, ^ shows len_dim
    // float[4][0][0][0] should match     float[-1][0][0][0]
    //       ^                                   ^
    // float[4][4][0][0] should match     float[-1][4][0][0]
    //       ^                                   ^
    // float[4][4][0][0] should not match float[-1][3][0][0]
    //       ^                                   ^
    //
    // float[1][3][0][0] should match float[3][0][0][0]
    // i.e. an array of size 1 with type float3 should be compatible with float3

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
        case TYPE_ATOMSET: return true;
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

    if (type_info_equal(from, to)) return true;

    if (from.base_type != to.base_type) {
        return false;
    }

    dim_prune_leading_ones(from.dim);
    dim_prune_leading_ones(to.dim);

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

    // This is essentially a logical XOR, we only want to support this operation if we have one unspecified array dimension.
    if (!is_variable_length(from) && is_variable_length(to)) {
        return is_type_dim_compatible(from, to);
    }

    return false;
}

// An absent argument, in place of an omitted PARAM_NULLABLE parameter: matches whatever the procedure takes there
static inline bool is_absent_type(type_info_t ti) {
    return ti.base_type == TYPE_UNDEFINED;
}

static bool compare_type_info_array(const type_info_t a[], const type_info_t b[], size_t num_arg_types) {
    for (size_t i = 0; i < num_arg_types; ++i) {
        if (is_absent_type(a[i])) continue;
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

static bool allocate_data(data_t* data, type_info_t type, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);
    
    // The convention should be that we only allocate data for known lengths: a length of -1 would make for a
    // negative size, taken for an enormous one
    const int64_t array_len = type_info_array_len(type);
    if (is_undefined_type(type) || array_len < 0) {
        MD_LOG_ERROR("Script: internal error, data allocated for a type %s", is_undefined_type(type) ? "which is undefined" : "of unknown length");
        return false;
    }

    // Do the base type allocation (array)
    const int64_t bytes = type_info_element_byte_stride(type) * array_len;

    if (bytes > 0) {
        data->ptr = md_alloc(alloc, bytes);
        if (!data->ptr) {
            MD_LOG_ERROR("Failed to allocate data in script!");
            return false;
        }
        // @NOTE (Robin): This can be a VERY wasteful operation most times as the memory will be overwritten anyways.
        MEMSET(data->ptr, 0, bytes);
    }
    data->size = bytes;
    data->type = type;
    data->value_range = (frange_t){0, 0};
    data->period = (period_t){0};

    data->unit[0] = md_unit_none();
    data->unit[1] = md_unit_none();

    // Bitsets need no initialization: zeroed memory is the empty set
    data->alloc = alloc;

    data->weight = NULL;
    if (type.base_type == TYPE_DISTRIBUTION && bytes > 0) {
        // Weights of 1: a plain histogram, unless the producer says otherwise
        const size_t count = (size_t)bytes / sizeof(float);
        data->weight = md_alloc(alloc, (size_t)bytes);
        if (!data->weight) {
            MD_LOG_ERROR("Failed to allocate data in script!");
            return false;
        }
        for (size_t i = 0; i < count; ++i) {
            data->weight[i] = 1.0f;
        }
    }

    return true;
}

// The weights of the elements of d from byte_offset on: at the same offset as their values (data_t.weight)
static inline float* data_weight_at(const data_t* d, size_t byte_offset) {
    return d->weight ? (float*)((char*)d->weight + byte_offset) : NULL;
}

// Copies bytes of values from src at src_offset to dst at dst_offset, and the weights beside them
static void copy_data_block(data_t* dst, size_t dst_offset, const data_t* src, size_t src_offset, size_t bytes) {
    MEMCPY((char*)dst->ptr + dst_offset, (const char*)src->ptr + src_offset, bytes);
    if (dst->weight) {
        float* w = data_weight_at(dst, dst_offset);
        if (src->weight) {
            MEMCPY(w, data_weight_at(src, src_offset), bytes);
        } else {
            for (size_t i = 0; i < bytes / sizeof(float); ++i) w[i] = 1.0f;
        }
    }
}

static bondset_t bondset_copy(bondset_t s, md_allocator_i* alloc) {
    bondset_t r = { .bonds = md_bitset_copy(s.bonds, alloc), .num_extra = s.num_extra };
    if (s.num_extra) {
        md_atom_pair_t*  pair  = md_alloc(alloc, sizeof(md_atom_pair_t)  * s.num_extra);
        md_bond_flags_t* flags = md_alloc(alloc, sizeof(md_bond_flags_t) * s.num_extra);
        MEMCPY(pair,  s.extra_pair,  sizeof(md_atom_pair_t)  * s.num_extra);
        MEMCPY(flags, s.extra_flags, sizeof(md_bond_flags_t) * s.num_extra);
        r.extra_pair  = pair;
        r.extra_flags = flags;
    }
    return r;
}

static void bondset_free(bondset_t* s, md_allocator_i* alloc) {
    md_bitset_free(&s->bonds, alloc);
    if (s->num_extra) {
        md_free(alloc, (void*)s->extra_pair,  sizeof(md_atom_pair_t)  * s->num_extra);
        md_free(alloc, (void*)s->extra_flags, sizeof(md_bond_flags_t) * s->num_extra);
    }
    MEMSET(s, 0, sizeof(bondset_t));
}

// Whether the values of a type carry a payload of their own (copy_values_deep)
static inline bool base_type_has_payload(base_type_t type) {
    return type == TYPE_ATOMSET || type == TYPE_BONDSET || type == TYPE_BOND_DELTA;
}

// Copies count values of a type whose values carry a payload of their own (sets of atoms and bonds, bond deltas), each
// with its payload allocated with alloc: a value copied as bytes would point into the payload of its source, which
// goes with whatever arena the source is in. False, and nothing copied, for a type of plain values.
static bool copy_values_deep(void* dst, const void* src, base_type_t type, size_t count, md_allocator_i* alloc) {
    switch (type) {
    case TYPE_ATOMSET:
        for (size_t i = 0; i < count; ++i) ((md_bitset_t*)dst)[i] = md_bitset_copy(((const md_bitset_t*)src)[i], alloc);
        return true;
    case TYPE_BONDSET:
        for (size_t i = 0; i < count; ++i) ((bondset_t*)dst)[i] = bondset_copy(((const bondset_t*)src)[i], alloc);
        return true;
    case TYPE_BOND_DELTA:
        for (size_t i = 0; i < count; ++i) {
            const bond_delta_value_t* s = (const bond_delta_value_t*)src + i;
            bond_delta_value_t* d = (bond_delta_value_t*)dst + i;
            d->touched = md_bitset_copy(s->touched, alloc);
            d->broken  = md_bitset_copy(s->broken,  alloc);
            d->formed  = bondset_copy(s->formed, alloc);
        }
        return true;
    default:
        return false;
    }
}

static void free_data(data_t* data, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);    
    if (data->ptr && data->size) {
        const int64_t num_elem = element_count(*data);
        if (data->type.base_type == TYPE_ATOMSET) {
            md_bitset_t* bs = data->ptr;
            for (int64_t i = 0; i < num_elem; ++i) {
                md_bitset_free(&bs[i], alloc);
            }
        } else if (data->type.base_type == TYPE_BONDSET) {
            bondset_t* bs = data->ptr;
            for (int64_t i = 0; i < num_elem; ++i) {
                bondset_free(&bs[i], alloc);
            }
        } else if (data->type.base_type == TYPE_BOND_DELTA) {
            bond_delta_value_t* d = data->ptr;
            for (int64_t i = 0; i < num_elem; ++i) {
                md_bitset_free(&d[i].touched, alloc);
                md_bitset_free(&d[i].broken, alloc);
                bondset_free(&d[i].formed, alloc);
            }
        }
        md_free(alloc, data->ptr, data->size);
        if (data->weight) {
            md_free(alloc, data->weight, data->size);
        }
    }
    data->ptr = 0;
    data->weight = 0;
    data->size = 0;
}

// The end of the domain of atom indices: masks built during an evaluation drop nothing below it
static inline uint32_t atom_domain_end(const eval_context_t* ctx) {
    return ctx->sys ? (uint32_t)ctx->sys->atom.count : MD_BITSET_MAX_INDEX;
}

// The allocator for the words of the bitsets written into d: the one d itself lives in, which allocate_data records
// (and anything that makes a data_t of bitsets otherwise sets). Script data lives in arenas (the IR's, the
// evaluation's, temp), so a value written over is not released: its words go with the arena.
// Every destination is expected to carry it. Should one not, a debug build stops here; a release build takes the
// allocator of what is produced (the IR's arena in a compilation, the evaluation's own otherwise), never the IR's in
// an evaluation.
static inline md_allocator_i* data_bitset_alloc(const data_t* d, const eval_context_t* ctx) {
    ASSERT(d && d->alloc && "a data_t holding bitsets without the allocator it lives in");
    if (d && d->alloc) return d->alloc;
    return ctx->alloc ? ctx->alloc : ctx->temp_alloc;
}

// False, with an error, if dst cannot hold src
static bool copy_data(data_t* dst, const data_t* src, const eval_context_t* ctx) {
    ASSERT(dst);
    ASSERT(src);
    // A smaller dst would be written past its end
    if (!INTERNAL_CHECK(ctx, ctx->op_token, dst->size >= src->size)) return false;

    dst->type = src->type;
    if (src->size == 0) {
        // An empty array, 'residue(x)' of an empty x: nothing to copy, and possibly no memory behind either
        return true;
    }
    if (!INTERNAL_CHECK(ctx, ctx->op_token, src->ptr && dst->ptr)) return false;

    // As many values as the bytes of src hold, whatever its type says: they are not read past
    const base_type_t base_type = dst->type.base_type;
    const size_t elem_size = base_type_element_byte_size(base_type);
    const size_t num_elem = elem_size ? MIN((size_t)element_count(*src), src->size / elem_size) : 0;
    if (!base_type_has_payload(base_type) || !copy_values_deep(dst->ptr, src->ptr, base_type, num_elem, data_bitset_alloc(dst, ctx))) {
        copy_data_block(dst, 0, src, 0, src->size);
    }
    return true;
}



static bool is_operator(ast_type_t type) {
    return (type == AST_ADD || type == AST_SUB || type == AST_MUL || type == AST_DIV || type == AST_INT_DIV || type == AST_UNARY_NEG ||
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
    case TYPE_ATOMSET: return "atomset";
    case TYPE_CONTACT: return "contact";
    case TYPE_DISTRIBUTION: return "distribution";
    case TYPE_PLANE: return "plane";
    case TYPE_PUCKER: return "pucker";
    case TYPE_BONDSET: return "bondset";
    case TYPE_BOND_DELTA: return "bond_delta";
    default: return "type out of range";
    }
}

#define LOG_WARNING(ir, tok, fmt, ...)  create_log_token_fmt(LOG_LEVEL_WARNING, ir, tok, fmt, ##__VA_ARGS__)
#define LOG_ERROR(ir, tok, fmt, ...)    create_log_token_fmt(LOG_LEVEL_ERROR, ir, tok, fmt, ##__VA_ARGS__)

#define LOG_WARNING_STR(ir, tok, str)   create_log_token_str(LOG_LEVEL_WARNING, ir, tok, str)
#define LOG_ERROR_STR(ir, tok, str)     create_log_token_str(LOG_LEVEL_ERROR, ir, tok, str)

static void create_log_token_str(log_level_t level, md_script_ir_t* ir, token_t token, str_t str) {
    ASSERT(ir);
    if (!ir->record_log) return;

    if (level >= LOG_LEVEL_ERROR) {
        // Stop after first error
        ir->record_log = false;
    }

    if (level == LOG_LEVEL_WARNING) {
        // The arguments of a procedure are checked again once it is resolved: a warning about one is made twice
        for (size_t i = 0; i < md_array_size(ir->warnings); ++i) {
            const md_log_token_t* w = &ir->warnings[i];
            if (w->range.beg == token.beg && w->range.end == token.end && str_eq(w->text, str)) return;
        }
    }

    md_log_token_t log_tok = {
        .range = {
            .beg = token.beg,
            .end = token.end,
        },
        .text = str_copy(str, ir->arena),
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
}

static void create_log_token_fmt(log_level_t level, md_script_ir_t* ir, token_t token, const char* format, ...) {
    char buf[1024];
    va_list args;
    va_start(args, format);
    int len = md_vsnprintf(buf, ARRAY_SIZE(buf), format, args);
    va_end(args);

    len = MAX(len, 0);

    str_t str = {buf, len};
    create_log_token_str(level, ir, token, str);
}

// Diagnostics raised where an eval_context_t is at hand: in the static checks of a compilation, and in the
// procedures, which run in compilations and evaluations alike. A compilation records them in its IR,
// which is how they reach the editor. An evaluation does not write to the IR: they go to its log.
static void eval_log_str(log_level_t level, const eval_context_t* ctx, token_t token, str_t str) {
    ASSERT(ctx);
    if (ctx->compile_ir) {
        create_log_token_str(level, ctx->compile_ir, token, str);
        return;
    }

    eval_log_t* log = ctx->log;
    if (!log) {
        MD_LOG_DEBUG("Script evaluation %s: "STR_FMT, level >= LOG_LEVEL_ERROR ? "error" : "warning", STR_ARG(str));
        return;
    }
    if (log->full) return;
    if (level >= LOG_LEVEL_ERROR) {
        log->full = true;
    }

    ASSERT(log->alloc);
    md_log_token_t log_tok = {
        .range = {
            .beg = token.beg,
            .end = token.end,
        },
        .text = str_copy(str, log->alloc),
    };
    if (level >= LOG_LEVEL_ERROR) {
        md_array_push(log->errors, log_tok, log->alloc);
    } else {
        md_array_push(log->warnings, log_tok, log->alloc);
    }
}

static void eval_log_fmt(log_level_t level, const eval_context_t* ctx, token_t token, const char* format, ...) {
    char buf[1024];
    va_list args;
    va_start(args, format);
    int len = md_vsnprintf(buf, ARRAY_SIZE(buf), format, args);
    va_end(args);

    len = CLAMP(len, 0, (int)ARRAY_SIZE(buf) - 1);

    str_t str = {buf, (size_t)len};
    eval_log_str(level, ctx, token, str);
}

#define EVAL_LOG_WARNING(ctx, tok, fmt, ...)  eval_log_fmt(LOG_LEVEL_WARNING, ctx, tok, fmt, ##__VA_ARGS__)
#define EVAL_LOG_ERROR(ctx, tok, fmt, ...)    eval_log_fmt(LOG_LEVEL_ERROR, ctx, tok, fmt, ##__VA_ARGS__)

#define EVAL_LOG_WARNING_STR(ctx, tok, str)   eval_log_str(LOG_LEVEL_WARNING, ctx, tok, str)
#define EVAL_LOG_ERROR_STR(ctx, tok, str)     eval_log_str(LOG_LEVEL_ERROR, ctx, tok, str)

// Reports a broken invariant (see INTERNAL_CHECK)
static void internal_error(const eval_context_t* ctx, token_t tok, const char* what, const char* file, int line) {
    // The name of the file is enough to find the line, the path it was built from is noise
    for (const char* c = file; *c; ++c) {
        if (*c == '/' || *c == '\\') file = c + 1;
    }
    MD_LOG_ERROR("Script: internal error, %s (%s:%d)", what, file, line);
    EVAL_LOG_ERROR(ctx, tok, "Internal error: %s (%s:%d). This is a bug, please report it", what, file, line);
}

// Whether str matches pattern, in which '*' stands for any run of characters (the empty one too) and every other
// character for itself, case and all. Several '*' anywhere: "C*", "*A", "*A*", "C*A", "C*1*". What a '*' takes is
// tried shortest first and widened when the rest does not match, which backtracks to the last '*' only: an earlier
// one never needs more, since the later one can take whatever it would have.
static bool glob_match(str_t pattern, str_t str) {
    size_t p = 0, s = 0;
    size_t star = SIZE_MAX;     // the last '*' seen in the pattern
    size_t star_s = 0;          // where in str what it takes ends, so far
    while (s < str.len) {
        if (p < pattern.len && pattern.ptr[p] == '*') {
            star = p++;
            star_s = s;
        } else if (p < pattern.len && pattern.ptr[p] == str.ptr[s]) {
            ++p;
            ++s;
        } else if (star != SIZE_MAX) {
            // The last '*' takes one more character
            p = star + 1;
            s = ++star_s;
        } else {
            return false;
        }
    }
    while (p < pattern.len && pattern.ptr[p] == '*') ++p;
    return p == pattern.len;
}

static inline bool str_has_glob(str_t str) {
    return str_find_char(NULL, str, '*');
}

// Whether path matches pattern segment by segment: as many segments, each matching its own. A '*' does not reach
// across a '/', so 'edr/*' is what the group edr holds and not what the groups below it hold.
static bool glob_match_path(str_t pattern, str_t path) {
    for (;;) {
        size_t pl = 0, sl = 0;
        const bool p_more = str_find_char(&pl, pattern, '/');
        const bool s_more = str_find_char(&sl, path, '/');
        if (p_more != s_more) return false;
        if (!p_more) return glob_match(pattern, path);
        if (!glob_match(str_substr(pattern, 0, pl), str_substr(path, 0, sl))) return false;
        pattern = str_substr(pattern, pl + 1, SIZE_MAX);
        path    = str_substr(path, sl + 1, SIZE_MAX);
    }
}

// The segments of path that a '*' of pattern (which it matches, glob_match_path) stood in, joined by '/': the part of
// what a pattern matches that differs from one match to the next. "coul_sr_sol" of "edr/coul_sr_*".
static str_t glob_varying_segments(str_t pattern, str_t path, md_allocator_i* alloc) {
    md_str_writer_t w = md_str_writer_create(alloc);
    for (;;) {
        size_t pl = 0, sl = 0;
        const bool more = str_find_char(&pl, pattern, '/') && str_find_char(&sl, path, '/');
        const str_t p_seg = more ? str_substr(pattern, 0, pl) : pattern;
        const str_t s_seg = more ? str_substr(path, 0, sl) : path;
        if (str_has_glob(p_seg)) {
            if (md_str_writer_len(&w)) md_str_writer_push_char(&w, '/');
            md_str_writer_push_str(&w, s_seg);
        }
        if (!more) break;
        pattern = str_substr(pattern, pl + 1, SIZE_MAX);
        path    = str_substr(path, sl + 1, SIZE_MAX);
    }
    return md_str_writer_to_str(&w);
}

// The last segment of a path
static str_t path_leaf(str_t path) {
    size_t loc;
    return str_rfind_char(&loc, path, '/') ? str_substr(path, loc + 1, SIZE_MAX) : path;
}

// Include all procedures, operators and defines
#ifdef MD_SCRIPT_TEST_HOOKS
// Lets the unittests count how often a procedure is evaluated, i.e. called with a destination
static struct {
    int (*proc_ptr)(data_t*, data_t[], eval_context_t*);
    size_t count;
} test_hook_proc_eval;

// Number of contact searches (md_contact_query_eval) contacts() has run
static size_t test_hook_contact_searches;
static size_t test_hook_contact_local_prepares;    // Evaluations which found no query prepared for their call site
#endif

#include "md_script_functions.inl"

// The volumes of a type: 1 for a volume, float[128][128][128] (or [1][128][128][128]), n for an array of them,
// float[n][128][128][128], which is what 'sdf(...) in residue(1:n)' or '{sdf(...), sdf(...)}' is. 0 for anything else.
static inline int volume_count(type_info_t ti) {
    if (ti.dim[0] == -1) return 0;
    if (is_type_directly_compatible(ti, (type_info_t)TI_VOLUME)) return 1;
    if (ti.base_type == TYPE_FLOAT && dim_ndims(ti.dim) == 4 && ti.dim[0] > 0 &&
        ti.dim[1] == MD_VOL_DIM && ti.dim[2] == MD_VOL_DIM && ti.dim[3] == MD_VOL_DIM) {
        return ti.dim[0];
    }
    return 0;
}


// ############################
// ###   HELPER FUNCTIONS   ###
// ############################

// Returns if the value type can be operated on using logical operators
static inline bool is_value_type_logical_operator_compatible(base_type_t type) {
    return type == TYPE_BOOL || type == TYPE_ATOMSET || type == TYPE_BONDSET;
}

static inline bool is_identifier_static_procedure(str_t ident) {
    if (str_eq(ident, STR_LIT("attr"))) {
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
    node->data.unit[0] = md_unit_none();
    node->data.unit[1] = md_unit_none();
    node->index = (uint32_t)md_array_size(ir->nodes);
    md_array_push(ir->nodes, node, ir->arena);
    return node;
}

// Copies the contents of src into dst, a node of its own: dst keeps its index
static void copy_node_contents(ast_node_t* dst, const ast_node_t* src) {
    const uint32_t index = dst->index;
    MEMCPY(dst, src, sizeof(ast_node_t));
    dst->index = index;
}

static void fix_precedence(ast_node_t** node) {
    // We parse from left to right thus unbalancing the tree in its right child node (child[1])
    // Therefore we only need to fix potential precedence issues down this path
    // This is unfortunately not the case with keyword in, since that can be unbalanced down in the left subtree...
    // Therefore we check both

    if (!node) return;
    ast_node_t* parent = *node;
    if (!parent) return;

    if (operator_binary((*node)->type)) {
        if (md_array_size((*node)->children) != 2) {
            MD_LOG_ERROR("Invalid number of children in binary node");
            return;
        }

        ast_node_t* child = parent->children[1];
        if (md_array_size(child->children) == 0) return;

        int prec_parent   = operator_precedence(parent->type);
        int prec_child    = operator_precedence(child->type);

        if (prec_child > prec_parent || (prec_child == prec_parent && operator_associativity(child->type) == ASSOCIATIVITY_LEFT_TO_RIGHT)) {
            ast_node_t* grand_child = child->children[0];
            child->children[0]  = parent;
            parent->children[1] = grand_child;
            *node = child;
            fix_precedence(&((*node)->children[0]));
        }
    } else {
        const size_t num_children = md_array_size(parent->children);
        for (size_t i = 0; i < num_children; ++i) {
            ast_node_t* child = parent->children[i];
            if (!child || md_array_size(child->children) == 0) continue;
            if (operator_precedence(child->type) > operator_precedence(parent->type)) {
                ast_node_t* grand_child = child->children[0];
                child->children[0]  = parent;
                parent->children[i] = grand_child;
                *node = child;
            }
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
        if (is_absent_type(arg_types[j]) || type_info_equal(arg_types[j], proc->arg_type[j])) {
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

static procedure_match_result_t find_procedure_supporting_arg_types_in_candidates(str_t name, const type_info_t arg_types[], size_t num_arg_types, procedure_t candidates[], size_t num_cantidates, bool allow_implicit_conversions) {
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
    case AST_INT_DIV: name = STR_LIT("//"); break;
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

// The same lookup for an evaluation, which only reads the IR
static const identifier_t* find_ir_identifier(const md_script_ir_t* ir, str_t name) {
    if (ir) {
        for (size_t i = 0; i < md_array_size(ir->identifiers); ++i) {
            if (str_eq(name, ir->identifiers[i].name)) {
                return &ir->identifiers[i];
            }
        }
    }
    return NULL;
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

    md_array_push(ir->identifiers, ident, ir->arena);
    return md_array_last(ir->identifiers);
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
            if (i == len) {
                // The comment runs to the end: there is nothing after it, and buf[len] is not ours to read
                tokenizer->cur = len;
                token.type = TOKEN_END;
                token.line_beg = line;
                token.line_end = line;
                token.col_beg  = len - tokenizer->line_offset;
                token.col_end  = len - tokenizer->line_offset;
                token.beg = len;
                token.end = len;
                break;
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
                else if (buf[j] == ';' || (buf[j] == '\r' && j + 1 < len && buf[j+1] == '\n')) { // We do not want to leak return carry '\r' into the text since that will malform any output
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
                {"==", TOKEN_EQ},
                {"//", TOKEN_INT_DIV},
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
    size_t len = md_snprintf(buf, cap, "%s", get_value_type_str(info.base_type));

    for (size_t i = 0; i < MAX_NUM_DIMS; ++i) {
        if (info.dim[i] == -1) {
            len += md_snprintf(buf + len, cap - MIN(len, cap), "[..]");
        }
        else if (info.dim[i] >= 1) {
            len += md_snprintf(buf + len, cap - MIN(len, cap), "[%i]", (int)info.dim[i]);
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

// NOTE: len accumulates snprintf return values, so it is the length the output WOULD have had and
// can exceed cap. buf + MIN(len, cap) keeps the pointer inside the buffer when that happens; the
// size argument is then zero and nothing is written. Callers that need the number of bytes
// actually present must clamp the result to cap - 1.
#define PRINT(...) len += md_snprintf(buf + MIN(len, cap), cap - MIN(len, cap), ##__VA_ARGS__)

static size_t print_bitfield(char* buf, size_t cap, const md_bitset_t* set) {
    ASSERT(set);
    if (cap <= 0) return 0;

    const size_t max_bits = 64;
    const size_t num_bits = MIN(set->end - set->beg, max_bits);
    size_t len = 0;

    PRINT("<%i,%i>[", set->beg, set->end);
    for (size_t i = 0; i < num_bits; ++i) {
        PRINT("%i", md_bitset_test(*set, set->beg + (uint32_t)i) ? 1 : 0);
    }
    PRINT("%s]", num_bits == max_bits ? "..." : "");

    return len;
}

static size_t print_data_value(char* buf, size_t cap, data_t data) {
    size_t len = 0;
    if (is_variable_length(data.type)) return 0;
    if (cap <= 0) return 0;

    if (data.ptr) {
        // An element of an array of arrays ([1][3]) is the inner array ([3]), not its first value
        while (data.type.dim[0] == 1 && data.type.dim[1] > 1) {
            dim_shift_left(data.type.dim);
        }
        if (is_array(data.type) && data.type.base_type != TYPE_ATOMSET) {
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
                    .weight = data_weight_at(&data, stride * i),
                };
                len += print_data_value(buf + len, cap - MIN(len, cap), elem_data);
                if (i < arr_len - 1) PRINT(",");
            }
            PRINT("}");
        } else {
            switch(data.type.base_type) {
            case TYPE_ATOMSET:
                //print_bitfield(file, (md_bitset_t*)data.ptr);
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
                if (rng.step > 1) {
                    PRINT(":%i:", rng.step);
                } else {
                    PRINT(":");
                }
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
            case TYPE_DISTRIBUTION:
            {
                const float v = *(const float*)data.ptr;
                const float w = data.weight ? *data.weight : 1.0f;
                PRINT("%.3f", w != 0.0f ? v / w : 0.0f);
                break;
            }
            case TYPE_PLANE:
            {
                const plane_t* p = (const plane_t*)data.ptr;
                PRINT("plane(centroid (%.2f, %.2f, %.2f), normal (%.2f, %.2f, %.2f), radius %.2f)",
                    p->centroid.x, p->centroid.y, p->centroid.z, p->normal.x, p->normal.y, p->normal.z, p->radius);
                break;
            }
            case TYPE_CONTACT:
            {
                const md_contact_set_t* set = (const md_contact_set_t*)data.ptr;
                PRINT("%zu pairs", set->count);
                break;
            }
            case TYPE_BONDSET:
            {
                const bondset_t* b = (const bondset_t*)data.ptr;
                PRINT("%zu bonds", (size_t)md_bitset_count(b->bonds) + b->num_extra);
                break;
            }
            case TYPE_BOND_DELTA:
            {
                const bond_delta_value_t* d = (const bond_delta_value_t*)data.ptr;
                PRINT("%u formed, %zu broken", d->formed.num_extra, (size_t)md_bitset_count(d->broken));
                break;
            }
            case TYPE_PUCKER:
            {
                const pucker_rec_t* p = (const pucker_rec_t*)data.ptr;
                PRINT("pucker(%s, Q %.2f)", p->name, p->cp.amplitude);
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
    if (data.ptr && !md_unit_is_none(data.unit[1])) {
        if (len + 1 < cap) buf[len] = ' ';
        len += 1;
        len += md_unit_print(buf + MIN(len, cap), cap - MIN(len, cap), data.unit[1]);
    }
    return len;
}

#undef PRINT

// ### VALUE SUMMARIES ###
// What hovering a value shows: a few lines which say what it is, rather than every element of it. A selection is its
// number of atoms, a contact set its pairs, a plane its centroid and normal, and a long array its first values.
// subidx picks one element of an array (-1 for all), as when stepping through them with the mouse wheel.

#define SUMMARY_MAX_VALUES 8

static void summarize_plane(md_writer_i* w, const plane_t* p) {
    md_writef(w, "centroid (%.2f, %.2f, %.2f) \xC3\x85\n", p->centroid.x, p->centroid.y, p->centroid.z);
    md_writef(w, "normal (%.3f, %.3f, %.3f)\n", p->normal.x, p->normal.y, p->normal.z);
    md_writef(w, "radius %.2f \xC3\x85, spread %.2f x %.2f \xC3\x85, deviation %.3f \xC3\x85 (%d atoms)", p->radius, p->extent[0], p->extent[1], p->deviation, p->count);
}

static void summarize_pucker(md_writer_i* w, const pucker_rec_t* p) {
    if (p->name[0]) md_writef(w, "%s\n", p->name);
    md_writef(w, "Q %.3f \xC3\x85", p->cp.amplitude);
    if (p->cp.size == 6) md_writef(w, ", theta %.1f\xC2\xB0", p->cp.theta);
    if (p->anchored || p->has_pseudo) {
        md_writef(w, ", phi %.1f\xC2\xB0", p->cp.phi2);
        if (p->cp.size == 7) md_writef(w, ", phi3 %.1f\xC2\xB0", p->cp.phi3);
    }
    if (p->has_pseudo) md_writef(w, "\nP %.1f\xC2\xB0, tau_m %.1f\xC2\xB0", p->pseudo_phase, p->pseudo_amplitude);
    if (!p->anchored && !p->has_pseudo) md_writef(w, "\n(the ring has no first atom, so no phase)");
}

static void summarize_data(md_writer_i* w, data_t data, int subidx) {
    if (!data.ptr) {
        return;
    }
    const int len = (int)element_count(data);
    const bool pick = subidx >= 0 && subidx < len;
    const size_t stride = type_info_element_byte_stride(data.type);

    switch (data.type.base_type) {
    case TYPE_ATOMSET: {
        const md_bitset_t* bf = (const md_bitset_t*)data.ptr;
        if (pick) {
            const size_t n = md_bitset_count(bf[subidx]);
            md_writef(w, "selection %d of %d: %zu atom%s", subidx + 1, len, n, n == 1 ? "" : "s");
        } else if (len == 1) {
            const size_t n = md_bitset_count(bf[0]);
            md_writef(w, "%zu atom%s", n, n == 1 ? "" : "s");
        } else {
            md_bitmask_t all = md_bitmask_create(0, MD_BITSET_MAX_INDEX, md_get_heap_allocator());
            size_t empty = 0;
            for (int i = 0; i < len; ++i) {
                md_bitmask_or(&all, bf[i]);
                empty += md_bitset_empty(bf[i]);
            }
            md_writef(w, "%d selections, %zu atoms in all", len, (size_t)md_bitmask_count(&all));
            if (empty) md_writef(w, " (%zu empty)", empty);
            md_bitmask_free(&all);
        }
        return;
    }
    case TYPE_BONDSET: {
        const bondset_t* bs = (const bondset_t*)data.ptr;
        if (pick || len == 1) {
            const bondset_t* b = &bs[pick ? subidx : 0];
            const size_t n = (size_t)md_bitset_count(b->bonds) + b->num_extra;
            if (pick) md_writef(w, "bond set %d of %d: ", subidx + 1, len);
            md_writef(w, "%zu bond%s", n, n == 1 ? "" : "s");
            if (b->num_extra) md_writef(w, " (%u computed for the frame)", b->num_extra);
        } else {
            size_t n = 0;
            for (int i = 0; i < len; ++i) n += (size_t)md_bitset_count(bs[i].bonds) + bs[i].num_extra;
            md_writef(w, "%d bond sets, %zu bonds", len, n);
        }
        return;
    }
    case TYPE_BOND_DELTA: {
        const bond_delta_value_t* d = (const bond_delta_value_t*)data.ptr;
        md_writef(w, "%u bond%s formed, %zu broken", d->formed.num_extra, d->formed.num_extra == 1 ? "" : "s", (size_t)md_bitset_count(d->broken));
        md_writef(w, "\nof %zu bonds in the scope", (size_t)md_bitset_count(d->touched));
        return;
    }
    case TYPE_CONTACT: {
        const md_contact_set_t* set = (const md_contact_set_t*)data.ptr;
        const bool self = set->flags & MD_CONTACT_FLAG_SELF;
        if (!set->i && set->count == 0) {
            // The shape only (the static check)
            if (self) md_writef(w, "pairs among %u groups", set->num_a);
            else      md_writef(w, "pairs between %u and %u groups", set->num_a, set->num_b);
            return;
        }
        double atom_pairs = 0;
        for (size_t k = 0; k < set->count; ++k) atom_pairs += set->atom_pairs ? set->atom_pairs[k] : 0;
        md_writef(w, "%zu pair%s of groups in contact", set->count, set->count == 1 ? "" : "s");
        if (self) md_writef(w, ", among %u groups", set->num_a);
        else      md_writef(w, ", between %u and %u groups", set->num_a, set->num_b);
        md_writef(w, "\n%.0f atom pairs", atom_pairs);
        if (set->count > 0 && set->d_min) {
            float d = set->d_min[0];
            for (size_t k = 1; k < set->count; ++k) d = MIN(d, set->d_min[k]);
            md_writef(w, ", closest %.2f \xC3\x85", d);
        }
        return;
    }
    case TYPE_PLANE: {
        const plane_t* p = (const plane_t*)data.ptr;
        if (pick) {
            md_writef(w, "plane %d of %d\n", subidx + 1, len);
            summarize_plane(w, &p[subidx]);
        } else if (len == 1) {
            summarize_plane(w, &p[0]);
        } else {
            float r = 0, dev = 0;
            for (int i = 0; i < len; ++i) {
                r   += p[i].radius;
                dev  = MAX(dev, p[i].deviation);
            }
            md_writef(w, "%d planes, mean radius %.2f \xC3\x85, largest deviation %.3f \xC3\x85\n", len, len ? r / len : 0.0f, dev);
            md_writef(w, "plane 1:\n");
            summarize_plane(w, &p[0]);
        }
        return;
    }
    case TYPE_PUCKER: {
        const pucker_rec_t* p = (const pucker_rec_t*)data.ptr;
        if (pick) {
            md_writef(w, "ring %d of %d\n", subidx + 1, len);
            summarize_pucker(w, &p[subidx]);
        } else if (len == 1) {
            summarize_pucker(w, &p[0]);
        } else {
            // How many rings are in each conformation, the most common first
            md_writef(w, "%d rings:", len);
            md_temp_scope_t temp = md_temp_begin_avoid(md_writer_allocator(w));
            bool* done = md_temp_alloc(temp, sizeof(bool) * len);
            MEMSET(done, 0, sizeof(bool) * len);
            int shown = 0;
            for (int i = 0; i < len && shown < 6; ++i) {
                int best = -1, best_count = 0;
                for (int k = 0; k < len; ++k) {
                    if (done[k]) continue;
                    int c = 0;
                    for (int m = k; m < len; ++m) c += !done[m] && !strcmp(p[m].name, p[k].name);
                    if (c > best_count) { best = k; best_count = c; }
                }
                if (best < 0) break;
                for (int m = best; m < len; ++m) if (!strcmp(p[m].name, p[best].name)) done[m] = true;
                md_writef(w, "%s %d %s", shown ? "," : "", best_count, p[best].name[0] ? p[best].name : "unnamed");
                shown += 1;
            }
            md_writef(w, "\nring 1:\n");
            summarize_pucker(w, &p[0]);
            md_temp_end(temp);
        }
        return;
    }
    case TYPE_DISTRIBUTION: {
        const int bins = data.type.dim[1];
        if (len == 1) md_writef(w, "distribution of %d bins", bins);
        else          md_writef(w, "%d distributions of %d bins", len, bins);
        if (data.value_range.end > data.value_range.beg) {
            md_writef(w, " over %.2f to %.2f", data.value_range.beg, data.value_range.end);
            char unit[64];
            const size_t n = md_unit_print(unit, sizeof(unit), data.unit[0]);
            if (n) md_writef(w, " %.*s", (int)MIN(n, sizeof(unit) - 1), unit);
        }
        return;
    }
    default:
        break;
    }

    // Plain values: printed, the first few of a long array. Angles are read in degrees.
    md_temp_scope_t temp = md_temp_begin_avoid(md_writer_allocator(w));
    data_t view = data;
    if (data.type.base_type == TYPE_FLOAT && md_unit_equal(data.unit[1], md_unit_radian())) {
        const size_t total = type_info_total_element_count(data.type);
        float* deg = md_temp_alloc(temp, sizeof(float) * MAX(total, 1));
        for (size_t i = 0; i < total; ++i) deg[i] = RAD_TO_DEG(((const float*)data.ptr)[i]);
        view.ptr = deg;
        view.unit[1] = md_unit_degree();
        data = view;
    }
    char buf[512];
    if (pick && is_array(data.type)) {
        view.ptr = (char*)data.ptr + stride * subidx;
        view.weight = data_weight_at(&data, stride * subidx);
        view.type.dim[0] = 1;
        md_writef(w, "element %d of %d: ", subidx + 1, len);
    } else if (is_array(data.type) && len > SUMMARY_MAX_VALUES) {
        view.type.dim[0] = SUMMARY_MAX_VALUES;
    }
    const size_t n = print_data_value(buf, sizeof(buf), view);
    md_write(w, buf, MIN(n, sizeof(buf) - 1));
    if (!pick && is_array(data.type) && len > SUMMARY_MAX_VALUES) {
        md_writef(w, " ... (%d values)", len);
    }
    md_temp_end(temp);
}

#undef SUMMARY_MAX_VALUES

#if DEBUG


static void print_label(md_file_t  file, const ast_node_t* node) {
    char buf[1024];
    switch(node->type) {
    case AST_CONSTANT_VALUE:
        print_data_value(buf, sizeof(buf), node->data);
        md_file_printf(file, "%s", buf);
        break;
    case AST_PROC_CALL:
        if (node->proc) {
            md_file_printf(file, ""STR_FMT"", (int)node->proc->name.len, node->proc->name.ptr);
        }
        else {
            md_file_printf(file, "NULL");
        }
        md_file_printf(file, "(proc call)");
        break;
    case AST_ASSIGNMENT:
        md_file_printf(file, "assignment =");
        break;
    case AST_CAST:
        ASSERT(node->children);
        print_type_info(buf, ARRAY_SIZE(buf), node->children[0]->data.type);
        md_file_printf(file, "cast (%s) ->", buf);
        break;
    case AST_IDENTIFIER:
        md_file_printf(file, ""STR_FMT"", (int)node->ident.len, node->ident.ptr);
        md_file_printf(file, " (identifier)");
        break;
    case AST_ARRAY:
        md_file_printf(file, "array");
        break;
    case AST_ARRAY_SUBSCRIPT:
        md_file_printf(file, "array subscript");
        break;
    case AST_AND:
    case AST_OR:
    case AST_NOT:
        md_file_printf(file, "logical op %s", get_token_type_str(node->token.type));
        break;
    case AST_ADD:
    case AST_SUB:
    case AST_MUL:
    case AST_DIV:
        md_file_printf(file, "binary op %s", get_token_type_str(node->token.type));
        break;
    default:
        md_file_printf(file, "%s", get_token_type_str(node->token.type));
        break;
    }

    print_type_info(buf, sizeof(buf), node->data.type);
    md_file_printf(file, " -> %s", buf);

    if (node->flags) {
        md_file_printf(file, " {");
        if (node->flags & FLAG_DYNAMIC) md_file_printf(file, "D");
        md_file_printf(file, "}");
    }

    if (node->type != AST_ASSIGNMENT) {
        if (node->data.ptr) {
            md_file_printf(file, " = ");
            print_data_value(buf, sizeof(buf), node->data);
            md_file_printf(file, "%s", buf);
        }
    }
}


static inline void indent(md_file_t  file, int amount) {
    for (int i = 0; i < amount; ++i) {
        md_file_printf(file, "\t");
    }
}

static void print_node(md_file_t  file, const ast_node_t* node, int depth) {
    if (!node || !md_file_valid(file)) return;

    indent(file, depth);
    md_file_printf(file, "{\"name\": \"");
    print_label(file, node);
    md_file_printf(file, "\",\n");

    const size_t num_children = md_array_size(node->children);
    if (num_children) {
        indent(file, depth);
        md_file_printf(file, "\"children\": [\n");
        for (size_t i = 0; i < num_children; ++i) {
            print_node(file, node->children[i], depth + 1);
        }
        indent(file, depth);
        md_file_printf(file, "]\n");
    }

    indent(file, depth);
    md_file_printf(file, "},\n");
}

static void print_expr(md_file_t  file, str_t str) {
    // We need to add escape character to quotation marks, since they are not representable within a json string...
    char buf[1024] = {0};
    int len = 0;
    for (const char* c = str.ptr; c != str.ptr + str.len; ++c) {
        if (*c == '\"') {
            buf[len++] = '\\';
            buf[len++] = *c;
        }
        else if (*c == '\r') {

        }
        else if (*c == '\n') {
            buf[len++] = '\\';
            buf[len++] = 'n';
        }
        else {
            buf[len++] = *c;
        }
    }
    md_file_printf(file, "\""STR_FMT"\"\n", len, buf);
}

static void ast_to_json(const ast_node_t* node, str_t filename) {
    if (!node) return;
    md_file_t file = {0};
    if (md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        md_file_printf(file, "var treeData = [\n");
        md_file_printf(file, "{\n\"node\" :\n");
        print_node(file, node, 1);
        md_file_printf(file, "\"flags\" : %u,\n", node->flags);
        md_file_printf(file, "\"expr\" : ");
        print_expr(file, node->token.str);
        md_file_printf(file, "},\n");
        md_file_printf(file, "];");
        md_file_close(&file);
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

// True if the next two tokens are 'identifier =', i.e. a named argument. Consumes nothing.
// '==' is a token of its own (TOKEN_EQ), so a comparison is never mistaken for one.
static bool tokenizer_next_is_named_argument(const tokenizer_t* tokenizer) {
    tokenizer_t peek = *tokenizer;
    if (tokenizer_consume_next(&peek).type != TOKEN_IDENT) return false;
    return tokenizer_peek_next(&peek).type == '=';
}

// out_named is optional. It is only supplied for the argument list of a procedure call, which is the one place
// where named arguments are valid. When supplied, it receives one entry per argument (parallel to the result) if
// any argument was named, and NULL otherwise.
static ast_node_t** parse_comma_separated_arguments_until_token(parse_context_t* ctx, token_type_t token_type, md_array(named_arg_t)* out_named) {
    ast_node_t** args = 0;
    md_array(named_arg_t) named = 0;
    bool any_named = false;
    token_t next = {0};

    if (out_named) *out_named = NULL;

    while (next = tokenizer_peek_next(ctx->tokenizer), next.type != TOKEN_END) {
        if (next.type == token_type) goto done;
        if (next.type == ',') {
            LOG_ERROR(ctx->ir, next, "Empty argument in argument list");
            return NULL;
        }

        named_arg_t arg_name = {0};
        // The name has to be recognized before the expression is parsed: a name which happens to match a
        // procedure (e.g. 'count') would otherwise be parsed as a call to it.
        if (tokenizer_next_is_named_argument(ctx->tokenizer)) {
            token_t name_tok = tokenizer_consume_next(ctx->tokenizer);
            tokenizer_consume_next(ctx->tokenizer); // '='
            if (!out_named) {
                LOG_ERROR(ctx->ir, name_tok, "Named arguments ('"STR_FMT"=') are only valid in procedure calls", STR_ARG(name_tok.str));
                return NULL;
            }
            for (size_t i = 0; i < md_array_size(named); ++i) {
                if (str_eq(named[i].name, name_tok.str)) {
                    LOG_ERROR(ctx->ir, name_tok, "The argument '"STR_FMT"' is given more than once", STR_ARG(name_tok.str));
                    return NULL;
                }
            }
            arg_name.name  = name_tok.str;
            arg_name.token = name_tok;
            any_named = true;
        } else if (any_named) {
            LOG_ERROR(ctx->ir, next, "A positional argument cannot follow a named argument");
            return NULL;
        }

        ctx->node = 0;
        ast_node_t* arg = parse_expression(ctx);
        next = tokenizer_peek_next(ctx->tokenizer);
        if (arg && (next.type == ',' || next.type == token_type)) {
            md_array_push(args, arg, ctx->temp_alloc);
            if (!str_empty(arg_name.name)) {
                arg_name.token = concat_tokens(arg_name.token, arg->token);
            }
            md_array_push(named, arg_name, ctx->temp_alloc);
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
    if (out_named && any_named) {
        *out_named = named;
    }
    return args;
}

static ast_node_t* parse_procedure_call(parse_context_t* ctx, token_t token) {
    // Here we know atleast that the marker refers to something which is listed as a procedure!
    token_t next = tokenizer_peek_next(ctx->tokenizer);
    ast_node_t* node = 0;
    
    if (next.type == '(') {
        tokenizer_consume_next(ctx->tokenizer); // '('

        md_array(named_arg_t) named = 0;
        ast_node_t **args = parse_comma_separated_arguments_until_token(ctx, ')', &named);
        next = tokenizer_consume_next(ctx->tokenizer);
        if (expect_token_type(ctx->ir, next, ')')) {
            node = create_node(ctx->ir, AST_PROC_CALL, token);
            node->ident = str_copy(token.str, ctx->ir->arena);
            const size_t num_args = md_array_size(args);
            if (num_args) {
                md_array_push_array(node->children, args, num_args, ctx->ir->arena);
            }
            if (named) {
                ASSERT(md_array_size(named) == num_args);
                md_array_push_array(node->named_args, named, num_args, ctx->ir->arena);
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

// attr, flatten and transpose are parsed as procedure calls but become nodes of their own, which are never
// bound against a signature. Named arguments would silently be ignored there, so they are rejected up front.
static bool reject_named_arguments(parse_context_t* ctx, const ast_node_t* node) {
    for (size_t i = 0; i < md_array_size(node->named_args); ++i) {
        if (!str_empty(node->named_args[i].name)) {
            LOG_ERROR(ctx->ir, node->named_args[i].token, "'"STR_FMT"' does not accept named arguments", STR_ARG(node->ident));
            return false;
        }
    }
    return true;
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
        // Gone: data from files is loaded into the system's attributes and read with attr().
        LOG_ERROR(ctx->ir, token, "import() has been removed: load the file alongside the trajectory, and read what it holds with attr(\"<path>\")");
        return NULL;
    } else if (str_eq(ident, STR_LIT("attr"))) {
        node = parse_procedure_call(ctx, token);
        if (node) {
            if (!reject_named_arguments(ctx, node)) return NULL;
            node->type = AST_ATTRIBUTE;
        }
    } else if (str_eq(ident, STR_LIT("flatten"))) {
        node = parse_procedure_call(ctx, token);
        if (node) {
            if (!reject_named_arguments(ctx, node)) return NULL;
            node->type = AST_FLATTEN;
        }
    } else if (str_eq(ident, STR_LIT("transpose"))) {
        node = parse_procedure_call(ctx, token);
        if (node) {
            if (!reject_named_arguments(ctx, node)) return NULL;
            node->type = AST_TRANSPOSE;
        }
    } else if (is_identifier_procedure(ident)) {
        // The identifier has matched some procedure name
        node = parse_procedure_call(ctx, token);
    } else if ((c = find_constant(ident)) != 0) {
        node = create_node(ctx->ir, AST_CONSTANT_VALUE, token);
        node->value._float = c->value;
        node->data.unit[0] = md_unit_none();
        node->data.unit[1] = c->unit;
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
    if (!lhs) {
        LOG_ERROR(ctx->ir, token, "Invalid comparison, left operand of '"STR_FMT"' is missing.", STR_ARG(token.str));
        return NULL;
    }

    ctx->node = 0; // The right hand side is parsed on its own, otherwise 'a < b' is mistaken for two adjacent identifiers
    ast_node_t* rhs = parse_expression(ctx);
    ast_node_t* node = 0;

    if (!rhs) {
        LOG_ERROR(ctx->ir, token, "Invalid comparison, right operand of '"STR_FMT"' is missing.", STR_ARG(token.str));
    } else {
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
    ASSERT(token.type == '-' || token.type == '+' || token.type == '*' || token.type == '/' || token.type == TOKEN_INT_DIV);

    ast_node_t* lhs = ctx->node;
    // Only '-' has a unary form. Any other operator without a left operand, e.g. the '*' of 'f(1, *2)', would get a
    // NULL child, which nothing after the parser expects
    if (!lhs && token.type != '-') {
        LOG_ERROR(ctx->ir, token, "Invalid arithmetic expression, left operand of '"STR_FMT"' is missing.", STR_ARG(token.str));
        return NULL;
    }

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
            case TOKEN_INT_DIV: type = AST_INT_DIV; break;
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
    int    step = 1;

        if (is_number(token.type)) {
            beg = parse_float(token.str);
            expect_token_type(ctx->ir, tokenizer_consume_next(ctx->tokenizer), ':');
        } else { // == ':'
            beg = -FLT_MAX;
        }
        
        next = tokenizer_peek_next(ctx->tokenizer);
        if (is_number(next.type)) {
            token = tokenizer_consume_next(ctx->tokenizer);
            double mid = parse_float(token.str);
            bool mid_is_float = (token.type == TOKEN_FLOAT);
            // Expand token span
            node->token = concat_tokens(node->token, next);

            // Peek to see if we have a second ':' indicating stride syntax start:step:end
            token_t maybe_colon = tokenizer_peek_next(ctx->tokenizer);
            if (maybe_colon.type == ':') {
                // Consume ':'
                (void)tokenizer_consume_next(ctx->tokenizer);
                node->token = concat_tokens(node->token, maybe_colon);
                // Next must be a number (end)
                token_t end_tok = tokenizer_peek_next(ctx->tokenizer);
                if (is_number(end_tok.type)) {
                    end = parse_float(tokenizer_consume_next(ctx->tokenizer).str);
                    node->token = concat_tokens(node->token, end_tok);
                } else {
                    // Missing end, treat as open range
                    end = FLT_MAX;
                }

                // Only integer stride is supported
                if (mid_is_float) {
                    // If stride is float, fall back to error later in static check by marking as float range
                    is_float = true; // Force FRANGE to trigger error or ignore stride
                } else {
                    step = (int)mid;
                }
            } else {
                // No second colon -> simple range beg:end (mid is end)
                end = mid;
                is_float |= mid_is_float;
            }
        } else {
            end = FLT_MAX;
        }

        if (is_float) {
            node->data.type = (type_info_t)TI_FRANGE;
            node->value._frange.beg = (float)beg;
            node->value._frange.end = (float)end;
        } else {
            node->data.type = (type_info_t)TI_IRANGE;
            node->value._irange.beg  = (beg == -FLT_MAX) ?  INT32_MIN : (int32_t)beg;
            node->value._irange.end  = (end ==  FLT_MAX) ?  INT32_MAX : (int32_t)end;
            // Preserve explicit stride as parsed; validation (step > 0) is handled in static_check_constant_value
            node->value._irange.step = step;
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

        ast_node_t** elements = parse_comma_separated_arguments_until_token(ctx, ']', NULL);
        token_t next = tokenizer_consume_next(ctx->tokenizer);
        if (expect_token_type(ctx->ir, next, ']')) {
            const size_t num_elements = md_array_size(elements);
            if (num_elements) {

                ast_node_t* node_copy = create_node(ctx->ir, node->type, token);
                copy_node_contents(node_copy, node);

                const uint32_t index = node->index;
                MEMSET(node, 0, sizeof(ast_node_t));
                node->index = index;
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

    ast_node_t** elements = parse_comma_separated_arguments_until_token(ctx, '}', NULL);
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


ast_node_t* parse_expression(parse_context_t* ctx) {
    token_t token = {0};
    while (token = tokenizer_peek_next(ctx->tokenizer), token.type != TOKEN_END) {
        // What begins an expression of its own ('3', 'x', '(', 'not', ...) right where one has ended: 'a = 3 4;',
        // 'f(x y)'. Parsed on, the new expression would silently take the place of the one before it.
        switch (token.type) {
        case TOKEN_IDENT:
        case TOKEN_FLOAT:
        case TOKEN_INT:
        case TOKEN_STRING:
        case ':':
        case '(':
        case '{':
        case TOKEN_NOT:
        case TOKEN_OUT:
            if (ctx->node) {
                LOG_ERROR(ctx->ir, concat_tokens(ctx->node->token, token), "Invalid syntax, expressions must be separated by an operator");
                return NULL;
            }
            break;
        default:
            break;
        }

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
            break;
        case '<':
        case '>':
        case TOKEN_GE:
        case TOKEN_LE:
        case TOKEN_EQ:
            ctx->node = parse_comparison(ctx);
            fix_precedence(&ctx->node);
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
        case TOKEN_INT_DIV:
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

static bool evaluate_node_alloc(data_t* out, const ast_node_t* node, eval_context_t* ctx, md_allocator_i* alloc);

// The evaluated arguments of a procedure call
typedef struct proc_args_t {
    ast_node_t* const* nodes;
    size_t  count;
    data_t  data  [MAX_SUPPORTED_PROC_ARGS];
    token_t tokens[MAX_SUPPORTED_PROC_ARGS];
    flags_t flags [MAX_SUPPORTED_PROC_ARGS];
    md_bitset_t empty[MAX_SUPPORTED_PROC_ARGS];   // What an argument which is an empty array of bitfields points at
} proc_args_t;

// Evaluates the arguments of a call to proc, each exactly once. An argument whose length is only known at run
// time is sized by evaluating it (see evaluate_node_alloc), not by a separate query ahead of its evaluation.
// Whatever succeeded is released by proc_args_free, also on failure.
static bool proc_args_evaluate(proc_args_t* pa, const procedure_t* proc, ast_node_t** const args, size_t num_args, eval_context_t* ctx) {
    ASSERT(pa);
    ASSERT(num_args < MAX_SUPPORTED_PROC_ARGS);
    MEMSET(pa, 0, sizeof(proc_args_t));
    pa->nodes = args;
    pa->count = num_args;

    md_script_vis_t* old_vis = ctx->vis;
    // If we are in visualization mode and the procedure is flagged with visualization, we do not want to propagate the visualization to the arguments
    if (ctx->vis && proc->flags & FLAG_VISUALIZE) {
        ctx->vis = 0;
    }

    bool result = true;
    for (size_t i = 0; i < num_args; ++i) {
        pa->tokens[i] = args[i]->token;
        pa->flags[i]  = args[i]->flags;

        if (args[i]->flags & FLAG_CONSTANT) {
            // Evaluated by the static check: a type (or absent, without data) and a length
            if (!INTERNAL_CHECK(ctx, args[i]->token, (args[i]->data.type.base_type != TYPE_UNDEFINED || args[i]->data.size == 0) && args[i]->data.type.dim[0] != -1)) {
                result = false;
                break;
            }
            MEMCPY(&pa->data[i], &args[i]->data, sizeof(data_t));
        } else if (!evaluate_node_alloc(&pa->data[i], args[i], ctx, ctx->temp_alloc)) {
            // Only used as arguments, so the data goes in the temp alloc
            result = false;
            break;
        }
        // An empty array of bitfields ('water' of a system without any) has no data at all, while a procedure which
        // takes a bitfield reads one through the pointer: it reads an empty one. Its length is still 0.
        if (pa->data[i].type.base_type == TYPE_ATOMSET && !pa->data[i].ptr) {
            pa->data[i].ptr = &pa->empty[i];
        }
        // A parameter of a fixed size (a float, a vec3, a range) of which the argument holds less. An empty array is
        // typed as a scalar (all its dimensions are 0): 'sqrt(coord_z(flatten(water)))' is matched with sqrt(float),
        // and has no float to read.
        const type_info_t param = proc->arg_type[i];
        if (param.base_type != TYPE_ATOMSET && param.base_type != TYPE_COORDINATE && !is_variable_length(param) &&
            !is_absent_type(pa->data[i].type) && pa->data[i].size < type_info_total_byte_size(param)) {
            EVAL_LOG_ERROR(ctx, args[i]->token, "The argument has no value: it is empty");
            result = false;
            break;
        }
    }

    ctx->vis = old_vis;
    return result;
}

static void proc_args_free(proc_args_t* pa, eval_context_t* ctx) {
    // @NOTE(Robin): We cannot simply reset the stack since atomsets are lazily allocated, meaning they allocate data on demand.
    // This means that they will most likely allocate data deeply nested within the calling scope of another procedure call with another reset point for its stack
    // Thus, the data for the atomsets will be reset with that reset point and we will overwrite it later.
    // To solve it I think we need another 'specific' allocator for just the atomsets which may outlive the other stack variables and can be reset separately.
    // Conceptually we could use the top end of the stack for this and have two different stack pointers.
    for (int64_t i = (int64_t)pa->count - 1; i >= 0; --i) {
        if (!(pa->nodes[i]->flags & FLAG_CONSTANT)) {
            free_data(&pa->data[i], ctx->temp_alloc);
        }
    }
}

// Calls proc on already evaluated arguments. dst == NULL is a query (or a static validation), which for a
// procedure flagged FLAG_QUERYABLE_LENGTH returns the length of its result given these arguments.
static int proc_invoke(data_t* dst, const procedure_t* proc, proc_args_t* pa, eval_context_t* ctx) {
    flags_t* old_arg_flags  = ctx->arg_flags;
    token_t* old_arg_tokens = ctx->arg_tokens;
    flags_t  old_proc_flags = ctx->proc_flags;

    ctx->arg_tokens = pa->tokens;
    ctx->arg_flags  = pa->flags;
    ctx->proc_flags = proc->flags;

#ifdef MD_SCRIPT_TEST_HOOKS
    if (dst && proc->proc_ptr == test_hook_proc_eval.proc_ptr) {
        test_hook_proc_eval.count += 1;
    }
#endif

    const int result = proc->proc_ptr(dst, pa->data, ctx);

    ctx->arg_flags  = old_arg_flags;
    ctx->arg_tokens = old_arg_tokens;
    ctx->proc_flags = old_proc_flags;
    return result;
}

static int do_proc_call(data_t* dst, const procedure_t* proc,  ast_node_t** const args, size_t num_args, eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(proc);
    ASSERT(num_args < MAX_SUPPORTED_PROC_ARGS);

    // Without a destination and without a visualization, a call can only be to validate the arguments or to ask
    // for the length of the result, which only some procedures do. The others are written for a destination
    // (they dereference it): for them there is nothing to do.
    if (!dst && !ctx->vis && !(proc->flags & (FLAG_STATIC_VALIDATION | FLAG_QUERYABLE_LENGTH))) {
        return 0;
    }

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
        return 0;
    }

    proc_args_t args_data;
    int result = -1;
    if (proc_args_evaluate(&args_data, proc, args, num_args, ctx)) {
        result = proc_invoke(dst, proc, &args_data, ctx);
    }
    proc_args_free(&args_data, ctx);
    return result;
}

static bool evaluate_proc_call(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    if (!INTERNAL_CHECK(ctx, node->token, node->proc && md_array_size(node->children) < MAX_SUPPORTED_PROC_ARGS)) return false;

    token_t old_token = ctx->op_token;
    const ast_node_t* old_node = ctx->op_node;
    ctx->op_token = node->token;
    ctx->op_node  = node;
    int result = do_proc_call(dst, node->proc, node->children, md_array_size(node->children), ctx);
    ctx->op_token = old_token;
    ctx->op_node  = old_node;

    return result >= 0;
}

static bool evaluate_constant_value(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_CONSTANT_VALUE);
    (void)ctx;

    if (dst) {
        // The static check set the data, and dst has room for it
        if (!INTERNAL_CHECK(ctx, node->token, node->data.ptr && node->data.size > 0)) return false;
        return copy_data(dst, &node->data, ctx);
    }

    return true;
}


static bool evaluate_attribute(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ATTRIBUTE);

    if (!dst) {
        return true;
    }
    if (!INTERNAL_CHECK(ctx, node->token, dst->ptr && dst->size >= type_info_total_byte_size(node->data.type))) return false;

    if (!ctx->cur_state || !md_state_has_frame(ctx->cur_state)) {
        MD_LOG_ERROR("attr: the evaluated state carries no frame");
        return false;
    }

    const md_attributes_t* attributes = &ctx->sys->attributes;
    // The frame is a frame of the run being evaluated, which need not be the run the attribute
    // was found in: another member of an ensemble is matched to it by time.
    const md_attribute_t* run_ax = md_attributes_get(attributes, ctx->run_axis_id ? ctx->run_axis_id : node->run_axis_id);
    const size_t frame = (size_t)md_state_frame_nearest(ctx->cur_state);

    // The members one after the other, each found on its own frame axis
    float* out = (float*)dst->ptr;
    size_t cap = dst->size / sizeof(float);
    const size_t num_members = md_array_size(node->attr_members);
    for (size_t m = 0; m < num_members; ++m) {
        const attr_member_t* member = &node->attr_members[m];
        const md_attribute_t* attr    = md_attributes_get(attributes, member->id);
        const md_attribute_t* attr_ax = md_attributes_get(attributes, member->axis_id);
        if (!attr || !attr_ax || !run_ax) {
            MD_LOG_ERROR("attr: '"STR_FMT"' is no longer in the system", STR_ARG(member->path));
            return false;
        }

        size_t row;
        if (!md_attribute_axis_map(&row, run_ax, frame, attr_ax)) {
            MD_LOG_ERROR("attr: '"STR_FMT"' has no value at frame %zu", STR_ARG(member->path), frame);
            return false;
        }

        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)row);
        const size_t count = md_attribute_slice_count(attr, slice);
        if (count > cap) {
            MD_LOG_ERROR("attr: '"STR_FMT"' changed shape since the script was compiled", STR_ARG(member->path));
            return false;
        }
        if (md_attribute_extract_f32(out, count, attr, slice, md_unit_none()) != count) {
            return false;
        }
        out += count;
        cap -= count;
    }
    return true;
}

static bool evaluate_flatten(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_FLATTEN);
    if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 1)) return false;
    const ast_node_t* arg = node->children[0];

    if (dst) {
        if (node->data.type.base_type == TYPE_ATOMSET) {
            data_t data = {0};
            const md_bitset_t* src_arr = 0;
            size_t src_len = 0;

            if (arg->flags & FLAG_CONSTANT && arg->data.ptr) {
                // Assume that it has data. No need to copy here...
                src_arr = arg->data.ptr;
                src_len = type_info_array_len(arg->data.type);
            } else {
                if (!allocate_data(&data, arg->data.type, ctx->temp_alloc)) {
                    EVAL_LOG_ERROR(ctx, node->token, "Failed to allocate data to evaluate flatten operation");
                    return false;
                }
                if (!evaluate_node(&data, arg, ctx)) {
                    return false;
                }
                src_arr = data.ptr;
                src_len = type_info_array_len(data.type);
            }

            md_bitmask_t b = md_bitmask_create(0, atom_domain_end(ctx), ctx->temp_alloc);
            md_bitmask_or(&b, *(md_bitset_t*)dst->ptr);
            for (size_t i = 0; i < src_len; ++i) {
                md_bitmask_or(&b, src_arr[i]);
            }
            *(md_bitset_t*)dst->ptr = md_bitmask_extract(&b, data_bitset_alloc(dst, ctx));
            md_bitmask_free(&b);
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
    if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 1)) return false;
    const ast_node_t* arg = node->children[0];
    int ndim = dim_ndims(arg->data.type.dim);

    if (dst && ndim == 2) {
        if (ndim != dim_ndims(dst->type.dim)) {
            EVAL_LOG_ERROR(ctx, node->token, "Incorrect dimensions in transpose operation");
            return false;
        }
        if (!INTERNAL_CHECK(ctx, node->token, dst->type.base_type == arg->data.type.base_type &&
            dst->type.dim[0] == arg->data.type.dim[1] && dst->type.dim[1] == arg->data.type.dim[0])) return false;

        data_t data = arg->data;
        data.ptr = 0;
        data.size = 0;
        /*
        // @NOTE: For now, we just support static length arrays as the transposition throws a bone into the machinery when it comes to determining length
        // That is, it shifts the dimension of the undetermined length
        if (arg->data.type.dim[0] == -1) {
            if (!finalize_type(&data.type, arg, ctx)) {
                EVAL_LOG_ERROR(ctx, node->token, "Failed to finalize type of subexpression in transpose operation");
                return false;
            }
            dst->type.dim[1] = arg->data.type.dim[0];
        }
        */
        if (!allocate_data(&data, data.type, ctx->temp_alloc)) {
            EVAL_LOG_ERROR(ctx, node->token, "Failed to allocate data to evaluate transpose operation");
            return false;
        }
        if (!evaluate_node(&data, arg, ctx)) {
            EVAL_LOG_ERROR(ctx, node->token, "Failed to evaluate subexpression in transpose operation");
            return false;
        }
        const size_t num_col = (size_t)dst->type.dim[0];
        const size_t num_row = (size_t)dst->type.dim[1];
        const size_t num_bytes = base_type_element_byte_size(dst->type.base_type);
        char* dst_base = dst->ptr;
        const char* src_base = data.ptr;
        for (size_t col = 0; col < num_col; ++col) {
            for (size_t row = 0; row < num_row; ++row) {
                MEMCPY(dst_base + (col * num_row + row) * num_bytes, src_base + (row * num_col + col) * num_bytes, num_bytes);
            }
        }
        return true;
    } else {
        return evaluate_node(dst, arg, ctx);
    }
}

static const identifier_t* find_static_identifier(str_t name, eval_context_t* ctx) {
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

// The value of a dynamic identifier which has already been evaluated for the current frame, or NULL.
// Identifiers are evaluated in the order they are declared and a reference can only follow its declaration,
// so within a frame the value is available by the time anything refers to it. It can only stand in for the
// expression when referenced in the same context it was evaluated in: 'd in residue(1:3)' evaluates the
// expression of d once per residue, which is a different value.
static const data_t* find_evaluated_identifier(str_t name, eval_context_t* ctx) {
    const identifier_t* ident = find_dynamic_identifier(name, ctx);
    if (ident && ident->data && ident->data->ptr && ident->mol_ctx == ctx->mol_ctx) {
        return ident->data;
    }
    return NULL;
}

// Evaluates an identifier declared by destructuring, '{x,y,z} = com(...)': the whole expression is evaluated
// and the element belonging to the identifier is copied out of it.
static bool evaluate_identifier_element(data_t* dst, const identifier_t* decl, const ast_node_t* expr, eval_context_t* ctx) {
    ASSERT(dst && decl && decl->is_element);

    type_info_t type = expr->data.type;
    if (is_variable_length(type) && !finalize_type(&type, expr, ctx)) {
        return false;
    }

    md_temp_scope_t temp = md_temp_begin_in(ctx->temp_alloc);
    data_t whole = {0};
    bool result = allocate_data(&whole, type, md_temp_allocator(temp)) && evaluate_node(&whole, expr, ctx);
    if (result) {
        const size_t stride = type_info_element_byte_stride(type);
        if ((decl->elem_idx + 1) * stride <= whole.size && stride == dst->size) {
            data_t elem = {
                .type = dst->type,
                .ptr  = (char*)whole.ptr + stride * decl->elem_idx,
                .size = stride,
                .weight = data_weight_at(&whole, stride * decl->elem_idx),
            };
            result = copy_data(dst, &elem, ctx);
        } else {
            MD_LOG_ERROR("Script: size mismatch when extracting element %u of '"STR_FMT"'", decl->elem_idx, STR_ARG(decl->name));
            result = false;
        }
    }
    md_temp_end(temp);
    return result;
}

static bool evaluate_identifier_reference(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_IDENTIFIER);
    if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 1 && node->children[0])) return false;

    // Only copy data if the node is constant
    if (dst && node->flags & FLAG_CONSTANT) {
        return copy_data(dst, &node->data, ctx);
    }

    if (dst) {
        // A dynamic identifier is evaluated once per frame, references reuse that value
        const data_t* value = find_evaluated_identifier(node->ident, ctx);
        if (value && value->size == dst->size) {
            if (!copy_data(dst, value, ctx)) return false;
            // Visualization is not memoized, it is produced by evaluating the expression
            return ctx->vis ? evaluate_node(NULL, node->children[0], ctx) : true;
        }

        // Not available (a different context, or no evaluation of the declaration in this pass):
        // an identifier declared by destructuring is one element of what its expression evaluates to
        const identifier_t* decl = find_ir_identifier(ctx->ir, node->ident);
        if (decl && decl->is_element) {
            if (!evaluate_identifier_element(dst, decl, node->children[0], ctx)) {
                return false;
            }
            return ctx->vis ? evaluate_node(NULL, node->children[0], ctx) : true;
        }
    }

    if (dst || ctx->vis) {
        return evaluate_node(dst, node->children[0], ctx);
    }

    return true;
}

static bool evaluate_assignment(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ASSIGNMENT);
    if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 2 && node->children[0] && node->children[1])) return false;
    ast_node_t* lhs = node->children[0];
    ast_node_t* rhs = node->children[1];
    if (!INTERNAL_CHECK(ctx, node->token, lhs->type == AST_IDENTIFIER || lhs->type == AST_ARRAY)) return false;

    // We know from the static check that the lhs (child[0]) is one or more identifiers and the rhs is a valid expression with a known data size
    // the data size could be zero if it is an array of length 0, in that case we don't need to allocate anything.
    // This should also be the conceptual entry point for all evaluations in this declarative language
    if (lhs->type == AST_IDENTIFIER) {
        if (!dst && ctx->vis) {
            return evaluate_node(NULL, rhs, ctx);
        }

        const identifier_t* static_ident = find_static_identifier(lhs->ident, ctx);
        if (static_ident) {
            if (dst && !copy_data(dst, static_ident->data, ctx)) {
                return false;
            }
            if (ctx->vis) {
                return evaluate_node(NULL, rhs, ctx);
            }
            return true;
        }

        const identifier_t* ident = find_dynamic_identifier(lhs->ident, ctx);
        if (!ident && ctx->alloc) {
            // Registered once evaluated, so a reference never picks up a value which is still being computed
            if (!evaluate_node(dst, rhs, ctx)) {
                return false;
            }
            identifier_t id = {
                .name = lhs->ident,
                .node = rhs,
                .data = dst,
                .mol_ctx = ctx->mol_ctx,
            };
            md_array_push(ctx->identifiers, id, ctx->alloc);
            return true;
        }
    } else if (lhs->type == AST_ARRAY) {
        // assignment of LHS which is an array of Identifiers
        // e.g {x,y,z} = coord(1);
        int count = (int)md_array_size(lhs->children);
        if (rhs->flags & FLAG_CONSTANT) {
            if (dst && !copy_data(dst, &rhs->data, ctx)) {
                return false;
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
            if (!INTERNAL_CHECK(ctx, child->token, child->type == AST_IDENTIFIER)) return false;
            identifier_t* ident = find_dynamic_identifier(child->ident, ctx);

            if (!ident && ctx->alloc) {
                identifier_t id = {
                    .name = child->ident,
                    .node = rhs,
                    .data = md_alloc(ctx->alloc, sizeof(data_t)),
                    .mol_ctx = ctx->mol_ctx,
                };
                MEMSET(id.data, 0, sizeof(data_t));
                if (dst && dst->ptr) {
                    *id.data = child->data;
                    id.data->ptr = (char*)dst->ptr + stride * i;
                    id.data->weight = data_weight_at(dst, stride * i);
                    id.data->size = stride;
                }
                md_array_push(ctx->identifiers, id, ctx->alloc);
                ident = md_array_last(ctx->identifiers);
            }
        }
        return true;
    }

    return false;
}

static bool index_in_range(int idx, irange_t range) {
    if (!(range.beg <= idx && idx < range.end)) return false;
    int st = range.step > 0 ? range.step : 1;
    return ((idx - range.beg) % st) == 0;
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
                .size = type_info_total_byte_size(type),
                .alloc = dst->alloc,
                .weight = data_weight_at(dst, byte_offset),
            };
            // Within a context dst has the type of the part, which need not be the type of the node: what matters is
            // that the elements fit
            if (!INTERNAL_CHECK(ctx, args[i]->token, byte_offset + data.size <= dst->size)) goto done;
            md_script_vis_t* vis = ctx->vis;
            if (ranges && !index_in_range((int)i, ranges[0])) {
                ctx->vis = 0;
            }
            const bool ok = evaluate_node(&data, args[i], ctx);
            ctx->vis = vis;
            if (!ok) {
                goto done;
            }
            byte_offset += data.size;
        }
    }
    else {
        // Without a destination there is nothing to evaluate but what is visualized. An element which is not, is
        // skipped: evaluated without either, a procedure call would be taken for a query of its length.
        for (size_t i = 0; i < num_args; ++i) {
            if (!ctx->vis || (ranges && !index_in_range((int)i, ranges[0]))) {
                continue;
            }
            if (!evaluate_node(NULL, args[i], ctx)) {
                goto done;
            }
        }
    }

    result = true;
done:
    ctx->subscript_ranges = ranges;
    return result;
}

static bool evaluate_array_subscript(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ARRAY_SUBSCRIPT);
    if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) > 0 && node->subscript_dim <= MAX_NUM_DIMS)) return false;

    const ast_node_t* arr = node->children[0];
    data_t arr_data = {0};

    md_array(irange_t) idx_ranges = md_array_create(irange_t, node->subscript_dim, ctx->temp_alloc);
    md_array(irange_t) ext_ranges = ctx->subscript_ranges;

    // @NOTE: This is a bit of a hack, we should not need to allocate the subscript ranges here
    // It is just to serve the poor design of the current propagation of array subscripts in the context
    MEMCPY(idx_ranges, node->subscript_ranges, sizeof(irange_t) * node->subscript_dim);

    if (ext_ranges) {
        // The ranges from outside may cover fewer dimensions: those beyond are taken as they are subscripted here
        const size_t num_ext = md_array_size(ext_ranges);
        for (size_t i = 0; i < node->subscript_dim; ++i) {
            if (i >= num_ext || (ext_ranges[i].beg == 0 && ext_ranges[i].end == 0)) {
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

    bool evaluated = false;
    if (dst) {
        evaluated = allocate_data(&arr_data, arr->data.type, ctx->temp_alloc) && evaluate_node(&arr_data, arr, ctx);
    } else {
        evaluated = evaluate_node(0, arr, ctx);
    }
    ctx->subscript_ranges = ext_ranges;
    if (!evaluated) return false;

    if (dst) {
        // Extract data from the underlying type given the subscript ranges. Within a context dst has the type of the
        // part (with a leading dimension of one), which is not the type of the node: the blocks are checked to fit.
        if (!INTERNAL_CHECK(ctx, node->token, arr_data.ptr && dst->ptr)) return false;

    const base_type_t base_type = arr_data.type.base_type;
        const size_t base_size = base_type_element_byte_size(base_type);
    // The payloads of the bitsets copied into dst go where dst lives (see data_bitset_alloc)
    md_allocator_i* dst_alloc = data_bitset_alloc(dst, ctx);
        size_t write_offset = 0;
        const irange_t* src_rng = node->subscript_ranges;
        const int* src_dim = arr_data.type.dim;

        switch (node->subscript_dim) {
            case 1: {
                size_t block_elems = 1;
                for (int d = 1; d < MAX_NUM_DIMS; ++d) {
                    if (src_dim[d] <= 0) break;
                    block_elems *= (size_t)src_dim[d];
                }
                const size_t block_bytes = block_elems * base_size;
                for (int i = src_rng[0].beg; i < src_rng[0].end; i += MAX(src_rng[0].step,1)) {
                    const size_t read_offset = (size_t)i * block_bytes;
                    if (!INTERNAL_CHECK(ctx, node->token, write_offset + block_bytes <= dst->size && read_offset + block_bytes <= arr_data.size)) return false;
                    if (!copy_values_deep((char*)dst->ptr + write_offset, (const char*)arr_data.ptr + read_offset, base_type, block_elems, dst_alloc)) {
                        copy_data_block(dst, write_offset, &arr_data, read_offset, block_bytes);
                    }
                    write_offset += block_bytes;
                }
            } break;
            case 2:
                {
                    size_t block_elems = 1;
                    for (int d = 2; d < MAX_NUM_DIMS; ++d) {
                        if (src_dim[d] <= 0) break;
                        block_elems *= (size_t)src_dim[d];
                    }
                    const size_t block_bytes = block_elems * base_size;
                    for (int i = src_rng[0].beg; i < src_rng[0].end; i += MAX(src_rng[0].step,1)) {
                        for (int j = src_rng[1].beg; j < src_rng[1].end; j += MAX(src_rng[1].step,1)) {
                            const size_t lin = (size_t)i * (size_t)src_dim[1] + (size_t)j;
                            const size_t read_offset = lin * block_bytes;
                            if (!INTERNAL_CHECK(ctx, node->token, write_offset + block_bytes <= dst->size && read_offset + block_bytes <= arr_data.size)) return false;
                            if (!copy_values_deep((char*)dst->ptr + write_offset, (const char*)arr_data.ptr + read_offset, base_type, block_elems, dst_alloc)) {
                                copy_data_block(dst, write_offset, &arr_data, read_offset, block_bytes);
                            }
                            write_offset += block_bytes;
                        }
                    }
                }
                break;
            case 3:
                // @NOTE: TO BE TESTED
                {
                    size_t block_elems = 1;
                    for (int d = 3; d < MAX_NUM_DIMS; ++d) {
                        if (src_dim[d] <= 0) break;
                        block_elems *= (size_t)src_dim[d];
                    }
                    const size_t block_bytes = block_elems * base_size;
                    for (int k = src_rng[0].beg; k < src_rng[0].end; k += MAX(src_rng[0].step,1)) {
                        for (int j = src_rng[1].beg; j < src_rng[1].end; j += MAX(src_rng[1].step,1)) {
                            for (int i = src_rng[2].beg; i < src_rng[2].end; i += MAX(src_rng[2].step,1)) {
                                const size_t lin = ((size_t)k * (size_t)src_dim[1] + (size_t)j) * (size_t)src_dim[2] + (size_t)i;
                                const size_t read_offset = lin * block_bytes;
                                if (!INTERNAL_CHECK(ctx, node->token, write_offset + block_bytes <= dst->size && read_offset + block_bytes <= arr_data.size)) return false;
                                if (!copy_values_deep((char*)dst->ptr + write_offset, (const char*)arr_data.ptr + read_offset, base_type, block_elems, dst_alloc)) {
                                    copy_data_block(dst, write_offset, &arr_data, read_offset, block_bytes);
                                }
                                write_offset += block_bytes;
                            }
                        }
                    }
                }
                break;
            case 4:
                // @NOTE: TO BE TESTED
                {
                    size_t block_elems = 1;
                    for (int d = 4; d < MAX_NUM_DIMS; ++d) {
                        if (src_dim[d] <= 0) break;
                        block_elems *= (size_t)src_dim[d];
                    }
                    const size_t block_bytes = block_elems * base_size;
                    for (int l = src_rng[0].beg; l < src_rng[0].end; l += MAX(src_rng[0].step,1)) {
                        for (int k = src_rng[1].beg; k < src_rng[1].end; k += MAX(src_rng[1].step,1)) {
                            for (int j = src_rng[2].beg; j < src_rng[2].end; j += MAX(src_rng[2].step,1)) {
                                for (int i = src_rng[3].beg; i < src_rng[3].end; i += MAX(src_rng[3].step,1)) {
                                    const size_t lin = (((size_t)l * (size_t)src_dim[1] + (size_t)k) * (size_t)src_dim[2] + (size_t)j) * (size_t)src_dim[3] + (size_t)i;
                                    const size_t read_offset = lin * block_bytes;
                                    if (!INTERNAL_CHECK(ctx, node->token, write_offset + block_bytes <= dst->size && read_offset + block_bytes <= arr_data.size)) return false;
                                    if (!copy_values_deep((char*)dst->ptr + write_offset, (const char*)arr_data.ptr + read_offset, base_type, block_elems, dst_alloc)) {
                                        copy_data_block(dst, write_offset, &arr_data, read_offset, block_bytes);
                                    }
                                    write_offset += block_bytes;
                                }
                            }
                        }
                    }
                }
                break;
            default:
                INTERNAL_ERROR(ctx, node->token, "an array subscript of more dimensions than there are");
                return false;
        }
    }

    return true;
}

static bool evaluate_out(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(node && node->type == AST_OUT);
    if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 1)) return false;

    const ast_node_t* expr = node->children[0];
    const md_bitset_t* old_ctx = ctx->mol_ctx;
    ctx->mol_ctx = 0;
    bool result = evaluate_node(dst, expr, ctx);
    ctx->mol_ctx = old_ctx;

    return result;
}

static bool evaluate_context(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(ctx);
    ASSERT(node && node->type == AST_CONTEXT);

    if (dst && node->flags & FLAG_CONSTANT) {
        return copy_data(dst, &node->data, ctx);
    }

    if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 2)) return false;

    const ast_node_t* expr_node = node->children[0];
    ast_node_t* ctx_node  = node->children[1];

    // No need to reevaluate CTX_NODE here since it has already been done during static type checking and its results should be stored in the data ptr

    if (!INTERNAL_CHECK(ctx, ctx_node->token, ctx_node->data.type.base_type == TYPE_ATOMSET)) return false;
    data_t ctx_data = ctx_node->data;
    if (!ctx_data.ptr) {
        // A compilation stores what it computes in the AST. An evaluation does not write to the IR, so it
        // computes the contexts for this call only.
        md_allocator_i* ctx_alloc = ctx->compile_ir ? ctx->compile_ir->arena : ctx->temp_alloc;
        if (!allocate_data(&ctx_data, ctx_node->data.type, ctx_alloc) || !evaluate_node(&ctx_data, ctx_node, ctx)) {
            return false;
        }
        if (ctx->compile_ir) {
            ctx_node->data = ctx_data;
        }
    }

    const int num_ctx = type_info_array_len(ctx_data.type);
    if (num_ctx < 0) {
        EVAL_LOG_ERROR(ctx, node->token, "Invalid number of contexts (%i) in context expression", num_ctx);
		return false;
    }
    const md_bitset_t* ctx_bf = (const md_bitset_t*)ctx_data.ptr;

    // A type for the left hand side within each context, from the static check
    const type_info_t* lhs_types = node->lhs_context_types;
    if (!INTERNAL_CHECK(ctx, node->token, lhs_types && (int)md_array_size(lhs_types) == num_ctx)) return false;

    data_t data = {0};
    data.unit[0] = expr_node->data.unit[0];
    data.unit[1] = expr_node->data.unit[1];
    data.value_range = expr_node->data.value_range;
    data.period = expr_node->data.period;
   
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
            data.weight = data_weight_at(dst, elem_size * dst_idx);
            data.size = type_info_total_byte_size(lhs_types[i]);
            data.alloc = dst->alloc;

            int64_t offset = type_info_array_len(lhs_types[i]);
            dst_idx += offset;
        }

        if (ctx->subscript_ranges) {
            if (!index_in_range((int)i, ctx->subscript_ranges[0])) {
                // Clear the visualization context if this index is not present in subscript ranges to prevent visualization
                sub_ctx.vis = 0;
            }
        }

        // Each context is a structure of its own: built while the expression is visualized within it, then
        // stored in the visualization
        md_bitmask_t* prev_vis_structure = 0;
        md_bitmask_t structure = {0};
        size_t structure_idx = SIZE_MAX;
        if (sub_ctx.vis) {
            prev_vis_structure = sub_ctx.vis_structure;
            if (sub_ctx.vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
                structure_idx = push_structure(sub_ctx.vis);
                structure = md_bitmask_create(0, atom_domain_end(ctx), ctx->temp_alloc);
                sub_ctx.vis_structure = &structure;
            }
        }

        if (!evaluate_node(sub_dst, expr_node, &sub_ctx)) {
            result = false;
        }

        if (sub_ctx.vis) {
            // Reset atom mask
            if (sub_ctx.vis_flags & MD_SCRIPT_VISUALIZE_ATOMS) {
                sub_ctx.vis->structures[structure_idx] = md_bitmask_extract(&structure, sub_ctx.vis->alloc);
                md_bitmask_free(&structure);
                sub_ctx.vis_structure = prev_vis_structure;
            }
        }

        sub_ctx.vis = ctx->vis;

        if (!result) break;
    }

    return result;
}

static bool evaluate_node(data_t* dst, const ast_node_t* node, eval_context_t* ctx) {
    if (!INTERNAL_CHECK(ctx, ctx->op_token, node)) return false;

    if (dst && node->flags & FLAG_CONSTANT) {
        return copy_data(dst, &node->data, ctx);
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
        case AST_ATTRIBUTE:
            return evaluate_attribute(dst, node, ctx);
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
            INTERNAL_ERROR(ctx, node->token, "a node which the static check should have resolved into a procedure call");
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

        if (node->flags & FLAG_CONSTANT && !(ctx->eval_flags & EVAL_FLAG_NO_STATIC_EVAL)) {
            if (to.dim[0] == -1) {
                if (res.procedure->flags & FLAG_DEDUCE_LENGTH_FROM_ARG) {
                    to.dim[0] = from.dim[0];
                } else {
                    if (!INTERNAL_CHECK(ctx, node->token, res.procedure->flags & FLAG_QUERYABLE_LENGTH)) return false;

                    token_t  old_op_token = ctx->op_token;
                    ctx->op_token = node->token;
                    // Perform the call to get length
                    int query_result = do_proc_call(NULL, res.procedure, &node, 1, ctx);
                    ctx->op_token = old_op_token;

                    if (query_result >= 0) { // Zero length is valid
                        to.dim[0] = query_result;
                    }
                    else {
                        EVAL_LOG_ERROR(ctx, node->token, "Unexpected return value (%i) when querying procedure for array length.", query_result);
                        return false;
                    }
                }
            }

            if (!INTERNAL_CHECK(ctx, node->token, to.dim[0] > -1)) return false;

            data_t new_data = {0};
            if (node->type == AST_CONSTANT_VALUE && to.base_type != TYPE_ATOMSET) {
                new_data = node->data;
                if (!INTERNAL_CHECK(ctx, node->token, is_scalar(to))) return false;
                new_data.type = to;
                new_data.size = type_info_total_byte_size(to);
                new_data.ptr  = &node->value;
            } else if (!allocate_data(&new_data, to, ctx->alloc)) {
                EVAL_LOG_ERROR(ctx, node->token, "Failed to allocate data for the conversion of static data");
                return false;
            }

            // Perform the data conversion
            if (do_proc_call(&new_data, res.procedure, &node, 1, ctx) < 0) {
                EVAL_LOG_ERROR(ctx, node->token, "Failed to perform data conversion of static data");
                return false;
            }
            // Finalize the node with the new converted data
            node->data = new_data;
            return true;
        } else {
            // We need to convert this node into a cast node and add the original node data as a child
            ast_node_t* node_copy = create_node(ctx->compile_ir, node->type, node->token);
            copy_node_contents(node_copy, node);
            node->data = (data_t){0};
            node->data.unit[0] = node_copy->data.unit[0];
            node->data.unit[1] = node_copy->data.unit[1];
            node->data.period = node_copy->data.period;
            node->type = AST_PROC_CALL;
            node->children = 0; // node_copy have taken over the children, we need to zero this to trigger a proper allocation in next step
            md_array_push(node->children, node_copy, ctx->compile_ir->arena);
            node->proc = res.procedure;
            node->proc_flags = res.flags;
            node->flags &= ~FLAG_CONSTANT;

            return finalize_proc_call(node, ctx);
        }
    }

    char from_buf[64], to_buf[64];
    print_type_info(from_buf, sizeof(from_buf), from);
    print_type_info(to_buf, sizeof(to_buf), to);
    EVAL_LOG_ERROR(ctx, node->token, "No conversion from '%s' to '%s'", from_buf, to_buf);
    return false;
}

static void deduce_type_dim_from_types(type_info_t* type, const type_info_t* arg_types, size_t num_args) {
    for (size_t i = 0; i < num_args; ++i) {
    // If the argument has leading ones in its type, we want to replicate that
        if (dim_ndims(arg_types[i].dim) > 1 && arg_types[i].dim[i] == 1) {
            dim_shift_right(type->dim);
            type->dim[1] = arg_types[i].dim[1];
        } else {
            type->dim[0] = MAX(type->dim[0], arg_types[i].dim[0]);
        }
    }
}

static bool deduce_type_dim_from_args(type_info_t* type, ast_node_t** args, size_t num_args, token_t token, eval_context_t* ctx) {
    ASSERT(ctx);
    (void)token;
    ASSERT(num_args <= MAX_SUPPORTED_PROC_ARGS);

    type_info_t arg_types[MAX_SUPPORTED_PROC_ARGS];
    for (size_t i = 0; i < num_args; ++i) {
        arg_types[i] = args[i]->data.type;
    }
    deduce_type_dim_from_types(type, arg_types, num_args);

    return true;
}

static bool finalize_type_proc(type_info_t* type, const ast_node_t* node, eval_context_t* ctx) {
    ASSERT(type);
    ASSERT(node);
    ASSERT(ctx);

    *type = node->data.type;
    if (!INTERNAL_CHECK(ctx, node->token, node->proc)) return false;

    ast_node_t** args = node->children;
    const size_t num_args = md_array_size(node->children);

    if (node->proc->flags & FLAG_DEDUCE_LENGTH_FROM_ARG) {
        if (!args || num_args == 0) {
            EVAL_LOG_ERROR(ctx, node->token, "Unexpected number of arguments (0), but procedure is marked to extract length from first argument");
            return false;
        }
        // We can deduce the dimensions from the input argument types
        if (!deduce_type_dim_from_args(type, args, num_args, node->token, ctx)) {
            return false;
        }
    } else {
        if (!INTERNAL_CHECK(ctx, node->token, node->proc->flags & FLAG_QUERYABLE_LENGTH)) return false;

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
            EVAL_LOG_ERROR(ctx, node->token, "Unexpected return value (%i) when querying procedure for array length.", query_result);
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
        if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 2)) return false;
        return finalize_type(type, node->children[1], ctx);
    case AST_ARRAY:
        return finalize_type_array(type, node, ctx);
    case AST_IDENTIFIER: {
        if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 1)) return false;
        // The length of an identifier already evaluated this frame is known, asking its expression would evaluate it again
        const data_t* value = find_evaluated_identifier(node->ident, ctx);
        if (value) {
            *type = value->type;
            return true;
        }
        return finalize_type(type, node->children[0], ctx);
    }
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

    case AST_OUT: {
        // The length of what is evaluated outside of the context
        if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 1)) return false;
        const md_bitset_t* old_ctx = ctx->mol_ctx;
        ctx->mol_ctx = NULL;
        const bool result = finalize_type(type, node->children[0], ctx);
        ctx->mol_ctx = old_ctx;
        return result;
    }
    case AST_ARRAY_SUBSCRIPT: // Array subscripts should not be variable length, since the input should be statically evaluated so no out of bounds accesses occur...
    case AST_CONSTANT_VALUE:
    case AST_CONTEXT:   // We don't allow dynamic length expressions for Contexts for now... Perhaps in the future...
    default:
        INTERNAL_ERROR(ctx, node->token, "the length of this expression is not known until it is evaluated, and it cannot be determined");
    }
    return false;
}

// ### VARIABLE LENGTH EVALUATION ###
// A value whose length is only known at run time used to be sized by finalize_type ahead of its evaluation,
// which for a procedure call meant evaluating its arguments once to ask for the length and once more to fill
// in the result. Nested variable length calls doubled the work per level. Here the arguments are evaluated
// once and the length follows from them, so every node is evaluated exactly once.

// Evaluates a procedure call whose result has a variable length into out, allocated from alloc.
static bool evaluate_proc_call_alloc(data_t* out, const ast_node_t* node, eval_context_t* ctx, md_allocator_i* alloc) {
    ASSERT(node && node->type == AST_PROC_CALL && node->proc);
    const procedure_t* proc = node->proc;

    token_t old_token = ctx->op_token;
    const ast_node_t* old_node = ctx->op_node;
    ctx->op_token = node->token;
    ctx->op_node  = node;

    proc_args_t args;
    bool result = proc_args_evaluate(&args, proc, node->children, md_array_size(node->children), ctx);

    type_info_t type = node->data.type;
    if (result) {
        if (proc->flags & FLAG_ARGS_EQUAL_LENGTH) {
            // Only checked statically for lengths known at compile time
            for (size_t i = 1; i < args.count; ++i) {
                if (type_info_array_len(args.data[i].type) != type_info_array_len(args.data[0].type)) {
                    MD_LOG_ERROR("Script: '"STR_FMT"': arguments of different length (%i and %i)", STR_ARG(node->token.str),
                        type_info_array_len(args.data[0].type), type_info_array_len(args.data[i].type));
                    result = false;
                }
            }
        }
        if (!result) {
            // Reported above
        } else if (proc->flags & FLAG_DEDUCE_LENGTH_FROM_ARG) {
            type_info_t arg_types[MAX_SUPPORTED_PROC_ARGS];
            for (size_t i = 0; i < args.count; ++i) {
                arg_types[i] = args.data[i].type;
            }
            deduce_type_dim_from_types(&type, arg_types, args.count);
        } else if (proc->flags & FLAG_QUERYABLE_LENGTH) {
            // Ask for the length given the evaluated arguments, then fill in the result from the same arguments
            md_script_vis_t* old_vis = ctx->vis;
            ctx->vis = NULL;
            const int length = proc_invoke(NULL, proc, &args, ctx);
            ctx->vis = old_vis;
            if (length < 0) {
                EVAL_LOG_ERROR(ctx, node->token, "Unexpected return value (%i) when querying procedure for array length.", length);
                result = false;
            }
            type.dim[0] = length;
        }
        if (result && is_variable_length(type)) {
            MD_LOG_ERROR("Script: failed to determine the length of '"STR_FMT"'", STR_ARG(node->token.str));
            result = false;
        }
    }

    result = result && allocate_data(out, type, alloc) && proc_invoke(out, proc, &args, ctx) >= 0;

    proc_args_free(&args, ctx);
    ctx->op_token = old_token;
    ctx->op_node  = old_node;
    return result;
}

// Evaluates node into out, allocating its storage from alloc. For a node of static length this is allocation
// followed by evaluate_node. A variable length is determined by the evaluation itself where that is possible
// (procedure calls, identifiers and their declarations), and by finalize_type otherwise.
static bool evaluate_node_alloc(data_t* out, const ast_node_t* node, eval_context_t* ctx, md_allocator_i* alloc) {
    ASSERT(out && node && ctx && alloc);

    type_info_t type = node->data.type;
    if (is_variable_length(type)) {
        switch (node->type) {
        case AST_PROC_CALL:
            if (node->proc) {
                return evaluate_proc_call_alloc(out, node, ctx, alloc);
            }
            break;
        case AST_IDENTIFIER:
            if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 1)) return false;
            if (!find_evaluated_identifier(node->ident, ctx)) {
                // Not evaluated in this context, the reference is its expression
                return evaluate_node_alloc(out, node->children[0], ctx, alloc);
            }
            break;  // Sized by its value in finalize_type, then copied by evaluate_node
        case AST_ASSIGNMENT: {
            // The declaration of a dynamic identifier: evaluate its expression and register the value
            const ast_node_t* lhs = node->children[0];
            ast_node_t* rhs = node->children[1];
            if (lhs->type == AST_IDENTIFIER && !ctx->vis && ctx->alloc &&
                !find_static_identifier(lhs->ident, ctx) && !find_dynamic_identifier(lhs->ident, ctx)) {
                if (!evaluate_node_alloc(out, rhs, ctx, alloc)) {
                    return false;
                }
                identifier_t id = {
                    .name = lhs->ident,
                    .node = rhs,
                    .data = out,
                    .mol_ctx = ctx->mol_ctx,
                };
                md_array_push(ctx->identifiers, id, ctx->alloc);
                return true;
            }
            break;
        }
        default:
            break;
        }

        if (!finalize_type(&type, node, ctx)) {
            MD_LOG_ERROR("Script: failed to determine the length of '"STR_FMT"'", STR_ARG(node->token.str));
            return false;
        }
    }

    return allocate_data(out, type, alloc) && evaluate_node(out, node, ctx);
}

static bool finalize_proc_call(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(ctx);
    if (!INTERNAL_CHECK(ctx, node->token, node->proc && node->type == AST_PROC_CALL)) return false;

    node->data.type = node->proc->return_type;

    const size_t num_args = md_array_size(node->children);
    ast_node_t** args = node->children;

    if (num_args > 0) {
        if (node->proc_flags & FLAG_SYMMETRIC_ARGS) {
            // If the symmetric flag is set, we swap the first two arguments
            if (!INTERNAL_CHECK(ctx, node->token, num_args >= 2)) return false;
            ast_node_t* tmp_child = node->children[0];
            node->children[0] = node->children[1];
            node->children[1] = tmp_child;
            node->proc_flags &= ~FLAG_SYMMETRIC_ARGS; // Clear flags so we don't swap again
        }

        // Make sure we cast the arguments into the expected types of the procedure.
        for (size_t i = 0; i < num_args; ++i) {
            if (!is_absent_type(args[i]->data.type) && !is_type_directly_compatible(args[i]->data.type, node->proc->arg_type[i])) {
                // Types are not directly compatible, but should be implicitly convertible (otherwise the match should never have been found)
                if (!INTERNAL_CHECK(ctx, args[i]->token, is_type_implicitly_convertible(args[i]->data.type, node->proc->arg_type[i]))) return false;
                if (!convert_node(args[i], node->proc->arg_type[i], ctx)) {
                    return false;
                }
                node->flags    |= args[i]->flags & FLAG_AST_PROPAGATION_MASK;
                ctx->compile_ir->flags |= args[i]->flags & FLAG_IR_PROPAGATION_MASK;
            }
        }

        if (node->proc->flags & FLAG_ARGS_EQUAL_LENGTH) {
            // This flag should only be set for procedures with multiple arguments
            if (!INTERNAL_CHECK(ctx, node->token, num_args > 1)) return false;
            // Make sure that the input arrays are of equal length. Otherwise create an error.
            const int expected_len = type_info_array_len(args[0]->data.type);
            for (size_t i = 1; i < num_args; ++i) {
                int len = type_info_array_len(args[i]->data.type);
                if (len != expected_len) {
                    EVAL_LOG_ERROR(ctx, node->token, "Expected array-length of arguments to match. arg 0 has length %i, arg %i has length %i.", expected_len, (int)i, len);
                    return false;
                }
            }
        }
    }

    // Propagate procedure flags
    node->flags    |= node->proc->flags & FLAG_AST_PROPAGATION_MASK;
    ctx->compile_ir->flags |= node->proc->flags & FLAG_IR_PROPAGATION_MASK;

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
                EVAL_LOG_ERROR(ctx, node->token, "Failed to determine length of procedure return type!");
                return false;
            }
        } else if (node->proc->flags & FLAG_DEDUCE_LENGTH_FROM_ARG) {
            // We can deduce the length of the array by using the type of the first argument
            if (!INTERNAL_CHECK(ctx, node->token, num_args > 0 && node->proc->arg_type[0].base_type == args[0]->data.type.base_type)) return false;

            // An argument whose length is only known at run time leaves the length of the result to run time as
            // well. Deducing from it here would take the -1 as a length and let any other argument decide it:
            // 'coord_x(within(3, x)) * 2' became a single value.
            for (size_t i = 0; i < num_args; ++i) {
                if (is_variable_length(args[i]->data.type)) {
                    node->flags |= FLAG_DYNAMIC_LENGTH;
                    return true;
                }
            }

            return deduce_type_dim_from_args(&node->data.type, args, num_args, node->token, ctx);
        } else {
            EVAL_LOG_ERROR(ctx, node->token, "Procedure returns variable length, but its length cannot be determined.");
            return false;
        }
    }

    return true;
}

static bool static_check_children(ast_node_t* node, eval_context_t* ctx) {
    const size_t num_children = md_array_size(node->children);
    for (size_t i = 0; i < num_children; ++i) {
        if (!INTERNAL_CHECK(ctx, node->token, node->children[i] && node->children[i] != node)) return false;
        if (static_check_node(node->children[i], ctx)) {
            node->flags |= node->children[i]->flags & FLAG_AST_PROPAGATION_MASK;
        } else {
            return false;
        }
    }
    return true;
}

// The value of a number known since the compilation: a literal, a constant, or what was evaluated of them ('180 / PI')
static bool constant_scalar(double* out, const ast_node_t* node) {
    if (!(node->flags & FLAG_CONSTANT) || !node->data.ptr || !is_scalar(node->data.type)) return false;
    switch (node->data.type.base_type) {
    case TYPE_FLOAT: *out = *(const float*)node->data.ptr; return true;
    case TYPE_INT:   *out = *(const int*)node->data.ptr;   return true;
    default: return false;
    }
}

// The period of what an arithmetic operator gives, from those of its operands (see PERIODIC VALUES). x is periodic and
// c a constant:
//   x + c, c + x, x - c, c - x     periodic, the period and its seam moved along by c (and turned around for c - x)
//   x * c, c * x, x / c            periodic, the period and its seam scaled by c (1 / c): dihedral(...) * 180 / PI
//   -x                             periodic, the period and its seam turned around
//   x + y, x - y                   both periodic in periods of the same length: periodic (period_combine)
//   x + v, x - v, v + x, v - x     v not periodic, of the unit of x: the period of x (turned around for v - x), which the
//                                  values are moved out of by v. They are still the same a period apart.
// Nothing else is: x * y, c / x and x // c, and whatever does not give a number.
static period_t period_of_operator(ast_type_t op, ast_node_t* const* arg, size_t num_args) {
    if (num_args == 1) {
        return op == AST_UNARY_NEG ? period_scale(arg[0]->data.period, -1.0) : (period_t){0};
    }
    if (num_args != 2) return (period_t){0};

    const period_t pa = arg[0]->data.period;
    const period_t pb = arg[1]->data.period;
    const bool a = period_valid(pa);
    const bool b = period_valid(pb);
    if (!a && !b) return (period_t){0};

    double c = 0;
    switch (op) {
    case AST_ADD:
    case AST_SUB: {
        const double sign = op == AST_ADD ? 1.0 : -1.0;
        if (a && b) {
            return period_same_length(pa, pb) ? period_combine(pa, pb, sign) : (period_t){0};
        }
        if (a) {
            if (constant_scalar(&c, arg[1]))                                   return period_shift(pa, sign * c);
            if (md_unit_equal(arg[0]->data.unit[1], arg[1]->data.unit[1]))     return pa;
            return (period_t){0};
        }
        // The periodic one is the second
        const period_t pr = op == AST_ADD ? pb : period_scale(pb, -1.0);
        if (constant_scalar(&c, arg[0]))                                       return period_shift(pr, c);
        if (md_unit_equal(arg[0]->data.unit[1], arg[1]->data.unit[1]))         return pr;
        return (period_t){0};
    }
    case AST_MUL:
        if (a && !b && constant_scalar(&c, arg[1])) return period_scale(pa, c);
        if (b && !a && constant_scalar(&c, arg[0])) return period_scale(pb, c);
        return (period_t){0};
    case AST_DIV:
        if (a && !b && constant_scalar(&c, arg[1]) && c != 0.0) return period_scale(pa, 1.0 / c);
        return (period_t){0};
    default:
        return (period_t){0};
    }
}

// The period of what a procedure gives which does not state one itself (see static_backchannel_t)
static period_t period_of_proc_call(const ast_node_t* node) {
    ASSERT(node && node->proc);
    const size_t num_args = md_array_size(node->children);
    int (*fn)(data_t*, data_t[], eval_context_t*) = node->proc->proc_ptr;

    // An angle in [-180, 180[, cut where atan2 cuts it: nothing is known of where its values are
    if (fn == _atan2f) {
        return period_make(-PI, PI, -PI);
    }
    if (num_args == 1) {
        const period_t p = node->children[0]->data.period;
        // The least and the greatest of angles are angles of the same period
        if (fn == _min_farr || fn == _max_farr) {
            return p;
        }
        // The sum of n angles is an angle, about n times the middle of their period
        if (fn == _sum_farr && period_valid(p)) {
            const int n = type_info_array_len(node->children[0]->data.type);
            return n > 0 ? period_shift(p, (n - 1) * 0.5 * ((double)p.beg + p.end)) : (period_t){0};
        }
    }
    return (period_t){0};
}

static bool static_check_operator(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(is_operator(node->type));
    ASSERT(ctx);
    // One operand or two: arg_type below has room for no more
    const size_t num_operands = md_array_size(node->children);
    if (!INTERNAL_CHECK(ctx, node->token, num_operands == 1 || num_operands == 2)) return false;

    // @NOTE: if we have a bitwise operator in the tree, we flag this
    // so that we can shortcut some ambigious operations such as 'protein' 'residue' etc.
    // That they can be concatenated into a single atomset instead of an array of atomsets.
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
                case AST_ADD:
                    node->data.unit[0] = md_unit_add(arg[0]->data.unit[0], arg[1]->data.unit[0]);
                    node->data.unit[1] = md_unit_add(arg[0]->data.unit[1], arg[1]->data.unit[1]);
                    break;
                case AST_SUB:
                    node->data.unit[0] = md_unit_sub(arg[0]->data.unit[0], arg[1]->data.unit[0]);
                    node->data.unit[1] = md_unit_sub(arg[0]->data.unit[1], arg[1]->data.unit[1]);
                    break;
                case AST_MUL:
                    node->data.unit[0] = md_unit_mul(arg[0]->data.unit[0], arg[1]->data.unit[0]);
                    node->data.unit[1] = md_unit_mul(arg[0]->data.unit[1], arg[1]->data.unit[1]);
                    break;
                case AST_DIV:
                case AST_INT_DIV:
                    node->data.unit[0] = md_unit_div(arg[0]->data.unit[0], arg[1]->data.unit[0]);
                    node->data.unit[1] = md_unit_div(arg[0]->data.unit[1], arg[1]->data.unit[1]);
                    break;
                default: break;
                }
            } else if (node->type == AST_UNARY_NEG) {
                // -x is in the unit of x
                node->data.unit[0] = arg[0]->data.unit[0];
                node->data.unit[1] = arg[0]->data.unit[1];
            }

            // From the operands as they are written, before finalize_proc_call swaps and converts them
            const period_t period = period_of_operator(node->type, arg, num_args);

            node->type = AST_PROC_CALL;
            node->proc = res.procedure;
            node->proc_flags = res.flags;
            result = finalize_proc_call(node, ctx);
            if (result) {
                node->data.period = period;
            }
        } else {
            char arg_type_str[2][64] = {"empty", "empty"};
            bool bonds = false;
            for (size_t i = 0; i < num_args; ++i) {
                if (arg[i]) print_type_info(arg_type_str[i], ARRAY_SIZE(arg_type_str[i]), arg[i]->data.type);
                bonds |= arg[i] && arg[i]->data.type.base_type == TYPE_BONDSET;
            }
            EVAL_LOG_ERROR(ctx, node->token,
                "Could not find supporting operator '%s' with left hand side type (%s) and right hand side type (%s)%s",
                get_token_type_str(node->token.type), arg_type_str[0], arg_type_str[1], bonds ? BONDS_DO_NOT_MIX_HINT : "");
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
            len += md_snprintf(buf + len, cap - len, ", ");
        }
    }
    return len;
}

// What attr() can read: a temporal attribute holding numbers, part of a run ('run/<name>/...') whose
// frames it follows, with no more dimensions per frame than a script value can hold. The static check
// reports why an attribute is not, and completion lists the ones that are.
typedef enum attr_readable_t {
    ATTR_READABLE = 0,
    ATTR_NOT_TEMPORAL,
    ATTR_NOT_NUMERIC,
    ATTR_NOT_IN_RUN,
    ATTR_RUN_WITHOUT_AXIS,
    ATTR_WITHOUT_AXIS,
    ATTR_TOO_MANY_DIMS,
} attr_readable_t;

typedef struct attr_read_info_t {
    str_t                 run;          // 'run/<name>', a view into the path of the attribute
    const md_attribute_t* run_axis;     // The frames of the run
    const md_attribute_t* attr_axis;    // The frames of the attribute
} attr_read_info_t;

static attr_readable_t attr_readable(const md_attributes_t* attributes, const md_attribute_t* attr, attr_read_info_t* info) {
    ASSERT(attributes && attr && info);
    MEMSET(info, 0, sizeof(*info));

    if (!(attr->flags & MD_ATTRIBUTE_FLAG_TEMPORAL)) {
        return ATTR_NOT_TEMPORAL;
    }
    if (attr->format.type == MD_ATTRIBUTE_TYPE_STR || attr->format.type == MD_ATTRIBUTE_TYPE_NONE) {
        return ATTR_NOT_NUMERIC;
    }

    // The run is the first two segments, and its frames are the ones the script steps through.
    const str_t full = attr->path;
    size_t name_len = 0;
    if (!str_begins_with(full, STR_LIT("run/")) || !str_find_char(&name_len, str_substr(full, 4, SIZE_MAX), '/')) {
        return ATTR_NOT_IN_RUN;
    }
    info->run = str_substr(full, 0, 4 + name_len);

    info->run_axis = md_attributes_find_in(attributes, info->run, STR_LIT("time"));
    if (!info->run_axis || md_attributes_axis(attributes, info->run_axis) != info->run_axis) {
        return ATTR_RUN_WITHOUT_AXIS;
    }
    info->attr_axis = md_attributes_axis(attributes, attr);
    if (!info->attr_axis) {
        return ATTR_WITHOUT_AXIS;
    }

    const uint32_t num_dims = (attr->format.rank - 1) + (attr->format.components > 1 ? 1 : 0);
    if (num_dims > MAX_NUM_DIMS) {
        return ATTR_TOO_MANY_DIMS;
    }
    return ATTR_READABLE;
}

// The path attr() reads attr by in its shortest form: relative to its run where that names it alone, in full where it
// does not. The path as given is looked for first and then in every run (static_check_attribute), so the relative one
// must name nothing else, neither outside the runs nor in another one.
static str_t attr_short_path(const md_attributes_t* attributes, const md_attribute_t* attr, str_t run) {
    const str_t relative = str_substr(attr->path, run.len + 1, SIZE_MAX);
    if (md_attributes_find(attributes, relative)) {
        return attr->path;
    }
    for (md_attribute_iter_t it = md_attributes_iter_children(attributes, STR_LIT("run")); md_attributes_next(&it);) {
        const md_attribute_t* other = md_attributes_find_in(attributes, it.child_path, relative);
        if (other && other != attr) {
            return attr->path;
        }
    }
    return relative;
}

str_t md_script_attr_path(const md_attributes_t* attributes, const md_attribute_t* attr) {
    if (!attributes || !attr) {
        return (str_t){0};
    }
    attr_read_info_t info;
    if (attr_readable(attributes, attr, &info) != ATTR_READABLE) {
        return (str_t){0};
    }
    return attr_short_path(attributes, attr, info.run);
}

// What one value of attr at a frame is: the index axes left after the frame axis, then the components of a value, from
// dimension d on. A 3x3 tensor is float[3][3], a vector float[3], a scalar float (nothing). Returns the dimension after.
static uint32_t attr_value_dims(int dim[MAX_NUM_DIMS], uint32_t d, const md_attribute_t* attr) {
    for (uint32_t i = 1; i < attr->format.rank && d < MAX_NUM_DIMS; ++i) {
        dim[d++] = (int)attr->format.shape[i];
    }
    if (attr->format.components > 1 && d < MAX_NUM_DIMS) {
        dim[d++] = (int)attr->format.components;
    }
    return d;
}

// Mixes what a value read from attr depends on into h: which attribute it is, its version and its frame axis's
static uint64_t attr_member_version(uint64_t h, const md_attribute_t* attr, const md_attribute_t* axis) {
    const uint64_t v[3] = { attr->id, attr->version, axis->version };
    return md_hash64(v, sizeof(v), h);
}

// What the members of a pattern's array are called: their labels (the names their file gives them, "Coul-SR:Protein-SOL")
// when every one has one and no two the same, otherwise the segments of their paths that the pattern's '*' stood in
// ("coul_sr_protein_sol" of "edr/coul_sr_*"), the part of them that differs. subjects are the paths the pattern was
// matched against (relative to the run where the pattern is).
static void attr_member_labels(attr_member_t* members, const str_t* subjects, size_t count, str_t pattern, const md_attributes_t* attributes, md_allocator_i* alloc) {
    bool own = true;
    for (size_t i = 0; i < count && own; ++i) {
        const md_attribute_t* a = md_attributes_get(attributes, members[i].id);
        own = a && !str_empty(a->label);
        for (size_t j = 0; j < i && own; ++j) {
            own = !str_eq(members[j].label, a->label);
        }
        if (own) members[i].label = str_copy(a->label, alloc);
    }
    if (!own) {
        for (size_t i = 0; i < count; ++i) {
            members[i].label = glob_varying_segments(pattern, subjects[i], alloc);
        }
    }
}

// attr("edr/coul_sr_*"): every attribute the pattern matches that attr() reads, one array over them in path order, of
// float[n] for n scalars, float[n][3] for n vectors. Absolute or relative to a run like a path, and all of the matches
// from one run. A frame axis is never one: it is the time the others are sampled at, not one of them. The matches have
// to agree on what one value is and on its unit to make one array, and each is read on its own frame axis.
static bool static_check_attribute_pattern(ast_node_t* node, token_t tok, str_t pattern, eval_context_t* ctx) {
    const md_attributes_t* attributes = &ctx->sys->attributes;
    md_allocator_i* arena = ctx->compile_ir->arena;
    const bool absolute = str_begins_with(pattern, STR_LIT("run/"));

    md_array(attr_member_t) members = 0;
    md_array(str_t) subjects = 0;           // what the pattern matched in each: the path, or the part below the run
    md_array(str_t) runs = 0;               // the runs they are from, each once
    const md_attribute_t* run_axis = NULL;  // of the first
    size_t num_unreadable = 0;
    size_t num_axes = 0;

    for (md_attribute_iter_t it = md_attributes_iter(attributes, STR_LIT("run")); md_attributes_next(&it);) {
        const md_attribute_t* a = it.attr;
        if (!str_begins_with(a->path, STR_LIT("run/"))) continue;
        str_t subject = a->path;
        if (!absolute) {
            size_t loc;
            const str_t below = str_substr(a->path, 4, SIZE_MAX);
            if (!str_find_char(&loc, below, '/')) continue;
            subject = str_substr(below, loc + 1, SIZE_MAX);
        }
        if (!glob_match_path(pattern, subject)) continue;

        attr_read_info_t info;
        if (attr_readable(attributes, a, &info) != ATTR_READABLE) {
            num_unreadable += 1;
            continue;
        }
        if (md_attributes_axis(attributes, a) == a) {
            num_axes += 1;
            continue;
        }

        if (!md_array_size(runs) || !str_eq(runs[md_array_size(runs) - 1], info.run)) {
            md_array_push(runs, info.run, ctx->temp_alloc);
        }
        if (!run_axis) run_axis = info.run_axis;
        attr_member_t m = { .path = str_copy(a->path, arena), .id = a->id, .axis_id = info.attr_axis->id };
        md_array_push(members, m, arena);
        md_array_push(subjects, subject, ctx->temp_alloc);
    }

    const size_t num = md_array_size(members);
    if (num == 0) {
        if (num_axes && !num_unreadable) {
            EVAL_LOG_ERROR(ctx, tok, "attr: '"STR_FMT"' matches only frame axes, the times the others are sampled at, which a pattern leaves out", STR_ARG(pattern));
        } else if (num_unreadable) {
            EVAL_LOG_ERROR(ctx, tok, "attr: '"STR_FMT"' matches %zu %s, but none that attr() reads: numbers that vary over the frames of a run",
                STR_ARG(pattern), num_unreadable, num_unreadable == 1 ? "attribute" : "attributes");
        } else {
            EVAL_LOG_ERROR(ctx, tok, "attr: nothing in any run matches '"STR_FMT"'", STR_ARG(pattern));
        }
        return false;
    }
    if (md_array_size(runs) > 1) {
        // The pattern below the run, which each run is suggested with
        str_t below = pattern;
        size_t loc;
        if (absolute && str_find_char(&loc, str_substr(pattern, 4, SIZE_MAX), '/')) {
            below = str_substr(pattern, 4 + loc + 1, SIZE_MAX);
        }
        md_str_writer_t w = md_str_writer_create(ctx->temp_alloc);
        for (size_t i = 0; i < md_array_size(runs); ++i) {
            md_str_writer_fmt(&w, "  '"STR_FMT"/"STR_FMT"'\n", STR_ARG(runs[i]), STR_ARG(below));
        }
        EVAL_LOG_ERROR(ctx, tok, "attr: '"STR_FMT"' matches in more than one run, name the run:\n"STR_FMT,
            STR_ARG(pattern), STR_ARG(md_str_writer_to_str(&w)));
        return false;
    }

    // One value of the array: the members' axis, then what one value of a member is, the same for all of them
    const md_attribute_t* first = md_attributes_get(attributes, members[0].id);
    const uint32_t member_dims = (first->format.rank - 1) + (first->format.components > 1 ? 1 : 0);
    if (1 + member_dims > MAX_NUM_DIMS) {
        EVAL_LOG_ERROR(ctx, tok, "attr: an array of what '"STR_FMT"' matches has more dimensions than a script value can hold (%d)", STR_ARG(pattern), MAX_NUM_DIMS);
        return false;
    }
    type_info_t type = { .base_type = TYPE_FLOAT, .dim = { (int)num } };
    attr_value_dims(type.dim, 1, first);
    uint64_t version = attr_member_version(0, first, md_attributes_get(attributes, members[0].axis_id));
    for (size_t i = 1; i < num; ++i) {
        const md_attribute_t* a = md_attributes_get(attributes, members[i].id);
        type_info_t t = { .base_type = TYPE_FLOAT, .dim = { (int)num } };
        attr_value_dims(t.dim, 1, a);
        if (MEMCMP(t.dim, type.dim, sizeof(type.dim)) != 0 || !md_unit_equal(a->unit, first->unit)) {
            char ta[64], tb[64], ua[32], ub[32];
            type_info_t one_a = type, one_b = t;
            // What one member is: a scalar prints as float
            dim_shift_left(one_a.dim);
            dim_shift_left(one_b.dim);
            print_type_info(ta, sizeof(ta), one_a);
            print_type_info(tb, sizeof(tb), one_b);
            md_unit_print(ua, sizeof(ua), first->unit);
            md_unit_print(ub, sizeof(ub), a->unit);
            EVAL_LOG_ERROR(ctx, tok, "attr: what '"STR_FMT"' matches does not make one array: '"STR_FMT"' is %s (%s), '"STR_FMT"' is %s (%s)",
                STR_ARG(pattern), STR_ARG(members[0].path), ta, ua[0] ? ua : "no unit", STR_ARG(members[i].path), tb, ub[0] ? ub : "no unit");
            return false;
        }
        version = attr_member_version(version, a, md_attributes_get(attributes, members[i].axis_id));
    }
    attr_member_labels(members, subjects, num, pattern, attributes, arena);

    node->attr_path    = str_copy(pattern, arena);
    node->attr_members = members;
    node->attr_pattern = true;
    node->run_axis_id  = run_axis->id;
    node->attr_version = md_hash64(&run_axis->version, sizeof(run_axis->version), version);

    node->data.type    = type;
    node->data.unit[0] = md_unit_none();
    node->data.unit[1] = first->unit;
    node->flags       |= FLAG_DYNAMIC;  // a different value at every frame

    return true;
}

// attr("<path>") names a temporal attribute of the system. The path is either absolute
// ("run/a/edr/potential") or relative to a run ("edr/potential"), in which case it has to exist in
// exactly one run: with a single trajectory that is always the case, and when it is not the error
// says which runs to choose between rather than picking one. A path with a '*' in it is a pattern
// (static_check_attribute_pattern), which reads every attribute it matches.
static bool static_check_attribute(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ATTRIBUTE);
    ASSERT(ctx);

    if (!static_check_children(node, ctx)) {
        return false;
    }

    ast_node_t** args = node->children;
    if (md_array_size(args) != 1 || args[0]->type != AST_CONSTANT_VALUE || !is_type_equivalent(args[0]->data.type, (type_info_t)TI_STRING)) {
        EVAL_LOG_ERROR(ctx, node->token, "attr: expected one constant string, the path of the attribute, e.g. attr(\"edr/potential\")");
        return false;
    }
    if (!ctx->sys) {
        EVAL_LOG_ERROR(ctx, node->token, "attr: there is no system to read attributes from");
        return false;
    }

    const md_attributes_t* attributes = &ctx->sys->attributes;
    const str_t path = args[0]->value._string;
    if (str_empty(path)) {
        EVAL_LOG_ERROR(ctx, args[0]->token, "attr: empty path");
        return false;
    }

    // Resolved anew: what was read before is no longer what is read
    node->attr_members = NULL;
    node->attr_pattern = false;

    if (str_has_glob(path)) {
        return static_check_attribute_pattern(node, args[0]->token, path, ctx);
    }

    const md_attribute_t* attr = md_attributes_find(attributes, path);
    if (!attr) {
        // Relative to a run: look in every one and insist on exactly one match.
        size_t num_matches = 0;
        md_str_writer_t candidates = md_str_writer_create(ctx->temp_alloc);
        for (md_attribute_iter_t run = md_attributes_iter_children(attributes, STR_LIT("run")); md_attributes_next(&run);) {
            const md_attribute_t* a = md_attributes_find_in(attributes, run.child_path, path);
            if (a) {
                attr = a;
                num_matches += 1;
                md_str_writer_fmt(&candidates, "  '"STR_FMT"'\n", STR_ARG(a->path));
            }
        }
        if (num_matches == 0) {
            EVAL_LOG_ERROR(ctx, args[0]->token, "attr: no attribute '"STR_FMT"' in the system or in any run", STR_ARG(path));
            return false;
        }
        if (num_matches > 1) {
            EVAL_LOG_ERROR(ctx, args[0]->token, "attr: '"STR_FMT"' exists in more than one run, name it in full:\n"STR_FMT,
                STR_ARG(path), STR_ARG(md_str_writer_to_str(&candidates)));
            return false;
        }
    }

    attr_read_info_t info;
    const str_t full = attr->path;
    switch (attr_readable(attributes, attr, &info)) {
    case ATTR_READABLE:
        break;
    case ATTR_NOT_TEMPORAL:
        EVAL_LOG_ERROR(ctx, args[0]->token, "attr: '"STR_FMT"' does not vary over the trajectory; attr() reads temporal attributes", STR_ARG(full));
        return false;
    case ATTR_NOT_NUMERIC:
        EVAL_LOG_ERROR(ctx, args[0]->token, "attr: '"STR_FMT"' does not hold numbers", STR_ARG(full));
        return false;
    case ATTR_NOT_IN_RUN:
        EVAL_LOG_ERROR(ctx, args[0]->token, "attr: '"STR_FMT"' is not part of a run ('run/<name>/...'), so it has no frames to follow", STR_ARG(full));
        return false;
    case ATTR_RUN_WITHOUT_AXIS:
        EVAL_LOG_ERROR(ctx, args[0]->token, "attr: the run '"STR_FMT"' has no frame axis ('"STR_FMT"/time')", STR_ARG(info.run), STR_ARG(info.run));
        return false;
    case ATTR_WITHOUT_AXIS:
        EVAL_LOG_ERROR(ctx, args[0]->token, "attr: '"STR_FMT"' has no frame axis", STR_ARG(full));
        return false;
    case ATTR_TOO_MANY_DIMS:
        EVAL_LOG_ERROR(ctx, args[0]->token, "attr: '"STR_FMT"' has more dimensions per frame than a script value can hold (%d)", STR_ARG(full), MAX_NUM_DIMS);
        return false;
    }
    const md_attribute_t* run_axis  = info.run_axis;
    const md_attribute_t* attr_axis = info.attr_axis;

    // One value per frame. A scalar is float.
    type_info_t type = { .base_type = TYPE_FLOAT, .dim = {1} };
    attr_value_dims(type.dim, 0, attr);

    attr_member_t member = { .path = str_copy(full, ctx->compile_ir->arena), .id = attr->id, .axis_id = attr_axis->id };
    md_array_push(node->attr_members, member, ctx->compile_ir->arena);
    node->attr_path    = member.path;
    node->run_axis_id  = run_axis->id;
    node->attr_version = md_hash64(&run_axis->version, sizeof(run_axis->version), attr_member_version(0, attr, attr_axis));

    node->data.type    = type;
    node->data.unit[0] = md_unit_none();
    node->data.unit[1] = attr->unit;
    node->flags       |= FLAG_DYNAMIC;  // a different value at every frame

    return true;
}

static bool static_check_flatten(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(node->type == AST_FLATTEN);
    ASSERT(ctx);

    size_t num_args = md_array_size(node->children);
    if (num_args != 1) {
        EVAL_LOG_ERROR(ctx, node->token, "Expected 1 argument to built in procedure flatten, got %zu", num_args);
        return false;
    }

    uint32_t flags = ctx->eval_flags;
    ctx->eval_flags |= EVAL_FLAG_FLATTEN;
    const bool children_ok = static_check_children(node, ctx);
    ctx->eval_flags = flags;
    if (!children_ok) {
        return false;
    }

    const ast_node_t* arg = node->children[0];
    if (arg->data.type.base_type == TYPE_UNDEFINED) {
        EVAL_LOG_ERROR(ctx, arg->token, "Argument to built in procedure 'flatten' has undefined base type");
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

    if (node->data.type.base_type == TYPE_ATOMSET) {
        // If we are dealing with atomsets, we flatten even further by collapsing it into a single atomset. That is not
        // a change of shape but their union (an atom in two of them is in the result once), which union() names.
        node->data.type.dim[0] = 1;
        EVAL_LOG_WARNING(ctx, node->token, "flatten() of selections is their union: use union(), flatten() is for the shape of arrays");
    }

    return true;
}

static bool static_check_transpose(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(node->type == AST_TRANSPOSE);
    ASSERT(ctx);

    size_t num_args = md_array_size(node->children);
    if (num_args != 1) {
        EVAL_LOG_ERROR(ctx, node->token, "Expected 1 argument to built in procedure 'transpose', got %zu", num_args);
        return false;
    }

    if (!static_check_children(node, ctx)) {
        return false;
    }

    const ast_node_t* arg = node->children[0];
    if (arg->data.type.base_type == TYPE_UNDEFINED) {
        EVAL_LOG_ERROR(ctx, arg->token, "Argument to built in procedure 'transpose' has undefined base type");
        return false;
    }

    int ndim = dim_ndims(arg->data.type.dim);
    if (ndim == 0) {
        EVAL_LOG_ERROR(ctx, arg->token, "Undefined number of dimensions of argument supplied to transpose");
        return false;
    }

    if (arg->data.type.dim[0] == -1) {
        EVAL_LOG_ERROR(ctx, arg->token, "Built in procedure 'transpose' does not support dynamic length types");
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
        EVAL_LOG_ERROR(ctx, arg->token, "objects of dim > 2 are currently not supported in transpose operation");
        return false;
    }
}

static const proc_sig_t* find_proc_signature(str_t name) {
    for (size_t i = 0; i < ARRAY_SIZE(signatures); ++i) {
        if (str_eq(signatures[i].proc, name)) {
            return &signatures[i];
        }
    }
    return NULL;
}

static int find_param_index(const proc_sig_t* sig, str_t name) {
    for (size_t i = 0; i < sig->num_params; ++i) {
        if (str_eq(sig->param[i].name, name)) {
            return (int)i;
        }
    }
    return -1;
}

static ast_node_t* create_absent_argument_node(md_script_ir_t* ir, token_t token) {
    ast_node_t* node = create_node(ir, AST_CONSTANT_VALUE, token);
    node->data.type = (type_info_t){0};
    node->flags |= FLAG_IMPLICIT;
    return node;
}

static ast_node_t* create_default_argument_node(md_script_ir_t* ir, const param_sig_t* param, token_t token) {
    ASSERT(param->flags & PARAM_DEFAULT);
    ASSERT(param->def_type.base_type != TYPE_UNDEFINED && is_scalar(param->def_type));
    ast_node_t* node = create_node(ir, AST_CONSTANT_VALUE, token);
    node->data.type = param->def_type;
    node->value     = param->def;
    node->flags    |= FLAG_IMPLICIT;
    return node;
}

// Rearranges the arguments of a procedure call into the positional order declared by sig, supplying defaults
// for omitted parameters, and clears node->named_args. Performed before overload resolution, which therefore
// only ever sees positional arguments. Idempotent: binding an already bound call leaves it as it is.
static bool bind_arguments(ast_node_t* node, const proc_sig_t* sig, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_PROC_CALL);
    ASSERT(sig);

    const size_t num_args = md_array_size(node->children);
    const named_arg_t* named = node->named_args;
    ASSERT(!named || md_array_size(named) == num_args);
    node->named_args = NULL;

    ast_node_t* bound[MAX_SUPPORTED_PROC_ARGS] = {0};
    ASSERT(sig->num_params <= MAX_SUPPORTED_PROC_ARGS);

    for (size_t i = 0; i < num_args; ++i) {
        size_t p = 0;
        if (!named || str_empty(named[i].name)) {
            // The parser guarantees that positional arguments precede named ones, so the position is the index
            p = i;
            if (p >= sig->num_params) {
                EVAL_LOG_ERROR(ctx, node->children[i]->token, "Too many arguments, '"STR_FMT"' takes at most %i", STR_ARG(sig->proc), (int)sig->num_params);
                return false;
            }
            if (sig->param[p].flags & PARAM_KW_ONLY) {
                EVAL_LOG_ERROR(ctx, node->children[i]->token, "The argument '"STR_FMT"' of '"STR_FMT"' can only be given by name", STR_ARG(sig->param[p].name), STR_ARG(sig->proc));
                return false;
            }
        } else {
            const int idx = find_param_index(sig, named[i].name);
            if (idx < 0) {
                md_str_writer_t sb = md_str_writer_create(ctx->temp_alloc);
                md_str_writer_fmt(&sb, "'"STR_FMT"' has no parameter named '"STR_FMT"', valid names are: ", STR_ARG(sig->proc), STR_ARG(named[i].name));
                for (size_t j = 0; j < sig->num_params; ++j) {
                    md_str_writer_push_str(&sb, sig->param[j].name);
                    if (j + 1 < sig->num_params) {
                        md_str_writer_push_str(&sb, STR_LIT(", "));
                    }
                }
                EVAL_LOG_ERROR_STR(ctx, named[i].token, md_str_writer_to_str(&sb));
                return false;
            }
            p = (size_t)idx;
            if (bound[p]) {
                // Duplicate names are rejected by the parser, so the earlier binding was positional
                EVAL_LOG_ERROR(ctx, named[i].token, "The argument '"STR_FMT"' is already given by position", STR_ARG(named[i].name));
                return false;
            }
        }
        bound[p] = node->children[i];
    }

    for (size_t p = 0; p < sig->num_params; ++p) {
        if (!bound[p] && (sig->param[p].flags & PARAM_DEFAULT)) {
            bound[p] = create_default_argument_node(ctx->compile_ir, &sig->param[p], node->token);
        } else if (!bound[p] && (sig->param[p].flags & PARAM_NULLABLE)) {
            bound[p] = create_absent_argument_node(ctx->compile_ir, node->token);
        }
    }

    // The bound arguments must be contiguous: anything after the last one is left to an overload taking fewer arguments
    size_t count = 0;
    for (size_t p = 0; p < sig->num_params; ++p) {
        if (bound[p]) count = p + 1;
    }
    for (size_t p = 0; p < sig->num_params; ++p) {
        if (bound[p]) continue;
        const param_sig_t* param = &sig->param[p];
        if (!(param->flags & PARAM_OPTIONAL)) {
            EVAL_LOG_ERROR(ctx, node->token, "Missing argument '"STR_FMT"' in call to '"STR_FMT"'", STR_ARG(param->name), STR_ARG(sig->proc));
            return false;
        }
        if (p < count) {
            EVAL_LOG_ERROR(ctx, node->token, "The argument '"STR_FMT"' of '"STR_FMT"' cannot be omitted when a later argument is given", STR_ARG(param->name), STR_ARG(sig->proc));
            return false;
        }
    }

    if (named || count != num_args) {
        md_array(ast_node_t*) children = 0;
        md_array_push_array(children, bound, count, ctx->compile_ir->arena);
        node->children = children;
    }

    return true;
}

static bool static_check_proc_call(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node->type == AST_PROC_CALL);
    uint32_t backup_flags = ctx->eval_flags;

    bool result = true;
    // @NOTE: We do not want to overwrite the procedure if already assigned, as it might have been assigned as an operator (which is a proc call)
    if (!node->proc) {
        // Bind the arguments into positional order before the procedure is resolved.
        // This runs for positional calls as well, since those may leave parameters to their defaults.
        const proc_sig_t* sig = find_proc_signature(node->ident);
        if (sig) {
            if (!bind_arguments(node, sig, ctx)) {
                return false;
            }
        } else if (node->named_args) {
            for (size_t i = 0; i < md_array_size(node->named_args); ++i) {
                if (!str_empty(node->named_args[i].name)) {
                    EVAL_LOG_ERROR(ctx, node->named_args[i].token, "'"STR_FMT"' does not accept named arguments", STR_ARG(node->ident));
                    break;
                }
            }
            node->named_args = NULL;
            return false;
        }

        // No point in letting expressions statically evaluate and store its data within the tree at this point
        // Conversions and other stuff may occur later
        uint32_t flags = ctx->eval_flags;
        ctx->eval_flags |= EVAL_FLAG_NO_STATIC_EVAL | 0;
        result = static_check_children(node, ctx);
        ctx->eval_flags = flags;
        if (result) {
            const size_t num_args = md_array_size(node->children);
            const str_t proc_name = node->ident;

            // No procedure takes this many, and the arguments would not fit below (nor where it is evaluated)
            if (num_args >= MAX_SUPPORTED_PROC_ARGS) {
                EVAL_LOG_ERROR(ctx, node->token, "Too many arguments (%zu) in call to '"STR_FMT"', no procedure takes more than %i",
                    num_args, STR_ARG(proc_name), MAX_SUPPORTED_PROC_ARGS - 1);
                ctx->eval_flags = backup_flags;
                return false;
            }
            type_info_t arg_type[MAX_SUPPORTED_PROC_ARGS] = { 0 };
            extract_argument_types(arg_type, ARRAY_SIZE(arg_type), node);

            procedure_match_result_t res = find_procedure_supporting_arg_types(proc_name, arg_type, num_args, true);
            result = result && res.success;

            if (res.success) {
                node->proc = res.procedure;
                node->proc_flags = res.flags;
            } else {
                if (num_args == 0) {
                    EVAL_LOG_ERROR(ctx, node->token,
                        "Could not find matching procedure '"STR_FMT"' which takes no arguments", STR_ARG(proc_name));
                } else {
                    char buf[512];
                    print_argument_list(buf, ARRAY_SIZE(buf), arg_type, num_args);
                    bool bonds = false;
                    for (size_t i = 0; i < num_args; ++i) bonds |= arg_type[i].base_type == TYPE_BONDSET;

                    EVAL_LOG_ERROR(ctx, node->token,
                        "Could not find matching procedure '"STR_FMT"' which takes the following argument(s): %s%s",
                        STR_ARG(proc_name), buf, bonds ? BONDS_DO_NOT_MIX_HINT : "");
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
        
        // Perform new child check here since children may have changed due to conversions etc.
        // Also allow for any static evaluation to occur and for length to be determined etc.
        result = result && static_check_children(node, ctx);
            
        const size_t num_args = md_array_size(node->children);
        ast_node_t** args = node->children;
        for (size_t i = 0; i < num_args; ++i) {
            if (node->proc->arg_type[i].base_type == TYPE_COORDINATE && !is_absent_type(args[i]->data.type)) {
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

        if (result && (node->proc->flags & FLAG_UNIT_OF_ARG) && md_array_size(node->children) > 0) {
            node->data.unit[0] = node->children[0]->data.unit[0];
            node->data.unit[1] = node->children[0]->data.unit[1];
        }

        // A procedure which states nothing about it during its static validation (below) is periodic by the rules of
        // period_of_proc_call
        if (result) {
            node->data.period = period_of_proc_call(node);
        }

        // Perform static validation by evaluating with NULL ptr
        if (result && node->proc->flags & FLAG_STATIC_VALIDATION) {
            //MD_LOG_DEBUG("Static validating procedure: '"STR_FMT"'", STR_ARG(node->ident));
            static_backchannel_t* prev_channel = ctx->backchannel;
            static_backchannel_t channel = {0};
            ctx->backchannel = &channel;
            result = result && evaluate_proc_call(NULL, node, ctx);
            node->data.unit[0] = channel.unit[0];
            node->data.unit[1] = channel.unit[1];
            node->data.value_range = channel.value_range;
            node->data.period = channel.period;
            ctx->backchannel = prev_channel;
        }
    }

    ctx->eval_flags = backup_flags;
    return result;
}

static bool static_check_constant_value(ast_node_t* node, eval_context_t* ctx) {
    (void)ctx; 
    ASSERT(node && node->type == AST_CONSTANT_VALUE);
    if (is_absent_type(node->data.type)) {
        // An absent argument (see PARAM_NULLABLE): no data
        node->data.ptr  = NULL;
        node->data.size = 0;
        node->flags |= FLAG_CONSTANT;
        return true;
    }
    if (!INTERNAL_CHECK(ctx, node->token, is_scalar(node->data.type))) return false;

    if (node->data.type.base_type == TYPE_IRANGE) {
        // Validate ascending format and positive step
        irange_t rng = node->value._irange;
        if (rng.beg > rng.end) {
            EVAL_LOG_ERROR(ctx, node->token, "The range is invalid, a range must have an ascending format, did you mean '%i:%i'?", rng.end, rng.beg);
            return false;
        }
        if (rng.step <= 0) {
            EVAL_LOG_ERROR(ctx, node->token, "The stride in a range must be a positive integer, got %i", rng.step);
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
        md_unit_t array_unit[2] = {md_unit_none(), md_unit_none()};
        
        // Pass 1: find a suitable base_type for the array
        for (size_t i = 0; i < num_elem; ++i) {
            int elem_len = type_info_array_len(elem[i]->data.type);
            type_info_t elem_type = type_info_element_type(elem[i]->data.type);
            if (is_variable_length(elem_type)) {
                array_len = -1;
            }

            /*
            // Normalize the type and dimensions here for the length to be computed correctly
            if (elem_type.base_type != TYPE_ATOMSET) {
                int ndim = dim_ndims(elem_type.dim);
                if (ndim == 1 || (ndim > 1 && elem_type.dim[0] != 1)) {
                    dim_shift_right(elem_type.dim);
                }
                // Turn the normalized form into a single entity
                elem_type.dim[0] = 1;
            }
            */
            if (elem_type.base_type == TYPE_ATOMSET) {
                // Atomset
                dim_flatten(elem_type.dim);
            }

            if (array_len >= 0) array_len += elem_len;

            if (array_type.base_type == TYPE_UNDEFINED) {
                // Assign type from first element
                array_type = elem_type;
            }
            else {
                if (array_type.base_type == TYPE_ATOMSET) {
                    // Atomsets are handled separately and concatenated
                    if (elem_type.base_type != TYPE_ATOMSET) {
                        EVAL_LOG_ERROR(ctx, elem[i]->token, "Incompatible types wihin array construct");
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
                                    EVAL_LOG_ERROR(ctx, elem[i]->token, "Incompatible types wihin array construct");
                                    return false;
                                }
                            }
                        }
                        else {
                            // Incompatible types...
                            EVAL_LOG_ERROR(ctx, elem[i]->token, "Incompatible types wihin array construct");
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

            if (!is_type_directly_compatible(elem[i]->data.type, converted_type)) {
                if (!is_type_implicitly_convertible(elem[i]->data.type, converted_type)) {
                    // An element accepted before the type of the array was promoted, which cannot follow it: the
                    // int[2] of '{{1,2},{3,4:5}}'. Taken as it is, it would not fill its part of the array.
                    EVAL_LOG_ERROR(ctx, elem[i]->token, "Incompatible types wihin array construct");
                    return false;
                }
                if (!convert_node(elem[i], converted_type, ctx)) {
                    return false;
                }
            }
        }

        // Pass 3: Deduce the unit, if defined and the same.
        if (num_elem > 0) {
            array_unit[0] = elem[0]->data.unit[0];
            array_unit[1] = elem[0]->data.unit[1];
            for (size_t i = 1; i < num_elem; ++i) {
                if (!md_unit_equal(elem[i]->data.unit[0], array_unit[0])) {
                    array_unit[0] = md_unit_none();
                    break;
                }
                if (!md_unit_equal(elem[i]->data.unit[1], array_unit[1])) {
                    array_unit[1] = md_unit_none();
                    break;
                }
            }
        }

        // Finalize the type of the array (dimensionality)
        array_type.dim[0] = (int32_t)array_len;

        // Volumes are stacked rather than concatenated: '{sdf(a), sdf(b)}' is float[2][128][128][128], two volumes,
        // not float[256][128][128], one twice as long along x. The data is the same, one volume after the other: only
        // the type says which it is.
        int num_volumes = 0;
        for (size_t i = 0; i < num_elem; ++i) {
            const int n = volume_count(elem[i]->data.type);
            if (n == 0) {
                num_volumes = 0;
                break;
            }
            num_volumes += n;
        }
        if (num_volumes > 0) {
            array_type = (type_info_t){TYPE_FLOAT, {num_volumes, MD_VOL_DIM, MD_VOL_DIM, MD_VOL_DIM}};
        }
        node->data.type = array_type;
        node->data.unit[0] = array_unit[0];
        node->data.unit[1] = array_unit[1];

        // Periodic if every element is, in periods of the same length: the period of the first. Elements given in
        // other periods of that length are still the same a period apart.
        node->data.period = num_elem > 0 ? elem[0]->data.period : (period_t){0};
        for (size_t i = 1; i < num_elem; ++i) {
            if (!period_valid(elem[i]->data.period) || !period_same_length(elem[i]->data.period, node->data.period)) {
                node->data.period = (period_t){0};
                break;
            }
        }

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
        EVAL_LOG_ERROR(ctx, node->token, "Missing arguments in array subscript");
        return false;
    }

    uint32_t eval_flags = ctx->eval_flags;
    ctx->eval_flags = ctx->eval_flags & ~(EVAL_FLAG_FLATTEN);
    const bool lhs_ok = static_check_node(elem[0], ctx);
    ctx->eval_flags = eval_flags;
    if (!lhs_ok) return false;

    for (size_t i = 1; i < num_elem; ++i) {
        if (!static_check_node(elem[i], ctx)) return false;
    }

    ast_node_t*  lhs  = elem[0];
    ast_node_t** args = elem + 1;
    const size_t num_args = num_elem - 1;

    if (is_variable_length(lhs->data.type)) {
        EVAL_LOG_ERROR(ctx, lhs->token, "Array subscript operator can only be applied to expressions which have a static length");
        return false;
    }

    size_t num_dim = (size_t)dim_ndims(lhs->data.type.dim);
    if (num_args > 1 && num_args != num_dim) {
    	EVAL_LOG_ERROR(ctx, elem[1]->token, "Invalid number of arguments (%i) in array subscript, expected number of arguments to match number of dimensions of lhs (%i)", (int)num_args, num_dim);
		return false;
    }

    type_info_t result_type = lhs->data.type;
    for (size_t i = 0; i < num_args; ++i) {
        int lhs_dim = lhs->data.type.dim[i];
        if (lhs_dim <= 0) {
            EVAL_LOG_ERROR(ctx, lhs->token, "Unexpected length (%i) in dimension (%i) of lhs", lhs_dim, (int)i);
            return false;
        }
        ast_node_t* arg = args[i];
        if (arg->flags & FLAG_DYNAMIC) {
            EVAL_LOG_ERROR(ctx, args[i]->token, "Only static expressions are allowed within array subscript");
            return false;
        }

        if (is_scalar(arg->data.type)) {
            if (arg->data.type.base_type == TYPE_INT || arg->data.type.base_type == TYPE_IRANGE) {
                irange_t range = {0};
                if (arg->data.type.base_type == TYPE_INT) {
                    range.beg = arg->value._int;
                    range.end = arg->value._int;
                    range.step = 1;
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
                    if (range.step <= 0) range.step = 1;
                }
                if (range.beg <= range.end && 1 <= range.beg && range.end <= lhs_dim) {
                    // Convert to 0-based half-open interval [beg0, end0)
                    irange_t sr = {0};
                    sr.beg  = range.beg - 1;
                    sr.end  = range.end;
                    sr.step = range.step;
                    node->subscript_ranges[i] = sr;
                    // Compute resulting length with stride: floor((end-beg)/step)+1
                    int cnt = (range.end - range.beg) / MAX(range.step,1) + 1;
                    result_type.dim[i] = cnt;
                } else {
                    EVAL_LOG_ERROR(ctx, arg->token, "Invalid Array subscript range");
                    return false;
                }
            } else {
                EVAL_LOG_ERROR(ctx, arg->token, "Only int and int-ranges are allowed inside array subscript");
                return false;
            }
        } else {
            EVAL_LOG_ERROR(ctx, arg->token, "No arrays are allowed inside array subscript");
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

    // Distributions keep their count: [1][bins] for one, never [bins], which would read as that many of them
    if (result_type.base_type == TYPE_DISTRIBUTION && dim_ndims(result_type.dim) == 1) {
        dim_shift_right(result_type.dim);
    }

    // SUCCESS!
    // @NOTE Subscript dim is not just the number of args supplied,
    // It is the number of dimensions
    node->subscript_dim = num_dim;
    node->flags = (lhs->flags & FLAG_AST_PROPAGATION_MASK);
    node->data.type = result_type;
    node->data.unit[0] = lhs->data.unit[0];
    node->data.unit[1] = lhs->data.unit[1];
    node->data.period = lhs->data.period;    // An element of a periodic array is periodic alike
    return true;
}

static bool static_check_identifier_reference(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_IDENTIFIER);
    
    if (str_empty(node->ident)) {
        node->ident = node->token.str;
    }

    identifier_t* ident = get_identifier(ctx->compile_ir, node->ident);
    if (ident && ident->node) {
        if (ident->node->data.type.base_type == TYPE_UNDEFINED) {
            EVAL_LOG_ERROR(ctx, node->token, "Identifier ("STR_FMT") has an unresolved type", ident->name.len, ident->name.ptr);
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
        EVAL_LOG_ERROR(ctx, node->token, "Unresolved reference to identifier ("STR_FMT")", node->ident.len, node->ident.ptr);
    }
    return false;
}

static bool static_check_assignment(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_ASSIGNMENT);
    if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 2 && node->children[0] && node->children[1])) return false;
    ast_node_t* lhs = node->children[0];
    ast_node_t* rhs = node->children[1];

    int num_idents = 1;
    ast_node_t** idents = &lhs;

    if (lhs->type == AST_ARRAY) {
        for (size_t i = 0; i < md_array_size(lhs->children); ++i) {
			if (lhs->children[i]->type != AST_IDENTIFIER) {
				EVAL_LOG_ERROR(ctx, node->token, "An element on the left hand side of assignment is not an identifier");
				return false;
			}
		}
        num_idents = (int)md_array_size(lhs->children);
        idents = lhs->children;
    } else if (lhs->type != AST_IDENTIFIER) {
        EVAL_LOG_ERROR(ctx, node->token, "Left hand side of assignment is not an identifier");
		return false;
    }

    for (int i = 0; i < num_idents; ++i) {
        str_t ident = idents[i]->ident;
        // Assignment, therefore a new identifier
        if (find_constant(ident)) {
            EVAL_LOG_ERROR(ctx, idents[i]->token, "The identifier is occupied by a constant and cannot be assigned.");
            return false;
        }
        else if (get_identifier(ctx->compile_ir, ident)) {
            EVAL_LOG_ERROR(ctx, idents[i]->token, "The identifier is already taken. Variables cannot be reassigned.");
            return false;
        }
    }

    if (static_check_node(rhs, ctx)) {
        // Right hand side type must be known
        if (!INTERNAL_CHECK(ctx, rhs->token, rhs->data.type.base_type != TYPE_UNDEFINED)) return false;
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
                EVAL_LOG_ERROR(ctx, node->token, "Assignment mismatch between left and right hand side: The right hand side has a variable length");
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
                EVAL_LOG_ERROR(ctx, node->token, "Assignment mismatch between left and right hand side: The number of identifiers on the left hand side must match the length of the right hand side");
                return false;
            }

            for (int i = 0; i < num_idents; ++i) {
                // Identifiers type should always be undefined until explicitly assigned.
                if (!INTERNAL_CHECK(ctx, idents[i]->token, idents[i]->data.type.base_type == TYPE_UNDEFINED)) return false;
                identifier_t* ident = create_identifier(ctx->compile_ir, idents[i]->ident);
                if (!ident) {
                    EVAL_LOG_ERROR(ctx, node->token, "Failed to create identifier. Is the identifier already taken?");
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
                    ident->is_element = true;
                    ident->elem_idx = (uint32_t)i;
                    idents[i]->data.type = type_info_element_type(rhs->data.type);
                    const int stride = (int)type_info_element_byte_stride(rhs->data.type);
                    ident->data = md_alloc(ctx->compile_ir->arena, sizeof(data_t));
                    *ident->data = rhs->data;
                    ident->data->type = type_info_element_type(rhs->data.type);
                    dim_prune_leading_ones(ident->data->type.dim);
                    ident->data->size = stride;
                    // A constant right hand side was evaluated at compile time: the identifier is its i-th element.
                    // Without this, references to a constant identifier would read through a null pointer.
                    ident->data->ptr = (rhs->data.ptr && (rhs->flags & FLAG_CONSTANT)) ? (uint8_t*)rhs->data.ptr + (size_t)i * stride : 0;
                    ident->data->weight = ident->data->ptr ? data_weight_at(&rhs->data, (size_t)i * stride) : NULL;
                }
                if (!static_check_node(idents[i], ctx)) {
                    return false;
                }
            }
        } else {
            // Identifiers type should always be undefined until explicitly assigned.
            if (!INTERNAL_CHECK(ctx, lhs->token, lhs->data.type.base_type == TYPE_UNDEFINED)) return false;

            identifier_t* ident = create_identifier(ctx->compile_ir, lhs->ident);
            if (!ident) {
                EVAL_LOG_ERROR(ctx, node->token, "Failed to create identifier. Is the identifier already taken?");
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

        // Mark identifier within assignment and rhs
        if (num_idents == 1) {
            node->ident = lhs->ident;
            rhs->ident  = lhs->ident;
        }
        return true;
    }

    return false;
}

static void propagate_contexts(ast_node_t* node, const md_bitset_t* contexts, int64_t num_contexts) {
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
    if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 1 && node->children[0])) return false;

    ast_node_t* expr = node->children[0];

    const md_bitset_t* old_ctx = ctx->mol_ctx;
    ctx->mol_ctx = 0;
    bool result = false;
    if (static_check_node(expr, ctx)) {
        node->flags |= expr->flags & FLAG_AST_PROPAGATION_MASK;
        node->data.type = expr->data.type;
        node->data.unit[0] = expr->data.unit[0];
        node->data.unit[1] = expr->data.unit[1];
        node->data.value_range = expr->data.value_range;
        node->data.period = expr->data.period;
        result = true;
    }
    ctx->mol_ctx = old_ctx;

    return result;
}

static bool static_check_context(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node && node->type == AST_CONTEXT);
    if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 2 && node->children[0] && node->children[1])) return false;

    // The rules are as follows:
    // A context RHS must be an atomset type.
    //     
    // If the the basetype of LHS is (float, int, float[4], bool or whatever), then the result is of the same type but with an extra dimension with the same length as the RHS.
    // I.e.
    // If RHS is atomset with 8 matching structures:
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
    
    if (result && rhs->data.type.base_type == TYPE_ATOMSET) {
        if (!(rhs->flags & FLAG_DYNAMIC)) {
            const int ctx_len = type_info_array_len(rhs->data.type);
            if (ctx_len > 0) {
                const size_t num_contexts = (size_t)ctx_len;
                // Store this persistently and set this as a context for all child nodes
                md_bitset_t* contexts = 0;

                if (rhs->flags & FLAG_CONSTANT) {
                    if (!INTERNAL_CHECK(ctx, rhs->token, rhs->data.ptr && rhs->data.size == num_contexts * sizeof(md_bitset_t))) return false;
                    contexts = (md_bitset_t*)rhs->data.ptr;
                } else {
                    md_array_resize(contexts, num_contexts, ctx->compile_ir->arena);
                    MEMSET(contexts, 0, num_contexts * sizeof(md_bitset_t));
                    rhs->data.ptr = contexts;
                    rhs->data.size = num_contexts * sizeof(md_bitset_t);
                    rhs->data.alloc = ctx->compile_ir->arena;
                    result = evaluate_node(&rhs->data, rhs, ctx);
                }

                //allocate_data(&rhs->data, ctx->compile_ir->arena);
                if (result) {

                    // We differentiate here if the LHS is an atomset or not
                    // if LHS is an atomset of length 1, the resulting atomset length is the same as RHS (M)
                    // if LHS is an atomset of length N, the resulting atomset length is (N*M), N may vary for each context.

                    const bool contexts_equivalent = are_bitfields_equivalent(contexts, num_contexts, ctx->sys); 

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
                            md_log_token_t* last_error = md_array_last(ctx->compile_ir->errors);
                            if (last_error) {
                                last_error->context = &contexts[i];
                            }
                            return false;
                        }

                        type_info_t local_type = lhs->data.type;

                        if (lhs->flags & FLAG_DYNAMIC_LENGTH) {
                            if (!finalize_type(&local_type, lhs, &local_ctx)) {
                                EVAL_LOG_ERROR(ctx, lhs->token, "Failed to deduce length of type which is required for determining type and size of context");
                                return false;
                            }
                        }

                        // Normalize the type and dimensions here for the length to be computed correctly
                        if (local_type.base_type != TYPE_ATOMSET) {
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
                            EVAL_LOG_ERROR(ctx, lhs->token, "The type of the left hand side of 'in' must have a static length");
                            return false;
                        }

                        md_array_push(node->lhs_context_types, local_type, ctx->compile_ir->arena);
                        arr_len += (size_t)len;
                    }

                    node->flags |= lhs->flags & FLAG_AST_PROPAGATION_MASK;
                    node->data.type = node->lhs_context_types[0];
                    if (node->data.type.base_type == TYPE_ATOMSET) {
                        // This is a flattening of the dimensions
                        MEMSET(node->data.type.dim, 0, sizeof(node->data.type.dim));
                        node->data.type.dim[0] = (int)arr_len;
                    } else {
                        node->data.type.dim[0] = (int)arr_len;
                    }
                    node->data.unit[0] = lhs->data.unit[0];
                    node->data.unit[1] = lhs->data.unit[1];
                    node->data.value_range = lhs->data.value_range;
                    node->data.period = lhs->data.period;

                    propagate_contexts(lhs, contexts, num_contexts);
                    if (contexts_equivalent) {
                        lhs->flags |= FLAG_CONTEXTS_EQUIVALENT;
                    }

                    result = true;

                } else {
                    EVAL_LOG_ERROR(ctx, node->token, "Right hand side of 'in' failed to evaluate at compile time.");
                    result = false;
                }
            } else {
                EVAL_LOG_ERROR(ctx, node->token, "The context is empty.");
                result = false;
            }
        } else {
            EVAL_LOG_ERROR(ctx, node->token, "Right hand side of 'in' must be known at compile time.");
            result = false;
        }
    } else {
        char buf[128];
        print_type_info(buf, ARRAY_SIZE(buf), rhs->data.type);
        EVAL_LOG_ERROR(ctx, node->token, "Right hand side of keyword 'in' has an incompatible type: Expected atomset, got '%s'%s.", buf,
                       rhs->data.type.base_type == TYPE_BONDSET ? BONDS_DO_NOT_MIX_HINT : "");
        result = false;
    }

    return result;
}

static bool static_check_node(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(ctx && ctx->compile_ir);     // Part of a compilation: an evaluation never type checks
    // This is probably more like a static check which encompass more than just checking the types...
    // The idea here is that we have an syntax tree given by node and we want to map the function calls into concrete functions and check for the type.
    // The type will propagate back throughout the tree and we will remap nodes into concrete function calls.
    // This is not done in the parsing step as we want the full context before attempting this operation.

    if (!INTERNAL_CHECK(ctx, ctx->op_token, node)) return false;

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
    case AST_ATTRIBUTE:
        result = static_check_attribute(node, ctx);
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
    case AST_INT_DIV:
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
        INTERNAL_ERROR(ctx, node->token, "a node of a kind the static check does not know");
        return false;
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
        hash = hash_node(node->children[i], hash);
    }
    return hash;
}

// What a node is on its own, apart from its children: what it does and what it produces. Two nodes over the
// same children differ in this alone, '1.0 + 2.0' and '1.0 - 2.0' in their type, two overloads in their procedure.
// The procedure counts by its address: a hash is only ever compared within the process which made it. Not by
// 'ident', which the static check overwrites with the name of the variable an expression is assigned to.
static uint64_t hash_node_self(const ast_node_t* node, uint64_t seed) {
    uint64_t hash = md_hash64(&node->type, sizeof(node->type), seed);
    hash = md_hash64(&node->data.type, sizeof(node->data.type), hash);
    hash = md_hash64(&node->proc, sizeof(node->proc), hash);
    hash = md_hash64(&node->proc_flags, sizeof(node->proc_flags), hash);
    if (node->type == AST_ARRAY_SUBSCRIPT && node->subscript_dim > 0) {
        hash = md_hash64(node->subscript_ranges, node->subscript_dim * sizeof(node->subscript_ranges[0]), hash);
    }
    return hash;
}

static uint64_t hash_node(const ast_node_t* node, uint64_t seed) {
    ASSERT(node);
    switch (node->type) {
    case AST_CONSTANT_VALUE:
        switch (node->data.type.base_type) {
        case TYPE_ATOMSET:
        {
            const md_bitset_t* bs = (const md_bitset_t*)node->data.ptr;
            return md_bitset_hash64(*bs, seed);
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
    case AST_ATTRIBUTE:
        // What was read and which version of it. The data itself can be large and lives in the
        // system, whose attribute versions already say when it changed.
        return md_hash64(node->attr_path.ptr, node->attr_path.len, md_hash64(&node->attr_version, sizeof(node->attr_version), seed));
    default:
        return hash_children(node, hash_node_self(node, seed));
    }
}

// The fingerprint of what a node evaluates to. As hash_node, except for identifiers: a reference counts by the
// definition it refers to, never by its name, so renaming a variable does not change what is computed from it,
// and redefining one changes everything computed from it, however indirectly. memo holds the fingerprint of each
// node of the IR by its index, 0 for not yet computed: a definition referred to many times is hashed once.
static uint64_t dependency_hash_node(const md_script_ir_t* ir, const ast_node_t* node, uint64_t* memo) {
    ASSERT(node);
    if (memo[node->index]) {
        return memo[node->index];
    }

    uint64_t hash = 0;
    switch (node->type) {
    case AST_IDENTIFIER: {
        // A reference is linked to the expression of its declaration by the static check. A constant one (a
        // selection added from outside) is declared by a node of its own, whose contents it is.
        const identifier_t* decl = find_ir_identifier(ir, node->ident);
        const ast_node_t* def = md_array_size(node->children) == 1 ? node->children[0] : (decl ? decl->node : NULL);
        hash = def ? dependency_hash_node(ir, def, memo) : md_hash64(node->ident.ptr, node->ident.len, 0);
        hash = md_hash64(&node->data.type, sizeof(node->data.type), hash);
        if (decl && decl->is_element) {
            // One element of what the expression evaluates to: '{x, y} = com(...)'
            hash = md_hash64(&decl->elem_idx, sizeof(decl->elem_idx), hash);
        }
        break;
    }
    case AST_CONSTANT_VALUE:
    case AST_ATTRIBUTE:
        hash = hash_node(node, 0);
        break;
    default:
        hash = hash_node_self(node, 0);
        for (size_t i = 0; i < md_array_size(node->children); ++i) {
            const uint64_t child = dependency_hash_node(ir, node->children[i], memo);
            hash = md_hash64(&child, sizeof(child), hash);
        }
        break;
    }

    hash |= 1;  // 0 marks a node not hashed yet
    memo[node->index] = hash;
    return hash;
}

// The fingerprint of a dynamic identifier: what its value is computed from, whatever it is called
static uint64_t identifier_fingerprint(const md_script_ir_t* ir, const identifier_t* decl, uint64_t seed, uint64_t* memo) {
    if (!decl || !decl->node) {
        return 0;
    }
    uint64_t hash = dependency_hash_node(ir, decl->node, memo);
    if (decl->is_element) {
        hash = md_hash64(&decl->elem_idx, sizeof(decl->elem_idx), hash);
    }
    return md_hash64(&hash, sizeof(hash), seed) | 1;
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

// An assignment below node, where the parser lets one through: 'f(1, (a = 2))'. An assignment declares an
// identifier, which only a statement does: as a part of an expression it would be checked as often as the
// expression is (the arguments of a call are checked twice), and declare its identifier more than once.
static const ast_node_t* find_nested_assignment(const ast_node_t* node) {
    for (size_t i = 0; i < md_array_size(node->children); ++i) {
        const ast_node_t* child = node->children[i];
        if (!child) continue;
        if (child->type == AST_ASSIGNMENT) return child;
        const ast_node_t* nested = find_nested_assignment(child);
        if (nested) return nested;
    }
    return NULL;
}

static bool parse_script(md_script_ir_t* ir, md_allocator_i* temp_arena) {
    ASSERT(ir);
    md_temp_scope_t temp = md_temp_begin_in(temp_arena);

    bool result = true;

    tokenizer_t tokenizer = tokenizer_init(ir->str);

    parse_context_t ctx = {
        .ir = ir,
        .tokenizer = &tokenizer,
        .node = 0,
        .temp_alloc = temp_arena,
    };

    ir->stage = "Parsing";

    // Parse statements until we have exhausted the tokenizer

    token_t tok = {0};
    while (tok = tokenizer_peek_next(ctx.tokenizer), tok.type != TOKEN_END) {
        ctx.node = 0;
        ir->record_log = true;   // We reset the error recording flag before each statement
        const size_t num_errors = md_array_size(ir->errors);
        ast_node_t* raw_node = parse_expression(&ctx);
        // A tree which comes back from a statement that logged an error is not what was written: a call whose
        // argument list failed to parse, say, comes back as a call without arguments. It is not checked or run.
        ast_node_t* pruned_node = (md_array_size(ir->errors) == num_errors) ? prune_expressions(raw_node) : NULL;
        const ast_node_t* nested = pruned_node ? find_nested_assignment(pruned_node) : NULL;
        if (nested) {
            LOG_ERROR(ir, nested->token, "An assignment is a statement of its own, it cannot be a part of an expression");
            pruned_node = NULL;
        }
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

    md_temp_end(temp);
    return result;
}

static bool static_type_check(md_script_ir_t* ir, const md_system_t* sys, md_allocator_i* temp_arena) {
    ASSERT(ir);
    ASSERT(sys);

    bool result = true;

    eval_context_t ctx = {
        .ir = ir,
        .compile_ir = ir,
        .temp_alloc = temp_arena,
        .alloc = ir->arena,
        .sys = sys,
        // Static checking has no frame; the reference configuration is the only state available.
        .cur_state = &sys->reference,
        .ref_state = &sys->reference,
    };

    ir->stage = "Static Type Checking";

    for (size_t i = 0; i < md_array_size(ir->parsed_expressions); ++i) {
        ir->record_log = true;   // We reset the error recording flag before each statement
        // A check which logged an error and still succeeded leaves a node which is not what was written, or has no
        // type: it is not evaluated
        const size_t num_errors = md_array_size(ir->errors);
        if (static_check_node(ir->parsed_expressions[i], &ctx) && md_array_size(ir->errors) == num_errors) {
            md_array_push(ir->type_checked_expressions, ir->parsed_expressions[i], ir->arena);
            ir->flags |= (ir->parsed_expressions[i]->flags & FLAG_IR_PROPAGATION_MASK);
        } else {
            MD_LOG_DEBUG("Static type checking failed for expression: '"STR_FMT"'", STR_ARG(ir->parsed_expressions[i]->token.str));
            result = false;
        }
    }

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
    return ti.base_type == TYPE_FLOAT && ti.dim[0] != -1 && (ndim == 1 || (ndim == 2 && ti.dim[1] == 1));
}

// One or more distributions: [n][bins]
static inline bool is_distribution_type(type_info_t ti) {
    return ti.base_type == TYPE_DISTRIBUTION && ti.dim[0] > 0 && ti.dim[1] > 0 && ti.dim[2] == 0;
}

static inline bool is_volume_type(type_info_t ti) {
    return volume_count(ti) > 0;
}

// An array of volumes evaluated as one property: each is 8 MiB (MD_VOL_DIM^3 floats), held by the evaluation, its
// periods and every thread evaluating a frame. 32 are 256 MiB.
#define MAX_PROPERTY_VOLUMES 32

static inline bool is_property_type(type_info_t ti) {
    return is_temporal_type(ti) || is_distribution_type(ti) || is_volume_type(ti);
}

static bool static_eval_node(ast_node_t* node, eval_context_t* ctx) {
    ASSERT(node);
    ASSERT(ctx);
    ASSERT(ctx->compile_ir);

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
    // An absent argument has no data to evaluate (see PARAM_NULLABLE)
    if (!(node->flags & FLAG_DYNAMIC) && !node->data.ptr && !is_absent_type(node->data.type)) {
        uint64_t hash = hash_node(node, 0);

        // Try to find in existing expressions
        const size_t num_expressions = md_array_size(ctx->compile_ir->static_expression_hash);
        for (size_t i = 0; i < num_expressions; ++i) {
            if (hash == ctx->compile_ir->static_expression_hash[i]) {
                // str_t str = ctx->compile_ir->static_expression_str[i];
                node->data = ctx->compile_ir->static_expression_data[i];
                node->flags |= FLAG_CONSTANT;
                return true;
            }
        }

        data_t data = {0};
        if (allocate_data(&data, node->data.type, ctx->compile_ir->arena)) {
            if (!evaluate_node(&data, node, ctx)) {
                EVAL_LOG_ERROR(ctx, node->token, "Failed to evaluate node during static evaluation");
                return false;
            }
            if (node->type == AST_CONTEXT) {
                ASSERT(node->children);
                // If its a context node, we copy the data to the child as well
                node->children[0]->data = data;
            }
            node->data = data;
            node->flags |= FLAG_CONSTANT;

            md_array_push(ctx->compile_ir->static_expression_hash, hash, ctx->alloc);
            md_array_push(ctx->compile_ir->static_expression_data, data, ctx->alloc);
            md_array_push(ctx->compile_ir->static_expression_str, node->token.str, ctx->alloc);

            return true;
        } else {
            EVAL_LOG_ERROR(ctx, node->token, "Could not allocate data for node during static evaluation");
            return false;
        }
    }

    return true;
}

static bool static_evaluation(md_script_ir_t* ir, const md_system_t* sys) {
    ASSERT(sys);

    md_temp_scope_t temp_scope = md_temp_begin_avoid(ir->arena);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
    const size_t temp_reset_pos = md_vm_arena_get_pos(temp_alloc);

    eval_context_t ctx = {
        .ir = ir,
        .compile_ir = ir,
        .sys = sys,
        .temp_alloc = temp_alloc,
        .alloc = ir->arena,
        .cur_state = &sys->reference,
        .ref_state = &sys->reference,
    };

    ir->stage = "Static Evaluation";
    bool result = true;

    // Evaluate every node which is not flagged with FLAG_DYNAMIC and store its value.
    for (size_t i = 0; i < md_array_size(ir->type_checked_expressions); ++i) {
        result &= static_eval_node(ir->type_checked_expressions[i], &ctx);
        md_vm_arena_set_pos_back(temp_alloc, temp_reset_pos);
    }

    md_temp_end(temp_scope);

    return result;
}

// ### PROPERTY ATTRIBUTES ###
//
// Every property an evaluation computes is published into the evaluation's own attribute table
// under 'script/', and the table is the ONLY storage there is: eval_property_t holds writable views
// into it and nothing else of substance.
//
// The property kind picks the shape and nothing else needs to be said about it:
//
//   temporal      script/<ident>          {num_frames, population}   MD_ATTRIBUTE_FLAG_TEMPORAL
//                 script/<ident>/range    rank 0, 2 components (min, max) of the values
//                 script/<ident>/mean     {num_frames}               only when a population exists
//                 script/<ident>/variance {num_frames}
//                 script/<ident>/extent   {num_frames}, 2 components (min, max)
//                 script/<ident>/labels   {population}, str: the names of the members, when they are what an
//                                         attr() pattern matched ('e = attr("edr/coul_sr_*")')
//   distribution  script/<ident>          {num_bins}
//                 script/<ident>/weight   {num_bins}
//                 script/<ident>/range    rank 0, 2 components (min, max) of the bin axis
//                 script/<ident>/bin      {num_bins}, virtual: the x coordinate of each bin
//   volume        script/<ident>          {x, y, z}
//
// A distribution's bin coordinates are a SIBLING attribute rather than a pair of numbers hung off
// the values, which is the same choice the coordinate axes elsewhere in the table make. A consumer
// that can already plot one attribute against another then needs no special case for a histogram.
//
// 'range' is the range the script DECLARED for the property when it declared one, and what the
// evaluation observed otherwise. For a temporal property that is the range of its values, which is
// what a histogram of them is binned over; for a distribution it is the range its bins cover, in
// the unit of /bin. Either way it is the domain a consumer bins or plots the property over.

// A script identifier is one path segment. Nothing in the grammar produces a '/', but a stray one
// would silently turn the leaf into a group, so it is folded rather than trusted.
static str_t script_attr_path(char* buf, size_t cap, str_t ident, const char* suffix) {
    size_t len = 0;
#define SCRIPT_ATTR_PUT(c) do { if (len + 1 < cap) buf[len++] = (c); } while (0)
    for (const char* c = "script/"; *c; ++c) SCRIPT_ATTR_PUT(*c);
    for (size_t i = 0; i < ident.len; ++i)   SCRIPT_ATTR_PUT(ident.ptr[i] == '/' ? '_' : ident.ptr[i]);
    for (const char* c = suffix ? suffix : ""; *c; ++c) SCRIPT_ATTR_PUT(*c);
#undef SCRIPT_ATTR_PUT
    buf[len] = '\0';
    return (str_t){buf, len};
}

// Reserves zeroed storage inside the table and hands back the writable view. A table that refuses
// the attribute costs the published NAME and nothing else: the caller still gets a buffer, so a
// rejected path can never break an evaluation.
static float* script_attr_storage(md_attributes_t* attributes, md_allocator_i* alloc, str_t ident, const char* suffix, md_unit_t unit,
                                  md_attribute_flags_t flags, uint32_t components, uint32_t rank, const uint32_t shape[])
{
    char buf[512];
    md_attribute_desc_t desc = {
        .path   = script_attr_path(buf, sizeof(buf), ident, suffix),
        .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = components, .rank = rank },
        .flags  = flags,
        .unit   = unit,
    };
    for (uint32_t i = 0; i < rank; ++i) {
        desc.format.shape[i] = shape[i];
    }

    md_attribute_id_t id = md_attributes_replace(attributes, &desc);
    float* ptr = id ? (float*)md_attributes_data(attributes, id, MD_ATTRIBUTE_TYPE_F32) : NULL;
    if (!ptr) {
        MD_LOG_DEBUG("Script eval: could not publish attribute '"STR_FMT"', falling back to plain storage", STR_ARG(desc.path));
        const size_t bytes = md_attribute_element_count(&desc.format) * sizeof(float);
        ptr = md_alloc(alloc, bytes);
        MEMSET(ptr, 0, bytes);
    }
    return ptr;
}

// The bin coordinates of a distribution are its range cut into equal bins, and that range is not
// known until frames have been evaluated - it may be widened by the data itself. Computing them on
// demand keeps them correct without anyone having to remember to rewrite them.
static size_t script_bin_coord_provide(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const eval_property_t* prop = (const eval_property_t*)user_data;
    if (!prop || !prop->range || attr->format.rank != 1) {
        return 0;
    }

    const size_t num_bins = attr->format.shape[0];
    size_t beg   = 0;
    size_t count = num_bins;
    if (slice && slice->num_idx > 0) {
        beg   = slice->idx[0];
        count = (beg < num_bins) ? 1 : 0;
    }
    if (count > cap) {
        count = cap;
    }

    const double x_min = prop->range->x;
    const double x_max = prop->range->y;
    const double x_scl = num_bins ? (x_max - x_min) / (double)num_bins : 0.0;

    float* out = (float*)dst;
    for (size_t i = 0; i < count; ++i) {
        out[i] = (float)(x_min + ((double)(beg + i) + 0.5) * x_scl);
    }
    return count;
}

// Whether two patterns name their members alike
static bool same_member_labels(const ast_node_t* a, const ast_node_t* b) {
    const size_t n = md_array_size(a->attr_members);
    if (n != md_array_size(b->attr_members)) return false;
    for (size_t i = 0; i < n; ++i) {
        if (!str_eq(a->attr_members[i].label, b->attr_members[i].label)) return false;
    }
    return true;
}

// The attr() pattern whose members the outermost axis of what node evaluates to holds, if any: the pattern itself, a
// name it is assigned to, and arithmetic on it that keeps that axis ('e * 2', 'e - e0', '-e', 'e + f' where f names its
// members alike). NULL for anything else, which leaves the axis unnamed: an element of it, a sum over it, and arithmetic
// with an array as long that names its members otherwise or not at all.
static const ast_node_t* labelled_pattern(const ast_node_t* node, int depth) {
    if (!node || depth > 64) return NULL;
    const size_t num_children = md_array_size(node->children);
    const ast_node_t* src = NULL;
    switch (node->type) {
    case AST_ATTRIBUTE:
        src = node->attr_pattern ? node : NULL;
        break;
    case AST_ASSIGNMENT:
        // To one name; one of several ('{a, b} = ...') is an element of what is assigned
        if (num_children == 2 && node->children[0]->type == AST_IDENTIFIER) {
            src = labelled_pattern(node->children[1], depth + 1);
        }
        break;
    case AST_IDENTIFIER:
        // A reference is linked to its definition by the static check
        if (num_children == 1) {
            src = labelled_pattern(node->children[0], depth + 1);
        }
        break;
    case AST_PROC_CALL:
        // An arithmetic operator, which the static check made a call of the operator's procedure
        if (node->proc && num_children >= 1 && num_children <= 2 &&
            (str_eq(node->proc->name, STR_LIT("+")) || str_eq(node->proc->name, STR_LIT("-")) ||
             str_eq(node->proc->name, STR_LIT("*")) || str_eq(node->proc->name, STR_LIT("/")))) {
            for (size_t i = 0; i < num_children; ++i) {
                const ast_node_t* l = labelled_pattern(node->children[i], depth + 1);
                if (!l) {
                    // An array as long without names of its own, whose members these are as much as the other's
                    if (node->children[i]->data.type.dim[0] > 1 && node->children[i]->data.type.dim[0] == node->data.type.dim[0]) return NULL;
                    continue;
                }
                if (src && src != l && !same_member_labels(src, l)) return NULL;
                src = l;
            }
        }
        break;
    default:
        break;
    }
    // Only while the outermost axis is still the members'
    if (src && node->data.type.dim[0] != (int)md_array_size(src->attr_members)) return NULL;
    return src;
}

// script/<ident>/labels: the names of the members of a temporal property that holds what an attr() pattern matched,
// one value per member per frame (float[n]), in the order of its population.
static void publish_property_labels(md_attributes_t* attributes, str_t ident, const ast_node_t* node, md_allocator_i* temp_alloc) {
    const ast_node_t* src = labelled_pattern(node, 0);
    if (!src || dim_ndims(node->data.type.dim) != 1) return;
    const size_t n = md_array_size(src->attr_members);
    str_t* labels = md_alloc(temp_alloc, n * sizeof(str_t));
    for (size_t i = 0; i < n; ++i) {
        labels[i] = src->attr_members[i].label;
    }
    char buf[512];
    md_attribute_desc_t desc = {
        .path        = script_attr_path(buf, sizeof(buf), ident, "/labels"),
        .format      = { .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 1, .shape = { (uint32_t)n } },
        .description = STR_INIT("The names of the members, of what attr() matched"),
        .data        = labels,
        .byte_size   = n * sizeof(str_t),
    };
    md_attributes_replace(attributes, &desc);
}

// Publishes the attributes of one property and points prop at them.
// unit[0] is the coordinate unit (a distribution's bin axis), unit[1] the unit of the values.
// flags are those of the property (md_script_ir_property_flags): its kind, and whether it is periodic, in period.
static void init_property(eval_property_t* prop, md_attributes_t* attributes, str_t ident, md_script_property_flags_t flags,
                          type_info_t type, const md_unit_t unit[2], period_t period, size_t num_frames, md_allocator_i* alloc) {
    ASSERT(prop);
    ASSERT(attributes);
    ASSERT(alloc);

    // The kind alone: what is evaluated, and how it is kept
    const md_script_property_flags_t kind = flags & (MD_SCRIPT_PROPERTY_FLAG_TEMPORAL | MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION | MD_SCRIPT_PROPERTY_FLAG_VOLUME);

    MEMSET(prop, 0, sizeof(eval_property_t));
    prop->ident = ident;
    prop->kind  = kind;
    md_mutex_init(&prop->accum_mutex);

    // @NOTE: We need to 'normalize' the dimensionality of the types in a consistent way, since the properties will be exposed
    // Therefore we make sure all types encode the array length in the first dimension, even if its a scalar, i.e. [1]

    // Step 1: Prune trailing ones
    for (int i = MAX_NUM_DIMS - 1; i > 0; --i) {
        if (type.dim[i] == 1) {
            type.dim[i] = 0;
        }
    }

    const md_unit_t coord_unit = unit[0];
    const md_unit_t value_unit = unit[1];

    // @NOTE: At this stage, the flags should only contain the property type
    switch (kind) {
        case MD_SCRIPT_PROPERTY_FLAG_TEMPORAL: {
            // For temporal data, we store and expose all values, this enables filtering to be performed afterwards to create distributions
            if (dim_ndims(type.dim) < 2) {
                dim_shift_right(type.dim);
            }
            const uint32_t shape[2] = { (uint32_t)num_frames, (uint32_t)type.dim[1] };
            prop->count      = shape[1];
            prop->num_values = (size_t)shape[0] * (size_t)shape[1];
            prop->values = script_attr_storage(attributes, alloc, ident, NULL, value_unit, MD_ATTRIBUTE_FLAG_TEMPORAL, 1, 2, shape);
            prop->range  = (vec2_t*)script_attr_storage(attributes, alloc, ident, "/range", value_unit, MD_ATTRIBUTE_FLAG_NONE, 2, 0, NULL);

            if (shape[1] > 1) {
                // The population has more than one member, so the per frame summary over it is a
                // quantity in its own right and gets published as one.
                const uint32_t frames[1] = { shape[0] };
                prop->mean = script_attr_storage(attributes, alloc, ident, "/mean", value_unit, MD_ATTRIBUTE_FLAG_TEMPORAL, 1, 1, frames);
                prop->var  = script_attr_storage(attributes, alloc, ident, "/variance", md_unit_mul(value_unit, value_unit), MD_ATTRIBUTE_FLAG_TEMPORAL, 1, 1, frames);
                // min and max are two COMPONENTS of one value, not two positions along an axis:
                // they are not interchangeable and nothing indexes between them.
                prop->ext  = (vec2_t*)script_attr_storage(attributes, alloc, ident, "/extent", value_unit, MD_ATTRIBUTE_FLAG_TEMPORAL, 2, 1, frames);
            }

            if ((flags & MD_SCRIPT_PROPERTY_FLAG_PERIODIC) && period_valid(period)) {
                // The period its values are given in, and where a plot of them is best cut (see PERIODIC VALUES).
                // Known since the compilation, so written once.
                prop->period = period;
                float* p = script_attr_storage(attributes, alloc, ident, "/period", value_unit, MD_ATTRIBUTE_FLAG_NONE, 2, 0, NULL);
                p[0] = period.beg;
                p[1] = period.end;
                float* s = script_attr_storage(attributes, alloc, ident, "/seam", value_unit, MD_ATTRIBUTE_FLAG_NONE, 1, 0, NULL);
                s[0] = period.seam;
            }
            break;
        }
        case MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION: {
            // The evaluated type is [n][bins]: n distributions over the same bins, their weights beside them.
            // Published as {bins} for one, {n, bins} for more; the weights in the same shape.
            const uint32_t num = (uint32_t)type.dim[0];
            const uint32_t bins[1] = { (uint32_t)type.dim[1] };
            const uint32_t shape[2] = { num, bins[0] };
            const uint32_t rank = num > 1 ? 2 : 1;
            const uint32_t* dims = num > 1 ? shape : bins;
            prop->count      = bins[0];
            prop->num_values = (size_t)num * bins[0];
            prop->values  = script_attr_storage(attributes, alloc, ident, NULL,      value_unit,    MD_ATTRIBUTE_FLAG_NONE, 1, rank, dims);
            prop->weights = script_attr_storage(attributes, alloc, ident, "/weight", md_unit_none(), MD_ATTRIBUTE_FLAG_NONE, 1, rank, dims);
            prop->range   = (vec2_t*)script_attr_storage(attributes, alloc, ident, "/range", coord_unit, MD_ATTRIBUTE_FLAG_NONE, 2, 0, NULL);
            prop->sum_values  = md_alloc(alloc, prop->num_values * sizeof(double));
            prop->sum_weights = md_alloc(alloc, prop->num_values * sizeof(double));
            for (size_t i = 0; i < prop->num_values; ++i) {
                prop->weights[i] = 1.0f;
            }

            char path[512];
            md_attribute_virtual_t virt = {
                .provider  = script_bin_coord_provide,
                // Borrowed: the property lives in the same arena as the table and outlives it by
                // construction, so there is nothing here for the table to release.
                .user_data = prop,
                .user_data_size = 0,
            };
            md_attribute_desc_t desc = {
                .path   = script_attr_path(path, sizeof(path), ident, "/bin"),
                .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 1, .shape = { bins[0] } },
                .unit   = coord_unit,
                .virt   = &virt,
            };
            md_attributes_replace(attributes, &desc);
            break;
        }
        case MD_SCRIPT_PROPERTY_FLAG_VOLUME: {
            // A volume is {x, y, z}, an array of n of them {n, x, y, z}: one more index axis, which picks the volume
            const int num_volumes = volume_count(type);
            if (dim_ndims(type.dim) < 4) {
                dim_shift_right(type.dim);
            }
            const uint32_t shape[4] = { (uint32_t)num_volumes, (uint32_t)type.dim[1], (uint32_t)type.dim[2], (uint32_t)type.dim[3] };
            prop->num_values = (size_t)shape[0] * (size_t)shape[1] * (size_t)shape[2] * (size_t)shape[3];
            prop->count      = (uint32_t)prop->num_values;
            prop->values = num_volumes == 1 ?
                script_attr_storage(attributes, alloc, ident, NULL, value_unit, MD_ATTRIBUTE_FLAG_NONE, 1, 3, shape + 1) :
                script_attr_storage(attributes, alloc, ident, NULL, value_unit, MD_ATTRIBUTE_FLAG_NONE, 1, 4, shape);
            break;
        }
        default:
            ASSERT(false);
            break;
    }
}

// Two passes rather than a running sum of squares, which cancels badly in single precision.
static void compute_min_max_mean_variance(float* out_min, float* out_max, float* out_mean, float* out_var, const float* data, size_t count) {
    ASSERT(out_min && out_max && out_mean && out_var);
    ASSERT(data && count > 0);

    float min = FLT_MAX;
    float max = -FLT_MAX;
    float sum = 0;
    for (size_t i = 0; i < count; ++i) {
        sum += data[i];
        min = MIN(min, data[i]);
        max = MAX(max, data[i]);
    }
    const float mean = sum / (float)count;

    float var = 0;
    for (size_t i = 0; i < count; ++i) {
        var += (data[i] - mean) * (data[i] - mean);
    }

    *out_min  = min;
    *out_max  = max;
    *out_mean = mean;
    *out_var  = var / (float)count;
}

static int compare_float(const void* a, const void* b) {
    const float x = *(const float*)a;
    const float y = *(const float*)b;
    return (x > y) - (x < y);
}

// The summary of periodic values, in period p, as circular statistics have it (Mardia & Jupp, Directional Statistics,
// 2000): the values as angles on a circle, of which the period is the whole way around.
//   mean       the direction of the mean of the angles as unit vectors (the circular mean), in [p.beg, p.end[
//   var        the square of the circular standard deviation sqrt(-2 ln R), R the length of that mean vector, in the
//              unit of the values. It is the standard deviation of values which are close together (R near 1), and
//              grows without bound as they spread around the circle (R near 0): it is capped at half the period, for
//              mean +- that covers the whole period already.
//   min, max   the shortest arc which holds every value: from where the widest gap between two of them ends, around
//              to where it begins. min is in [p.beg, p.end[ and max - min is the length of the arc, so max may be past
//              p.end: the arc crosses the end of the period. Values which do not cross it give their least and greatest.
// scratch holds count floats.
static void compute_circular_summary(float* out_min, float* out_max, float* out_mean, float* out_var, const float* data, size_t count,
                                     period_t p, float* scratch) {
    ASSERT(out_min && out_max && out_mean && out_var);
    ASSERT(data && count > 0 && scratch);
    ASSERT(period_valid(p));

    const double beg = p.beg;
    const double len = period_length(p);
    const double to_angle = TWO_PI / len;

    double sum_cos = 0, sum_sin = 0;
    for (size_t i = 0; i < count; ++i) {
        const double v = period_wrap(data[i], beg, len);
        // Rounded to float, a value just short of the end may land on it: the same value, at the start
        scratch[i] = (float)v < p.end ? (float)v : p.beg;
        const double a = (v - beg) * to_angle;
        sum_cos += cos(a);
        sum_sin += sin(a);
    }
    sum_cos /= (double)count;
    sum_sin /= (double)count;

    const double R = sqrt(sum_cos * sum_cos + sum_sin * sum_sin);
    const double mean = period_wrap(beg + atan2(sum_sin, sum_cos) / to_angle, beg, len);
    // R is at most 1, up to rounding, and the deviation is that of all directions alike at most
    const double var_angle = R > 0.0 ? CLAMP(-2.0 * log(R), 0.0, PI * PI) : PI * PI;
    const double sd = sqrt(var_angle) / to_angle;

    // The widest gap between neighbours around the circle. The one across the end of the period is taken unless
    // another is wider, so that values which do not cross it come out as their least and greatest.
    qsort(scratch, count, sizeof(float), compare_float);
    size_t first = 0;
    double gap = (double)scratch[0] + len - (double)scratch[count - 1];
    for (size_t i = 1; i < count; ++i) {
        const double g = (double)scratch[i] - (double)scratch[i - 1];
        if (g > gap) {
            gap = g;
            first = i;
        }
    }

    *out_min  = scratch[first];
    *out_max  = (float)((double)scratch[first] + (len - gap));
    *out_mean = (float)mean;
    *out_var  = (float)(sd * sd);
}

static void fill_float(float* dst, size_t count, float value) {
    for (size_t i = 0; i < count; ++i) {
        dst[i] = value;
    }
}

static void clear_property(eval_property_t* prop) {
    // A temporal property holds a value per frame, and a frame which is not evaluated (yet) holds NaN, so
    // it cannot be taken for a value. Distributions and volumes are running means over the evaluated frames
    // (see eval_properties) and start from zero: NaN would poison them for good, as NaN * 0 is NaN.
    // NaN is stored here, but not tested for: the library is built with -ffast-math (/fp:fast), under which
    // isnan() may be folded away. done is what says which frames are evaluated.
    const float clear_value = (prop->kind == MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) ? NAN : 0.0f;
    fill_float(prop->values, prop->num_values, clear_value);
    if (prop->mean) {
        // Temporal only
        const size_t num_frames = prop->num_values / prop->count;
        fill_float(prop->mean, num_frames, NAN);
        fill_float(prop->var,  num_frames, NAN);
        for (size_t i = 0; i < num_frames; ++i) {
            prop->ext[i] = (vec2_t){NAN, NAN};
        }
    }
    if (prop->sum_values) {
        // Distribution: nothing accumulated. The weights go back to their 1.0 default, not 0, so values / weights is 0
        MEMSET(prop->sum_values,  0, prop->num_values * sizeof(double));
        MEMSET(prop->sum_weights, 0, prop->num_values * sizeof(double));
        fill_float(prop->weights, prop->num_values, 1.0f);
    }
    if (prop->range) {
        *prop->range = (vec2_t){0, 0};
    }
    prop->min_value = +FLT_MAX;
    prop->max_value = -FLT_MAX;
    prop->accum_count = 0;
    md_bitmask_clear(&prop->done);
}

// A version of a slot: never repeats, also across md_script_eval_create_from, and is never 0
static uint64_t next_version(md_script_eval_t* eval) {
    return md_atomic_fetch_add64(&eval->version_counter, 1) + 1;
}

static inline bool is_aggregate(const eval_property_t* prop) {
    return prop->kind == MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION || prop->kind == MD_SCRIPT_PROPERTY_FLAG_VOLUME;
}

// Whether the slot has frame f to add or to take out. Only a slot added up as frames are evaluated has: one of a
// distribution kept per frame needs nothing evaluated. The caller holds the slot's lock.
static inline bool slot_needs_frame(const eval_slot_t* s, uint32_t f) {
    return s->active && s->sum_values && md_bitmask_test_bit(&s->wanted, f) != md_bitmask_test_bit(&s->included, f);
}

// Adds frame f to the slot or takes it out, whichever it needs, from the values (and weights) of the frame. The
// caller holds the slot's lock.
static void slot_add_frame(md_script_eval_t* eval, eval_slot_t* s, uint32_t f, const float* values, const float* weights) {
    if (!slot_needs_frame(s, f)) return;
    const bool add = md_bitmask_test_bit(&s->wanted, f);
    if (add) {
        for (size_t i = 0; i < s->num_values; ++i) s->sum_values[i] += values[i];
        if (s->sum_weights) {
            for (size_t i = 0; i < s->num_values; ++i) s->sum_weights[i] += weights[i];
        }
        md_bitmask_set_bit(&s->included, f);
        s->count += 1;
    } else {
        s->count -= 1;
        md_bitmask_clear_bit(&s->included, f);
        if (s->count == 0) {
            // What is taken out is never quite what was added, in floating point: empty is zero, exactly
            MEMSET(s->sum_values, 0, s->num_values * sizeof(double));
            if (s->sum_weights) MEMSET(s->sum_weights, 0, s->num_values * sizeof(double));
        } else {
            for (size_t i = 0; i < s->num_values; ++i) s->sum_values[i] -= values[i];
            if (s->sum_weights) {
                for (size_t i = 0; i < s->num_values; ++i) s->sum_weights[i] -= weights[i];
            }
        }
    }
    s->version = next_version(eval);
}

// What frame f is still needed for: the properties which do not hold it, and those with a slot which has it to add
// or to take out. todo[p] for property p, may be NULL. Takes the eval's frame_lock, then each slot's lock in turn.
// Returns whether there is anything.
static bool frame_todo(md_script_eval_t* eval, uint32_t f, bool* todo) {
    const size_t num_props = md_array_size(eval->props);
    bool any = false;

    md_mutex_lock(&eval->frame_lock);
    for (size_t p = 0; p < num_props; ++p) {
        const bool t = !md_bitmask_test_bit(&eval->props[p].done, f);
        if (todo) todo[p] = t;
        any |= t;
    }
    for (uint32_t i = 0; i < eval->num_slots_used; ++i) {
        eval_slot_t* s = &eval->slots[i];
        if (!s->active || !s->sum_values) continue;
        md_mutex_lock(&s->lock);
        if (slot_needs_frame(s, f)) {
            if (todo) todo[s->prop] = true;
            any = true;
        }
        md_mutex_unlock(&s->lock);
    }
    md_mutex_unlock(&eval->frame_lock);
    return any;
}

static bool eval_properties(md_script_eval_t* eval, const md_system_t* sys, str_t run, const md_script_ir_t* ir, const uint32_t* frames, size_t num_frames) {
    ASSERT(eval);
    ASSERT(sys);
    ASSERT(ir);

    // No properties to evaluate!
    if (md_array_size(eval->props) == 0) return true;
    
    const size_t num_expr = md_array_size(ir->eval_targets);
    ast_node_t** const expr = ir->eval_targets;
    
    md_temp_scope_t temp = md_temp_begin_avoid(eval->arena);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool result = true;

    const size_t num_props = md_array_size(eval->props);
    ASSERT(num_expr == eval->num_targets);
    // Per property: what the frame is still needed for, whether this added it to its published aggregate, and what
    // it evaluated to
    bool* todo    = md_temp_alloc_array(temp, bool, num_props);
    bool* claimed = md_temp_alloc_array(temp, bool, num_props);
    const data_t** prop_data = md_temp_alloc_array(temp, const data_t*, num_props);
    // The targets those need evaluated
    uint8_t* needed = md_temp_alloc_array(temp, uint8_t, MAX(num_expr, 1));

    // What the procedures report while evaluating. Not the IR's (it is read-only here): reported below, with the frame
    eval_log_t log = { .alloc = temp_alloc };

    // Per worker states. Each evaluation range owns its coordinates, which is why the whole system
    // had to be copied before: it was the only way to get a private coordinate buffer.
    md_system_state_t cur_state = { .alloc = temp_alloc };
    md_system_state_init(&cur_state, sys->atom.count);

    // No per frame attribute table on this state. The scratch it would live in is rewound after
    // every frame, so a table filled on frame f would be read on frame f+1 out of memory the
    // evaluation has since reused. The evaluation asks for nothing that would land there anyway.
    cur_state.attributes.alloc = NULL;

    // coordinate data for reading trajectory frames into
    const size_t stride = ALIGN_TO(sys->atom.count, 16);    // Round up allocation size to simd width to allow for vectorized operations

    float* atom_mass   = md_temp_alloc_array(temp, float, stride);
    float* atom_radius = md_temp_alloc_array(temp, float, stride);

    // One extraction context for the range: the files a run streams from stay open across its
    // frames. It is made on this thread and used by it alone, which is what a context asks for.
    const str_t coord_paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };
    md_system_extract_t* extract = md_system_extract_begin(sys, run, coord_paths, ARRAY_SIZE(coord_paths), md_get_heap_allocator());
    if (!extract) {
        MD_LOG_ERROR("Failed to begin extracting frames for evaluation");
        result = false;
        goto done;
    }
    
    // This data is meant to hold the evaluated expressions
    data_t* data = md_temp_alloc_array(temp, data_t, num_expr);

    md_atom_extract_masses(atom_mass,   0, sys->atom.count, &sys->atom);
    md_atom_extract_radii (atom_radius, 0, sys->atom.count, &sys->atom);

    const size_t STACK_RESET_POINT = md_vm_arena_get_pos(temp_alloc);

    //int thread_id = (md_thread_id() & 0xFFFF);
    //md_logf(MD_LOG_TYPE_DEBUG, "Starting evaluation on thread %i, range (%i,%i) arena size: %.2f MB", thread_id, (int)frame_beg, (int)frame_end, (double)vm_arena.commit_pos / (double)MEGABYTES(1));
    //uint64_t max_arena_pos = 0;

    eval_context_t ctx = {
        .ir = ir,   // Read-only, shared with the evaluations of the other frame ranges
        .log = &log,
        .sys = sys,
        .atom_mass = atom_mass,
        .atom_radius = atom_radius,
        .temp_alloc = temp_alloc,
        .alloc = temp_alloc,
        .cur_state = &cur_state,
        .ref_state = &sys->reference,
    };
    {
        const md_attribute_t* axis = md_attributes_find_in(&sys->attributes, run, STR_LIT("time"));
        ctx.run_axis_id = axis ? axis->id : MD_ATTRIBUTE_INVALID;
    }

    // We evaluate each frame, one at a time
    for (size_t list_idx = 0; list_idx < num_frames; ++list_idx) {
        const uint32_t f_idx = frames[list_idx];
        if (eval->interrupt) {
            goto done;
        }

        // Nothing left to do for it: every property holds it, and no slot has it to add or take out
        const bool pending = f_idx < eval->frame_count && frame_todo(eval, f_idx, todo);
        if (!pending) {
            continue;
        }

        // What the properties still to do need evaluated, and nothing else
        MEMSET(needed, 0, num_expr);
        for (size_t p = 0; p < num_props; ++p) {
            if (!todo[p]) continue;
            const uint8_t* row = eval->needs + p * num_expr;
            for (size_t i = 0; i < num_expr; ++i) {
                needed[i] |= row[i];
            }
        }
        
        result = md_system_extract_frame(extract, f_idx, &cur_state);

        if (!result) {
            MD_LOG_ERROR("Failed to load frame during evaluation");
            goto done;
        }
        
        md_vm_arena_set_pos_back(temp_alloc, STACK_RESET_POINT);
        ctx.identifiers = NULL;
        log = (eval_log_t){ .alloc = temp_alloc };     // its arrays went with the reset
		MEMSET(&ctx.spatial_acc, 0, sizeof(ctx.spatial_acc));
        ctx.spatial_acc_cutoff = 0.0;

        for (size_t i = 0; i < num_expr; ++i) {
            if (!needed[i]) {
                continue;
            }
            if (!evaluate_node_alloc(&data[i], expr[i], &ctx, temp_alloc)) {
                str_t str = expr[i]->token.str;
                MD_LOG_ERROR("Evaluation error when evaluating the following expression '"STR_FMT"' at frame %i", STR_ARG(str), (int)f_idx);
                for (size_t e = 0; e < md_array_size(log.errors); ++e) {
                    MD_LOG_ERROR("    "STR_FMT, STR_ARG(log.errors[e].text));
                }
                result = false;
                goto done;
            }
            // Nothing reads these during the evaluation, the properties read them afterwards
            data[i].unit[0] = expr[i]->data.unit[0];
            data[i].unit[1] = expr[i]->data.unit[1];
            data[i].value_range = expr[i]->data.value_range;
            data[i].period = expr[i]->data.period;
            if (data[i].value_range.beg == 0 && data[i].value_range.end == 0) {
                data[i].value_range.beg = -FLT_MAX;
                data[i].value_range.end = +FLT_MAX;
            }
        }

        // What each property still to do evaluated to
        for (size_t p = 0; p < num_props; ++p) {
            prop_data[p] = NULL;
            if (!todo[p]) continue;
            eval_property_t* prop = &eval->props[p];
            ASSERT(str_eq(prop->ident, ir->property_names[p]));
            const identifier_t* ident = find_dynamic_identifier(prop->ident, &ctx);
            if (!ident || !ident->data) {
                MD_LOG_ERROR("Script eval: no evaluated data for property '"STR_FMT"'", STR_ARG(prop->ident));
                result = false;
                goto done;
            }
            prop_data[p] = ident->data;
            ASSERT(prop->values);
        }

        // In three steps. First what is kept per frame: the temporal values, and the distributions kept per frame.
        // Then, under the lock, the properties are marked as holding the frame, the published aggregates claim it,
        // and the slots of distributions kept per frame take it. Then the aggregates which claimed it, and the slots
        // added up as frames are evaluated, add it. In that order, a frame a property holds always has its values
        // written: done is the mask a reader on another thread takes them over.
        for (size_t p = 0; p < num_props; ++p) {
            if (!prop_data[p]) continue;
            eval_property_t* prop = &eval->props[p];
            const float* values = (const float*)prop_data[p]->ptr;
            const frange_t value_range = prop_data[p]->value_range;

            if (prop->kind == MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
                const uint32_t size = (uint32_t)prop_data[p]->type.dim[0];
                ASSERT(size == prop->count);
                MEMCPY(prop->values + (size_t)f_idx * size, values, size * sizeof(float));

                float min, max, mean, var;
                compute_min_max_mean_variance(&min, &max, &mean, &var, values, size);

                if (prop->mean) {
                    if (period_valid(prop->period)) {
                        // Around the circle the values lie on: the mean of 179 and -179 degrees is 180, not 0
                        float c_min, c_max, c_mean, c_var;
                        float* scratch = md_alloc(temp_alloc, size * sizeof(float));
                        compute_circular_summary(&c_min, &c_max, &c_mean, &c_var, values, size, prop->period, scratch);
                        prop->mean[f_idx] = c_mean;
                        prop->var [f_idx] = c_var;
                        prop->ext [f_idx] = (vec2_t){c_min, c_max};
                    } else {
                        prop->mean[f_idx] = mean;
                        prop->var [f_idx] = var;
                        prop->ext [f_idx] = (vec2_t){min, max};
                    }
                }

                // Frame ranges are evaluated concurrently, and the extremes are shared by all of them
                md_mutex_lock(&prop->accum_mutex);
                prop->min_value = MIN(prop->min_value, min);
                prop->max_value = MAX(prop->max_value, max);
                prop->range->x = (value_range.beg == -FLT_MAX) ? prop->min_value : value_range.beg;
                prop->range->y = (value_range.end == +FLT_MAX) ? prop->max_value : value_range.end;
                md_mutex_unlock(&prop->accum_mutex);
            } else if (prop->frame_values) {
                ASSERT(prop_data[p]->weight);
                MEMCPY(prop->frame_values  + (size_t)f_idx * prop->num_values, values,                prop->num_values * sizeof(float));
                MEMCPY(prop->frame_weights + (size_t)f_idx * prop->num_values, prop_data[p]->weight, prop->num_values * sizeof(float));
            }
        }

        // Claimed only now that it is evaluated, and under the lock: whatever else evaluates the same frame, each
        // aggregate adds it exactly once
        md_mutex_lock(&eval->frame_lock);
        for (size_t p = 0; p < num_props; ++p) {
            eval_property_t* prop = &eval->props[p];
            claimed[p] = prop_data[p] && !md_bitmask_test_bit(&prop->done, f_idx);
            if (claimed[p]) {
                md_bitmask_set_bit(&prop->done, f_idx);
            }
        }
        for (uint32_t i = 0; i < eval->num_slots_used; ++i) {
            eval_slot_t* s = &eval->slots[i];
            if (!s->active || s->sum_values || !claimed[s->prop]) continue;
            md_mutex_lock(&s->lock);
            if (md_bitmask_test_bit(&s->wanted, f_idx) && !md_bitmask_test_bit(&s->included, f_idx)) {
                md_bitmask_set_bit(&s->included, f_idx);
                s->count += 1;
                s->version = next_version(eval);
            }
            md_mutex_unlock(&s->lock);
        }
        md_mutex_unlock(&eval->frame_lock);

        for (size_t p = 0; p < num_props; ++p) {
            if (!claimed[p]) continue;
            eval_property_t* prop = &eval->props[p];
            const float* values = (const float*)prop_data[p]->ptr;
            const frange_t value_range = prop_data[p]->value_range;

            switch (prop->kind) {
            case MD_SCRIPT_PROPERTY_FLAG_TEMPORAL:
                break;
            case MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION: {
                // Every bin of every distribution in it
                const size_t num_bins = prop->num_values;
                ASSERT(prop_data[p]->type.base_type == TYPE_DISTRIBUTION);
                ASSERT(type_info_total_element_count(prop_data[p]->type) == num_bins);
                float min, max, mean, var;
                compute_min_max_mean_variance(&min, &max, &mean, &var, values, num_bins);

                const float* weights = prop_data[p]->weight;
                ASSERT(weights);

                md_mutex_lock(&prop->accum_mutex);
                {
                    // Both the values and their weights are summed, and what is published is the mean of each.
                    // Their ratio, which is what is shown (also after merging bins), is then the ratio of the sums:
                    // the estimate over the frames evaluated so far, whichever they are and in whatever order.
                    // Weights differ between frames (an rdf normalizes by the pair count of its frame), so taking
                    // those of any one frame is not the same.
                    const uint32_t count = ++prop->accum_count;
                    const double scl = 1.0 / (double)count;
                    for (size_t i = 0; i < num_bins; ++i) {
                        prop->sum_values [i] += values [i];
                        prop->sum_weights[i] += weights[i];
                        prop->values [i] = (float)(prop->sum_values [i] * scl);
                        prop->weights[i] = (float)(prop->sum_weights[i] * scl);
                    }

                    prop->min_value = MIN(prop->min_value, min);
                    prop->max_value = MAX(prop->max_value, max);
                    prop->range->x = (value_range.beg == -FLT_MAX) ? prop->min_value : value_range.beg;
                    prop->range->y = (value_range.end == +FLT_MAX) ? prop->max_value : value_range.end;
                }
                md_mutex_unlock(&prop->accum_mutex);
                break;
            }
            case MD_SCRIPT_PROPERTY_FLAG_VOLUME: {
                ASSERT(prop->num_values % 8 == 0); // This should always be the case if we nice powers of 2 for our volumes

                md_mutex_lock(&prop->accum_mutex);
                {
                    // Cumulative moving average
                    const uint32_t count = prop->accum_count++;
                    const md_256 N = md_mm256_set1_ps((float)(count));
                    const md_256 scl = md_mm256_set1_ps(1.0f / (float)(count + 1));

                    for (size_t i = 0; i < prop->num_values; i += 8) {
                        md_256 old_val = md_mm256_mul_ps(md_mm256_loadu_ps(prop->values + i), N);
                        md_256 new_val = md_mm256_loadu_ps(values + i);
                        md_mm256_storeu_ps(prop->values + i, md_mm256_mul_ps(md_mm256_add_ps(new_val, old_val), scl));
                    }
                }
                md_mutex_unlock(&prop->accum_mutex);
                break;
            }
            default:
                ASSERT(false);
                break;
            }
        }

        // The slots added up as frames are evaluated, of every property evaluated for the frame: whichever needed
        // it, a slot which has it to add or take out does so now. Read under nothing but the slot's own lock.
        md_mutex_lock(&eval->frame_lock);
        const uint32_t num_slots = eval->num_slots_used;
        md_mutex_unlock(&eval->frame_lock);
        for (uint32_t i = 0; i < num_slots; ++i) {
            eval_slot_t* s = &eval->slots[i];
            md_mutex_lock(&s->lock);
            if (s->active && s->sum_values && prop_data[s->prop]) {
                const data_t* d = prop_data[s->prop];
                slot_add_frame(eval, s, f_idx, (const float*)d->ptr, d->weight);
            }
            md_mutex_unlock(&s->lock);
        }
        
        //max_arena_pos = MAX(max_arena_pos, vm_arena.commit_pos);
    }
done:
    //md_logf(MD_LOG_TYPE_DEBUG, "Finished evaluation on thread %i, max arena size: %.2f MB", thread_id, (double)max_arena_pos / (double)MEGABYTES(1));
    md_temp_end(temp);
    md_system_extract_end(extract);
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
    ir->flags |= ctx_ir->flags;

    return true;
}

// A process wide counter: an IR freed and another allocated in its place must not share an id
static uint64_t next_ir_id(void) {
    static volatile int64_t counter = 0;
    return (uint64_t)atomic_add_i64(&counter, 1);
}

static md_script_vis_ref_t make_vis_ref(const md_script_ir_t* ir, const ast_node_t* node) {
    md_script_vis_ref_t ref = {0};
    // Only nodes of this IR: a node reached through an identifier of a context IR belongs to that one
    if (ir && ir->id && node && node->index < md_array_size(ir->nodes) && ir->nodes[node->index] == node) {
        ref.ir_id    = ir->id;
        ref.node_idx = node->index;
    }
    return ref;
}

static const ast_node_t* resolve_vis_ref(const md_script_ir_t* ir, md_script_vis_ref_t ref) {
    if (validate_ir(ir) && ref.ir_id != 0 && ref.ir_id == ir->id && ref.node_idx < md_array_size(ir->nodes)) {
        return ir->nodes[ref.node_idx];
    }
    return NULL;
}

static void create_vis_tokens(md_script_ir_t* ir, const ast_node_t* node, const ast_node_t* node_override, int32_t depth) {
    md_temp_scope_t temp = md_temp_begin_avoid(ir->arena);
    md_allocator_i* temp_arena = md_temp_allocator(temp);
    md_script_vis_token_t vis = {0};

    ASSERT(node->type != AST_UNDEFINED);

    // An argument the call did not write has nothing in the source to hover, and it carries the token of the call:
    // a marker of its own would take over the text of the call's
    if (node->flags & FLAG_IMPLICIT) {
        md_temp_end(temp);
        return;
    }

    md_str_writer_t sb = {0};
    md_str_writer_init(&sb, temp_arena);

    char type_buf[128];
    size_t type_len = print_type_info(type_buf, (int)sizeof(type_buf), node->data.type);
    md_str_writer_push_cstrl(&sb, type_buf, type_len);

    char unit_buf[128];
    int unit_len = (int)md_unit_print(unit_buf, (int)sizeof(unit_buf), node->data.unit[1]);
    if (unit_len) {
        md_str_writer_fmt(&sb, " ("STR_FMT")", unit_len, unit_buf);
    }

    if (node->type == AST_ATTRIBUTE) {
        if (node->attr_pattern) {
            // What the pattern matched, in the order of the array
            const size_t num = md_array_size(node->attr_members);
            const size_t shown = MIN(num, 16);
            md_str_writer_fmt(&sb, "\n"STR_FMT": %zu %s", STR_ARG(node->attr_path), num, num == 1 ? "match" : "matches");
            for (size_t i = 0; i < shown; ++i) {
                const attr_member_t* m = &node->attr_members[i];
                md_str_writer_fmt(&sb, "\n  [%zu] "STR_FMT, i + 1, STR_ARG(m->path));
                if (!str_eq(m->label, path_leaf(m->path))) {
                    md_str_writer_fmt(&sb, "  "STR_FMT, STR_ARG(m->label));
                }
            }
            if (shown < num) {
                md_str_writer_fmt(&sb, "\n  ... and %zu more", num - shown);
            }
        } else {
            md_str_writer_fmt(&sb, "\n"STR_FMT, STR_ARG(node->attr_path));
        }
    } else if (!(node->flags & FLAG_DYNAMIC)) {
        // A name, or what is assigned to it, holds the value of its expression
        const ast_node_t* vnode = node;
        if (!vnode->data.ptr && vnode->type == AST_IDENTIFIER && md_array_size(vnode->children) > 0) vnode = vnode->children[0];
        if (!vnode->data.ptr && vnode->type == AST_ASSIGNMENT && md_array_size(vnode->children) > 1) vnode = vnode->children[1];
        if (vnode->data.ptr) {
            md_str_writer_t val = md_str_writer_create(temp_arena);
            summarize_data(&val.w, vnode->data, -1);
            if (md_str_writer_len(&val)) {
                md_str_writer_push_char(&sb, '\n');
                md_str_writer_push_str(&sb, md_str_writer_to_str(&val));
            }
        }
    } else {
        md_str_writer_push_str(&sb, STR_LIT(" [d]"));
    }

    vis.range.beg = node->token.beg;
    vis.range.end = node->token.end;
    vis.depth = depth;
    vis.text = str_copy(md_str_writer_to_str(&sb), ir->arena);
    vis.ref = make_vis_ref(ir, node_override ? node_override : node);

    // Parent marker should be last marker added (unless empty)
    md_script_vis_token_t* last = md_array_last(ir->vis_tokens);
    if (last &&
        vis.range.beg == last->range.beg &&
        vis.range.end == last->range.end)
    {
        // If the current marker has the exact same line, col_beg and col_end, we only want to update the text, not the reference (the node)
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
    md_temp_end(temp);
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
        md_array_push(ir->identifier_names, ir->identifiers[i].name, ir->arena);
    }
    return true;
}

static bool create_property(md_script_ir_t* ir, str_t ident, const ast_node_t* node) {
    const type_info_t type = node->data.type;
    md_script_property_flags_t flags =
        is_temporal_type(type)     ? MD_SCRIPT_PROPERTY_FLAG_TEMPORAL :
        is_distribution_type(type) ? MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION :
        is_volume_type(type)       ? MD_SCRIPT_PROPERTY_FLAG_VOLUME :
                                     MD_SCRIPT_PROPERTY_FLAG_NONE;
    ASSERT(flags != MD_SCRIPT_PROPERTY_FLAG_NONE);
    // Only values over time are summarized and shown as periodic
    if (flags == MD_SCRIPT_PROPERTY_FLAG_TEMPORAL && period_valid(node->data.period)) {
        flags |= MD_SCRIPT_PROPERTY_FLAG_PERIODIC;
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
            const int num_volumes = volume_count(expr->data.type);
            if (num_volumes > MAX_PROPERTY_VOLUMES) {
                LOG_ERROR(ir, expr->token, "An array of %i volumes is %i MiB to evaluate, no more than %i (%i MiB) are evaluated as one property",
                    num_volumes, num_volumes * 8, MAX_PROPERTY_VOLUMES, MAX_PROPERTY_VOLUMES * 8);
                return false;
            }
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

bool md_script_ir_compile_from_source(md_script_ir_t* ir, str_t src, const md_system_t* sys, const md_script_ir_t* ctx_ir) {
    if (!validate_ir(ir)) {
        MD_LOG_ERROR("Script Compile: IR is not valid");
        return false;
    }

    if (str_empty(src)) {
        MD_LOG_ERROR("Script Compile: Source string was empty");
        return false;
    }

    if (!sys) {
        MD_LOG_ERROR("Script Compile: Molecule was not supplied");
        return false;
    }
    
    ir->str = str_copy(src, ir->arena);
    ir->id  = next_ir_id();

    if (ctx_ir) {
        add_ir_ctx(ir, ctx_ir);
    }

    md_temp_scope_t temp_scope = md_temp_begin_avoid(ir->arena);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);

    ir->compile_success =
        parse_script(ir, temp_alloc) &&
        static_type_check(ir, sys, temp_alloc) &&
        // optimize ast?
        //static_evaluation(ir, sys) &&
        extract_dynamic_evaluation_targets(ir) &&
        extract_identifiers(ir) &&
        extract_properties(ir);

    extract_vis_tokens(ir);

    // What each property is computed from. The topology seeds them: procedures read it without it appearing in
    // their arguments (residue(), the backbone of ramachandran()), so a property over another topology is another.
    if (ir->compile_success) {
        uint64_t* memo = md_temp_alloc_array(temp_scope, uint64_t, MAX(md_array_size(ir->nodes), 1));
        MEMSET(memo, 0, md_array_size(ir->nodes) * sizeof(uint64_t));
        const uint64_t seed = md_hash64(&sys->topology_version, sizeof(sys->topology_version), sys->atom.count);
        for (size_t i = 0; i < md_array_size(ir->property_names); ++i) {
            const identifier_t* decl = find_ir_identifier(ir, ir->property_names[i]);
            const uint64_t fingerprint = identifier_fingerprint(ir, decl, seed ^ (uint64_t)ir->property_flags[i], memo);
            md_array_push(ir->property_fingerprints, fingerprint, ir->arena);
        }
    }

    md_temp_end(temp_scope);

    if (ir->compile_success) {
        extract_vis_tokens(ir);

        uint64_t hash = 0;
        const size_t num_expr = md_array_size(ir->type_checked_expressions);
        for (size_t i = 0; i < num_expr; ++i) {
            hash = hash_node(ir->type_checked_expressions[i], hash);
        }
        ir->fingerprint = hash;
    }

    return ir->compile_success;
}

bool md_script_ir_add_identifier_bitset(md_script_ir_t* ir, str_t name, md_bitset_t set) {
    if (!validate_ir(ir)) {
        MD_LOG_ERROR("IR is not valid");
        return false;
    }

    if (!md_script_identifier_name_valid(name)) {
        MD_LOG_ERROR("'"STR_FMT"' is not a valid identifier", STR_ARG(name));
        return false;
    }


    identifier_t* ident = create_identifier(ir, name);
    if (!ident) {
        MD_LOG_ERROR("Failed to create identifier with name '"STR_FMT"', it is perhaps already taken.", STR_ARG(name));
        return false;
    }

    ast_node_t* node = create_node(ir, AST_CONSTANT_VALUE, (token_t){0});
    node->data.ptr = &node->value._bitset;
    node->data.size = sizeof(md_bitset_t);
    node->data.type = (type_info_t)TI_ATOMSET;
    node->data.alloc = ir->arena;
    node->value._bitset = md_bitset_copy(set, ir->arena);

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
        md_allocator_i* alloc = md_arena_allocator_backing(ir->arena);
        md_arena_allocator_destroy(ir->arena);
        MEMSET(ir, 0, sizeof(md_script_ir_t));
        md_free(alloc, ir, sizeof(md_script_ir_t));
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

// ### DIAGNOSTICS ###
// The errors and warnings of a compilation as text, in the format of GCC and clang.

#define DIAG_SOURCE_NAME "<script>"  // Where a file name would be (as <stdin> is for GCC and clang)
#define DIAG_TAB_WIDTH  4   // A tab is echoed as this many spaces, in the source line and the marker line alike
#define DIAG_EDGE_LINES 2   // A range over more lines than twice this shows only its first and last ones

typedef struct diag_entry_t {
    const md_log_token_t* tok;
    bool is_error;
} diag_entry_t;

// The trailing bytes of a UTF-8 sequence: they do not start a character
static inline bool diag_utf8_cont(char c) {
    return ((unsigned char)c & 0xC0) == 0x80;
}

// The columns the byte c takes when echoed
static inline int diag_char_width(char c) {
    return c == '\t' ? DIAG_TAB_WIDTH : diag_utf8_cont(c) ? 0 : 1;
}

static int diag_num_digits(int x) {
    int n = 1;
    while (x >= 10) {
        x /= 10;
        n += 1;
    }
    return n;
}

// One diagnostic as GCC and clang print theirs: its location and message, then each line of the source its range
// covers, with the range marked beneath it ('^' where it begins, '~' for the rest of it):
//
// <script>:2:5: error: Undefined identifier 'y'
//     2 | b = y + 2;
//       |     ^
static void diag_append(md_writer_i* w, str_t src, const md_log_token_t* tok, bool is_error) {
    const char* level = is_error ? "error" : "warning";
    const str_t text  = str_trim_end(tok->text);
    const int   len   = (int)MIN(src.len, (size_t)INT32_MAX);

    if (len == 0) {
        md_writef(w, "%s: %s: " STR_FMT "\n", DIAG_SOURCE_NAME, level, STR_ARG(text));
        return;
    }

    // A range from elsewhere (a token of a context IR) is clamped to the source rather than trusted
    const int beg = CLAMP(tok->range.beg, 0, len);
    const int end = CLAMP(tok->range.end, beg, len);

    // Line and column of beg, 1-based. The column counts characters, not bytes
    int line = 1;
    int line_beg = 0;
    for (int i = 0; i < beg; ++i) {
        if (src.ptr[i] == '\n') {
            line += 1;
            line_beg = i + 1;
        }
    }
    int col = 1;
    for (int i = line_beg; i < beg; ++i) {
        col += !diag_utf8_cont(src.ptr[i]);
    }

    md_writef(w, "%s:%d:%d: %s: " STR_FMT "\n", DIAG_SOURCE_NAME, line, col, level, STR_ARG(text));

    // The lines covered: through the one holding the last character of the range, the newline ending it excluded
    int num_lines = 1;
    for (int i = beg; i < end - 1; ++i) {
        num_lines += src.ptr[i] == '\n';
    }
    const int  gutter = MAX(5, diag_num_digits(line + num_lines - 1));
    const bool elide  = num_lines > 2 * DIAG_EDGE_LINES + 1;

    int cur = line_beg;
    for (int k = 0; k < num_lines; ++k) {
        int eol = cur;
        while (eol < len && src.ptr[eol] != '\n') ++eol;
        const int next = eol < len ? eol + 1 : len;
        if (eol > cur && src.ptr[eol - 1] == '\r') --eol;

        if (elide && DIAG_EDGE_LINES <= k && k < num_lines - DIAG_EDGE_LINES) {
            if (k == DIAG_EDGE_LINES) {
                md_writef(w, "%*s | ...\n", gutter, "");
            }
            cur = next;
            continue;
        }

        md_writef(w, "%*d | ", gutter, line + k);
        for (int i = cur; i < eol; ++i) {
            if (src.ptr[i] == '\t') {
                md_write_fill(w, ' ', DIAG_TAB_WIDTH);
            } else {
                md_write_char(w, src.ptr[i]);
            }
        }
        md_write_char(w, '\n');

        // What of the line the range marks [mb, me[: on the first line from beg, on the following ones from the
        // first character which is not blank, so that indentation is not marked
        int mb = MIN(MAX(beg, cur), eol);
        int me = MAX(mb, MIN(end, eol));
        if (k > 0) {
            while (mb < me && (src.ptr[mb] == ' ' || src.ptr[mb] == '\t')) ++mb;
        }

        // The first line always gets a marker, if only a '^' where the range begins (it may be empty)
        if (k == 0 || mb < me) {
            md_writef(w, "%*s | ", gutter, "");
            for (int i = cur; i < mb; ++i) {
                md_write_fill(w, ' ', (size_t)diag_char_width(src.ptr[i]));
            }
            if (mb == me) {
                md_write_char(w, '^');
            } else {
                for (int i = mb; i < me; ++i) {
                    const int cw = diag_char_width(src.ptr[i]);
                    if (cw > 0) {
                        md_write_char(w, i == beg ? '^' : '~');
                        md_write_fill(w, '~', (size_t)(cw - 1));
                    }
                }
            }
            md_write_char(w, '\n');
        }

        cur = next;
    }
}

void md_script_ir_write_diagnostics(md_writer_i* w, const md_script_ir_t* ir) {
    ASSERT(w);
    if (!validate_ir(ir)) {
        MD_LOG_ERROR("IR is not valid");
        return;
    }

    const size_t num_errors   = md_array_size(ir->errors);
    const size_t num_warnings = md_array_size(ir->warnings);
    const size_t num_entries  = num_errors + num_warnings;

    // A compilation which got under way (it has an id) and failed without saying why: that is worth knowing too
    const bool silent_failure = ir->id != 0 && !ir->compile_success && num_errors == 0;

    if (num_entries == 0 && !silent_failure) {
        return;
    }

    // The writer may grow in a temp arena: the scope must not be opened in the same one
    md_temp_scope_t temp = md_temp_begin_avoid(md_writer_allocator(w));

    // In the order they appear in the source, an error before a warning at the same place
    diag_entry_t* entries = md_temp_alloc_array(temp, diag_entry_t, MAX(num_entries, 1));
    size_t n = 0;
    for (size_t i = 0; i < num_errors;   ++i) entries[n++] = (diag_entry_t){ &ir->errors[i],   true  };
    for (size_t i = 0; i < num_warnings; ++i) entries[n++] = (diag_entry_t){ &ir->warnings[i], false };
    for (size_t i = 1; i < n; ++i) {
        const diag_entry_t e = entries[i];
        size_t j = i;
        while (j > 0 && (entries[j-1].tok->range.beg > e.tok->range.beg ||
                        (entries[j-1].tok->range.beg == e.tok->range.beg && !entries[j-1].is_error && e.is_error))) {
            entries[j] = entries[j-1];
            --j;
        }
        entries[j] = e;
    }

    for (size_t i = 0; i < n; ++i) {
        diag_append(w, ir->str, entries[i].tok, entries[i].is_error);
    }

    if (silent_failure) {
        md_write_cstr(w, DIAG_SOURCE_NAME ": error: compilation failed without reporting an error\n");
    }

    // As clang sums up
    if (num_warnings && num_errors) {
        md_writef(w, "%d warning%s and %d error%s generated.\n", (int)num_warnings, num_warnings == 1 ? "" : "s", (int)num_errors, num_errors == 1 ? "" : "s");
    } else if (num_warnings) {
        md_writef(w, "%d warning%s generated.\n", (int)num_warnings, num_warnings == 1 ? "" : "s");
    } else if (num_errors) {
        md_writef(w, "%d error%s generated.\n", (int)num_errors, num_errors == 1 ? "" : "s");
    }

    md_temp_end(temp);
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

// Index of the named property, -1 if there is none
static int ir_property_index(const md_script_ir_t* ir, str_t name) {
    if (validate_ir(ir)) {
        for (size_t i = 0; i < md_array_size(ir->property_names); ++i) {
            if (str_eq(ir->property_names[i], name)) {
                return (int)i;
            }
        }
    }
    return -1;
}

md_script_property_flags_t md_script_ir_property_flags(const md_script_ir_t* ir, str_t name) {
    const int idx = ir_property_index(ir, name);
    return idx < 0 ? MD_SCRIPT_PROPERTY_FLAG_NONE : ir->property_flags[idx];
}

uint64_t md_script_ir_property_fingerprint(const md_script_ir_t* ir, str_t name) {
    const int idx = ir_property_index(ir, name);
    return (idx < 0 || (size_t)idx >= md_array_size(ir->property_fingerprints)) ? 0 : ir->property_fingerprints[idx];
}

md_script_range_marker_t md_script_ir_property_range(const md_script_ir_t* ir, str_t name) {
    const int idx = ir_property_index(ir, name);
    if (idx < 0) return (md_script_range_marker_t){0};
    // The assignment, or the name within a list of them: its name either way
    const ast_node_t* node = ir->property_nodes[idx];
    if (node->type == AST_ASSIGNMENT && md_array_size(node->children) > 0) {
        node = node->children[0];
    }
    return (md_script_range_marker_t){ node->token.beg, node->token.end };
}

md_script_vis_ref_t md_script_ir_property_vis_ref(const md_script_ir_t* ir, str_t name) {
    const int idx = ir_property_index(ir, name);
    return idx < 0 ? (md_script_vis_ref_t){0} : make_vis_ref(ir, ir->property_nodes[idx]);
}

bool md_script_vis_ref_valid(const md_script_ir_t* ir, md_script_vis_ref_t ref) {
    return resolve_vis_ref(ir, ref) != NULL;
}

str_t md_script_vis_ref_ident(const md_script_ir_t* ir, md_script_vis_ref_t ref) {
    const ast_node_t* node = resolve_vis_ref(ir, ref);
    return node ? node->ident : (str_t){0};
}

int md_script_vis_ref_dim(const md_script_ir_t* ir, md_script_vis_ref_t ref) {
    const ast_node_t* node = resolve_vis_ref(ir, ref);
    return node ? node->data.type.dim[0] : 0;
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
    const md_script_eval_desc_t desc = { .num_frames = num_frames };
    return md_script_eval_create_desc(ir, &desc, alloc);
}

// Takes slot s into use for property p under key, holding nothing and asked for nothing. The caller holds the eval's
// frame_lock and the slot's lock, or is the only one with the eval.
static void slot_take(md_script_eval_t* eval, eval_slot_t* s, uint64_t key, uint32_t p) {
    const eval_property_t* prop = &eval->props[p];
    if (!s->wanted.alloc) {
        // Over every frame up front, so that marking a frame on an evaluating thread never allocates
        md_bitmask_init(&s->wanted,   0, (uint32_t)eval->frame_count, eval->arena);
        md_bitmask_init(&s->included, 0, (uint32_t)eval->frame_count, eval->arena);
        md_bitmask_reserve(&s->wanted,   0, (uint32_t)eval->frame_count);
        md_bitmask_reserve(&s->included, 0, (uint32_t)eval->frame_count);
    }
    md_bitmask_clear(&s->wanted);
    md_bitmask_clear(&s->included);
    s->active = true;
    s->key    = key;
    s->prop   = p;
    s->count  = 0;
    s->num_values  = prop->num_values;
    s->sum_values  = NULL;
    s->sum_weights = NULL;
    if (!prop->frame_values) {
        md_allocator_i* heap = md_get_heap_allocator();
        s->sum_values = md_alloc(heap, s->num_values * sizeof(double));
        MEMSET(s->sum_values, 0, s->num_values * sizeof(double));
        if (prop->kind == MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) {
            s->sum_weights = md_alloc(heap, s->num_values * sizeof(double));
            MEMSET(s->sum_weights, 0, s->num_values * sizeof(double));
        }
    }
    s->version = next_version(eval);
}

// Lets go of what slot s holds and takes it out of use. The caller holds its lock (and the eval's frame_lock, unless
// it is the only one with the eval).
static void slot_let_go(eval_slot_t* s) {
    md_allocator_i* heap = md_get_heap_allocator();
    if (s->sum_values)  md_free(heap, s->sum_values,  s->num_values * sizeof(double));
    if (s->sum_weights) md_free(heap, s->sum_weights, s->num_values * sizeof(double));
    s->sum_values  = NULL;
    s->sum_weights = NULL;
    if (s->wanted.alloc) {
        md_bitmask_clear(&s->wanted);
        md_bitmask_clear(&s->included);
    }
    s->active = false;
    s->key    = 0;
    s->prop   = 0;
    s->count  = 0;
}

static int find_property(const md_script_eval_t* eval, str_t ident) {
    for (size_t p = 0; p < md_array_size(eval->props); ++p) {
        if (str_eq(eval->props[p].ident, ident)) return (int)p;
    }
    return -1;
}

// The slot under key for property p, NULL if there is none. The caller holds frame_lock.
static eval_slot_t* find_slot(const md_script_eval_t* eval, uint64_t key, uint32_t p) {
    for (uint32_t i = 0; i < eval->num_slots_used; ++i) {
        eval_slot_t* s = &eval->slots[i];
        if (s->active && s->key == key && s->prop == p) return s;
    }
    return NULL;
}

// The evaluation target which assigns the variable, -1 if none does (a constant, or a name the script does not assign)
static int target_assigning(const md_script_ir_t* ir, str_t name) {
    for (size_t t = 0; t < md_array_size(ir->eval_targets); ++t) {
        const ast_node_t* lhs = ir->eval_targets[t]->children[0];
        if (lhs->type == AST_IDENTIFIER) {
            if (str_eq(lhs->ident, name)) return (int)t;
        } else if (lhs->type == AST_ARRAY) {
            for (size_t j = 0; j < md_array_size(lhs->children); ++j) {
                if (lhs->children[j]->type == AST_IDENTIFIER && str_eq(lhs->children[j]->ident, name)) return (int)t;
            }
        }
    }
    return -1;
}

// Adds to the row of target t the targets of the variables node refers to, and theirs in turn. A reference is not
// descended into: it stands for the target which assigns it, whose row is complete already, as a variable can only
// refer to the ones declared before it.
static void mark_references(uint8_t* row, size_t t, const ast_node_t* node, const md_script_ir_t* ir, const uint8_t* target_needs, size_t num_targets) {
    if (node->type == AST_IDENTIFIER) {
        const int u = target_assigning(ir, node->ident);
        if (u >= 0 && (size_t)u < t) {
            const uint8_t* dep = target_needs + (size_t)u * num_targets;
            for (size_t i = 0; i < num_targets; ++i) {
                row[i] |= dep[i];
            }
        }
        return;
    }
    for (size_t i = 0; i < md_array_size(node->children); ++i) {
        mark_references(row, t, node->children[i], ir, target_needs, num_targets);
    }
}

// Which targets each property needs evaluated: eval->needs
static void compute_needs(md_script_eval_t* eval, const md_script_ir_t* ir) {
    const size_t num_targets = md_array_size(ir->eval_targets);
    const size_t num_props   = md_array_size(ir->property_names);

    md_temp_scope_t temp = md_temp_begin_avoid(eval->arena);
    uint8_t* target_needs = md_temp_alloc_array(temp, uint8_t, MAX(num_targets * num_targets, 1));
    MEMSET(target_needs, 0, num_targets * num_targets);
    for (size_t t = 0; t < num_targets; ++t) {
        uint8_t* row = target_needs + t * num_targets;
        row[t] = 1;
        const ast_node_t* expr = ir->eval_targets[t];
        if (md_array_size(expr->children) == 2) {
            mark_references(row, t, expr->children[1], ir, target_needs, num_targets);
        }
    }

    eval->num_targets = num_targets;
    eval->needs = md_alloc(eval->arena, MAX(num_props * num_targets, 1));
    for (size_t p = 0; p < num_props; ++p) {
        uint8_t* row = eval->needs + p * num_targets;
        const int t = target_assigning(ir, ir->property_names[p]);
        if (t >= 0) {
            MEMCPY(row, target_needs + (size_t)t * num_targets, num_targets);
        } else {
            // Cannot happen for a property of the IR; if it did, everything is evaluated as before
            MEMSET(row, 1, num_targets);
        }
    }
    md_temp_end(temp);
}

// The property of prev computed from the same thing (md_script_ir_property_fingerprint), the first of them: what a
// property of that fingerprint takes over (md_script_eval_create_from). -1 if there is none, and over another number
// of frames, which is another run.
static int taken_over_from(const md_script_eval_t* prev, uint64_t fingerprint, size_t num_frames) {
    if (!prev || !fingerprint || prev->frame_count != num_frames) return -1;
    for (size_t j = 0; j < md_array_size(prev->props); ++j) {
        if (prev->props[j].fingerprint == fingerprint) return (int)j;
    }
    return -1;
}

// An evaluation for ir. Of what prev holds nothing is taken over here, but a distribution which will be is kept the way
// it is kept in prev, per frame or not, so that it comes along as it is (see md_script_eval_create_from).
static md_script_eval_t* eval_create(const md_script_ir_t* ir, const md_script_eval_desc_t* desc, md_allocator_i* alloc, const md_script_eval_t* prev) {
    if (!desc) {
        MD_LOG_ERROR("Script eval: Descriptor was null");
        return NULL;
    }
    const size_t num_frames = desc->num_frames;
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
    compute_needs(eval, ir);

    md_attributes_t* attributes = &eval->attributes;
    attributes->alloc = eval->arena;

    // The frame axis every temporal property is checked against. It lives at the ROOT, outside
    // 'script/', so that a property the user happens to call "time" can never become the axis of
    // its neighbours. The evaluation knows how many frames it covers but not when they were taken,
    // so the coordinates are ordinals and the unit is empty, which is the convention for exactly
    // that. The time of a frame is the trajectory's business.
    if (num_frames > 0) {
        md_attribute_desc_t axis = {
            .path   = STR_INIT("time"),
            .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)num_frames } },
            .flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit   = md_unit_none(),
            .label  = STR_INIT("Frame"),
        };
        md_attribute_id_t axis_id = md_attributes_create(attributes, &axis);
        double* ordinals = axis_id ? (double*)md_attributes_data(attributes, axis_id, MD_ATTRIBUTE_TYPE_F64) : NULL;
        if (ordinals) {
            for (size_t i = 0; i < num_frames; ++i) {
                ordinals[i] = (double)i;
            }
        }
    }

    // Sized up front rather than pushed one at a time: a published attribute may hold the address
    // of the property it describes, and a growing array would move it out from under one.
    const size_t num_props = md_array_size(ir->property_names);
    md_array_resize(eval->props, num_props, eval->arena);
    for (size_t i = 0; i < num_props; ++i) {
        eval_property_t* prop = &eval->props[i];
        const ast_node_t* node = ir->property_nodes[i];
        init_property(prop, attributes, ir->property_names[i], ir->property_flags[i], node->data.type, node->data.unit, node->data.period, num_frames, eval->arena);
        if (prop->kind == MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
            publish_property_labels(attributes, ir->property_names[i], node, eval->arena);
        }
        prop->fingerprint = i < md_array_size(ir->property_fingerprints) ? ir->property_fingerprints[i] : 0;
        // Reserved up front, so marking a frame never moves the words while another thread reads them
        md_bitmask_init(&prop->done, 0, (uint32_t)num_frames, eval->arena);
        md_bitmask_reserve(&prop->done, 0, (uint32_t)num_frames);
    }

    // Distributions are kept per frame while there is room for them in the budget, in script order. One taken over
    // from prev is kept the way it was kept there, and has its room first.
    const size_t budget = desc->frame_budget ? desc->frame_budget : MD_SCRIPT_EVAL_FRAME_BUDGET_DEFAULT;
    size_t kept = 0;
    for (int pass = 0; pass < 2; ++pass) {
        for (size_t i = 0; i < num_props; ++i) {
            eval_property_t* prop = &eval->props[i];
            if (prop->kind != MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) continue;
            const int j = taken_over_from(prev, prop->fingerprint, num_frames);
            if ((pass == 0) != (j >= 0)) continue;
            const size_t bytes = num_frames * prop->num_values * sizeof(float);
            const bool per_frame = j >= 0 ? prev->props[j].frame_values != NULL : kept <= budget && 2 * bytes <= budget - kept;
            if (per_frame) {
                prop->frame_values  = md_alloc(eval->arena, bytes);
                prop->frame_weights = md_alloc(eval->arena, bytes);
                kept += 2 * bytes;
            }
        }
    }

    for (size_t i = 0; i < num_props; ++i) {
        clear_property(&eval->props[i]);
    }

    eval->slots = md_alloc(eval->arena, MD_SCRIPT_EVAL_MAX_SLOTS * sizeof(eval_slot_t));
    MEMSET(eval->slots, 0, MD_SCRIPT_EVAL_MAX_SLOTS * sizeof(eval_slot_t));
    for (size_t i = 0; i < MD_SCRIPT_EVAL_MAX_SLOTS; ++i) {
        md_mutex_init(&eval->slots[i].lock);
    }

    md_mutex_init(&eval->frame_lock);

    return eval;
}

md_script_eval_t* md_script_eval_create_desc(const md_script_ir_t* ir, const md_script_eval_desc_t* desc, md_allocator_i* alloc) {
    return eval_create(ir, desc, alloc, NULL);
}

// Copies what src holds into dst, a property of the same kind and shape: its values and the frames they are of
static void adopt_property(eval_property_t* dst, const eval_property_t* src, size_t num_frames) {
    MEMCPY(dst->values, src->values, dst->num_values * sizeof(float));
    if (dst->weights && src->weights) {
        MEMCPY(dst->weights, src->weights, dst->num_values * sizeof(float));
    }
    if (dst->sum_values && src->sum_values) {
        MEMCPY(dst->sum_values,  src->sum_values,  dst->num_values * sizeof(double));
        MEMCPY(dst->sum_weights, src->sum_weights, dst->num_values * sizeof(double));
    }
    if (dst->mean && src->mean) {
        MEMCPY(dst->mean, src->mean, num_frames * sizeof(float));
        MEMCPY(dst->var,  src->var,  num_frames * sizeof(float));
        MEMCPY(dst->ext,  src->ext,  num_frames * sizeof(vec2_t));
    }
    if (dst->frame_values && src->frame_values) {
        MEMCPY(dst->frame_values,  src->frame_values,  num_frames * dst->num_values * sizeof(float));
        MEMCPY(dst->frame_weights, src->frame_weights, num_frames * dst->num_values * sizeof(float));
    }
    if (dst->range && src->range) {
        *dst->range = *src->range;
    }
    dst->min_value   = src->min_value;
    dst->max_value   = src->max_value;
    dst->accum_count = src->accum_count;
    md_bitmask_assign(&dst->done, md_bitmask_view(&src->done));
}

// Of the same kind and shape, and kept the same way: per frame in both or in neither
static bool same_shape(const eval_property_t* a, const eval_property_t* b) {
    return a->kind == b->kind && a->count == b->count && a->num_values == b->num_values && a->values && b->values &&
        (a->frame_values != NULL) == (b->frame_values != NULL);
}

// The property of eval a slot of prev's property j goes to: the one which took over j's values (of the same name
// first, as two properties computed alike hold the same values), whose slots it holds as they are. Otherwise an
// aggregate of j's name, which holds nothing yet. -1 for neither. from[i] is the property of prev that property i
// took over, -1 for none.
static int slot_destination(bool* out_held, const md_script_eval_t* eval, const md_script_eval_t* prev, uint32_t j, const int* from) {
    const eval_property_t* src = &prev->props[j];
    int named = -1, alike = -1;
    for (size_t i = 0; i < md_array_size(eval->props); ++i) {
        const eval_property_t* dst = &eval->props[i];
        const bool same_name = str_eq(dst->ident, src->ident);
        if (from[i] >= 0 && dst->fingerprint == src->fingerprint) {
            if (same_name) { *out_held = true; return (int)i; }
            if (alike < 0) alike = (int)i;
        }
        if (same_name && is_aggregate(dst) && named < 0) named = (int)i;
    }
    *out_held = alike >= 0;
    return alike >= 0 ? alike : named;
}

md_script_eval_t* md_script_eval_create_from(const md_script_eval_t* prev, const md_script_ir_t* ir, const md_script_eval_desc_t* desc, md_allocator_i* alloc) {
    if (!validate_eval(prev)) {
        return md_script_eval_create_desc(ir, desc, alloc);
    }
    md_script_eval_t* eval = eval_create(ir, desc, alloc, prev);
    if (!eval) {
        return eval;
    }

    // The versions carry on from the previous evaluation's: a version of this one never equals one of that one, so
    // nothing derived from a previous property or slot can be taken for what is published now, taken over or not.
    // Before any slot is taken into use below, which takes a version.
    eval->version_counter = MAX(eval->version_counter, md_atomic_load64(&prev->version_counter));

    const size_t num_props  = md_array_size(eval->props);
    md_temp_scope_t temp = md_temp_begin();
    int* from = md_temp_alloc_array(temp, int, MAX(num_props, 1));
    for (size_t i = 0; i < num_props; ++i) from[i] = -1;

    // A property is taken over by what it is computed from, whatever it is called now (taken_over_from). Nothing over
    // another number of frames, which is another run.
    const bool same_run = prev->frame_count == eval->frame_count;
    md_mutex_t* prev_lock = (md_mutex_t*)&prev->frame_lock;
    md_mutex_lock(prev_lock);
    for (size_t i = 0; i < num_props; ++i) {
        eval_property_t* dst = &eval->props[i];
        const int j = taken_over_from(prev, dst->fingerprint, eval->frame_count);
        // Computed alike, it is alike: of the same kind and shape, and kept the same way (eval_create)
        if (j >= 0 && same_shape(dst, &prev->props[j])) {
            adopt_property(dst, &prev->props[j], eval->frame_count);
            from[i] = j;
        }
    }

    // And the slots, asked for what they were asked for: with what they hold where their property was taken over
    for (uint32_t k = 0; k < prev->num_slots_used && same_run; ++k) {
        eval_slot_t* src = &prev->slots[k];
        md_mutex_lock(&src->lock);
        bool held = false;
        const int i = src->active ? slot_destination(&held, eval, prev, src->prop, from) : -1;
        if (i >= 0 && !find_slot(eval, src->key, (uint32_t)i) && eval->num_slots_used < MD_SCRIPT_EVAL_MAX_SLOTS) {
            eval_slot_t* dst = &eval->slots[eval->num_slots_used++];
            slot_take(eval, dst, src->key, (uint32_t)i);
            md_bitmask_assign(&dst->wanted, md_bitmask_view(&src->wanted));
            if (dst->sum_values) {
                if (held) {
                    md_bitmask_assign(&dst->included, md_bitmask_view(&src->included));
                    dst->count = src->count;
                    MEMCPY(dst->sum_values, src->sum_values, dst->num_values * sizeof(double));
                    if (dst->sum_weights) MEMCPY(dst->sum_weights, src->sum_weights, dst->num_values * sizeof(double));
                }
            } else {
                // Kept per frame: it holds what it wants of the frames its property holds
                md_bitmask_assign(&dst->included, md_bitmask_view(&dst->wanted));
                md_bitmask_and(&dst->included, md_bitmask_view(&eval->props[i].done));
                dst->count = (uint32_t)md_bitmask_count(&dst->included);
            }
        }
        md_mutex_unlock(&src->lock);
    }
    md_mutex_unlock(prev_lock);
    md_temp_end(temp);

    md_attributes_t* table = &eval->attributes;
    table->version_counter = MAX(table->version_counter, prev->attributes.version_counter);
    for (md_attribute_iter_t it = md_attributes_iter(table, (str_t){0}); md_attributes_next(&it);) {
        md_attributes_touch(table, it.attr->id);
    }

    return eval;
}

void md_script_eval_clear_data(md_script_eval_t* eval) {
    ASSERT(eval);
    ASSERT(eval->magic == SCRIPT_EVAL_MAGIC);
    md_mutex_lock(&eval->frame_lock);
    for (size_t p = 0; p < md_array_size(eval->props); ++p) {
        clear_property(&eval->props[p]);
    }
    // The slots hold nothing, and are asked for what they were asked for
    for (uint32_t i = 0; i < eval->num_slots_used; ++i) {
        eval_slot_t* s = &eval->slots[i];
        md_mutex_lock(&s->lock);
        if (s->active) {
            md_bitmask_clear(&s->included);
            s->count = 0;
            if (s->sum_values)  MEMSET(s->sum_values,  0, s->num_values * sizeof(double));
            if (s->sum_weights) MEMSET(s->sum_weights, 0, s->num_values * sizeof(double));
            s->version = next_version(eval);
        }
        md_mutex_unlock(&s->lock);
    }
    md_mutex_unlock(&eval->frame_lock);
    eval->interrupt = false;
}

// The frames which every property of 'kinds' holds, or property p alone when p >= 0. Taken under the eval's
// frame_lock, allocated with alloc.
static md_bitset_t frames_held(const md_script_eval_t* eval, uint32_t kinds, int p, md_allocator_i* alloc) {
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_bitmask_t m = md_bitmask_create(0, (uint32_t)eval->frame_count, md_temp_allocator(temp));
    md_bitmask_set_range(&m, 0, (uint32_t)eval->frame_count);
    md_mutex_t* lock = (md_mutex_t*)&eval->frame_lock;
    md_mutex_lock(lock);
    for (size_t i = 0; i < md_array_size(eval->props); ++i) {
        if (p >= 0 && (size_t)p != i) continue;
        const eval_property_t* prop = &eval->props[i];
        if (!(prop->kind & kinds)) continue;
        md_bitmask_and(&m, md_bitmask_view(&prop->done));
    }
    md_mutex_unlock(lock);
    md_bitset_t result = md_bitmask_extract(&m, alloc);
    md_temp_end(temp);
    return result;
}

static const uint32_t all_kinds = MD_SCRIPT_PROPERTY_FLAG_TEMPORAL | MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION | MD_SCRIPT_PROPERTY_FLAG_VOLUME;

md_bitset_t md_script_eval_property_frames(const md_script_eval_t* eval, str_t ident, md_allocator_i* alloc) {
    if (!validate_eval(eval) || !alloc) return md_bitset_empty_set();
    const int p = find_property(eval, ident);
    return p >= 0 ? frames_held(eval, all_kinds, p, alloc) : md_bitset_empty_set();
}

uint64_t md_script_eval_property_fingerprint(const md_script_eval_t* eval, str_t ident) {
    if (!validate_eval(eval)) return 0;
    const int p = find_property(eval, ident);
    return p >= 0 ? eval->props[p].fingerprint : 0;
}

bool md_script_eval_plan(md_script_eval_plan_t* out_plans, size_t cap, md_bitset_t* out_pending, const md_script_eval_t* eval, const md_script_ir_t* ir, size_t num_frames, md_allocator_i* alloc) {
    if (!out_plans || !md_script_ir_valid(ir) || !alloc || cap < md_array_size(ir->property_names)) return false;
    const bool have = validate_eval(eval);
    const uint32_t n = (uint32_t)num_frames;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_bitmask_t pending = md_bitmask_create(0, n, temp_alloc);
    md_bitmask_t mine    = md_bitmask_create(0, n, temp_alloc);
    md_bitmask_t moving  = md_bitmask_create(0, n, temp_alloc);

    // As md_script_eval_create_from would have it: what is computed alike is taken over (taken_over_from), with the
    // frames it holds and its slots as they are, and kept the way it is kept. Anything else starts over.
    md_mutex_t* lock = have ? (md_mutex_t*)&eval->frame_lock : NULL;
    if (lock) md_mutex_lock(lock);
    for (size_t i = 0; i < md_array_size(ir->property_names); ++i) {
        md_script_eval_plan_t* plan = &out_plans[i];
        MEMSET(plan, 0, sizeof(*plan));
        const uint64_t fingerprint = i < md_array_size(ir->property_fingerprints) ? ir->property_fingerprints[i] : 0;
        const int j = have ? taken_over_from(eval, fingerprint, num_frames) : -1;
        if (j < 0) {
            plan->status = (have && find_property(eval, ir->property_names[i]) >= 0) ? MD_SCRIPT_EVAL_STATUS_CHANGED : MD_SCRIPT_EVAL_STATUS_NEW;
            plan->num_pending = n;
            md_bitmask_set_range(&pending, 0, n);
            continue;
        }
        const eval_property_t* src = &eval->props[j];
        plan->status     = MD_SCRIPT_EVAL_STATUS_HELD;
        plan->from       = str_copy(src->ident, alloc);
        plan->num_frames = (uint32_t)md_bitmask_count(&src->done);

        // The frames it does not hold, and those its slots have to add or take out
        md_bitmask_clear(&mine);
        md_bitmask_set_range(&mine, 0, n);
        md_bitmask_andnot(&mine, md_bitmask_view(&src->done));
        for (uint32_t k = 0; k < eval->num_slots_used; ++k) {
            eval_slot_t* s = &eval->slots[k];
            if (!s->active || !s->sum_values || s->prop != (uint32_t)j) continue;
            md_mutex_lock(&s->lock);
            md_bitmask_assign(&moving, md_bitmask_view(&s->wanted));
            md_bitmask_xor(&moving, md_bitmask_view(&s->included));
            md_bitmask_or(&mine, md_bitmask_view(&moving));
            md_mutex_unlock(&s->lock);
        }
        plan->num_pending = (uint32_t)md_bitmask_count(&mine);
        md_bitmask_or(&pending, md_bitmask_view(&mine));
    }
    if (lock) md_mutex_unlock(lock);

    if (out_pending) {
        *out_pending = md_bitmask_extract(&pending, alloc);
    }
    md_temp_end(temp);
    return true;
}

md_bitset_t md_script_eval_pending_frames(md_script_eval_t* eval, md_allocator_i* alloc) {
    if (!validate_eval(eval) || !alloc) return md_bitset_empty_set();
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    const uint32_t num_frames = (uint32_t)eval->frame_count;
    md_bitmask_t held    = md_bitmask_create(0, num_frames, temp_alloc);
    md_bitmask_t pending = md_bitmask_create(0, num_frames, temp_alloc);
    md_bitmask_t moving  = md_bitmask_create(0, num_frames, temp_alloc);
    md_bitmask_set_range(&held,    0, num_frames);
    md_bitmask_set_range(&pending, 0, num_frames);

    md_mutex_lock(&eval->frame_lock);
    // The frames some property does not hold
    for (size_t p = 0; p < md_array_size(eval->props); ++p) {
        md_bitmask_and(&held, md_bitmask_view(&eval->props[p].done));
    }
    md_bitmask_andnot(&pending, md_bitmask_view(&held));
    // and those a slot has to add or take out
    for (uint32_t i = 0; i < eval->num_slots_used; ++i) {
        eval_slot_t* s = &eval->slots[i];
        if (!s->active || !s->sum_values) continue;
        md_mutex_lock(&s->lock);
        {
            md_bitmask_assign(&moving, md_bitmask_view(&s->wanted));
            md_bitmask_xor(&moving, md_bitmask_view(&s->included));
            md_bitmask_or(&pending, md_bitmask_view(&moving));
        }
        md_mutex_unlock(&s->lock);
    }
    md_mutex_unlock(&eval->frame_lock);

    md_bitset_t result = md_bitmask_extract(&pending, alloc);
    md_temp_end(temp);
    return result;
}

// ### SLOTS ###

bool md_script_eval_slot_request(md_script_eval_t* eval, uint64_t key, str_t ident, const md_bitset_t* frames) {
    if (!validate_eval(eval) || !frames) return false;
    const int p = find_property(eval, ident);
    if (p < 0 || !is_aggregate(&eval->props[p])) return false;
    const eval_property_t* prop = &eval->props[p];

    md_temp_scope_t temp = md_temp_begin();
    const md_bitset_t want = md_bitset_and(*frames, md_bitset_range(0, (uint32_t)eval->frame_count), md_temp_allocator(temp));

    // Under frame_lock: a slot is taken into use under it, and one kept per frame holds what its property does
    md_mutex_lock(&eval->frame_lock);
    eval_slot_t* s = find_slot(eval, key, (uint32_t)p);
    if (s) {
        md_mutex_lock(&s->lock);
    } else {
        for (uint32_t i = 0; i < MD_SCRIPT_EVAL_MAX_SLOTS && !s; ++i) {
            if (!eval->slots[i].active) s = &eval->slots[i];
        }
        if (s) {
            md_mutex_lock(&s->lock);
            slot_take(eval, s, key, (uint32_t)p);
            eval->num_slots_used = MAX(eval->num_slots_used, (uint32_t)(s - eval->slots) + 1);
        }
    }
    if (s) {
        if (!md_bitset_equal(md_bitmask_view(&s->wanted), want)) {
            md_bitmask_assign(&s->wanted, want);
            if (!s->sum_values) {
                // Kept per frame: it holds at once what it wants of the frames evaluated
                md_bitmask_assign(&s->included, md_bitmask_view(&s->wanted));
                md_bitmask_and(&s->included, md_bitmask_view(&prop->done));
                s->count = (uint32_t)md_bitmask_count(&s->included);
            }
            s->version = next_version(eval);
        }
        md_mutex_unlock(&s->lock);
    }
    md_mutex_unlock(&eval->frame_lock);
    md_temp_end(temp);

    if (!s) {
        MD_LOG_ERROR("Script eval: no room for another slot, all %d are in use", MD_SCRIPT_EVAL_MAX_SLOTS);
    }
    return s != NULL;
}

void md_script_eval_slot_release(md_script_eval_t* eval, uint64_t key, str_t ident) {
    if (!validate_eval(eval)) return;
    const int p = find_property(eval, ident);
    if (p < 0) return;
    md_mutex_lock(&eval->frame_lock);
    eval_slot_t* s = find_slot(eval, key, (uint32_t)p);
    if (s) {
        md_mutex_lock(&s->lock);
        slot_let_go(s);
        md_mutex_unlock(&s->lock);
    }
    md_mutex_unlock(&eval->frame_lock);
}

// Finds the slot and returns it locked, NULL if there is none
static eval_slot_t* lock_slot(const md_script_eval_t* eval, uint64_t key, str_t ident) {
    const int p = find_property(eval, ident);
    if (p < 0) return NULL;
    md_mutex_t* frame_lock = (md_mutex_t*)&eval->frame_lock;
    md_mutex_lock(frame_lock);
    eval_slot_t* s = find_slot(eval, key, (uint32_t)p);
    if (s) md_mutex_lock(&s->lock);
    md_mutex_unlock(frame_lock);
    return s;
}

// What slot s holds. The caller holds its lock.
static void slot_info(md_script_eval_slot_info_t* out, const eval_slot_t* s) {
    out->version    = s->version;
    out->num_frames = s->count;
    out->num_values = s->num_values;
    const uint32_t wanted = (uint32_t)md_bitmask_count(&s->wanted);
    if (s->sum_values) {
        // Those it wants and does not hold, and those it holds and does not want
        md_temp_scope_t temp = md_temp_begin();
        md_bitmask_t both = md_bitmask_create(0, s->wanted.end, md_temp_allocator(temp));
        md_bitmask_assign(&both, md_bitmask_view(&s->wanted));
        md_bitmask_and(&both, md_bitmask_view(&s->included));
        out->num_pending = wanted + s->count - 2 * (uint32_t)md_bitmask_count(&both);
        md_temp_end(temp);
    } else {
        // It holds only frames it wants
        out->num_pending = wanted - s->count;
    }
}

bool md_script_eval_slot_info(md_script_eval_slot_info_t* out_info, const md_script_eval_t* eval, uint64_t key, str_t ident) {
    if (!validate_eval(eval) || !out_info) return false;
    eval_slot_t* s = lock_slot(eval, key, ident);
    if (!s) return false;
    slot_info(out_info, s);
    md_mutex_unlock(&s->lock);
    return true;
}

// Means of the sums over count frames; with none, values of 0 and weights of 1, as a published aggregate starts out
static void write_means(float* out_values, float* out_weights, const double* sum_values, const double* sum_weights, size_t num_values, uint32_t count) {
    const double scl = count ? 1.0 / (double)count : 0.0;
    for (size_t i = 0; i < num_values; ++i) {
        out_values[i] = (float)(sum_values[i] * scl);
    }
    if (out_weights && sum_weights) {
        for (size_t i = 0; i < num_values; ++i) {
            out_weights[i] = count ? (float)(sum_weights[i] * scl) : 1.0f;
        }
    }
}

bool md_script_eval_slot_extract(md_script_eval_slot_info_t* out_info, float* out_values, float* out_weights, size_t num_values,
                                 const md_script_eval_t* eval, uint64_t key, str_t ident) {
    if (!validate_eval(eval)) return false;
    eval_slot_t* s = lock_slot(eval, key, ident);
    if (!s) return false;

    md_script_eval_slot_info_t info = {0};
    slot_info(&info, s);
    if (out_info) *out_info = info;
    if (!out_values || num_values < s->num_values) {
        md_mutex_unlock(&s->lock);
        return false;
    }
    const size_t n = s->num_values;

    if (s->sum_values) {
        // Added up as frames are evaluated: copied under the lock, which holds back only adding a frame to it
        write_means(out_values, out_weights, s->sum_values, s->sum_weights, n, s->count);
        md_mutex_unlock(&s->lock);
        return true;
    }

    // Kept per frame: summed outside the lock, over the frames it holds as they were when it was taken. What is kept
    // of a frame is written before the frame is held, and not written again.
    md_temp_scope_t temp = md_temp_begin();
    const md_bitset_t frames = md_bitmask_extract(&s->included, md_temp_allocator(temp));
    const uint32_t count = s->count;
    const eval_property_t* prop = &eval->props[s->prop];
    md_mutex_unlock(&s->lock);

    double* sum_values  = md_temp_alloc_array(temp, double, n);
    double* sum_weights = md_temp_alloc_array(temp, double, n);
    MEMSET(sum_values,  0, n * sizeof(double));
    MEMSET(sum_weights, 0, n * sizeof(double));
    md_bitset_iter_t it = md_bitset_iter(frames);
    while (md_bitset_iter_next(&it)) {
        const float* values  = prop->frame_values  + (size_t)it.idx * n;
        const float* weights = prop->frame_weights + (size_t)it.idx * n;
        for (size_t i = 0; i < n; ++i) {
            sum_values [i] += values [i];
            sum_weights[i] += weights[i];
        }
    }
    write_means(out_values, out_weights, sum_values, sum_weights, n, count);
    md_temp_end(temp);
    return true;
}

void md_script_eval_reset_interrupt(md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        eval->interrupt = false;
    }
}

bool md_script_eval_frame_range(md_script_eval_t* eval, const struct md_script_ir_t* ir, const struct md_system_t* sys, str_t run, uint32_t frame_beg, uint32_t frame_end) {
    if (frame_beg > frame_end) {
        MD_LOG_ERROR("Script eval: Invalid frame range");
        return false;
    }
    md_temp_scope_t temp = md_temp_begin();
    const uint32_t count = frame_end - frame_beg;
    uint32_t* frames = md_temp_alloc_array(temp, uint32_t, MAX(count, 1));
    for (uint32_t i = 0; i < count; ++i) {
        frames[i] = frame_beg + i;
    }
    const bool result = md_script_eval_frames(eval, ir, sys, run, frames, count);
    md_temp_end(temp);
    return result;
}

bool md_script_eval_frames(md_script_eval_t* eval, const struct md_script_ir_t* ir, const struct md_system_t* sys, str_t run, const uint32_t* frames, size_t num_frames_to_eval) {
    ASSERT(eval);

    if (!ir) {
        MD_LOG_ERROR("Script eval: Immediate representation was null");
        return false;
    }
    if (!sys) {
        MD_LOG_ERROR("Script eval: Molecule was null");
        return false;
    }

    const md_attribute_t* run_axis = str_empty(run) ? NULL : md_attributes_find_in(&sys->attributes, run, STR_LIT("time"));
    const uint32_t num_frames = run_axis ? run_axis->format.shape[0] : 0;
    if (num_frames == 0) {
        MD_LOG_ERROR("Script eval: no frames in the run '" STR_FMT "'", STR_ARG(run));
        return false;
    }
    for (size_t i = 0; i < num_frames_to_eval; ++i) {
        if (frames[i] >= num_frames || frames[i] >= eval->frame_count) {
            MD_LOG_ERROR("Script eval: Frame %u is out of range", frames[i]);
            return false;
        }
    }

    if (md_array_size(eval->props) == 0) {
        MD_LOG_INFO("Script eval: No properties present, nothing to evaluate");
        return false;
    }
    
    bool result = eval_properties(eval, sys, run, ir, frames, num_frames_to_eval);

    // The buffers were written through md_attributes_data, which deliberately does not bump a
    // version. This is the producer saying it is done. Under the lock, as several threads evaluating at once finish
    // at once, and the table's counter is not for them to race on.
    md_mutex_lock(&eval->frame_lock);
    for (md_attribute_iter_t it = md_attributes_iter(&eval->attributes, (str_t){0}); md_attributes_next(&it);) {
        md_attributes_touch(&eval->attributes, it.attr->id);
    }
    md_mutex_unlock(&eval->frame_lock);

    return result;
}

const md_attributes_t* md_script_eval_attributes(const md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        return &eval->attributes;
    }
    return NULL;
}

uint64_t md_script_eval_ir_fingerprint(const md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        return eval->ir_fingerprint;
    }
    return 0;
}

void md_script_eval_free(md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        for (size_t i = 0; i < md_array_size(eval->props); ++i) {
            md_mutex_destroy(&eval->props[i].accum_mutex);
        }
        for (size_t i = 0; i < MD_SCRIPT_EVAL_MAX_SLOTS; ++i) {
            slot_let_go(&eval->slots[i]);
            md_mutex_destroy(&eval->slots[i].lock);
        }
        md_mutex_destroy(&eval->frame_lock);
        md_arena_allocator_destroy(eval->arena);
    }
}

size_t md_script_eval_property_count(const md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        return md_array_size(eval->props);
    }
    return 0;
}

size_t md_script_eval_frame_count(const md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        return eval->frame_count;
    }
    return 0;
}

md_bitset_t md_script_eval_frame_mask(const md_script_eval_t* eval, md_allocator_i* alloc) {
    if (!validate_eval(eval) || !alloc) return md_bitset_empty_set();
    bool any_temporal = false;
    for (size_t i = 0; i < md_array_size(eval->props); ++i) {
        any_temporal |= eval->props[i].kind == MD_SCRIPT_PROPERTY_FLAG_TEMPORAL;
    }
    return frames_held(eval, any_temporal ? MD_SCRIPT_PROPERTY_FLAG_TEMPORAL : all_kinds, -1, alloc);
}

void md_script_eval_interrupt(md_script_eval_t* eval) {
    if (validate_eval(eval)) {
        eval->interrupt = true;
    }
}

static bool eval_expression(data_t* dst, str_t expr, md_system_t* sys, md_allocator_i* alloc) {
    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);

    md_script_ir_t* ir = create_ir(temp_alloc);
    ir->str = str_copy(expr, ir->arena);

    tokenizer_t tokenizer = tokenizer_init(ir->str);
    bool result = false;

    eval_log_t log = { .alloc = temp_alloc };

    ast_node_t* node = prune_expressions(parse_expression(&(parse_context_t){ .ir = ir, .tokenizer = &tokenizer, .temp_alloc = temp_alloc}));
    if (node) {
        eval_context_t ctx = {
            .ir = ir,
            .compile_ir = ir,
            .sys = sys,
            .temp_alloc = temp_alloc,
            .alloc = temp_alloc,
            .cur_state = &sys->reference,
            .ref_state = &sys->reference,
        };

        if (static_check_node(node, &ctx)) {
            ctx.compile_ir = NULL;
            ctx.log = &log;
            allocate_data(dst, node->data.type, alloc);
            if (evaluate_node(dst, node, &ctx)) {
                if (dst->type.base_type == TYPE_STRING) {
                    // @HACK: We reallocate strings here here: If the data type is a str_t, then it gets a shallow copy
                    // Which means that the actual string data is contained within the ir->arena => temp_alloc
                    const size_t num_elem = element_count(*dst);
                    str_t* dst_str = dst->ptr;
                    for (size_t i = 0; i < num_elem; ++i) {
                        dst_str[i] = str_copy(dst_str[i], alloc);
                    }
                }
                result = true;
            }
        }
    }

    for (size_t i = 0; i < md_array_size(ir->errors); ++i) {
        MD_LOG_ERROR(""STR_FMT"", ir->errors[i].text.len, ir->errors[i].text.ptr);
    }
    for (size_t i = 0; i < md_array_size(log.errors); ++i) {
        MD_LOG_ERROR(""STR_FMT"", log.errors[i].text.len, log.errors[i].text.ptr);
    }

    md_temp_end(temp_scope);

    return result;
}

#if DEBUG
static void parse_type_check_and_print_expression_to_json(str_t expr, const md_system_t* sys, str_t filename) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
    md_script_ir_t* ir = create_ir(temp_alloc);
    ir->str = str_copy(expr, ir->arena);

    tokenizer_t tokenizer = tokenizer_init(ir->str);
    ast_node_t* node = prune_expressions(parse_expression(&(parse_context_t){ .ir = ir, .tokenizer = &tokenizer, .temp_alloc = temp_alloc}));
    if (node) {
        eval_context_t ctx = {
            .ir = ir,
            .compile_ir = ir,
            .sys = sys,
            .temp_alloc = temp_alloc,
            .alloc = temp_alloc,
            .cur_state = &sys->reference,
            .ref_state = &sys->reference,
            .eval_flags = EVAL_FLAG_NO_STATIC_EVAL,
        };

        if (static_check_node(node, &ctx)) {
            ast_to_json(node, filename);
        }
    }

    md_temp_end(temp_scope);
}
#endif

// The front half shared by md_filter and md_filter_evaluate: parse expr and type check it in ctx.
// Returns the checked node, or NULL if it did not compile. A variable length is left to the evaluation,
// which determines it (evaluate_node_alloc). What went wrong in the compilation is recorded in ir->errors.
// ctx is left ready for the evaluation: reading ir, and logging to ctx->log.
static ast_node_t* filter_compile(md_script_ir_t* ir, eval_context_t* ctx, str_t expr, const md_script_ir_t* ctx_ir) {
    ir->str = str_copy(expr, ir->arena);
    if (ctx_ir) {
        add_ir_ctx(ir, ctx_ir);
    }
    ir->stage = "Filter evaluate";
    ir->record_log = true;

    tokenizer_t tokenizer = tokenizer_init(ir->str);
    parse_context_t parse_ctx = {
        .ir = ir,
        .tokenizer = &tokenizer,
        .temp_alloc = ctx->temp_alloc,
    };

    const size_t num_errors = md_array_size(ir->errors);
    ast_node_t* node = prune_expressions(parse_expression(&parse_ctx));
    if (md_array_size(ir->errors) != num_errors) {
        node = NULL;    // As for a statement of a script: the tree of a parse which logged an error is not what was written
    }
    const ast_node_t* nested = node ? find_nested_assignment(node) : NULL;
    if (nested) {
        LOG_ERROR(ir, nested->token, "An assignment is a statement of its own, it cannot be a part of an expression");
        node = NULL;
    }
    ctx->compile_ir = ir;
    const bool checked = node && static_check_node(node, ctx) && md_array_size(ir->errors) == num_errors;
    ctx->compile_ir = NULL;
    return checked ? node : NULL;
}

// The errors of the compilation, then those of the evaluation, one line each, and where the first one is. What
// went wrong otherwise (msg) only when nothing more specific did.
static void filter_write_errors(md_filter_result_t* res, const md_script_ir_t* ir, const eval_log_t* log, const char* msg) {
    const md_log_token_t* lists[2] = { ir->errors, log->errors };
    const size_t cap = sizeof(res->error);
    size_t len = 0;
    bool first = true;
    for (size_t l = 0; l < ARRAY_SIZE(lists); ++l) {
        for (size_t i = 0; i < md_array_size(lists[l]) && len + 1 < cap; ++i) {
            const int n = md_snprintf(res->error + len, cap - len, "%s" STR_FMT, len ? "\n" : "", STR_ARG(lists[l][i].text));
            if (n < 0) return;
            len = MIN(len + (size_t)n, cap - 1);
            if (first) {
                res->error_range = lists[l][i].range;
                first = false;
            }
        }
    }
    if (len == 0 && msg) {
        md_snprintf(res->error, cap, "%s", msg);
    }
}

// What both filters share: the context of the evaluation, and the expression compiled in it
typedef struct filter_eval_t {
    md_script_ir_t* ir;
    eval_log_t      log;
    eval_context_t  ctx;
    ast_node_t*     node;   // NULL if it did not compile
} filter_eval_t;

static bool filter_begin(filter_eval_t* f, md_filter_result_t* res, str_t expr, const md_filter_context_t* ctx, uint32_t eval_flags, md_allocator_i* temp_alloc) {
    MEMSET(res, 0, sizeof(*res));
    if (!ctx || !ctx->sys) {
        md_snprintf(res->error, sizeof(res->error), "No system to pick atoms from");
        return false;
    }
    if (ctx->sys->atom.count == 0) {
        md_snprintf(res->error, sizeof(res->error), "The system has no atoms");
        return false;
    }
    // Nothing written picks everything, as 'all' does
    if (str_empty(str_trim(expr))) {
        expr = STR_LIT("all");
    }

    f->ir = create_ir(temp_alloc);
    f->log = (eval_log_t){ .alloc = temp_alloc };
    f->ctx = (eval_context_t){
        .ir = f->ir,
        .log = &f->log,
        .sys = ctx->sys,
        // The context the left hand side of 'in' is compiled and evaluated in: the static validation ("did not match
        // any residue within the context") sees it as much as the evaluation does
        .mol_ctx = ctx->atoms,
        .temp_alloc = temp_alloc,
        .alloc = temp_alloc,
        .cur_state = ctx->state ? ctx->state : &ctx->sys->reference,
        .ref_state = &ctx->sys->reference,
        .eval_flags = eval_flags,
    };
    f->node = filter_compile(f->ir, &f->ctx, expr, ctx->ir);
    if (f->node) {
        res->is_dynamic = (f->node->flags & FLAG_DYNAMIC) != 0;
    }
    return true;
}

// A filter of atoms which gives bonds picks the atoms of those bonds, as atom(b) does, within the context. The one place
// where a bond set becomes atoms without being asked to: within an expression, bonds and atoms do not mix.
static md_bitset_t filter_bond_atoms(const bondset_t* sets, size_t count, const md_filter_context_t* ctx, md_allocator_i* alloc, md_allocator_i* temp) {
    md_bitmask_t b = md_bitmask_create(0, (uint32_t)ctx->sys->atom.count, temp);
    for (size_t i = 0; i < count; ++i) {
        bondset_atoms_into(&b, &sets[i], ctx->sys);
    }
    if (ctx->atoms) md_bitmask_and(&b, *ctx->atoms);
    return md_bitmask_extract(&b, alloc);
}

md_filter_result_t md_filter_sets(md_array(md_bitset_t)* dst, str_t expr, const md_filter_context_t* ctx, md_allocator_i* alloc) {
    ASSERT(dst);
    ASSERT(alloc);
    md_filter_result_t res;

    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);

    filter_eval_t f;
    if (filter_begin(&f, &res, expr, ctx, EVAL_FLAG_NO_STATIC_EVAL, temp_alloc)) {
        const char* msg = NULL;
        if (f.node && f.node->data.type.base_type == TYPE_BONDSET) {
            // The atoms of each set of bonds
            data_t data = {0};
            if (evaluate_node_alloc(&data, f.node, &f.ctx, temp_alloc)) {
                const int64_t len = data.ptr ? type_info_array_len(data.type) : 0;
                const bondset_t* sets = data.ptr;
                for (int64_t i = 0; i < len; ++i) {
                    md_array_push(*dst, filter_bond_atoms(&sets[i], 1, ctx, alloc, temp_alloc), alloc);
                }
                res.success = true;
            }
        } else if (f.node && f.node->data.type.base_type == TYPE_ATOMSET) {
            data_t data = {0};
            if (evaluate_node_alloc(&data, f.node, &f.ctx, temp_alloc)) {
                const int64_t len = type_info_array_len(data.type);
                const md_bitset_t* bs_arr = data.ptr;
                for (int64_t i = 0; bs_arr && i < len; ++i) {
                    md_array_push(*dst, md_bitset_copy(bs_arr[i], alloc), alloc);
                }
                res.success = true;
            }
        } else if (f.node) {
            msg = "The expression picks neither atoms nor bonds";
        }
        if (!res.success) {
            filter_write_errors(&res, f.ir, &f.log, msg);
        }
    }

    md_temp_end(temp_scope);
    return res;
}

md_filter_result_t md_filter(md_bitset_t* dst, str_t expr, const md_filter_context_t* ctx, md_allocator_i* alloc) {
    ASSERT(dst);
    ASSERT(alloc);
    md_filter_result_t res;

    md_bitset_free(dst, alloc);

    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);

    filter_eval_t f;
    if (filter_begin(&f, &res, expr, ctx, EVAL_FLAG_FLATTEN | EVAL_FLAG_NO_STATIC_EVAL, temp_alloc)) {
        const char* msg = NULL;
        ast_node_t* node = f.node;
        if (node && node->data.type.base_type == TYPE_BONDSET) {
            // The atoms of the bonds
            data_t data = {0};
            if (evaluate_node_alloc(&data, node, &f.ctx, temp_alloc)) {
                const size_t len = data.ptr ? (size_t)type_info_array_len(data.type) : 0;
                *dst = filter_bond_atoms((const bondset_t*)data.ptr, len, ctx, alloc, temp_alloc);
                res.success = true;
            }
        } else if (node) {
            const bool is_bitfield = node->data.type.base_type == TYPE_ATOMSET;
            const bool variable = is_bitfield && is_variable_length(node->data.type);

            // A variable length is known once evaluated
            data_t data = {0};
            bool evaluated = false;
            if (variable) {
                evaluated = evaluate_node_alloc(&data, node, &f.ctx, temp_alloc);
            }
            const int len = !is_bitfield ? 0 : variable ? (evaluated ? (int)type_info_array_len(data.type) : -1) : (int)type_info_array_len(node->data.type);

            if (variable && !evaluated) {
                res.success = false;
            } else if (len == 1 && !evaluated) {
                // A single set needs no array. Evaluated in temp like everything in a script (whose data lives in
                // arenas), then copied out to the caller's allocator, which may well be a heap.
                md_bitset_t result = {0};
                data = (data_t){
                    .type = node->data.type,
                    .ptr  = &result,
                    .size = sizeof(md_bitset_t),
                    .alloc = temp_alloc,
                };
                res.success = evaluate_node(&data, node, &f.ctx);
                *dst = res.success ? md_bitset_copy(result, alloc) : md_bitset_empty_set();
            } else if (len >= 1) {
                res.success = evaluated || (allocate_data(&data, node->data.type, temp_alloc) && evaluate_node(&data, node, &f.ctx));
                if (res.success) {
                    const md_bitset_t* src = (const md_bitset_t*)data.ptr;
                    md_bitmask_t b = md_bitmask_create(0, (uint32_t)ctx->sys->atom.count, temp_alloc);
                    for (int i = 0; i < len; ++i) {
                        md_bitmask_or(&b, src[i]);
                    }
                    *dst = md_bitmask_extract(&b, alloc);
                }
            } else {
                msg = is_bitfield ? "The expression does not pick a valid set of atoms" : "The expression picks neither atoms nor bonds";
            }
        }
        if (!res.success) {
            filter_write_errors(&res, f.ir, &f.log, msg);
        }
    }

    md_temp_end(temp_scope);
    return res;
}

md_filter_result_t md_filter_bonds(md_array(md_atom_pair_t)* dst, str_t expr, const md_filter_context_t* ctx, md_allocator_i* alloc) {
    ASSERT(dst);
    ASSERT(alloc);
    md_filter_result_t res;

    md_array_shrink(*dst, 0);

    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);

    filter_eval_t f;
    if (filter_begin(&f, &res, expr, ctx, EVAL_FLAG_FLATTEN | EVAL_FLAG_NO_STATIC_EVAL, temp_alloc)) {
        const char* msg = NULL;
        const base_type_t type = f.node ? f.node->data.type.base_type : TYPE_UNDEFINED;
        if (type == TYPE_BONDSET || type == TYPE_ATOMSET) {
            data_t data = {0};
            if (evaluate_node_alloc(&data, f.node, &f.ctx, temp_alloc)) {
                const size_t len = data.ptr ? (size_t)type_info_array_len(data.type) : 0;
                if (type == TYPE_BONDSET) {
                    md_atom_pair_t* pairs = NULL;
                    const size_t n = len ? bondset_union_pairs(&pairs, (const bondset_t*)data.ptr, len, &f.ctx, temp_alloc) : 0;
                    for (size_t i = 0; i < n; ++i) {
                        // Within the context, which bond() keeps to but the bonds of bond_delta() need not: those of
                        // its scope reach out of it
                        if (ctx->atoms && !(md_bitset_test(*ctx->atoms, (uint32_t)pairs[i].idx[0]) && md_bitset_test(*ctx->atoms, (uint32_t)pairs[i].idx[1]))) continue;
                        md_array_push(*dst, pairs[i], alloc);
                    }
                } else if (len) {
                    // The persistent bonds among the atoms, of both ranges and of every kind, as bond() gives them: every
                    // atom ('all', nothing written) is every bond. The elastic bonds too, which licorice and ball and
                    // stick leave out unless asked to: a representation of bonds shows what its filter picks.
                    const md_system_t* sys = ctx->sys;
                    md_bitmask_t b = md_bitmask_create(0, (uint32_t)sys->atom.count, temp_alloc);
                    const md_bitset_t* sets = (const md_bitset_t*)data.ptr;
                    for (size_t i = 0; i < len; ++i) {
                        md_bitmask_or(&b, sets[i]);
                    }
                    const md_bitset_t atoms = md_bitmask_extract(&b, temp_alloc);
                    const size_t total = sys->bond.count + sys->bond.noncovalent_count;
                    for (size_t i = 0; i < total; ++i) {
                        const md_atom_pair_t p = sys->bond.pairs[i];
                        if (md_bitset_test(atoms, (uint32_t)p.idx[0]) && md_bitset_test(atoms, (uint32_t)p.idx[1])) {
                            md_array_push(*dst, p, alloc);
                        }
                    }
                }
                res.success = true;
            }
        } else if (f.node) {
            msg = "The expression picks neither bonds nor atoms";
        }
        if (!res.success) {
            filter_write_errors(&res, f.ir, &f.log, msg);
        }
    }

    md_temp_end(temp_scope);
    return res;
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
        if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 1 && node->children[0])) return;
        node = node->children[0];
        break;
    case AST_ASSIGNMENT:
        if (!INTERNAL_CHECK(ctx, node->token, md_array_size(node->children) == 2 && node->children[1])) return;
        node = node->children[1];
        break;
    default:
        break;
    }

    if (node->data.type.base_type == TYPE_ATOMSET) {
        data_t data = {0};

        if (!evaluate_node_alloc(&data, node, ctx, ctx->temp_alloc)) {
            md_log(MD_LOG_TYPE_DEBUG, "Vis Eval: Failed to evaluate expression");
            return;
        }
        if (data.ptr) {
            const md_bitset_t* bs_arr = data.ptr;
            ASSERT(ctx->vis_atom_mask);
            // Each element visualized is a group of its own
            md_script_vis_t* vis = ctx->vis;
            if (ctx->subscript_ranges) {
                int count = (int)element_count(data);
                for (size_t i = 0; i < md_array_size(ctx->subscript_ranges); ++i) {
                    irange_t rng = ctx->subscript_ranges[i];
                    for (int j = rng.beg; j < rng.end; ++j) {
                        if (j >= count) break;
                        md_bitmask_or(ctx->vis_atom_mask, bs_arr[j]);
                        if (vis) md_array_push(vis->groups, md_bitset_copy(bs_arr[j], vis->alloc), vis->alloc);
                    }
                }
            } else {
                for (size_t i = 0; i < element_count(data); ++i) {
                    md_bitmask_or(ctx->vis_atom_mask, bs_arr[i]);
                    if (vis) md_array_push(vis->groups, md_bitset_copy(bs_arr[i], vis->alloc), vis->alloc);
                }
            }
            free_data(&data, ctx->temp_alloc);
        } else {
            md_log(MD_LOG_TYPE_DEBUG,"Failed to allocate data for atomset visualization");
        }
    } else if (node->data.type.base_type == TYPE_BONDSET) {
        // The atoms of the bonds, and a line along each
        data_t data = {0};
        if (!evaluate_node_alloc(&data, node, ctx, ctx->temp_alloc)) {
            md_log(MD_LOG_TYPE_DEBUG, "Vis Eval: Failed to evaluate expression");
            return;
        }
        const bondset_t* bs_arr = data.ptr;
        const int count = (int)element_count(data);
        for (int i = 0; bs_arr && i < count; ++i) {
            if (ctx->subscript_ranges) {
                bool picked = false;
                for (size_t r = 0; r < md_array_size(ctx->subscript_ranges) && !picked; ++r) {
                    picked = ctx->subscript_ranges[r].beg <= i && i < ctx->subscript_ranges[r].end;
                }
                if (!picked) continue;
            }
            visualize_bondset(&bs_arr[i], ctx);
        }
        free_data(&data, ctx->temp_alloc);
    } else if (node->flags & FLAG_COORDINATE) {
        data_t tmp_data = {0};
        const data_t* data = &node->data;
        if (!node->data.ptr) {
            // Evaluate the non constant value expression here into a data_t item and visualize it.
            if (!evaluate_node_alloc(&tmp_data, node, ctx, ctx->temp_alloc)) {
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

}

bool md_script_vis_eval_ref(md_script_vis_t* vis, md_script_vis_ref_t ref, int subidx, const md_script_vis_ctx_t* vis_ctx, md_script_vis_flags_t flags) {
    ASSERT(vis);
    ASSERT(vis_ctx);

    // Resolved against the IR it is visualized in, so it can only be the node it was made for
    const ast_node_t* node = resolve_vis_ref(vis_ctx->ir, ref);
    if (!node) {
        return false;
    }
    
    if (vis->magic != VIS_MAGIC) {
        MD_LOG_ERROR("Visualize: vis object not initialized.");
        return false;
    }

    if (subidx < -1) {
        MD_LOG_ERROR("Visualize: invalid subidx");
        return false;
    }
    // An element past the end resolves to nothing, like a stale reference: the procedures index by it unchecked
    if (subidx >= 0 && node->data.type.dim[0] >= 0 && subidx >= node->data.type.dim[0]) {
        return false;
    }

    if (flags == 0) flags = 0xFFFFFFFFU;
    
    md_bitset_free(&vis->atom_mask, vis->alloc);

    md_temp_scope_t temp = md_temp_begin_avoid(vis->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_bitmask_t atom_mask = md_bitmask_create(0, (uint32_t)vis_ctx->sys->atom.count, temp_alloc);


    md_array(irange_t) ranges = 0;
    if (subidx > -1) {
        irange_t range = {subidx, subidx+1};
        md_array_push(ranges, range, temp_alloc);
    }

    eval_context_t ctx = {
        .ir = vis_ctx->ir,  // Read-only, possibly shared with evaluations running on other threads
        .sys = vis_ctx->sys,
        .temp_alloc = temp_alloc,
        .vis = vis,
        .vis_flags = flags,
        .vis_structure = &atom_mask,
        .vis_atom_mask = &atom_mask,
        .cur_state = vis_ctx->state,
        .ref_state = &vis_ctx->sys->reference,
        .subscript_ranges = ranges,
    };

    visualize_node(node, &ctx);

    // This just to see that we conceptually did not make any errors when evaluating sub_contexts
    ASSERT(ctx.vis_structure == &atom_mask);

    // Append all structures into the 'global' atom mask
    for (size_t i = 0; i < md_array_size(vis->structures); ++i) {
        md_bitmask_or(&atom_mask, vis->structures[i]);
    }
    vis->atom_mask = md_bitmask_extract(&atom_mask, vis->alloc);

    md_temp_end(temp);

    return true;
}

size_t md_script_vis_value_text(char* buf, size_t cap, md_script_vis_ref_t ref, int subidx, const md_script_vis_ctx_t* vis_ctx) {
    ASSERT(vis_ctx);
    if (!buf || cap == 0) return 0;
    buf[0] = '\0';

    const ast_node_t* node = resolve_vis_ref(vis_ctx->ir, ref);
    if (!node || !vis_ctx->sys || !vis_ctx->state) return 0;
    switch (node->type) {
    case AST_IDENTIFIER:
    case AST_ASSIGNMENT:
        if (md_array_size(node->children) == 0) return 0;
        node = node->type == AST_IDENTIFIER ? node->children[0] : node->children[md_array_size(node->children) - 1];
        break;
    default:
        break;
    }
    if (subidx < -1 || (subidx >= 0 && node->data.type.dim[0] >= 0 && subidx >= node->data.type.dim[0])) return 0;

    // Written straight into buf: what does not fit is dropped
    md_buf_writer_t bw = md_buf_writer_create(buf, cap);

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    if (!(node->flags & FLAG_DYNAMIC) && node->data.ptr) {
        // Known since the compilation
        summarize_data(&bw.w, node->data, subidx);
    } else {
        // Within a context ('x in y') the value is that of the first: the one the hover shows first
        eval_context_t ctx = {
            .ir = vis_ctx->ir,
            .sys = vis_ctx->sys,
            .mol_ctx = (node->type != AST_CONTEXT && node->num_contexts) ? &node->contexts[0] : NULL,
            .temp_alloc = temp_alloc,
            .alloc = temp_alloc,
            .cur_state = vis_ctx->state,
            .ref_state = &vis_ctx->sys->reference,
        };
        data_t data = {0};
        if (evaluate_node_alloc(&data, node, &ctx, temp_alloc)) {
            data.unit[0] = node->data.unit[0];
            data.unit[1] = node->data.unit[1];
            data.value_range = node->data.value_range;
            data.period = node->data.period;
            summarize_data(&bw.w, data, subidx);
            if (node->type != AST_CONTEXT && node->num_contexts > 1) {
                md_writef(&bw.w, "\n(the first of %zu contexts)", (size_t)node->num_contexts);
            }
        }
    }

    md_temp_end(temp);
    return md_buf_writer_to_str(&bw).len;
}

bool md_script_vis_eval_string(md_script_vis_t* vis, str_t expr, const md_script_vis_ctx_t* vis_ctx, md_script_vis_flags_t flags) {
    ASSERT(vis);
    ASSERT(vis_ctx);

    if (vis->magic != VIS_MAGIC) {
        MD_LOG_ERROR("Visualize: vis object not initialized.");
        return false;
    }

    if (flags == 0) flags = 0xFFFFFFFFU;

    bool success = false;

    // Not the arena of the visualization, whose results must outlive this scope
    md_temp_scope_t temp = md_temp_begin_avoid(vis->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    md_script_ir_t* ir = create_ir(temp_alloc);
    ir->str = str_copy(expr, ir->arena);

    if (vis_ctx->ir) {
        add_ir_ctx(ir, vis_ctx->ir);
    }

    if (parse_script(ir, temp_alloc)) {
        if (static_type_check(ir, vis_ctx->sys, temp_alloc)) {
            size_t num_nodes = md_array_size(ir->type_checked_expressions);
            if (num_nodes) {
                md_bitset_free(&vis->atom_mask, vis->alloc);
                md_bitmask_t atom_mask = md_bitmask_create(0, (uint32_t)vis_ctx->sys->atom.count, temp_alloc);

                eval_context_t ctx = {
                    .ir = ir,
                    .sys = vis_ctx->sys,
                    .alloc = temp_alloc,
                    .temp_alloc = temp_alloc,
                    .vis = vis,
                    .vis_flags = flags,
                    .vis_structure = &atom_mask,
                    .vis_atom_mask = &atom_mask,
                    .cur_state = vis_ctx->state,
                    .ref_state = &vis_ctx->sys->reference,
                    .eval_flags = EVAL_FLAG_FLATTEN | EVAL_FLAG_NO_STATIC_EVAL
                };

                for (size_t i = 0; i < num_nodes; ++i) {
                    ast_node_t* node = ir->type_checked_expressions[i];
                    if (node) {
                         visualize_node(node, &ctx);
                    }
                }

                // Append all structures into the 'global' atom mask
                for (size_t i = 0; i < md_array_size(vis->structures); ++i) {
                    md_bitmask_or(&atom_mask, vis->structures[i]);
                }
                vis->atom_mask = md_bitmask_extract(&atom_mask, vis->alloc);
                success = true;
            }
        }
    }

    md_temp_end(temp);
    return success;
}

void md_script_vis_init(md_script_vis_t* vis, md_allocator_i* alloc) {
    ASSERT(vis);
    ASSERT(alloc);
    MEMSET(vis, 0, sizeof(md_script_vis_t));
    vis->magic = VIS_MAGIC;
    vis->alloc = alloc;
}

bool md_script_vis_free(md_script_vis_t* vis) {
    ASSERT(vis);
    if (vis->magic == VIS_MAGIC) {
        ASSERT(vis->alloc);

        for (size_t i = 0; i < md_array_size(vis->structures); ++i) {
            md_bitset_free(&vis->structures[i], vis->alloc);
        }

        md_array_free(vis->structures, vis->alloc);

        for (size_t i = 0; i < md_array_size(vis->groups); ++i) {
            md_bitset_free(&vis->groups[i], vis->alloc);
        }
        md_array_free(vis->groups, vis->alloc);

        md_array_free(vis->points, vis->alloc);
        md_array_free(vis->lines, vis->alloc);
        md_array_free(vis->triangles, vis->alloc);
        md_array_free(vis->spheres, vis->alloc);
        md_array_free(vis->text, vis->alloc);
        md_array_free(vis->structures, vis->alloc);

        for (size_t i = 0; i < md_array_size(vis->sdf.structures); ++i) {
            md_bitset_free(&vis->sdf.structures[i], vis->alloc);
        }

        md_array_free(vis->sdf.matrices, vis->alloc);
        md_array_free(vis->sdf.structures, vis->alloc);

        md_bitset_free(&vis->atom_mask, vis->alloc);
    }
    MEMSET(vis, 0, sizeof(md_script_vis_t));
    return true;
}

bool md_script_vis_clear(md_script_vis_t* vis) {
    if (vis->magic != VIS_MAGIC) {
        MD_LOG_ERROR("Vis: Failed to clear object: not initialized");
        return false;
    }

    for (size_t i = 0; i < md_array_size(vis->structures); ++i) {
        md_bitset_free(&vis->structures[i], vis->alloc);
    }

    for (size_t i = 0; i < md_array_size(vis->sdf.structures); ++i) {
        md_bitset_free(&vis->sdf.structures[i], vis->alloc);
    }

    for (size_t i = 0; i < md_array_size(vis->groups); ++i) {
        md_bitset_free(&vis->groups[i], vis->alloc);
    }
    md_array_shrink(vis->groups, 0);

    md_array_shrink(vis->structures, 0);
    md_array_shrink(vis->points, 0);
    md_array_shrink(vis->lines, 0);
    md_array_shrink(vis->triangles, 0);
    md_array_shrink(vis->spheres, 0);
    md_array_shrink(vis->text, 0);
    md_array_shrink(vis->structures, 0);

    md_array_shrink(vis->sdf.matrices, 0);
    md_array_shrink(vis->sdf.structures, 0);

    md_bitset_free(&vis->atom_mask, vis->alloc);

    md_allocator_i* alloc = vis->alloc;
    MEMSET(vis, 0, sizeof(md_script_vis_t));
    vis->magic = VIS_MAGIC;
    vis->alloc = alloc;
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


// Keep in step with the tokenizer (tokenizer_get_next_from_buffer)
static const str_t script_keywords[] = {
    STR_INIT("and"), STR_INIT("or"), STR_INIT("xor"), STR_INIT("not"), STR_INIT("in"), STR_INIT("of"), STR_INIT("out"),
};

// Procedures that the parser handles itself instead of looking them up in the procedure table
static const str_t script_intrinsics[] = {
    STR_INIT("attr"), STR_INIT("flatten"), STR_INIT("transpose"),
};


size_t md_script_num_keywords(void) {
    return ARRAY_SIZE(script_keywords);
}

const str_t* md_script_keywords(void) {
    return script_keywords;
}

static bool procedure_name_before(size_t idx, str_t name) {
    for (size_t i = 0; i < idx; ++i) {
        if (str_eq(procedures[i].name, name)) return true;
    }
    return false;
}

size_t md_script_builtin_identifiers(str_t* out, size_t cap) {
    size_t count = 0;
    #define EMIT(name) do { if (out && count < cap) out[count] = (name); ++count; } while (0)

    for (size_t i = 0; i < ARRAY_SIZE(procedures); ++i) {
        // Overloads share a name, list it once
        if (!procedure_name_before(i, procedures[i].name)) EMIT(procedures[i].name);
    }
    for (size_t i = 0; i < ARRAY_SIZE(script_intrinsics); ++i) {
        if (!procedure_name_before(ARRAY_SIZE(procedures), script_intrinsics[i])) EMIT(script_intrinsics[i]);
    }
    for (size_t i = 0; i < ARRAY_SIZE(constants); ++i) {
        EMIT(constants[i].name);
    }

    #undef EMIT
    return count;
}

#include "md_script_complete.inl"
