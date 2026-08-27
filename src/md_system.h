#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <md_types.h>
#include <core/md_unit.h>

// Forward declare trajectory type to avoid including trajectory header here
typedef struct md_trajectory_i md_trajectory_i;

typedef struct md_atom_type_data_t {
    size_t count;

    md_label_t*     name;
    md_atomic_number_t* z;
    float*          mass;
    float*          radius;
    uint32_t*       color;
    md_flags_t*     flags;
} md_atom_type_data_t;

typedef struct md_atom_data_t {
    size_t count;

    // Coordinates live in md_system_state_t, not here. A system is time invariant; where the atoms
    // are is not.
    md_atom_type_idx_t* type_idx;
    md_flags_t* flags;

    md_atom_type_data_t type;
} md_atom_data_t;

// Component (Residue)
typedef struct md_component_data_t {
    size_t count;
    md_label_t* name;
    md_sequence_id_t* seq_id;
    uint32_t* atom_offset;
    md_flags_t* flags;
} md_component_data_t;

// Instance (e.g. Chains + other)
typedef struct md_instance_data_t {
    size_t count;
    md_label_t* id;
    md_label_t* auth_id;
    uint32_t* comp_offset;
    md_entity_idx_t* entity_idx;
} md_instance_data_t;

// Entities are used to describe the types found within a system
typedef struct md_entity_data_t {
    size_t count;
    md_label_t* id;
    md_flags_t* flags;
    str_t* description;
} md_entity_data_t;

// @TODO: Split this into two or more structucomp,
// One for nucleic backbone and one for protein backbone
typedef struct md_protein_backbone_data_t {
    // This holds the consecutive ranges which form the backbones
    struct {
        size_t count;
        uint32_t* offset; // Offsets into the segments
        md_instance_idx_t* inst_idx; // Reference to the instance in which the backbone is located
    } range;

    // These fields share the same length 'count'
    struct {
        size_t count;
        md_amino_acid_atoms_t* atoms;
        md_backbone_angles_t* angle;
        md_secondary_structure_t* secondary_structure;
        md_ramachandran_type_t* rama_type;
        md_component_idx_t* comp_idx;                  // Index to the component which contains the backbone
    } segment;
} md_protein_backbone_data_t;

typedef struct md_nucleic_backbone_data_t {
    // This holds the consecutive ranges which form the backbones
    struct {
        size_t count;
        uint32_t* offset; // Offsets into the backbone fields stored bellow
        md_instance_idx_t* inst_idx; // Reference to the instance in which the backbone is located
    } range;

    // These fields share the same length 'count'
    struct {
        size_t count;
        md_nucleic_acid_atoms_t* atoms;
        md_component_idx_t* comp_idx;                  // Index to the component which contains the backbone segment
    } segment;
} md_nucleic_backbone_data_t;

// This represents symmetries which are instanced, commonly found
// in PDB and mmcif data. It is up to the renderer to properly render this instanced data.
typedef struct md_assembly_data_t {
    size_t count;
    md_urange_t* atom_range;
    md_label_t* label;
    mat4_t* transform;
} md_assembly_data_t;

// Atom centric representation of bonds
typedef struct md_bond_conn_data_t {
    size_t count;
    md_atom_idx_t* atom_idx; // Indices to the 'other' atoms
    md_bond_idx_t* bond_idx; // Indices to the bonds
    // The offsets into the atom_idx and bond_idx for each atom.
    // Consequently offset_count should be atom count + 1
    size_t offset_count;
    uint32_t* offset;
} md_bond_conn_data_t;

// Bond centric representation
typedef struct md_bond_data_t {
    size_t count;
    md_atom_pair_t*  pairs;
    md_bond_flags_t* flags;
    md_bond_conn_data_t  conn;   // Connectivity
} md_bond_data_t;

typedef struct md_bond_iter_t {
    const md_bond_data_t* data;
	uint32_t i;
    uint32_t end_idx;
} md_bond_iter_t;

// Structure represents one connected component within the system
// And contains the iteration order of atoms and parent-child relationships for traversing the structure as a tree
// This is used for unwrapping for example.
typedef struct md_structure_t {
    const int32_t* atom_idx;
    const int32_t* parent_idx;
    size_t   count;
} md_structure_t;

// This is the contiguous storage of all structures, the offsets field is used to index into the atom_idx and parent_idx arrays for each structure
// From this individual structures (md_structure_t) can be extracted
//
// The index into atom_idx / parent_idx is referred to as a SLOT. Both arrays are indexed by slot
// and both hold GLOBAL atom indices:
//   atom_idx[slot]   - the atom
//   parent_idx[slot] - the atom it was reached from during the traversal
// The root of each structure is its own parent, so a root is parent_idx[slot] == atom_idx[slot] and
// consumers need no sentinel branch. There is exactly one root per structure.
//
// atom_slot is the reverse map, indexed by GLOBAL atom index and covering every atom in the system:
//   atom_slot[atom] -> the slot which holds that atom
// It turns "what is the parent of this atom" into an O(1) lookup, which is what lets a caller walk
// the hierarchy upwards from an arbitrary subset of atoms rather than only downwards from a root
// (see md_util_unwrap_structure).
//
// ROOT SELECTION:
//   The root is the most topologically central atom of the structure - the atom whose greatest edge
//   distance to any other atom in the structure is smallest, i.e. the graph center. This bounds the
//   depth of the traversal, which bounds both the cost of walking the hierarchy upwards and the
//   number of minimum image steps accumulated when unwrapping. It also makes the root a property of
//   the molecule rather than of the atom ordering in the source file, so two identical molecules
//   stored in different orders get the same traversal.
//
//   The center is located with the standard double sweep (farthest atom from an arbitrary seed,
//   then farthest atom from that, then the midpoint of the path between them). That is exact for
//   acyclic structures and a close approximation in the presence of rings, which are local and
//   small in practice.
//
// INVARIANT - load bearing, md_util_unwrap_structure depends on it:
//   Within a structure, atoms are stored in BFS order from the root, so every atom appears after
//   the one it was reached from: slot(parent) < slot(child). A single forward pass therefore always
//   sees a placed parent, and sorting an arbitrary set of slots ascending yields a valid
//   topological order.
//   Do not reorder atom_idx (for locality or otherwise) without reordering parent_idx and updating
//   atom_slot to match.
typedef struct md_structure_data_t {
    size_t count;
    uint32_t* offset;    // Offsets into the structure fields stored bellow, includes sentinel at the end (so length is count + 1)
    int32_t* atom_idx;   // [slot] -> global atom index
    int32_t* parent_idx; // [slot] -> global atom index of the atom it was reached from (a root is its own parent)
    int32_t* atom_slot;  // [global atom index] -> slot. Length is the system atom count
} md_structure_data_t;

typedef struct md_hydrogen_bond_candidates_t {
    struct {
        size_t count;
        md_atom_idx_t* idx;
        int* num_lone_pairs;
    } acceptor;
   
    struct {
        size_t count;
        md_atom_idx_t* d_idx;
        md_atom_idx_t* h_idx;
    } donor;
} md_hydrogen_bond_candidates_t;

typedef struct md_hydrogen_bond_pair_t {
	uint32_t don_idx; // Index into hydrogen bond candidate donors
	uint32_t acc_idx; // Index into hydrogen bond candidate acceptors
} md_hydrogen_bond_pair_t;

typedef struct md_hydrogen_bond_data_t {
    md_hydrogen_bond_candidates_t candidate;

    size_t num_bonds;
    md_hydrogen_bond_pair_t* bonds;
} md_hydrogen_bond_data_t;

// ATTRIBUTES
//
// Auxiliary data attached to a system, keyed by an HDF5 style path. The path is the whole
// key: the leading segments name the category the data applies to, the last segment names
// the field.
//
//     atom/charge/mulliken    per atom Mulliken charges
//     atom/velocity           per atom velocities
//     dipole/magnetic/origin  the magnetic dipole origin
//     dipole/magnetic/vector  the magnetic dipole vector
//
// Categories are NOT predeclared. A producer writes whatever path fits its data; a consumer
// asks for it by name, or asks for everything below a prefix with md_attributes_query.
// Prefixes match at segment boundaries, so "atom" matches "atom/charge" but never
// "atomic_number/z".
//
// The path carries the meaning, the format carries only the layout. There is deliberately
// no semantic tag on the format: atom/velocity is a vector because of what it is called,
// and a second field saying so could only ever disagree with the name.
//
// LAYOUT. The trailing axis of shape is the components of one value; the leading axes index
// whatever the category is over.
//
//     rank 0     a single value                 SCF energy
//     {N}        N scalars                      atom/charge/mulliken
//     {N,3}      N three component values       atom/position, atom/velocity
//     {N,6}      N six component values         atom/adp
//     {M,N,3}    M sets of N                    atom/normal_mode
//
// {N,3} and not {3,N}: the components of one value are contiguous, so a row can be handed
// straight to something expecting a vec3.
//
// Independent quantities over the same category are SIBLING PATHS, never an extra axis.
// Three charge schemes are atom/charge/{mulliken,hirshfeld,lowdin}, not one {N,3}
// attribute, which would be indistinguishable from a velocity.
//
// ANCHORING. A vector quantity whose category has no implicit anchor publishes its origin as a
// sibling in the same group, with the SAME shape as the vector so no consumer has to branch:
//
//     dipole/ground_state/vector        {1,3}    the moment
//     dipole/ground_state/origin        {1,3}    where to draw it from
//     dipole/electric_transition/vector {N,3}    one per excited state
//     dipole/electric_transition/origin {N,3}    one per excited state
//
// atom/velocity needs no origin: atom i is anchored at atom i's position, and that is what makes
// the per atom case work without any machinery. A dipole has no index space behind it, so its
// anchor has to be said out loud. Note that the two are commonly in DIFFERENT units, which is
// the reason unit is carried per attribute rather than per group.
//
// UNIT. Dimensionless is (md_unit_t){0}, so a unit is never invalid and a caller with nothing to
// say passes md_unit_none(). It is descriptive only: nothing here converts, and two attributes in
// the same group may legitimately disagree, as vector and origin above do.
//
// EXTENT IS NOT CHECKED. Nothing here knows that "atom/..." ought to have
// shape[0] == sys->atom.count; that is the price of not predeclaring categories. A consumer
// that indexes by atom index checks the extent itself before it does.
//
// The table is dataset scoped: it lives on md_system_t and a path is unique within one table
// only. Two datasets both holding atom/charge/mulliken is the normal case, so do not hoist
// this into anything global.
//
// Set alloc before the first create, exactly as md_system_state_t and md_index_data_t
// require. md_attributes_free clears it along with everything else.

#define MD_ATTRIBUTE_INVALID  ((md_attribute_id_t)0)
#define MD_ATTRIBUTE_MAX_RANK 4

// Hash of the full path. Zero is the invalid value, so a zeroed field is not a valid
// attribute. Stable across a reload, so a consumer may store one and re-find its attribute.
typedef uint64_t md_attribute_id_t;

typedef enum md_attribute_type_t {
    MD_ATTRIBUTE_TYPE_NONE = 0,
    MD_ATTRIBUTE_TYPE_F32,
    MD_ATTRIBUTE_TYPE_F64,
    MD_ATTRIBUTE_TYPE_I8,
    MD_ATTRIBUTE_TYPE_U8,
    MD_ATTRIBUTE_TYPE_I16,
    MD_ATTRIBUTE_TYPE_U16,
    MD_ATTRIBUTE_TYPE_I32,
    MD_ATTRIBUTE_TYPE_U32,
    MD_ATTRIBUTE_TYPE_I64,
    MD_ATTRIBUTE_TYPE_U64,
    MD_ATTRIBUTE_TYPE_COUNT,
} md_attribute_type_t;

typedef struct md_attribute_format_t {
    md_attribute_type_t type;
    uint32_t            rank;                          // 0 means a single value
    uint32_t            shape[MD_ATTRIBUTE_MAX_RANK];  // trailing axis is components per value
} md_attribute_format_t;

typedef struct md_attribute_t {
    md_attribute_id_t     id;
    str_t                 name;    // full path, owned by the table
    md_attribute_format_t format;
    md_unit_t             unit;    // (md_unit_t){0} is dimensionless, a valid and common value
    void*                 data;    // owned by the table, md_attribute_byte_size() bytes, zeroed on create
} md_attribute_t;

typedef struct md_attributes_t {
    struct md_allocator_i*   alloc;
    md_array(md_attribute_t) attr;  // kept sorted by name
} md_attributes_t;

// A snapshot of the geometric state of a system: where the atoms are and what box they are in.
//
// This holds exactly the fields which share one interpolation contract - same type in and out,
// periodic boundary aware, handled as a unit by md_util_interpolate_*. Quantities which vary over
// time but do not obey that contract (per atom charge, backbone angles, secondary structure) do
// not belong here; they interpolate differently and belong with whoever produces them.
//
// Two presence bits, both self describing:
//   num_atoms == 0        -> no coordinates
//   unitcell.flags == 0   -> no cell
// num_atoms is used rather than testing x != NULL because md_array_ensure allocates capacity
// without setting size, so a non NULL x does not imply the coordinates were populated.
// Ownership: alloc non-NULL means this state owns x/y/z and must be freed with
// md_system_state_free. alloc NULL means the state is a non owning view over coordinates somebody
// else owns - a scratch arena during script evaluation, a GPU mapped buffer during interpolation.
// Both forms are load bearing, and this is the only thing that distinguishes them.
//
// Set alloc before handing the state to a producer, exactly as md_system_t requires sys->alloc to
// be set before loading. The two then read symmetrically at the call site, and the state's lifetime
// is free to differ from the system's - a temp allocator for the state, a persistent one for the
// system, is a legitimate and useful combination.
// frame is the trajectory ordinal the coordinates came from, as a continuous quantity: the integer
// part selects the frame, the fractional part is how far between that frame and the next the state
// has been interpolated. A NEGATIVE value means the state did not come from a trajectory at all
// (a topology's own coordinates, a scratch buffer), which is why 0.0 cannot serve as that marker.
// Use md_state_has_frame / md_state_frame_floor / md_state_frame_nearest / md_state_frame_frac
// rather than reading the field: a raw cast of the absent value lands on frame 0, which is in range
// and plausible and therefore the kind of mistake that survives.
//
// CAVEAT, unlike the two presence bits above: -1 does not survive zero initialisation. A state
// built as {0} or with designated initialisers that omit frame reads as frame 0, not as absent.
// md_system_state_init stamps -1, so any state that went through it reads as absent until something
// writes a real frame; a hand rolled {0} does not. Only md_trajectory_reader_load_frame and the
// interpolation which produces a state write a non negative value.
//
// @NOTE: physical time is deliberately NOT stored here. It is derivable from frame and the
// trajectory's frame times (md_trajectory_time_at_frame), and storing both would reintroduce the
// very thing this struct exists to prevent - two fields which must agree, with nothing enforcing it.
typedef struct md_system_state_t {
    size_t num_atoms;
    float* x;
    float* y;
    float* z;
    md_unitcell_t unitcell;
    double frame;
    md_allocator_i* alloc;
} md_system_state_t;

// This represents the persistent portion (topology) of a system which does not change over time, such as the atom types, bonds, components, etc.
// It may of course be modified through some special operations, though it is not expected to change frequently.
typedef struct md_system_t {
    md_allocator_i*             alloc;
    md_trajectory_i*            trajectory;

    // The state from which the derived topology below (bonds, rings, structures, backbones) was
    // inferred. Written by md_util_system_infer as part of performing the inference, so it is by
    // construction the input which produced that topology and cannot go stale.
    //
    // RULE: only inference and load time completion read this. Every other operation takes the
    // state it works on as an explicit parameter. Reaching for sys->reference to avoid threading
    // a state through is how the previous coupling arose.
    md_system_state_t           reference;

    md_atom_data_t              atom;
    md_component_data_t         component;
    md_instance_data_t          instance;
    md_entity_data_t            entity;

    md_protein_backbone_data_t  protein_backbone;
    md_nucleic_backbone_data_t  nucleic_backbone;
    
    md_bond_data_t              bond;               // Persistent covalent bonds
    md_hydrogen_bond_data_t     hydrogen_bond;      // Hydrogen bonds
    
    md_index_data_t             ring;               // Ring structures formed by persistent bonds
    md_structure_data_t         structure;          // Isolated structures connected by persistent bonds

    md_assembly_data_t          assembly;           // Assemblies of  (duplications of ranges with new transforms)
    
    md_attributes_t             attributes;         // Custom attributes

    str_t                       description;
} md_system_t;

#ifdef __cplusplus
extern "C" {
#endif

// ATTRIBUTE FORMAT
// These operate on the format alone, so a producer can size a buffer before it creates anything.

// Bytes of one component. Zero for MD_ATTRIBUTE_TYPE_NONE.
size_t md_attribute_type_size(md_attribute_type_t type);

// Components of one value: the trailing axis, or 1 when rank < 2.
size_t md_attribute_components(const md_attribute_format_t* format);

// Number of values: the product of the leading axes, or 1 when rank < 2.
size_t md_attribute_value_count(const md_attribute_format_t* format);

// Total components: value_count * components == product(shape), or 1 when rank == 0.
size_t md_attribute_element_count(const md_attribute_format_t* format);

// element_count * type_size.
size_t md_attribute_byte_size(const md_attribute_format_t* format);

// ATTRIBUTE PATH
// Views into the attribute's own name, nothing is allocated.

// "atom/charge/mulliken" -> "atom/charge". Empty if the path has no separator.
str_t md_attribute_category(const md_attribute_t* attr);

// "atom/charge/mulliken" -> "mulliken". The whole path if it has no separator.
str_t md_attribute_field(const md_attribute_t* attr);

// Extracts the whole attribute into dst, converting the stored type to f32 and the stored unit to
// dst_unit. Returns the number of floats written, 0 on failure. cap is in floats and must be at
// least md_attribute_element_count().
//
// dst_unit of md_unit_none() means "as stored": no conversion is attempted, whatever the attribute
// carries. Otherwise the two units must share dimensions or the call FAILS - asking for e a0 in
// Angstrom is a refusal, not a rescale, because an arrow drawn at 1.9x the right length looks
// entirely plausible on screen and nothing downstream will catch it.
//
// The conversion is applied in double and narrowed once, at the end. Narrowing first loses
// precision against the small factors (1e-10 for Angstrom in SI) that make up most real conversions.
//
// Integral attributes convert numerically, which is right for a charge stored as an integer and
// meaningless for an index or a label column. Nothing here can tell those apart; the path is what
// tells the consumer, as designed.
size_t md_attribute_extract_f32(float dst[], size_t cap, const md_attribute_t* attr, md_unit_t dst_unit);

// ATTRIBUTE TABLE

void md_attributes_free(md_attributes_t* attributes);

size_t md_attributes_count(const md_attributes_t* attributes);

// Creates an attribute and returns its id, MD_ATTRIBUTE_INVALID on failure.
//
// data NULL reserves the storage and zeroes it, for a producer that computes straight into the
// attribute through md_attributes_data; byte_size must then be 0. Otherwise byte_size must equal
// md_attribute_byte_size(&format) exactly and the buffer is copied in.
//
// The size is not decoration. Nothing here can inspect what data points at, so it is the only
// guard against a buffer whose element type or count disagrees with the declared format: handing
// a float array to a format declaring F64 reads twice past the end and the size catches it. A
// mismatch is a rejection, never a truncated copy.
//
// Populating here rather than afterwards also keeps creation atomic. A create followed by a
// separate write leaves a fully zeroed, valid looking attribute in the table in between, which
// any early return in that window makes permanent.
//
// Rejects an unset allocator, an empty path, a leading or trailing '/', an empty segment, a
// duplicate path, a type of NONE, a rank above MD_ATTRIBUTE_MAX_RANK, a zero extent and a
// byte_size disagreeing with the format.
md_attribute_id_t md_attributes_create(md_attributes_t* attributes, str_t name, md_attribute_format_t format, md_unit_t unit, const void* data, size_t byte_size);

bool md_attributes_remove(md_attributes_t* attributes, md_attribute_id_t id);

// NULL if absent.
// INVALIDATED BY THE NEXT create OR remove: the array reallocates and is kept sorted by name.
// Hold the id, not the pointer.
const md_attribute_t* md_attributes_get (const md_attributes_t* attributes, md_attribute_id_t id);
const md_attribute_t* md_attributes_find(const md_attributes_t* attributes, str_t name);

// Writable view for whoever fills the attribute in. NULL if the id is unknown or the stored
// type is not expected_type, so nobody memcpys f64 into an f32 attribute and finds out at
// render time. Same invalidation rule as above.

void* md_attributes_data(md_attributes_t* attributes, md_attribute_id_t id, md_attribute_type_t expected_type);

// Ids of every attribute at or below prefix, in name order. A prefix matches at segment
// boundaries and an optional trailing '/' is ignored, so "atom" and "atom/" both match
// "atom/charge" and "atom/charge/mulliken", and neither matches "atomic_number/z". An empty
// prefix matches everything.
// Returns the total number of matches and writes at most cap of them, so pass cap 0 to count.
size_t md_attributes_query(md_attribute_id_t out_ids[], size_t cap, const md_attributes_t* attributes, str_t prefix);

// Immediate child segments below prefix, deduplicated, in name order: "atom/" yields
// "charge", "velocity". For walking the paths as a tree. The returned strings are views into
// the stored names and follow the same invalidation rule.
// Returns the total number of children and writes at most cap of them.
size_t md_attributes_query_children(str_t out_names[], size_t cap, const md_attributes_t* attributes, str_t prefix);

// Atom type table helper functions
static inline size_t md_atom_type_count(const md_atom_type_data_t* atom_type) {
    ASSERT(atom_type);
    return atom_type->count;
}

static inline md_atom_type_idx_t md_atom_type_find(const md_atom_type_data_t* atom_type, str_t name, md_atomic_number_t z) {
    ASSERT(atom_type);
	md_atom_type_idx_t type_idx = 0; // Zero is sentinel for "not found"
    for (size_t i = 0; i < atom_type->count; ++i) {
        str_t atom_type_name = LBL_TO_STR(atom_type->name[i]);
        if (str_eq(atom_type_name, name) && atom_type->z[i] == z) {
            type_idx = (md_atom_type_idx_t)i;
            break;
        }
    }
    return type_idx;
}

static inline md_atom_type_idx_t md_atom_type_find_or_add(md_atom_type_data_t* atom_type, str_t name, md_atomic_number_t z, float mass, float radius, uint32_t color, md_flags_t flags, struct md_allocator_i* alloc) {
    ASSERT(atom_type);
    ASSERT(alloc);
    
    // First try to find existing atom type
    md_atom_type_idx_t type_idx = md_atom_type_find(atom_type, name, z);
    if (type_idx != 0) {
        return type_idx;
    }
    
    // Add new atom type
    md_array_push(atom_type->name, make_label(name), alloc);
    md_array_push(atom_type->z, z, alloc);
    md_array_push(atom_type->mass, mass, alloc);
    md_array_push(atom_type->radius, radius, alloc);
    md_array_push(atom_type->color, color, alloc);
    md_array_push(atom_type->flags, flags, alloc);
    atom_type->count++;
    
    return (md_atom_type_idx_t)(atom_type->count - 1);
}

static inline md_atomic_number_t md_atom_type_atomic_number(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->z[type_idx];
    }
    return 0;
}

static inline float md_atom_type_mass(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->mass[type_idx];
    }
    return 0;
}

static inline md_flags_t md_atom_type_flags(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->flags[type_idx];
    }
    return MD_FLAG_NONE;
}

static inline float md_atom_type_radius(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->radius[type_idx];
    }
    return 0;
}

static inline uint32_t md_atom_type_color(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->color[type_idx];
    }
    return 0;
}

static inline str_t md_atom_type_name(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return LBL_TO_STR(type_data->name[type_idx]);
    }
    return STR_LIT("");
}

// State helpers

// A state which did not come from a trajectory carries a negative frame.
// @NOTE: NaN would pass this test, so nothing may ever write one into the field.
static inline bool md_state_has_frame(const md_system_state_t* state) {
    ASSERT(state);
    return state->frame >= 0.0;
}

// The frame the state sits on or just after: the one to interpolate FORWARD from.
// @NOTE: a cast truncates toward zero, which is floor for the non negative values md_state_has_frame
// guarantees. Deliberately not floor() - keeping <math.h> out of this header matters, see below.
static inline int64_t md_state_frame_floor(const md_system_state_t* state) {
    ASSERT(state);
    ASSERT(md_state_has_frame(state));
    return (int64_t)state->frame;
}

// The frame the state most closely corresponds to. Use this where a single frame must be named -
// indexing a per frame array, for instance - not md_state_frame_floor, which would truncate 3.99 to 3.
static inline int64_t md_state_frame_nearest(const md_system_state_t* state) {
    ASSERT(state);
    ASSERT(md_state_has_frame(state));
    return (int64_t)(state->frame + 0.5);
}

// How far the state lies between md_state_frame_floor() and the frame after it, in [0,1).
static inline double md_state_frame_frac(const md_system_state_t* state) {
    ASSERT(state);
    ASSERT(md_state_has_frame(state));
    return state->frame - (double)(int64_t)state->frame;
}

// @NOTE: do NOT add <math.h> to this header. simde-math.h keys off HUGE_VAL to detect "a math header
// was already included", and when that fires under C++ it assumes the header was <cmath> and starts
// emitting std::trunc. On glibc that is harmless (math.h defines an isnan MACRO, which simde checks
// for first) and on gcc/clang __has_builtin wins before it matters - but MSVC has neither, so any
// C++ translation unit that reaches simde through this header fails with "trunc is not a member of
// std". md_system.h is included nearly everywhere, so it is the worst possible place to trip that.

// Atom helpers

static inline md_atom_type_idx_t md_atom_type_idx(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    if (atom->type_idx && atom_idx < atom->count) {
        return atom->type_idx[atom_idx];
	}
    return 0;
}

static inline md_atomic_number_t md_atom_atomic_number(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    
    // Try atom type table first if type_idx is available
    if (atom->type_idx && (size_t)atom->type_idx[atom_idx] < atom->type.count) {
        return atom->type.z[atom->type_idx[atom_idx]];
    }
    
    return 0;
}

static inline size_t md_atom_count(const md_atom_data_t* atom_data) {
    ASSERT(atom_data);
    return atom_data->count;
}

static inline vec3_t md_state_coord(const md_system_state_t* state, size_t atom_idx) {
    ASSERT(state);
    if (atom_idx < state->num_atoms) {
        return vec3_set(state->x[atom_idx], state->y[atom_idx], state->z[atom_idx]);
    }
    return vec3_zero();
}

static inline float md_atom_mass(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    
    if (atom_idx < atom->count) {
		return md_atom_type_mass(&atom->type, atom->type_idx[atom_idx]);
    }
    
    return 0.0f;
}

static inline float md_atom_radius(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    
    if (atom_idx < atom->count) {
		return md_atom_type_radius(&atom->type, atom->type_idx[atom_idx]);
    }
    
    return 0.0f;
}

static inline str_t md_atom_name(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    if (atom_idx < atom->count) {
		return md_atom_type_name(&atom->type, atom->type_idx[atom_idx]);
    }

    return STR_LIT("");
}

static inline md_flags_t md_atom_flags(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    if (atom_idx < atom->count && atom->flags) {
        return atom->flags[atom_idx];
	}
    return MD_FLAG_NONE;
}

// Component

static inline size_t md_component_count(const md_component_data_t* comp) {
    ASSERT(comp);
    return comp->count;
}

static inline str_t md_component_name(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
    str_t name = STR_LIT("");
    if (comp->name && comp_idx < comp->count) {
        name = LBL_TO_STR(comp->name[comp_idx]);
    }
    return name;
}

static inline md_sequence_id_t md_component_seq_id(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
    md_sequence_id_t id = 0;
    if (comp->seq_id && comp_idx < comp->count) {
        id = comp->seq_id[comp_idx];
    }
    return id;
}

static inline md_urange_t md_component_atom_range(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
	md_urange_t range = {0};
	if (comp->atom_offset && comp_idx < comp->count) {
		range.beg = comp->atom_offset[comp_idx];
		range.end = comp->atom_offset[comp_idx + 1];
	}
	return range;
}

static inline md_flags_t md_component_flags(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
    md_flags_t flags = MD_FLAG_NONE;
    if (comp->flags && comp_idx < comp->count) {
        flags = comp->flags[comp_idx];
    }
    return flags;
}

static inline md_component_idx_t md_component_find_by_atom_idx(const md_component_data_t* comp, size_t atom_idx) {
    ASSERT(comp);

    md_component_idx_t comp_idx = -1;
    if (comp->atom_offset) {
        for (size_t i = 0; i < comp->count; ++i) {
            size_t comp_beg = comp->atom_offset[i];
            size_t comp_end = comp->atom_offset[i + 1];
            if (comp_beg <= atom_idx && atom_idx < comp_end) {
                comp_idx = (md_component_idx_t)i;
                break;
            }
            if (comp_beg > atom_idx) {
                break;
            }
        }
    }

    return comp_idx;
}

static inline size_t md_component_atom_count(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
    size_t count = 0;

    if (comp->atom_offset && comp_idx < comp->count) {
        count = comp->atom_offset[comp_idx + 1] - comp->atom_offset[comp_idx];
    }
    return count;
}

// Instance

static inline size_t md_instance_count(const md_instance_data_t* inst) {
    ASSERT(inst);
    return inst->count;
}

static inline md_urange_t md_instance_component_range(const md_instance_data_t* inst, size_t inst_idx) {
    ASSERT(inst);

    md_urange_t range = {0};
    if (inst->comp_offset && inst_idx < inst->count) {
        range.beg = inst->comp_offset[inst_idx];
        range.end = inst->comp_offset[inst_idx + 1];
    }
    return range;
}

static inline md_instance_idx_t md_instance_find_by_comp_idx(const md_instance_data_t* inst, size_t comp_idx) {
    ASSERT(inst);

    md_instance_idx_t inst_idx = -1;
    if (inst->comp_offset) {
        for (size_t i = 0; i < inst->count; ++i) {
            size_t inst_beg = inst->comp_offset[i];
            size_t inst_end = inst->comp_offset[i + 1];
            if (inst_beg <= comp_idx && comp_idx < inst_end) {
                inst_idx = (md_instance_idx_t)i;
                break;
            }
            if (inst_beg > comp_idx) {
                break;
            }
        }
    }
    return inst_idx;
}

/*
static inline md_instance_idx_t md_inst_find_by_atom_idx(const md_instance_data_t* inst, size_t atom_idx) {
    ASSERT(inst);

    md_instance_idx_t inst_idx = -1;
    if (inst->atom_range) {
        int ai = (int)atom_idx;
        for (size_t i = 0; i < inst->count; ++i) {
            md_range_t range = inst->atom_range[i];
            if (range.beg <= ai && ai < range.end) {
                inst_idx = (md_instance_idx_t)i;
                break;
            }
            if (range.beg > ai) {
                break;
            }
        }
    }
    return inst_idx;
}
*/

static inline size_t md_instance_comp_count(const md_instance_data_t* inst, size_t inst_idx) {
    ASSERT(inst);

    size_t count = 0;
    if (inst->comp_offset && inst_idx < inst->count) {
        count = inst->comp_offset[inst_idx + 1] - inst->comp_offset[inst_idx];
    }
    return count;
}

/*
static inline md_urange_t md_inst_atom_range(const md_instance_data_t* inst, size_t inst_idx) {
	ASSERT(inst);

    md_range_t range = {0};
    if (inst->comp_range && inst_idx < inst->count) {
        uint32_t cbeg = inst->comp_range[inst_idx].beg;
        uint32_t cend = inst->comp_range[inst_idx].end;
        if (inst->atom_range && cbeg < cend) {
            range.beg = inst->atom_range[cbeg].beg;
            range.end = inst->atom_range[cend - 1].end;
        }
        range = inst->atom_range[inst_idx];
    }
    return range;
}

static inline size_t md_inst_atom_count(const md_instance_data_t* inst, size_t inst_idx) {
    size_t count = 0;
    if (inst->atom_range && inst_idx < inst->count) {
        md_range_t range = inst->atom_range[inst_idx];
        count = range.end - range.beg;
    }
    return count;
}
*/

static inline str_t md_instance_id(const md_instance_data_t* inst, size_t inst_idx) {
    ASSERT(inst);
    str_t id = STR_LIT("");
    if (inst->id && inst_idx < inst->count) {
        id = LBL_TO_STR(inst->id[inst_idx]);
    }
    return id;
}

static inline str_t md_instance_auth_id(const md_instance_data_t* inst, size_t inst_idx) {
    ASSERT(inst);
    str_t auth_id = STR_LIT("");
    if (inst->auth_id && inst_idx < inst->count) {
        auth_id = LBL_TO_STR(inst->auth_id[inst_idx]);
    }
    return auth_id;
}

static inline md_entity_idx_t md_instance_entity_idx(const md_instance_data_t* inst, size_t inst_idx) {
    ASSERT(inst);
    md_entity_idx_t entity_idx = -1;
    if (inst->entity_idx && inst_idx < inst->count) {
        entity_idx = inst->entity_idx[inst_idx];
    }
    return entity_idx;
}

static inline size_t md_entity_count(const md_entity_data_t* entity) {
    ASSERT(entity);
    return entity->count;
}

static inline md_entity_idx_t md_entity_find_by_id(const md_entity_data_t* entity, str_t id) {
    ASSERT(entity);
    md_entity_idx_t entity_idx = -1;
    if (entity->id) {
        for (size_t i = 0; i < entity->count; ++i) {
            str_t entity_id = LBL_TO_STR(entity->id[i]);
            if (str_eq(entity_id, id)) {
                entity_idx = (md_entity_idx_t)i;
                break;
            }
        }
    }
    return entity_idx;
}

static inline str_t md_entity_id(const md_entity_data_t* entity, size_t entity_idx) {
    ASSERT(entity);
    str_t label = STR_LIT("");
    if (entity->id && entity_idx < entity->count) {
        label = LBL_TO_STR(entity->id[entity_idx]);
    }
    return label;
}

static inline str_t md_entity_description(const md_entity_data_t* entity, size_t entity_idx) {
    ASSERT(entity);
    str_t desc = STR_LIT("");
    if (entity->description && entity_idx < entity->count) {
        desc = entity->description[entity_idx];
    }
    return desc;
}

static inline md_flags_t md_entity_flags(const md_entity_data_t* entity, size_t entity_idx) {
    ASSERT(entity);
    md_flags_t flags = MD_FLAG_NONE;
    if (entity->flags && entity_idx < entity->count) {
        flags = entity->flags[entity_idx];
    }
    return flags;
}

static inline size_t md_structure_count(const md_structure_data_t* structure) {
    ASSERT(structure);
    return structure->count;
}

static inline bool md_structure_extract(md_structure_t* out_structure, const md_structure_data_t* structure_data, size_t struct_idx) {
    ASSERT(out_structure);
    ASSERT(structure_data);
    if (struct_idx < structure_data->count) {
        uint32_t offset = structure_data->offset[struct_idx];
        uint32_t next_offset = structure_data->offset[struct_idx + 1];
        out_structure->count = next_offset - offset;
        out_structure->atom_idx = &structure_data->atom_idx[offset];
        out_structure->parent_idx = &structure_data->parent_idx[offset];
        return true;
    }
    return false;
}

// Slot of an atom within the flat structure arrays.
// Requires md_util_system_infer_structures to have been run.
static inline int32_t md_structure_atom_slot(const md_structure_data_t* structure_data, int32_t atom_idx) {
    ASSERT(structure_data);
    ASSERT(structure_data->atom_slot);
    return structure_data->atom_slot[atom_idx];
}

// Global index of the atom the supplied atom was reached from during the traversal.
// A root atom is its own parent, so parent == atom_idx identifies a root.
// Requires md_util_system_infer_structures to have been run.
static inline int32_t md_structure_atom_parent(const md_structure_data_t* structure_data, int32_t atom_idx) {
    ASSERT(structure_data);
    ASSERT(structure_data->atom_slot);
    return structure_data->parent_idx[structure_data->atom_slot[atom_idx]];
}

// SYSTEM
// System level convenience accessors

static inline size_t md_system_atom_count(const md_system_t* sys) {
    ASSERT(sys);
    return md_atom_count(&sys->atom);
}

static inline md_flags_t md_system_atom_flags(const md_system_t* sys, size_t atom_idx) {
    ASSERT(sys);
    return md_atom_flags(&sys->atom, atom_idx);
}

static inline size_t md_system_atom_type_count(const md_system_t* sys) {
    ASSERT(sys);
    return md_atom_type_count(&sys->atom.type);
}

static inline md_flags_t md_system_atom_type_flags(const md_system_t* sys, size_t type_idx) {
    ASSERT(sys);
    return md_atom_type_flags(&sys->atom.type, type_idx);
}

static inline size_t md_system_component_count(const md_system_t* sys) {
    ASSERT(sys);
    return md_component_count(&sys->component);
}

static inline md_flags_t md_system_component_flags(const md_system_t* sys, size_t comp_idx) {
    ASSERT(sys);
    return md_component_flags(&sys->component, comp_idx);
}

static inline size_t md_system_instance_count(const md_system_t* sys) {
    ASSERT(sys);
    return md_instance_count(&sys->instance);
}

static inline size_t md_system_bond_count(const md_system_t* sys) {
    ASSERT(sys);
    return sys->bond.count;
}

static inline size_t md_system_entity_count(const md_system_t* sys) {
    ASSERT(sys);
    return sys->entity.count;
}

static inline md_flags_t md_system_entity_flags(const md_system_t* sys, size_t ent_idx) {
    ASSERT(sys);
    return md_entity_flags(&sys->entity, ent_idx);
}

static inline md_flags_t md_system_instance_flags(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    if (sys->instance.entity_idx && inst_idx < sys->instance.count) {
        return md_entity_flags(&sys->entity, sys->instance.entity_idx[inst_idx]);
    }
    return MD_FLAG_NONE;
}

static inline str_t md_system_instance_id(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    str_t id = STR_LIT("");
    if (sys->instance.id && inst_idx < sys->instance.count) {
        id = md_instance_id(&sys->instance, inst_idx);
    }
    return id;
}

static inline str_t md_system_instance_auth_id(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    str_t id = STR_LIT("");
    if (sys->instance.id && inst_idx < sys->instance.count) {
        id = md_instance_auth_id(&sys->instance, inst_idx);
    }
    return id;
}

static inline size_t md_system_instance_comp_count(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    return md_instance_comp_count(&sys->instance, inst_idx);
}

static inline md_urange_t md_system_instance_comp_range(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    return md_instance_component_range(&sys->instance, inst_idx);
}

static inline md_urange_t md_system_instance_atom_range(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    md_urange_t range = {0};
    if (inst_idx < sys->instance.count) {
        md_urange_t comp_range = md_instance_component_range(&sys->instance, inst_idx);
        if (comp_range.beg != comp_range.end) {
            range.beg = md_component_atom_range(&sys->component, comp_range.beg).beg;
            range.end = md_component_atom_range(&sys->component, comp_range.end - 1).end;
        }
    }
    return range;
}

static inline size_t md_system_instance_atom_count(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    md_urange_t atom_range = md_system_instance_atom_range(sys, inst_idx);
    return atom_range.end - atom_range.beg;
}

static inline size_t md_system_instance_entity_idx(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    return md_instance_entity_idx(&sys->instance, inst_idx);
}

static inline md_sequence_id_t md_system_component_seq_id(const md_system_t* sys, size_t comp_idx) {
    ASSERT(sys);
    return md_component_seq_id(&sys->component, comp_idx);
}

static inline md_urange_t md_system_component_atom_range(const md_system_t* sys, size_t comp_idx) {
    ASSERT(sys);
    return md_component_atom_range(&sys->component, comp_idx);
}

static inline size_t md_system_component_atom_count(const md_system_t* sys, size_t comp_idx) {
    ASSERT(sys);
    return md_component_atom_count(&sys->component, comp_idx);
}

static inline md_component_idx_t md_system_component_find_by_atom_idx(const md_system_t* sys, size_t atom_idx) {
    ASSERT(sys);
    return md_component_find_by_atom_idx(&sys->component, atom_idx);
}

static inline md_instance_idx_t md_system_instance_find_by_atom_idx(const md_system_t* sys, size_t atom_idx) {
    ASSERT(sys);
    md_instance_idx_t inst_idx = -1;
    md_component_idx_t comp_idx = md_system_component_find_by_atom_idx(sys, atom_idx);
    if (comp_idx >= 0) {
        inst_idx = md_instance_find_by_comp_idx(&sys->instance, comp_idx);
    }
    return inst_idx;
}

// Convenience functions to extract atom properties into arrays
static inline void md_atom_extract_radii(float out_radii[], size_t offset, size_t length, const md_atom_data_t* atom_data) {
    ASSERT(out_radii);
    ASSERT(atom_data);
    ASSERT(offset + length <= atom_data->count);
    
    for (size_t i = 0; i < length; ++i) {
        out_radii[i] = md_atom_radius(atom_data, offset + i);
    }
}

static inline void md_atom_extract_masses(float out_masses[], size_t offset, size_t length, const md_atom_data_t* atom_data) {
    ASSERT(out_masses);
    ASSERT(atom_data);
    ASSERT(offset + length <= atom_data->count);
    
    for (size_t i = 0; i < length; ++i) {
        out_masses[i] = md_atom_mass(atom_data, offset + i);
    }
}

static inline void md_atom_extract_atomic_numbers(md_atomic_number_t out_z[], size_t offset, size_t length, const md_atom_data_t* atom_data) {
    ASSERT(out_z);
    ASSERT(atom_data);
    ASSERT(offset + length <= atom_data->count);
    
    for (size_t i = 0; i < length; ++i) {
        out_z[i] = md_atom_atomic_number(atom_data, offset + i);
    }
}

static inline md_bond_iter_t md_bond_iter(const md_bond_data_t* bond_data, size_t atom_idx) {
    md_bond_iter_t it = {0};
    if (bond_data && bond_data->conn.offset && atom_idx < bond_data->conn.offset_count) {
        it.data = bond_data;
        it.i = bond_data->conn.offset[atom_idx];
		it.end_idx = bond_data->conn.offset[atom_idx + 1];
    }
    return it;
}

// This is not something which should be done frequently
void md_bond_build_connectivity(md_bond_data_t* in_out_bond, size_t atom_count, md_allocator_i* alloc);
void md_system_bond_build_connectivity(md_system_t* sys);

static inline size_t md_bond_conn_count(const md_bond_data_t* bond_data, size_t atom_idx) {
    ASSERT(bond_data);
    return bond_data->conn.offset[atom_idx + 1] - bond_data->conn.offset[atom_idx];
}

static inline md_atom_idx_t md_bond_conn_atom_idx(const md_bond_data_t* bond_data, uint32_t atom_conn_idx, uint32_t idx) {
    ASSERT(bond_data);
    return bond_data->conn.atom_idx[atom_conn_idx + idx];
}

static inline md_bond_idx_t md_bond_conn_bond_idx(const md_bond_data_t* bond_data, uint32_t atom_conn_idx, uint32_t idx) {
    ASSERT(bond_data);
    return bond_data->conn.bond_idx[atom_conn_idx + idx];
}

static inline bool md_bond_iter_has_next(const md_bond_iter_t* it) {
    ASSERT(it);
    return it->i < it->end_idx;
}

static inline void md_bond_iter_next(md_bond_iter_t* it) {
    ASSERT(it);
    it->i += 1;
}

static inline md_atom_idx_t md_bond_iter_atom_index(const md_bond_iter_t* it) {
    ASSERT(it);
	return it->data->conn.atom_idx[it->i];
}

static inline md_atom_idx_t md_bond_iter_bond_index(const md_bond_iter_t* it) {
    ASSERT(it);
    return it->data->conn.bond_idx[it->i];
}

static inline uint32_t md_bond_iter_bond_flags(const md_bond_iter_t* it) {
    ASSERT(it);
    return it->data->flags[it->data->conn.bond_idx[it->i]];
}

static inline md_bond_idx_t md_bond_find(const md_bond_data_t* bond_data, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b) {
    ASSERT(bond_data);
	md_bond_iter_t it = md_bond_iter(bond_data, atom_idx_a);
    while (md_bond_iter_has_next(&it)) {
        md_atom_idx_t other_atom_idx = md_bond_iter_atom_index(&it);
        if (other_atom_idx == atom_idx_b) {
            return md_bond_iter_bond_index(&it);
        }
        md_bond_iter_next(&it);
	}
	return -1;
}

static inline void md_bond_insert(md_bond_data_t* bond_data, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b, md_bond_flags_t flags, md_allocator_i* alloc) {
    ASSERT(bond_data);
    ASSERT(alloc);

    // Ensure that the bond does not exist
	md_bond_idx_t bond_idx = md_bond_find(bond_data, atom_idx_a, atom_idx_b);
    if (bond_idx != -1) {
        return;
    }

    // Add bond
    md_atom_pair_t pair = {atom_idx_a, atom_idx_b};
    md_array_push(bond_data->pairs, pair, alloc);
    md_array_push(bond_data->flags, flags, alloc);
    bond_data->count++;
}

// This requires a rebuild of connectivity to be valid again
static inline void md_bond_remove(md_bond_data_t* bond_data, md_bond_idx_t bond_idx) {
    ASSERT(bond_data);
    ASSERT(bond_idx < (md_bond_idx_t)bond_data->count);

    // Swap and pop
    size_t last_idx = bond_data->count - 1;
    if ((size_t)bond_idx != last_idx) {
        bond_data->pairs[bond_idx] = bond_data->pairs[last_idx];
        bond_data->flags[bond_idx] = bond_data->flags[last_idx];
    }
    bond_data->count--;
    md_array_shrink(bond_data->pairs, bond_data->count);
    md_array_shrink(bond_data->flags, bond_data->count);
}

static inline md_atom_pair_t md_bond_pair(const md_bond_data_t* bond_data, md_bond_idx_t bond_idx) {
    ASSERT(bond_data);
    if (bond_idx < (md_bond_idx_t)bond_data->count) {
        return bond_data->pairs[bond_idx];
    }
    md_atom_pair_t invalid_pair = { -1, -1 };
    return invalid_pair;
}

static inline md_bond_idx_t md_system_bond_find(const md_system_t* sys, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b) {
    ASSERT(sys);
    return md_bond_find(&sys->bond, atom_idx_a, atom_idx_b);
}

static inline void md_system_bond_insert(md_system_t* sys, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b, md_bond_flags_t flags) {
    ASSERT(sys);
    md_bond_insert(&sys->bond, atom_idx_a,  atom_idx_b, flags, sys->alloc);
}

static inline void md_system_bond_remove(md_system_t* sys, md_bond_idx_t bond_idx) {
    ASSERT(sys);
    md_bond_remove(&sys->bond, bond_idx);
}

static inline md_bond_flags_t md_system_bond_flags(const md_system_t* sys, md_bond_idx_t bond_idx) {
    ASSERT(sys);
    if (bond_idx < (md_bond_idx_t)sys->bond.count) {
        return sys->bond.flags[bond_idx];
    }
    return MD_BOND_FLAG_NONE;
}

static inline void md_bond_conn_clear(md_bond_conn_data_t* conn_data) {
    ASSERT(conn_data);
    conn_data->count = 0;
    md_array_shrink(conn_data->atom_idx, 0);
    md_array_shrink(conn_data->bond_idx, 0);

    conn_data->offset_count = 0;
    md_array_shrink(conn_data->offset, 0);
}

static inline void md_bond_data_clear(md_bond_data_t* bond_data) {
    ASSERT(bond_data);

    bond_data->count = 0;
    md_array_shrink(bond_data->pairs, 0);
    md_array_shrink(bond_data->flags, 0);
    
    bond_data->conn.count = 0;
    md_array_shrink(bond_data->conn.atom_idx, 0);
    md_array_shrink(bond_data->conn.bond_idx, 0);

    md_bond_conn_clear(&bond_data->conn);
}

static inline bool md_atom_is_connected_to_atomic_numbers(const md_atom_data_t* atom_data, const md_bond_data_t* bond_data, size_t atom_idx, const md_atomic_number_t z_list[], size_t z_count) {
    ASSERT(bond_data);
    ASSERT(atom_data);
    ASSERT(atom_idx < atom_data->count);
    bool found = false;
    md_bond_iter_t it = md_bond_iter(bond_data, atom_idx);
    while (md_bond_iter_has_next(&it) && !found) {
        md_atom_idx_t other_atom_idx = md_bond_iter_atom_index(&it);
        md_atomic_number_t other_z = md_atom_atomic_number(atom_data, other_atom_idx);
        for (size_t i = 0; i < z_count; ++i) {
            if (other_z == z_list[i]) {
                found = true;
                break;
            }
        }
        md_bond_iter_next(&it);
    }
    return found;
}

static inline md_atom_idx_t md_hydrogen_bond_donor_atom_idx(const md_hydrogen_bond_candidates_t* hbond_cand, size_t don_idx) {
	ASSERT(hbond_cand);
    if (don_idx < hbond_cand->donor.count) {
        return hbond_cand->donor.d_idx[don_idx];
	}
    return -1;
}

static inline md_atom_idx_t md_hydrogen_bond_donor_hydrogen_atom_idx(const md_hydrogen_bond_candidates_t* hbond_cand, size_t don_idx) {
    ASSERT(hbond_cand);
    if (don_idx < hbond_cand->donor.count) {
        return hbond_cand->donor.h_idx[don_idx];
    }
    return -1;
}

static inline md_atom_idx_t md_hydrogen_bond_acceptor_atom_idx(const md_hydrogen_bond_candidates_t* hbond_cand, size_t acc_idx) {
    ASSERT(hbond_cand);
    if (acc_idx < hbond_cand->acceptor.count) {
        return hbond_cand->acceptor.idx[acc_idx];
    }
    return -1;
}

static inline int md_hydrogen_bond_acceptor_num_lone_pairs(const md_hydrogen_bond_candidates_t* hbond_cand, size_t acc_idx) {
    ASSERT(hbond_cand);
    if (acc_idx < hbond_cand->acceptor.count) {
        return hbond_cand->acceptor.num_lone_pairs[acc_idx];
    }
    return 0;
}


// @NOTE(Robin): This is just to be thorough,
// I would recommend using an explicit arena allocator for the molecule and just clearing that in one go instead of calling this.
// Initialize a system with an allocator. This helper records the allocator on the system so
// subsequent calls that accept a NULL allocator may fall back to `sys->alloc`.

void md_system_reset(md_system_t* sys); // Reset to empty state, maintain allocator
void md_system_free(md_system_t* sys); // Free all memory associated with the system, including the allocator if set.

// Attach helpers: set or create-and-attach trajectories to a system.
void md_system_attach_trajectory(md_system_t* sys, struct md_trajectory_i* traj);

// STATE
// A state owns its coordinate arrays and must be freed with the same allocator it was created with.

// True if the state carries coordinates. A state may legitimately carry a cell but no coordinates
// (a topology only format such as PSF), or coordinates but no cell (xyz without a cell).
static inline bool md_system_state_has_coords(const md_system_state_t* state) {
    return state && state->num_atoms > 0 && state->x && state->y && state->z;
}

static inline bool md_system_state_has_unitcell(const md_system_state_t* state) {
    return state && state->unitcell.flags != 0;
}

// Allocate coordinate storage for num_atoms using state->alloc, which must be set.
// Existing storage is freed first, so this doubles as the reset for a state being reused, and the
// padding up to the simd aligned capacity is zeroed. Returns false if the state has no allocator.
//
// CONTRACT for anything that produces a system and a state together (the format parsers):
// validate sys->alloc and state->alloc, then call md_system_reset(sys) and md_system_state_init(state, N)
// as a pair before touching either. Pass the exact atom count for N when it is known up front and
// write coordinates by index; pass 0 when atoms are filtered while parsing, then reserve with
// md_array_ensure(state->x, capacity, state->alloc) and push. Either way finish with
// state->num_atoms = sys->atom.count, and grow the coordinate arrays with state->alloc and never
// with sys->alloc - md_system_state_free releases them with state->alloc, and the two allocators
// are routinely different.
bool md_system_state_init(md_system_state_t* state, size_t num_atoms);

// Free the coordinate storage and zero the state, preserving the allocator so the state can be
// reused. A view (alloc == NULL) is simply zeroed. Takes no allocator argument by design: the state
// records the one it was allocated with, so it cannot be freed with the wrong one.
void md_system_state_free(md_system_state_t* state);

// Copy src into dst, reallocating dst as needed. dst->alloc must be set.
bool md_system_state_copy(md_system_state_t* dst, const md_system_state_t* src);

#ifdef __cplusplus
}
#endif
