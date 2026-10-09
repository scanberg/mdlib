#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <md_types.h>
#include <md_unitcell.h>
#include <core/md_unit.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#include <md_attributes.h>

typedef struct md_atom_type_data_t {
    size_t count;

    md_label_t*     name;
    // The force field's own name for the type ('opls_135', 'CT', or a Martini bead type such as
    // 'Q5'), when the source has one; empty otherwise. It is what separates two types with the same
    // name that are different particles - a Martini SC1 is a different bead in every residue - and
    // it is too long for a label, so these are owned strings.
    str_t*          ff_type;
    md_atomic_number_t* z;
    float*          mass;
    float*          radius;
    uint32_t*       color;
    md_atom_type_flags_t* flags;    // The particle kind, see md_types.h
} md_atom_type_data_t;

typedef struct md_atom_data_t {
    size_t count;

    // Coordinates live in md_system_state_t, not here. A system is time invariant; where the atoms
    // are is not.
    md_atom_type_idx_t* type_idx;
    md_atom_flags_t* flags;     // Role in the component and chemistry, see md_types.h

    // Chemistry of each atom, as given by the file or perceived by md_chem_perceive. NULL when unknown. Per atom
    // and not per type, as atom types are shared by name (the OD1 of Asp and of Asn).
    int8_t*  formal_charge;     // [count]
    uint8_t* hydrogen_count;    // [count] Hydrogens attached, explicit (bonded H atoms) and implicit

    md_atom_type_data_t type;
} md_atom_data_t;

// Component (Residue): a contiguous range of atoms
typedef struct md_component_data_t {
    size_t count;
    md_label_t* name;
    md_sequence_id_t* seq_id;
    uint32_t* atom_offset;          // [count + 1]
    md_component_flags_t* flags;    // Kind (amino acid, nucleotide, water, ion) and place in the chain, see md_types.h
} md_component_data_t;

// Instance: ONE molecule, or ONE polymer chain, as a contiguous range of components.
//   - A polymer chain is one instance, also when it is broken (residues missing in a crystal structure).
//   - Every other molecule is an instance of its own: each water, each ion, each ligand, each lipid.
//   - A molecule of several components is one instance when they are bonded (a lipid split into head and tails).
// Instances are what a loader gives (mmCIF asyms, the molecules of a GROMACS topology), else they are inferred
// (md_util_system_infer_entity_and_instance).
//
// id is the instance's label (mmCIF label_asym_id) and auth_id the author's chain id (auth_asym_id, the chain id of a
// PDB file), "" when there is none. The id is unique for polymer chains but SHARED by the small molecules of one
// asym: mmCIF puts all the waters of a chain in one asym, and inference gives a run of molecules of the same
// non-polymer entity one id in the same way. Look up instances by index, not by id.
typedef struct md_instance_data_t {
    size_t count;
    md_label_t* id;
    md_label_t* auth_id;
    uint32_t* comp_offset;          // [count + 1]
    md_entity_idx_t* entity_idx;
} md_instance_data_t;

// Entity: a molecule type, or a polymer sequence, of which the system holds one or more instances.
// What it is lives in its flags (md_entity_kind_t), which its instances share.
typedef struct md_entity_data_t {
    size_t count;
    md_label_t* id;
    md_entity_flags_t* flags;       // Kind (peptide, DNA, water, ...) and whether it was inferred, see md_types.h
    str_t* description;
} md_entity_data_t;

// The backbones hold only what is static with regard to the topology: which atoms of which components form them, in
// consecutive ranges along a chain. What depends on a frame (the backbone angles, the secondary structure) or serves a
// single consumer (the ramachandran classification) is computed from them by whoever needs it, into arrays of its own
// with one entry per segment (md_util_backbone_angles_compute, md_util_backbone_secondary_structure_infer,
// md_util_backbone_ramachandran_classify).
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

// A snapshot of the geometric state of a system: where the atoms are and what box they are in.
//
// The FIELDS hold exactly what shares one interpolation contract - same type in and out, periodic
// boundary aware, handled as a unit by md_util_interpolate_*. Nothing which interpolates differently
// becomes a field here, because a field in this struct is a promise that it does.
//
// 'attributes' is how a frame carries everything else. A TRR frame has velocities and forces beside
// its coordinates; they belong to that frame and to no other, and load_frame hands back one object
// precisely so a caller cannot pair one frame's velocities with another frame's positions. They are
// NOT part of the interpolation contract and nothing pretends otherwise: md_util_interpolate_* takes
// raw coordinate arrays and never sees this table, so an interpolated state carries whatever its
// producer chose to put there and no quantity is silently blended.
//
// The table's allocator is the state's own, set by md_system_state_init, so a view state
// (alloc NULL) has no table - the same ownership rule the coordinates follow.
//
// What is computed from a frame's coordinates belongs there too: the backbone angles and the secondary structure of
// the frame (md_util_state_backbone_compute and the accessors beside it in md_util.h).
//
// Two presence bits, both self describing:
//   num_atoms == 0        -> no coordinates
//   unitcell.flags == 0   -> no cell
// num_atoms is used rather than testing xyz != NULL because md_array_ensure allocates capacity
// without setting size, so a non NULL xyz does not imply the coordinates were populated.
//
// xyz is packed, one vec3_t of 12 bytes per atom - the layout every source and every consumer
// already has: the files, a run's atom/position attribute, the GPU buffers. Kernels that want x, y
// and z apart load four or eight atoms and de-interleave in registers (md_mm_load_xyz_packed_ps).
// The array is padded to a multiple of 16 atoms, zeroed past num_atoms, so such loads never need a
// scalar tail for reading.
//
// Ownership: alloc non-NULL means this state owns xyz and must be freed with
// md_system_state_free. alloc NULL means the state is a non owning view over coordinates somebody
// else owns - a scratch arena during script evaluation, a GPU mapped buffer during interpolation.
// Both forms are load bearing, and this is the only thing that distinguishes them.
//
// Set alloc before handing the state to a producer, exactly as md_system_t requires sys->alloc to
// be set before loading. The two then read symmetrically at the call site, and the state's lifetime
// is free to differ from the system's - a temp allocator for the state, a persistent one for the
// system, is a legitimate and useful combination.
// frame is the ordinal of the run frame the coordinates came from, as a continuous quantity: the
// integer part selects the frame, the fractional part is how far between that frame and the next the
// state has been interpolated. A NEGATIVE value means the state did not come from a run at all
// (a topology's own coordinates, a scratch buffer), which is why 0.0 cannot serve as that marker.
// Use md_state_has_frame / md_state_frame_floor / md_state_frame_nearest / md_state_frame_frac
// rather than reading the field: a raw cast of the absent value lands on frame 0, which is in range
// and plausible and therefore the kind of mistake that survives.
//
// CAVEAT, unlike the two presence bits above: -1 does not survive zero initialisation. A state
// built as {0} or with designated initialisers that omit frame reads as frame 0, not as absent.
// md_system_state_init stamps -1, so any state that went through it reads as absent until something
// writes a real frame; a hand rolled {0} does not. Only md_system_extract_frame and the
// interpolation which produces a state write a non negative value.
//
// @NOTE: physical time is deliberately NOT stored here. It is derivable from frame and the run's
// frame axis ("<run>/time"), and storing both would reintroduce the
// very thing this struct exists to prevent - two fields which must agree, with nothing enforcing it.
typedef struct md_system_state_t {
    size_t num_atoms;
    vec3_t* xyz;
    md_unitcell_t unitcell;
    double frame;
    md_attributes_t attributes;   // per frame quantities beyond the interpolation contract above
    md_allocator_i* alloc;
} md_system_state_t;

// This represents the persistent portion (topology) of a system which does not change over time, such as the atom types, bonds, components, etc.
// It may of course be modified through some special operations, though it is not expected to change frequently.
typedef struct md_system_t {
    md_allocator_i*             alloc;

    // TOPOLOGY VERSION. Changes whenever the topology does: the atoms, their types (what the particles are), their
    // flags, the components, instances and entities, the bonds, and what is derived from them (rings, structures,
    // backbones, the perceived chemistry). The properties of the atom types (mass, radius, color) are not topology.
    //
    // A consumer which keeps something derived from the topology (a query, a lookup, a classification) stores the
    // version it was built from and rebuilds when it differs. Versions come from one counter shared by every system and
    // are never handed out twice, so a system freed and loaded anew has a version no consumer has seen. Compare for
    // equality only: it is not a count of changes. 0 is a system no mdlib function has built.
    //
    // Every mdlib function which changes the topology bumps it: md_system_reset (and so every loader),
    // md_util_system_infer and the md_util_system_infer_* functions, md_chem_perceive, md_system_bond_insert and
    // md_system_bond_remove, md_itp_system_supplement. Code which changes the arrays by hand calls
    // md_system_topology_changed afterwards.
    uint64_t                    topology_version;

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
    
    md_index_data_t             ring;               // Ring structures formed by persistent bonds
    md_structure_data_t         structure;          // Isolated structures connected by persistent bonds (plus hierarchy links for coarse grained systems)

    md_assembly_data_t          assembly;           // Assemblies of  (duplications of ranges with new transforms)
    
    // ONE table, holding everything this system carries that is not one of the fields above.
    //
    // Whether a quantity varies over the trajectory is a FLAG on the attribute
    // (MD_ATTRIBUTE_FLAG_TEMPORAL), not a separate table, because a single producer emits both
    // kinds into one namespace: a script yields temporal series AND distributions, and both are
    // "script/...". Splitting by kind would cut a group in half and make enumerating or removing
    // that namespace two operations on two containers.
    //
    // "What varies over time in this dataset" is an md_attributes_iter over a prefix, testing the
    // temporal bit.
    md_attributes_t             attributes;

    // The non-bonded force field (md_nonbonded.h), when the system was loaded with one (a tpr): particle types,
    // charges, pair parameters, exclusions and how the simulation cut the interactions off. NULL otherwise.
    struct md_nb_forcefield_t*  nonbonded;

    str_t                       description;
} md_system_t;

#ifdef __cplusplus
extern "C" {
#endif

// Records that the topology of the system changed: gives it a new topology_version, which it returns
uint64_t md_system_topology_changed(md_system_t* sys);

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

// Adds a type unconditionally, for a loader that decides itself which particles share a type.
// ff_type may be empty.
static inline md_atom_type_idx_t md_atom_type_add(md_atom_type_data_t* atom_type, str_t name, str_t ff_type, md_atomic_number_t z, float mass, float radius, uint32_t color, md_atom_type_flags_t flags, struct md_allocator_i* alloc) {
    ASSERT(atom_type);
    ASSERT(alloc);

    md_array_push(atom_type->name, make_label(name), alloc);
    str_t ff_copy = {0};
    if (!str_empty(ff_type)) {
        ff_copy = str_copy(ff_type, alloc);
    }
    md_array_push(atom_type->ff_type, ff_copy, alloc);
    md_array_push(atom_type->z, z, alloc);
    md_array_push(atom_type->mass, mass, alloc);
    md_array_push(atom_type->radius, radius, alloc);
    md_array_push(atom_type->color, color, alloc);
    md_array_push(atom_type->flags, flags, alloc);
    atom_type->count++;

    return (md_atom_type_idx_t)(atom_type->count - 1);
}

static inline md_atom_type_idx_t md_atom_type_find_or_add(md_atom_type_data_t* atom_type, str_t name, md_atomic_number_t z, float mass, float radius, uint32_t color, md_atom_type_flags_t flags, struct md_allocator_i* alloc) {
    ASSERT(atom_type);
    ASSERT(alloc);
    
    // First try to find existing atom type
    md_atom_type_idx_t type_idx = md_atom_type_find(atom_type, name, z);
    if (type_idx != 0) {
        return type_idx;
    }
    
    const str_t no_ff_type = {0};
    return md_atom_type_add(atom_type, name, no_ff_type, z, mass, radius, color, flags, alloc);
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

static inline md_atom_type_flags_t md_atom_type_flags(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count && type_data->flags) {
        return type_data->flags[type_idx];
    }
    return MD_ATOM_TYPE_FLAG_NONE;
}

static inline md_particle_kind_t md_atom_type_particle_kind(const md_atom_type_data_t* type_data, size_t type_idx) {
    return md_atom_type_flags_particle_kind(md_atom_type_flags(type_data, type_idx));
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

// Empty when the source had no force field type for it
static inline str_t md_atom_type_ff_type(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_data->ff_type && type_idx < type_data->count) {
        return type_data->ff_type[type_idx];
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
        return state->xyz[atom_idx];
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

// Formal charge of the atom, 0 when unknown
static inline int md_atom_formal_charge(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    return (atom->formal_charge && atom_idx < atom->count) ? atom->formal_charge[atom_idx] : 0;
}

// Hydrogens attached to the atom, explicit and implicit; -1 when unknown
static inline int md_atom_hydrogen_count(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    return (atom->hydrogen_count && atom_idx < atom->count) ? atom->hydrogen_count[atom_idx] : -1;
}

static inline md_atom_flags_t md_atom_flags(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    if (atom_idx < atom->count && atom->flags) {
        return atom->flags[atom_idx];
    }
    return MD_ATOM_FLAG_NONE;
}

static inline md_hybridization_t md_atom_hybridization(const md_atom_data_t* atom, size_t atom_idx) {
    return md_atom_flags_hybridization(md_atom_flags(atom, atom_idx));
}

// The kind of particle, read through the atom's type
static inline md_particle_kind_t md_atom_particle_kind(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    if (atom->type_idx && atom_idx < atom->count) {
        return md_atom_type_particle_kind(&atom->type, atom->type_idx[atom_idx]);
    }
    return MD_PARTICLE_ATOM;
}

// Component

static inline size_t md_component_count(const md_component_data_t* comp) {
    ASSERT(comp);
    return comp->count;
}

static inline str_t md_component_name(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
    str_t name = STR_INIT("");
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

static inline md_component_flags_t md_component_flags(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
    md_component_flags_t flags = MD_COMPONENT_FLAG_NONE;
    if (comp->flags && comp_idx < comp->count) {
        flags = comp->flags[comp_idx];
    }
    return flags;
}

static inline md_component_kind_t md_component_kind(const md_component_data_t* comp, size_t comp_idx) {
    return md_component_flags_kind(md_component_flags(comp, comp_idx));
}

// The index i of the range [offset[i], offset[i + 1]) which holds value, -1 when there is none.
// offset is non decreasing and count + 1 long. A binary search.
static inline int32_t md_offset_range_find(const uint32_t* offset, size_t count, size_t value) {
    if (!offset || count == 0 || value < offset[0] || value >= offset[count]) return -1;
    size_t lo = 0, hi = count;  // offset[lo] <= value < offset[hi]
    while (hi - lo > 1) {
        const size_t mid = lo + (hi - lo) / 2;
        if (offset[mid] <= value) lo = mid;
        else hi = mid;
    }
    return (int32_t)lo;
}

static inline md_component_idx_t md_component_find_by_atom_idx(const md_component_data_t* comp, size_t atom_idx) {
    ASSERT(comp);
    return (md_component_idx_t)md_offset_range_find(comp->atom_offset, comp->count, atom_idx);
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
    return (md_instance_idx_t)md_offset_range_find(inst->comp_offset, inst->count, comp_idx);
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
    str_t id = STR_INIT("");
    if (inst->id && inst_idx < inst->count) {
        id = LBL_TO_STR(inst->id[inst_idx]);
    }
    return id;
}

static inline str_t md_instance_auth_id(const md_instance_data_t* inst, size_t inst_idx) {
    ASSERT(inst);
    str_t auth_id = STR_INIT("");
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
    str_t label = STR_INIT("");
    if (entity->id && entity_idx < entity->count) {
        label = LBL_TO_STR(entity->id[entity_idx]);
    }
    return label;
}

static inline str_t md_entity_description(const md_entity_data_t* entity, size_t entity_idx) {
    ASSERT(entity);
    str_t desc = STR_INIT("");
    if (entity->description && entity_idx < entity->count) {
        desc = entity->description[entity_idx];
    }
    return desc;
}

static inline md_entity_flags_t md_entity_flags(const md_entity_data_t* entity, size_t entity_idx) {
    ASSERT(entity);
    md_entity_flags_t flags = MD_ENTITY_FLAG_NONE;
    if (entity->flags && entity_idx < entity->count) {
        flags = entity->flags[entity_idx];
    }
    return flags;
}

static inline md_entity_kind_t md_entity_kind(const md_entity_data_t* entity, size_t entity_idx) {
    return md_entity_flags_kind(md_entity_flags(entity, entity_idx));
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

static inline md_atom_flags_t md_system_atom_flags(const md_system_t* sys, size_t atom_idx) {
    ASSERT(sys);
    return md_atom_flags(&sys->atom, atom_idx);
}

static inline md_particle_kind_t md_system_atom_particle_kind(const md_system_t* sys, size_t atom_idx) {
    ASSERT(sys);
    return md_atom_particle_kind(&sys->atom, atom_idx);
}

static inline size_t md_system_atom_type_count(const md_system_t* sys) {
    ASSERT(sys);
    return md_atom_type_count(&sys->atom.type);
}

static inline md_atom_type_flags_t md_system_atom_type_flags(const md_system_t* sys, size_t type_idx) {
    ASSERT(sys);
    return md_atom_type_flags(&sys->atom.type, type_idx);
}

// A system is coarse grained when any of its particles is a bead
static inline bool md_system_is_coarse_grained(const md_system_t* sys) {
    ASSERT(sys);
    for (size_t i = 0; i < sys->atom.type.count; ++i) {
        if (md_atom_type_particle_kind(&sys->atom.type, i) == MD_PARTICLE_BEAD) return true;
    }
    return false;
}

static inline size_t md_system_component_count(const md_system_t* sys) {
    ASSERT(sys);
    return md_component_count(&sys->component);
}

static inline md_component_flags_t md_system_component_flags(const md_system_t* sys, size_t comp_idx) {
    ASSERT(sys);
    return md_component_flags(&sys->component, comp_idx);
}

static inline md_component_kind_t md_system_component_kind(const md_system_t* sys, size_t comp_idx) {
    ASSERT(sys);
    return md_component_kind(&sys->component, comp_idx);
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

static inline md_entity_flags_t md_system_entity_flags(const md_system_t* sys, size_t ent_idx) {
    ASSERT(sys);
    return md_entity_flags(&sys->entity, ent_idx);
}

static inline md_entity_kind_t md_system_entity_kind(const md_system_t* sys, size_t ent_idx) {
    ASSERT(sys);
    return md_entity_kind(&sys->entity, ent_idx);
}

// The kind of the instance's entity, MD_ENTITY_KIND_UNKNOWN when it has none
static inline md_entity_kind_t md_system_instance_entity_kind(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    if (sys->instance.entity_idx && inst_idx < sys->instance.count) {
        const md_entity_idx_t ent_idx = sys->instance.entity_idx[inst_idx];
        if (ent_idx >= 0) return md_entity_kind(&sys->entity, (size_t)ent_idx);
    }
    return MD_ENTITY_KIND_UNKNOWN;
}

static inline str_t md_system_instance_id(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    str_t id = STR_INIT("");
    if (sys->instance.id && inst_idx < sys->instance.count) {
        id = md_instance_id(&sys->instance, inst_idx);
    }
    return id;
}

static inline str_t md_system_instance_auth_id(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    str_t id = STR_INIT("");
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

// The kind of the component which holds the atom, MD_COMPONENT_KIND_OTHER when it is in none.
// A binary search: a loop over all atoms is better served walking the components.
static inline md_component_kind_t md_system_atom_component_kind(const md_system_t* sys, size_t atom_idx) {
    ASSERT(sys);
    const md_component_idx_t comp_idx = md_system_component_find_by_atom_idx(sys, atom_idx);
    return comp_idx >= 0 ? md_component_kind(&sys->component, (size_t)comp_idx) : MD_COMPONENT_KIND_OTHER;
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
    md_system_topology_changed(sys);
}

static inline void md_system_bond_remove(md_system_t* sys, md_bond_idx_t bond_idx) {
    ASSERT(sys);
    md_bond_remove(&sys->bond, bond_idx);
    md_system_topology_changed(sys);
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


// @NOTE(Robin): This is just to be thorough,
// I would recommend using an explicit arena allocator for the molecule and just clearing that in one go instead of calling this.
// Initialize a system with an allocator. This helper records the allocator on the system so
// subsequent calls that accept a NULL allocator may fall back to `sys->alloc`.

void md_system_reset(md_system_t* sys); // Reset to empty state, maintain allocator
void md_system_free(md_system_t* sys); // Free all memory associated with the system, including the allocator if set.

// STATE
// A state owns its coordinate arrays and must be freed with the same allocator it was created with.

// True if the state carries coordinates. A state may legitimately carry a cell but no coordinates
// (a topology only format such as PSF), or coordinates but no cell (xyz without a cell).
static inline bool md_system_state_has_coords(const md_system_state_t* state) {
    return state && state->num_atoms > 0 && state->xyz;
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
// md_array_ensure(state->xyz, capacity, state->alloc) and push. Either way finish with
// state->num_atoms = sys->atom.count, and grow the coordinate array with state->alloc and never
// with sys->alloc - md_system_state_free releases them with state->alloc, and the two allocators
// are routinely different.
bool md_system_state_init(md_system_state_t* state, size_t num_atoms);

// Free the coordinate storage and zero the state, preserving the allocator so the state can be
// reused. A view (alloc == NULL) is simply zeroed. Takes no allocator argument by design: the state
// records the one it was allocated with, so it cannot be freed with the wrong one.
void md_system_state_free(md_system_state_t* state);

// Copy src into dst, reallocating dst as needed. dst->alloc must be set.
bool md_system_state_copy(md_system_state_t* dst, const md_system_state_t* src);

// EXTRACTION
//
// A state IS a snapshot of a run's temporal attributes at one frame, and these take it. What to take
// is said once, up front, and the context that comes back is then asked for frames:
//
//     const str_t paths[] = { STR_LIT("atom/position"), STR_LIT("unitcell") };
//     md_system_extract_t* ex = md_system_extract_begin(sys, run, paths, 2, alloc);
//     for (f = beg; f < end; ++f) md_system_extract_frame(ex, f, &state);
//     md_system_extract_end(ex);
//
// Saying it once is what lets the context be worth having. It resolves the paths a single time, and
// it keeps what a source is expensive to reopen - the files a trajectory streams from - open for as
// long as it lives. Opening a file is cheap on a local disk and not on the file servers of a
// cluster, where every open is a round trip to a metadata server that every rank shares.
//
// Paths are relative to the run ("<run>/<path>") and each lands in the state:
//
//   atom/position   into xyz, straight from the source: both are packed.
//   unitcell        into unitcell, a {F,3,3} box per frame with row i box vector i.
//   anything else   into out->attributes under the same relative path, as the value at that frame -
//                   a {F,N} c3 attribute arrives as {N} c3 - and without the temporal flag, since a
//                   snapshot has no frame axis. An attribute along an axis of its own (an energy file
//                   written more often than the coordinates, velocities written less often) is read
//                   at the row whose time matches the frame. A frame with no such row leaves it OUT
//                   of the state: its absence is the answer, and a state reused from frame to frame
//                   never shows another frame's value in its place.
//
// Every extract stamps out->frame. out->num_atoms must be 0 or the run's atom count, and a state
// asking for positions must have storage for that many; nothing is allocated for coordinates here.
//
// OWNERSHIP AND LIFETIME. The context belongs to the caller and to ONE thread at a time; a thread
// pool keeps one per thread. It holds attribute ids rather than pointers and looks each up again on
// every frame, so an attribute removed underneath it fails the extract instead of being read. It must
// still be ENDED before the run it reads is removed, since the files it holds are that run's.
typedef struct md_system_extract_t md_system_extract_t;

// NULL when the run, or any of the paths in it, does not exist. alloc is what the context lives in;
// the thread that begins it need not be the one that uses it, but only one may use it at a time.
md_system_extract_t* md_system_extract_begin(const md_system_t* sys, str_t run, const str_t paths[], size_t num_paths, struct md_allocator_i* alloc);
bool                 md_system_extract_frame(md_system_extract_t* ex, int64_t frame, md_system_state_t* out);
void                 md_system_extract_end(md_system_extract_t* ex);

// PUBLISHING A RUN
//
// What a trajectory format publishes is the same shape whatever the format (see RUNS):
//
//     <run>/time            F64 {F}        the frame axis; without a unit it holds frame ordinals
//     <run>/step            I64 {F}        the simulation step of each frame, when the file says
//     <run>/unitcell        F32 {F,3,3}    Angstrom, row i box vector i, zero where a frame has none
//     <run>/atom/position   F32 {F,N} c3   Angstrom, VIRTUAL: read from the file per frame
//     <run>/source/path     STR            the file the frames are read from
//     <run>/source/offset   I64 {F}        where each frame starts in it
//     <run>/source/size     I64 {F}        and how many bytes it has
//
// md_run_publish puts all of it in place from one description, or none of it: a publish that fails
// part way removes the run's prefix. What a format adds beyond this (velocities, how its frames are
// laid out) it publishes itself afterwards, under the same run, removing the prefix on failure too.
// Flags for a format's run publisher
typedef enum md_run_flag_t {
    MD_RUN_FLAG_NONE                = 0,
    MD_RUN_FLAG_DISABLE_CACHE_WRITE = 1,    // leave no '<file>.cache' index beside the file
} md_run_flag_t;

typedef uint32_t md_run_flags_t;

// INDEX CACHE. A format that has to scan its file to know where the frames are keeps what the scan
// learned beside it, as '<file>.cache', so the scan is paid once per file. The cache starts with
// this header and the format appends its own blocks after it.
//
// A cache is only as good as its match with the file: md_run_cache_open accepts it when the magic
// and version are the format's and the file's size and modification time, as the operating system
// reports them now, are the ones the cache was made from. Anything else - a run still being written,
// a file replaced by one of the same size - is a cache to make again.
typedef struct md_run_cache_header_t {
    uint64_t       magic;
    uint64_t       version;
    uint64_t       source_size;        // bytes
    md_file_time_t source_modified;    // the OS modification time, nanoseconds since the epoch
    uint64_t       num_atoms;
    uint64_t       num_frames;
} md_run_cache_header_t;

// Opens '<source_path>.cache' and reads its header. True, with the file left open just past the
// header for the format's own blocks, when the cache matches the source file as it is now.
bool md_run_cache_open(md_file_t* out_file, md_run_cache_header_t* out_header, str_t source_path, uint64_t magic, uint64_t version);

// Creates '<source_path>.cache' and writes the header, stamped with the source file's size and
// modification time as they were when it was scanned: take scanned with md_file_info_extract before
// the scan, so a file that grows during it leaves a cache that does not match it. True, with the
// file left open for the format's own blocks.
bool md_run_cache_create(md_file_t* out_file, str_t source_path, const md_file_info_t* scanned, uint64_t magic, uint64_t version, size_t num_atoms, size_t num_frames);

typedef struct md_run_desc_t {
    size_t          num_frames;
    size_t          num_atoms;

    const double*   time;               // num_frames values, required
    md_unit_t       time_unit;          // none when the file does not know time: time is then ordinals
    const int64_t*  step;               // num_frames values, or NULL when the file has none

    // The cell: RESIDENT from num_frames * 9 floats, or VIRTUAL from a provider (a format that
    // keeps it in the frames), or neither for a run without one.
    const float*                  unitcell;
    const md_attribute_virtual_t* unitcell_virt;

    str_t           source_path;        // copied
    const int64_t*  source_offset;      // num_frames values
    const int64_t*  source_size;        // num_frames values

    const md_attribute_virtual_t* position_virt;   // required
} md_run_desc_t;

bool md_run_publish(md_system_t* sys, str_t run, const md_run_desc_t* desc);

// "<run>/<leaf>" into buf. Empty when it does not fit.
str_t md_run_path(char* buf, size_t cap, str_t run, str_t leaf);

// For a provider of a run's attribute: the run it belongs to - its path minus "/<leaf>" - and the
// file and frame table the run reads from. False, and logged, when the attribute is not at <leaf>
// or the run has lost its source attributes.
typedef struct md_run_source_t {
    str_t          run;
    str_t          path;
    const int64_t* offset;
    const int64_t* size;
    size_t         num_frames;
} md_run_source_t;

bool md_run_source(md_run_source_t* out, const md_attributes_t* attributes, const md_attribute_t* attr, str_t leaf);

// SERIES ALONG A RUN
//
// Columns of numbers sampled along a run - an .xvg or .csv of per frame quantities - published as a
// group below it, "<run>/<group>", so that attr() and extraction read them like anything else in the
// run. Each column becomes "<run>/<group>/<name>", the name folded to lower case letters, digits and
// '_' (the original is the label); a name taken twice gets a suffix.
//
// With a time column the group has a frame axis of its own, "<run>/<group>/time", and every frame of
// the run must find its time in it, as an energy file's must; the columns are read at the matching
// row. Time without a unit is taken to be in the run's. Without a time column the rows ARE the run's
// frames, and there must be exactly as many.
//
// The group replaces one of the same name, all or nothing. source_path, when given, is published as
// "<run>/<group>/source" so a session can load it again.
typedef struct md_run_series_desc_t {
    str_t               group;          // below the run, e.g. "xvg/energy"
    size_t              num_rows;
    size_t              num_columns;
    const double*       time;           // num_rows values, or NULL: the rows are the run's frames
    md_unit_t           time_unit;
    const str_t*        names;          // num_columns
    const md_unit_t*    units;          // num_columns, or NULL for none
    const float* const* columns;        // num_columns arrays of num_rows values
    str_t               source_path;    // optional
} md_run_series_desc_t;

bool md_run_publish_series(md_system_t* sys, str_t run, const md_run_series_desc_t* desc);

#ifdef __cplusplus
}
#endif
