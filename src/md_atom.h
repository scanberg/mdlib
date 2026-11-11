#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_types.h>
#include <core/md_array.h>
#include <core/md_label.h>

struct md_system_t;
struct md_executor_t;

// Atoms and their particle types (md_atom_data_t, md_particle_type_data_t).

// ### PARTICLE TYPE ### (see CLASSIFICATION in md_types.h)

// What a particle is. A property of its type: every particle of a type is the same kind of particle.
typedef enum md_particle_kind_t {
    MD_PARTICLE_ATOM         = 0,   // An atom. Its element is the type's atomic number, 0 when it is not known.
    MD_PARTICLE_BEAD         = 1,   // A coarse grained bead which stands for a group of atoms. A system with beads is
                                    // coarse grained (md_system_is_coarse_grained): bonds and hydrogen bonds are not
                                    // inferred from the geometry, as atomic radii do not apply.
    MD_PARTICLE_VIRTUAL_SITE = 2,   // A massless interaction site without an element: the M site of 4 site water models
                                    // (TIP4P, OPC), the lone pairs of 5 site models (TIP5P). It belongs to its molecule
                                    // but forms no bonds.
} md_particle_kind_t;

typedef enum {
    MD_PARTICLE_TYPE_FLAG_NONE              = 0,
    MD_PARTICLE_TYPE_FLAG_PARTICLE_KIND_MASK = 0x3,     // md_particle_kind_t, see md_particle_type_flags_particle_kind
} md_particle_type_flags_t;

MD_ENUM_FLAGS(md_particle_type_flags_t)

// ### ATOM ###

// Hybridization of a heavy atom, see md_chem.h. SP2 includes the N and O whose lone pair takes part in a pi system
// (amide and aniline N, pyrrole N, furan O).
typedef enum md_hybridization_t {
    MD_HYBRIDIZATION_UNKNOWN = 0,   // Not perceived, or not applicable (hydrogens, metals)
    MD_HYBRIDIZATION_SP      = 1,
    MD_HYBRIDIZATION_SP2     = 2,
    MD_HYBRIDIZATION_SP3     = 3,
} md_hybridization_t;

typedef enum {
    MD_ATOM_FLAG_NONE                   = 0,

    // The atom's role in an amino acid or nucleotide of a chain, found from the covalent bonds whatever the atoms are
    // called, also without residues (POLYMER ATOMS in md_util.c, md_util_system_infer_residue_flags), and for the
    // beads of the predefined coarse grained types. A match which is no residue of a chain (a ligand which holds the
    // atoms of one) gets no roles.
    MD_ATOM_FLAG_BACKBONE               = 0x1,      // N, CA, C of an amino acid; P, O5', C5', C4', C3', O3' of a nucleotide
    MD_ATOM_FLAG_SIDE_CHAIN             = 0x2,      // Amino acid side chain, from CB outwards
    MD_ATOM_FLAG_NUCLEOSIDE             = 0x4,      // Sugar and base of a nucleotide
    MD_ATOM_FLAG_NUCLEOBASE             = 0x8,      // Base of a nucleotide
    MD_ATOM_FLAG_TERMINAL_BEG           = 0x10,     // In the terminal group which begins a chain (N-terminus, 5' end)
    MD_ATOM_FLAG_TERMINAL_END           = 0x20,     // In the terminal group which ends a chain (C-terminus, 3' end)
    // 0x40 unused: the atoms of the tails of the lipids are in md_lipid_data_t (md_lipid_tail_atom_mask)
    // 0x80 unused: the atoms of a metallic phase are those of its metallic bonds (md_bond_metallic_atoms)

    // Chemistry, see md_chem.h (md_chem_perceive)
    MD_ATOM_FLAG_HYBRIDIZATION_MASK     = 0x300,    // md_hybridization_t, see md_atom_flags_hybridization
    MD_ATOM_FLAG_AROMATIC               = 0x400,    // In an aromatic ring

    // Treated quantum mechanically: an atom of the QM calculation the system was loaded from, set by its loader.
    // A system from a QM calculation alone is QM throughout; one with an embedding or a QM/MM partition has its
    // QM region flagged, and the atoms without the flag are its environment. Nothing infers it: a system no
    // loader said it of has no QM region. Selected in scripts by 'qm' and 'environment'.
    MD_ATOM_FLAG_QM                     = 0x800,

    // Hydrogen bond roles are not flags: they depend on the hydrogen bond model and its options, and a donor is a
    // D-H pair rather than an atom. See md_hbond_perceive_roles and md_hbond_query_t.
} md_atom_flags_t;

MD_ENUM_FLAGS(md_atom_flags_t)

static inline md_particle_kind_t md_particle_type_flags_particle_kind(md_particle_type_flags_t flags) {
    return (md_particle_kind_t)((unsigned int)flags & (unsigned int)MD_PARTICLE_TYPE_FLAG_PARTICLE_KIND_MASK);
}

static inline md_particle_type_flags_t md_particle_type_flags_set_particle_kind(md_particle_type_flags_t flags, md_particle_kind_t kind) {
    const unsigned int mask = (unsigned int)MD_PARTICLE_TYPE_FLAG_PARTICLE_KIND_MASK;
    return (md_particle_type_flags_t)(((unsigned int)flags & ~mask) | ((unsigned int)kind & mask));
}

#define MD_ATOM_HYBRIDIZATION_SHIFT 8

static inline md_hybridization_t md_atom_flags_hybridization(md_atom_flags_t flags) {
    return (md_hybridization_t)(((unsigned int)flags & (unsigned int)MD_ATOM_FLAG_HYBRIDIZATION_MASK) >> MD_ATOM_HYBRIDIZATION_SHIFT);
}

static inline md_atom_flags_t md_atom_flags_set_hybridization(md_atom_flags_t flags, md_hybridization_t hyb) {
    const unsigned int mask = (unsigned int)MD_ATOM_FLAG_HYBRIDIZATION_MASK;
    return (md_atom_flags_t)(((unsigned int)flags & ~mask) | (((unsigned int)hyb << MD_ATOM_HYBRIDIZATION_SHIFT) & mask));
}

#ifdef __cplusplus
extern "C" {
#endif

// Names for printing ("atom", "bead", "sp2", ...), "" for a value out of range
const char* md_particle_kind_name(md_particle_kind_t kind);
const char* md_hybridization_name(md_hybridization_t hyb);

#ifdef __cplusplus
}
#endif

typedef struct md_particle_type_data_t {
    size_t count;

    md_label_t*     name;
    // The force field's own name for the type ('opls_135', 'CT', or a Martini bead type such as
    // 'Q5'), when the source has one; empty otherwise. It is what separates two types with the same
    // name that are different particles - a Martini SC1 is a different bead in every residue. Up to 16
    // characters (md_label16_t), which every naming scheme of a force field fits.
    md_label16_t*   ff_type;
    md_atomic_number_t* z;
    float*          mass;
    float*          radius;
    uint32_t*       color;
    md_particle_type_flags_t* flags; // The particle kind, see md_particle_type_flags_t above
} md_particle_type_data_t;

typedef struct md_atom_data_t {
    size_t count;

    // Coordinates live in md_system_state_t, not here. A system is time invariant; where the atoms
    // are is not.
    md_particle_type_idx_t* type_idx;
    md_atom_flags_t* flags;     // Role in the residue and chemistry, see md_atom_flags_t above

    // Chemistry of each atom, as given by the file or perceived by md_chem_perceive. NULL when unknown. Per atom
    // and not per type, as particle types are shared by name (the OD1 of Asp and of Asn).
    int8_t*  formal_charge;     // [count]
    uint8_t* hydrogen_count;    // [count] Hydrogens attached, explicit (bonded H atoms) and implicit

    md_particle_type_data_t type;
} md_atom_data_t;

#ifdef __cplusplus
extern "C" {
#endif

// Particle type table helper functions
static inline size_t md_particle_type_count(const md_particle_type_data_t* particle_type) {
    MD_ASSERT(particle_type);
    return particle_type->count;
}

static inline md_particle_type_idx_t md_particle_type_find(const md_particle_type_data_t* particle_type, md_str_t name, md_atomic_number_t z) {
    MD_ASSERT(particle_type);
    // The name as a type keeps it, its first 8 characters, compared as one word
    const uint64_t key = md_label_key(md_label_make(name));
    for (size_t i = 1; i < particle_type->count; ++i) {
        if (md_label_key(particle_type->name[i]) == key && particle_type->z[i] == z) {
            return (md_particle_type_idx_t)i;
        }
    }
    return 0; // Zero is sentinel for "not found"
}

// Adds a type unconditionally, for a loader that decides itself which particles share a type.
// ff_type may be empty.
static inline md_particle_type_idx_t md_particle_type_add(md_particle_type_data_t* particle_type, md_str_t name, md_str_t ff_type, md_atomic_number_t z, float mass, float radius, uint32_t color, md_particle_type_flags_t flags, struct md_allocator_i* alloc) {
    MD_ASSERT(particle_type);
    MD_ASSERT(alloc);

    md_array_push(particle_type->name, md_label_make(name), alloc);
    md_array_push(particle_type->ff_type, md_label16_make(ff_type), alloc);
    md_array_push(particle_type->z, z, alloc);
    md_array_push(particle_type->mass, mass, alloc);
    md_array_push(particle_type->radius, radius, alloc);
    md_array_push(particle_type->color, color, alloc);
    md_array_push(particle_type->flags, flags, alloc);
    particle_type->count++;

    return (md_particle_type_idx_t)(particle_type->count - 1);
}

static inline md_particle_type_idx_t md_particle_type_find_or_add(md_particle_type_data_t* particle_type, md_str_t name, md_atomic_number_t z, float mass, float radius, uint32_t color, md_particle_type_flags_t flags, struct md_allocator_i* alloc) {
    MD_ASSERT(particle_type);
    MD_ASSERT(alloc);
    
    // First try to find existing particle type
    md_particle_type_idx_t type_idx = md_particle_type_find(particle_type, name, z);
    if (type_idx != 0) {
        return type_idx;
    }
    
    const md_str_t no_ff_type = {0};
    return md_particle_type_add(particle_type, name, no_ff_type, z, mass, radius, color, flags, alloc);
}

static inline md_atomic_number_t md_particle_type_atomic_number(const md_particle_type_data_t* type_data, size_t type_idx) {
    MD_ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->z[type_idx];
    }
    return 0;
}

static inline float md_particle_type_mass(const md_particle_type_data_t* type_data, size_t type_idx) {
    MD_ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->mass[type_idx];
    }
    return 0;
}

static inline md_particle_type_flags_t md_particle_type_flags(const md_particle_type_data_t* type_data, size_t type_idx) {
    MD_ASSERT(type_data);
    if (type_idx < type_data->count && type_data->flags) {
        return type_data->flags[type_idx];
    }
    return MD_PARTICLE_TYPE_FLAG_NONE;
}

static inline md_particle_kind_t md_particle_type_particle_kind(const md_particle_type_data_t* type_data, size_t type_idx) {
    return md_particle_type_flags_particle_kind(md_particle_type_flags(type_data, type_idx));
}

static inline float md_particle_type_radius(const md_particle_type_data_t* type_data, size_t type_idx) {
    MD_ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->radius[type_idx];
    }
    return 0;
}

static inline uint32_t md_particle_type_color(const md_particle_type_data_t* type_data, size_t type_idx) {
    MD_ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->color[type_idx];
    }
    return 0;
}

static inline md_str_t md_particle_type_name(const md_particle_type_data_t* type_data, size_t type_idx) {
    MD_ASSERT(type_data);
    if (type_idx < type_data->count) {
        return MD_LABEL_STR(type_data->name[type_idx]);
    }
    return MD_STR_LIT("");
}

// Empty when the source had no force field type for it
static inline md_str_t md_particle_type_ff_type(const md_particle_type_data_t* type_data, size_t type_idx) {
    MD_ASSERT(type_data);
    if (type_data->ff_type && type_idx < type_data->count) {
        return md_label16_str(&type_data->ff_type[type_idx]);
    }
    return MD_STR_LIT("");
}

// Atom helpers

static inline md_particle_type_idx_t md_atom_particle_type_idx(const md_atom_data_t* atom, size_t atom_idx) {
    MD_ASSERT(atom);
    if (atom->type_idx && atom_idx < atom->count) {
        return atom->type_idx[atom_idx];
	}
    return 0;
}

static inline md_atomic_number_t md_atom_atomic_number(const md_atom_data_t* atom, size_t atom_idx) {
    MD_ASSERT(atom);
    
    // Try particle type table first if type_idx is available
    if (atom->type_idx && (size_t)atom->type_idx[atom_idx] < atom->type.count) {
        return atom->type.z[atom->type_idx[atom_idx]];
    }
    
    return 0;
}

static inline size_t md_atom_count(const md_atom_data_t* atom_data) {
    MD_ASSERT(atom_data);
    return atom_data->count;
}

static inline float md_atom_mass(const md_atom_data_t* atom, size_t atom_idx) {
    MD_ASSERT(atom);
    
    if (atom_idx < atom->count) {
		return md_particle_type_mass(&atom->type, atom->type_idx[atom_idx]);
    }
    
    return 0.0f;
}

static inline float md_atom_radius(const md_atom_data_t* atom, size_t atom_idx) {
    MD_ASSERT(atom);
    
    if (atom_idx < atom->count) {
        return md_particle_type_radius(&atom->type, atom->type_idx[atom_idx]);
    }
    return 0.0f;
}

static inline md_str_t md_atom_name(const md_atom_data_t* atom, size_t atom_idx) {
    MD_ASSERT(atom);
    if (atom_idx < atom->count) {
        return md_particle_type_name(&atom->type, atom->type_idx[atom_idx]);
    }
    return MD_STR_LIT("");
}

// Formal charge of the atom, 0 when unknown
static inline int md_atom_formal_charge(const md_atom_data_t* atom, size_t atom_idx) {
    MD_ASSERT(atom);
    return (atom->formal_charge && atom_idx < atom->count) ? atom->formal_charge[atom_idx] : 0;
}

// Hydrogens attached to the atom, explicit and implicit; -1 when unknown
static inline int md_atom_hydrogen_count(const md_atom_data_t* atom, size_t atom_idx) {
    MD_ASSERT(atom);
    return (atom->hydrogen_count && atom_idx < atom->count) ? atom->hydrogen_count[atom_idx] : -1;
}

static inline md_atom_flags_t md_atom_flags(const md_atom_data_t* atom, size_t atom_idx) {
    MD_ASSERT(atom);
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
    MD_ASSERT(atom);
    if (atom->type_idx && atom_idx < atom->count) {
        return md_particle_type_particle_kind(&atom->type, atom->type_idx[atom_idx]);
    }
    return MD_PARTICLE_ATOM;
}

// Convenience functions to extract atom properties into arrays
static inline void md_atom_extract_radii(float out_radii[], size_t offset, size_t length, const md_atom_data_t* atom_data) {
    MD_ASSERT(out_radii);
    MD_ASSERT(atom_data);
    MD_ASSERT(offset + length <= atom_data->count);
    
    for (size_t i = 0; i < length; ++i) {
        out_radii[i] = md_atom_radius(atom_data, offset + i);
    }
}

static inline void md_atom_extract_masses(float out_masses[], size_t offset, size_t length, const md_atom_data_t* atom_data) {
    MD_ASSERT(out_masses);
    MD_ASSERT(atom_data);
    MD_ASSERT(offset + length <= atom_data->count);
    
    for (size_t i = 0; i < length; ++i) {
        out_masses[i] = md_atom_mass(atom_data, offset + i);
    }
}

static inline void md_atom_extract_atomic_numbers(md_atomic_number_t out_z[], size_t offset, size_t length, const md_atom_data_t* atom_data) {
    MD_ASSERT(out_z);
    MD_ASSERT(atom_data);
    MD_ASSERT(offset + length <= atom_data->count);
    
    for (size_t i = 0; i < length; ++i) {
        out_z[i] = md_atom_atomic_number(atom_data, offset + i);
    }
}

// ### PARTICLE TYPES OF A SYSTEM, IN PLACE ###

// Identify particle types within the system from the atom labels and the names of their residues. Coarse grained beads
// come from predefined tables. Martini beads of amino acids (BB, SC#) and nucleotides (BB#, SC#) are recognized in any
// residue named as one, and get the masses and sizes of the version of Martini (2, polarizable 2 or 3) the names of
// the system's residues and beads tell (Martini 3 when they tell nothing); a massless bead is a virtual site. Atoms of
// one label with a different element, mass, size or particle kind are of different types.
void md_system_infer_particle_types(struct md_system_t* sys, const md_str_t atom_labels[], const struct md_executor_t* ex);
// The same, with the labels as labels (8 characters): what a loader keeps per atom at half the size of a string.
void md_system_infer_particle_types_from_labels(struct md_system_t* sys, const md_label_t atom_labels[], const struct md_executor_t* ex);

// Applies the predefined particle type tables (coarse grained beads) to atoms that already have a type,
// for loaders that know more about their atoms than a name does (a tpr carries the mass and the LJ
// parameters of every bead). Only types without an element (z == 0) are considered. The table gives
// the atoms their role (backbone, side chain), the residues their kind (amino acid, water) where
// they have none, and the types their particle kind where the loader did not say. Mass and radius
// are left alone: they are the loader's.
void md_system_augment_particle_types(struct md_system_t* sys);

#ifdef __cplusplus
}
#endif
