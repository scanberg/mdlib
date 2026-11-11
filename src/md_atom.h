#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_types.h>

// Atoms and their particle types (md_atom_data_t, md_particle_type_data_t).

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
    md_particle_type_flags_t* flags; // The particle kind, see md_types.h
} md_particle_type_data_t;

typedef struct md_atom_data_t {
    size_t count;

    // Coordinates live in md_system_state_t, not here. A system is time invariant; where the atoms
    // are is not.
    md_particle_type_idx_t* type_idx;
    md_atom_flags_t* flags;     // Role in the residue and chemistry, see md_types.h

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
    const uint64_t key = md_label_key(make_label(name));
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

    md_array_push(particle_type->name, make_label(name), alloc);
    md_array_push(particle_type->ff_type, make_label16(ff_type), alloc);
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
        return LBL_TO_STR(type_data->name[type_idx]);
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

#ifdef __cplusplus
}
#endif
