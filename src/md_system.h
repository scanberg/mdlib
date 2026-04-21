#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <md_types.h>

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

    // Coordinates
    float* x;
    float* y;
    float* z;

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
typedef struct md_structure_data_t {
    size_t count;
    uint32_t* offset; // Offsets into the structure fields stored bellow, includes sentinel at the end (so length is count + 1)
    int32_t* atom_idx;
    int32_t* parent_idx;
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

typedef struct md_system_state_t {
    float* atom_x;
    float* atom_y;
    float* atom_z;
    md_unitcell_t unitcell;
} md_system_state_t;

// This represents the persistent portion (topology) of a system which does not change over time, such as the atom types, bonds, components, etc.
// It may of course be modified through some special operations, though it is not expected to change frequently.
typedef struct md_system_t {
    md_allocator_i*             alloc;
    md_trajectory_i*            trajectory;

    md_unitcell_t               initial_unitcell;
    md_unitcell_t               unitcell; // compatibility alias for consumers still using sys->unitcell

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

    str_t                       description;
} md_system_t;

#ifdef __cplusplus
extern "C" {
#endif



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

static inline vec3_t md_atom_coord(const md_atom_data_t* atom_data, size_t atom_idx) {
    ASSERT(atom_data);
    if (atom_idx < atom_data->count) {
        return vec3_set(atom_data->x[atom_idx], atom_data->y[atom_idx], atom_data->z[atom_idx]);
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

// Ensure that the dst_sys has the allocator set
bool md_system_copy(md_system_t* dst_sys, const md_system_t* src_sys);

// Attach helpers: set or create-and-attach trajectories to a system.
void md_system_attach_trajectory(md_system_t* sys, struct md_trajectory_i* traj);

#ifdef __cplusplus
}
#endif
