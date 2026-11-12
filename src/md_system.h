#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <md_types.h>
#include <md_ring.h>
#include <md_fragment.h>

#include <md_atom.h>
#include <md_hierarchy.h>
#include <md_backbone.h>
#include <md_lipid.h>
#include <md_bond.h>
#include <md_assembly.h>
#include <md_group.h>
#include <md_system_state.h>

// @NOTE: do NOT add <math.h> to this header. simde-math.h keys off HUGE_VAL to detect "a math header
// was already included", and when that fires under C++ it assumes the header was <cmath> and starts
// emitting std::trunc. On glibc that is harmless (math.h defines an isnan MACRO, which simde checks
// for first) and on gcc/clang __has_builtin wins before it matters - but MSVC has neither, so any
// C++ translation unit that reaches simde through this header fails with "trunc is not a member of
// std". md_system.h is included nearly everywhere, so it is the worst possible place to trip that.

// This represents the persistent portion (topology) of a system which does not change over time, such as the particle types, bonds, residues, etc.
// It may of course be modified through some special operations, though it is not expected to change frequently.
typedef struct md_system_t {
    md_allocator_i*             alloc;

    // TOPOLOGY VERSION. Changes whenever the topology does: the atoms, their types (what the particles are), their
    // flags, the residues, the instances, the parts and their types, the bonds, and what is derived from them (rings,
    // molecules and their types, backbones, the perceived chemistry). The properties of the particle types (mass, radius, color) are not topology.
    //
    // A consumer which keeps something derived from the topology (a query, a lookup, a classification) stores the
    // version it was built from and rebuilds when it differs. Versions come from one counter shared by every system and
    // are never handed out twice, so a system freed and loaded anew has a version no consumer has seen. Compare for
    // equality only: it is not a count of changes. 0 is a system no mdlib function has built.
    //
    // Every mdlib function which changes the topology bumps it: md_system_reset (and so every loader),
    // md_system_infer and the md_util_system_infer_* functions, md_chem_perceive, md_system_bond_insert and
    // md_system_bond_remove, md_itp_system_supplement. Code which changes the arrays by hand calls
    // md_system_topology_changed afterwards.
    uint64_t                    topology_version;

    // The state from which the derived topology below (bonds, rings, molecules, backbones) was
    // inferred. Written by md_system_infer as part of performing the inference, so it is by
    // construction the input which produced that topology and cannot go stale.
    //
    // RULE: only inference and load time completion read this. Every other operation takes the
    // state it works on as an explicit parameter. Reaching for sys->reference to avoid threading
    // a state through is how the previous coupling arose.
    md_system_state_t           reference;

    md_atom_data_t              atom;
    md_residue_data_t           residue;
    md_instance_data_t          instance;           // What the data gives, the input of inference, see md_instance_data_t
    md_chain_data_t             chain;              // The instances by the chain id the data gives them, see md_chain_data_t
    md_part_data_t              part;               // Strands, glycan trees, ligands, ions, waters, see md_part_data_t

    md_protein_backbone_data_t  protein_backbone;
    md_nucleic_backbone_data_t  nucleic_backbone;
    md_lipid_data_t             lipid;              // The heads and tails of the lipids, see md_lipid_data_t
    
    md_bond_data_t              bond;               // Persistent bonds: the bond graph, then the non-covalent bonds (md_bond_data_t)
    
    md_ring_data_t              ring;               // Rings formed by persistent bonds: 3 to 7 ordered atoms and their flags
    md_molecule_data_t          molecule;           // The molecules: connected by persistent bonds, within residues and along polymer chains (md_molecule_t)
    md_fragment_data_t          fragment;           // What the molecules are made of: ring systems, functional groups, backbones (md_fragment.h)

    md_assembly_data_t          assembly;           // The biological assemblies of a deposited structure, see md_assembly_data_t
    md_group_data_t             group;              // Named sets of atoms the data gives (GROMACS index groups), see md_group_data_t
    
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

    md_str_t                    description;
} md_system_t;

#ifdef __cplusplus
extern "C" {
#endif

// Records that the topology of the system changed: gives it a new topology_version, which it returns
uint64_t md_system_topology_changed(md_system_t* sys);

// ASSEMBLIES (md_assembly_data_t). Sets the assemblies of the system from the rows of a file and its operators: the
// rows of an assembly id are merged into its groups in canonical form, the instances of each operator the union of those
// of the rows which name it, and the operators with the same instances one group. The assemblies are in the order their
// ids first appear, an assembly's groups in the order of their first operator, and the operators of a group in the order
// they first appear. Instances out of range are left out; a row without instances or operators adds nothing.
// An operator which is not rigid is logged, and kept. Replaces the previous assemblies; false (and none) on bad input.
bool md_system_assemblies_init(md_system_t* sys, const md_assembly_gen_t* gen, size_t gen_count, const md_f32x4x4_t* oper_transform, const md_label_t* oper_id, size_t oper_count);
void md_system_assemblies_free(md_system_t* sys);

// GROUPS (md_group_data_t)

// Adds a group of atoms, a copy of atoms, with the roles it has (bits of md_group_role_t, 0 for none) and returns its
// index. A group of the same name and atoms is not added again: it gains the roles, and its index is returned.
// Atoms past the system's are left out. Records a change of the topology (scripts select groups by name).
size_t md_system_group_add(md_system_t* sys, md_str_t name, md_bitset_t atoms, uint32_t roles);

// Sets the groups of a role, in its order, replacing what it had: indices of groups, which gain the role
void md_system_group_set_role(md_system_t* sys, md_group_role_t role, const uint32_t* group_idx, size_t count);

void md_system_groups_free(md_system_t* sys);

// RESIDUES (md_residue_data_t)

static inline size_t md_system_residue_count(const md_system_t* sys) {
    MD_ASSERT(sys);
    return sys->residue.count;
}

static inline md_str_t md_system_residue_name(const md_system_t* sys, size_t res_idx) {
    MD_ASSERT(sys);
    if (sys->residue.name && res_idx < sys->residue.count) return MD_LABEL_STR(sys->residue.name[res_idx]);
    return MD_STR_LIT("");
}

static inline md_sequence_id_t md_system_residue_seq_id(const md_system_t* sys, size_t res_idx) {
    MD_ASSERT(sys);
    return (sys->residue.seq_id && res_idx < sys->residue.count) ? sys->residue.seq_id[res_idx] : 0;
}

static inline md_urange_t md_system_residue_atom_range(const md_system_t* sys, size_t res_idx) {
    MD_ASSERT(sys);
    md_urange_t range = {0, 0};
    if (sys->residue.atom_offset && res_idx < sys->residue.count) {
        range.beg = sys->residue.atom_offset[res_idx];
        range.end = sys->residue.atom_offset[res_idx + 1];
    }
    return range;
}

static inline size_t md_system_residue_atom_count(const md_system_t* sys, size_t res_idx) {
    const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
    return range.end - range.beg;
}

static inline md_residue_flags_t md_system_residue_flags(const md_system_t* sys, size_t res_idx) {
    MD_ASSERT(sys);
    return (sys->residue.flags && res_idx < sys->residue.count) ? sys->residue.flags[res_idx] : MD_RESIDUE_FLAG_NONE;
}

// The kind of the residue: that of its type (md_residue_type_data_t). OTHER for a residue of no type yet (one never
// classified, before md_system_infer_residue_types)
static inline md_residue_kind_t md_system_residue_kind(const md_system_t* sys, size_t res_idx) {
    MD_ASSERT(sys);
    if (!sys->residue.type_idx || res_idx >= sys->residue.count || res_idx >= md_array_size(sys->residue.type_idx)) return MD_RESIDUE_KIND_OTHER;
    const md_residue_type_idx_t t = sys->residue.type_idx[res_idx];
    return (t >= 0 && (size_t)t < sys->residue.type.count) ? sys->residue.type.kind[t] : MD_RESIDUE_KIND_OTHER;
}

// The residue of the atom, -1 for an atom in none
static inline md_residue_idx_t md_system_residue_find_by_atom_idx(const md_system_t* sys, size_t atom_idx) {
    MD_ASSERT(sys);
    return (md_residue_idx_t)md_offset_range_find(sys->residue.atom_offset, sys->residue.count, atom_idx);
}

// The type of the residue, -1 for one of none yet
static inline md_residue_type_idx_t md_system_residue_type_idx(const md_system_t* sys, size_t res_idx) {
    MD_ASSERT(sys);
    return (sys->residue.type_idx && res_idx < sys->residue.count && res_idx < md_array_size(sys->residue.type_idx)) ? sys->residue.type_idx[res_idx] : -1;
}

// PART TYPES (md_part_type_data_t)

static inline size_t md_system_part_type_count(const md_system_t* sys) {
    MD_ASSERT(sys);
    return sys->part.type.count;
}

// The part type of the id, -1 for none
static inline md_part_type_idx_t md_system_part_type_find_by_id(const md_system_t* sys, md_str_t id) {
    MD_ASSERT(sys);
    if (sys->part.type.id) {
        for (size_t i = 0; i < sys->part.type.count; ++i) {
            if (md_str_eq(MD_LABEL_STR(sys->part.type.id[i]), id)) return (md_part_type_idx_t)i;
        }
    }
    return -1;
}

static inline md_str_t md_system_part_type_id(const md_system_t* sys, size_t type_idx) {
    MD_ASSERT(sys);
    if (sys->part.type.id && type_idx < sys->part.type.count) return MD_LABEL_STR(sys->part.type.id[type_idx]);
    return MD_STR_LIT("");
}

static inline md_str_t md_system_part_type_description(const md_system_t* sys, size_t type_idx) {
    MD_ASSERT(sys);
    if (sys->part.type.description && type_idx < sys->part.type.count) return sys->part.type.description[type_idx];
    return MD_STR_LIT("");
}

static inline md_part_type_flags_t md_system_part_type_flags(const md_system_t* sys, size_t type_idx) {
    MD_ASSERT(sys);
    return (sys->part.type.flags && type_idx < sys->part.type.count) ? sys->part.type.flags[type_idx] : MD_PART_TYPE_FLAG_NONE;
}

static inline md_part_kind_t md_system_part_type_kind(const md_system_t* sys, size_t type_idx) {
    return md_part_type_flags_kind(md_system_part_type_flags(sys, type_idx));
}

// INSTANCES (md_instance_data_t)

static inline size_t md_system_instance_count(const md_system_t* sys) {
    MD_ASSERT(sys);
    return sys->instance.count;
}

static inline md_urange_t md_system_instance_residue_range(const md_system_t* sys, size_t inst_idx) {
    MD_ASSERT(sys);
    md_urange_t range = {0, 0};
    if (sys->instance.res_offset && inst_idx < sys->instance.count) {
        range.beg = sys->instance.res_offset[inst_idx];
        range.end = sys->instance.res_offset[inst_idx + 1];
    }
    return range;
}

static inline size_t md_system_instance_residue_count(const md_system_t* sys, size_t inst_idx) {
    const md_urange_t range = md_system_instance_residue_range(sys, inst_idx);
    return range.end - range.beg;
}

// The atoms of the instance: one range, as its residues are consecutive
static inline md_urange_t md_system_instance_atom_range(const md_system_t* sys, size_t inst_idx) {
    const md_urange_t res = md_system_instance_residue_range(sys, inst_idx);
    md_urange_t range = {0, 0};
    if (res.beg < res.end) {
        range.beg = md_system_residue_atom_range(sys, res.beg).beg;
        range.end = md_system_residue_atom_range(sys, res.end - 1).end;
    }
    return range;
}

static inline size_t md_system_instance_atom_count(const md_system_t* sys, size_t inst_idx) {
    const md_urange_t range = md_system_instance_atom_range(sys, inst_idx);
    return range.end - range.beg;
}

static inline md_str_t md_system_instance_id(const md_system_t* sys, size_t inst_idx) {
    MD_ASSERT(sys);
    if (sys->instance.id && inst_idx < sys->instance.count) return MD_LABEL_STR(sys->instance.id[inst_idx]);
    return MD_STR_LIT("");
}

static inline md_str_t md_system_instance_auth_id(const md_system_t* sys, size_t inst_idx) {
    MD_ASSERT(sys);
    if (sys->instance.auth_id && inst_idx < sys->instance.count) return MD_LABEL_STR(sys->instance.auth_id[inst_idx]);
    return MD_STR_LIT("");
}

// The part type the data gives the instance (an mmCIF entity), -1 for none
static inline md_part_type_idx_t md_system_instance_part_type_idx(const md_system_t* sys, size_t inst_idx) {
    MD_ASSERT(sys);
    return (sys->instance.part_type_idx && inst_idx < sys->instance.count) ? sys->instance.part_type_idx[inst_idx] : -1;
}

// The molecule type the data gives the instance (a molecule type of a topology), -1 for none
static inline md_molecule_type_idx_t md_system_instance_molecule_type_idx(const md_system_t* sys, size_t inst_idx) {
    MD_ASSERT(sys);
    return (sys->instance.mol_type_idx && inst_idx < sys->instance.count) ? sys->instance.mol_type_idx[inst_idx] : -1;
}

// The instance of the residue, -1 for none
static inline md_instance_idx_t md_system_instance_find_by_residue_idx(const md_system_t* sys, size_t res_idx) {
    MD_ASSERT(sys);
    return (md_instance_idx_t)md_offset_range_find(sys->instance.res_offset, sys->instance.count, res_idx);
}

// The instance of the atom, -1 for none
static inline md_instance_idx_t md_system_instance_find_by_atom_idx(const md_system_t* sys, size_t atom_idx) {
    const md_residue_idx_t res_idx = md_system_residue_find_by_atom_idx(sys, atom_idx);
    return res_idx >= 0 ? md_system_instance_find_by_residue_idx(sys, (size_t)res_idx) : -1;
}

// SYSTEM
// System level convenience accessors

static inline size_t md_system_atom_count(const md_system_t* sys) {
    MD_ASSERT(sys);
    return md_atom_count(&sys->atom);
}

static inline md_atom_flags_t md_system_atom_flags(const md_system_t* sys, size_t atom_idx) {
    MD_ASSERT(sys);
    return md_atom_flags(&sys->atom, atom_idx);
}

static inline md_particle_kind_t md_system_atom_particle_kind(const md_system_t* sys, size_t atom_idx) {
    MD_ASSERT(sys);
    return md_atom_particle_kind(&sys->atom, atom_idx);
}

static inline size_t md_system_particle_type_count(const md_system_t* sys) {
    MD_ASSERT(sys);
    return md_particle_type_count(&sys->atom.type);
}

static inline md_particle_type_flags_t md_system_particle_type_flags(const md_system_t* sys, size_t type_idx) {
    MD_ASSERT(sys);
    return md_particle_type_flags(&sys->atom.type, type_idx);
}

// A system is coarse grained when any of its particles is a bead
static inline bool md_system_is_coarse_grained(const md_system_t* sys) {
    MD_ASSERT(sys);
    for (size_t i = 0; i < sys->atom.type.count; ++i) {
        if (md_particle_type_particle_kind(&sys->atom.type, i) == MD_PARTICLE_BEAD) return true;
    }
    return false;
}

// The bonds of the bond graph. The persistent non-covalent bonds come after them (md_bond_data_t).
static inline size_t md_system_bond_count(const md_system_t* sys) {
    MD_ASSERT(sys);
    return sys->bond.count;
}

static inline size_t md_system_chain_count(const md_system_t* sys) {
    MD_ASSERT(sys);
    return sys->chain.count;
}

// The author's identifier of the chain
static inline md_str_t md_system_chain_id(const md_system_t* sys, size_t chain_idx) {
    MD_ASSERT(sys);
    md_str_t id = MD_STR_INIT("");
    if (sys->chain.id && chain_idx < sys->chain.count) {
        id = MD_LABEL_STR(sys->chain.id[chain_idx]);
    }
    return id;
}

// The runs of consecutive residues of the chain, ascending: [beg, end)
static inline const md_urange_t* md_system_chain_run_beg(const md_system_t* sys, size_t chain_idx) {
    MD_ASSERT(sys);
    return (sys->chain.run && chain_idx < sys->chain.count) ? sys->chain.run + sys->chain.run_offset[chain_idx] : NULL;
}

static inline const md_urange_t* md_system_chain_run_end(const md_system_t* sys, size_t chain_idx) {
    MD_ASSERT(sys);
    return (sys->chain.run && chain_idx < sys->chain.count) ? sys->chain.run + sys->chain.run_offset[chain_idx + 1] : NULL;
}

// The atoms of a run of residues
static inline md_urange_t md_system_chain_run_atom_range(const md_system_t* sys, md_urange_t run) {
    MD_ASSERT(sys);
    md_urange_t r = {0, 0};
    if (sys->residue.atom_offset && run.beg < run.end && run.end <= sys->residue.count) {
        r.beg = sys->residue.atom_offset[run.beg];
        r.end = sys->residue.atom_offset[run.end];
    }
    return r;
}

// The chain of the residue, -1 for a residue in none
static inline md_chain_idx_t md_system_chain_find_by_residue_idx(const md_system_t* sys, size_t res_idx) {
    MD_ASSERT(sys);
    return (sys->chain.res_chain && res_idx < sys->residue.count) ? sys->chain.res_chain[res_idx] : -1;
}

// The chain of the atom, -1 for an atom in none
static inline md_chain_idx_t md_system_chain_find_by_atom_idx(const md_system_t* sys, size_t atom_idx) {
    MD_ASSERT(sys);
    const md_residue_idx_t res_idx = md_system_residue_find_by_atom_idx(sys, atom_idx);
    return res_idx >= 0 ? md_system_chain_find_by_residue_idx(sys, (size_t)res_idx) : -1;
}

// The atoms of the chain: the atom ranges of its runs. Mostly one, as a chain is mostly one strand, but the ligands and
// waters a file files under a chain id come after all of the polymers, and they are of the chain too.
static inline size_t md_system_chain_atom_count(const md_system_t* sys, size_t chain_idx) {
    MD_ASSERT(sys);
    size_t n = 0;
    const md_urange_t* end = md_system_chain_run_end(sys, chain_idx);
    for (const md_urange_t* it = md_system_chain_run_beg(sys, chain_idx); it && it < end; ++it) {
        const md_urange_t r = md_system_chain_run_atom_range(sys, *it);
        n += r.end - r.beg;
    }
    return n;
}

static inline size_t md_system_part_count(const md_system_t* sys) {
    MD_ASSERT(sys);
    return sys->part.count;
}

// The residues of the part
static inline md_urange_t md_system_part_residue_range(const md_system_t* sys, size_t part_idx) {
    MD_ASSERT(sys);
    md_urange_t r = {0, 0};
    if (sys->part.res_offset && part_idx < sys->part.count) {
        r.beg = sys->part.res_offset[part_idx];
        r.end = sys->part.res_offset[part_idx + 1];
    }
    return r;
}

// The atoms of the part: one range, as its residues are consecutive
static inline md_urange_t md_system_part_atom_range(const md_system_t* sys, size_t part_idx) {
    MD_ASSERT(sys);
    const md_urange_t residues = md_system_part_residue_range(sys, part_idx);
    md_urange_t r = {0, 0};
    if (residues.beg < residues.end && sys->residue.atom_offset && residues.end <= sys->residue.count) {
        r.beg = sys->residue.atom_offset[residues.beg];
        r.end = sys->residue.atom_offset[residues.end];
    }
    return r;
}

// The part of the residue, -1 for none (a system whose parts were not inferred)
static inline md_part_idx_t md_system_part_find_by_residue_idx(const md_system_t* sys, size_t res_idx) {
    MD_ASSERT(sys);
    return (md_part_idx_t)md_offset_range_find(sys->part.res_offset, sys->part.count, res_idx);
}

// The part of the atom, -1 for an atom in none
static inline md_part_idx_t md_system_part_find_by_atom_idx(const md_system_t* sys, size_t atom_idx) {
    MD_ASSERT(sys);
    const md_residue_idx_t res_idx = md_system_residue_find_by_atom_idx(sys, atom_idx);
    return res_idx >= 0 ? md_system_part_find_by_residue_idx(sys, (size_t)res_idx) : -1;
}

// The type of the part, -1 for none
static inline md_part_type_idx_t md_system_part_type_idx(const md_system_t* sys, size_t part_idx) {
    MD_ASSERT(sys);
    return (sys->part.type_idx && part_idx < sys->part.count) ? sys->part.type_idx[part_idx] : -1;
}

// What the part is, the kind of its type: MD_PART_KIND_UNKNOWN when it has none
static inline md_part_kind_t md_system_part_kind(const md_system_t* sys, size_t part_idx) {
    const md_part_type_idx_t type_idx = md_system_part_type_idx(sys, part_idx);
    return type_idx >= 0 ? md_system_part_type_kind(sys, (size_t)type_idx) : MD_PART_KIND_UNKNOWN;
}

static inline size_t md_system_molecule_count(const md_system_t* sys) {
    MD_ASSERT(sys);
    return sys->molecule.count;
}

// The molecule of the atom, -1 when the molecules have not been inferred (md_molecule_infer). The atoms
// of a molecule need not be consecutive: it is found by the atom's slot (md_molecule_data_t).
static inline md_molecule_idx_t md_system_molecule_find_by_atom_idx(const md_system_t* sys, size_t atom_idx) {
    MD_ASSERT(sys);
    const md_molecule_data_t* st = &sys->molecule;
    if (!st->atom_slot || !st->offset || atom_idx >= sys->atom.count) return -1;
    return (md_molecule_idx_t)md_offset_range_find(st->offset, st->count, (size_t)st->atom_slot[atom_idx]);
}

static inline size_t md_system_molecule_atom_count(const md_system_t* sys, size_t mol_idx) {
    MD_ASSERT(sys);
    const md_molecule_data_t* st = &sys->molecule;
    return (st->offset && mol_idx < st->count) ? st->offset[mol_idx + 1] - st->offset[mol_idx] : 0;
}

// The type of the molecule, -1 for none
static inline md_molecule_type_idx_t md_system_molecule_type_idx(const md_system_t* sys, size_t mol_idx) {
    MD_ASSERT(sys);
    return (sys->molecule.type_idx && mol_idx < sys->molecule.count) ? sys->molecule.type_idx[mol_idx] : -1;
}

static inline size_t md_system_molecule_type_count(const md_system_t* sys) {
    MD_ASSERT(sys);
    return sys->molecule.type.count;
}

static inline md_str_t md_system_molecule_type_name(const md_system_t* sys, size_t type_idx) {
    MD_ASSERT(sys);
    if (sys->molecule.type.name && type_idx < sys->molecule.type.count) return sys->molecule.type.name[type_idx];
    return MD_STR_LIT("");
}

static inline md_molecule_type_flags_t md_system_molecule_type_flags(const md_system_t* sys, size_t type_idx) {
    MD_ASSERT(sys);
    return (sys->molecule.type.flags && type_idx < sys->molecule.type.count) ? sys->molecule.type.flags[type_idx] : MD_MOLECULE_TYPE_FLAG_NONE;
}

// The atoms of the molecule in the order of its tree, with the atom each was reached from (md_molecule_t). Of no atoms
// for an index out of range.
static inline md_molecule_t md_system_molecule(const md_system_t* sys, size_t mol_idx) {
    MD_ASSERT(sys);
    const md_molecule_data_t* st = &sys->molecule;
    md_molecule_t mol = {0};
    if (st->offset && mol_idx < st->count) {
        const uint32_t beg = st->offset[mol_idx];
        mol.count      = st->offset[mol_idx + 1] - beg;
        mol.atom_idx   = st->atom_idx + beg;
        mol.parent_idx = st->parent_idx + beg;
    }
    return mol;
}

// The slot of the atom in the arrays of the molecules (md_molecule_data_t). Requires the molecules
// (md_molecule_infer).
static inline int32_t md_system_molecule_atom_slot(const md_system_t* sys, size_t atom_idx) {
    MD_ASSERT(sys);
    MD_ASSERT(sys->molecule.atom_slot);
    return sys->molecule.atom_slot[atom_idx];
}

// The atom the atom was reached from in the tree of its molecule; a root is its own parent. Requires the molecules
// (md_molecule_infer).
static inline int32_t md_system_molecule_atom_parent(const md_system_t* sys, size_t atom_idx) {
    MD_ASSERT(sys);
    MD_ASSERT(sys->molecule.atom_slot);
    return sys->molecule.parent_idx[sys->molecule.atom_slot[atom_idx]];
}

// The kind of the residue which holds the atom, MD_RESIDUE_KIND_OTHER when it is in none.
// A binary search: a loop over all atoms is better served walking the residues.
static inline md_residue_kind_t md_system_atom_residue_kind(const md_system_t* sys, size_t atom_idx) {
    MD_ASSERT(sys);
    const md_residue_idx_t res_idx = md_system_residue_find_by_atom_idx(sys, atom_idx);
    return res_idx >= 0 ? md_system_residue_kind(sys, (size_t)res_idx) : MD_RESIDUE_KIND_OTHER;
}

void md_system_bond_build_connectivity(md_system_t* sys);

static inline md_bond_idx_t md_system_bond_find(const md_system_t* sys, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b) {
    MD_ASSERT(sys);
    return md_bond_find(&sys->bond, atom_idx_a, atom_idx_b);
}

static inline void md_system_bond_insert(md_system_t* sys, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b, md_bond_flags_t flags) {
    MD_ASSERT(sys);
    md_bond_insert(&sys->bond, atom_idx_a,  atom_idx_b, flags, sys->alloc);
    md_system_topology_changed(sys);
}

static inline void md_system_bond_remove(md_system_t* sys, md_bond_idx_t bond_idx) {
    MD_ASSERT(sys);
    md_bond_remove(&sys->bond, bond_idx);
    md_system_topology_changed(sys);
}

// Flags of ring ring_idx of sys->ring, NONE when they have not been inferred
static inline md_ring_flags_t md_system_ring_flags(const md_system_t* sys, size_t ring_idx) {
    MD_ASSERT(sys);
    return md_ring_flags(&sys->ring, ring_idx);
}

static inline md_bond_flags_t md_system_bond_flags(const md_system_t* sys, md_bond_idx_t bond_idx) {
    MD_ASSERT(sys);
    if (bond_idx >= 0 && (size_t)bond_idx < sys->bond.count + sys->bond.noncovalent_count) {
        return sys->bond.flags[bond_idx];
    }
    return MD_BOND_FLAG_NONE;
}

// @NOTE(Robin): This is just to be thorough,
// I would recommend using an explicit arena allocator for the molecule and just clearing that in one go instead of calling this.
// Initialize a system with an allocator. This helper records the allocator on the system so
// subsequent calls that accept a NULL allocator may fall back to `sys->alloc`.

void md_system_reset(md_system_t* sys); // Reset to empty state, maintain allocator
void md_system_free(md_system_t* sys); // Free all memory associated with the system, including the allocator if set.

// A copy of src in dst's allocator, which dst keeps: every array of the system, its reference state and description,
// and its topology version (it is the same topology). Not the attribute table, which dst gets empty, nor the force
// field (nonbonded): what is computed on a copy (inference, in the background) reads neither. What dst had is freed.
bool md_system_copy(md_system_t* dst, const md_system_t* src);

#ifdef __cplusplus
}
#endif
