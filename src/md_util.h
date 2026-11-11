#pragma once

#include <md_system.h>
#include <md_element.h>

#include <core/md_str.h>
#include <core/md_vec.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>

struct md_allocator_i;
struct md_executor_t;
struct md_bitset_t;
struct md_bitmask_t;

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    MD_UTIL_INFER_NONE                = 0,
    MD_UTIL_INFER_COLOR_BIT           = 1u << 1,
    MD_UTIL_INFER_BOND_BIT            = 1u << 2,
    MD_UTIL_INFER_PART_BIT            = 1u << 3,     // The chains and the parts with their types (md_chain_infer, _parts)
    MD_UTIL_INFER_BACKBONE_BIT        = 1u << 4,
    MD_UTIL_INFER_MOLECULE_BIT        = 1u << 5,
    MD_UTIL_INFER_RESIDUE_BIT         = 1u << 6,     // The kinds of the residues and the roles of their atoms (md_system_infer_residue_flags), after the bonds which tell them
	//MD_UTIL_INFER_SECONDARY_STRUCTURE_BIT = 0x0400,
	MD_UTIL_INFER_UNWRAP_MOLECULE_BIT = 1u << 7,
    MD_UTIL_INFER_CHEMISTRY_BIT       = 1u << 8,     // Bond orders, aromaticity, formal charges, hydrogen counts and hybridization, see md_chem.h
    MD_UTIL_INFER_LIPID_BIT           = 1u << 9,     // The heads and tails of the lipids (md_lipid_infer), and the kinds of their residues
    MD_UTIL_INFER_FRAGMENT_BIT        = 1u << 10,    // Ring systems, functional groups, backbones (md_util_system_infer_fragments, md_fragment.h), last

    MD_UTIL_INFER_ALL                 = -1,
} md_infer_flags_t;

MD_ENUM_FLAGS(md_infer_flags_t)

bool md_util_resname_dna(md_str_t str);
bool md_util_resname_rna(md_str_t str);
bool md_util_resname_acidic(md_str_t str);
bool md_util_resname_basic(md_str_t str);
bool md_util_resname_neutral(md_str_t str);
bool md_util_resname_water(md_str_t str);
bool md_util_resname_ion(md_str_t str);         // The common monatomic ions (NA, CL, ...)
bool md_util_resname_lipid(md_str_t str);      // The lipids of the Martini lipidomes (their names are mostly the atomistic ones)
bool md_util_resname_hydrophobic(md_str_t str);
bool md_util_resname_amino_acid(md_str_t str);
bool md_util_resname_nucleotide(md_str_t str);

void md_util_system_extract_xyzw_from_mask(md_f32x4_t* out_xyzw, const struct md_bitset_t* mask, const md_system_t* sys, const md_system_state_t* state);

// Grow a mask by bonds up to a certain extent (counted as number of bonds from the original mask)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_mask_grow_by_bonds(struct md_bitmask_t* mask, const struct md_system_t* sys, size_t extent, const struct md_bitset_t* viable_mask);

// Grow a mask by radius (in Angstrom)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_mask_grow_by_radius(struct md_bitmask_t* mask, const md_system_state_t* state, double radius, const struct md_bitset_t* viable_mask);

// The type of an atom of the name and element, added when the system has none, with the mass, size and color of the
// element. map, which the caller keeps for a load (zeroed, with an allocator), finds a type it has given before in one
// lookup rather than a search of the types.
md_particle_type_idx_t md_util_particle_type_find_or_add(md_system_t* sys, md_hashmap32_t* map, md_str_t name, md_atomic_number_t z);

// Attempts to generate missing data such as covalent bonds, chains, secondary structures, backbone angles etc.
// Infers the derivable parts of a system (covalent bonds, rings, molecules, backbones, hydrogen
// bond roles) from the supplied state, and records that state as sys->reference. A system without residues has its
// backbones found from the covalent bonds alone (md_protein_backbone_data_t).
//
// The state is an explicit parameter rather than read off the system so that the recorded
// reference is by construction the input which produced the topology. Requires coordinates for
// anything geometric; a state without them (a topology only format such as PSF) still yields the
// purely graph derived parts, provided the bonds were supplied by the format.
//
// The phases run one after the other, each on the threads of the executor where it spreads its work (NULL: the calling
// thread); the result is the same whatever the threads.
bool md_util_system_infer(struct md_system_t* sys, const md_system_state_t* state, md_infer_flags_t flags, const struct md_executor_t* executor);

// INFERENCE WHILE THE SYSTEM IS IN USE. Inference reads the CORE of a system (the atoms, their types, the residues, the
// instances, the assemblies) and writes the rest, the DERIVED (md_system_fields.h). So a system can be shown and used
// while it is inferred on another thread, without a copy of it:
//
//   md_util_system_infer_prepare     The work: of the system, what inference writes, copied into alloc; the coordinates
//                                    in state, which become its reference; the per atom arrays of the core shared. On
//                                    the thread which owns the system.
//   md_util_system_infer_run         Inference on the work: on any thread, while the system is read (and its core
//                                    is not written: what changes the topology waits for the work, or discards it).
//                                    All but MD_UTIL_INFER_COLOR_BIT, which is of the types of the atoms: the system
//                                    is given its colors before, md_util_system_infer with that bit alone.
//   md_util_system_infer_publish     The derived arrays of the work become the system's, which keeps its core, its
//                                    attributes and its force field. On the thread which owns the system, nothing
//                                    else reading it. False, and nothing done, when the topology of the system changed
//                                    after the work was prepared: the work is of a system which is no more.
//
// The arrays of the work are in alloc, which the caller keeps with the system once published (they become the system's)
// and frees otherwise (the work is never md_system_free'd: it shares arrays of the system). The types of the atoms are
// copied too (a few per system), so that the system's may be changed meanwhile (their colors, their radii): they stay
// the system's when published.
typedef struct md_util_infer_work_t {
    struct md_system_t sys;     // Inferred on: md_util_system_infer_run
    uint64_t base_version;      // The topology version of the system it was prepared of
} md_util_infer_work_t;

bool md_util_system_infer_prepare(md_util_infer_work_t* work, const struct md_system_t* sys, const md_system_state_t* state, struct md_allocator_i* alloc);
bool md_util_system_infer_run(md_util_infer_work_t* work, md_infer_flags_t flags, const struct md_executor_t* executor);
bool md_util_system_infer_publish(struct md_system_t* sys, md_util_infer_work_t* work);

#ifdef __cplusplus
}
#endif
