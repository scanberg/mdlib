#pragma once

// The inference of a system: what follows from its atoms, residues and coordinates (bonds, rings, chains, parts,
// molecules, backbones, lipids, chemistry, fragments), each step declared with the data it infers (md_bond.h,
// md_ring.h, md_hierarchy.h, md_backbone.h, md_lipid.h, md_atom.h, md_chem.h, md_fragment.h). md_system_infer runs them
// in order and gives the system what they infer.

#include <md_system.h>

struct md_allocator_i;
struct md_executor_t;
struct md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    MD_SYSTEM_INFER_NONE                = 0,
    MD_SYSTEM_INFER_COLOR_BIT           = 1u << 1,
    MD_SYSTEM_INFER_BOND_BIT            = 1u << 2,
    MD_SYSTEM_INFER_PART_BIT            = 1u << 3,     // The chains and the parts with their types (md_chain_infer, _parts)
    MD_SYSTEM_INFER_BACKBONE_BIT        = 1u << 4,
    MD_SYSTEM_INFER_MOLECULE_BIT        = 1u << 5,
    MD_SYSTEM_INFER_RESIDUE_BIT         = 1u << 6,     // The kinds of the residues and the roles of their atoms (md_system_infer_residue_flags), after the bonds which tell them
	//MD_SYSTEM_INFER_SECONDARY_STRUCTURE_BIT = 0x0400,
	MD_SYSTEM_INFER_UNWRAP_MOLECULE_BIT = 1u << 7,
    MD_SYSTEM_INFER_CHEMISTRY_BIT       = 1u << 8,     // Bond orders, aromaticity, formal charges, hydrogen counts and hybridization, see md_chem.h
    MD_SYSTEM_INFER_LIPID_BIT           = 1u << 9,     // The heads and tails of the lipids (md_lipid_infer), and the kinds of their residues
    MD_SYSTEM_INFER_FRAGMENT_BIT        = 1u << 10,    // Ring systems, functional groups, backbones (md_util_system_infer_fragments, md_fragment.h), last

    MD_SYSTEM_INFER_ALL                 = -1,
} md_infer_flags_t;

MD_ENUM_FLAGS(md_infer_flags_t)

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
bool md_system_infer(struct md_system_t* sys, const struct md_system_state_t* state, md_infer_flags_t flags, const struct md_executor_t* executor);

// INFERENCE WHILE THE SYSTEM IS IN USE. Inference reads the CORE of a system (the atoms, their types, the residues, the
// instances, the assemblies) and writes the rest, the DERIVED (md_system_fields.h). So a system can be shown and used
// while it is inferred on another thread, without a copy of it:
//
//   md_system_infer_prepare   The work: of the system, what inference writes, copied into alloc; the coordinates
//                                  in state, which become its reference; the per atom arrays of the core shared. On
//                                  the thread which owns the system.
//   md_system_infer_run       Inference on the work: on any thread, while the system is read (and its core
//                                  is not written: what changes the topology waits for the work, or discards it).
//                                  All but MD_SYSTEM_INFER_COLOR_BIT, which is of the types of the atoms: the system
//                                  is given its colors before, md_system_infer with that bit alone.
//   md_system_infer_publish   The derived arrays of the work become the system's, which keeps its core, its
//                                  attributes and its force field. On the thread which owns the system, nothing
//                                  else reading it. False, and nothing done, when the topology of the system changed
//                                  after the work was prepared: the work is of a system which is no more.
//
// The arrays of the work are in alloc, which the caller keeps with the system once published (they become the system's)
// and frees otherwise (the work is never md_system_free'd: it shares arrays of the system). The types of the atoms are
// copied too (a few per system), so that the system's may be changed meanwhile (their colors, their radii): they stay
// the system's when published.
typedef struct md_system_infer_work_t {
    struct md_system_t sys;     // Inferred on: md_system_infer_run
    uint64_t base_version;      // The topology version of the system it was prepared of
} md_system_infer_work_t;

bool md_system_infer_prepare(md_system_infer_work_t* work, const struct md_system_t* sys, const struct md_system_state_t* state, struct md_allocator_i* alloc);
bool md_system_infer_run(md_system_infer_work_t* work, md_infer_flags_t flags, const struct md_executor_t* executor);
bool md_system_infer_publish(struct md_system_t* sys, md_system_infer_work_t* work);

#ifdef __cplusplus
}
#endif
