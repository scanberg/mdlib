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

// Infers the derived of a system (md_system_fields.h), all of it anew, from its core and the supplied state, and
// records that state as sys->reference: the bonds, the coordination, the kinds of the residues and the roles of their
// atoms, the chains and the parts, the molecules, the rings, the backbones, the chemistry, the lipids, the kinds of the
// residues and the fragments, in that order.
//
// The bonds are those given (md_system_t::given_bond) with those inferred from the coordinates between the atoms whose
// bonds are not given complete (md_atom_data_t::bonds_complete, md_bond_merge): of a system a topology describes whole,
// none is inferred, and of one without bonds given (most structure files) all of them are.
//
// The state is an explicit parameter rather than read off the system so that the recorded reference is by
// construction the input which produced the topology. Requires coordinates for anything geometric; a state without
// them (a topology only format such as PSF) still yields the purely graph derived parts, from the bonds the format
// gives.
//
// The phases run one after the other, each on the threads of the executor where it spreads its work (NULL: the calling
// thread); the result is the same whatever the threads.
//
// The chemistry is perceived for the molecules it infers, and a molecule which is a copy of another
// (md_molecule_data_t, STRUCTURES) with all of its hydrogens explicit takes the chemistry of its reference: where the
// geometry chooses between equivalent forms (the double bond of a carboxylate, the charge of a guanidinium, the Kekule
// structure of a ring), every copy has the form of its reference. md_chem_perceive does the same while the structures
// are of the system (md_system_molecule_structures_valid).
bool md_system_infer(struct md_system_t* sys, const struct md_system_state_t* state, const struct md_executor_t* executor);

// The colors of the particle types which have none: that of the element (CPK), or one of the index of the type. Of
// the core (the types are what the data gives), so not of md_system_infer: on the thread which owns the system, once
// it is loaded.
void md_system_infer_type_colors(struct md_system_t* sys);

// INFERENCE WHILE THE SYSTEM IS IN USE. Inference reads the CORE of a system (the atoms, their types, the residues, the
// instances, the assemblies) and writes the rest, the DERIVED (md_system_fields.h). So a system can be shown and used
// while it is inferred on another thread, without a copy of it:
//
//   md_system_infer_prepare   The work: of the system, what inference writes, copied into alloc; the coordinates
//                                  in state, which become its reference; the per atom arrays of the core shared. On
//                                  the thread which owns the system.
//   md_system_infer_run       Inference on the work: on any thread, while the system is read (and its core
//                                  is not written: what changes the topology waits for the work, or discards it).
//   md_system_infer_publish   The derived arrays of the work become the system's, which keeps its core, its
//                                  attributes and its force field. On the thread which owns the system, nothing
//                                  else reading it. False, and nothing done, when the topology of the system changed
//                                  after the work was prepared: the work is of a system which is no more.
//
// What the work infers (its DERIVED arrays and its reference state) is allocated with the system's allocator of the
// derived (md_system_t, ALLOCATORS), which must allow several threads (the heap does), and becomes the system's when
// published, or is freed with md_system_infer_discard. The rest of the work, the types of the atoms (copied, a few per
// system, so that the system's may be changed meanwhile: their colors, their radii), is in alloc, which the caller
// frees after either. Of a system without an allocator of the derived, all of the work is in alloc, which the caller
// then keeps with the system once published (its arrays become the system's). The work is never md_system_free'd: it
// shares arrays of the system.
typedef struct md_system_infer_work_t {
    struct md_system_t sys;     // Inferred on: md_system_infer_run
    uint64_t base_version;      // The topology version of the system it was prepared of
} md_system_infer_work_t;

bool md_system_infer_prepare(md_system_infer_work_t* work, const struct md_system_t* sys, const struct md_system_state_t* state, struct md_allocator_i* alloc);
bool md_system_infer_run(md_system_infer_work_t* work, const struct md_executor_t* executor);
bool md_system_infer_publish(struct md_system_t* sys, md_system_infer_work_t* work);
// A work not published (cancelled, or false from md_system_infer_publish): what it inferred freed
void md_system_infer_discard(md_system_infer_work_t* work);

#ifdef __cplusplus
}
#endif
