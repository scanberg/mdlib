#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <md_types.h>
#include <core/md_common.h>

struct md_system_t;
struct md_allocator_i;

struct md_allocator_i;

// Lipids: the heads and tails of each lipid molecule.

// What a tail of a lipid is (md_lipid_data_t)
typedef enum md_lipid_tail_kind_t {
    MD_LIPID_TAIL_CHAIN     = 0,    // A hydrocarbon chain: acyl, alkyl, alkenyl or sphingoid, branched or not (phytanyl).
                                    // It hangs off a head, or spans two (the chains of a bolalipid).
    MD_LIPID_TAIL_SKELETON  = 1,    // The fused ring system of a sterol, with its side chain
} md_lipid_tail_kind_t;

// LIPIDS (md_lipid_infer): the anatomy of each lipid molecule, its heads and its tails.
//
// A tail is a hydrocarbon part of the molecule which hangs off the rest of it at one point: an acyl, alkyl or sphingoid
// chain (MD_LIPID_TAIL_CHAIN), or the fused ring skeleton of a sterol with its side chain (MD_LIPID_TAIL_SKELETON). A
// chain which hangs between two points spans them: the chains of a bolalipid, which join its two heads. A head is a
// connected part of what is not a tail, so a lipid has one head, a bolalipid two. For an atomistic lipid the tails are
// found from the bonds alone (carbons bonded to carbons and hydrogens only), whatever the atoms are named; for a
// coarse grained (Martini) lipid from the names of its beads.
//
// Each tail has a SPINE: its atoms (heavy atoms or beads) in order along it, from where it hangs off its head out to its
// tip (or, for a spanning chain, on to where it joins its other head). Methyl branches and the atoms of the rings off the
// shortest way through a skeleton are not on it. The spine is what a spline of the tail follows.
//
// A lipid is a molecule, residues bonded to each other wherever they are in the system (a lipid split into head and
// tails, as the AMBER lipids are, is one; so is a lipid A whose acyl chains are ligands of their own), all of them
// MD_RESIDUE_KIND_LIPID. Of a molecule of several residues which holds a lipid and more (a glycolipid or a
// lipopolysaccharide whose sugars are residues of their own, a PEG lipid whose monomers are), the lipid is the
// residues which hold the tails and those between them, and the residues bonded to these: the head group of a
// residue of its own, the first sugar of a glycolipid, the first Kdo of a lipid A. The rest of the molecule is attached
// to the lipid but no part of it, nor of its heads, and keeps its kind. Every lipid molecule of the system is here, also
// one whose heads and tails could not be found (a lipid by name in a system without bonds): it has none. The atoms of
// the tails: md_lipid_tail_atom_mask.
typedef struct md_lipid_data_t {
    size_t count;                       // Lipid molecules
    uint32_t* residue_offset;           // [count + 1] The residues of lipid i are residue_idx[residue_offset[i] .. residue_offset[i + 1]]
    int32_t*  residue_idx;              // Ascending within a lipid (not necessarily consecutive)
    uint32_t* head_offset;              // [count + 1] The heads of lipid i are head_offset[i] .. head_offset[i + 1]
    uint32_t* tail_offset;              // [count + 1] Its tails, tail_offset[i] .. tail_offset[i + 1]

    struct {
        size_t count;
        uint32_t* offset;               // [count + 1] The atoms of head h are atom_idx[offset[h] .. offset[h + 1]]
        md_atom_idx_t* atom_idx;        // The heavy atoms (beads) of each head, ascending
    } head;

    struct {
        size_t count;
        uint32_t* spine_offset;         // [count + 1] The spine of tail t is spine_idx[spine_offset[t] .. spine_offset[t + 1]]
        md_atom_idx_t* spine_idx;       // The spine of each tail, in order from its head
        uint32_t* atom_offset;          // [count + 1] The atoms of tail t are atom_idx[atom_offset[t] .. atom_offset[t + 1]]
        md_atom_idx_t* atom_idx;        // All atoms of each tail (branches and hydrogens too), ascending
        md_atom_idx_t* anchor_beg;      // [count] The atom of the head the tail hangs off (the carbonyl C of an ester, a GL bead)
        md_atom_idx_t* anchor_end;      // [count] The atom of the other head a spanning chain joins, -1 for a free tip
        uint32_t* head_beg;             // [count] The head (index into head) the tail hangs off
        uint32_t* head_end;             // [count] The head a spanning chain joins, UINT32_MAX for a free tip
        md_lipid_tail_kind_t* kind;     // [count]
    } tail;
} md_lipid_data_t;

#ifdef __cplusplus
extern "C" {
#endif

// Frees the lipids of a system (allocated from alloc, the system's) and leaves them empty
void md_lipid_data_free(md_lipid_data_t* lipid, struct md_allocator_i* alloc);

// The ranges of a lipid: NULL and 0 when the index is out of range
static inline const int32_t* md_lipid_residues(const md_lipid_data_t* lipid, size_t lipid_idx) {
    MD_ASSERT(lipid);
    return lipid_idx < lipid->count ? lipid->residue_idx + lipid->residue_offset[lipid_idx] : NULL;
}

static inline size_t md_lipid_residue_count(const md_lipid_data_t* lipid, size_t lipid_idx) {
    MD_ASSERT(lipid);
    return lipid_idx < lipid->count ? lipid->residue_offset[lipid_idx + 1] - lipid->residue_offset[lipid_idx] : 0;
}

static inline const md_atom_idx_t* md_lipid_head_atoms(const md_lipid_data_t* lipid, size_t head_idx) {
    MD_ASSERT(lipid);
    return head_idx < lipid->head.count ? lipid->head.atom_idx + lipid->head.offset[head_idx] : NULL;
}

static inline size_t md_lipid_head_atom_count(const md_lipid_data_t* lipid, size_t head_idx) {
    MD_ASSERT(lipid);
    return head_idx < lipid->head.count ? lipid->head.offset[head_idx + 1] - lipid->head.offset[head_idx] : 0;
}

static inline const md_atom_idx_t* md_lipid_tail_spine(const md_lipid_data_t* lipid, size_t tail_idx) {
    MD_ASSERT(lipid);
    return tail_idx < lipid->tail.count ? lipid->tail.spine_idx + lipid->tail.spine_offset[tail_idx] : NULL;
}

static inline size_t md_lipid_tail_spine_size(const md_lipid_data_t* lipid, size_t tail_idx) {
    MD_ASSERT(lipid);
    return tail_idx < lipid->tail.count ? lipid->tail.spine_offset[tail_idx + 1] - lipid->tail.spine_offset[tail_idx] : 0;
}

static inline const md_atom_idx_t* md_lipid_tail_atoms(const md_lipid_data_t* lipid, size_t tail_idx) {
    MD_ASSERT(lipid);
    return tail_idx < lipid->tail.count ? lipid->tail.atom_idx + lipid->tail.atom_offset[tail_idx] : NULL;
}

static inline size_t md_lipid_tail_atom_count(const md_lipid_data_t* lipid, size_t tail_idx) {
    MD_ASSERT(lipid);
    return tail_idx < lipid->tail.count ? lipid->tail.atom_offset[tail_idx + 1] - lipid->tail.atom_offset[tail_idx] : 0;
}

// ### INFERENCE ###

// The lipids of the system (md_lipid_data_t): every lipid molecule, its heads and tails. A molecule (a run of consecutive residues bonded to each other, none of them an amino
// acid, a nucleotide, water or an ion) whose residues are lipids by name gets the tails found in it. One whose
// residues are not is a lipid by its structure when it has a polar head and a hydrocarbon chain of at least 8 carbons
// or the skeleton of a sterol (atomistic molecules only): lipids of
// names nothing knows (the 3 letter names of PDB files, the codes of the PDB, lipids split into several residues).
// Tails are found from the bonds alone for atoms, and from the names of the beads (the conventions of the Martini
// lipidomes) for coarse grained lipids. Requires the classified residues, the bonds for the heads and tails and for
// the molecules of several residues (without them a lipid is its residue), and for sterols the rings.
// Into out (allocated with alloc, never &sys->lipid); the system is read and never written: the kinds of the residues
// follow from the lipids, md_system_infer_lipid_kinds.
bool md_lipid_infer(md_lipid_data_t* out, const struct md_system_t* sys, struct md_allocator_i* alloc);

// The kinds of the residues from the lipids of the system (sys->lipid): the residues of a lipid are
// MD_RESIDUE_KIND_LIPID (lipids of names nothing knows, found by their structure), and a residue of kind LIPID which is
// of no lipid now is MD_RESIDUE_KIND_OTHER, unless it is a lipid by its name (md_resname_lipid) or a coarse grained
// one (its kind is the force field's): a lipid by its structure alone stays one only while it is found to be one.
void md_system_infer_lipid_kinds(struct md_system_t* sys);

// The atoms of the tails of the lipids: out[i] is 1 for an atom of a tail, 0 for any other, of atom_count atoms.
// Returns how many there are.
size_t md_lipid_tail_atom_mask(uint8_t* out, const md_lipid_data_t* lipid, size_t atom_count);

#ifdef __cplusplus
}
#endif
