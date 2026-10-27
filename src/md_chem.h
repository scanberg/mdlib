#pragma once

#include <md_types.h>

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_system_t;
struct md_system_state_t;
struct md_bond_delta_t;
struct md_allocator_i;

// ### CHEMISTRY PERCEPTION ###
//
// From the elements, the covalent bonds and the geometry of a reference state, md_chem_perceive works out what the
// connectivity leaves implicit:
//   - the hydrogens of each heavy atom, explicit and implicit (md_atom_data_t.hydrogen_count)
//   - the order of each covalent bond, one Kekule structure (md_bond_order)
//   - aromatic rings (MD_ATOM_FLAG_AROMATIC on atoms, MD_BOND_FLAG_AROMATIC on bonds)
//   - delocalized groups outside rings, whose formal double bond and charge could sit on any of their atoms:
//     carboxylate, nitro, phosphate, sulfonate, guanidinium and amidinium, azide (MD_BOND_FLAG_DELOCALIZED)
//   - formal charges (md_atom_data_t.formal_charge)
//   - hybridization (md_atom_hybridization)
// Bonds whose order it sets get MD_BOND_FLAG_ORDER_PERCEIVED. Orders a file gives (known and without that flag) are
// kept and the rest is fitted around them, and so are the nonzero formal charges the file gives: the
// atom/formal_charge column of sys->attributes, which the loaders publish (mmCIF).
// Atoms outside of organic chemistry take no part: metals, coarse grained beads and virtual sites; bonds to metals
// and coordination bonds (MD_BOND_KIND_COORDINATE) do not count towards valences, and elastic bonds
// (MD_BOND_KIND_ELASTIC) are no bonds here at all. Running it again gives the same result.
//
// HYDROGENS. Each residue (or molecule, without residues) is handled by what it carries:
//   - all hydrogens (some on carbon): every hydrogen is explicit, and orders and charges follow from the valences.
//     This is exact up to resonance.
//   - polar hydrogens only (united atom force fields): the hydrogens of carbons are implicit, from the geometry
//   - no hydrogens (most crystal structures): every hydrogen is implicit. The geometry decides between single and
//     double bonds (lengths relative to the covalent radii, angles), the tautomer where it does not: amide, lactam
//     and thione over imidic acid and thiol, amino over imino. Atoms bound to a metal are deprotonated instead (the
//     cysteines of a zinc finger, the N of a porphyrin). Protonation follows MD_CHEM_FLAG_PROTONATE_PH7.
//
// METHOD. Each atom wants the pi bonds its valence leaves after its sigma bonds; these are placed by a maximum
// matching over the bonds between such atoms (Edmonds), which prefers the shorter bonds. What the matching cannot
// place is resolved the way chemistry does: a cation (the N+ of nitro, N-oxides, pyridinium, guanidinium), an anion
// (carboxylate, phenolate) or, without hydrogens, a proton.
//
// AROMATICITY. Rings whose pi electrons count 4n+2, alone or fused with a neighbour, are aromatic, counted as RDKit
// counts them: the model SMILES are written in, so that a lowercase atom of a SMILES pattern (md_match.h) finds its
// atom. A ring atom with an exocyclic double bond to O, N or S (the C=O of a pyridone, uracil, guanine or flavin)
// gives the ring no electron but does not break it; any other exocyclic double bond does. The PDB chemical component
// dictionary differs there, and leaves those rings non-aromatic.
//
// Validated against the PDB chemical component dictionary (8500 components with their ideal coordinates, organic
// elements), with the aromaticity RDKit perceives for each: orders, aromaticity and charges agree for 99.4% of them
// with hydrogens and 98% from the heavy atoms alone. Most of the rest are macrocycles and 7 membered rings (rings are
// perceived up to 6 atoms, so a porphyrin is aromatic as pyrroles), charges missing in the dictionary (boron,
// quaternary atoms) and protonation the geometry cannot tell.

typedef enum md_chem_flags_t {
    MD_CHEM_FLAG_NONE           = 0,
    // Residues without hydrogens: the protonation states at pH 7. Aliphatic amines and acyclic amidines and
    // guanidines cationic, carboxylic, phosphoric and sulfonic acids anionic, histidine neutral. Without it they are
    // left neutral.
    MD_CHEM_FLAG_PROTONATE_PH7  = 1u << 0,
} md_chem_flags_t;

#ifdef __cplusplus
extern "C" {
#endif

// reference is optional but needed for anything without hydrogens: the geometry tells sp2 from sp3 and single from
// double bonds. The rings of the system are used (md_util_system_infer_rings) and inferred if missing.
bool md_chem_perceive(struct md_system_t* sys, const struct md_system_state_t* reference, uint32_t flags);

// ### A FRAME WITH OTHER BONDS ###
//
// A state whose bonds are not those of the system (bonds formed and broken in it: md_bond_delta_t, md_util_bond_delta)
// has a chemistry of its own where they changed. md_chem_perceive_delta works it out for the REGION the changes
// reach, as md_chem_perceive would for the bonds of the state, and leaves the system as it is: what lies outside of
// the region keeps what md_chem_perceive gave it, as nothing there can change.
//
// THE REGION. The orders come from a matching of each conjugated network on its own, and a bond formed or broken
// changes the sigma bonds of its two atoms alone. So the region is the atoms of the bonds that changed, grown in the
// bonds of the state through every atom which can take part in a network: an atom with a pi bond, a charge or a
// valence it does not fill in the system, and those a leftover can be resolved on (an N with three sigma bonds, an O
// or S in a ring: N-oxides, guanidinium, pyridinium). Then by every atom of the rings of 4 to 7 atoms it reaches,
// whose aromaticity it decides, and so on until it grows no more. The atoms bonded to it are its edge: they keep
// their state in the system and want no pi bond, which they did not want there either. Each atom keeps the hydrogen
// mode of its residue (or molecule) in the system (HYDROGENS), as that is how the file was written. An atom without
// any bond in the state (a bromide leaving) is a monatomic ion.
// A region of more than MD_CHEM_DELTA_MAX_ATOMS atoms (a reaction on a graphene sheet) is not perceived.
//
// The state's coordinates are the geometry (implicit hydrogens and the evidence of the bond lengths, as for
// md_chem_perceive). A state between two shapes (a transition, while a bond formed and the one it replaces are both
// there) can give an atom more sigma bonds than any valence: such an atom has no pi bond, and is charged only where
// its valence says so (N+ with four, O+ with three).

#define MD_CHEM_DELTA_MAX_ATOMS 4096

typedef struct md_chem_delta_t {
    size_t           num_atoms;
    int32_t*         atom_idx;          // [num_atoms] The atoms of the region, ascending
    md_atom_flags_t* atom_flags;        // [num_atoms] Their flags in the system with the hybridization and MD_ATOM_FLAG_AROMATIC of the state
    int8_t*          formal_charge;     // [num_atoms]
    uint8_t*         hydrogen_count;    // [num_atoms] Explicit and implicit
    size_t           num_bonds;
    uint32_t*        bond_idx;          // [num_bonds] The covalent bonds with an atom in the region, ascending: one of the system by its index,
                                        // one formed by sys->bond.count + its index in the delta
    md_bond_flags_t* bond_flags;        // [num_bonds] Their flags with the order, aromaticity and delocalization of the state
    md_ring_data_t   ring;              // The rings with an atom in the region, in the state (md_util_infer_rings)
} md_chem_delta_t;

// The chemistry of state (the coordinates of sys's atoms) where its bonds are those of sys changed by delta, see A
// FRAME WITH OTHER BONDS, with the flags sys was perceived with (md_chem_flags_t). out is freed and filled, allocated
// with alloc; empty when nothing changed. Returns false when it is not perceived: the region is too large, or sys has
// no perceived chemistry (md_chem_perceive).
bool md_chem_perceive_delta(md_chem_delta_t* out, const struct md_system_t* sys, const struct md_system_state_t* state, const struct md_bond_delta_t* delta, uint32_t flags, struct md_allocator_i* alloc);
void md_chem_delta_free(md_chem_delta_t* chem, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
