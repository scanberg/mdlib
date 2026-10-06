#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_system_t;
struct md_system_state_t;

// ### CHEMISTRY PERCEPTION ###
//
// From the elements, the covalent bonds and the geometry of a reference state, md_chem_perceive works out what the
// connectivity leaves implicit:
//   - the hydrogens of each heavy atom, explicit and implicit (md_atom_data_t.hydrogen_count)
//   - the order of each covalent bond, one Kekule structure (md_bond_order)
//   - aromatic rings (MD_FLAG_AROMATIC on atoms, MD_BOND_FLAG_AROMATIC on bonds)
//   - delocalized groups outside rings, whose formal double bond and charge could sit on any of their atoms:
//     carboxylate, nitro, phosphate, sulfonate, guanidinium and amidinium, azide (MD_BOND_FLAG_DELOCALIZED)
//   - formal charges (md_atom_data_t.formal_charge)
//   - hybridization (MD_FLAG_SP, MD_FLAG_SP2, MD_FLAG_SP3)
// Bonds whose order it sets get MD_BOND_FLAG_ORDER_PERCEIVED. Orders a file gives (known and without that flag) are
// kept and the rest is fitted around them, and so are the nonzero formal charges the file gives: the
// atom/formal_charge column of sys->attributes, which the loaders publish (mmCIF).
// Atoms outside of organic chemistry take no part: metals, coarse grained beads and virtual sites; bonds to metals
// and coordination bonds (MD_BOND_FLAG_COORDINATE) do not count towards valences. Running it again gives the same result.
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

#ifdef __cplusplus
}
#endif
