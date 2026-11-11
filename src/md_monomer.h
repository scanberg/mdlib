#pragma once

#include <stdint.h>

#include <core/md_str.h>

// The standard monomers of proteins and nucleic acids: the predefined table residue types are matched against, as
// particle types are matched against the elements. A monomer has its id in the Chemical Component Dictionary of the
// PDB, a one letter code and a default color (md_monomer_*). The values are not stable across versions: store a monomer
// by its id (md_monomer_id, md_monomer_from_id), not by its value.
//
// The functions which take or give a kind of residue are declared with it in md_hierarchy.h: md_monomer_kind and
// md_monomer_from_residue_name, which gives the monomer a residue type is a form of.
typedef enum md_monomer_t {
    MD_MONOMER_NONE = 0,    // Not a form of a standard monomer: a ligand, water, an ion, a lipid, ...

    // Amino acids
    MD_MONOMER_ALA, MD_MONOMER_ARG, MD_MONOMER_ASN, MD_MONOMER_ASP, MD_MONOMER_CYS,
    MD_MONOMER_GLN, MD_MONOMER_GLU, MD_MONOMER_GLY, MD_MONOMER_HIS, MD_MONOMER_ILE,
    MD_MONOMER_LEU, MD_MONOMER_LYS, MD_MONOMER_MET, MD_MONOMER_PHE, MD_MONOMER_PRO,
    MD_MONOMER_SER, MD_MONOMER_THR, MD_MONOMER_TRP, MD_MONOMER_TYR, MD_MONOMER_VAL,
    MD_MONOMER_SEC,         // Selenocysteine (U)
    MD_MONOMER_PYL,         // Pyrrolysine (O)
    MD_MONOMER_UNK,         // An amino acid which is none of the above (a modified residue of no known parent), or a
                            // cap of a chain (ACE, NME)

    // Deoxyribonucleotides
    MD_MONOMER_DA, MD_MONOMER_DC, MD_MONOMER_DG, MD_MONOMER_DT, MD_MONOMER_DU, MD_MONOMER_DI,
    // Ribonucleotides
    MD_MONOMER_A, MD_MONOMER_C, MD_MONOMER_G, MD_MONOMER_U, MD_MONOMER_I,
    MD_MONOMER_N,           // A nucleotide which is none of the above

    MD_MONOMER_COUNT
} md_monomer_t;

#ifdef __cplusplus
extern "C" {
#endif

// Standard monomers (md_monomer_t), "" / 0 for MD_MONOMER_NONE or a value out of range
md_str_t    md_monomer_id(md_monomer_t monomer);          // Id in the Chemical Component Dictionary: "ALA", "DA", "A"
md_str_t    md_monomer_name(md_monomer_t monomer);        // "Alanine", "Deoxyadenosine", "Adenosine"
char     md_monomer_letter(md_monomer_t monomer);      // One letter code: 'A' for ALA, DA and A
uint32_t md_monomer_color(md_monomer_t monomer);       // Default color (RGBA, R in the lowest byte, as the CPK colors)

// The monomer of the given id ("HIS"), exactly and ignoring case. MD_MONOMER_NONE if there is none.
md_monomer_t md_monomer_from_id(md_str_t id);

#ifdef __cplusplus
}
#endif
