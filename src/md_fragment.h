#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <core/md_str.h>

// Chemically meaningful fragments of a system, as a named tree: the level at which a chemist talks about
// where an excitation, a charge or a contact sits. "The amine", "the phenyl ring", "TYR 66".
//
// THE TREE
//
//   molecule set   identical small molecules together ("water x24"), when there are several of them
//   molecule       a connected part of the COVALENT bond graph (coordination bonds do not join molecules,
//                  so the metal and each ligand of a complex are molecules of their own)
//   unit           what a molecule is built of (below); a molecule of one unit is that unit and has no
//                  children
//
// Every atom of the system with an element lies in exactly one leaf (md_fragments_t.atom_fragment). The
// fragments are in pre-order, a parent before its children, and their atoms are ordered so that every
// fragment is one contiguous range of md_fragments_t.atom_idx, its children's ranges inside its own.
//
// UNITS
//
//   residue         An amino acid or a nucleotide of the system's residues (md_residue_kind_t) is one unit,
//                   named by its residue name and number ("TYR 66").
//   metal           Every metal atom is a unit of its own.
//   ring system     Rings (sys->ring, 3 to 7 atoms, covalent bonds only) sharing an atom are one ring system:
//                   naphthalene and indole are one unit, biphenyl is two. Terminal atoms double bonded to a
//                   ring atom (the O of a quinone or a lactam, the =CH2 of a methylene ring) belong to it.
//   functional      Ertl's algorithm (J. Cheminform. 9:36, 2017), which needs no list of known groups: the
//   group           heteroatoms, the carbons double or triple bonded to a heteroatom, the carbons of C=C and
//                   C#C outside of rings, the acetal carbons (sp3, single bonds to two or more of O, N, S) and
//                   the carbons carrying two or more halogens (with those halogens) are marked, and each
//                   connected set of marked atoms is a group.
//   alkyl, linker   What is left of the carbon skeleton, as connected pieces. A piece of one or two carbons
//                   which touches exactly one unit is part of that unit (the methyls of a dimethylamine, the
//                   ethyl of an ethyl ester); larger ones are units of their own, a LINKER between two or more
//                   units and an ALKYL otherwise.
//   atom            Anything else on its own (a halide ion, a lone hydrogen).
//
// A single halogen belongs to the unit of the atom it is bonded to (chlorobenzene is a phenyl), and every
// hydrogen to the unit of its heavy atom. Conjugated groups are NOT merged with the ring they hang from:
// nitrobenzene is a phenyl and a nitro, 4-(dimethylamino)benzonitrile an amine, a phenyl and a nitrile -
// which is what lets a donor - pi - acceptor read as three things.
//
// NAMES
//
// A unit is named by its core (the ring system or the marked atoms, without the alkyl pieces it absorbed),
// identified as a WHOLE against a library of SMILES (md_match_identify): phenyl, pyridyl, naphthyl, amide,
// ester, nitro, nitrile, ... A core of one atom is named by its element and its bonds (amine, hydroxyl,
// ether, thiol); one the library does not know by its kind and composition ("aryl C9N", "group C2NO").
// Names repeated within a molecule are numbered ("phenyl 1", "phenyl 2"). A molecule of several units is
// named "protein" or "nucleic acid" when it holds residues of one, by its formula otherwise ("C10H15N").
//
// What it needs: the elements and the covalent bonds; for the functional groups the bond orders and the
// aromaticity (md_chem_perceive), without which every C=C and C=O reads as single; and the rings
// (md_util_system_infer_rings), which are computed here when the system has none.

struct md_system_t;
struct md_allocator_i;

#define MD_FRAGMENT_NONE     0xFFFFFFFFu
#define MD_FRAGMENT_NAME_CAP 32

typedef enum md_fragment_kind_t {
    MD_FRAGMENT_MOLECULE_SET = 0,
    MD_FRAGMENT_MOLECULE,
    MD_FRAGMENT_RESIDUE,
    MD_FRAGMENT_METAL,
    MD_FRAGMENT_RING_SYSTEM,
    MD_FRAGMENT_GROUP,
    MD_FRAGMENT_LINKER,
    MD_FRAGMENT_ALKYL,
    MD_FRAGMENT_ATOM,
    MD_FRAGMENT_KIND_COUNT,
} md_fragment_kind_t;

typedef struct md_fragments_t {
    size_t    count;
    uint32_t* parent;           // [count] MD_FRAGMENT_NONE for a root; a parent comes before its children
    uint8_t*  kind;             // [count] md_fragment_kind_t. A molecule of one unit has the kind of that unit
    char*     name;             // [count * MD_FRAGMENT_NAME_CAP] zero terminated
    uint32_t* atom_beg;         // [count] the fragment's atoms are atom_idx[atom_beg .. atom_end)
    uint32_t* atom_end;         // [count]
    int32_t*  atom_idx;         // [num_atom_idx]
    size_t    num_atom_idx;
    uint32_t* atom_fragment;    // [num_atoms] the leaf holding each atom, MD_FRAGMENT_NONE for none (no element)
    size_t    num_atoms;
    struct md_allocator_i* alloc;
} md_fragments_t;

#ifdef __cplusplus
extern "C" {
#endif

// Computes the fragments of the system. Returns false only for invalid input; a system without bonds is
// a set of single atoms.
bool md_fragments_compute(md_fragments_t* out, const struct md_system_t* sys, struct md_allocator_i* alloc);
void md_fragments_free(md_fragments_t* frags);

static inline md_str_t md_fragment_name(const md_fragments_t* frags, size_t idx) {
    return idx < frags->count ? md_str_from_cstr(frags->name + idx * MD_FRAGMENT_NAME_CAP) : (md_str_t){0};
}

static inline bool md_fragment_is_leaf(const md_fragments_t* frags, size_t idx) {
    // Children follow their parent directly, so a fragment is a leaf when the next one is not its child
    return idx < frags->count && (idx + 1 >= frags->count || frags->parent[idx + 1] != (uint32_t)idx);
}

const char* md_fragment_kind_str(md_fragment_kind_t kind);   // "ring system", "functional group", ...

#ifdef __cplusplus
}
#endif
