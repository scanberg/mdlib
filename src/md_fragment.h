#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <core/md_str.h>
#include <md_types.h>

// FRAGMENTS: the chemically meaningful pieces a molecule is made of - a ring system, a functional group, a backbone,
// an alkyl chain. The level at which a chemist says where an excitation, a charge or a contact sits ("the phenyl of
// TYR 66", "the nitrile"), and selects ("every carboxylate").
//
// A LEVEL OF THE SYSTEM
//
// Every atom with an element is in exactly one fragment, and every fragment is within one molecule. The fragments are
// inferred with the rest of the topology (md_util_system_infer, MD_UTIL_INFER_FRAGMENT_BIT, last: after the molecules,
// the residue flags and atom roles, the rings, the chemistry and the lipids, which they are made from), are DERIVED
// (md_system_fields.h) and covered by sys->topology_version, and are what a selection of fragment granularity grows to.
//
// Fragments and residues are both ways of cutting a molecule, and neither nests in the other: the peptide group joins
// two amino acids, the phosphate two nucleotides, an ester the residues of a lipid. So a fragment is of every residue
// it has atoms in, and a residue's fragments are those of its atoms - both found from the atoms
// (md_system_residue_fragments, md_system_fragment_residues), neither stored. Quantities aggregate by atoms (a
// residue's mass, charge or share of an excitation is over its atoms), never by adding up a residue's fragments,
// which counts such a fragment in each of its residues.
//
// WHAT A FRAGMENT IS
//
//   backbone        Of a polymer whose residues were resolved (their atoms have roles, MD_RESIDUE_FLAG_RESOLVED). Of a
//                   protein: the alpha carbon of each amino acid (CA alone, with its hydrogen), the peptide group
//                   between two (C, O and the next N: "peptide"),
//                   and the termini (an ammonium, a carboxylate). Of a nucleic acid: the sugar of each nucleotide (its
//                   ring and C5': "deoxyribose", "ribose") and the phosphate between two. Any group with an atom of
//                   MD_ATOM_FLAG_BACKBONE is of the backbone; the side chains and bases are fragmented as below.
//   lipid tail      The atoms of MD_ATOM_FLAG_LIPID_TAIL, one fragment per tail.
//   ring system     Rings (sys->ring) sharing an atom, with the terminal atoms double bonded to them (the O of a
//                   quinone or a lactam).
//   functional      Ertl's algorithm (J. Cheminform. 9:36, 2017): the heteroatoms, the carbons double or triple
//   group           bonded to a heteroatom, the carbons of C=C and C#C outside of rings, the acetal carbons and the
//                   carbons carrying two or more halogens are marked; each connected set of marked atoms is a group.
//   alkyl, linker   What is left of the carbon skeleton, as connected pieces. A piece of one or two carbons on one
//                   functional group, and nothing else but the backbone, is part of the group, as a chemist names them
//                   together (the methyls of a dimethylamino, the methyl of a methoxy, the CB of a serine's hydroxyl).
//                   Any other piece is an ALKYL (a methyl on a ring; the methyl, isopropyl, isobutyl and sec-butyl of
//                   alanine, valine, leucine, isoleucine) or a LINKER between two or more fragments (the methylene
//                   between a phenylalanine's alpha carbon and its ring). Neither a ring system nor the backbone takes
//                   a piece of the skeleton: what is drawn as a phenyl is the six carbons of its ring.
//   atom            Anything else on its own: a metal, an ion, a lone atom.
//
// A single halogen belongs to the fragment of the atom it is bonded to, every hydrogen to that of its heavy atom.
// Conjugated groups are not merged with the ring they hang from: 4-(dimethylamino)benzonitrile is an amine, a phenyl
// and a nitrile, which is what lets a donor - pi - acceptor read as three things. Coordination bonds join nothing.
//
// TYPES: what a fragment is called
//
// A fragment type is a kind and a name ("functional group", "amide"), shared by every fragment so called, as a residue
// type is a name and a kind: what a selection by name takes and a legend lists, and it carries a colour, which
// inferring again keeps. Names: the core (the ring system, or the marked atoms) identified as a whole against a
// library of SMILES (md_match_identify: phenyl, pyridyl, amide, nitro, ...); a core of one atom by its element and
// bonds (amine, hydroxyl, ether, water); the backbone by its part in it (alpha carbon, peptide, deoxyribose); else by
// kind and composition ("aryl C9N", "group C2NO", "C16 tail"). Two phenyls of one molecule are both "phenyl": a view
// that needs them apart numbers them.
//
// The kind is what a fragment is to its molecule, not its chemistry: a peptide group is an amide, but of kind backbone
// and named "peptide", where the amide of an asparagine is a functional group named "amide". So fragment("backbone")
// is the whole backbone and fragment("amide") the amides of the side chains.
//
// TEMPLATES: analysed once per kind of molecule, looked up for every molecule of that kind
//
// A molecule's fragments depend on nothing outside of it, so two molecules alike have the same fragments, atom for
// atom, computed once as a TEMPLATE. Alike is exact: the molecule in its slot order (md_molecule_data_t, a traversal
// from the molecule's graph centre) encoded as each atom's element, formal charge, hydrogen count, role flags and
// residue kind and flags, and each bond's two slots, order and aromaticity. The molecules are grouped by a hash of the
// encoding and the encodings compared in full, so a collision costs a comparison, never a wrong template. Two molecules
// alike but traversed differently are two templates: less sharing, never a wrong lookup. A solvent box is a template
// per solvent; a protein is one, shared by its copies.
//
// RINGS are an input, not a part of the result: a ring's fragment is that of any of its atoms
// (md_system_fragment_find_by_ring), and a consumer that needs fragment -> rings in bulk builds it in one pass.
//
// SELECTION
//
// A fragment is not a range of atoms (hydrogens are interleaved with heavy atoms, and slots are not atom order), so it
// is selected atom by atom: md_system_fragment_find_by_atom_idx, then md_fragment_atom_iter, each a lookup per atom.
// "The fragments of X, whole" - of a residue, of a selection - is the fragments of X's atoms: what a selection grown
// to fragment granularity is.
//
// INDICES
//
//   fragment            md_fragment_idx_t: the fragments of molecule m are [frag_offset[m], frag_offset[m + 1]), in
//                       the order of its template
//   fragment type       md_fragment_type_idx_t: a kind and a name, with a colour
//   local slot          An atom's position in its molecule: molecule.atom_slot[atom] - molecule.offset[m]
//
// COST
//
// Per molecule: its template and its first fragment (8 bytes). Per template: 8 bytes per atom of its representative
// (the fragment of each slot, the slots of each fragment) and 8 per fragment. A water box is one template and 8 bytes
// a water; a protein of 100k atoms 0.8 MB once, whatever its number of copies.

#define MD_FRAGMENT_NONE 0xFFFFFFFFu

typedef uint32_t md_fragment_idx_t;
typedef uint32_t md_fragment_type_idx_t;

typedef enum md_fragment_kind_t {
    MD_FRAGMENT_BACKBONE = 0,
    MD_FRAGMENT_LIPID_TAIL,
    MD_FRAGMENT_RING_SYSTEM,
    MD_FRAGMENT_GROUP,
    MD_FRAGMENT_LINKER,
    MD_FRAGMENT_ALKYL,
    MD_FRAGMENT_ATOM,
    MD_FRAGMENT_KIND_COUNT,
} md_fragment_kind_t;

// The types: one per kind and name in the system
typedef struct md_fragment_type_data_t {
    size_t    count;
    uint8_t*  kind;                     // [count] md_fragment_kind_t
    md_str_t* name;                     // [count]
    uint32_t* color;                    // [count] RGBA, as the residue types'
} md_fragment_type_data_t;

// sys->fragment. Every array an md_array of the system's allocator.
typedef struct md_fragment_data_t {
    size_t    count;                    // Fragments: frag_offset[molecule count]
    uint32_t* mol_template;             // [molecule count]
    uint32_t* frag_offset;              // [molecule count + 1]

    // TEMPLATES
    size_t    template_count;
    uint32_t* template_molecule;        // [template_count] Its representative
    uint32_t* template_slot_offset;     // [template_count + 1] Its slots in slot_frag
    uint32_t* slot_frag;                // The fragment of each slot of each template, local to it; MD_FRAGMENT_NONE for
                                        // an atom of none (no element)
    uint32_t* template_frag_offset;     // [template_count + 1] Its fragments in frag_type and frag_slot_offset
    md_fragment_type_idx_t* frag_type;  // The type of each fragment of each template
    uint32_t* frag_slot_offset;         // [+ 1] The slots of each fragment of each template in frag_slot
    uint32_t* frag_slot;                // Local slots, ascending

    md_fragment_type_data_t type;
} md_fragment_data_t;

typedef struct md_fragment_atom_iter_t {
    const uint32_t* slot;               // The fragment's local slots, still to come
    const uint32_t* end;
    const int32_t*  mol_atom;           // molecule.atom_idx from its molecule's first slot
    md_atom_idx_t   atom;               // Current atom
} md_fragment_atom_iter_t;

struct md_system_t;
struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

// Infers sys->fragment anew (replacing what it held, keeping the colours of the types it had) and bumps the topology
// version. Needs the molecules (md_util_system_infer_molecules, which clears the fragments: infer them again after
// it); without them the system has no fragments. Reads the bond orders and aromaticity (md_chem_perceive), without
// which every C=O reads as single, the rings, the residue flags and atom roles (md_util_system_infer_residue_flags)
// and the lipid tails. Done by md_util_system_infer with MD_UTIL_INFER_FRAGMENT_BIT.
bool md_util_system_infer_fragments(struct md_system_t* sys);
void md_fragment_data_free(md_fragment_data_t* data, struct md_allocator_i* alloc);

size_t md_system_fragment_count(const struct md_system_t* sys);

// The fragment of an atom, MD_FRAGMENT_NONE for none
md_fragment_idx_t md_system_fragment_find_by_atom_idx(const struct md_system_t* sys, size_t atom_idx);
// The fragment of a ring of sys->ring: that of its first atom
md_fragment_idx_t md_system_fragment_find_by_ring(const struct md_system_t* sys, size_t ring_idx);

size_t                 md_system_fragment_molecule(const struct md_system_t* sys, md_fragment_idx_t frag);
md_fragment_type_idx_t md_system_fragment_type(const struct md_system_t* sys, md_fragment_idx_t frag);
md_fragment_kind_t     md_system_fragment_kind(const struct md_system_t* sys, md_fragment_idx_t frag);
md_str_t               md_system_fragment_name(const struct md_system_t* sys, md_fragment_idx_t frag);
size_t                 md_system_fragment_atom_count(const struct md_system_t* sys, md_fragment_idx_t frag);

// The atoms of a fragment, in its molecule's slot order:
//   for (md_fragment_atom_iter_t it = md_fragment_atom_iter(sys, f); md_fragment_atom_iter_next(&it);) use(it.atom);
md_fragment_atom_iter_t md_fragment_atom_iter(const struct md_system_t* sys, md_fragment_idx_t frag);
bool md_fragment_atom_iter_next(md_fragment_atom_iter_t* it);

// Membership with the residues, found from the atoms: each writes at most cap, ascending and distinct, and returns how
// many there are
size_t md_system_residue_fragments(md_fragment_idx_t* out, size_t cap, const struct md_system_t* sys, size_t res_idx);
size_t md_system_fragment_residues(md_residue_idx_t* out, size_t cap, const struct md_system_t* sys, md_fragment_idx_t frag);

// Every fragment of a type, ascending: writes at most cap, returns how many there are
size_t md_system_fragment_find_by_type(md_fragment_idx_t* out, size_t cap, const struct md_system_t* sys, md_fragment_type_idx_t type);

size_t md_system_fragment_type_count(const struct md_system_t* sys);
md_fragment_kind_t md_fragment_type_kind(const struct md_system_t* sys, md_fragment_type_idx_t type);
md_str_t           md_fragment_type_name(const struct md_system_t* sys, md_fragment_type_idx_t type);
uint32_t           md_fragment_type_color(const struct md_system_t* sys, md_fragment_type_idx_t type);
// The type of a name (and kind, MD_FRAGMENT_KIND_COUNT for any), MD_FRAGMENT_NONE for none
md_fragment_type_idx_t md_fragment_type_find(const struct md_system_t* sys, md_str_t name, md_fragment_kind_t kind);

const char* md_fragment_kind_str(md_fragment_kind_t kind);   // "ring system", "functional group", ...

#ifdef __cplusplus
}
#endif
