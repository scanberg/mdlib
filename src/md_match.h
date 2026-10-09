#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <md_types.h>

struct md_allocator_i;
struct md_bitfield_t;
struct md_system_t;
struct md_smiles_error_t;

// ### STRUCTURE MATCHING ###
// Finds where a query structure occurs in the covalent bond graph of a system (substructure search).
//
// QUERY    md_match_query_t. A small labelled graph: atoms carry constraints (element, atom name, hydrogen count),
//          bonds an order. Plain data, made from
//            - a SMILES string                 md_match_query_init_smiles
//            - atoms of a system (reference)   md_match_query_init_atoms
//            - by hand, or SMARTS later
//          A query keeps nothing of the system it was made from, so a reference taken in one system can be looked
//          for in another.
//
// SEARCH   md_match_desc_t. Where to look (level, mask), what to report (mode) and how strictly (resolution).
//          md_match_for_each streams the matches to a callback, md_match_find collects them.
//
// RESULT   md_match_result_t. One row per match holding one atom index per query atom, in query atom order, so
//          column j of every row is the atom which plays the part of query atom j.
//
// LIBRARY  md_match_library_t. Many queries, each of a known structure, and md_match_identify, which tells for each
//          residue (molecule, chain) which of them it is. See IDENTIFICATION.
//
// ## WHAT A MATCH IS
//
// A match maps every query atom to its own atom of the system such that each satisfies the constraints of the query
// atom mapped onto it, and every query bond is a bond of the system (of a compatible order, where orders are tested).
// The system may have more: further atoms bonded to the matched ones, and bonds between matched atoms that the query
// does not state. A chain in the query therefore also matches along a ring. With MD_MATCH_FLAG_WHOLE the system may
// not have more: the match is then all of its unit (see LEVELS).
//
// ## THE GRAPH
//
// What is searched is the covalent bond graph of the system, which is not quite all of sys->bond:
//   - A coordination bond (MD_BOND_FLAG_COORDINATE, between a metal and a non-metal) between different components
//     is not part of either molecule: the Mg between the phosphates of an ATP, the Zn of a zinc finger. One within
//     a residue is part of it (the Fe of a heme, an iron-sulfur cluster). In a system without components every
//     coordination bond is left out.
//   - Virtual sites (MD_PARTICLE_VIRTUAL_SITE: the M site of a TIP4P water) are not atoms of the graph.
// A MOLECULE is a connected part of this graph. Molecules are numbered in the order of their first atom, as
// sys->structure numbers structures: in a system without coordination bonds and virtual sites, which is not coarse
// grained, they are its structures, index for index. Otherwise a structure may hold several (an ATP and its Mg are two
// molecules).
//
// ## WHAT THE SYSTEM CAN TELL
//
// A system describes its atoms less completely than a SMILES does. A crystal structure without hydrogens, a united
// atom model (hydrogens on carbon folded into the carbon) and an all atom simulation of the same molecule differ in
// what they state, and the same query should find the molecule in all three. md_chem_perceive (md_chem.h) works out
// most of what they leave implicit. What is still unknown is not tested, and what is known up to something is tested
// up to that something:
//
//   HYDROGENS    With hydrogen counts in the system (md_atom_data_t.hydrogen_count, explicit and implicit), tested on
//                every atom except those whose protonation is unknown (below). Without them, tested on an atom when
//                its molecule holds at least one atom of the same element with a hydrogen bonded to it, and against
//                the hydrogens bonded to it. All atom: tested everywhere. United atom: on N and O, not on C. No
//                hydrogens: nowhere. Per molecule, so a ligand without hydrogens in an all atom protein is still
//                found, and per element, which is what tells a united atom model apart.
//   PROTONATION  A crystal structure does not show where its protons are. The hydrogens and charge of N, O, P, S and
//                Se in a molecule without hydrogen atoms are md_chem's guess (MD_CHEM_FLAG_PROTONATE_PH7), and a
//                proton more or less is as likely. They are compared up to that: what is tested is the hydrogens minus
//                the charge, which a proton leaves as it is and the bonds alone decide (the valence the bonds leave).
//                C(=O)[OH] and C(=O)[O-] both find every Asp and Glu of a crystal structure, and neither finds an
//                ester; [O-] finds the hydroxyls as well, and [NH3+] the primary amides with the amines: by the heavy
//                atoms alone they are the same. Within a resonance group md_chem also chose where the double bonds
//                are, which may make the difference one. A query atom which states only one of the two is not tested
//                on either. With hydrogen atoms in the molecule the protonation is known and both are tested as they
//                are. To hold a crystal structure to md_chem's guess, test both ALWAYS (see below).
//   CHARGES      Tested when the system has formal charges (md_atom_data_t.formal_charge). Within a resonance group
//                (atoms joined by aromatic or delocalized bonds: a carboxylate, a guanidinium, an imidazolium) a
//                charge belongs to the group rather than to the atom md_chem put it on, and the charges a query
//                states for atoms of a group are compared as a sum: with the sum of the charges of the atoms they are
//                matched to, or, as the group's charge may sit on any of its atoms, with any share of it (from 0 to
//                all of it), but all of it when the query states a charge for every atom of the group.
//                NC(=[NH2+])N finds every arginine whichever N carries the charge, [O-] both oxygens of a
//                carboxylate, [O-]C=[O-] neither.
//   BOND ORDERS  Tested on every bond whose order is known: md_bond_order is not unknown, or the bond is aromatic or
//                delocalized.
//   AROMATICITY  Tested on atoms with at least one bond of known order. An atom is aromatic with MD_ATOM_FLAG_AROMATIC or
//                an aromatic bond.
//
// md_match_desc_t.hydrogens, .charges and .bond_orders override this per search. ALWAYS tests regardless and as the
// system has it: an atom without a hydrogen count has the hydrogens bonded to it, an atom without a charge has none, a
// bond of unknown order reads as single, an atom without bonds of known order as not aromatic. NEVER does not test
// (bond_orders covers aromaticity too). Resonance groups are compared as groups either way.
//
// ## HYDROGENS IN A QUERY
//
// Hydrogens are not searched for. A query hydrogen (MD_MATCH_ATOM_ELEMENT, z = 1) bonded to exactly one atom which is
// not a hydrogen is a property of that atom rather than an atom of its own: it requires the atom to carry at least
// that many hydrogens, and once the other atoms are matched it is mapped onto one of the system's hydrogens on the
// matched atom which satisfies its constraints (a name), lowest index first. Searching for them would multiply every
// match by the permutations of the hydrogens on each atom (3! = 6 per methyl group) and report the same atoms over and
// over. Where the matched atom has no hydrogen atom left to map onto (its hydrogens are implicit, or not resolved),
// such a query hydrogen is mapped to -1; where it has one but none satisfies the constraints, the match fails.
// Any other query hydrogen ([H][H], [H+], a bridging hydrogen) is searched for like any other atom.
// Hydrogen counts, and the hydrogens mapped this way, are those bonded to the atom whether or not they are within the
// unit and the mask: they belong to the atom. The count of a query atom (h_min, h_max) is all of its hydrogens, those
// given as query atoms included.
//
// ## LEVELS (md_match_level_t)
//
// Every match lies within one unit: a molecule (see THE GRAPH), a component (residue) or an instance (a polymer chain
// or a molecule, see md_instance_data_t). At the component and instance levels the bonds leaving the unit are not part
// of the graph that is searched: a residue pattern ends at the peptide bonds, and MD_MATCH_FLAG_WHOLE holds a residue
// against its own atoms only. A query made of disconnected parts (SMILES '.') has all of its parts matched within the
// same unit.
// The unit of a match (md_match_result_t.unit) is the index of its molecule, component or instance.
//
// MD_MATCH_FLAG_WHOLE: the match covers every atom of its unit, hydrogens aside, and the unit has no bonds between
// those atoms beyond those of the query. "Which residues are glycine" is NCC=O at the component level with this flag;
// without it every amino acid contains glycine.
// Hydrogens are left out of the count because whether a unit has hydrogen atoms is what the system may not tell (see
// above); a query which states hydrogen counts still has them tested. With a mask, a unit is its atoms within the mask.
//
// ## MODES (md_match_mode_t)
//
//   UNIQUE        One match per distinct set of matched atoms. A symmetric query matches the same atoms in several ways
//                 (benzene in 12), which are one match (see SYMMETRY). The default.
//   ALL           Every mapping, the symmetric ones included. Each mapping is reported once.
//   ONE_PER_UNIT  The first match found in each unit. "Is it in there", one row per residue, molecule or chain.
//   DISJOINT      Matches share no atoms; of overlapping ones, the one found first is kept. Counting repeat units of a
//                 polymer (glucose rings along a cellulose chain at the structure level) is this mode.
//
// Matches are reported unit by unit, in the order found, which is fixed for a given system and query. Only the
// searched atoms count towards distinct and overlapping, not the hydrogens mapped onto them.
//
// ## SMILES AS A QUERY
//
//   - Organic subset atoms (C, N, O, c, n, ...) leave the hydrogen count and the charge open, bracket atoms fix both:
//     [CH2] has exactly two hydrogens and no charge, [NH3+] is an ammonium. SMARTS reads them the same way, and it is
//     what lets a SMILES describe a fragment: CC(=O)O finds acids, carboxylates and esters alike, C(=O)[OH] only the
//     acids, [CH3]C(=O)[OH] acetic acid itself. Hydrogens written as atoms add to the count of a bracket atom:
//     N[C@@H](C)C=O and N[C@@]([H])(C)C=O are both one hydrogen.
//   - Lowercase atoms have to be aromatic, where aromaticity is known. Uppercase ones in a ring are not held to be
//     aliphatic, so a Kekule SMILES (C1=CC=CC=C1) finds an aromatic ring too. Uppercase atoms outside of the rings of
//     a SMILES which has lowercase atoms are held to be aliphatic: whoever wrote c1ccccc1 for the ring meant the C of
//     Cc1ccccc1 not to be aromatic, and it does not find a fused ring (Cc1cncn1 is a histidine side chain, not a
//     purine). A bond without a symbol is single, or aromatic between two aromatic atoms; see md_match_bond_order_t
//     for how orders compare.
//   - The atom class (:n) becomes the tag of the atom (md_match_atom_t.tag), which the search ignores: it names the
//     atoms of a pattern for whoever reads the result (OC(=O)[CH:1](C)[NH2:2]: the CA is tagged 1, the N 2).
//   - Isotope and chirality are read but not tested: md_system_t has no isotopes, and chirality would take
//     coordinates.
//   - '*' matches any atom which is not a hydrogen.
//
// ## REFERENCE ATOMS AS A QUERY
//
//   - The atoms are labelled by element, and with MD_MATCH_LABEL_NAME also by atom (type) name. Names only mean
//     something between systems which name their atoms the same way, and make the match far more selective there.
//   - The bonds are those of the graph (see THE GRAPH) between the given atoms, with their order where it is known; a
//     delocalized bond becomes MD_MATCH_BOND_AROMATIC, which takes either order. Aromatic atoms are required to be
//     aromatic.
//   - An atom whose hydrogens are all among the given atoms has its hydrogen count fixed, and its charge where the
//     system has charges; with some or none of them given both are open. Selecting the hydrogens of an atom pins its
//     protonation, leaving them out leaves it open. This holds where the reference system resolves hydrogens on the
//     atom (see WHAT THE SYSTEM CAN TELL), and only when the given atoms include hydrogens at all: atoms given without
//     any hydrogens are a skeleton, and a carbonyl carbon in it does not insist on having none.
//   - Atoms which are not bonded to each other make a query of several parts, matched within one unit: a chain with a
//     gap is two parts, each looked for anywhere in the chain.
//
// ## IDENTIFICATION
//
// md_match_identify answers "which of these is each residue" for a library of known structures (amino acids,
// nucleotides, ligands, as SMILES or references): for every unit, the first entry of the library, in the order added,
// which matches it. With MD_MATCH_FLAG_WHOLE (what identification usually means) the unit is that structure; without,
// the unit contains it. It is one pass over the units, and an entry is only tried on a unit whose heavy atom
// composition it fits (with WHOLE: the same elements, as many of each, and as many bonds), so the cost grows with the
// units rather than the entries: 49 entries (the amino acids and nucleotides) identify the 10687 residues of a system
// of 161742 atoms in about 10 ms, where a search per entry takes 15 times as long. Each identified unit comes with the
// mapping of the entry's atoms onto its own, which is what carries the names, orders or charges of a template over to
// the system.
//
// ## HOW THE SEARCH RUNS
//
// The query atoms are put in an order in which every atom after the first of its part is bonded to an earlier one,
// preferring the atom with the most bonds back into those already ordered (the RI ordering, Bonnici et al. 2013), so
// that ring closures are tested as early as possible, and among those breadth first, so that a wrong choice between
// the neighbours of an atom shows a few steps later rather than at the far end of a chain. The first atom is the one
// whose element is rarest among the candidate atoms of the system. Candidates for each later atom are the unmatched
// neighbours of the system atom its ordered parent was matched to, filtered by its constraints, by having at least as
// many bonds to heavy atoms (exactly as many with WHOLE), and by the bonds back to earlier atoms. A step costs in the
// order of the degree of the atoms involved, never the size of the structure, and the search is iterative, so
// references of thousands of atoms (a whole chain: 2 ms for 253 chains of 319 atoms) are no different from small
// ones. Each mapping is reached exactly once, from one start atom. The search runs unit by unit, and what is decided
// per unit (its size, for WHOLE) or per molecule (resolution) is worked out once.
// The search runs on one thread. Units are independent of each other, which is where it splits when it needs to.
//
// SYMMETRY. A query with symmetries fits the same atoms in several ways: a phenyl ring flipped over, the oxygens of a
// sulfonate in any order, and a reference of a chain with 40 such groups in 2^40 and more. Outside of ALL they are not
// all searched. Atoms with a single bond which hang from the same atom (the O of a sulfonate, the methyls of a valine)
// are searched last and mapped as a set onto its neighbours, and the symmetries of the rest are broken by conditions
// on the order of the atoms matched (Grochow and Kellis 2007), worked out from the query when the search begins: what
// is left to search is close to one mapping per set of atoms. Atoms which are alike take atoms of the system in their
// own order, so a reference given in the order of its atoms (as a script gives them) maps onto its own atoms as given.
//
// ## SMARTS
//
// The atom and bond constraints here are a conjunction of primitives, the part of SMARTS with only the '&' and ';'
// operators. The search tests an atom or a bond through one function each; SMARTS adds an expression per atom and
// per bond (or, not, recursion) behind those two functions, and a parser. Most of the primitives (element, aromatic,
// hydrogen count, charge, degree, bond order) can be answered from the system now; ring membership and ring sizes
// (R, r, x, @) need a per atom view of the rings of the system first.

// ### QUERY ###

typedef enum md_match_atom_flags_t {
    MD_MATCH_ATOM_NONE      = 0,
    MD_MATCH_ATOM_ELEMENT   = 0x1,      // z equal
    MD_MATCH_ATOM_NAME      = 0x2,      // Atom (type) name equal
    MD_MATCH_ATOM_HCOUNT    = 0x4,      // Hydrogens in [h_min, h_max]
    MD_MATCH_ATOM_CHARGE    = 0x8,      // Formal charge equal
    MD_MATCH_ATOM_AROMATIC  = 0x10,     // Aromatic
    MD_MATCH_ATOM_ALIPHATIC = 0x20,     // Not aromatic
} md_match_atom_flags_t;

typedef struct md_match_atom_t {
    uint32_t           flags;           // md_match_atom_flags_t: which of the fields below constrain the match
    md_atomic_number_t z;
    uint8_t            h_min;           // All hydrogens of the atom, those given as query atoms included
    uint8_t            h_max;
    int8_t             charge;
    md_label_t         name;
    int32_t            source;          // Where the atom came from: its character offset in a SMILES string, its index
                                        // in the reference system. -1 for none. Not used by the search.
    uint16_t           tag;             // The atom class of a SMILES, 0 for none. Not used by the search.
} md_match_atom_t;

// Where orders are tested: ANY matches every bond. SINGLE and DOUBLE match their own order and aromatic and
// delocalized bonds, AROMATIC matches single, double, aromatic and delocalized bonds, so Kekule, aromatic and
// resonance forms find each other (C(=O)O finds a carboxylate either way round). TRIPLE and QUADRUPLE match their own
// order only.
typedef enum md_match_bond_order_t {
    MD_MATCH_BOND_ANY       = 0,
    MD_MATCH_BOND_SINGLE    = 1,
    MD_MATCH_BOND_DOUBLE    = 2,
    MD_MATCH_BOND_TRIPLE    = 3,
    MD_MATCH_BOND_QUADRUPLE = 4,
    MD_MATCH_BOND_AROMATIC  = 5,
} md_match_bond_order_t;

typedef struct md_match_bond_t {
    uint32_t a;
    uint32_t b;
    uint32_t order;                     // md_match_bond_order_t
} md_match_bond_t;

// Plain data, may be filled by hand. The search validates it (indices in range, no bond from an atom to itself, no
// bond given twice) and derives everything else it needs per search.
typedef struct md_match_query_t {
    size_t           num_atoms;
    md_match_atom_t* atoms;
    size_t           num_bonds;
    md_match_bond_t* bonds;
    struct md_allocator_i* alloc;
} md_match_query_t;

// How md_match_query_init_atoms labels the atoms of a reference.
typedef enum md_match_label_t {
    MD_MATCH_LABEL_ELEMENT = 0,         // By element
    MD_MATCH_LABEL_NAME    = 1,         // By element and atom (type) name
} md_match_label_t;

// ### SEARCH ###

typedef enum md_match_level_t {
    MD_MATCH_LEVEL_STRUCTURE = 0,       // Within one molecule (see THE GRAPH)
    MD_MATCH_LEVEL_COMPONENT = 1,       // Within one component (residue)
    MD_MATCH_LEVEL_INSTANCE  = 2,       // Within one instance (polymer chain or molecule)
} md_match_level_t;

typedef enum md_match_mode_t {
    MD_MATCH_MODE_UNIQUE       = 0,
    MD_MATCH_MODE_ALL          = 1,
    MD_MATCH_MODE_ONE_PER_UNIT = 2,
    MD_MATCH_MODE_DISJOINT     = 3,
} md_match_mode_t;

typedef enum md_match_flags_t {
    MD_MATCH_FLAG_NONE  = 0,
    MD_MATCH_FLAG_WHOLE = 0x1,          // The match is all of its unit (see LEVELS)
} md_match_flags_t;

typedef enum md_match_resolve_t {
    MD_MATCH_RESOLVE_AUTO   = 0,        // Tested where the system carries the information (see WHAT THE SYSTEM CAN TELL)
    MD_MATCH_RESOLVE_ALWAYS = 1,        // Always tested, as stated: a system without the information then never matches
    MD_MATCH_RESOLVE_NEVER  = 2,        // Never tested
} md_match_resolve_t;

typedef struct md_match_desc_t {
    const md_match_query_t* query;

    md_match_level_t   level;
    md_match_mode_t    mode;
    uint32_t           flags;           // md_match_flags_t

    md_match_resolve_t hydrogens;
    md_match_resolve_t bond_orders;
    md_match_resolve_t charges;

    // Optional. Only these atoms are searched: the graph is the one the bonds form between them, as for a unit.
    // This is what a script's 'in' context is.
    const struct md_bitfield_t* mask;

    // Optional, 0 for no limit. The search stops after this many matches.
    size_t limit;
} md_match_desc_t;

// ### RESULT ###

typedef struct md_match_result_t {
    size_t    count;                    // Number of matches
    size_t    width;                    // Atoms per match: the number of atoms of the query
    int32_t*  atom_idx;                 // [count * width]: the atom matched to query atom j in match i is
                                        // atom_idx[i * width + j]; -1 for a query hydrogen the system does not resolve
    uint32_t* unit;                     // [count]: index of the molecule, component or instance of each match
    struct md_allocator_i* alloc;
} md_match_result_t;

// Receives each match as it is found, atom_idx as one row of md_match_result_t. Valid during the call only.
// Return false to stop the search.
typedef bool (*md_match_callback_t)(const int32_t* atom_idx, size_t width, uint32_t unit, void* user_param);

typedef enum md_match_select_flags_t {
    MD_MATCH_SELECT_NONE      = 0,
    MD_MATCH_SELECT_HYDROGENS = 0x1,    // Also the hydrogens bonded to the matched atoms which the query did not map
} md_match_select_flags_t;

// ### IDENTIFICATION ###

// Entries are copies of the queries given, kept in the order added
typedef struct md_match_library_t md_match_library_t;

typedef struct md_match_identify_desc_t {
    const md_match_library_t* library;

    md_match_level_t   level;
    uint32_t           flags;           // md_match_flags_t, MD_MATCH_FLAG_WHOLE to identify rather than to find within

    md_match_resolve_t hydrogens;
    md_match_resolve_t bond_orders;
    md_match_resolve_t charges;

    const struct md_bitfield_t* mask;   // Optional, as for md_match_desc_t
} md_match_identify_desc_t;

typedef struct md_match_identify_result_t {
    size_t    count;                    // Units identified
    uint32_t* unit;                     // [count] Index of each unit identified, ascending
    uint32_t* entry;                    // [count] The library entry it was identified as
    uint32_t* offset;                   // [count + 1] The match of unit i is atom_idx[offset[i] .. offset[i + 1]),
    int32_t*  atom_idx;                 // in the atom order of its entry's query, as a row of md_match_result_t
    struct md_allocator_i* alloc;
} md_match_identify_result_t;

#ifdef __cplusplus
extern "C" {
#endif

// Query from SMILES. Returns false on a syntax error, described in err (optional).
bool md_match_query_init_smiles(md_match_query_t* query, str_t smiles, struct md_allocator_i* alloc, struct md_smiles_error_t* err);

// Query from atoms of a system, in the order given: the order of the columns of the result.
// The atoms need not be connected; a disconnected set is a query of several parts.
bool md_match_query_init_atoms(md_match_query_t* query, const int32_t* atom_idx, size_t count, md_match_label_t label, const struct md_system_t* sys, struct md_allocator_i* alloc);

void md_match_query_free(md_match_query_t* query);

// Streams the matches to callback, which may be NULL to only count them (UNIQUE and DISJOINT still keep the atoms
// seen so far). The number of matches goes to out_count (optional). Returns false if the query or the description is
// invalid (logged), true otherwise, also when nothing is found.
bool md_match_for_each(size_t* out_count, const md_match_desc_t* desc, const struct md_system_t* sys, md_match_callback_t callback, void* user_param);

// Collects the matches. Returns false if the query or the description is invalid (logged), true otherwise, also when
// nothing is found.
bool md_match_find(md_match_result_t* out, const md_match_desc_t* desc, const struct md_system_t* sys, struct md_allocator_i* alloc);

void md_match_result_free(md_match_result_t* result);

// One selection per match, written to out[0 .. min(count, cap)), which must be initialized bitfields.
// Returns the number written. flags is a combination of md_match_select_flags_t.
size_t md_match_result_select(struct md_bitfield_t* out, size_t cap, const md_match_result_t* result, const struct md_system_t* sys, uint32_t flags);

// All of the matches in one selection: out, an initialized bitfield, is cleared and set to their union.
void md_match_result_select_all(struct md_bitfield_t* out, const md_match_result_t* result, const struct md_system_t* sys, uint32_t flags);

// ### LIBRARY ###

md_match_library_t* md_match_library_create(struct md_allocator_i* alloc);
void md_match_library_destroy(md_match_library_t* library);

// Adds a copy of the query under a name (optional, copied). Returns the index of the entry, or -1 if the query is
// invalid (logged).
int32_t md_match_library_add(md_match_library_t* library, const md_match_query_t* query, str_t name);

// Adds the query of a SMILES. Returns the index of the entry, or -1 on a syntax error, described in err (optional).
int32_t md_match_library_add_smiles(md_match_library_t* library, str_t smiles, str_t name, struct md_smiles_error_t* err);

size_t md_match_library_count(const md_match_library_t* library);
const md_match_query_t* md_match_library_query(const md_match_library_t* library, size_t entry);
str_t md_match_library_name(const md_match_library_t* library, size_t entry);

// For every unit, the first entry of the library which matches it (see IDENTIFICATION). Returns false if the
// description is invalid (logged), true otherwise, also when nothing is identified.
bool md_match_identify(md_match_identify_result_t* out, const md_match_identify_desc_t* desc, const struct md_system_t* sys, struct md_allocator_i* alloc);

void md_match_identify_result_free(md_match_identify_result_t* result);

#ifdef __cplusplus
}
#endif
