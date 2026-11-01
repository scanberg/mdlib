#pragma once

// ### MAXIMUM COMMON SUBSTRUCTURE ###
// Which atoms of one structure play the part of which atoms of another, for structures that are not the same: a
// mutant and its wild type, GTP and GDP, the rings of PHE and TYR. What it gives is a correspondence, pairs of atoms,
// which a superposition, an RMSD or an SDF over different kinds of structure (VIAL's sdf) can be computed over.
//
// WHAT IS MATCHED. The heavy atoms of the two sets of atoms, each atom labelled by its element (or by its element and
// name, MD_MCS_FLAG_BY_NAME), and the bonds between them (sys->bond, the graph of the persistent bonds), whatever
// their order: inferred orders differ between structures that are alike. Hydrogens are left out: they add nothing to
// a superposition, and every methyl group would multiply the ways of matching by six.
//
// TWO WAYS, chosen by what the sets are:
//
//   GRAPH        The largest connected set of heavy atoms of a and of b whose bonds are the same both ways round (a
//                maximum common connected induced subgraph), found exactly by McSplit (McCreesh, Prosser and
//                Trimble, IJCAI 2017). What it costs grows with the size of what is common and how symmetric it is:
//                microseconds for residues, a fraction of a millisecond for nucleotides and ligands of 60 heavy atoms
//                (GTP and GDP: 0.2 ms), and beyond a few hundred atoms the search may not finish, so it gives up after
//                max_steps (incomplete) with the largest it found. More than MD_MCS_GRAPH_MAX heavy atoms is not
//                searched at all.
//   BY RESIDUE   Polymers: when both sets span two residues or more, all of them amino acids or nucleotides. The
//                residues are aligned by name, in order (Needleman-Wunsch: a residue of the same name scores 4, of
//                another -1, a gap -3), and each aligned pair is matched by atom name when the two are of one name,
//                by GRAPH when not (once per pair of names). Linear in the residues: α against β tubulin, 433
//                residues each, in 2.4 ms. MD_MCS_FLAG_GRAPH matches polymers as one graph instead.
//
// SYMMETRY. The bonds alone do not tell the oxygens of a phosphate apart, or the two sides of a phenyl ring, and the
// mapping the search gives first is any of them. With positions for both (md_mcs_side_t.xyz), every mapping of a's
// common atoms onto b's that the bonds allow is tried (md_match, ALL), and the one that superposes with the least
// RMSD is kept: GTP onto GDP fits at 0.29 Å, the first mapping found at 1.00 Å. Without positions it is the first.
// BY RESIDUE takes the atom names as they are.
//
// HOW MUCH IS COMMON. Two structures always have something in common, an atom of carbon at least: GTP and vinblastine
// share 13 heavy atoms. Whether that is enough is the caller's to say; md_mcs_result_t has the heavy atoms of both.

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <core/md_vec_math.h>

struct md_system_t;
struct md_allocator_i;

#define MD_MCS_GRAPH_MAX 1024

typedef enum md_mcs_flags_t {
    MD_MCS_FLAG_NONE    = 0,
    MD_MCS_FLAG_BY_NAME = 0x1,      // Atoms alike by element and name, not by element alone
    MD_MCS_FLAG_GRAPH   = 0x2,      // Polymers as one graph, not residue by residue
} md_mcs_flags_t;

typedef struct md_mcs_side_t {
    const struct md_system_t* sys;
    const int32_t*  atoms;          // Atoms of sys, in any order; hydrogens among them are not matched
    size_t          count;
    const vec3_t*   xyz;            // Optional: the positions of every atom of sys, to choose among symmetric mappings
} md_mcs_side_t;

typedef struct md_mcs_desc_t {
    md_mcs_side_t   a;
    md_mcs_side_t   b;
    uint32_t        flags;          // md_mcs_flags_t
    uint64_t        max_steps;      // Steps of the graph search before it gives up; 0 for 10^7
} md_mcs_desc_t;

typedef struct md_mcs_result_t {
    size_t      count;              // Pairs
    int32_t*    a;                  // [count] Atoms of a, ascending
    int32_t*    b;                  // [count] The atom of b that plays the part of each
    size_t      a_heavy;            // Heavy atoms of each set: what count is out of
    size_t      b_heavy;
    bool        by_residue;         // Matched residue by residue
    bool        incomplete;         // The graph search gave up at max_steps: count is the largest found by then
    double      rmsd;               // Å, of the pairs superposed, with positions for both; NaN without
    struct md_allocator_i* alloc;
} md_mcs_result_t;

#ifdef __cplusplus
extern "C" {
#endif

// The correspondence of the atoms of desc->a and desc->b. False when the description is not valid or a set is too
// large for the graph search (logged); true otherwise, also when nothing is common.
bool md_mcs(md_mcs_result_t* out, const md_mcs_desc_t* desc, struct md_allocator_i* alloc);
void md_mcs_result_free(md_mcs_result_t* result);

#ifdef __cplusplus
}
#endif
