#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_bitset.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_system_t;

// GROMACS index files (.ndx): named groups of atoms, what grompp and the gmx tools are given with -n. A topology never
// refers to one; an mdp names its groups (energygrps = Protein SOL) and grompp finds them in the index file it is given.
//
//     [ Protein ]
//        1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
//       16   17 ...
//     [ SOL ]
//      ...
//
// A group is a header '[ name ]' and the atom numbers after it (1 based, as many to a line as there are), up to the
// next header. A group may hold an atom more than once and in any order; the groups may overlap. Text after ';' is a
// comment. An index file names atoms by number only: it is of whatever system has that many atoms or more.

typedef struct md_ndx_group_t {
    md_str_t  name;
    size_t    count;
    uint32_t* atom_idx;     // [count] Zero based, in the order of the file
} md_ndx_group_t;

typedef struct md_ndx_data_t {
    size_t          num_groups;
    md_ndx_group_t* groups;     // In the order of the file
    uint32_t        max_atom;   // One past the largest atom index of any group: the atoms a system must have
} md_ndx_data_t;

bool md_ndx_data_parse_str (md_ndx_data_t* data, md_str_t str, struct md_allocator_i* alloc);
bool md_ndx_data_parse_file(md_ndx_data_t* data, md_str_t filename, struct md_allocator_i* alloc);
void md_ndx_data_free(md_ndx_data_t* data, struct md_allocator_i* alloc);

// Adds the groups of the index file to the system (md_system_group_add), without roles: an index file says what groups
// are, not what a simulation used them for. The energy file of the run tells that (md_edr_system_infer_group_roles). Two groups of the same name and atoms are one; a group already in the
// system (from the tpr of the same run) gains nothing. False, and nothing added, when the file names an atom the system
// does not have: it is of another system.
bool md_ndx_system_supplement(struct md_system_t* sys, const md_ndx_data_t* data);
bool md_ndx_system_supplement_from_file(struct md_system_t* sys, md_str_t filename);

// GROMACS' DEFAULT GROUPS. Without an index file, grompp (and gmx make_ndx, which starts from them) makes groups of the
// topology by its residue types (share/top/residuetypes.dat) and atom names, and an mdp names those: tc-grps = Protein
// non-Protein, comm-grps = System. In GROMACS' order:
//
//     System                     every atom
//     Protein, Protein-H, C-alpha, Backbone, MainChain, MainChain+Cb, MainChain+H, SideChain, SideChain-H, Prot-Masses,
//     non-Protein                where there are residues of protein
//     Water, SOL, non-Water      where there is water
//     Ion                        where there are ions
//     Other (or DNA, RNA)        the first other residue type there is, and then a group for each residue name that is
//                                not of protein, DNA, RNA or water, in the order they first appear (ions too: CA, NA)
//     Water_and_ions             where there are both
//
// The atoms of the default group of a name, matched without regard to case, as grompp matches the names of an mdp: false,
// and out untouched, when there is no such group. The words of out are allocated with alloc.
bool md_ndx_default_group(md_bitset_t* out, const struct md_system_t* sys, md_str_t name, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
