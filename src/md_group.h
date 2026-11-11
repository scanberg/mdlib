#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_types.h>
#include <core/md_bitset.h>

// Named groups of atoms the data gives (GROMACS index groups) and the roles they had.

// GROUPS. Named sets of atoms the data gives, the way GROMACS has them: the index groups grompp was given (a tpr keeps
// those the simulation used) or an index file (.ndx). They are not of the hierarchy: a group holds any atoms, whole
// molecules or parts of them, and groups may overlap. Inference neither reads nor writes them.
//
// A ROLE is what the simulation used groups for, an mdp option naming them: tc-grps, energygrps and so on. The groups of
// a role are a partition of the atoms, in the order the option names them, the atoms named by none of them in a group of
// their own named 'rest'. The order is load bearing: an energy file names the pairs of energy groups by it
// (md_edr_term_groups). A group is stored once whatever roles it has; two groups may share a name ('rest' of two roles
// is two sets of atoms, as two index files may give one name to two sets).
typedef enum md_group_role_t {
    MD_GROUP_ROLE_TC = 0,               // tc-grps: temperature coupling
    MD_GROUP_ROLE_ENERGY,               // energygrps: energy output, the pairs of the non-bonded terms
    MD_GROUP_ROLE_ACC,                  // acc-grps: acceleration
    MD_GROUP_ROLE_FREEZE,               // freezegrps
    MD_GROUP_ROLE_USER1,                // user1-grps
    MD_GROUP_ROLE_USER2,                // user2-grps
    MD_GROUP_ROLE_COMM,                 // comm-grps: removal of the center of mass motion
    MD_GROUP_ROLE_COMPRESSED_X,         // compressed-x-grps: the atoms written to the .xtc
    MD_GROUP_ROLE_ORIRE_FIT,            // orire-fitgrp: fit of the orientation restraints
    MD_GROUP_ROLE_QMMM,                 // qmmm-grps: the QM atoms
    MD_GROUP_ROLE_COUNT,
} md_group_role_t;

typedef struct md_group_data_t {
    size_t count;
    md_str_t*    name;                  // [count]
    md_bitset_t* atoms;                 // [count] Each its own words, in the system's allocator
    uint32_t*    roles;                 // [count] Bit r: the group is one of role r's (md_group_role_t)
    // The groups of role r, in its order: role_group[role_offset[r] .. role_offset[r + 1]]. A role the data does not
    // give, or gives as one group of every atom it never named ('rest' of all), has none.
    uint32_t     role_offset[MD_GROUP_ROLE_COUNT + 1];
    uint32_t*    role_group;
} md_group_data_t;

#ifdef __cplusplus
extern "C" {
#endif

static inline md_urange_t md_system_group_role_range(const md_group_data_t* g, md_group_role_t role) {
    md_urange_t r = {0, 0};
    if (g && (uint32_t)role < MD_GROUP_ROLE_COUNT) {
        r.beg = g->role_offset[role];
        r.end = g->role_offset[role + 1];
    }
    return r;
}

// The name of a role as GROMACS' mdp option has it ("energygrps"), and as a script names it ("energy"): NULL for none
const char* md_group_role_mdp_option(md_group_role_t role);
const char* md_group_role_name(md_group_role_t role);

#ifdef __cplusplus
}
#endif
