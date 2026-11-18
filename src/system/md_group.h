#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_types.h>
#include <core/md_str.h>

struct md_bitset_t;

// Named groups of atoms the data gives (GROMACS index groups) and the roles they had.

// GROUPS. Named sets of atoms the data gives, the way GROMACS has them: the index groups grompp was given (a tpr keeps
// those the simulation used) or an index file (.ndx). They are not of the hierarchy: a group holds any atoms, whole
// molecules or parts of them, and groups may overlap. Inference neither reads nor writes them.
//
// A ROLE is what the simulation used groups for, an mdp option naming them: tc-grps, energygrps and so on. The groups of
// a role are a partition of the atoms, in the order the option names them, the atoms named by none of them in a group of
// their own named 'rest'. The order is kept: a script takes a role's groups in it, and it tells apart the pairs of an
// energy file whose names hold '-' (md_edr_term_groups). A group is stored once whatever roles it has; two groups may
// share a name ('rest' of two roles is two sets of atoms, as two index files may give one name to two sets).
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

// A consumer reads a group's name, atoms and roles: whether a group is one of a role's is a bit of its roles, and a
// role no group has is one the data does not give. The order of a role's groups is read through
// md_system_group_role_order; md_system_group_set_role keeps both, so they always agree.
typedef struct md_group_data_t {
    size_t count;
    md_str_t*    name;                  // [count]
    struct md_bitset_t* atoms;          // [count] Each its own words, in the system's allocator
    uint32_t*    roles;                 // [count] Bit r: the group is one of role r's (md_group_role_t)

    // The groups of each role in its order, which md_system_group_role_order gives: the groups of role r are
    // group[offset[r] .. offset[r + 1]]. A role the data does not give, or gives as one group of every atom it never
    // named ('rest' of all), has none.
    struct {
        uint32_t  offset[MD_GROUP_ROLE_COUNT + 1];
        uint32_t* group;
    } role_order;
} md_group_data_t;

#ifdef __cplusplus
extern "C" {
#endif

// The groups of a role, as indices of groups in the role's order (the order of its mdp option): count receives how many,
// 0 (and NULL returned) for a role without groups
static inline const uint32_t* md_system_group_role_order(const md_group_data_t* g, md_group_role_t role, size_t* count) {
    *count = 0;
    if (!g || (uint32_t)role >= MD_GROUP_ROLE_COUNT) return NULL;
    const uint32_t beg = g->role_order.offset[role];
    const uint32_t end = g->role_order.offset[role + 1];
    *count = end - beg;
    return end > beg ? g->role_order.group + beg : NULL;
}

// The name of a role as GROMACS' mdp option has it ("energygrps"), and as a script names it ("energy"): NULL for none
const char* md_group_role_mdp_option(md_group_role_t role);
const char* md_group_role_name(md_group_role_t role);

#ifdef __cplusplus
}
#endif
