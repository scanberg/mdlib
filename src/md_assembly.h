#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_types.h>
#include <core/md_label.h>
#include <core/md_vec.h>

// The biological assemblies of a deposited structure.

// ASSEMBLIES. The biological assemblies of a deposited structure (mmCIF pdbx_struct_assembly, PDB REMARK 350): an
// assembly is made of groups, and a group is a set of instances of the system copied by each of a set of operators,
// rigid transforms of the coordinates of the system. The atoms are not copied: an assembly is a way of showing those of
// the system, with the copies of its instances where the operators take them. The coordinates of the system are an
// assembly only where one says so: one which keeps them as they are has the identity among its operators, and of the
// instances only those it names. The assemblies name instances (the asyms and chains of the file), which inference
// leaves as they are, rather than parts, which it may make anew (md_instance_data_t).
//
// The form is canonical: within an assembly every operator is in one group only (the rows of a file which share an
// operator are merged, md_system_assemblies_init), so that the operator names one copy, and what is bonded across
// instances of a copy (a glycan, a ligand, a disulfide between chains) is in the same group. The operators are the
// system's, shared by its assemblies; a product of operators (mmCIF "(1-60)(61-88)") is an operator of its own.
typedef struct md_assembly_data_t {
    size_t              count;              // Assemblies
    md_label_t*         id;                 // [count] As the file names them ("1", "2", "PAU")
    uint32_t*           group_offset;       // [count + 1] The groups of assembly i: [group_offset[i], group_offset[i + 1])

    size_t              group_count;
    uint32_t*           group_inst_offset;  // [group_count + 1] The instances of each group in inst
    md_instance_idx_t*  inst;               // In increasing order within a group
    uint32_t*           group_oper_offset;  // [group_count + 1] The operators of each group in oper
    uint32_t*           oper;               // Operator indices, in the order the file gives them

    size_t              oper_count;
    md_label_t*         oper_id;            // [oper_count] As the file names it ("1", "X0"); a product its factors ("1x61")
    md_float4x4_t*             oper_transform;     // [oper_count] Rigid: a rotation (possibly improper) and a translation
} md_assembly_data_t;

// One row of the description of an assembly, as the files give them: the instances which each of the operators copies
// (mmCIF pdbx_struct_assembly_gen, a group of chains of REMARK 350 with its BIOMT)
typedef struct md_assembly_gen_t {
    md_label_t                  assembly_id;
    const md_instance_idx_t*    inst;
    size_t                      inst_count;
    const uint32_t*             oper;       // Indices of operators
    size_t                      oper_count;
} md_assembly_gen_t;

#ifdef __cplusplus
extern "C" {
#endif

static inline size_t md_assembly_count(const md_assembly_data_t* a) {
    MD_ASSERT(a);
    return a->count;
}

// The groups of an assembly
static inline md_urange_t md_assembly_group_range(const md_assembly_data_t* a, size_t asm_idx) {
    MD_ASSERT(a);
    md_urange_t r = {0};
    if (a->group_offset && asm_idx < a->count) {
        r.beg = a->group_offset[asm_idx];
        r.end = a->group_offset[asm_idx + 1];
    }
    return r;
}

// The instances of a group: inst[beg, end)
static inline md_urange_t md_assembly_group_inst_range(const md_assembly_data_t* a, size_t group_idx) {
    MD_ASSERT(a);
    md_urange_t r = {0};
    if (a->group_inst_offset && group_idx < a->group_count) {
        r.beg = a->group_inst_offset[group_idx];
        r.end = a->group_inst_offset[group_idx + 1];
    }
    return r;
}

// The operators of a group: oper[beg, end)
static inline md_urange_t md_assembly_group_oper_range(const md_assembly_data_t* a, size_t group_idx) {
    MD_ASSERT(a);
    md_urange_t r = {0};
    if (a->group_oper_offset && group_idx < a->group_count) {
        r.beg = a->group_oper_offset[group_idx];
        r.end = a->group_oper_offset[group_idx + 1];
    }
    return r;
}

static inline md_str_t md_assembly_id(const md_assembly_data_t* a, size_t asm_idx) {
    MD_ASSERT(a);
    md_str_t id = MD_STR_INIT("");
    if (a->id && asm_idx < a->count) {
        id = MD_LABEL_STR(a->id[asm_idx]);
    }
    return id;
}

// The copies an assembly is made of: the operators of all of its groups
static inline size_t md_assembly_copy_count(const md_assembly_data_t* a, size_t asm_idx) {
    MD_ASSERT(a);
    size_t n = 0;
    const md_urange_t groups = md_assembly_group_range(a, asm_idx);
    for (uint32_t p = groups.beg; p < groups.end; ++p) {
        const md_urange_t ops = md_assembly_group_oper_range(a, p);
        n += ops.end - ops.beg;
    }
    return n;
}

#ifdef __cplusplus
}
#endif
