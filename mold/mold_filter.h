#ifndef _MOLD_FILTER_H_
#define _MOLD_FILTER_H_

#include <stdint.h>
#include <stdbool.h>
#include "mold_molecule.h"

#ifdef __cplusplus
extern "C" {
#endif

//typedef struct mold_filter mold_filter;
typedef struct mold_filter_context mold_filter_context;
typedef struct mold_filter_stored_selection mold_filter_stored_selection;
typedef struct mold_filter_molecule mold_filter_molecule;

// mold filter molecule represents a const subset of the fields within mold_molecule that are used internally in the filter
struct mold_filter_molecule {
    struct {
        uint32_t count;
        const float*            x;
        const float*            y;
        const float*            z;
        const float*            radius;
        const float*            mass;
        const uint8_t*          element;
        const char**            name;
        const float*            bfactor;
        const float*            occupancy;
        const mold_bond_entry*  bond;
    } atom;

    struct {
        uint32_t count;
        const char**        name;
        const uint32_t*     id;
        const mold_range*   atom_range;
        const uint32_t*     secondary_structure;
    } residue;

    struct {
        uint32_t count;
        const char**      id;
        const mold_range* atom_range;
    } chain;
};

struct mold_filter_stored_selection {
    const char*     ident;
    const uint64_t* bits;
};

struct mold_filter_context {
    const mold_filter_molecule*         mol;
    uint64_t                            bit_count;
    uint64_t*                           bits;
    uint32_t                            sel_count;
    const mold_filter_stored_selection* sel;
};

//typedef struct mold_filter mold_filter;

void mold_filter_molecule_extract(mold_filter_molecule* filt_mol, const mold_molecule* mol);
bool mold_filter_parse(const char* expression, const mold_filter_context* context);
//mold_filter* mold_filter_create(const char* expr_str, uint32_t expr_len);
//void         mold_filter_destroy(mold_filter* filter);
//bool         mold_filter_apply(uint64_t* bits, uint64_t bit_count, const mold_filter* filter, const mold_filter_context* context);

#ifdef __cplusplus
}
#endif

#endif