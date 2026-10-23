#pragma once

// Shared structural invariants for md_system_t.
//
// Several fields of md_system_t are sentinel-terminated offset arrays: an array of
// 'count' offsets plus one trailing end-offset at index [count]. The range accessors
// in md_system.h read offset[i] and offset[i+1] unconditionally, so an array that is
// missing its sentinel, or whose offsets were recorded against a stale counter, is
// silently wrong (empty/inverted ranges) or reads out of bounds.
//
// This is easy to get wrong in a new loader and invisible in tests that only assert
// count > 0, so every code path that produces or transforms an md_system_t should run
// these checks.

#include "utest.h"

#include <md_system.h>
#include <core/md_array.h>

#define EXPECT_VALID_SYSTEM(sys) expect_valid_system(utest_result, __FILE__, __LINE__, (sys))

static void expect_valid_system(int* utest_result, const char* file, int line, const md_system_t* sys) {
#define FAIL(fmt, ...) do { \
        printf("%s:%i: failure: " fmt "\n", file, line, ##__VA_ARGS__); \
        *utest_result = UTEST_TEST_FAILURE; \
    } while (0)

    if (sys->residue.count) {
        if (!sys->residue.atom_offset) {
            FAIL("component.atom_offset is NULL with component.count = %zu", sys->residue.count);
        } else {
            if (md_array_size(sys->residue.atom_offset) != sys->residue.count + 1) {
                FAIL("component.atom_offset holds %zu entries, expected count+1 = %zu (missing sentinel?)",
                     md_array_size(sys->residue.atom_offset), sys->residue.count + 1);
            }
            if (sys->residue.atom_offset[0] != 0) {
                FAIL("component.atom_offset[0] = %u, expected 0", sys->residue.atom_offset[0]);
            }
            if (sys->residue.atom_offset[sys->residue.count] != sys->atom.count) {
                FAIL("component sentinel = %u, expected atom.count = %zu",
                     sys->residue.atom_offset[sys->residue.count], sys->atom.count);
            }
            for (size_t i = 0; i < sys->residue.count; ++i) {
                uint32_t beg = sys->residue.atom_offset[i];
                uint32_t end = sys->residue.atom_offset[i + 1];
                if (beg > end) {
                    FAIL("component[%zu] range is inverted: [%u, %u)", i, beg, end);
                    break;
                }
                if (beg == end) {
                    FAIL("component[%zu] range is empty: [%u, %u)", i, beg, end);
                    break;
                }
            }
        }
    }

    // The component types, once inferred: a type for every component, of its name and kind
    {
        const md_residue_type_data_t* type = &sys->residue.type;
        if (md_array_size(type->name) != type->count || md_array_size(type->kind) != type->count ||
            md_array_size(type->monomer) != type->count || md_array_size(type->color) != type->count) {
            FAIL("component.type arrays do not hold count = %zu entries", type->count);
        }
        if (sys->residue.type_idx) {
            if (md_array_size(sys->residue.type_idx) != sys->residue.count) {
                FAIL("component.type_idx holds %zu entries, expected component.count = %zu",
                     md_array_size(sys->residue.type_idx), sys->residue.count);
            } else {
                for (size_t i = 0; i < sys->residue.count; ++i) {
                    const md_residue_type_idx_t t = sys->residue.type_idx[i];
                    if (t < 0 || (size_t)t >= type->count) {
                        FAIL("component[%zu] has type %d, out of [0, %zu)", i, t, type->count);
                        break;
                    }
                    if (type->kind[t] != md_system_residue_kind(sys, i) ||
                        !str_eq(LBL_TO_STR(type->name[t]), md_system_residue_name(sys, i))) {
                        FAIL("component[%zu] is not of the name and kind of its type %d", i, t);
                        break;
                    }
                }
            }
        }
    }

    if (sys->instance.count) {
        if (!sys->instance.res_offset) {
            FAIL("instance.comp_offset is NULL with instance.count = %zu", sys->instance.count);
        } else {
            if (md_array_size(sys->instance.res_offset) != sys->instance.count + 1) {
                FAIL("instance.comp_offset holds %zu entries, expected count+1 = %zu (missing sentinel?)",
                     md_array_size(sys->instance.res_offset), sys->instance.count + 1);
            }
            if (sys->instance.res_offset[sys->instance.count] != sys->residue.count) {
                FAIL("instance sentinel = %u, expected component.count = %zu",
                     sys->instance.res_offset[sys->instance.count], sys->residue.count);
            }
            for (size_t i = 0; i < sys->instance.count; ++i) {
                if (sys->instance.res_offset[i] > sys->instance.res_offset[i + 1]) {
                    FAIL("instance[%zu] range is inverted: [%u, %u)", i,
                         sys->instance.res_offset[i], sys->instance.res_offset[i + 1]);
                    break;
                }
            }
        }
    }

    if (sys->protein_backbone.range.count) {
        if (md_array_size(sys->protein_backbone.range.offset) != sys->protein_backbone.range.count + 1) {
            FAIL("protein_backbone.range.offset holds %zu entries, expected count+1 = %zu",
                 md_array_size(sys->protein_backbone.range.offset), sys->protein_backbone.range.count + 1);
        } else if (sys->protein_backbone.range.offset[sys->protein_backbone.range.count] != sys->protein_backbone.segment.count) {
            FAIL("protein_backbone sentinel = %u, expected segment.count = %zu",
                 sys->protein_backbone.range.offset[sys->protein_backbone.range.count], sys->protein_backbone.segment.count);
        }
        if (md_array_size(sys->protein_backbone.range.inst_idx) != sys->protein_backbone.range.count) {
            FAIL("protein_backbone.range.inst_idx holds %zu entries, expected count = %zu",
                 md_array_size(sys->protein_backbone.range.inst_idx), sys->protein_backbone.range.count);
        }
    }

    if (sys->nucleic_backbone.range.count) {
        if (md_array_size(sys->nucleic_backbone.range.offset) != sys->nucleic_backbone.range.count + 1) {
            FAIL("nucleic_backbone.range.offset holds %zu entries, expected count+1 = %zu",
                 md_array_size(sys->nucleic_backbone.range.offset), sys->nucleic_backbone.range.count + 1);
        } else if (sys->nucleic_backbone.range.offset[sys->nucleic_backbone.range.count] != sys->nucleic_backbone.segment.count) {
            FAIL("nucleic_backbone sentinel = %u, expected segment.count = %zu",
                 sys->nucleic_backbone.range.offset[sys->nucleic_backbone.range.count], sys->nucleic_backbone.segment.count);
        }
        if (md_array_size(sys->nucleic_backbone.range.inst_idx) != sys->nucleic_backbone.range.count) {
            FAIL("nucleic_backbone.range.inst_idx holds %zu entries, expected count = %zu",
                 md_array_size(sys->nucleic_backbone.range.inst_idx), sys->nucleic_backbone.range.count);
        }
    }

    if (sys->molecule.count) {
        if (md_array_size(sys->molecule.offset) != sys->molecule.count + 1) {
            FAIL("structure.offset holds %zu entries, expected count+1 = %zu",
                 md_array_size(sys->molecule.offset), sys->molecule.count + 1);
        } else if (sys->molecule.offset[sys->molecule.count] != md_array_size(sys->molecule.atom_idx)) {
            FAIL("structure sentinel = %u, expected atom_idx length = %zu",
                 sys->molecule.offset[sys->molecule.count], md_array_size(sys->molecule.atom_idx));
        }
    }

#undef FAIL
}
