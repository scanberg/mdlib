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

    if (sys->component.count) {
        if (!sys->component.atom_offset) {
            FAIL("component.atom_offset is NULL with component.count = %zu", sys->component.count);
        } else {
            if (md_array_size(sys->component.atom_offset) != sys->component.count + 1) {
                FAIL("component.atom_offset holds %zu entries, expected count+1 = %zu (missing sentinel?)",
                     md_array_size(sys->component.atom_offset), sys->component.count + 1);
            }
            if (sys->component.atom_offset[0] != 0) {
                FAIL("component.atom_offset[0] = %u, expected 0", sys->component.atom_offset[0]);
            }
            if (sys->component.atom_offset[sys->component.count] != sys->atom.count) {
                FAIL("component sentinel = %u, expected atom.count = %zu",
                     sys->component.atom_offset[sys->component.count], sys->atom.count);
            }
            for (size_t i = 0; i < sys->component.count; ++i) {
                uint32_t beg = sys->component.atom_offset[i];
                uint32_t end = sys->component.atom_offset[i + 1];
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

    if (sys->instance.count) {
        if (!sys->instance.comp_offset) {
            FAIL("instance.comp_offset is NULL with instance.count = %zu", sys->instance.count);
        } else {
            if (md_array_size(sys->instance.comp_offset) != sys->instance.count + 1) {
                FAIL("instance.comp_offset holds %zu entries, expected count+1 = %zu (missing sentinel?)",
                     md_array_size(sys->instance.comp_offset), sys->instance.count + 1);
            }
            if (sys->instance.comp_offset[sys->instance.count] != sys->component.count) {
                FAIL("instance sentinel = %u, expected component.count = %zu",
                     sys->instance.comp_offset[sys->instance.count], sys->component.count);
            }
            for (size_t i = 0; i < sys->instance.count; ++i) {
                if (sys->instance.comp_offset[i] > sys->instance.comp_offset[i + 1]) {
                    FAIL("instance[%zu] range is inverted: [%u, %u)", i,
                         sys->instance.comp_offset[i], sys->instance.comp_offset[i + 1]);
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

    if (sys->structure.count) {
        if (md_array_size(sys->structure.offset) != sys->structure.count + 1) {
            FAIL("structure.offset holds %zu entries, expected count+1 = %zu",
                 md_array_size(sys->structure.offset), sys->structure.count + 1);
        } else if (sys->structure.offset[sys->structure.count] != md_array_size(sys->structure.atom_idx)) {
            FAIL("structure sentinel = %u, expected atom_idx length = %zu",
                 sys->structure.offset[sys->structure.count], md_array_size(sys->structure.atom_idx));
        }
    }

#undef FAIL
}
