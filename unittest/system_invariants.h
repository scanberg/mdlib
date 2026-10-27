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
            if (md_array_size(sys->instance.id) != sys->instance.count || md_array_size(sys->instance.auth_id) != sys->instance.count ||
                md_array_size(sys->instance.part_type_idx) != sys->instance.count || md_array_size(sys->instance.mol_type_idx) != sys->instance.count) {
                FAIL("instance arrays do not hold count = %zu entries", sys->instance.count);
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

    // THE HIERARCHY: every part of a type, the given types before the inferred, and every part within one molecule
    if (sys->part.count) {
        const md_part_type_data_t* type = &sys->part.type;
        if (md_array_size(sys->part.res_offset) != sys->part.count + 1 || md_array_size(sys->part.type_idx) != sys->part.count) {
            FAIL("part arrays do not hold count = %zu entries", sys->part.count);
        } else if (sys->part.res_offset[sys->part.count] != sys->residue.count) {
            FAIL("part sentinel = %u, expected residue.count = %zu", sys->part.res_offset[sys->part.count], sys->residue.count);
        } else {
            for (size_t p = 0; p < sys->part.count; ++p) {
                const md_part_type_idx_t t = sys->part.type_idx[p];
                if (t < 0 || (size_t)t >= type->count) {
                    FAIL("part[%zu] has type %d, out of [0, %zu)", p, t, type->count);
                    break;
                }
            }
        }
        if (md_array_size(type->id) != type->count || md_array_size(type->flags) != type->count || md_array_size(type->description) != type->count) {
            FAIL("part.type arrays do not hold count = %zu entries", type->count);
        } else {
            bool inferred = false;
            for (size_t t = 0; t < type->count; ++t) {
                const bool f = (type->flags[t] & MD_PART_TYPE_FLAG_INFERRED) != 0;
                if (inferred && !f) {
                    FAIL("part type %zu is given but comes after an inferred one", t);
                    break;
                }
                inferred |= f;
            }
        }
        if (sys->molecule.count && sys->molecule.atom_slot) {
            for (size_t p = 0; p < sys->part.count; ++p) {
                const md_urange_t r = md_system_part_atom_range(sys, p);
                if (r.beg >= r.end) continue;
                const md_molecule_idx_t m = md_system_molecule_find_by_atom_idx(sys, r.beg);
                bool whole = true;
                for (uint32_t a = r.beg + 1; a < r.end && whole; ++a) whole = md_system_molecule_find_by_atom_idx(sys, a) == m;
                if (!whole) {
                    FAIL("part[%zu] is not within one molecule", p);
                    break;
                }
            }
        }
    }

    if (sys->molecule.count && sys->molecule.type_idx) {
        const md_molecule_type_data_t* type = &sys->molecule.type;
        if (md_array_size(sys->molecule.type_idx) != sys->molecule.count) {
            FAIL("molecule.type_idx holds %zu entries, expected molecule.count = %zu", md_array_size(sys->molecule.type_idx), sys->molecule.count);
        } else {
            for (size_t m = 0; m < sys->molecule.count; ++m) {
                const md_molecule_type_idx_t t = sys->molecule.type_idx[m];
                if (t < 0 || (size_t)t >= type->count) {
                    FAIL("molecule[%zu] has type %d, out of [0, %zu)", m, t, type->count);
                    break;
                }
            }
        }
        bool inferred = false;
        for (size_t t = 0; t < type->count; ++t) {
            const bool f = (type->flags[t] & MD_MOLECULE_TYPE_FLAG_INFERRED) != 0;
            if (inferred && !f) {
                FAIL("molecule type %zu is given but comes after an inferred one", t);
                break;
            }
            inferred |= f;
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
        if (md_array_size(sys->protein_backbone.range.part_idx) != sys->protein_backbone.range.count) {
            FAIL("protein_backbone.range.part_idx holds %zu entries, expected count = %zu",
                 md_array_size(sys->protein_backbone.range.part_idx), sys->protein_backbone.range.count);
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
        if (md_array_size(sys->nucleic_backbone.range.part_idx) != sys->nucleic_backbone.range.count) {
            FAIL("nucleic_backbone.range.part_idx holds %zu entries, expected count = %zu",
                 md_array_size(sys->nucleic_backbone.range.part_idx), sys->nucleic_backbone.range.count);
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
