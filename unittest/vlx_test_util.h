#pragma once

// Reading a VeloxChem file the way every consumer now does.
//
// md_vlx.h exposes three entry points and nothing else: no reader object, no accessors, no types
// describing what the file held. A file is parsed into an md_system_t, and everything it carried is
// read back out of that system - its atoms, its state, and its ATTRIBUTE TABLE.
//
// Everything that reads such a system is in qm_test_util.h, because none of it is VeloxChem
// specific - the paths it reaches for are the format neutral ones. What is left here is the two
// things that ARE: calling this reader, and the frontier orbital rule its tests use.

#include "qm_test_util.h"

#include <md_vlx.h>

#define VLX_TEST_ANGSTROM_TO_BOHR QM_TEST_ANGSTROM_TO_BOHR

typedef qm_test_t vlx_test_t;

// Parses 'filename' into a fresh system. Returns false exactly when the reader declines the file.
static inline bool vlx_test_load(vlx_test_t* out, str_t filename, size_t arena_bytes) {
    qm_test_init(out, arena_bytes);
    return md_vlx_system_init_from_file(&out->sys, &out->state, filename);
}

// The frontier orbital, derived from the published occupations by the SAME rule the reader used
// before the occupations were the only thing published: the first orbital with zero occupancy is
// the LUMO. Returns the molecular orbital count when every orbital is occupied.
static inline size_t vlx_test_lumo_idx(const vlx_test_t* t, str_t occupation_path) {
    const md_attribute_t* a = md_attributes_find(&t->sys.attributes, occupation_path);
    if (!a) {
        return 0;
    }
    const size_t num_mo = md_attribute_element_count(&a->format);
    double* occ = (double*)md_alloc(t->alloc, sizeof(double) * num_mo);
    size_t lumo = num_mo;
    if (md_attribute_extract_f64(occ, num_mo, a, md_unit_none()) == num_mo) {
        for (size_t i = 0; i < num_mo; ++i) {
            if (occ[i] == 0.0) {
                lumo = i;
                break;
            }
        }
    }
    md_free(t->alloc, occ, sizeof(double) * num_mo);
    return lumo;
}

// The QM atom geometry, converted to the bohr the GTO evaluator works in. Returns the atom count.
static inline size_t vlx_test_atom_xyz_bohr(float* dst, size_t cap, const vlx_test_t* t) {
    const md_attribute_t* a = md_attributes_find(&t->sys.attributes, STR_LIT("qm/atom/coordinate"));
    if (!a) {
        return 0;
    }
    const size_t num_values = md_attribute_element_count(&a->format);   // value_count * 3
    if (num_values > cap) {
        return 0;
    }
    double* xyz = (double*)md_alloc(t->alloc, sizeof(double) * num_values);
    size_t num_atoms = 0;
    if (md_attribute_extract_f64(xyz, num_values, a, md_unit_none()) == num_values) {
        for (size_t i = 0; i < num_values; ++i) {
            dst[i] = (float)(xyz[i] * VLX_TEST_ANGSTROM_TO_BOHR);
        }
        num_atoms = num_values / 3;
    }
    md_free(t->alloc, xyz, sizeof(double) * num_values);
    return num_atoms;
}
