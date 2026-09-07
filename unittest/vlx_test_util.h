#pragma once

// Reading a VeloxChem file the way every consumer now does.
//
// md_vlx.h exposes three entry points and nothing else: no reader object, no accessors, no types
// describing what the file held. A file is parsed into an md_system_t, and everything it carried is
// read back out of that system - its atoms, its state, and its ATTRIBUTE TABLE. So the vlx tests
// speak the same vocabulary a real consumer does, and these helpers are deliberately thin: each one
// is the two or three calls a consumer would otherwise make itself, named once so a test reads as
// the assertion it is making rather than as attribute plumbing.

#include <string.h>

#include <md_vlx.h>
#include <md_system.h>
#include <md_gto.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#define VLX_TEST_ANGSTROM_TO_BOHR 1.8897261246257702

typedef struct vlx_test_t {
	md_allocator_i*   alloc;   // Owns everything below; destroying it frees the system
	md_system_t       sys;
	md_system_state_t state;
} vlx_test_t;

// Parses 'filename' into a fresh system. Returns false exactly when the reader declines the file.
static inline bool vlx_test_load(vlx_test_t* out, str_t filename, size_t arena_bytes) {
	memset(out, 0, sizeof(*out));
	out->alloc = md_arena_allocator_create(md_get_heap_allocator(), arena_bytes);
	out->sys   = (md_system_t){ .alloc = out->alloc };
	out->state = (md_system_state_t){ .alloc = out->alloc };
	return md_vlx_system_init_from_file(&out->sys, &out->state, filename);
}

static inline void vlx_test_free(vlx_test_t* t) {
	if (t->alloc) {
		md_arena_allocator_destroy(t->alloc);
	}
	memset(t, 0, sizeof(*t));
}

static inline const md_attribute_t* vlx_test_attr(const vlx_test_t* t, str_t path) {
	return md_attributes_find(&t->sys.attributes, path);
}

static inline bool vlx_test_has(const vlx_test_t* t, str_t path) {
	return md_attributes_find(&t->sys.attributes, path) != NULL;
}

// A rank 0 attribute. 'fallback' is what an absent path reads as, which is the caller's way of
// saying whether absence is a failure or simply a block the file did not carry.
static inline double vlx_test_scalar(const vlx_test_t* t, str_t path, double fallback) {
	const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
	double value = fallback;
	if (a && md_attribute_extract_f64(&value, 1, a, md_unit_none()) == 1) {
		return value;
	}
	return fallback;
}

static inline str_t vlx_test_string(const vlx_test_t* t, str_t path) {
	const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
	str_t value = {0};
	if (a && md_attribute_extract_str(&value, 1, &t->sys.attributes, a) == 1) {
		return value;
	}
	return (str_t){0};
}

// Length of a rank 1 series, 0 when the path is absent.
static inline size_t vlx_test_count(const vlx_test_t* t, str_t path) {
	const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
	return a ? md_attribute_element_count(&a->format) : 0;
}

static inline size_t vlx_test_series(double* dst, size_t cap, const vlx_test_t* t, str_t path) {
	const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
	return a ? md_attribute_extract_f64(dst, cap, a, md_unit_none()) : 0;
}

// One row of a rank 2 {rows, cols} attribute, or one plane of a rank 3 {n, r, c}. The slice is the
// point: a computed attribute reconstructs only the plane that was asked for.
static inline size_t vlx_test_row(double* dst, size_t cap, const vlx_test_t* t, str_t path, size_t row) {
	const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
	if (!a) {
		return 0;
	}
	const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)row);
	return md_attribute_extract_slice_f64(dst, cap, a, &slice, md_unit_none());
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

// The basis, rebuilt from basis/shell/* and basis/primitive/* - the format neutral paths, so this
// is the same call a consumer of any other QM reader makes.
static inline bool vlx_test_basis(md_gto_basis_t* out, const vlx_test_t* t) {
	return md_gto_basis_extract_attributes(out, &t->sys.attributes, t->alloc);
}
