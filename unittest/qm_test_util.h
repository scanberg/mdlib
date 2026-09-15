#pragma once

// Reading a quantum chemistry file the way every consumer now does.
//
// md_molden.h and md_trexio.h expose entry points and nothing else: no reader object, no accessors,
// no types describing what the file held. A file is parsed into an md_system_t, and everything it
// carried is read back out of that system - its atoms, its state, and its ATTRIBUTE TABLE. So these
// tests speak the same vocabulary a real consumer does, and a value that cannot be reached this way
// is a value no consumer can reach either, which is the property they exist to hold.
//
// The helpers are deliberately thin: each one is the two or three calls a consumer would otherwise
// make itself, named once so a test reads as the assertion it is making rather than as attribute
// plumbing. vlx_test_util.h is the same idea for the VeloxChem reader; this one is format neutral
// because the paths it reaches for are.

#include <string.h>
#include <math.h>

#include <md_system.h>
#include <md_gto.h>
#include <md_qm.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#define QM_TEST_ANGSTROM_TO_BOHR 1.8897261246257702

typedef struct qm_test_t {
    md_allocator_i*   alloc;   // Owns everything below; destroying it frees the system
    md_system_t       sys;
    md_system_state_t state;
} qm_test_t;

// Prepares a fresh system for a reader to fill in. The reader is then called by the test itself,
// because which one to call is the one thing these helpers must not decide.
static inline void qm_test_init(qm_test_t* out, size_t arena_bytes) {
    memset(out, 0, sizeof(*out));
    out->alloc = md_arena_allocator_create(md_get_heap_allocator(), arena_bytes);
    out->sys   = (md_system_t){ .alloc = out->alloc };
    out->state = (md_system_state_t){ .alloc = out->alloc };
}

static inline void qm_test_free(qm_test_t* t) {
    if (t->alloc) {
        md_arena_allocator_destroy(t->alloc);
    }
    memset(t, 0, sizeof(*t));
}

static inline const md_attribute_t* qm_test_attr(const qm_test_t* t, str_t path) {
    return md_attributes_find(&t->sys.attributes, path);
}

static inline bool qm_test_has(const qm_test_t* t, str_t path) {
    return md_attributes_find(&t->sys.attributes, path) != NULL;
}

// A rank 0 attribute. 'fallback' is what an absent path reads as, which is the caller's way of
// saying whether absence is a failure or simply a block the file did not carry.
static inline double qm_test_scalar(const qm_test_t* t, str_t path, double fallback) {
    const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
    double value = fallback;
    if (a && md_attribute_extract_f64(&value, 1, a, md_unit_none()) == 1) {
        return value;
    }
    return fallback;
}

// One entry of a string attribute; index 0 is the whole value for the rank 1 {1} case a single
// string is published as.
static inline str_t qm_test_string_at(const qm_test_t* t, str_t path, size_t index) {
    const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
    return a ? md_attribute_str(&t->sys.attributes, a, index) : (str_t){0};
}

static inline str_t qm_test_string(const qm_test_t* t, str_t path) {
    const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
    str_t value = {0};
    if (a && md_attribute_extract_str(&value, 1, &t->sys.attributes, a) == 1) {
        return value;
    }
    return (str_t){0};
}

// Length of a rank 1 series, 0 when the path is absent.
static inline size_t qm_test_count(const qm_test_t* t, str_t path) {
    const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
    return a ? md_attribute_element_count(&a->format) : 0;
}

static inline size_t qm_test_series(double* dst, size_t cap, const qm_test_t* t, str_t path) {
    const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
    return a ? md_attribute_extract_f64(dst, cap, a, md_unit_none()) : 0;
}

// One row of a rank 2 {rows, cols} attribute. The slice is the point: a computed attribute
// reconstructs only the plane that was asked for.
static inline size_t qm_test_row(double* dst, size_t cap, const qm_test_t* t, str_t path, size_t row) {
    const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
    if (!a) {
        return 0;
    }
    const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)row);
    return md_attribute_extract_slice_f64(dst, cap, a, &slice, md_unit_none());
}

// A whole rank 2 attribute, allocated from the test's own arena. NULL when the path is absent or
// the extraction comes up short, so a caller asserts on the pointer rather than on a count.
static inline double* qm_test_matrix(const qm_test_t* t, str_t path, size_t* out_dim) {
    const md_attribute_t* a = md_attributes_find(&t->sys.attributes, path);
    if (!a || a->format.rank != 2) {
        return NULL;
    }
    const size_t n = md_attribute_element_count(&a->format);
    double* dst = (double*)md_alloc(t->alloc, sizeof(double) * n);
    if (!dst || md_attribute_extract_f64(dst, n, a, md_unit_none()) != n) {
        return NULL;
    }
    if (out_dim) {
        *out_dim = a->format.shape[1];
    }
    return dst;
}

// The basis, rebuilt from basis/shell/* and basis/primitive/* - the format neutral paths, so this
// is the same call a consumer of any QM reader makes.
static inline bool qm_test_basis(md_gto_basis_t* out, const qm_test_t* t) {
    return md_gto_basis_extract_attributes(out, &t->sys.attributes, t->alloc);
}

// max |C S C^T - I| over the molecular orbitals of one spin channel.
//
// This is THE test of a QM reader, and it is one number. It fails if the atomic orbitals were put
// in the wrong order, if a shell's contraction was normalised against the wrong convention, if the
// spherical to Cartesian expansion lost a factor, or if the overlap does not describe the same
// basis the coefficients are stated against - and it cannot pass by accident, because an
// orthonormal set of orbitals in the wrong basis is not orthonormal in this one.
static inline double qm_test_orthonormality(const qm_test_t* t, str_t coefficient_path) {
    size_t num_ao = 0, dim = 0;
    const md_attribute_t* ca = md_attributes_find(&t->sys.attributes, coefficient_path);
    if (!ca || ca->format.rank != 2) {
        return INFINITY;
    }
    const size_t num_mo = ca->format.shape[0];
    double* C = qm_test_matrix(t, coefficient_path, &num_ao);
    double* S = qm_test_matrix(t, STR_LIT("basis/overlap"), &dim);
    if (!C || !S || dim != num_ao) {
        return INFINITY;
    }

    double* SC = (double*)md_alloc(t->alloc, sizeof(double) * num_mo * num_ao);
    if (!SC) {
        return INFINITY;
    }
    for (size_t m = 0; m < num_mo; ++m) {
        for (size_t a = 0; a < num_ao; ++a) {
            double v = 0.0;
            for (size_t b = 0; b < num_ao; ++b) {
                v += S[a * num_ao + b] * C[m * num_ao + b];
            }
            SC[m * num_ao + a] = v;
        }
    }

    double worst = 0.0;
    for (size_t i = 0; i < num_mo; ++i) {
        for (size_t j = 0; j < num_mo; ++j) {
            double v = 0.0;
            for (size_t a = 0; a < num_ao; ++a) {
                v += C[i * num_ao + a] * SC[j * num_ao + a];
            }
            const double err = fabs(v - ((i == j) ? 1.0 : 0.0));
            worst = (err > worst) ? err : worst;
        }
    }
    return worst;
}

// tr(D S), the electron count the density and the overlap agree on. The one cheap check that a
// density built from the published coefficients and occupations describes the right number of
// electrons - and it is sensitive to the occupation convention, which is what makes it worth
// asserting separately from orthonormality.
static inline double qm_test_electron_count(const qm_test_t* t, str_t density_path) {
    size_t dn = 0, sn = 0;
    double* D = qm_test_matrix(t, density_path, &dn);
    double* S = qm_test_matrix(t, STR_LIT("basis/overlap"), &sn);
    if (!D || !S || dn != sn) {
        return -1.0;
    }
    double tr = 0.0;
    for (size_t i = 0; i < dn; ++i) {
        for (size_t j = 0; j < dn; ++j) {
            tr += D[i * dn + j] * S[j * dn + i];
        }
    }
    return tr;
}
