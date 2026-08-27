#pragma once

// Reader for the polarizable-embedding potential files (.pot) written by PyFraME
// (https://github.com/VeloxChem/PyFraME) and consumed by VeloxChem.
//
// The format is a flat sequence of sections, each opened by a '@'-prefixed keyword
// and closed by '@end'. Blank lines between sections are ignored.
//
//   @environment
//   units: angstrom
//   xyz:
//   O     -5.672000     2.390000    -4.911000  HOH_pe  1  OW
//   H     -6.099000     2.789000    -4.153000  HOH_pe  1  H1
//   ...
//   @end
//
//   @charges
//   O   -0.67444000  HOH_pe
//   H    0.33722000  HOH_pe
//   H    0.33722000  HOH_pe
//   @end
//
//   @polarizabilities
//   O    5.73935000   0.00000000   0.00000000   5.73935000   0.00000000   5.73935000  HOH_pe
//   ...
//   @end
//
// NOTE ON GRANULARITY:
// '@environment' lists every site in the system. '@charges' and '@polarizabilities'
// do NOT — they are per fragment *type*: one row per atom of a single representative
// fragment, tagged with that fragment's name. A site therefore picks up its parameters
// by (fragment name, index within its fragment). A fragment type may be absent from a
// section entirely (in the reference file 'HOH_npe' has charges but no polarizabilities,
// i.e. it is embedded non-polarizably).
//
// This reader deliberately performs no such resolution and no unit conversion — it
// reproduces the file contents. Mapping sites to parameters is the caller's job.

#include <core/md_str.h>

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;

typedef enum md_pot_unit_t {
    MD_POT_UNIT_UNKNOWN  = 0,   // No 'units:' directive was present, or it was not recognized
    MD_POT_UNIT_ANGSTROM = 1,   // 'angstrom' / 'aa' / 'a'
    MD_POT_UNIT_BOHR     = 2,   // 'au' / 'bohr' / 'atomic'
} md_pot_unit_t;

// Component order of md_pot_polarizability_t::alpha.
// The six values pack the upper triangle of the symmetric dipole-dipole polarizability
// row by row. This matches the Dalton/MOLCAS LoProp ordering that PyFraME reads its
// polarizabilities from.
//
// WARNING: every polarizability in the reference file is isotropic (all off-diagonals
// are exactly zero), so that file cannot discriminate between this order and the
// alternative column-wise packing (xx, xy, yy, xz, yz, zz). If you ever get an
// anisotropic .pot, verify against it before trusting an off-diagonal component.
enum {
    MD_POT_ALPHA_XX = 0,
    MD_POT_ALPHA_XY = 1,
    MD_POT_ALPHA_XZ = 2,
    MD_POT_ALPHA_YY = 3,
    MD_POT_ALPHA_YZ = 4,
    MD_POT_ALPHA_ZZ = 5,
};

// One row of the '@environment' xyz: block.
typedef struct md_pot_site_t {
    char   element[4];      // Element symbol as written, NUL-terminated
    double coord[3];        // In md_pot_t::unit, exactly as written
    str_t  fragment_name;   // e.g. "HOH_pe". Empty if the column was absent
    int    fragment_id;     // Fragment number as written. -1 if the column was absent
    str_t  atom_name;       // e.g. "OW". Empty if the column was absent
} md_pot_site_t;

// One row of the '@charges' block.
typedef struct md_pot_charge_t {
    char   element[4];
    double charge;
    str_t  fragment_name;   // Empty if the column was absent
} md_pot_charge_t;

// One row of the '@polarizabilities' block.
typedef struct md_pot_polarizability_t {
    char   element[4];
    double alpha[6];        // Indexed by MD_POT_ALPHA_*
    str_t  fragment_name;   // Empty if the column was absent
} md_pot_polarizability_t;

typedef struct md_pot_t {
    md_pot_unit_t unit;     // Unit of md_pot_site_t::coord, from the 'units:' directive

    size_t num_sites;
    md_pot_site_t* sites;

    size_t num_charges;
    md_pot_charge_t* charges;

    size_t num_polarizabilities;
    md_pot_polarizability_t* polarizabilities;

    // Interned strings backing every str_t above. Unique strings only, so identical
    // fragment names compare equal by pointer. Owned by this object; do not modify.
    size_t num_strings;
    str_t* strings;
} md_pot_t;

// Parse a .pot text blob. The returned object owns its memory and does not reference
// 'str', so 'str' may be freed immediately after.
bool md_pot_parse_str (md_pot_t* pot, str_t str, struct md_allocator_i* alloc);
bool md_pot_parse_file(md_pot_t* pot, str_t filename, struct md_allocator_i* alloc);

// Only free objects that were filled in by one of the parse functions.
void md_pot_free(md_pot_t* pot, struct md_allocator_i* alloc);

// Scale factor to convert md_pot_site_t::coord into Ångström.
// Returns 1.0 for MD_POT_UNIT_UNKNOWN (the format's default is angstrom).
double md_pot_unit_to_angstrom(md_pot_unit_t unit);

// Expand the packed upper triangle into a full symmetric 3x3 matrix.
static inline void md_pot_alpha_to_mat3(double dst[3][3], const double alpha[6]) {
    dst[0][0] = alpha[MD_POT_ALPHA_XX];
    dst[0][1] = dst[1][0] = alpha[MD_POT_ALPHA_XY];
    dst[0][2] = dst[2][0] = alpha[MD_POT_ALPHA_XZ];
    dst[1][1] = alpha[MD_POT_ALPHA_YY];
    dst[1][2] = dst[2][1] = alpha[MD_POT_ALPHA_YZ];
    dst[2][2] = alpha[MD_POT_ALPHA_ZZ];
}

#ifdef __cplusplus
}
#endif
