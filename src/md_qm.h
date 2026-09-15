#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <md_system.h>
#include <core/md_str.h>
#include <core/md_unit.h>
#include <core/md_vec_math.h>

// Shared conventions for QUANTUM CHEMISTRY READERS.
//
// A QM reader parses a file into an md_system_t and everything it carried is read back out of that
// system's ATTRIBUTE TABLE - see the ATTRIBUTES section of md_system.h, and md_vlx.h for the
// worked example. Several of those paths are FORMAT NEUTRAL: basis/, orbital/ and qm/ mean the
// same thing whichever program wrote the file, and md_gto_basis_extract_attributes reads the basis
// back without knowing which reader filled it in. This is the WRITE side of that contract, so a
// second reader does not re-derive a layout a consumer already depends on.
//
// It also holds the NORMALISATION a reader has to apply. md_gto.h's AO CONVENTION block fixes what
// a coefficient means to this library; a file's coefficients mean something slightly different in
// every format, and getting from one to the other is arithmetic that must exist exactly once.
//
// md_vlx.c, md_molden.c and md_trexio.c all publish through this. Nothing here is VeloxChem,
// Molden or TREXIO specific: if a third of a thing lives here and two thirds in a reader, the
// reader is doing something its format actually requires.
//
// The PUBLISHING VOCABULARY below is attribute plumbing rather than quantum chemistry, and it sits
// here because the three QM readers are its only callers. If a non-QM producer ever wants the same
// thing it graduates to md_system.h; until then, keeping it here is what stops a fourth reader
// inventing a thirteenth spelling of "publish a series of doubles".

struct md_gto_basis_t;

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// UNITS
// ---------------------------------------------------------------------------
// The ones md_unit.h does not already name and more than one reader needs. Everything else -
// hartree, bohr, e a0, the Bohr magneton - is already there and is called directly.

md_unit_t md_qm_unit_wavenumber(void);   // cm^-1, what a vibrational frequency is stated in
md_unit_t md_qm_unit_km_per_mol(void);   // what an IR intensity is stated in
md_unit_t md_qm_unit_amu(void);          // u (Da) as a scaled kilogram, for a reduced mass

// ---------------------------------------------------------------------------
// PUBLISHING VOCABULARY
// ---------------------------------------------------------------------------
// One spelling each for the shapes a QM reader actually publishes. Every one of them REPLACES what
// is already at that path - see md_attributes_replace on why that is what a producer wants - and
// every one returns the id, so an alias or a follow up write has something to name.
//
// A NULL or empty input publishes NOTHING and returns MD_ATTRIBUTE_INVALID rather than an empty
// attribute: an absent path is how a consumer learns a block is missing, and a zero length series
// would read as a block that was there and held nothing.

// The general form. 'format' says what the layout is; see the ATTRIBUTES section of md_system.h.
md_attribute_id_t md_qm_publish(struct md_system_t* sys, str_t path, str_t label, md_unit_t unit,
                                md_attribute_format_t format, const void* data, size_t byte_size);

// The same, for an attribute computed through a provider instead of one copied in.
md_attribute_id_t md_qm_publish_virtual(struct md_system_t* sys, str_t path, str_t label, md_unit_t unit,
                                        md_attribute_format_t format, const md_attribute_virtual_t* virt);

md_attribute_id_t md_qm_publish_scalar     (struct md_system_t* sys, str_t path, str_t label, md_unit_t unit, double value);                                  // rank 0
md_attribute_id_t md_qm_publish_series     (struct md_system_t* sys, str_t path, str_t label, md_unit_t unit, const double* values, size_t count);            // rank 1 {N}
md_attribute_id_t md_qm_publish_vec3_series(struct md_system_t* sys, str_t path, str_t label, md_unit_t unit, const dvec3_t* values, size_t count);           // rank 1 {N} x 3
md_attribute_id_t md_qm_publish_matrix     (struct md_system_t* sys, str_t path, str_t label, md_unit_t unit, const double* values, size_t rows, size_t cols); // rank 2 {R,C}
md_attribute_id_t md_qm_publish_str        (struct md_system_t* sys, str_t path, str_t label, str_t value);                                                   // rank 1 {1} string
md_attribute_id_t md_qm_publish_strings    (struct md_system_t* sys, str_t path, str_t label, const str_t* values, size_t count);                             // rank 1 {N} strings

// One COLUMN of an array of structs, gathered into the table's own storage. A record of six fields
// is six sibling paths over one index space, because a value has ONE type - so the transposition
// from a file's row layout happens once, here, rather than in every consumer.
md_attribute_id_t md_qm_publish_column(struct md_system_t* sys, str_t path, str_t label, md_unit_t unit,
                                       md_attribute_type_t type, const void* base, size_t stride, size_t count);

// The anchor of a vector group: rank 0, one 3 component value, constant over whatever index space
// the group's vector has. Angstrom, because it is a point in system space. NOT replicated to match
// the vector's shape - group members share an index space, not a shape.
md_attribute_id_t md_qm_publish_origin(struct md_system_t* sys, str_t path, dvec3_t origin);

// A second name for an attribute already published. Both names then read one datum - no copy, and a
// consumer of either is unaffected when the other appears or goes.
md_attribute_id_t md_qm_alias(struct md_system_t* sys, md_attribute_id_t target, str_t path);

// Publishes beta's copy of a per orbital series, or a SECOND NAME for alpha's when the two share
// storage. Comparing the POINTERS is what tells the cases apart, and it is the only test that gets
// restricted open shell right, where the orbitals are shared and the occupations are not.
md_attribute_id_t md_qm_publish_or_alias(struct md_system_t* sys, md_attribute_id_t alpha_id, str_t path, str_t label,
                                         md_unit_t unit, const double* alpha_values, const double* beta_values, size_t count);

// Builds an attribute path from a fixed group prefix and a name taken from a file. A '/' inside the
// name would silently introduce a group level in a namespace where the separator is the only
// structure there is, so it is folded to '_'. Returns an empty str_t when the name does not fit,
// which a caller treats as "skip this one" rather than as a reason to stop publishing.
str_t md_qm_attribute_path(char* buf, size_t cap, str_t group, str_t name);

// ---------------------------------------------------------------------------
// NORMALISATION
// ---------------------------------------------------------------------------
// md_gto evaluates AO number n of a shell as
//
//     phi_n(r) = f(i,j,k) * x^i y^j z^k * sum_p coeff[p] * exp(-alpha[p] * r^2)
//
// so the shell's coeff[] carries the whole radial normalisation and f(i,j,k) the per monomial part.
// A file, by contrast, states a contraction over PRIMITIVES THAT ARE ALREADY NORMALISED and leaves
// the rest implicit. The two functions below are that gap, and applying them in order - scale each
// primitive, then scale the shell - lands a reader exactly on md_gto's convention.

// The coefficient of a single normalised primitive of exponent 'alpha' in a shell of angular
// momentum l, expressed against md_gto's Cartesian AO. Multiply a file's contraction coefficient by
// this. The value depends on l because the radial normalisation does.
double md_qm_primitive_norm_factor(uint32_t l, double alpha);

// The factor a whole contraction must be DIVIDED by so that the AO md_gto evaluates comes out with
// unit norm, with the primitive factor above already folded into coeff[].
//
// 'spherical' picks WHICH AO is normalised, and the distinction is not cosmetic: md_gto's spherical
// functions are the real solid harmonics as its expansion tables spell them, which are mutually
// orthogonal but NOT unit normalised - the d functions come out sqrt(12) too large, the f functions
// sqrt(60). A reader of a file that stores spherical coefficients therefore has to normalise the
// SPHERICAL function, or every d coefficient it publishes is off by a constant that looks entirely
// plausible on screen. Pass false for a file that stores Cartesian coefficients, where md_gto's own
// f(i,j,k) already equalises the monomials within a shell.
//
// Returns 0 for a degenerate shell (no primitives, or a contraction that sums to nothing), which a
// caller treats as "skip this shell" rather than dividing by it.
double md_qm_shell_norm_factor(uint32_t l, const double alpha[], const double coeff[], size_t num_primitives, bool spherical);

// The factor a CARTESIAN AO coefficient read from a file must be multiplied by, for Cartesian AO
// 'cart_idx' of a shell of angular momentum l.
//
// Gaussian, Molden's [6D] and friends normalise a Cartesian shell once, against the (l,0,0)
// monomial, and use that one constant for all (l+1)(l+2)/2 functions - so their xy function is not
// unit normalised and their xx one is. md_gto normalises every monomial in the shell. This is the
// ratio between the two, sqrt((2i-1)!!(2j-1)!!(2k-1)!! / (2l-1)!!), and it is 1 for every AO of an
// s or p shell, which is why a reader that forgets it still looks correct until the first d.
double md_qm_cart_coeff_factor(uint32_t l, uint32_t cart_idx);

// Converts a whole [num_mo][n_sph] coefficient matrix into [num_mo][n_cart], row by row, where
// n_sph is md_gto_basis_num_sph_ao(basis) and n_cart is md_gto_basis_num_ao(basis). 'src' must
// already be in the basis's own shell order and in md_gto's m-ascending order within each shell;
// getting a file's order onto that is the reader's job and the one part of this that cannot be
// shared. Returns num_mo on success, 0 on failure. In place is not supported: the row stride grows.
size_t md_qm_sph_to_cart_coefficients(double* dst, const double* src, size_t num_mo, const struct md_gto_basis_t* basis);

// ---------------------------------------------------------------------------
// PUBLISHING
// ---------------------------------------------------------------------------

// basis/shell/{atom_index,primitive_offset,primitive_count,angular_momentum} and
// basis/primitive/{exponent,coefficient} - the format neutral layout
// md_gto_basis_extract_attributes reads back. Published as COLUMNS and not as md_gto_shell_t
// records, so nothing stored is a struct layout contract.
//
// basis/shell/atom_index indexes the QM ATOM DOMAIN (qm/atom/*), not the system's atoms.
bool md_qm_publish_basis(struct md_system_t* sys, const struct md_gto_basis_t* basis);

// qm/atom/{atomic_number,coordinate} - the atoms the calculation covered, in ITS order and at ITS
// geometry, which is not necessarily the system's atom set. Coordinates are taken in ANGSTROM, to
// match the system's own, and are published as such.
bool md_qm_publish_atoms(struct md_system_t* sys, const uint8_t atomic_number[], const dvec3_t coord_angstrom[], size_t count);

// basis/overlap - the AO overlap S[a][b] = <phi_a|phi_b> over the Cartesian AOs md_gto evaluates, as
// a VIRTUAL attribute computed from basis/shell/*, basis/primitive/* and qm/atom/coordinate, which
// must already be published.
//
// COMPUTED AND NOT READ, even from a file that carries its own. A file that stored spherical AO
// data states its overlap in that basis, and the Cartesian embedding of a spherical basis is not
// invertible - S_sph = T S_cart T^T loses exactly the contaminant directions - so there is no
// converting one into the other. What can be done instead is to integrate it, which is exact, needs
// nothing from the file, and comes out in the SAME order and convention as the coefficients beside
// it. Note that it is SINGULAR for a file that stored spherical data, which is harmless for
// Mulliken partitioning and tr(DS) and fatal for anything that inverts or factorises it.
bool md_qm_publish_overlap(struct md_system_t* sys);

// S[a][b] over the basis's Cartesian AOs, row major, md_gto_basis_num_ao(basis) on a side.
// 'atom_coord_bohr' holds one position per atom the shell list indexes, in BOHR - the unit the
// exponents are stated in, not the Angstrom the system's own coordinates use.
// Returns the dimension written, 0 on failure.
size_t md_qm_compute_overlap(double* out, const struct md_gto_basis_t* basis, const dvec3_t atom_coord_bohr[], size_t num_atoms);

// orbital/{alpha,beta}/density, and orbital/{total,difference}/density when both spins are present,
// as VIRTUAL attributes computed from orbital/{alpha,beta}/{coefficient,occupation} - which must
// already be published, since this reads them back rather than taking them as arguments.
//
// Virtual rather than resident because a density is an exact function of the coefficients and the
// occupations beside it: a resident copy would be a second {A,A} matrix with nothing keeping it in
// step. See the provider notes in md_system.h; 'sys' is borrowed and must outlive the table.
bool md_qm_publish_orbital_densities(struct md_system_t* sys);

#ifdef __cplusplus
}
#endif
