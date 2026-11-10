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
struct md_allocator_i;

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
md_attribute_id_t md_qm_publish(struct md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit,
                                md_attribute_format_t format, const void* data, size_t byte_size);

// The same, for an attribute computed through a provider instead of one copied in.
md_attribute_id_t md_qm_publish_virtual(struct md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit,
                                        md_attribute_format_t format, const md_attribute_virtual_t* virt);

md_attribute_id_t md_qm_publish_scalar     (struct md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit, double value);                                  // rank 0
md_attribute_id_t md_qm_publish_series     (struct md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit, const double* values, size_t count);            // rank 1 {N}
md_attribute_id_t md_qm_publish_vec3_series(struct md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit, const dvec3_t* values, size_t count);           // rank 1 {N} x 3
md_attribute_id_t md_qm_publish_matrix     (struct md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit, const double* values, size_t rows, size_t cols); // rank 2 {R,C}
md_attribute_id_t md_qm_publish_str        (struct md_system_t* sys, md_str_t path, md_str_t label, md_str_t value);                                                   // rank 1 {1} string
md_attribute_id_t md_qm_publish_strings    (struct md_system_t* sys, md_str_t path, md_str_t label, const md_str_t* values, size_t count);                             // rank 1 {N} strings

// One COLUMN of an array of structs, gathered into the table's own storage. A record of six fields
// is six sibling paths over one index space, because a value has ONE type - so the transposition
// from a file's row layout happens once, here, rather than in every consumer.
md_attribute_id_t md_qm_publish_column(struct md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit,
                                       md_attribute_type_t type, const void* base, size_t stride, size_t count);

// The anchor of a vector group: rank 0, one 3 component value, constant over whatever index space
// the group's vector has. Angstrom, because it is a point in system space. NOT replicated to match
// the vector's shape - group members share an index space, not a shape.
md_attribute_id_t md_qm_publish_origin(struct md_system_t* sys, md_str_t path, dvec3_t origin);

// A second name for an attribute already published. Both names then read one datum - no copy, and a
// consumer of either is unaffected when the other appears or goes.
md_attribute_id_t md_qm_alias(struct md_system_t* sys, md_attribute_id_t target, md_str_t path);

// Publishes beta's copy of a per orbital series, or a SECOND NAME for alpha's when the two share
// storage. Comparing the POINTERS is what tells the cases apart, and it is the only test that gets
// restricted open shell right, where the orbitals are shared and the occupations are not.
md_attribute_id_t md_qm_publish_or_alias(struct md_system_t* sys, md_attribute_id_t alpha_id, md_str_t path, md_str_t label,
                                         md_unit_t unit, const double* alpha_values, const double* beta_values, size_t count);

// Builds an attribute path from a fixed group prefix and a name taken from a file. A '/' inside the
// name would silently introduce a group level in a namespace where the separator is the only
// structure there is, so it is folded to '_'. Returns an empty md_str_t when the name does not fit,
// which a caller treats as "skip this one" rather than as a reason to stop publishing.
md_str_t md_qm_attribute_path(char* buf, size_t cap, md_str_t group, md_str_t name);

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
// For a format that normalises a Cartesian shell once, against the (l,0,0) monomial, and uses that one
// constant for all (l+1)(l+2)/2 functions - so its xy function is not unit normalised and its xx one is
// (TREXIO, before its per function ao_normalization; Molden files from Turbomole, CFOUR 2.1 and PSI4
// <= 1.3.2). md_gto normalises every monomial in the shell. This is the ratio between the two,
// sqrt((2i-1)!!(2j-1)!!(2k-1)!! / (2l-1)!!), and it is 1 for every AO of an s or p shell, which is why
// a reader that gets the convention wrong still looks correct until the first d. NOT the Molden format's
// own convention: there each Cartesian function is normalised on its own, as in md_gto, and the
// coefficients are taken as they are.
double md_qm_cart_coeff_factor(uint32_t l, uint32_t cart_idx);

// Converts a whole [num_mo][n_sph] coefficient matrix into [num_mo][n_cart] (md_gto_sph_to_cart_orbitals), where
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

// qm/atom/{atomic_number,nuclear_charge,coordinate} - the atoms the calculation covered, in ITS order
// and at ITS geometry, which is not necessarily the system's atom set. Coordinates are taken in
// ANGSTROM, to match the system's own, and are published as such.
//
// nuclear_charge (e, f64) is the charge of each nucleus as the electrons of the calculation see it:
// the atomic number, less the core electrons an effective core potential replaces. It is what the
// density was solved against, so it is what an electrostatic potential, a dipole or a nuclear
// repulsion has to be computed from - with Z instead, the potential of an ECP calculation is off by
// n_core / r. Optional: NULL publishes the atomic numbers, which is right for an all electron
// calculation and the most a reader can say when its format does not state the charges. The
// attribute is published either way, so a consumer reads it and never chooses.
bool md_qm_publish_atoms(struct md_system_t* sys, const uint8_t atomic_number[], const double nuclear_charge[], const dvec3_t coord_angstrom[], size_t count);

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

// SYMMETRIC AO MATRICES. A density - or an overlap - is a symmetric N x N matrix in the AO basis, and is
// published one of two ways: square ({A,A}, or {S,A,A} indexed by state), or as its packed upper
// triangle (MD_ATTRIBUTE_FLAG_PACKED_SYMMETRIC, {A(A+1)/2} or {S,A(A+1)/2}) when it is large and read
// from disk on demand. These two read either one, so a consumer never has to know which it was given.
// 'slice' narrows 'attr' to ONE matrix. Both return N, or 0 on failure; with dst NULL they return N
// and extract nothing. cap is in elements.

// The packed upper triangle as float - row i holding columns i..N-1, N(N+1)/2 values - which is the
// layout the GTO density kernels take. A packed attribute is extracted straight into dst, so nothing
// larger than the triangle is ever made; a square one is extracted whole into scratch first.
size_t md_qm_extract_packed_symmetric_f32(float* dst, size_t cap, const struct md_attribute_t* attr, md_attribute_slice_t slice, md_attribute_io_t* io);

// The whole N x N matrix as double, row major: a packed attribute is extracted into scratch and its
// lower half mirrored from the upper; a square one is extracted as it is.
//
// Both through io (may be NULL): a density computed on demand (orbital/*/density, from the orbitals) is built on
// its executor (md_attribute_io_create).
size_t md_qm_extract_symmetric_f64(double* dst, size_t cap, const struct md_attribute_t* attr, md_attribute_slice_t slice, md_attribute_io_t* io);

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

// The orbitals an SCF density is a sum of squares of: D = sum_k w_k c_k c_k^T, for 'density_path' one of
// orbital/{alpha,beta,total,difference}/density, from orbital/{alpha,beta}/{coefficient,occupation}.
// What md_gto_iso_field_create_products takes (left = rows, weight = weights): exact, where a matrix
// would have to be factorized again. Orbitals of zero weight are left out, and an orbital both spins
// share appears once with its weights summed (a restricted total: twice the alpha occupation; a
// restricted difference: nothing). rows: count x num_ao, weights: count, both allocated from 'alloc'
// (free with sizeof(double) * count * num_ao and sizeof(double) * count). Returns count; 0 for another
// path, missing orbitals, or a density that is zero.
size_t md_qm_density_orbitals(double** out_rows, double** out_weights, size_t* out_num_ao, const struct md_system_t* sys,
                              md_str_t density_path, struct md_allocator_i* alloc);

// The orbitals of the AVERAGE LOCAL IONIZATION ENERGY (Sjoberg, Murray, Brinck and Politzer, Can. J. Chem. 68,
// 1440 (1990)),
//
//     I(r) = sum_k n_k |e_k| |phi_k(r)|^2 / sum_k n_k |phi_k(r)|^2
//
// over the occupied orbitals of both spins, n_k the occupation and e_k the orbital energy: where the least
// tightly held electrons are, the lowest values (on an isodensity surface, the sites an electrophile or a
// radical attacks). Both sums are of squares of the same orbitals, so this hands out the rows once and two
// sets of weights: weights, n_k (exactly md_qm_density_orbitals' for orbital/total/density: the denominator
// is the electron density), and energy_weights, n_k |e_k| in hartree (the numerator). What
// md_gto_iso_field_create_products takes, once with each.
//
// An orbital both spins share appears once with its weights summed: n_a + n_b, and n_a |e_a| + n_b |e_b| -
// a restricted open shell's two spins may state different energies for the same orbital. Reads
// orbital/{alpha,beta}/{coefficient,occupation,energy}.
//
// Canonical orbitals only: every occupation 0, 1 or 2 (within 1e-3). Natural orbitals (fractional
// occupations) have no orbital energies to weight by, and a file that gives them a column of energies
// anyway would give a plausible looking I(r) that means nothing; 0 is returned for them, and for a file
// without orbital energies (TREXIO may leave them out) or with an energy that is not finite. The energies
// of Kohn-Sham orbitals are not ionization energies, and I(r) from them lies below the Hartree-Fock one,
// as is usual.
//
// Where it means something: I(r) >= |e_HOMO| everywhere, and on an isodensity surface of 0.001-0.002 au it
// is the usual 0.3-0.6 hartree. Far out it is not |e_HOMO|, as with exact orbitals: a gaussian basis' tails
// are its most diffuse primitives', and the ratio there is of two vanishing sums (water 8.7 bohr out, rho
// 7e-12: 0.49 against |e_HOMO| 0.29). A consumer evaluates it where rho is substantial.
//
// Every array from 'alloc' (sizeof(double) * count * num_ao, and sizeof(double) * count each). Returns
// count; 0 on failure.
size_t md_qm_alie_orbitals(double** out_rows, double** out_weights, double** out_energy_weights, size_t* out_num_ao,
                           const struct md_system_t* sys, struct md_allocator_i* alloc);

// Whether md_qm_alie_orbitals can make them: the alpha orbitals' coefficients, occupations and energies are
// there and the occupations are those of canonical orbitals. Reads the occupations only.
bool md_qm_alie_available(const struct md_system_t* sys);

// The orbitals of the LOCAL ELECTRON ATTACHMENT ENERGY (Brinck, Carlqvist and Stenlid, J. Phys. Chem. A 120,
// 10023 (2016)),
//
//     E(r) = sum_{i >= LUMO, e_i < 0} e_i rho_i(r) / rho(r)
//
// the counterpart of the average local ionization energy for nucleophiles: the bound virtual orbitals (only
// they bind an added electron, by Janak's theorem) weighted by their energies, against the ground state
// density. Never positive; the lowest values mark where a nucleophile attacks or a halogen bond forms (a
// halogen's sigma hole). The numerator is handed out here, as rows and weights for
// md_gto_iso_field_create_products; the denominator is the total density (md_qm_density_orbitals,
// orbital/total/density).
//
// rho_i is the density the virtual would hold filled: per SPIN orbital, the weight is (1 - n) e for one with
// e < 0, summed over both spins. A closed shell's virtual so counts twice, 2 e |phi|^2 - the paper's eq. 10,
// E ~ 2 rho_LUMO e_LUMO / rho, and its Table S1, whose numbers come out only with the 2 - and a restricted and an
// unrestricted calculation of one closed shell give the same E(r). Unrestricted: each spin's virtuals once.
// Restricted open shell: a singly occupied orbital contributes its empty spin, when that spin's energy is below
// zero. An orbital both spins share is there once with the two shares summed.
//
// Canonical orbitals only: per spin occupations of 0 or 1 (within 1e-3), as the readers publish them. Nothing
// (0) when no virtual is bound: then E(r) is 0 everywhere, which is what Hartree-Fock virtuals, nearly always
// positive, give. The paper's level: B3LYP/6-31+G(d,p) (hybrid GKS-DFT), on the 0.004 au isodensity surface.
// Strongly method and basis dependent, and the virtuals are more diffuse than the occupied orbitals: far out
// the ratio grows without bound. A consumer evaluates it on an isodensity surface.
//
// Arrays from 'alloc' (sizeof(double) * count * num_ao and sizeof(double) * count). Returns count; 0 on failure
// or when no virtual is bound.
size_t md_qm_attachment_orbitals(double** out_rows, double** out_weights, size_t* out_num_ao, const struct md_system_t* sys,
                                 struct md_allocator_i* alloc);

// Whether md_qm_attachment_orbitals has orbitals to give: canonical orbitals with energies, at least one virtual
// bound. Reads the occupations and energies.
bool md_qm_attachment_available(const struct md_system_t* sys);

// The KINETIC ENERGY DENSITY of the density at density_path, one of orbital/{alpha,beta,total}/density,
//
//     tau(r) = 1/2 sum_k n_k |grad phi_k(r)|^2
//
// (the positive definite form, which the electron localization function and meta-GGA functionals take) as a sum
// of squares of orbitals over the DERIVATIVE BASIS of the system's basis (md_gto_basis_derivative): each orbital
// of md_qm_density_orbitals three times, its d/dx, d/dy and d/dz (md_gto_orbital_gradients), with weight n_k / 2.
// What md_gto_iso_field_create_products takes over that basis - in double with its gradient, or in a GPU set -
// and out_basis receives it (from 'alloc': md_gto_basis_free). rows: count x num_ao of the derivative basis, and
// weights: count, from 'alloc' (sizeof(double) * count * num_ao and sizeof(double) * count). Returns count, three
// per orbital; 0 for another path, missing orbitals or basis, or a basis with g functions, whose derivatives need
// h functions, more than is evaluated (MD_GTO_MAX_ANGULAR_MOMENTUM).
size_t md_qm_kinetic_orbitals(double** out_rows, double** out_weights, size_t* out_num_ao, struct md_gto_basis_t* out_basis,
                              const struct md_system_t* sys, md_str_t density_path, struct md_allocator_i* alloc);

// Whether md_qm_kinetic_orbitals can make them: the alpha orbitals and the basis are there, and the basis has
// no shell of the highest angular momentum evaluated.
bool md_qm_kinetic_available(const struct md_system_t* sys);

// The ELECTRON LOCALIZATION FUNCTION (Becke and Edgecombe, J. Chem. Phys. 92, 5397 (1990), as Savin et al. read it,
// Angew. Chem. Int. Ed. 31, 187 (1992)) at one point, from each spin's density rho_s, the square of its gradient
// g2_s = |grad rho_s|^2 and its kinetic energy density tau_s (md_qm_kinetic_orbitals):
//
//     D_s  = tau_s - g2_s / (8 rho_s)                  spin s's excess kinetic energy density, due to the Pauli
//                                                      principle: 0 where a single orbital makes its density
//     D0_s = 3/10 (6 pi^2)^(2/3) rho_s^(5/3)           the same for the uniform electron gas of that density
//     ELF  = 1 / (1 + chi^2),   chi = (D_a + D_b) / (D0_a + D0_b)
//
// The spin polarized form of Kohout and Savin (Int. J. Quantum Chem. 60, 875 (1996)), which Multiwfn uses. For a
// closed shell it is the usual ELF of the total density (md_qm_elf_closed_shell): D = tau - |grad rho|^2 / (8 rho),
// D0 = 3/10 (3 pi^2)^(2/3) rho^(5/3). In [0, 1]: 1 where one orbital or one pair alone makes the density (cores,
// bonds, lone pairs), 1/2 where it is like the electron gas, low between shells and pairs. A spin of no density
// adds nothing; D_s is clamped at 0 (tau_s >= g2_s / (8 rho_s) holds exactly, rounding aside).
//
// In the tails: D (both spins) plus d_min, as Multiwfn has it (function.f90, ELF_LOL: D = D + 1e-5 by default,
// settings.ini ELF_addminimal = 1). Far from the atoms D and D0 are a difference and a power of vanishing sums: with a
// gaussian basis' tails their ratio takes any value, and as the most diffuse orbital comes to make the density alone
// D goes to 0 and ELF to 1 in the vacuum around a molecule. With d_min it goes to 0 instead, smoothly, as D0 falls
// below d_min: MD_QM_ELF_D_MIN (1e-5) has ELF follow the exact one to near the 1e-3 isodensity (lower by 0.11 to 0.15
// there: d_min / D0 is 0.35 at rho 1e-3) and gone by rho 1e-4, as Multiwfn's (and the figures made with it). d_min also keeps ELF of float data from its noise
// (md_qm_elf_noise), which is far below it wherever ELF is not already 0. d_min 0: the function as defined, out to
// where the density vanishes. 0 where the density (both spins) is below rho_min.
#define MD_QM_ELF_D_MIN 1.0e-5
double md_qm_elf(double rho_a, double g2_a, double tau_a, double rho_b, double g2_b, double tau_b, double rho_min, double d_min);
double md_qm_elf_closed_shell(double rho, double g2, double tau, double rho_min, double d_min);

// How far ELF (d_min 0) can be trusted at a point, from float inputs: chi = D / D0 is a small difference of large terms
// over D0, and an error of float's epsilon in them (the basis, the coordinates, the sums), as any evaluation from
// float data has, moves chi by about epsilon * (tau_a + tau_b) / (D0_a + D0_b), which this returns (inf where D0 is 0).
// It grows smoothly outward as rho^(-2/3); measured against double, ELF's error is mostly within a few times it (mean
// a fifth of it, worst about 20 times): about 1e-2 by rho 1e-7.
double md_qm_elf_noise(double rho_a, double tau_a, double rho_b, double tau_b);

// The LOCALIZED ORBITAL LOCATOR (Schmider and Becke, J. Mol. Struct. THEOCHEM 527, 51 (2000)) from the same:
// t = (D0_a + D0_b) / (tau_a + tau_b), LOL = t / (1 + t), in [0, 1]. It reads tau without the gradient term, so it
// is smoother than ELF and its maxima are where ELF's are; high where the orbitals' gradients are small (cores,
// bonds, lone pairs). 0 below rho_min, as ELF.
double md_qm_lol(double rho_a, double tau_a, double rho_b, double tau_b, double rho_min);
double md_qm_lol_closed_shell(double rho, double tau, double rho_min);

// The DENSITY SCALE LENGTH L = rho / |grad rho| = 1 / |grad ln rho| (bohr), from the total density and the square of its
// gradient g2: the distance over which the density would change by a factor e at its rate of change there (the density
// gradient scale length of plasma physics). On an isodensity surface rho = c it is c / |grad rho|, and it is how far
// the surface moves along its normal per relative change of the density there: a change d rho / rho moves the surface
// by L d rho / rho. So mapped onto the surface it shows where the electron cloud's boundary is soft (large L: the
// surface moves far for a given compression or expansion of the density) and where it is stiff (small L, close to
// the nuclei). Of the exact density it goes to 1 / (2 sqrt(2 I)) in the far tails, for the first ionization energy I;
// of a gaussian basis' density it shrinks there instead, as 1 / (4 alpha r) of the most diffuse function. Infinite
// where the gradient vanishes (the critical points: on a surface only where it passes through one): capped at
// length_max (0: FLT_MAX). 0 where rho is below rho_min.
double md_qm_density_scale_length(double rho, double g2, double rho_min, double length_max);

// Whether the ground state density can be made of orbitals (md_qm_density_orbitals): the alpha orbitals and their
// occupations are there. Cheap: the attributes are looked up.
bool md_qm_density_available(const struct md_system_t* sys);

// ELF and LOL at points, in double, on the CPU, from the system's orbitals (md_qm_density_orbitals and
// md_qm_kinetic_orbitals of each spin, once for a closed shell): the reference the faster paths are checked
// against, and what a consumer without them takes. On the GPU nothing of ELF's own is needed: rho and its
// gradient from a set of the density over the basis, tau from a set of md_qm_kinetic_orbitals over the derivative
// basis, both by md_gto_iso_gpu_evaluate at the same points, then md_qm_elf per point; with both sets at the AO
// cutoff 1e-8 (md_gto_iso_gpu_set_cutoff), within 1e-5 of this where rho > 1e-3 (float, against double). xyz: count points, stride bytes apart (0 = packed float3), in
// bohr, as atom_xyz (one per atom of the basis, bohr, atom_xyz_stride bytes apart, 0 = packed). Either output may
// be NULL. rho_min, d_min: as md_qm_elf. Builds the fields on every call: for repeated evaluation, hold them. Returns
// count, 0 on failure.
size_t md_qm_elf_evaluate(double* out_elf, double* out_lol, const float* xyz, size_t count, size_t stride,
                          const struct md_system_t* sys, const float* atom_xyz, size_t atom_xyz_stride, double rho_min, double d_min);

// orbital/{total,difference}/density (alpha + beta, alpha - beta) as virtual attributes over whatever
// orbital/{alpha,beta}/density are: for a producer that publishes those two itself (read from a file
// that stores them, say). Where beta is an alias of alpha, alpha is read once for the total and not
// at all for the difference. md_qm_publish_orbital_densities ends with this.
bool md_qm_publish_orbital_density_combinations(struct md_system_t* sys);

#ifdef __cplusplus
}
#endif
