#pragma once

#include <stdbool.h>
#include <core/md_str.h>

struct md_system_t;
struct md_system_state_t;

// TREXIO (.trexio / .h5) reader, for the HDF5 back end.
//
// THE ENTIRE SURFACE IS THE ENTRY POINTS BELOW. There is no trexio object, no accessor for one, and
// no type describing what the file held: a file is parsed into an md_system_t and everything it
// carried is read back out of that system - its atoms, its state, and its ATTRIBUTE TABLE. See the
// ATTRIBUTES section of md_system.h for how the table is read, md_qm.h for the format neutral paths
// this shares with every other QM reader, and md_vlx.h for the same contract stated at length.
//
// TREXIO's HDF5 back end is a group per section, scalars as group attributes and arrays as
// datasets, which is a thin enough layout to read with HDF5 directly. This reader therefore does
// NOT link against libtrexio: mdlib already depends on HDF5 for the VeloxChem reader, and a second
// dependency buys only the TEXT back end, which is a directory of files rather than a file and is
// not what anybody hands to a viewer. A file this reader cannot open is declined, not guessed at.
//
// What a parsed file publishes, and where:
//
//   ATOMS AND STATE          the system's own atom list and coordinates, converted from the bohr
//                            TREXIO stores to the Angstrom a system uses
//
//   trexio/metadata/*        rank 1 {1} strings, whatever the metadata group carried
//   trexio/basis_type,       rank 1 {1} strings
//   trexio/mo_type
//   trexio/ao_convention     rank 1 {1} string, "spherical" or "cartesian": WHICH angular
//                            convention the file's coefficients were in. The published basis is
//                            Cartesian either way, so nothing else in the table records it
//   trexio/electron_count/   rank 0
//     {total,up,down}
//   trexio/nuclear_repulsion_energy   rank 0, when the file states it
//   trexio/mo/{occupation,   {M} as the FILE states them, before the split into spin channels
//     spin}                  described below - so nothing the file said is lost
//
//   qm/atom/{atomic_number,  the QM ATOM DOMAIN - the atoms this calculation covered, in its order
//            coordinate}     and at its geometry. basis/shell/atom_index indexes THIS space
//
//   basis/shell/*,           the GTO basis in the format neutral layout
//   basis/primitive/*        md_gto_basis_extract_attributes reads
//   basis/overlap            {A,A} the AO overlap, INTEGRATED from the published basis rather than
//                            read from the file's own ao_1e_int - see md_qm_publish_overlap for why
//                            a spherical overlap cannot be converted into a Cartesian one
//
//   orbital/{alpha,beta}/    {M} per molecular orbital
//     {energy,occupation,
//      symmetry}
//   orbital/{alpha,beta}/    {M,A} MO coefficients, Cartesian AO order
//     coefficient
//   orbital/{alpha,beta,     {A,A}, computed on demand from the coefficients and occupations
//     total,difference}/
//     density
//
// SPIN CHANNELS. TREXIO states one occupation per molecular orbital and an optional mo_spin column.
// Where that column is absent or all zero the calculation is restricted, the occupation is the
// TOTAL over both spins, and what is published is HALF of it in each channel with beta an alias of
// alpha - which is the per channel convention md_vlx.h uses and what makes the total and difference
// densities come out right. The file's own numbers are published unaltered under trexio/mo/, so a
// consumer that wants them as stated has them. A restricted OPEN shell calculation written with no
// mo_spin column cannot be told apart from a closed shell one and is split the same way, which is
// wrong for its singly occupied orbital; a file that fills in mo_spin is read exactly.
//
// NORMALISATION. basis_prim_factor, basis_shell_factor and ao_normalization are applied as TREXIO
// defines them, and the contraction is then rescaled onto md_gto's convention - see md_qm.h. A file
// whose radial functions are not Gaussian (basis_type other than "Gaussian", or a non zero
// basis_r_power) is DECLINED rather than read as if they were.

#ifdef __cplusplus
extern "C" {
#endif

// Parse a TREXIO file into a system, replacing its atoms, its state and its attribute table with
// what the file carried.
bool md_trexio_system_init_from_file(struct md_system_t* sys, struct md_system_state_t* state, str_t filename);

// Whether the file is an HDF5 file holding a TREXIO nucleus group. Cheap, and meant for a caller
// deciding which reader to hand a file to: .h5 is shared with several other formats, so the
// extension does not settle it.
bool md_trexio_file_is_trexio(str_t filename);

#ifdef __cplusplus
}
#endif
