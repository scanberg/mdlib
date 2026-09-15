#pragma once

#include <stdbool.h>
#include <core/md_str.h>

struct md_system_t;
struct md_system_state_t;

// Molden (.molden / .molden.input / .mold) reader.
//
// THE ENTIRE SURFACE IS THE ENTRY POINTS BELOW. There is no molden object, no accessor for one, and
// no type describing what the file held: a file is parsed into an md_system_t and everything it
// carried is read back out of that system - its atoms, its state, and its ATTRIBUTE TABLE. The
// reader's own representation is an implementation detail of md_molden.c which does not outlive the
// call. See the ATTRIBUTES section of md_system.h for how the table is read, md_qm.h for the
// format neutral paths this shares with every other QM reader, and md_vlx.h for the same contract
// stated at length.
//
// What a parsed file publishes, and where:
//
//   ATOMS AND STATE          the system's own atom list and coordinates (md_system_t /
//                            md_system_state_t), in Angstrom whichever unit [Atoms] was written in
//
//   molden/title             rank 1 {1} strings: [Title] and, when the file names one, the program
//   molden/program           that wrote it
//   molden/ao_convention     rank 1 {1} string, "6D10F15G" style: WHICH angular convention the
//                            coefficients below were read under. The published basis is Cartesian
//                            either way - this records what the file said, because nothing else in
//                            the table can be used to work it out afterwards
//
//   qm/atom/{atomic_number,  the QM ATOM DOMAIN - the atoms this calculation covered, in its order
//            coordinate}     and at its geometry. basis/shell/atom_index indexes THIS space
//
//   basis/shell/*,           the GTO basis in the format neutral layout
//   basis/primitive/*        md_gto_basis_extract_attributes reads
//   basis/overlap            {A,A} the AO overlap, integrated from that basis - Molden carries no
//                            integrals of its own
//
//   orbital/{alpha,beta}/    {M} per molecular orbital. 'symmetry' is the [MO] Sym= label and is
//     {energy,occupation,    text, so it is a string attribute and not an index into anything
//      symmetry}
//   orbital/{alpha,beta}/    {M,A} MO coefficients, Cartesian AO order - see the AO CONVENTION
//     coefficient            block in md_gto.h. A file which stored spherical coefficients is
//                            converted once, here
//   orbital/{alpha,beta,     {A,A}, computed on demand from the coefficients and occupations
//     total,difference}/
//     density
//
//   molden/vib/frequency,    {D} per normal mode, from [FREQ] and [INT]
//   molden/vib/ir_intensity
//   qm/atom/normal_mode      {D,N} x 3, one displacement per QM atom per mode, from
//                            [FR-NORM-COORD]
//
// A block the file does not contain publishes nothing, so a consumer asks the attribute table what
// is there rather than asking a reader what it parsed.
//
// The molden/ prefix is for what is specific to this format; qm/, basis/ and orbital/ are the
// format neutral conventions md_system.h and md_qm.h document, so those land there directly.
//
// WHAT THE FORMAT DOES NOT SAY, and why the reader does not guess:
//
//   - The ANGULAR CONVENTION defaults to CARTESIAN (6D, 10F, 15G) and is made spherical by the
//     [5D], [5D7F], [5D10F], [7F] and [9G] markers, which is what the format specifies. A writer
//     which emits spherical coefficients and forgets the marker produces a file nothing can read
//     correctly; the convention actually used is published as molden/ao_convention so a consumer
//     can at least see what was assumed.
//   - A Molden file carries NO CHARGE and NO SPIN MULTIPLICITY, so neither is published. The
//     occupations are what a consumer counts electrons from.
//   - AO coefficients are stated against NORMALISED contracted functions, which is what every
//     writer means and what the primitive coefficients are scaled to here - see md_qm.h. A file
//     whose contractions are deliberately unnormalised is not distinguishable from one whose are,
//     and is read as if they were.

#ifdef __cplusplus
extern "C" {
#endif

// Parse a Molden file into a system, replacing its atoms, its state and its attribute table with
// what the file carried.
bool md_molden_system_init_from_file(struct md_system_t* sys, struct md_system_state_t* state, str_t filename);

// The same, for a file already in memory.
bool md_molden_system_init_from_str(struct md_system_t* sys, struct md_system_state_t* state, str_t str);

// Whether the file looks like Molden: it opens, and the first non-empty line is the [Molden Format]
// header. Cheap, and meant for a caller deciding which reader to hand a file to - the extension
// alone does not settle it, since .molden.input and .mold are both in use and .input is not.
bool md_molden_file_is_molden(str_t filename);

#ifdef __cplusplus
}
#endif
