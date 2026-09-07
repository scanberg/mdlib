#pragma once

#include <stdbool.h>
#include <core/md_str.h>

struct md_system_t;
struct md_system_state_t;

// VeloxChem (.h5 / .out) reader.
//
// THE ENTIRE SURFACE IS THE THREE ENTRY POINTS BELOW. There is no vlx object, no accessor for one,
// and no type describing what the file held: a file is parsed into an md_system_t and everything it
// carried is read back out of that system - its atoms, its state, and its ATTRIBUTE TABLE. The
// reader's own representation is an implementation detail of md_vlx.c which does not outlive the
// call, so nothing downstream can depend on it and no consumer has to keep a parser alive to ask a
// question. See the ATTRIBUTES section of md_system.h for how the table is read.
//
// What a parsed file publishes, and where:
//
//   ATOMS AND STATE      the system's own atom list and coordinates (md_system_t / md_system_state_t)
//
//   vlx/molecular_charge, vlx/nuclear_repulsion_energy,      rank 0, single values
//   vlx/spin_multiplicity, vlx/electron_count/{alpha,beta}
//   vlx/basis_set, vlx/dft_functional, vlx/potfile           rank 1 {1} strings
//   vlx/scf/type, vlx/rsp/type, vlx/opt/type                 rank 1 {1} strings: WHICH kind of run
//                                                            this was, lowercased enumerator names.
//                                                            Absent when the file does not say,
//                                                            which is what a consumer treats as
//                                                            unknown. Text rather than an integer
//                                                            so nothing stored depends on an
//                                                            internal enum ordering.
//   vlx/scf/history/*                                        {I} per SCF iteration
//   vlx/scf/orbital/{alpha,beta}/*                           {M} per molecular orbital, aliased to
//                                                            orbital/{alpha,beta}/* as well
//   vlx/rsp/{oscillator,rotatory}_strength                   {S} per excited state
//   vlx/rsp/frequency, vlx/rsp/{cpp,tpa}/*                   {F} over the response frequencies
//   vlx/rsp/rixs/*                                           {P}, {C}, and {F,P} for the 2D maps,
//                                                            plus scattering_amplitude_{re,im}
//                                                            {F,P,3,3} - the complex tensor as
//                                                            two real attributes, because a
//                                                            value has one type
//   vlx/rsp/nto/lambda                                       {S,Lmax}, ragged rows zero padded
//   vlx/rsp/nto/{particle,hole}/coefficient                  {S,Lmax,A}
//   vlx/rsp/transition_density/*                             {S,A,A}, computed on demand
//   vlx/rsp/solution_matrix, vlx/rsp/num_{core,valence,       the raw response eigenvectors and the
//   virtual}                                                  occupied/virtual split they are
//                                                             indexed by. Not meant for direct
//                                                             consumption - they are what lets the
//                                                             transition densities above be rebuilt
//                                                             on demand with no reader in reach
//   vlx/xps/*                                                {C} per core-hole state, one path per
//                                                            field of the record
//   vlx/vib/*                                                {D} per normal mode, plus
//                                                            vlx/vib/raman_activity {E,D} - one row
//                                                            per external frequency
//   vlx/opt/{energy,coordinate}                              {P} per optimisation step;
//                                                            vlx/opt/state_index and
//                                                            vlx/opt/irc_ts_index are rank 0
//   vlx/density_property/<dataset name>                      {A,A} density properties as the file
//                                                            carried them
//   atom/<dataset name>                                      per atom properties from the file
//   orbital/{alpha,beta}/coefficient                         {M,A} MO coefficients, Cartesian AO
//                                                            order - see the AO CONVENTION block in
//                                                            md_gto.h
//   orbital/{alpha,beta,total,difference}/density            {A,A}, computed on demand from the
//                                                            coefficients and occupations
//   basis/shell/*, basis/primitive/*                         the GTO basis, in the format neutral
//                                                            layout md_gto_basis_extract_attributes
//                                                            reads
//   basis/overlap                                            {A,A} the AO overlap, in the same
//                                                            Cartesian order and convention as the
//                                                            coefficients, and singular with it
//   qm/atom/{atomic_number,coordinate}                       the QM ATOM DOMAIN - the atoms this
//                                                            calculation covered, in its order and
//                                                            at its geometry. NOT the system's
//                                                            atoms; basis/shell/atom_index indexes
//                                                            THIS space, and qm/atom/system_index
//                                                            is the only bridge across
//   qm/atom/normal_mode                                      {M,N} x 3, one displacement per QM
//                                                            atom per mode
//   dipole/{ground_state,electric_transition,                each a vector and an origin of the same
//           magnetic_transition,velocity_transition}/*       shape, anchored at the centre of charge
//
// Most are read with md_attribute_extract_f32; the multi axis ones by building an
// md_attribute_slice_t, asking md_attribute_slice_count how big it is and extracting into that.
//
// A block the file does not contain publishes nothing, so a consumer asks the attribute table what
// is there rather than asking a reader what it parsed.
//
// The vlx/ prefix is deliberate for the format specific tree: these are one program's output, and a
// path is a promise, so a quantity moves to a format neutral name once a second loader produces the
// same thing rather than in anticipation of one. atom/, orbital/, basis/, qm/ and dipole/ are
// already format neutral conventions which md_system.h documents, so those land there directly.
//
// Publishing is idempotent: a path already present is replaced, so loading another file into the
// same system leaves no stale series behind.

#ifdef __cplusplus
extern "C" {
#endif

// Parse a VeloxChem .h5 or .out file into a system, replacing its atoms, its state and its
// attribute table with what the file carried.
bool md_vlx_system_init_from_file(struct md_system_t* sys, struct md_system_state_t* state, str_t filename);

// For a file whose atoms are a subset of an ALREADY loaded system - see
// md_vlx_system_is_file_supplemental, which is how a caller decides between the two. Parses and
// publishes onto that system, leaving its atoms and its state untouched where
// md_vlx_system_init_from_file would replace them.
bool md_vlx_system_supplement_from_file(struct md_system_t* sys, str_t filename);

bool md_vlx_system_is_file_supplemental(const struct md_system_t* sys, str_t filename);

#ifdef __cplusplus
}
#endif
