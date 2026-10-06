#pragma once

#include <core/md_str.h>

struct md_system_t;
struct md_system_state_t;

// H5MD (.h5md), HDF5 for molecular data: https://www.nongnu.org/h5md/h5md.html, version 1.x.
//
// One file holds a system and its trajectory, so there are two entry points and nothing else: the
// system is read the way a structure file is, and the trajectory is published as a run the way a
// trajectory file is. There is no h5md object and no type describing the file; everything it
// carries is read back out of the system.
//
// WHAT IS READ. The H5MD core, as the specification defines it, and on top of it what GROMACS
// (2026, mdrun -o traj.h5md) writes beyond the core:
//
//   /h5md                           required, version 1.x; its author and creator are published
//   /h5md/modules/units             the 'unit' attribute on a dataset, parsed as the module defines it
//                                   ("nm", "nm ps-1", "kJ mol-1 nm-1", "10+3 m"). A position or box
//                                   without one is taken to be in Angstrom; anything else without one
//                                   is published without a unit.
//   /particles/<group>              ONE group is read: the one with the most particles. A particle
//                                   group is a subset of the simulation and they may overlap, so they
//                                   are not concatenated. GROMACS writes one, "system".
//     box                           dimension 3, boundary, edges: fixed or per frame, cuboid [3] or
//                                   triclinic [3][3] with the box vectors as rows. A dimension whose
//                                   boundary is "none" is not periodic in the cell.
//     position                      [F][N][3] per frame, or [N][3] for a structure without frames
//     velocity, force, image,       any other per particle element; see the run below
//     mass, charge, species, id     the topology when there is no GROMACS module (see SYSTEM)
//   /connectivity/<name>            [B][2] pairs of the group, as bonds. A dataset with a
//                                   particles_group reference is used when it names the group read;
//                                   one without is used when the file has a single particle group,
//                                   which is how GROMACS writes /connectivity/bonds.
//   /observables/<name>             published with the run
//
//   /h5md/modules/gromacs_topology  GROMACS' molecule types (particle names, residues, masses,
//                                   charges, atomic numbers) and molecule blocks
//
// A time dependent element is a group holding value, step and optionally time, and both ways the
// specification stores step and time are read: a dataset per frame, or a scalar increment with an
// offset attribute. Frames are identified by step, which is what the specification makes exact.
//
// WHAT IS DECLINED, not guessed at: an H5MD root below the file root, a major version other than 1,
// a particle group with an id that changes over time (particles reordered, inserted or removed
// between frames), and positions in other than three dimensions.
//
// SYSTEM. With the GROMACS module the system is the molecule types laid out by the blocks, built
// the way md_tpr.h builds one - same atom types, same residue numbering, the bonds complete (of origin
// MD_BOND_ORIGIN_TOPOLOGY) - so a system read from the trajectory and one read from the run input it
// was simulated from agree. What the module does not carry is missing: force field type names,
// non-bonded parameters and particle types (an atom without mass is taken as a virtual site).
// Without the module the core elements are all there is: species names the atom types (the names
// of an enumeration, the number of an integer species), mass their mass, and an element is taken
// from the name, or from the mass for an integer species. Neither case infers bonds; a caller
// wanting them for a file without connectivity infers them itself.
//
// The coordinates and the cell are those of the first frame. Published:
//
//     atom/charge                   F32 {N}      when the file carries charges
//     atom/velocity, atom/force,    F32 {N} cD   when the file carries them WITHOUT frames; with
//     atom/image                                 frames they belong to the run
//     h5md/version                  I32 {2}
//     h5md/{author,creator}/...     STR          name, email; name, version
//
// RUN. The group's positions as a run, in the layout md_run_publish describes, and every other
// per particle element beside them:
//
//     <run>/time, step              the position frames. Without a time dataset the time is the
//                                   step itself, without a unit
//     <run>/unitcell                F32 {F,3,3}  Angstrom, the box at each position frame
//     <run>/atom/position           F32 {F,N} c3 Angstrom, VIRTUAL
//     <run>/source/path             STR          the file
//     <run>/source/offset, size     I64 {F}      where each frame of position/value is in the file,
//                                                -1 and 0 where it is not one contiguous run of bytes
//
//     <run>/atom/<E>                F32 {F,N} cK VIRTUAL, E sampled at the position steps: velocity,
//                                                force, image, ... in the file's own unit
//     <run>/h5md/particles/<E>/...  E sampled at other steps (GROMACS' nstvout and nstfout against
//                                   nstxout): its own time and step, and its values at atom/<E>. An
//                                   element without a time of its own is given one from its step,
//                                   through the run's time step when the run has a constant one
//     <run>/h5md/particles/<E>/source/offset     I64, where each frame of E is, as above
//     <run>/h5md/observables/<path>/value        F64, the observable as the file has it, with its
//     <run>/h5md/observables/<path>/{time,step}  own frame axis when it is time dependent
//
// READING. GROMACS writes each frame as one uncompressed chunk, and a frame that is one contiguous
// run of bytes - contiguous storage, or chunks spanning whole frames without filters - is read like
// any other trajectory: through the io it is handed, straight from the file, with no HDF5 call and
// safely from any number of threads. Anything else (compressed chunks, integer images, a type that
// is not IEEE float) is read through HDF5, one reader at a time, since HDF5 is not thread safe. The
// providers' user_data is their own; they read sys, borrowed, so the run goes before the system.

#ifdef __cplusplus
extern "C" {
#endif

// Read the system: its atoms, its bonds, the first frame's coordinates and cell, and what is
// published above, replacing what the system held.
bool md_h5md_system_init_from_file(struct md_system_t* sys, struct md_system_state_t* state, str_t filename);

// Publish the trajectory as the run 'run'. The particle group read must hold as many particles as
// the system has atoms. flags are md_run_flag_t; none applies, as HDF5 keeps its own index and
// there is no '<file>.cache' to write.
bool md_h5md_system_publish_run(struct md_system_t* sys, str_t filename, str_t run, uint32_t flags);

#ifdef __cplusplus
}
#endif
