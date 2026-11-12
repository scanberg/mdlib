#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_str.h>

// Reading files into a system: what every reader takes, and the table of the formats.
//
// A reader of a format has some of these, all with the same signature, all taking an optional md_load_desc_t:
//
//   md_<fmt>_system_init_from_file(sys, state, path, desc)    the system (atoms, bonds, ...) and its coordinates
//   md_<fmt>_system_init_from_str (sys, state, str,  desc)    the same of text in memory
//   md_<fmt>_system_supplement_from_file(sys, path, desc)     added to a loaded system: its topology, groups, or a
//                                                              series along a run (desc->run)
//   md_<fmt>_system_publish_run(sys, path, run, flags)        frames streamed from the file as a run (md_run.h)
//
// The table below says which of them a format has, what it is called, what it is told by (its extension, and where that
// is not enough its contents) and what it holds, for an application to dispatch on, rather than every application
// keeping a table of its own.

struct md_system_t;
struct md_system_state_t;
struct md_executor_t;

enum {
    MD_LOAD_FLAG_NONE       = 0,
    MD_LOAD_FLAG_ALL_MODELS = 1,    // A file of several models (pdb): the atoms of all of them, rather than of the first
};

typedef struct md_load_desc_t {
    const struct md_executor_t* executor;   // Reads in parallel where a reader can (gro, pdb, mmcif, vlx); NULL: serially
    uint32_t flags;                         // MD_LOAD_FLAG_*
    md_str_t format;                        // What the file does not say itself: the atom style of a LAMMPS data file
                                            // ("full", "atomic", ... or the columns); empty: told from the file
    md_str_t run;                           // The run a series sampled along a trajectory joins ("run/<name>"): a
                                            // supplement which is MD_IO_FLAG_TEMPORAL fails without one
} md_load_desc_t;

// ### FORMATS ###

typedef enum md_io_format_t {
    MD_IO_FORMAT_UNKNOWN = 0,
    MD_IO_FORMAT_PDB,
    MD_IO_FORMAT_GRO,
    MD_IO_FORMAT_XYZ,
    MD_IO_FORMAT_XMOL,
    MD_IO_FORMAT_ARC,
    MD_IO_FORMAT_MMCIF,
    MD_IO_FORMAT_LAMMPS_DATA,
    MD_IO_FORMAT_LAMMPS_TRJ,
    MD_IO_FORMAT_XTC,
    MD_IO_FORMAT_TRR,
    MD_IO_FORMAT_DCD,
    MD_IO_FORMAT_VLX,               // VeloxChem (.h5)
    MD_IO_FORMAT_MOLDEN,
    MD_IO_FORMAT_TREXIO,
    MD_IO_FORMAT_H5MD,
    MD_IO_FORMAT_ITP,               // GROMACS topology (.itp, .top)
    MD_IO_FORMAT_TPR,
    MD_IO_FORMAT_NDX,
    MD_IO_FORMAT_EDR,
    MD_IO_FORMAT_XVG,
    MD_IO_FORMAT_CSV,
    MD_IO_FORMAT_ASE_TRAJ,
    MD_IO_FORMAT_AMBER_PRMTOP,
    MD_IO_FORMAT_AMBER_RESTART,     // ASCII (.rst7, .inpcrd) or NetCDF (.ncrst)
    MD_IO_FORMAT_AMBER_NETCDF,      // A trajectory
    MD_IO_FORMAT_AMBER_ENERGY,      // mdout or mden
    MD_IO_FORMAT_PSF,
    MD_IO_FORMAT_NAMD_COOR,
    MD_IO_FORMAT_COUNT
} md_io_format_t;

// What a format holds and how it joins what is loaded
enum {
    MD_IO_FLAG_SYSTEM            = 1u << 1,     // A system of its own (system_init)
    MD_IO_FLAG_RUN               = 1u << 2,     // Frames, published as a run (publish_run)
    MD_IO_FLAG_MM                = 1u << 5,     // Molecular mechanics data
    MD_IO_FLAG_QM                = 1u << 6,     // Quantum mechanics data
    MD_IO_FLAG_SUPPLEMENT        = 1u << 7,     // Adds to a loaded system (supplement)
    MD_IO_FLAG_TOPOLOGY          = 1u << 8,     // A complete topology: its bonds are not to be inferred
    MD_IO_FLAG_TEMPORAL          = 1u << 9,     // A supplement sampled along a trajectory: it joins a run (desc->run)
    MD_IO_FLAG_NEEDS_COORDINATES = 1u << 11,    // A topology without coordinates (prmtop, psf): they come from another file
    MD_IO_FLAG_COORDINATES       = 1u << 12,    // Coordinates without atoms (a restart): for a topology (state_init)
};

typedef struct md_io_format_info_t {
    md_str_t name;          // "Gromacs Structure (gro)"
    md_str_t ext;           // The extension it is known by ("gro"), also what a format is stored as
    uint32_t flags;         // MD_IO_FLAG_*

    // The readers of the format; NULL where it has none, all of them where the format is not built in (VeloxChem,
    // TREXIO and H5MD without HDF5)
    bool   (*system_init)(struct md_system_t* sys, struct md_system_state_t* state, md_str_t path, const md_load_desc_t* desc);
    bool   (*supplement) (struct md_system_t* sys, md_str_t path, const md_load_desc_t* desc);
    bool   (*publish_run)(struct md_system_t* sys, md_str_t path, md_str_t run, uint32_t flags);
    // The coordinates of a file of coordinates (MD_IO_FLAG_COORDINATES, or the first frame of a trajectory) for a
    // system of num_atoms atoms
    bool   (*state_init) (struct md_system_state_t* state, md_str_t path, size_t num_atoms);
    // The atom count its header gives, without reading the rest; NULL where a format does not say up front
    size_t (*atom_count) (md_str_t path);
} md_io_format_info_t;

// The format's entry. MD_IO_FORMAT_UNKNOWN, and anything past the end, is an entry of nothing.
const md_io_format_info_t* md_io_format_info(md_io_format_t format);

// Whether the format is built in: false for the formats of HDF5 in a build without it
bool md_io_format_available(md_io_format_t format);

// The format of an extension, without the dot and of any case (also .top, .parm7, .inpcrd, .ncrst, .ncdf, .mden, ...):
// where several share it (.top, .nc, .coor) the one it mostly is. UNKNOWN for none.
md_io_format_t md_io_format_from_ext(md_str_t ext);

// The format of a file: by its extension, and by its contents where the extension does not tell or is shared - Amber
// names its files without one, .top is GROMACS or Amber, .nc a trajectory or a restart, .coor binary or PDB, .h5
// VeloxChem or TREXIO, and a Molden file is often .molden.input. UNKNOWN for none.
md_io_format_t md_io_format_from_file(md_str_t path);
