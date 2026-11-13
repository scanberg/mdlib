#include <md_io.h>

#include <core/md_common.h>
#include <core/md_path.h>

#include <md_system.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_xyz.h>
#include <md_mmcif.h>
#include <md_lammps.h>
#include <md_xtc.h>
#include <md_trr.h>
#include <md_dcd.h>
#include <md_molden.h>
#include <md_itp.h>
#include <md_tpr.h>
#include <md_ndx.h>
#include <md_edr.h>
#include <md_xvg.h>
#include <md_csv.h>
#include <md_ase_traj.h>
#include <md_amber_prmtop.h>
#include <md_amber_rst7.h>
#include <md_amber_netcdf.h>
#include <md_amber_energy.h>
#include <md_psf.h>
#include <md_namdbin.h>

#if MD_HDF5
#include <md_vlx.h>
#include <md_trexio.h>
#include <md_h5md.h>
#endif

// Where a format needs more than its reader's function

// An index file names the groups an energy file of the run, loaded before it, used for its energy and T-coupling groups
static bool ndx_supplement(md_system_t* sys, md_str_t path, const md_load_desc_t* desc) {
    if (!md_ndx_system_supplement_from_file(sys, path, desc)) return false;
    if (desc && !md_str_empty(desc->run)) md_edr_system_infer_group_roles_of_run(sys, desc->run);
    return true;
}

// An Amber restart is ASCII or NetCDF
static bool amber_restart_state_init(md_system_state_t* state, md_str_t path, size_t num_atoms) {
    if (md_amber_netcdf_file_kind(path) == MD_AMBER_NETCDF_RESTART) return md_amber_netcdf_state_init_from_file(state, path, num_atoms);
    return md_amber_rst7_state_init_from_file(state, path, num_atoms);
}

static size_t amber_restart_atom_count(md_str_t path) {
    const size_t n = md_amber_netcdf_file_atom_count(path);
    return n ? n : md_amber_rst7_file_atom_count(path);
}

#define S(x) MD_STR_INIT(x)

static const md_io_format_info_t formats[MD_IO_FORMAT_COUNT] = {
    [MD_IO_FORMAT_UNKNOWN] = { S("Undefined"), S(""), 0 },
    [MD_IO_FORMAT_PDB] = { S("Standard Protein Data Bank (pdb)"), S("pdb"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_RUN | MD_IO_FLAG_MM,
        .system_init = md_pdb_system_init_from_file, .publish_run = md_pdb_system_publish_run },
    [MD_IO_FORMAT_GRO] = { S("Gromacs Structure (gro)"), S("gro"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_MM,
        .system_init = md_gro_system_init_from_file },
    [MD_IO_FORMAT_XYZ] = { S("xyz (xyz)"), S("xyz"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_RUN | MD_IO_FLAG_MM,
        .system_init = md_xyz_system_init_from_file, .publish_run = md_xyz_system_publish_run },
    [MD_IO_FORMAT_XMOL] = { S("xyz (xmol)"), S("xmol"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_RUN | MD_IO_FLAG_MM,
        .system_init = md_xyz_system_init_from_file, .publish_run = md_xyz_system_publish_run },
    [MD_IO_FORMAT_ARC] = { S("xyz (arc)"), S("arc"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_RUN | MD_IO_FLAG_MM,
        .system_init = md_xyz_system_init_from_file, .publish_run = md_xyz_system_publish_run },
    [MD_IO_FORMAT_MMCIF] = { S("PDBx/mmCIF (cif)"), S("cif"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_MM,
        .system_init = md_mmcif_system_init_from_file },
    [MD_IO_FORMAT_LAMMPS_DATA] = { S("LAMMPS (data)"), S("data"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_MM,
        .system_init = md_lammps_system_init_from_file },
    [MD_IO_FORMAT_LAMMPS_TRJ] = { S("LAMMPS Trajectory (lammpstrj)"), S("lammpstrj"), MD_IO_FLAG_RUN | MD_IO_FLAG_MM,
        .publish_run = md_lammps_system_publish_run },
    [MD_IO_FORMAT_XTC] = { S("Gromacs Compressed Trajectory (xtc)"), S("xtc"), MD_IO_FLAG_RUN | MD_IO_FLAG_MM,
        .publish_run = md_xtc_system_publish_run, .atom_count = md_xtc_file_atom_count },
    [MD_IO_FORMAT_TRR] = { S("Gromacs Lossless Trajectory (trr)"), S("trr"), MD_IO_FLAG_RUN | MD_IO_FLAG_MM,
        .publish_run = md_trr_system_publish_run, .atom_count = md_trr_file_atom_count },
    [MD_IO_FORMAT_DCD] = { S("DCD Trajectory (dcd)"), S("dcd"), MD_IO_FLAG_RUN | MD_IO_FLAG_MM,
        .publish_run = md_dcd_system_publish_run, .atom_count = md_dcd_file_atom_count },
    // VeloxChem: a system, or the QM data of the atoms of one loaded (md_vlx_system_is_file_supplemental). RUN as it is
    // read along with the trajectories of a system; it publishes none of its own.
    [MD_IO_FORMAT_VLX] = { S("VeloxChem (h5)"), S("h5"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_RUN | MD_IO_FLAG_MM | MD_IO_FLAG_QM,
#if MD_HDF5
        .system_init = md_vlx_system_init_from_file, .supplement = md_vlx_system_supplement_from_file,
#endif
    },
    [MD_IO_FORMAT_MOLDEN] = { S("Molden (molden)"), S("molden"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_QM,
        .system_init = md_molden_system_init_from_file },
    [MD_IO_FORMAT_TREXIO] = { S("TREXIO (trexio)"), S("trexio"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_QM,
#if MD_HDF5
        .system_init = md_trexio_system_init_from_file,
#endif
    },
    // The system and its positions as a run in one file. TOPOLOGY: GROMACS writes its bonds (a file without any has them
    // inferred)
    [MD_IO_FORMAT_H5MD] = { S("H5MD (h5md)"), S("h5md"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_RUN | MD_IO_FLAG_MM | MD_IO_FLAG_TOPOLOGY,
#if MD_HDF5
        .system_init = md_h5md_system_init_from_file, .publish_run = md_h5md_system_publish_run,
#endif
    },
    [MD_IO_FORMAT_ITP] = { S("Gromacs Topology (itp/top)"), S("itp"), MD_IO_FLAG_SUPPLEMENT | MD_IO_FLAG_MM,
        .supplement = md_itp_system_supplement_from_file },
    [MD_IO_FORMAT_TPR] = { S("Gromacs Run Input (tpr)"), S("tpr"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_MM | MD_IO_FLAG_TOPOLOGY,
        .system_init = md_tpr_system_init_from_file },
    [MD_IO_FORMAT_NDX] = { S("Gromacs Index Groups (ndx)"), S("ndx"), MD_IO_FLAG_SUPPLEMENT | MD_IO_FLAG_MM,
        .supplement = ndx_supplement },
    [MD_IO_FORMAT_EDR] = { S("Gromacs Energy (edr)"), S("edr"), MD_IO_FLAG_SUPPLEMENT | MD_IO_FLAG_TEMPORAL | MD_IO_FLAG_MM,
        .supplement = md_edr_system_supplement_from_file },
    [MD_IO_FORMAT_XVG] = { S("xmgrace columns (xvg)"), S("xvg"), MD_IO_FLAG_SUPPLEMENT | MD_IO_FLAG_TEMPORAL,
        .supplement = md_xvg_system_supplement_from_file },
    [MD_IO_FORMAT_CSV] = { S("Comma separated values (csv)"), S("csv"), MD_IO_FLAG_SUPPLEMENT | MD_IO_FLAG_TEMPORAL,
        .supplement = md_csv_system_supplement_from_file },
    [MD_IO_FORMAT_ASE_TRAJ] = { S("ASE Trajectory (traj)"), S("traj"), MD_IO_FLAG_SYSTEM | MD_IO_FLAG_RUN | MD_IO_FLAG_MM,
        .system_init = md_ase_traj_system_init_from_file, .publish_run = md_ase_traj_system_publish_run },
    // A topology: a system once its coordinates are there (another file of them, as state_init of its format reads), or
    // the bonds and types of a loaded system of its atoms
    [MD_IO_FORMAT_AMBER_PRMTOP] = { S("Amber Topology (prmtop/parm7)"), S("prmtop"),
        MD_IO_FLAG_SYSTEM | MD_IO_FLAG_MM | MD_IO_FLAG_TOPOLOGY | MD_IO_FLAG_NEEDS_COORDINATES,
        .system_init = md_amber_prmtop_system_init_from_file, .supplement = md_amber_prmtop_system_supplement_from_file,
        .atom_count = md_amber_prmtop_file_atom_count },
    [MD_IO_FORMAT_AMBER_RESTART] = { S("Amber Restart (rst7/inpcrd/ncrst)"), S("rst7"), MD_IO_FLAG_COORDINATES | MD_IO_FLAG_MM,
        .state_init = amber_restart_state_init, .atom_count = amber_restart_atom_count },
    [MD_IO_FORMAT_AMBER_NETCDF] = { S("Amber NetCDF Trajectory (nc)"), S("nc"), MD_IO_FLAG_RUN | MD_IO_FLAG_MM,
        .publish_run = md_amber_netcdf_system_publish_run, .state_init = md_amber_netcdf_state_init_from_file,
        .atom_count = md_amber_netcdf_file_atom_count },
    [MD_IO_FORMAT_AMBER_ENERGY] = { S("Amber Energies (mdout/mden)"), S("mdout"), MD_IO_FLAG_SUPPLEMENT | MD_IO_FLAG_TEMPORAL | MD_IO_FLAG_MM,
        .supplement = md_amber_energy_system_supplement_from_file },
    [MD_IO_FORMAT_PSF] = { S("CHARMM/NAMD Topology (psf)"), S("psf"),
        MD_IO_FLAG_SYSTEM | MD_IO_FLAG_MM | MD_IO_FLAG_TOPOLOGY | MD_IO_FLAG_NEEDS_COORDINATES,
        .system_init = md_psf_system_init_from_file, .supplement = md_psf_system_supplement_from_file,
        .atom_count = md_psf_file_atom_count },
    [MD_IO_FORMAT_NAMD_COOR] = { S("NAMD Binary Coordinates (coor)"), S("coor"), MD_IO_FLAG_COORDINATES | MD_IO_FLAG_MM,
        .state_init = md_namdbin_state_init_from_file, .atom_count = md_namdbin_file_atom_count },
};

#undef S

const md_io_format_info_t* md_io_format_info(md_io_format_t format) {
    return &formats[(uint32_t)format < MD_IO_FORMAT_COUNT ? format : MD_IO_FORMAT_UNKNOWN];
}

bool md_io_format_available(md_io_format_t format) {
    const md_io_format_info_t* f = md_io_format_info(format);
    return f->system_init || f->supplement || f->publish_run || f->state_init;
}

md_io_format_t md_io_format_from_ext(md_str_t ext) {
    if (md_str_empty(ext)) return MD_IO_FORMAT_UNKNOWN;
    // Extensions other than the one a format is known by: a .top is an .itp with a [ molecules ] section (or an Amber
    // topology, which md_io_format_from_file tells from the contents)
    static const struct { const char* ext; md_io_format_t format; } aliases[] = {
        { "top", MD_IO_FORMAT_ITP },
        { "parm7", MD_IO_FORMAT_AMBER_PRMTOP },
        { "inpcrd", MD_IO_FORMAT_AMBER_RESTART }, { "restrt", MD_IO_FORMAT_AMBER_RESTART }, { "rst", MD_IO_FORMAT_AMBER_RESTART },
        { "ncrst", MD_IO_FORMAT_AMBER_RESTART },
        { "ncdf", MD_IO_FORMAT_AMBER_NETCDF }, { "netcdf", MD_IO_FORMAT_AMBER_NETCDF },
        { "mden", MD_IO_FORMAT_AMBER_ENERGY },
    };
    for (size_t i = 0; i < MD_COUNTOF(aliases); ++i) {
        if (md_str_eq_cstr_ignore_case(ext, aliases[i].ext)) return aliases[i].format;
    }
    for (int i = 1; i < MD_IO_FORMAT_COUNT; ++i) {
        if (md_str_eq_ignore_case(ext, formats[i].ext)) return (md_io_format_t)i;
    }
    return MD_IO_FORMAT_UNKNOWN;
}

// Where the extension does not tell, or is shared by several formats: Amber names its files without one (prmtop, inpcrd,
// mdout, mden), .top is a GROMACS topology or an Amber one, .out and .log are anything, .coor is binary or a PDB file
// and .nc is NetCDF of any convention. The contents decide.
static md_io_format_t format_from_content(md_str_t path, md_str_t ext, md_io_format_t from_ext) {
    switch (from_ext) {
    case MD_IO_FORMAT_ITP:
        return md_amber_prmtop_file_is_topology(path) ? MD_IO_FORMAT_AMBER_PRMTOP : MD_IO_FORMAT_ITP;
    case MD_IO_FORMAT_AMBER_NETCDF:
    case MD_IO_FORMAT_AMBER_RESTART:
        switch (md_amber_netcdf_file_kind(path)) {
        case MD_AMBER_NETCDF_TRAJECTORY: return MD_IO_FORMAT_AMBER_NETCDF;
        case MD_AMBER_NETCDF_RESTART:    return MD_IO_FORMAT_AMBER_RESTART;
        default:                         return from_ext;   // the reader says what is wrong with it
        }
    case MD_IO_FORMAT_NAMD_COOR:
        // binaryoutput off: NAMD writes its coordinates as a PDB file of the same name
        return md_namdbin_file_is_coor(path) ? MD_IO_FORMAT_NAMD_COOR : MD_IO_FORMAT_PDB;
#if MD_HDF5
    case MD_IO_FORMAT_VLX:
        // .h5 is not one format: a TREXIO file has a nucleus group and a VeloxChem one does not
        return md_trexio_file_is_trexio(path) ? MD_IO_FORMAT_TREXIO : MD_IO_FORMAT_VLX;
#endif
    case MD_IO_FORMAT_UNKNOWN:
        break;
    default:
        return from_ext;
    }
    if (!md_str_empty(ext) && !md_str_eq_ignore_case(ext, MD_STR_LIT("out")) && !md_str_eq_ignore_case(ext, MD_STR_LIT("log")) &&
        !md_str_eq_ignore_case(ext, MD_STR_LIT("txt")) && !md_str_eq_ignore_case(ext, MD_STR_LIT("dat"))) {
        return MD_IO_FORMAT_UNKNOWN;
    }
    if (md_amber_energy_file_is_mdout(path) || md_amber_energy_file_is_mden(path)) return MD_IO_FORMAT_AMBER_ENERGY;
    if (!md_str_empty(ext)) return MD_IO_FORMAT_UNKNOWN;
    if (md_amber_prmtop_file_is_topology(path)) return MD_IO_FORMAT_AMBER_PRMTOP;
    switch (md_amber_netcdf_file_kind(path)) {
    case MD_AMBER_NETCDF_TRAJECTORY: return MD_IO_FORMAT_AMBER_NETCDF;
    case MD_AMBER_NETCDF_RESTART:    return MD_IO_FORMAT_AMBER_RESTART;
    default: break;
    }
    if (md_amber_rst7_file_is_restart(path)) return MD_IO_FORMAT_AMBER_RESTART;
    if (md_psf_file_is_psf(path)) return MD_IO_FORMAT_PSF;
    return MD_IO_FORMAT_UNKNOWN;
}

md_io_format_t md_io_format_from_file(md_str_t path) {
    if (md_str_empty(path)) return MD_IO_FORMAT_UNKNOWN;
    const md_str_t ext = md_path_ext(path);
    const md_io_format_t format = format_from_content(path, ext, md_io_format_from_ext(ext));
    if (format != MD_IO_FORMAT_UNKNOWN) return format;
    // A Molden file is commonly called .molden.input or .mold, and .input is not an extension worth claiming
    if (md_molden_file_is_molden(path)) return MD_IO_FORMAT_MOLDEN;
    return MD_IO_FORMAT_UNKNOWN;
}
