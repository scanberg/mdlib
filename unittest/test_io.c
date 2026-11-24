#include "utest.h"
#include <string.h>
#include <stdio.h>

#include <md_io.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_file.h>
#include <core/md_str.h>

// Every format has a name and the extension it is known by, which gives the format back; what its flags say it holds,
// it has a reader of where it is built in.
UTEST(io, format_table) {
    EXPECT_FALSE(md_io_format_available(MD_IO_FORMAT_UNKNOWN));
    EXPECT_EQ(md_io_format_info(MD_IO_FORMAT_UNKNOWN), md_io_format_info(MD_IO_FORMAT_COUNT));
    for (int i = 1; i < MD_IO_FORMAT_COUNT; ++i) {
        const md_io_format_t format = (md_io_format_t)i;
        const md_io_format_info_t* f = md_io_format_info(format);
        EXPECT_FALSE(md_str_empty(f->name));
        EXPECT_FALSE(md_str_empty(f->ext));
        EXPECT_EQ(format, md_io_format_from_ext(f->ext));
        if (!md_io_format_available(format)) continue;
        // Every format loads; data of a system's atoms is checked as well, and only it is
        EXPECT_TRUE(f->load != NULL);
        EXPECT_EQ((f->flags & MD_IO_FLAG_SUPPLEMENT) != 0, f->check != NULL);
        if (f->flags & MD_IO_FLAG_TEMPORAL)     EXPECT_TRUE((f->flags & MD_IO_FLAG_SUPPLEMENT) != 0);
        if (f->publish_run)                     EXPECT_TRUE((f->flags & MD_IO_FLAG_RUN) != 0);
    }
#if MD_HDF5
    EXPECT_TRUE(md_io_format_available(MD_IO_FORMAT_H5MD));
#else
    EXPECT_FALSE(md_io_format_available(MD_IO_FORMAT_VLX));
    EXPECT_FALSE(md_io_format_available(MD_IO_FORMAT_TREXIO));
    EXPECT_FALSE(md_io_format_available(MD_IO_FORMAT_H5MD));
#endif
}

UTEST(io, format_from_ext) {
    EXPECT_EQ(MD_IO_FORMAT_PDB, md_io_format_from_ext(MD_STR_LIT("PDB")));
    EXPECT_EQ(MD_IO_FORMAT_ITP, md_io_format_from_ext(MD_STR_LIT("top")));
    EXPECT_EQ(MD_IO_FORMAT_AMBER_PRMTOP, md_io_format_from_ext(MD_STR_LIT("parm7")));
    EXPECT_EQ(MD_IO_FORMAT_AMBER_RESTART, md_io_format_from_ext(MD_STR_LIT("inpcrd")));
    EXPECT_EQ(MD_IO_FORMAT_AMBER_RESTART, md_io_format_from_ext(MD_STR_LIT("ncrst")));
    EXPECT_EQ(MD_IO_FORMAT_AMBER_NETCDF, md_io_format_from_ext(MD_STR_LIT("ncdf")));
    EXPECT_EQ(MD_IO_FORMAT_AMBER_ENERGY, md_io_format_from_ext(MD_STR_LIT("mden")));
    EXPECT_EQ(MD_IO_FORMAT_UNKNOWN, md_io_format_from_ext(MD_STR_LIT("")));
    EXPECT_EQ(MD_IO_FORMAT_UNKNOWN, md_io_format_from_ext(MD_STR_LIT("docx")));
}

UTEST(io, format_from_file) {
    static const struct { const char* file; md_io_format_t format; } files[] = {
        { MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb",                MD_IO_FORMAT_PDB },
        { MD_UNITTEST_DATA_DIR "/centered.gro",                  MD_IO_FORMAT_GRO },
        { MD_UNITTEST_DATA_DIR "/1fez.cif",                      MD_IO_FORMAT_MMCIF },
        { MD_UNITTEST_DATA_DIR "/catalyst.xtc",                  MD_IO_FORMAT_XTC },
        { MD_UNITTEST_DATA_DIR "/tryptophan-md.trr",             MD_IO_FORMAT_TRR },
        { MD_UNITTEST_DATA_DIR "/C60-Ih.xyz",                    MD_IO_FORMAT_XYZ },
        { MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol",        MD_IO_FORMAT_XMOL },
        { MD_UNITTEST_DATA_DIR "/ch4.arc",                       MD_IO_FORMAT_ARC },
        { MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data",  MD_IO_FORMAT_LAMMPS_DATA },
        { MD_UNITTEST_DATA_DIR "/cubic_standardASCII.lammpstrj", MD_IO_FORMAT_LAMMPS_TRJ },
        { MD_UNITTEST_DATA_DIR "/ener.edr",                      MD_IO_FORMAT_EDR },
        { MD_UNITTEST_DATA_DIR "/energy.xvg",                    MD_IO_FORMAT_XVG },
        { MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden",      MD_IO_FORMAT_MOLDEN },
        { MD_UNITTEST_DATA_DIR "/h5md/peptide_tip3p.tpr",        MD_IO_FORMAT_TPR },
    };
    for (size_t i = 0; i < MD_COUNTOF(files); ++i) {
        EXPECT_EQ(files[i].format, md_io_format_from_file(md_str_from_cstr(files[i].file)));
    }
    EXPECT_EQ(MD_IO_FORMAT_UNKNOWN, md_io_format_from_file(MD_STR_LIT("")));
}

static bool io_write_text(md_str_t path, const char* text) {
    md_file_t f = {0};
    if (!md_file_open(&f, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
    md_file_write(f, text, strlen(text));
    md_file_close(&f);
    return true;
}

static const char* io_psf_text =
    "PSF EXT\n"
    "\n"
    "         1 !NTITLE\n"
    " REMARKS a water\n"
    "\n"
    "         3 !NATOM\n"
    "         1 WAT      1        TIP3     OH2      OT      -0.834000       15.9994           0\n"
    "         2 WAT      1        TIP3     H1       HT       0.417000        1.0080           0\n"
    "         3 WAT      1        TIP3     H2       HT       0.417000        1.0080           0\n"
    "\n"
    "         2 !NBOND: bonds\n"
    "         1         2         1         3\n"
    "\n";

// The contents tell what the name does not: a PSF file without an extension
UTEST(io, format_from_contents) {
    const md_str_t psf = MD_STR_LIT("md_unittest_io_topology");
    ASSERT_TRUE(io_write_text(psf, io_psf_text));
    EXPECT_EQ(MD_IO_FORMAT_PSF, md_io_format_from_file(psf));
    const md_io_format_info_t* f = md_io_format_info(MD_IO_FORMAT_PSF);
    EXPECT_EQ(3u, f->atom_count(psf));

    // The table's reader of it, as an application would call it
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = { .alloc = arena };
    md_system_state_t state = { .alloc = arena };
    EXPECT_EQ(MD_LOAD_OK, f->load(&sys, &state, psf, NULL));
    EXPECT_EQ(3u, sys.atom.count);
    EXPECT_EQ(2u, md_bond_list_count(&sys.given_bond));
    md_arena_allocator_destroy(arena);
    remove(psf.ptr);
}

// A trajectory alone, through the table: a system of its atoms by their count, at the first frame, which the run of
// the same file then joins
UTEST(io, trajectory_alone) {
    static const md_io_format_t formats[] = { MD_IO_FORMAT_XTC, MD_IO_FORMAT_TRR, MD_IO_FORMAT_LAMMPS_TRJ };
    static const char* files[] = {
        MD_UNITTEST_DATA_DIR "/catalyst.xtc",
        MD_UNITTEST_DATA_DIR "/tryptophan-md.trr",
        MD_UNITTEST_DATA_DIR "/cubic_standardASCII.lammpstrj",
    };
    for (size_t i = 0; i < MD_COUNTOF(formats); ++i) {
        const md_io_format_info_t* f = md_io_format_info(formats[i]);
        const md_str_t path = md_str_from_cstr(files[i]);
        md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
        md_system_t sys = { .alloc = arena };
        md_system_state_t state = { .alloc = arena };
        ASSERT_EQ(MD_LOAD_OK, f->load(&sys, &state, path, NULL));
        EXPECT_TRUE(sys.atom.count > 0);
        EXPECT_EQ(sys.atom.count, state.num_atoms);
        EXPECT_TRUE(md_system_state_has_coords(&state));
        // Every atom of a type, of no element, which has a size: drawn, and read through its type
        ASSERT_TRUE(sys.atom.type_idx != NULL);
        EXPECT_EQ(0, md_atom_atomic_number(&sys.atom, 0));
        EXPECT_TRUE(md_atom_radius(&sys.atom, sys.atom.count - 1) > 0.0f);
        if (f->atom_count) EXPECT_EQ(f->atom_count(path), sys.atom.count);
        EXPECT_TRUE(f->publish_run(&sys, path, MD_STR_LIT("run/a"), 0));
        // A system of its own: not onto atoms
        EXPECT_EQ(MD_LOAD_FAILED, f->load(&sys, &state, path, NULL));
        md_system_free(&sys);
        md_arena_allocator_destroy(arena);
    }
}
