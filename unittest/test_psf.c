#include "utest.h"
#include <string.h>
#include <stdio.h>

#include <md_psf.h>
#include <md_namdbin.h>
#include <md_system.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_str.h>

static bool psf_write_text(md_str_t path, const char* text) {
    md_file_t f = {0};
    if (!md_file_open(&f, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
    md_file_write(f, text, strlen(text));
    md_file_close(&f);
    return true;
}

// A water and a potassium of CHARMM, in two segments, the second residue numbered with an insertion code. No
// elements: they come from names and masses (the potassium is POT, with CHARMM's 39.102).
static const char* psf_text =
    "PSF EXT\n"
    "\n"
    "         1 !NTITLE\n"
    " REMARKS a water and an ion\n"
    "\n"
    "         4 !NATOM\n"
    "         1 WAT      1        TIP3     OH2      OT      -0.834000       15.9994           0\n"
    "         2 WAT      1        TIP3     H1       HT       0.417000        1.0080           0\n"
    "         3 WAT      1        TIP3     H2       HT       0.417000        1.0080           0\n"
    "         4 ION      1A       POT      POT      POT      1.000000       39.1020           0\n"
    "\n"
    "         3 !NBOND: bonds\n"
    "         1         2         1         3         2         3\n"
    "\n"
    "         0 !NTHETA: angles\n";

UTEST(psf, builds_a_system_without_coordinates) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const md_str_t path = MD_STR_LIT("md_unittest.psf");
    ASSERT_TRUE(psf_write_text(path, psf_text));
    EXPECT_TRUE(md_psf_file_is_psf(path));
    EXPECT_EQ(4u, md_psf_file_atom_count(path));

    md_system_t sys = { .alloc = arena };
    md_system_state_t state = { .alloc = arena };
    ASSERT_TRUE(md_psf_system_init_from_file(&sys, &state, path));
    EXPECT_EQ(4u, sys.atom.count);
    EXPECT_EQ(2u, sys.residue.count);
    EXPECT_EQ(3u, sys.bond.count);
    EXPECT_FALSE(md_system_state_has_coords(&state));
    EXPECT_EQ(8,  sys.atom.type.z[sys.atom.type_idx[0]]);
    EXPECT_EQ(1,  sys.atom.type.z[sys.atom.type_idx[1]]);
    EXPECT_EQ(19, sys.atom.type.z[sys.atom.type_idx[3]]);
    EXPECT_EQ(1, sys.residue.seq_id[1]);
    EXPECT_TRUE(md_str_eq_cstr(md_system_residue_name(&sys, 1), "POT"));
    const md_attribute_t* charge = md_attributes_find(&sys.attributes, MD_STR_LIT("atom/charge"));
    ASSERT_TRUE(charge != NULL);
    EXPECT_NEAR(-0.834f, ((const float*)charge->data)[0], 1e-5f);

    // Applied to the same atoms: the bonds of the file
    md_bond_data_clear(&sys.bond);
    md_bond_build_connectivity(&sys.bond, sys.atom.count, arena);
    ASSERT_TRUE(md_psf_system_supplement_from_file(&sys, path));
    EXPECT_EQ(3u, sys.bond.count);

    ASSERT_TRUE(psf_write_text(path, "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N\n"));
    EXPECT_FALSE(md_psf_file_is_psf(path));
    md_system_free(&sys);
    remove(path.ptr);
    md_arena_allocator_destroy(arena);
}

UTEST(namdbin, coordinates_in_either_byte_order) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const md_str_t path = MD_STR_LIT("md_unittest.coor");
    for (int swap = 0; swap < 2; ++swap) {
        md_file_t f = {0};
        ASSERT_TRUE(md_file_open(&f, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
        uint32_t n = 3;
        if (swap) n = md_bswap32(n);
        md_file_write(f, &n, 4);
        for (int i = 0; i < 9; ++i) {
            double v = 1.5 * i - 4.0;
            uint64_t u;
            md_memcpy(&u, &v, 8);
            if (swap) u = md_bswap64(u);
            md_file_write(f, &u, 8);
        }
        md_file_close(&f);
        EXPECT_TRUE(md_namdbin_file_is_coor(path));
        EXPECT_EQ(3u, md_namdbin_file_atom_count(path));
        md_system_state_t s = { .alloc = arena };
        ASSERT_TRUE(md_namdbin_state_init_from_file(&s, path, 3));
        EXPECT_EQ(-4.0f, s.xyz[0].x);
        EXPECT_EQ(8.0f, s.xyz[2].z);
        EXPECT_FALSE(md_namdbin_state_init_from_file(&s, path, 4));
    }
    // A PDB file of the same name is not one
    ASSERT_TRUE(psf_write_text(path, "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N\n"));
    EXPECT_FALSE(md_namdbin_file_is_coor(path));
    remove(path.ptr);
    md_arena_allocator_destroy(arena);
}

// A blank segment id, as VMD and others write for a system without segments: the columns of the record are fixed
// (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4 and I10,1X,A8,... in EXT), and read as words a blank field shifts every field
// after it - the residue number read as the segment, the residue name as the number.
UTEST(psf, blank_segment_ids) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const md_str_t path = MD_STR_LIT("md_unittest_blank_segid.psf");
    static const char* texts[2] = {
        "PSF\n"
        "\n"
        "       1 !NTITLE\n"
        " REMARKS no segments\n"
        "\n"
        "       2 !NATOM\n"
        "       1      7    ALA  N    NH3   -0.300000       14.0070           0\n"
        "       2      7    ALA  CA   CT1    0.210000       12.0110           0\n"
        "\n"
        "       1 !NBOND: bonds\n"
        "       1       2\n",
        "PSF EXT\n"
        "\n"
        "         1 !NTITLE\n"
        " REMARKS no segments\n"
        "\n"
        "         2 !NATOM\n"
        "         1          7        ALA      N        NH3     -0.300000       14.0070           0\n"
        "         2          7        ALA      CA       CT1      0.210000       12.0110           0\n"
        "\n"
        "         1 !NBOND: bonds\n"
        "         1         2\n",
    };
    for (int k = 0; k < 2; ++k) {
        ASSERT_TRUE(psf_write_text(path, texts[k]));
        md_system_t sys = { .alloc = arena };
        md_system_state_t state = { .alloc = arena };
        ASSERT_TRUE(md_psf_system_init_from_file(&sys, &state, path));
        ASSERT_EQ(2u, sys.atom.count);
        ASSERT_EQ(1u, sys.residue.count);
        EXPECT_TRUE(md_str_eq_cstr(md_system_residue_name(&sys, 0), "ALA"));
        EXPECT_EQ(7, sys.residue.seq_id[0]);
        EXPECT_TRUE(md_str_eq_cstr(md_atom_name(&sys.atom, 1), "CA"));
        const md_attribute_t* charge = md_attributes_find(&sys.attributes, MD_STR_LIT("atom/charge"));
        ASSERT_TRUE(charge != NULL);
        EXPECT_NEAR(0.21f, ((const float*)charge->data)[1], 1e-5f);
        EXPECT_EQ(1u, sys.bond.count);
        md_system_free(&sys);
    }
    remove(path.ptr);
    md_arena_allocator_destroy(arena);
}
