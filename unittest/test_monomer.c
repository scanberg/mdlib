#include "utest.h"
#include "system_invariants.h"

#include <md_types.h>
#include <md_system.h>
#include <md_infer.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_element.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <md_io.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <md_xyz.h>
#include <md_mmcif.h>
#include <md_tpr.h>
#include <md_lammps.h>
#include <md_run.h>
#include <core/md_geom.h>
#include <core/md_sort.h>
#include <md_dssp.h>
#include <md_itp.h>
#include <core/md_writer.h>
#include <core/md_os.h>
#include <core/md_path.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>
#include <core/md_vec_math.h>
#include "rmsd.h"
#include "system_util.h"

// MONOMERS: the predefined table of standard monomers, the names of the residues it resolves, and the residue types
// of md_system_infer_residue_types which link to it.

UTEST(monomer, table) {
    for (int m = MD_MONOMER_NONE + 1; m < MD_MONOMER_COUNT; ++m) {
        const md_monomer_t monomer = (md_monomer_t)m;
        const md_str_t id = md_monomer_id(monomer);
        EXPECT_FALSE(md_str_empty(id));
        EXPECT_FALSE(md_str_empty(md_monomer_name(monomer)));
        EXPECT_NE(0, md_monomer_letter(monomer));
        EXPECT_NE(0u, md_monomer_color(monomer) >> 24);    // Opaque
        const md_residue_kind_t kind = md_monomer_kind(monomer);
        EXPECT_TRUE(kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE);
        // Ids are unique: the id finds the monomer back, and the id of a residue is its own monomer
        EXPECT_EQ(monomer, md_monomer_from_id(id));
        EXPECT_EQ(monomer, md_monomer_from_residue_name(id, kind));
    }

    EXPECT_TRUE(md_str_empty(md_monomer_id(MD_MONOMER_NONE)));
    EXPECT_TRUE(md_str_empty(md_monomer_id(MD_MONOMER_COUNT)));
    EXPECT_EQ(0, md_monomer_letter(MD_MONOMER_NONE));
    EXPECT_EQ(0u, md_monomer_color(MD_MONOMER_NONE));
    EXPECT_EQ(MD_MONOMER_NONE, md_monomer_from_id(MD_STR_LIT("")));
    EXPECT_EQ(MD_MONOMER_NONE, md_monomer_from_id(MD_STR_LIT("HIE")));   // An alias is not an id
    EXPECT_EQ(MD_MONOMER_HIS, md_monomer_from_id(MD_STR_LIT("his")));

    EXPECT_EQ('A', md_monomer_letter(MD_MONOMER_ALA));
    EXPECT_EQ('W', md_monomer_letter(MD_MONOMER_TRP));
    EXPECT_EQ('A', md_monomer_letter(MD_MONOMER_DA));
    EXPECT_EQ('U', md_monomer_letter(MD_MONOMER_U));
    EXPECT_EQ('X', md_monomer_letter(MD_MONOMER_UNK));
}

typedef struct name_case_t {
    const char* name;
    md_residue_kind_t kind;
    md_monomer_t expected;
} name_case_t;

UTEST(monomer, from_res_name) {
    const md_residue_kind_t AA = MD_RESIDUE_KIND_AMINO_ACID;
    const md_residue_kind_t NT = MD_RESIDUE_KIND_NUCLEOTIDE;
    const name_case_t cases[] = {
        // Histidine of Amber, CHARMM, GROMACS, GROMOS, Martini
        {"HIE", AA, MD_MONOMER_HIS}, {"HID", AA, MD_MONOMER_HIS}, {"HIP", AA, MD_MONOMER_HIS},
        {"HSD", AA, MD_MONOMER_HIS}, {"HSP", AA, MD_MONOMER_HIS}, {"HISH", AA, MD_MONOMER_HIS},
        {"HISA", AA, MD_MONOMER_HIS}, {"HIS1", AA, MD_MONOMER_HIS}, {"HIH", AA, MD_MONOMER_HIS},
        // Other protonation states and bonding
        {"CYX", AA, MD_MONOMER_CYS}, {"CYM", AA, MD_MONOMER_CYS}, {"CYS2", AA, MD_MONOMER_CYS},
        {"ASH", AA, MD_MONOMER_ASP}, {"ASPP", AA, MD_MONOMER_ASP}, {"GLH", AA, MD_MONOMER_GLU},
        {"LYN", AA, MD_MONOMER_LYS}, {"LSN", AA, MD_MONOMER_LYS}, {"ARG0", AA, MD_MONOMER_ARG},
        // Termini of Amber, which do not shadow names of their own
        {"NALA", AA, MD_MONOMER_ALA}, {"CALA", AA, MD_MONOMER_ALA}, {"NHIE", AA, MD_MONOMER_HIS},
        {"CYSH", AA, MD_MONOMER_CYS},
        // Modified residues
        {"MSE", AA, MD_MONOMER_MET}, {"SEP", AA, MD_MONOMER_SER}, {"TPO", AA, MD_MONOMER_THR}, {"PTR", AA, MD_MONOMER_TYR},
        // Unknown amino acids, and the caps of a chain whatever they were classified as
        {"XYZ", AA, MD_MONOMER_UNK}, {"ACE", AA, MD_MONOMER_UNK}, {"NME", MD_RESIDUE_KIND_OTHER, MD_MONOMER_UNK},
        // Case and padding
        {"hie", AA, MD_MONOMER_HIS}, {" ALA ", AA, MD_MONOMER_ALA},

        // Nucleotides: Amber termini and RNA, CHARMM, old PDB thymidine, modified
        {"DA5", NT, MD_MONOMER_DA}, {"DA3", NT, MD_MONOMER_DA}, {"DAN", NT, MD_MONOMER_DA},
        {"RA5", NT, MD_MONOMER_A}, {"RU3", NT, MD_MONOMER_U}, {"RG", NT, MD_MONOMER_G}, {"C5", NT, MD_MONOMER_C},
        {"ADE", NT, MD_MONOMER_A}, {"THY", NT, MD_MONOMER_DT}, {"URA", NT, MD_MONOMER_U}, {"T", NT, MD_MONOMER_DT},
        {"PSU", NT, MD_MONOMER_U}, {"5MC", NT, MD_MONOMER_C}, {"7MG", NT, MD_MONOMER_G},
        {"DN", NT, MD_MONOMER_N}, {"QQQ", NT, MD_MONOMER_N},

        // The kind decides: a name alone makes no monomer
        {"ALA", MD_RESIDUE_KIND_OTHER, MD_MONOMER_NONE},
        {"ALA", NT, MD_MONOMER_N},
        {"DA", AA, MD_MONOMER_UNK},
        {"HOH", MD_RESIDUE_KIND_WATER, MD_MONOMER_NONE},
        {"NA", MD_RESIDUE_KIND_ION, MD_MONOMER_NONE},
        {"POPC", MD_RESIDUE_KIND_LIPID, MD_MONOMER_NONE},
        {"LIG", MD_RESIDUE_KIND_OTHER, MD_MONOMER_NONE},
    };

    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        const md_monomer_t m = md_monomer_from_residue_name(md_str_from_cstr(cases[i].name), cases[i].kind);
        if (m != cases[i].expected) {
            printf("'%s' (kind %d): expected " MD_STR_FMT ", got " MD_STR_FMT "\n", cases[i].name, (int)cases[i].kind,
                   MD_STR_ARG(md_monomer_id(cases[i].expected)), MD_STR_ARG(md_monomer_id(m)));
        }
        EXPECT_EQ(cases[i].expected, m);
    }
}

// The type of every component with the given name, -1 if there is none or they differ
static md_residue_type_idx_t type_of_name(const md_system_t* sys, const char* name) {
    md_residue_type_idx_t t = -1;
    for (size_t i = 0; i < sys->residue.count; ++i) {
        if (!md_str_eq_cstr(md_system_residue_name(sys, i), name)) continue;
        const md_residue_type_idx_t ti = md_system_residue_type_idx(sys, i);
        if (t >= 0 && ti != t) return -1;
        t = ti;
    }
    return t;
}

UTEST(monomer, residue_types_nucleosome) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"), NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    EXPECT_VALID_SYSTEM(&sys);

    const md_residue_type_data_t* type = &sys.residue.type;
    ASSERT_TRUE(sys.residue.type_idx != NULL);
    EXPECT_GT(type->count, (size_t)20);
    EXPECT_LT(type->count, sys.residue.count);

    // Each name and kind once
    for (size_t a = 0; a < type->count; ++a) {
        for (size_t b = a + 1; b < type->count; ++b) {
            EXPECT_FALSE(type->kind[a] == type->kind[b] && md_str_eq(MD_LABEL_STR(type->name[a]), MD_LABEL_STR(type->name[b])));
        }
    }

    const struct { const char* name; md_monomer_t monomer; } expected[] = {
        {"DA", MD_MONOMER_DA}, {"DT", MD_MONOMER_DT}, {"DG", MD_MONOMER_DG}, {"DC", MD_MONOMER_DC},
        {"ARG", MD_MONOMER_ARG}, {"LYS", MD_MONOMER_LYS}, {"ALA", MD_MONOMER_ALA},
    };
    for (size_t i = 0; i < MD_COUNTOF(expected); ++i) {
        const md_residue_type_idx_t t = type_of_name(&sys, expected[i].name);
        ASSERT_GE(t, 0);
        EXPECT_EQ(expected[i].monomer, md_residue_type_monomer(type, t));
        EXPECT_EQ(md_monomer_color(expected[i].monomer), md_residue_type_color(type, t));
    }

    // Inferring again keeps the types, and the colors they were given since
    const size_t count = type->count;
    const md_residue_type_idx_t t_da = type_of_name(&sys, "DA");
    type->color[t_da] = 0xFF123456u;
    ASSERT_TRUE(md_system_infer_residue_types(&sys));
    EXPECT_EQ(count, sys.residue.type.count);
    EXPECT_EQ(t_da, type_of_name(&sys, "DA"));
    EXPECT_EQ(0xFF123456u, md_residue_type_color(&sys.residue.type, t_da));
    EXPECT_VALID_SYSTEM(&sys);

    md_vm_arena_destroy(alloc);
}

static const char* ligand_water_ion_pdb =
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
    "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
    "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
    "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
    "ATOM      5  CB  ALA A   1       1.988  -0.773  -1.199  1.00  0.00           C\n"
    "HETATM    6  O   HOH A 101      10.000  10.000  10.000  1.00  0.00           O\n"
    "HETATM    7  O   HOH A 102      13.000  10.000  10.000  1.00  0.00           O\n"
    "HETATM    8 NA    NA A 201      20.000  20.000  20.000  1.00  0.00          NA\n"
    "HETATM    9  C1  LIG A 301      30.000  30.000  30.000  1.00  0.00           C\n"
    "HETATM   10  O1  LIG A 301      31.200  30.000  30.000  1.00  0.00           O\n"
    "END\n";

UTEST(monomer, residue_types_without_monomer) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(ligand_water_ion_pdb), NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    EXPECT_VALID_SYSTEM(&sys);

    const md_residue_type_data_t* type = &sys.residue.type;
    EXPECT_EQ((size_t)4, type->count);  // ALA, HOH (both waters), NA, LIG

    const md_residue_type_idx_t ala = type_of_name(&sys, "ALA");
    const md_residue_type_idx_t hoh = type_of_name(&sys, "HOH");
    const md_residue_type_idx_t na    = type_of_name(&sys, "NA");
    const md_residue_type_idx_t lig = type_of_name(&sys, "LIG");
    ASSERT_GE(ala, 0);
    ASSERT_GE(hoh, 0);
    ASSERT_GE(na, 0);
    ASSERT_GE(lig, 0);

    EXPECT_EQ(MD_MONOMER_ALA,  md_residue_type_monomer(type, ala));
    EXPECT_EQ(MD_MONOMER_NONE, md_residue_type_monomer(type, hoh));
    EXPECT_EQ(MD_MONOMER_NONE, md_residue_type_monomer(type, na));
    EXPECT_EQ(MD_MONOMER_NONE, md_residue_type_monomer(type, lig));
    EXPECT_EQ(MD_RESIDUE_KIND_WATER, md_residue_type_kind(type, hoh));

    // A monatomic ion takes the color of its element, the others one of their name, opaque and the same every time
    EXPECT_EQ(md_atomic_number_cpk_color(MD_Z_Na), md_residue_type_color(type, na));
    EXPECT_NE(0u, md_residue_type_color(type, lig) >> 24);
    EXPECT_NE(0u, md_residue_type_color(type, hoh) >> 24);
    EXPECT_NE(md_residue_type_color(type, hoh), md_residue_type_color(type, lig));

    md_system_t sys2 = { .alloc = alloc };
    md_system_state_t state2 = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys2, &state2, md_str_from_cstr(ligand_water_ion_pdb), NULL));
    ASSERT_TRUE(md_system_infer(&sys2, &state2, MD_SYSTEM_INFER_ALL, NULL));
    EXPECT_EQ(md_residue_type_color(type, lig), md_residue_type_color(&sys2.residue.type, type_of_name(&sys2, "LIG")));

    md_vm_arena_destroy(alloc);
}

// A modified residue larger than any standard one: the palmitoylated cysteine of a CHARMM-GUI spike protein (6VSB),
// named CYS, 58 atoms with its hydrogens. Its backbone is found by its atom names whatever its size, as a nucleotide's is.
static const char palmitoyl_cys_pdb[] =
    "ATOM  19115  N   CYS A1236     120.160 136.980 189.740  0.00  0.00           N\n"
    "ATOM  19116  HN  CYS A1236     120.310 136.990 190.720  0.00  0.00           H\n"
    "ATOM  19117  CA  CYS A1236     121.240 136.510 188.940  0.00  0.00           C\n"
    "ATOM  19118  HA  CYS A1236     121.650 137.310 188.340  0.00  0.00           H\n"
    "ATOM  19119  C   CYS A1236     120.710 135.530 187.950  0.00  0.00           C\n"
    "ATOM  19120  O   CYS A1236     120.880 135.590 186.730  0.00  0.00           O\n"
    "ATOM  19121  CB  CYS A1236     122.260 135.850 189.880  0.00  0.00           C\n"
    "ATOM  19122  HB1 CYS A1236     121.720 135.070 190.470  0.00  0.00           H\n"
    "ATOM  19123  HB2 CYS A1236     123.050 135.500 189.190  0.00  0.00           H\n"
    "ATOM  19124  SG  CYS A1236     123.200 136.950 191.060  0.00  0.00           S\n"
    "ATOM  19125  C1  CYS A1236     124.530 135.950 191.790  0.00  0.00           C\n"
    "ATOM  19126  O1  CYS A1236     124.750 134.930 191.150  0.00  0.00           O\n"
    "ATOM  19127  C2  CYS A1236     125.130 136.410 193.130  0.00  0.00           C\n"
    "ATOM  19128  H2A CYS A1236     124.670 137.340 193.530  0.00  0.00           H\n"
    "ATOM  19129  H2B CYS A1236     126.240 136.490 193.130  0.00  0.00           H\n"
    "ATOM  19130  C3  CYS A1236     124.810 135.380 194.180  0.00  0.00           C\n"
    "ATOM  19131  H3A CYS A1236     125.220 134.400 193.880  0.00  0.00           H\n"
    "ATOM  19132  H3B CYS A1236     123.720 135.180 194.310  0.00  0.00           H\n"
    "ATOM  19133  C4  CYS A1236     125.540 135.780 195.480  0.00  0.00           C\n"
    "ATOM  19134  H4A CYS A1236     125.330 136.760 195.940  0.00  0.00           H\n"
    "ATOM  19135  H4B CYS A1236     126.620 135.670 195.230  0.00  0.00           H\n"
    "ATOM  19136  C5  CYS A1236     125.150 134.810 196.570  0.00  0.00           C\n"
    "ATOM  19137  H5A CYS A1236     125.320 133.760 196.230  0.00  0.00           H\n"
    "ATOM  19138  H5B CYS A1236     124.110 134.990 196.910  0.00  0.00           H\n"
    "ATOM  19139  C6  CYS A1236     126.030 135.050 197.800  0.00  0.00           C\n"
    "ATOM  19140  H6A CYS A1236     125.820 136.010 198.320  0.00  0.00           H\n"
    "ATOM  19141  H6B CYS A1236     127.100 134.980 197.490  0.00  0.00           H\n"
    "ATOM  19142  C7  CYS A1236     125.820 133.990 198.900  0.00  0.00           C\n"
    "ATOM  19143  H7A CYS A1236     126.090 132.980 198.520  0.00  0.00           H\n"
    "ATOM  19144  H7B CYS A1236     124.760 133.900 199.240  0.00  0.00           H\n"
    "ATOM  19145  C8  CYS A1236     126.660 134.290 200.100  0.00  0.00           C\n"
    "ATOM  19146  H8A CYS A1236     126.210 135.060 200.780  0.00  0.00           H\n"
    "ATOM  19147  H8B CYS A1236     127.610 134.760 199.770  0.00  0.00           H\n"
    "ATOM  19148  C9  CYS A1236     127.010 133.060 200.990  0.00  0.00           C\n"
    "ATOM  19149  H9A CYS A1236     127.370 132.230 200.350  0.00  0.00           H\n"
    "ATOM  19150  H9B CYS A1236     125.990 132.790 201.330  0.00  0.00           H\n"
    "ATOM  19151  C10 CYS A1236     128.130 133.330 202.020  0.00  0.00           C\n"
    "ATOM  19152 H10A CYS A1236     128.000 134.290 202.560  0.00  0.00           H\n"
    "ATOM  19153 H10B CYS A1236     129.160 133.250 201.630  0.00  0.00           H\n"
    "ATOM  19154  C11 CYS A1236     127.990 132.280 203.120  0.00  0.00           C\n"
    "ATOM  19155 H11A CYS A1236     128.250 131.250 202.800  0.00  0.00           H\n"
    "ATOM  19156 H11B CYS A1236     126.970 132.430 203.520  0.00  0.00           H\n"
    "ATOM  19157  C12 CYS A1236     129.020 132.650 204.310  0.00  0.00           C\n"
    "ATOM  19158 H12A CYS A1236     128.830 133.620 204.800  0.00  0.00           H\n"
    "ATOM  19159 H12B CYS A1236     129.970 132.610 203.740  0.00  0.00           H\n"
    "ATOM  19160  C13 CYS A1236     128.870 131.530 205.400  0.00  0.00           C\n"
    "ATOM  19161 H13A CYS A1236     128.760 130.550 204.900  0.00  0.00           H\n"
    "ATOM  19162 H13B CYS A1236     127.910 131.720 205.930  0.00  0.00           H\n"
    "ATOM  19163  C14 CYS A1236     130.040 131.660 206.440  0.00  0.00           C\n"
    "ATOM  19164 H14A CYS A1236     131.060 131.450 206.040  0.00  0.00           H\n"
    "ATOM  19165 H14B CYS A1236     129.820 130.980 207.280  0.00  0.00           H\n"
    "ATOM  19166  C15 CYS A1236     130.200 133.070 207.240  0.00  0.00           C\n"
    "ATOM  19167 H15A CYS A1236     129.130 133.390 207.220  0.00  0.00           H\n"
    "ATOM  19168 H15B CYS A1236     130.700 133.870 206.650  0.00  0.00           H\n"
    "ATOM  19169  C16 CYS A1236     130.740 133.110 208.620  0.00  0.00           C\n"
    "ATOM  19170 H16A CYS A1236     130.110 132.440 209.250  0.00  0.00           H\n"
    "ATOM  19171 H16B CYS A1236     131.740 132.620 208.620  0.00  0.00           H\n"
    "ATOM  19172 H16C CYS A1236     130.750 134.090 209.140  0.00  0.00           H\n"
    "END\n";

// A nucleotide by name is one of the standard names, padded or with the terminal suffix of Amber and GROMACS (DA5, A3),
// and nothing else: not a ligand of the extended codes of the chemical component dictionary (A1IK6, A1JB9), whose letters
// before the digit are a standard name, nor U5P, nor TIP3P water (T3P).
UTEST(monomer, nucleotide_resname) {
    const char* yes_rna[] = {"A", "U", " G ", "C5", "A3"};
    const char* yes_dna[] = {"DA", "DT", "DG5", "DC3", " DA "};
    const char* no[]      = {"A1IK6", "A1JB9", "U5P", "T3P", "T4P", "T5P", "C14", "A2M", "N1", "5DA", "XA", "DA53", "DA7", "3", "", "ALA", "ADE"};
    for (size_t i = 0; i < MD_COUNTOF(yes_rna); ++i) {
        EXPECT_TRUE(md_resname_rna(md_str_from_cstr(yes_rna[i])));
        EXPECT_TRUE(md_resname_nucleotide(md_str_from_cstr(yes_rna[i])));
    }
    for (size_t i = 0; i < MD_COUNTOF(yes_dna); ++i) {
        EXPECT_TRUE(md_resname_dna(md_str_from_cstr(yes_dna[i])));
        EXPECT_TRUE(md_resname_nucleotide(md_str_from_cstr(yes_dna[i])));
    }
    for (size_t i = 0; i < MD_COUNTOF(no); ++i) {
        EXPECT_FALSE(md_resname_nucleotide(md_str_from_cstr(no[i])));
    }
}

// The ligand A1JB9 of 9qzl, a covalent inhibitor of NonO, has no atoms named as the backbone of a nucleotide: its name
// alone could make it one, and must not
static const char a1jb9_cif[] =
    "data_A1JB9\n"
    "loop_\n"
    "_atom_site.group_PDB\n"
    "_atom_site.id\n"
    "_atom_site.type_symbol\n"
    "_atom_site.label_atom_id\n"
    "_atom_site.label_alt_id\n"
    "_atom_site.label_comp_id\n"
    "_atom_site.label_asym_id\n"
    "_atom_site.label_entity_id\n"
    "_atom_site.label_seq_id\n"
    "_atom_site.Cartn_x\n"
    "_atom_site.Cartn_y\n"
    "_atom_site.Cartn_z\n"
    "_atom_site.auth_asym_id\n"
    "_atom_site.auth_seq_id\n"
    "HETATM 1 C C21 . A1JB9 C 2 . -14.89088 34.98645 8.41739 C 401\n"
    "HETATM 2 C C22 . A1JB9 C 2 . -13.54556 35.10987 8.12241 C 401\n"
    "HETATM 3 C C25 . A1JB9 C 2 . -10.62193 32.96442 7.63592 C 401\n"
    "HETATM 4 C C02 . A1JB9 C 2 . -13.84536 33.61134 13.30844 C 401\n"
    "HETATM 5 C C04 . A1JB9 C 2 . -15.65973 34.74598 12.08074 C 401\n"
    "HETATM 6 C C05 . A1JB9 C 2 . -16.95876 34.26957 11.68238 C 401\n"
    "HETATM 7 C C06 . A1JB9 C 2 . -18.12849 34.80314 12.38640 C 401\n"
    "HETATM 8 C C26 . A1JB9 C 2 . -13.34836 32.71988 7.99413 C 401\n"
    "HETATM 9 C C27 . A1JB9 C 2 . -14.69327 32.58878 8.28672 C 401\n"
    "HETATM 10 C C30 . A1JB9 C 2 . -16.32890 32.04194 11.23282 C 401\n"
    "HETATM 11 N N03 . A1JB9 C 2 . -14.79462 33.60739 12.37260 C 401\n"
    "HETATM 12 C C20 . A1JB9 C 2 . -15.45832 33.72549 8.50526 C 401\n"
    "HETATM 13 C C23 . A1JB9 C 2 . -12.77472 33.97918 7.91633 C 401\n"
    "HETATM 14 C C31 . A1JB9 C 2 . -14.88078 32.41807 11.51333 C 401\n"
    "HETATM 15 C C32 . A1JB9 C 2 . -13.05443 34.83487 13.50164 C 401\n"
    "HETATM 16 N N07 . A1JB9 C 2 . -17.93222 35.07418 13.69030 C 401\n"
    "HETATM 17 N N18 . A1JB9 C 2 . -17.05650 33.23376 10.72339 C 401\n"
    "HETATM 18 O O01 . A1JB9 C 2 . -13.57555 32.59012 13.96340 C 401\n"
    "HETATM 19 O O17 . A1JB9 C 2 . -19.21791 34.95509 11.83548 C 401\n"
    "HETATM 20 O O24 . A1JB9 C 2 . -11.42294 34.14170 7.75883 C 401\n"
    "HETATM 21 O O28 . A1JB9 C 2 . -17.75356 34.83010 8.95684 C 401\n"
    "HETATM 22 O O29 . A1JB9 C 2 . -17.68353 32.41143 8.45239 C 401\n"
    "HETATM 23 S S19 . A1JB9 C 2 . -17.11222 33.56123 9.08767 C 401\n"
    "HETATM 24 C C1 . A1JB9 C 2 . -19.05794 35.27664 14.59092 C 401\n"
    "HETATM 25 C C3 . A1JB9 C 2 . -19.76980 37.42108 13.48350 C 401\n"
    "HETATM 26 C C4 . A1JB9 C 2 . -20.15111 38.74873 13.52344 C 401\n"
    "HETATM 27 C C5 . A1JB9 C 2 . -20.27179 39.39254 14.74358 C 401\n"
    "HETATM 28 C C6 . A1JB9 C 2 . -20.00894 38.70341 15.91567 C 401\n"
    "HETATM 29 C C2 . A1JB9 C 2 . -19.50314 36.71711 14.64668 C 401\n"
    "HETATM 30 C C8 . A1JB9 C 2 . -20.91030 41.38337 13.57517 C 401\n"
    "HETATM 31 C C7 . A1JB9 C 2 . -19.62866 37.37564 15.86020 C 401\n"
    "HETATM 32 O O1 . A1JB9 C 2 . -20.65529 40.70556 14.80583 C 401\n"
    "#\n";

UTEST(monomer, extended_ccd_ligand_is_not_nucleotide) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_mmcif_system_init_from_str(&sys, &state, (md_str_t){a1jb9_cif, sizeof(a1jb9_cif) - 1}, NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ(1u, md_system_residue_count(&sys));
    EXPECT_TRUE(md_str_eq_cstr(md_system_residue_name(&sys, 0), "A1JB9"));
    EXPECT_EQ((int)MD_RESIDUE_KIND_OTHER, (int)md_system_residue_kind(&sys, 0));
    md_vm_arena_destroy(alloc);
}

// TIP3P water named T3P is water, not the nucleotide T
UTEST(monomer, t3p_is_water) {
    static const char gro[] =
        "T3P\n"
        "    3\n"
        "    1T3P     OW    1   0.000   0.000   0.000\n"
        "    1T3P    HW1    2   0.096   0.000   0.000\n"
        "    1T3P    HW2    3  -0.024   0.093   0.000\n"
        "   3.00000   3.00000   3.00000\n";
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, (md_str_t){gro, sizeof(gro) - 1}, NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ(1u, md_system_residue_count(&sys));
    EXPECT_EQ((int)MD_RESIDUE_KIND_WATER, (int)md_system_residue_kind(&sys, 0));
    md_vm_arena_destroy(alloc);
}

UTEST(monomer, amino_acid_of_any_size) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(palmitoyl_cys_pdb), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ(1u, md_system_residue_count(&sys));
    ASSERT_EQ(58u, (unsigned)sys.atom.count);
    EXPECT_EQ((int)MD_RESIDUE_KIND_AMINO_ACID, (int)md_system_residue_kind(&sys, 0));
    size_t backbone = 0, side = 0;
    for (size_t i = 0; i < sys.atom.count; ++i) {
        const md_atom_flags_t f = md_atom_flags(&sys.atom, i);
        backbone += (f & MD_ATOM_FLAG_BACKBONE) != 0;
        side += (f & MD_ATOM_FLAG_SIDE_CHAIN) != 0;
        if (md_str_eq_cstr(md_atom_name(&sys.atom, i), "C16")) EXPECT_TRUE((f & MD_ATOM_FLAG_SIDE_CHAIN) != 0);
    }
    EXPECT_EQ(3u, backbone);        // N, CA, C
    EXPECT_GT(side, 40u);           // CB, SG and the palmitoyl, with their hydrogens
    md_vm_arena_destroy(alloc);
}

// The roles of the atoms of a residue stay within it: a side chain does not run on through the N of a proline into the
// backbone before it, nor through a coordination bond into a metal, and a base is found from C1' whatever atom is
// listed first among its bonds.
UTEST(monomer, atom_roles_stay_in_their_residue) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1a64.pdb"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        size_t num_pro = 0, num_with_side_chain = 0, num_aa = 0;
        for (size_t c = 0; c < sys.residue.count; ++c) {
            if (md_system_residue_kind(&sys, c) != MD_RESIDUE_KIND_AMINO_ACID) continue;
            num_aa += 1;
            const md_str_t resname = md_system_residue_name(&sys, c);
            const bool pro = md_str_eq_cstr(resname, "PRO");
            num_pro += pro;
            size_t side = 0;
            const md_urange_t r = md_system_residue_atom_range(&sys, c);
            for (uint32_t i = r.beg; i < r.end; ++i) {
                const md_atom_flags_t f = md_atom_flags(&sys.atom, i);
                EXPECT_FALSE((f & MD_ATOM_FLAG_SIDE_CHAIN) && (f & MD_ATOM_FLAG_BACKBONE));
                if (!(f & MD_ATOM_FLAG_SIDE_CHAIN)) continue;
                side += 1;
                const md_str_t name = md_atom_name(&sys.atom, i);
                if (pro) EXPECT_TRUE(name.len >= 2 && (name.ptr[0] == 'C' || name.ptr[0] == 'H') && strchr("BGD", name.ptr[1]));
            }
            num_with_side_chain += side > 0;
        }
        EXPECT_GT(num_pro, 0u);
        EXPECT_GT(num_with_side_chain, num_aa / 2);
    }
    {
        // The Mg of 1fez is coordinated by an aspartate
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        size_t num_mg = 0;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            if (md_atom_atomic_number(&sys.atom, i) != MD_Z_Mg) continue;
            num_mg += 1;
            EXPECT_EQ(0, (int)(md_atom_flags(&sys.atom, i) & (MD_ATOM_FLAG_SIDE_CHAIN | MD_ATOM_FLAG_BACKBONE)));
        }
        EXPECT_EQ(2u, num_mg);
    }
    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        size_t num_nuc = 0, num_with_base = 0;
        for (size_t c = 0; c < sys.residue.count; ++c) {
            if (md_system_residue_kind(&sys, c) != MD_RESIDUE_KIND_NUCLEOTIDE) continue;
            num_nuc += 1;
            size_t base = 0;
            const md_urange_t r = md_system_residue_atom_range(&sys, c);
            for (uint32_t i = r.beg; i < r.end; ++i) {
                const md_atom_flags_t f = md_atom_flags(&sys.atom, i);
                EXPECT_FALSE((f & MD_ATOM_FLAG_NUCLEOBASE) && (f & MD_ATOM_FLAG_BACKBONE));
                base += (f & MD_ATOM_FLAG_NUCLEOBASE) != 0;
            }
            num_with_base += base > 0;
        }
        EXPECT_GT(num_nuc, 0u);
        EXPECT_EQ(num_nuc, num_with_base);
        // Nor does a base reach another residue
        for (size_t i = 0; i < sys.atom.count; ++i) {
            const md_residue_idx_t c = md_system_residue_find_by_atom_idx(&sys, i);
            if (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_NUCLEOBASE) EXPECT_EQ((int)MD_RESIDUE_KIND_NUCLEOTIDE, (int)md_system_residue_kind(&sys, c));
        }
    }
    md_vm_arena_destroy(alloc);
}

// A residue of one atom of an element which forms monatomic ions (O, H, Na, ...) is an ion where its name says so, or
// where it is bonded to nothing but by coordination. Orthosilicic acid, Si(OH)4, and a sodium, each atom a residue
// "MOL" as a LAMMPS data file gives them: classified first without bonds (as the loader does), every O and H is taken
// for an ion; with the bonds inferred they are not, and the sodium still is.
UTEST(monomer, lone_atom_ion_only_without_covalent_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(MD_MEGABYTES(64));
    const float d = 1.63f / 1.7320508f;     // Si-O 1.63 along the diagonals of a cube
    const float h = 2.59f / 1.7320508f;     // O-H 0.96 further out
    const float dir[4][3] = { {1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1} };
    char text[2048];
    int len = 0, serial = 1;
    #define ATOM(name, elem, x, y, z) len += snprintf(text + len, sizeof(text) - len, "HETATM%5d %-4s MOL A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n", serial, name, serial, (double)(x), (double)(y), (double)(z), elem), serial += 1
    ATOM("SI", "SI", 0, 0, 0);
    for (int k = 0; k < 4; ++k) ATOM("O", " O", d * dir[k][0], d * dir[k][1], d * dir[k][2]);
    for (int k = 0; k < 4; ++k) ATOM("H", " H", h * dir[k][0], h * dir[k][1], h * dir[k][2]);
    ATOM("NA", "NA", 8, 0, 0);
    #undef ATOM
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, (md_str_t){ text, (size_t)len }, NULL));
    ASSERT_EQ(10u, sys.residue.count);
    ASSERT_EQ(0u, sys.bond.count);

    ASSERT_TRUE(md_system_infer_residue_flags(&sys, &state));
    for (size_t r = 1; r < 10; ++r) EXPECT_EQ(MD_RESIDUE_KIND_ION, md_system_residue_kind(&sys, r));

    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    EXPECT_EQ(8u, sys.bond.count);
    for (size_t r = 1; r < 9; ++r) EXPECT_EQ(MD_RESIDUE_KIND_OTHER, md_system_residue_kind(&sys, r));
    EXPECT_EQ(MD_RESIDUE_KIND_ION, md_system_residue_kind(&sys, 9));
    md_vm_arena_destroy(alloc);
}

// Inferring again classifies the residues of atoms anew, from what the system is now: nothing found before is kept.
// The free tryptophan of tryptophan-md.gro (atoms named C15, N4, ...) is an amino acid by its bonds and its name;
// renamed LIG it is a ligand which holds the atoms of one, linked to nothing: no amino acid, its flags gone.
UTEST(monomer, residue_classified_anew) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"), NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    size_t trp = SIZE_MAX;
    for (size_t r = 0; r < sys.residue.count && trp == SIZE_MAX; ++r) {
        if (md_str_eq_cstr(md_system_residue_name(&sys, r), "TRP")) trp = r;
    }
    ASSERT_NE(SIZE_MAX, trp);
    EXPECT_EQ(MD_RESIDUE_KIND_AMINO_ACID, md_system_residue_kind(&sys, trp));
    EXPECT_EQ((unsigned)(MD_RESIDUE_FLAG_RESOLVED | MD_RESIDUE_FLAG_TERMINAL_BEG | MD_RESIDUE_FLAG_TERMINAL_END),
              (unsigned)md_system_residue_flags(&sys, trp) & (MD_RESIDUE_FLAG_RESOLVED | MD_RESIDUE_FLAG_TERMINAL_BEG | MD_RESIDUE_FLAG_TERMINAL_END));

    sys.residue.name[trp] = md_label_make(MD_STR_LIT("LIG"));
    md_system_topology_changed(&sys);
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    EXPECT_EQ(MD_RESIDUE_KIND_OTHER, md_system_residue_kind(&sys, trp));
    EXPECT_EQ(0u, (unsigned)md_system_residue_flags(&sys, trp) & (MD_RESIDUE_FLAG_RESOLVED | MD_RESIDUE_FLAG_TERMINAL_BEG | MD_RESIDUE_FLAG_TERMINAL_END));
    // The waters are classified anew too, the same
    size_t num_water = 0;
    for (size_t r = 0; r < sys.residue.count; ++r) num_water += md_system_residue_kind(&sys, r) == MD_RESIDUE_KIND_WATER;
    EXPECT_EQ(sys.residue.count - 1, num_water);
    md_vm_arena_destroy(alloc);
}
