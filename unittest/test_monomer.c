#include "utest.h"
#include "system_invariants.h"

#include <md_types.h>
#include <md_system.h>
#include <md_util.h>
#include <md_gro.h>
#include <md_pdb.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

// MONOMERS: the predefined table of standard monomers, the names of the residues it resolves, and the component types
// of md_util_system_infer_component_types which link to it.

UTEST(monomer, table) {
    for (int m = MD_MONOMER_NONE + 1; m < MD_MONOMER_COUNT; ++m) {
        const md_monomer_t monomer = (md_monomer_t)m;
        const str_t id = md_monomer_id(monomer);
        EXPECT_FALSE(str_empty(id));
        EXPECT_FALSE(str_empty(md_monomer_name(monomer)));
        EXPECT_NE(0, md_monomer_letter(monomer));
        EXPECT_NE(0u, md_monomer_color(monomer) >> 24);    // Opaque
        const md_component_kind_t kind = md_monomer_kind(monomer);
        EXPECT_TRUE(kind == MD_COMPONENT_KIND_AMINO_ACID || kind == MD_COMPONENT_KIND_NUCLEOTIDE);
        // Ids are unique: the id finds the monomer back, and the id of a residue is its own monomer
        EXPECT_EQ(monomer, md_monomer_from_id(id));
        EXPECT_EQ(monomer, md_monomer_from_comp_name(id, kind));
    }

    EXPECT_TRUE(str_empty(md_monomer_id(MD_MONOMER_NONE)));
    EXPECT_TRUE(str_empty(md_monomer_id(MD_MONOMER_COUNT)));
    EXPECT_EQ(0, md_monomer_letter(MD_MONOMER_NONE));
    EXPECT_EQ(0u, md_monomer_color(MD_MONOMER_NONE));
    EXPECT_EQ(MD_MONOMER_NONE, md_monomer_from_id(STR_LIT("")));
    EXPECT_EQ(MD_MONOMER_NONE, md_monomer_from_id(STR_LIT("HIE")));   // An alias is not an id
    EXPECT_EQ(MD_MONOMER_HIS, md_monomer_from_id(STR_LIT("his")));

    EXPECT_EQ('A', md_monomer_letter(MD_MONOMER_ALA));
    EXPECT_EQ('W', md_monomer_letter(MD_MONOMER_TRP));
    EXPECT_EQ('A', md_monomer_letter(MD_MONOMER_DA));
    EXPECT_EQ('U', md_monomer_letter(MD_MONOMER_U));
    EXPECT_EQ('X', md_monomer_letter(MD_MONOMER_UNK));
}

typedef struct name_case_t {
    const char* name;
    md_component_kind_t kind;
    md_monomer_t expected;
} name_case_t;

UTEST(monomer, from_comp_name) {
    const md_component_kind_t AA = MD_COMPONENT_KIND_AMINO_ACID;
    const md_component_kind_t NT = MD_COMPONENT_KIND_NUCLEOTIDE;
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
        {"XYZ", AA, MD_MONOMER_UNK}, {"ACE", AA, MD_MONOMER_UNK}, {"NME", MD_COMPONENT_KIND_OTHER, MD_MONOMER_UNK},
        // Case and padding
        {"hie", AA, MD_MONOMER_HIS}, {" ALA ", AA, MD_MONOMER_ALA},

        // Nucleotides: Amber termini and RNA, CHARMM, old PDB thymidine, modified
        {"DA5", NT, MD_MONOMER_DA}, {"DA3", NT, MD_MONOMER_DA}, {"DAN", NT, MD_MONOMER_DA},
        {"RA5", NT, MD_MONOMER_A}, {"RU3", NT, MD_MONOMER_U}, {"RG", NT, MD_MONOMER_G}, {"C5", NT, MD_MONOMER_C},
        {"ADE", NT, MD_MONOMER_A}, {"THY", NT, MD_MONOMER_DT}, {"URA", NT, MD_MONOMER_U}, {"T", NT, MD_MONOMER_DT},
        {"PSU", NT, MD_MONOMER_U}, {"5MC", NT, MD_MONOMER_C}, {"7MG", NT, MD_MONOMER_G},
        {"DN", NT, MD_MONOMER_N}, {"QQQ", NT, MD_MONOMER_N},

        // The kind decides: a name alone makes no monomer
        {"ALA", MD_COMPONENT_KIND_OTHER, MD_MONOMER_NONE},
        {"ALA", NT, MD_MONOMER_N},
        {"DA", AA, MD_MONOMER_UNK},
        {"HOH", MD_COMPONENT_KIND_WATER, MD_MONOMER_NONE},
        {"NA", MD_COMPONENT_KIND_ION, MD_MONOMER_NONE},
        {"POPC", MD_COMPONENT_KIND_LIPID, MD_MONOMER_NONE},
        {"LIG", MD_COMPONENT_KIND_OTHER, MD_MONOMER_NONE},
    };

    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        const md_monomer_t m = md_monomer_from_comp_name(str_from_cstr(cases[i].name), cases[i].kind);
        if (m != cases[i].expected) {
            printf("'%s' (kind %d): expected " STR_FMT ", got " STR_FMT "\n", cases[i].name, (int)cases[i].kind,
                   STR_ARG(md_monomer_id(cases[i].expected)), STR_ARG(md_monomer_id(m)));
        }
        EXPECT_EQ(cases[i].expected, m);
    }
}

// The type of every component with the given name, -1 if there is none or they differ
static md_component_type_idx_t type_of_name(const md_system_t* sys, const char* name) {
    md_component_type_idx_t t = -1;
    for (size_t i = 0; i < sys->component.count; ++i) {
        if (!str_eq_cstr(md_component_name(&sys->component, i), name)) continue;
        const md_component_type_idx_t ti = md_component_type_idx(&sys->component, i);
        if (t >= 0 && ti != t) return -1;
        t = ti;
    }
    return t;
}

UTEST(monomer, component_types_nucleosome) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro")));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    EXPECT_VALID_SYSTEM(&sys);

    const md_component_type_data_t* type = &sys.component.type;
    ASSERT_TRUE(sys.component.type_idx != NULL);
    EXPECT_GT(type->count, (size_t)20);
    EXPECT_LT(type->count, sys.component.count);

    // Each name and kind once
    for (size_t a = 0; a < type->count; ++a) {
        for (size_t b = a + 1; b < type->count; ++b) {
            EXPECT_FALSE(type->kind[a] == type->kind[b] && str_eq(LBL_TO_STR(type->name[a]), LBL_TO_STR(type->name[b])));
        }
    }

    const struct { const char* name; md_monomer_t monomer; } expected[] = {
        {"DA", MD_MONOMER_DA}, {"DT", MD_MONOMER_DT}, {"DG", MD_MONOMER_DG}, {"DC", MD_MONOMER_DC},
        {"ARG", MD_MONOMER_ARG}, {"LYS", MD_MONOMER_LYS}, {"ALA", MD_MONOMER_ALA},
    };
    for (size_t i = 0; i < ARRAY_SIZE(expected); ++i) {
        const md_component_type_idx_t t = type_of_name(&sys, expected[i].name);
        ASSERT_GE(t, 0);
        EXPECT_EQ(expected[i].monomer, md_component_type_monomer(type, t));
        EXPECT_EQ(md_monomer_color(expected[i].monomer), md_component_type_color(type, t));
    }

    // Inferring again keeps the types, and the colors they were given since
    const size_t count = type->count;
    const md_component_type_idx_t t_da = type_of_name(&sys, "DA");
    type->color[t_da] = 0xFF123456u;
    ASSERT_TRUE(md_util_system_infer_component_types(&sys));
    EXPECT_EQ(count, sys.component.type.count);
    EXPECT_EQ(t_da, type_of_name(&sys, "DA"));
    EXPECT_EQ(0xFF123456u, md_component_type_color(&sys.component.type, t_da));
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

UTEST(monomer, component_types_without_monomer) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(ligand_water_ion_pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    EXPECT_VALID_SYSTEM(&sys);

    const md_component_type_data_t* type = &sys.component.type;
    EXPECT_EQ((size_t)4, type->count);  // ALA, HOH (both waters), NA, LIG

    const md_component_type_idx_t ala = type_of_name(&sys, "ALA");
    const md_component_type_idx_t hoh = type_of_name(&sys, "HOH");
    const md_component_type_idx_t na  = type_of_name(&sys, "NA");
    const md_component_type_idx_t lig = type_of_name(&sys, "LIG");
    ASSERT_GE(ala, 0);
    ASSERT_GE(hoh, 0);
    ASSERT_GE(na, 0);
    ASSERT_GE(lig, 0);

    EXPECT_EQ(MD_MONOMER_ALA,  md_component_type_monomer(type, ala));
    EXPECT_EQ(MD_MONOMER_NONE, md_component_type_monomer(type, hoh));
    EXPECT_EQ(MD_MONOMER_NONE, md_component_type_monomer(type, na));
    EXPECT_EQ(MD_MONOMER_NONE, md_component_type_monomer(type, lig));
    EXPECT_EQ(MD_COMPONENT_KIND_WATER, md_component_type_kind(type, hoh));

    // A monatomic ion takes the color of its element, the others one of their name, opaque and the same every time
    EXPECT_EQ(md_atomic_number_cpk_color(MD_Z_Na), md_component_type_color(type, na));
    EXPECT_NE(0u, md_component_type_color(type, lig) >> 24);
    EXPECT_NE(0u, md_component_type_color(type, hoh) >> 24);
    EXPECT_NE(md_component_type_color(type, hoh), md_component_type_color(type, lig));

    md_system_t sys2 = { .alloc = alloc };
    md_system_state_t state2 = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys2, &state2, str_from_cstr(ligand_water_ion_pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys2, &state2, MD_UTIL_INFER_ALL));
    EXPECT_EQ(md_component_type_color(type, lig), md_component_type_color(&sys2.component.type, type_of_name(&sys2, "LIG")));

    md_vm_arena_destroy(alloc);
}
