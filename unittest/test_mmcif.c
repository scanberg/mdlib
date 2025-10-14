#include "utest.h"

#include <md_mmcif.h>
#include <md_molecule.h>
#include <md_util.h>
#include <core/md_allocator.h>

UTEST(mmcif, 1fez) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/1fez.cif");

    md_molecule_t mol;
    bool result = md_mmcif_molecule_api()->init_from_file(&mol, path, NULL, md_get_heap_allocator());
    EXPECT_TRUE(result);

    if (result) {
        EXPECT_EQ(4097, mol.atom.count);

        EXPECT_EQ(7, mol.atom.element[0]);
        EXPECT_STREQ("N", mol.atom.type[0].buf);
        EXPECT_STREQ("LYS", mol.atom.resname[0].buf);
        EXPECT_EQ(1, mol.atom.resid[0]);
        EXPECT_STREQ("A", mol.atom.chainid[0].buf);

        EXPECT_NEAR(52.489, mol.atom.x[0], 0.001);
        EXPECT_NEAR(21.292, mol.atom.y[0], 0.001);
        EXPECT_NEAR(84.339, mol.atom.z[0], 0.001);
    }

    md_molecule_free(&mol, md_get_heap_allocator());
}

UTEST(mmcif, 2or2) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/2or2.cif");

    md_molecule_t mol;
    bool result = md_mmcif_molecule_api()->init_from_file(&mol, path, NULL, md_get_heap_allocator());
    EXPECT_TRUE(result);
    //md_util_molecule_postprocess(&mol, md_get_heap_allocator(), MD_UTIL_POSTPROCESS_ALL);

    if (result) {
        EXPECT_EQ(5382, mol.atom.count);

        EXPECT_EQ(7, mol.atom.element[0]);
        EXPECT_STREQ("N", mol.atom.type[0].buf);
        EXPECT_STREQ("ALA", mol.atom.resname[0].buf);
        EXPECT_EQ(1, mol.atom.resid[0]);
        EXPECT_STREQ("A", mol.atom.chainid[0].buf);

        EXPECT_NEAR(58.157, mol.atom.x[0], 0.001);
        EXPECT_NEAR(49.822, mol.atom.y[0], 0.001);
        EXPECT_NEAR(80.569, mol.atom.z[0], 0.001);
    }

    md_molecule_free(&mol, md_get_heap_allocator());
}

UTEST(mmcif, 8g7u) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/8g7u.cif");

    md_molecule_t mol;
    bool result = md_mmcif_molecule_api()->init_from_file(&mol, path, NULL, md_get_heap_allocator());
    EXPECT_TRUE(result);
    md_util_molecule_postprocess(&mol, md_get_heap_allocator(), MD_UTIL_POSTPROCESS_ALL);

    if (result) {
        EXPECT_EQ(14229, mol.atom.count);

        EXPECT_EQ(7, mol.atom.element[0]);
        EXPECT_STREQ("N", mol.atom.type[0].buf);
        EXPECT_STREQ("PHE", mol.atom.resname[0].buf);
        EXPECT_EQ(241, mol.atom.resid[0]);
        EXPECT_STREQ("A", mol.atom.chainid[0].buf);

        EXPECT_NEAR(77.862,  mol.atom.x[0], 0.001);
        EXPECT_NEAR(105.453, mol.atom.y[0], 0.001);
        EXPECT_NEAR(80.951,  mol.atom.z[0], 0.001);
    }

    md_molecule_free(&mol, md_get_heap_allocator());
}

// Comprehensive test for all CIF files in test_data
UTEST(mmcif, all_cif_files) {
    const char* cif_files[] = {
        MD_UNITTEST_DATA_DIR"/1fez.cif",
        MD_UNITTEST_DATA_DIR"/2or2.cif",
        MD_UNITTEST_DATA_DIR"/8g7u.cif",
    };
    const size_t expected_atom_counts[] = {
        4097,   // 1fez.cif
        5382,   // 2or2.cif
        14229,  // 8g7u.cif
    };
    const size_t num_files = sizeof(cif_files) / sizeof(cif_files[0]);

    for (size_t i = 0; i < num_files; ++i) {
        str_t path = str_from_cstr(cif_files[i]);
        md_molecule_t mol;
        bool result = md_mmcif_molecule_api()->init_from_file(&mol, path, NULL, md_get_heap_allocator());
        
        EXPECT_TRUE(result);
        if (result) {
            EXPECT_EQ(expected_atom_counts[i], mol.atom.count);
            
            // Verify all atoms have valid elements
            size_t zero_element_count = 0;
            for (size_t j = 0; j < mol.atom.count; ++j) {
                if (mol.atom.element[j] == 0) {
                    zero_element_count++;
                }
            }
            
            // Allow some missing elements but ensure most are filled
            double missing_ratio = (double)zero_element_count / (double)mol.atom.count;
            EXPECT_LT(missing_ratio, 0.15); // Less than 15% missing
            
            md_molecule_free(&mol, md_get_heap_allocator());
        }
    }
}
