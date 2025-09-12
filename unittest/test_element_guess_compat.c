#include "utest.h"
#include <md_util.h>
#include <md_molecule.h>
#include <md_pdb.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <string.h>

// Test backward compatibility with md_util_element_guess
UTEST(element_guess_compat, basic_inference) {
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(1));
    
    // Create a simple molecule structure
    md_molecule_t mol = {0};
    
    // Atom data
    const size_t atom_count = 5;
    mol.atom.count = atom_count;
    
    // Allocate arrays
    mol.atom.type = md_alloc(alloc, sizeof(md_label_t) * atom_count);
    mol.atom.resname = md_alloc(alloc, sizeof(md_label_t) * atom_count);
    mol.atom.element = md_alloc(alloc, sizeof(md_element_t) * atom_count);
    
    // Set up atom types and residue names
    // HOH water oxygen
    strncpy(mol.atom.type[0].buf, "O", sizeof(mol.atom.type[0].buf));
    mol.atom.type[0].len = 1;
    strncpy(mol.atom.resname[0].buf, "HOH", sizeof(mol.atom.resname[0].buf));
    mol.atom.resname[0].len = 3;
    mol.atom.element[0] = 0; // Start unknown
    
    // HOH water hydrogen
    strncpy(mol.atom.type[1].buf, "H1", sizeof(mol.atom.type[1].buf));
    mol.atom.type[1].len = 2;
    strncpy(mol.atom.resname[1].buf, "HOH", sizeof(mol.atom.resname[1].buf));
    mol.atom.resname[1].len = 3;
    mol.atom.element[1] = 0; // Start unknown
    
    // Alanine alpha carbon
    strncpy(mol.atom.type[2].buf, "CA", sizeof(mol.atom.type[2].buf));
    mol.atom.type[2].len = 2;
    strncpy(mol.atom.resname[2].buf, "ALA", sizeof(mol.atom.resname[2].buf));
    mol.atom.resname[2].len = 3;
    mol.atom.element[2] = 0; // Start unknown
    
    // Sodium ion
    strncpy(mol.atom.type[3].buf, "NA", sizeof(mol.atom.type[3].buf));
    mol.atom.type[3].len = 2;
    strncpy(mol.atom.resname[3].buf, "NA", sizeof(mol.atom.resname[3].buf));
    mol.atom.resname[3].len = 2;
    mol.atom.element[3] = 0; // Start unknown
    
    // Generic carbon
    strncpy(mol.atom.type[4].buf, "C1", sizeof(mol.atom.type[4].buf));
    mol.atom.type[4].len = 2;
    strncpy(mol.atom.resname[4].buf, "", sizeof(mol.atom.resname[4].buf));
    mol.atom.resname[4].len = 0;
    mol.atom.element[4] = 0; // Start unknown
    
    // Call the element guess function
    bool result = md_util_element_guess(mol.atom.element, atom_count, &mol);
    
    // Verify results
    EXPECT_TRUE(result);
    EXPECT_EQ(mol.atom.element[0], MD_Z_O);  // Water oxygen
    EXPECT_EQ(mol.atom.element[1], MD_Z_H);  // Water hydrogen
    EXPECT_EQ(mol.atom.element[2], MD_Z_C);  // Alanine alpha carbon (not calcium!)
    EXPECT_EQ(mol.atom.element[3], MD_Z_Na); // Sodium ion
    EXPECT_EQ(mol.atom.element[4], MD_Z_C);  // Generic carbon from C1
    
    md_vm_arena_destroy(alloc);
}

// Test element inference against all PDB files with explicit element symbols
UTEST(element_guess_compat, real_pdb_validation) {
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(16));
    
    // Test each PDB file individually
    
    // Test 1a64.pdb
    {
        str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/1a64.pdb");
        md_pdb_data_t pdb_data = {0};
        bool parse_result = md_pdb_data_parse_file(&pdb_data, path, alloc);
        EXPECT_TRUE(parse_result);
        
        if (parse_result) {
            md_molecule_t mol = {0};
            bool mol_result = md_pdb_molecule_init(&mol, &pdb_data, MD_PDB_OPTION_NONE, alloc);
            EXPECT_TRUE(mol_result);
            
            if (mol_result && mol.atom.count > 0) {
                size_t total_atoms = 0;
                size_t explicit_element_atoms = 0;
                size_t correct_inferences = 0;
                size_t cd_carbon_cases = 0;
                size_t cd_carbon_correct = 0;
                
                for (size_t i = 0; i < mol.atom.count && i < pdb_data.num_atom_coordinates; ++i) {
                    total_atoms++;
                    
                    const char* explicit_element = pdb_data.atom_coordinates[i].element;
                    if (explicit_element[0] != '\0' && explicit_element[0] != ' ') {
                        explicit_element_atoms++;
                        
                        str_t explicit_str = {explicit_element, strlen(explicit_element)};
                        explicit_str = str_trim(explicit_str);
                        md_element_t expected_element = md_util_element_lookup_ignore_case(explicit_str);
                        md_element_t inferred_element = mol.atom.element[i];
                        
                        if (expected_element != 0 && inferred_element == expected_element) {
                            correct_inferences++;
                        }
                        
                        // Special check for CD atoms in amino acids
                        str_t atom_name = LBL_TO_STR(mol.atom.type[i]);
                        str_t res_name = mol.atom.resname ? LBL_TO_STR(mol.atom.resname[i]) : (str_t){0};
                        if (str_eq_ignore_case(atom_name, STR_LIT("CD"))) {
                            cd_carbon_cases++;
                            if (res_name.len > 0 && md_util_resname_amino_acid(res_name)) {
                                if (inferred_element == MD_Z_C && expected_element == MD_Z_C) {
                                    cd_carbon_correct++;
                                }
                            }
                        }
                    }
                }
                
                // Validate accuracy is reasonable
                if (explicit_element_atoms > 0) {
                    double accuracy = (double)correct_inferences / explicit_element_atoms;
                    EXPECT_GT(accuracy, 0.85); // At least 85% accuracy for real files
                }
                
                // All CD atoms in amino acids should be correctly identified as Carbon
                if (cd_carbon_cases > 0) {
                    EXPECT_EQ(cd_carbon_correct, cd_carbon_cases);
                }
            }
            
            md_molecule_free(&mol, alloc);
        }
        md_pdb_data_free(&pdb_data, alloc);
    }
    
    // Test 1k4r.pdb
    {
        str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/1k4r.pdb");
        md_pdb_data_t pdb_data = {0};
        bool parse_result = md_pdb_data_parse_file(&pdb_data, path, alloc);
        EXPECT_TRUE(parse_result);
        
        if (parse_result) {
            md_molecule_t mol = {0};
            bool mol_result = md_pdb_molecule_init(&mol, &pdb_data, MD_PDB_OPTION_NONE, alloc);
            EXPECT_TRUE(mol_result);
            
            if (mol_result && mol.atom.count > 0) {
                size_t correct_inferences = 0;
                size_t explicit_element_atoms = 0;
                
                for (size_t i = 0; i < mol.atom.count && i < pdb_data.num_atom_coordinates; ++i) {
                    const char* explicit_element = pdb_data.atom_coordinates[i].element;
                    if (explicit_element[0] != '\0' && explicit_element[0] != ' ') {
                        explicit_element_atoms++;
                        
                        str_t explicit_str = {explicit_element, strlen(explicit_element)};
                        explicit_str = str_trim(explicit_str);
                        md_element_t expected_element = md_util_element_lookup_ignore_case(explicit_str);
                        md_element_t inferred_element = mol.atom.element[i];
                        
                        if (expected_element != 0 && inferred_element == expected_element) {
                            correct_inferences++;
                        }
                    }
                }
                
                if (explicit_element_atoms > 0) {
                    double accuracy = (double)correct_inferences / explicit_element_atoms;
                    EXPECT_GT(accuracy, 0.85);
                }
            }
            
            md_molecule_free(&mol, alloc);
        }
        md_pdb_data_free(&pdb_data, alloc);
    }
    
    // Test c60.pdb 
    {
        str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/c60.pdb");
        md_pdb_data_t pdb_data = {0};
        bool parse_result = md_pdb_data_parse_file(&pdb_data, path, alloc);
        EXPECT_TRUE(parse_result);
        
        if (parse_result) {
            md_molecule_t mol = {0};
            bool mol_result = md_pdb_molecule_init(&mol, &pdb_data, MD_PDB_OPTION_NONE, alloc);
            EXPECT_TRUE(mol_result);
            
            if (mol_result && mol.atom.count > 0) {
                // C60 should have all carbon atoms
                for (size_t i = 0; i < mol.atom.count; ++i) {
                    EXPECT_EQ(mol.atom.element[i], MD_Z_C);
                }
            }
            
            md_molecule_free(&mol, alloc);
        }
        md_pdb_data_free(&pdb_data, alloc);
    }
    
    // Test ciprofloxacin.pdb - small molecule with F, O atoms
    {
        str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/ciprofloxacin.pdb");
        md_pdb_data_t pdb_data = {0};
        bool parse_result = md_pdb_data_parse_file(&pdb_data, path, alloc);
        EXPECT_TRUE(parse_result);
        
        if (parse_result) {
            md_molecule_t mol = {0};
            bool mol_result = md_pdb_molecule_init(&mol, &pdb_data, MD_PDB_OPTION_NONE, alloc);
            EXPECT_TRUE(mol_result);
            
            if (mol_result && mol.atom.count > 0) {
                size_t correct_inferences = 0;
                size_t explicit_element_atoms = 0;
                
                for (size_t i = 0; i < mol.atom.count && i < pdb_data.num_atom_coordinates; ++i) {
                    const char* explicit_element = pdb_data.atom_coordinates[i].element;
                    if (explicit_element[0] != '\0' && explicit_element[0] != ' ') {
                        explicit_element_atoms++;
                        
                        str_t explicit_str = {explicit_element, strlen(explicit_element)};
                        explicit_str = str_trim(explicit_str);
                        md_element_t expected_element = md_util_element_lookup_ignore_case(explicit_str);
                        md_element_t inferred_element = mol.atom.element[i];
                        
                        if (expected_element != 0 && inferred_element == expected_element) {
                            correct_inferences++;
                        }
                    }
                }
                
                if (explicit_element_atoms > 0) {
                    double accuracy = (double)correct_inferences / explicit_element_atoms;
                    EXPECT_GT(accuracy, 0.85);
                }
            }
            
            md_molecule_free(&mol, alloc);
        }
        md_pdb_data_free(&pdb_data, alloc);
    }
    
    // Test tryptophan.pdb - amino acid
    {
        str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/tryptophan.pdb");
        md_pdb_data_t pdb_data = {0};
        bool parse_result = md_pdb_data_parse_file(&pdb_data, path, alloc);
        EXPECT_TRUE(parse_result);
        
        if (parse_result) {
            md_molecule_t mol = {0};
            bool mol_result = md_pdb_molecule_init(&mol, &pdb_data, MD_PDB_OPTION_NONE, alloc);
            EXPECT_TRUE(mol_result);
            
            if (mol_result && mol.atom.count > 0) {
                size_t correct_inferences = 0;
                size_t explicit_element_atoms = 0;
                
                for (size_t i = 0; i < mol.atom.count && i < pdb_data.num_atom_coordinates; ++i) {
                    const char* explicit_element = pdb_data.atom_coordinates[i].element;
                    if (explicit_element[0] != '\0' && explicit_element[0] != ' ') {
                        explicit_element_atoms++;
                        
                        str_t explicit_str = {explicit_element, strlen(explicit_element)};
                        explicit_str = str_trim(explicit_str);
                        md_element_t expected_element = md_util_element_lookup_ignore_case(explicit_str);
                        md_element_t inferred_element = mol.atom.element[i];
                        
                        if (expected_element != 0 && inferred_element == expected_element) {
                            correct_inferences++;
                        }
                    }
                }
                
                if (explicit_element_atoms > 0) {
                    double accuracy = (double)correct_inferences / explicit_element_atoms;
                    EXPECT_GT(accuracy, 0.85);
                }
            }
            
            md_molecule_free(&mol, alloc);
        }
        md_pdb_data_free(&pdb_data, alloc);
    }
    
    // Test 1ALA-560ns.pdb - trajectory file with alanine
    {
        str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/1ALA-560ns.pdb");
        md_pdb_data_t pdb_data = {0};
        bool parse_result = md_pdb_data_parse_file(&pdb_data, path, alloc);
        EXPECT_TRUE(parse_result);
        
        if (parse_result) {
            md_molecule_t mol = {0};
            bool mol_result = md_pdb_molecule_init(&mol, &pdb_data, MD_PDB_OPTION_NONE, alloc);
            EXPECT_TRUE(mol_result);
            
            if (mol_result && mol.atom.count > 0) {
                size_t correct_inferences = 0;
                size_t explicit_element_atoms = 0;
                
                for (size_t i = 0; i < mol.atom.count && i < pdb_data.num_atom_coordinates; ++i) {
                    const char* explicit_element = pdb_data.atom_coordinates[i].element;
                    if (explicit_element[0] != '\0' && explicit_element[0] != ' ') {
                        explicit_element_atoms++;
                        
                        str_t explicit_str = {explicit_element, strlen(explicit_element)};
                        explicit_str = str_trim(explicit_str);
                        md_element_t expected_element = md_util_element_lookup_ignore_case(explicit_str);
                        md_element_t inferred_element = mol.atom.element[i];
                        
                        if (expected_element != 0 && inferred_element == expected_element) {
                            correct_inferences++;
                        }
                    }
                }
                
                if (explicit_element_atoms > 0) {
                    double accuracy = (double)correct_inferences / explicit_element_atoms;
                    EXPECT_GT(accuracy, 0.85);
                }
            }
            
            md_molecule_free(&mol, alloc);
        }
        md_pdb_data_free(&pdb_data, alloc);
    }
    
    md_vm_arena_destroy(alloc);
}

// Comprehensive test that reports detailed statistics about element inference
UTEST(element_guess_compat, detailed_pdb_inference_stats) {
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(16));
    
    // Test 1a64.pdb and provide detailed statistics
    {
        str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/1a64.pdb");
        md_pdb_data_t pdb_data = {0};
        bool parse_result = md_pdb_data_parse_file(&pdb_data, path, alloc);
        
        if (parse_result) {
            md_molecule_t mol = {0};
            bool mol_result = md_pdb_molecule_init(&mol, &pdb_data, MD_PDB_OPTION_NONE, alloc);
            
            if (mol_result && mol.atom.count > 0) {
                size_t total_atoms = 0;
                size_t explicit_element_atoms = 0; 
                size_t correct_inferences = 0;
                size_t cd_cases = 0;
                size_t cd_correct = 0;
                size_t ca_cases = 0;
                size_t ca_correct = 0;
                
                for (size_t i = 0; i < mol.atom.count && i < pdb_data.num_atom_coordinates; ++i) {
                    total_atoms++;
                    
                    const char* explicit_element = pdb_data.atom_coordinates[i].element;
                    if (explicit_element[0] != '\0' && explicit_element[0] != ' ') {
                        explicit_element_atoms++;
                        
                        str_t explicit_str = {explicit_element, strlen(explicit_element)};
                        explicit_str = str_trim(explicit_str);
                        md_element_t expected_element = md_util_element_lookup_ignore_case(explicit_str);
                        md_element_t inferred_element = mol.atom.element[i];
                        
                        if (expected_element != 0 && inferred_element == expected_element) {
                            correct_inferences++;
                        }
                        
                        str_t atom_name = LBL_TO_STR(mol.atom.type[i]);
                        str_t res_name = mol.atom.resname ? LBL_TO_STR(mol.atom.resname[i]) : (str_t){0};
                        
                        // Track CD cases (the main issue we're fixing)
                        if (str_eq_ignore_case(atom_name, STR_LIT("CD"))) {
                            cd_cases++;
                            if (res_name.len > 0 && md_util_resname_amino_acid(res_name)) {
                                if (inferred_element == MD_Z_C && expected_element == MD_Z_C) {
                                    cd_correct++;
                                }
                            }
                        }
                        
                        // Track CA cases (should be Carbon in amino acids, not Calcium)
                        if (str_eq_ignore_case(atom_name, STR_LIT("CA"))) {
                            ca_cases++;
                            if (res_name.len > 0 && md_util_resname_amino_acid(res_name)) {
                                if (inferred_element == MD_Z_C && expected_element == MD_Z_C) {
                                    ca_correct++;
                                }
                            }
                        }
                    }
                }
                
                // Expected high accuracy for 1a64.pdb
                if (explicit_element_atoms > 0) {
                    double accuracy = (double)correct_inferences / explicit_element_atoms;
                    EXPECT_GT(accuracy, 0.85);
                }
                
                // All CD cases should be correctly identified as Carbon
                if (cd_cases > 0) {
                    EXPECT_EQ(cd_correct, cd_cases);
                }
                
                // All CA cases should be correctly identified as Carbon in amino acids
                if (ca_cases > 0) {
                    EXPECT_EQ(ca_correct, ca_cases);
                }
            }
            
            md_molecule_free(&mol, alloc);
        }
        md_pdb_data_free(&pdb_data, alloc);
    }
    
    md_vm_arena_destroy(alloc);
}