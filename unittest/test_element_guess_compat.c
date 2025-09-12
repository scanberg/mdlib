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
UTEST(element_guess_compat, all_pdb_validation) {
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(16));
    
    // List of PDB files to test
    const char* pdb_files[] = {
        MD_UNITTEST_DATA_DIR"/1a64.pdb",
        MD_UNITTEST_DATA_DIR"/1k4r.pdb", 
        MD_UNITTEST_DATA_DIR"/c60.pdb",
        MD_UNITTEST_DATA_DIR"/ciprofloxacin.pdb",
        MD_UNITTEST_DATA_DIR"/tryptophan.pdb",
        MD_UNITTEST_DATA_DIR"/1ALA-560ns.pdb",
        MD_UNITTEST_DATA_DIR"/dppc64.pdb"
    };
    const size_t num_files = sizeof(pdb_files) / sizeof(pdb_files[0]);
    
    // Aggregate statistics across all files
    size_t total_explicit_elements = 0;
    size_t total_correct_inferences = 0;
    size_t files_processed = 0;
    
    for (size_t file_idx = 0; file_idx < num_files; ++file_idx) {
        str_t path = {pdb_files[file_idx], strlen(pdb_files[file_idx])};
        md_pdb_data_t pdb_data = {0};
        bool parse_result = md_pdb_data_parse_file(&pdb_data, path, alloc);
        
        if (!parse_result) {
            // Some files might not exist, skip gracefully
            continue;
        }
        
        md_molecule_t mol = {0};
        bool mol_result = md_pdb_molecule_init(&mol, &pdb_data, MD_PDB_OPTION_NONE, alloc);
        
        if (mol_result && mol.atom.count > 0) {
            files_processed++;
            size_t file_explicit_elements = 0;
            size_t file_correct_inferences = 0;
            
            for (size_t i = 0; i < mol.atom.count && i < pdb_data.num_atom_coordinates; ++i) {
                const char* explicit_element = pdb_data.atom_coordinates[i].element;
                if (explicit_element[0] != '\0' && explicit_element[0] != ' ') {
                    file_explicit_elements++;
                    total_explicit_elements++;
                    
                    str_t explicit_str = {explicit_element, strlen(explicit_element)};
                    explicit_str = str_trim(explicit_str);
                    md_element_t expected_element = md_util_element_lookup_ignore_case(explicit_str);
                    md_element_t inferred_element = mol.atom.element[i];
                    
                    if (expected_element != 0 && inferred_element == expected_element) {
                        file_correct_inferences++;
                        total_correct_inferences++;
                    }
                }
            }
            
            md_molecule_free(&mol, alloc);
        }
        md_pdb_data_free(&pdb_data, alloc);
    }
    
    // Validate that we processed at least some files
    EXPECT_GT(files_processed, 0);
    
    // Calculate overall failure percentage across all matches
    if (total_explicit_elements > 0) {
        double overall_accuracy = (double)total_correct_inferences / total_explicit_elements;
        double failure_percentage = (1.0 - overall_accuracy) * 100.0;
        
        // Element inference should maintain > 85% accuracy across all test files
        EXPECT_GT(overall_accuracy, 0.85);
        EXPECT_LT(failure_percentage, 15.0);
    }
    
    md_vm_arena_destroy(alloc);
}