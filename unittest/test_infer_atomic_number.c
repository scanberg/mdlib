#include "utest.h"
#include <md_util.h>
#include <md_molecule.h>
#include <md_pdb.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <string.h>

// Test element inference against all PDB files with explicit element symbols
UTEST(element_guess_compat, all_pdb_validation) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    
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

    size_t temp_pos = md_vm_arena_get_pos(arena);
    
    for (size_t file_idx = 0; file_idx < num_files; ++file_idx) {

        str_t path = {pdb_files[file_idx], strlen(pdb_files[file_idx])};
        md_pdb_data_t pdb_data = {0};
        bool parse_result = md_pdb_data_parse_file(&pdb_data, path, arena);
        
        if (!parse_result) {
            // Some files might not exist, skip gracefully
            continue;
        }
        
        if (pdb_data.num_atom_coordinates > 0) {
            files_processed++;
            
            for (size_t i = 0; i < pdb_data.num_atom_coordinates; ++i) {
                str_t explicit_symbol = str_trim(str_from_cstr(pdb_data.atom_coordinates[i].element));
                if (str_empty(explicit_symbol)) {
                    continue;
                }
                total_explicit_elements++;

                str_t atom_name    = str_trim(str_from_cstr(pdb_data.atom_coordinates[i].atom_name));
                str_t atom_resname = str_trim(str_from_cstr(pdb_data.atom_coordinates[i].res_name));
                    
                md_element_t expected_element = md_util_element_lookup_ignore_case(explicit_symbol);
                md_element_t inferred_element = md_atom_infer_atomic_number(atom_name, atom_resname);
                    
                if (expected_element != 0 && inferred_element == expected_element) {
                    total_correct_inferences++;
                }
            }
        }

        md_vm_arena_set_pos_back(arena, temp_pos);
    }
    
    // Validate that we processed at least some files
    EXPECT_GT(files_processed, 0);    
    EXPECT_EQ(total_correct_inferences, total_explicit_elements);
    
    md_vm_arena_destroy(arena);
}