#include "utest.h"
#include <md_util.h>
#include <md_system.h>
#include <core/md_os.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <string.h>

// Test element inference against all PDB files with explicit element symbols
UTEST(element_guess_compat, all_pdb_validation) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp_scope);
    
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
        // The columns of the atom lines themselves: the element (77 to 78), the atom name (13 to 16), the residue name
        // (18 to 20). Not through a loader: one which reads the element column names its types by it.
        str_t text = load_textfile(path, arena);
        if (str_empty(text)) {
            // Some files might not exist, skip gracefully
            continue;
        }
        files_processed++;
        str_t rest = text, line;
        while (str_extract_line(&line, &rest)) {
            if (!(str_eq_cstr_n(line, "ATOM", 4) || str_eq_cstr_n(line, "HETATM", 6)) || line.len < 78) continue;
            str_t explicit_symbol = str_trim(str_substr(line, 76, 2));
            if (str_empty(explicit_symbol)) {
                continue;
            }
            total_explicit_elements++;

            str_t atom_name    = str_trim(str_substr(line, 12, 4));
            str_t atom_resname = str_trim(str_substr(line, 17, 3));

            md_element_t expected_element = md_atomic_number_from_symbol(explicit_symbol, true);
            md_element_t inferred_element = md_atomic_number_infer_from_label(atom_name, atom_resname, 0);

            if (expected_element != 0 && inferred_element == expected_element) {
                total_correct_inferences++;
            }
        }

        md_vm_arena_set_pos_back(arena, temp_pos);
    }
    
    // Validate that we processed at least some files
    EXPECT_GT(files_processed, 0);    
    EXPECT_EQ(total_correct_inferences, total_explicit_elements);
    
    md_temp_end(temp_scope);
}