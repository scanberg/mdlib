#include "utest.h"
#include <md_script.h>
#include <md_allocator.h>
#include <md_molecule.h>
#include <core/str_util.h>
#include <core/common.h>

UTEST(script, test) {

    // Create molecule for evaulation
    uint32_t atom_count = 16;
    float x[] = {1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8};
    float y[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    float z[] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
    float r[] = {1,2,3,4,4,4,5,1,1,2,3,4,4,4,5,1};
    float m[] = {1,2,2,2,2,4,4,4,1,2,2,2,2,4,4,4};
    uint8_t e[] = {1,8,1,8,5,8,8,8,1,8,1,8,5,8,8,8};
    md_label n[] = {"H", "O", "H", "He", "C", "N", "CA", "O"};
    float b[] = {0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0};
    float o[] = {1,1,1,1,1,2,2,2,1,1,1,1,1,2,2,2};
    uint32_t r_idx[] = {0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,3};
    uint32_t c_idx[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    uint32_t res_count = 4;
    md_label r_name[] = {"SOL", "LYS", "PFT", "PFT"};
    uint32_t r_id[] = {1, 2, 3, 4};
    md_range r_range[] = {{0, 3}, {3, 8}, {8,12}, {12,16}};

    uint32_t chain_count = 1;
    md_label c_id[] = {"A"};
    md_range c_arange[] = {0,16};
    md_range c_rrange[] = {0,4};

    md_molecule mol = {
        .atom = {
            .count = atom_count,
            .x = x,
            .y = y,
            .z = z,
            .radius = r,
            .mass = m,
            .element = e,
            .name = n,
            .bfactor = b,
            .occupancy = o,
            .residue_idx = r_idx,
            .chain_idx = c_idx
        },
        .residue = {
            .count = res_count,
            .name = r_name,
            .id = r_id,
            .atom_range = r_range
        },
        .chain = {
            .count = chain_count,
            .id = c_id,
            .atom_range = c_arange,
            .residue_range = c_rrange
        }
    };

    const char script_file[] = MD_UNITTEST_DATA_DIR "/script.txt";
    str_t script_text = load_textfile(script_file, ARRAY_SIZE(script_file), default_allocator);

    struct md_script_ir* ir = md_script_compile(script_text.ptr, script_text.len, &mol, default_allocator);

    free_str(script_text, default_allocator);
}