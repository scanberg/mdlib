#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_molecule_t;
struct md_allocator_i;
struct md_trajectory_i;
struct md_trajectory_loader_i;
struct md_molecule_loader_i;
struct md_mat4_t;

enum {
	MD_PDB_COORD_FLAG_HETATM        = 1,
    MD_PDB_COORD_FLAG_TERMINATOR    = 2,   // Terminator for chain
};

enum {
    MD_PDB_OPTION_NONE = 0,
    MD_PDB_OPTION_CONCAT_MODELS = 1,
};

typedef uint32_t md_pdb_options_t;

typedef enum md_pdb_helix_class_t {
    Helix_Unknown           = 0,
    Helix_RH_alpha          = 1,
    Helix_RH_omega          = 2,
    Helix_RH_pi             = 3,
    Helix_RH_gamma          = 4,
    Helix_RH_3_10           = 5,
    Helix_LH_alpha          = 6,
    Helix_LH_omega          = 7,
    Helix_LH_gamma          = 8,
    Helix_2_7_ribbon_helix  = 9,
    Helix_Polyproline       = 10
} md_pdb_helix_class_t;

typedef enum md_pdb_sheet_sense_t {
    Sheet_AntiParallel = -1,
    Sheet_First = 0,
    Sheet_Parallel = 1,
} md_pdb_sheet_sense_t;

typedef struct md_pdb_model_t {
    int32_t serial;
    int32_t beg_atom_serial;
    int32_t end_atom_serial;

    // These are just meta data used for extracting trajectory information
    int32_t beg_atom_index; // Inclusive
    int32_t end_atom_index; // Exclusive
    size_t byte_offset;    
} md_pdb_model_t;

typedef struct md_pdb_helix_t {
    int32_t serial_number;
    char    id[4];
    int32_t init_res_seq;
    int32_t end_res_seq;
    char    init_res_name[4];
    char    end_res_name[4];
    char    init_res_i_code;
    char    end_res_i_code;
    char    init_chain_id;
    char    end_chain_id;
    md_pdb_helix_class_t helix_class;
    char    comment[30];
    int32_t length;
} md_pdb_helix_t;

typedef struct md_pdb_sheet_t {
    char    id[4];
    int32_t strand;
    int32_t num_strands;
    md_pdb_sheet_sense_t sense;

    char    init_res_name[4];
    char    end_res_name[4];
    char    cur_res_name[4];
    char    prev_res_name[4];
    char    cur_atom_name[4];
    char    prev_atom_name[4];
    
    int32_t init_res_seq;
    int32_t end_res_seq;
    int32_t cur_res_seq;
    int32_t prev_res_seq;

    char    init_chain_id;
    char    end_chain_id;
    char    cur_chain_id;
    char    prev_chain_id;

    char    init_res_i_code;
    char    end_res_i_code;
    char    cur_res_i_code;
    char    prev_res_i_code;
} md_pdb_sheet_t;

typedef struct md_pdb_coordinate_t {
    int32_t atom_serial;
    char atom_name[4];
    char res_name[5];
    char alt_loc;
    char chain_id;
    char icode;
    int32_t res_seq;
    float x;
    float y;
    float z;
    float occupancy;
    float temp_factor;
    char element[4];
    char charge[4];
	uint32_t flags; // MD_PDB_COORD_FLAG_
} md_pdb_coordinate_t;

typedef struct md_pdb_connect_t {
    int32_t atom_serial[8];
} md_pdb_connect_t;

typedef struct md_pdb_cryst1_t {
    float a;
    float b;
    float c;
    float alpha;
    float beta;
    float gamma;
    char space_group[16];
    int32_t z;
} md_pdb_cryst1_t;

// This structure represents the Biomolecule Assembly information found in REMARK 350
// In order to fully reconstruct the molecular structure, atoms within certain chains should be duplicated and transformed
typedef struct md_pdb_assembly_t {
    int32_t id;
    char    apply_to_chains[32];
    // These offset and count refer into the separate transform array
    uint32_t transform_offset;
    uint32_t transform_count;
} md_pdb_assembly_t;

typedef struct md_pdb_data_t {
    size_t num_cryst1;
    md_pdb_cryst1_t* cryst1;
    size_t num_models;
    md_pdb_model_t* models;
    size_t num_atom_coordinates;
    md_pdb_coordinate_t* atom_coordinates;
    size_t num_helices;
    md_pdb_helix_t* helices;
    size_t num_sheets;
    md_pdb_sheet_t* sheets;
    size_t num_connections;
    md_pdb_connect_t* connections;
    size_t num_assemblies;
    md_pdb_assembly_t* assemblies;
    size_t num_transforms;
    struct mat4_t* transforms;
} md_pdb_data_t;

// RAW FUNCTIONS
// Parse a text-blob as PDB
// This will append to the supplied pdb_data so make sure the data is in a valid state or zero
bool md_pdb_data_parse_str(md_pdb_data_t* data, str_t str, struct md_allocator_i* alloc);
bool md_pdb_data_parse_file(md_pdb_data_t* data, str_t filename, struct md_allocator_i* alloc);
void md_pdb_data_free(md_pdb_data_t* data, struct md_allocator_i* alloc);

// MOLECULE
bool md_pdb_molecule_init(struct md_molecule_t* mol, const md_pdb_data_t* data, md_pdb_options_t options, struct md_allocator_i* alloc);
struct md_molecule_loader_i* md_pdb_molecule_api(void);

// TRAJECTORY
struct md_trajectory_i* md_pdb_trajectory_create(str_t filename, struct md_allocator_i* alloc, uint32_t flags);
void md_pdb_trajectory_free(struct md_trajectory_i* traj);

struct md_trajectory_loader_i* md_pdb_trajectory_loader(void);

#ifdef __cplusplus
}
#endif
