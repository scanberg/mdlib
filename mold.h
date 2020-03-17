/*
 *  MOLD - Molecule Drawcall Generator
 *  Copyright 2020 Robin Skånberg
 *  
 *  LICENSE see bottom of file
 *  
 *  Example use case:
 *  @TODO: Fill in
 *  
 */

#pragma once

#include <stdint.h>
#include <stdbool.h>

typedef uint8_t  mold_element_t;
typedef uint8_t  mold_secondary_structure_t;
typedef uint16_t mold_texture_t;
typedef uint32_t mold_atom_idx_t;

typedef struct mold_result_t {
    bool status;
    const char* error;
} mold_result_t;

/*
 *  MOLECULE
 *  Interface to create and modify molecule data
 */

typedef enum {
    MOLD_SECONDARY_STRUCTURE_UNKNOWN,
    MOLD_SECONDARY_STRUCTURE_COIL = MOLD_SECONDARY_STRUCTURE_UNKNOWN,
    MOLD_SECONDARY_STRUCTURE_HELIX,
    MOLD_SECONDARY_STRUCTURE_SHEET
} mold_secondary_structure;

typedef struct mold_mol_t {
    uint32_t id;
} mold_mol_t;

typedef struct mold_position_t {
    float x, y, z;
} mold_position_t;

typedef struct mold_bond_t {
    mold_atom_idx_t idx[2];
} mold_bond_t;

typedef struct mold_color_t {
    uint8_t r, g, b, a;
} mold_color_t;

typedef struct mold_backbone_segment_indices_t {
    mold_atom_idx_t ca_idx;
    mold_atom_idx_t n_idx;
    mold_atom_idx_t c_idx;
    mold_atom_idx_t o_idx;
} mold_backbone_segment_t;

typedef struct mold_backbone_sequence_t {
    uint32_t offset;
    uint32_t length;
} mold_backbone_sequence_t;

typedef struct {
    struct {
        uint32_t count;
        mold_position_t* position_data;
        mold_element_t*  element_data;
    } atom;

    struct {
        uint32_t count;
        mold_bond_t* data;        
    } bond;

    struct {
        struct {
            uint32_t count;
            mold_backbone_segment_indices_t*    index_data;
            mold_secondary_structure_t*         secondary_structure_data;
        } segment;
        struct {
            uint32_t count;
            mold_backbone_sequence_t* data;
        } sequence;
    } backbone;
} mold_mol_desc_t;

result_t mold_init_molecule(mold_mol_t* mol, const mold_mol_desc_t* desc);
result_t mold_free_molecule(mold_mol_t* mol);

result_t mold_set_atom_position_data(mold_mol_t molecule, const mold_position_t* ptr, uint32_t count, uint32_t offset);
result_t mold_set_atom_position_data(mold_mol_t molecule, const float* pos_x, const float* pos_y, const float* pos_z, uint32_t count, uint32_t offset);

result_t mold_set_atom_element_data(mold_mol_t molecule, const mold_element_t* ptr, uint32_t count, uint32_t offset);

result_t mold_set_covalent_bond_data(mold_mol_t molecule, const mold_bond_t* ptr, uint32_t count, uint32_t offset);
result_t mold_set_hydrogen_bond_data(mold_mol_t molecule, const mold_bond_t* ptr, uint32_t count, uint32_t offset);

result_t mold_set_backbone_segment_index_data(mold_mol_t molecule, const mold_backbone_segment_indices_t* ptr, uint32_t count, uint32_t offset);
result_t mold_set_backbone_secondary_structure_data(mold_mol_t molecule, const mold_secondary_structure_t* ptr, uint32_t count, uint32_t offset);
result_t mold_set_backbone_sequence_data(mold_mol_t molecule, const mold_backbone_sequence_t* ptr, uint32_t count, uint32_t offset);

/*
 *  REPRESENTATIONS
 *  Interface for creating visual representations
 */

typedef struct {
    uint32_t id;
} mold_rep_t;

typedef enum mold_rep_type {
    MOLD_REP_TYPE_SPACE_FILL,
    MOLD_REP_TYPE_LICORICE,
    MOLD_REP_TYPE_RIBBONS,
    MOLD_REP_TYPE_CARTOON,
    MOLD_REP_TYPE_SOLVENT_EXCLUDED_SURFACE
} mold_rep_type;

typedef union mold_rep_args_t {
    struct {
        float radius_scale;
    } space_fill;

    struct {
        float radius_scale;
    } licorice;

    struct {
        float width_scale;
        float thickness_scale;
    } ribbons;

    struct {
        float width_scale;
        float thickness_scale;
    } cartoon;

    struct {
        float probe_radius;
    } solvent_excluded_surface;
} mold_rep_args_t;

typedef struct {
    mold_rep_type_t type;

    struct {
        mold_color_t* color_data;
    } color;
} mold_rep_desc_t;

result_t mold_init_representation(mold_rep_t* rep, const mold_mol_desc_t* desc);
result_t mold_free_representation(mold_rep_t* rep);

result_t mold_set_representation(mold_rep_t rep, const mold_rep_desc_t* desc);
result_t mold_set_color_data(mold_rep_t rep, const mold_color_t* color_data, uint32_t count, uint32_t offset);

/*
 *  DRAW
 *  Interface for drawing representations of molecules
 */

typedef enum mold_gl_version {
    MOLD_GL_330_CORE,
    MOLD_GL_430_CORE
} mold_gl_version;

typedef uint32_t mold_gl_version_t;

typedef struct mold_gl_desc_t {
    mold_gl_version_t gl_version;
    
    struct {
        uint32_t width;
        uint32_t height;
        mold_texture_t depth;
        mold_texture_t color;
        mold_texture_t view_normal;
        mold_texture_t view_velocity;
        mold_texture_t atom_idx;
    } render_targets;
} mold_gl_desc_t;

mold_result_t mold_init_gl(const mold_gl_desc_t* desc);
mold_result_t mold_shutdown();

typedef struct mold_draw_desc_t {
    struct {
        float model_matrix[16];
        float view_matrix[16];
        float proj_matrix[16];
    } view_param;
} mold_draw_desc_t;

mold_result_t mold_draw(mold_mol_t mol, const mold_draw_desc_t* desc, const mold_rep_t representations[], uint32_t num_representations);
