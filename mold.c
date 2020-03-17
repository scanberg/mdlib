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

typedef uint8_t  mold_renderer_t;
typedef uint8_t  mold_element_t;
typedef uint8_t  mold_sstruct_t;
typedef uint8_t  mold_rep_type_t;
typedef uint16_t mold_texture_t;
typedef uint32_t mold_atom_idx_t;

typedef struct mold_result_t {
    bool status;
    const char* error;
} mold_result_t;

/*
 *  LIBRARY
 *  
 *  Interface to initialize and shutdown the library as a whole
 *  Will use the supplied allocator if it is declared, otherwise it will use malloc/free
 *  It will try to minimize the amount of allocations required.
 */

typedef enum {
    mold_renderer_opengl_33,
    mold_renderer_opengl_43,
/*  mold_renderer_vulkan  */    // We wish!
/*  mold_renderer_metal   */    // We wish we didn't have to!
} mold_renderer;

// @TODO: We will see if we possibly can avoid allocation anything at all, thus not needing this
typedef struct mold_allocator_t {
    void* (*alloc)(size_t num_bytes);
    void  (*free)(void* ptr);
} mold_allocator_t;

typedef struct mold_desc_t {
    mold_allocator_t allocator;
    mold_renderer_t renderer;
    struct {
        mold_texture_t depth;
        mold_texture_t color;
        mold_texture_t view_normal;
        mold_texture_t view_velocity;
        mold_texture_t atom_idx;
    } render_targets;
} mold_desc_t;

mold_result_t mold_initialize(const mold_desc_t* desc);
mold_result_t mold_shutdown();

/*
 *  MOLECULE
 *  Interface to create and modify molecule data
 */

typedef enum {
    mold_sstruct_unknown,
    mold_sstruct_coil,
    mold_sstruct_helix,
    mold_sstruct_sheet
} mold_sstruct;

typedef struct {
    uint32_t id;
} mold_mol_t;

typedef struct {
    float x, y, z;
} mold_position_t;

typedef struct {
    mold_atom_idx_t idx[2];
} mold_bond_t;

typedef struct {
    uint8_t r, g, b, a;
} mold_color_t;

typedef struct {
    mold_atom_idx_t ca_idx;
    mold_atom_idx_t n_idx;
    mold_atom_idx_t c_idx;
    mold_atom_idx_t o_idx;
    mold_sstruct_t  ss_type;
} mold_bb_segment_t;

typedef struct {
    uint32_t offset;
    uint32_t length;
} mold_bb_sequence_t;

typedef struct {
    mold_bb_segment_t* segments;
    mold_bb_sequence_t* seequences;
    uint32_t num_segments;
    uint32_t num_sequences;
} mold_backbone_t;

typedef struct {
    struct {
        uint32_t count;
        mold_position_t* position_data;
        mold_element_t*  element_data;
    } atom;

    struct {
        struct {
            uint32_t count;
            mold_bond_t* data;
        } covalent;
        struct {
            uint32_t count;
            mold_bond_t* data;
        } hydrogen;
    }

    struct {
        struct {
            uint32_t count;
            mold_backbone_segment_t* data;
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

result_t mold_set_bb_segment_data(mold_mol_t molecule, const mold_backbone_segment_t* ptr, uint32_t count, uint32_t offset);
result_t mold_set_bb_sequence_data(mold_mol_t molecule, const mold_backbone_sequence_t* ptr, uint32_t count, uint32_t offset);

/*
 *  REPRESENTATIONS
 *  Interface for creating visual representations
 */

typedef struct {
    uint32_t id;
} mold_rep_t;

typedef enum {
    mold_rep_type_space_fill,
    mold_rep_type_licorice,
    mold_rep_type_ribbons,
    mold_rep_type_cartoon,
    mold_rep_type_solvent_excluded_surface
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
typedef struct {
    struct {
        float model_mat[16];
        float view_mat[16];
        float proj_mat[16];
        float uv_jitter[2];     /* Screen-space jitter in UV-texture coordinates  */
    } view_param;
} mold_draw_t;

result_t draw_molecule(mold_mol_t mol, mold_rep_t* representations, uint32_t num_representations);