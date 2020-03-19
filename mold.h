#pragma once

/*
 *  MOLD - Molecule Draw
 *  Copyright 2020 Robin Skånberg
 *  
 *  LICENSE see bottom of file
 *  
 *  Example use case:
 *  @TODO: Fill in
 *  
 */

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef uint32_t mold_error_code_t;
typedef enum mold_error_code{
    MOLD_SUCCESS = 0,
    MOLD_UNKNOWN_ERROR,
    MOLD_REQUIRED_ARGUMENT_IS_NULL,
    MOLD_OUT_OF_MOLECULE_SLOTS,
    MOLD_OUT_OF_REPRESENTATION_SLOTS,
    MOLD_GL_VERSION_NOT_SUPPORTED,
} mold_error_code;

typedef struct mold_result_t {
    mold_error_code_t error_code;
    const char* error_str;
} mold_result_t;

inline int mold_success(mold_result_t res) {
    return res.error_code == MOLD_SUCCESS;
}

// Range defining a continous range of indices from [beg, end[ where beg is included, but end is excluded
// Example range{0, 5} would represent [0,1,2,3,4]
typedef struct mold_range_t {
    uint32_t beg;
    uint32_t end;
} mold_range_t;

inline uint32_t mold_extent(mold_range_t range) {
    return range.end - range.beg;
}

typedef uint32_t mold_texture_t;
typedef struct mold_context_desc_t {
    struct {
        void* (*allocate)(size_t num_bytes);
        void  (*free)(void* mem);
    } allocator;

    // GL specifics
    struct {
        struct {
            uint32_t width;
            uint32_t height;
            mold_texture_t depth;
            mold_texture_t color;
            mold_texture_t view_normal;
            mold_texture_t view_velocity;
            mold_texture_t atom_idx;
        } render_target;
    } gl;
} mold_context_desc_t;

typedef struct mold_context_t mold_context_t;

// Initializes the context, memory is allocated using the supplied allocator in desc
// if allocator is not specified, malloc and free will be used.
mold_result_t mold_context_init(mold_context_t* ctx, const mold_context_desc_t* desc);
mold_result_t mold_context_free(mold_context_t* ctx);

/*
 *  MOLECULE
 *  Interface to create and modify molecule data
 */
typedef uint8_t mold_secondary_structure_t;
typedef enum {
    MOLD_SECONDARY_STRUCTURE_UNKNOWN = 0,
    MOLD_SECONDARY_STRUCTURE_COIL,
    MOLD_SECONDARY_STRUCTURE_HELIX,
    MOLD_SECONDARY_STRUCTURE_SHEET
} mold_secondary_structure;

typedef struct mold_mol_t {
    uint32_t id;
} mold_mol_t;

typedef uint32_t mold_atom_idx_t;

typedef struct mold_position_t {
    float x, y, z;
} mold_position_t;

typedef float mold_radius_t;

typedef struct mold_bond_t {
    mold_atom_idx_t idx[2];
} mold_bond_t;

typedef struct mold_color_t {
    uint8_t r, g, b, a;
} mold_color_t;

typedef struct mold_segment_indices_t {
    mold_atom_idx_t ca_idx;
    mold_atom_idx_t c_idx;
    mold_atom_idx_t o_idx;
} mold_segment_indices_t;

typedef struct {
    struct {
        uint32_t count;
        mold_position_t* position_data;
        mold_radius_t*   radius_data;
    } atom;

    struct {
        uint32_t count;
        mold_range_t* atom_range_data;
    } residue;

    struct {
        uint32_t count;
        mold_bond_t* data;        
    } covalent_bond;

    struct {
        struct {
            uint32_t count;
            mold_segment_indices_t*     index_data;
            mold_secondary_structure_t* secondary_structure_data;
        } segment;
        struct {
            uint32_t count;
            mold_range_t* segment_range_data;
        } sequence;
    } backbone;
} mold_mol_desc_t;

mold_result_t mold_generate_molecule(mold_context_t* ctx, mold_mol_t* mol);
mold_result_t mold_init_molecule(mold_context_t* ctx, mold_mol_t mol, const mold_mol_desc_t* desc);
mold_result_t mold_free_molecule(mold_context_t* ctx, mold_mol_t mol);

mold_result_t mold_set_atom_position_data(mold_context_t* ctx, const mold_position_t pos_xyz[], uint32_t count, uint32_t offset);
mold_result_t mold_set_atom_position_data_soa(mold_context_t* ctx, const float pos_x[], const float pos_y[], const float pos_z[], uint32_t count, uint32_t offset);
mold_result_t mold_set_atom_radius_data(mold_context_t* ctx, const float radius[], uint32_t count, uint32_t offset);

mold_result_t mold_set_covalent_bond_data(mold_context_t* ctx, const mold_bond_t* bond[], uint32_t count, uint32_t offset);

mold_result_t mold_set_backbone_segment_index_data(mold_context_t* ctx, const mold_segment_indices_t segment_indices[], uint32_t count, uint32_t offset);
mold_result_t mold_set_backbone_secondary_structure_data(mold_context_t* ctx, const mold_secondary_structure_t secondary_structure[], uint32_t count, uint32_t offset);
mold_result_t mold_set_backbone_sequence_data(mold_context_t* ctx, const mold_range_t segment_sequence[], uint32_t count, uint32_t offset);

/*
 *  REPRESENTATIONS
 *  Interface for creating visual representations for molecules
 */

typedef struct {
    uint32_t id;
} mold_rep_t;

typedef uint32_t mold_rep_type_t;
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
    mold_rep_args_t args;

    // The range of atoms which the representation applies to
    // If supplied as (0,0) then its 
    mold_range_t atom_range;

    // [Optional] Atom color data, if supplied, it is expected to have the size of the extent of the atom_range 
    mold_color_t* atom_color_data;
} mold_rep_desc_t;

mold_result_t mold_gen_representation(mold_context_t* ctx, mold_rep_t* rep);
mold_result_t mold_init_representation(mold_context_t* ctx, mold_rep_t rep, const mold_mol_desc_t* desc);
mold_result_t mold_free_representation(mold_context_t* ctx, mold_rep_t rep);

mold_result_t mold_set_representation(mold_context_t* ctx, mold_rep_t rep, const mold_rep_desc_t* desc);
mold_result_t mold_set_color_data(mold_context_t* ctx, mold_rep_t rep, const mold_color_t color_data[], uint32_t count, uint32_t offset);

/*
 *  DRAW
 *  Interface for drawing representations of molecules
 */

typedef uint32_t mold_draw_flags_t;
typedef enum mold_draw_flags {
    MOLD_RESIDUE_OCCLUSION_CULLING_BIT = 0x1,
} mold_draw_flags;

typedef struct mold_draw_desc_t {
    mold_draw_flags_t flags;

    mold_mol_t molecule;
    struct {
        uint32_t count;
        const mold_rep_t* data;
    } representation;

    struct {
        float model[4][4];
        float view[4][4];
        float proj[4][4];
    } matrix;
} mold_draw_desc_t;

mold_result_t mold_draw(const mold_draw_desc_t* desc);

#ifdef __cplusplus
}
#endif