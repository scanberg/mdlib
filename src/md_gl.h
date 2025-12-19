#pragma once

//  MD GL - Molecule Dynamics GL
//  Copyright 2020 Robin Skånberg
// 
//  LICENSE: see license file

/*
@NOTE(Robin): This is a clusterf*ck API and I'm too lazy to change it by now since it will at some point be superseded by something that supports next gen rendering APIs.
*/

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <md_types.h>

// Forward declarations
struct md_system_t;

#if 0
// Predefined Palettes
enum {
    MD_GL_PALETTE_CPK,
    MD_GL_PALETTE_SECONDARY_STRUCTURE,
    MD_GL_PALETTE_ORDINAL_DIVERGING_HSV,
};
#endif

#ifdef __cplusplus
extern "C" {
#endif

// ### MD GL ###
// Initialize the GL module and its resources

void md_gl_initialize(void);
void md_gl_shutdown(void);

/*

In order to interface with your application, the shaders have a customizable output.
A function called 'write_fragment' is called in the fragment shader to finalize and write the fragments to the respective outputs.
This part is customizable by providing a snipped which controls what parameters to write and to which channels of the framebuffer.

If no custom_shader_snippet is supplied to md_gl_context_init, then a default_shader_snippet is used which is listed below:

layout(location = 0) out vec4 out_color;

void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint element_index) {
   out_color  = color;
}
*/

// The shaders are customizable in their output as mentioned in the comment above.
// You could have multiple versions of shaders each compiled for a specific output.
typedef struct {
    uint32_t id;
} md_gl_shaders_t;

md_gl_shaders_t md_gl_shaders_create(str_t str);
void md_gl_shaders_destroy(md_gl_shaders_t handle);


// ### GL MOLECULE ###
//
// The gl molecule represents the state of a single molecule which can be drawn.
// The required fields depends on the visual representation applied to the molecule, but to cover all cases the following fields are used:
//  - atom
//      - x      : f32
//      - y      : f32
//      - z      : f32
//      - radius : f32
//      - flags  : u8
//  - residue
//      - atom_range : u32[2]
//  - bond
//      - bond : md_bond_t = u32[2]
//  - backbone
//      - backbone_atoms      : md_backbone_atoms_t = u32[4] (N, CA, C, O)
//      - secondary_structure : md_secondary_structure_t = u8[4]

typedef struct {
    uint32_t id;
} md_gl_mol_t;

// Secondary Structure Type
// This is different from md_secondary_structure_t to allow interpolated values
// It essentially controls the blend (barycentric) between the geometrical representations of coil, helix, and sheet
// coil is implied when both helix and sheet are zero
typedef struct {
    float helix;  // 1.0 = full helix
    float sheet;  // 1.0 = full sheet
} md_gl_secondary_structure_t;

static inline md_gl_secondary_structure_t md_gl_secondary_structure_convert(md_secondary_structure_t ss) {
    md_gl_secondary_structure_t gl_ss = {0};
    switch (ss) {
        case MD_SECONDARY_STRUCTURE_HELIX_310:
        case MD_SECONDARY_STRUCTURE_HELIX_ALPHA:
        case MD_SECONDARY_STRUCTURE_HELIX_PI:
            gl_ss.helix = 1.0f;
            gl_ss.sheet = 0.0f;
            break;
        case MD_SECONDARY_STRUCTURE_BETA_SHEET:
        case MD_SECONDARY_STRUCTURE_BETA_BRIDGE:
            gl_ss.helix = 0.0f;
            gl_ss.sheet = 1.0f;
            break;
        case MD_SECONDARY_STRUCTURE_COIL:
        case MD_SECONDARY_STRUCTURE_TURN:
        case MD_SECONDARY_STRUCTURE_BEND:
        case MD_SECONDARY_STRUCTURE_UNKNOWN:
        default:
            gl_ss.helix = 0.0f;
            gl_ss.sheet = 0.0f;
            break;
    }
    return gl_ss;
}


md_gl_mol_t md_gl_mol_create(const struct md_system_t* mol);
void md_gl_mol_destroy(md_gl_mol_t handle);

// ### Set molecule data fields ###
void md_gl_mol_set_index_base   (md_gl_mol_t handle, uint32_t atom_index_base, uint32_t bond_index_base);
void md_gl_mol_set_atom_position(md_gl_mol_t handle, uint32_t atom_offset, uint32_t atom_count, const float* x, const float* y, const float* z, uint32_t byte_stride);
void md_gl_mol_set_atom_velocity(md_gl_mol_t handle, uint32_t atom_offset, uint32_t atom_count, const float* x, const float* y, const float* z, uint32_t byte_stride);
void md_gl_mol_set_atom_radius  (md_gl_mol_t handle, uint32_t atom_offset, uint32_t atom_count, const float* radius,  uint32_t byte_stride);
void md_gl_mol_set_atom_flags   (md_gl_mol_t handle, uint32_t atom_offset, uint32_t atom_count, const uint8_t* flags, uint32_t byte_stride);

// This is a simpler version which assumes packed xyz data (identical to the internal representation and therefore faster to copy)
void md_gl_mol_set_atom_position_xyz(md_gl_mol_t handle, uint32_t atom_offset, uint32_t atom_count, const float* xyz);

// Call this function after setting new atomic positions to update velocities
// It will compute a new velocity as the difference between new and old atomic positions
// 
// pbc_ext is an optional parameter which corresponds to the extent of the periodic boundry conditions.
// If supplied, it will make sure to deperiodize the atomic positions before computing the velocity,
// which is crucial in order to resolve an accrucate velocity in a periodic domain.
void md_gl_mol_compute_velocity(md_gl_mol_t handle, const float pbc_ext[3]);

// Clear the velocity to zero.
void md_gl_mol_zero_velocity(md_gl_mol_t handle);

void md_gl_mol_set_bonds(md_gl_mol_t mol, uint32_t offset, uint32_t count, const struct md_atom_pair_t* bond_pairs, uint32_t byte_stride);
void md_gl_mol_set_backbone_secondary_structure(md_gl_mol_t mol, uint32_t offset, uint32_t count, const md_gl_secondary_structure_t* secondary_structure, uint32_t byte_stride);

#if 0
typedef struct {
    uint32_t id;
} md_gl_palette_t;

md_gl_palette_t md_gl_palette_create();
void md_gl_palette_destroy(md_gl_palette_t pal);

// The color palette 
void md_gl_palette_set_colors(md_gl_palette_t pal, uint32_t offset, uint32_t count, const uint32_t* colors);

#endif

/*
 *  REPRESENTATIONS
 *  Interface for creating visual representations for molecules
 *  A representations is just a set of colors for each atom in the molecule
 */

typedef struct {
    uint32_t id;
} md_gl_rep_t;

md_gl_rep_t md_gl_rep_create(md_gl_mol_t mol);
void        md_gl_rep_destroy(md_gl_rep_t rep);

void        md_gl_rep_set_color(md_gl_rep_t rep, uint32_t offset, uint32_t count, const uint32_t* color, uint32_t byte_stride);

/*
 *  DRAW
 *  Interface for defining representations of molecules
 *
 */

typedef enum {
    MD_GL_REP_SPACE_FILL,
    MD_GL_REP_LICORICE,
    MD_GL_REP_BALL_AND_STICK,
    MD_GL_REP_RIBBONS,
    MD_GL_REP_CARTOON,
} md_gl_rep_type_t;

typedef struct md_gl_draw_op_t {
    md_gl_rep_type_t type;

    union {
        struct {
            float radius_scale;
        } space_fill;

        struct {
            float radius;
        } licorice;

        struct {
            float ball_scale;
            float stick_radius;
        } ball_and_stick;

        struct {
            float width_scale;
            float thickness_scale;
        } ribbons;

        struct {
            float coil_scale;
            float sheet_scale;
            float helix_scale;
        } cartoon;

        struct {
            float probe_radius;
        } solvent_excluded_surface;
    } args;

    md_gl_rep_t rep;
    const float* model_matrix;              // Column major float[4][4]
} md_gl_draw_op_t;

typedef struct md_gl_draw_args_t {
    md_gl_shaders_t shaders;

    struct {
        size_t count;
        md_gl_draw_op_t* ops;
    } draw_operations;

    struct {
        // @NOTE: Matrices are column major float[4][4]
        const float* view_matrix;
        const float* proj_matrix;

        // [Optional] These fields are only used to compute screen space velocity
        const float* prev_view_matrix; 
        const float* prev_proj_matrix;
    } view_transform;

    uint32_t atom_mask;
} md_gl_draw_args_t;

bool md_gl_draw(const md_gl_draw_args_t* args);

#ifdef __cplusplus
}
#endif
