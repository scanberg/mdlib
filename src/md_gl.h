#pragma once

//  MD GL - Molecule Dynamics GL
//  Copyright 2020 Robin Skånberg
// 
//  LICENSE: see license file

#include <stdint.h>
#include <stdbool.h>
#include "md_molecule.h"

#ifdef __cplusplus
extern "C" {
#endif

// ### MD GL ###
// Initialize the GL module and its resources

bool md_gl_initialize();
void md_gl_shutdown();

/*

In order to interface with your application, the shaders have a customizable output.
A function called 'write_fragment' is called in the fragment shader to finalize and write the fragments to the respective outputs.
This part is customizable by providing a snipped which controls what parameters to write and to which channels of the framebuffer.

If no custom_shader_snippet is supplied to md_gl_context_init, then a default_shader_snippet is used which is listed bellow:

layout(location = 0) out vec4 out_color;

void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint atom_index) {
   out_color  = color;
}
*/

// The shaders are customizable in their output as mentioned in the comment above.
// You could have multiple versions of shaders each compiled for a specific output.
typedef struct md_gl_shaders_t {
    uint64_t _opaque[4];
} md_gl_shaders_t;

bool md_gl_shaders_init(md_gl_shaders_t* shaders, const char* custom_shader_snippet);
bool md_gl_shaders_free(md_gl_shaders_t* shaders);


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

typedef struct md_gl_molecule_t {
    uint64_t _mem[12];
} md_gl_molecule_t;

bool md_gl_molecule_init(md_gl_molecule_t* gl_mol, const md_molecule_t* mol);
bool md_gl_molecule_free(md_gl_molecule_t* gl_mol);

// ### Set molecule data fields ###
bool md_gl_molecule_set_atom_position(md_gl_molecule_t* mol, uint32_t atom_offset, uint32_t atom_count, const float* x, const float* y, const float* z, uint32_t byte_stride);
bool md_gl_molecule_set_atom_velocity(md_gl_molecule_t* mol, uint32_t atom_offset, uint32_t atom_count, const float* x, const float* y, const float* z, uint32_t byte_stride);
bool md_gl_molecule_set_atom_radius  (md_gl_molecule_t* mol, uint32_t atom_offset, uint32_t atom_count, const float* radius, uint32_t byte_stride);
bool md_gl_molecule_set_atom_flags   (md_gl_molecule_t* mol, uint32_t atom_offset, uint32_t atom_count, const uint8_t* flags, uint32_t byte_stride);

// Call this function after setting new atomic positions to update velocities
// It will compute a new velocity as the difference between new and old atomic positions
// 
// pbc_ext is an optional parameter which corresponds to the extent of the periodic boundry conditions.
// If supplied, it will make sure to deperiodize the atomic positions before computing the velocity,
// which is crucial in order to resolve an accrucate velocity in a periodic domain.
bool md_gl_molecule_compute_velocity(md_gl_molecule_t* mol, const float pbc_ext[3]);

// Clear the velocity to zero.
bool md_gl_molecule_zero_velocity(md_gl_molecule_t* mol);

bool md_gl_molecule_set_covalent_bonds(md_gl_molecule_t* mol, uint32_t offset, uint32_t count, const md_bond_t* bonds, uint32_t byte_stride);
bool md_gl_molecule_set_backbone_secondary_structure(md_gl_molecule_t* mol, uint32_t offset, uint32_t count, const md_secondary_structure_t* secondary_structure, uint32_t byte_stride);

/*
 *  REPRESENTATIONS
 *  Interface for creating visual representations for molecules
 */
typedef struct md_gl_representation_t md_gl_representation_t;

struct md_gl_representation_t { // This is an opaque blob which matches the size of internal memory used by the structure
    uint64_t _mem[4];
};

typedef uint16_t md_gl_representation_type_t;
enum {
    MD_GL_REP_DEFAULT                    = 0,
    MD_GL_REP_SPACE_FILL                 = 1,
    MD_GL_REP_SOLVENT_EXCLUDED_SURFACE   = 2,
    MD_GL_REP_RIBBONS                    = 3,
    MD_GL_REP_CARTOON                    = 4,
    MD_GL_REP_LICORICE                   = 5
};

typedef union md_gl_representation_args_t {
    struct {
        float radius_scale;
    } space_fill;

    struct {
        float radius;
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
} md_gl_representation_args_t;

bool md_gl_representation_init(md_gl_representation_t* rep, const md_gl_molecule_t* mol);
bool md_gl_representation_free(md_gl_representation_t* rep);

bool md_gl_representation_set_type_and_args(md_gl_representation_t* rep, md_gl_representation_type_t type, md_gl_representation_args_t args);
bool md_gl_representation_set_color(md_gl_representation_t* rep, uint32_t offset, uint32_t count, const uint32_t* color, uint32_t byte_stride);

/*
 *  DRAW
 *  Interface for defining representations of molecules
 *
 */

typedef uint32_t md_gl_options_t;
enum {
    MD_GL_OPTION_NONE                      = 0,
    MD_GL_OPTION_RESIDUE_OCCLUSION_CULLING = 1,
};

typedef struct md_gl_draw_args_t {
    md_gl_shaders_t* shaders;

    struct {
        uint32_t count;
        const md_gl_representation_t** data;
        const float** model_matrix;
    } representation;

    struct {
        // @NOTE: Matrices are column major float[4][4]
        const float* view_matrix;
        const float* projection_matrix;

        // [Optional] These fields are only required if md_gl_rendertarget_t is used and a texture_view_velocity is present
        const float* prev_view_matrix; 
        const float* prev_projection_matrix;
    } view_transform;

    // If atom_mask is non-zero, the value is used as a mask and is ANDed with the atom flag of the molecule,
    // if the result after the operation is non-zero, the atom will be drawn.
    // Some representations (such as ribbons, cartoon) use spline segments derived from CA atoms, hide the CA atom => hide the segment
    uint32_t atom_mask;

    md_gl_options_t options;
} md_gl_draw_args_t;

bool md_gl_draw(const md_gl_draw_args_t* args);

#ifdef __cplusplus
}
#endif
