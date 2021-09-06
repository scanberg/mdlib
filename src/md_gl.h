#ifndef _MD_GL_H_
#define _MD_GL_H_

//  MD GL - Molecule Dynamics GL
//  Copyright 2020 Robin Skånberg
// 
//  LICENSE: see license file
// 
//  Example use case:
//  @TODO: Fill in

#include <stdint.h>
#include <stdbool.h>
#include "md_molecule.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MD_SPLINE_MAX_SUBDIVISION_COUNT 8

// ### CONTEXT ###
//
// The context holds the shared resources such as buffers, shaders and textures needed to draw molecules
typedef struct md_gl_context md_gl_context;
struct md_gl_context {
    uint64_t _opaque[72];
};

bool md_gl_context_init(md_gl_context* ctx);
bool md_gl_context_free(md_gl_context* ctx);


// ### GL MOLECULE ###
//
// The gl molecule represents the state of a single molecule which can be drawn
// Its data can either be initialized from a md_molecule or the data fields can be set individually
// The required fields depends on the visual representation applied to the molecule, but to cover all cases the following fields are used:
//  - atom
//      - x      : f32
//      - y      : f32
//      - z      : f32
//      - radius : f32
//      - flags  : u8
//  - residue
//      - atom_range          : u32[2]
//      - backbone_atoms      : u32[4] (N, CA, C, O)
//      - secondary_structure : u32
//  - bond
//      - bond : u32[2]

typedef struct md_gl_molecule md_gl_molecule;
struct md_gl_molecule {
    uint64_t _mem[12];
};

typedef struct md_gl_molecule_desc md_gl_molecule_desc;
struct md_gl_molecule_desc {
    struct {
        uint32_t count;
        const float*   x;
        const float*   y;
        const float*   z;
        const float*   radius;
        const uint8_t* flags;
    } atom;

    struct {
        uint32_t count;
        const md_bond_t* atom_bond;
    } covalent_bond;

    struct {
        uint32_t count;
        const md_range_t* atom_range; // This is just for building AABBs for culling
    } residue;

    struct {
        uint32_t count;
        //const md_range_t*               atom_range;
        const md_backbone_atoms_t*      atoms;
        const md_secondary_structure_t* secondary_structure;
    } backbone;

    struct {
        uint32_t count;
        const md_range_t* backbone_range;
    } chain;
};

bool md_gl_molecule_init(md_gl_molecule* mol, const md_gl_molecule_desc* desc);
bool md_gl_molecule_free(md_gl_molecule* mol);

// ### Update by modifying the dynamic data fields ###
bool md_gl_molecule_set_atom_position(md_gl_molecule* mol, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t byte_stride);
bool md_gl_molecule_set_atom_radius(md_gl_molecule* mol, uint32_t offset, uint32_t count, const float* radius, uint32_t byte_stride);
bool md_gl_molecule_set_atom_flags(md_gl_molecule* mol, uint32_t offset, uint32_t count, const uint8_t* flags, uint32_t byte_stride);

bool md_gl_molecule_set_backbone_secondary_structure(md_gl_molecule* mol, uint32_t offset, uint32_t count, const md_secondary_structure_t* secondary_structure, uint32_t byte_stride);

// This is called to copy the atom position buffer to previous atom position
// Usually just before calling md_gl_molecule_set_atom_position() to set a new current position
bool md_gl_molecule_update_atom_previous_position(md_gl_molecule* mol);
/*
 *  REPRESENTATIONS
 *  Interface for creating visual representations for molecules
 */
typedef struct md_gl_representation md_gl_representation;

struct md_gl_representation { // This is an opaque blob which matches the size of internal memory used by the structure
    uint64_t _mem[4];
};

typedef uint16_t md_gl_representation_type;
enum {
    MD_GL_REP_DEFAULT                    = 0,
    MD_GL_REP_SPACE_FILL                 = 1,
    MD_GL_REP_SOLVENT_EXCLUDED_SURFACE   = 2,
    MD_GL_REP_RIBBONS                    = 3,
    MD_GL_REP_CARTOON                    = 4,
    MD_GL_REP_LICORICE                   = 5
};

typedef union md_gl_representation_args {
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
} md_gl_representation_args;

bool md_gl_representation_init(md_gl_representation* rep, const md_gl_molecule* mol);
bool md_gl_representation_free(md_gl_representation* rep);

bool md_gl_representation_set_type_and_args(md_gl_representation* rep, md_gl_representation_type type, md_gl_representation_args args);
bool md_gl_representation_set_color(md_gl_representation* rep, uint32_t offset, uint32_t count, const uint32_t* color, uint32_t byte_stride);

/*
 *  DRAW
 *  Interface for defining representations of molecules
 *
 */

typedef struct md_gl_rendertarget {
    uint32_t width;
    uint32_t height;
    uint32_t texture_depth;
    uint32_t texture_color;             // [Optional] (0 = ignore)
    uint32_t texture_atom_index;        // [Optional] (0 = ignore)
    uint32_t texture_view_normal;       // [Optional] (0 = ignore)
    uint32_t texture_view_velocity;     // [Optional] (0 = ignore)
} md_gl_rendertarget;

typedef uint32_t md_gl_options;
enum {
    md_gl_OPTION_NONE                      = 0,
    md_gl_OPTION_RESIDUE_OCCLUSION_CULLING = 1,
};

typedef struct md_gl_draw_args {
    struct {
        uint32_t count;
        const md_gl_representation** data;
    } representation;

    struct {
        // @NOTE: Matrices are column major float[4][4]
        const float* model_view_matrix;
        const float* projection_matrix;

        // [Optional] These fields are only required if md_gl_rendertarget is used and a texture_view_velocity is present
        const float* prev_model_view_matrix; 
        const float* prev_projection_matrix;
    } view_transform;

    const md_gl_rendertarget* render_target;     // [Optional] If NULL, then the global gl state of the viewport and framebuffer will be used

    // If atom_mask is non-zero, the value is used as a mask and is ANDed with the atom flag of the molecule,
    // if the result after the operation is non-zero, the atom will be drawn.
    // Some representations (such as ribbons, cartoon) use spline segments derived from CA atoms, hide the CA atom => hide the segment
    uint32_t atom_mask;

    md_gl_options options;
} md_gl_draw_args_t;

bool md_gl_draw(md_gl_context* ctx, const md_gl_draw_args_t* args);

#ifdef __cplusplus
}
#endif

#endif