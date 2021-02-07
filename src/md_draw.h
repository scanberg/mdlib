#ifndef _MD_DRAW_H_
#define _MD_DRAW_H_

//  MOLD - Molecule Device Draw
//  Copyright 2020 Robin Skånberg
// 
//  LICENSE: see license file
// 
//  Example use case:
//  @TODO: Fill in

#include <stdint.h>
#include <md_molecule.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MD_SPLINE_MAX_SUBDIVISION_COUNT 8

#define MD_DRAW_SUCCESS                                       0
#define MD_DRAW_UNKNOWN_ERROR                                -1
#define MD_DRAW_ARGUMENT_IS_NULL                             -2
#define MD_DRAW_CONTEXT_INVALID                              -3
#define MD_DRAW_MOLECULE_INVALID                             -4
#define MD_DRAW_REPRESENTATION_INVALID                       -5
#define MD_DRAW_FILE_NOT_FOUND                               -6
#define MD_DRAW_GL_VERSION_NOT_SUPPORTED                     -7
#define MD_DRAW_GL_INVALID_BUFFER                            -8
#define MD_DRAW_GL_ATTEMPTING_TO_WRITE_OUT_OF_BUFFER_RANGE   -9
#define MD_DRAW_GL_COULD_NOT_MAP_BUFFER                     -10
#define MD_DRAW_GL_SHADER_COMPILE_ERROR                     -11
#define MD_DRAW_GL_PROGRAM_LINK_ERROR                       -12
#define MD_DRAW_GL_FRAMEBUFFER_ERROR                        -13

typedef int md_draw_error;

// ### DRAW CONTEXT ###
//
// The draw context holds the shared resources such as buffers, shaders and textures needed to draw molecules
typedef struct md_draw_context md_draw_context;
struct md_draw_context {
    uint64_t _opaque[72];
};

md_draw_error md_draw_context_init(md_draw_context* ctx);
md_draw_error md_draw_context_free(md_draw_context* ctx);


// ### DRAW MOLECULE ###
//
// The draw molecule represents the state of a single molecule which can be drawn
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

typedef struct md_draw_molecule md_draw_molecule;
struct md_draw_molecule {
    uint64_t _mem[12];
};

typedef struct md_draw_molecule_desc md_draw_molecule_desc;
struct md_draw_molecule_desc {
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
        const md_bond* atom_bond;
    } bond;

    struct {
        uint32_t count;
        const md_range*               atom_range;
        const md_backbone_atoms*      backbone_atoms;
        const md_secondary_structure* secondary_structure;
    } residue;

    struct {
        uint32_t count;
        const md_range* residue_range;
    } chain;
};

md_draw_error md_draw_molecule_init(md_draw_molecule* mol, const md_draw_molecule_desc* desc);
md_draw_error md_draw_molecule_free(md_draw_molecule* mol);

// ### Update by modifying the dynamic data fields ###
md_draw_error md_draw_molecule_set_atom_position(md_draw_molecule* mol, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t byte_stride);
md_draw_error md_draw_molecule_set_atom_radius(md_draw_molecule* mol, uint32_t offset, uint32_t count, const float* radius, uint32_t byte_stride);
md_draw_error md_draw_molecule_set_atom_flags(md_draw_molecule* mol, uint32_t offset, uint32_t count, const uint8_t* flags, uint32_t byte_stride);

md_draw_error md_draw_molecule_set_residue_secondary_structure(md_draw_molecule* mol, uint32_t offset, uint32_t count, const md_secondary_structure* secondary_structure, uint32_t byte_stride);

// This is called to copy the atom position buffer to previous atom position
// Usually just before calling md_draw_molecule_set_atom_position() to set a new current position
md_draw_error md_draw_molecule_update_atom_previous_position(md_draw_molecule* mol);
/*
 *  REPRESENTATIONS
 *  Interface for creating visual representations for molecules
 */
typedef struct md_draw_representation md_draw_representation;

struct md_draw_representation { // This is an opaque blob which matches the size of internal memory used by the structure
    uint64_t _mem[4];
};

typedef uint16_t md_draw_representation_type;
enum {
    MD_DRAW_REP_DEFAULT                    = 0,
    MD_DRAW_REP_SPACE_FILL                 = 1,
    MD_DRAW_REP_SOLVENT_EXCLUDED_SURFACE   = 2,
    MD_DRAW_REP_RIBBONS                    = 3,
    MD_DRAW_REP_CARTOON                    = 4,
    MD_DRAW_REP_LICORICE                   = 5
};

typedef union md_draw_representation_args {
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
} md_draw_representation_args;

md_draw_error md_draw_representation_init(md_draw_representation* rep, const md_draw_molecule* mol);
md_draw_error md_draw_representation_free(md_draw_representation* rep);

md_draw_error md_draw_representation_set_type_and_args(md_draw_representation* rep, md_draw_representation_type type, md_draw_representation_args args);
md_draw_error md_draw_representation_set_color(md_draw_representation* rep, uint32_t offset, uint32_t count, const uint32_t* color, uint32_t byte_stride);

/*
 *  DRAW
 *  Interface for defining representations of molecules
 *
 */

typedef struct md_draw_rendertarget {
    uint32_t width;
    uint32_t height;
    uint32_t texture_depth;
    uint32_t texture_color;             // [Optional] (0 = ignore)
    uint32_t texture_atom_index;        // [Optional] (0 = ignore)
    uint32_t texture_view_normal;       // [Optional] (0 = ignore)
    uint32_t texture_view_velocity;     // [Optional] (0 = ignore)
} md_draw_rendertarget;

typedef uint32_t md_draw_options;
enum {
    MD_DRAW_OPTION_NONE                      = 0,
    MD_DRAW_OPTION_RESIDUE_OCCLUSION_CULLING = 1,
};

typedef struct md_draw_desc {
    struct {
        uint32_t count;
        const md_draw_representation** data;
    } representation;

    struct {
        // @NOTE: Matrices are column major float[4][4]
        const float* model_view_matrix;
        const float* projection_matrix;

        // [Optional] These fields are only required if md_draw_rendertarget is used and a texture_view_velocity is present
        const float* prev_model_view_matrix; 
        const float* prev_projection_matrix;
    } view_transform;

    const md_draw_rendertarget* render_target;     // [Optional] If NULL, then the global gl state of the viewport and framebuffer will be used

    // If mol_mask is non-zero, the value is used as a mask and is ANDed with the atom flag of the molecule,
    // if the result after the operation is non-zero, the atom will be drawn.
    // Some representations (such as ribbons, cartoon) use spline segments derived from CA atoms, hide the CA atom => hide the segment
    uint32_t mol_mask;

    md_draw_options options;
} md_draw_desc;

md_draw_error md_draw(md_draw_context* ctx, const md_draw_desc* desc);

#ifdef __cplusplus
}
#endif

#endif