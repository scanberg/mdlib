#ifndef _MOLD_DRAW_H_
#define _MOLD_DRAW_H_

//  MOLD - Molecule Device Draw
//  Copyright 2020 Robin Skånberg
// 
//  LICENSE: see license file
// 
//  Example use case:
//  @TODO: Fill in


#include <stdint.h>
#include "mold_molecule.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MOLD_SPLINE_SUBDIVISION_COUNT 8

typedef uint32_t mold_error;
enum {
    MOLD_DRAW_SUCCESS = 0,
    MOLD_DRAW_UNKNOWN_ERROR,
    MOLD_DRAW_ARGUMENT_IS_NULL,
    MOLD_DRAW_CONTEXT_NOT_INITIALIZED,
    MOLD_DRAW_MOL_INVALID,
    MOLD_DRAW_REP_INVALID,
    MOLD_DRAW_FILE_NOT_FOUND,
    MOLD_DRAW_GL_VERSION_NOT_SUPPORTED,
    MOLD_DRAW_GL_INVALID_BUFFER,
    MOLD_DRAW_GL_ATTEMPTING_TO_WRITE_OUT_OF_BUFFER_RANGE,
    MOLD_DRAW_GL_COULD_NOT_MAP_BUFFER,
    MOLD_DRAW_GL_SHADER_COMPILE_ERROR,
    MOLD_DRAW_GL_PROGRAM_LINK_ERROR,
    MOLD_DRAW_GL_FRAMEBUFFER_ERROR,
};

const char* mold_draw_get_error_str();

// ### DRAW CONTEXT ###
//
// The draw context holds the shared resources such as buffers, shaders and textures needed to draw molecules
typedef struct mold_draw_ctx mold_draw_ctx;
struct mold_draw_ctx {
    uint64_t _mem[72];
};

typedef struct mold_draw_ctx_desc mold_draw_ctx_desc;
struct mold_draw_ctx_desc {
    // This function snipped will be injected into the shaders and will control how the atom index is encoded.
    // If supplied, use the following function declaration and name:
    //      vec4 encode_atom_idx(in uint);
    //
    // If left null, the default version will be used which packs the index in a straight forward way:
    //      vec4 encode_atom_idx(in uint atom_idx) {
    //      return vec4(
    //          (data & 0x000000FFU) >> 0U,
    //          (data & 0x0000FF00U) >> 8U,
    //          (data & 0x00FF0000U) >> 16U,
    //          (data & 0xFF000000U) >> 24U) / 255.0;
    //      }
    const char* encode_atom_idx_src;

    // This function snipped will be injected into the shaders and will control how the view normal is encoded.
    // If supplied, use the following function declaration and name:
    //      vec4 encode_view_normal(in vec3);
    // If left null, the default version will be used:
    //      vec4 encode_view_normal (in vec3 n) {
    //          float p = sqrt(n.z*8+8);
    //          return vec4(n.xy/p + 0.5,0,0);
    //      }
    // This is spherical mapping version, which is listed here https://aras-p.info/texts/CompactNormalStorage.html
    // Note that this version only uses the first two channels
    // The decode function looks like this:
    //      vec3 decode_view_normal(vec2 enc) {
    //          vec2 fenc = enc*4-2;
    //          float f = dot(fenc,fenc);
    //          float g = sqrt(1-f/4.0);
    //          vec3 n;
    //          n.xy = fenc*g;
    //          n.z = 1-f/2.0;
    //          return n;
    //      }
    const char* encode_view_normal_src;
};

mold_error mold_draw_ctx_init(mold_draw_ctx* ctx);
mold_error mold_draw_ctx_free(mold_draw_ctx* ctx);


// ### DRAW MOLECULE ###
//
// The draw molecule represents the state of a single molecule which can be drawn
// Its data can either be initialized from a mold_molecule or the data fields can be set individually
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

typedef struct mold_draw_mol mold_draw_mol;
struct mold_draw_mol {
    uint64_t _mem[12];
};

typedef struct mold_draw_mol_desc mold_draw_mol_desc;

struct mold_draw_mol_desc {
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
        const mold_bond* atom_bond;
    } bond;

    struct {
        uint32_t count;             
        const mold_range*               atom_range;
        const mold_backbone_atoms*      backbone_atoms;
        const mold_secondary_structure* secondary_structure;
    } residue;

    struct {
        uint32_t count;
        const mold_range* residue_range;
    } chain;
};

void       mold_draw_mol_desc_extract(mold_draw_mol_desc* desc, const mold_molecule* mold_mol);
mold_error mold_draw_mol_init(mold_draw_mol* mol, const mold_draw_mol_desc* desc);
mold_error mold_draw_mol_free(mold_draw_mol* mol);

// ### Update by modifying the dynamic data fields ###
mold_error mold_draw_mol_set_atom_position(mold_draw_mol* mol, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t byte_stride);
mold_error mold_draw_mol_set_atom_radius(mold_draw_mol* mol, uint32_t offset, uint32_t count, const float* radius, uint32_t byte_stride);
mold_error mold_draw_mol_set_atom_flags(mold_draw_mol* mol, uint32_t offset, uint32_t count, const uint8_t* flags, uint32_t byte_stride);

mold_error mold_draw_mol_set_residue_secondary_structure(mold_draw_mol* mol, uint32_t offset, uint32_t count, const mold_secondary_structure* secondary_structure, uint32_t byte_stride);

// This is called to copy the atom position buffer to previous atom position
// Usually just before calling mold_draw_mol_set_atom_position() to set a new current position
mold_error mold_draw_mol_update_atom_previous_position(mold_draw_mol* mol);
/*
 *  REPRESENTATIONS
 *  Interface for creating visual representations for molecules
 */
typedef struct mold_draw_rep mold_draw_rep;

struct mold_draw_rep { // This is an opaque blob which matches the size of internal memory used by the structure
    uint64_t _mem[4];
};

typedef uint16_t mold_draw_rep_type;
enum {
    MOLD_DRAW_REP_DEFAULT                    = 0,
    MOLD_DRAW_REP_SPACE_FILL                 = 1,
    MOLD_DRAW_REP_SOLVENT_EXCLUDED_SURFACE   = 2,
    MOLD_DRAW_REP_RIBBONS                    = 3,
    MOLD_DRAW_REP_CARTOON                    = 4,
    MOLD_DRAW_REP_LICORICE                   = 5
};

typedef union mold_draw_rep_args {
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
} mold_draw_rep_args;

mold_error mold_draw_rep_init(mold_draw_rep* rep, const mold_draw_mol* mol);
mold_error mold_draw_rep_free(mold_draw_rep* rep);

mold_error mold_draw_rep_set_type_and_args(mold_draw_rep* rep, mold_draw_rep_type type, mold_draw_rep_args args);
mold_error mold_draw_rep_set_color(mold_draw_rep* rep, uint32_t offset, uint32_t count, const uint32_t* color, uint32_t byte_stride);

/*
 *  DRAW
 *  Interface for drawing representations of molecules
 *
 */

typedef struct mold_draw_rendertarget {
    uint32_t width;
    uint32_t height;
    uint32_t texture_depth;
    uint32_t texture_color;             // [Optional] (0 = ignore)
    uint32_t texture_atom_index;        // [Optional] (0 = ignore)
    uint32_t texture_view_normal;       // [Optional] (0 = ignore)
    uint32_t texture_view_velocity;     // [Optional] (0 = ignore)
} mold_draw_rendertarget;

typedef uint32_t mold_draw_options;
enum {
    MOLD_DRAW_OPTION_NONE                      = 0,
    MOLD_DRAW_OPTION_RESIDUE_OCCLUSION_CULLING = 1,
};

typedef struct mold_draw_desc {
    struct {
        uint32_t count;
        const mold_draw_rep** data;
    } representation;

    struct {
        // @NOTE: Matrices are column major float[4][4]
        const float* model_view_matrix;
        const float* projection_matrix;

        // [Optional] These fields are only required if mold_draw_rendertarget is used and a texture_view_velocity is present
        const float* prev_model_view_matrix; 
        const float* prev_projection_matrix;
    } view_transform;

    const mold_draw_rendertarget* render_target;     // [Optional] If NULL, then the global gl state of the viewport and framebuffer will be used

    // If mol_mask is non-zero, the value is used as a mask and is ANDed with the atom flag of the molecule,
    // if the result after the operation is non-zero, the atom will be drawn.
    // Some representations (such as ribbons, cartoon) use spline segments derived from CA atoms, hide the CA atom => hide the segment
    uint32_t mol_mask;

    mold_draw_options options;
} mold_draw_desc;

mold_error mold_draw(const mold_draw_ctx* ctx, const mold_draw_desc* desc);

#ifdef __cplusplus
}
#endif

#endif