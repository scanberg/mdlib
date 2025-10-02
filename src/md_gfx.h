#pragma once
#include <stdint.h>
#include <stdbool.h>

struct md_system_t;
struct vec3_t;
struct mat4_t;

typedef uint32_t md_gfx_config_flags_t;
typedef uint32_t md_gfx_draw_flags_t;

enum {
    MD_GFX_DRAW_DISABLE_CULLING = 1,
};

typedef enum md_gfx_secondary_structure_t {
    MD_GFX_SECONDARY_STRUCTURE_UNKNOWN = 0x00000000,
    MD_GFX_SECONDARY_STRUCTURE_COIL    = 0x000000FF,
    MD_GFX_SECONDARY_STRUCTURE_HELIX   = 0x0000FF00,
    MD_GFX_SECONDARY_STRUCTURE_SHEET   = 0x00FF0000,
} md_gfx_secondary_structure_t;

typedef enum md_gfx_rep_type_t {
    MD_GFX_REP_TYPE_SPACEFILL,
    MD_GFX_REP_TYPE_LICORICE,
    MD_GFX_REP_TYPE_VDW,
    MD_GFX_REP_TYPE_CARTOON,
    MD_GFX_REP_TYPE_RIBBONS,
    MD_GFX_REP_TYPE_GAUSSIAN_SURFACE,
    MD_GFX_REP_TYPE_SOLVENT_EXCLUDED_SURFACE,
} md_gfx_rep_type_t;

typedef struct md_gfx_handle_t {
    uint16_t gen;
    uint16_t idx;
} md_gfx_handle_t;

typedef struct md_gfx_color_t {
    uint8_t r,g,b,a;
} md_gfx_color_t;

typedef struct md_gfx_range_t {
    uint32_t beg_idx;
    uint32_t end_idx;
} md_gfx_range_t;

typedef struct md_gfx_backbone_segment_t {
    uint32_t ca_idx;
    uint32_t c_idx;
    uint32_t o_idx;
} md_gfx_backbone_segment_t;

typedef struct md_gfx_bond_t {
    uint32_t idx[2];
} md_gfx_bond_t;

typedef struct md_gfx_draw_op_t {
    md_gfx_handle_t structure;
    md_gfx_handle_t representation;
    uint32_t mask;                   // Optional: 0 implies no masking is done, otherwise each atom is masked with an bitwise AND against this flag
    const struct mat4_t* model_mat;  // Optional: NULL implies identity matrix
} md_gfx_draw_op_t;

typedef struct md_gfx_rep_attr_t {
    union {
        struct {
            float radius_scale;
        } spacefill;

        struct {
            float radius_scale;
        } licorice;

        struct {
            float ball_radius_scale;
            float stick_radius_scale;
        } vdw;

        struct {
            float profile_width;
            float profile_height;
        } cartoon;

        struct {
            float profile_width;
            float profile_height;
        } ribbons;

        struct {
            float sigma;
        } gaussian_surface;

        struct {
            float probe_radius;
        } solvent_excluded_surface;
    };
} md_gfx_rep_attr_t;

#ifdef __cplusplus
extern "C" {
#endif

// SYSTEM
bool md_gfx_initialize(const char* shader_base_dir, uint32_t fbo_width, uint32_t fbo_height, md_gfx_config_flags_t flags);
void md_gfx_shutdown(void);

// STRUCTURE
md_gfx_handle_t md_gfx_structure_create_from_mol(const struct md_system_t* mol);
// Explicit initialization of data for a structure
md_gfx_handle_t md_gfx_structure_create(uint32_t atom_count, uint32_t bond_count, uint32_t backbone_segment_count, uint32_t backbone_range_count, uint32_t group_count, uint32_t instance_count);
bool md_gfx_structure_destroy(md_gfx_handle_t id);

// ATOM FIELDS
bool md_gfx_structure_set_atom_position_soa(md_gfx_handle_t id, const float* x, const float* y, const float* z, uint32_t count);
bool md_gfx_structure_set_atom_position(md_gfx_handle_t id, const struct vec3_t* xyz, uint32_t count, uint32_t byte_stride);
bool md_gfx_structure_set_atom_radius(md_gfx_handle_t id, const float* radius, uint32_t count, uint32_t byte_stride);
bool md_gfx_structure_set_atom_flags(md_gfx_handle_t id, const uint32_t* flags, uint32_t count, uint32_t byte_stride);

bool md_gfx_structure_set_group_atom_ranges(md_gfx_handle_t id, const md_gfx_range_t* atom_ranges, uint32_t count, uint32_t byte_stride);

// Parts of a structure can be replicated and transformed as instances
// An instance is defined as an atom range and a transform
bool md_gfx_structure_set_instance_atom_ranges(md_gfx_handle_t id, const md_gfx_range_t* atom_ranges, uint32_t count, uint32_t byte_stride);
bool md_gfx_structure_set_instance_transforms(md_gfx_handle_t id, const struct mat4_t* transforms, uint32_t count, uint32_t byte_stride);

// Required for licorice and vdw
void md_gfx_structure_set_bonds(md_gfx_handle_t id, const md_gfx_bond_t* bonds, uint32_t count, uint32_t byte_stride);

// Required for cartoon and ribbons
void md_gfx_structure_set_backbone_segment(md_gfx_handle_t id, const md_gfx_backbone_segment_t* segments, uint32_t count, uint32_t byte_stride);
void md_gfx_structure_set_backbone_secondary_structure(md_gfx_handle_t id, const md_gfx_secondary_structure_t* secondary_structures, uint32_t count, uint32_t byte_stride);
void md_gfx_structure_set_backbone_range(md_gfx_handle_t id, const md_gfx_range_t* backbone_ranges, uint32_t count, uint32_t byte_stride);

// Representations (visual representations of structures)
md_gfx_handle_t md_gfx_rep_create(uint32_t color_count);
bool md_gfx_rep_destroy(md_gfx_handle_t id);

bool md_gfx_rep_set_data(md_gfx_handle_t id, md_gfx_rep_type_t type, md_gfx_rep_attr_t attr, const md_gfx_color_t* color, uint32_t count, uint32_t byte_stride);

// Draw (Commit draw operations to render structures with representations)
bool md_gfx_draw(uint32_t draw_op_count, const md_gfx_draw_op_t* draw_ops, const struct mat4_t* proj_mat, const struct mat4_t* view_mat, const struct mat4_t* inv_proj_mat, const struct mat4_t* inv_view_mat);

// Essentially a render pass

// Call before the submitting draw commands for the render pass
bool md_gfx_draw_begin(const struct mat4_t* proj_mat, const struct mat4_t* view_mat, const struct mat4_t* inv_proj_mat, const struct mat4_t* inv_view_mat, uint32_t draw_flags);

// Submit draw command
// struct_id: Handle to a valid structure
// rep_id: Handle to a valid representation
// model_mat: pointer to model matrix (column major float[4][4]), [OPTIONAL] NULL implies identity matrix.
// atom_mask: mask atoms by a bitwise AND operation, [OPTIONAL] 0 implies all atoms are drawn.
bool md_gfx_draw_submit(md_gfx_handle_t struct_id, md_gfx_handle_t rep_id, const struct mat4_t* model_mat, uint32_t atom_mask);

// Call after draw commands have been submitted for the render pass
bool md_gfx_draw_end(void);

// Call after all render passes are done
bool md_gfx_compose_to_fbo(void);

// Perform this after the draw operation
bool md_gfx_query_picking(uint32_t mouse_x, uint32_t mouse_y);

uint32_t md_gfx_get_picking_idx(void);
float    md_gfx_get_picking_depth(void);

#ifdef __cplusplus
}
#endif
