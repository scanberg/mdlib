#pragma once
#include <stdint.h>
#include <stdbool.h>

struct md_molecule_t;
struct md_mat4_t;

typedef uint32_t md_gfx_config_flags_t;

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

typedef struct md_gfx_instance_t {
    md_gfx_range_t group_range;
    const struct mat4_t* transform;
} md_gfx_instance_t;

typedef struct md_gfx_draw_op_t {
    md_gfx_handle_t structure;
    md_gfx_handle_t representation;
    const struct mat4_t* model_mat;  // Optional: NULL implies identity matrix
} md_gfx_draw_op_t;

// This could conceptually be a union, but there are benefits in keeping state of other representation types in case the user wants to swap back and forth.
typedef struct md_gfx_rep_attr_t {
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
} md_gfx_rep_attr_t;

#ifdef __cplusplus
extern "C" {
#endif

// System
bool md_gfx_initialize(uint32_t width, uint32_t height, md_gfx_config_flags_t flags);
void md_gfx_shutdown();

// Structures (holds common data)
md_gfx_handle_t md_gfx_structure_create(uint32_t atom_count, uint32_t bond_count, uint32_t backbone_segment_count, uint32_t backbone_range_count, uint32_t group_count, uint32_t instance_count);
bool md_gfx_structure_destroy(md_gfx_handle_t id);

bool md_gfx_structure_init_from_mol(md_gfx_handle_t id, const struct md_molecule_t* mol);

bool md_gfx_structure_set_atom_position(md_gfx_handle_t id, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t stride);
bool md_gfx_structure_set_atom_velocity(md_gfx_handle_t id, uint32_t offset, uint32_t count, const float* vx, const float* vy, const float* vz, uint32_t stride);

bool md_gfx_structure_zero_atom_velocity(md_gfx_handle_t id);

// Required for spacefill and surfaces
bool md_gfx_structure_set_atom_radius(md_gfx_handle_t id, uint32_t offset, uint32_t count, const float* radius, uint32_t stride);

// Groups are spatially coherent ranges of atoms which are suitable for culling
// Residues generally map well to this, but larger ranges could also be used as
// long as they remain spatially 'packed' and form a tight bounding box.
bool md_gfx_structure_set_group_ranges(md_gfx_handle_t id, uint32_t offset, uint32_t count, const md_gfx_range_t* ranges, uint32_t stride);

// Parts of a structure can be replicated and transformed as instances
// This is defined through an atom_range and a transform
void md_gfx_structure_set_instances(md_gfx_handle_t id, uint32_t offset, uint32_t count, const md_gfx_instance_t* instances, uint32_t stride);

// Required for licorice and vdw
void md_gfx_structure_set_bonds(md_gfx_handle_t id, uint32_t offset, uint32_t count, const md_gfx_bond_t* bonds, uint32_t stride);

// Required for cartoon and ribbons
void md_gfx_structure_set_backbone_segment(md_gfx_handle_t id, uint32_t offset, uint32_t count, const md_gfx_backbone_segment_t* segments, uint32_t stride);
void md_gfx_structure_set_backbone_secondary_structure(md_gfx_handle_t id, uint32_t offset, uint32_t count, const md_gfx_secondary_structure_t* secondary_structures, uint32_t stride);
void md_gfx_structure_set_backbone_range(md_gfx_handle_t id, uint32_t offset, uint32_t count, const md_gfx_range_t* backbone_ranges, uint32_t stride);

// Representations (visual representations of structures)
md_gfx_handle_t md_gfx_rep_create(uint32_t color_count);
bool md_gfx_rep_destroy(md_gfx_handle_t id);

bool md_gfx_rep_set_type_and_attr(md_gfx_handle_t id, md_gfx_rep_type_t type, const md_gfx_rep_attr_t* attr);
bool md_gfx_rep_set_color(md_gfx_handle_t id, uint32_t offset, uint32_t count, const md_gfx_color_t* color, uint32_t stride);

// Draw (Commit draw operations to render structures with representations)
bool md_gfx_draw(uint32_t draw_op_count, const md_gfx_draw_op_t* draw_ops, const struct mat4_t* proj_mat, const struct mat4_t* view_mat);

#ifdef __cplusplus
}
#endif
