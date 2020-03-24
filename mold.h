
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

typedef uint32_t mold_error_code;
enum mold_error_code_ {
    MOLD_SUCCESS = 0,
    MOLD_ERROR_UNKNOWN,
    MOLD_ERROR_ARGUMENT_IS_NULL,
    MOLD_ERROR_CONTEXT_NOT_INITIALIZED,
    MOLD_ERROR_MOLECULE_SLOT_OUT_OF_RANGE,
    MOLD_ERROR_REPRESENTATION_SLOT_OUT_OF_RANGE,
    MOLD_ERROR_INVALID_HANDLE,
    MOLD_ERROR_GL_VERSION_NOT_SUPPORTED
};

typedef struct mold_result {
    mold_error_code error_code;
    const char* error_str;
} mold_result;

inline int mold_success(mold_result res) {
    return res.error_code == MOLD_SUCCESS;
}

// Range defining a continous range of indices from [beg, end[ where beg is included, but end is excluded
// Example range{0, 5} would represent [0,1,2,3,4]
typedef struct mold_range {
    uint32_t beg;
    uint32_t end;
} mold_range;

inline uint32_t mold_extent(mold_range range) {
    return range.end - range.beg;
}

/*
 * CONTEXT
 * Interface to initialize and free shared resources within the context
 */
mold_result mold_initialize_context();
mold_result mold_shutdown_context();

/*
 *  MOLECULE
 *  Interface to create and modify molecule data
 */
typedef struct mold_molecule_t* mold_molecule;

typedef uint8_t mold_secondary_structure;
typedef enum {
    MOLD_SECONDARY_STRUCTURE_UNKNOWN = 0,
    MOLD_SECONDARY_STRUCTURE_COIL,
    MOLD_SECONDARY_STRUCTURE_HELIX,
    MOLD_SECONDARY_STRUCTURE_SHEET
} mold_secondary_structure_;

typedef uint32_t mold_atom_idx;

typedef struct mold_position {
    float x, y, z;
} mold_position;

typedef float mold_radius;

typedef struct mold_bond {
    mold_atom_idx idx[2];
} mold_bond;

typedef struct mold_color {
    uint8_t r, g, b, a;
} mold_color;

typedef struct mold_segment_indices {
    mold_atom_idx ca_idx;
    mold_atom_idx c_idx;
    mold_atom_idx o_idx;
} mold_segment_indices;

typedef struct mold_molecule_desc {
    struct {
        uint32_t count;
        const mold_position* position_data;
        const mold_radius*   radius_data;
    } atom;

    struct {
        uint32_t count;
        const mold_range* atom_range_data;
    } residue;

    struct {
        uint32_t count;
        const mold_bond* data;        
    } covalent_bond;

    struct {
        struct {
            uint32_t count;
            const mold_segment_indices*     index_data;
            const mold_secondary_structure* secondary_structure_data;
        } segment;
        struct {
            uint32_t count;
            const mold_range* segment_range_data;
        } sequence;
    } backbone;
} mold_molecule_desc;

mold_result mold_create_molecule(mold_molecule* mol, const mold_molecule_desc* desc);
mold_result mold_destroy_molecule(mold_molecule mol);

mold_result mold_set_atom_position_data(mold_molecule mol, const mold_position pos_xyz[], uint32_t count, uint32_t offset);
mold_result mold_set_atom_position_data_soa(mold_molecule mol, const float pos_x[], const float pos_y[], const float pos_z[], uint32_t count, uint32_t offset);
mold_result mold_set_atom_radius_data(mold_molecule mol, const float radius[], uint32_t count, uint32_t offset);

mold_result mold_set_covalent_bond_data(mold_molecule mol, const mold_bond* bond[], uint32_t count, uint32_t offset);

mold_result mold_set_backbone_segment_index_data(mold_molecule mol, const mold_segment_indices segment_indices[], uint32_t count, uint32_t offset);
mold_result mold_set_backbone_secondary_structure_data(mold_molecule mol, const mold_secondary_structure secondary_structure[], uint32_t count, uint32_t offset);
mold_result mold_set_backbone_sequence_data(mold_molecule mol, const mold_range segment_sequence[], uint32_t count, uint32_t offset);

/*
 *  REPRESENTATIONS
 *  Interface for creating visual representations for molecules
 */
typedef struct mold_representation_t* mold_representation;

typedef uint32_t mold_rep_type;
enum mold_rep_type_ {
    MOLD_REPYPE_SPACE_FILL,
    MOLD_REPYPE_LICORICE,
    MOLD_REPYPE_RIBBONS,
    MOLD_REPYPE_CARTOON,
    MOLD_REPYPE_SOLVENT_EXCLUDED_SURFACE
};

typedef union mold_rep_args {
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
} mold_rep_args;

typedef struct mold_representation_desc {
    mold_rep_type type;
    mold_rep_args args;

    // [Optional] Atom color data, if supplied, it is expected to have the size of the extent of the molecule
    mold_color* atom_color_data;
} mold_representation_desc;

mold_result mold_create_molecule_representation(mold_representation* rep, mold_molecule mol, const mold_representation_desc* desc);
mold_result mold_destroy_molecule_representation(mold_representation rep);

mold_result mold_set_representation_type_and_args(mold_representation rep, mold_rep_type type, mold_rep_args args);
mold_result mold_set_representation_color_data(mold_representation rep, const mold_color color_data[], uint32_t count, uint32_t offset);

/*
 *  DRAW
 *  Interface for drawing representations of molecules
 */

typedef uint32_t mold_texture;

typedef uint32_t mold_draw_flags;
enum mold_draw_flags_ {
    MOLD_DRAW_FLAGS_RESIDUE_OCCLUSION_CULLING = 0x1,
};

typedef struct mold_rendertarget {
    uint32_t width;
    uint32_t height;
    mold_texture depth;
    mold_texture color;
    mold_texture atom_idx;
    mold_texture view_normal;
    mold_texture view_velocity;
} mold_rendertarget;

typedef struct mold_draw_desc {
    mold_molecule molecule;
    struct {
        uint32_t count;
        const mold_representation* data;
    } representation;

    struct {
        float model[4][4];
        float view[4][4];
        float proj[4][4];
    } matrix;

    // (Optional) NULL for forward rendering
    const mold_rendertarget* render_target;
    mold_draw_flags flags;
} mold_draw_desc;

mold_result mold_draw(const mold_draw_desc* desc);

#ifdef __cplusplus
}
#endif