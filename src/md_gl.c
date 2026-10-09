#include <md_gl.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_compiler.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_os.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_handle.h>
#include <core/md_gl_util.h>

#include <md_system.h>

#include <stdbool.h>
#include <string.h>
#include <stdio.h>      // printf etc
#include <stddef.h>     // offsetof

#include <GL/gl3w.h>

// Baked shaders
#include <gl_shaders.inl>

#define PUSH_GPU_SECTION(lbl)                                                                   \
{                                                                                               \
    if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
}
#define POP_GPU_SECTION()                       \
{                                           \
    if (glPopDebugGroup) glPopDebugGroup(); \
}

#define MAGIC 0xfacb8172U

// Resolution of the backbone representations (cartoon, ribbons)
// Segments along the spline per residue (even)
#ifndef MD_GL_BACKBONE_SEGMENT_COUNT
#define MD_GL_BACKBONE_SEGMENT_COUNT 12
#endif
// Vertices around the elliptic cross section of the cartoon (even)
#ifndef MD_GL_BACKBONE_PROFILE_COUNT
#define MD_GL_BACKBONE_PROFILE_COUNT 16
#endif
// Segments across the wide faces of the ribbons
#define BACKBONE_RIBBONS_FACE_SUBDIVISIONS 4
// Levels of detail, each halves the resolution of the previous (or less, the segments must divide the finest level)
#define BACKBONE_LOD_COUNT 3
// A chain gets the coarsest level of detail whose segments are at most this many pixels long where it is closest to
// the camera
#define BACKBONE_LOD_SEGMENT_PIXELS 2.5f
// Length of the backbone spline per residue (CA - CA distance, Ångström)
#define BACKBONE_RESIDUE_LENGTH 3.8f

#define UBO_SIZE (1 << 10)
#define SHADER_BUF_SIZE KILOBYTES(14)


static const str_t default_shader_output = STR_INIT(
"layout(location = 0) out vec4 out_color;\n"
"void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint atom_index) {\n"
"    out_color = color;\n"
"}\n");


enum {
    GL_VERSION_UNKNOWN = 0,
    GL_VERSION_330
};

enum {
    GL_PROGRAM_EXTRACT_CONTROL_POINTS,
    GL_PROGRAM_ORIENT_CONTROL_POINTS,
    GL_PROGRAM_SUBDIVIDE_SPLINE,
    GL_PROGRAM_COMPUTE_VELOCITY,
    GL_PROGRAM_COUNT
};

enum {
    GL_TEXTURE_BUFFER_0,
    GL_TEXTURE_BUFFER_1,
    GL_TEXTURE_BUFFER_2,
    GL_TEXTURE_BUFFER_3,
    GL_TEXTURE_BUFFER_4,
    GL_TEXTURE_MAX_DEPTH,
    GL_TEXTURE_COUNT
};

enum {
    GL_BUFFER_ATOM_POSITION,                   // vec3
    GL_BUFFER_ATOM_POSITION_PREV,              // vec3
    GL_BUFFER_ATOM_RADIUS,                     // float
    GL_BUFFER_ATOM_VELOCITY,                   // vec3
    GL_BUFFER_ATOM_FLAGS,                      // u8
    GL_BUFFER_BOND_ATOM_INDICES,               // u32[2]
    GL_BUFFER_BACKBONE_DATA,                   // u32: residue index, u32: residue atom offset, u8: CA index, C index and O Index, u8: flags
    GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE,    // u8[4]  (0: Unknown, 1: Coil, 2: Helix, 3: Sheet)
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA,     // Extracted control points
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT,      // Oriented control points (support vector + relation to the next control point)
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT_PREV, // The previous computation's oriented control points, swapped with the one above each computation
    GL_BUFFER_BACKBONE_NEIGHBOR,               // u32[4] per control point: prev, next, prev2, next2 (clamped to the chain)
    GL_BUFFER_BACKBONE_RING_DATA,              // Evaluated cross section frames of the spline, MD_GL_BACKBONE_SEGMENT_COUNT per control point (see compute_spline_subdivide.vert)
    GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX,    // u32, LINE_STRIP_ADJACENCY Indices for legacy and debugging paths
    GL_BUFFER_INSTANCE_TRANSFORM,              // mat4 instance transformation matrices
    GL_BUFFER_COUNT
};

enum {
    PERMUTATION_BIT_ORTHO = 1
};

enum {
    MOL_FLAG_HAS_BACKBONE = 1
};

enum {
    DRAW_FLAG_COMPUTE_BACKBONE_SPLINE = 1
};

#define MAX_SHADERS 32
#define MAX_MOLECULES 32
#define MAX_PALETTES 32
#define MAX_REPRESENTATIONS 256
#define MAX_SHADER_PERMUTATIONS 2

static inline bool is_ortho_proj_matrix(const mat4_t M) { return M.elem[2][3] == 0.0f; }

static inline void extract_jitter_uv(float jitter[2], const mat4_t  M) {
    if (is_ortho_proj_matrix(M)) {
        jitter[0] = -M.elem[3][0] * 0.5f;
        jitter[1] = -M.elem[3][1] * 0.5f;
    }
    else {
        jitter[0] = M.elem[2][0] * 0.5f;
        jitter[1] = M.elem[2][1] * 0.5f;
    }
}

typedef struct {
    GLuint id;
    GLenum usage_hint;
    size_t size;
} gl_buffer_t;

typedef struct {
    GLuint id;
    uint32_t width;
    uint32_t height;
} gl_texture_t;

typedef struct {
    GLuint id;
} gl_program_t;

typedef struct {
    uint32_t bb_seg_idx;            // Backbone Segment index (To access secondary structure)
    uint32_t atom_base_idx;         // Residue atom base offset
    uint8_t  ca_idx, c_idx, o_idx;  // Local residue indices, Add to base offset of residue to get global atom index
    uint8_t  flags;
} gl_backbone_data_t;

typedef struct {
    float position[3];
    uint32_t atom_idx;
    float velocity[3];
    float segment_t;                    // @NOTE: stores the segment index (integer part) and the fraction within the segment (fractional part)
    uint8_t secondary_structure[3];     // @NOTE: Secondary structure as fractions within the components (0 = coil, 1 = helix, 2 = sheet) so we later can smoothly blend between them when subdividing.
    uint8_t flags;                      // @NOTE: Bitfield for setting flags, bits: [0] beg_chain, [1] end_chain, [2] beg_secondary_structure, [3] end_secondary_structure, [4] flip_next (support of the next control point is used flipped, set by the orient pass)
    int16_t support_vector[3];
    int16_t tangent_vector[3];
} gl_control_point_t;

typedef struct {
    uint32_t idx[4];
} gl_neighbor_idx_t;

// Written by compute_spline_subdivide.vert
typedef struct {
    uint32_t data[12];
} gl_backbone_ring_t;

typedef struct {
    mat4_t world_to_view;
    mat4_t world_to_view_normal;
    mat4_t world_to_clip;
    mat4_t view_to_clip;
    mat4_t view_to_world;
    mat4_t clip_to_view;
    mat4_t prev_world_to_clip;
    mat4_t curr_view_to_prev_clip;
} gl_view_transform_t;

// Shared ubo data for all shaders
typedef struct {
    gl_view_transform_t view_transform;
    vec4_t jitter_uv;
    uint32_t atom_mask;
    uint32_t atom_index_base;
    uint32_t bond_index_base;
    uint32_t _pad[1];
} gl_ubo_base_t;

typedef struct {
    uint32_t atom_idx[2];
} gl_bond_t;

typedef struct {
    uint32_t id;
    gl_program_t spacefill[MAX_SHADER_PERMUTATIONS];
    gl_program_t licorice[MAX_SHADER_PERMUTATIONS];
    gl_program_t ribbons[MAX_SHADER_PERMUTATIONS];
    gl_program_t cartoon[MAX_SHADER_PERMUTATIONS];
} shaders_t;

// CPU side bounds of the backbone chains, for culling and level of detail
typedef struct {
    uint32_t  chain_count;
    uint32_t* chain_offset;     // [chain_count + 1] Control point ranges of the chains
    uint32_t* ca_atom_idx;      // [control point count]
    vec3_t*   ca_xyz;           // [control point count] Last positions of the control points (CA)
    vec4_t*   chain_sphere;     // [chain_count] Bounding sphere of the control points of the chain (xyz center, w radius)
    bool      valid;            // The positions have been set
} backbone_bounds_t;

typedef struct {
    uint32_t id;
    uint32_t flags;

    uint32_t atom_count;
    uint32_t comp_count;
    uint32_t bond_count;
    uint32_t backbone_count;
    uint32_t backbone_control_point_index_count;

    // GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT holds a previous computation that the next one can
    // continue from. Cleared on creation and by md_gl_mol_reset_backbone_history.
    bool backbone_orient_history;

    backbone_bounds_t backbone_bounds;

    gl_buffer_t buffer[GL_BUFFER_COUNT];
} molecule_t;

typedef struct {
    uint32_t id;
    gl_texture_t tex;
} palette_t;

enum {
    BACKBONE_PROFILE_ELLIPSE,   // Cartoon
    BACKBONE_PROFILE_BOX,       // Ribbons
    BACKBONE_PROFILE_COUNT
};

// Static triangulation of one backbone instance (the part of the spline that belongs to one residue)
typedef struct {
    uint32_t segments;          // S: segments along the spline, S + 1 rings
    uint32_t profile_count;     // P: vertices per ring
    uint32_t outline_count;     // O: vertices on the outline of a cap
    uint32_t face_subdivisions; // K: segments across the wide faces of the box profile
    uint32_t index_count;
    gl_buffer_t index_buffer;   // u16, GL_TRIANGLES
} backbone_mesh_t;

typedef struct {
    uint32_t id;
    uint32_t mol_id; // molecule
    uint32_t pal_id; // palette
    gl_buffer_t atom_color;
} representation_t;

typedef struct {
    // Shared resources
    GLuint vao;
    GLuint fbo;
    gl_buffer_t ubo;
    gl_buffer_t instance_ubo;
    gl_texture_t texture[GL_TEXTURE_COUNT];
    gl_program_t program[GL_PROGRAM_COUNT];
    backbone_mesh_t backbone_mesh[BACKBONE_PROFILE_COUNT][BACKBONE_LOD_COUNT];
    uint32_t version;

    // Handle pools
    md_handle_pool_t shader_pool;
    md_handle_pool_t molecule_pool;
    md_handle_pool_t palette_pool;
    md_handle_pool_t representation_pool;

    // Internal structures which are referenced by handles
    shaders_t shaders[MAX_SHADERS];
    molecule_t molecules[MAX_MOLECULES];
    palette_t palettes[MAX_PALETTES];
    representation_t representations[MAX_REPRESENTATIONS];

    md_allocator_i* arena;
} context_t;

static context_t ctx = {0};

static inline shaders_t* shad_lookup(uint32_t id) {
    if (id == MD_HANDLE_INVALID_ID) {
        return NULL;
    }
    int index = md_handle_index(id);
    if (index >= ctx.shader_pool.size || ctx.shaders[index].id != id) {
        MD_LOG_ERROR("Invalid id");
        return NULL;
    }
    return ctx.shaders + index;
}

static inline molecule_t* mol_lookup(uint32_t id) {
    if (id == MD_HANDLE_INVALID_ID) {
        return NULL;
    }
    int index = md_handle_index(id);
    if (index >= ctx.molecule_pool.size || ctx.molecules[index].id != id) {
        MD_LOG_ERROR("Invalid id");
        return NULL;
    }
    return ctx.molecules + index;
}

static inline palette_t* pal_lookup(uint32_t id) {
    if (id == MD_HANDLE_INVALID_ID) {
        return NULL;
    }
    int index = md_handle_index(id);
    if (index >= ctx.palette_pool.size || ctx.palettes[index].id != id) {
        MD_LOG_ERROR("Invalid id");
        return NULL;
    }
    return ctx.palettes + index;
}

static inline representation_t* rep_lookup(uint32_t id) {
    if (id == MD_HANDLE_INVALID_ID) {
        return NULL;
    }
    int index = md_handle_index(id);
    if (index >= ctx.representation_pool.size || ctx.representations[index].id != id) {
        MD_LOG_ERROR("Invalid id");
        return NULL;
    }
    return ctx.representations + index;
}

static inline bool validate_context(void) {
    if (ctx.version == 0) {
        MD_LOG_ERROR("MD GL module has not been initialized");
        return false;
    }
    return true;
}


static inline gl_buffer_t gl_buffer_create(uint32_t num_bytes, const void* data, GLenum usage_hint) {
    gl_buffer_t buf = {0};
    glGenBuffers(1, &buf.id);
    glBindBuffer(GL_ARRAY_BUFFER, buf.id);
    glBufferData(GL_ARRAY_BUFFER, num_bytes, data, usage_hint);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    buf.usage_hint = usage_hint;
    buf.size = num_bytes;
    return buf;
}

static inline uint32_t gl_buffer_size(gl_buffer_t buf) {
    GLint size;
    glBindBuffer(GL_ARRAY_BUFFER, buf.id);
    glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    return (uint32_t)size;
}

static inline void gl_buffer_conditional_delete(gl_buffer_t* buf) {
    if (buf->id) {
        GLuint id = buf->id;
        glDeleteBuffers(1, &id);
        buf->id = 0;
    }
}

static inline void gl_buffer_set_data(gl_buffer_t buf, uint32_t num_bytes, const void* data) {
    glBindBuffer(GL_ARRAY_BUFFER, buf.id);
    glBufferData(GL_ARRAY_BUFFER, num_bytes, data, buf.usage_hint);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    buf.size = num_bytes;
}

static inline void gl_buffer_set_sub_data(gl_buffer_t buf, uint32_t byte_offset, uint32_t byte_size, const void* data) {
    glBindBuffer(GL_ARRAY_BUFFER, buf.id);
    glBufferSubData(GL_ARRAY_BUFFER, byte_offset, byte_size, data);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static inline void gl_buffer_clear(gl_buffer_t buf) {
    glBindBuffer(GL_ARRAY_BUFFER, buf.id);
    //    if (ctx.version >= 430) {
    //        uint8_t data = 0;
    //        glClearBufferSubData(GL_ARRAY_BUFFER, GL_R8UI, 0, buf.size, GL_RED, GL_UNSIGNED_BYTE, &data);
    //    } else {
    char* ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    if (!ptr) {
        MD_LOG_ERROR("Failed to map buffer");
        return;
    }
    MEMSET(ptr, 0, buf.size);
    glUnmapBuffer(GL_ARRAY_BUFFER);
    //    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static void backbone_bounds_free(backbone_bounds_t* bounds) {
    md_allocator_i* alloc = md_get_heap_allocator();
    if (bounds->chain_offset) {
        const uint32_t cp_count = bounds->chain_offset[bounds->chain_count];
        md_free(alloc, bounds->ca_atom_idx,  cp_count * sizeof(uint32_t));
        md_free(alloc, bounds->ca_xyz,       cp_count * sizeof(vec3_t));
        md_free(alloc, bounds->chain_sphere, bounds->chain_count * sizeof(vec4_t));
        md_free(alloc, bounds->chain_offset, (bounds->chain_count + 1) * sizeof(uint32_t));
    }
    MEMSET(bounds, 0, sizeof(backbone_bounds_t));
}

// xyz holds the positions of the atoms [offset, offset + count)
static void backbone_bounds_update(backbone_bounds_t* bounds, uint32_t offset, uint32_t count, const vec3_t* xyz) {
    if (!bounds->chain_count) return;
    const uint32_t cp_count = bounds->chain_offset[bounds->chain_count];
    for (uint32_t i = 0; i < cp_count; ++i) {
        const uint32_t atom_idx = bounds->ca_atom_idx[i];
        if (offset <= atom_idx && atom_idx < offset + count) {
            bounds->ca_xyz[i] = xyz[atom_idx - offset];
        }
    }
    for (uint32_t c = 0; c < bounds->chain_count; ++c) {
        const uint32_t beg = bounds->chain_offset[c];
        const uint32_t end = bounds->chain_offset[c + 1];
        if (beg == end) {
            bounds->chain_sphere[c] = (vec4_t){0};
            continue;
        }
        vec3_t min_box = bounds->ca_xyz[beg];
        vec3_t max_box = bounds->ca_xyz[beg];
        for (uint32_t i = beg + 1; i < end; ++i) {
            min_box = vec3_min(min_box, bounds->ca_xyz[i]);
            max_box = vec3_max(max_box, bounds->ca_xyz[i]);
        }
        const vec3_t center = vec3_mul1(vec3_add(min_box, max_box), 0.5f);
        float r2 = 0.0f;
        for (uint32_t i = beg; i < end; ++i) {
            r2 = MAX(r2, vec3_distance_squared(center, bounds->ca_xyz[i]));
        }
        bounds->chain_sphere[c] = (vec4_t){center.x, center.y, center.z, sqrtf(r2)};
    }
    bounds->valid = true;
}

// The buffers hold packed xyz, as the positions come: one upload, after the current positions are
// kept as the previous ones (the renderer derives motion from the two).
void md_gl_mol_set_atom_position(md_gl_mol_t handle, uint32_t offset, uint32_t count, const vec3_t* xyz) {
    if (xyz == NULL) {
        MD_LOG_ERROR("Missing argument: xyz");
        return;
    }
    molecule_t* mol = mol_lookup(handle.id);
    if (mol) {
        if (!mol->buffer[GL_BUFFER_ATOM_POSITION].id) {
            MD_LOG_ERROR("Molecule position buffer missing");
            return;
        }
        if (!mol->buffer[GL_BUFFER_ATOM_POSITION_PREV].id) {
            MD_LOG_ERROR("Molecule previous position buffer missing");
            return;
        }
        if (offset + count > mol->atom_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return;
        }

        // Copy position buffer to previous position buffer
        uint32_t buffer_size = gl_buffer_size(mol->buffer[GL_BUFFER_ATOM_POSITION]);
        glBindBuffer(GL_COPY_READ_BUFFER, mol->buffer[GL_BUFFER_ATOM_POSITION].id);
        glBindBuffer(GL_COPY_WRITE_BUFFER, mol->buffer[GL_BUFFER_ATOM_POSITION_PREV].id);
        glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, buffer_size);
        glBindBuffer(GL_COPY_READ_BUFFER, 0);
        glBindBuffer(GL_COPY_WRITE_BUFFER, 0);

        gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_ATOM_POSITION], offset * sizeof(vec3_t), count * sizeof(vec3_t), xyz);
        backbone_bounds_update(&mol->backbone_bounds, offset, count, xyz);
    }
}

void md_gl_mol_set_atom_velocity(md_gl_mol_t handle, uint32_t offset, uint32_t count, const vec3_t* xyz) {
    if (xyz == NULL) {
        MD_LOG_ERROR("Missing argument: xyz");
        return;
    }
    molecule_t* mol = mol_lookup(handle.id);
    if (mol) {
        if (!mol->buffer[GL_BUFFER_ATOM_VELOCITY].id) {
            MD_LOG_ERROR("Molecule velocity buffer missing");
            return;
        }
        if (offset + count > mol->atom_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return;
        }
        gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_ATOM_VELOCITY], offset * sizeof(vec3_t), count * sizeof(vec3_t), xyz);
    }
}

void md_gl_mol_set_atom_radius(md_gl_mol_t handle, uint32_t offset, uint32_t count, const float* radius, uint32_t byte_stride) {
    if (!radius) {
        MD_LOG_ERROR("One or more arguments are missing");
        return;
    }
    molecule_t* mol = mol_lookup(handle.id);
    if (mol) {
        if (!mol->buffer[GL_BUFFER_ATOM_RADIUS].id) {
            MD_LOG_ERROR("Molecule radius buffer missing");
            return;
        }
        if (offset + count > mol->atom_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return;
        }
        byte_stride = MAX(sizeof(float), byte_stride);
        if (byte_stride > sizeof(float)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_RADIUS].id);
            float* radius_data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (radius_data == NULL) {
                MD_LOG_ERROR("Failed to map molecule radius buffer");
                return;
            }
            for (uint32_t i = offset; i < count; ++i) {
                radius_data[i] = *(const float*)((const uint8_t*)radius + byte_stride * i);
            }
            glUnmapBuffer(GL_ARRAY_BUFFER);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_ATOM_RADIUS], offset * sizeof(float), count * sizeof(float), radius);
        }
    }
}

void md_gl_mol_set_atom_flags(md_gl_mol_t handle, uint32_t offset, uint32_t count, const uint8_t* flag_data, uint32_t byte_stride) {
    if (!flag_data) {
        MD_LOG_ERROR("One or more arguments are missing");
        return;
    }
    molecule_t* mol = mol_lookup(handle.id);
    if (mol) {
        if (!mol->buffer[GL_BUFFER_ATOM_FLAGS].id) {
            MD_LOG_ERROR("Molecule flags buffer missing");
            return;
        }
        if (offset + count > mol->atom_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return;
        }
        byte_stride = MAX(sizeof(uint8_t), byte_stride);
        if (byte_stride > sizeof(uint8_t)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
            uint8_t* data = (uint8_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (data == NULL) {
                MD_LOG_ERROR("Failed to map molecule flags buffer");
                return;
            }
            for (uint32_t i = offset; i < count; ++i) {
                data[i] = *(flag_data + byte_stride * i);
            }
            glUnmapBuffer(GL_ARRAY_BUFFER);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_ATOM_FLAGS], offset * sizeof(uint8_t), count * sizeof(uint8_t), flag_data);
        }
    }
}

void md_gl_mol_set_bonds(md_gl_mol_t handle, uint32_t offset, uint32_t count, const md_atom_pair_t* bonds, uint32_t byte_stride) {
    if (bonds == NULL) {
        MD_LOG_ERROR("One or more arguments are missing");
        return;
    }

    molecule_t* mol = mol_lookup(handle.id);
    if (mol) {
        if (!mol->buffer[GL_BUFFER_BOND_ATOM_INDICES].id) {
            MD_LOG_ERROR("Molecule bond buffer buffer missing");
            return;
        }

        if (offset > 0 && offset + count > mol->bond_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return;
        }

        if (offset == 0) {
            gl_buffer_set_data(mol->buffer[GL_BUFFER_BOND_ATOM_INDICES], count * sizeof(md_atom_pair_t), bonds);
            mol->bond_count = count;
            return;
        }

        byte_stride = MAX(sizeof(md_atom_pair_t), byte_stride);
        glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_BOND_ATOM_INDICES].id);
        gl_bond_t* bond_data = (gl_bond_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (bond_data == NULL) {
            MD_LOG_ERROR("Failed to map molecule bond buffer");
            return;
        }
        for (uint32_t i = offset; i < count; ++i) {
            const md_atom_pair_t* bond = (const md_atom_pair_t*)((const uint8_t*)bonds + byte_stride * i);
			uint32_t dst = offset + i;
            bond_data[dst].atom_idx[0] = bond->idx[0];
            bond_data[dst].atom_idx[1] = bond->idx[1];
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

// Essentially packing normalized secondary structure floats into a uint32_t
static inline uint32_t pack_gl_secondary_structure(md_gl_secondary_structure_t ss) {
    float helix = CLAMP(ss.helix, 0.0f, 1.0f);
    float sheet = CLAMP(ss.sheet, 0.0f, 1.0f);
    float coil  = CLAMP(1.0f - (helix + sheet), 0.0f, 1.0f);

    uint32_t packed = 0;
    packed |= (uint8_t)(coil  * 255.0f) << 0;
    packed |= (uint8_t)(helix * 255.0f) << 8;
    packed |= (uint8_t)(sheet * 255.0f) << 16;
    return packed;
}

void md_gl_mol_set_backbone_secondary_structure(md_gl_mol_t handle, uint32_t offset, uint32_t count, const md_gl_secondary_structure_t* secondary_structure, uint32_t byte_stride) {
    if (secondary_structure == NULL) {
        MD_LOG_ERROR("One or more arguments are missing");
        return;
    }
    molecule_t* mol = mol_lookup(handle.id);
    if (mol) {
        if (!mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE].id) {
            MD_LOG_ERROR("Molecule secondary structure buffer missing");
            return;
        }
        if (offset + count > mol->backbone_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return;
        }
        byte_stride = MAX(sizeof(md_gl_secondary_structure_t), byte_stride);
        glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE].id);
        uint32_t* buffer_data = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (buffer_data == NULL) {
            MD_LOG_ERROR("Failed to map molecule secondary structure buffer");
            return;
        }
        for (uint32_t i = offset; i < count; ++i) {
			const uint8_t* secondary_structure_raw = (const uint8_t*)secondary_structure + byte_stride * i;
            md_gl_secondary_structure_t ss = *(const md_gl_secondary_structure_t*)secondary_structure_raw;
			uint32_t dst = offset + i;
            buffer_data[dst] = pack_gl_secondary_structure(ss);
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

void md_gl_mol_reset_backbone_history(md_gl_mol_t handle) {
    molecule_t* mol = mol_lookup(handle.id);
    if (mol) {
        mol->backbone_orient_history = false;
    }
}

void md_gl_mol_compute_velocity(md_gl_mol_t handle, const float pbc_ext[3]) {
    if (!validate_context()) {
        return;
    }

    if (pbc_ext == NULL) {
        MD_LOG_ERROR("Missing one or more arguments");
        return;
    }

    molecule_t* mol = mol_lookup(handle.id);
    if (!mol) {
        return;
    }

    if (mol->buffer[GL_BUFFER_ATOM_VELOCITY].id == 0) {
        MD_LOG_ERROR("Velocity buffer is zero");
        return;
    }

    glEnable(GL_RASTERIZER_DISCARD);
    glBindVertexArray(ctx.vao);

    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_POSITION_PREV].id);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    GLuint program = ctx.program[GL_PROGRAM_COMPUTE_VELOCITY].id;

    const GLint pbc_ext_loc = glGetUniformLocation(program, "u_pbc_ext");

    glUseProgram(program);
    glUniform3fv(pbc_ext_loc, 1, pbc_ext);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    glBeginTransformFeedback(GL_POINTS);
    glDrawArrays(GL_POINTS, 0, mol->atom_count);
    glEndTransformFeedback();
    glUseProgram(0);

    glBindVertexArray(0);
    glDisable(GL_RASTERIZER_DISCARD);
}

void md_gl_mol_zero_velocity(md_gl_mol_t handle) {
    molecule_t* mol = mol_lookup(handle.id);
    if (!mol) {
        return;
    }
    if (!mol->buffer[GL_BUFFER_ATOM_VELOCITY].id) {
        MD_LOG_ERROR("Molecule position buffer missing");
        return;
    }
    gl_buffer_clear(mol->buffer[GL_BUFFER_ATOM_VELOCITY]);
}

// geom_src is optional. defines is injected after the version line of every stage.
bool create_permuted_program(str_t identifier, gl_program_t* program_permutations, str_t vert_src, str_t geom_src, str_t frag_src, str_t frag_output_src, str_t defines) {
    const bool has_geom = !str_empty(geom_src);
    GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint geom_shader = has_geom ? glCreateShader(GL_GEOMETRY_SHADER) : 0;
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

    const str_t perm_str[] = {
        STR_INIT("#define ORTHO 0"),
        STR_INIT("#define ORTHO 1"),
    };
    
    ASSERT(ARRAY_SIZE(perm_str) <= MAX_SHADER_PERMUTATIONS);

    for (uint32_t perm = 0; perm < MAX_SHADER_PERMUTATIONS; ++perm) {
        md_gl_shader_src_injection_t injections[] = { {perm_str[perm], {0}}, {defines, {0}}, {frag_output_src, STR_INIT("EXTRA_SRC")} };

        if (!str_empty(vert_src) && !md_gl_shader_compile(vert_shader, vert_src, injections, 2)) {
            MD_LOG_ERROR("Error occured when compiling vertex shader for: '%.*s'", STR_ARG(identifier));
            return false;
        }
            
        if (has_geom && !md_gl_shader_compile(geom_shader, geom_src, injections, 2)) {
            MD_LOG_ERROR("Error occured when compiling geometry shader for: '%.*s'", STR_ARG(identifier));
            return false;
        }
        if (!str_empty(frag_src) && !md_gl_shader_compile(frag_shader, frag_src, injections, 3)) {
            MD_LOG_ERROR("Error occured when compiling fragment shader for: '%.*s'", STR_ARG(identifier));
            return false;
        }

        program_permutations[perm].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, frag_shader, geom_shader};
        if (!md_gl_program_attach_and_link(program_permutations[perm].id, shaders, has_geom ? 3 : 2)) return false;
    }

    glDeleteShader(vert_shader);
    if (geom_shader) glDeleteShader(geom_shader);
    glDeleteShader(frag_shader);

    return true;
}

// Triangulation of one backbone instance, matching the vertex decoding in backbone.vert:
// Tube vertices ring * P + k for ring in [0, S], k in [0, P), then a cap (center, O outline vertices) at the beginning and one at the end.
// Triangles are counter clockwise seen from the outside.
// segments: even, profile_count: even (ellipse), face_subdivisions: >= 1 (box)
static bool create_backbone_mesh(backbone_mesh_t* mesh, uint32_t profile, uint32_t segments, uint32_t profile_count, uint32_t face_subdivisions) {
    ASSERT(mesh);
    ASSERT(segments >= 2 && segments % 2 == 0);
    const uint32_t S = segments;
    uint32_t P, O, K = 0;
    uint32_t edges[64][2];
    uint32_t edge_count = 0;

    if (profile == BACKBONE_PROFILE_ELLIPSE) {
        ASSERT(profile_count >= 4 && profile_count % 2 == 0 && profile_count <= 64);
        P = profile_count;
        O = P;
        for (uint32_t k = 0; k < P; ++k) {
            edges[edge_count][0] = k;
            edges[edge_count][1] = (k + 1) % P;
            ++edge_count;
        }
    } else {
        ASSERT(face_subdivisions >= 1 && face_subdivisions <= 29);
        // Bottom face K + 1 vertices, right face 2 vertices, then their point reflection (top face, left face)
        K = face_subdivisions;
        const uint32_t half = K + 3;
        P = 2 * half;
        O = 2 * (K + 1);
        for (uint32_t h = 0; h < 2; ++h) {
            const uint32_t base = h * half;
            for (uint32_t m = 0; m < K; ++m) {
                edges[edge_count][0] = base + m;
                edges[edge_count][1] = base + m + 1;
                ++edge_count;
            }
            edges[edge_count][0] = base + K + 1;
            edges[edge_count][1] = base + K + 2;
            ++edge_count;
        }
    }

    const uint32_t vertex_count = (S + 1) * P + 2 * (O + 1);
    const uint32_t index_count  = S * edge_count * 6 + 2 * O * 3;
    if (vertex_count > 0xFFFF) {
        MD_LOG_ERROR("Backbone mesh resolution is too high");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    uint16_t* idx = md_temp_alloc(temp, index_count * sizeof(uint16_t));
    uint32_t len = 0;

    for (uint32_t j = 0; j < S; ++j) {
        for (uint32_t e = 0; e < edge_count; ++e) {
            const uint16_t a = (uint16_t)( j      * P + edges[e][0]);
            const uint16_t b = (uint16_t)( j      * P + edges[e][1]);
            const uint16_t c = (uint16_t)((j + 1) * P + edges[e][0]);
            const uint16_t d = (uint16_t)((j + 1) * P + edges[e][1]);
            idx[len++] = a; idx[len++] = b; idx[len++] = c;
            idx[len++] = b; idx[len++] = d; idx[len++] = c;
        }
    }

    // Caps: the beginning faces backwards, the end forwards
    const uint32_t cap_base[2] = { (S + 1) * P, (S + 1) * P + O + 1 };
    for (uint32_t cap = 0; cap < 2; ++cap) {
        const uint16_t center = (uint16_t)cap_base[cap];
        for (uint32_t m = 0; m < O; ++m) {
            const uint16_t r0 = (uint16_t)(cap_base[cap] + 1 + m);
            const uint16_t r1 = (uint16_t)(cap_base[cap] + 1 + (m + 1) % O);
            idx[len++] = center;
            idx[len++] = cap == 0 ? r1 : r0;
            idx[len++] = cap == 0 ? r0 : r1;
        }
    }
    ASSERT(len == index_count);

    gl_buffer_conditional_delete(&mesh->index_buffer);
    mesh->index_buffer = gl_buffer_create(index_count * sizeof(uint16_t), idx, GL_STATIC_DRAW);
    mesh->segments = S;
    mesh->profile_count = P;
    mesh->outline_count = O;
    mesh->face_subdivisions = K;
    mesh->index_count = index_count;

    md_temp_end(temp);
    return true;
}

void md_gl_initialize(void) {
    if (gl3wInit() != GL3W_OK) {
        MD_LOG_ERROR("Could not load OpenGL extensions");
        return;
    }

    GLint major, minor;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);

    if (major < 3 && minor < 3) {
        MD_LOG_ERROR("OpenGL version %i.%i is not supported", major, minor);
        return;
    }

    ctx.version = major * 100 + minor * 10;

    if (!ctx.vao) {
        glGenVertexArrays(1, &ctx.vao);
    }

    glGenFramebuffers(1, &ctx.fbo);
    ctx.ubo = gl_buffer_create(UBO_SIZE, NULL, GL_DYNAMIC_DRAW);

    for (uint32_t i = 0; i < GL_TEXTURE_COUNT; ++i) {
        glGenTextures(1, &ctx.texture[i].id);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);

        if (!md_gl_shader_compile(vert_shader, (str_t){(const char*)compute_velocity_vert, compute_velocity_vert_size}, 0, 0)) {
            return;
        }
        ctx.program[GL_PROGRAM_COMPUTE_VELOCITY].id = glCreateProgram();
        const GLuint shaders[] = { vert_shader };
        const GLchar* varyings[] = { "out_velocity" };
        if (!md_gl_program_attach_and_link_transform_feedback(ctx.program[GL_PROGRAM_COMPUTE_VELOCITY].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) {
            return;
        }

        glDeleteShader(vert_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);

        if (!md_gl_shader_compile(vert_shader, (str_t){(const char*)compute_spline_extract_vert, compute_spline_extract_vert_size}, 0, 0)) {
            return;
        }

        ctx.program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id = glCreateProgram();
        const GLuint shaders[] = { vert_shader };
        const GLchar* varyings[] = { "out_position", "out_atom_idx", "out_velocity", "out_segment_t", "out_secondary_structure_and_flags", "out_support_and_tangent_vector" };
        if (!md_gl_program_attach_and_link_transform_feedback(ctx.program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) {
            return;
        }

        glDeleteShader(vert_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);

        if (!md_gl_shader_compile(vert_shader, (str_t){(const char*)compute_spline_orient_vert, compute_spline_orient_vert_size}, 0, 0)) {
            return;
        }

        ctx.program[GL_PROGRAM_ORIENT_CONTROL_POINTS].id = glCreateProgram();
        const GLuint shaders[] = { vert_shader };
        const GLchar* varyings[] = { "out_position", "out_atom_idx", "out_velocity", "out_segment_t", "out_secondary_structure_and_flags", "out_support_and_tangent_vector" };
        if (!md_gl_program_attach_and_link_transform_feedback(ctx.program[GL_PROGRAM_ORIENT_CONTROL_POINTS].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) {
            return;
        }

        glDeleteShader(vert_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);

        if (!md_gl_shader_compile(vert_shader, (str_t){(const char*)compute_spline_subdivide_vert, compute_spline_subdivide_vert_size}, 0, 0)) {
            return;
        }

        ctx.program[GL_PROGRAM_SUBDIVIDE_SPLINE].id = glCreateProgram();
        const GLuint shaders[] = { vert_shader };
        const GLchar* varyings[] = { "out_ring_0", "out_ring_1", "out_ring_2" };
        if (!md_gl_program_attach_and_link_transform_feedback(ctx.program[GL_PROGRAM_SUBDIVIDE_SPLINE].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) {
            return;
        }

        glDeleteShader(vert_shader);
    }

    {
        STATIC_ASSERT(MD_GL_BACKBONE_SEGMENT_COUNT >= 2 && MD_GL_BACKBONE_SEGMENT_COUNT % 2 == 0, "Invalid backbone segment count");
        STATIC_ASSERT(MD_GL_BACKBONE_PROFILE_COUNT >= 4 && MD_GL_BACKBONE_PROFILE_COUNT % 2 == 0 && MD_GL_BACKBONE_PROFILE_COUNT <= 64, "Invalid cartoon profile count");
        STATIC_ASSERT(BACKBONE_RIBBONS_FACE_SUBDIVISIONS >= 1 && BACKBONE_RIBBONS_FACE_SUBDIVISIONS <= 29, "Invalid ribbons face subdivisions");

        uint32_t segments = MD_GL_BACKBONE_SEGMENT_COUNT;
        uint32_t profile_count = MD_GL_BACKBONE_PROFILE_COUNT;
        uint32_t face_subdivisions = BACKBONE_RIBBONS_FACE_SUBDIVISIONS;
        for (uint32_t lod = 0; lod < BACKBONE_LOD_COUNT; ++lod) {
            if (!create_backbone_mesh(&ctx.backbone_mesh[BACKBONE_PROFILE_ELLIPSE][lod], BACKBONE_PROFILE_ELLIPSE, segments, profile_count, 0) ||
                !create_backbone_mesh(&ctx.backbone_mesh[BACKBONE_PROFILE_BOX][lod],     BACKBONE_PROFILE_BOX,     segments, 0, face_subdivisions)) {
                MD_LOG_ERROR("Failed to create backbone mesh");
                return;
            }
            // Halve the resolution: the segments must be an even divisor of the finest level (they pick its rings)
            uint32_t next = MAX(2, segments / 2);
            while (next > 2 && (next % 2 != 0 || MD_GL_BACKBONE_SEGMENT_COUNT % next != 0)) --next;
            segments = next;
            profile_count = MAX(4, (profile_count / 2) & ~1u);
            face_subdivisions = MAX(1, face_subdivisions / 2);
        }
    }

    if (!ctx.arena) {
        ctx.arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    }

    md_handle_pool_init(&ctx.shader_pool,           MAX_SHADERS,            ctx.arena);
    md_handle_pool_init(&ctx.molecule_pool,         MAX_MOLECULES,          ctx.arena);
    md_handle_pool_init(&ctx.palette_pool,          MAX_PALETTES,           ctx.arena);
    md_handle_pool_init(&ctx.representation_pool,   MAX_REPRESENTATIONS,    ctx.arena);
}

void md_gl_shutdown(void) {
    if (ctx.vao) glDeleteVertexArrays(1, &ctx.vao);
    if (ctx.fbo) glDeleteFramebuffers(1, &ctx.fbo);
    if (ctx.ubo.id) glDeleteBuffers(1, &ctx.ubo.id);
    for (uint32_t i = 0; i < GL_TEXTURE_COUNT; ++i) {
        if (ctx.texture[i].id) glDeleteTextures(1, &ctx.texture[i].id);
    }
    for (uint32_t i = 0; i < GL_PROGRAM_COUNT; ++i) {
        if (ctx.program[i].id) glDeleteProgram(ctx.program[i].id);
    }
    for (uint32_t i = 0; i < BACKBONE_PROFILE_COUNT; ++i) {
        for (uint32_t lod = 0; lod < BACKBONE_LOD_COUNT; ++lod) {
            gl_buffer_conditional_delete(&ctx.backbone_mesh[i][lod].index_buffer);
        }
    }

    md_arena_allocator_destroy(ctx.arena);
}

md_gl_shaders_t md_gl_shaders_create(str_t str) {
    md_gl_shaders_t handle = {0};

    uint32_t id = {md_handle_pool_alloc_slot(&ctx.shader_pool)};
    if (id == MD_HANDLE_INVALID_ID) {
        MD_LOG_ERROR("Fatal error: ran out of slots!");
        return handle;
    }

    int idx = md_handle_index(id);
    ctx.shaders[idx].id = id;
    shaders_t* shaders = ctx.shaders + idx;

    if (str_empty(str)) {
        str = default_shader_output;
    }
    const str_t backbone_vert_src = {(const char*)backbone_vert, backbone_vert_size};
    const str_t backbone_frag_src = {(const char*)backbone_frag, backbone_frag_size};
    if (!create_permuted_program(STR_LIT("SpaceFill"), shaders->spacefill,  (str_t){(const char*)spacefill_vert, spacefill_vert_size}, (str_t){(const char*)spacefill_geom, spacefill_geom_size},   (str_t){(const char*)spacefill_frag, spacefill_frag_size},    str, (str_t){0})) return handle;
    if (!create_permuted_program(STR_LIT("Licorice"),  shaders->licorice,   (str_t){(const char*)licorice_vert, licorice_vert_size},   (str_t){(const char*)licorice_geom, licorice_geom_size},     (str_t){(const char*)licorice_frag, licorice_frag_size},      str, (str_t){0})) return handle;
    if (!create_permuted_program(STR_LIT("Ribbons"),   shaders->ribbons,    backbone_vert_src, (str_t){0}, backbone_frag_src, str, STR_LIT("#define REP_RIBBONS 1"))) return handle;
    if (!create_permuted_program(STR_LIT("Cartoon"),   shaders->cartoon,    backbone_vert_src, (str_t){0}, backbone_frag_src, str, STR_LIT("#define REP_RIBBONS 0"))) return handle;

    handle.id = id;
    return handle;
}

void md_gl_shaders_destroy(md_gl_shaders_t handle) {
    if (handle.id == MD_HANDLE_INVALID_ID) return;
    shaders_t* shaders = shad_lookup(handle.id);
    if (shaders) {
        for (int i = 0; i < MAX_SHADER_PERMUTATIONS; ++i) {
            if (glIsProgram(shaders->spacefill[i].id))  glDeleteProgram(shaders->spacefill[i].id);
            if (glIsProgram(shaders->licorice[i].id))   glDeleteProgram(shaders->licorice[i].id);
            if (glIsProgram(shaders->ribbons[i].id))    glDeleteProgram(shaders->ribbons[i].id);
            if (glIsProgram(shaders->cartoon[i].id))    glDeleteProgram(shaders->cartoon[i].id);
        }
        md_handle_pool_free_slot(&ctx.shader_pool, handle.id);
    }
}

md_gl_mol_t md_gl_mol_create(const md_system_t* sys) {
    md_gl_mol_t handle = {0};
    md_temp_scope_t temp_scope = md_temp_begin();

    if (sys) {
        if (sys->atom.count == 0) {
            MD_LOG_ERROR("The supplied molecule has no atoms.");
            goto done;
        }

        uint32_t id = {md_handle_pool_alloc_slot(&ctx.molecule_pool)};
        if (id == MD_HANDLE_INVALID_ID) {
            MD_LOG_ERROR("Fatal error: ran out of slots!");
            goto done;
        }

        int index = md_handle_index(id);
        ctx.molecules[index].id = id;
        handle.id = id;
        molecule_t* gl_mol = ctx.molecules + index;
        gl_mol->backbone_orient_history = false;

        gl_mol->atom_count = (uint32_t)sys->atom.count;
        gl_mol->buffer[GL_BUFFER_ATOM_POSITION]        = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 3,   NULL, GL_DYNAMIC_DRAW);
        gl_mol->buffer[GL_BUFFER_ATOM_POSITION_PREV]   = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 3,   NULL, GL_DYNAMIC_COPY);
        gl_mol->buffer[GL_BUFFER_ATOM_VELOCITY]        = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 3,   NULL, GL_DYNAMIC_COPY);
        gl_mol->buffer[GL_BUFFER_ATOM_RADIUS]          = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 1,   NULL, GL_STATIC_DRAW);
        gl_mol->buffer[GL_BUFFER_ATOM_FLAGS]           = gl_buffer_create(gl_mol->atom_count * sizeof(uint8_t),     NULL, GL_STATIC_DRAW);

        // Seed with the reference configuration; callers push per frame positions afterwards.
        if (md_system_state_has_coords(&sys->reference)) {
            md_gl_mol_set_atom_position(handle, 0, gl_mol->atom_count, sys->reference.xyz);
        }
        md_gl_mol_zero_velocity(handle);

        float* radii = md_temp_alloc(temp_scope, sizeof(float) * sys->atom.count);
        md_atom_extract_radii(radii, 0, sys->atom.count, &sys->atom);

        md_gl_mol_set_atom_radius(handle, 0, gl_mol->atom_count, radii, 0);
        //if (mol->atom.flags)  md_gl_molecule_set_atom_flags(ext_mol,  0, gl_mol->atom_count, mol->atom.flags, 0);

        gl_mol->comp_count = (uint32_t)sys->component.count;
        //gl_mol->buffer[GL_BUFFER_RESIDUE_ATOM_RANGE]          = gl_buffer_create(gl_mol->comp_count * sizeof(md_range_t),   NULL, GL_STATIC_DRAW);
        //gl_mol->buffer[GL_BUFFER_RESIDUE_AABB]                = gl_buffer_create(gl_mol->comp_count * sizeof(float) * 6,    NULL, GL_DYNAMIC_COPY);
        //gl_mol->buffer[GL_BUFFER_RESIDUE_VISIBLE]             = gl_buffer_create(gl_mol->comp_count * sizeof(int),          NULL, GL_DYNAMIC_COPY);

        //if (sys->comp.atom_offset)           gl_buffer_set_sub_data(gl_mol->buffer[GL_BUFFER_RESIDUE_ATOM_RANGE], 0, gl_mol->comp_count * sizeof(uint32_t) * 2, sys->comp.atom_range);
        //if (desc->residue.backbone_atoms)       gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_RESIDUE_BACKBONE_ATOMS], 0, sys->comp_count * sizeof(uint8_t) * 4, desc->residue.backbone_atoms);

        if (sys->protein_backbone.range.count > 0 && sys->protein_backbone.range.offset && sys->protein_backbone.segment.atoms) {
            uint32_t backbone_residue_count = 0;
            for (uint32_t i = 0; i < (uint32_t)sys->protein_backbone.range.count; ++i) {
                uint32_t res_count = sys->protein_backbone.range.offset[i+1] - sys->protein_backbone.range.offset[i];
                backbone_residue_count += res_count;
            }

            const uint32_t backbone_count                     = backbone_residue_count;
            const uint32_t backbone_control_point_data_count  = backbone_residue_count;
            const uint32_t backbone_control_point_index_count = backbone_residue_count + (uint32_t)sys->protein_backbone.range.count * (2 + 1); // Duplicate pair first and last in each chain for adjacency + primitive restart between

            gl_mol->buffer[GL_BUFFER_BACKBONE_DATA]                = gl_buffer_create(backbone_count                     * sizeof(gl_backbone_data_t),         NULL, GL_STATIC_DRAW);
            gl_mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE] = gl_buffer_create(backbone_count                     * sizeof(md_secondary_structure_t),   NULL, GL_DYNAMIC_DRAW);
            gl_mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA]  = gl_buffer_create(backbone_control_point_data_count  * sizeof(gl_control_point_t),         NULL, GL_DYNAMIC_COPY);
            gl_mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT]      = gl_buffer_create(backbone_control_point_data_count * sizeof(gl_control_point_t), NULL, GL_DYNAMIC_COPY);
            gl_mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT_PREV] = gl_buffer_create(backbone_control_point_data_count * sizeof(gl_control_point_t), NULL, GL_DYNAMIC_COPY);
            gl_mol->buffer[GL_BUFFER_BACKBONE_NEIGHBOR]            = gl_buffer_create(backbone_control_point_data_count  * sizeof(gl_neighbor_idx_t),         NULL, GL_STATIC_DRAW);
            gl_mol->buffer[GL_BUFFER_BACKBONE_RING_DATA]           = gl_buffer_create(backbone_control_point_data_count  * MD_GL_BACKBONE_SEGMENT_COUNT * sizeof(gl_backbone_ring_t), NULL, GL_DYNAMIC_COPY);
            gl_mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX] = gl_buffer_create(backbone_control_point_index_count * sizeof(uint32_t),                   NULL, GL_STATIC_DRAW);

            //gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE], 0, desc->backbone.count * sizeof(uint8_t) * 4, desc->backbone.secondary_structure);

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_BACKBONE_DATA].id);
            gl_backbone_data_t* backbone_data = (gl_backbone_data_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (backbone_data) {
                uint32_t idx = 0;
                for (size_t i = 0; i < sys->protein_backbone.range.count; ++i) {
                    uint32_t beg = sys->protein_backbone.range.offset[i];
                    uint32_t end = sys->protein_backbone.range.offset[i+1];
                    for (uint32_t j = beg; j < end; ++j) {
                        md_component_idx_t comp_idx = sys->protein_backbone.segment.comp_idx[j];
                        uint32_t comp_atom_offset = sys->component.atom_offset[comp_idx];

                        backbone_data[idx].bb_seg_idx = j;
                        backbone_data[idx].atom_base_idx = comp_atom_offset;
                        backbone_data[idx].ca_idx = (uint8_t)(sys->protein_backbone.segment.atoms[j].ca - comp_atom_offset);
                        backbone_data[idx].c_idx  = (uint8_t)(sys->protein_backbone.segment.atoms[j].c  - comp_atom_offset);
                        backbone_data[idx].o_idx  = (uint8_t)(sys->protein_backbone.segment.atoms[j].o  - comp_atom_offset);
                        backbone_data[idx].flags  = (uint8_t)((j == beg ? 1 : 0) | (j == end - 1 ? 2 : 0));
                        ++idx;
                    }
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                goto done;
            }

            // The secondary structure depends on the frame, and is not the system's: coil until the caller sets it
            // (md_gl_mol_set_backbone_secondary_structure)
            const uint32_t ss_coil = pack_gl_secondary_structure((md_gl_secondary_structure_t){0});

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE].id);
            uint32_t* secondary_structure = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (secondary_structure) {
                for (size_t i = 0; i < sys->protein_backbone.segment.count; ++i) {
                    secondary_structure[i] = ss_coil;
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                goto done;
            }

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX].id);
            uint32_t* control_point_index = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (control_point_index) {
                uint32_t len = 0;
                for (size_t i = 0; i < sys->protein_backbone.range.count; ++i) {
                    uint32_t bb_beg = sys->protein_backbone.range.offset[i];
                    uint32_t bb_end = sys->protein_backbone.range.offset[i+1];
                    control_point_index[len++] = bb_beg;
                    for (uint32_t j = bb_beg; j < bb_end; ++j) {
                        control_point_index[len++] = j;
                    }
                    control_point_index[len++] = bb_end-1;
                    control_point_index[len++] = 0xFFFFFFFFU;
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            }
            else {
                goto done;
            }

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_BACKBONE_NEIGHBOR].id);
            gl_neighbor_idx_t* neighbor = (gl_neighbor_idx_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (neighbor) {
                for (uint32_t i = 0; i < backbone_control_point_data_count; ++i) {
                    neighbor[i].idx[0] = i;
                    neighbor[i].idx[1] = i;
                    neighbor[i].idx[2] = i;
                    neighbor[i].idx[3] = i;
                }

                for (uint32_t i = 0; i < (uint32_t)sys->protein_backbone.range.count; ++i) {
                    uint32_t beg = sys->protein_backbone.range.offset[i];
                    uint32_t end = sys->protein_backbone.range.offset[i + 1];
                    for (uint32_t j = beg; j < end; ++j) {
                        const uint32_t p1 = (j > beg) ? (j - 1) : j;
                        const uint32_t n1 = (j + 1 < end) ? (j + 1) : j;
                        neighbor[j].idx[0] = p1;
                        neighbor[j].idx[1] = n1;
                        neighbor[j].idx[2] = (p1 > beg) ? (p1 - 1) : p1;
                        neighbor[j].idx[3] = (n1 + 1 < end) ? (n1 + 1) : n1;
                    }
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                goto done;
            }

            glBindBuffer(GL_ARRAY_BUFFER, 0);

            gl_mol->backbone_control_point_index_count = backbone_control_point_index_count;
            gl_mol->backbone_count = backbone_count;
            gl_mol->flags |= MOL_FLAG_HAS_BACKBONE;

            {
                backbone_bounds_t* bounds = &gl_mol->backbone_bounds;
                md_allocator_i* alloc = md_get_heap_allocator();
                const uint32_t chain_count = (uint32_t)sys->protein_backbone.range.count;
                bounds->chain_count  = chain_count;
                bounds->chain_offset = md_alloc(alloc, (chain_count + 1) * sizeof(uint32_t));
                bounds->chain_sphere = md_alloc(alloc, chain_count * sizeof(vec4_t));
                bounds->ca_atom_idx  = md_alloc(alloc, backbone_count * sizeof(uint32_t));
                bounds->ca_xyz       = md_alloc(alloc, backbone_count * sizeof(vec3_t));
                MEMSET(bounds->chain_sphere, 0, chain_count * sizeof(vec4_t));
                MEMSET(bounds->ca_xyz, 0, backbone_count * sizeof(vec3_t));
                for (uint32_t i = 0; i <= chain_count; ++i) {
                    bounds->chain_offset[i] = sys->protein_backbone.range.offset[i];
                }
                for (uint32_t i = 0; i < backbone_count; ++i) {
                    bounds->ca_atom_idx[i] = (uint32_t)sys->protein_backbone.segment.atoms[i].ca;
                }
                if (md_system_state_has_coords(&sys->reference)) {
                    backbone_bounds_update(bounds, 0, gl_mol->atom_count, sys->reference.xyz);
                }
            }
        }

        gl_mol->bond_count = (uint32_t)sys->bond.count;
        gl_mol->buffer[GL_BUFFER_BOND_ATOM_INDICES] = gl_buffer_create(gl_mol->bond_count * sizeof(uint32_t) * 2, NULL, GL_DYNAMIC_COPY);

        if (sys->bond.pairs) {
            md_gl_mol_set_bonds(handle, 0, gl_mol->bond_count, sys->bond.pairs, sizeof(md_atom_pair_t));
        }
    }

done:
    md_temp_end(temp_scope);
    return handle;
}

void md_gl_mol_destroy(md_gl_mol_t handle) {
    if (handle.id == MD_HANDLE_INVALID_ID) return;
    molecule_t* mol = mol_lookup(handle.id);
    if (mol) {
        for (uint32_t i = 0; i < GL_BUFFER_COUNT; ++i) {
            gl_buffer_conditional_delete(&mol->buffer[i]);
        }
        backbone_bounds_free(&mol->backbone_bounds);
        MEMSET(mol, 0, sizeof(molecule_t));
        md_handle_pool_free_slot(&ctx.molecule_pool, handle.id);
    }
}

md_gl_rep_t md_gl_rep_create(md_gl_mol_t mol_handle) {
    md_gl_rep_t handle = {0};
    molecule_t* mol = mol_lookup(mol_handle.id);
    if (mol) {
        if (mol->atom_count == 0) {
            MD_LOG_ERROR("Supplied molecule has no atoms.");
            return handle;
        }

        uint32_t id = {md_handle_pool_alloc_slot(&ctx.representation_pool)};
        if (id == MD_HANDLE_INVALID_ID) {
            MD_LOG_ERROR("Fatal error: ran out of slots!");
            return handle;
        }

        int idx = md_handle_index(id);
        handle.id = id;
        representation_t* rep = ctx.representations + idx;

        rep->id = id;
        rep->mol_id = mol_handle.id;
        rep->pal_id = 0;
        rep->atom_color = gl_buffer_create(mol->atom_count * sizeof(uint32_t), NULL, GL_STATIC_DRAW);
        const uint32_t color = (uint32_t)((255 << 24) | (127 << 16) | (127 << 8) | (127 << 0));
        glBindBuffer(GL_ARRAY_BUFFER, rep->atom_color.id);
        uint32_t* data = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (data) {
            for (uint32_t i = 0; i < mol->atom_count; ++i) {
                data[i] = color;
            }
            glUnmapBuffer(GL_ARRAY_BUFFER);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            //md_update_visible_atom_color_range(rep);
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    return handle;
}

void md_gl_rep_destroy(md_gl_rep_t handle) {
    if (handle.id == MD_HANDLE_INVALID_ID) return;
    representation_t* rep = rep_lookup(handle.id);
    if (rep) {
        gl_buffer_conditional_delete(&rep->atom_color);
        MEMSET(rep, 0, sizeof(representation_t));
        md_handle_pool_free_slot(&ctx.representation_pool, handle.id);
    }
}

void md_gl_rep_set_atom_colors(md_gl_rep_t handle, uint32_t offset, uint32_t count, const uint32_t* color_data, uint32_t byte_stride) {
    if (color_data == NULL) {
        MD_LOG_ERROR("color_data ptr was NULL");
        return;
    }
    representation_t* rep = rep_lookup(handle.id);
    if (rep == NULL) {
        MD_LOG_ERROR("representation was invalid");
        return;
    }
    molecule_t* mol = mol_lookup(rep->mol_id);
    if (mol == NULL) {
        MD_LOG_ERROR("Representation's molecule is invalid");
        return;
    }
    if (offset + count > mol->atom_count) {
        MD_LOG_ERROR("Attempting to write out of bounds");
        return;
    }
    ASSERT(glIsBuffer(rep->atom_color.id));
    if (byte_stride) {
        glBindBuffer(GL_ARRAY_BUFFER, rep->atom_color.id);
        uint32_t* data = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (data) {
            for (uint32_t i = offset; i < offset + count; ++i) {
                data[i] = *(uint32_t*)((uint8_t*)color_data + i * byte_stride);
            }
            glUnmapBuffer(GL_ARRAY_BUFFER);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    else {
        gl_buffer_set_sub_data(rep->atom_color, offset * sizeof(uint32_t), count * sizeof(uint32_t), color_data);
    }
    //md_update_visible_atom_color_range(rep);
}

static bool compute_spline(molecule_t* mol);

static bool draw_space_fill(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, float scale);
static bool draw_licorice  (gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, float radius, float max_length, md_gl_bond_mode_t mode, float sharpness, uint32_t uniform_color);
static bool draw_backbone  (gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, uint32_t profile, const float scale[4], float max_extent, const mat4_t* world_to_clip, float viewport_half_height);

static inline void init_ubo_base_data(gl_ubo_base_t* ubo_data, const md_gl_draw_args_t* args, const mat4_t* model_matrix) {
    ASSERT(ubo_data);
    ASSERT(args);

    MEMCPY(&ubo_data->view_transform.world_to_view, args->view_transform.view_matrix, sizeof(mat4_t));
    if (model_matrix) {
        ubo_data->view_transform.world_to_view = mat4_mul(ubo_data->view_transform.world_to_view, *model_matrix);
    }
    MEMCPY(&ubo_data->view_transform.view_to_clip,  args->view_transform.proj_matrix, sizeof(mat4_t));

    ubo_data->view_transform.world_to_clip = mat4_mul(ubo_data->view_transform.view_to_clip, ubo_data->view_transform.world_to_view);
    ubo_data->view_transform.world_to_view_normal = mat4_transpose(mat4_inverse(ubo_data->view_transform.world_to_view));
    ubo_data->view_transform.view_to_world = mat4_inverse(ubo_data->view_transform.world_to_view);
    ubo_data->view_transform.clip_to_view = mat4_inverse(ubo_data->view_transform.view_to_clip);

    if (args->view_transform.prev_view_matrix && args->view_transform.proj_matrix) {
        const mat4_t* prev_world_to_view = (const mat4_t*)args->view_transform.prev_view_matrix;
        const mat4_t* prev_view_to_clip  = (const mat4_t*)args->view_transform.prev_proj_matrix;
        ubo_data->view_transform.prev_world_to_clip = mat4_mul(*prev_view_to_clip, *prev_world_to_view);
        ubo_data->view_transform.curr_view_to_prev_clip = mat4_mul(ubo_data->view_transform.prev_world_to_clip, ubo_data->view_transform.view_to_world);
        extract_jitter_uv(ubo_data->jitter_uv.elem + 0, ubo_data->view_transform.view_to_clip);
        extract_jitter_uv(ubo_data->jitter_uv.elem + 2, *prev_view_to_clip);
    }
    ubo_data->atom_mask = args->atom_mask;
    ubo_data->atom_index_base = args->picking_offset.atom_base;
    ubo_data->bond_index_base = args->picking_offset.bond_base;
}

static inline bool is_backbone_representation_type(md_gl_rep_type_t type) {
    return type == MD_GL_REP_RIBBONS || type == MD_GL_REP_CARTOON;
}

bool md_gl_draw(const md_gl_draw_args_t* args) {
    if (!args) {
        MD_LOG_ERROR("draw args object was NULL");
        return false;
    }
    if (!validate_context()) return false;

    shaders_t* shaders = shad_lookup(args->shaders.id);
    if (!shaders) {
        MD_LOG_ERROR("shaders object was invalid");
        return false;
    }

    PUSH_GPU_SECTION("MOLD DRAW")
            
    gl_ubo_base_t ubo_data = {0};
    init_ubo_base_data(&ubo_data, args, NULL);

    gl_buffer_set_sub_data(ctx.ubo, 0, sizeof(ubo_data), &ubo_data);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ctx.ubo.id);

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    bool result = false;

    // Valid draw operations to issue
    md_array(md_gl_draw_op_t const*) draw_ops = 0;
    md_array(molecule_t*) bb_mols = 0;
        
    // Validate and extract backbone molecules
    for (size_t i = 0; i < args->draw_operations.count; ++i) {
        const md_gl_draw_op_t* draw_op = &args->draw_operations.ops[i];
        const representation_t* rep = rep_lookup(draw_op->rep.id);
        if (!rep) {
            MD_LOG_ERROR("Invalid representation");
            continue;
        }
        molecule_t* mol = mol_lookup(rep->mol_id);
        if (!mol) {
            MD_LOG_ERROR("Invalid molecule");
            continue;
        }
        bool bb_rep = is_backbone_representation_type(draw_op->type);
        
        if (bb_rep && !(mol->flags & MOL_FLAG_HAS_BACKBONE)) {
            MD_LOG_ERROR("Incompatible representation type for molecule: Molecule is missing a protein_backbone");
            continue;
        }

        md_array_push(draw_ops, draw_op, alloc);

        if (bb_rep) {
            int64_t mol_idx = -1;
            for (size_t j = 0; j < md_array_size(bb_mols); j++) {
                if (bb_mols[j] == mol) {
                    mol_idx = (int64_t)j;
                    break;
                }
            }
            // Not found
            if (mol_idx == -1) {
                md_array_push(bb_mols, mol, alloc);
            }
        }
    }
    
    //qsort((void*)draw_ent, draw_ent_count, sizeof(draw_entity_t), compare_draw_ent);
            
    PUSH_GPU_SECTION("COMPUTE SPLINE")
    for (size_t i = 0; i < md_array_size(bb_mols); i++) {
        compute_spline(bb_mols[i]);
    }
    POP_GPU_SECTION()
        
    //bool using_internal_depth = false;

    int program_permutation = 0;
    if (is_ortho_proj_matrix(ubo_data.view_transform.view_to_clip)) program_permutation |= PERMUTATION_BIT_ORTHO;

    // Maximum bond length in units (Ångström assumed)
    const float max_length = args->max_bond_length > 0.0f ? args->max_bond_length : 5.0f;

    // For the level of detail of the backbone representations
    GLint viewport[4] = {0};
    glGetIntegerv(GL_VIEWPORT, viewport);
    const float viewport_half_height = (float)viewport[3] * 0.5f;
        
    PUSH_GPU_SECTION("DRAW REPRESENTATIONS")
    for (size_t i = 0; i < md_array_size(draw_ops); i++) {
        const md_gl_draw_op_t* draw_op = draw_ops[i];
        const representation_t* rep = rep_lookup(draw_op->rep.id);
        const molecule_t* mol       = mol_lookup(rep->mol_id); 
        const mat4_t* model_matrix  = (const mat4_t*)draw_op->model_matrix;
        float scale = 1.0f;
        mat4_t world_to_clip = ubo_data.view_transform.world_to_clip;

        if (model_matrix) {
            // If we have a model matrix, we need to recompute the entire matrix stack...
            gl_ubo_base_t ubo_tmp = {0};
            init_ubo_base_data(&ubo_tmp, args, model_matrix);
            gl_buffer_set_sub_data(ctx.ubo, 0, sizeof(gl_view_transform_t), &ubo_tmp);
            world_to_clip = ubo_tmp.view_transform.world_to_clip;
            const vec3_t model_scale = {
                vec3_length(vec3_from_vec4(ubo_tmp.view_transform.world_to_view.col[0])),
                vec3_length(vec3_from_vec4(ubo_tmp.view_transform.world_to_view.col[1])),
                vec3_length(vec3_from_vec4(ubo_tmp.view_transform.world_to_view.col[2])),
            };
            const float mean = (model_scale.x + model_scale.y + model_scale.z) / 3.0f;

            // Non uniform scale is not supported, it will cause rendering artifact for radius scaling parameters which are assumed to be uniform.
            scale = mean;
        }

        switch (draw_op->type) {
        case MD_GL_REP_SPACE_FILL:
            draw_space_fill(shaders->spacefill[program_permutation], mol, rep->atom_color, scale * draw_op->args.space_fill.radius_scale);
            break;
        case MD_GL_REP_LICORICE:
            draw_licorice(shaders->licorice[program_permutation],    mol, rep->atom_color, 0.2f * scale * draw_op->args.licorice.radius, max_length, draw_op->args.licorice.color_mode, draw_op->args.licorice.sharpness, draw_op->args.licorice.uniform_color);
            break;
        case MD_GL_REP_BALL_AND_STICK:
            draw_licorice(shaders->licorice[program_permutation],    mol, rep->atom_color, 0.2f * scale * draw_op->args.ball_and_stick.stick_radius, max_length, draw_op->args.ball_and_stick.color_mode, draw_op->args.ball_and_stick.sharpness, draw_op->args.ball_and_stick.uniform_color);
            draw_space_fill(shaders->spacefill[program_permutation], mol, rep->atom_color, 0.2f * scale * draw_op->args.ball_and_stick.ball_scale);
            break;
        case MD_GL_REP_RIBBONS: {
            // Half width and half thickness of the box profile
            const float profile_scale[4] = { scale * draw_op->args.ribbons.width_scale, scale * draw_op->args.ribbons.thickness_scale * 0.1f, 0.0f, 0.0f };
            const float max_extent = sqrtf(profile_scale[0] * profile_scale[0] + profile_scale[1] * profile_scale[1]);
            draw_backbone(shaders->ribbons[program_permutation], mol, rep->atom_color, BACKBONE_PROFILE_BOX, profile_scale, max_extent, &world_to_clip, viewport_half_height);
            break;
        }
        case MD_GL_REP_CARTOON: {
            // Scales of the coil, helix and sheet profiles
            const float profile_scale[4] = { scale * draw_op->args.cartoon.coil_scale, draw_op->args.cartoon.helix_scale, draw_op->args.cartoon.sheet_scale, 0.0f };
            // Largest semi axis of the profiles in backbone.vert
            const float max_extent = MAX(MAX(0.2f * profile_scale[0], 1.2f * profile_scale[1]), 1.5f * profile_scale[2]);
            draw_backbone(shaders->cartoon[program_permutation], mol, rep->atom_color, BACKBONE_PROFILE_ELLIPSE, profile_scale, max_extent, &world_to_clip, viewport_half_height);
            break;
        }
        default:
            MD_LOG_ERROR("Representation had unexpected type");
            goto done;
        }

        if (model_matrix) {
            // Reset matrix stack
            gl_buffer_set_sub_data(ctx.ubo, 0, sizeof(gl_view_transform_t), &ubo_data);
        }
    }

    result = true;

done:
    POP_GPU_SECTION()

    POP_GPU_SECTION()

    md_temp_end(temp);
    
    return result;
}

static bool draw_space_fill(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, float scale) {
    ASSERT(mol);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_RADIUS].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(atom_color.id);

    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(scale), &scale);

    glBindVertexArray(ctx.vao);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_RADIUS].id);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_BYTE, 0, 0);

    glEnableVertexAttribArray(4);
    glBindBuffer(GL_ARRAY_BUFFER, atom_color.id);
    glVertexAttribPointer(4, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
    
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glUseProgram(program.id);
    glDrawArrays(GL_POINTS, 0, mol->atom_count);
    glUseProgram(0);
    
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);
    
    glBindVertexArray(0);
    
    return true;
}

static bool draw_licorice(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, float radius, float max_length, md_gl_bond_mode_t mode, float sharpness, uint32_t uniform_color) {
    ASSERT(mol);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(mol->buffer[GL_BUFFER_BOND_ATOM_INDICES].id);
    ASSERT(atom_color.id);

    if (max_length == 0) {
        max_length = 1000.0f;
    }

    struct {
		float radius;
		float max_d2;
        int   mode;
        float sharpness;
		uint32_t uniform_color;
	} params = { radius, max_length * max_length, mode, sharpness, uniform_color };

    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(params), &params);

    glBindVertexArray(ctx.vao);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    glVertexAttribIPointer(2, 1, GL_UNSIGNED_BYTE, 0, 0);

    glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, atom_color.id);
    glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mol->buffer[GL_BUFFER_BOND_ATOM_INDICES].id);

    glUseProgram(program.id);
    glDrawElements(GL_LINES, mol->bond_count * 2, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);

    glBindVertexArray(0);

    return true;
}

enum {
    BACKBONE_CHAIN_CULLED = 0xFF
};

// Level of detail of each chain, or BACKBONE_CHAIN_CULLED: the chains are culled against the view frustum and get the
// coarsest level whose segments are at most BACKBONE_LOD_SEGMENT_PIXELS long where the chain is closest to the camera.
static void backbone_chain_lod(uint8_t* chain_lod, const backbone_bounds_t* bounds, const backbone_mesh_t meshes[BACKBONE_LOD_COUNT], float max_extent, const mat4_t* world_to_clip, float viewport_half_height) {
    // Rows of the matrix, the clip planes are combinations of them
    vec4_t row[4];
    for (int r = 0; r < 4; ++r) {
        row[r] = (vec4_t){world_to_clip->elem[0][r], world_to_clip->elem[1][r], world_to_clip->elem[2][r], world_to_clip->elem[3][r]};
    }
    vec4_t plane[6] = {
        vec4_add(row[3], row[0]), vec4_sub(row[3], row[0]),
        vec4_add(row[3], row[1]), vec4_sub(row[3], row[1]),
        vec4_add(row[3], row[2]), vec4_sub(row[3], row[2]),
    };
    for (int i = 0; i < 6; ++i) {
        const float len = vec3_length(vec3_from_vec4(plane[i]));
        plane[i] = len > 0.0f ? vec4_mul1(plane[i], 1.0f / len) : (vec4_t){0, 0, 0, 1};
    }
    const float w_scale = vec3_length(vec3_from_vec4(row[3]));
    // Pixels per unit of length at clip w = 1
    const float pixel_scale = vec3_length(vec3_from_vec4(row[1])) * viewport_half_height;

    for (uint32_t c = 0; c < bounds->chain_count; ++c) {
        const uint32_t cp_count = bounds->chain_offset[c + 1] - bounds->chain_offset[c];
        if (cp_count < 2) {
            chain_lod[c] = BACKBONE_CHAIN_CULLED;
            continue;
        }
        chain_lod[c] = 0;
        if (!bounds->valid) {
            continue;
        }

        const vec4_t sphere = bounds->chain_sphere[c];
        const vec4_t center = {sphere.x, sphere.y, sphere.z, 1.0f};
        const float radius = sphere.w + max_extent;

        bool culled = false;
        for (int i = 0; i < 6; ++i) {
            if (vec4_dot(plane[i], center) < -radius) {
                culled = true;
                break;
            }
        }
        if (culled) {
            chain_lod[c] = BACKBONE_CHAIN_CULLED;
            continue;
        }

        const float w_near = vec4_dot(row[3], center) - radius * w_scale;
        if (w_near <= 0.0f || pixel_scale <= 0.0f) {
            continue;
        }
        const float pixels_per_unit = pixel_scale / w_near;
        for (uint32_t lod = BACKBONE_LOD_COUNT - 1; lod > 0; --lod) {
            if (BACKBONE_RESIDUE_LENGTH / (float)meshes[lod].segments * pixels_per_unit <= BACKBONE_LOD_SEGMENT_PIXELS) {
                chain_lod[c] = (uint8_t)lod;
                break;
            }
        }
    }
}

// One instance per control point (residue), see backbone.vert. Consecutive chains with the same level of detail are drawn together.
static bool draw_backbone(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, uint32_t profile, const float scale[4], float max_extent, const mat4_t* world_to_clip, float viewport_half_height) {
    ASSERT(mol);
    ASSERT(profile < BACKBONE_PROFILE_COUNT);
    ASSERT(world_to_clip);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_NEIGHBOR].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_RING_DATA].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(atom_color.id);

    const backbone_mesh_t* meshes = ctx.backbone_mesh[profile];
    const backbone_bounds_t* bounds = &mol->backbone_bounds;
    if (!meshes[0].index_buffer.id || mol->backbone_count == 0 || bounds->chain_count == 0) {
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    uint8_t* chain_lod = md_temp_alloc(temp, bounds->chain_count);
    backbone_chain_lod(chain_lod, bounds, meshes, max_extent, world_to_clip, viewport_half_height);

    glBindVertexArray(ctx.vao);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color.id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, mol->buffer[GL_BUFFER_ATOM_FLAGS].id);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_2].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT].id);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_3].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, mol->buffer[GL_BUFFER_BACKBONE_NEIGHBOR].id);

    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_4].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, mol->buffer[GL_BUFFER_BACKBONE_RING_DATA].id);

    glUseProgram(program.id);
    glUniform1i(glGetUniformLocation(program.id, "u_atom_color_buffer"), 0);
    glUniform1i(glGetUniformLocation(program.id, "u_atom_flags_buffer"), 1);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_control_points"), 2);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_neighbors"), 3);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_rings"), 4);
    const GLint instance_offset_loc = glGetUniformLocation(program.id, "u_instance_offset");

    for (uint32_t lod = 0; lod < BACKBONE_LOD_COUNT; ++lod) {
        const backbone_mesh_t* mesh = &meshes[lod];
        bool bound = false;

        uint32_t c = 0;
        while (c < bounds->chain_count) {
            if (chain_lod[c] != lod) {
                ++c;
                continue;
            }
            const uint32_t beg = bounds->chain_offset[c];
            while (c < bounds->chain_count && chain_lod[c] == lod) ++c;
            const uint32_t end = bounds->chain_offset[c];

            if (!bound) {
                struct {
                    float    scale[4];
                    uint32_t res[4];
                    uint32_t rings[4];
                } params = {
                    {scale[0], scale[1], scale[2], scale[3]},
                    {mesh->segments, mesh->profile_count, mesh->outline_count, mesh->face_subdivisions},
                    {MD_GL_BACKBONE_SEGMENT_COUNT, MD_GL_BACKBONE_SEGMENT_COUNT / mesh->segments, 0, 0},
                };
                gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(params), &params);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->index_buffer.id);
                bound = true;
            }

            glUniform1ui(instance_offset_loc, beg);
            glDrawElementsInstanced(GL_TRIANGLES, (GLsizei)mesh->index_count, GL_UNSIGNED_SHORT, 0, (GLsizei)(end - beg));
        }
    }

    glUseProgram(0);

    glActiveTexture(GL_TEXTURE0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    md_temp_end(temp);
    return true;
}

static bool compute_spline(molecule_t* mol) {
    ASSERT(ctx.program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id);
    ASSERT(ctx.program[GL_PROGRAM_ORIENT_CONTROL_POINTS].id);
    ASSERT(ctx.program[GL_PROGRAM_SUBDIVIDE_SPLINE].id);
    ASSERT(mol);

    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_DATA].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT_PREV].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_NEIGHBOR].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_RING_DATA].id);

    if (mol->buffer[GL_BUFFER_BACKBONE_DATA].id == 0) {
        MD_LOG_ERROR("Backbone data buffer is zero, which is required to compute the protein_backbone. Is the molecule missing a protein_backbone?");
        return false;
    }

    glEnable(GL_RASTERIZER_DISCARD);
    glBindVertexArray(ctx.vao);

    // Pass 1: Extract control points
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_BACKBONE_DATA].id);

    glEnableVertexAttribArray(0);
    glVertexAttribIPointer(0, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, bb_seg_idx));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, atom_base_idx));

    glEnableVertexAttribArray(2);
    glVertexAttribIPointer(2, 1, GL_UNSIGNED_BYTE, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, ca_idx));

    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_BYTE, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, c_idx));

    glEnableVertexAttribArray(4);
    glVertexAttribIPointer(4, 1, GL_UNSIGNED_BYTE, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, o_idx));

    glEnableVertexAttribArray(5);
    glVertexAttribIPointer(5, 1, GL_UNSIGNED_BYTE, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, flags));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_POSITION].id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_2].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE].id);

    {
        GLuint program = ctx.program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id;
        const GLint buf_atom_pos_loc            = glGetUniformLocation(program, "u_buf_atom_pos");
        const GLint buf_atom_vel_loc            = glGetUniformLocation(program, "u_buf_atom_vel");
        const GLint buf_secondary_structure_loc = glGetUniformLocation(program, "u_buf_secondary_structure");

        glUseProgram(program);
        glUniform1i(buf_atom_pos_loc, 0);
        glUniform1i(buf_atom_vel_loc, 1);
        glUniform1i(buf_secondary_structure_loc, 2);
        glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA].id);
        glBeginTransformFeedback(GL_POINTS);
        glDrawArrays(GL_POINTS, 0, mol->backbone_count);
        glEndTransformFeedback();
        glUseProgram(0);
    }

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);
    glDisableVertexAttribArray(5);

    // Pass 2: Orient control points
    // The last oriented control points become the previous ones, which the orient pass continues from
    // (the relation between neighbouring support vectors is kept with hysteresis).
    {
        gl_buffer_t tmp = mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT];
        mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT] = mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT_PREV];
        mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT_PREV] = tmp;
    }

    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA].id);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, position));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, atom_idx));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, velocity));

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, segment_t));

    glEnableVertexAttribArray(4);
    glVertexAttribIPointer(4, 1, GL_UNSIGNED_INT, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, secondary_structure));

    glEnableVertexAttribArray(5);
    glVertexAttribIPointer(5, 3, GL_UNSIGNED_INT, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, support_vector));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA].id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, mol->buffer[GL_BUFFER_BACKBONE_NEIGHBOR].id);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_2].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT_PREV].id);

    {
        // A relation is kept until the angle between neighbouring support vectors under it exceeds
        // 90 degrees + this margin.
        const float hysteresis_deg = 40.0f;
        GLuint program = ctx.program[GL_PROGRAM_ORIENT_CONTROL_POINTS].id;
        glUseProgram(program);
        glUniform1i(glGetUniformLocation(program, "u_buf_control_point_words"), 0);
        glUniform1i(glGetUniformLocation(program, "u_buf_neighbors"), 1);
        glUniform1i(glGetUniformLocation(program, "u_buf_prev_control_point_words"), 2);
        glUniform1i(glGetUniformLocation(program, "u_has_history"), mol->backbone_orient_history ? 1 : 0);
        glUniform1f(glGetUniformLocation(program, "u_flip_threshold"), -sinf((float)DEG_TO_RAD(hysteresis_deg)));
        glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT].id);
        glBeginTransformFeedback(GL_POINTS);
        glDrawArrays(GL_POINTS, 0, mol->backbone_count);
        glEndTransformFeedback();
        glUseProgram(0);
        mol->backbone_orient_history = true;
    }

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);
    glDisableVertexAttribArray(5);

    // Pass 3: Evaluate the rings (cross section frames) of the spline, MD_GL_BACKBONE_SEGMENT_COUNT per control point
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT].id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, mol->buffer[GL_BUFFER_BACKBONE_NEIGHBOR].id);

    {
        GLuint program = ctx.program[GL_PROGRAM_SUBDIVIDE_SPLINE].id;
        glUseProgram(program);
        glUniform1i(glGetUniformLocation(program, "u_buf_control_points"), 0);
        glUniform1i(glGetUniformLocation(program, "u_buf_neighbors"), 1);
        glUniform1i(glGetUniformLocation(program, "u_segments"), MD_GL_BACKBONE_SEGMENT_COUNT);
        glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->buffer[GL_BUFFER_BACKBONE_RING_DATA].id);
        glBeginTransformFeedback(GL_POINTS);
        glDrawArraysInstanced(GL_POINTS, 0, MD_GL_BACKBONE_SEGMENT_COUNT, mol->backbone_count);
        glEndTransformFeedback();
        glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, 0);
        glUseProgram(0);
    }

    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(0);
    glDisable(GL_RASTERIZER_DISCARD);

    return true;
}
