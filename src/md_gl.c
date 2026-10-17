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
// Largest radius of a tube (MD_GL_REP_TUBES, Ångström): the radius of a control point is packed as a fraction of it in
// the w of its secondary structure. Must match TUBE_RADIUS_MAX of compute_spline_extract.vert,
// compute_spline_subdivide.vert and backbone.vert
#define TUBE_RADIUS_MAX 8.0f

// Bonds are drawn as one quad each (bond.vert), from a static index buffer of this many quads (u16 indices: at most
// 16384), as instances of the whole buffer followed by a draw of the rest
#define BOND_QUADS_PER_INSTANCE 16384
// Pattern of the dashed partners of delocalized bonds (and of aromatic bonds outside of a ring drawn as a circle)
#define BOND_PARTNER_DASHES 3
#define BOND_PARTNER_FILL   0.5f

// The gaps of the chains (md_gl_mol_set_chain_gaps): the radius of the dashed line relative to the coil of the cartoon
// (and to the thickness of the ribbons), and to the sticks of licorice and ball and stick, and its dashes, of a length of
// their own along it (Ångström, the period of the pattern) whatever the length of the gap
#define CHAIN_GAP_RADIUS_TRACE  0.7f
#define CHAIN_GAP_RADIUS_LINK   0.5f
#define CHAIN_GAP_DASH_PERIOD   1.2f
#define CHAIN_GAP_DASH_FILL     0.5f

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
    GL_PROGRAM_COMPUTE_DERIVED,
    GL_PROGRAM_COUNT
};

enum {
    GL_TEXTURE_BUFFER_0,
    GL_TEXTURE_BUFFER_1,
    GL_TEXTURE_BUFFER_2,
    GL_TEXTURE_BUFFER_3,
    GL_TEXTURE_BUFFER_4,
    GL_TEXTURE_BUFFER_5,
    GL_TEXTURE_BUFFER_6,
    GL_TEXTURE_BUFFER_7,
    GL_TEXTURE_BUFFER_8,
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
    GL_BUFFER_BOND_ORDER,                      // u32 per bond: kind | index into GL_BUFFER_BOND_ORDER_REF << 3 (md_gl_mol_set_bond_orders)
    GL_BUFFER_BOND_ORDER_REF,                  // u32[2] per double or triple bond: atoms bonded to its first and its second atom
    GL_BUFFER_BOND_PARTNER,                    // u32[4] per dashed partner of a delocalized bond, or of an aromatic bond in no circle: bond, atoms bonded to its first and its second atom, 0
    GL_BUFFER_AROMATIC_RING,                   // u32[8] per aromatic ring drawn as a circle: ring index | atom count << 29, its atoms in order around it (RING_ATOM_NONE past the count)
    GL_BUFFER_BACKBONE_DATA,                   // u32: residue index, u32: residue atom offset, u8: CA index, C index and O Index, u8: flags
    GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE,    // u8[4]  (0: Unknown, 1: Coil, 2: Helix, 3: Sheet)
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA,     // Extracted control points
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT,      // Oriented control points (support vector + relation to the next control point)
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT_PREV, // The previous computation's oriented control points, swapped with the one above each computation
    GL_BUFFER_BACKBONE_NEIGHBOR,               // u32[4] per control point: prev, next, prev2, next2 (clamped to the chain)
    GL_BUFFER_BACKBONE_RING_DATA,              // Evaluated cross section frames of the spline, MD_GL_BACKBONE_SEGMENT_COUNT per control point (see compute_spline_subdivide.vert)
    GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX,    // u32, LINE_STRIP_ADJACENCY Indices for legacy and debugging paths
    GL_BUFFER_INSTANCE_TRANSFORM,              // mat4 instance transformation matrices
    GL_BUFFER_DERIVED_DATA,                    // u32[3] per particle derived from a source (md_gl_mol_set_derived): offset of its atoms in the two below, their count, its reference atom
    GL_BUFFER_DERIVED_ATOM,                    // u32 Atoms of the source
    GL_BUFFER_DERIVED_WEIGHT,                  // f32 Their weights
    GL_BUFFER_SLAB,                            // u32[8] per slab (md_gl_mol_set_slabs): its corner count, then up to SLAB_MAX_CORNERS particles in order around it
    GL_BUFFER_CHAIN_GAP,                       // u32[4] per gap of a chain (md_gl_mol_set_chain_gaps): atom a, its neighbour, atom b, its neighbour
    GL_BUFFER_CHAIN_GAP_COMP,                  // u32 per gap of a chain: the component before its break, what it picks as
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

// What the extraction of a control point reads (compute_spline_extract.vert)
typedef struct {
    uint32_t bb_seg_idx;            // Backbone Segment index (To access secondary structure)
    uint32_t ca;                    // The atom (particle) of the control point
    uint32_t c, o;                  // The fallback orientation O - C (a protein: its carbonyl), or with BACKBONE_FLAG_FIXED_FACE the normal of the flat faces
    uint32_t flags;                 // BACKBONE_FLAG_*
} gl_backbone_data_t;

// Flags of a control point (gl_backbone_data_t, and the top byte of the secondary structure of gl_control_point_t)
enum {
    BACKBONE_FLAG_CHAIN_BEG  = 1,
    BACKBONE_FLAG_CHAIN_END  = 2,
    BACKBONE_FLAG_FLIP_NEXT  = 16,  // Set by the orient pass
    BACKBONE_FLAG_FIXED_FACE = 32,  // The flat faces of the profile face O - C, rather than following the curve (compute_spline_orient.vert)
    BACKBONE_FLAG_THROUGH    = 64,  // The spline passes through the control point rather than near it (compute_spline_orient.vert)
};

typedef struct {
    float position[3];
    uint32_t atom_idx;
    float velocity[3];
    float radius;                       // @NOTE: Of a tube at the control point (MD_GL_REP_TUBES), Ångström: TUBE_RADIUS_MAX times the w of its secondary structure, 0 for a protein backbone
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
    uint32_t backbone_index_base;
} gl_ubo_base_t;

typedef struct {
    uint32_t atom_idx[2];
} gl_bond_t;

typedef struct {
    uint32_t id;
    gl_program_t spacefill[MAX_SHADER_PERMUTATIONS];
    gl_program_t licorice[MAX_SHADER_PERMUTATIONS];
    gl_program_t weak_bonds[MAX_SHADER_PERMUTATIONS];
    gl_program_t bond_partners[MAX_SHADER_PERMUTATIONS];
    gl_program_t chain_gaps[MAX_SHADER_PERMUTATIONS];
    gl_program_t aromatic_rings[MAX_SHADER_PERMUTATIONS];
    gl_program_t ribbons[MAX_SHADER_PERMUTATIONS];
    gl_program_t cartoon[MAX_SHADER_PERMUTATIONS];
    gl_program_t tubes[MAX_SHADER_PERMUTATIONS];
    gl_program_t slabs[MAX_SHADER_PERMUTATIONS];
    gl_program_t cartoon_rc[MAX_SHADER_PERMUTATIONS];   // Ray cast (backbone_rc.vert)
    gl_program_t tubes_rc[MAX_SHADER_PERMUTATIONS];
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
    uint32_t bond_order_count;      // Bonds GL_BUFFER_BOND_ORDER was computed for, 0 when the bonds changed since
    uint32_t bond_partner_count;    // Entries of GL_BUFFER_BOND_PARTNER
    uint32_t aromatic_ring_count;   // Entries of GL_BUFFER_AROMATIC_RING
    uint32_t gap_trace_count;       // Entries of GL_BUFFER_CHAIN_GAP for the cartoon and ribbons, first
    uint32_t gap_link_count;        // Entries of GL_BUFFER_CHAIN_GAP for licorice and ball and stick, after them
    uint32_t backbone_count;
    uint32_t backbone_control_point_index_count;

    // GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT holds a previous computation that the next one can
    // continue from. Cleared on creation and by md_gl_mol_reset_backbone_history.
    bool backbone_orient_history;

    backbone_bounds_t backbone_bounds;

    // The backbone is tubes through the particles (md_gl_mol_set_tubes), the largest half width (l + r) of their profiles
    bool  tubes;
    float tube_max_extent;

    // Entries of GL_BUFFER_SLAB (md_gl_mol_set_slabs)
    uint32_t slab_count;

    // Particles derived from the atoms of a source molecule (md_gl_mol_set_derived), the arrays of their buffers kept on
    // the CPU too: where the control points of the tubes are, for their culling and level of detail
    struct {
        uint32_t  source_id;
        uint32_t  source_atom_count;
        uint32_t  atom_count;       // Entries of atom and weight
        uint32_t* data;             // [3 * particle count] (GL_BUFFER_DERIVED_DATA)
        uint32_t* atom;
        float*    weight;
    } derived;

    gl_buffer_t buffer[GL_BUFFER_COUNT];
} molecule_t;

typedef struct {
    uint32_t id;
    gl_texture_t tex;
} palette_t;

enum {
    BACKBONE_PROFILE_STADIUM,   // Cartoon: a stadium (a 2D capsule), flat caps
    BACKBONE_PROFILE_BOX,       // Ribbons
    BACKBONE_PROFILE_TUBE,      // Tubes: the stadium with round caps
    BACKBONE_PROFILE_COUNT
};

// Static triangulation of one backbone instance (the part of the spline that belongs to one residue)
typedef struct {
    uint32_t segments;          // S: segments along the spline, S + 1 rings
    uint32_t profile_count;     // P: vertices per ring
    uint32_t outline_count;     // O: vertices on the outline of a cap
    uint32_t face_subdivisions; // K: segments across the wide faces of the box profile
    uint32_t dome_rings;        // D: segments of a round cap of the tube profile, D + 1 rings from its rim to its top
    uint32_t index_count;
    gl_buffer_t index_buffer;   // u16, GL_TRIANGLES
} backbone_mesh_t;

typedef struct {
    uint32_t id;
    uint32_t mol_id; // molecule
    uint32_t pal_id; // palette
    gl_buffer_t atom_color;

    // Weak bonds of the representation (md_gl_rep_set_weak_bonds), capacity in bonds of the buffers
    uint32_t weak_bond_count;
    uint32_t weak_bond_capacity;
    gl_buffer_t weak_bond_atoms;    // u32[2]
    gl_buffer_t weak_bond_weight;   // f32
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
    gl_buffer_t bond_quad_index;    // u16, BOND_QUADS_PER_INSTANCE quads of two triangles
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

// The spheres of the chains from the positions of their control points (ca_xyz)
static void backbone_bounds_update_spheres(backbone_bounds_t* bounds) {
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
    backbone_bounds_update_spheres(bounds);
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

        // The bond orders refer to the bonds by index: set them again (md_gl_mol_set_bond_orders)
        mol->bond_order_count = 0;
        mol->bond_partner_count = 0;
        mol->aromatic_ring_count = 0;

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
        for (uint32_t i = 0; i < count; ++i) {
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
// Tube vertices ring * P + k for ring in [0, S], k in [0, P), then a cap at the beginning and one at the end. A cap is
// flat, its center and its O outline vertices, or round for the tube profile: a dome of D + 1 rings of O vertices,
// d * O + m for ring d in [0, D] (ring 0 is the rim, ring D the top, a point or a segment).
// Triangles are counter clockwise seen from the outside.
// segments: even, profile_count: even (stadium, tube), face_subdivisions: >= 1 (box)
static bool create_backbone_mesh(backbone_mesh_t* mesh, uint32_t profile, uint32_t segments, uint32_t profile_count, uint32_t face_subdivisions) {
    ASSERT(mesh);
    ASSERT(segments >= 2 && segments % 2 == 0);
    const uint32_t S = segments;
    uint32_t P, O, K = 0, D = 0;
    uint32_t edges[64][2];
    uint32_t edge_count = 0;

    if (profile == BACKBONE_PROFILE_STADIUM || profile == BACKBONE_PROFILE_TUBE) {
        ASSERT(profile_count >= 4 && profile_count % 2 == 0 && profile_count <= 64);
        P = profile_count;
        O = P;
        // The quarter circle of a dome as finely as the profile
        D = profile == BACKBONE_PROFILE_TUBE ? MAX(1, P / 4) : 0;
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

    const uint32_t cap_vertex_count = D ? (D + 1) * O : O + 1;
    const uint32_t cap_index_count  = D ? D * O * 6 : O * 3;
    const uint32_t vertex_count = (S + 1) * P + 2 * cap_vertex_count;
    const uint32_t index_count  = S * edge_count * 6 + 2 * cap_index_count;
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
    const uint32_t cap_base[2] = { (S + 1) * P, (S + 1) * P + cap_vertex_count };
    for (uint32_t cap = 0; cap < 2; ++cap) {
        if (D == 0) {
            const uint16_t center = (uint16_t)cap_base[cap];
            for (uint32_t m = 0; m < O; ++m) {
                const uint16_t r0 = (uint16_t)(cap_base[cap] + 1 + m);
                const uint16_t r1 = (uint16_t)(cap_base[cap] + 1 + (m + 1) % O);
                idx[len++] = center;
                idx[len++] = cap == 0 ? r1 : r0;
                idx[len++] = cap == 0 ? r0 : r1;
            }
            continue;
        }
        // A dome: ring d + 1 lies further out than ring d, beyond the end of the tube (backwards for the beginning). The
        // end continues the tube as its next ring would, the beginning as its previous ring would. The top ring has no
        // area (all of its vertices are on a point or a segment), so nothing closes it.
        for (uint32_t d = 0; d < D; ++d) {
            for (uint32_t m = 0; m < O; ++m) {
                const uint16_t a = (uint16_t)(cap_base[cap] + d * O + m);
                const uint16_t b = (uint16_t)(cap_base[cap] + d * O + (m + 1) % O);
                const uint16_t c = (uint16_t)(a + O);
                const uint16_t e = (uint16_t)(b + O);
                if (cap == 1) {
                    idx[len++] = a; idx[len++] = b; idx[len++] = c;
                    idx[len++] = b; idx[len++] = e; idx[len++] = c;
                } else {
                    idx[len++] = c; idx[len++] = e; idx[len++] = a;
                    idx[len++] = e; idx[len++] = b; idx[len++] = a;
                }
            }
        }
    }
    ASSERT(len == index_count);

    gl_buffer_conditional_delete(&mesh->index_buffer);
    mesh->index_buffer = gl_buffer_create(index_count * sizeof(uint16_t), idx, GL_STATIC_DRAW);
    mesh->segments = S;
    mesh->profile_count = P;
    mesh->outline_count = O;
    mesh->face_subdivisions = K;
    mesh->dome_rings = D;
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

        if (!md_gl_shader_compile(vert_shader, (str_t){(const char*)compute_derived_vert, compute_derived_vert_size}, 0, 0)) {
            return;
        }
        ctx.program[GL_PROGRAM_COMPUTE_DERIVED].id = glCreateProgram();
        const GLuint shaders[] = { vert_shader };
        // Into the positions and the velocities of the particles, two buffers
        const GLchar* varyings[] = { "out_position", "out_velocity" };
        if (!md_gl_program_attach_and_link_transform_feedback(ctx.program[GL_PROGRAM_COMPUTE_DERIVED].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_SEPARATE_ATTRIBS)) {
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
        const GLchar* varyings[] = { "out_position", "out_atom_idx", "out_velocity", "out_radius", "out_secondary_structure_and_flags", "out_support_and_tangent_vector" };
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
        const GLchar* varyings[] = { "out_position", "out_atom_idx", "out_velocity", "out_radius", "out_secondary_structure_and_flags", "out_support_and_tangent_vector" };
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
            if (!create_backbone_mesh(&ctx.backbone_mesh[BACKBONE_PROFILE_STADIUM][lod], BACKBONE_PROFILE_STADIUM, segments, profile_count, 0) ||
                !create_backbone_mesh(&ctx.backbone_mesh[BACKBONE_PROFILE_BOX][lod],     BACKBONE_PROFILE_BOX,     segments, 0, face_subdivisions) ||
                !create_backbone_mesh(&ctx.backbone_mesh[BACKBONE_PROFILE_TUBE][lod],    BACKBONE_PROFILE_TUBE,    segments, profile_count, 0)) {
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

    {
        // Quad q has the vertices 4q .. 4q + 3 (bond.vert): triangles (0, 1, 2) and (2, 1, 3), counter clockwise
        STATIC_ASSERT(BOND_QUADS_PER_INSTANCE * 4 <= 65536, "Bond quad indices must fit in u16");
        const uint32_t count = BOND_QUADS_PER_INSTANCE * 6;
        uint16_t* indices = md_alloc(md_get_heap_allocator(), count * sizeof(uint16_t));
        for (uint32_t q = 0; q < BOND_QUADS_PER_INSTANCE; ++q) {
            const uint16_t base = (uint16_t)(q * 4);
            indices[q * 6 + 0] = base + 0;
            indices[q * 6 + 1] = base + 1;
            indices[q * 6 + 2] = base + 2;
            indices[q * 6 + 3] = base + 2;
            indices[q * 6 + 4] = base + 1;
            indices[q * 6 + 5] = base + 3;
        }
        ctx.bond_quad_index = gl_buffer_create(count * sizeof(uint16_t), indices, GL_STATIC_DRAW);
        md_free(md_get_heap_allocator(), indices, count * sizeof(uint16_t));
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
    gl_buffer_conditional_delete(&ctx.bond_quad_index);

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
    const str_t bond_vert_src = {(const char*)bond_vert, bond_vert_size};
    const str_t bond_frag_src = {(const char*)bond_frag, bond_frag_size};
    const str_t ring_vert_src = {(const char*)ring_vert, ring_vert_size};
    const str_t ring_frag_src = {(const char*)ring_frag, ring_frag_size};
    if (!create_permuted_program(STR_LIT("SpaceFill"), shaders->spacefill,  (str_t){(const char*)spacefill_vert, spacefill_vert_size}, (str_t){(const char*)spacefill_geom, spacefill_geom_size},   (str_t){(const char*)spacefill_frag, spacefill_frag_size},    str, (str_t){0})) return handle;
    if (!create_permuted_program(STR_LIT("Licorice"),  shaders->licorice,   bond_vert_src, (str_t){0}, bond_frag_src, str, STR_LIT("#define WEAK_BONDS 0"))) return handle;
    if (!create_permuted_program(STR_LIT("WeakBonds"), shaders->weak_bonds, bond_vert_src, (str_t){0}, bond_frag_src, str, STR_LIT("#define WEAK_BONDS 1"))) return handle;
    if (!create_permuted_program(STR_LIT("BondPartners"), shaders->bond_partners, bond_vert_src, (str_t){0}, bond_frag_src, str, STR_LIT("#define BOND_PARTNERS 1"))) return handle;
    if (!create_permuted_program(STR_LIT("ChainGaps"), shaders->chain_gaps, bond_vert_src, (str_t){0}, bond_frag_src, str, STR_LIT("#define CHAIN_GAPS 1"))) return handle;
    if (!create_permuted_program(STR_LIT("AromaticRings"), shaders->aromatic_rings, ring_vert_src, (str_t){0}, ring_frag_src, str, (str_t){0})) return handle;
    if (!create_permuted_program(STR_LIT("Ribbons"),   shaders->ribbons,    backbone_vert_src, (str_t){0}, backbone_frag_src, str, STR_LIT("#define REP_RIBBONS 1"))) return handle;
    if (!create_permuted_program(STR_LIT("Cartoon"),   shaders->cartoon,    backbone_vert_src, (str_t){0}, backbone_frag_src, str, STR_LIT("#define REP_RIBBONS 0"))) return handle;
    if (!create_permuted_program(STR_LIT("Tubes"),     shaders->tubes,      backbone_vert_src, (str_t){0}, backbone_frag_src, str, STR_LIT("#define REP_TUBES 1"))) return handle;
    const str_t backbone_rc_vert_src = {(const char*)backbone_rc_vert, backbone_rc_vert_size};
    const str_t backbone_rc_frag_src = {(const char*)backbone_rc_frag, backbone_rc_frag_size};
    if (!create_permuted_program(STR_LIT("CartoonRaycast"), shaders->cartoon_rc, backbone_rc_vert_src, (str_t){0}, backbone_rc_frag_src, str, STR_LIT("#define REP_TUBES 0"))) return handle;
    if (!create_permuted_program(STR_LIT("TubesRaycast"),   shaders->tubes_rc,   backbone_rc_vert_src, (str_t){0}, backbone_rc_frag_src, str, STR_LIT("#define REP_TUBES 1"))) return handle;
    if (!create_permuted_program(STR_LIT("Slabs"),     shaders->slabs,      (str_t){(const char*)slab_vert, slab_vert_size}, (str_t){0}, (str_t){(const char*)slab_frag, slab_frag_size}, str, (str_t){0})) return handle;

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
            if (glIsProgram(shaders->weak_bonds[i].id)) glDeleteProgram(shaders->weak_bonds[i].id);
            if (glIsProgram(shaders->bond_partners[i].id)) glDeleteProgram(shaders->bond_partners[i].id);
            if (glIsProgram(shaders->chain_gaps[i].id)) glDeleteProgram(shaders->chain_gaps[i].id);
            if (glIsProgram(shaders->aromatic_rings[i].id)) glDeleteProgram(shaders->aromatic_rings[i].id);
            if (glIsProgram(shaders->ribbons[i].id))    glDeleteProgram(shaders->ribbons[i].id);
            if (glIsProgram(shaders->cartoon[i].id))    glDeleteProgram(shaders->cartoon[i].id);
            if (glIsProgram(shaders->tubes[i].id))      glDeleteProgram(shaders->tubes[i].id);
            if (glIsProgram(shaders->slabs[i].id))      glDeleteProgram(shaders->slabs[i].id);
            if (glIsProgram(shaders->cartoon_rc[i].id)) glDeleteProgram(shaders->cartoon_rc[i].id);
            if (glIsProgram(shaders->tubes_rc[i].id))   glDeleteProgram(shaders->tubes_rc[i].id);
        }
        md_handle_pool_free_slot(&ctx.shader_pool, handle.id);
    }
}

// ### BOND ORDERS ###
//
// What bond.vert needs to show the order of a bond: the plane of the bond, given by its surroundings, an atom bonded to
// either end which lies in it. The plane is that of the sp2 centers: the substituents of a double bond, the ring of an
// aromatic bond, the group of a delocalized one.
//   - Double and triple bonds have a kind per bond, which replaces the stick, and their surroundings.
//   - Aromatic rings are a circle inside the ring (ring.vert), an entry per ring: its index and its atoms. Their bonds
//     keep the stick. Two fused rings have a circle each.
//   - Delocalized bonds keep the stick and get dashed partners, a list of entries of their own: a bond and the
//     surroundings which give the side of the partner. So does an aromatic bond in no circle (in a ring larger than
//     MD_RING_MAX_SIZE, or of a system without ring flags): no aromatic bond goes unmarked.
//
// BEGIN BOND ORDER DERIVATION (plain arrays in, the contents of the buffers out; tested on its own)

// Kinds of bonds in GL_BUFFER_BOND_ORDER (bond.vert)
enum {
    BOND_KIND_PLAIN   = 0,
    BOND_KIND_DOUBLE  = 1,
    BOND_KIND_TRIPLE  = 2,
};
#define BOND_REF_NONE 0xFFFFFFFFu

// Entries of GL_BUFFER_AROMATIC_RING (ring.vert): u32[8], the index of the ring | its atom count << 29, then its atoms
// in order around it and RING_ATOM_NONE past the count. Seven atoms at most, MD_RING_MAX_SIZE.
#define RING_ENTRY_U32  8
#define RING_ATOM_NONE  0xFFFFFFFFu
#define RING_COUNT_SHIFT 29

typedef struct {
    uint32_t        bond_count;
    const uint32_t* pairs;          // [bond_count * 2]
    const uint32_t* flags;          // [bond_count] md_bond_flags_t
    uint32_t        atom_count;
    const uint32_t* adj_offset;     // [atom_count + 1] Bonds of each atom in adj_atom and adj_bond
    const uint32_t* adj_atom;       // The other atom
    const uint32_t* adj_bond;       // The bond
    const uint8_t*  hybridization;  // [atom_count] md_hybridization_t (optional)
    const float*    xyz;            // [atom_count * 3] Reference coordinates (optional)
    uint32_t        ring_count;
    const uint32_t* ring_offset;    // [ring_count + 1]
    const int32_t*  ring_atom;      // Atoms of each ring, in order around it
    const md_ring_flags_t* ring_flags;  // [ring_count] (optional, no circles without)
} bond_order_input_t;

static uint32_t bo_find_bond(const bond_order_input_t* in, uint32_t a, uint32_t b) {
    for (uint32_t j = in->adj_offset[a]; j < in->adj_offset[a + 1]; ++j) {
        if (in->adj_atom[j] == b) return in->adj_bond[j];
    }
    return BOND_REF_NONE;
}

static bool bo_valid_bond(const bond_order_input_t* in, uint32_t i) {
    const uint32_t a = in->pairs[2 * i + 0];
    const uint32_t b = in->pairs[2 * i + 1];
    return a < in->atom_count && b < in->atom_count && a != b && !(in->flags[i] & MD_BOND_FLAG_COORDINATE);
}

// c is off the axis e - o by more than 10 degrees in the reference coordinates (or there are none)
static bool bo_off_axis(const bond_order_input_t* in, uint32_t e, uint32_t o, uint32_t c) {
    if (!in->xyz) return true;
    const float* pe = in->xyz + 3 * e;
    const float* po = in->xyz + 3 * o;
    const float* pc = in->xyz + 3 * c;
    const float u[3] = {po[0] - pe[0], po[1] - pe[1], po[2] - pe[2]};
    const float v[3] = {pc[0] - pe[0], pc[1] - pe[1], pc[2] - pe[2]};
    const float x[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
    const float uu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    const float vv = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    const float xx = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    return xx > 0.03f * uu * vv;
}

static bool bo_is_multiple(uint32_t flags) {
    return md_bond_order((md_bond_flags_t)flags) >= 2 || (flags & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED));
}

// An atom bonded to end e of the bond to o which gives the plane of the bond: off the axis, through a bond that is not
// a coordination bond. An end which is not sp2 (the S of a sulfone, the P of a phosphine oxide) prefers one through a
// multiple bond, so the bond lies in the plane they share (O=S=O). A linear end (sp) gives none.
static uint32_t bo_plane_ref(const bond_order_input_t* in, uint32_t e, uint32_t o) {
    const uint8_t hyb = in->hybridization ? in->hybridization[e] : (uint8_t)MD_HYBRIDIZATION_UNKNOWN;
    if (hyb == MD_HYBRIDIZATION_SP) return BOND_REF_NONE;
    uint32_t ref = BOND_REF_NONE;
    for (uint32_t j = in->adj_offset[e]; j < in->adj_offset[e + 1]; ++j) {
        const uint32_t c = in->adj_atom[j];
        const uint32_t f = in->flags[in->adj_bond[j]];
        if (c == o || (f & MD_BOND_FLAG_COORDINATE) || !bo_off_axis(in, e, o, c)) continue;
        if (hyb != MD_HYBRIDIZATION_SP2 && bo_is_multiple(f)) return c;
        if (ref == BOND_REF_NONE) ref = c;
    }
    return ref;
}

// A delocalized bond has its partner within the group. The center of the group is the end with more delocalized bonds
// (the C of a carboxylate or guanidinium, the N of a nitro group); the partner lies towards the next of them in order of
// atom index, so the partners of a carboxylate are both between its O, and those of a guanidinium one in each angle.
static void bo_delocalized_refs(const bond_order_input_t* in, uint32_t a, uint32_t b, uint32_t* ref_a, uint32_t* ref_b) {
    uint32_t na = 0, nb = 0;
    for (uint32_t j = in->adj_offset[a]; j < in->adj_offset[a + 1]; ++j) na += (in->flags[in->adj_bond[j]] & MD_BOND_FLAG_DELOCALIZED) ? 1 : 0;
    for (uint32_t j = in->adj_offset[b]; j < in->adj_offset[b + 1]; ++j) nb += (in->flags[in->adj_bond[j]] & MD_BOND_FLAG_DELOCALIZED) ? 1 : 0;

    *ref_a = BOND_REF_NONE;
    *ref_b = BOND_REF_NONE;
    if (na < 2 && nb < 2) {
        // No group to speak of: the plane only
        *ref_a = bo_plane_ref(in, a, b);
        *ref_b = bo_plane_ref(in, b, a);
        return;
    }
    const uint32_t center = na >= nb ? a : b;
    const uint32_t other  = center == a ? b : a;

    // The delocalized neighbors of the center, sorted
    uint32_t nbr[8];
    uint32_t n = 0;
    for (uint32_t j = in->adj_offset[center]; j < in->adj_offset[center + 1] && n < 8; ++j) {
        if (in->flags[in->adj_bond[j]] & MD_BOND_FLAG_DELOCALIZED) {
            uint32_t k = n++;
            const uint32_t c = in->adj_atom[j];
            while (k > 0 && nbr[k - 1] > c) {
                nbr[k] = nbr[k - 1];
                --k;
            }
            nbr[k] = c;
        }
    }
    uint32_t ref = BOND_REF_NONE;
    for (uint32_t k = 0; k < n; ++k) {
        if (nbr[k] == other) {
            ref = nbr[(k + 1) % n];
            break;
        }
    }
    if (ref == other) ref = BOND_REF_NONE;
    if (center == a) *ref_a = ref;
    else             *ref_b = ref;
}

// An aromatic bond in no circle has its partner on the side of the aromatic bonds which continue from its ends (along its
// ring, when that is too large for a circle), failing those in its plane
static uint32_t bo_aromatic_ref(const bond_order_input_t* in, uint32_t e, uint32_t o) {
    for (uint32_t j = in->adj_offset[e]; j < in->adj_offset[e + 1]; ++j) {
        const uint32_t c = in->adj_atom[j];
        const uint32_t f = in->flags[in->adj_bond[j]];
        if (c != o && (f & MD_BOND_FLAG_AROMATIC) && !(f & MD_BOND_FLAG_COORDINATE) && bo_off_axis(in, e, o, c)) return c;
    }
    return bo_plane_ref(in, e, o);
}

// The ring r is drawn as a circle: aromatic (MD_RING_FLAG_AROMATIC, which says that every bond of it is), of at most
// RING_ENTRY_U32 - 1 atoms, each bonded to the next by a bond the molecule draws
static bool bo_circle_ring(const bond_order_input_t* in, uint32_t r) {
    if (!in->ring_flags || !(in->ring_flags[r] & MD_RING_FLAG_AROMATIC) || r >= (1u << RING_COUNT_SHIFT)) return false;
    const uint32_t beg = in->ring_offset[r];
    const uint32_t n   = in->ring_offset[r + 1] - beg;
    if (n < 3 || n > RING_ENTRY_U32 - 1) return false;
    const int32_t* ring = in->ring_atom + beg;
    for (uint32_t k = 0; k < n; ++k) {
        const uint32_t a = (uint32_t)ring[k];
        const uint32_t b = (uint32_t)ring[(k + 1) % n];
        const uint32_t bond = (a < in->atom_count && b < in->atom_count) ? bo_find_bond(in, a, b) : BOND_REF_NONE;
        if (bond == BOND_REF_NONE || !bo_valid_bond(in, bond)) return false;
    }
    return true;
}

// The aromatic rings drawn as circles, RING_ENTRY_U32 u32 per entry (see GL_BUFFER_AROMATIC_RING). Marks the bonds of
// those rings in out_covered [bond_count] (optional): they need no dashed partner (compute_bond_partners).
// Returns the number of entries, and writes them when out_ring is not NULL: call once to count.
static uint32_t compute_aromatic_rings(uint32_t* out_ring, uint8_t* out_covered, const bond_order_input_t* in) {
    uint32_t count = 0;
    for (uint32_t r = 0; r < in->ring_count; ++r) {
        if (!bo_circle_ring(in, r)) continue;
        const uint32_t beg = in->ring_offset[r];
        const uint32_t n   = in->ring_offset[r + 1] - beg;
        const int32_t* ring = in->ring_atom + beg;
        if (out_ring) {
            uint32_t* e = out_ring + RING_ENTRY_U32 * count;
            e[0] = r | (n << RING_COUNT_SHIFT);
            for (uint32_t k = 0; k < RING_ENTRY_U32 - 1; ++k) {
                e[1 + k] = k < n ? (uint32_t)ring[k] : RING_ATOM_NONE;
            }
        }
        if (out_covered) {
            for (uint32_t k = 0; k < n; ++k) {
                out_covered[bo_find_bond(in, (uint32_t)ring[k], (uint32_t)ring[(k + 1) % n])] = 1;
            }
        }
        ++count;
    }
    return count;
}

// Double and triple bonds. out_order [bond_count]: kind | entry << 3, where entry indexes out_ref for double and triple
// bonds. out_ref [2 * bond_count]: the atoms bonded to the first and to the second atom of each such bond, or
// BOND_REF_NONE. Returns the number of entries of out_ref. Aromatic and delocalized bonds are plain here: they keep their
// stick (compute_bond_partners).
static uint32_t compute_bond_orders(uint32_t* out_order, uint32_t* out_ref, const bond_order_input_t* in) {
    uint32_t entries = 0;
    for (uint32_t i = 0; i < in->bond_count; ++i) {
        const uint32_t f = in->flags[i];
        out_order[i] = BOND_KIND_PLAIN;
        if (!bo_valid_bond(in, i) || (f & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED))) continue;
        const uint32_t a = in->pairs[2 * i + 0];
        const uint32_t b = in->pairs[2 * i + 1];
        const int order = md_bond_order((md_bond_flags_t)f);
        if (order == 2) {
            out_ref[2 * entries + 0] = bo_plane_ref(in, a, b);
            out_ref[2 * entries + 1] = bo_plane_ref(in, b, a);
            out_order[i] = BOND_KIND_DOUBLE | (entries << 3);
            ++entries;
        } else if (order >= 3) {
            // No plane, perpendicular to the view
            out_ref[2 * entries + 0] = BOND_REF_NONE;
            out_ref[2 * entries + 1] = BOND_REF_NONE;
            out_order[i] = BOND_KIND_TRIPLE | (entries << 3);
            ++entries;
        }
    }
    return entries;
}

// The dashed partners, u32[4] per entry: the bond, the atoms bonded to its first and to its second atom which give the
// side of the partner, and 0.
//   - A delocalized bond has one, within its group (bo_delocalized_refs).
//   - An aromatic bond which no circle covers (covered [bond_count], from compute_aromatic_rings, optional) has one on
//     the side of the aromatic bonds that continue from it (bo_aromatic_ref). Those of the circles have none.
// Returns the number of entries, and writes them when out_partner is not NULL: call once to count.
static uint32_t compute_bond_partners(uint32_t* out_partner, const uint8_t* covered, const bond_order_input_t* in) {
    uint32_t count = 0;
    for (uint32_t i = 0; i < in->bond_count; ++i) {
        const uint32_t f = in->flags[i];
        const bool aromatic    = (f & MD_BOND_FLAG_AROMATIC) && !(covered && covered[i]);
        const bool delocalized = (f & MD_BOND_FLAG_DELOCALIZED) && !(f & MD_BOND_FLAG_AROMATIC);
        if (!(aromatic || delocalized) || !bo_valid_bond(in, i)) continue;
        if (out_partner) {
            const uint32_t a = in->pairs[2 * i + 0];
            const uint32_t b = in->pairs[2 * i + 1];
            uint32_t* e = out_partner + 4 * count;
            e[0] = i;
            if (aromatic) {
                e[1] = bo_aromatic_ref(in, a, b);
                e[2] = bo_aromatic_ref(in, b, a);
            } else {
                bo_delocalized_refs(in, a, b, &e[1], &e[2]);
            }
            e[3] = 0;
        }
        ++count;
    }
    return count;
}

// END BOND ORDER DERIVATION

void md_gl_mol_set_bond_orders(md_gl_mol_t handle, const md_system_t* sys) {
    molecule_t* mol = mol_lookup(handle.id);
    if (!mol) {
        MD_LOG_ERROR("Invalid molecule");
        return;
    }
    mol->bond_order_count = 0;
    mol->bond_partner_count = 0;
    mol->aromatic_ring_count = 0;
    if (!sys) {
        MD_LOG_ERROR("The system is missing");
        return;
    }
    const uint32_t bond_count = (uint32_t)sys->bond.count;
    if (bond_count != mol->bond_count) {
        MD_LOG_ERROR("The system has %u bonds and the molecule %u: set the bonds first (md_gl_mol_set_bonds)", bond_count, mol->bond_count);
        return;
    }
    if (bond_count == 0 || !sys->bond.pairs || !sys->bond.flags) {
        return;
    }
    if (bond_count >= (1u << 29)) {
        MD_LOG_ERROR("Too many bonds for bond orders");
        return;
    }
    STATIC_ASSERT(sizeof(md_atom_pair_t) == 2 * sizeof(uint32_t), "Bond pairs are read as u32[2]");

    md_allocator_i* alloc = md_get_heap_allocator();
    const uint32_t atom_count = mol->atom_count;

    // Connectivity from the pairs, which are what the molecule draws
    const uint32_t* pairs = (const uint32_t*)sys->bond.pairs;
    uint32_t* adj_offset = md_alloc(alloc, (atom_count + 1) * sizeof(uint32_t));
    uint32_t* adj_atom   = md_alloc(alloc, 2 * bond_count * sizeof(uint32_t));
    uint32_t* adj_bond   = md_alloc(alloc, 2 * bond_count * sizeof(uint32_t));
    uint32_t* flags      = md_alloc(alloc, bond_count * sizeof(uint32_t));
    uint8_t*  hyb        = md_alloc(alloc, atom_count * sizeof(uint8_t));
    uint32_t* order      = md_alloc(alloc, bond_count * sizeof(uint32_t));
    uint32_t* ref        = md_alloc(alloc, 2 * bond_count * sizeof(uint32_t));

    MEMSET(adj_offset, 0, (atom_count + 1) * sizeof(uint32_t));
    for (uint32_t i = 0; i < bond_count; ++i) {
        const uint32_t a = pairs[2 * i + 0];
        const uint32_t b = pairs[2 * i + 1];
        flags[i] = (uint32_t)sys->bond.flags[i];
        if (a < atom_count && b < atom_count) {
            adj_offset[a + 1] += 1;
            adj_offset[b + 1] += 1;
        }
    }
    for (uint32_t i = 0; i < atom_count; ++i) {
        adj_offset[i + 1] += adj_offset[i];
    }
    uint32_t* cursor = md_alloc(alloc, atom_count * sizeof(uint32_t));
    MEMCPY(cursor, adj_offset, atom_count * sizeof(uint32_t));
    for (uint32_t i = 0; i < bond_count; ++i) {
        const uint32_t a = pairs[2 * i + 0];
        const uint32_t b = pairs[2 * i + 1];
        if (a < atom_count && b < atom_count) {
            adj_atom[cursor[a]] = b;
            adj_bond[cursor[a]++] = i;
            adj_atom[cursor[b]] = a;
            adj_bond[cursor[b]++] = i;
        }
    }
    md_free(alloc, cursor, atom_count * sizeof(uint32_t));

    for (uint32_t i = 0; i < atom_count; ++i) {
        hyb[i] = (uint8_t)md_atom_hybridization(&sys->atom, i);
    }

    const bool has_xyz = md_system_state_has_coords(&sys->reference) && sys->reference.num_atoms >= atom_count;
    const uint32_t ring_count = (uint32_t)md_index_data_num_ranges(&sys->ring);
    // Without flags for every ring there are no circles, and the aromatic bonds get dashed partners
    const bool has_ring_flags = ring_count && sys->ring_flags && md_array_size(sys->ring_flags) == ring_count;

    const bond_order_input_t in = {
        .bond_count    = bond_count,
        .pairs         = pairs,
        .flags         = flags,
        .atom_count    = atom_count,
        .adj_offset    = adj_offset,
        .adj_atom      = adj_atom,
        .adj_bond      = adj_bond,
        .hybridization = hyb,
        .xyz           = has_xyz ? (const float*)sys->reference.xyz : NULL,
        .ring_count    = ring_count,
        .ring_offset   = ring_count ? sys->ring.offsets : NULL,
        .ring_atom     = ring_count ? sys->ring.indices : NULL,
        .ring_flags    = has_ring_flags ? sys->ring_flags : NULL,
    };
    const uint32_t entries = compute_bond_orders(order, ref, &in);

    // The circles first: the bonds they cover need no dashed partner
    uint8_t* covered = md_alloc(alloc, bond_count * sizeof(uint8_t));
    MEMSET(covered, 0, bond_count * sizeof(uint8_t));
    const uint32_t ring_entry_count = compute_aromatic_rings(NULL, NULL, &in);
    const size_t   ring_entry_bytes = MAX(ring_entry_count, 1) * RING_ENTRY_U32 * sizeof(uint32_t);
    uint32_t* ring_entry = md_alloc(alloc, ring_entry_bytes);
    MEMSET(ring_entry, 0, ring_entry_bytes);
    compute_aromatic_rings(ring_entry, covered, &in);

    const uint32_t partner_count = compute_bond_partners(NULL, covered, &in);
    uint32_t* partner = md_alloc(alloc, MAX(partner_count, 1) * 4 * sizeof(uint32_t));
    MEMSET(partner, 0, MAX(partner_count, 1) * 4 * sizeof(uint32_t));
    compute_bond_partners(partner, covered, &in);

    gl_buffer_conditional_delete(&mol->buffer[GL_BUFFER_BOND_ORDER]);
    gl_buffer_conditional_delete(&mol->buffer[GL_BUFFER_BOND_ORDER_REF]);
    gl_buffer_conditional_delete(&mol->buffer[GL_BUFFER_BOND_PARTNER]);
    gl_buffer_conditional_delete(&mol->buffer[GL_BUFFER_AROMATIC_RING]);
    mol->buffer[GL_BUFFER_BOND_ORDER]     = gl_buffer_create(bond_count * sizeof(uint32_t), order, GL_STATIC_DRAW);
    mol->buffer[GL_BUFFER_BOND_ORDER_REF] = gl_buffer_create(MAX(entries, 1) * 2 * sizeof(uint32_t), ref, GL_STATIC_DRAW);
    mol->buffer[GL_BUFFER_BOND_PARTNER]   = gl_buffer_create(MAX(partner_count, 1) * 4 * sizeof(uint32_t), partner, GL_STATIC_DRAW);
    mol->buffer[GL_BUFFER_AROMATIC_RING]  = gl_buffer_create((uint32_t)ring_entry_bytes, ring_entry, GL_STATIC_DRAW);
    mol->bond_order_count    = bond_count;
    mol->bond_partner_count  = partner_count;
    mol->aromatic_ring_count = ring_entry_count;
    md_free(alloc, partner, MAX(partner_count, 1) * 4 * sizeof(uint32_t));
    md_free(alloc, ring_entry, ring_entry_bytes);
    md_free(alloc, covered, bond_count * sizeof(uint8_t));

    md_free(alloc, adj_offset, (atom_count + 1) * sizeof(uint32_t));
    md_free(alloc, adj_atom,   2 * bond_count * sizeof(uint32_t));
    md_free(alloc, adj_bond,   2 * bond_count * sizeof(uint32_t));
    md_free(alloc, flags,      bond_count * sizeof(uint32_t));
    md_free(alloc, hyb,        atom_count * sizeof(uint8_t));
    md_free(alloc, order,      bond_count * sizeof(uint32_t));
    md_free(alloc, ref,        2 * bond_count * sizeof(uint32_t));
}

// ### CHAIN GAPS ###
//
// One entry per break (md_gl_mol_set_chain_gaps), u32[4] in GL_BUFFER_CHAIN_GAP: the atoms of its two ends and a
// neighbour of each, (a, neighbour of a, b, neighbour of b), and in GL_BUFFER_CHAIN_GAP_COMP the component before the
// break, which it picks as. bond.vert (CHAIN_GAPS) draws the gap from 5/6 of the way from the neighbour of a to a, to
// the same point at b: where the uniform cubic B-spline of the cartoon ends at the control point which ends a chain,
// with its neighbours clamped to the chain (compute_spline_subdivide.vert), so that the line starts at the end of the
// cartoon rather than a little past it. An end which is its own neighbour is the atom itself.
// The entries of the cartoon and ribbons (trace) come first: from the control point (CA, or the bead of a trace) of the
// protein backbone on either side, or where a side has no cartoon (a piece of chain too short for one), from the atom
// of the link there. Then those of licorice and ball and stick (link): between the atoms a link would join
// (md_util_chain_break_atoms), each its own neighbour.
void md_gl_mol_set_chain_gaps(md_gl_mol_t handle, const md_system_t* sys) {
    molecule_t* mol = mol_lookup(handle.id);
    if (!mol) {
        MD_LOG_ERROR("Invalid molecule");
        return;
    }
    gl_buffer_conditional_delete(&mol->buffer[GL_BUFFER_CHAIN_GAP]);
    gl_buffer_conditional_delete(&mol->buffer[GL_BUFFER_CHAIN_GAP_COMP]);
    mol->gap_trace_count = 0;
    mol->gap_link_count  = 0;
    if (!sys || !sys->component.flags || sys->component.count < 2) return;
    if (sys->atom.count != mol->atom_count) {
        MD_LOG_ERROR("The system is not the one of the molecule: it has %zu atoms, the molecule %u", sys->atom.count, mol->atom_count);
        return;
    }

    const size_t comp_count = MIN(sys->component.count, md_array_size(sys->component.flags));
    uint32_t break_count = 0;
    for (size_t i = 0; i + 1 < comp_count; ++i) {
        if (sys->component.flags[i] & MD_COMPONENT_FLAG_BREAK_AFTER) ++break_count;
    }
    if (break_count == 0) return;

    md_temp_scope_t temp = md_temp_begin();

    // The segment of each component in the protein backbone (-1 for none) and the range of each segment
    const md_protein_backbone_data_t* bb = &sys->protein_backbone;
    int32_t* comp_seg = md_temp_alloc(temp, comp_count * sizeof(int32_t));
    MEMSET(comp_seg, -1, comp_count * sizeof(int32_t));
    uint32_t* seg_beg = NULL;
    uint32_t* seg_end = NULL;
    if (bb->range.count > 0 && bb->range.offset && bb->segment.count > 0 && bb->segment.atoms && bb->segment.comp_idx) {
        seg_beg = md_temp_alloc(temp, bb->segment.count * sizeof(uint32_t));
        seg_end = md_temp_alloc(temp, bb->segment.count * sizeof(uint32_t));
        for (size_t r = 0; r < bb->range.count; ++r) {
            const uint32_t beg = bb->range.offset[r];
            const uint32_t end = (uint32_t)MIN(bb->range.offset[r + 1], bb->segment.count);
            for (uint32_t j = beg; j < end; ++j) {
                seg_beg[j] = beg;
                seg_end[j] = end;
                const int32_t c = bb->segment.comp_idx[j];
                if (0 <= c && (size_t)c < comp_count) comp_seg[c] = (int32_t)j;
            }
        }
    }

    uint32_t* trace      = md_temp_alloc(temp, 2 * (size_t)break_count * 4 * sizeof(uint32_t));
    uint32_t* trace_comp = md_temp_alloc(temp, 2 * (size_t)break_count * sizeof(uint32_t));
    uint32_t* link       = md_temp_alloc(temp, (size_t)break_count * 4 * sizeof(uint32_t));
    uint32_t* link_comp  = md_temp_alloc(temp, (size_t)break_count * sizeof(uint32_t));
    uint32_t trace_count = 0;
    uint32_t link_count  = 0;

    for (size_t i = 0; i + 1 < comp_count; ++i) {
        if (!(sys->component.flags[i] & MD_COMPONENT_FLAG_BREAK_AFTER)) continue;
        md_atom_idx_t atoms[2];
        if (!md_util_chain_break_atoms(atoms, sys, i)) continue;
        const uint32_t la = (uint32_t)atoms[0];
        const uint32_t lb = (uint32_t)atoms[1];
        if (la >= mol->atom_count || lb >= mol->atom_count) continue;

        uint32_t* e = link + 4 * link_count;
        e[0] = la;
        e[1] = la;
        e[2] = lb;
        e[3] = lb;
        link_comp[link_count++] = (uint32_t)i;

        // The cartoon: a protein chain with a cartoon on at least one side
        const int32_t sa = comp_seg[i];
        const int32_t sb = comp_seg[i + 1];
        if (sa < 0 && sb < 0) continue;
        uint32_t a  = la;
        uint32_t na = la;
        uint32_t b  = lb;
        uint32_t nb = lb;
        if (sa >= 0) {
            // The end of a range has the control point before it as its neighbour. A segment the range goes on past (a
            // break the backbone did not end the range at) or one alone is its own.
            a  = (uint32_t)bb->segment.atoms[sa].ca;
            na = ((uint32_t)sa + 1 == seg_end[sa] && (uint32_t)sa > seg_beg[sa]) ? (uint32_t)bb->segment.atoms[sa - 1].ca : a;
        }
        if (sb >= 0) {
            b  = (uint32_t)bb->segment.atoms[sb].ca;
            nb = ((uint32_t)sb == seg_beg[sb] && (uint32_t)sb + 1 < seg_end[sb]) ? (uint32_t)bb->segment.atoms[sb + 1].ca : b;
        }
        if (a >= mol->atom_count || na >= mol->atom_count || b >= mol->atom_count || nb >= mol->atom_count) continue;
        e = trace + 4 * trace_count;
        e[0] = a;
        e[1] = na;
        e[2] = b;
        e[3] = nb;
        trace_comp[trace_count++] = (uint32_t)i;
    }

    const uint32_t total = trace_count + link_count;
    if (total > 0) {
        // The links after the traces
        MEMCPY(trace + 4 * trace_count, link, link_count * 4 * sizeof(uint32_t));
        MEMCPY(trace_comp + trace_count, link_comp, link_count * sizeof(uint32_t));
        mol->buffer[GL_BUFFER_CHAIN_GAP]      = gl_buffer_create(total * 4 * sizeof(uint32_t), trace, GL_STATIC_DRAW);
        mol->buffer[GL_BUFFER_CHAIN_GAP_COMP] = gl_buffer_create(total * sizeof(uint32_t), trace_comp, GL_STATIC_DRAW);
        mol->gap_trace_count = trace_count;
        mol->gap_link_count  = link_count;
    }

    md_temp_end(temp);
}

static const uint32_t backbone_buffers[] = {
    GL_BUFFER_BACKBONE_DATA,
    GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE,
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA,
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT,
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT_PREV,
    GL_BUFFER_BACKBONE_NEIGHBOR,
    GL_BUFFER_BACKBONE_RING_DATA,
    GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX,
};

static void mol_clear_backbone(molecule_t* mol) {
    for (size_t i = 0; i < ARRAY_SIZE(backbone_buffers); ++i) {
        gl_buffer_conditional_delete(&mol->buffer[backbone_buffers[i]]);
    }
    backbone_bounds_free(&mol->backbone_bounds);
    mol->backbone_count = 0;
    mol->backbone_control_point_index_count = 0;
    mol->backbone_orient_history = false;
    mol->flags &= ~MOL_FLAG_HAS_BACKBONE;
    mol->tubes = false;
    mol->tube_max_extent = 0.0f;
}

// The backbone of a molecule (ribbons, cartoon, tubes): cp_count control points in chain_count chains, chain c the control
// points [chain_offset[c], chain_offset[c + 1]) (from chain_offset[0] = 0 to chain_offset[chain_count] = cp_count).
//  data:        [cp_count] what the extraction of each control point reads (compute_spline_extract.vert), its bb_seg_idx
//               the control point itself (which is what it picks as)
//  ss:          [cp_count] its secondary structure (pack_gl_secondary_structure), and the radius of a tube in the w
//  bounds_atom: [cp_count] the atom of the molecule where it is, for the culling and level of detail of the chains
//               (backbone_bounds_t): set positions after (md_gl_mol_set_atom_position, backbone_bounds_update)
// Replaces the previous backbone. The orientation starts over (md_gl_mol_reset_backbone_history).
static void mol_set_backbone(molecule_t* mol, uint32_t cp_count, uint32_t chain_count, const uint32_t* chain_offset, const gl_backbone_data_t* data, const uint32_t* ss, const uint32_t* bounds_atom) {
    mol_clear_backbone(mol);
    if (cp_count == 0 || chain_count == 0) return;
    ASSERT(chain_offset[0] == 0 && chain_offset[chain_count] == cp_count);

    md_temp_scope_t temp = md_temp_begin();
    // The control points of each chain with the first and the last duplicated for adjacency, and a primitive restart
    // after (LINE_STRIP_ADJACENCY, legacy and debugging paths)
    uint32_t* cp_index = md_temp_alloc(temp, (cp_count + 3 * chain_count) * sizeof(uint32_t));
    gl_neighbor_idx_t* neighbor = md_temp_alloc(temp, cp_count * sizeof(gl_neighbor_idx_t));
    uint32_t index_count = 0;
    for (uint32_t c = 0; c < chain_count; ++c) {
        const uint32_t beg = chain_offset[c];
        const uint32_t end = chain_offset[c + 1];
        if (beg == end) continue;
        cp_index[index_count++] = beg;
        for (uint32_t j = beg; j < end; ++j) {
            cp_index[index_count++] = j;
            const uint32_t p1 = (j > beg) ? (j - 1) : j;
            const uint32_t n1 = (j + 1 < end) ? (j + 1) : j;
            neighbor[j].idx[0] = p1;
            neighbor[j].idx[1] = n1;
            neighbor[j].idx[2] = (p1 > beg) ? (p1 - 1) : p1;
            neighbor[j].idx[3] = (n1 + 1 < end) ? (n1 + 1) : n1;
        }
        cp_index[index_count++] = end - 1;
        cp_index[index_count++] = 0xFFFFFFFFU;
    }

    mol->buffer[GL_BUFFER_BACKBONE_DATA]                         = gl_buffer_create(cp_count * sizeof(gl_backbone_data_t), data, GL_STATIC_DRAW);
    mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE]          = gl_buffer_create(cp_count * sizeof(uint32_t), ss, GL_DYNAMIC_DRAW);
    mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA]           = gl_buffer_create(cp_count * sizeof(gl_control_point_t), NULL, GL_DYNAMIC_COPY);
    mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT]      = gl_buffer_create(cp_count * sizeof(gl_control_point_t), NULL, GL_DYNAMIC_COPY);
    mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT_PREV] = gl_buffer_create(cp_count * sizeof(gl_control_point_t), NULL, GL_DYNAMIC_COPY);
    mol->buffer[GL_BUFFER_BACKBONE_NEIGHBOR]                     = gl_buffer_create(cp_count * sizeof(gl_neighbor_idx_t), neighbor, GL_STATIC_DRAW);
    mol->buffer[GL_BUFFER_BACKBONE_RING_DATA]                    = gl_buffer_create(cp_count * MD_GL_BACKBONE_SEGMENT_COUNT * sizeof(gl_backbone_ring_t), NULL, GL_DYNAMIC_COPY);
    mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX]          = gl_buffer_create(index_count * sizeof(uint32_t), cp_index, GL_STATIC_DRAW);

    mol->backbone_count = cp_count;
    mol->backbone_control_point_index_count = index_count;
    mol->flags |= MOL_FLAG_HAS_BACKBONE;

    backbone_bounds_t* bounds = &mol->backbone_bounds;
    md_allocator_i* alloc = md_get_heap_allocator();
    bounds->chain_count  = chain_count;
    bounds->chain_offset = md_alloc(alloc, (chain_count + 1) * sizeof(uint32_t));
    bounds->chain_sphere = md_alloc(alloc, chain_count * sizeof(vec4_t));
    bounds->ca_atom_idx  = md_alloc(alloc, cp_count * sizeof(uint32_t));
    bounds->ca_xyz       = md_alloc(alloc, cp_count * sizeof(vec3_t));
    MEMCPY(bounds->chain_offset, chain_offset, (chain_count + 1) * sizeof(uint32_t));
    MEMCPY(bounds->ca_atom_idx, bounds_atom, cp_count * sizeof(uint32_t));
    MEMSET(bounds->chain_sphere, 0, chain_count * sizeof(vec4_t));
    MEMSET(bounds->ca_xyz, 0, cp_count * sizeof(vec3_t));
    bounds->valid = false;

    md_temp_end(temp);
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
            const uint32_t chain_count = (uint32_t)sys->protein_backbone.range.count;
            const uint32_t backbone_count = sys->protein_backbone.range.offset[chain_count];
            // Control point idx is what the backbone shaders pick by (backbone_base + idx), which is only the segment
            // index while the ranges cover the segments in order from 0
            ASSERT(sys->protein_backbone.range.offset[0] == 0);

            gl_backbone_data_t* backbone_data = md_temp_alloc(temp_scope, backbone_count * sizeof(gl_backbone_data_t));
            uint32_t* secondary_structure     = md_temp_alloc(temp_scope, backbone_count * sizeof(uint32_t));
            uint32_t* ca_atom_idx             = md_temp_alloc(temp_scope, backbone_count * sizeof(uint32_t));

            // The secondary structure depends on the frame, and is not the system's: coil until the caller sets it
            // (md_gl_mol_set_backbone_secondary_structure)
            const uint32_t ss_coil = pack_gl_secondary_structure((md_gl_secondary_structure_t){0});

            for (uint32_t i = 0; i < chain_count; ++i) {
                uint32_t beg = sys->protein_backbone.range.offset[i];
                uint32_t end = sys->protein_backbone.range.offset[i+1];
                for (uint32_t j = beg; j < end; ++j) {
                    backbone_data[j].bb_seg_idx = j;
                    // A trace (coarse grained) has its bead in ca and no C or O: the bead stands in for both. The
                    // orientation comes from the frame of the trace itself; only where the trace runs straight does
                    // it fall back on the carbonyl, which is then a fixed direction (compute_spline_extract.vert)
                    const md_amino_acid_atoms_t atoms = sys->protein_backbone.segment.atoms[j];
                    backbone_data[j].ca     = (uint32_t)atoms.ca;
                    backbone_data[j].c      = (uint32_t)(atoms.c >= 0 ? atoms.c : atoms.ca);
                    backbone_data[j].o      = (uint32_t)(atoms.o >= 0 ? atoms.o : atoms.ca);
                    backbone_data[j].flags  = (j == beg ? BACKBONE_FLAG_CHAIN_BEG : 0) | (j == end - 1 ? BACKBONE_FLAG_CHAIN_END : 0);
                    secondary_structure[j] = ss_coil;
                    ca_atom_idx[j] = (uint32_t)atoms.ca;
                }
            }

            mol_set_backbone(gl_mol, backbone_count, chain_count, sys->protein_backbone.range.offset, backbone_data, secondary_structure, ca_atom_idx);
            if (md_system_state_has_coords(&sys->reference)) {
                backbone_bounds_update(&gl_mol->backbone_bounds, 0, gl_mol->atom_count, sys->reference.xyz);
            }
        }

        gl_mol->bond_count = (uint32_t)sys->bond.count;
        gl_mol->buffer[GL_BUFFER_BOND_ATOM_INDICES] = gl_buffer_create(gl_mol->bond_count * sizeof(uint32_t) * 2, NULL, GL_DYNAMIC_COPY);

        if (sys->bond.pairs) {
            md_gl_mol_set_bonds(handle, 0, gl_mol->bond_count, sys->bond.pairs, sizeof(md_atom_pair_t));
            md_gl_mol_set_bond_orders(handle, sys);
        }
        md_gl_mol_set_chain_gaps(handle, sys);
    }

done:
    md_temp_end(temp_scope);
    return handle;
}

md_gl_mol_t md_gl_mol_create_particles(uint32_t count) {
    md_gl_mol_t handle = {0};
    if (count == 0) {
        MD_LOG_ERROR("A molecule of no particles");
        return handle;
    }
    uint32_t id = md_handle_pool_alloc_slot(&ctx.molecule_pool);
    if (id == MD_HANDLE_INVALID_ID) {
        MD_LOG_ERROR("Fatal error: ran out of slots!");
        return handle;
    }
    molecule_t* gl_mol = ctx.molecules + md_handle_index(id);
    MEMSET(gl_mol, 0, sizeof(molecule_t));
    gl_mol->id = id;
    handle.id = id;

    gl_mol->atom_count = count;
    gl_mol->buffer[GL_BUFFER_ATOM_POSITION]      = gl_buffer_create(count * sizeof(float) * 3, NULL, GL_DYNAMIC_DRAW);
    gl_mol->buffer[GL_BUFFER_ATOM_POSITION_PREV] = gl_buffer_create(count * sizeof(float) * 3, NULL, GL_DYNAMIC_COPY);
    gl_mol->buffer[GL_BUFFER_ATOM_VELOCITY]      = gl_buffer_create(count * sizeof(float) * 3, NULL, GL_DYNAMIC_COPY);
    gl_mol->buffer[GL_BUFFER_ATOM_RADIUS]        = gl_buffer_create(count * sizeof(float) * 1, NULL, GL_STATIC_DRAW);
    gl_mol->buffer[GL_BUFFER_ATOM_FLAGS]         = gl_buffer_create(count * sizeof(uint8_t),   NULL, GL_STATIC_DRAW);

    md_temp_scope_t temp = md_temp_begin();
    void* zero = md_temp_alloc_zero(temp, count * sizeof(float) * 3);
    gl_buffer_set_sub_data(gl_mol->buffer[GL_BUFFER_ATOM_POSITION],      0, count * sizeof(float) * 3, zero);
    gl_buffer_set_sub_data(gl_mol->buffer[GL_BUFFER_ATOM_POSITION_PREV], 0, count * sizeof(float) * 3, zero);
    gl_buffer_set_sub_data(gl_mol->buffer[GL_BUFFER_ATOM_RADIUS],        0, count * sizeof(float),     zero);
    gl_buffer_set_sub_data(gl_mol->buffer[GL_BUFFER_ATOM_FLAGS],         0, count * sizeof(uint8_t),   zero);
    md_temp_end(temp);
    md_gl_mol_zero_velocity(handle);
    return handle;
}

// ### DERIVED PARTICLES ###

static void derived_free(molecule_t* mol) {
    md_allocator_i* alloc = md_get_heap_allocator();
    if (mol->derived.data) {
        md_free(alloc, mol->derived.data,   3 * mol->atom_count * sizeof(uint32_t));
        md_free(alloc, mol->derived.atom,   MAX(mol->derived.atom_count, 1) * sizeof(uint32_t));
        md_free(alloc, mol->derived.weight, MAX(mol->derived.atom_count, 1) * sizeof(float));
    }
    MEMSET(&mol->derived, 0, sizeof(mol->derived));
    gl_buffer_conditional_delete(&mol->buffer[GL_BUFFER_DERIVED_DATA]);
    gl_buffer_conditional_delete(&mol->buffer[GL_BUFFER_DERIVED_ATOM]);
    gl_buffer_conditional_delete(&mol->buffer[GL_BUFFER_DERIVED_WEIGHT]);
}

void md_gl_mol_set_derived(md_gl_mol_t handle, md_gl_mol_t source_handle, const uint32_t* offset, const uint32_t* atom, const float* weight, const uint32_t* ref) {
    molecule_t* mol = mol_lookup(handle.id);
    const molecule_t* src = mol_lookup(source_handle.id);
    if (!mol || !src) {
        MD_LOG_ERROR("Invalid molecule");
        return;
    }
    if (mol == src) {
        MD_LOG_ERROR("A molecule cannot be derived from itself");
        return;
    }
    const uint32_t n = mol->atom_count;
    if (!offset || !ref || (offset[n] > offset[0] && (!atom || !weight))) {
        MD_LOG_ERROR("Missing one or more arguments");
        return;
    }
    for (uint32_t i = 0; i < n; ++i) {
        if (offset[i + 1] < offset[i]) {
            MD_LOG_ERROR("The offsets of the atoms of the derived particles are not in order");
            return;
        }
    }

    derived_free(mol);

    // Atoms outside of the source are left out: the shader reads them without bounds checks
    uint32_t total = 0;
    for (uint32_t j = offset[0]; j < offset[n]; ++j) {
        total += atom[j] < src->atom_count ? 1 : 0;
    }
    const uint32_t dropped = (offset[n] - offset[0]) - total;

    md_allocator_i* alloc = md_get_heap_allocator();
    uint32_t* data       = md_alloc(alloc, 3 * n * sizeof(uint32_t));
    uint32_t* atom_out   = md_alloc(alloc, MAX(total, 1) * sizeof(uint32_t));
    float*    weight_out = md_alloc(alloc, MAX(total, 1) * sizeof(float));
    atom_out[0] = 0;
    weight_out[0] = 0.0f;

    uint32_t len = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t beg = len;
        for (uint32_t j = offset[i]; j < offset[i + 1]; ++j) {
            if (atom[j] < src->atom_count) {
                atom_out[len]   = atom[j];
                weight_out[len] = weight[j];
                ++len;
            }
        }
        data[3 * i + 0] = beg;
        data[3 * i + 1] = len - beg;
        // Without a reference of its own, its first atom: a particle of no atoms is at the origin whatever it is
        data[3 * i + 2] = ref[i] < src->atom_count ? ref[i] : (len > beg ? atom_out[beg] : 0);
    }
    ASSERT(len == total);
    if (dropped) {
        MD_LOG_ERROR("%u atoms of derived particles were outside of the source (%u atoms) and were left out", dropped, src->atom_count);
    }

    mol->buffer[GL_BUFFER_DERIVED_DATA]   = gl_buffer_create(3 * n * sizeof(uint32_t),      data,       GL_STATIC_DRAW);
    mol->buffer[GL_BUFFER_DERIVED_ATOM]   = gl_buffer_create(MAX(total, 1) * sizeof(uint32_t), atom_out,   GL_STATIC_DRAW);
    mol->buffer[GL_BUFFER_DERIVED_WEIGHT] = gl_buffer_create(MAX(total, 1) * sizeof(float),    weight_out, GL_STATIC_DRAW);
    mol->derived.source_id         = source_handle.id;
    mol->derived.source_atom_count = src->atom_count;
    mol->derived.atom_count        = total;
    mol->derived.data              = data;
    mol->derived.atom              = atom_out;
    mol->derived.weight            = weight_out;
}

// The periodic images of compute_derived.vert: the basis vectors of the cell as the columns of A, its inverse I and which
// of them are periodic (mask). Without a cell, or a cell of no volume, there are none.
static void derived_cell(mat3_t* A, mat3_t* I, vec3_t* mask, const md_unitcell_t* cell) {
    *A = mat3_ident();
    *I = mat3_ident();
    *mask = (vec3_t){0, 0, 0};
    if (!cell || !(cell->flags & MD_UNITCELL_PBC_ALL)) return;
    mat3_t M;
    md_unitcell_A_extract_float(M.elem, cell);
    const float det = mat3_determinant(M);
    if (!(fabsf(det) > 1.0e-6f)) return;
    *A = M;
    *I = mat3_inverse(M);
    *mask = (vec3_t){
        (cell->flags & MD_UNITCELL_PBC_X) ? 1.0f : 0.0f,
        (cell->flags & MD_UNITCELL_PBC_Y) ? 1.0f : 0.0f,
        (cell->flags & MD_UNITCELL_PBC_Z) ? 1.0f : 0.0f,
    };
}

// Where the control points of the tubes are, from the positions of the atoms of the source on the CPU: the particles
// they are, as compute_derived.vert places them
static void derived_bounds_update(molecule_t* mol, const vec3_t* xyz, const mat3_t* A, const mat3_t* I, vec3_t mask) {
    backbone_bounds_t* bounds = &mol->backbone_bounds;
    if (!bounds->chain_count || !mol->derived.data) return;
    const uint32_t cp_count = bounds->chain_offset[bounds->chain_count];
    for (uint32_t i = 0; i < cp_count; ++i) {
        const uint32_t p = bounds->ca_atom_idx[i];
        if (p >= mol->atom_count) continue;
        const uint32_t* d = mol->derived.data + 3 * p;
        const vec3_t ref = xyz[d[2]];
        vec3_t sum = {0, 0, 0};
        float sum_w = 0.0f;
        for (uint32_t k = d[0]; k < d[0] + d[1]; ++k) {
            const float w = mol->derived.weight[k];
            vec3_t dx = vec3_sub(xyz[mol->derived.atom[k]], ref);
            vec3_t f = mat3_mul_vec3(*I, dx);
            f = (vec3_t){floorf(f.x + 0.5f) * mask.x, floorf(f.y + 0.5f) * mask.y, floorf(f.z + 0.5f) * mask.z};
            dx = vec3_sub(dx, mat3_mul_vec3(*A, f));
            sum = vec3_add(sum, vec3_mul1(dx, w));
            sum_w += w;
        }
        bounds->ca_xyz[i] = vec3_add(vec3_mul1(ref, sum_w), sum);
    }
    backbone_bounds_update_spheres(bounds);
}

void md_gl_mol_compute_derived(md_gl_mol_t handle, const md_unitcell_t* cell, const vec3_t* source_xyz) {
    if (!validate_context()) {
        return;
    }
    molecule_t* mol = mol_lookup(handle.id);
    if (!mol) {
        return;
    }
    if (!mol->buffer[GL_BUFFER_DERIVED_DATA].id) {
        MD_LOG_ERROR("The particles are not derived from a source (md_gl_mol_set_derived)");
        return;
    }
    const molecule_t* src = mol_lookup(mol->derived.source_id);
    if (!src || src->atom_count != mol->derived.source_atom_count) {
        MD_LOG_ERROR("The source of the derived particles is gone");
        return;
    }

    mat3_t A, I;
    vec3_t mask;
    derived_cell(&A, &I, &mask, cell);

    glEnable(GL_RASTERIZER_DISCARD);
    glBindVertexArray(ctx.vao);

    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_DERIVED_DATA].id);
    glEnableVertexAttribArray(0);
    glVertexAttribIPointer(0, 3, GL_UNSIGNED_INT, 3 * sizeof(uint32_t), 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, src->buffer[GL_BUFFER_ATOM_POSITION].id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, src->buffer[GL_BUFFER_ATOM_VELOCITY].id);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_2].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, mol->buffer[GL_BUFFER_DERIVED_ATOM].id);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_3].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, mol->buffer[GL_BUFFER_DERIVED_WEIGHT].id);

    const GLuint program = ctx.program[GL_PROGRAM_COMPUTE_DERIVED].id;
    glUseProgram(program);
    glUniform1i(glGetUniformLocation(program, "u_buf_atom_pos"), 0);
    glUniform1i(glGetUniformLocation(program, "u_buf_atom_vel"), 1);
    glUniform1i(glGetUniformLocation(program, "u_buf_derived_atom"), 2);
    glUniform1i(glGetUniformLocation(program, "u_buf_derived_weight"), 3);
    glUniformMatrix3fv(glGetUniformLocation(program, "u_cell"), 1, GL_FALSE, &A.elem[0][0]);
    glUniformMatrix3fv(glGetUniformLocation(program, "u_inv_cell"), 1, GL_FALSE, &I.elem[0][0]);
    glUniform3fv(glGetUniformLocation(program, "u_pbc_mask"), 1, mask.elem);

    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 1, mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    glBeginTransformFeedback(GL_POINTS);
    glDrawArrays(GL_POINTS, 0, mol->atom_count);
    glEndTransformFeedback();
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, 0);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 1, 0);
    glUseProgram(0);

    glDisableVertexAttribArray(0);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(0);
    glDisable(GL_RASTERIZER_DISCARD);

    if (source_xyz) {
        derived_bounds_update(mol, source_xyz, &A, &I, mask);
    }
}

// ### TUBES ###

static inline uint32_t pack_tube_length(float v) {
    return (uint32_t)(CLAMP(v, 0.0f, TUBE_RADIUS_MAX) / TUBE_RADIUS_MAX * 255.0f + 0.5f);
}

void md_gl_mol_set_tubes(md_gl_mol_t handle, const md_gl_tubes_t* tubes) {
    molecule_t* mol = mol_lookup(handle.id);
    if (!mol) {
        MD_LOG_ERROR("Invalid molecule");
        return;
    }
    const uint32_t count = tubes ? tubes->count : 0;
    const uint32_t* offset   = count ? tubes->offset : NULL;
    const uint32_t* particle = count ? tubes->particle : NULL;
    const float*    radius   = count ? tubes->radius : NULL;
    const float*    flat     = count ? tubes->flat : NULL;
    const uint32_t* face     = count ? tubes->face : NULL;
    const uint8_t*  through  = count ? tubes->through : NULL;
    if (count > 0 && (!offset || !particle || !radius)) {
        MD_LOG_ERROR("Missing one or more arguments");
        return;
    }
    if (mol->buffer[GL_BUFFER_BACKBONE_DATA].id && !mol->tubes) {
        MD_LOG_ERROR("The molecule has a protein backbone: tubes are for molecules of particles (md_gl_mol_create_particles)");
        return;
    }
    for (uint32_t t = 0; t < count; ++t) {
        if (offset[t + 1] < offset[t]) {
            MD_LOG_ERROR("The offsets of the tubes are not in order");
            return;
        }
    }
    const uint32_t base = count ? offset[0] : 0;
    const uint32_t cp_count = count ? offset[count] - base : 0;
    for (uint32_t j = 0; j < cp_count; ++j) {
        if (particle[base + j] >= mol->atom_count || (face && face[base + j] != UINT32_MAX && face[base + j] >= mol->atom_count)) {
            MD_LOG_ERROR("A control point of a tube is not a particle of the molecule (%u particles)", mol->atom_count);
            return;
        }
    }

    md_temp_scope_t temp = md_temp_begin();
    uint32_t* chain_offset    = md_temp_alloc(temp, (count + 1) * sizeof(uint32_t));
    gl_backbone_data_t* data  = md_temp_alloc(temp, MAX(cp_count, 1) * sizeof(gl_backbone_data_t));
    uint32_t* ss              = md_temp_alloc(temp, MAX(cp_count, 1) * sizeof(uint32_t));
    uint32_t* bounds_particle = md_temp_alloc(temp, MAX(cp_count, 1) * sizeof(uint32_t));

    // Coil, with the half length of the flat sides in the y and the radius in the w (compute_spline_extract.vert)
    const uint32_t ss_coil = pack_gl_secondary_structure((md_gl_secondary_structure_t){0}) & 0xFFu;
    float max_extent = 0.0f;
    for (uint32_t t = 0; t < count; ++t) {
        const uint32_t beg = offset[t] - base;
        const uint32_t end = offset[t + 1] - base;
        chain_offset[t] = beg;
        for (uint32_t j = beg; j < end; ++j) {
            const uint32_t p = particle[base + j];
            const uint32_t f = face ? face[base + j] : UINT32_MAX;
            const float    r = CLAMP(radius[base + j], 0.0f, TUBE_RADIUS_MAX);
            const float    l = flat ? CLAMP(flat[base + j], 0.0f, TUBE_RADIUS_MAX) : 0.0f;
            data[j] = (gl_backbone_data_t){
                .bb_seg_idx = j,
                .ca         = p,
                .c          = p,
                .o          = f != UINT32_MAX ? f : p,
                .flags      = (j == beg ? BACKBONE_FLAG_CHAIN_BEG : 0) | (j == end - 1 ? BACKBONE_FLAG_CHAIN_END : 0) |
                              (f != UINT32_MAX ? BACKBONE_FLAG_FIXED_FACE : 0) | ((through && through[t]) ? BACKBONE_FLAG_THROUGH : 0),
            };
            ss[j] = ss_coil | (pack_tube_length(l) << 8) | (pack_tube_length(r) << 24);
            bounds_particle[j] = p;
            max_extent = MAX(max_extent, l + r);
        }
    }
    chain_offset[count] = cp_count;

    mol_set_backbone(mol, cp_count, count, chain_offset, data, ss, bounds_particle);
    mol->tubes = cp_count > 0;
    mol->tube_max_extent = max_extent;
    md_temp_end(temp);
}

// ### SLABS ###

#define SLAB_MAX_CORNERS 7
#define SLAB_ENTRY_U32   8

void md_gl_mol_set_slabs(md_gl_mol_t handle, uint32_t count, const uint32_t* offset, const uint32_t* particle) {
    molecule_t* mol = mol_lookup(handle.id);
    if (!mol) {
        MD_LOG_ERROR("Invalid molecule");
        return;
    }
    if (count > 0 && (!offset || !particle)) {
        MD_LOG_ERROR("Missing one or more arguments");
        return;
    }
    for (uint32_t i = 0; i < count; ++i) {
        const uint32_t n = offset[i + 1] - offset[i];
        if (offset[i + 1] < offset[i] || n < 2 || n > SLAB_MAX_CORNERS) {
            MD_LOG_ERROR("A slab has %u corners, it takes 2 to %u", offset[i + 1] < offset[i] ? 0 : n, SLAB_MAX_CORNERS);
            return;
        }
        for (uint32_t j = offset[i]; j < offset[i + 1]; ++j) {
            if (particle[j] >= mol->atom_count) {
                MD_LOG_ERROR("A corner of a slab is not a particle of the molecule (%u particles)", mol->atom_count);
                return;
            }
        }
    }

    gl_buffer_conditional_delete(&mol->buffer[GL_BUFFER_SLAB]);
    mol->slab_count = 0;
    if (count == 0) return;

    md_temp_scope_t temp = md_temp_begin();
    uint32_t* entry = md_temp_alloc_zero(temp, count * SLAB_ENTRY_U32 * sizeof(uint32_t));
    for (uint32_t i = 0; i < count; ++i) {
        uint32_t* e = entry + i * SLAB_ENTRY_U32;
        e[0] = offset[i + 1] - offset[i];
        for (uint32_t j = offset[i]; j < offset[i + 1]; ++j) {
            e[1 + j - offset[i]] = particle[j];
        }
    }
    mol->buffer[GL_BUFFER_SLAB] = gl_buffer_create(count * SLAB_ENTRY_U32 * sizeof(uint32_t), entry, GL_STATIC_DRAW);
    mol->slab_count = count;
    md_temp_end(temp);
}

void md_gl_mol_destroy(md_gl_mol_t handle) {
    if (handle.id == MD_HANDLE_INVALID_ID) return;
    molecule_t* mol = mol_lookup(handle.id);
    if (mol) {
        derived_free(mol);
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
        const uint32_t color = (255u << 24) | (127u << 16) | (127u << 8) | (127u << 0);
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
        gl_buffer_conditional_delete(&rep->weak_bond_atoms);
        gl_buffer_conditional_delete(&rep->weak_bond_weight);
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

void md_gl_rep_set_weak_bonds(md_gl_rep_t handle, uint32_t count, const uint32_t* atom_a, const uint32_t* atom_b, const float* weight) {
    representation_t* rep = rep_lookup(handle.id);
    if (rep == NULL) {
        MD_LOG_ERROR("representation was invalid");
        return;
    }
    const molecule_t* mol = mol_lookup(rep->mol_id);
    if (mol == NULL) {
        MD_LOG_ERROR("Representation's molecule is invalid");
        return;
    }
    if (count > 0 && (atom_a == NULL || atom_b == NULL)) {
        MD_LOG_ERROR("Missing atom indices of the weak bonds");
        return;
    }

    rep->weak_bond_count = 0;
    if (count == 0) {
        return;
    }

    if (count > rep->weak_bond_capacity) {
        // The bonds typically change every frame: grow with some slack
        const uint32_t capacity = MAX(count, rep->weak_bond_capacity + rep->weak_bond_capacity / 2);
        gl_buffer_conditional_delete(&rep->weak_bond_atoms);
        gl_buffer_conditional_delete(&rep->weak_bond_weight);
        rep->weak_bond_atoms  = gl_buffer_create(capacity * sizeof(uint32_t) * 2, NULL, GL_DYNAMIC_DRAW);
        rep->weak_bond_weight = gl_buffer_create(capacity * sizeof(float),        NULL, GL_DYNAMIC_DRAW);
        rep->weak_bond_capacity = capacity;
    }

    // Bonds with an atom outside of the molecule are left out (the shaders fetch the atoms without bounds checks)
    uint32_t written = 0;
    glBindBuffer(GL_ARRAY_BUFFER, rep->weak_bond_atoms.id);
    uint32_t* atoms = (uint32_t*)glMapBufferRange(GL_ARRAY_BUFFER, 0, count * sizeof(uint32_t) * 2, GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);
    if (atoms) {
        for (uint32_t i = 0; i < count; ++i) {
            if (atom_a[i] < mol->atom_count && atom_b[i] < mol->atom_count) {
                atoms[written * 2 + 0] = atom_a[i];
                atoms[written * 2 + 1] = atom_b[i];
                ++written;
            }
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }

    glBindBuffer(GL_ARRAY_BUFFER, rep->weak_bond_weight.id);
    float* w = (float*)glMapBufferRange(GL_ARRAY_BUFFER, 0, count * sizeof(float), GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);
    if (w) {
        uint32_t j = 0;
        for (uint32_t i = 0; i < count; ++i) {
            if (atom_a[i] < mol->atom_count && atom_b[i] < mol->atom_count) {
                w[j++] = weight ? weight[i] : 1.0f;
            }
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    if (!atoms || !w) {
        MD_LOG_ERROR("Failed to map the buffers of the weak bonds");
        return;
    }
    if (written < count) {
        MD_LOG_ERROR("%u weak bonds had atom indices outside of the molecule (%u atoms) and were left out", count - written, mol->atom_count);
    }
    rep->weak_bond_count = written;
}

static bool compute_spline(molecule_t* mol);

typedef struct {
    gl_buffer_t atoms;      // u32[2] per bond
    gl_buffer_t weight;     // f32 per bond, weak bonds only
    gl_buffer_t order;      // u32 per bond, covalent bonds with bond orders only (GL_BUFFER_BOND_ORDER)
    gl_buffer_t order_ref;  // u32[2] per entry (GL_BUFFER_BOND_ORDER_REF)
    gl_buffer_t partner;    // u32[4] per dashed partner (GL_BUFFER_BOND_PARTNER): the primitives are its entries
    uint32_t    count;      // Primitives: bonds, or the entries of partner
} bond_set_t;

// Parameters of bond.vert and bond.frag, after gl_ubo_base_t
typedef struct {
    float    radius;
    float    max_d2;
    int32_t  mode;
    float    sharpness;
    uint32_t uniform_color;
    uint32_t dash_count;
    float    dash_fill;
    uint32_t bond_orders;   // Show the orders of covalent bonds (bond_set_t order, order_ref)
    float    spread;        // Envelope of the thinner sticks of multiple bonds, relative to the radius
    float    dash_period;   // The gaps of the chains (CHAIN_GAPS): the period of the dashes (Ångström), dash_count unused
    uint32_t gap_index_base;// The gaps of the chains: their picking base (md_gl_picking_offset_t gap_base)
} bond_params_t;

// The parameters of the gaps of the chains (draw_chain_gaps): dashes of the radius, colored as the bonds are, drawn
// at any length
static inline bond_params_t chain_gap_params(float radius, int32_t mode, float sharpness, uint32_t uniform_color, uint32_t gap_index_base) {
    bond_params_t params = {
        .radius = radius,
        .max_d2 = 3.4e38f,
        .mode = mode,
        .sharpness = sharpness,
        .uniform_color = uniform_color,
        .dash_count = 2,
        .dash_fill = CHAIN_GAP_DASH_FILL,
        .dash_period = CHAIN_GAP_DASH_PERIOD,
        .gap_index_base = gap_index_base,
    };
    return params;
}

// Parameters of ring.vert and ring.frag, after gl_ubo_base_t. The picking base of the rings is here rather than in
// gl_ubo_base_t, which every shader declares: it is 16 byte aligned (vec4_t), and a fifth index would pad it.
typedef struct {
    float    radius;        // Of the sticks
    float    max_d2;
    int32_t  mode;
    float    sharpness;
    uint32_t uniform_color;
    float    spread;
    uint32_t ring_index_base;
} ring_params_t;

static bool draw_space_fill(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, float scale);
static bool draw_bonds     (gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, bond_set_t bonds, const bond_params_t* params, const float viewport_size[2]);
static bool draw_bond_partners(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, const bond_params_t* bond_params, const float viewport_size[2]);
static bool draw_chain_gaps(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, bool trace, const bond_params_t* params, const float viewport_size[2]);
static bool draw_aromatic_rings(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, const bond_params_t* bond_params, uint32_t ring_index_base, const float viewport_size[2]);
static bool draw_backbone  (gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, uint32_t profile, const float scale[4], float max_extent, const mat4_t* world_to_clip, float viewport_half_height, bool raycast, const float viewport_size[2]);
static bool draw_slabs     (gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, float radius, float plate, float max_length, uint32_t slab_index_base, const float viewport_size[2]);

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
        if (model_matrix) {
            // The model matrix as it was in the previous frame too (the images of a periodic cell, say). Without it the
            // velocities took view_to_world, which holds its inverse, against the previous view without it: a copy moved
            // by the matrix showed as moving by it, every frame
            ubo_data->view_transform.prev_world_to_clip = mat4_mul(ubo_data->view_transform.prev_world_to_clip, *model_matrix);
        }
        ubo_data->view_transform.curr_view_to_prev_clip = mat4_mul(ubo_data->view_transform.prev_world_to_clip, ubo_data->view_transform.view_to_world);
        extract_jitter_uv(ubo_data->jitter_uv.elem + 0, ubo_data->view_transform.view_to_clip);
        extract_jitter_uv(ubo_data->jitter_uv.elem + 2, *prev_view_to_clip);
    }
    ubo_data->atom_mask = args->atom_mask;
    ubo_data->atom_index_base = args->picking_offset.atom_base;
    ubo_data->bond_index_base = args->picking_offset.bond_base;
    ubo_data->backbone_index_base = args->picking_offset.backbone_base;
}

static inline bool is_backbone_representation_type(md_gl_rep_type_t type) {
    return type == MD_GL_REP_RIBBONS || type == MD_GL_REP_CARTOON || type == MD_GL_REP_TUBES;
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
            MD_LOG_ERROR("Incompatible representation type for molecule: Molecule is missing a backbone (a protein backbone, or tubes)");
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
    const float viewport_size[2] = {(float)viewport[2], (float)viewport[3]};
        
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
        case MD_GL_REP_LICORICE: {
            const bool orders = draw_op->args.licorice.bond_orders && mol->bond_order_count == mol->bond_count;
            const bond_params_t params = {
                .radius = 0.2f * scale * draw_op->args.licorice.radius,
                .max_d2 = max_length * max_length,
                .mode = draw_op->args.licorice.color_mode,
                .sharpness = draw_op->args.licorice.sharpness,
                .uniform_color = draw_op->args.licorice.uniform_color,
                .bond_orders = orders,
                // Within the stick: their ends stay within the sphere at the atom
                .spread = 1.0f,
            };
            const bond_set_t bonds = {mol->buffer[GL_BUFFER_BOND_ATOM_INDICES], {0}, mol->buffer[GL_BUFFER_BOND_ORDER], mol->buffer[GL_BUFFER_BOND_ORDER_REF], {0}, mol->bond_count};
            draw_bonds(shaders->licorice[program_permutation], mol, rep->atom_color, bonds, &params, viewport_size);
            if (orders) {
                draw_bond_partners(shaders->bond_partners[program_permutation], mol, rep->atom_color, &params, viewport_size);
                draw_aromatic_rings(shaders->aromatic_rings[program_permutation], mol, rep->atom_color, &params, args->picking_offset.ring_base, viewport_size);
            }
            if (draw_op->args.licorice.chain_gaps) {
                const bond_params_t gap = chain_gap_params(CHAIN_GAP_RADIUS_LINK * params.radius, params.mode, params.sharpness, params.uniform_color, args->picking_offset.gap_base);
                draw_chain_gaps(shaders->chain_gaps[program_permutation], mol, rep->atom_color, false, &gap, viewport_size);
            }
            break;
        }
        case MD_GL_REP_BALL_AND_STICK: {
            const bool orders = draw_op->args.ball_and_stick.bond_orders && mol->bond_order_count == mol->bond_count;
            const bond_params_t params = {
                .radius = 0.2f * scale * draw_op->args.ball_and_stick.stick_radius,
                .max_d2 = max_length * max_length,
                .mode = draw_op->args.ball_and_stick.color_mode,
                .sharpness = draw_op->args.ball_and_stick.sharpness,
                .uniform_color = draw_op->args.ball_and_stick.uniform_color,
                .bond_orders = orders,
                // The balls are larger than the sticks
                .spread = 1.3f,
            };
            const bond_set_t bonds = {mol->buffer[GL_BUFFER_BOND_ATOM_INDICES], {0}, mol->buffer[GL_BUFFER_BOND_ORDER], mol->buffer[GL_BUFFER_BOND_ORDER_REF], {0}, mol->bond_count};
            draw_bonds(shaders->licorice[program_permutation], mol, rep->atom_color, bonds, &params, viewport_size);
            if (orders) {
                draw_bond_partners(shaders->bond_partners[program_permutation], mol, rep->atom_color, &params, viewport_size);
                draw_aromatic_rings(shaders->aromatic_rings[program_permutation], mol, rep->atom_color, &params, args->picking_offset.ring_base, viewport_size);
            }
            if (draw_op->args.ball_and_stick.chain_gaps) {
                const bond_params_t gap = chain_gap_params(CHAIN_GAP_RADIUS_LINK * params.radius, params.mode, params.sharpness, params.uniform_color, args->picking_offset.gap_base);
                draw_chain_gaps(shaders->chain_gaps[program_permutation], mol, rep->atom_color, false, &gap, viewport_size);
            }
            draw_space_fill(shaders->spacefill[program_permutation], mol, rep->atom_color, 0.2f * scale * draw_op->args.ball_and_stick.ball_scale);
            break;
        }
        case MD_GL_REP_WEAK_BONDS: {
            const bond_params_t params = {
                .radius = scale * draw_op->args.weak_bonds.radius,
                .max_d2 = max_length * max_length,
                .mode = draw_op->args.weak_bonds.color_mode,
                .sharpness = draw_op->args.weak_bonds.sharpness,
                .uniform_color = draw_op->args.weak_bonds.uniform_color,
                .dash_count = MAX(1, draw_op->args.weak_bonds.dash_count),
                .dash_fill = CLAMP(draw_op->args.weak_bonds.dash_fill, 0.0f, 1.0f),
            };
            const bond_set_t bonds = {rep->weak_bond_atoms, rep->weak_bond_weight, {0}, {0}, {0}, rep->weak_bond_count};
            draw_bonds(shaders->weak_bonds[program_permutation], mol, rep->atom_color, bonds, &params, viewport_size);
            break;
        }
        case MD_GL_REP_RIBBONS: {
            // Half width and half thickness of the box profile
            const float profile_scale[4] = { scale * draw_op->args.ribbons.width_scale, scale * draw_op->args.ribbons.thickness_scale * 0.1f, 0.0f, 0.0f };
            const float max_extent = sqrtf(profile_scale[0] * profile_scale[0] + profile_scale[1] * profile_scale[1]);
            draw_backbone(shaders->ribbons[program_permutation], mol, rep->atom_color, BACKBONE_PROFILE_BOX, profile_scale, max_extent, &world_to_clip, viewport_half_height, false, viewport_size);
            if (draw_op->args.ribbons.chain_gaps) {
                // As thick as the gap of a cartoon of the same thickness: the box is 0.1 thick at a scale of 1, the coil 0.2
                const float radius = CHAIN_GAP_RADIUS_TRACE * 0.2f * scale * draw_op->args.ribbons.thickness_scale;
                const bond_params_t gap = chain_gap_params(radius, MD_GL_BOND_MODE_SMOOTH, 0.0f, 0, args->picking_offset.gap_base);
                draw_chain_gaps(shaders->chain_gaps[program_permutation], mol, rep->atom_color, true, &gap, viewport_size);
            }
            break;
        }
        case MD_GL_REP_CARTOON: {
            // Scales of the coil, helix and sheet profiles
            const float profile_scale[4] = { scale * draw_op->args.cartoon.coil_scale, draw_op->args.cartoon.helix_scale, draw_op->args.cartoon.sheet_scale, 0.0f };
            // Largest semi axis of the profiles in backbone.vert
            const float max_extent = MAX(MAX(0.2f * profile_scale[0], 1.2f * profile_scale[1]), 1.5f * profile_scale[2]);
            const bool rc = draw_op->args.cartoon.raycast;
            draw_backbone(rc ? shaders->cartoon_rc[program_permutation] : shaders->cartoon[program_permutation], mol, rep->atom_color, BACKBONE_PROFILE_STADIUM, profile_scale, max_extent, &world_to_clip, viewport_half_height, rc, viewport_size);
            if (draw_op->args.cartoon.chain_gaps) {
                // Thinner than the coil (radius 0.2 at a scale of 1), its colors blended smoothly from end to end
                const float radius = CHAIN_GAP_RADIUS_TRACE * 0.2f * profile_scale[0];
                const bond_params_t gap = chain_gap_params(radius, MD_GL_BOND_MODE_SMOOTH, 0.0f, 0, args->picking_offset.gap_base);
                draw_chain_gaps(shaders->chain_gaps[program_permutation], mol, rep->atom_color, true, &gap, viewport_size);
            }
            break;
        }
        case MD_GL_REP_SLABS: {
            const float radius = scale * draw_op->args.slabs.radius;
            const float plate  = radius * CLAMP(draw_op->args.slabs.plate, 0.0f, 1.0f);
            draw_slabs(shaders->slabs[program_permutation], mol, rep->atom_color, radius, plate, max_length, args->picking_offset.ring_base, viewport_size);
            break;
        }
        case MD_GL_REP_TUBES: {
            // The (l, r) of each control point times the scale, the r of a band (l > 0) times the band scale, and round caps
            const float band = draw_op->args.tubes.band_scale > 0.0f ? draw_op->args.tubes.band_scale : 1.0f;
            const float profile_scale[4] = { scale * draw_op->args.tubes.radius_scale, band, 0.0f, 0.0f };
            const float max_extent = mol->tube_max_extent * profile_scale[0] * MAX(band, 1.0f);
            const bool rc = draw_op->args.tubes.raycast;
            draw_backbone(rc ? shaders->tubes_rc[program_permutation] : shaders->tubes[program_permutation], mol, rep->atom_color, BACKBONE_PROFILE_TUBE, profile_scale, max_extent, &world_to_clip, viewport_half_height, rc, viewport_size);
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

// Bonds (licorice, weak bonds) by vertex pulling (bond.vert): one quad per bond, the atoms and bonds are read from
// texture buffers. Full instances of the quad index buffer, then the rest.
static bool draw_bonds(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, bond_set_t bonds, const bond_params_t* params, const float viewport_size[2]) {
    ASSERT(mol);
    ASSERT(params);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(atom_color.id);

    if (bonds.count == 0 || !bonds.atoms.id || !ctx.bond_quad_index.id) {
        return false;
    }

    bond_params_t p = *params;
    p.bond_orders = (p.bond_orders && bonds.order.id && bonds.order_ref.id) ? 1 : 0;
    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(bond_params_t), &p);

    glBindVertexArray(ctx.vao);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_POSITION].id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_2].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, mol->buffer[GL_BUFFER_ATOM_FLAGS].id);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_3].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color.id);

    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_4].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RG32UI, bonds.atoms.id);

    if (bonds.weight.id) {
        glActiveTexture(GL_TEXTURE5);
        glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_5].id);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, bonds.weight.id);
    }

    if (bonds.partner.id) {
        glActiveTexture(GL_TEXTURE8);
        glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_8].id);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, bonds.partner.id);
    }

    if (p.bond_orders) {
        glActiveTexture(GL_TEXTURE6);
        glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_6].id);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, bonds.order.id);

        glActiveTexture(GL_TEXTURE7);
        glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_7].id);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RG32UI, bonds.order_ref.id);
    }

    glUseProgram(program.id);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_pos"),   0);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_vel"),   1);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_flags"), 2);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_color"), 3);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_bond"),       4);
    if (bonds.weight.id) {
        glUniform1i(glGetUniformLocation(program.id, "u_buf_bond_weight"), 5);
    }
    // Always on units of their own, also when unused: samplers of different types must not share a unit
    glUniform1i(glGetUniformLocation(program.id, "u_buf_bond_order"), 6);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_bond_ref"),   7);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_bond_partner"), 8);
    glUniform1ui(glGetUniformLocation(program.id, "u_quads_per_instance"), BOND_QUADS_PER_INSTANCE);
    glUniform2f(glGetUniformLocation(program.id, "u_viewport_size"), viewport_size[0], viewport_size[1]);
    const GLint bond_offset_loc = glGetUniformLocation(program.id, "u_bond_offset");

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ctx.bond_quad_index.id);

    const uint32_t full = bonds.count / BOND_QUADS_PER_INSTANCE;
    const uint32_t rest = bonds.count % BOND_QUADS_PER_INSTANCE;
    if (full > 0) {
        glUniform1ui(bond_offset_loc, 0);
        glDrawElementsInstanced(GL_TRIANGLES, BOND_QUADS_PER_INSTANCE * 6, GL_UNSIGNED_SHORT, 0, (GLsizei)full);
    }
    if (rest > 0) {
        glUniform1ui(bond_offset_loc, full * BOND_QUADS_PER_INSTANCE);
        glDrawElements(GL_TRIANGLES, (GLsizei)(rest * 6), GL_UNSIGNED_SHORT, 0);
    }

    glUseProgram(0);

    glActiveTexture(GL_TEXTURE0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return true;
}

// The dashed partners of the aromatic and delocalized bonds (md_gl_mol_set_bond_orders), drawn with the bonds of the
// licorice and ball and stick representations: thin dashed cylinders, one primitive per entry
static bool draw_bond_partners(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, const bond_params_t* bond_params, const float viewport_size[2]) {
    ASSERT(mol);
    ASSERT(bond_params);
    if (mol->bond_partner_count == 0 || mol->bond_order_count != mol->bond_count) {
        return false;
    }
    bond_params_t params = *bond_params;
    params.dash_count  = BOND_PARTNER_DASHES;
    params.dash_fill   = BOND_PARTNER_FILL;
    params.bond_orders = 0;
    const bond_set_t partners = {mol->buffer[GL_BUFFER_BOND_ATOM_INDICES], {0}, {0}, {0}, mol->buffer[GL_BUFFER_BOND_PARTNER], mol->bond_partner_count};
    return draw_bonds(program, mol, atom_color, partners, &params, viewport_size);
}

// The gaps of the chains (md_gl_mol_set_chain_gaps) as dashed cylinders (bond.vert and bond.frag, CHAIN_GAPS): the
// entries of the cartoon and ribbons (trace), or those of licorice and ball and stick, one primitive per entry
static bool draw_chain_gaps(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, bool trace, const bond_params_t* params, const float viewport_size[2]) {
    ASSERT(mol);
    ASSERT(params);
    ASSERT(atom_color.id);

    const uint32_t offset = trace ? 0 : mol->gap_trace_count;
    const uint32_t count  = trace ? mol->gap_trace_count : mol->gap_link_count;
    if (count == 0 || !mol->buffer[GL_BUFFER_CHAIN_GAP].id || !mol->buffer[GL_BUFFER_CHAIN_GAP_COMP].id || !ctx.bond_quad_index.id) {
        return false;
    }

    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(bond_params_t), params);

    glBindVertexArray(ctx.vao);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_POSITION].id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_2].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, mol->buffer[GL_BUFFER_ATOM_FLAGS].id);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_3].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color.id);

    glActiveTexture(GL_TEXTURE6);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_6].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, mol->buffer[GL_BUFFER_CHAIN_GAP_COMP].id);

    glActiveTexture(GL_TEXTURE8);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_8].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, mol->buffer[GL_BUFFER_CHAIN_GAP].id);

    glUseProgram(program.id);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_pos"),   0);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_vel"),   1);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_flags"), 2);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_color"), 3);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_gap_comp"),   6);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_gap"),        8);
    glUniform1ui(glGetUniformLocation(program.id, "u_quads_per_instance"), BOND_QUADS_PER_INSTANCE);
    glUniform2f(glGetUniformLocation(program.id, "u_viewport_size"), viewport_size[0], viewport_size[1]);
    const GLint offset_loc = glGetUniformLocation(program.id, "u_bond_offset");

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ctx.bond_quad_index.id);

    const uint32_t full = count / BOND_QUADS_PER_INSTANCE;
    const uint32_t rest = count % BOND_QUADS_PER_INSTANCE;
    if (full > 0) {
        glUniform1ui(offset_loc, offset);
        glDrawElementsInstanced(GL_TRIANGLES, BOND_QUADS_PER_INSTANCE * 6, GL_UNSIGNED_SHORT, 0, (GLsizei)full);
    }
    if (rest > 0) {
        glUniform1ui(offset_loc, offset + full * BOND_QUADS_PER_INSTANCE);
        glDrawElements(GL_TRIANGLES, (GLsizei)(rest * 6), GL_UNSIGNED_SHORT, 0);
    }

    glUseProgram(0);

    glActiveTexture(GL_TEXTURE0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return true;
}

// The circles of the aromatic rings (md_gl_mol_set_bond_orders), drawn with the bonds of the licorice and ball and stick
// representations, by vertex pulling (ring.vert): one quad per ring, from the index buffer of the bonds. They take the
// radius, the colors and the reach of the bonds (bond_params), and write ring_index_base + the index of the ring to the
// picking target.
static bool draw_aromatic_rings(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, const bond_params_t* bond_params, uint32_t ring_index_base, const float viewport_size[2]) {
    ASSERT(mol);
    ASSERT(bond_params);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(atom_color.id);

    if (mol->aromatic_ring_count == 0 || mol->bond_order_count != mol->bond_count || !mol->buffer[GL_BUFFER_AROMATIC_RING].id || !ctx.bond_quad_index.id) {
        return false;
    }

    const ring_params_t params = {
        .radius          = bond_params->radius,
        .max_d2          = bond_params->max_d2,
        .mode            = bond_params->mode,
        .sharpness       = bond_params->sharpness,
        .uniform_color   = bond_params->uniform_color,
        .spread          = bond_params->spread,
        .ring_index_base = ring_index_base,
    };
    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(ring_params_t), &params);

    glBindVertexArray(ctx.vao);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_POSITION].id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_2].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, mol->buffer[GL_BUFFER_ATOM_FLAGS].id);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_3].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color.id);

    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_4].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, mol->buffer[GL_BUFFER_AROMATIC_RING].id);

    glUseProgram(program.id);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_pos"),   0);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_vel"),   1);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_flags"), 2);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_color"), 3);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_ring"),       4);
    glUniform1ui(glGetUniformLocation(program.id, "u_quads_per_instance"), BOND_QUADS_PER_INSTANCE);
    glUniform2f(glGetUniformLocation(program.id, "u_viewport_size"), viewport_size[0], viewport_size[1]);
    const GLint ring_offset_loc = glGetUniformLocation(program.id, "u_ring_offset");

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ctx.bond_quad_index.id);

    const uint32_t count = mol->aromatic_ring_count;
    const uint32_t full  = count / BOND_QUADS_PER_INSTANCE;
    const uint32_t rest  = count % BOND_QUADS_PER_INSTANCE;
    if (full > 0) {
        glUniform1ui(ring_offset_loc, 0);
        glDrawElementsInstanced(GL_TRIANGLES, BOND_QUADS_PER_INSTANCE * 6, GL_UNSIGNED_SHORT, 0, (GLsizei)full);
    }
    if (rest > 0) {
        glUniform1ui(ring_offset_loc, full * BOND_QUADS_PER_INSTANCE);
        glDrawElements(GL_TRIANGLES, (GLsizei)(rest * 6), GL_UNSIGNED_SHORT, 0);
    }

    glUseProgram(0);

    glActiveTexture(GL_TEXTURE0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
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

// Slabs by vertex pulling (slab.vert): one quad per slab, as the aromatic rings
static bool draw_slabs(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, float radius, float plate, float max_length, uint32_t slab_index_base, const float viewport_size[2]) {
    ASSERT(mol);
    ASSERT(atom_color.id);
    if (mol->slab_count == 0 || !mol->buffer[GL_BUFFER_SLAB].id || !ctx.bond_quad_index.id || !(radius > 0.0f)) {
        return false;
    }

    // After gl_ubo_base_t, as in slab.vert
    const struct {
        float    radius;
        float    plate;
        float    max_d2;
        uint32_t slab_index_base;
    } params = { radius, plate, max_length * max_length, slab_index_base };
    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(params), &params);

    glBindVertexArray(ctx.vao);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_POSITION].id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_2].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, mol->buffer[GL_BUFFER_ATOM_FLAGS].id);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_3].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, atom_color.id);

    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_4].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, mol->buffer[GL_BUFFER_SLAB].id);

    glUseProgram(program.id);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_pos"),   0);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_vel"),   1);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_flags"), 2);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_color"), 3);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_slab"),       4);
    glUniform1ui(glGetUniformLocation(program.id, "u_quads_per_instance"), BOND_QUADS_PER_INSTANCE);
    glUniform2f(glGetUniformLocation(program.id, "u_viewport_size"), viewport_size[0], viewport_size[1]);
    const GLint slab_offset_loc = glGetUniformLocation(program.id, "u_slab_offset");

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ctx.bond_quad_index.id);

    const uint32_t count = mol->slab_count;
    const uint32_t full  = count / BOND_QUADS_PER_INSTANCE;
    const uint32_t rest  = count % BOND_QUADS_PER_INSTANCE;
    if (full > 0) {
        glUniform1ui(slab_offset_loc, 0);
        glDrawElementsInstanced(GL_TRIANGLES, BOND_QUADS_PER_INSTANCE * 6, GL_UNSIGNED_SHORT, 0, (GLsizei)full);
    }
    if (rest > 0) {
        glUniform1ui(slab_offset_loc, full * BOND_QUADS_PER_INSTANCE);
        glDrawElements(GL_TRIANGLES, (GLsizei)(rest * 6), GL_UNSIGNED_SHORT, 0);
    }

    glUseProgram(0);

    glActiveTexture(GL_TEXTURE0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return true;
}

// One instance per control point (residue), see backbone.vert. Consecutive chains with the same level of detail are drawn together.
// Ray cast (backbone_rc.vert): the most segments per control point at each level of detail, even divisors of
// MD_GL_BACKBONE_SEGMENT_COUNT (those of the mesh), and the largest angle a segment may turn by, which decides how many
// of them each control point uses (a helix all, a strand few). The segments are straight: where they meet the surface
// creases by that angle.
static const uint32_t backbone_rc_segments[BACKBONE_LOD_COUNT] = { 12, 6, 2 };
static const float    backbone_rc_turn[BACKBONE_LOD_COUNT]     = { (float)DEG_TO_RAD(6.0), (float)DEG_TO_RAD(12.0), 0.0f };

// With raycast, the stadium profile is ray cast (backbone_rc.vert) with a quad per segment, the meshes of profile only
// decide the level of detail
static bool draw_backbone(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, uint32_t profile, const float scale[4], float max_extent, const mat4_t* world_to_clip, float viewport_half_height, bool raycast, const float viewport_size[2]) {
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
    if (raycast) {
        glUniform2f(glGetUniformLocation(program.id, "u_viewport_size"), viewport_size[0], viewport_size[1]);
    }

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

            const uint32_t segments = raycast ? backbone_rc_segments[lod] : mesh->segments;
            if (!bound) {
                struct {
                    float    scale[4];
                    uint32_t res[4];
                    uint32_t rings[4];
                } params = {
                    {scale[0], scale[1], scale[2], raycast ? backbone_rc_turn[lod] : scale[3]},
                    {segments, mesh->profile_count, mesh->outline_count, profile == BACKBONE_PROFILE_BOX ? mesh->face_subdivisions : mesh->dome_rings},
                    {MD_GL_BACKBONE_SEGMENT_COUNT, MD_GL_BACKBONE_SEGMENT_COUNT / segments, 0, 0},
                };
                gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(params), &params);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, raycast ? ctx.bond_quad_index.id : mesh->index_buffer.id);
                bound = true;
            }

            glUniform1ui(instance_offset_loc, beg);
            // A quad (two triangles of 6 indices) per segment, or the mesh
            const GLsizei index_count = raycast ? (GLsizei)(segments * 6) : (GLsizei)mesh->index_count;
            glDrawElementsInstanced(GL_TRIANGLES, index_count, GL_UNSIGNED_SHORT, 0, (GLsizei)(end - beg));
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
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, ca));

    glEnableVertexAttribArray(2);
    glVertexAttribIPointer(2, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, c));

    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, o));

    glEnableVertexAttribArray(4);
    glVertexAttribIPointer(4, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, flags));

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
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, radius));

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
