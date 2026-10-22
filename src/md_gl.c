#include <md_gl.h>
#include <md_util.h>
#include <md_chem.h>

#include <core/md_common.h>
#include <core/md_compiler.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_os.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_handle.h>
#include <core/md_gl_util.h>
#include <core/md_intrinsics.h>
#include <core/md_bitset.h>

#include <md_system.h>

#include <stdbool.h>
#include <string.h>
#include <stdio.h>      // printf etc
#include <stdlib.h>     // qsort
#include <stddef.h>     // offsetof
#include <float.h>      // FLT_MAX

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

// Spheres, bonds, rings, slabs and gaps are drawn as one quad each (spacefill.vert, bond.vert, ring.vert, slab.vert), in
// blocks of this many quads: an instance of the draw is one block of one copy of the molecule, from a static index
// buffer of one block (QUADS in instance_shader_src). A partial block at the end culls the quads past the count in the
// vertex shader. A draw with an atom mask draws a list of blocks, those which hold marked atoms (MARKS); of spheres the
// list of the marked atoms themselves.
#define QUAD_BLOCK 128

// INSTANCES (md_gl_copies_t). Each draw of md_gl_draw uploads the copies of its ops which are in view to a buffer read
// as a texture buffer, INSTANCE_TEXELS RGBA32F texels each: the rows of the 3x4 transform, then those of its transform
// in the previous frame. A draw names the first of its copies (u_instance_base) and the shader reads that plus the
// instance of the vertex (instance_shader_src).
#define INSTANCE_TEXELS 6
// Texture units of the instances and of the lists of the draw (MARKS), set in every program once, after those the
// draws bind and before MD_GL_USER_TEXTURE_UNIT
#define INSTANCE_TEXTURE_UNIT 9
#define LIST_TEXTURE_UNIT 10
// DASHED LINES. Every dashed line md_gl draws (the weak bonds, the partners of delocalized bonds, the gaps of chains)
// is given by the length of its dashes and of the gaps between them, never by a number of dashes: as many dashes fit
// as come closest to the lengths, at least two, the first starting at one end and the last ending at the other, the
// dashes and gaps stretched or shrunk alike to fit (bond.vert). So no dash is cut off at an end, and each kind of line
// keeps a look of its own at any length. The lines md_gl styles itself have their lengths in proportion to their
// radius, so they keep that look at any thickness too.
// The one exception are the dashed circles of fused aromatic rings (ring.vert): a closed loop tied to the atoms of its
// ring, so a number of dashes per bond, the gaps centred on the atoms, which follows from the size of the ring and the
// radius of the sticks alone and so does not change from frame to frame as the atoms move.

// The dashed partners of delocalized bonds (and of aromatic bonds outside of a ring drawn as a circle): the lengths of
// their dashes and gaps relative to the radius of the stick (the partner is a third of it thick): short and tight,
// about three dashes along the partner of a bond in a ring
#define BOND_PARTNER_DASH_LENGTH 0.75f
#define BOND_PARTNER_GAP_LENGTH  0.75f

// The gaps of the chains (md_gl_mol_set_chain_gaps): the radius of the dashed line relative to the coil of the cartoon
// (and to the thickness of the ribbons), and to the sticks of licorice and ball and stick, and the lengths of its
// dashes and gaps relative to that radius: long dashes, which tell a gap from a bond
#define CHAIN_GAP_RADIUS_TRACE    0.7f
#define CHAIN_GAP_RADIUS_LINK     0.5f
#define CHAIN_GAP_DASH_LENGTH     4.0f
#define CHAIN_GAP_GAP_LENGTH      4.0f

#define UBO_SIZE (1 << 10)

// The uniform buffers of the draws (ubo_write): a ring of buffer objects, each write of the uniforms to the next one.
// One buffer written before every draw (glBufferSubData) is a buffer that the draws queued before it still read, and a
// GL that does not copy the data waits for them to finish (Apple documents this of its GL): the CPU then waits for the
// GPU at every md_gl_draw and every op in it, and work drawn in several passes over the same geometry (transparency) runs
// in lockstep, the CPU time of each pass the GPU time of the one before (as an M4 Pro profiled). A buffer of the ring is
// written again UBO_RING writes later, when what read it is long done (a frame writes some tens), and respecified whole
// (glBufferData), which lets a driver give it new storage should it still be read. The instances and the lists of a
// draw are respecified whole already (glBufferData), not written into.
#define UBO_RING 256
#define SHADER_BUF_SIZE KILOBYTES(14)


static const str_t default_shader_output = STR_INIT(
"layout(location = 0) out vec4 out_color;\n"
"void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint atom_index, uint atom_flags) {\n"
"    out_color = color;\n"
"}\n");

// Injected at '#pragma INSTANCES' of the vertex shaders of the representations, after their uniform block.
// INSTANCES: instance_load reads the transforms of instance i of the draw, and the positions and directions of the
// molecule are taken into world space through them. The velocity of a point is its displacement over the last frame in
// world space: the motion of the atom (v, in the molecule) carried by the transform of the previous frame, and the
// motion of the instance itself, which is 0 when the transform did not change.
// QUADS: a draw of quads is u_block_count blocks of QUAD_BLOCK quads per instance, instance-major in gl_InstanceID: all
// blocks of the primitives [u_prim_offset, u_prim_offset + u_prim_count) in order, or (u_block_list >= 0) the blocks in
// the list there in u_buf_list, or (u_prim_list >= 0) the u_prim_count primitives in the list there, in blocks of
// QUAD_BLOCK entries of it (a subset, SUBSETS), each u_prim_offset + its entry. quad_primitive returns the primitive of
// the vertex, QUAD_NONE past the count.
static const str_t instance_shader_src = STR_INIT(
"#define QUAD_BLOCK " STRINGIFY_VAL(QUAD_BLOCK) "u\n"
"#define QUAD_NONE 0xFFFFFFFFu\n"
"uniform samplerBuffer  u_buf_instance;  // RGBA32F, " STRINGIFY_VAL(INSTANCE_TEXELS) " per instance\n"
"uniform usamplerBuffer u_buf_list;      // R32UI, the lists of the draw\n"
"uniform int  u_instance_base;\n"
"uniform uint u_prim_offset;\n"
"uniform uint u_prim_count;\n"
"uniform uint u_block_count;\n"
"uniform int  u_block_list;\n"
"uniform int  u_prim_list;\n"
"vec4 instance_row[6];\n"
"void instance_load(uint i) {\n"
"    int base = (u_instance_base + int(i)) * 6;\n"
"    for (int k = 0; k < 6; ++k) instance_row[k] = texelFetch(u_buf_instance, base + k);\n"
"}\n"
"vec3 instance_point(vec3 p) {\n"
"    vec4 q = vec4(p, 1.0);\n"
"    return vec3(dot(instance_row[0], q), dot(instance_row[1], q), dot(instance_row[2], q));\n"
"}\n"
"vec3 instance_vector(vec3 v) {\n"
"    return vec3(dot(instance_row[0].xyz, v), dot(instance_row[1].xyz, v), dot(instance_row[2].xyz, v));\n"
"}\n"
"vec3 instance_velocity(vec3 p, vec3 v) {\n"
"    vec4 q = vec4(p - v, 1.0);\n"
"    return instance_point(p) - vec3(dot(instance_row[3], q), dot(instance_row[4], q), dot(instance_row[5], q));\n"
"}\n"
"uint quad_primitive() {\n"
"    uint i = uint(gl_InstanceID) / u_block_count;\n"
"    uint b = uint(gl_InstanceID) - i * u_block_count;\n"
"    instance_load(i);\n"
"    uint q = b * QUAD_BLOCK + (uint(gl_VertexID) >> 2u);\n"
"    if (u_prim_list >= 0) {\n"
"        return q < u_prim_count ? u_prim_offset + texelFetch(u_buf_list, u_prim_list + int(q)).x : QUAD_NONE;\n"
"    }\n"
"    if (u_block_list >= 0) q = texelFetch(u_buf_list, u_block_list + int(b)).x * QUAD_BLOCK + (uint(gl_VertexID) >> 2u);\n"
"    return q < u_prim_count ? u_prim_offset + q : QUAD_NONE;\n"
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
    GL_TEXTURE_BUFFER_INSTANCE,     // ctx.instance_buf, on INSTANCE_TEXTURE_UNIT
    GL_TEXTURE_BUFFER_LIST,         // ctx.list_buf, on LIST_TEXTURE_UNIT
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
    GL_BUFFER_BOND_ORDER,                      // u32 per bond: kind (plain, double, triple, elastic) | index into GL_BUFFER_BOND_ORDER_REF << 3 (md_gl_mol_set_bond_orders)
    GL_BUFFER_BOND_ORDER_REF,                  // u32[4] per double or triple bond: two atoms bonded to its first atom, two to its second
    GL_BUFFER_BOND_PARTNER,                    // u32[4] per dashed partner of a delocalized bond, or of an aromatic bond in no circle: bond, atoms bonded to its first and its second atom, 0
    GL_BUFFER_AROMATIC_RING,                   // u32[8] per aromatic ring drawn as a circle: ring index | dashed << 28 | atom count << 29, its atoms in order around it (RING_ATOM_NONE past the count)
    GL_BUFFER_BACKBONE_DATA,                   // u32: residue index, u32: residue atom offset, u8: CA index, C index and O Index, u8: flags
    GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE,    // u8[4]  (0: Unknown, 1: Coil, 2: Helix, 3: Sheet)
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA,     // Extracted control points
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT,      // Oriented control points (support vector + relation to the next control point)
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT_PREV, // The previous computation's oriented control points, swapped with the one above each computation
    GL_BUFFER_BACKBONE_NEIGHBOR,               // u32[4] per control point: prev, next, prev2, next2 (clamped to the chain)
    GL_BUFFER_BACKBONE_RING_DATA,              // Evaluated cross section frames of the spline, MD_GL_BACKBONE_SEGMENT_COUNT per control point (see compute_spline_subdivide.vert)
    GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX,    // u32, LINE_STRIP_ADJACENCY Indices for legacy and debugging paths
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
#define MAX_SUBSETS 4096
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
    uint32_t c, o;                  // The fallback orientation O - C (a protein: its carbonyl), or with BACKBONE_FLAG_FIXED_FACE the normal of the flat faces, with BACKBONE_FLAG_FIXED_SIDE the long axis of the profile
    uint32_t flags;                 // BACKBONE_FLAG_*
} gl_backbone_data_t;

// Flags of a control point (gl_backbone_data_t, and the top byte of the secondary structure of gl_control_point_t)
enum {
    BACKBONE_FLAG_CHAIN_BEG  = 1,
    BACKBONE_FLAG_CHAIN_END  = 2,
    BACKBONE_FLAG_FLIP_NEXT  = 16,  // Set by the orient pass
    BACKBONE_FLAG_FIXED_FACE = 32,  // The flat faces of the profile face O - C, rather than following the curve (compute_spline_orient.vert)
    BACKBONE_FLAG_THROUGH    = 64,  // The spline passes through the control point rather than near it (compute_spline_orient.vert)
    BACKBONE_FLAG_FIXED_SIDE = 128, // The long axis of the profile points along O - C, rather than following the curve (compute_spline_orient.vert)
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
} shaders_t;

// CPU side bounds of the backbone chains, for culling and level of detail
typedef struct {
    uint32_t  chain_count;
    uint32_t* chain_offset;     // [chain_count + 1] Control point ranges of the chains
    uint32_t* ca_atom_idx;      // [control point count]
    vec3_t*   ca_xyz;           // [control point count] Last positions of the control points (CA)
    vec4_t*   chain_sphere;     // [chain_count] Bounding sphere of the control points of the chain (xyz center, w radius)
    vec4_t    all_sphere;       // A sphere which holds those of all chains: what is within it holds the levels of detail of all of them
    bool      valid;            // The positions have been set
} backbone_bounds_t;

// The primitives of a kind (bonds, rings, control points, ...). In blocks of QUAD_BLOCK, with the range of the atoms
// each block reads the flags of: a block holds an atom with the mask of a draw only if one in its range has it (MARKS).
// And the atoms each primitive is shown with, PRIM_ATOM_NONE past those of one: a primitive is in a subset when all of
// them are (SUBSETS). Those of the control points and of the weak bonds are kept as they are set, the others are read
// back from the buffers of the molecule when a subset first needs them (mol_kind_atoms), so that a molecule drawn
// without subsets keeps no copy of its bonds. Which primitives each atom is in (adj) is found then too.
typedef struct {
    uint32_t  count;            // Blocks
    uint32_t* range;            // [2 * count] (first, last) atom, first > last for none
    uint32_t  prim_count;
    uint32_t  stride;           // Atoms per primitive
    uint32_t* atoms;            // [prim_count * stride] NULL until read back, where they are not kept
    uint64_t  version;          // Changes with the atoms of the primitives, never 0 once initialized
    uint32_t  adj_n;            // Atoms adj was found for, 0 for not found
    uint64_t  adj_version;      // Of the atoms it was found from
    uint32_t* adj_offset;       // [adj_n + 1] The primitives of atom a: adj_prim[adj_offset[a] .. adj_offset[a + 1])
    uint32_t* adj_prim;
} prim_blocks_t;

#define PRIM_ATOM_NONE 0xFFFFFFFFu

// MARKS. The atoms with all bits of the mask of a draw (md_gl_draw_args_t atom_mask), found from the flags as last set on
// the CPU, and kept until those or the mask change: as bits with the count of the bits before each word, for whether a
// range of atoms holds one, as the atoms in order (the spheres) and as the control points of the backbone whose atom
// is one, in order.
typedef struct {
    bool      valid;
    uint32_t  mask;
    uint64_t  flags_version;
    uint32_t  n;                // Atoms of the molecule it was found for
    uint32_t  word_count;
    uint64_t* bits;             // [word_count]
    uint32_t* prefix;           // [word_count + 1]
    uint32_t  atom_count;
    uint32_t* atoms;            // [n]
    uint32_t  cp_count;
    uint32_t* cps;              // [cp_capacity]
    uint32_t  cp_capacity;
} marked_t;

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

    // The box of the positions of the atoms, as last set or derived on the CPU, for the culling of the copies of a
    // draw (md_gl_copies_t). Not valid: never culled.
    vec3_t aabb_min;
    vec3_t aabb_max;
    bool   aabb_valid;
    float  max_radius;              // The largest radius of an atom set (md_gl_mol_set_atom_radius)

    // What the spline of the backbone is computed from changed (compute_spline): the positions, velocities, secondary
    // structure or the backbone itself. Computed in the first draw after, once.
    bool spline_dirty;

    // The flags as last set (md_gl_mol_set_atom_flags), on the CPU as well, for the draws with an atom mask (MARKS)
    uint8_t* atom_flags;            // [atom_count], NULL until set: all 0
    uint64_t flags_version;
    marked_t marked;

    // The blocks of the primitives of each kind (MARKS, SUBSETS)
    prim_blocks_t bond_blocks;
    prim_blocks_t partner_blocks;
    prim_blocks_t ring_blocks;
    prim_blocks_t gap_trace_blocks;
    prim_blocks_t gap_link_blocks;
    prim_blocks_t slab_blocks;
    prim_blocks_t cp_blocks;        // The control points of the backbone, each with its atom (backbone_bounds_t ca_atom_idx)

    // A bit per atom, all clear between uses: the atoms of a subset while what it holds is found (SUBSETS)
    uint64_t* scratch_bits;

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
    float weak_bond_max_weight;     // The largest weight, which scales the radius
    prim_blocks_t weak_bond_blocks; // MARKS, SUBSETS
} representation_t;

// SUBSETS. The kinds of primitives a subset holds lists of, found from its atoms and kept until they or the primitives
// change (the weak bonds of a representation change every frame and are found in each draw instead)
enum {
    SUBSET_CPS,
    SUBSET_BONDS,
    SUBSET_PARTNERS,
    SUBSET_RINGS,
    SUBSET_GAP_TRACE,
    SUBSET_GAP_LINK,
    SUBSET_SLABS,
    SUBSET_KIND_COUNT
};

#define SUBSET_LIST_STALE UINT64_MAX

// The control points of a subset as drawn (backbone.vert): with whether the one before and the one after are out of it
#define CP_PREV_OUT 0x80000000u
#define CP_NEXT_OUT 0x40000000u

typedef struct {
    uint64_t  version;          // Of the primitives it was found from (prim_blocks_t version), SUBSET_LIST_STALE for none
    uint32_t  count;
    uint32_t  capacity;
    uint32_t* prims;            // [count] The primitives within the subset, in order
} subset_list_t;

// The lists a draw of a subset reads (LIST_* bits of mol_lists_t): its atoms, control points and the primitives of each
// kind, one after another in a buffer of its own on the GPU, uploaded when one of them changes: a draw of a subset
// which did not change uploads nothing
#define SUBSET_GPU_LISTS 8

typedef struct {
    uint32_t id;
    uint32_t mol_id;
    uint32_t atom_count;
    uint32_t atom_capacity;
    uint32_t* atoms;            // [atom_count] In order, each once, within the molecule
    subset_list_t list[SUBSET_KIND_COUNT];

    gl_buffer_t gpu;
    uint32_t    gpu_lists;      // LIST_* bits of the lists in it
    uint32_t    gpu_offset[SUBSET_GPU_LISTS];
    uint32_t    gpu_count[SUBSET_GPU_LISTS];
    bool        gpu_stale;      // Its atoms changed since

    // The box of the positions of its atoms, as last set (md_gl_mol_set_atom_position) or derived, for the culling of
    // the copies of an op of it. Not valid: found from the positions on the GPU when a draw of copies needs it.
    vec3_t aabb_min;
    vec3_t aabb_max;
    bool   aabb_valid;
} subset_t;

typedef struct {
    // Shared resources
    GLuint vao;
    GLuint fbo;
    GLuint ubo[UBO_RING];           // The ring of uniform buffers (ubo_write)
    uint32_t ubo_next;              // The one written next
    uint8_t ubo_data[UBO_SIZE];     // What the draws read: the shared part (gl_ubo_base_t), a draw's parameters after it
    gl_buffer_t instance_buf;       // INSTANCES of the draw, RGBA32F
    gl_buffer_t list_buf;           // Lists of the draw (MARKS), u32: atoms (indices of the spheres), control points, blocks
    gl_texture_t texture[GL_TEXTURE_COUNT];
    gl_program_t program[GL_PROGRAM_COUNT];
    backbone_mesh_t backbone_mesh[BACKBONE_PROFILE_COUNT][BACKBONE_LOD_COUNT];
    gl_buffer_t quad_index;         // u16, QUAD_BLOCK quads of two triangles
    uint32_t version;

    // Handle pools
    md_handle_pool_t shader_pool;
    md_handle_pool_t molecule_pool;
    md_handle_pool_t palette_pool;
    md_handle_pool_t representation_pool;
    md_handle_pool_t subset_pool;

    // Internal structures which are referenced by handles
    shaders_t shaders[MAX_SHADERS];
    molecule_t molecules[MAX_MOLECULES];
    palette_t palettes[MAX_PALETTES];
    representation_t representations[MAX_REPRESENTATIONS];
    subset_t subsets[MAX_SUBSETS];

    // The last version given to the primitives of a kind (prim_blocks_t version)
    uint64_t prim_version;

    // The buffer the list texture of a draw is on (LIST_TEXTURE_UNIT): the lists of the draw, or those of a subset
    GLuint list_bound;

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

static inline subset_t* subset_lookup(uint32_t id) {
    if (id == MD_HANDLE_INVALID_ID) {
        return NULL;
    }
    int index = md_handle_index(id);
    if (index >= ctx.subset_pool.size || ctx.subsets[index].id != id) {
        MD_LOG_ERROR("Invalid id");
        return NULL;
    }
    return ctx.subsets + index;
}

// SUBSETS: the boxes of the subsets of a molecule, from positions of its atoms [offset, offset + count) as they are set:
// of all of them, or grown by those of a part. Not known: found again when needed.
static void subsets_aabb_update(const molecule_t* mol, uint32_t offset, uint32_t count, const vec3_t* xyz);
static void subsets_aabb_invalidate(const molecule_t* mol);

static inline bool validate_context(void) {
    if (ctx.version == 0) {
        MD_LOG_ERROR("MD GL module has not been initialized");
        return false;
    }
    return true;
}

// The uniforms of the draws from offset on (the shared part at 0, a draw's parameters after it), with the rest as last
// written, into the next buffer of the ring (UBO_RING), which is then the one bound
static void ubo_write(uint32_t offset, uint32_t size, const void* data) {
    ASSERT(offset + size <= UBO_SIZE);
    MEMCPY(ctx.ubo_data + offset, data, size);
    const GLuint buf = ctx.ubo[ctx.ubo_next];
    ctx.ubo_next = (ctx.ubo_next + 1) % UBO_RING;
    glBindBuffer(GL_UNIFORM_BUFFER, buf);
    glBufferData(GL_UNIFORM_BUFFER, UBO_SIZE, ctx.ubo_data, GL_STREAM_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, buf);
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

// ### BLOCKS AND MARKS ###

static void blocks_free_adjacency(prim_blocks_t* blocks) {
    md_allocator_i* alloc = md_get_heap_allocator();
    if (blocks->adj_offset) {
        md_free(alloc, blocks->adj_prim,   MAX(blocks->adj_offset[blocks->adj_n], 1) * sizeof(uint32_t));
        md_free(alloc, blocks->adj_offset, (blocks->adj_n + 1) * sizeof(uint32_t));
    }
    blocks->adj_offset  = NULL;
    blocks->adj_prim    = NULL;
    blocks->adj_n       = 0;
    blocks->adj_version = 0;
}

static void blocks_free(prim_blocks_t* blocks) {
    md_allocator_i* alloc = md_get_heap_allocator();
    if (blocks->range) {
        md_free(alloc, blocks->range, 2 * blocks->count * sizeof(uint32_t));
    }
    if (blocks->atoms) {
        md_free(alloc, blocks->atoms, (size_t)blocks->prim_count * blocks->stride * sizeof(uint32_t));
    }
    blocks_free_adjacency(blocks);
    blocks->range = NULL;
    blocks->count = 0;
    blocks->atoms = NULL;
    blocks->prim_count = 0;
    blocks->stride = 0;
    // Other than any it had: a subset found from them finds again
    blocks->version = ++ctx.prim_version;
}

// The blocks of prim_count primitives of at most stride atoms each, none of which holds an atom yet: add each atom of
// each primitive after. keep_atoms: their atoms are kept as they are added, otherwise read back when needed.
static void blocks_init(prim_blocks_t* blocks, uint32_t prim_count, uint32_t stride, bool keep_atoms) {
    blocks_free(blocks);
    blocks->count = DIV_UP(prim_count, QUAD_BLOCK);
    if (blocks->count == 0) return;
    md_allocator_i* alloc = md_get_heap_allocator();
    blocks->range = md_alloc(alloc, 2 * blocks->count * sizeof(uint32_t));
    for (uint32_t b = 0; b < blocks->count; ++b) {
        blocks->range[2 * b + 0] = UINT32_MAX;
        blocks->range[2 * b + 1] = 0;
    }
    blocks->prim_count = prim_count;
    blocks->stride = stride;
    if (keep_atoms) {
        blocks->atoms = md_alloc(alloc, (size_t)prim_count * stride * sizeof(uint32_t));
        MEMSET(blocks->atoms, 0xFF, (size_t)prim_count * stride * sizeof(uint32_t));    // PRIM_ATOM_NONE
    }
}

static inline void blocks_add(prim_blocks_t* blocks, uint32_t prim, uint32_t atom) {
    const uint32_t b = prim / QUAD_BLOCK;
    if (b >= blocks->count || prim >= blocks->prim_count) return;
    uint32_t* r = blocks->range + 2 * b;
    r[0] = MIN(r[0], atom);
    r[1] = MAX(r[1], atom);
    if (blocks->atoms) {
        uint32_t* a = blocks->atoms + (size_t)prim * blocks->stride;
        for (uint32_t k = 0; k < blocks->stride; ++k) {
            if (a[k] == PRIM_ATOM_NONE) {
                a[k] = atom;
                break;
            }
        }
    }
    blocks->version = ++ctx.prim_version;
}

// The primitive holds no atoms: set its atoms again (blocks_add). The range of its block keeps those it had.
static inline void blocks_clear_prim(prim_blocks_t* blocks, uint32_t prim) {
    if (prim >= blocks->prim_count) return;
    if (blocks->atoms) {
        uint32_t* a = blocks->atoms + (size_t)prim * blocks->stride;
        for (uint32_t k = 0; k < blocks->stride; ++k) a[k] = PRIM_ATOM_NONE;
    }
    blocks->version = ++ctx.prim_version;
}

// Which primitives each of the n atoms is in, found again if the primitives changed since (SUBSETS)
static void blocks_adjacency(prim_blocks_t* blocks, uint32_t n) {
    if (blocks->adj_offset && blocks->adj_n == n && blocks->adj_version == blocks->version) return;
    blocks_free_adjacency(blocks);
    if (n == 0) return;
    md_allocator_i* alloc = md_get_heap_allocator();
    uint32_t* offset = md_alloc(alloc, (n + 1) * sizeof(uint32_t));
    MEMSET(offset, 0, (n + 1) * sizeof(uint32_t));
    const size_t slots = (size_t)blocks->prim_count * blocks->stride;
    for (size_t i = 0; i < slots; ++i) {
        const uint32_t a = blocks->atoms[i];
        if (a < n) offset[a + 1] += 1;
    }
    for (uint32_t a = 0; a < n; ++a) offset[a + 1] += offset[a];
    const uint32_t total = offset[n];
    uint32_t* prim = md_alloc(alloc, MAX(total, 1) * sizeof(uint32_t));
    // Filled from the start of each atom's range, which moves along: shifted back after
    for (size_t i = 0; i < slots; ++i) {
        const uint32_t a = blocks->atoms[i];
        if (a < n) prim[offset[a]++] = (uint32_t)(i / blocks->stride);
    }
    for (uint32_t a = n; a > 0; --a) offset[a] = offset[a - 1];
    offset[0] = 0;
    blocks->adj_offset  = offset;
    blocks->adj_prim    = prim;
    blocks->adj_n       = n;
    blocks->adj_version = blocks->version;
}

static void marked_free(marked_t* m) {
    md_allocator_i* alloc = md_get_heap_allocator();
    if (m->bits)   md_free(alloc, m->bits,   m->word_count * sizeof(uint64_t));
    if (m->prefix) md_free(alloc, m->prefix, (m->word_count + 1) * sizeof(uint32_t));
    if (m->atoms)  md_free(alloc, m->atoms,  m->n * sizeof(uint32_t));
    if (m->cps)    md_free(alloc, m->cps,    m->cp_capacity * sizeof(uint32_t));
    MEMSET(m, 0, sizeof(marked_t));
}

// Atoms before atom i which are marked
static inline uint32_t marked_before(const marked_t* m, uint32_t i) {
    const uint32_t w = i >> 6;
    const uint32_t b = i & 63;
    uint32_t count = m->prefix[w];
    if (b) count += (uint32_t)popcnt64(m->bits[w] & ((1ull << b) - 1));
    return count;
}

// Whether an atom in [first, last] is marked
static inline bool marked_any(const marked_t* m, uint32_t first, uint32_t last) {
    if (first > last || first >= m->n) return false;
    const uint32_t end = (last >= m->n - 1) ? m->n : last + 1;
    return marked_before(m, end) > marked_before(m, first);
}

// The atoms of the molecule with all bits of the mask (MARKS)
static const marked_t* mol_marked(molecule_t* mol, uint32_t mask) {
    marked_t* m = &mol->marked;
    if (m->valid && m->mask == mask && m->flags_version == mol->flags_version && m->n == mol->atom_count) {
        return m;
    }
    const backbone_bounds_t* bounds = &mol->backbone_bounds;
    const uint32_t n = mol->atom_count;
    const uint32_t words = DIV_UP(n, 64);
    const uint32_t cp_count = bounds->chain_count ? bounds->chain_offset[bounds->chain_count] : 0;
    if (m->n != n || m->cp_capacity != MAX(cp_count, 1) || !m->bits) {
        marked_free(m);
        md_allocator_i* alloc = md_get_heap_allocator();
        // A molecule has atoms (md_gl_mol_create, md_gl_mol_create_particles): n and words are at least 1
        m->n = n;
        m->word_count  = words;
        m->cp_capacity = MAX(cp_count, 1);
        m->bits   = md_alloc(alloc, words * sizeof(uint64_t));
        m->prefix = md_alloc(alloc, (words + 1) * sizeof(uint32_t));
        m->atoms  = md_alloc(alloc, n * sizeof(uint32_t));
        m->cps    = md_alloc(alloc, m->cp_capacity * sizeof(uint32_t));
    }

    const uint8_t* flags = mol->atom_flags;
    uint32_t count = 0;
    m->prefix[0] = 0;
    for (uint32_t w = 0; w < words; ++w) {
        const uint32_t beg = w * 64;
        const uint32_t end = MIN(beg + 64, n);
        uint64_t bits = 0;
        if (flags) {
            for (uint32_t i = beg; i < end; ++i) {
                bits |= (uint64_t)((flags[i] & mask) == mask) << (i - beg);
            }
        } else if (mask == 0) {
            bits = (end - beg == 64) ? ~0ull : ((1ull << (end - beg)) - 1);
        }
        m->bits[w] = bits;
        while (bits) {
            m->atoms[count++] = beg + (uint32_t)ctz64(bits);
            bits &= bits - 1;
        }
        m->prefix[w + 1] = count;
    }
    m->atom_count = count;

    uint32_t cps = 0;
    for (uint32_t c = 0; c < cp_count; ++c) {
        const uint32_t atom = bounds->ca_atom_idx[c];
        if (atom < n && ((m->bits[atom >> 6] >> (atom & 63)) & 1)) {
            m->cps[cps++] = c;
        }
    }
    m->cp_count = cps;

    m->mask = mask;
    m->flags_version = mol->flags_version;
    m->valid = true;
    return m;
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

    // The sphere of all chains: about the center of their boxes, out to the farthest point of any of them
    vec3_t min_all = { FLT_MAX,  FLT_MAX,  FLT_MAX};
    vec3_t max_all = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    for (uint32_t c = 0; c < bounds->chain_count; ++c) {
        const vec4_t s = bounds->chain_sphere[c];
        const vec3_t r = {s.w, s.w, s.w};
        min_all = vec3_min(min_all, vec3_sub(vec3_from_vec4(s), r));
        max_all = vec3_max(max_all, vec3_add(vec3_from_vec4(s), r));
    }
    const vec3_t center = vec3_mul1(vec3_add(min_all, max_all), 0.5f);
    float radius = 0.0f;
    for (uint32_t c = 0; c < bounds->chain_count; ++c) {
        const vec4_t s = bounds->chain_sphere[c];
        radius = MAX(radius, vec3_distance(center, vec3_from_vec4(s)) + s.w);
    }
    // Rounding: a little more, so that every chain is within it as computed
    bounds->all_sphere = (vec4_t){center.x, center.y, center.z, radius * 1.0001f + 1.0e-3f};
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
        mol->spline_dirty = true;

        // The box of the positions, for the culling of copies: of all of them, or grown by a part of them
        if (count > 0) {
            vec3_t mn, mx;
            md_util_aabb_compute(mn.elem, mx.elem, xyz, NULL, NULL, count);
            if (offset == 0 && count == mol->atom_count) {
                mol->aabb_min = mn;
                mol->aabb_max = mx;
                mol->aabb_valid = true;
            } else if (mol->aabb_valid) {
                mol->aabb_min = vec3_min(mol->aabb_min, mn);
                mol->aabb_max = vec3_max(mol->aabb_max, mx);
            }
            // And of the subsets of the molecule (SUBSETS)
            subsets_aabb_update(mol, offset, count, xyz);
        }
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
        mol->spline_dirty = true;
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
        // The largest radius, for the culling of copies: it only grows, as the radii of a part are set
        float max_radius = offset == 0 && count == mol->atom_count ? 0.0f : mol->max_radius;
        for (uint32_t i = 0; i < count; ++i) {
            max_radius = MAX(max_radius, *(const float*)((const uint8_t*)radius + byte_stride * i));
        }
        mol->max_radius = max_radius;
        if (byte_stride > sizeof(float)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_RADIUS].id);
            float* radius_data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (radius_data == NULL) {
                MD_LOG_ERROR("Failed to map molecule radius buffer");
                return;
            }
            for (uint32_t i = 0; i < count; ++i) {
                radius_data[offset + i] = *(const float*)((const uint8_t*)radius + byte_stride * i);
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
        // Kept on the CPU too, for the draws with an atom mask (MARKS)
        if (!mol->atom_flags) {
            mol->atom_flags = md_alloc(md_get_heap_allocator(), mol->atom_count * sizeof(uint8_t));
            MEMSET(mol->atom_flags, 0, mol->atom_count * sizeof(uint8_t));
        }
        for (uint32_t i = 0; i < count; ++i) {
            mol->atom_flags[offset + i] = flag_data[byte_stride * i];
        }
        mol->flags_version += 1;
        if (byte_stride > sizeof(uint8_t)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
            uint8_t* data = (uint8_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (data == NULL) {
                MD_LOG_ERROR("Failed to map molecule flags buffer");
                return;
            }
            for (uint32_t i = 0; i < count; ++i) {
                data[offset + i] = *(flag_data + byte_stride * i);
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

        byte_stride = MAX(sizeof(md_atom_pair_t), byte_stride);
        if (offset == 0) {
            gl_buffer_set_data(mol->buffer[GL_BUFFER_BOND_ATOM_INDICES], count * sizeof(md_atom_pair_t), bonds);
            mol->bond_count = count;
            blocks_init(&mol->bond_blocks, count, 2, false);
        }
        // The blocks: of all bonds, or grown by those of a part
        for (uint32_t i = 0; i < count; ++i) {
            const md_atom_pair_t* bond = (const md_atom_pair_t*)((const uint8_t*)bonds + byte_stride * i);
            blocks_clear_prim(&mol->bond_blocks, offset + i);
            blocks_add(&mol->bond_blocks, offset + i, (uint32_t)bond->idx[0]);
            blocks_add(&mol->bond_blocks, offset + i, (uint32_t)bond->idx[1]);
        }
        if (offset == 0) {
            return;
        }

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
        mol->spline_dirty = true;
    }
}

void md_gl_mol_reset_backbone_history(md_gl_mol_t handle) {
    molecule_t* mol = mol_lookup(handle.id);
    if (mol) {
        mol->backbone_orient_history = false;
        mol->spline_dirty = true;
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
    mol->spline_dirty = true;
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
    mol->spline_dirty = true;
}

// A vertex and a fragment shader: no geometry shaders, which cost more than they save, the more so where they are
// emulated (Apple's GL over Metal). Primitives are made in the vertex shader instead (QUADS). defines is injected after
// the version line of both stages.
static bool create_permuted_program(str_t identifier, gl_program_t* program_permutations, str_t vert_src, str_t frag_src, str_t frag_output_src, str_t defines) {
    GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

    const str_t perm_str[] = {
        STR_INIT("#define ORTHO 0"),
        STR_INIT("#define ORTHO 1"),
    };
    
    ASSERT(ARRAY_SIZE(perm_str) <= MAX_SHADER_PERMUTATIONS);

    for (uint32_t perm = 0; perm < MAX_SHADER_PERMUTATIONS; ++perm) {
        // The vertex shader takes the instances (INSTANCES), the fragment shader the output snippet
        md_gl_shader_src_injection_t injections[] = { {perm_str[perm], {0}}, {defines, {0}}, {instance_shader_src, STR_INIT("INSTANCES")}, {frag_output_src, STR_INIT("EXTRA_SRC")} };

        if (!str_empty(vert_src) && !md_gl_shader_compile(vert_shader, vert_src, injections, 3)) {
            MD_LOG_ERROR("Error occured when compiling vertex shader for: '%.*s'", STR_ARG(identifier));
            return false;
        }

        if (!str_empty(frag_src) && !md_gl_shader_compile(frag_shader, frag_src, injections, 4)) {
            MD_LOG_ERROR("Error occured when compiling fragment shader for: '%.*s'", STR_ARG(identifier));
            return false;
        }

        program_permutations[perm].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, frag_shader};
        if (!md_gl_program_attach_and_link(program_permutations[perm].id, shaders, 2)) return false;

        // The textures of the caller's that the output snippet samples: their units are set once, here, as no draw sets them.
        // So are those of the instances and the lists of a draw, which every draw binds to the same units.
        glUseProgram(program_permutations[perm].id);
        for (int i = 0; i < MD_GL_USER_TEXTURE_COUNT; ++i) {
            char name[32];
            snprintf(name, sizeof(name), "md_user_texture_%d", i);
            const GLint loc = glGetUniformLocation(program_permutations[perm].id, name);
            if (loc != -1) glUniform1i(loc, MD_GL_USER_TEXTURE_UNIT + i);
        }
        {
            const GLint loc_instance = glGetUniformLocation(program_permutations[perm].id, "u_buf_instance");
            const GLint loc_list     = glGetUniformLocation(program_permutations[perm].id, "u_buf_list");
            if (loc_instance != -1) glUniform1i(loc_instance, INSTANCE_TEXTURE_UNIT);
            if (loc_list != -1)     glUniform1i(loc_list, LIST_TEXTURE_UNIT);
        }
        glUseProgram(0);
    }

    glDeleteShader(vert_shader);
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
    glGenBuffers(UBO_RING, ctx.ubo);
    for (uint32_t i = 0; i < UBO_RING; ++i) {
        glBindBuffer(GL_UNIFORM_BUFFER, ctx.ubo[i]);
        glBufferData(GL_UNIFORM_BUFFER, UBO_SIZE, NULL, GL_STREAM_DRAW);
    }
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    ctx.ubo_next = 0;

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
        STATIC_ASSERT(QUAD_BLOCK * 4 <= 65536, "Quad indices must fit in u16");
        STATIC_ASSERT(LIST_TEXTURE_UNIT < MD_GL_USER_TEXTURE_UNIT && INSTANCE_TEXTURE_UNIT < MD_GL_USER_TEXTURE_UNIT, "The units of the instances and the lists must not be the caller's");
        const uint32_t count = QUAD_BLOCK * 6;
        uint16_t* indices = md_alloc(md_get_heap_allocator(), count * sizeof(uint16_t));
        for (uint32_t q = 0; q < QUAD_BLOCK; ++q) {
            const uint16_t base = (uint16_t)(q * 4);
            indices[q * 6 + 0] = base + 0;
            indices[q * 6 + 1] = base + 1;
            indices[q * 6 + 2] = base + 2;
            indices[q * 6 + 3] = base + 2;
            indices[q * 6 + 4] = base + 1;
            indices[q * 6 + 5] = base + 3;
        }
        ctx.quad_index = gl_buffer_create(count * sizeof(uint16_t), indices, GL_STATIC_DRAW);
        md_free(md_get_heap_allocator(), indices, count * sizeof(uint16_t));
    }

    // Filled by each draw (md_gl_draw)
    ctx.instance_buf = gl_buffer_create(INSTANCE_TEXELS * sizeof(vec4_t), NULL, GL_STREAM_DRAW);
    ctx.list_buf     = gl_buffer_create(sizeof(uint32_t), NULL, GL_STREAM_DRAW);

    if (!ctx.arena) {
        ctx.arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    }

    md_handle_pool_init(&ctx.shader_pool,           MAX_SHADERS,            ctx.arena);
    md_handle_pool_init(&ctx.molecule_pool,         MAX_MOLECULES,          ctx.arena);
    md_handle_pool_init(&ctx.palette_pool,          MAX_PALETTES,           ctx.arena);
    md_handle_pool_init(&ctx.representation_pool,   MAX_REPRESENTATIONS,    ctx.arena);
    md_handle_pool_init(&ctx.subset_pool,           MAX_SUBSETS,            ctx.arena);
}

void md_gl_shutdown(void) {
    if (ctx.vao) glDeleteVertexArrays(1, &ctx.vao);
    if (ctx.fbo) glDeleteFramebuffers(1, &ctx.fbo);
    if (ctx.ubo[0]) glDeleteBuffers(UBO_RING, ctx.ubo);
    MEMSET(ctx.ubo, 0, sizeof(ctx.ubo));
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
    gl_buffer_conditional_delete(&ctx.quad_index);
    gl_buffer_conditional_delete(&ctx.instance_buf);
    gl_buffer_conditional_delete(&ctx.list_buf);

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
    if (!create_permuted_program(STR_LIT("SpaceFill"), shaders->spacefill,  (str_t){(const char*)spacefill_vert, spacefill_vert_size}, (str_t){(const char*)spacefill_frag, spacefill_frag_size}, str, (str_t){0})) return handle;
    if (!create_permuted_program(STR_LIT("Licorice"),  shaders->licorice,   bond_vert_src, bond_frag_src, str, STR_LIT("#define WEAK_BONDS 0"))) return handle;
    if (!create_permuted_program(STR_LIT("WeakBonds"), shaders->weak_bonds, bond_vert_src, bond_frag_src, str, STR_LIT("#define WEAK_BONDS 1"))) return handle;
    if (!create_permuted_program(STR_LIT("BondPartners"), shaders->bond_partners, bond_vert_src, bond_frag_src, str, STR_LIT("#define BOND_PARTNERS 1"))) return handle;
    if (!create_permuted_program(STR_LIT("ChainGaps"), shaders->chain_gaps, bond_vert_src, bond_frag_src, str, STR_LIT("#define CHAIN_GAPS 1"))) return handle;
    if (!create_permuted_program(STR_LIT("AromaticRings"), shaders->aromatic_rings, ring_vert_src, ring_frag_src, str, (str_t){0})) return handle;
    if (!create_permuted_program(STR_LIT("Ribbons"),   shaders->ribbons,    backbone_vert_src, backbone_frag_src, str, STR_LIT("#define REP_RIBBONS 1"))) return handle;
    if (!create_permuted_program(STR_LIT("Cartoon"),   shaders->cartoon,    backbone_vert_src, backbone_frag_src, str, STR_LIT("#define REP_RIBBONS 0"))) return handle;
    if (!create_permuted_program(STR_LIT("Tubes"),     shaders->tubes,      backbone_vert_src, backbone_frag_src, str, STR_LIT("#define REP_TUBES 1"))) return handle;
    if (!create_permuted_program(STR_LIT("Slabs"),     shaders->slabs,      (str_t){(const char*)slab_vert, slab_vert_size}, (str_t){(const char*)slab_frag, slab_frag_size}, str, (str_t){0})) return handle;

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
        }
        md_handle_pool_free_slot(&ctx.shader_pool, handle.id);
    }
}

// ### BOND ORDERS ###
//
// What bond.vert needs to show the order of a bond: the plane of the bond, given by its surroundings, atoms bonded to
// either end which lie in it. The plane is that of the sp2 centers: the substituents of a double bond, the ring of an
// aromatic bond, the group of a delocalized one.
//   - Double and triple bonds have a kind per bond, which replaces the stick, and their surroundings: up to two atoms
//     per end, so that bond.vert can tell how planar the end is in the frame shown (bo_plane_refs).
//   - Aromatic rings are a circle inside the ring (ring.vert), an entry per ring: its index and its atoms. Their bonds
//     keep the stick. Two fused rings have a circle each. The circle is solid where the ring holds a sextet of its
//     own (benzene, the benzene ring of tetralin), dashed in the rings of a fused aromatic system, which share their
//     sextets (MD_RING_FLAG_FUSED_AROMATIC: naphthalene, indole, a purine): it says that they are aromatic, and claims
//     no more about which ring holds which electrons.
//   - Delocalized bonds keep the stick and get dashed partners, a list of entries of their own: a bond and the
//     surroundings which give the side of the partner. So does an aromatic bond in no circle (in a ring larger than
//     MD_RING_MAX_SIZE, or of a system without ring flags): no aromatic bond goes unmarked.
//
// BEGIN BOND ORDER DERIVATION (plain arrays in, the contents of the buffers out; tested on its own)

// Kinds of bonds in GL_BUFFER_BOND_ORDER (bond.vert). An elastic bond (MD_BOND_KIND_ELASTIC) is drawn as a plain one
// where the representation shows the elastic bonds, and not at all elsewhere. A metallic bond (MD_BOND_KIND_METALLIC)
// is never drawn by licorice and ball and stick: a Bonds representation draws it.
enum {
    BOND_KIND_PLAIN    = 0,
    BOND_KIND_DOUBLE   = 1,
    BOND_KIND_TRIPLE   = 2,
    BOND_KIND_ELASTIC  = 3,
    BOND_KIND_METALLIC = 4,
};
#define BOND_REF_NONE 0xFFFFFFFFu

// Entries of GL_BUFFER_AROMATIC_RING (ring.vert): u32[8], the index of the ring | dashed << 28 | its atom count << 29,
// then its atoms in order around it and RING_ATOM_NONE past the count. Seven atoms at most, MD_RING_MAX_SIZE. Dashed:
// the circle of a ring of a fused aromatic system (MD_RING_FLAG_FUSED_AROMATIC).
#define RING_ENTRY_U32  8
#define RING_ATOM_NONE  0xFFFFFFFFu
#define RING_DASHED_SHIFT 28
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
    const uint8_t*  hydrogen;       // [atom_count] 1 for a hydrogen (optional): the last choice of surroundings
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
    return a < in->atom_count && b < in->atom_count && a != b && md_bond_kind(in->flags[i]) == MD_BOND_KIND_COVALENT;
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

// The atoms bonded to end e of the bond to o which give the plane of the bond, up to two (out, BOND_REF_NONE past those
// found): off the axis, through bonds that are not coordination bonds, heavy atoms before hydrogens (whose motion out of
// the plane is the largest). An end which is not sp2 (the S of a sulfone, the P of a phosphine oxide) takes those
// through a multiple bond when it has any, so the bond lies in the plane they share (O=S=O). A linear end (sp) gives
// none. bond.vert weighs an end of two by how planar it is in the frame shown: the hybridization is of the bonds, the
// plane of the coordinates, and an end that is pyramidal there gives no plane, whatever its order says.
static void bo_plane_refs(const bond_order_input_t* in, uint32_t e, uint32_t o, uint32_t out[2]) {
    out[0] = BOND_REF_NONE;
    out[1] = BOND_REF_NONE;
    const uint8_t hyb = in->hybridization ? in->hybridization[e] : (uint8_t)MD_HYBRIDIZATION_UNKNOWN;
    if (hyb == MD_HYBRIDIZATION_SP) return;
    // In order of preference: 0 through a multiple bond (an end which is not sp2, and then only those), 1 a heavy atom,
    // 2 a hydrogen
    uint32_t n = 0;
    for (int rank = 0; rank < 3 && n < 2; ++rank) {
        for (uint32_t j = in->adj_offset[e]; j < in->adj_offset[e + 1] && n < 2; ++j) {
            const uint32_t c = in->adj_atom[j];
            const uint32_t f = in->flags[in->adj_bond[j]];
            if (c == o || md_bond_kind(f) == MD_BOND_KIND_COORDINATE || !bo_off_axis(in, e, o, c)) continue;
            const int r = (hyb != MD_HYBRIDIZATION_SP2 && bo_is_multiple(f)) ? 0 : (in->hydrogen && in->hydrogen[c]) ? 2 : 1;
            if (r == rank) out[n++] = c;
        }
        if (rank == 0 && n > 0) return;
    }
}

// The first of bo_plane_refs, for the side of a dashed partner
static uint32_t bo_plane_ref(const bond_order_input_t* in, uint32_t e, uint32_t o) {
    uint32_t ref[2];
    bo_plane_refs(in, e, o, ref);
    return ref[0];
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
        if (c != o && (f & MD_BOND_FLAG_AROMATIC) && md_bond_kind(f) != MD_BOND_KIND_COORDINATE && bo_off_axis(in, e, o, c)) return c;
    }
    return bo_plane_ref(in, e, o);
}

// The ring r is drawn as a circle: aromatic (MD_RING_FLAG_AROMATIC, which says that every bond of it is), of at most
// RING_ENTRY_U32 - 1 atoms, each bonded to the next by a bond the molecule draws, and its index fits in its entry
static bool bo_circle_ring(const bond_order_input_t* in, uint32_t r) {
    if (!in->ring_flags || !(in->ring_flags[r] & MD_RING_FLAG_AROMATIC) || r >= (1u << RING_DASHED_SHIFT)) return false;
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
            const uint32_t dashed = (in->ring_flags[r] & MD_RING_FLAG_FUSED_AROMATIC) ? 1u : 0u;
            uint32_t* e = out_ring + RING_ENTRY_U32 * count;
            e[0] = r | (dashed << RING_DASHED_SHIFT) | (n << RING_COUNT_SHIFT);
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
// bonds. out_ref [4 * bond_count]: per such bond the atoms bonded to its first atom which give its plane, then those
// bonded to its second atom, two each (bo_plane_refs), BOND_REF_NONE where there are fewer. Returns the number of
// entries of out_ref. Aromatic and delocalized bonds are plain here: they keep their stick (compute_bond_partners).
static uint32_t compute_bond_orders(uint32_t* out_order, uint32_t* out_ref, const bond_order_input_t* in) {
    uint32_t entries = 0;
    for (uint32_t i = 0; i < in->bond_count; ++i) {
        const uint32_t f = in->flags[i];
        const md_bond_kind_t kind = md_bond_kind(f);
        out_order[i] = kind == MD_BOND_KIND_ELASTIC ? BOND_KIND_ELASTIC : kind == MD_BOND_KIND_METALLIC ? BOND_KIND_METALLIC : BOND_KIND_PLAIN;
        if (!bo_valid_bond(in, i) || (f & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED))) continue;
        const uint32_t a = in->pairs[2 * i + 0];
        const uint32_t b = in->pairs[2 * i + 1];
        const int order = md_bond_order((md_bond_flags_t)f);
        uint32_t* ref = out_ref + 4 * entries;
        if (order == 2) {
            bo_plane_refs(in, a, b, ref + 0);
            bo_plane_refs(in, b, a, ref + 2);
            out_order[i] = BOND_KIND_DOUBLE | (entries << 3);
            ++entries;
        } else if (order >= 3) {
            // No plane, perpendicular to the view
            ref[0] = ref[1] = ref[2] = ref[3] = BOND_REF_NONE;
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

void md_gl_mol_set_bond_orders(md_gl_mol_t handle, const md_system_t* sys, const md_bond_delta_t* delta, const md_chem_delta_t* chem) {
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
    // The bonds of the graph come first in the molecule, then those formed in a frame (delta), as many as it has
    const uint32_t sys_bond_count = (uint32_t)sys->bond.count;
    if (sys_bond_count > mol->bond_count) {
        MD_LOG_ERROR("The system has %u bonds and the molecule %u: set the bonds first (md_gl_mol_set_bonds)", sys_bond_count, mol->bond_count);
        return;
    }
    const uint32_t mol_bond_count = mol->bond_count;
    if (sys_bond_count == 0 || !sys->bond.pairs || !sys->bond.flags) {
        return;
    }
    const uint32_t num_formed = delta ? (uint32_t)MIN(delta->num_formed, (size_t)(mol_bond_count - sys_bond_count)) : 0;
    const uint32_t bond_count = sys_bond_count + num_formed;
    if (bond_count >= (1u << 29)) {
        MD_LOG_ERROR("Too many bonds for bond orders");
        return;
    }
    STATIC_ASSERT(sizeof(md_atom_pair_t) == 2 * sizeof(uint32_t), "Bond pairs are read as u32[2]");

    md_allocator_i* alloc = md_get_heap_allocator();
    const uint32_t atom_count = mol->atom_count;

    // The bonds as the frame has them: a broken one is set to one atom twice (no bond), the formed ones follow, and the
    // region of the frame's chemistry has its own flags
    uint32_t* pairs = md_alloc(alloc, 2 * MAX(bond_count, 1) * sizeof(uint32_t));
    uint32_t* flags = md_alloc(alloc, MAX(bond_count, 1) * sizeof(uint32_t));
    MEMCPY(pairs, sys->bond.pairs, 2 * sys_bond_count * sizeof(uint32_t));
    for (uint32_t i = 0; i < sys_bond_count; ++i) flags[i] = (uint32_t)sys->bond.flags[i];
    if (delta) {
        md_bitset_iter_t it = md_bitset_iter(delta->broken);
        while (md_bitset_iter_next(&it)) {
            if (it.idx < sys_bond_count) pairs[2 * it.idx + 1] = pairs[2 * it.idx + 0];
        }
        for (uint32_t j = 0; j < num_formed; ++j) {
            pairs[2 * (sys_bond_count + j) + 0] = (uint32_t)delta->formed_pair[j].idx[0];
            pairs[2 * (sys_bond_count + j) + 1] = (uint32_t)delta->formed_pair[j].idx[1];
            flags[sys_bond_count + j] = (uint32_t)delta->formed_flags[j];
        }
    }
    uint8_t* hyb      = md_alloc(alloc, MAX(atom_count, 1) * sizeof(uint8_t));
    uint8_t* hydrogen = md_alloc(alloc, MAX(atom_count, 1) * sizeof(uint8_t));
    for (uint32_t i = 0; i < atom_count; ++i) {
        hyb[i]      = (uint8_t)md_atom_hybridization(&sys->atom, i);
        hydrogen[i] = md_atom_atomic_number(&sys->atom, i) == 1 ? 1 : 0;
    }
    uint8_t* region = NULL;     // [atom_count] Atoms of the region of the frame's chemistry
    if (chem && chem->num_atoms) {
        region = md_alloc(alloc, MAX(atom_count, 1));
        MEMSET(region, 0, MAX(atom_count, 1));
        for (size_t k = 0; k < chem->num_atoms; ++k) {
            const int32_t a = chem->atom_idx[k];
            if (a < 0 || (uint32_t)a >= atom_count) continue;
            region[a] = 1;
            hyb[a] = (uint8_t)md_atom_flags_hybridization(chem->atom_flags[k]);
        }
        for (size_t k = 0; k < chem->num_bonds; ++k) {
            if (chem->bond_idx[k] < bond_count) flags[chem->bond_idx[k]] = (uint32_t)chem->bond_flags[k];
        }
    }

    // Connectivity from the pairs, which are what the molecule draws. The elastic bonds are no part of it: they give
    // no plane to a double bond and close no ring.
    uint32_t* adj_offset = md_alloc(alloc, (atom_count + 1) * sizeof(uint32_t));
    uint32_t* adj_atom   = md_alloc(alloc, 2 * MAX(bond_count, 1) * sizeof(uint32_t));
    uint32_t* adj_bond   = md_alloc(alloc, 2 * MAX(bond_count, 1) * sizeof(uint32_t));
    uint32_t* order      = md_alloc(alloc, MAX(mol_bond_count, 1) * sizeof(uint32_t));
    uint32_t* ref        = md_alloc(alloc, 4 * MAX(bond_count, 1) * sizeof(uint32_t));

    MEMSET(adj_offset, 0, (atom_count + 1) * sizeof(uint32_t));
    for (uint32_t i = 0; i < bond_count; ++i) {
        const uint32_t a = pairs[2 * i + 0];
        const uint32_t b = pairs[2 * i + 1];
        if (a < atom_count && b < atom_count && a != b && md_bond_is_graph(flags[i])) {
            adj_offset[a + 1] += 1;
            adj_offset[b + 1] += 1;
        }
    }
    for (uint32_t i = 0; i < atom_count; ++i) {
        adj_offset[i + 1] += adj_offset[i];
    }
    uint32_t* cursor = md_alloc(alloc, MAX(atom_count, 1) * sizeof(uint32_t));
    MEMCPY(cursor, adj_offset, atom_count * sizeof(uint32_t));
    for (uint32_t i = 0; i < bond_count; ++i) {
        const uint32_t a = pairs[2 * i + 0];
        const uint32_t b = pairs[2 * i + 1];
        if (a < atom_count && b < atom_count && a != b && md_bond_is_graph(flags[i])) {
            adj_atom[cursor[a]] = b;
            adj_bond[cursor[a]++] = i;
            adj_atom[cursor[b]] = a;
            adj_bond[cursor[b]++] = i;
        }
    }
    md_free(alloc, cursor, MAX(atom_count, 1) * sizeof(uint32_t));

    // The rings of the system keep their indices, those with an atom in the region drawn no more; the rings of the
    // region follow them. Without flags for every ring there are no circles, and the aromatic bonds get dashed partners.
    const bool has_xyz = md_system_state_has_coords(&sys->reference) && sys->reference.num_atoms >= atom_count;
    const uint32_t sys_ring_count = (uint32_t)md_index_data_num_ranges(&sys->ring);
    const bool has_ring_flags = sys_ring_count && sys->ring_flags && md_array_size(sys->ring_flags) == sys_ring_count;
    const uint32_t chem_ring_count = region ? (uint32_t)md_index_data_num_ranges(&chem->ring) : 0;
    const uint32_t ring_count = sys_ring_count + chem_ring_count;
    const uint32_t* ring_offset = ring_count ? sys->ring.offsets : NULL;
    const int32_t*  ring_atom   = ring_count ? sys->ring.indices : NULL;
    const md_ring_flags_t* ring_flags = has_ring_flags ? sys->ring_flags : NULL;
    uint32_t* own_offset = NULL;
    int32_t*  own_atom = NULL;
    md_ring_flags_t* own_flags = NULL;
    size_t own_atom_count = 0;
    if (chem_ring_count) {
        const uint32_t sys_ring_atoms = sys_ring_count ? sys->ring.offsets[sys_ring_count] : 0;
        own_atom_count = sys_ring_atoms + md_array_size(chem->ring.indices);
        own_offset = md_alloc(alloc, (ring_count + 1) * sizeof(uint32_t));
        own_atom   = md_alloc(alloc, MAX(own_atom_count, 1) * sizeof(int32_t));
        own_flags  = md_alloc(alloc, ring_count * sizeof(md_ring_flags_t));
        if (sys_ring_count) {
            MEMCPY(own_offset, sys->ring.offsets, (sys_ring_count + 1) * sizeof(uint32_t));
            MEMCPY(own_atom, sys->ring.indices, sys_ring_atoms * sizeof(int32_t));
        } else {
            own_offset[0] = 0;
        }
        for (uint32_t r = 0; r < sys_ring_count; ++r) {
            md_ring_flags_t f = has_ring_flags ? sys->ring_flags[r] : MD_RING_FLAG_NONE;
            for (uint32_t k = own_offset[r]; k < own_offset[r + 1]; ++k) {
                if (own_atom[k] >= 0 && (uint32_t)own_atom[k] < atom_count && region[own_atom[k]]) f = MD_RING_FLAG_NONE;
            }
            own_flags[r] = f;
        }
        for (uint32_t q = 0; q < chem_ring_count; ++q) {
            const int32_t* rq = md_index_range_beg(&chem->ring, q);
            const uint32_t n = (uint32_t)md_index_range_size(&chem->ring, q);
            const uint32_t at = own_offset[sys_ring_count + q];
            MEMCPY(own_atom + at, rq, n * sizeof(int32_t));
            own_offset[sys_ring_count + q + 1] = at + n;
            own_flags[sys_ring_count + q] = chem->ring_flags[q];
        }
        ring_offset = own_offset;
        ring_atom   = own_atom;
        ring_flags  = own_flags;
    }

    const bond_order_input_t in = {
        .bond_count    = bond_count,
        .pairs         = pairs,
        .flags         = flags,
        .atom_count    = atom_count,
        .adj_offset    = adj_offset,
        .adj_atom      = adj_atom,
        .adj_bond      = adj_bond,
        .hybridization = hyb,
        .hydrogen      = hydrogen,
        .xyz           = has_xyz ? (const float*)sys->reference.xyz : NULL,
        .ring_count    = ring_count,
        .ring_offset   = ring_offset,
        .ring_atom     = ring_atom,
        .ring_flags    = ring_flags,
    };
    const uint32_t entries = compute_bond_orders(order, ref, &in);

    // The circles first: the bonds they cover need no dashed partner
    uint8_t* covered = md_alloc(alloc, MAX(bond_count, 1) * sizeof(uint8_t));
    MEMSET(covered, 0, MAX(bond_count, 1) * sizeof(uint8_t));
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
    // The bonds past those given: plain (KIND_PLAIN, 0)
    if (mol_bond_count > bond_count) MEMSET(order + bond_count, 0, (mol_bond_count - bond_count) * sizeof(uint32_t));
    mol->buffer[GL_BUFFER_BOND_ORDER]     = gl_buffer_create(mol_bond_count * sizeof(uint32_t), order, GL_STATIC_DRAW);
    mol->buffer[GL_BUFFER_BOND_ORDER_REF] = gl_buffer_create(MAX(entries, 1) * 4 * sizeof(uint32_t), ref, GL_STATIC_DRAW);
    mol->buffer[GL_BUFFER_BOND_PARTNER]   = gl_buffer_create(MAX(partner_count, 1) * 4 * sizeof(uint32_t), partner, GL_STATIC_DRAW);
    mol->buffer[GL_BUFFER_AROMATIC_RING]  = gl_buffer_create((uint32_t)ring_entry_bytes, ring_entry, GL_STATIC_DRAW);
    mol->bond_order_count    = mol_bond_count;
    mol->bond_partner_count  = partner_count;
    mol->aromatic_ring_count = ring_entry_count;

    // The blocks (MARKS): a partner is shown when the atoms of its bond are, a circle when all atoms of its ring are
    blocks_init(&mol->partner_blocks, partner_count, 2, false);
    for (uint32_t i = 0; i < partner_count; ++i) {
        const uint32_t bond = partner[4 * i + 0];
        if (bond < bond_count) {
            blocks_add(&mol->partner_blocks, i, pairs[2 * bond + 0]);
            blocks_add(&mol->partner_blocks, i, pairs[2 * bond + 1]);
        }
    }
    blocks_init(&mol->ring_blocks, ring_entry_count, RING_ENTRY_U32 - 1, false);
    for (uint32_t i = 0; i < ring_entry_count; ++i) {
        const uint32_t* e = ring_entry + i * RING_ENTRY_U32;
        const uint32_t n = MIN(e[0] >> RING_COUNT_SHIFT, RING_ENTRY_U32 - 1);
        for (uint32_t k = 0; k < n; ++k) {
            blocks_add(&mol->ring_blocks, i, e[1 + k]);
        }
    }
    md_free(alloc, partner, MAX(partner_count, 1) * 4 * sizeof(uint32_t));
    md_free(alloc, ring_entry, ring_entry_bytes);
    md_free(alloc, covered, MAX(bond_count, 1) * sizeof(uint8_t));
    if (own_offset) {
        md_free(alloc, own_offset, (ring_count + 1) * sizeof(uint32_t));
        md_free(alloc, own_atom, MAX(own_atom_count, 1) * sizeof(int32_t));
        md_free(alloc, own_flags, ring_count * sizeof(md_ring_flags_t));
    }
    if (region) md_free(alloc, region, MAX(atom_count, 1));

    md_free(alloc, adj_offset, (atom_count + 1) * sizeof(uint32_t));
    md_free(alloc, adj_atom,   2 * MAX(bond_count, 1) * sizeof(uint32_t));
    md_free(alloc, adj_bond,   2 * MAX(bond_count, 1) * sizeof(uint32_t));
    md_free(alloc, pairs,      2 * MAX(bond_count, 1) * sizeof(uint32_t));
    md_free(alloc, flags,      MAX(bond_count, 1) * sizeof(uint32_t));
    md_free(alloc, hyb,        MAX(atom_count, 1) * sizeof(uint8_t));
    md_free(alloc, hydrogen,   MAX(atom_count, 1) * sizeof(uint8_t));
    md_free(alloc, order,      MAX(mol_bond_count, 1) * sizeof(uint32_t));
    md_free(alloc, ref,        4 * MAX(bond_count, 1) * sizeof(uint32_t));
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
    blocks_free(&mol->gap_trace_blocks);
    blocks_free(&mol->gap_link_blocks);
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

        // The blocks of each kind (MARKS): a gap is shown when its two ends are
        blocks_init(&mol->gap_trace_blocks, trace_count, 2, false);
        blocks_init(&mol->gap_link_blocks, link_count, 2, false);
        for (uint32_t i = 0; i < total; ++i) {
            const uint32_t* e = trace + 4 * i;
            prim_blocks_t* blocks = i < trace_count ? &mol->gap_trace_blocks : &mol->gap_link_blocks;
            const uint32_t prim = i < trace_count ? i : i - trace_count;
            blocks_add(blocks, prim, e[0]);
            blocks_add(blocks, prim, e[2]);
        }
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
    blocks_free(&mol->cp_blocks);
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
    mol->spline_dirty = true;
    // The control points of the marks are of the previous backbone
    mol->marked.valid = false;

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

    // The control points as primitives of their atoms (SUBSETS)
    blocks_init(&mol->cp_blocks, cp_count, 1, true);
    for (uint32_t c = 0; c < cp_count; ++c) {
        blocks_add(&mol->cp_blocks, c, bounds_atom[c]);
    }

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
            md_gl_mol_set_bond_orders(handle, sys, NULL, NULL);
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
    // All at the origin
    gl_mol->aabb_valid = true;
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

// Where particle p is, from the positions of the atoms of the source on the CPU, as compute_derived.vert places it
static vec3_t derived_position(const molecule_t* mol, uint32_t p, const vec3_t* xyz, const mat3_t* A, const mat3_t* I, vec3_t mask) {
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
    return vec3_add(vec3_mul1(ref, sum_w), sum);
}

// Where the particles are, from the positions of the atoms of the source on the CPU: the box of all of them, for the
// culling of copies, and the control points of the tubes, for their culling and level of detail. Each particle once.
static void derived_bounds_update(molecule_t* mol, const vec3_t* xyz, const mat3_t* A, const mat3_t* I, vec3_t mask) {
    if (!mol->derived.data || mol->atom_count == 0) return;
    md_temp_scope_t temp = md_temp_begin();
    vec3_t* pos = md_temp_alloc(temp, mol->atom_count * sizeof(vec3_t));
    vec3_t mn = { FLT_MAX,  FLT_MAX,  FLT_MAX};
    vec3_t mx = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    for (uint32_t p = 0; p < mol->atom_count; ++p) {
        pos[p] = derived_position(mol, p, xyz, A, I, mask);
        mn = vec3_min(mn, pos[p]);
        mx = vec3_max(mx, pos[p]);
    }
    mol->aabb_min = mn;
    mol->aabb_max = mx;
    mol->aabb_valid = true;
    subsets_aabb_update(mol, 0, mol->atom_count, pos);

    backbone_bounds_t* bounds = &mol->backbone_bounds;
    if (bounds->chain_count) {
        const uint32_t cp_count = bounds->chain_offset[bounds->chain_count];
        for (uint32_t i = 0; i < cp_count; ++i) {
            const uint32_t p = bounds->ca_atom_idx[i];
            if (p < mol->atom_count) bounds->ca_xyz[i] = pos[p];
        }
        backbone_bounds_update_spheres(bounds);
    }
    md_temp_end(temp);
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
    mol->spline_dirty = true;

    if (source_xyz) {
        derived_bounds_update(mol, source_xyz, &A, &I, mask);
    } else {
        // Moved, to where is not known on the CPU: not culled
        mol->aabb_valid = false;
        subsets_aabb_invalidate(mol);
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
    const uint32_t* side     = count ? tubes->side : NULL;
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
        if (particle[base + j] >= mol->atom_count || (face && face[base + j] != UINT32_MAX && face[base + j] >= mol->atom_count) ||
            (side && side[base + j] != UINT32_MAX && side[base + j] >= mol->atom_count)) {
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
            const uint32_t sd = (f == UINT32_MAX && side) ? side[base + j] : UINT32_MAX;
            const float    r = CLAMP(radius[base + j], 0.0f, TUBE_RADIUS_MAX);
            const float    l = flat ? CLAMP(flat[base + j], 0.0f, TUBE_RADIUS_MAX) : 0.0f;
            data[j] = (gl_backbone_data_t){
                .bb_seg_idx = j,
                .ca         = p,
                .c          = p,
                .o          = f != UINT32_MAX ? f : (sd != UINT32_MAX ? sd : p),
                .flags      = (j == beg ? BACKBONE_FLAG_CHAIN_BEG : 0) | (j == end - 1 ? BACKBONE_FLAG_CHAIN_END : 0) |
                              (f != UINT32_MAX ? BACKBONE_FLAG_FIXED_FACE : 0) | (sd != UINT32_MAX ? BACKBONE_FLAG_FIXED_SIDE : 0) |
                              ((through && through[t]) ? BACKBONE_FLAG_THROUGH : 0),
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
    blocks_free(&mol->slab_blocks);
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
    // The blocks (MARKS): a slab is shown when all its corners are
    blocks_init(&mol->slab_blocks, count, SLAB_MAX_CORNERS, false);
    for (uint32_t i = 0; i < count; ++i) {
        for (uint32_t j = offset[i]; j < offset[i + 1]; ++j) {
            blocks_add(&mol->slab_blocks, i, particle[j]);
        }
    }
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
        if (mol->atom_flags) md_free(md_get_heap_allocator(), mol->atom_flags, mol->atom_count * sizeof(uint8_t));
        marked_free(&mol->marked);
        blocks_free(&mol->bond_blocks);
        blocks_free(&mol->partner_blocks);
        blocks_free(&mol->ring_blocks);
        blocks_free(&mol->gap_trace_blocks);
        blocks_free(&mol->gap_link_blocks);
        blocks_free(&mol->slab_blocks);
        blocks_free(&mol->cp_blocks);
        if (mol->scratch_bits) md_free(md_get_heap_allocator(), mol->scratch_bits, DIV_UP(mol->atom_count, 64) * sizeof(uint64_t));
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
        blocks_free(&rep->weak_bond_blocks);
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
    rep->weak_bond_max_weight = 0.0f;
    blocks_free(&rep->weak_bond_blocks);
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

    // The reach of the bonds beyond their atoms, for the culling of copies, and their blocks (MARKS)
    blocks_init(&rep->weak_bond_blocks, written, 2, true);
    uint32_t j = 0;
    for (uint32_t i = 0; i < count; ++i) {
        if (atom_a[i] < mol->atom_count && atom_b[i] < mol->atom_count) {
            rep->weak_bond_max_weight = MAX(rep->weak_bond_max_weight, weight ? weight[i] : 1.0f);
            blocks_add(&rep->weak_bond_blocks, j, atom_a[i]);
            blocks_add(&rep->weak_bond_blocks, j, atom_b[i]);
            ++j;
        }
    }
}

// ### SUBSETS ###

static void u32_reserve(uint32_t** arr, uint32_t* capacity, uint32_t need) {
    if (need <= *capacity) return;
    md_allocator_i* alloc = md_get_heap_allocator();
    const uint32_t cap = MAX(need, *capacity + *capacity / 2);
    uint32_t* mem = md_alloc(alloc, cap * sizeof(uint32_t));
    if (*arr) {
        MEMCPY(mem, *arr, *capacity * sizeof(uint32_t));
        md_free(alloc, *arr, *capacity * sizeof(uint32_t));
    }
    *arr = mem;
    *capacity = cap;
}

static int cmp_u32(const void* a, const void* b) {
    const uint32_t x = *(const uint32_t*)a;
    const uint32_t y = *(const uint32_t*)b;
    return (x > y) - (x < y);
}

// In order and each once
static uint32_t sort_unique_u32(uint32_t* items, uint32_t count) {
    bool sorted = true;
    for (uint32_t i = 1; i < count && sorted; ++i) sorted = items[i - 1] < items[i];
    if (sorted) return count;
    qsort(items, count, sizeof(uint32_t), cmp_u32);
    uint32_t n = 0;
    for (uint32_t i = 0; i < count; ++i) {
        if (n == 0 || items[n - 1] != items[i]) items[n++] = items[i];
    }
    return n;
}

static uint64_t* subset_bits_set(molecule_t* mol, const subset_t* s);
static void subset_bits_clear(molecule_t* mol, const subset_t* s);
static inline bool bit_test(const uint64_t* bits, uint32_t i) {
    return (bits[i >> 6] >> (i & 63)) & 1;
}

static void subset_free(subset_t* s) {
    md_allocator_i* alloc = md_get_heap_allocator();
    if (s->atoms) md_free(alloc, s->atoms, s->atom_capacity * sizeof(uint32_t));
    for (int k = 0; k < SUBSET_KIND_COUNT; ++k) {
        if (s->list[k].prims) md_free(alloc, s->list[k].prims, s->list[k].capacity * sizeof(uint32_t));
    }
    gl_buffer_conditional_delete(&s->gpu);
    MEMSET(s, 0, sizeof(subset_t));
}

md_gl_subset_t md_gl_subset_create(md_gl_mol_t mol_handle) {
    md_gl_subset_t handle = {0};
    if (!validate_context()) return handle;
    const molecule_t* mol = mol_lookup(mol_handle.id);
    if (!mol) {
        MD_LOG_ERROR("Invalid molecule");
        return handle;
    }
    const uint32_t id = md_handle_pool_alloc_slot(&ctx.subset_pool);
    if (id == MD_HANDLE_INVALID_ID) {
        MD_LOG_ERROR("Fatal error: ran out of slots for subsets!");
        return handle;
    }
    subset_t* s = ctx.subsets + md_handle_index(id);
    MEMSET(s, 0, sizeof(subset_t));
    s->id = id;
    s->mol_id = mol_handle.id;
    for (int k = 0; k < SUBSET_KIND_COUNT; ++k) {
        s->list[k].version = SUBSET_LIST_STALE;
    }
    handle.id = id;
    return handle;
}

void md_gl_subset_destroy(md_gl_subset_t handle) {
    if (handle.id == MD_HANDLE_INVALID_ID) return;
    subset_t* s = subset_lookup(handle.id);
    if (s) {
        subset_free(s);
        md_handle_pool_free_slot(&ctx.subset_pool, handle.id);
    }
}

void md_gl_subset_set_atoms(md_gl_subset_t handle, const uint32_t* atoms, uint32_t count) {
    subset_t* s = subset_lookup(handle.id);
    if (!s) {
        MD_LOG_ERROR("Invalid subset");
        return;
    }
    const molecule_t* mol = mol_lookup(s->mol_id);
    if (!mol) {
        MD_LOG_ERROR("The molecule of the subset is invalid");
        return;
    }
    if (count > 0 && !atoms) {
        MD_LOG_ERROR("Missing the atoms of the subset");
        return;
    }
    u32_reserve(&s->atoms, &s->atom_capacity, count);
    uint32_t n = 0;
    for (uint32_t i = 0; i < count; ++i) {
        if (atoms[i] < mol->atom_count) s->atoms[n++] = atoms[i];
    }
    s->atom_count = sort_unique_u32(s->atoms, n);
    for (int k = 0; k < SUBSET_KIND_COUNT; ++k) {
        s->list[k].version = SUBSET_LIST_STALE;
    }
    s->gpu_stale = true;
    s->aabb_valid = false;
}

void md_gl_subset_set_derived(md_gl_subset_t handle, md_gl_subset_t source_handle) {
    subset_t* s = subset_lookup(handle.id);
    const subset_t* src = subset_lookup(source_handle.id);
    if (!s || !src) {
        MD_LOG_ERROR("Invalid subset");
        return;
    }
    const molecule_t* mol = mol_lookup(s->mol_id);
    molecule_t* src_mol = mol_lookup(src->mol_id);
    if (!mol || !src_mol) {
        MD_LOG_ERROR("The molecule of a subset is invalid");
        return;
    }
    if (!mol->derived.data || mol->derived.source_id != src_mol->id) {
        MD_LOG_ERROR("The molecule of the subset is not derived from that of the source (md_gl_mol_set_derived)");
        return;
    }
    // The particles whose reference atom is in the source
    md_temp_scope_t temp = md_temp_begin();
    uint32_t* particles = md_temp_alloc(temp, MAX(mol->atom_count, 1) * sizeof(uint32_t));
    uint32_t count = 0;
    const uint64_t* bits = subset_bits_set(src_mol, src);
    for (uint32_t p = 0; p < mol->atom_count; ++p) {
        const uint32_t ref = mol->derived.data[3 * p + 2];
        if (ref < src_mol->atom_count && bit_test(bits, ref)) particles[count++] = p;
    }
    subset_bits_clear(src_mol, src);
    md_gl_subset_set_atoms(handle, particles, count);
    md_temp_end(temp);
}

static prim_blocks_t* mol_kind_blocks(molecule_t* mol, int kind) {
    switch (kind) {
    case SUBSET_CPS:        return &mol->cp_blocks;
    case SUBSET_BONDS:      return &mol->bond_blocks;
    case SUBSET_PARTNERS:   return &mol->partner_blocks;
    case SUBSET_RINGS:      return &mol->ring_blocks;
    case SUBSET_GAP_TRACE:  return &mol->gap_trace_blocks;
    case SUBSET_GAP_LINK:   return &mol->gap_link_blocks;
    case SUBSET_SLABS:      return &mol->slab_blocks;
    default: ASSERT(false); return NULL;
    }
}

// The bits of the atoms of the subset set in the scratch bits of its molecule, or cleared again
static uint64_t* subset_bits_set(molecule_t* mol, const subset_t* s) {
    if (!mol->scratch_bits) {
        const size_t bytes = DIV_UP(mol->atom_count, 64) * sizeof(uint64_t);
        mol->scratch_bits = md_alloc(md_get_heap_allocator(), MAX(bytes, sizeof(uint64_t)));
        MEMSET(mol->scratch_bits, 0, MAX(bytes, sizeof(uint64_t)));
    }
    for (uint32_t i = 0; i < s->atom_count; ++i) {
        const uint32_t a = s->atoms[i];
        mol->scratch_bits[a >> 6] |= 1ull << (a & 63);
    }
    return mol->scratch_bits;
}

static void subset_bits_clear(molecule_t* mol, const subset_t* s) {
    for (uint32_t i = 0; i < s->atom_count; ++i) {
        mol->scratch_bits[s->atoms[i] >> 6] = 0;
    }
}

// Whether all atoms of the primitive are in bits, and the first of them
static inline bool prim_within(const prim_blocks_t* b, uint32_t prim, const uint64_t* bits, uint32_t n, uint32_t* first) {
    const uint32_t* pa = b->atoms + (size_t)prim * b->stride;
    uint32_t lo = PRIM_ATOM_NONE;
    for (uint32_t k = 0; k < b->stride; ++k) {
        const uint32_t x = pa[k];
        if (x == PRIM_ATOM_NONE) break;
        if (x >= n || !bit_test(bits, x)) return false;
        lo = MIN(lo, x);
    }
    *first = lo;
    return lo != PRIM_ATOM_NONE;
}

// The primitives all of whose atoms are in the subset, its atoms in bits: from the primitives each of its atoms is in,
// each primitive at the first of its atoms, in order
static void subset_find(subset_list_t* l, const subset_t* s, prim_blocks_t* b, const uint64_t* bits, uint32_t n) {
    l->count = 0;
    if (b->prim_count == 0 || s->atom_count == 0) return;
    blocks_adjacency(b, n);
    for (uint32_t i = 0; i < s->atom_count; ++i) {
        const uint32_t a = s->atoms[i];
        for (uint32_t j = b->adj_offset[a]; j < b->adj_offset[a + 1]; ++j) {
            const uint32_t p = b->adj_prim[j];
            uint32_t first;
            // The same atom twice in a primitive lists it twice in a row
            if (l->count > 0 && l->prims[l->count - 1] == p) continue;
            if (prim_within(b, p, bits, n, &first) && first == a) {
                u32_reserve(&l->prims, &l->capacity, l->count + 1);
                l->prims[l->count++] = p;
            }
        }
    }
    qsort(l->prims, l->count, sizeof(uint32_t), cmp_u32);
}

// The u32 of a buffer, read back from the GPU
static void buffer_read_u32(uint32_t* out, gl_buffer_t buf, size_t count) {
    if (count == 0 || !buf.id) return;
    glBindBuffer(GL_COPY_READ_BUFFER, buf.id);
    glGetBufferSubData(GL_COPY_READ_BUFFER, 0, (GLsizeiptr)(count * sizeof(uint32_t)), out);
    glBindBuffer(GL_COPY_READ_BUFFER, 0);
}

// The atoms of the primitives of a kind, read back from the buffers of the molecule where they are not kept on the CPU
// (prim_blocks_t): once a subset needs them, and again after the primitives change (blocks_init frees them)
static const prim_blocks_t* mol_kind_atoms(molecule_t* mol, int kind) {
    prim_blocks_t* b = mol_kind_blocks(mol, kind);
    if (b->atoms || b->prim_count == 0) return b;
    md_allocator_i* heap = md_get_heap_allocator();
    const uint32_t n = b->prim_count;
    const size_t slots = (size_t)n * b->stride;
    b->atoms = md_alloc(heap, slots * sizeof(uint32_t));
    MEMSET(b->atoms, 0xFF, slots * sizeof(uint32_t));   // PRIM_ATOM_NONE
    md_temp_scope_t temp = md_temp_begin();
    switch (kind) {
    case SUBSET_BONDS:
        // The pairs as they are, two atoms each
        if (n == mol->bond_count) buffer_read_u32(b->atoms, mol->buffer[GL_BUFFER_BOND_ATOM_INDICES], 2 * (size_t)n);
        break;
    case SUBSET_PARTNERS:
        // The atoms of the bond of each (GL_BUFFER_BOND_PARTNER)
        if (n == mol->bond_partner_count && mol->bond_count > 0) {
            uint32_t* partner = md_temp_alloc(temp, 4 * (size_t)n * sizeof(uint32_t));
            uint32_t* pairs   = md_temp_alloc(temp, 2 * (size_t)mol->bond_count * sizeof(uint32_t));
            buffer_read_u32(partner, mol->buffer[GL_BUFFER_BOND_PARTNER], 4 * (size_t)n);
            buffer_read_u32(pairs, mol->buffer[GL_BUFFER_BOND_ATOM_INDICES], 2 * (size_t)mol->bond_count);
            for (uint32_t i = 0; i < n; ++i) {
                const uint32_t bond = partner[4 * i];
                if (bond >= mol->bond_count) continue;
                b->atoms[2 * i + 0] = pairs[2 * bond + 0];
                b->atoms[2 * i + 1] = pairs[2 * bond + 1];
            }
        }
        break;
    case SUBSET_RINGS:
        // Their atoms in order around them (GL_BUFFER_AROMATIC_RING)
        if (n == mol->aromatic_ring_count) {
            uint32_t* entry = md_temp_alloc(temp, RING_ENTRY_U32 * (size_t)n * sizeof(uint32_t));
            buffer_read_u32(entry, mol->buffer[GL_BUFFER_AROMATIC_RING], RING_ENTRY_U32 * (size_t)n);
            for (uint32_t i = 0; i < n; ++i) {
                const uint32_t* e = entry + RING_ENTRY_U32 * i;
                const uint32_t c = MIN(e[0] >> RING_COUNT_SHIFT, (uint32_t)(RING_ENTRY_U32 - 1));
                for (uint32_t k = 0; k < c; ++k) b->atoms[(RING_ENTRY_U32 - 1) * i + k] = e[1 + k];
            }
        }
        break;
    case SUBSET_GAP_TRACE:
    case SUBSET_GAP_LINK: {
        // Their two ends (GL_BUFFER_CHAIN_GAP), the traces first
        const uint32_t total = mol->gap_trace_count + mol->gap_link_count;
        const uint32_t base  = kind == SUBSET_GAP_TRACE ? 0 : mol->gap_trace_count;
        const uint32_t count = kind == SUBSET_GAP_TRACE ? mol->gap_trace_count : mol->gap_link_count;
        if (n == count) {
            uint32_t* entry = md_temp_alloc(temp, 4 * (size_t)total * sizeof(uint32_t));
            buffer_read_u32(entry, mol->buffer[GL_BUFFER_CHAIN_GAP], 4 * (size_t)total);
            for (uint32_t i = 0; i < n; ++i) {
                b->atoms[2 * i + 0] = entry[4 * (base + i) + 0];
                b->atoms[2 * i + 1] = entry[4 * (base + i) + 2];
            }
        }
        break;
    }
    case SUBSET_SLABS:
        // Their corners (GL_BUFFER_SLAB)
        if (n == mol->slab_count) {
            uint32_t* entry = md_temp_alloc(temp, SLAB_ENTRY_U32 * (size_t)n * sizeof(uint32_t));
            buffer_read_u32(entry, mol->buffer[GL_BUFFER_SLAB], SLAB_ENTRY_U32 * (size_t)n);
            for (uint32_t i = 0; i < n; ++i) {
                const uint32_t* e = entry + SLAB_ENTRY_U32 * i;
                const uint32_t c = MIN(e[0], (uint32_t)SLAB_MAX_CORNERS);
                for (uint32_t k = 0; k < c; ++k) b->atoms[SLAB_MAX_CORNERS * i + k] = e[1 + k];
            }
        }
        break;
    default:
        // The control points keep their atoms
        ASSERT(false);
        break;
    }
    md_temp_end(temp);
    return b;
}

// The lists of the kinds (SUBSET_* bits) found again where the subset or the primitives changed since. Returns whether
// any was.
static bool subset_refresh(subset_t* s, molecule_t* mol, uint32_t kinds) {
    uint32_t stale = 0;
    for (int k = 0; k < SUBSET_KIND_COUNT; ++k) {
        if ((kinds & (1u << k)) && s->list[k].version != mol_kind_blocks(mol, k)->version) stale |= 1u << k;
    }
    if (!stale) return false;
    const uint64_t* bits = subset_bits_set(mol, s);
    for (int k = 0; k < SUBSET_KIND_COUNT; ++k) {
        if (!(stale & (1u << k))) continue;
        prim_blocks_t* b = (prim_blocks_t*)mol_kind_atoms(mol, k);
        subset_find(s->list + k, s, b, bits, mol->atom_count);
        s->list[k].version = b->version;
    }
    subset_bits_clear(mol, s);
    return true;
}

// The box of the atoms of a subset in positions xyz of the atoms [offset, offset + count) of its molecule, merged into
// mn, mx; whether it has any there
static bool subset_aabb_of(const subset_t* s, uint32_t offset, uint32_t count, const vec3_t* xyz, vec3_t* mn, vec3_t* mx) {
    // The first of its atoms at or past offset: they are in order
    uint32_t lo = 0, hi = s->atom_count;
    while (lo < hi) { const uint32_t mid = (lo + hi) / 2; if (s->atoms[mid] < offset) lo = mid + 1; else hi = mid; }
    bool any = false;
    for (uint32_t i = lo; i < s->atom_count && s->atoms[i] < offset + count; ++i) {
        const vec3_t p = xyz[s->atoms[i] - offset];
        *mn = vec3_min(*mn, p);
        *mx = vec3_max(*mx, p);
        any = true;
    }
    return any;
}

static void subsets_aabb_update(const molecule_t* mol, uint32_t offset, uint32_t count, const vec3_t* xyz) {
    const bool all = offset == 0 && count == mol->atom_count;
    for (int i = 1; i < ctx.subset_pool.size; ++i) {
        subset_t* s = ctx.subsets + i;
        if (s->id == MD_HANDLE_INVALID_ID || s->mol_id != mol->id) continue;
        if (all) {
            s->aabb_min = (vec3_t){ FLT_MAX,  FLT_MAX,  FLT_MAX};
            s->aabb_max = (vec3_t){-FLT_MAX, -FLT_MAX, -FLT_MAX};
            subset_aabb_of(s, 0, count, xyz, &s->aabb_min, &s->aabb_max);
            s->aabb_valid = true;
        } else if (s->aabb_valid) {
            subset_aabb_of(s, offset, count, xyz, &s->aabb_min, &s->aabb_max);
        }
    }
}

static void subsets_aabb_invalidate(const molecule_t* mol) {
    for (int i = 1; i < ctx.subset_pool.size; ++i) {
        subset_t* s = ctx.subsets + i;
        if (s->id != MD_HANDLE_INVALID_ID && s->mol_id == mol->id) s->aabb_valid = false;
    }
}

// The box of a subset for the culling of the copies of an op, found from the positions on the GPU if it is not known
// (read back once per draw and molecule, into *xyz). False when there is none to cull with: the box of the molecule then.
static bool subset_aabb(subset_t* s, const molecule_t* mol, const vec3_t** xyz, md_allocator_i* alloc) {
    if (s->aabb_valid) return true;
    if (!mol->aabb_valid) return false;
    if (!*xyz) {
        vec3_t* pos = md_alloc(alloc, MAX(mol->atom_count, 1) * sizeof(vec3_t));
        STATIC_ASSERT(sizeof(vec3_t) == 3 * sizeof(uint32_t), "The positions are packed xyz");
        buffer_read_u32((uint32_t*)pos, mol->buffer[GL_BUFFER_ATOM_POSITION], 3 * (size_t)mol->atom_count);
        *xyz = pos;
    }
    s->aabb_min = (vec3_t){ FLT_MAX,  FLT_MAX,  FLT_MAX};
    s->aabb_max = (vec3_t){-FLT_MAX, -FLT_MAX, -FLT_MAX};
    subset_aabb_of(s, 0, mol->atom_count, *xyz, &s->aabb_min, &s->aabb_max);
    s->aabb_valid = true;
    return true;
}

static bool compute_spline(molecule_t* mol);

// A range of the list buffer of a draw (ctx.list_buf), in u32
typedef struct {
    uint32_t offset;
    uint32_t count;
} list_range_t;

// What a draw covers: the copies of its op in the instance buffer of the draw (INSTANCES), and with an atom mask
// or a subset (masked) a list: with a mask what holds marked atoms (MARKS), the blocks of a kind of quads; of a subset
// (exact) the primitives within it (SUBSETS); the atoms of the spheres either way
typedef struct {
    uint32_t     inst_base;
    uint32_t     inst_count;
    bool         masked;
    bool         exact;
    list_range_t list;
    GLuint       buf;           // The buffer the list is in: a subset's own, 0 for the lists of the draw (ctx.list_buf)
} draw_range_t;

typedef struct {
    gl_buffer_t atoms;      // u32[2] per bond
    gl_buffer_t weight;     // f32 per bond, weak bonds only
    gl_buffer_t order;      // u32 per bond, covalent bonds with bond orders only (GL_BUFFER_BOND_ORDER)
    gl_buffer_t order_ref;  // u32[4] per entry (GL_BUFFER_BOND_ORDER_REF)
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
    float    dash_length;   // Dashed lines (weak bonds, partners, gaps): the length of the dashes (Ångström), see DASHED LINES
    float    gap_length;    // And of the gaps between them
    uint32_t bond_orders;   // Show the orders of covalent bonds (bond_set_t order, order_ref)
    float    spread;        // Envelope of the thinner sticks of multiple bonds, relative to the radius
    uint32_t gap_index_base;// The gaps of the chains: their picking base (md_gl_picking_offset_t gap_base)
    uint32_t bond_kinds;    // The kind of each bond is known (bond_set_t order, which tells the elastic bonds)
    uint32_t elastic_bonds; // Show the elastic bonds (BOND_KIND_ELASTIC), which are left out otherwise
} bond_params_t;

// The parameters of the gaps of the chains (draw_chain_gaps): a line of the radius with its dashes in proportion,
// colored as the bonds are, drawn at any length
static inline bond_params_t chain_gap_params(float radius, int32_t mode, float sharpness, uint32_t uniform_color, uint32_t gap_index_base) {
    bond_params_t params = {
        .radius = radius,
        .max_d2 = 3.4e38f,
        .mode = mode,
        .sharpness = sharpness,
        .uniform_color = uniform_color,
        .dash_length = CHAIN_GAP_DASH_LENGTH * radius,
        .gap_length = CHAIN_GAP_GAP_LENGTH * radius,
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

// ### INSTANCES OF A DRAW ###

enum {
    BACKBONE_CHAIN_CULLED = 0xFF
};

// Of the backbone representations, how far the spline may reach beyond its control points (Ångström): a spline through
// them (BACKBONE_FLAG_THROUGH) overshoots a little where it bends
#define BACKBONE_SPLINE_SLACK 1.0f

// The instances of a backbone draw with the same levels of detail of the chains, which are drawn together, and whether
// they are mirror images: their meshes face the other way
typedef struct {
    uint32_t       inst_base;
    uint32_t       inst_count;
    const uint8_t* chain_lod;   // [chain_count] Level of detail of each chain, or BACKBONE_CHAIN_CULLED
    uint64_t       key;         // Of chain_lod
    bool           mirrored;
} bb_group_t;

// The lists of a molecule in a draw with an atom mask (MARKS), or of a subset (SUBSETS), pushed to the list buffer as
// an op needs them
typedef struct {
    molecule_t*     mol;
    const marked_t* marked;     // With an atom mask
    subset_t*       subset;     // Or of a subset: its lists are exact
    uint32_t        pushed;     // LIST_* bits of those pushed
    list_range_t    atoms;      // Marked atoms, or of the subset
    list_range_t    cps;        // Control points of those atoms, in order
    list_range_t    bonds, partners, rings, gap_trace, gap_link, slabs;    // Blocks with marked atoms, or the primitives of the subset
    const uint32_t* cp_items;   // The control points of cps, on the CPU
    uint32_t        cp_item_count;
    GLuint          buf;        // The buffer the lists are in: of a subset its own, otherwise 0 for those of the draw
} mol_lists_t;

enum {
    LIST_ATOMS     = 1,
    LIST_CPS       = 2,
    LIST_BONDS     = 4,
    LIST_PARTNERS  = 8,
    LIST_RINGS     = 16,
    LIST_GAP_TRACE = 32,
    LIST_GAP_LINK  = 64,
    LIST_SLABS     = 128,
};

// What md_gl_draw found for an op before drawing it
typedef struct {
    uint32_t     inst_base;     // Its copies in view, in the instance buffer
    uint32_t     inst_count;
    uint32_t     group_count;   // Backbone: the copies in groups
    bb_group_t*  groups;
    mol_lists_t* lists;         // With an atom mask, the lists of its molecule
    list_range_t weak_blocks;   // With an atom mask, the blocks of the weak bonds of its representation, of a subset those within it
} op_prep_t;

// The planes of the view volume of a matrix to clip space, inward and of unit normal (xyz), and what the level of detail
// of the backbone takes from it: the row of w, the length of its xyz and the pixels per unit of length at w = 1
typedef struct {
    vec4_t plane[6];
    vec4_t row_w;
    float  w_scale;
    float  pixel_scale;
} view_frame_t;

static view_frame_t view_frame_make(const mat4_t* to_clip, float viewport_half_height) {
    vec4_t row[4];
    for (int r = 0; r < 4; ++r) {
        row[r] = (vec4_t){to_clip->elem[0][r], to_clip->elem[1][r], to_clip->elem[2][r], to_clip->elem[3][r]};
    }
    view_frame_t f = {
        .plane = {
            vec4_add(row[3], row[0]), vec4_sub(row[3], row[0]),
            vec4_add(row[3], row[1]), vec4_sub(row[3], row[1]),
            vec4_add(row[3], row[2]), vec4_sub(row[3], row[2]),
        },
        .row_w = row[3],
        .w_scale = vec3_length(vec3_from_vec4(row[3])),
        .pixel_scale = vec3_length(vec3_from_vec4(row[1])) * viewport_half_height,
    };
    for (int i = 0; i < 6; ++i) {
        const float len = vec3_length(vec3_from_vec4(f.plane[i]));
        f.plane[i] = len > 0.0f ? vec4_mul1(f.plane[i], 1.0f / len) : (vec4_t){0, 0, 0, 1};
    }
    return f;
}

// The coarsest level of detail whose segments are at most BACKBONE_LOD_SEGMENT_PIXELS long at the depth w (clip space)
static uint8_t backbone_lod_at(const view_frame_t* f, const backbone_mesh_t meshes[BACKBONE_LOD_COUNT], float w) {
    if (w <= 0.0f || f->pixel_scale <= 0.0f) {
        return 0;
    }
    const float pixels_per_unit = f->pixel_scale / w;
    for (uint32_t lod = BACKBONE_LOD_COUNT - 1; lod > 0; --lod) {
        if (BACKBONE_RESIDUE_LENGTH / (float)meshes[lod].segments * pixels_per_unit <= BACKBONE_LOD_SEGMENT_PIXELS) {
            return (uint8_t)lod;
        }
    }
    return 0;
}

// Level of detail of each chain, or BACKBONE_CHAIN_CULLED: the chains are culled against the view frustum and get the
// coarsest level whose segments are at most BACKBONE_LOD_SEGMENT_PIXELS long where the chain is closest to the camera.
static void backbone_chain_lod(uint8_t* chain_lod, const backbone_bounds_t* bounds, const backbone_mesh_t meshes[BACKBONE_LOD_COUNT], float max_extent, const view_frame_t* f) {
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
            if (vec4_dot(f->plane[i], center) < -radius) {
                culled = true;
                break;
            }
        }
        if (culled) {
            chain_lod[c] = BACKBONE_CHAIN_CULLED;
            continue;
        }
        chain_lod[c] = backbone_lod_at(f, meshes, vec4_dot(f->row_w, center) - radius * f->w_scale);
    }
}

// The level of detail all chains get when they are all in view and all at the same level, which the sphere of all
// chains tells without looking at each (no chain is closer or farther than it): -1 when it does not
static int backbone_uniform_lod(const backbone_bounds_t* bounds, const backbone_mesh_t meshes[BACKBONE_LOD_COUNT], float max_extent, const view_frame_t* f) {
    if (!bounds->valid) return -1;
    const vec4_t s = bounds->all_sphere;
    const vec4_t center = {s.x, s.y, s.z, 1.0f};
    const float radius = s.w + max_extent;
    for (int i = 0; i < 6; ++i) {
        if (vec4_dot(f->plane[i], center) < radius) return -1;
    }
    const float w = vec4_dot(f->row_w, center);
    const uint8_t lod_near = backbone_lod_at(f, meshes, w - radius * f->w_scale);
    const uint8_t lod_far  = backbone_lod_at(f, meshes, w + radius * f->w_scale);
    return lod_near == lod_far ? lod_near : -1;
}

// Whether the box [lo, hi] grown by margin, moved by the rigid transform T, is out of view
static bool box_out_of_view(const view_frame_t* f, const mat4_t* T, vec3_t lo, vec3_t hi, float margin) {
    const vec3_t c_local = vec3_mul1(vec3_add(lo, hi), 0.5f);
    const vec3_t e = vec3_mul1(vec3_sub(hi, lo), 0.5f);
    const vec4_t c = mat4_mul_vec4(*T, vec4_from_vec3(c_local, 1.0f));
    const vec3_t ax = vec3_from_vec4(T->col[0]);
    const vec3_t ay = vec3_from_vec4(T->col[1]);
    const vec3_t az = vec3_from_vec4(T->col[2]);
    for (int i = 0; i < 6; ++i) {
        const vec3_t n = vec3_from_vec4(f->plane[i]);
        const float dist = vec3_dot(n, vec3_from_vec4(c)) + f->plane[i].w;
        const float reach = fabsf(vec3_dot(n, ax)) * e.x + fabsf(vec3_dot(n, ay)) * e.y + fabsf(vec3_dot(n, az)) * e.z + margin;
        if (dist < -reach) return true;
    }
    return false;
}

static inline mat4_t instance_transform(const float* transforms, uint32_t i) {
    mat4_t T;
    MEMCPY(&T, transforms + 16 * i, sizeof(mat4_t));
    return T;
}

static inline bool is_backbone_representation_type(md_gl_rep_type_t type) {
    return type == MD_GL_REP_RIBBONS || type == MD_GL_REP_CARTOON || type == MD_GL_REP_TUBES;
}

// The backbone of an op: its profile, the scales of the profile (u_scale of backbone.vert) and the largest half width
// of the profiles, how far the surface reaches from the spline
static void backbone_params(uint32_t* profile, float scale[4], float* max_extent, const md_gl_draw_op_t* op, const molecule_t* mol) {
    switch (op->type) {
    case MD_GL_REP_RIBBONS:
        // Half width and half thickness of the box profile
        *profile = BACKBONE_PROFILE_BOX;
        scale[0] = op->args.ribbons.width_scale;
        scale[1] = op->args.ribbons.thickness_scale * 0.1f;
        scale[2] = scale[3] = 0.0f;
        *max_extent = sqrtf(scale[0] * scale[0] + scale[1] * scale[1]);
        break;
    case MD_GL_REP_CARTOON:
        // Scales of the coil, helix and sheet profiles, and the largest semi axis of the profiles in backbone.vert
        *profile = BACKBONE_PROFILE_STADIUM;
        scale[0] = op->args.cartoon.coil_scale;
        scale[1] = op->args.cartoon.helix_scale;
        scale[2] = op->args.cartoon.sheet_scale;
        scale[3] = 0.0f;
        *max_extent = MAX(MAX(0.2f * scale[0], 1.2f * scale[1]), 1.5f * scale[2]);
        break;
    case MD_GL_REP_TUBES: {
        // The (l, r) of each control point times the scale, the r of a band (l > 0) times the band scale, and round caps
        const float band = op->args.tubes.band_scale > 0.0f ? op->args.tubes.band_scale : 1.0f;
        *profile = BACKBONE_PROFILE_TUBE;
        scale[0] = op->args.tubes.radius_scale;
        scale[1] = band;
        scale[2] = scale[3] = 0.0f;
        *max_extent = mol->tube_max_extent * scale[0] * MAX(band, 1.0f);
        break;
    }
    default:
        ASSERT(false);
        break;
    }
}

// How far what an op draws reaches beyond the atoms of its molecule (Ångström), which its copies are culled with
static float op_reach(const md_gl_draw_op_t* op, const molecule_t* mol, const representation_t* rep) {
    // The sticks of licorice and ball and stick at radius r: their multiple bonds within u_spread (1.3) r, the dashed
    // partners of delocalized bonds offset by less than 1.9 * 1.3 r and a third of r thick (bond.vert)
    const float stick_reach = 3.0f;
    float reach = 0.0f;
    switch (op->type) {
    case MD_GL_REP_SPACE_FILL:
        reach = mol->max_radius * op->args.space_fill.radius_scale;
        break;
    case MD_GL_REP_LICORICE:
        reach = stick_reach * 0.2f * op->args.licorice.radius;
        break;
    case MD_GL_REP_BALL_AND_STICK:
        reach = MAX(mol->max_radius * 0.2f * op->args.ball_and_stick.ball_scale, stick_reach * 0.2f * op->args.ball_and_stick.stick_radius);
        break;
    case MD_GL_REP_WEAK_BONDS:
        reach = op->args.weak_bonds.radius * MAX(rep->weak_bond_max_weight, 1.0f);
        break;
    case MD_GL_REP_RIBBONS:
    case MD_GL_REP_CARTOON:
    case MD_GL_REP_TUBES: {
        uint32_t profile;
        float scale[4];
        float max_extent;
        backbone_params(&profile, scale, &max_extent, op, mol);
        reach = max_extent + BACKBONE_SPLINE_SLACK;
        break;
    }
    case MD_GL_REP_SLABS:
        reach = op->args.slabs.radius;
        break;
    default:
        break;
    }
    return fabsf(reach) * 1.01f + 1.0e-2f;
}

static uint64_t hash_bytes(const uint8_t* data, size_t size) {
    uint64_t h = 1469598103934665603ull;
    for (size_t i = 0; i < size; ++i) {
        h = (h ^ data[i]) * 1099511628211ull;
    }
    return h;
}

// The rows of a transform and of its previous one, as instance_load reads them (INSTANCES)
static void push_instance(md_array(vec4_t)* data, const mat4_t* T, const mat4_t* T_prev, md_allocator_i* alloc) {
    for (int k = 0; k < 2; ++k) {
        const mat4_t* M = k == 0 ? T : T_prev;
        for (int r = 0; r < 3; ++r) {
            const vec4_t row = {M->elem[0][r], M->elem[1][r], M->elem[2][r], M->elem[3][r]};
            md_array_push(*data, row, alloc);
        }
    }
}

static list_range_t push_list(md_array(uint32_t)* data, const uint32_t* items, uint32_t count, md_allocator_i* alloc) {
    const list_range_t range = {(uint32_t)md_array_size(*data), count};
    md_array_push_array(*data, items, count, alloc);
    return range;
}

// The blocks which hold a marked atom
static list_range_t push_blocks(md_array(uint32_t)* data, const prim_blocks_t* blocks, const marked_t* marked, md_allocator_i* alloc) {
    list_range_t range = {(uint32_t)md_array_size(*data), 0};
    for (uint32_t b = 0; b < blocks->count; ++b) {
        if (marked_any(marked, blocks->range[2 * b + 0], blocks->range[2 * b + 1])) {
            md_array_push(*data, b, alloc);
            range.count += 1;
        }
    }
    return range;
}

// The kind of primitives of a subset each list is of (LIST_* bit i is SUBSET_* list_subset_kind[i]), -1 for the atoms
static const int list_subset_kind[] = { -1, SUBSET_CPS, SUBSET_BONDS, SUBSET_PARTNERS, SUBSET_RINGS, SUBSET_GAP_TRACE, SUBSET_GAP_LINK, SUBSET_SLABS };

// The lists of a subset which an op needs (LIST_* bits), pushed once per draw
// The lists of a subset, those of lists (LIST_* bits) which are found, uploaded to its buffer on the GPU one after another
static void subset_upload(subset_t* s, const molecule_t* mol, uint32_t lists) {
    uint32_t count[SUBSET_GPU_LISTS] = {0};
    uint32_t total = 0;
    uint32_t uploaded = 0;
    for (int i = 0; i < SUBSET_GPU_LISTS; ++i) {
        if (!(lists & (1u << i))) continue;
        const int k = list_subset_kind[i];
        if (k >= 0 && s->list[k].version != mol_kind_blocks((molecule_t*)mol, k)->version) continue;   // Not found: not drawn
        count[i] = k < 0 ? s->atom_count : s->list[k].count;
        total += count[i];
        uploaded |= 1u << i;
    }
    md_temp_scope_t temp = md_temp_begin();
    uint32_t* data = md_temp_alloc(temp, MAX(total, 1) * sizeof(uint32_t));
    uint32_t at = 0;
    for (int i = 0; i < SUBSET_GPU_LISTS; ++i) {
        s->gpu_offset[i] = at;
        s->gpu_count[i] = count[i];
        if (!(uploaded & (1u << i))) continue;
        const int k = list_subset_kind[i];
        if (k < 0) {
            MEMCPY(data + at, s->atoms, count[i] * sizeof(uint32_t));
        } else if (k == SUBSET_CPS) {
            // Each with whether its neighbours are out of the subset, where the backbone is closed (backbone.vert)
            const subset_list_t* cps = s->list + SUBSET_CPS;
            for (uint32_t j = 0; j < cps->count; ++j) {
                const uint32_t c = cps->prims[j];
                const bool prev_in = j > 0 && cps->prims[j - 1] + 1 == c;
                const bool next_in = j + 1 < cps->count && cps->prims[j + 1] == c + 1;
                data[at + j] = c | (prev_in ? 0 : CP_PREV_OUT) | (next_in ? 0 : CP_NEXT_OUT);
            }
        } else {
            MEMCPY(data + at, s->list[k].prims, count[i] * sizeof(uint32_t));
        }
        at += count[i];
    }
    const uint32_t bytes = MAX(total, 1) * sizeof(uint32_t);
    if (!s->gpu.id) {
        s->gpu = gl_buffer_create(bytes, data, GL_STATIC_DRAW);
    } else {
        gl_buffer_set_data(s->gpu, bytes, data);
    }
    s->gpu_lists = uploaded;
    s->gpu_stale = false;
    md_temp_end(temp);
}

// The lists of a subset which an op needs (LIST_* bits): found again where the subset or the primitives changed, and
// then uploaded again, with those it had, to the buffer of the subset, which the draws of the op read (mol_lists_t buf).
// A subset which did not change uploads nothing.
static void subset_lists_need(mol_lists_t* l, uint32_t need) {
    const uint32_t todo = need & ~l->pushed;
    subset_t* s = l->subset;
    uint32_t kinds = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(list_subset_kind); ++i) {
        if ((todo & (1u << i)) && list_subset_kind[i] >= 0) kinds |= 1u << list_subset_kind[i];
    }
    const bool found = subset_refresh(s, l->mol, kinds);
    if (found || s->gpu_stale || !s->gpu.id || (todo & ~s->gpu_lists)) {
        subset_upload(s, l->mol, s->gpu_lists | todo);
    }
    list_range_t* range[SUBSET_GPU_LISTS] = {&l->atoms, &l->cps, &l->bonds, &l->partners, &l->rings, &l->gap_trace, &l->gap_link, &l->slabs};
    for (int i = 0; i < SUBSET_GPU_LISTS; ++i) {
        if (todo & (1u << i)) *range[i] = (list_range_t){s->gpu_offset[i], s->gpu_count[i]};
    }
    if (todo & LIST_CPS) {
        l->cp_items = s->list[SUBSET_CPS].prims;
        l->cp_item_count = s->list[SUBSET_CPS].count;
    }
    l->buf = s->gpu.id;
    l->pushed |= need;
}

// The weak bonds of a representation within a subset: found in each draw, as they change every frame
static list_range_t push_weak_subset(md_array(uint32_t)* data, const prim_blocks_t* weak, molecule_t* mol, const subset_t* s, md_allocator_i* alloc) {
    list_range_t range = {(uint32_t)md_array_size(*data), 0};
    if (weak->prim_count == 0 || s->atom_count == 0) return range;
    const uint64_t* bits = subset_bits_set(mol, s);
    for (uint32_t p = 0; p < weak->prim_count; ++p) {
        uint32_t first;
        if (prim_within(weak, p, bits, mol->atom_count, &first)) {
            md_array_push(*data, p, alloc);
            range.count += 1;
        }
    }
    subset_bits_clear(mol, s);
    return range;
}

// The lists of a molecule which an op needs (LIST_* bits), pushed once per draw
static void mol_lists_need(mol_lists_t* l, uint32_t need, md_array(uint32_t)* data, md_allocator_i* alloc) {
    if (l->subset) {
        subset_lists_need(l, need);
        return;
    }
    const uint32_t todo = need & ~l->pushed;
    const molecule_t* mol = l->mol;
    const marked_t* m = l->marked;
    if (todo & LIST_ATOMS)     l->atoms     = push_list(data, m->atoms, m->atom_count, alloc);
    if (todo & LIST_CPS) {
        l->cps = push_list(data, m->cps, m->cp_count, alloc);
        l->cp_items = m->cps;
        l->cp_item_count = m->cp_count;
    }
    if (todo & LIST_BONDS)     l->bonds     = push_blocks(data, &mol->bond_blocks, m, alloc);
    if (todo & LIST_PARTNERS)  l->partners  = push_blocks(data, &mol->partner_blocks, m, alloc);
    if (todo & LIST_RINGS)     l->rings     = push_blocks(data, &mol->ring_blocks, m, alloc);
    if (todo & LIST_GAP_TRACE) l->gap_trace = push_blocks(data, &mol->gap_trace_blocks, m, alloc);
    if (todo & LIST_GAP_LINK)  l->gap_link  = push_blocks(data, &mol->gap_link_blocks, m, alloc);
    if (todo & LIST_SLABS)     l->slabs     = push_blocks(data, &mol->slab_blocks, m, alloc);
    l->pushed |= need;
}

// The lists an op of a type draws from
static uint32_t op_lists(const md_gl_draw_op_t* op) {
    switch (op->type) {
    case MD_GL_REP_SPACE_FILL:      return LIST_ATOMS;
    case MD_GL_REP_LICORICE:        return LIST_BONDS | LIST_PARTNERS | LIST_RINGS | LIST_GAP_LINK;
    case MD_GL_REP_BALL_AND_STICK:  return LIST_ATOMS | LIST_BONDS | LIST_PARTNERS | LIST_RINGS | LIST_GAP_LINK;
    case MD_GL_REP_RIBBONS:
    case MD_GL_REP_CARTOON:         return LIST_CPS | LIST_GAP_TRACE;
    case MD_GL_REP_TUBES:           return LIST_CPS;
    case MD_GL_REP_SLABS:           return LIST_SLABS;
    default:                        return 0;
    }
}

// buf: the buffer of the list, 0 for the lists of the draw
static inline draw_range_t draw_range(const op_prep_t* prep, const list_range_t* list, GLuint buf) {
    draw_range_t r = {
        .inst_base  = prep->inst_base,
        .inst_count = prep->inst_count,
        .masked     = list != NULL,
        .exact      = list != NULL && prep->lists && prep->lists->subset != NULL,
        .buf        = buf,
    };
    if (list) r.list = *list;
    return r;
}

// The list texture (LIST_TEXTURE_UNIT) on the buffer of a list: a subset's own, or (0) the lists of the draw
static void bind_list(GLuint buf) {
    if (!buf) buf = ctx.list_buf.id;
    if (buf == ctx.list_bound) return;
    glActiveTexture(GL_TEXTURE0 + LIST_TEXTURE_UNIT);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, buf);
    glActiveTexture(GL_TEXTURE0);
    ctx.list_bound = buf;
}

typedef struct {
    const view_frame_t* frame;          // Of world space
    const mat4_t*       world_to_clip;
    float               viewport_half_height;
    md_array(vec4_t)*   inst_data;
    md_array(uint32_t)* list_data;
    md_array(mol_lists_t)* mol_lists;   // With an atom mask or subsets
    uint32_t            atom_mask;
    md_allocator_i*     alloc;
    const molecule_t*   xyz_mol;        // The positions of this molecule read back for the boxes of its subsets, if any
    const vec3_t*       xyz;
} prep_ctx_t;

// The copies of an op in view, pushed to the instance data (in groups of the levels of detail of the chains for the
// backbone), and with an atom mask the lists it draws from
static void prepare_op(op_prep_t* prep, const md_gl_draw_op_t* op, molecule_t* mol, const representation_t* rep, prep_ctx_t* pc) {
    MEMSET(prep, 0, sizeof(op_prep_t));
    md_allocator_i* alloc = pc->alloc;
    const md_gl_copies_t* inst = &op->copies;
    const uint32_t count = (inst->count && inst->transforms) ? inst->count : 1;
    const float reach = op_reach(op, mol, rep);

    // The box of what the op shows, which its copies are culled with: of its subset, or of the molecule
    subset_t* subset = op->subset.id ? subset_lookup(op->subset.id) : NULL;
    bool   culled  = mol->aabb_valid;
    vec3_t box_min = mol->aabb_min;
    vec3_t box_max = mol->aabb_max;
    if (subset) {
        if (subset->atom_count == 0) return;     // Nothing to show
        if (pc->xyz_mol != mol) {
            pc->xyz_mol = mol;
            pc->xyz = NULL;
        }
        if (subset_aabb(subset, mol, &pc->xyz, alloc)) {
            culled  = true;
            box_min = subset->aabb_min;
            box_max = subset->aabb_max;
        }
    }

    // In view
    uint32_t* visible = md_alloc(alloc, count * sizeof(uint32_t));
    uint32_t visible_count = 0;
    for (uint32_t i = 0; i < count; ++i) {
        const mat4_t T = (inst->count && inst->transforms) ? instance_transform(inst->transforms, i) : mat4_ident();
        if (culled && box_out_of_view(pc->frame, &T, box_min, box_max, reach)) {
            continue;
        }
        visible[visible_count++] = i;
    }

    // The order the copies go in: the backbone in groups of the same levels of detail and facing
    uint32_t* order = visible;
    if (is_backbone_representation_type(op->type) && visible_count > 0) {
        const backbone_bounds_t* bounds = &mol->backbone_bounds;
        uint32_t profile;
        float scale[4];
        float max_extent;
        backbone_params(&profile, scale, &max_extent, op, mol);
        const backbone_mesh_t* meshes = ctx.backbone_mesh[profile];
        const uint32_t chain_count = bounds->chain_count;

        // The same level of detail for all chains, shared by the copies which have it
        uint8_t* uniform[BACKBONE_LOD_COUNT] = {0};
        bb_group_t* groups = md_alloc(alloc, visible_count * sizeof(bb_group_t));
        uint32_t*   member = md_alloc(alloc, visible_count * sizeof(uint32_t));   // Group of each visible copy
        uint32_t group_count = 0;

        for (uint32_t v = 0; v < visible_count; ++v) {
            const uint32_t i = visible[v];
            const mat4_t T = (inst->count && inst->transforms) ? instance_transform(inst->transforms, i) : mat4_ident();
            const mat4_t to_clip = mat4_mul(*pc->world_to_clip, T);
            const view_frame_t f = view_frame_make(&to_clip, pc->viewport_half_height);

            const uint8_t* lod = NULL;
            const int u = backbone_uniform_lod(bounds, meshes, max_extent, &f);
            if (u >= 0) {
                if (!uniform[u]) {
                    uint8_t* arr = md_alloc(alloc, MAX(chain_count, 1));
                    for (uint32_t c = 0; c < chain_count; ++c) {
                        arr[c] = (bounds->chain_offset[c + 1] - bounds->chain_offset[c] < 2) ? BACKBONE_CHAIN_CULLED : (uint8_t)u;
                    }
                    uniform[u] = arr;
                }
                lod = uniform[u];
            } else {
                uint8_t* arr = md_alloc(alloc, MAX(chain_count, 1));
                backbone_chain_lod(arr, bounds, meshes, max_extent, &f);
                lod = arr;
            }
            const uint64_t key = hash_bytes(lod, chain_count);
            const bool mirrored = mat3_determinant(mat3_from_mat4(T)) < 0.0f;

            uint32_t g = 0;
            for (; g < group_count; ++g) {
                if (groups[g].key == key && groups[g].mirrored == mirrored && (groups[g].chain_lod == lod || MEMCMP(groups[g].chain_lod, lod, chain_count) == 0)) break;
            }
            if (g == group_count) {
                groups[g] = (bb_group_t){.chain_lod = lod, .key = key, .mirrored = mirrored};
                group_count += 1;
            }
            groups[g].inst_count += 1;
            member[v] = g;
        }

        // The copies by group
        uint32_t base = 0;
        for (uint32_t g = 0; g < group_count; ++g) {
            groups[g].inst_base = base;
            base += groups[g].inst_count;
            groups[g].inst_count = 0;
        }
        order = md_alloc(alloc, visible_count * sizeof(uint32_t));
        for (uint32_t v = 0; v < visible_count; ++v) {
            bb_group_t* grp = groups + member[v];
            order[grp->inst_base + grp->inst_count++] = visible[v];
        }
        prep->groups = groups;
        prep->group_count = group_count;
    }

    prep->inst_base  = (uint32_t)(md_array_size(*pc->inst_data) / INSTANCE_TEXELS);
    prep->inst_count = visible_count;
    for (uint32_t v = 0; v < visible_count; ++v) {
        const uint32_t i = order[v];
        if (inst->count && inst->transforms) {
            const mat4_t T = instance_transform(inst->transforms, i);
            const mat4_t T_prev = inst->prev_transforms ? instance_transform(inst->prev_transforms, i) : T;
            push_instance(pc->inst_data, &T, &T_prev, alloc);
        } else {
            const mat4_t I = mat4_ident();
            push_instance(pc->inst_data, &I, &I, alloc);
        }
    }
    for (uint32_t g = 0; g < prep->group_count; ++g) {
        prep->groups[g].inst_base += prep->inst_base;
    }

    // MARKS, SUBSETS: of a subset its lists, which the atom mask applies within (in the shaders); otherwise with an atom
    // mask those of the marked atoms
    if ((pc->atom_mask || subset) && visible_count > 0) {
        mol_lists_t* lists = NULL;
        for (size_t k = 0; k < md_array_size(*pc->mol_lists); ++k) {
            mol_lists_t* l = *pc->mol_lists + k;
            if (subset ? l->subset == subset : (l->mol == mol && !l->subset)) {
                lists = l;
                break;
            }
        }
        if (!lists) {
            const mol_lists_t entry = subset ? (mol_lists_t){.mol = mol, .subset = subset}
                                             : (mol_lists_t){.mol = mol, .marked = mol_marked(mol, pc->atom_mask)};
            md_array_push(*pc->mol_lists, entry, alloc);
            lists = md_array_last(*pc->mol_lists);
        }
        mol_lists_need(lists, op_lists(op), pc->list_data, alloc);
        prep->lists = lists;
        if (op->type == MD_GL_REP_WEAK_BONDS) {
            prep->weak_blocks = subset ? push_weak_subset(pc->list_data, &rep->weak_bond_blocks, mol, subset, alloc)
                                       : push_blocks(pc->list_data, &rep->weak_bond_blocks, lists->marked, alloc);
        }
    }
}

static bool draw_space_fill(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, float scale, const draw_range_t* range);
static bool draw_bonds     (gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, bond_set_t bonds, const bond_params_t* params, const float viewport_size[2], const draw_range_t* range);
static bool draw_bond_partners(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, const bond_params_t* bond_params, const float viewport_size[2], const draw_range_t* range);
static bool draw_chain_gaps(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, bool trace, const bond_params_t* params, const float viewport_size[2], const draw_range_t* range);
static bool draw_aromatic_rings(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, const bond_params_t* bond_params, uint32_t ring_index_base, const float viewport_size[2], const draw_range_t* range);
static bool draw_backbone  (gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, uint32_t profile, const float scale[4], const op_prep_t* prep, GLint front_face);
static bool draw_slabs     (gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, float radius, float plate, float max_length, uint32_t slab_index_base, const float viewport_size[2], const draw_range_t* range);

static inline void init_ubo_base_data(gl_ubo_base_t* ubo_data, const md_gl_draw_args_t* args) {
    ASSERT(ubo_data);
    ASSERT(args);

    MEMCPY(&ubo_data->view_transform.world_to_view, args->view_transform.view_matrix, sizeof(mat4_t));
    MEMCPY(&ubo_data->view_transform.view_to_clip,  args->view_transform.proj_matrix, sizeof(mat4_t));

    ubo_data->view_transform.world_to_clip = mat4_mul(ubo_data->view_transform.view_to_clip, ubo_data->view_transform.world_to_view);
    ubo_data->view_transform.world_to_view_normal = mat4_transpose(mat4_inverse(ubo_data->view_transform.world_to_view));
    ubo_data->view_transform.view_to_world = mat4_inverse(ubo_data->view_transform.world_to_view);
    ubo_data->view_transform.clip_to_view = mat4_inverse(ubo_data->view_transform.view_to_clip);

    if (args->view_transform.prev_view_matrix && args->view_transform.proj_matrix) {
        const mat4_t* prev_world_to_view = (const mat4_t*)args->view_transform.prev_view_matrix;
        const mat4_t* prev_view_to_clip  = (const mat4_t*)args->view_transform.prev_proj_matrix;
        // The copies are taken into world space by the vertex shaders (INSTANCES), with the motion of each in the velocity
        ubo_data->view_transform.prev_world_to_clip = mat4_mul(*prev_view_to_clip, *prev_world_to_view);
        ubo_data->view_transform.curr_view_to_prev_clip = mat4_mul(ubo_data->view_transform.prev_world_to_clip, ubo_data->view_transform.view_to_world);
        extract_jitter_uv(ubo_data->jitter_uv.elem + 0, ubo_data->view_transform.view_to_clip);
        extract_jitter_uv(ubo_data->jitter_uv.elem + 2, *prev_view_to_clip);
    }
    ubo_data->atom_mask = args->atom_mask;
    ubo_data->atom_index_base = args->picking_offset.atom_base;
    ubo_data->bond_index_base = args->picking_offset.bond_base;
    ubo_data->backbone_index_base = args->picking_offset.backbone_base;
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
    init_ubo_base_data(&ubo_data, args);

    ubo_write(0, sizeof(ubo_data), &ubo_data);

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    bool result = false;

    // Valid draw operations to issue
    md_array(md_gl_draw_op_t const*) draw_ops = 0;

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
        const bool bb_rep = is_backbone_representation_type(draw_op->type);
        if (bb_rep && !(mol->flags & MOL_FLAG_HAS_BACKBONE)) {
            MD_LOG_ERROR("Incompatible representation type for molecule: Molecule is missing a backbone (a protein backbone, or tubes)");
            continue;
        }
        if (draw_op->subset.id) {
            const subset_t* subset = subset_lookup(draw_op->subset.id);
            if (!subset) {
                MD_LOG_ERROR("Invalid subset");
                continue;
            }
            if (subset->mol_id != rep->mol_id) {
                MD_LOG_ERROR("The subset is of another molecule than the representation");
                continue;
            }
        }
        md_array_push(draw_ops, draw_op, alloc);
    }

    // The splines of the backbones drawn, once after what they are computed from changed: not once per draw, which
    // would also move their orientation history on (md_gl_mol_reset_backbone_history) more than once per frame
    PUSH_GPU_SECTION("COMPUTE SPLINE")
    for (size_t i = 0; i < md_array_size(draw_ops); i++) {
        if (!is_backbone_representation_type(draw_ops[i]->type)) continue;
        molecule_t* mol = mol_lookup(rep_lookup(draw_ops[i]->rep.id)->mol_id);
        if (mol->spline_dirty) {
            compute_spline(mol);
            mol->spline_dirty = false;
        }
    }
    POP_GPU_SECTION()

    int program_permutation = 0;
    if (is_ortho_proj_matrix(ubo_data.view_transform.view_to_clip)) program_permutation |= PERMUTATION_BIT_ORTHO;

    // Maximum bond length in units (Ångström assumed)
    const float max_length = args->max_bond_length > 0.0f ? args->max_bond_length : 5.0f;

    // For the level of detail of the backbone representations
    GLint viewport[4] = {0};
    glGetIntegerv(GL_VIEWPORT, viewport);
    const float viewport_half_height = (float)viewport[3] * 0.5f;
    const float viewport_size[2] = {(float)viewport[2], (float)viewport[3]};

    // INSTANCES and MARKS: what each op draws, culled on the CPU, uploaded once for all ops
    const size_t op_count = md_array_size(draw_ops);
    op_prep_t* preps = md_alloc(alloc, MAX(op_count, 1) * sizeof(op_prep_t));
    md_array(vec4_t)      inst_data = 0;
    md_array(uint32_t)    list_data = 0;
    md_array(mol_lists_t) mol_lists = 0;
    // Room for the lists of the molecules: their entries are pointed to and must not move
    md_array_ensure(mol_lists, op_count, alloc);
    {
        const view_frame_t frame = view_frame_make(&ubo_data.view_transform.world_to_clip, viewport_half_height);
        prep_ctx_t pc = {
            .frame = &frame,
            .world_to_clip = &ubo_data.view_transform.world_to_clip,
            .viewport_half_height = viewport_half_height,
            .inst_data = &inst_data,
            .list_data = &list_data,
            .mol_lists = &mol_lists,
            .atom_mask = args->atom_mask,
            .alloc = alloc,
        };
        for (size_t i = 0; i < op_count; ++i) {
            const representation_t* rep = rep_lookup(draw_ops[i]->rep.id);
            prepare_op(&preps[i], draw_ops[i], mol_lookup(rep->mol_id), rep, &pc);
        }
    }

    // Uploaded whole: a new store for each draw, which the driver need not wait for
    glBindBuffer(GL_ARRAY_BUFFER, ctx.instance_buf.id);
    glBufferData(GL_ARRAY_BUFFER, MAX(md_array_bytes(inst_data), sizeof(vec4_t)), inst_data, GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, ctx.list_buf.id);
    glBufferData(GL_ARRAY_BUFFER, MAX(md_array_bytes(list_data), sizeof(uint32_t)), list_data, GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glActiveTexture(GL_TEXTURE0 + INSTANCE_TEXTURE_UNIT);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_INSTANCE].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32F, ctx.instance_buf.id);
    glActiveTexture(GL_TEXTURE0 + LIST_TEXTURE_UNIT);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_LIST].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, ctx.list_buf.id);
    ctx.list_bound = ctx.list_buf.id;
    glActiveTexture(GL_TEXTURE0);

    // The meshes of mirrored copies face the other way (draw_backbone): the caller's front face, asked for only then,
    // as a query may wait for a driver of its own thread
    GLint front_face = GL_CCW;
    bool mirrored = false;
    for (size_t i = 0; i < op_count && !mirrored; ++i) {
        for (uint32_t g = 0; g < preps[i].group_count; ++g) {
            mirrored |= preps[i].groups[g].mirrored;
        }
    }
    if (mirrored) {
        glGetIntegerv(GL_FRONT_FACE, &front_face);
    }

    PUSH_GPU_SECTION("DRAW REPRESENTATIONS")
    for (size_t i = 0; i < op_count; i++) {
        const md_gl_draw_op_t* draw_op = draw_ops[i];
        const op_prep_t* prep = &preps[i];
        if (prep->inst_count == 0) continue;

        const representation_t* rep = rep_lookup(draw_op->rep.id);
        const molecule_t* mol       = mol_lookup(rep->mol_id);
        const mol_lists_t* lists    = prep->lists;
        // With an atom mask, what holds marked atoms
#define RANGE(member) draw_range(prep, lists ? &lists->member : NULL, lists ? lists->buf : 0)
        switch (draw_op->type) {
        case MD_GL_REP_SPACE_FILL: {
            const draw_range_t range = RANGE(atoms);
            draw_space_fill(shaders->spacefill[program_permutation], mol, rep->atom_color, draw_op->args.space_fill.radius_scale, &range);
            break;
        }
        case MD_GL_REP_LICORICE: {
            const bool kinds  = mol->bond_order_count == mol->bond_count;
            const bool orders = draw_op->args.licorice.bond_orders && kinds;
            const bond_params_t params = {
                .radius = 0.2f * draw_op->args.licorice.radius,
                .max_d2 = max_length * max_length,
                .mode = draw_op->args.licorice.color_mode,
                .sharpness = draw_op->args.licorice.sharpness,
                .uniform_color = draw_op->args.licorice.uniform_color,
                .bond_orders = orders,
                // Within the stick: their ends stay within the sphere at the atom
                .spread = 1.0f,
                .bond_kinds = kinds,
                .elastic_bonds = draw_op->args.licorice.elastic_bonds,
            };
            const bond_set_t bonds = {mol->buffer[GL_BUFFER_BOND_ATOM_INDICES], {0}, mol->buffer[GL_BUFFER_BOND_ORDER], mol->buffer[GL_BUFFER_BOND_ORDER_REF], {0}, mol->bond_count};
            const draw_range_t bond_range = RANGE(bonds);
            draw_bonds(shaders->licorice[program_permutation], mol, rep->atom_color, bonds, &params, viewport_size, &bond_range);
            if (orders) {
                const draw_range_t partner_range = RANGE(partners);
                const draw_range_t ring_range = RANGE(rings);
                draw_bond_partners(shaders->bond_partners[program_permutation], mol, rep->atom_color, &params, viewport_size, &partner_range);
                draw_aromatic_rings(shaders->aromatic_rings[program_permutation], mol, rep->atom_color, &params, args->picking_offset.ring_base, viewport_size, &ring_range);
            }
            if (draw_op->args.licorice.chain_gaps) {
                const bond_params_t gap = chain_gap_params(CHAIN_GAP_RADIUS_LINK * params.radius, params.mode, params.sharpness, params.uniform_color, args->picking_offset.gap_base);
                const draw_range_t gap_range = RANGE(gap_link);
                draw_chain_gaps(shaders->chain_gaps[program_permutation], mol, rep->atom_color, false, &gap, viewport_size, &gap_range);
            }
            break;
        }
        case MD_GL_REP_BALL_AND_STICK: {
            const bool kinds  = mol->bond_order_count == mol->bond_count;
            const bool orders = draw_op->args.ball_and_stick.bond_orders && kinds;
            const bond_params_t params = {
                .radius = 0.2f * draw_op->args.ball_and_stick.stick_radius,
                .max_d2 = max_length * max_length,
                .mode = draw_op->args.ball_and_stick.color_mode,
                .sharpness = draw_op->args.ball_and_stick.sharpness,
                .uniform_color = draw_op->args.ball_and_stick.uniform_color,
                .bond_orders = orders,
                // The balls are larger than the sticks
                .spread = 1.3f,
                .bond_kinds = kinds,
                .elastic_bonds = draw_op->args.ball_and_stick.elastic_bonds,
            };
            const bond_set_t bonds = {mol->buffer[GL_BUFFER_BOND_ATOM_INDICES], {0}, mol->buffer[GL_BUFFER_BOND_ORDER], mol->buffer[GL_BUFFER_BOND_ORDER_REF], {0}, mol->bond_count};
            const draw_range_t bond_range = RANGE(bonds);
            draw_bonds(shaders->licorice[program_permutation], mol, rep->atom_color, bonds, &params, viewport_size, &bond_range);
            if (orders) {
                const draw_range_t partner_range = RANGE(partners);
                const draw_range_t ring_range = RANGE(rings);
                draw_bond_partners(shaders->bond_partners[program_permutation], mol, rep->atom_color, &params, viewport_size, &partner_range);
                draw_aromatic_rings(shaders->aromatic_rings[program_permutation], mol, rep->atom_color, &params, args->picking_offset.ring_base, viewport_size, &ring_range);
            }
            if (draw_op->args.ball_and_stick.chain_gaps) {
                const bond_params_t gap = chain_gap_params(CHAIN_GAP_RADIUS_LINK * params.radius, params.mode, params.sharpness, params.uniform_color, args->picking_offset.gap_base);
                const draw_range_t gap_range = RANGE(gap_link);
                draw_chain_gaps(shaders->chain_gaps[program_permutation], mol, rep->atom_color, false, &gap, viewport_size, &gap_range);
            }
            const draw_range_t atom_range = RANGE(atoms);
            draw_space_fill(shaders->spacefill[program_permutation], mol, rep->atom_color, 0.2f * draw_op->args.ball_and_stick.ball_scale, &atom_range);
            break;
        }
        case MD_GL_REP_WEAK_BONDS: {
            const float dash = draw_op->args.weak_bonds.dash_length;
            const float gap  = draw_op->args.weak_bonds.gap_length > 0.0f ? draw_op->args.weak_bonds.gap_length : dash;
            if (!(dash > 0.0f)) {
                MD_LOG_ERROR("Weak bonds with a dash length of %f: it must be positive", dash);
                break;
            }
            const bond_params_t params = {
                .radius = draw_op->args.weak_bonds.radius,
                .max_d2 = max_length * max_length,
                .mode = draw_op->args.weak_bonds.color_mode,
                .sharpness = draw_op->args.weak_bonds.sharpness,
                .uniform_color = draw_op->args.weak_bonds.uniform_color,
                .dash_length = dash,
                .gap_length = gap,
            };
            const bond_set_t bonds = {rep->weak_bond_atoms, rep->weak_bond_weight, {0}, {0}, {0}, rep->weak_bond_count};
            // Found in this draw, in its lists
            const draw_range_t range = draw_range(prep, lists ? &prep->weak_blocks : NULL, 0);
            draw_bonds(shaders->weak_bonds[program_permutation], mol, rep->atom_color, bonds, &params, viewport_size, &range);
            break;
        }
        case MD_GL_REP_RIBBONS:
        case MD_GL_REP_CARTOON:
        case MD_GL_REP_TUBES: {
            uint32_t profile;
            float profile_scale[4];
            float max_extent;
            backbone_params(&profile, profile_scale, &max_extent, draw_op, mol);
            const gl_program_t program = draw_op->type == MD_GL_REP_RIBBONS ? shaders->ribbons[program_permutation] :
                                         draw_op->type == MD_GL_REP_CARTOON ? shaders->cartoon[program_permutation] : shaders->tubes[program_permutation];
            draw_backbone(program, mol, rep->atom_color, profile, profile_scale, prep, front_face);

            const bool gaps = draw_op->type == MD_GL_REP_RIBBONS ? draw_op->args.ribbons.chain_gaps :
                              draw_op->type == MD_GL_REP_CARTOON ? draw_op->args.cartoon.chain_gaps : false;
            if (gaps) {
                // As thick as the gap of a cartoon of the same thickness: the box is 0.1 thick at a scale of 1, the coil
                // 0.2. The gap of a cartoon is thinner than the coil (radius 0.2 at a scale of 1). Its colors are
                // blended smoothly from end to end.
                const float radius = draw_op->type == MD_GL_REP_RIBBONS ? CHAIN_GAP_RADIUS_TRACE * 0.2f * draw_op->args.ribbons.thickness_scale :
                                                                          CHAIN_GAP_RADIUS_TRACE * 0.2f * profile_scale[0];
                const bond_params_t gap = chain_gap_params(radius, MD_GL_BOND_MODE_SMOOTH, 0.0f, 0, args->picking_offset.gap_base);
                const draw_range_t gap_range = RANGE(gap_trace);
                draw_chain_gaps(shaders->chain_gaps[program_permutation], mol, rep->atom_color, true, &gap, viewport_size, &gap_range);
            }
            break;
        }
        case MD_GL_REP_SLABS: {
            const float radius = draw_op->args.slabs.radius;
            const float plate  = radius * CLAMP(draw_op->args.slabs.plate, 0.0f, 1.0f);
            const draw_range_t range = RANGE(slabs);
            draw_slabs(shaders->slabs[program_permutation], mol, rep->atom_color, radius, plate, max_length, args->picking_offset.ring_base, viewport_size, &range);
            break;
        }
        default:
            MD_LOG_ERROR("Representation had unexpected type");
            goto done;
        }
#undef RANGE
    }

    result = true;

done:
    POP_GPU_SECTION()

    if (mirrored) {
        glFrontFace(front_face);
    }
    glActiveTexture(GL_TEXTURE0 + INSTANCE_TEXTURE_UNIT);
    glBindTexture(GL_TEXTURE_BUFFER, 0);
    glActiveTexture(GL_TEXTURE0 + LIST_TEXTURE_UNIT);
    glBindTexture(GL_TEXTURE_BUFFER, 0);
    glActiveTexture(GL_TEXTURE0);

    POP_GPU_SECTION()

    md_temp_end(temp);

    return result;
}

// Draws the quads [offset, offset + count) of a kind (QUADS) with the program in use, its textures bound, for the
// copies of the range: all their blocks, or with an atom mask the blocks of the list
static void issue_quads(GLuint program, uint32_t offset, uint32_t count, const draw_range_t* range) {
    // Of a subset the primitives of its list, in blocks; with an atom mask the blocks of the list; otherwise all blocks
    const uint32_t blocks = range->exact ? DIV_UP(range->list.count, QUAD_BLOCK) : range->masked ? range->list.count : DIV_UP(count, QUAD_BLOCK);
    if (blocks == 0 || range->inst_count == 0 || count == 0) return;

    glUniform1ui(glGetUniformLocation(program, "u_prim_offset"), offset);
    glUniform1ui(glGetUniformLocation(program, "u_prim_count"),  range->exact ? range->list.count : count);
    glUniform1ui(glGetUniformLocation(program, "u_block_count"), blocks);
    glUniform1i (glGetUniformLocation(program, "u_block_list"),  range->masked && !range->exact ? (GLint)range->list.offset : -1);
    glUniform1i (glGetUniformLocation(program, "u_prim_list"),   range->exact ? (GLint)range->list.offset : -1);
    const GLint base_loc = glGetUniformLocation(program, "u_instance_base");
    if (range->masked) bind_list(range->buf);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ctx.quad_index.id);
    // gl_InstanceID is a signed int: the copies in as few draws as keep blocks * copies within it
    const uint32_t per_draw = MAX(1u, (uint32_t)(INT32_MAX / blocks));
    for (uint32_t i = 0; i < range->inst_count; i += per_draw) {
        const uint32_t n = MIN(per_draw, range->inst_count - i);
        glUniform1i(base_loc, (GLint)(range->inst_base + i));
        glDrawElementsInstanced(GL_TRIANGLES, QUAD_BLOCK * 6, GL_UNSIGNED_SHORT, 0, (GLsizei)(blocks * n));
    }
}

// Spheres by vertex pulling (spacefill.vert): one quad per atom, the atoms read from texture buffers, in blocks (QUADS).
// With an atom mask or a subset the list of the range is of the atoms themselves, drawn as the primitives of a subset.
static bool draw_space_fill(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, float scale, const draw_range_t* range) {
    ASSERT(mol);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_RADIUS].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(atom_color.id);

    if (!ctx.quad_index.id || range->inst_count == 0 || (range->masked && range->list.count == 0)) {
        return false;
    }

    ubo_write(sizeof(gl_ubo_base_t), sizeof(scale), &scale);

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
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, mol->buffer[GL_BUFFER_ATOM_RADIUS].id);

    glUseProgram(program.id);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_pos"),   0);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_vel"),   1);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_flags"), 2);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_color"), 3);
    glUniform1i(glGetUniformLocation(program.id, "u_buf_atom_rad"),   4);

    // The list of a masked range holds the atoms, as a subset's holds its primitives
    draw_range_t atoms = *range;
    atoms.exact = range->masked;
    issue_quads(program.id, 0, mol->atom_count, &atoms);

    glUseProgram(0);

    glActiveTexture(GL_TEXTURE0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return true;
}

// Bonds (licorice, weak bonds) by vertex pulling (bond.vert): one quad per bond, the atoms and bonds are read from
// texture buffers, in blocks (QUADS)
static bool draw_bonds(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, bond_set_t bonds, const bond_params_t* params, const float viewport_size[2], const draw_range_t* range) {
    ASSERT(mol);
    ASSERT(params);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(atom_color.id);

    if (bonds.count == 0 || !bonds.atoms.id || !ctx.quad_index.id || range->inst_count == 0 || (range->masked && range->list.count == 0)) {
        return false;
    }

    bond_params_t p = *params;
    p.bond_orders = (p.bond_orders && bonds.order.id && bonds.order_ref.id) ? 1 : 0;
    p.bond_kinds  = (p.bond_kinds && bonds.order.id) ? 1 : 0;
    p.elastic_bonds = p.elastic_bonds ? 1 : 0;
    ubo_write(sizeof(gl_ubo_base_t), sizeof(bond_params_t), &p);

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

    if (p.bond_kinds || p.bond_orders) {
        glActiveTexture(GL_TEXTURE6);
        glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_6].id);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, bonds.order.id);
    }
    if (p.bond_orders) {
        glActiveTexture(GL_TEXTURE7);
        glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_7].id);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, bonds.order_ref.id);
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
    glUniform2f(glGetUniformLocation(program.id, "u_viewport_size"), viewport_size[0], viewport_size[1]);

    issue_quads(program.id, 0, bonds.count, range);

    glUseProgram(0);

    glActiveTexture(GL_TEXTURE0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return true;
}

// The dashed partners of the aromatic and delocalized bonds (md_gl_mol_set_bond_orders), drawn with the bonds of the
// licorice and ball and stick representations: thin dashed cylinders, one primitive per entry
static bool draw_bond_partners(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, const bond_params_t* bond_params, const float viewport_size[2], const draw_range_t* range) {
    ASSERT(mol);
    ASSERT(bond_params);
    if (mol->bond_partner_count == 0 || mol->bond_order_count != mol->bond_count) {
        return false;
    }
    bond_params_t params = *bond_params;
    params.dash_length = BOND_PARTNER_DASH_LENGTH * bond_params->radius;
    params.gap_length  = BOND_PARTNER_GAP_LENGTH  * bond_params->radius;
    params.bond_orders = 0;
    const bond_set_t partners = {mol->buffer[GL_BUFFER_BOND_ATOM_INDICES], {0}, {0}, {0}, mol->buffer[GL_BUFFER_BOND_PARTNER], mol->bond_partner_count};
    return draw_bonds(program, mol, atom_color, partners, &params, viewport_size, range);
}

// The gaps of the chains (md_gl_mol_set_chain_gaps) as dashed cylinders (bond.vert and bond.frag, CHAIN_GAPS): the
// entries of the cartoon and ribbons (trace), or those of licorice and ball and stick, one primitive per entry
static bool draw_chain_gaps(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, bool trace, const bond_params_t* params, const float viewport_size[2], const draw_range_t* range) {
    ASSERT(mol);
    ASSERT(params);
    ASSERT(atom_color.id);

    const uint32_t offset = trace ? 0 : mol->gap_trace_count;
    const uint32_t count  = trace ? mol->gap_trace_count : mol->gap_link_count;
    if (count == 0 || !mol->buffer[GL_BUFFER_CHAIN_GAP].id || !mol->buffer[GL_BUFFER_CHAIN_GAP_COMP].id || !ctx.quad_index.id || range->inst_count == 0 || (range->masked && range->list.count == 0)) {
        return false;
    }

    ubo_write(sizeof(gl_ubo_base_t), sizeof(bond_params_t), params);

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
    glUniform2f(glGetUniformLocation(program.id, "u_viewport_size"), viewport_size[0], viewport_size[1]);

    issue_quads(program.id, offset, count, range);

    glUseProgram(0);

    glActiveTexture(GL_TEXTURE0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return true;
}

// The circles of the aromatic rings (md_gl_mol_set_bond_orders), drawn with the bonds of the licorice and ball and stick
// representations, by vertex pulling (ring.vert): one quad per ring, in blocks as the bonds. They take the radius, the
// colors and the reach of the bonds (bond_params), and write ring_index_base + the index of the ring to the picking
// target.
static bool draw_aromatic_rings(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, const bond_params_t* bond_params, uint32_t ring_index_base, const float viewport_size[2], const draw_range_t* range) {
    ASSERT(mol);
    ASSERT(bond_params);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(atom_color.id);

    if (mol->aromatic_ring_count == 0 || mol->bond_order_count != mol->bond_count || !mol->buffer[GL_BUFFER_AROMATIC_RING].id || !ctx.quad_index.id || range->inst_count == 0 || (range->masked && range->list.count == 0)) {
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
    ubo_write(sizeof(gl_ubo_base_t), sizeof(ring_params_t), &params);

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
    glUniform2f(glGetUniformLocation(program.id, "u_viewport_size"), viewport_size[0], viewport_size[1]);

    issue_quads(program.id, 0, mol->aromatic_ring_count, range);

    glUseProgram(0);

    glActiveTexture(GL_TEXTURE0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return true;
}

// Slabs by vertex pulling (slab.vert): one quad per slab, as the aromatic rings
static bool draw_slabs(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, float radius, float plate, float max_length, uint32_t slab_index_base, const float viewport_size[2], const draw_range_t* range) {
    ASSERT(mol);
    ASSERT(atom_color.id);
    if (mol->slab_count == 0 || !mol->buffer[GL_BUFFER_SLAB].id || !ctx.quad_index.id || !(radius > 0.0f) || range->inst_count == 0 || (range->masked && range->list.count == 0)) {
        return false;
    }

    // After gl_ubo_base_t, as in slab.vert
    const struct {
        float    radius;
        float    plate;
        float    max_d2;
        uint32_t slab_index_base;
    } params = { radius, plate, max_length * max_length, slab_index_base };
    ubo_write(sizeof(gl_ubo_base_t), sizeof(params), &params);

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
    glUniform2f(glGetUniformLocation(program.id, "u_viewport_size"), viewport_size[0], viewport_size[1]);

    issue_quads(program.id, 0, mol->slab_count, range);

    glUseProgram(0);

    glActiveTexture(GL_TEXTURE0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return true;
}

// One instance of the mesh per control point (residue) and copy of the molecule, see backbone.vert: per group of
// copies of the same levels of detail (op_prep_t), consecutive chains of the same level together, for all
// copies of the group in one draw. With an atom mask, of the control points of the list in the chains, those of
// marked atoms. A group of mirror images is drawn with the other front face, as its triangles turn the other way.
static bool draw_backbone(gl_program_t program, const molecule_t* mol, gl_buffer_t atom_color, uint32_t profile, const float scale[4], const op_prep_t* prep, GLint front_face) {
    ASSERT(mol);
    ASSERT(profile < BACKBONE_PROFILE_COUNT);
    ASSERT(prep);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA_ORIENT].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_NEIGHBOR].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_RING_DATA].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(atom_color.id);

    const backbone_mesh_t* meshes = ctx.backbone_mesh[profile];
    const backbone_bounds_t* bounds = &mol->backbone_bounds;
    if (!meshes[0].index_buffer.id || mol->backbone_count == 0 || bounds->chain_count == 0 || prep->group_count == 0) {
        return false;
    }
    // With an atom mask the control points of marked atoms, of a subset those within it, in order
    const bool masked = prep->lists != NULL;
    const uint32_t* cps = masked ? prep->lists->cp_items : NULL;
    const uint32_t  cp_count = masked ? prep->lists->cp_item_count : 0;
    const uint32_t  cp_list = masked ? prep->lists->cps.offset : 0;
    if (masked && cp_count == 0) {
        return false;
    }
    if (masked) bind_list(prep->lists->buf);

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
    const GLint cp_offset_loc  = glGetUniformLocation(program.id, "u_cp_offset");
    const GLint run_length_loc = glGetUniformLocation(program.id, "u_run_length");
    const GLint cp_list_loc    = glGetUniformLocation(program.id, "u_cp_list");
    const GLint inst_base_loc  = glGetUniformLocation(program.id, "u_instance_base");

    for (uint32_t g = 0; g < prep->group_count; ++g) {
        const bb_group_t* group = prep->groups + g;
        const uint8_t* chain_lod = group->chain_lod;
        if (group->mirrored) {
            glFrontFace(front_face == GL_CCW ? GL_CW : GL_CCW);
        }
        glUniform1i(inst_base_loc, (GLint)group->inst_base);

        for (uint32_t lod = 0; lod < BACKBONE_LOD_COUNT; ++lod) {
            const backbone_mesh_t* mesh = &meshes[lod];
            bool bound = false;

            uint32_t c = 0;
            while (c < bounds->chain_count) {
                if (chain_lod[c] != lod) {
                    ++c;
                    continue;
                }
                uint32_t beg = bounds->chain_offset[c];
                while (c < bounds->chain_count && chain_lod[c] == lod) ++c;
                uint32_t end = bounds->chain_offset[c];

                // With an atom mask, the control points of the list within [beg, end)
                uint32_t list_offset = 0;
                if (masked) {
                    uint32_t lo = 0, hi = cp_count;
                    while (lo < hi) { const uint32_t mid = (lo + hi) / 2; if (cps[mid] < beg) lo = mid + 1; else hi = mid; }
                    const uint32_t first = lo;
                    hi = cp_count;
                    while (lo < hi) { const uint32_t mid = (lo + hi) / 2; if (cps[mid] < end) lo = mid + 1; else hi = mid; }
                    list_offset = cp_list + first;
                    beg = first;
                    end = lo;
                }
                const uint32_t run = end - beg;
                if (run == 0) continue;

                const uint32_t segments = mesh->segments;
                if (!bound) {
                    struct {
                        float    scale[4];
                        uint32_t res[4];
                        uint32_t rings[4];
                    } params = {
                        {scale[0], scale[1], scale[2], scale[3]},
                        {segments, mesh->profile_count, mesh->outline_count, profile == BACKBONE_PROFILE_BOX ? mesh->face_subdivisions : mesh->dome_rings},
                        {MD_GL_BACKBONE_SEGMENT_COUNT, MD_GL_BACKBONE_SEGMENT_COUNT / segments, 0, 0},
                    };
                    ubo_write(sizeof(gl_ubo_base_t), sizeof(params), &params);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->index_buffer.id);
                    bound = true;
                }

                glUniform1ui(cp_offset_loc, masked ? 0 : beg);
                glUniform1i(cp_list_loc, masked ? (GLint)list_offset : -1);
                glUniform1ui(run_length_loc, run);
                // gl_InstanceID is a signed int: the copies in as few draws as keep run * copies within it
                const uint32_t per_draw = MAX(1u, (uint32_t)(INT32_MAX / run));
                for (uint32_t i = 0; i < group->inst_count; i += per_draw) {
                    const uint32_t n = MIN(per_draw, group->inst_count - i);
                    if (i > 0) glUniform1i(inst_base_loc, (GLint)(group->inst_base + i));
                    glDrawElementsInstanced(GL_TRIANGLES, (GLsizei)mesh->index_count, GL_UNSIGNED_SHORT, 0, (GLsizei)(run * n));
                }
                if (group->inst_count > per_draw) glUniform1i(inst_base_loc, (GLint)group->inst_base);
            }
        }
        if (group->mirrored) {
            glFrontFace(front_face);
        }
    }

    glUseProgram(0);

    glActiveTexture(GL_TEXTURE0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

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
