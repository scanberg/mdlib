#include <md_gfx.h>

#include <md_molecule.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>

#include <GL/gl3w.h>

#include <string.h>     // memset, memcpy, strstr
#include <stddef.h>
#include <float.h>

typedef mat4_t mat4;
typedef vec3_t vec3;
typedef vec4_t vec4;
typedef uint32_t uint;
typedef struct UniformData UniformData;
typedef struct DebugData DebugData;

typedef struct DrawArraysIndirectCommand DrawArraysIndirectCommand;
typedef struct DrawElementsIndirectCommand DrawElementsIndirectCommand;
typedef struct DispatchIndirectCommand DispatchIndirectCommand;
typedef struct DrawSpheresIndirectCommand DrawSpheresIndirectCommand;
typedef struct RastSpheresIndirectCommand RastSpheresIndirectCommand;
typedef struct DrawParameters DrawParameters;
typedef struct InstanceData InstanceData;
typedef struct ClusterBounds ClusterBounds;

// Common shared header between
#include <shaders/gfx/common.h>

#define CULL_ATOMS 1
#define CULL_GROUPS 0

// For subdividing backbone segments
#define MIN_SPLINE_SUBDIVISION_COUNT 1
#define MAX_SPLINE_SUBDIVISION_COUNT 32

#define MAX_INSTANCE_COUNT (1024*1024)

#define MIN_COUNT_PER_CLUSTER 128
#define MAX_COUNT_PER_CLUSTER 255

#define MAX_STRUCTURE_COUNT 64
#define MAX_REPRESENTATION_COUNT 64
#define MAX_HANDLE_COUNT (MAX_STRUCTURE_COUNT + MAX_REPRESENTATION_COUNT)

#define INVALID_INDEX (~0U)

#ifndef GL_PUSH_GPU_SECTION
#define GL_PUSH_GPU_SECTION(lbl)                                                                \
{                                                                                               \
    if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
}
#endif

#ifndef GL_POP_GPU_SECTION
#define GL_POP_GPU_SECTION()                \
{                                           \
    if (glPopDebugGroup) glPopDebugGroup(); \
}
#endif

static inline uint32_t previous_pow2(uint32_t v) {
    uint32_t r = 1;
    while (r * 2 < v)
        r *= 2;
    return r;
}

static inline bool is_ortho_proj_matrix(const mat4_t P) {
    return P.elem[2][3] == 0.0f;
}

static inline float extract_znear(const mat4_t Pinv) {
    return 1.f / (Pinv.elem[3][3] - Pinv.elem[2][3]);
}

static inline float extract_zfar(const mat4_t Pinv) {
    return 1.f / (Pinv.elem[3][3] + Pinv.elem[2][3]);
}

static inline vec2_t extract_jitter_uv(const mat4_t P) {
    if (is_ortho_proj_matrix(P)) {
        return (vec2_t) {
            -P.elem[3][0] * 0.5f,
            -P.elem[3][1] * 0.5f,
        };
    }
    else {
        return (vec2_t) {
            P.elem[2][0] * 0.5f,
            P.elem[2][1] * 0.5f,
        };
    }
}

static inline uint32_t compute_mip_count(uint32_t width, uint32_t height) {
    uint32_t mips = 1;
    uint32_t max_dim = MAX(width, height);
    while (max_dim >>= 1) ++mips;
    return mips;
}

typedef struct gl_buffer_t {
    uint32_t id;
    uint32_t size;
} gl_buffer_t;

typedef struct gl_texture_t {
    uint32_t id;
    uint32_t width;
    uint32_t height;
    uint32_t levels;
} gl_texture_t;

typedef struct gl_program_t {
    uint32_t id;
    uint32_t flags;     // Compute_bit, vertex_bit, fragment_bit etc. whatever stages this program contains
} gl_program_t;

typedef uint32_t gl_version_t;

// Internal packed structures which maps to the layout used in GPU-memory
// Pack the first index as an absolute offset and second as a relative to the first
typedef struct bond_t {
    uint32_t first : 20;
    uint32_t other : 12;
} bond_t;

typedef struct range_t {
    uint32_t offset;
    uint32_t count;
} range_t;

typedef struct segment_t {
    uint32_t offset : 20;
    uint32_t ca_idx : 4;
    uint32_t c_idx  : 4;
    uint32_t o_idx  : 4;
} segment_t;

typedef struct instance_t {
    mat4_t  transform;
    vec4_t  min_xyzr;
    vec4_t  max_xyzr;
    range_t atom_range;
    range_t cluster_range;
} instance_t;

typedef struct allocation_field_t {
    uint32_t offset;
    uint32_t capacity;
} allocation_field_t;

typedef struct allocation_range_t {
    uint32_t offset;
    uint32_t count;
    uint32_t capacity;
} allocation_range_t;

typedef struct structure_t {
    allocation_range_t atom;
    allocation_range_t bond;
    allocation_range_t bb_seg;
    allocation_range_t bb_range;
    allocation_range_t cluster;

    // Dynamic array
    range_t*    groups;

    instance_t  base_instance;
    instance_t* instances;

    uint16_t h_idx; // Handle index
} structure_t;

typedef struct representation_t {
    md_gfx_rep_type_t type;
    md_gfx_rep_attr_t attr;

    allocation_range_t color;

    uint16_t h_idx; // Handle index
} representation_t;

typedef struct gl_context_t {
    gl_version_t version;

    // For fixed pipeline framebuffer
    uint32_t vao_id;
    uint32_t fbo_id;
    uint32_t fbo_width;
    uint32_t fbo_height;

    gl_texture_t depth_tex;         // D32F   32-BITS depth buffer for depth testing in fixed pipeline
    gl_texture_t color_tex;         // RGBA8  32-BITS
    gl_texture_t index_tex;         // R32UI  32-BITS
    gl_texture_t normal_tex;        // RG8    16-BITS
    gl_texture_t ss_vel_tex;        // RG16F  32-BITS

    gl_texture_t depth_pyramid_tex; // R32F   32-BITS depth pyramid for occlusion culling (power of two, not conforming in size to others)

    // This should ideally be a texture for better locality when operating with in shaders, but no vendor seems to support that extension yet...
    gl_buffer_t visibility_buf;     // R64UI  64-BITS hi: 32-bit depth + lo: 32-bit payload, for compute rasterization of spheres

    // Structure Buffers
    gl_buffer_t position_buf;
    gl_buffer_t radius_buf;
    gl_buffer_t flags_buf;
    gl_buffer_t bond_buf;
    gl_buffer_t segment_buf;
    gl_buffer_t sec_str_buf;
    gl_buffer_t bb_range_buf;

    // Cluster Count in length
    gl_buffer_t cluster_bounds_buf;
    gl_buffer_t cluster_range_buf;

    // Atom count in length
    gl_buffer_t cluster_data_pos_rad_buf;
    gl_buffer_t cluster_data_atom_idx_buf;

    allocation_field_t atom;
    allocation_field_t bond;
    allocation_field_t bb_seg;
    allocation_field_t bb_range;
    allocation_field_t cluster;
    allocation_field_t transform;

    // Representation
    gl_buffer_t color_buf;

    allocation_field_t color;

    gl_program_t spacefill_prog;
    gl_program_t licorice_prog;
    gl_program_t ribbons_prog;
    gl_program_t cartoon_prog;
    gl_program_t gaussian_surface_prog;
    gl_program_t solvent_excluded_surface_prog;

    // Computes cluster boundaries and writes compressed coordinates within that space
    gl_program_t cluster_compute_data_prog;

    // This culls instances supplied by the CPU, stores them into cluster instances and adds them to visible or occluded list
    gl_program_t instance_cull_prog;
    gl_program_t instance_cull_late_prog;
    gl_program_t instance_cluster_cull_prog;
    gl_program_t cluster_cull_late_prog;

    gl_program_t cluster_raster_prog;

    // This sets up a draw element write command by writing element indices
    gl_program_t cluster_write_element_indices_prog;

    // Computes depth pyramid
    gl_program_t vis_reduce_prog;
    gl_program_t depth_reduce_prog;

    // Write back picking data
    gl_program_t write_picking_prog;

    // Composes final image and outputs
    gl_program_t compose_prog;

    // Transient GPU resident buffers
    gl_buffer_t ubo_buf;
    gl_buffer_t vertex_buf;
    gl_buffer_t index_buf;

    gl_buffer_t draw_ind_param_buf;
    
    gl_buffer_t draw_sphere_idx_buf;

    gl_buffer_t instance_vis_buf; // Visible instance indices
    gl_buffer_t instance_occ_buf; // Occluded instance indices

    gl_buffer_t cluster_instance_data_buf;
    gl_buffer_t cluster_instance_vis_idx_buf;   // Visible cluster instance indices
    gl_buffer_t cluster_instance_occ_idx_buf;   // Occluded cluster instance indices

    // List of cluster instances to rasterize
    gl_buffer_t cluster_instance_rast_idx_buf;

    // Buffers which we fill from CPU side to kick of work
    gl_buffer_t transform_buf;
    gl_buffer_t instance_data_buf;

    gl_buffer_t DEBUG_BUF;
    DebugData* debug_data;

    mat4_t prev_view_mat;
    mat4_t prev_proj_mat;
    uint32_t frame_idx;

    structure_t structures[MAX_STRUCTURE_COUNT];
    representation_t representations[MAX_REPRESENTATION_COUNT];

    md_gfx_handle_t handles[MAX_HANDLE_COUNT];

    uint16_t representation_count;
    uint16_t structure_count;
    uint16_t handle_next_idx;
} gl_context_t;



// These are absolute for now, but may be compacted later on
typedef struct position_t {
    float x, y, z;
} position_t;

typedef struct velocity_t {
    float x, y, z;
} velocity_t;

STATIC_ASSERT(sizeof(bond_t) == sizeof(uint32_t), "Invalid size of bond");
STATIC_ASSERT(sizeof(segment_t) == sizeof(uint32_t), "Invalid size of segment");

static inline gl_texture_t gl_texture_create(uint32_t width, uint32_t height, uint32_t levels, uint32_t format) {
    gl_texture_t tex = {0};
    glCreateTextures(GL_TEXTURE_2D, 1, &tex.id);
    glTextureStorage2D(tex.id, levels, format, width, height);
    tex.width = width;
    tex.height = height;
    tex.levels = levels;
    return tex;
}

static inline void gl_texture_destroy(gl_texture_t* tex) {
    if (tex->id) {
        glDeleteTextures(1, &tex->id);
    }
    *tex = (gl_texture_t){0};
}

static inline gl_buffer_t gl_buffer_create(uint32_t size, const void* data, uint32_t flags) {
    gl_buffer_t buf = {0};
    glCreateBuffers(1, &buf.id);
    glNamedBufferStorage(buf.id, size, data, flags);
    buf.size = size;
    return buf;
}

static inline void gl_buffer_destroy(gl_buffer_t* buf) {
    if (buf->id) {
        glDeleteBuffers(1, &buf->id);
    }
    *buf = (gl_buffer_t){0};
}

static inline void gl_buffer_set_sub_data(gl_buffer_t buf, uint32_t byte_offset, uint32_t byte_size, const void* data) {
    glNamedBufferSubData(buf.id, byte_offset, byte_size, data);
}

static bool extract_version_string(char version_str[], uint32_t version_cap, const char** src) {
    ASSERT(version_str);
    ASSERT(version_cap > 0);
    ASSERT(src);
    if (strncmp(*src, "#version", 8) == 0) {
        const char* endline;
        if ((endline = strchr(*src, '\n')) != NULL) {
            ++endline;
            const size_t length = endline - *src;
            strncpy(version_str, *src, MIN(version_cap, length));
        }
        *src = endline;
        return true;
    }
    return false;
}

typedef struct include_file_t {
    const char* filename;
    const char* source;
} include_file_t;

#define APPEND_STR(a, str, alloc)   md_array_push_array(a, str.ptr, str.len, alloc)
#define APPEND_CSTR(a, cstr, alloc) md_array_push_array(a, cstr, (int64_t)strlen(cstr), alloc)
#define APPEND_CHAR(a, ch, alloc) md_array_push(a, ch, alloc);
#define APPEND_LINE(a, line, alloc) APPEND_STR(a, str_printf(alloc, "#line %u\n", line), alloc)

static bool compile_shader_from_source(GLuint shader, const char* source, const char* defines, const include_file_t* include_files, uint32_t include_file_count) {
    if (!source) {
        MD_LOG_ERROR("Missing shader source");
        return false;
    }

    str_t str = str_from_cstr(source);

    str_t version_str;
    if (!str_extract_line(&version_str, &str) || !str_eq_cstr_n(version_str, "#version", 8)) {
        MD_LOG_ERROR("Missing version as first line in shader!");
        return false;
    }

    md_allocator_i* alloc = md_get_temp_allocator();

    char* buf = 0;
    md_array_ensure(buf, KILOBYTES(4), alloc);

    uint32_t line_count = 0;
    APPEND_STR(buf, version_str, alloc);
    line_count += 1;

    str_t glsl_define = STR_LIT("#ifndef GLSL\n#define GLSL\n#endif");
    APPEND_STR(buf, glsl_define, alloc);
    line_count += 3;

    if (defines) {
        str_t def = str_from_cstr(defines);
        line_count += (uint32_t)str_count_occur_char(def, '\n');        
        APPEND_LINE(buf, line_count, alloc);
    }

    str_t line;
    while (str_extract_line(&line, &str)) {
        if (str_eq_cstr_n(line, "#include", 8)) {
            if (include_file_count == 0 || !include_files) {
                MD_LOG_ERROR("Failed to parse include in shader file: no include files have been supplied");
                return false;
            }
            str_t file = str_trim(str_substr(line, 8, SIZE_MAX));
            if (!file.len) {
                MD_LOG_ERROR("Failed to parse include in shader file: file is missing");
                return false;
            }
            char beg_delim = file.ptr[0];
            char end_delim = file.ptr[file.len-1];
            if (!(beg_delim == '"' || beg_delim == '<') ||
                !(end_delim == '"' || end_delim == '>') ||
                 (beg_delim == '"' && end_delim != '"') ||
                 (beg_delim == '<' && end_delim != '>')) {
                MD_LOG_ERROR("Failed to parse include in shader file: missing or mismatched delimiters");
                return false;
            }
            file = str_substr(file, 1, file.len-2);
            str_t src = {0};
            for (uint32_t i = 0; i < include_file_count; ++i) {
                if (str_eq_cstr(file, include_files[i].filename)) {
                    src = str_from_cstr(include_files[i].source);
                    break;
                }
            }
            if (str_empty(src)) {
                MD_LOG_ERROR("Failed to parse include in shader file: could not find include file '%.*s'", (int)file.len, file.ptr);
                return false;
            }
            APPEND_STR(buf, src, alloc);
            APPEND_LINE(buf, line_count, alloc);
        }
        else {
            APPEND_STR(buf, line, alloc);
            line_count += 1;
        }
    }
    APPEND_CHAR(buf, '\0', alloc);

    const char* c_src[] = {buf};
    glShaderSource(shader, ARRAY_SIZE(c_src), c_src, 0);
    
    GLint success;
    glCompileShader(shader);
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char err_buf[1024];
        glGetShaderInfoLog(shader, ARRAY_SIZE(err_buf), NULL, err_buf);
        MD_LOG_ERROR("Shader compile error:\n%s\n", err_buf);
        return false;
    }
    
    return true;
}

static bool compile_shader_from_file(GLuint shader, const char* filename, const char* defines, const include_file_t* include_files, uint32_t include_file_count) {
    str_t src = load_textfile(str_from_cstr(filename), md_get_temp_allocator());
    if (!str_empty(src)) {
        const bool success = compile_shader_from_source(shader, src.ptr, defines, include_files, include_file_count);
        const md_log_type_t log_type = success ? MD_LOG_TYPE_INFO : MD_LOG_TYPE_ERROR;
        const char* res_str = success ? "Success" : "Fail";
        md_logf(log_type,  "Compiling shader %-40s %s", filename, res_str);
        return success;
    } else {
        MD_LOG_ERROR("Could not open file file '%s'", filename);
        return false;
    }
}

static bool link_program(GLuint program, const GLuint shader[], uint32_t count) {
    ASSERT(program);
     
    for (uint32_t i = 0; i < count; i++) {
        GLint compile_status;
        glGetShaderiv(shader[i], GL_COMPILE_STATUS, &compile_status);
        if (glIsShader(shader[i]) && compile_status) {
            glAttachShader(program, shader[i]);
        } else {
            MD_LOG_ERROR("Program link error: One or more shaders are invalid\n");
            return false;
        }
    }
    
    GLint success;
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        char err_buf[1024];
        glGetProgramInfoLog(program, ARRAY_SIZE(err_buf), NULL, err_buf);
        MD_LOG_ERROR("Program link error:\n%s\n", err_buf);
        return false;
    }
    
    for (uint32_t i = 0; i < count; i++) {
        glDetachShader(program, shader[i]);
    }
    
    return true;
}

static inline uint32_t gl_extract_shader_type_from_file_extension(const char* filename) {
    if (strstr(filename, ".comp")) {
        return GL_COMPUTE_SHADER;
    } else if (strstr(filename, ".vert")) {
        return GL_VERTEX_SHADER;
    } else if (strstr(filename, ".geom")) {
        return GL_GEOMETRY_SHADER;
    } else if (strstr(filename, ".frag")) {
        return GL_FRAGMENT_SHADER;
    }
    return 0;
}

static inline uint32_t gl_get_shader_type_bit(uint32_t shader_type) {
    switch (shader_type) {
    case GL_VERTEX_SHADER:   return GL_VERTEX_SHADER_BIT;
    case GL_FRAGMENT_SHADER: return GL_FRAGMENT_SHADER_BIT;
    case GL_COMPUTE_SHADER:  return GL_COMPUTE_SHADER_BIT;
    }
    return 0;
}

static gl_program_t gl_program_create_from_files(const char* shader_files[], uint32_t shader_file_count, const char* defines, const include_file_t* include_files, uint32_t include_file_count) {
    uint32_t shaders[8] = {0};
    ASSERT(shader_file_count < ARRAY_SIZE(shaders));
    gl_program_t prog = {0};

    uint32_t flags = 0;
    for (uint32_t i = 0; i < shader_file_count; ++i) {
        uint32_t type = gl_extract_shader_type_from_file_extension(shader_files[i]);
        flags |= gl_get_shader_type_bit(type);
        shaders[i] = glCreateShader(type);
        compile_shader_from_file(shaders[i], shader_files[i], defines, include_files, include_file_count);
    }
    
    uint32_t id = glCreateProgram();
    if (link_program(id, shaders, shader_file_count)) {
        prog.id = id;
        prog.flags = flags;  
    } else {
        glDeleteProgram(id);
    }

    for (uint32_t i = 0; i < shader_file_count; ++i) {
        glDeleteShader(shaders[i]);
    }

    return prog;
}

static gl_context_t ctx = {0};

static inline bool validate_context() {
    if (ctx.version == 0) {
        MD_LOG_ERROR("gfx context has not been initialized");
        return false;
    }
    return true;
}

void init_handles() {
    for (uint32_t i = 0; i < ARRAY_SIZE(ctx.handles); ++i) {
        ctx.handles[i].idx = (uint16_t)(i + 1);
        // Scramble generations a bit
        ctx.handles[i].gen = (uint16_t)((i + 1) * (1236285 ^ 23645) + 123745);
    }
    ctx.handle_next_idx = 0;
}

static inline void free_handle(uint16_t idx) {
    ctx.handles[idx].gen += 1;
    ctx.handles[idx].idx = (uint16_t)ctx.handle_next_idx;
    ctx.handle_next_idx = idx;
}

static inline uint16_t alloc_handle() {
    if (ctx.handle_next_idx >= ARRAY_SIZE(ctx.handles)) {
        ASSERT(false);
    }
    uint16_t idx = ctx.handle_next_idx;
    ctx.handle_next_idx = ctx.handles[idx].idx;
    return idx;
}

static inline bool validate_handle(md_gfx_handle_t id) {
    if (id.idx >= ARRAY_SIZE(ctx.handles)) return false;
    return ctx.handles[id.idx].gen == id.gen;
}

static inline representation_t* get_representation(md_gfx_handle_t id) {
    if (id.idx >= ARRAY_SIZE(ctx.handles)) return NULL;
    return (ctx.handles[id.idx].gen == id.gen) ? &ctx.representations[ctx.handles[id.idx].idx] : NULL;
}

static inline structure_t* get_structure(md_gfx_handle_t id) {
    if (id.idx >= ARRAY_SIZE(ctx.handles)) return NULL;
    return (ctx.handles[id.idx].gen == id.gen) ? &ctx.structures[ctx.handles[id.idx].idx] : NULL;
}

bool md_gfx_initialize(const char* shader_base_dir, uint32_t width, uint32_t height, md_gfx_config_flags_t flags) {
    if (!ctx.version) {
        if (gl3wInit() != GL3W_OK) {
            MD_LOG_ERROR("Could not load OpenGL extensions");
            return false;
        }

        GLint major, minor;
        glGetIntegerv(GL_MAJOR_VERSION, &major);
        glGetIntegerv(GL_MINOR_VERSION, &minor);

        if (major < 4 && minor < 6) {
            MD_LOG_ERROR("OpenGL version %i.%i is not supported, this renderer requires version 4.6", major, minor);
            return false;
        }

        glCreateFramebuffers(1, &ctx.fbo_id);
        glCreateVertexArrays(1, &ctx.vao_id);

        ctx.ubo_buf = gl_buffer_create(KILOBYTES(1), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

        ctx.instance_data_buf = gl_buffer_create(sizeof(InstanceData) * MAX_INSTANCE_COUNT, NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.instance_vis_buf  = gl_buffer_create(sizeof(uint) * MAX_INSTANCE_COUNT, NULL, 0);
        ctx.instance_occ_buf  = gl_buffer_create(sizeof(uint) * MAX_INSTANCE_COUNT, NULL, 0);

        // Fully GPU Resident Buffers
        ctx.draw_ind_param_buf      = gl_buffer_create(sizeof(DrawParameters), NULL, GL_DYNAMIC_STORAGE_BIT);

        ctx.cluster_instance_data_buf       = gl_buffer_create(MEGABYTES(256), NULL, 0);
        ctx.cluster_instance_vis_idx_buf    = gl_buffer_create(MEGABYTES(64), NULL, 0);
        ctx.cluster_instance_occ_idx_buf    = gl_buffer_create(MEGABYTES(64), NULL, 0);

        ctx.vertex_buf                      = gl_buffer_create(MEGABYTES(32), NULL, 0);
        ctx.index_buf                       = gl_buffer_create(MEGABYTES(32), NULL, 0);

        // @TODO: This can be improved with drawMultiIndirectArraysCount
        // The current approach is suboptimal using individual elements for cluster_idx + sphere_idx
        ctx.draw_sphere_idx_buf             = gl_buffer_create(MEGABYTES(128), NULL, 0);

        ctx.cluster_instance_rast_idx_buf   = gl_buffer_create(MEGABYTES(256),  NULL, 0);

        ctx.DEBUG_BUF = gl_buffer_create(KILOBYTES(1), NULL, GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT);
        ctx.debug_data = glMapNamedBuffer(ctx.DEBUG_BUF.id, GL_READ_ONLY);

        // Storage buffers for structure data
        ctx.atom.capacity = 4 * 1000 * 1000;
        ctx.position_buf        = gl_buffer_create(ctx.atom.capacity * sizeof(vec3_t),   NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.radius_buf          = gl_buffer_create(ctx.atom.capacity * sizeof(float),    NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.flags_buf           = gl_buffer_create(ctx.atom.capacity * sizeof(uint32_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

        ctx.cluster_data_atom_idx_buf = gl_buffer_create(ctx.atom.capacity * sizeof(uint32_t), NULL, GL_DYNAMIC_STORAGE_BIT);
        ctx.cluster_data_pos_rad_buf  = gl_buffer_create(ctx.atom.capacity * sizeof(struct CompressedPosRad), NULL, 0);
        
        ctx.bond.capacity = 2 * 1000 * 1000;
        ctx.bond_buf     = gl_buffer_create(ctx.bond.capacity * sizeof(bond_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

        ctx.bb_seg.capacity = 1 * 1000 * 1000;
        ctx.segment_buf  = gl_buffer_create(ctx.bb_seg.capacity * sizeof(md_gfx_backbone_segment_t),    NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.sec_str_buf  = gl_buffer_create(ctx.bb_seg.capacity * sizeof(md_gfx_secondary_structure_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

        ctx.bb_range.capacity = 10 * 1000;
        ctx.bb_range_buf = gl_buffer_create(ctx.bb_range.capacity * sizeof(range_t),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

        ctx.cluster.capacity = 1 * 1000 * 1000;
        ctx.cluster_range_buf  = gl_buffer_create(ctx.cluster.capacity * sizeof(uint32_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.cluster_bounds_buf = gl_buffer_create(ctx.cluster.capacity * sizeof(uint32_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

        ctx.transform.capacity = 100000;
        ctx.transform_buf = gl_buffer_create(ctx.transform.capacity * sizeof(mat4_t),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

        // Storage buffers for representation data
        ctx.color.capacity = 4 * 1000 * 1000;
        ctx.color_buf = gl_buffer_create(ctx.color.capacity * sizeof(md_gfx_color_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

#define CONCAT_STR(STR_A, STR_B) (str_printf(md_get_temp_allocator(), "%s%s", STR_A, STR_B).ptr)
#define BASE shader_base_dir

        str_t common_src = load_textfile(str_from_cstr(CONCAT_STR(BASE, "/gfx/common.h")), md_get_temp_allocator());
        if (str_empty(common_src)) {
            MD_LOG_ERROR("Failed to read common.h required for shaders");
            return false;
        }

        str_t culling_src = load_textfile(str_from_cstr(CONCAT_STR(BASE, "/gfx/culling.glsl")), md_get_temp_allocator());
        if (str_empty(culling_src)) {
            MD_LOG_ERROR("Failed to read culling.glsl required for shaders");
            return false;
        }

        include_file_t include_files[] = {
            {
                .filename = "common.h",
                .source = common_src.ptr,
            },
            {
                .filename = "culling.glsl",
                .source = culling_src.ptr,
            },
        };

        // Compile shader programs
        {
            const char* files[] = {
                CONCAT_STR(BASE, "/gfx/spacefill.vert"),
                CONCAT_STR(BASE, "/gfx/spacefill.geom"),
                CONCAT_STR(BASE, "/gfx/spacefill.frag"),
            };
            ctx.spacefill_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        {
            const char* files[] = { CONCAT_STR(BASE, "/gfx/cluster_compute_data.comp") };
            ctx.cluster_compute_data_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        {
            const char* files[] = { CONCAT_STR(BASE, "/gfx/instance_cull.comp") };
            ctx.instance_cull_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        {
            const char* files[] = { CONCAT_STR(BASE, "/gfx/instance_cull_late.comp") };
            ctx.instance_cull_late_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        {
            const char* files[] = { CONCAT_STR(BASE, "/gfx/instance_cluster_cull.comp") };
            ctx.instance_cluster_cull_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        {
            const char* files[] = { CONCAT_STR(BASE, "/gfx/cluster_cull_late.comp") };
            ctx.cluster_cull_late_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        {
            const char* files[] = { CONCAT_STR(BASE, "/gfx/cluster_write_elem_idx.comp") };
            ctx.cluster_write_element_indices_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        {
            const char* files[] = { CONCAT_STR(BASE, "/gfx/cluster_rasterize.comp") };
            ctx.cluster_raster_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        {
            const char* files[] = { CONCAT_STR(BASE, "/gfx/vis_reduce.comp") };
            ctx.vis_reduce_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        {
            const char* files[] = { CONCAT_STR(BASE, "/gfx/depth_reduce.comp") };
            ctx.depth_reduce_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        {
            const char* files[] = { CONCAT_STR(BASE, "/gfx/write_picking_data.comp") };
            ctx.write_picking_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        {
            const char* files[] = { CONCAT_STR(BASE, "/gfx/fs_quad.vert"), CONCAT_STR(BASE, "/gfx/compose.frag"), };
            ctx.compose_prog = gl_program_create_from_files(files, ARRAY_SIZE(files), 0, include_files, ARRAY_SIZE(include_files));
        }

        if (!ctx.spacefill_prog.id ||
            !ctx.cluster_compute_data_prog.id ||
            !ctx.instance_cull_prog.id ||
            !ctx.instance_cull_late_prog.id ||
            !ctx.instance_cluster_cull_prog.id ||
            !ctx.cluster_cull_late_prog.id ||
            !ctx.cluster_write_element_indices_prog.id ||
            !ctx.vis_reduce_prog.id ||
            !ctx.depth_reduce_prog.id ||
            !ctx.write_picking_prog.id ||
            !ctx.compose_prog.id)
        {
            MD_LOG_ERROR("Failed to compile one or more shader programs.");
            md_gfx_shutdown();
            return false;
        }

        init_handles();

        // We only set version when everything is successful
        ctx.version = major * 100 + minor * 10;
    }

    if (ctx.version) {
        if (ctx.fbo_width != width || ctx.fbo_height != height) {
            gl_texture_destroy(&ctx.depth_tex);
            ctx.depth_tex = gl_texture_create(width, height, 1, GL_DEPTH_COMPONENT32F);
            glTextureParameteri(ctx.depth_tex.id, GL_TEXTURE_REDUCTION_MODE_ARB, GL_MAX);

            gl_texture_destroy(&ctx.color_tex);
            ctx.color_tex = gl_texture_create(width, height, 1, GL_RGBA8);

            gl_texture_destroy(&ctx.normal_tex);
            ctx.normal_tex = gl_texture_create(width, height, 1, GL_RG8);

            gl_texture_destroy(&ctx.index_tex);
            ctx.index_tex = gl_texture_create(width, height, 1, GL_RGBA8);

            gl_texture_destroy(&ctx.ss_vel_tex);
            ctx.ss_vel_tex = gl_texture_create(width, height, 1, GL_RG16F);

            gl_buffer_destroy(&ctx.visibility_buf);
            ctx.visibility_buf = gl_buffer_create(width * height * (uint32_t)sizeof(uint64_t), 0, 0);

            glNamedFramebufferTexture(ctx.fbo_id, GL_DEPTH_ATTACHMENT, ctx.depth_tex.id,  0);
            glNamedFramebufferTexture(ctx.fbo_id, GL_COLOR_ATTACHMENT0 + COLOR_TEX_ATTACHMENT,  ctx.color_tex.id,  0);
            glNamedFramebufferTexture(ctx.fbo_id, GL_COLOR_ATTACHMENT0 + NORMAL_TEX_ATTACHMENT, ctx.normal_tex.id, 0);
            glNamedFramebufferTexture(ctx.fbo_id, GL_COLOR_ATTACHMENT0 + SS_VEL_TEX_ATTACHMENT, ctx.ss_vel_tex.id, 0);
            glNamedFramebufferTexture(ctx.fbo_id, GL_COLOR_ATTACHMENT0 + INDEX_TEX_ATTACHMENT,  ctx.index_tex.id,  0);

            ctx.fbo_width = width;
            ctx.fbo_height = height;
        }
        uint32_t w = previous_pow2(width);
        uint32_t h = previous_pow2(height);
        uint32_t l = compute_mip_count(w, h);
        if (ctx.depth_pyramid_tex.width != w || ctx.depth_pyramid_tex.height != h || ctx.depth_pyramid_tex.levels != l) {
            gl_texture_destroy(&ctx.depth_pyramid_tex);
            ctx.depth_pyramid_tex = gl_texture_create(w, h, l, GL_R32F);
            glTextureParameteri(ctx.depth_pyramid_tex.id, GL_TEXTURE_REDUCTION_MODE_ARB, GL_MAX);
            glTextureParameteri(ctx.depth_pyramid_tex.id, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
            glTextureParameteri(ctx.depth_pyramid_tex.id, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTextureParameteri(ctx.depth_pyramid_tex.id, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTextureParameteri(ctx.depth_pyramid_tex.id, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        }

        ctx.prev_proj_mat = mat4_ident();
        ctx.prev_view_mat = mat4_ident();
        ctx.frame_idx = 0;
    }

    return true;
}

void md_gfx_shutdown() {
    if (ctx.fbo_id) glDeleteFramebuffers(1, &ctx.fbo_id);
    gl_texture_destroy(&ctx.color_tex);
    gl_texture_destroy(&ctx.normal_tex);
    gl_texture_destroy(&ctx.ss_vel_tex);
    gl_texture_destroy(&ctx.index_tex);

    gl_texture_destroy(&ctx.depth_pyramid_tex);
    gl_texture_destroy(&ctx.depth_tex);

    gl_buffer_destroy(&ctx.visibility_buf);

    gl_buffer_destroy(&ctx.ubo_buf);
    gl_buffer_destroy(&ctx.vertex_buf);
    gl_buffer_destroy(&ctx.index_buf);
    //gl_buffer_destroy(&ctx.draw_op_buf);

    gl_buffer_destroy(&ctx.position_buf);
    gl_buffer_destroy(&ctx.cluster_data_atom_idx_buf);
    gl_buffer_destroy(&ctx.radius_buf);
    gl_buffer_destroy(&ctx.bond_buf);
    gl_buffer_destroy(&ctx.segment_buf);
    gl_buffer_destroy(&ctx.sec_str_buf);
    gl_buffer_destroy(&ctx.bb_range_buf);
    gl_buffer_destroy(&ctx.cluster_range_buf);
    gl_buffer_destroy(&ctx.cluster_data_pos_rad_buf);

    ctx = (gl_context_t){0};
}

allocation_range_t field_alloc(allocation_field_t* field, uint32_t size) {
    uint32_t offset = field->offset + size < field->capacity ? field->offset : 0;
    field->offset += size;
    return (allocation_range_t){offset, 0, size};
}

static void field_free(allocation_field_t* field, allocation_range_t range) {
    if (field->offset == range.offset + range.capacity) {
        field->offset -= range.capacity;
    }
}

md_gfx_handle_t md_gfx_structure_create_from_mol(const md_molecule_t* mol) {
    ASSERT(mol);
    return md_gfx_structure_create((uint32_t)mol->atom.count, (uint32_t)md_array_size(mol->bond.pairs), (uint32_t)mol->protein_backbone.count, (uint32_t)mol->protein_backbone.range.count, (uint32_t)mol->residue.count, (uint32_t)mol->instance.count);
}

md_gfx_handle_t md_gfx_structure_create(uint32_t atom_count, uint32_t bond_count, uint32_t backbone_segment_count, uint32_t backbone_range_count, uint32_t group_count, uint32_t instance_count) {
    if (!validate_context()) {
        return (md_gfx_handle_t){0};
    }

    if (ctx.structure_count >= MAX_STRUCTURE_COUNT) {
        MD_LOG_ERROR("Unable to create structure, pool is full!");
        return (md_gfx_handle_t){0};
    }
    uint16_t s_idx = (uint16_t)ctx.structure_count++;

    structure_t* s = &ctx.structures[s_idx];
    s->h_idx = alloc_handle();

    // Map handle index to resource index
    ctx.handles[s->h_idx].idx = s_idx;

    // Allocate internal fields
    s->atom     = field_alloc(&ctx.atom, atom_count);
    s->bond     = field_alloc(&ctx.bond, bond_count);
    s->bb_seg   = field_alloc(&ctx.bb_seg, backbone_segment_count);
    s->bb_range = field_alloc(&ctx.bb_range, backbone_range_count);
    s->cluster  = field_alloc(&ctx.cluster, DIV_UP(atom_count, MIN_COUNT_PER_CLUSTER)); // Reserve the maximum amount of clusters we need (this changes dynamically when recomputing clusters)

    md_array_ensure(s->groups, group_count, md_get_heap_allocator());
    md_array_ensure(s->instances, instance_count, md_get_heap_allocator());

    // Create external handle for user
    return (md_gfx_handle_t) {ctx.handles[s->h_idx].gen, s->h_idx};
}

bool md_gfx_structure_destroy(md_gfx_handle_t id) {
    if (!validate_handle(id)) {
        MD_LOG_ERROR("Failed to destroy structure: Handle is invalid");
        return false;
    }
    uint16_t h_idx = id.idx;
    uint16_t s_idx = ctx.handles[h_idx].idx;
    free_handle(h_idx);

    ASSERT(ctx.structure_count > 0);
    ctx.structure_count -= 1;

    structure_t* s = &ctx.structures[s_idx];
    
    field_free(&ctx.atom,       s->atom);
    field_free(&ctx.bond,       s->bond);
    field_free(&ctx.bb_seg,     s->bb_seg);
    field_free(&ctx.bb_range,   s->bb_range);

    md_array_free(s->instances, md_get_heap_allocator());

    if (ctx.structure_count > 1) {
        // Swap back and pop resource array to keep it packed
        ctx.structures[s_idx] = ctx.structures[ctx.structure_count];
        
        // Reassign handle indices to new location
        ctx.handles[ctx.structures[s_idx].h_idx].idx = s_idx;
    }

    return true;
}

typedef struct aabb_t {
    vec3_t min_box;
    vec3_t max_box;
} aabb_t;

uint32_t* create_cluster_ranges(uint32_t atom_range_beg, uint32_t atom_range_end, md_allocator_i* alloc) {
    uint32_t* ranges = 0;
    for (uint32_t i = atom_range_beg; i < atom_range_end;) {
        uint32_t cluster_offset = i;
        uint32_t cluster_size   = MIN(MAX_COUNT_PER_CLUSTER, atom_range_end - i);

        uint32_t cluster_range = (cluster_offset << 8) | cluster_size;
        md_array_push(ranges, cluster_range, alloc);
        i += cluster_size;
    }
    return ranges;
}

uint32_t* create_cluster_ranges_from_groups(range_t* groups, uint32_t group_count, md_allocator_i* alloc) {
    uint32_t* ranges = 0;
    uint32_t cluster_offset = 0;
    uint32_t cluster_size   = 0;

    for (uint32_t i = 0; i < group_count; ++i) {
        if (cluster_size + groups[i].count < MAX_COUNT_PER_CLUSTER) {
            cluster_size += groups[i].count;
        }
        else {
            uint32_t range = (cluster_offset << 8) | cluster_size;
            md_array_push(ranges, range, alloc);
            cluster_offset = groups[i].offset;
            cluster_size   = groups[i].count;
        }
    }

    const uint32_t range = (cluster_offset << 8) | cluster_size;
    md_array_push(ranges, range, alloc);

    return ranges;
}

void recompute_clusters2(structure_t* s, const float* x, const float* y, const float* z, uint32_t count, uint32_t byte_stride) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    // This is a bit fiddly to get the mapping right, especially if you have instances defined.
    // We make the assumption that if you have groups and instances defined, the instances are defined on the boarders of the defined groups.
    // And that all atoms belong to a group.

    const uint32_t grp_count = (uint32_t)md_array_size(s->groups);
    const uint32_t inst_count = (uint32_t)md_array_size(s->instances);

    uint32_t* grp_indices  = md_vm_arena_push(arena, sizeof(uint32_t) * count);

    aabb_t* grp_aabb = md_vm_arena_push(arena, sizeof(aabb_t) * grp_count);
    vec3_t* grp_cent = md_vm_arena_push(arena, sizeof(vec3_t) * grp_count);
    int32_t* grp_inst_idx = md_vm_arena_push(arena, sizeof(int32_t) * grp_count);

    for (uint32_t g_idx = 0; g_idx < grp_count; ++g_idx) {
        vec3_t aabb_min = {+FLT_MAX, +FLT_MAX, +FLT_MAX};
        vec3_t aabb_max = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
        uint32_t beg = s->groups[g_idx].offset;
        uint32_t end = s->groups[g_idx].offset + s->groups[g_idx].count;
        for (uint32_t i = beg; i < end; ++i) {
            vec3_t p = {x[i], y[i], z[i]};
            aabb_min = vec3_min(aabb_min, p);
            aabb_max = vec3_max(aabb_max, p);
            grp_indices[i] = g_idx;
            grp_inst_idx[i] = -1;
        }
        aabb_t aabb = {aabb_min, aabb_max};
        vec3_t cent = vec3_mul_f(vec3_add(aabb_min, aabb_max), 0.5f);
        grp_aabb[g_idx] = aabb;
        grp_cent[g_idx] = cent;
    }

    range_t* inst_grp_ranges = md_vm_arena_push(arena, sizeof(range_t) * md_array_size(s->instances));

    // Mark any atoms part of an instance, these are special cases
    for (uint32_t i_idx = 0; i_idx < md_array_size(s->instances); ++i_idx) {
        uint32_t first = s->instances[i_idx].atom_range.offset;
        uint32_t last  = s->instances[i_idx].atom_range.offset + s->instances[i_idx].atom_range.count - 1;
        range_t range = {grp_indices[first], grp_indices[last]};
        inst_grp_ranges[i_idx] = range;
        for (uint32_t i = range.offset; i < range.offset + range.count; ++i) {
            grp_inst_idx[i] = i_idx;
        }
    }

    uint32_t* grp_src_indices = md_vm_arena_push(arena, sizeof(uint32_t) * grp_count);
    md_util_sort_spatial_vec3(grp_src_indices, grp_cent, grp_count);

    uint32_t* src_indices = NULL;
    range_t* base_groups  = NULL;
    
    // Add according to group and ignore groups part of instances
    for (uint32_t j = 0; j < grp_count; ++j) {
        uint32_t g_idx = grp_src_indices[j];
        if (grp_inst_idx[g_idx] > -1) continue;
        md_array_push(base_groups, s->groups[g_idx], arena);

        uint32_t beg = s->groups[g_idx].offset;
        uint32_t end = s->groups[g_idx].offset + s->groups[g_idx].count;
        for (uint32_t i = beg; i < end; ++i) {
            md_array_push(src_indices, i, arena);
        }
    }

    uint32_t* cluster_base_ranges = create_cluster_ranges_from_groups(base_groups, (uint32_t)md_array_size(base_groups), arena);

    uint32_t* cluster_ranges = NULL;
    md_array_push_array(cluster_ranges, cluster_base_ranges, md_array_size(cluster_base_ranges), arena);

    // Add according to group and part of instances
    for (uint32_t i = 0; i < inst_count; ++i) {
        range_t* groups = 0;
        uint32_t beg = inst_grp_ranges[i].offset;
        uint32_t end = inst_grp_ranges[i].offset + inst_grp_ranges[i].count;

        for (uint32_t j = beg; j < end; ++j) {
            md_array_push(groups, s->groups[j], arena);
        }
        uint32_t* inst_ranges = create_cluster_ranges_from_groups(groups, (uint32_t)md_array_size(groups), arena);
        s->instances[i].cluster_range = (range_t){(uint32_t)md_array_size(cluster_ranges), (uint32_t)md_array_size(inst_ranges)};
        md_array_push_array(cluster_ranges, inst_ranges, md_array_size(inst_ranges), arena);
    }

    uint32_t cluster_count = (uint32_t)md_array_size(cluster_ranges);

    s->atom.count = count;
    s->cluster.count = cluster_count;

    gl_buffer_set_sub_data(ctx.cluster_data_atom_idx_buf, s->atom.offset * sizeof(uint32_t),     s->atom.count * sizeof(uint32_t),     src_indices);
    gl_buffer_set_sub_data(ctx.cluster_range_buf,         s->cluster.offset * sizeof(uint32_t),  s->cluster.count * sizeof(uint32_t),  cluster_ranges);

    md_vm_arena_destroy(arena);
}

// Great resource and explanation of some of the techiques used when building BVH trees
// https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/


void recompute_clusters(structure_t* s, const float* x, const float* y, const float* z, uint32_t count, uint32_t byte_stride) {
    md_allocator_i* alloc = md_get_heap_allocator();

    // Compute clusters
    // This is a tricky problem to solve in real-time
    // We go for a greedy approach where we spatially sort the positions using spatial binning in morton order.
    // Then we try to pick 'good' continous ranges which make up clusters.
    // The constraints we have for each cluster is that it should have a size which falls between MIN_COUNT_PER_CLUSTER and MAX_COUNT_PER_CLUSTER
    // This problem is similar to BVH computation, but we are disregarding the Hierarchy and going straight for the nodes which are equivalent to the clusters.
    // Resorting to Surface Area Heuristics, might be a good metric in order to determine the boundary between clusters.

    // The optimal solution albeit possibly costly would be to compute this using a high quality BVH split approach with residue groups as leafs.
    // With a stopping condition when a node contains less than 256 entries.

    // @TODO: Implement this in a compute shader

    uint32_t* cluster_ranges = NULL;
    uint32_t* src_indices = md_alloc(alloc, sizeof(uint32_t) * count);

    md_util_sort_spatial(src_indices, x, y, z, count);
    md_array_ensure(cluster_ranges, DIV_UP(count, MIN_COUNT_PER_CLUSTER), alloc);

    vec3_t min_xyz = {+FLT_MAX, +FLT_MAX, +FLT_MAX};
    vec3_t max_xyz = {-FLT_MAX, -FLT_MAX, -FLT_MAX};

    // Create the clusters from consecutive ranges
    // We basically start at the minimum size and keep on adding points until we get a substantial increase in the surface area of the bounding box
    for (uint32_t i = 0; i < count;) {
        uint32_t cluster_offset = i;
        uint32_t cluster_size   = MIN(MIN_COUNT_PER_CLUSTER, count - i);
        vec3_t cluster_aabb_min = {+FLT_MAX, +FLT_MAX, +FLT_MAX};
        vec3_t cluster_aabb_max = {-FLT_MAX, -FLT_MAX, -FLT_MAX};

        for (; i < cluster_offset + cluster_size; ++i) {
            uint32_t idx = src_indices[i];
            vec3_t p = (vec3_t){x[idx], y[idx], z[idx]};
            cluster_aabb_min = vec3_min(cluster_aabb_min, p);
            cluster_aabb_max = vec3_max(cluster_aabb_max, p);
        }

        if (i == count) break;

        vec3_t ext = vec3_sub(cluster_aabb_max, cluster_aabb_min);
        float area = 2.0f * (ext.x*ext.y + ext.x*ext.z + ext.y*ext.z);
        for (; i < MIN(cluster_offset + MAX_COUNT_PER_CLUSTER, count); ++i) {
            uint32_t idx = src_indices[i];
            vec3_t p = (vec3_t){x[idx], y[idx], z[idx]};

            vec3_t aabb_min = vec3_min(cluster_aabb_min, p);
            vec3_t aabb_max = vec3_max(cluster_aabb_max, p);
            vec3_t aabb_ext = vec3_sub(aabb_max, aabb_min);
            float  new_area = 2.0f * (aabb_ext.x*aabb_ext.y + aabb_ext.x*aabb_ext.z + aabb_ext.y*aabb_ext.z);

            // If adding the next point increases the surface area of the bounding box by more than X%, then it goes into a new cluster
            float ratio = new_area / area;
            if (ratio > 2.0f) {
                break;
            }

            // Commit
            cluster_size += 1;
            cluster_aabb_min = aabb_min;
            cluster_aabb_max = aabb_max;
            //area = new_area;
        }

        uint32_t cluster_range = (cluster_offset << 8) | cluster_size;
        md_array_push(cluster_ranges, cluster_range, alloc);

        min_xyz = vec3_min(min_xyz, cluster_aabb_min);
        max_xyz = vec3_max(max_xyz, cluster_aabb_max);
    }

    uint32_t cluster_count = (uint32_t)md_array_size(cluster_ranges);

    s->base_instance.min_xyzr = vec4_from_vec3(min_xyz, s->base_instance.min_xyzr.w);
    s->base_instance.min_xyzr = vec4_from_vec3(max_xyz, s->base_instance.max_xyzr.w);
    s->base_instance.cluster_range = (range_t){0, cluster_count};

    s->atom.count = count;
    s->cluster.count = cluster_count;

    gl_buffer_set_sub_data(ctx.cluster_data_atom_idx_buf, s->atom.offset * sizeof(uint32_t),     s->atom.count * sizeof(uint32_t),     src_indices);
    gl_buffer_set_sub_data(ctx.cluster_range_buf,         s->cluster.offset * sizeof(uint32_t),  s->cluster.count * sizeof(uint32_t),  cluster_ranges);

    md_array_free(cluster_ranges, alloc);
    md_free(alloc, src_indices, sizeof(uint32_t) * count);
}

void update_cluster_data(const structure_t* s) {
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, POSITION_BINDING,            ctx.position_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, RADIUS_BINDING,              ctx.radius_buf.id);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_RANGE_BINDING,       ctx.cluster_range_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_BOUNDS_BINDING,      ctx.cluster_bounds_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_POS_RAD_BINDING,     ctx.cluster_data_pos_rad_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_ATOM_INDEX_BINDING,  ctx.cluster_data_atom_idx_buf.id);

    glUseProgram(ctx.cluster_compute_data_prog.id);
    glDispatchCompute(s->cluster.count, 1, 1);
    glUseProgram(0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, POSITION_BINDING,            0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, RADIUS_BINDING,              0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_RANGE_BINDING,       0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_BOUNDS_BINDING,      0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_POS_RAD_BINDING,     0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_ATOM_INDEX_BINDING,  0);

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
}

bool md_gfx_structure_set_atom_position_soa(md_gfx_handle_t id, const float* x, const float* y, const float* z, uint32_t count) {
    structure_t* s = get_structure(id);
    if (!s) {
        MD_LOG_ERROR("Failed to set atom position data: Handle is invalid");
        return false;
    }

    if (!x || !y || ! z) {
        MD_LOG_ERROR("One or more arguments are missing, must pass x, y and z for position.");
        return false;
    }

    if (count > s->atom.capacity) {
        MD_LOG_ERROR("Attempting to write out of bounds");
        return false;
    }
    
    vec3_t* pos = glMapNamedBufferRange(ctx.position_buf.id, s->atom.offset * sizeof(vec3_t), s->atom.capacity * sizeof(vec3_t), GL_MAP_WRITE_BIT);
    if (!pos) {
        MD_LOG_ERROR("Failed to set atom position data: Could not map buffer");
        return false;
    }
   
    for (uint32_t i = 0; i < count; ++i) {
        pos[i].x = x[i];
        pos[i].y = y[i];
        pos[i].z = z[i];
    }
    glUnmapNamedBuffer(ctx.position_buf.id);

    // Recompute clusters based on new positional data
    recompute_clusters(s, x, y, z, count, 0);
    update_cluster_data(s);
    
    return true;
}

bool md_gfx_structure_set_atom_position(md_gfx_handle_t id, const vec3_t* xyz, uint32_t count, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        MD_LOG_ERROR("Failed to set atom position data: Handle is invalid");
        return false;
    }

    if (!xyz) {
        MD_LOG_ERROR("Argument is missing");
        return false;
    }

    if (count > s->atom.capacity) {
        MD_LOG_ERROR("Attempting to write out of bounds");
        return false;
    }

    byte_stride = MAX(sizeof(vec3_t), byte_stride);
    if (byte_stride > sizeof(vec3_t)) {
        vec3_t* pos = glMapNamedBufferRange(ctx.position_buf.id, s->atom.offset * sizeof(vec3_t), s->atom.capacity * sizeof(vec3_t), GL_MAP_READ_BIT | GL_MAP_WRITE_BIT);
        if (!pos) {
            MD_LOG_ERROR("Failed to set atom position data: Could not map buffer");
            return false;
        }

        for (uint32_t i = 0; i < count; ++i) {
            pos[i] = *(const vec3_t*)((const uint8_t*)xyz + byte_stride * i);
        }
        glUnmapNamedBuffer(ctx.position_buf.id);
    } else {
        glNamedBufferSubData(ctx.position_buf.id, s->atom.offset * sizeof(vec3_t), s->atom.capacity * sizeof(vec3_t), xyz);
    }

    // Recompute clusters based on new positional data
    recompute_clusters(s, &xyz->x, &xyz->y, &xyz->z, count, byte_stride);
    update_cluster_data(s);

    return true;
}

bool md_gfx_structure_set_atom_radius(md_gfx_handle_t id, const float* radius, uint32_t count, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        MD_LOG_ERROR("Failed to set atom radius data: Handle is invalid");
        return false;
    }

    if (!radius) {
        MD_LOG_ERROR("Failed to set atom radius data: Radius argument is missing.");
        return false;
    }

    if (count > s->atom.capacity) {
        MD_LOG_ERROR("Attempting to write out of bounds");
        return false;
    }

    const uint32_t element_size = sizeof(float);
    if (byte_stride > element_size) {
        float* ptr = glMapNamedBufferRange(ctx.radius_buf.id, s->atom.offset * element_size, s->atom.capacity * element_size, GL_MAP_WRITE_BIT);
        if (!ptr) {
            MD_LOG_ERROR("Failed to set atom radius data: Could not map buffer");
            return false;
        }
        byte_stride = MAX(element_size, byte_stride);
        for (uint32_t i = 0; i < count; ++i) {
            ptr[i] = *(const float*)((const uint8_t*)radius + byte_stride * i);
        }
        glUnmapNamedBuffer(ctx.radius_buf.id);
    } else {
        gl_buffer_set_sub_data(ctx.radius_buf, s->atom.offset * element_size, count * element_size, radius);
    }

    float min_r = +FLT_MAX;
    float max_r = -FLT_MAX;
    for (uint32_t i = 0; i < count; ++i) {
        min_r = MIN(min_r, radius[i]);
        max_r = MAX(max_r, radius[i]);
    }

    s->base_instance.min_xyzr.w = min_r;
    s->base_instance.max_xyzr.w = max_r;

    update_cluster_data(s);

    return true;
}

bool md_gfx_structure_set_atom_flags(md_gfx_handle_t id, const uint32_t* flags, uint32_t count, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        MD_LOG_ERROR("Failed to set atom flags data: Handle is invalid");
        return false;
    }

    if (!flags) {
        MD_LOG_ERROR("Failed to set atom flags data: Flags argument is missing.");
        return false;
    }

    if (count > s->atom.capacity) {
        MD_LOG_ERROR("Attempting to write out of bounds");
        return false;
    }

    const uint32_t element_size = sizeof(uint32_t);
    if (byte_stride > element_size) {
        uint32_t* ptr = glMapNamedBufferRange(ctx.flags_buf.id, s->atom.offset * element_size, s->atom.capacity * element_size, GL_MAP_WRITE_BIT);
        if (!ptr) {
            MD_LOG_ERROR("Failed to set atom flags data: Could not map buffer");
            return false;
        }
        byte_stride = MAX(element_size, byte_stride);
        for (uint32_t i = 0; i < count; ++i) {
            ptr[i] = *(const uint32_t*)((const uint8_t*)flags + byte_stride * i);
        }
        glUnmapNamedBuffer(ctx.flags_buf.id);
    } else {
        gl_buffer_set_sub_data(ctx.flags_buf, s->atom.offset * element_size, count * element_size, flags);
    }
    return true;
}

bool md_gfx_structure_set_group_atom_ranges(md_gfx_handle_t id, const md_gfx_range_t* atom_ranges, uint32_t count, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        MD_LOG_ERROR("Failed to set group atom range data: Handle is invalid");
        return false;
    }

    if (!atom_ranges) {
        MD_LOG_ERROR("Failed to set group atom range data: Argument is missing.");
        return false;
    }

    if (count > md_array_capacity(s->groups)) {
        MD_LOG_ERROR("Attempting to write out of bounds");
        return false;
    }

    byte_stride = MAX(byte_stride, sizeof(md_gfx_range_t));
    for (uint32_t i = 0; i < count; ++i) {
        md_gfx_range_t range = *(const md_gfx_range_t*)((const uint8_t*)atom_ranges + byte_stride * i);
        s->groups[i].offset = range.beg_idx;
        s->groups[i].count  = range.end_idx - range.beg_idx;
    }

    return true;
}


bool md_gfx_structure_set_instance_atom_ranges(md_gfx_handle_t id, const md_gfx_range_t* atom_ranges, uint32_t count, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        MD_LOG_ERROR("Failed to set instance atom range data: Handle is invalid");
        return false;
    }

    if (!atom_ranges) {
        MD_LOG_ERROR("Failed to set instance atom range data: Argument is missing.");
        return false;
    }

    if (count > md_array_capacity(s->instances)) {
        MD_LOG_ERROR("Attempting to write out of bounds");
        return false;
    }

    byte_stride = MAX(byte_stride, sizeof(md_gfx_range_t));
    for (uint32_t i = 0; i < count; ++i) {
        md_gfx_range_t range = *(const md_gfx_range_t*)((const uint8_t*)atom_ranges + byte_stride * i);
        s->instances[i].atom_range.offset = range.beg_idx;
        s->instances[i].atom_range.count  = range.end_idx - range.beg_idx;
    }

    return true;
}

bool md_gfx_structure_set_instance_transforms(md_gfx_handle_t id, const struct mat4_t* transforms, uint32_t count, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        MD_LOG_ERROR("Failed to set instance group range data: Handle is invalid");
        return false;
    }

    if (!transforms) {
        MD_LOG_ERROR("Failed to set instance group range data: Argument is missing.");
        return false;
    }

    if (count > md_array_capacity(s->instances)) {
        MD_LOG_ERROR("Attempting to write out of bounds");
        return false;
    }

    md_array_resize(s->instances, count, md_get_heap_allocator());

    byte_stride = MAX(byte_stride, sizeof(mat4_t));
    for (uint32_t i = 0; i < count; ++i) {
        mat4_t mat = *(const mat4_t*)((const uint8_t*)transforms + byte_stride * i);
        s->instances[i].transform = mat;
    }

    return true;
}

md_gfx_handle_t md_gfx_rep_create(uint32_t color_count) {
    if (!validate_context()) {
        return (md_gfx_handle_t){0};
    }

    if (ctx.representation_count >= MAX_REPRESENTATION_COUNT) {
        MD_LOG_ERROR("Unable to create representation, pool is full!");
        return (md_gfx_handle_t){0};
    }
    uint16_t r_idx = ctx.representation_count++;

    representation_t* r = &ctx.representations[r_idx];
    r->color = field_alloc(&ctx.color, color_count); 

    r->h_idx = alloc_handle();

    // Map handle index to resource index
    ctx.handles[r->h_idx].idx = r_idx;

    // Create external handle for user
    return (md_gfx_handle_t) {ctx.handles[r->h_idx].gen, r->h_idx};
}

bool md_gfx_rep_destroy(md_gfx_handle_t id) {
    if (!validate_handle(id)) {
        MD_LOG_ERROR("Failed to destroy representation: Handle is invalid");
        return false;
    }
    uint16_t h_idx = id.idx;
    uint16_t r_idx = ctx.handles[h_idx].idx;
    free_handle(h_idx);

    ASSERT(ctx.representation_count > 0);
    ctx.representation_count -= 1;

    representation_t* r = &ctx.representations[r_idx];
    field_free(&ctx.color, r->color);

    if (ctx.representation_count > 1) {
        // Swap back and pop resource array to keep it packed
        ctx.representations[r_idx] = ctx.representations[ctx.representation_count];

        // Reassign handle indices to new location
        ctx.handles[ctx.representations[r_idx].h_idx].idx = r_idx;
    }

    return true;
}

bool md_gfx_rep_set_data(md_gfx_handle_t id, md_gfx_rep_type_t type, md_gfx_rep_attr_t attr, const md_gfx_color_t* color, uint32_t count, uint32_t byte_stride) {
    representation_t* r = get_representation(id);
    if (!r) {
        MD_LOG_ERROR("Failed to set representation data: Handle is invalid");
        return false;
    }

    if (!color) {
        MD_LOG_ERROR("Failed to set color data: Color argument is missing.");
        return false;
    }

    if (count > r->color.capacity) {
        MD_LOG_ERROR("Attempting to write out of bounds.");
        return false;
    }

    const uint32_t element_size = sizeof(md_gfx_color_t);
    if (byte_stride > element_size) {
        md_gfx_color_t* ptr = glMapNamedBufferRange(ctx.color_buf.id, r->color.offset * element_size, r->color.capacity * element_size, GL_WRITE_ONLY);
        if (!ptr) {
            MD_LOG_ERROR("Failed to set color data: Could not map buffer");
            return false;
        }
        byte_stride = MAX(element_size, byte_stride);
        for (uint32_t i = 0; i < count; ++i) {
            ptr[i] = *(const md_gfx_color_t*)((const uint8_t*)color + byte_stride * i);
        }
        glUnmapNamedBuffer(ctx.color_buf.id);
    } else {
        gl_buffer_set_sub_data(ctx.color_buf, r->color.offset * element_size, count * element_size, color);
    }

    r->color.count = count;

    r->type = type;
    r->attr = attr;
    return true;
}

bool md_gfx_draw(uint32_t in_draw_op_count, const md_gfx_draw_op_t* in_draw_ops, const mat4_t* proj_mat, const mat4_t* view_mat, const mat4_t* inv_proj_mat, const mat4_t* inv_view_mat) {
    if (!validate_context()) return false;

    if (in_draw_op_count > 0 && !in_draw_ops) {
        MD_LOG_ERROR("Failed to Draw: draw ops was NULL");
        return false;
    }

    if (!proj_mat || !view_mat || !inv_proj_mat || !inv_view_mat) {
        MD_LOG_ERROR("Failed to Draw: one or more missing matrix arguments");
        return false;
    }

    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    // Store old state
    uint32_t  old_draw_buffers[16];
    uint32_t  old_viewport[4];
    uint32_t  old_draw_buffer_count;
    uint32_t  old_fbo;
    uint32_t  old_vao;
    uint32_t  old_cull_mode;
    uint32_t  old_depth_func;
    float     old_clear_color[4];
    float     old_clear_depth;
    bool      old_color_mask[4];
    bool      old_depth_mask;
    bool      old_depth_test;
    bool      old_cull_face;

    glGetBooleanv(GL_COLOR_WRITEMASK, (GLboolean*)(old_color_mask));
    glGetBooleanv(GL_DEPTH_TEST, (GLboolean*)(&old_depth_test));
    glGetBooleanv(GL_DEPTH_WRITEMASK, (GLboolean*)(&old_depth_mask));
    glGetBooleanv(GL_CULL_FACE, (GLboolean*)(&old_cull_face));

    glGetIntegerv(GL_DEPTH_FUNC, (GLint*)(&old_depth_func));
    glGetIntegerv(GL_CULL_FACE_MODE, (GLint*)(&old_cull_mode));
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, (GLint*)(&old_fbo));
    glGetIntegerv(GL_VERTEX_ARRAY_BINDING, (GLint*)(&old_vao));
    glGetIntegerv(GL_VIEWPORT, (GLint*)old_viewport);
    glGetFloatv(GL_COLOR_CLEAR_VALUE, old_clear_color);
    glGetFloatv(GL_DEPTH_CLEAR_VALUE, &old_clear_depth);

    for (old_draw_buffer_count = 0; old_draw_buffer_count < ARRAY_SIZE(old_draw_buffers); ++old_draw_buffer_count) {
        glGetIntegerv(GL_DRAW_BUFFER0 + old_draw_buffer_count, (GLint*)(&old_draw_buffers[old_draw_buffer_count]));
        if (old_draw_buffers[old_draw_buffer_count] == 0) {
            break;
        }
    }

    // Set state
    const uint32_t draw_buffers[] = { GL_COLOR_ATTACHMENT0 + COLOR_TEX_ATTACHMENT, GL_COLOR_ATTACHMENT0 + NORMAL_TEX_ATTACHMENT, GL_COLOR_ATTACHMENT0 + SS_VEL_TEX_ATTACHMENT, GL_COLOR_ATTACHMENT0 + INDEX_TEX_ATTACHMENT };
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ctx.fbo_id);
    glDrawBuffers(ARRAY_SIZE(draw_buffers), draw_buffers);
    glViewport(0, 0, ctx.fbo_width, ctx.fbo_height);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);    
    glDepthFunc(GL_LESS);
    glDepthMask(1);
    glColorMask(1, 1, 1, 1);

    mat4_t curr_transform = mat4_mul(*proj_mat, *view_mat);
    mat4_t prev_transform = mat4_mul(ctx.prev_proj_mat, ctx.prev_view_mat);
    mat4_t inv_transform = mat4_mul(*inv_view_mat, *inv_proj_mat);
    mat4_t curr_clip_to_prev_clip_mat = mat4_mul(prev_transform, inv_transform);

    UniformData ubo_data = {
        .world_to_view = *view_mat,
        .world_to_view_normal = mat4_transpose(*inv_view_mat),
        .view_to_clip = *proj_mat,
        .clip_to_view = *inv_proj_mat,
        .world_to_clip = curr_transform,
        .curr_clip_to_prev_clip = curr_clip_to_prev_clip_mat,
        .prev_world_to_clip = prev_transform,
        .prev_world_to_view = ctx.prev_view_mat,
        .prev_view_to_clip = ctx.prev_proj_mat,
        .frustum_culling = 1,
        .depth_culling = 1,
        .depth_pyramid_width = ctx.depth_pyramid_tex.width,
        .depth_pyramid_height = ctx.depth_pyramid_tex.height,
        .render_width = ctx.fbo_width,
        .render_height = ctx.fbo_height,
        .znear = extract_znear(*inv_proj_mat),
        .zfar = extract_zfar(*inv_proj_mat),
        .late = 0,
    };

    ctx.prev_proj_mat = *proj_mat;
    ctx.prev_view_mat = *view_mat;

    gl_buffer_set_sub_data(ctx.ubo_buf, 0, sizeof(ubo_data), &ubo_data);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ctx.ubo_buf.id);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, POSITION_BINDING,                ctx.position_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, RADIUS_BINDING,                  ctx.radius_buf.id);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_RANGE_BINDING,           ctx.cluster_range_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_BOUNDS_BINDING,          ctx.cluster_bounds_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_POS_RAD_BINDING,         ctx.cluster_data_pos_rad_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_ATOM_INDEX_BINDING,      ctx.cluster_data_atom_idx_buf.id);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, COLOR_BINDING,                   ctx.color_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, TRANSFORM_BINDING,               ctx.transform_buf.id);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_INDIRECT_PARAM_BINDING,     ctx.draw_ind_param_buf.id);
    //glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_INDIRECT_CMD_BINDING,     ctx.draw_ind_cmd_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_SPHERE_INDEX_BINDING,       ctx.draw_sphere_idx_buf.id);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_INST_DATA_BINDING,       ctx.cluster_instance_data_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_INST_DRAW_IDX_BINDING,   ctx.cluster_instance_vis_idx_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_INST_RAST_IDX_BINDING,   ctx.cluster_instance_rast_idx_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_INST_OCC_IDX_BINDING,    ctx.cluster_instance_occ_idx_buf.id);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, INSTANCE_DATA_BINDING,           ctx.instance_data_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, INSTANCE_VIS_IDX_BINDING,        ctx.instance_vis_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, INSTANCE_OCC_IDX_BINDING,        ctx.instance_occ_buf.id);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DEBUG_BINDING,                   ctx.DEBUG_BUF.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, VIS_BINDING,                     ctx.visibility_buf.id);

    glBindBuffer(GL_DISPATCH_INDIRECT_BUFFER,   ctx.draw_ind_param_buf.id);
    glBindBuffer(GL_DRAW_INDIRECT_BUFFER,       ctx.draw_ind_param_buf.id);

    glBindTextureUnit(0, ctx.depth_pyramid_tex.id);
    
    InstanceData* inst_data = 0;

    mat4_t* transforms = 0;
    md_array_push(transforms, mat4_ident(), arena);

    for (uint32_t in_idx = 0; in_idx < in_draw_op_count; ++in_idx) {
        const md_gfx_draw_op_t* in_op = in_draw_ops + in_idx;

        // Validate draw op
        structure_t* s = get_structure(in_op->structure);
        if (!s) {
            MD_LOG_ERROR("Invalid structure supplied in draw operation");
            continue;   
        }
        representation_t* r = get_representation(in_op->representation);
        if (!r) {
            MD_LOG_ERROR("Invalid representation supplied in draw operation");
            continue;
        }

        if (s->atom.count != r->color.count) {
            MD_LOG_ERROR("Number of colors in representation does not match number of atoms in structure");
            continue;
        }

        // Create draw op for the base
        {
            uint32_t transform_idx = 0; // This is the default and represents identity matrix
            if (in_op->model_mat && !mat4_equal(*in_op->model_mat, mat4_ident())) {
                transform_idx = (uint32_t)md_array_size(transforms);
                md_array_push(transforms, *in_op->model_mat, arena);
            }

            InstanceData out = {0};
            out.min_xyzr = s->base_instance.min_xyzr;
            out.max_xyzr = s->base_instance.max_xyzr;
            out.cluster_offset = s->base_instance.cluster_range.offset;
            out.cluster_count  = s->base_instance.cluster_range.count;
            out.transform_idx = ctx.transform.offset + transform_idx;
            md_array_push(inst_data, out, arena);
        }

        // Create draw ops for instanced clusters
        for (uint32_t i = 0; i < (uint32_t)md_array_size(s->instances); ++i) {
            const instance_t* inst = s->instances + i;
            mat4_t M = inst->transform;
            if (in_op->model_mat && !mat4_equal(*in_op->model_mat, mat4_ident())) {
                M = mat4_mul(*in_op->model_mat, M);
            }
            uint32_t transform_idx = (uint32_t)md_array_size(transforms);
            md_array_push(transforms, M, arena);

            InstanceData out = {0};
            out.min_xyzr = inst->min_xyzr;
            out.max_xyzr = inst->max_xyzr;
            out.cluster_offset = inst->cluster_range.offset;
            out.cluster_count  = inst->cluster_range.count;
            out.transform_idx = ctx.transform.offset + transform_idx;
            md_array_push(inst_data, out, arena);
        }
    }

    uint32_t instance_count  = (uint32_t)md_array_size(inst_data);
    uint32_t transform_count = (uint32_t)md_array_size(transforms);

{
    GL_PUSH_GPU_SECTION("CLEAR BUFFERS");

    DispatchIndirectCommand     disp_ind_cmd = { .num_groups_x = 0, .num_groups_y = 1, .num_groups_z = 1};
    DrawElementsIndirectCommand draw_ind_cmd = { .instance_count = 1 };

    // this is what we clear with
    struct DrawParameters draw_param = {0};
    draw_param.write_vis_elem_cmd       = disp_ind_cmd;
    draw_param.inst_clust_cull_cmd      = disp_ind_cmd;
    draw_param.clust_raster_cmd         = disp_ind_cmd;
    draw_param.draw_sphere_cmd          = draw_ind_cmd;

    draw_param.inst_cull_late_cmd       = disp_ind_cmd;
    draw_param.clust_cull_late_cmd      = disp_ind_cmd;

    draw_param.instance_count           = instance_count;

    const uint64_t vis_clear_val = (0x3f800000ULL << 32) | 0;   // Hex value of float (1.0f)
    glClearNamedBufferData(ctx.visibility_buf.id,     GL_RG32UI, GL_RG,  GL_UNSIGNED_INT, &vis_clear_val);
    glNamedBufferSubData(ctx.draw_ind_param_buf.id, 0, sizeof(DrawParameters), &draw_param);

    MEMSET(ctx.debug_data, 0, sizeof(DebugData));

    GL_POP_GPU_SECTION();
}

{
    GL_PUSH_GPU_SECTION("UPLOAD INSTANCE DATA + TRANSFORMS");
    gl_buffer_set_sub_data(ctx.instance_data_buf, 0, sizeof(InstanceData) * instance_count, inst_data);
    gl_buffer_set_sub_data(ctx.transform_buf, ctx.transform.offset * sizeof(mat4_t), sizeof(mat4_t) * transform_count, transforms);
    GL_POP_GPU_SECTION();
}

    glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);

{
    GL_PUSH_GPU_SECTION("INSTANCE CULL");
    glUseProgram(ctx.instance_cull_prog.id);
    uint wg_size = DIV_UP(instance_count, INSTANCE_CULL_GROUP_SIZE);
    glDispatchCompute(wg_size, 1, 1);
    GL_POP_GPU_SECTION();
}

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);

{
    GL_PUSH_GPU_SECTION("INSTANCE CLUSTER CULL");
    glUseProgram(ctx.instance_cluster_cull_prog.id);
    glDispatchComputeIndirect((GLintptr)offsetof(DrawParameters, inst_clust_cull_cmd));
    GL_POP_GPU_SECTION();
}

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);

{
    GL_PUSH_GPU_SECTION("WRITE VISIBLE ELEMENT INDICES")
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_WRITE_ELEM_BINDING, ctx.cluster_instance_vis_idx_buf.id);
    glUseProgram(ctx.cluster_write_element_indices_prog.id);
    glDispatchComputeIndirect((GLintptr)offsetof(DrawParameters, write_vis_elem_cmd));
    GL_POP_GPU_SECTION();
}

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);
 
{
    GL_PUSH_GPU_SECTION("SPACEFILL");
    glUseProgram(ctx.spacefill_prog.id);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ctx.draw_sphere_idx_buf.id);
    glDrawElementsIndirect(GL_POINTS, GL_UNSIGNED_INT, (const void*)offsetof(DrawParameters, draw_sphere_cmd));
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    GL_POP_GPU_SECTION();
}
{
    GL_PUSH_GPU_SECTION("CLUSTER RASTER");
    glUseProgram(ctx.cluster_raster_prog.id);
    glDispatchComputeIndirect((GLintptr)offsetof(DrawParameters, clust_raster_cmd));
    GL_POP_GPU_SECTION();
}

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    // Reset the some counters for next round of retesting previously occluded clusters
    {
        const uint zero = 0;
        const uint one  = 1;
        glNamedBufferSubData(ctx.draw_ind_param_buf.id, offsetof(DrawParameters, write_vis_elem_cmd.num_groups_x),  sizeof(uint), &zero);
        glNamedBufferSubData(ctx.draw_ind_param_buf.id, offsetof(DrawParameters, clust_raster_cmd.num_groups_x),    sizeof(uint), &zero);
        glNamedBufferSubData(ctx.draw_ind_param_buf.id, offsetof(DrawParameters, inst_clust_cull_cmd.num_groups_x), sizeof(uint), &zero);
        glNamedBufferSubData(ctx.draw_ind_param_buf.id, offsetof(DrawParameters, draw_sphere_cmd.count),            sizeof(uint), &zero);
        glNamedBufferSubData(ctx.draw_ind_param_buf.id, offsetof(DrawParameters, cluster_draw_count),               sizeof(uint), &zero);
        glNamedBufferSubData(ctx.draw_ind_param_buf.id, offsetof(DrawParameters, cluster_rast_count),               sizeof(uint), &zero);
        glNamedBufferSubData(ctx.draw_ind_param_buf.id, offsetof(DrawParameters, instance_vis_count),               sizeof(uint), &zero);
        glNamedBufferSubData(ctx.ubo_buf.id, offsetof(UniformData, late), sizeof(uint), &one);
    }
    
{
    GL_PUSH_GPU_SECTION("DEPTH REDUCE");
    glUseProgram(ctx.vis_reduce_prog.id);
    glBindImageTexture(0, ctx.depth_pyramid_tex.id, 0, 0, 0, GL_WRITE_ONLY, GL_R32F);
    glDispatchCompute(DIV_UP(ctx.fbo_width, DEPTH_REDUCE_GROUP_SIZE), DIV_UP(ctx.fbo_height, DEPTH_REDUCE_GROUP_SIZE), 1);
    glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

    glUseProgram(ctx.depth_reduce_prog.id);
    glBindTextureUnit(1, ctx.depth_pyramid_tex.id);
    for (uint32_t i = 1; i < ctx.depth_pyramid_tex.levels; ++i) {
        // Output
        uint width  = MAX(1, ctx.depth_pyramid_tex.width >> i);
        uint height = MAX(1, ctx.depth_pyramid_tex.height >> i);
        float lod = (float)(i-1);
        glUniform2ui(0, width, height);
        glUniform1f(1, lod);
        glBindImageTexture(0, ctx.depth_pyramid_tex.id, i, 0, 0, GL_WRITE_ONLY, GL_R32F);
        glDispatchCompute(DIV_UP(width, DEPTH_REDUCE_GROUP_SIZE), DIV_UP(height, DEPTH_REDUCE_GROUP_SIZE), 1);
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    }
    GL_POP_GPU_SECTION();
}

glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);

{
    GL_PUSH_GPU_SECTION("INSTANCE CULL LATE");
    glUseProgram(ctx.instance_cull_late_prog.id);
    glDispatchComputeIndirect((GLintptr)offsetof(DrawParameters, inst_cull_late_cmd));
    GL_POP_GPU_SECTION();
}

glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);

{
    GL_PUSH_GPU_SECTION("INSTANCE CLUSTER CULL LATE");
    glUseProgram(ctx.instance_cluster_cull_prog.id);
    glDispatchComputeIndirect((GLintptr)offsetof(DrawParameters, inst_clust_cull_cmd));
    GL_POP_GPU_SECTION();
}

{
    GL_PUSH_GPU_SECTION("CLUSTER CULL LATE");
    glUseProgram(ctx.cluster_cull_late_prog.id);
    glDispatchComputeIndirect((GLintptr)offsetof(DrawParameters, clust_cull_late_cmd));
    GL_POP_GPU_SECTION();
}

glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);

{
    GL_PUSH_GPU_SECTION("WRITE VISIBLE ELEMENT INDICES LATE")
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_WRITE_ELEM_BINDING, ctx.cluster_instance_vis_idx_buf.id);
    glUseProgram(ctx.cluster_write_element_indices_prog.id);
    glDispatchComputeIndirect((GLintptr)offsetof(DrawParameters, write_vis_elem_cmd));
    GL_POP_GPU_SECTION();
}

glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);

{
    GL_PUSH_GPU_SECTION("SPACEFILL LATE");
    glUseProgram(ctx.spacefill_prog.id);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ctx.draw_sphere_idx_buf.id);
    glDrawElementsIndirect(GL_POINTS, GL_UNSIGNED_INT, (const void*)offsetof(DrawParameters, draw_sphere_cmd));
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    GL_POP_GPU_SECTION();
}

{
    GL_PUSH_GPU_SECTION("CLUSTER RASTER LATE");
    glUseProgram(ctx.cluster_raster_prog.id);
    glDispatchComputeIndirect((GLintptr)offsetof(DrawParameters, clust_raster_cmd));
    GL_POP_GPU_SECTION();
}

glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);


    // DRAW REST
    // TODO

    // COMPOSE & OUTPUT
GL_PUSH_GPU_SECTION("COMPOSE");
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, old_fbo);
    glDrawBuffers(old_draw_buffer_count, old_draw_buffers);
    glViewport(old_viewport[0], old_viewport[1], old_viewport[2], old_viewport[3]);

    glUseProgram(ctx.compose_prog.id);
    glBindTextureUnit(DEPTH_TEX_BINDING,  ctx.depth_tex.id);
    glBindTextureUnit(COLOR_TEX_BINDING,  ctx.color_tex.id);
    glBindTextureUnit(NORMAL_TEX_BINDING, ctx.normal_tex.id);
    glBindTextureUnit(SS_VEL_TEX_BINDING, ctx.ss_vel_tex.id);
    glBindTextureUnit(INDEX_TEX_BINDING,  ctx.index_tex.id);

    glDrawArrays(GL_TRIANGLES, 0, 3);
GL_POP_GPU_SECTION();

{
    // COMPUTE DEPTH PYRAMID LATE
    GL_PUSH_GPU_SECTION("DEPTH REDUCE LATE");
    glUseProgram(ctx.vis_reduce_prog.id);
    glBindImageTexture(0, ctx.depth_pyramid_tex.id, 0, 0, 0, GL_WRITE_ONLY, GL_R32F);
    glDispatchCompute(DIV_UP(ctx.fbo_width, DEPTH_REDUCE_GROUP_SIZE), DIV_UP(ctx.fbo_height, DEPTH_REDUCE_GROUP_SIZE), 1);
    glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

    glUseProgram(ctx.depth_reduce_prog.id);
    glBindTextureUnit(1, ctx.depth_pyramid_tex.id);
    for (uint32_t i = 1; i < ctx.depth_pyramid_tex.levels; ++i) {
        // Output
        uint width  = MAX(1, ctx.depth_pyramid_tex.width >> i);
        uint height = MAX(1, ctx.depth_pyramid_tex.height >> i);
        float lod = (float)(i-1);
        glUniform2ui(0, width, height);
        glUniform1f(1, lod);
        glBindImageTexture(0, ctx.depth_pyramid_tex.id, i, 0, 0, GL_WRITE_ONLY, GL_R32F);
        glDispatchCompute(DIV_UP(width, DEPTH_REDUCE_GROUP_SIZE), DIV_UP(height, DEPTH_REDUCE_GROUP_SIZE), 1);
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    }
    GL_POP_GPU_SECTION();
}

//md_logf(MD_LOG_TYPE_DEBUG, "Instance late count %i, cluster late count %i", ctx.debug_data->instance_late_count, ctx.debug_data->instance_clust_late_count);

GL_PUSH_GPU_SECTION("RESET STATE");
    glUseProgram(0);
    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0);
    glBindBuffer(GL_DISPATCH_INDIRECT_BUFFER, 0);
    glBindBuffer(GL_PARAMETER_BUFFER, 0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, POSITION_BINDING,            0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, RADIUS_BINDING,              0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_RANGE_BINDING,       0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_BOUNDS_BINDING,      0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_POS_RAD_BINDING,     0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_ATOM_INDEX_BINDING,  0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, COLOR_BINDING,               0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, TRANSFORM_BINDING,           0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, INSTANCE_DATA_BINDING,       0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, GPU_DRAW_OP_BINDING,         0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_INDIRECT_PARAM_BINDING, 0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_INDIRECT_CMD_BINDING,   0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_SPHERE_INDEX_BINDING,   0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_INST_DATA_BINDING,       0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_INST_DRAW_IDX_BINDING,   0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_INST_RAST_IDX_BINDING,   0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_INST_OCC_IDX_BINDING,    0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CLUSTER_WRITE_ELEM_BINDING,      0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DEBUG_BINDING,               0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, VIS_BINDING,                 0);

    // reset state
    glCullFace(old_cull_mode);    
    glDepthFunc(old_depth_func);
    glDepthMask(old_depth_mask);
    glColorMask(old_color_mask[0], old_color_mask[1], old_color_mask[2], old_color_mask[3]);
    glBindVertexArray(old_vao);
    
    if (old_cull_face) {
        glEnable(GL_CULL_FACE);
    } else {
        glDisable(GL_CULL_FACE);
    }

    if (old_depth_test) {
        glEnable(GL_DEPTH_TEST);
    } else {
        glDisable(GL_DEPTH_TEST);
    }
GL_POP_GPU_SECTION();

    md_vm_arena_destroy(arena);

    return true;
}

bool md_gfx_query_picking(uint32_t mouse_x, uint32_t mouse_y) {
    if (!validate_context()) return false;

GL_PUSH_GPU_SECTION("WRITE PICKING DATA");
    glUseProgram(ctx.write_picking_prog.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DEBUG_BINDING, ctx.DEBUG_BUF.id);
    glBindTextureUnit(DEPTH_TEX_BINDING,  ctx.depth_tex.id);
    glBindImageTexture(0, ctx.index_tex.id, 0, GL_FALSE, 0, GL_READ_ONLY, GL_R32UI);
    glUniform2i (0, (int)mouse_x, (int)mouse_y);
    glUniform2f (1, (float)ctx.index_tex.width, (float)ctx.index_tex.height);
    glDispatchCompute(1,1,1);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DEBUG_BINDING, 0);
    glUseProgram(0);
GL_POP_GPU_SECTION();

    return true;
}

uint32_t md_gfx_get_picking_idx() {
    if (!validate_context()) return INVALID_INDEX;
    uint32_t idx = ctx.debug_data->pick_index;
    return idx == 0 ? INVALID_INDEX : idx - 1;
}

float md_gfx_get_picking_depth() {
    if (!validate_context()) return 1.0f;
    return ctx.debug_data->pick_depth;
}
