#include "md_gfx.h"

#include "md_util.h"
#include <GL/gl3w.h>

#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_file.h>
#include <core/md_vec_math.h>
#include <core/md_compiler.h>
#include <core/md_allocator.h>
#include <core/md_array.inl>

#include <string.h>     // memset, memcpy, strstr
#include <stdlib.h>     // malloc
#include <stddef.h>

typedef mat4_t mat4;
typedef vec4_t vec4;
typedef uint32_t uint;
typedef struct UniformData UniformData;
typedef struct DebugData DebugData;

typedef struct DrawArraysIndirectCommand DrawArraysIndirectCommand;
typedef struct DrawElementsIndirectCommand DrawElementsIndirectCommand;
typedef struct DispatchIndirectCommand DispatchIndirectCommand;
typedef struct DrawSpheresIndirectCommand DrawSpheresIndirectCommand;
typedef struct DrawIndirect DrawIndirect;
typedef struct DrawOp DrawOp;

// Common shared header between
#include <shaders/gfx/common.h>

// For subdividing backbone segments
#define MIN_SPLINE_SUBDIVISION_COUNT 1
#define MAX_SPLINE_SUBDIVISION_COUNT 32

#define MAX_STRUCTURE_COUNT 64
#define MAX_REPRESENTATION_COUNT 64
#define MAX_HANDLE_COUNT (MAX_STRUCTURE_COUNT + MAX_REPRESENTATION_COUNT)

#define MAX_INSTANCE_PER_STRUCTURE_COUNT 32

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
    vec2_t jitter;
    if (is_ortho_proj_matrix(P)) {
        jitter.x = -P.elem[3][0] * 0.5f;
        jitter.y = -P.elem[3][1] * 0.5f;
    }
    else {
        jitter.x = P.elem[2][0] * 0.5f;
        jitter.y = P.elem[2][1] * 0.5f;
    }
    return jitter;
}

static inline uint32_t compute_mip_count(uint32_t width, uint32_t height) {
    uint32_t mips = 1;
    uint32_t max_dim = MAX(width, height);
    while (max_dim >>= 1) ++mips;
    return mips;
}

typedef struct gl_draw_elements_indirect_command_t {
    uint32_t count;
    uint32_t instance_count;
    uint32_t first_index;
    uint32_t base_vertex;
    uint32_t base_instance;
} gl_draw_elements_indirect_command_t;

typedef struct gl_dispatch_indirect_command_t {
    uint32_t num_groups_x;
    uint32_t num_groups_y;
    uint32_t num_groups_z;
} gl_dispatch_indirect_command_t;

// This is the structure where we record keep track of indirect command counts
// Which feed the indirect rendering
typedef struct gl_command_buf_t {
    uint32_t sphere_drawcall_count;
    gl_dispatch_indirect_command_t sphere_software_raster_command;

} gl_command_buf_t;

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
    mat4_t   transform;
    range_t  group_range;
} instance_t;

typedef struct structure_t {
    uint32_t atom_offset;
    uint32_t atom_count;

    uint32_t bond_offset;
    uint32_t bond_count;

    uint32_t bb_seg_offset;
    uint32_t bb_seg_count;

    uint32_t bb_range_offset;
    uint32_t bb_range_count;

    uint32_t group_offset;
    uint32_t group_count;

    /*
    uint32_t instance_offset;
    uint32_t instance_count;
    */
    /*
    uint32_t transform_offset;
    uint32_t transform_count;
    */

    uint32_t instance_count;
    instance_t instance[MAX_INSTANCE_PER_STRUCTURE_COUNT];

    uint16_t h_idx; // Handle index
} structure_t;

typedef struct representation_t {
    md_gfx_rep_type_t type;
    md_gfx_rep_attr_t attr;

    uint32_t color_offset;
    uint32_t color_count;

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
    gl_texture_t visibility_tex;    // R64UI  64-BITS hi: 32-bit depth + lo: 32-bit payload, for compute rasterization of spheres

    // Structure Buffers
    gl_buffer_t position_buf;
    gl_buffer_t velocity_buf;
    gl_buffer_t radius_buf;
    gl_buffer_t bond_buf;
    gl_buffer_t segment_buf;
    gl_buffer_t sec_str_buf;
    gl_buffer_t bb_range_buf;
    gl_buffer_t group_buf;
//    gl_buffer_t instance_buf;

    // Allocation Data
    uint32_t atom_offset;
    uint32_t atom_capacity;

    uint32_t bond_offset;
    uint32_t bond_capacity;

    uint32_t bb_seg_offset;
    uint32_t bb_seg_capacity;

    uint32_t bb_range_offset;
    uint32_t bb_range_capacity;

    uint32_t group_offset;
    uint32_t group_capacity;

    /*
    uint32_t instance_offset;
    uint32_t instance_capacity;
    */

    uint32_t transform_offset;
    uint32_t transform_capacity;

    // Representation Buffers
    gl_buffer_t color_buf;

    uint32_t color_offset;
    uint32_t color_capacity;

    gl_program_t spacefill_prog;
    gl_program_t licorice_prog;
    gl_program_t ribbons_prog;
    gl_program_t cartoon_prog;
    gl_program_t gaussian_surface_prog;
    gl_program_t solvent_excluded_surface_prog;

    // This culls structures and generates indirect draw and compute commands
    gl_program_t cull_prog;

    // Computes depth pyramid
    gl_program_t depth_reduce_prog;

    // Write back picking data
    gl_program_t write_picking_prog;

    // Composes final image and outputs
    gl_program_t compose_prog;

    // Transient GPU resident buffers
    gl_buffer_t ubo_buf;
    gl_buffer_t vertex_buf;
    gl_buffer_t index_buf;
    gl_buffer_t transform_buf;

    gl_buffer_t draw_op_buf;
    gl_buffer_t draw_ind_buf;
    gl_buffer_t draw_sphere_idx_buf;

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
#define APPEND_LINE(a, line, alloc) APPEND_STR(a, alloc_printf(alloc, "#line %u", line), alloc)

static bool compile_shader_from_source(GLuint shader, const char* source, const char* defines, const include_file_t* include_files, uint32_t include_file_count) {
    if (!source) {
        md_print(MD_LOG_TYPE_ERROR, "Missing shader source");
        return false;
    }

    str_t str = str_from_cstr(source);

    str_t version_str;
    if (!extract_line(&version_str, &str) || !str_equal_cstr_n(version_str, "#version", 8)) {
        md_print(MD_LOG_TYPE_ERROR, "Missing version as first line in shader!");
        return false;
    }

    md_allocator_i* alloc = default_temp_allocator;

    char* buf = 0;
    md_array_ensure(buf, KILOBYTES(4), alloc);

    uint32_t line_count = 0;
    APPEND_STR(buf, version_str, alloc);
    line_count += 1;

    if (defines) {
        str_t def = str_from_cstr(defines);
        line_count += (uint32_t)str_count_char_occur(def, '\n');        
        APPEND_LINE(buf, line_count, alloc);
    }

    str_t line;
    while (extract_line(&line, &str)) {
        if (str_equal_cstr_n(line, "#include", 8)) {
            if (include_file_count == 0 || !include_files) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to parse include in shader file: no include files have been supplied");
                return false;
            }
            str_t file = str_trim_whitespace(substr(line, 8, -1));
            if (!file.len) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to parse include in shader file: file is missing");
                return false;
            }
            char beg_delim = file.ptr[0];
            char end_delim = file.ptr[file.len-1];
            if (!(beg_delim == '"' || beg_delim == '<') ||
                !(end_delim == '"' || end_delim == '>') ||
                 (beg_delim == '"' && end_delim != '"') ||
                 (beg_delim == '<' && end_delim != '>')) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to parse include in shader file: missing or mismatched delimiters");
                return false;
            }
            file = substr(file, 1, file.len-2);
            str_t src = {0};
            for (uint32_t i = 0; i < include_file_count; ++i) {
                if (str_equal_cstr(file, include_files[i].filename)) {
                    src = str_from_cstr(include_files[i].source);
                    break;
                }
            }
            if (str_empty(src)) {
                md_printf(MD_LOG_TYPE_ERROR, "Failed to parse include in shader file: could not find include file '%.*s'", (int)file.len, file.ptr);
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
        md_printf(MD_LOG_TYPE_ERROR, "Shader compile error:\n%s\n", err_buf);
        return false;
    }
    
    return true;
}

static bool compile_shader_from_file(GLuint shader, const char* filename, const char* defines, const include_file_t* include_files, uint32_t include_file_count) {
    str_t src = load_textfile(str_from_cstr(filename), default_temp_allocator);
    if (!str_empty(src)) {
        const bool success = compile_shader_from_source(shader, src.ptr, defines, include_files, include_file_count);
        const md_log_type_t log_type = success ? MD_LOG_TYPE_INFO : MD_LOG_TYPE_ERROR;
        const char* res_str = success ? "Success" : "Fail";
        md_printf(log_type,  "Compiling shader %-40s %s", filename, res_str);
        return success;
    } else {
        md_printf(MD_LOG_TYPE_ERROR, "Could not open file file '%s'", filename);
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
            md_print(MD_LOG_TYPE_ERROR, "Program link error: One or more shaders are invalid\n");
            return false;
        }
    }
    
    GLint success;
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        char err_buf[1024];
        glGetProgramInfoLog(program, ARRAY_SIZE(err_buf), NULL, err_buf);
        md_printf(MD_LOG_TYPE_ERROR, "Program link error:\n%s\n", err_buf);
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
    uint32_t shaders[8];
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
        md_print(MD_LOG_TYPE_ERROR, "gfx context has not been initialized");
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

bool md_gfx_initialize(uint32_t width, uint32_t height, md_gfx_config_flags_t flags) {
    if (!ctx.version) {
        if (gl3wInit() != GL3W_OK) {
            md_print(MD_LOG_TYPE_ERROR, "Could not load OpenGL extensions");
            return false;
        }

        GLint major, minor;
        glGetIntegerv(GL_MAJOR_VERSION, &major);
        glGetIntegerv(GL_MINOR_VERSION, &minor);

        if (major < 4 && minor < 6) {
            md_printf(MD_LOG_TYPE_ERROR, "OpenGL version %i.%i is not supported, this renderer requires version 4.6", major, minor);
            return false;
        }

        glCreateFramebuffers(1, &ctx.fbo_id);
        glCreateVertexArrays(1, &ctx.vao_id);

        ctx.ubo_buf             = gl_buffer_create(KILOBYTES(1),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.draw_op_buf         = gl_buffer_create(MEGABYTES(1),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

        // Fully GPU Resident Buffers
        ctx.draw_ind_buf        = gl_buffer_create(KILOBYTES(32), NULL, 0);
        ctx.vertex_buf          = gl_buffer_create(MEGABYTES(64), NULL, 0);
        ctx.index_buf           = gl_buffer_create(MEGABYTES(64), NULL, 0);
        ctx.draw_sphere_idx_buf = gl_buffer_create(MEGABYTES(4),  NULL, 0);

        ctx.DEBUG_BUF = gl_buffer_create(KILOBYTES(1), NULL, GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT);
        ctx.debug_data = glMapNamedBuffer(ctx.DEBUG_BUF.id, GL_READ_ONLY);

        // Storage buffers for structure data
        const uint32_t atom_cap = 4 * 1000 * 1000;
        ctx.position_buf = gl_buffer_create(atom_cap * sizeof(vec3_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.velocity_buf = gl_buffer_create(atom_cap * sizeof(vec3_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.radius_buf   = gl_buffer_create(atom_cap * sizeof(float),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.atom_capacity = atom_cap;
        
        const uint32_t bond_cap = 2 * 1000 * 1000;
        ctx.bond_buf     = gl_buffer_create(bond_cap * sizeof(bond_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.bond_capacity = bond_cap;

        const uint32_t bb_cap = 1 * 1000 * 1000;
        ctx.segment_buf  = gl_buffer_create(bb_cap * sizeof(md_gfx_backbone_segment_t),    NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.sec_str_buf  = gl_buffer_create(bb_cap * sizeof(md_gfx_secondary_structure_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.bb_seg_capacity = bb_cap;

        const uint32_t bb_range_cap = 10 * 1000;
        ctx.bb_range_buf = gl_buffer_create(bb_range_cap * sizeof(range_t),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.bb_range_capacity = bb_range_cap;

        const uint32_t group_cap = 1 * 1000 * 1000;
        ctx.group_buf    = gl_buffer_create(group_cap * sizeof(range_t),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.group_capacity = group_cap;

        /*
        const uint32_t inst_cap = 10 * 1000;
        ctx.instance_buf = gl_buffer_create(inst_cap * sizeof(instance_t),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.instance_capacity = inst_cap;
        */

        const uint32_t transform_cap = 100000;
        ctx.transform_buf = gl_buffer_create(transform_cap * sizeof(mat4_t),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.transform_capacity = transform_cap;

        // Storage buffers for representation data
        ctx.color_buf = gl_buffer_create(MEGABYTES(64), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

        str_t common_src = load_textfile(str_from_cstr(MD_SHADER_DIR "/gfx/common.h"), default_temp_allocator);
        if (str_empty(common_src)) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to read common header for shaders");
            return false;
        }

        include_file_t include_files[] = {
            {
                .filename = "common.h",
                .source = common_src.ptr,
            },
        };

        // Compile shader programs
        const char* spacefill_files[] = {
            MD_SHADER_DIR "/gfx/spacefill.vert",
            MD_SHADER_DIR "/gfx/spacefill.geom",
            MD_SHADER_DIR "/gfx/spacefill.frag",
        };
        ctx.spacefill_prog = gl_program_create_from_files(spacefill_files, ARRAY_SIZE(spacefill_files), 0, include_files, ARRAY_SIZE(include_files));

        const char* cull_files[] = {
            MD_SHADER_DIR "/gfx/cull.comp"
        };
        ctx.cull_prog = gl_program_create_from_files(cull_files, ARRAY_SIZE(cull_files), 0, include_files, ARRAY_SIZE(include_files));

        const char* depth_reduce_files[] = {
            MD_SHADER_DIR "/gfx/depth_reduce.comp",
        };
        ctx.depth_reduce_prog = gl_program_create_from_files(depth_reduce_files, ARRAY_SIZE(depth_reduce_files), 0, include_files, ARRAY_SIZE(include_files));

        const char* picking_files[] = {
            MD_SHADER_DIR "/gfx/write_picking_data.comp",
        };
        ctx.write_picking_prog = gl_program_create_from_files(picking_files, ARRAY_SIZE(picking_files), 0, include_files, ARRAY_SIZE(include_files));

        const char* compose_files[] = {
            MD_SHADER_DIR "/gfx/fs_quad.vert",
            MD_SHADER_DIR "/gfx/compose.frag",
        };
        ctx.compose_prog = gl_program_create_from_files(compose_files, ARRAY_SIZE(compose_files), 0, include_files, ARRAY_SIZE(include_files));

        if (!ctx.spacefill_prog.id || !ctx.cull_prog.id || !ctx.depth_reduce_prog.id || !ctx.write_picking_prog.id || !ctx.compose_prog.id) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to compile one or more shader programs.");
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

            gl_texture_destroy(&ctx.visibility_tex);
            ctx.visibility_tex = gl_texture_create(width, height, 1, GL_RG32UI);

            glNamedFramebufferTexture(ctx.fbo_id, GL_DEPTH_ATTACHMENT,  ctx.depth_tex.id,  0);
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
    gl_texture_destroy(&ctx.visibility_tex);

    gl_buffer_destroy(&ctx.ubo_buf);
    gl_buffer_destroy(&ctx.vertex_buf);
    gl_buffer_destroy(&ctx.index_buf);
    gl_buffer_destroy(&ctx.draw_op_buf);

    gl_buffer_destroy(&ctx.position_buf);
    gl_buffer_destroy(&ctx.velocity_buf);
    gl_buffer_destroy(&ctx.radius_buf);
    gl_buffer_destroy(&ctx.bond_buf);
    gl_buffer_destroy(&ctx.segment_buf);
    gl_buffer_destroy(&ctx.sec_str_buf);
    gl_buffer_destroy(&ctx.bb_range_buf);
    gl_buffer_destroy(&ctx.group_buf);
    //gl_buffer_destroy(&ctx.instance_buf);

    ctx = (gl_context_t){0};
}

// Allocate ranges for structure
static bool alloc_structure_data(structure_t* s, uint32_t atom_count, uint32_t bond_count, uint32_t bb_seg_count, uint32_t bb_range_count, uint32_t group_count, uint32_t instance_count) {
    // Find range (for now we consider this as a ring-buffer and will replace old data when full)
    if (ctx.atom_offset + atom_count > ctx.atom_capacity) {
        ctx.atom_offset = 0;
    }

    if (ctx.bond_offset + bond_count > ctx.bond_capacity) {
        ctx.bond_offset = 0;
    }

    if (ctx.bb_seg_offset + bb_seg_count > ctx.bb_seg_capacity) {
        ctx.bb_seg_offset = 0;
    }

    if (ctx.group_offset + group_count > ctx.group_capacity) {
        ctx.group_offset = 0;
    }

    s->atom_offset = ctx.atom_offset;
    ctx.atom_offset += atom_count;

    s->bond_offset = ctx.bond_offset;
    ctx.bond_offset += bond_count;

    s->bb_seg_offset = ctx.bb_seg_offset;
    ctx.bb_seg_offset += bb_seg_count;

    s->bb_range_offset = ctx.bb_range_offset;
    ctx.bb_range_offset += bb_range_count;

    s->group_offset = ctx.group_offset;
    ctx.group_offset += group_count;

    ASSERT(instance_count < ARRAY_SIZE(s->instance));

    s->atom_count = atom_count;
    s->bond_count = bond_count;
    s->bb_seg_count = bb_seg_count;
    s->bb_range_count = bb_range_count;
    s->group_count = group_count;
    s->instance_count = instance_count;

    return true;
}

static bool free_structure_data(structure_t* s) {
    (void)s;
    return true;
}

md_gfx_handle_t md_gfx_structure_create(uint32_t atom_count, uint32_t bond_count, uint32_t backbone_segment_count, uint32_t backbone_range_count, uint32_t group_count, uint32_t instance_count) {
    if (!validate_context()) {
        return (md_gfx_handle_t){0};
    }

    if (ctx.structure_count >= MAX_STRUCTURE_COUNT) {
        md_print(MD_LOG_TYPE_ERROR, "Unable to create structure, pool is full!");
        return (md_gfx_handle_t){0};
    }
    uint16_t s_idx = (uint16_t)ctx.structure_count++;

    structure_t* s = &ctx.structures[s_idx];
    if (!alloc_structure_data(s, atom_count, bond_count, backbone_segment_count, backbone_range_count, group_count, instance_count)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to allocate storage required for structure");
        return (md_gfx_handle_t){0};
    }

    s->h_idx = alloc_handle();

    // Map handle index to resource index
    ctx.handles[s->h_idx].idx = s_idx;

    // Create external handle for user
    return (md_gfx_handle_t) {ctx.handles[s->h_idx].gen, s->h_idx};
}

bool md_gfx_structure_destroy(md_gfx_handle_t id) {
    if (!validate_handle(id)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to destroy structure: Handle is invalid");
        return false;
    }
    uint16_t h_idx = id.idx;
    uint16_t s_idx = ctx.handles[h_idx].idx;
    free_handle(h_idx);

    ASSERT(ctx.structure_count > 0);
    ctx.structure_count -= 1;

    free_structure_data(&ctx.structures[s_idx]);
    if (ctx.structure_count > 1) {
        // Swap back and pop resource array to keep it packed
        ctx.structures[s_idx] = ctx.structures[ctx.structure_count];
        
        // Reassign handle indices to new location
        ctx.handles[ctx.structures[s_idx].h_idx].idx = s_idx;
    }

    return true;
}

bool md_gfx_structure_set_atom_position(md_gfx_handle_t id, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom position data: Handle is invalid");
        return false;
    }

    if (!x || !y || ! z) {
        md_print(MD_LOG_TYPE_ERROR, "One or more arguments are missing, must pass x, y and z for position.");
        return false;
    }

    if (offset + count > s->atom_count) {
        md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds.");
        return false;
    }

    float* ptr = glMapNamedBufferRange(ctx.position_buf.id, s->atom_offset * sizeof(vec3_t), s->atom_count * sizeof(vec3_t), GL_MAP_WRITE_BIT);
    if (!ptr) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom position data: Could not map buffer");
        return false;
    }
   
    byte_stride = MAX(sizeof(float), byte_stride);
    for (uint32_t i = offset; i < count; ++i) {
        ptr[i * 3 + 0] = *(const float*)((const uint8_t*)x + byte_stride * i);
        ptr[i * 3 + 1] = *(const float*)((const uint8_t*)y + byte_stride * i);
        ptr[i * 3 + 2] = *(const float*)((const uint8_t*)z + byte_stride * i);        
    }
    glUnmapNamedBuffer(ctx.position_buf.id);

    return true;
}

bool md_gfx_structure_set_atom_velocity(md_gfx_handle_t id, uint32_t offset, uint32_t count, const float* vx, const float* vy, const float* vz, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom velocity data: Handle is invalid");
        return false;
    }

    if (!vx || !vy || ! vz) {
        md_print(MD_LOG_TYPE_ERROR, "One or more arguments are missing, must pass vx, vy and vz for velocity.");
        return false;
    }

    if (offset + count > s->atom_count) {
        md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds.");
        return false;
    }

    float* ptr = glMapNamedBufferRange(ctx.velocity_buf.id, s->atom_offset * sizeof(vec3_t), s->atom_count * sizeof(vec3_t), GL_MAP_WRITE_BIT);
    if (!ptr) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom velocity data: Could not map buffer");
        return false;
    }

    byte_stride = MAX(sizeof(float), byte_stride);
    for (uint32_t i = offset; i < count; ++i) {
        ptr[i * 3 + 0] = *(const float*)((const uint8_t*)vx + byte_stride * i);
        ptr[i * 3 + 1] = *(const float*)((const uint8_t*)vy + byte_stride * i);
        ptr[i * 3 + 2] = *(const float*)((const uint8_t*)vz + byte_stride * i);
    }
    glUnmapNamedBuffer(ctx.velocity_buf.id);

    return true;
}

bool md_gfx_structure_set_atom_radius(md_gfx_handle_t id, uint32_t offset, uint32_t count, const float* radius, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom radius data: Handle is invalid");
        return false;
    }

    if (!radius) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom radius data: Radius argument is missing.");
        return false;
    }

    if (offset + count > s->atom_count) {
        md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds.");
        return false;
    }

    const uint32_t element_size = sizeof(float);
    if (byte_stride > element_size) {
        float* ptr = glMapNamedBufferRange(ctx.radius_buf.id, s->atom_offset * element_size, s->atom_count * element_size, GL_MAP_WRITE_BIT);
        if (!ptr) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to set atom radius data: Could not map buffer");
            return false;
        }
        byte_stride = MAX(element_size, byte_stride);
        for (uint32_t i = offset; i < count; ++i) {
            ptr[i] = *(const float*)((const uint8_t*)radius + byte_stride * i);
        }
        glUnmapNamedBuffer(ctx.radius_buf.id);
    } else {
        gl_buffer_set_sub_data(ctx.radius_buf, (s->atom_offset + offset) * element_size, count * element_size, radius);
    }
    return true;
}

bool md_gfx_structure_zero_atom_velocity(md_gfx_handle_t id) {
    structure_t* s = get_structure(id);
    if (!s) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set clear atom velocity data: Handle is invalid");
        return false;
    }

    const uint32_t element_size = sizeof(vec3_t);
    glClearNamedBufferSubData(ctx.velocity_buf.id, GL_RGB32F, s->atom_offset * element_size, s->atom_count * element_size, GL_RGB, GL_FLOAT, 0);

    return true;
}

bool md_gfx_structure_set_group_ranges(md_gfx_handle_t id, uint32_t offset, uint32_t count, const md_gfx_range_t* ranges, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom group data: Handle is invalid");
        return false;
    }

    if (!ranges) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom group data: ranges argument is missing.");
        return false;
    }

    if (offset + count > s->atom_count) {
        md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds.");
        return false;
    }

    const uint32_t element_size = sizeof(range_t);
    range_t* ptr = glMapNamedBufferRange(ctx.group_buf.id, s->group_offset * element_size, s->group_count * element_size, GL_MAP_WRITE_BIT);
    if (!ptr) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom radius data: Could not map buffer");
        return false;
    }

    byte_stride = MAX(element_size, byte_stride);
    for (uint32_t i = offset; i < count; ++i) {
        const md_gfx_range_t* in_range = (const md_gfx_range_t*)((const uint8_t*)ranges + byte_stride * i);
        ptr[i] = (range_t){in_range->beg_idx, in_range->end_idx - in_range->beg_idx};
    }
    glUnmapNamedBuffer(ctx.group_buf.id);
    return true;
}

bool md_gfx_structure_set_instance_group_ranges(md_gfx_handle_t id, uint32_t offset, uint32_t count, const md_gfx_range_t* group_ranges, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set instance group range data: Handle is invalid");
        return false;
    }

    if (!group_ranges) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set instance group range data: Argument is missing.");
        return false;
    }

    if (offset + count > s->instance_count) {
        md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds.");
        return false;
    }

    byte_stride = MAX(byte_stride, sizeof(md_gfx_range_t));
    for (uint32_t i = offset; i < offset + count; ++i) {
        md_gfx_range_t range = *(const md_gfx_range_t*)((const uint8_t*)group_ranges + byte_stride * i);
        s->instance[i].group_range.offset = range.beg_idx;
        s->instance[i].group_range.count  = range.end_idx - range.beg_idx;
    }

    return true;
}

bool md_gfx_structure_set_instance_transforms(md_gfx_handle_t id, uint32_t offset, uint32_t count, const struct mat4_t* transforms, uint32_t byte_stride) {
    structure_t* s = get_structure(id);
    if (!s) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set instance group range data: Handle is invalid");
        return false;
    }

    if (!transforms) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set instance group range data: Argument is missing.");
        return false;
    }

    if (offset + count > s->instance_count) {
        md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds.");
        return false;
    }

    byte_stride = MAX(byte_stride, sizeof(mat4_t));
    for (uint32_t i = offset; i < offset + count; ++i) {
        mat4_t mat = *(const mat4_t*)((const uint8_t*)transforms + byte_stride * i);
        s->instance[i].transform = mat;
    }

    return true;
}

static bool alloc_representation_data(representation_t* rep, uint32_t color_count) {
    // Find range (for now we consider this as a ring-buffer and will replace old data when full)
    if (ctx.color_offset + color_count > ctx.color_capacity) {
        ctx.color_offset = 0;
    }
    rep->color_offset = ctx.color_offset;
    ctx.color_offset += color_count;

    rep->color_count = color_count;

    return true;
}

static bool free_representation_data(representation_t* r) {
    (void)r;
    return true;
}

md_gfx_handle_t md_gfx_rep_create(uint32_t color_count) {
    if (!validate_context()) {
        return (md_gfx_handle_t){0};
    }

    if (ctx.representation_count >= MAX_REPRESENTATION_COUNT) {
        md_print(MD_LOG_TYPE_ERROR, "Unable to create representation, pool is full!");
        return (md_gfx_handle_t){0};
    }
    uint16_t r_idx = ctx.representation_count++;

    representation_t* r = &ctx.representations[r_idx];
    if (!alloc_representation_data(r, color_count)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to allocate storage required for representation");
        return (md_gfx_handle_t){0};
    }

    r->h_idx = alloc_handle();

    // Map handle index to resource index
    ctx.handles[r->h_idx].idx = r_idx;

    // Create external handle for user
    return (md_gfx_handle_t) {ctx.handles[r->h_idx].gen, r->h_idx};
}

bool md_gfx_rep_destroy(md_gfx_handle_t id) {
    if (!validate_handle(id)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to destroy representation: Handle is invalid");
        return false;
    }
    uint16_t h_idx = id.idx;
    uint16_t r_idx = ctx.handles[h_idx].idx;
    free_handle(h_idx);

    ASSERT(ctx.representation_count > 0);
    ctx.representation_count -= 1;

    free_representation_data(&ctx.representations[r_idx]);

    if (ctx.representation_count > 1) {
        // Swap back and pop resource array to keep it packed
        ctx.representations[r_idx] = ctx.representations[ctx.representation_count];

        // Reassign handle indices to new location
        ctx.handles[ctx.representations[r_idx].h_idx].idx = r_idx;
    }

    return true;
}

bool md_gfx_rep_set_type_and_attr(md_gfx_handle_t id, md_gfx_rep_type_t type, const md_gfx_rep_attr_t* attr) {
    representation_t* r = get_representation(id);
    if (!r) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set representation type and atrributes: Handle is invalid");
        return false;
    }

    if (!attr) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set representation type and atrributes: Attribute is missing");
        return false;
    }

    r->type = type;
    r->attr = *attr;
    return true;
}

bool md_gfx_rep_set_color(md_gfx_handle_t id, uint32_t offset, uint32_t count, const md_gfx_color_t* color, uint32_t byte_stride) {
    representation_t* r = get_representation(id);
    if (!r) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set color data: Handle is invalid");
        return false;
    }

    if (!color) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set color data: Color argument is missing.");
        return false;
    }

    if (offset + count > r->color_count) {
        md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds.");
        return false;
    }

    const uint32_t element_size = sizeof(md_gfx_color_t);
    if (byte_stride > element_size) {
        md_gfx_color_t* ptr = glMapNamedBufferRange(ctx.color_buf.id, r->color_offset * element_size, r->color_count * element_size, GL_WRITE_ONLY);
        if (!ptr) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to set color data: Could not map buffer");
            return false;
        }
        byte_stride = MAX(element_size, byte_stride);
        for (uint32_t i = offset; i < count; ++i) {
            ptr[i] = *(const md_gfx_color_t*)((const uint8_t*)color + byte_stride * i);
        }
        glUnmapNamedBuffer(ctx.color_buf.id);
    } else {
        gl_buffer_set_sub_data(ctx.color_buf, (r->color_offset + offset) * element_size, count * element_size, color);
    }
    return true;
}

bool md_gfx_draw(uint32_t in_draw_op_count, const md_gfx_draw_op_t* in_draw_ops, const mat4_t* proj_mat, const mat4_t* view_mat, const mat4_t* inv_proj_mat, const mat4_t* inv_view_mat) {
    if (!validate_context()) return false;

    if (in_draw_op_count > 0 && !in_draw_ops) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to Draw: draw ops was NULL");
        return false;
    }

    if (!proj_mat || !view_mat || !inv_proj_mat || !inv_view_mat) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to Draw: one or more missing matrix arguments");
        return false;
    }

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

    glClearColor(0,0,0,0);
    glClearDepth(1.0);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    mat4_t prev_transform = mat4_mul(ctx.prev_proj_mat, ctx.prev_view_mat);
    mat4_t inv_transform = mat4_mul(*inv_view_mat, *inv_proj_mat);
    mat4_t curr_clip_to_prev_clip_mat = mat4_mul(prev_transform, inv_transform);

    ctx.prev_proj_mat = *proj_mat;
    ctx.prev_view_mat = *view_mat;

    UniformData ubo_data = {
        .world_to_view = *view_mat,
        .world_to_view_normal = mat4_transpose(*inv_view_mat),
        .view_to_clip = *proj_mat,
        .clip_to_view = *inv_proj_mat,
        .curr_clip_to_prev_clip = curr_clip_to_prev_clip_mat,
        .prev_world_to_clip = prev_transform,
        .depth_pyramid_width = ctx.depth_pyramid_tex.width,
        .depth_pyramid_height = ctx.depth_pyramid_tex.height,
        .frustum_culling = 1,
        .depth_culling = 1,
        .znear = extract_znear(*inv_proj_mat),
        .zfar = extract_zfar(*inv_proj_mat),
    };

    gl_buffer_set_sub_data(ctx.ubo_buf, 0, sizeof(ubo_data), &ubo_data);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ctx.ubo_buf.id);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, POSITION_BINDING,            ctx.position_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, RADIUS_BINDING,              ctx.radius_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, GROUP_BINDING,               ctx.group_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, COLOR_BINDING,               ctx.color_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, TRANSFORM_BINDING,           ctx.transform_buf.id);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_OP_BINDING,             ctx.draw_op_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_INDIRECT_BINDING,       ctx.draw_ind_buf.id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_SPHERE_INDEX_BINDING,   ctx.draw_sphere_idx_buf.id);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DEBUG_BINDING,               ctx.DEBUG_BUF.id);

    // Clear buffer
    glClearNamedBufferData(ctx.draw_ind_buf.id, GL_R32UI, GL_RED, GL_UNSIGNED_INT, 0);

    uint32_t draw_op_count = 0;
    DrawOp   draw_ops[512];

    uint32_t transform_count = 0;
    mat4_t transforms[128];
    transforms[transform_count++] = mat4_ident();

    for (uint32_t in_idx = 0; in_idx < in_draw_op_count; ++in_idx) {
        const md_gfx_draw_op_t* in_op = in_draw_ops + in_idx;

        // Validate draw op
        structure_t* s = get_structure(in_op->structure);
        if (!s) {
            md_print(MD_LOG_TYPE_ERROR, "Invalid structure supplied in draw operation");
            continue;   
        }
        representation_t* r = get_representation(in_op->representation);
        if (!r) {
            md_print(MD_LOG_TYPE_ERROR, "Invalid representation supplied in draw operation");
            continue;
        }

        if (s->atom_count > r->color_count) {
            md_print(MD_LOG_TYPE_ERROR, "Representation does not hold a sufficient quantity of colors to match the structure in draw operation");
            continue;
        }

        if (s->instance_count > 0) {
            DrawOp out_op;
            out_op.group_offset  = s->group_offset;
            out_op.group_count   = s->group_count;
            out_op.atom_offset   = s->atom_offset;
            out_op.color_offset  = r->color_offset;
            out_op.rep_type      = MD_GFX_REP_TYPE_SPACEFILL;
            out_op.rep_args[0]   = r->attr.spacefill.radius_scale;

            for (uint32_t inst_idx = 0; inst_idx < s->instance_count; ++inst_idx) {
                uint32_t transform_idx = 0;
                const instance_t* inst = s->instance + inst_idx;
                mat4_t model_mat = in_op->model_mat ? mat4_mul(*in_op->model_mat, inst->transform) : inst->transform;
                if (!mat4_equal(model_mat, mat4_ident())) {
                    transform_idx = transform_count;
                    transforms[transform_count++] = model_mat;
                }
                out_op.transform_idx = ctx.transform_offset + transform_idx;

                // We only process 1024 groups at a time in the compute shader
                // If the structure contains more than this, we split it into several draw ops
                uint32_t out_count = DIV_UP(inst->group_range.count, CULL_GROUP_SIZE);
                for (uint32_t i = 0; i < out_count; ++i) {
                    out_op.group_offset = inst->group_range.offset + i * CULL_GROUP_SIZE;
                    out_op.group_count  = (i < (out_count - 1)) ? CULL_GROUP_SIZE : inst->group_range.count % CULL_GROUP_SIZE;
                    ASSERT(draw_op_count < ARRAY_SIZE(draw_ops));
                    draw_ops[draw_op_count++] = out_op;
                }
            }
        } else {
            uint32_t transform_idx = 0; // This is the default and represents identity matrix
            if (in_op->model_mat && !mat4_equal(*in_op->model_mat, mat4_ident())) {
                ASSERT(transform_count < ARRAY_SIZE(transforms));
                transform_idx = transform_count;
                transforms[transform_count++] = *in_op->model_mat;
            }

            // Create internal draw op and fill buffer
            DrawOp out_op;
            out_op.group_offset  = s->group_offset;
            out_op.group_count   = s->group_count;
            out_op.atom_offset   = s->atom_offset;
            out_op.color_offset  = r->color_offset;
            out_op.rep_type      = MD_GFX_REP_TYPE_SPACEFILL;
            out_op.rep_args[0]   = r->attr.spacefill.radius_scale;
            out_op.transform_idx = ctx.transform_offset + transform_idx;

            // We only process 1024 groups at a time in the compute shader
            // If the structure contains more than this, we split it into several draw ops
            uint32_t out_count = DIV_UP(s->group_count, CULL_GROUP_SIZE);
            for (uint32_t i = 0; i < out_count; ++i) {
                out_op.group_offset = s->group_offset + i * CULL_GROUP_SIZE;
                out_op.group_count  = (i < (out_count - 1)) ? CULL_GROUP_SIZE : s->group_count % CULL_GROUP_SIZE;
                ASSERT(draw_op_count < ARRAY_SIZE(draw_ops));
                draw_ops[draw_op_count++] = out_op;
            }
        }
    }

GL_PUSH_GPU_SECTION("UPLOAD DRAW OPS + TRANSFORMS");
    gl_buffer_set_sub_data(ctx.draw_op_buf, 0, sizeof(DrawOp) * draw_op_count, draw_ops);
    gl_buffer_set_sub_data(ctx.transform_buf, ctx.transform_offset * sizeof(mat4_t), sizeof(mat4_t) * transform_count, transforms);
GL_POP_GPU_SECTION();

GL_PUSH_GPU_SECTION("CULL + GENERATE DRAW COMMANDS");
    // CULL + GENERATE DRAW DATA
    glUseProgram(ctx.cull_prog.id);
    glBindTextureUnit(0, ctx.depth_pyramid_tex.id);
    glDispatchCompute(draw_op_count, 1, 1);
GL_POP_GPU_SECTION();

    // BARRIER
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_COMMAND_BARRIER_BIT);

GL_PUSH_GPU_SECTION("SPACEFILL");
    // DRAW SPACEFILL
    glBindBuffer(GL_PARAMETER_BUFFER, ctx.draw_ind_buf.id);
    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, ctx.draw_ind_buf.id);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ctx.draw_sphere_idx_buf.id);
    glUseProgram(ctx.spacefill_prog.id);
    glMultiDrawElementsIndirectCount(GL_POINTS, GL_UNSIGNED_INT, (const void*)offsetof(DrawIndirect, draw_sphere_cmd), (GLintptr)offsetof(DrawIndirect, draw_sphere_cmd_count), draw_op_count, sizeof(DrawSpheresIndirectCommand));
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0);
    glBindBuffer(GL_PARAMETER_BUFFER, 0);
GL_POP_GPU_SECTION();

    // DRAW REST
    // TODO

    glMemoryBarrier(GL_FRAMEBUFFER_BARRIER_BIT);

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

    // COMPUTE DEPTH PYRAMID FOR NEXT ITERATION
GL_PUSH_GPU_SECTION("DEPTH REDUCE");
    glUseProgram(ctx.depth_reduce_prog.id);
    for (uint32_t i = 0; i < ctx.depth_pyramid_tex.levels; ++i) {
        // Output
        uint width  = MAX(1, ctx.depth_pyramid_tex.width >> i);
        uint height = MAX(1, ctx.depth_pyramid_tex.height >> i);
        glUniform2ui(0, width, height);
        glUniform1f (1, (float)MAX(0, (int)i-1));
        glBindTextureUnit(1, i == 0 ? ctx.depth_tex.id : ctx.depth_pyramid_tex.id);
        glBindImageTexture(0, ctx.depth_pyramid_tex.id, i, 0, 0, GL_WRITE_ONLY, GL_R32F);
        glDispatchCompute(DIV_UP(width, DEPTH_REDUCE_GROUP_SIZE), DIV_UP(height, DEPTH_REDUCE_GROUP_SIZE), 1);
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    }
GL_POP_GPU_SECTION();

GL_PUSH_GPU_SECTION("RESET STATE");
    glUseProgram(0);
    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0);
    glBindBuffer(GL_PARAMETER_BUFFER, 0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, POSITION_BINDING,            0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, RADIUS_BINDING,              0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, GROUP_BINDING,               0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, COLOR_BINDING,               0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, TRANSFORM_BINDING,           0);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_OP_BINDING,             0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_INDIRECT_BINDING,       0);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, DRAW_SPHERE_INDEX_BINDING,   0);

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
    return ctx.debug_data->pick_index;
}

float md_gfx_get_picking_depth() {
    if (!validate_context()) return 1.0f;
    return ctx.debug_data->pick_depth;
}