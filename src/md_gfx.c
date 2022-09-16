#include "md_gfx.h"

#include "md_util.h"
#include <gl3w/gl3w.h>

#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_file.h>
#include <core/md_vec_math.h>
#include <core/md_compiler.h>
#include <core/md_allocator.h>
#include <core/md_array.inl>

#include <stdbool.h>
#include <string.h>     // memset, memcpy

// For subdividing backbone segments
#define MIN_SPLINE_SUBDIVISION_COUNT 1
#define MAX_SPLINE_SUBDIVISION_COUNT 32

#define MAX_STRUCTURE_COUNT 64
#define MAX_REPRESENTATION_COUNT 64
#define MAX_HANDLE_COUNT (MAX_STRUCTURE_COUNT + MAX_REPRESENTATION_COUNT)

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

static inline bool is_ortho_proj_matrix(const mat4_t M) {
    return M.elem[2][3] == 0.0f;
}

static inline void extract_jitter_uv(float jitter[2], const mat4_t M) {
    if (is_ortho_proj_matrix(M)) {
        jitter[0] = -M.elem[3][0] * 0.5f;
        jitter[1] = -M.elem[3][1] * 0.5f;
    }
    else {
        jitter[0] = M.elem[2][0] * 0.5f;
        jitter[1] = M.elem[2][1] * 0.5f;
    }
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

// Shared ubo data for all shaders
typedef struct gl_ubo_base_t {
    mat4_t world_to_view;
    mat4_t world_to_view_normal;
    mat4_t world_to_clip;
    mat4_t view_to_clip;
} gl_ubo_base_t;

typedef uint32_t gl_version_t;

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

    uint32_t instance_offset;
    uint32_t instance_count;

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

    uint32_t fbo;
    gl_texture_t depth_tex;      // depth buffer for depth testing in fixed pipeline
    gl_texture_t visibility_tex; // 64-bit buffer with: hi: 32-bit depth + lo: 32-bit payload
    gl_texture_t hi_z_tex;       // depth pyramid for occlusion culling (smaller power of two)

    // Transient GPU resident
    gl_buffer_t ubo_buf;
    gl_buffer_t vertex_buf;
    gl_buffer_t index_buf;
    gl_buffer_t drawcall_buf;
    gl_buffer_t dispatch_buf;
    gl_buffer_t command_buf;

    // Structure Buffers
    gl_buffer_t position_buf;
    gl_buffer_t velocity_buf;
    gl_buffer_t radius_buf;
    gl_buffer_t bond_buf;
    gl_buffer_t segment_buf;
    gl_buffer_t sec_str_buf;
    gl_buffer_t bb_range_buf;
    gl_buffer_t group_buf;
    gl_buffer_t instance_buf;

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

    uint32_t instance_offset;
    uint32_t instance_capacity;

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

    structure_t structures[MAX_STRUCTURE_COUNT];
    representation_t representations[MAX_REPRESENTATION_COUNT];

    md_gfx_handle_t handles[MAX_HANDLE_COUNT];

    uint16_t representation_count;
    uint16_t structure_count;
    uint16_t handle_next_idx;
} gl_context_t;



// Internal packed structures which maps to the layout used in GPU-memory
// Pack the first index as an absolute offset and second as a relative to the first
typedef struct bond_t {
    uint32_t first : 20;
    uint32_t other : 12;
} bond_t;

typedef struct range_t {
    uint32_t offset : 20;
    uint32_t length : 12;
} range_t;

typedef struct segment_t {
    uint32_t offset : 20;
    uint32_t ca_idx : 4;
    uint32_t c_idx  : 4;
    uint32_t o_idx  : 4;
} segment_t;

typedef struct instance_t {
    range_t  atom_range;
    uint32_t transform_idx;
} instance_t;

// These are absolute for now, but may be compacted later on
typedef struct position_t {
    float x, y, z;
} position_t;

typedef struct velocity_t {
    float x, y, z;
} velocity_t;

STATIC_ASSERT(sizeof(bond_t) == sizeof(uint32_t), "Invalid size of bond");
STATIC_ASSERT(sizeof(range_t) == sizeof(uint32_t), "Invalid size of range");
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

static bool compile_shader_from_source(GLuint shader, char* shader_src, const char* defines) {
    const char* sources[8] = {0};
    uint32_t num_sources = 0;

    if (defines) {
        const char* src = shader_src;
        char version_str[32] = {0};
        if (!extract_version_string(version_str, sizeof(version_str), &src)) {
            md_print(MD_LOG_TYPE_ERROR, "Missing version string as first line in shader!");
            return false;
        }

        sources[num_sources++] = version_str;
        sources[num_sources++] = "\n";

        if (defines) {
            sources[num_sources++] = defines;
            sources[num_sources++] = "\n";
        }

        sources[num_sources++] = src;

        glShaderSource(shader, num_sources, sources, 0);
    } else {
        const char* c_src = shader_src;
        glShaderSource(shader, 1, &c_src, 0);
    }
    
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

static bool compile_shader_from_file(GLuint shader, const char* filename, const char* defines) {
    md_file_o* file = md_file_open((str_t){.ptr = filename, .len = strlen(filename)}, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        char buffer[KILOBYTES(32)] = {0};
        uint64_t file_size = md_file_size(file);
        const uint64_t size = MIN(file_size, ARRAY_SIZE(buffer) - 1);
        md_file_read(file, buffer, size);
        const bool success = compile_shader_from_source(shader, buffer, defines);
        const md_log_type_t log_type = success ? MD_LOG_TYPE_INFO : MD_LOG_TYPE_ERROR;
        const char* res_str = success ? "Success" : "Fail";
        md_printf(log_type,  "Compiling shader %-40s %s", filename, res_str);
        md_file_close(file);
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

static gl_program_t gl_program_create_from_files(const char* filenames[], uint32_t count, const char* defines) {
    uint32_t shaders[8];
    ASSERT(count < ARRAY_SIZE(shaders));
    gl_program_t prog = {0};

    uint32_t flags = 0;
    for (uint32_t i = 0; i < count; ++i) {
        uint32_t type = gl_extract_shader_type_from_file_extension(filenames[i]);
        flags |= gl_get_shader_type_bit(type);
        shaders[i] = glCreateShader(type);
        compile_shader_from_file(shaders[i], filenames[i], defines);
    }
    
    uint32_t id = glCreateProgram();
    if (link_program(id, shaders, count)) {
        prog.id = id;
        prog.flags = flags;  
    } else {
        glDeleteProgram(id);
    }

    for (uint32_t i = 0; i < count; ++i) {
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
        ctx.handles[i].gen = (uint16_t)((i + 1) * 1236285 ^ 23645 + 123745);
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

        glGenFramebuffers(1, &ctx.fbo);

        ctx.hi_z_tex        = gl_texture_create(1024, 1024, 10, GL_DEPTH_COMPONENT32F);

        ctx.ubo_buf      = gl_buffer_create(KILOBYTES(1), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);

        ctx.vertex_buf   = gl_buffer_create(MEGABYTES(64), NULL, 0);
        ctx.index_buf    = gl_buffer_create(MEGABYTES(64), NULL, 0);
        ctx.drawcall_buf = gl_buffer_create(MEGABYTES(8),  NULL, 0);
        ctx.dispatch_buf = gl_buffer_create(MEGABYTES(8),  NULL, 0);
        ctx.command_buf  = gl_buffer_create(sizeof(gl_command_buf_t), NULL, 0);

        const uint32_t atom_cap = 4 * 1000 * 1000;
        ctx.position_buf = gl_buffer_create(atom_cap * sizeof(vec3_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.velocity_buf = gl_buffer_create(atom_cap * sizeof(vec3_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.radius_buf   = gl_buffer_create(atom_cap * sizeof(float),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.atom_capacity = atom_cap;
        
        const uint32_t bond_cap = 2 * 1000 * 1000;
        ctx.bond_buf     = gl_buffer_create(bond_cap * sizeof(md_gfx_bond_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.bond_capacity = bond_cap;

        const uint32_t bb_cap = 1 * 1000 * 1000;
        ctx.segment_buf  = gl_buffer_create(bb_cap * sizeof(md_gfx_backbone_segment_t),    NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.sec_str_buf  = gl_buffer_create(bb_cap * sizeof(md_gfx_secondary_structure_t), NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.bb_seg_capacity = bb_cap;

        const uint32_t bb_range_cap = 10 * 1000;
        ctx.bb_range_buf = gl_buffer_create(bb_range_cap * sizeof(md_gfx_range_t),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.bb_range_capacity = bb_range_cap;

        const uint32_t group_cap = 1 * 1000 * 1000;
        ctx.group_buf    = gl_buffer_create(group_cap * sizeof(md_gfx_range_t),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.group_capacity = group_cap;

        const uint32_t inst_cap = 10 * 1000;
        ctx.instance_buf = gl_buffer_create(inst_cap * sizeof(md_gfx_instance_t),  NULL, GL_DYNAMIC_STORAGE_BIT | GL_MAP_WRITE_BIT);
        ctx.instance_capacity = inst_cap;

        // Compile shader programs
        const char* spacefill_files[] = {
            MD_SHADER_DIR "/gfx/spacefill.vert",
            MD_SHADER_DIR "/gfx/spacefill.frag",
        };
        ctx.spacefill_prog = gl_program_create_from_files(spacefill_files, ARRAY_SIZE(spacefill_files), "");

        if (!ctx.spacefill_prog.id) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to compile one or more shader programs.");
            md_gfx_shutdown();
            return false;
        }

        init_handles();

        // We only set version when everything is successful
        ctx.version = major * 100 + minor * 10;
    }

    if (ctx.version) {
        gl_texture_destroy(&ctx.depth_tex);
        ctx.depth_tex       = gl_texture_create(width, height, 1, GL_DEPTH_COMPONENT32F);

        ctx.visibility_tex  = gl_texture_create(width, height, 1, GL_RG32UI);
        gl_texture_destroy(&ctx.visibility_tex);
    }

    return true;
}

void md_gfx_shutdown() {
    if (ctx.fbo) glDeleteFramebuffers(1, &ctx.fbo);

    gl_texture_destroy(&ctx.hi_z_tex);
    gl_texture_destroy(&ctx.depth_tex);
    gl_texture_destroy(&ctx.visibility_tex);

    gl_buffer_destroy(&ctx.ubo_buf);
    gl_buffer_destroy(&ctx.vertex_buf);
    gl_buffer_destroy(&ctx.index_buf);
    gl_buffer_destroy(&ctx.drawcall_buf);
    gl_buffer_destroy(&ctx.dispatch_buf);
    gl_buffer_destroy(&ctx.command_buf);

    gl_buffer_destroy(&ctx.position_buf);
    gl_buffer_destroy(&ctx.velocity_buf);
    gl_buffer_destroy(&ctx.radius_buf);
    gl_buffer_destroy(&ctx.bond_buf);
    gl_buffer_destroy(&ctx.segment_buf);
    gl_buffer_destroy(&ctx.sec_str_buf);
    gl_buffer_destroy(&ctx.bb_range_buf);
    gl_buffer_destroy(&ctx.group_buf);
    gl_buffer_destroy(&ctx.instance_buf);

    ctx = (gl_context_t){0};
}

// Allocate ranges for structure
bool allocate_structure(structure_t* s, uint32_t atom_count, uint32_t bond_count, uint32_t bb_seg_count, uint32_t bb_range_count, uint32_t group_count, uint32_t instance_count) {
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

    if (ctx.instance_offset + instance_count > ctx.instance_capacity) {
        ctx.instance_offset = 0;
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

    s->instance_offset = ctx.instance_offset;
    ctx.instance_offset += instance_count;

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
    if (!allocate_structure(s, atom_count, bond_count, backbone_segment_count, backbone_range_count, group_count, instance_count)) {
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

    if (ctx.structure_count > 1) {
        // Swap back and pop resource array to keep it packed
        ctx.structures[s_idx] = ctx.structures[ctx.structure_count];
        
        // Reassign handle indices to new location
        ctx.handles[ctx.structures[s_idx].h_idx].idx = s_idx;
    }

    return true;
}

bool md_gfx_structure_set_atom_position(md_gfx_handle_t id, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t byte_stride) {
    if (!validate_handle(id)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom position data: Handle is invalid");
        return false;
    }

    if (!x || !y || ! z) {
        md_print(MD_LOG_TYPE_ERROR, "One or more arguments are missing, must pass x, y and z for position.");
        return false;
    }

    structure_t* s = &ctx.structures[id.idx];

    if (offset + count > s->atom_count) {
        md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds.");
        return false;
    }

    float* ptr = glMapNamedBufferRange(ctx.position_buf.id, s->atom_offset * sizeof(vec3_t), s->atom_count * sizeof(vec3_t), GL_WRITE_ONLY);
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
    if (!validate_handle(id)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom velocity data: Handle is invalid");
        return false;
    }

    if (!vx || !vy || ! vz) {
        md_print(MD_LOG_TYPE_ERROR, "One or more arguments are missing, must pass vx, vy and vz for velocity.");
        return false;
    }

    structure_t* s = &ctx.structures[id.idx];

    if (offset + count > s->atom_count) {
        md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds.");
        return false;
    }

    float* ptr = glMapNamedBufferRange(ctx.velocity_buf.id, s->atom_offset * sizeof(vec3_t), s->atom_count * sizeof(vec3_t), GL_WRITE_ONLY);
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
    if (!validate_handle(id)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom radius data: Handle is invalid");
        return false;
    }

    if (!radius) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set atom radius data: Radius argument is missing.");
        return false;
    }

    structure_t* s = &ctx.structures[id.idx];
    const uint32_t element_size = sizeof(float);

    if (offset + count > s->atom_count) {
        md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds.");
        return false;
    }

    if (byte_stride > element_size) {
        float* ptr = glMapNamedBufferRange(ctx.radius_buf.id, s->atom_offset * element_size, s->atom_count * element_size, GL_WRITE_ONLY);
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
        gl_buffer_set_sub_data(ctx.radius_buf, offset * element_size, count * element_size, radius);
    }
    return true;
}

bool md_gfx_structure_zero_atom_velocity(md_gfx_handle_t id) {
    if (!validate_handle(id)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to zero atom velocity data: Handle is invalid");
        return false;
    }

    structure_t* s = &ctx.structures[id.idx];
    const uint32_t element_size = sizeof(vec3_t);

    glClearNamedBufferSubData(ctx.velocity_buf.id, GL_RGB32F, s->atom_offset * element_size, s->atom_count * element_size, GL_RGB, GL_FLOAT, 0);

    return true;
}

bool allocate_representation(representation_t* rep, uint32_t color_count) {
    // Find range (for now we consider this as a ring-buffer and will replace old data when full)
    if (ctx.color_offset + color_count > ctx.color_capacity) {
        ctx.color_offset = 0;
    }
    rep->color_offset = ctx.color_offset;
    ctx.color_offset += color_count;

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
    if (!allocate_representation(r, color_count)) {
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

    if (ctx.representation_count > 1) {
        // Swap back and pop resource array to keep it packed
        ctx.representations[r_idx] = ctx.representations[ctx.representation_count];

        // Reassign handle indices to new location
        ctx.handles[ctx.representations[r_idx].h_idx].idx = r_idx;
    }

    return true;
}

bool md_gfx_rep_set_type_and_attr(md_gfx_handle_t id, md_gfx_rep_type_t type, const md_gfx_rep_attr_t* attr) {
    if (!validate_handle(id)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set representation type and atrributes: Handle is invalid");
        return false;
    }

    if (!attr) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set representation type and atrributes: Attribute is missing");
        return false;
    }

    representation_t* r = &ctx.representations[id.idx];
    r->type = type;
    r->attr = *attr;
    return true;
}

bool md_gfx_rep_set_color(md_gfx_handle_t id, uint32_t offset, uint32_t count, const md_gfx_color_t* color, uint32_t byte_stride) {
    if (!validate_handle(id)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set color data: Handle is invalid");
        return false;
    }

    if (!color) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to set color data: Color argument is missing.");
        return false;
    }

    representation_t* r = &ctx.representations[id.idx];
    const uint32_t element_size = sizeof(md_gfx_color_t);

    if (offset + count > r->color_count) {
        md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds.");
        return false;
    }

    if (byte_stride > element_size) {
        float* ptr = glMapNamedBufferRange(ctx.color_buf.id, r->color_offset * element_size, r->color_count * element_size, GL_WRITE_ONLY);
        if (!ptr) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to set color data: Could not map buffer");
            return false;
        }
        byte_stride = MAX(element_size, byte_stride);
        for (uint32_t i = offset; i < count; ++i) {
            ptr[i] = *(const float*)((const uint8_t*)color + byte_stride * i);
        }
        glUnmapNamedBuffer(ctx.color_buf.id);
    } else {
        gl_buffer_set_sub_data(ctx.color_buf, offset * element_size, count * element_size, color);
    }
    return true;
}

bool md_gfx_draw(uint32_t draw_op_count, const md_gfx_draw_op_t* draw_ops, const mat4_t* proj_mat, const mat4_t* view_mat) {
    if (!validate_context()) return false;

    if (!draw_ops) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to Draw: draw ops was NULL");
        return false;
    }
            
    gl_ubo_base_t ubo_data = {
        .world_to_view = *view_mat,
        .world_to_view_normal = mat4_transpose(mat4_inverse(*view_mat)),
        .world_to_clip = mat4_mul(*proj_mat, *view_mat),
        .view_to_clip = *proj_mat,
    };
    gl_buffer_set_sub_data(ctx.ubo_buf, 0, sizeof(ubo_data), &ubo_data);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ctx.ubo_buf.id);

    return true;
}
