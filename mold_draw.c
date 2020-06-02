#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#pragma warning( disable : 4204)            // non-constant aggregate initialization (fully supported in proper c compilers)
#define _Thread_local __declspec(thread)    // _Thread_local is standardized in C11
#endif

#include "mold_draw.h"
#include "mold_util.h"
#include "ext/gl3w.h"

#include <string.h>     // memset, memcpy
#include <stdbool.h>
#include <stdlib.h>     // qsort
#include <stdio.h>      // printf

#ifndef ASSERT
#include <assert.h>
#define ASSERT assert
#endif

#define internal static

static _Thread_local char error_buf[256];

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

#define ARRAY_SIZE(x) (sizeof(x) / sizeof(x[0]))
#define KILOBYTE(x) (x * 1024)

#define PUSH_GPU_SECTION(lbl)                                                                       \
{                                                                                               \
    if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
}
#define POP_GPU_SECTION()                       \
{                                           \
    if (glPopDebugGroup) glPopDebugGroup(); \
}

#define MOLD_UBO_SIZE (1 << 10)
#define MOLD_SHADER_BUF_SIZE KILOBYTE(14)

enum {
    GL_VERSION_UNKNOWN = 0,
    GL_VERSION_330,
    GL_VERSION_430
};

enum {
    GL_PROGRAM_DOWNSAMPLE_DEPTH,
    GL_PROGRAM_COMPUTE_AABB_330,
    GL_PROGRAM_COMPUTE_AABB_430,
    GL_PROGRAM_CULL_AABB,
    GL_PROGRAM_EXTRACT_CONTROL_POINTS,
    GL_PROGRAM_COMPUTE_SPLINE,
    GL_PROGRAM_SPACE_FILL,
    GL_PROGRAM_LICORICE,
    GL_PROGRAM_RIBBONS,
    GL_PROGRAM_CARTOON,
    GL_PROGRAM_COUNT
};

enum {
    GL_TEXTURE_BUFFER_0,
    GL_TEXTURE_BUFFER_1,
    GL_TEXTURE_BUFFER_2,
    GL_TEXTURE_BUFFER_3,
    GL_TEXTURE_MAX_DEPTH,
    GL_TEXTURE_COUNT
};

enum {
    GL_BUFFER_MOL_ATOM_POSITION,                   // vec3
    GL_BUFFER_MOL_ATOM_PREV_POSITION,              // vec3
    GL_BUFFER_MOL_ATOM_RADIUS,                     // float
    GL_BUFFER_MOL_ATOM_FLAGS,                      // u8
    GL_BUFFER_MOL_RESIDUE_ATOM_RANGE,              // u32[2]   (beg, end)
    GL_BUFFER_MOL_RESIDUE_AABB,                    // vec3[2]  (aabb_min, aabb_max)
    GL_BUFFER_MOL_RESIDUE_VISIBLE,                 // int
    GL_BUFFER_MOL_RESIDUE_SECONDARY_STRUCTURE,     // u8       (0: Unknown, 1: Coil, 2: Helix, 3: Sheet)
    GL_BUFFER_MOL_BOND_ATOM_INDICES,               // u32[2]
    GL_BUFFER_MOL_BACKBONE_DATA,                   // u32[5] residue index, segment index, CA index, C index and O Index
    GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_DATA,     // Extracted control points before spline subdivision
    GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_INDEX,    // u32, LINE_STRIP_ADJACENCY Indices for spline subdivision, seperated by primitive restart index 0xFFFFFFFF
    GL_BUFFER_MOL_BACKBONE_SPLINE_DATA,            // Subdivided control points of spline
    GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX,           // u32, LINE_STRIP indices for rendering, seperated by primitive restart index 0xFFFFFFFF
    GL_BUFFER_MOL_COUNT
};

typedef union vec4 {
    struct {
        float x, y, z, w;
    };
    float elem[4];
} vec4;

typedef union mat4 {
    float elem[4][4];
    vec4 col[4];
} mat4;

internal inline vec4 vec4_mul(const vec4 a, const vec4 b) {
    vec4 c;
    c.x = a.x * b.x;
    c.y = a.y * b.y;
    c.z = a.z * b.z;
    c.w = a.w * b.w;
    return c;
}

internal inline vec4 vec4_add(const vec4 a, const vec4 b) {
    vec4 c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    c.w = a.w + b.w;
    return c;
}

internal inline vec4 vec4_sub(const vec4 a, const vec4 b) {
    vec4 c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    c.w = a.w - b.w;
    return c;
}

internal inline mat4 mat4_mul(const mat4 A, const mat4 B) {
    mat4 C;
#define MULT(col, row) \
    A.elem[0][row] * B.elem[col][0] + A.elem[1][row] * B.elem[col][1] + A.elem[2][row] * B.elem[col][2] + A.elem[3][row] * B.elem[col][3]
        C.elem[0][0] = MULT(0, 0);
    C.elem[0][1] = MULT(0, 1);
    C.elem[0][2] = MULT(0, 2);
    C.elem[0][3] = MULT(0, 3);
    
    C.elem[1][0] = MULT(1, 0);
    C.elem[1][1] = MULT(1, 1);
    C.elem[1][2] = MULT(1, 2);
    C.elem[1][3] = MULT(1, 3);
    
    C.elem[2][0] = MULT(2, 0);
    C.elem[2][1] = MULT(2, 1);
    C.elem[2][2] = MULT(2, 2);
    C.elem[2][3] = MULT(2, 3);
    
    C.elem[3][0] = MULT(3, 0);
    C.elem[3][1] = MULT(3, 1);
    C.elem[3][2] = MULT(3, 2);
    C.elem[3][3] = MULT(3, 3);
#undef MULT
    return C;
}

internal inline mat4 mat4_mul_f(const mat4 M, float s) {
    mat4 C;
    
    C.elem[0][0] = M.elem[0][0] * s;
    C.elem[0][1] = M.elem[0][1] * s;
    C.elem[0][2] = M.elem[0][2] * s;
    C.elem[0][3] = M.elem[0][3] * s;
    
    C.elem[1][0] = M.elem[1][0] * s;
    C.elem[1][1] = M.elem[1][1] * s;
    C.elem[1][2] = M.elem[1][2] * s;
    C.elem[1][3] = M.elem[1][3] * s;
    
    C.elem[2][0] = M.elem[2][0] * s;
    C.elem[2][1] = M.elem[2][1] * s;
    C.elem[2][2] = M.elem[2][2] * s;
    C.elem[2][3] = M.elem[2][3] * s;
    
    C.elem[3][0] = M.elem[3][0] * s;
    C.elem[3][1] = M.elem[3][1] * s;
    C.elem[3][2] = M.elem[3][2] * s;
    C.elem[3][3] = M.elem[3][3] * s;
    
    return C;
}

internal inline mat4 mat4_identity() {
    mat4 M = {0};
    M.elem[0][0] = 1.0f;
    M.elem[1][1] = 1.0f;
    M.elem[2][2] = 1.0f;
    M.elem[3][3] = 1.0f;
    return M;
}

internal inline mat4 mat4_inverse(const mat4 M) {
    const float c00 = M.elem[2][2] * M.elem[3][3] - M.elem[3][2] * M.elem[2][3];
    const float c02 = M.elem[1][2] * M.elem[3][3] - M.elem[3][2] * M.elem[1][3];
    const float c03 = M.elem[1][2] * M.elem[2][3] - M.elem[2][2] * M.elem[1][3];
    
    const float c04 = M.elem[2][1] * M.elem[3][3] - M.elem[3][1] * M.elem[2][3];
    const float c06 = M.elem[1][1] * M.elem[3][3] - M.elem[3][1] * M.elem[1][3];
    const float c07 = M.elem[1][1] * M.elem[2][3] - M.elem[2][1] * M.elem[1][3];
    
    const float c08 = M.elem[2][1] * M.elem[3][2] - M.elem[3][1] * M.elem[2][2];
    const float c10 = M.elem[1][1] * M.elem[3][2] - M.elem[3][1] * M.elem[1][2];
    const float c11 = M.elem[1][1] * M.elem[2][2] - M.elem[2][1] * M.elem[1][2];
    
    const float c12 = M.elem[2][0] * M.elem[3][3] - M.elem[3][0] * M.elem[2][3];
    const float c14 = M.elem[1][0] * M.elem[3][3] - M.elem[3][0] * M.elem[1][3];
    const float c15 = M.elem[1][0] * M.elem[2][3] - M.elem[2][0] * M.elem[1][3];
    
    const float c16 = M.elem[2][0] * M.elem[3][2] - M.elem[3][0] * M.elem[2][2];
    const float c18 = M.elem[1][0] * M.elem[3][2] - M.elem[3][0] * M.elem[1][2];
    const float c19 = M.elem[1][0] * M.elem[2][2] - M.elem[2][0] * M.elem[1][2];
    
    const float c20 = M.elem[2][0] * M.elem[3][1] - M.elem[3][0] * M.elem[2][1];
    const float c22 = M.elem[1][0] * M.elem[3][1] - M.elem[3][0] * M.elem[1][1];
    const float c23 = M.elem[1][0] * M.elem[2][1] - M.elem[2][0] * M.elem[1][1];
    
    const vec4 f0 = {c00, c00, c02, c03};
    const vec4 f1 = {c04, c04, c06, c07};
    const vec4 f2 = {c08, c08, c10, c11};
    const vec4 f3 = {c12, c12, c14, c15};
    const vec4 f4 = {c16, c16, c18, c19};
    const vec4 f5 = {c20, c20, c22, c23};
    
    const vec4 v0 = {M.elem[1][0], M.elem[0][0], M.elem[0][0], M.elem[0][0]};
    const vec4 v1 = {M.elem[1][1], M.elem[0][1], M.elem[0][1], M.elem[0][1]};
    const vec4 v2 = {M.elem[1][2], M.elem[0][2], M.elem[0][2], M.elem[0][2]};
    const vec4 v3 = {M.elem[1][3], M.elem[0][3], M.elem[0][3], M.elem[0][3]};
    
    const vec4 i0 = vec4_add(vec4_sub(vec4_mul(v1, f0), vec4_mul(v2, f1)), vec4_mul(v3, f2));
    const vec4 i1 = vec4_add(vec4_sub(vec4_mul(v0, f0), vec4_mul(v2, f3)), vec4_mul(v3, f4));
    const vec4 i2 = vec4_add(vec4_sub(vec4_mul(v0, f1), vec4_mul(v1, f3)), vec4_mul(v3, f5));
    const vec4 i3 = vec4_add(vec4_sub(vec4_mul(v0, f2), vec4_mul(v1, f4)), vec4_mul(v2, f5));
    
    const vec4 sign_a = {+1, -1, +1, -1};
    const vec4 sign_b = {-1, +1, -1, +1};
    
    mat4 I;
    I.col[0] = vec4_mul(i0, sign_a);
    I.col[1] = vec4_mul(i1, sign_b);
    I.col[2] = vec4_mul(i2, sign_a);
    I.col[3] = vec4_mul(i3, sign_b);
    
    const vec4 row0 = {I.elem[0][0], I.elem[1][0], I.elem[2][0], I.elem[3][0]};
    const vec4 dot0 = vec4_mul(M.col[0], row0);
    
    return mat4_mul_f(I, 1.0f / (dot0.x + dot0.y + dot0.z + dot0.w));
}

internal inline mat4 mat4_transpose(const mat4 M) {
    mat4 T;
    T.elem[0][0] = M.elem[0][0];
    T.elem[0][1] = M.elem[1][0];
    T.elem[0][2] = M.elem[2][0];
    T.elem[0][3] = M.elem[3][0];

    T.elem[1][0] = M.elem[0][1];
    T.elem[1][1] = M.elem[1][1];
    T.elem[1][2] = M.elem[2][1];
    T.elem[1][3] = M.elem[3][1];

    T.elem[2][0] = M.elem[0][2];
    T.elem[2][1] = M.elem[1][2];
    T.elem[2][2] = M.elem[2][2];
    T.elem[2][3] = M.elem[3][2];

    T.elem[3][0] = M.elem[0][3];
    T.elem[3][1] = M.elem[1][3];
    T.elem[3][2] = M.elem[2][3];
    T.elem[3][3] = M.elem[3][3];
    
    return T;
}

internal inline bool is_ortho_proj_matrix(const mat4 M) { return M.elem[2][3] == 0.0f; }

internal inline void extract_jitter_uv(float jitter[2], const mat4 M) {
    if (is_ortho_proj_matrix(M)) {
        jitter[0] = -M.elem[3][0] * 0.5f;
        jitter[1] = -M.elem[3][1] * 0.5f;
    }
    else {
        jitter[0] = M.elem[2][0] * 0.5f;
        jitter[1] = M.elem[2][1] * 0.5f;
    }
}

internal inline uint32_t compute_mip_count(uint32_t width, uint32_t height) {
    uint32_t mips = 1;
    uint32_t max_dim = MAX(width, height);
    while (max_dim >>= 1) ++mips;
    return mips;
}

typedef struct gl_aabb {
    float min[3];
    float max[3];
} gl_aabb;

typedef struct gl_buffer {
    GLuint id;
} gl_buffer;

typedef struct gl_texture {
    GLuint id;
} gl_texture;

typedef struct gl_program {
    GLuint id;
} gl_program;

typedef struct gl_backbone_data {
    uint32_t residue_idx;               // @NOTE: the residue index which the backbone segment originates from
    uint32_t segment_idx;               // @NOTE: the segment index for the backbone, local to each chain -> 0 to N-1 where N is the number of residues per chain
    uint32_t ca_idx, c_idx, o_idx;
} gl_backbone_data;

typedef struct gl_control_point {
    float position[3];
    uint32_t atom_idx;
    float velocity[3];
    float segment_t;                    // @NOTE: stores the segment index (integer part) and the fraction within the segment (fractional part)
    uint8_t secondary_structure[3];     // @NOTE: Secondary structure as fractions within the components (0 = coil, 1 = helix, 2 = sheet) so we later can smoothly blend between them when subdividing.
    uint8_t flags;                      // @NOTE: Bitfield for setting flags, bits: [0] beg_chain, [1] end_chain, [2] beg_secondary_structure, [3] end_secondary_structure
    int16_t support_vector[3];
    int16_t tangent_vector[3];
} gl_control_point;

typedef struct gl_view_transform {
    mat4 world_to_view;
    mat4 world_to_view_normal;
    mat4 world_to_clip;
    mat4 view_to_clip;
    mat4 view_to_world;
    mat4 clip_to_view;
    mat4 prev_world_to_clip;
    mat4 curr_view_to_prev_clip;
} gl_view_transform;

// Shared ubo data for all shaders
typedef struct gl_ubo_base {
    gl_view_transform view_transform;
    vec4 jitter_uv;
    uint32_t atom_mask;
} gl_ubo_base;

static const int ubo_base_size = sizeof(gl_ubo_base);

typedef uint32_t gl_version;

static const int control_point_size = sizeof(gl_control_point);

typedef struct internal_rep internal_rep;
typedef struct internal_mol internal_mol;
typedef struct internal_ctx internal_ctx;

struct internal_rep {
    const internal_mol* mol;
    mold_draw_rep_args  args;
    mold_draw_rep_type  type;
    uint16_t            mol_version;
    
    //mold_range visible_range;      // Minimal spanning visible range of atoms to draw
    gl_buffer color;
};

struct internal_mol {
    gl_buffer buffer[GL_BUFFER_MOL_COUNT];
    uint16_t version;

    uint32_t atom_count;
    uint32_t residue_count;
    uint32_t bond_count;
    uint32_t backbone_data_count;    
    uint32_t backbone_control_point_index_count;
    uint32_t backbone_spline_index_count;
};

struct internal_ctx {
   //  gl_view_transform view_transform[GL_VIEW_TRANSFORM_COUNT];

    gl_program program[GL_PROGRAM_COUNT];
    gl_texture texture[GL_TEXTURE_COUNT];
    gl_buffer  ubo;

    GLuint vao;
    GLuint fbo[2];

    gl_version version;
};

static const int ctx_size = sizeof(struct internal_ctx);
static const int mol_size = sizeof(struct internal_mol);
static const int rep_size = sizeof(struct internal_rep);

static_assert(sizeof(internal_ctx) <= sizeof(mold_draw_ctx), "Internal draw ctx does not fit into containing structure");
static_assert(sizeof(internal_mol) <= sizeof(mold_draw_mol), "Internal draw mol does not fit into containing structure");
static_assert(sizeof(internal_rep) <= sizeof(mold_draw_rep), "Internal draw rep does not fit into containing structure");

internal inline gl_buffer gl_buffer_create(uint32_t num_bytes, const void* data, GLenum usage_hint) {
    gl_buffer buf = {0};
    glGenBuffers(1, &buf.id);
    glBindBuffer(GL_ARRAY_BUFFER, buf.id);
    glBufferData(GL_ARRAY_BUFFER, num_bytes, data, usage_hint);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    return buf;
}

internal inline uint32_t gl_buffer_size(gl_buffer buf) {
    GLint size;
    glBindBuffer          (GL_ARRAY_BUFFER, buf.id);
    glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
    glBindBuffer          (GL_ARRAY_BUFFER, 0);
    return (uint32_t)size;
}

internal inline void gl_buffer_conditional_delete(gl_buffer* buf) {
    if (buf->id) {
        GLuint id = buf->id;
        glDeleteBuffers(1, &id);
        buf->id = 0;
    }
}

internal inline void gl_buffer_set_sub_data(gl_buffer buf, uint32_t byte_offset, uint32_t byte_size, const void* data) {
    glBindBuffer(GL_COPY_WRITE_BUFFER, buf.id);
    glBufferSubData(GL_COPY_WRITE_BUFFER, byte_offset, byte_size, data);
    glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
}

internal inline uint32_t gl_texture_width(gl_texture tex, uint32_t level) {
    GLint width;
    glBindTexture(GL_TEXTURE_2D, tex.id);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, level, GL_TEXTURE_WIDTH, &width);
    glBindTexture(GL_TEXTURE_2D, 0);
    return (uint32_t)width;
}

internal inline uint32_t gl_texture_height(gl_texture tex, uint32_t level) {
    GLint height;
    glBindTexture(GL_TEXTURE_2D, tex.id);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, level, GL_TEXTURE_HEIGHT, &height);
    glBindTexture(GL_TEXTURE_2D, 0);
    return (uint32_t)height;
}

internal mold_error compile_shader_from_source(GLuint shader, const char* shader_src, const char* defines, char* err_buf, uint32_t err_len) {
    if (defines) {
        const char* src = shader_src;
        char version_str[32] = {0};
        
        if (strncmp(shader_src, "#version", 8) == 0) {
            const char* endline;
            if ((endline = strchr(shader_src, '\n')) != NULL) {
                ++endline;
                const size_t length = endline - shader_src;
                strncpy(version_str, shader_src, length);
            }
            src = endline;
        }
        const char* sources[5];
        sources[0] = version_str;
        sources[1] = "\n";
        sources[2] = defines;
        sources[3] = "\n";
        sources[4] = src;
        glShaderSource(shader, 5, sources, 0);
    } else {
        const char* c_src = shader_src;
        glShaderSource(shader, 1, &c_src, 0);
    }
    
    GLint success;
    glCompileShader(shader);
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(shader, err_len, NULL, err_buf);
        return MOLD_DRAW_GL_SHADER_COMPILE_ERROR;
    }
    
    return MOLD_DRAW_SUCCESS;
}

internal mold_error compile_shader_from_file(GLuint shader, const char* filename, const char* defines) {
    FILE* file = fopen(filename, "rb");
    if (file) {
        char buffer[MOLD_SHADER_BUF_SIZE] = {0};
        fseek(file, 0, SEEK_END);
        int size = MIN(ftell(file), ARRAY_SIZE(buffer) - 1);
        rewind(file);
        fread(buffer, 1, size, file);
        char err_buf[256];
        mold_error err = compile_shader_from_source(shader, buffer, defines, err_buf, ARRAY_SIZE(err_buf));
        if (err != MOLD_DRAW_SUCCESS) {
            snprintf(error_buf, ARRAY_SIZE(error_buf), "Error compiling file: '%s':\n%s", filename, err_buf);
        }
        fclose(file);
        return err;
    } else {
        snprintf(error_buf, ARRAY_SIZE(error_buf), "Could not open file file '%s'", filename);
        return MOLD_DRAW_FILE_NOT_FOUND;
    }
}

internal mold_error link_program(GLuint program, const GLuint shader[], uint32_t count) {
    ASSERT(program);
     
    for (uint32_t i = 0; i < count; i++) {
        GLint compile_status;
        glGetShaderiv(shader[i], GL_COMPILE_STATUS, &compile_status);
        if (glIsShader(shader[i]) && compile_status) {
            glAttachShader(program, shader[i]);
        } else {
            snprintf(error_buf, ARRAY_SIZE(error_buf), "One or more shaders are invalid");
            return MOLD_DRAW_GL_PROGRAM_LINK_ERROR;
        }
    }
    
    GLint success;
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(program, ARRAY_SIZE(error_buf), NULL, error_buf);
        return MOLD_DRAW_GL_PROGRAM_LINK_ERROR;
    }
    
    for (uint32_t i = 0; i < count; i++) {
        glDetachShader(program, shader[i]);
    }
    
    return MOLD_DRAW_SUCCESS;
}

internal mold_error link_program_transform_feedback(GLuint program, const GLuint shader[], uint32_t shader_count, const char* varying[],
                                            uint32_t varying_count, GLenum capture_mode) {
    ASSERT(program);
    ASSERT(capture_mode == GL_INTERLEAVED_ATTRIBS || capture_mode == GL_SEPARATE_ATTRIBS);
    
    for (uint32_t i = 0; i < shader_count; i++) {
        GLint compile_status;
        glGetShaderiv(shader[i], GL_COMPILE_STATUS, &compile_status);
        if (glIsShader(shader[i]) && compile_status) {
            glAttachShader(program, shader[i]);
        } else {
            snprintf(error_buf, ARRAY_SIZE(error_buf), "One or more shaders are invalid");
            return MOLD_DRAW_GL_PROGRAM_LINK_ERROR;
        }
    }
    
    GLint link_status;
    glTransformFeedbackVaryings(program, varying_count, varying, capture_mode);
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &link_status);
    if (!link_status) {
        glGetProgramInfoLog(program, ARRAY_SIZE(error_buf), NULL, error_buf);
        return MOLD_DRAW_GL_PROGRAM_LINK_ERROR;
    }
    
    for (uint32_t i = 0; i < shader_count; i++) {
        glDetachShader(program, shader[i]);
    }
    
    return MOLD_DRAW_SUCCESS;
}

const char* mold_draw_get_error_str() {
    return error_buf;
}

mold_error validate_context(const internal_ctx* ctx) {
    if (!ctx || (ctx->version != GL_VERSION_330 && ctx->version != GL_VERSION_430)) {
        snprintf(error_buf, ARRAY_SIZE(error_buf), "Internal Context has not been initialized");
        return MOLD_DRAW_CONTEXT_NOT_INITIALIZED;
    }
    return MOLD_DRAW_SUCCESS;
}

mold_error mold_draw_mol_set_atom_position(mold_draw_mol* ext_mol, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t byte_stride) {
    if (ext_mol && x && y && z) {
        internal_mol* mol = (internal_mol*)ext_mol;
        if (offset + count > mol->atom_count) {
            return MOLD_DRAW_GL_ATTEMPTING_TO_WRITE_OUT_OF_BUFFER_RANGE;
        }
        byte_stride = MAX(sizeof(float), byte_stride);
        glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
        float* data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (data) {
            for (uint32_t i = offset; i < count; ++i) {
                data[i * 3 + 0] = *(const float*)((const uint8_t*)x + byte_stride * i);
                data[i * 3 + 1] = *(const float*)((const uint8_t*)y + byte_stride * i);
                data[i * 3 + 2] = *(const float*)((const uint8_t*)z + byte_stride * i);
            }
            glUnmapBuffer(GL_ARRAY_BUFFER);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return MOLD_DRAW_SUCCESS;
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        return MOLD_DRAW_GL_COULD_NOT_MAP_BUFFER;
    }
    return MOLD_DRAW_ARGUMENT_IS_NULL;
}

mold_error mold_draw_mol_set_atom_radius(mold_draw_mol* ext_mol, uint32_t offset, uint32_t count, const float* radius, uint32_t byte_stride) {
    if (ext_mol && radius) {
        internal_mol* mol = (internal_mol*)ext_mol;
        byte_stride = MAX(sizeof(float), byte_stride);
        if (byte_stride > sizeof(float)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS].id);
            float* radius_data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (radius_data) {
                for (uint32_t i = 0; i < mol->atom_count; ++i) {
                    radius_data[i] = *(const float*)((const uint8_t*)radius + byte_stride * i);
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return MOLD_DRAW_SUCCESS;
            }
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return MOLD_DRAW_GL_COULD_NOT_MAP_BUFFER;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS], offset * sizeof(float), count * sizeof(float), radius);
            return MOLD_DRAW_SUCCESS;
        }
    }
    return MOLD_DRAW_ARGUMENT_IS_NULL;
}

mold_error mold_draw_mol_set_atom_flags(mold_draw_mol* ext_mol, uint32_t offset, uint32_t count, const uint8_t* flags, uint32_t byte_stride) {
    if (ext_mol && flags) {
        internal_mol* mol = (internal_mol*)ext_mol;
        byte_stride = MAX(sizeof(uint8_t), byte_stride);
        if (byte_stride > sizeof(uint8_t)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);
            uint8_t* data = (uint8_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (data) {
                for (uint32_t i = 0; i < mol->atom_count; ++i) {
                    data[i] = flags[byte_stride * i];
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return MOLD_DRAW_SUCCESS;
            }
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return MOLD_DRAW_GL_COULD_NOT_MAP_BUFFER;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS], offset * sizeof(uint8_t), count * sizeof(uint8_t), flags);
            return MOLD_DRAW_SUCCESS;
        }
    }
    return MOLD_DRAW_ARGUMENT_IS_NULL;
}

mold_error mold_draw_mol_set_residue_secondary_structure(mold_draw_mol* ext_mol, uint32_t offset, uint32_t count, const mold_secondary_structure* secondary_structure, uint32_t byte_stride) {
    if (ext_mol && secondary_structure) {
        internal_mol* mol = (internal_mol*)ext_mol;
        byte_stride = MAX(sizeof(mold_secondary_structure), byte_stride);
        if (byte_stride > sizeof(mold_secondary_structure)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_RESIDUE_SECONDARY_STRUCTURE].id);
            mold_secondary_structure* data = (mold_secondary_structure*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (data) {
                for (uint32_t i = 0; i < mol->atom_count; ++i) {
                    data[i] = secondary_structure[byte_stride * i];
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return MOLD_DRAW_SUCCESS;
            }
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return MOLD_DRAW_GL_COULD_NOT_MAP_BUFFER;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_RESIDUE_SECONDARY_STRUCTURE], offset * sizeof(mold_secondary_structure), count * sizeof(mold_secondary_structure), secondary_structure);
            return MOLD_DRAW_SUCCESS;
        }
    }
    return MOLD_DRAW_ARGUMENT_IS_NULL;
}

/*
internal mold_error mold_update_visible_atom_color_range(internal_rep* rep) {
    if (!rep) return MOLD_DRAW_REP_INVALID;
    if (!rep->mol) return MOLD_DRAW_MOL_INVALID;

    glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
    void* ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);
    if (ptr) {
        uint32_t* col = (uint32_t*)ptr;
        mold_range range = {0xFFFFFFFF, 0};
            
        for (uint32_t i = 0; i < rep->mol->atom_count; i++) {
            if (col[i] & 0xFF000000) {
                range.beg = MIN(range.beg, i);
                range.end = MAX(range.end, i);
            }
        }
            
        if (range.beg == 0xFFFFFFFFU) {
            range.beg = range.end = 0;
        } else {
            range.end++;    // @NOTE: Increment end one past last found index
        }
            
        rep->visible_range = range;
        glUnmapBuffer(GL_ARRAY_BUFFER);
    } else {
        return MOLD_DRAW_GL_COULD_NOT_MAP_BUFFER;
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    return MOLD_DRAW_SUCCESS;
}
*/

mold_error mold_draw_rep_set_type_and_args(mold_draw_rep* ext_rep, mold_draw_rep_type type, mold_draw_rep_args args) {
    internal_rep* rep = (internal_rep*)ext_rep;
    rep->type = type;
    rep->args = args;
    return MOLD_DRAW_SUCCESS;
}

mold_error mold_draw_rep_set_color(mold_draw_rep* ext_rep, uint32_t offset, uint32_t count, const uint32_t* color_data, uint32_t byte_stride) {
    if (ext_rep && color_data) {
        internal_rep* rep = (internal_rep*)ext_rep;
        if (rep->mol && glIsBuffer(rep->color.id)) {
            if (offset + count > rep->mol->atom_count) {
                return MOLD_DRAW_GL_ATTEMPTING_TO_WRITE_OUT_OF_BUFFER_RANGE;
            }

            if (byte_stride) {
                glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
                uint32_t* data = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
                if (data) {
                    for (uint32_t i = offset; i < offset + count; ++i) {
                        data[i] = *(uint32_t*)((uint8_t*)color_data + i * byte_stride);
                    }
                    glUnmapBuffer(GL_ARRAY_BUFFER);
                    glBindBuffer(GL_ARRAY_BUFFER, 0);
                    return MOLD_DRAW_SUCCESS;
                }
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return MOLD_DRAW_GL_COULD_NOT_MAP_BUFFER;
            }
            else {
                gl_buffer_set_sub_data(rep->color, offset * sizeof(uint32_t), count * sizeof(uint32_t), color_data);
            }
            //mold_update_visible_atom_color_range(rep);
        } else {
            return MOLD_DRAW_GL_INVALID_BUFFER;
        }
    }
    return MOLD_DRAW_ARGUMENT_IS_NULL;
}

mold_error mold_draw_ctx_init(mold_draw_ctx* ext_ctx) {
    if (!ext_ctx) {
        snprintf(error_buf, ARRAY_SIZE(error_buf), "Context pointer is NULL");
        return MOLD_DRAW_ARGUMENT_IS_NULL;
    }

    internal_ctx* ctx = (internal_ctx*)ext_ctx;
    memset(ctx, 0, sizeof(internal_ctx));

    gl3wInit();

    GLint major, minor;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);

    if (major >= 4 && minor >= 3) {
        ctx->version = GL_VERSION_430;
    } else if (major >= 3 && minor >= 3) {
        ctx->version = GL_VERSION_330;
    } else {
        ctx->version = GL_VERSION_UNKNOWN;
        return MOLD_DRAW_CONTEXT_NOT_INITIALIZED;
    }

    glGenVertexArrays(1, &ctx->vao);
    glGenFramebuffers(2, ctx->fbo);
    ctx->ubo = gl_buffer_create(MOLD_UBO_SIZE, NULL, GL_DYNAMIC_DRAW);

    for (uint32_t i = 0; i < GL_TEXTURE_COUNT; ++i) {
        glGenTextures(1, &ctx->texture[i].id);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
        GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);
        GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

        mold_error err;
        if (((err = compile_shader_from_file(vert_shader, MOLD_SHADER_DIR "/space_fill.vert", NULL)) != MOLD_DRAW_SUCCESS) ||
            ((err = compile_shader_from_file(geom_shader, MOLD_SHADER_DIR "/space_fill.geom", NULL)) != MOLD_DRAW_SUCCESS) ||
            ((err = compile_shader_from_file(frag_shader, MOLD_SHADER_DIR "/space_fill.frag", NULL)) != MOLD_DRAW_SUCCESS)) {
            return err;
        }

        ctx->program[GL_PROGRAM_SPACE_FILL].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, geom_shader, frag_shader};
        if ((err = link_program(ctx->program[GL_PROGRAM_SPACE_FILL].id, shaders, ARRAY_SIZE(shaders))) != MOLD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
        glDeleteShader(geom_shader);
        glDeleteShader(frag_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
        GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);
        GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

        mold_error err;
        if (((err = compile_shader_from_file(vert_shader, MOLD_SHADER_DIR "/licorice.vert", NULL)) != MOLD_DRAW_SUCCESS) ||
            ((err = compile_shader_from_file(geom_shader, MOLD_SHADER_DIR "/licorice.geom", NULL)) != MOLD_DRAW_SUCCESS) ||
            ((err = compile_shader_from_file(frag_shader, MOLD_SHADER_DIR "/licorice.frag", NULL)) != MOLD_DRAW_SUCCESS)) {
            return err;
        }

        ctx->program[GL_PROGRAM_LICORICE].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, geom_shader, frag_shader};
        if ((err = link_program(ctx->program[GL_PROGRAM_LICORICE].id, shaders, ARRAY_SIZE(shaders))) != MOLD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
        glDeleteShader(geom_shader);
        glDeleteShader(frag_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
        GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);
        GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

        mold_error err;
        if (((err = compile_shader_from_file(vert_shader, MOLD_SHADER_DIR "/ribbons.vert", NULL)) != MOLD_DRAW_SUCCESS) ||
            ((err = compile_shader_from_file(geom_shader, MOLD_SHADER_DIR "/ribbons.geom", NULL)) != MOLD_DRAW_SUCCESS) ||
            ((err = compile_shader_from_file(frag_shader, MOLD_SHADER_DIR "/ribbons.frag", NULL)) != MOLD_DRAW_SUCCESS)) {
            return err;
        }

        ctx->program[GL_PROGRAM_RIBBONS].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, geom_shader, frag_shader};
        if ((err = link_program(ctx->program[GL_PROGRAM_RIBBONS].id, shaders, ARRAY_SIZE(shaders))) != MOLD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
        glDeleteShader(geom_shader);
        glDeleteShader(frag_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
        GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);
        GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

        mold_error err;
        if (((err = compile_shader_from_file(vert_shader, MOLD_SHADER_DIR "/cartoon.vert", NULL)) != MOLD_DRAW_SUCCESS) ||
            ((err = compile_shader_from_file(geom_shader, MOLD_SHADER_DIR "/cartoon.geom", NULL)) != MOLD_DRAW_SUCCESS) ||
            ((err = compile_shader_from_file(frag_shader, MOLD_SHADER_DIR "/cartoon.frag", NULL)) != MOLD_DRAW_SUCCESS)) {
            return err;
        }

        ctx->program[GL_PROGRAM_CARTOON].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, geom_shader, frag_shader};
        if ((err = link_program(ctx->program[GL_PROGRAM_CARTOON].id, shaders, ARRAY_SIZE(shaders))) != MOLD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
        glDeleteShader(geom_shader);
        glDeleteShader(frag_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);

        mold_error err;
        if ((err = compile_shader_from_file(vert_shader, MOLD_SHADER_DIR "/compute_aabb.vert", NULL)) != MOLD_DRAW_SUCCESS) {
            return err;
        }
        ctx->program[GL_PROGRAM_COMPUTE_AABB_330].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader};
        const GLchar* varyings[] = {"out_aabb_min", "out_aabb_max"};
        if ((err = link_program_transform_feedback(ctx->program[GL_PROGRAM_COMPUTE_AABB_330].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) != MOLD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);

        mold_error err;
        if ((err = compile_shader_from_file(vert_shader, MOLD_SHADER_DIR "/cull_aabb.vert", NULL)) != MOLD_DRAW_SUCCESS) {
            return err;
        }
        ctx->program[GL_PROGRAM_CULL_AABB].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader};
        const GLchar* varyings[] = {"out_visible"};
        if ((err = link_program_transform_feedback(ctx->program[GL_PROGRAM_CULL_AABB].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) != MOLD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);

        mold_error err;
        if ((err = compile_shader_from_file(vert_shader, MOLD_SHADER_DIR "/extract_control_points.vert", NULL)) != MOLD_DRAW_SUCCESS) {
            return err;
        }
        ctx->program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader};
        const GLchar* varyings[] = {"out_position", "out_atom_idx", "out_velocity", "out_segment_t", "out_secondary_structure_and_flags", "out_support_and_tangent_vector"};
        if ((err = link_program_transform_feedback(ctx->program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) != MOLD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
        GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);

        char defines[32];
        snprintf(defines, ARRAY_SIZE(defines), "#define NUM_SUBDIVISIONS %i", MOLD_SPLINE_SUBDIVISION_COUNT);

        mold_error err;
        if ((err = compile_shader_from_file(vert_shader, MOLD_SHADER_DIR "/compute_spline.vert", defines)) != MOLD_DRAW_SUCCESS ||
            (err = compile_shader_from_file(geom_shader, MOLD_SHADER_DIR "/compute_spline.geom", defines)) != MOLD_DRAW_SUCCESS) {
            return err;
        }
        ctx->program[GL_PROGRAM_COMPUTE_SPLINE].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, geom_shader};
        const GLchar* varyings[] = {"out_position", "out_atom_idx", "out_velocity", "out_segment_t", "out_secondary_structure_and_flags", "out_support_and_tangent_vector"};
        if ((err = link_program_transform_feedback(ctx->program[GL_PROGRAM_COMPUTE_SPLINE].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) != MOLD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
        glDeleteShader(geom_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
        GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

        mold_error err;
        if ((err = compile_shader_from_file(vert_shader, MOLD_SHADER_DIR "/fs_quad.vert", NULL)) != MOLD_DRAW_SUCCESS || 
            (err = compile_shader_from_file(frag_shader, MOLD_SHADER_DIR "/downsample_depth.frag", NULL)) != MOLD_DRAW_SUCCESS) {
            return err;
        }
        ctx->program[GL_PROGRAM_DOWNSAMPLE_DEPTH].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, frag_shader};
        if ((err = link_program(ctx->program[GL_PROGRAM_DOWNSAMPLE_DEPTH].id, shaders, ARRAY_SIZE(shaders))) != MOLD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
        glDeleteShader(frag_shader);
    }

    return MOLD_DRAW_SUCCESS;
}

mold_error mold_draw_ctx_free(mold_draw_ctx* ext_ctx) {
    if (ext_ctx) {
        internal_ctx* ctx = (internal_ctx*)ext_ctx;
        if (ctx->vao) glDeleteVertexArrays(1, &ctx->vao);
        if (ctx->fbo) glDeleteFramebuffers(2, ctx->fbo);
        if (ctx->ubo.id) glDeleteBuffers(1, &ctx->ubo.id);
        for (uint32_t i = 0; i < GL_TEXTURE_COUNT; ++i) {
            if (ctx->texture[i].id) glDeleteTextures(1, &ctx->texture[i].id);
        }
        for (uint32_t i = 0; i < GL_PROGRAM_COUNT; ++i) {
            if (ctx->program[i].id) glDeleteProgram(ctx->program[i].id);
        }
        return MOLD_DRAW_SUCCESS;
    }
    return MOLD_DRAW_ARGUMENT_IS_NULL;
}

void mold_draw_mol_desc_extract(mold_draw_mol_desc* desc, const mold_molecule* mold_mol) {
    desc->atom.count     = mold_mol->atom.count;
    desc->atom.x         = mold_mol->atom.x;
    desc->atom.y         = mold_mol->atom.y;
    desc->atom.z         = mold_mol->atom.z;
    desc->atom.radius    = mold_mol->atom.radius;
    
    desc->residue.count               = mold_mol->residue.count;
    desc->residue.atom_range          = mold_mol->residue.atom_range;
    desc->residue.secondary_structure = mold_mol->residue.secondary_structure;
    
    desc->bond.count     = mold_mol->bond.count;
    desc->bond.atom_bond = mold_mol->bond.bond;
}

mold_error mold_draw_mol_init(mold_draw_mol* ext_mol, const mold_draw_mol_desc* desc) {
    if (ext_mol && desc) {
        mold_draw_mol_free(ext_mol);
        internal_mol* mol = (internal_mol*)ext_mol;

        mol->atom_count = desc->atom.count;
        mol->buffer[GL_BUFFER_MOL_ATOM_POSITION]      = gl_buffer_create(mol->atom_count * sizeof(float) * 3, NULL, GL_DYNAMIC_DRAW);
        mol->buffer[GL_BUFFER_MOL_ATOM_PREV_POSITION] = gl_buffer_create(mol->atom_count * sizeof(float) * 3, NULL, GL_DYNAMIC_COPY);
        mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS]        = gl_buffer_create(mol->atom_count * sizeof(float) * 1, NULL, GL_STATIC_DRAW);
        mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS]         = gl_buffer_create(mol->atom_count * sizeof(uint8_t),   NULL, GL_STATIC_DRAW);

        if (desc->atom.x && desc->atom.y && desc->atom.z)   mold_draw_mol_set_atom_position(ext_mol, 0, mol->atom_count, desc->atom.x, desc->atom.y, desc->atom.z, 0);
        if (desc->atom.radius)                              mold_draw_mol_set_atom_radius(ext_mol, 0, mol->atom_count, desc->atom.radius, 0);
        if (desc->atom.flags)                               mold_draw_mol_set_atom_flags(ext_mol, 0, mol->atom_count, desc->atom.flags, 0);

        mol->residue_count = desc->residue.count;
        mol->buffer[GL_BUFFER_MOL_RESIDUE_ATOM_RANGE]          = gl_buffer_create(mol->residue_count * sizeof(uint32_t) * 2, NULL, GL_STATIC_DRAW);
        mol->buffer[GL_BUFFER_MOL_RESIDUE_AABB]                = gl_buffer_create(mol->residue_count * sizeof(float) * 6,    NULL, GL_DYNAMIC_COPY);
        mol->buffer[GL_BUFFER_MOL_RESIDUE_VISIBLE]             = gl_buffer_create(mol->residue_count * sizeof(int),          NULL, GL_DYNAMIC_COPY);
        mol->buffer[GL_BUFFER_MOL_RESIDUE_SECONDARY_STRUCTURE] = gl_buffer_create(mol->residue_count * sizeof(uint8_t) * 4,  NULL, GL_DYNAMIC_DRAW);

        if (desc->residue.atom_range)           gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_RESIDUE_ATOM_RANGE], 0, mol->residue_count * sizeof(uint32_t) * 2, desc->residue.atom_range);
        //if (desc->residue.backbone_atoms)       gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_RESIDUE_BACKBONE_ATOMS], 0, mol->residue_count * sizeof(uint8_t) * 4, desc->residue.backbone_atoms);
        if (desc->residue.secondary_structure)  gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_RESIDUE_SECONDARY_STRUCTURE], 0, mol->residue_count * sizeof(uint8_t) * 4, desc->residue.secondary_structure);

        if (desc->chain.count > 0 && desc->chain.residue_range && desc->residue.backbone_atoms && desc->residue.atom_range && desc->residue.secondary_structure) {
            uint32_t backbone_residue_count = 0;
            uint32_t backbone_spline_count = 0;
            for (uint32_t i = 0; i < desc->chain.count; ++i) {
                uint32_t res_count = desc->chain.residue_range[i].end - desc->chain.residue_range[i].beg;
                backbone_residue_count += res_count;
                backbone_spline_count += (res_count - 1) * MOLD_SPLINE_SUBDIVISION_COUNT;
            }

            const uint32_t backbone_data_count                = backbone_residue_count;
            const uint32_t backbone_control_point_data_count  = backbone_residue_count;
            const uint32_t backbone_control_point_index_count = backbone_residue_count + desc->chain.count * (2 + 1); // Duplicate pair first and last in each chain for adjacency + primitive restart between
            const uint32_t backbone_spline_data_count         = backbone_spline_count;
            const uint32_t backbone_spline_index_count        = backbone_spline_count + desc->chain.count * (1); // primitive restart between each chain

            mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA]                = gl_buffer_create(backbone_data_count                * sizeof(gl_backbone_data), NULL, GL_STATIC_DRAW);
            mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_DATA]  = gl_buffer_create(backbone_control_point_data_count  * sizeof(gl_control_point), NULL, GL_DYNAMIC_COPY);
            mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_INDEX] = gl_buffer_create(backbone_control_point_index_count * sizeof(uint32_t),         NULL, GL_STATIC_DRAW);
            mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA]         = gl_buffer_create(backbone_spline_data_count         * sizeof(gl_control_point), NULL, GL_DYNAMIC_COPY);
            mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX]        = gl_buffer_create(backbone_spline_index_count        * sizeof(uint32_t),         NULL, GL_STATIC_DRAW);

            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA].id);
            gl_backbone_data* backbone_data = (gl_backbone_data*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (backbone_data) {
                uint32_t idx = 0;
                for (uint32_t i = 0; i < desc->chain.count; ++i) {
                    for (uint32_t j = desc->chain.residue_range[i].beg; j < desc->chain.residue_range[i].end; ++j) {
                        backbone_data[idx].residue_idx = j;
                        backbone_data[idx].segment_idx = j - desc->chain.residue_range[i].beg;
                        backbone_data[idx].ca_idx = desc->residue.backbone_atoms[j].ca;
                        backbone_data[idx].c_idx  = desc->residue.backbone_atoms[j].c;
                        backbone_data[idx].o_idx  = desc->residue.backbone_atoms[j].o;
                        ++idx;
                    }
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                return MOLD_DRAW_GL_COULD_NOT_MAP_BUFFER;
            }

            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_INDEX].id);
            uint32_t* control_point_index = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (control_point_index) {
                uint32_t idx = 0;
                uint32_t len = 0;
                for (uint32_t i = 0; i < desc->chain.count; ++i) {
                    control_point_index[len++] = idx;
                    for (uint32_t j = desc->chain.residue_range[i].beg; j < desc->chain.residue_range[i].end; ++j) {
                        control_point_index[len++] = idx++;
                    }
                    control_point_index[len++] = idx-1;
                    control_point_index[len++] = 0xFFFFFFFF;
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            }
            else {
                return MOLD_DRAW_GL_COULD_NOT_MAP_BUFFER;
            }

            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX].id);
            uint32_t* spline_index = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (spline_index) {
                uint32_t idx = 0;
                uint32_t len = 0;
                for (uint32_t i = 0; i < desc->chain.count; ++i) {
                    uint32_t res_count = desc->chain.residue_range[i].end - desc->chain.residue_range[i].beg;
                    if (res_count > 0) {
                        for (uint32_t j = 0; j < (res_count - 1) * MOLD_SPLINE_SUBDIVISION_COUNT; ++j) {
                            spline_index[len++] = idx++;
                        }
                    }
                    spline_index[len++] = 0xFFFFFFFF;
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                return MOLD_DRAW_GL_COULD_NOT_MAP_BUFFER;
            }

            glBindBuffer(GL_ARRAY_BUFFER, 0);

            mol->backbone_control_point_index_count = backbone_control_point_index_count;
            mol->backbone_spline_index_count = backbone_spline_index_count;
            mol->backbone_data_count = backbone_data_count;
        }

        mol->bond_count = desc->bond.count;
        mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES] = gl_buffer_create(mol->bond_count * sizeof(uint32_t) * 2, NULL, GL_DYNAMIC_COPY);

        if (desc->bond.atom_bond) gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES], 0, mol->bond_count * sizeof(uint32_t) * 2, desc->bond.atom_bond);
       
        return MOLD_DRAW_SUCCESS;
    }
    return MOLD_DRAW_ARGUMENT_IS_NULL;
}

mold_error mold_draw_mol_free(mold_draw_mol* ext_mol) {
    if (ext_mol) {
        internal_mol* mol = (internal_mol*)ext_mol;
        for (uint32_t i = 0; i < GL_BUFFER_MOL_COUNT; ++i) {
            gl_buffer_conditional_delete(&mol->buffer[i]);
        }
        uint16_t new_version = mol->version + 1;
        memset(mol, 0, sizeof(mold_draw_mol));
        mol->version = new_version;
        return MOLD_DRAW_SUCCESS;
    }
    return MOLD_DRAW_ARGUMENT_IS_NULL;
}

mold_error mold_draw_rep_init(mold_draw_rep* ext_rep, const mold_draw_mol* ext_mol) {   
    if (ext_rep && ext_mol) {
        mold_draw_rep_free(ext_rep);
        internal_rep* rep = (internal_rep*)ext_rep;
        memset(rep, 0, sizeof(internal_rep));
        rep->mol = (const internal_mol*)ext_mol;
        rep->mol_version = rep->mol->version;
        rep->color = gl_buffer_create(rep->mol->atom_count * sizeof(uint32_t), NULL, GL_STATIC_DRAW);
        const uint32_t color = (uint32_t)((255 << 24) | (127 << 16) | (127 << 8) | (127 << 0));
        glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
        uint32_t* data = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (data) {
            for (uint32_t i = 0; i < rep->mol->atom_count; ++i) {
                data[i] = color;
            }
            glUnmapBuffer(GL_ARRAY_BUFFER);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            //mold_update_visible_atom_color_range(rep);
            return MOLD_DRAW_SUCCESS;
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        return MOLD_DRAW_GL_COULD_NOT_MAP_BUFFER;
    }
    return MOLD_DRAW_ARGUMENT_IS_NULL;
}

mold_error mold_draw_rep_free(mold_draw_rep* ext_rep) {
    if (ext_rep) {
        internal_rep* rep = (internal_rep*)ext_rep;
        gl_buffer_conditional_delete(&rep->color);
        return MOLD_DRAW_SUCCESS;
    }
    return MOLD_DRAW_ARGUMENT_IS_NULL;
}

mold_error gl_build_depth_mipmap(const internal_ctx* ctx, gl_texture dst, gl_texture src);
mold_error gl_compute_residue_aabb(const internal_ctx* ctx, const internal_mol* mol);
mold_error gl_cull_aabb(const internal_ctx* ctx, const internal_mol* mol);
mold_error gl_compute_spline(const internal_ctx* ctx, const internal_mol* mol);

mold_error gl_draw_space_fill(const internal_ctx* ctx, const internal_rep* rep);
mold_error gl_draw_licorice(const internal_ctx* ctx, const internal_rep* rep);
mold_error gl_draw_ribbons(const internal_ctx* ctx, const internal_rep* rep);
mold_error gl_draw_cartoon(const internal_ctx* ctx, const internal_rep* rep);

int compare_draw_rep(const void* elem1, const void* elem2) {
    const internal_rep* r1 = (const internal_rep*)elem1;
    const internal_rep* r2 = (const internal_rep*)elem2;
    
    return (r1->type - r2->type) * 2 + (int)(r1->mol - r2->mol);
}

mold_error mold_draw(const mold_draw_ctx* ext_ctx, const mold_draw_desc* desc) {
    const internal_ctx* ctx = (const internal_ctx*)ext_ctx;
    if (!desc) return MOLD_DRAW_ARGUMENT_IS_NULL;
    mold_error err;
    if ((err = validate_context(ctx)) != MOLD_DRAW_SUCCESS) return err;

    PUSH_GPU_SECTION("MOLD DRAW")
            
    GLint bound_fbo;
    GLint bound_viewport[4];
    GLint bound_draw_buffer[8] = {0};
    GLint bound_draw_buffer_count = 0;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &bound_fbo);
    glGetIntegerv(GL_VIEWPORT, bound_viewport);
    for (int i = 0; i < 8; ++i) {
        glGetIntegerv(GL_DRAW_BUFFER0 + i, &bound_draw_buffer[i]);
        // @NOTE: Assume that its tightly packed and if we stumple upon a zero draw buffer index, we enterpret that as the 'end'
        if (bound_draw_buffer[i] == GL_NONE) {
            bound_draw_buffer_count = i;
            break;
        }
    }

    gl_ubo_base ubo_data = {0};
    memcpy(&ubo_data.view_transform.world_to_view, desc->view_transform.model_view_matrix, sizeof(mat4));
    memcpy(&ubo_data.view_transform.view_to_clip,  desc->view_transform.projection_matrix, sizeof(mat4));
    ubo_data.view_transform.world_to_clip = mat4_mul(ubo_data.view_transform.view_to_clip, ubo_data.view_transform.world_to_view);
    ubo_data.view_transform.world_to_view_normal = mat4_transpose(mat4_inverse(ubo_data.view_transform.world_to_view));
    ubo_data.view_transform.view_to_world = mat4_inverse(ubo_data.view_transform.world_to_view);
    ubo_data.view_transform.clip_to_view = mat4_inverse(ubo_data.view_transform.view_to_clip);

    if (desc->view_transform.prev_model_view_matrix && desc->view_transform.projection_matrix) {
        const mat4* prev_world_to_view = (const mat4*)desc->view_transform.prev_model_view_matrix;
        const mat4* prev_view_to_clip  = (const mat4*)desc->view_transform.prev_projection_matrix;
        ubo_data.view_transform.prev_world_to_clip = mat4_mul(*prev_view_to_clip, *prev_world_to_view);
        ubo_data.view_transform.curr_view_to_prev_clip = mat4_mul(ubo_data.view_transform.prev_world_to_clip, ubo_data.view_transform.view_to_world);
        extract_jitter_uv(&ubo_data.jitter_uv.x, ubo_data.view_transform.view_to_clip);
        extract_jitter_uv(&ubo_data.jitter_uv.z, *prev_view_to_clip);
    }
    ubo_data.atom_mask = desc->mol_mask;

    gl_buffer_set_sub_data(ctx->ubo, 0, sizeof(ubo_data), &ubo_data);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ctx->ubo.id);
        
    // Extract most potentially occluding structures first
    // Surface > Space Fill > Ribbons > Cartoon > Licorice
    const internal_rep* draw_rep[64];
    uint32_t            draw_rep_count = 0;
        
    const internal_mol* unique_mol_ptr[64];
    uint32_t            unique_mol_flags[64];
    uint32_t            unique_mol_count = 0;

    const uint32_t MOL_FLAG_SPLINE = 1;
        
    ASSERT(desc->representation.count < ARRAY_SIZE(draw_rep));
    for (uint32_t i = 0; i < desc->representation.count; i++) {
        const internal_rep* rep = (internal_rep*)desc->representation.data[i];
        if (rep && rep->mol) {
            draw_rep[draw_rep_count] = rep;
            draw_rep_count++;
                
            uint32_t unique_mol_idx = (uint32_t)-1;
            for (uint32_t j = 0; j < unique_mol_count; j++) {
                if (unique_mol_ptr[j] == rep->mol) {
                    unique_mol_idx = j;
                    break;
                }
            }
            if (unique_mol_idx == -1) {
                unique_mol_idx = unique_mol_count++;
                unique_mol_ptr[unique_mol_idx] = rep->mol;
                unique_mol_flags[unique_mol_idx] = 0;
            }
            if (rep->type == MOLD_DRAW_REP_RIBBONS || rep->type == MOLD_DRAW_REP_CARTOON) {
                unique_mol_flags[unique_mol_idx] |= MOL_FLAG_SPLINE;
            }
        }
    }
    
    qsort((void*)draw_rep, draw_rep_count, sizeof(mold_draw_rep*), compare_draw_rep);
        
    // Downsample previous frames depth into hi_z texture
        
    // Compute residue AABBs for culling
    /*
    PUSH_GPU_SECTION("COMPUTE RESIDUE AABB")
    for (uint32_t i = 0; i < unique_mol_count; i++) {
        gl_compute_residue_aabb(ctx, unique_mol_ptr[i]);
    }
    POP_GPU_SECTION()
            
    PUSH_GPU_SECTION("CULL RESIDUE AABB")
    for (uint32_t i = 0; i < unique_mol_count; i++) {
        gl_cull_aabb(ctx, unique_mol_ptr[i]);
    }
    POP_GPU_SECTION()
    */
            
    PUSH_GPU_SECTION("COMPUTE SPLINE")
    for (uint32_t i = 0; i < unique_mol_count; i++) {
        if (unique_mol_flags[i] & MOL_FLAG_SPLINE) {
            gl_compute_spline(ctx, unique_mol_ptr[i]);
        }
    }
    POP_GPU_SECTION()
            
    // Possibly resize internal depth buffer to match render_target size
    const gl_texture* max_tex = &ctx->texture[GL_TEXTURE_MAX_DEPTH];
    uint32_t max_tex_w = gl_texture_width(*max_tex, 0);
    uint32_t max_tex_h = gl_texture_height(*max_tex, 0);
    if (max_tex_w != desc->render_target->width || max_tex_h != desc->render_target->height) {
        max_tex_w = desc->render_target->width;
        max_tex_h = desc->render_target->height;
        int mips = compute_mip_count(max_tex_w, max_tex_h);
        if (max_tex->id) glDeleteTextures(1, &max_tex->id);
        glGenTextures(1, &max_tex->id);
        glBindTexture(GL_TEXTURE_2D, max_tex->id);
        glTexStorage2D(GL_TEXTURE_2D, mips, GL_DEPTH24_STENCIL8, max_tex_w, max_tex_h);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
        
    bool using_internal_depth = false;
        
    // If we have rendertargets -> setup and enable fbo
    if (desc->render_target) {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ctx->fbo[0]);
        GLenum draw_buffers[4];
        uint32_t count = 0;
        if (desc->render_target->texture_depth) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, desc->render_target->texture_depth, 0);
        } else {
            using_internal_depth = true;
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, ctx->texture[GL_TEXTURE_MAX_DEPTH].id, 0);
        }
        if (desc->render_target->texture_color) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, desc->render_target->texture_color, 0);
            draw_buffers[count++] = GL_COLOR_ATTACHMENT0;
        }
        if (desc->render_target->texture_view_normal) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, desc->render_target->texture_view_normal, 0);
            draw_buffers[count++] = GL_COLOR_ATTACHMENT1;
        }
        if (desc->render_target->texture_view_velocity) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, desc->render_target->texture_view_velocity, 0);
            draw_buffers[count++] = GL_COLOR_ATTACHMENT2;
        }
        if (desc->render_target->texture_atom_index) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_2D, desc->render_target->texture_atom_index, 0);
            draw_buffers[count++] = GL_COLOR_ATTACHMENT3;
        }
        glDrawBuffers(count, draw_buffers);
        glViewport(0, 0, desc->render_target->width, desc->render_target->height);
            
        GLenum fbo_status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (fbo_status != GL_FRAMEBUFFER_COMPLETE) {
            switch (fbo_status) {
            case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
                snprintf(error_buf, ARRAY_SIZE(error_buf), "Incomplete fbo attachment");
                return MOLD_DRAW_GL_FRAMEBUFFER_ERROR;
            case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
                snprintf(error_buf, ARRAY_SIZE(error_buf), "Missing fbo attachment");
                return MOLD_DRAW_GL_FRAMEBUFFER_ERROR;
            case GL_FRAMEBUFFER_UNSUPPORTED:

                snprintf(error_buf, ARRAY_SIZE(error_buf), "Unsupported fbo attachment format");
                return MOLD_DRAW_GL_FRAMEBUFFER_ERROR;
            default:
                break;
            };
        }
            
        if (using_internal_depth) {
            glClear(GL_DEPTH_BUFFER_BIT);
        }
            
    }
        
    PUSH_GPU_SECTION("DRAW REP")
    {
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glColorMask(1,1,1,1);
        glDepthMask(1);
        
        for (uint32_t i = 0; i < draw_rep_count; i++) {
            switch (draw_rep[i]->type) {
            case MOLD_DRAW_REP_SOLVENT_EXCLUDED_SURFACE:
                break;
            case MOLD_DRAW_REP_DEFAULT:
            case MOLD_DRAW_REP_SPACE_FILL:
                gl_draw_space_fill(ctx, draw_rep[i]);
                break;
            case MOLD_DRAW_REP_RIBBONS:
                gl_draw_ribbons(ctx, draw_rep[i]);
                break;
            case MOLD_DRAW_REP_CARTOON:
                gl_draw_cartoon(ctx, draw_rep[i]);
                break;
            case MOLD_DRAW_REP_LICORICE:
                gl_draw_licorice(ctx, draw_rep[i]);
                break;
            default:
                snprintf(error_buf, ARRAY_SIZE(error_buf), "Representation had unexpected type");
                return MOLD_DRAW_UNKNOWN_ERROR;
                break;
            }
        }
    }
    POP_GPU_SECTION()
       
#if 1
    PUSH_GPU_SECTION("DOWNSAMPLE MAX DEPTH")
    {
        gl_texture dst = ctx->texture[GL_TEXTURE_MAX_DEPTH];
        if (using_internal_depth) {
            gl_build_depth_mipmap(ctx, dst, dst);
        } else {
            ASSERT(desc->render_target);
            gl_texture src = {desc->render_target->texture_depth};
            gl_build_depth_mipmap(ctx, dst, src);
        }
    }
    POP_GPU_SECTION()
#endif
#if 1
    PUSH_GPU_SECTION("COPY POSITION DATA")    
    for (uint32_t i = 0; i < unique_mol_count; i++) {
        glBindBuffer(GL_COPY_READ_BUFFER,  unique_mol_ptr[i]->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
        glBindBuffer(GL_COPY_WRITE_BUFFER, unique_mol_ptr[i]->buffer[GL_BUFFER_MOL_ATOM_PREV_POSITION].id);
        glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, gl_buffer_size(unique_mol_ptr[i]->buffer[GL_BUFFER_MOL_ATOM_POSITION]));
        glBindBuffer(GL_COPY_READ_BUFFER, 0);
        glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
    }
    POP_GPU_SECTION()
#endif
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, bound_fbo);
    glViewport(bound_viewport[0], bound_viewport[1], bound_viewport[2], bound_viewport[3]);
    glDrawBuffers(bound_draw_buffer_count, (GLenum*)bound_draw_buffer);
        
    POP_GPU_SECTION()
    
    return MOLD_DRAW_SUCCESS;
}

mold_error gl_draw_space_fill(const internal_ctx* ctx, const internal_rep* rep) {
    const float radius_scale = rep->args.space_fill.radius_scale;
    gl_buffer_set_sub_data(ctx->ubo, sizeof(gl_ubo_base), sizeof(radius_scale), &radius_scale);

    glBindVertexArray(ctx->vao);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, 0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_PREV_POSITION].id);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, 0);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS].id);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(float), 0);

    glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
    glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(uint32_t), 0);
    
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glUseProgram(ctx->program[GL_PROGRAM_SPACE_FILL].id);
    glDrawArrays(GL_POINTS, 0, rep->mol->atom_count);
    glUseProgram(0);
    
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    
    glBindVertexArray(0);
    
    return MOLD_DRAW_SUCCESS;
}

mold_error gl_draw_licorice(const internal_ctx* ctx, const internal_rep* rep) {
    const float radius = rep->args.licorice.radius;
    gl_buffer_set_sub_data(ctx->ubo, sizeof(gl_ubo_base), sizeof(radius), &radius);

    glBindVertexArray(ctx->vao);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, 0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_PREV_POSITION].id);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, 0);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
    glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(uint32_t), 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES].id);

    glUseProgram(ctx->program[GL_PROGRAM_LICORICE].id);
    glDrawElements(GL_LINES, rep->mol->bond_count * 2, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);

    glBindVertexArray(0);

    return MOLD_DRAW_SUCCESS;
}

mold_error gl_draw_ribbons(const internal_ctx* ctx, const internal_rep* rep) {
    ASSERT(rep->mol);
    if (!rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id) {
        snprintf(error_buf, ARRAY_SIZE(error_buf), "No spline present for molecule, which is required for rendering ribbons representation");
        return MOLD_DRAW_GL_INVALID_BUFFER;
    }

    float scale[2] = {
        rep->args.ribbons.width_scale,
        rep->args.ribbons.thickness_scale * 0.1f
    };
    gl_buffer_set_sub_data(ctx->ubo, sizeof(gl_ubo_base), sizeof(scale), &scale);

    glBindVertexArray(ctx->vao);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, position));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, atom_idx));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, velocity));

    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_BYTE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, flags));

    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 3, GL_SHORT, GL_TRUE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, support_vector));

    glEnableVertexAttribArray(5);
    glVertexAttribPointer(5, 3, GL_SHORT, GL_TRUE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, tangent_vector));

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX].id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx->texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, rep->color.id);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);

    glUseProgram(ctx->program[GL_PROGRAM_RIBBONS].id);
    glDrawElements(GL_LINE_STRIP, rep->mol->backbone_spline_index_count, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);

    glBindVertexArray(0);

    return MOLD_DRAW_SUCCESS;
}

mold_error gl_draw_cartoon(const internal_ctx* ctx, const internal_rep* rep) {
    ASSERT(rep->mol);
    if (!rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id) {
        snprintf(error_buf, ARRAY_SIZE(error_buf), "No spline present for molecule, which is required for rendering ribbons representation");
        return MOLD_DRAW_GL_INVALID_BUFFER;
    }

    float scale[2] = {
        rep->args.ribbons.width_scale,
        rep->args.ribbons.thickness_scale
    };
    gl_buffer_set_sub_data(ctx->ubo, sizeof(gl_ubo_base), sizeof(scale), &scale);

    glBindVertexArray(ctx->vao);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, position));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, atom_idx));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, velocity));

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, segment_t));

    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 3, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, secondary_structure));

    glEnableVertexAttribArray(5);
    glVertexAttribIPointer(5, 1, GL_UNSIGNED_BYTE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, flags));

    glEnableVertexAttribArray(6);
    glVertexAttribPointer(6, 3, GL_SHORT, GL_TRUE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, support_vector));

    glEnableVertexAttribArray(7);
    glVertexAttribPointer(7, 3, GL_SHORT, GL_TRUE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, tangent_vector));

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX].id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx->texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, rep->color.id);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);

    glUseProgram(ctx->program[GL_PROGRAM_CARTOON].id);
    glDrawElements(GL_LINE_STRIP, rep->mol->backbone_spline_index_count, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);
    glDisableVertexAttribArray(5);

    glBindVertexArray(0);

    return MOLD_DRAW_SUCCESS;
}


mold_error gl_build_depth_mipmap(const internal_ctx* ctx, gl_texture dst, gl_texture src) {
    ASSERT(ctx);
    const GLint uniform_loc_level = glGetUniformLocation(ctx->program[GL_PROGRAM_DOWNSAMPLE_DEPTH].id, "u_level");
    const GLint uniform_loc_even  = glGetUniformLocation(ctx->program[GL_PROGRAM_DOWNSAMPLE_DEPTH].id, "u_even");
    
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ctx->fbo[1]);
    glDepthFunc(GL_ALWAYS);
    glActiveTexture(GL_TEXTURE0);
    glDepthMask(1);
    glColorMask(0,0,0,0);
    glUseProgram(ctx->program[GL_PROGRAM_DOWNSAMPLE_DEPTH].id);
    glBindTexture(GL_TEXTURE_2D, dst.id);

    uint32_t src_w = gl_texture_width(src, 0);
    uint32_t src_h = gl_texture_height(src, 0);
    uint32_t dst_w = gl_texture_width(dst, 0);
    uint32_t dst_h = gl_texture_height(dst, 0);
    
    if (dst.id != src.id) {
        // Copy data into src
        glBindFramebuffer(GL_READ_FRAMEBUFFER, ctx->fbo[0]);
        glFramebufferTexture2D(GL_READ_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, src.id, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, dst.id, 0);
        glBlitFramebuffer(0, 0, src_w, src_h, 0, 0, dst_w, dst_h, GL_DEPTH_BUFFER_BIT, GL_NEAREST);
        glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    }
    
    uint32_t level = 0;
    uint32_t dim = MAX(dst_w, dst_h);
    uint32_t width = dst_w;
    uint32_t height = dst_h;
    uint32_t even = 0;
    
    while (dim) {
        if (level) {
            width = MAX(1, width);
            height = MAX(1, height);
            glViewport(0, 0, width, height);
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, dst.id, level);
            glUniform1i(uniform_loc_level, level - 1);
            glUniform1i(uniform_loc_even, even);
            glDrawArrays(GL_TRIANGLES, 0, 3);
        }
        
        even = (width % 2 == 0) && (height % 2 == 0);
        dim /= 2;
        width /= 2;
        height /= 2;
        level++;
    }
    
    glUseProgram(0);
    glDepthMask(1);
    glColorMask(1,1,1,1);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glDepthFunc(GL_LEQUAL);
    glViewport(0, 0, width, height);
    
    return MOLD_DRAW_SUCCESS;
};

mold_error gl_compute_residue_aabb(const internal_ctx* ctx, const internal_mol* mol) {   
    //    if (ctx.version == GL_VERSION_330) {
    // Use transform feedback
    glEnable(GL_RASTERIZER_DISCARD);
    
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx->texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
    
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx->texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS].id);
    
    glBindVertexArray(ctx->vao);
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_RESIDUE_ATOM_RANGE].id);
    glEnableVertexAttribArray(0);
    glVertexAttribIPointer(0, 2, GL_UNSIGNED_INT, sizeof(uint32_t) * 2, 0);
    
    const GLuint prog = ctx->program[GL_PROGRAM_COMPUTE_AABB_330].id;
    glUseProgram(prog);
    
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->buffer[GL_BUFFER_MOL_RESIDUE_VISIBLE].id);
    glUniform1i(glGetUniformLocation(prog, "u_atom_pos_buf"), 0);
    glUniform1i(glGetUniformLocation(prog, "u_atom_rad_buf"), 1);
    
    glBeginTransformFeedback(GL_POINTS);
    glDrawArrays(GL_POINTS, 0, mol->residue_count);
    glEndTransformFeedback();
    
    glUseProgram(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDisableVertexAttribArray(0);
    glBindVertexArray(0);
    glDisable(GL_RASTERIZER_DISCARD);
    /*
        }
        else if (ctx.version == GL_VERSION_430) {
            // Use compute shader
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, out_residue_aabb.id);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, in_residue_atom_range.id);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, in_atom_position.id);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, in_atom_radius.id);

            glUseProgram(ctx.program.compute_aabb_430.id);
            glDispatchCompute(residue_count, 1, 1);
            glUseProgram(0);
        }
      */
    return MOLD_DRAW_SUCCESS;
}

mold_error gl_cull_aabb(const internal_ctx* ctx, const internal_mol* mol) {

#if 1
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_RESIDUE_VISIBLE].id);
    void* ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    if (ptr) {
        uint32_t* data = (uint32_t*)ptr;
        for (uint32_t i = 0; i < mol->residue_count; i++) {
            //data[i] = (i % 2U) ? 1U : 0U;
            data[i] = 1U;
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }
#else
    glEnable(GL_RASTERIZER_DISCARD);
    
    glBindVertexArray(ctx.vao);
    glBindBuffer(GL_ARRAY_BUFFER, mol->residue.aabb.id);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(gl_aabb), offsetof(gl_aabb, min));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(gl_aabb), offsetof(gl_aabb, max));
    
    glBindBuffer(GL_UNIFORM_BUFFER, ctx.ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(float) * 16, mvp_mat);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    
    glUseProgram(ctx.program.cull_aabb.id);
    
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ctx.ubo);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->residue.visible.id);
    
    glBeginTransformFeedback(GL_POINTS);
    glDrawArrays(GL_POINTS, 0, mol->residue_count);
    glEndTransformFeedback();
    
    glUseProgram(0);
    
    glDisable(GL_RASTERIZER_DISCARD);
#endif
    
    return MOLD_DRAW_SUCCESS;
}

mold_error gl_compute_spline(const internal_ctx* ctx, const internal_mol* mol) {
    if (mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA].id == 0) {
        snprintf(error_buf, ARRAY_SIZE(error_buf), "Backbone data buffer is zero, which is required to compute the backbone. Is the molecule missing a backbone?");
        return MOLD_DRAW_GL_INVALID_BUFFER;
    }

    glEnable(GL_RASTERIZER_DISCARD);
    glBindVertexArray(ctx->vao);

    // Step 1: Extract control points
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA].id);

    glEnableVertexAttribArray(0);
    glVertexAttribIPointer(0, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data), (const void*)offsetof(gl_backbone_data, residue_idx));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data), (const void*)offsetof(gl_backbone_data, segment_idx));

    glEnableVertexAttribArray(2);
    glVertexAttribIPointer(2, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data), (const void*)offsetof(gl_backbone_data, ca_idx));

    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data), (const void*)offsetof(gl_backbone_data, c_idx));

    glEnableVertexAttribArray(4);
    glVertexAttribIPointer(4, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data), (const void*)offsetof(gl_backbone_data, o_idx));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx->texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx->texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_MOL_ATOM_PREV_POSITION].id);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, ctx->texture[GL_TEXTURE_BUFFER_2].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, mol->buffer[GL_BUFFER_MOL_RESIDUE_SECONDARY_STRUCTURE].id);

    GLuint program = ctx->program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id;
    const GLint buf_atom_pos_loc            = glGetUniformLocation(program, "u_buf_atom_pos");
    const GLint buf_atom_prev_pos_loc       = glGetUniformLocation(program, "u_buf_atom_prev_pos");
    const GLint buf_secondary_structure_loc = glGetUniformLocation(program, "u_buf_secondary_structure");

    glUseProgram(program);
    glUniform1i(buf_atom_pos_loc, 0);
    glUniform1i(buf_atom_prev_pos_loc, 1);
    glUniform1i(buf_secondary_structure_loc, 2);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_DATA].id);
    glBeginTransformFeedback(GL_POINTS);
    glDrawArrays(GL_POINTS, 0, mol->backbone_data_count);
    glEndTransformFeedback();
    glUseProgram(0);

    // Step 2: Compute splines using control points
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_DATA].id);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,        sizeof(gl_control_point), (const void*)offsetof(gl_control_point, position));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT,          sizeof(gl_control_point), (const void*)offsetof(gl_control_point, atom_idx));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE,        sizeof(gl_control_point), (const void*)offsetof(gl_control_point, velocity));

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE,        sizeof(gl_control_point), (const void*)offsetof(gl_control_point, segment_t));

    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(gl_control_point), (const void*)offsetof(gl_control_point, secondary_structure));

    glEnableVertexAttribArray(5);
    glVertexAttribIPointer(5, 1, GL_UNSIGNED_BYTE,         sizeof(gl_control_point), (const void*)offsetof(gl_control_point, flags));

    glEnableVertexAttribArray(6);
    glVertexAttribPointer(6, 3, GL_SHORT, GL_TRUE,         sizeof(gl_control_point), (const void*)offsetof(gl_control_point, support_vector));
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_INDEX].id);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);

    glUseProgram(ctx->program[GL_PROGRAM_COMPUTE_SPLINE].id);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id);
    glBeginTransformFeedback(GL_POINTS);
    glDrawElements(GL_LINE_STRIP_ADJACENCY, mol->backbone_control_point_index_count, GL_UNSIGNED_INT, 0);
    glEndTransformFeedback();
    glUseProgram(0);

    glDisable(GL_PRIMITIVE_RESTART);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);

    glBindVertexArray(0);
    glDisable(GL_RASTERIZER_DISCARD);

    return MOLD_DRAW_SUCCESS;
}