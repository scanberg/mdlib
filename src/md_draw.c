
#include "md_draw.h"
#include "md_util.h"
#include "ext/gl3w/gl3w.h"
#include <core/md_common.h>
#include <core/md_compiler.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_file.h>

#include <stdbool.h>
#include <string.h>     // memset, memcpy
#include <stdlib.h>     // qsort
#include <stdio.h>      // printf etc
//#include <stdarg.h>     // ...

#define PUSH_GPU_SECTION(lbl)                                                                   \
{                                                                                               \
    if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
}
#define POP_GPU_SECTION()                       \
{                                           \
    if (glPopDebugGroup) glPopDebugGroup(); \
}

#define MD_UBO_SIZE (1 << 10)
#define MD_SHADER_BUF_SIZE KILOBYTES(14)

enum {
    GL_VERSION_UNKNOWN = 0,
    GL_VERSION_330,
    GL_VERSION_430
};

enum {
    GL_PROGRAM_COMPUTE_AABB,
    GL_PROGRAM_CULL_AABB,
    GL_PROGRAM_EXTRACT_CONTROL_POINTS,
    GL_PROGRAM_COMPUTE_SPLINE,
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
    GL_BUFFER_MOL_BOND_ATOM_INDICES,               // u32[2]
    GL_BUFFER_MOL_BACKBONE_DATA,                   // u32[5] residue index, segment index, CA index, C index and O Index
    GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE,    // u8[4]  (0: Unknown, 1: Coil, 2: Helix, 3: Sheet)
    GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_DATA,     // Extracted control points before spline subdivision
    GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_INDEX,    // u32, LINE_STRIP_ADJACENCY Indices for spline subdivision, seperated by primitive restart index 0xFFFFFFFF
    GL_BUFFER_MOL_BACKBONE_SPLINE_DATA,            // Subdivided control points of spline
    GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX,           // u32, LINE_STRIP indices for rendering, seperated by primitive restart index 0xFFFFFFFF
    GL_BUFFER_MOL_COUNT
};

enum {
    MOL_FLAG_POSITION_CHANGED = (1 << 0)
};

enum {
    PERMUTATION_BIT_ORTHO = 1,
    PERMUTATION_BIT_COLOR = 2,
    PERMUTATION_BIT_NORMAL = 4,
    PERMUTATION_BIT_VELOCITY = 8,
    PERMUTATION_BIT_INDEX = 16,
};

#define MAX_SHADER_PERMUTATIONS 32

internal inline bool is_ortho_proj_matrix(const mat4_t M) { return M.elem[2][3] == 0.0f; }

internal inline void extract_jitter_uv(float jitter[2], const mat4_t  M) {
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
    mat4_t world_to_view;
    mat4_t world_to_view_normal;
    mat4_t world_to_clip;
    mat4_t view_to_clip;
    mat4_t view_to_world;
    mat4_t clip_to_view;
    mat4_t prev_world_to_clip;
    mat4_t curr_view_to_prev_clip;
} gl_view_transform;

// Shared ubo data for all shaders
typedef struct gl_ubo_base {
    gl_view_transform view_transform;
    vec4_t jitter_uv;
    uint32_t atom_mask;
    uint32_t _pad[3];
} gl_ubo_base;

static const int ubo_base_size = sizeof(gl_ubo_base);

typedef uint32_t gl_version;

static const int control_point_size = sizeof(gl_control_point);

typedef struct internal_rep internal_rep;
typedef struct internal_mol internal_mol;
typedef struct internal_ctx internal_ctx;

struct internal_rep {
    const internal_mol* mol;
    md_draw_representation_args  args;
    md_draw_representation_type  type;
    uint32_t            mol_version;
    
    //md_range visible_range;      // Minimal spanning visible range of atoms to draw
    gl_buffer color;
};

struct internal_mol {
    gl_buffer buffer[GL_BUFFER_MOL_COUNT];
    uint32_t version;
    uint32_t flags;

    uint32_t atom_count;
    uint32_t residue_count;
    uint32_t bond_count;
    uint32_t backbone_count;    
    uint32_t backbone_control_point_index_count;
    uint32_t backbone_spline_index_count;
};

struct internal_ctx {
    struct {
        gl_program space_fill[MAX_SHADER_PERMUTATIONS];
        gl_program licorice[MAX_SHADER_PERMUTATIONS];
        gl_program ribbons[MAX_SHADER_PERMUTATIONS];
        gl_program cartoon[MAX_SHADER_PERMUTATIONS];
    } permuted_program;

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

STATIC_ASSERT(sizeof(internal_ctx) <= sizeof(md_draw_context), "Internal draw ctx does not fit into containing structure");
STATIC_ASSERT(sizeof(internal_mol) <= sizeof(md_draw_molecule), "Internal draw mol does not fit into containing structure");
STATIC_ASSERT(sizeof(internal_rep) <= sizeof(md_draw_representation), "Internal draw rep does not fit into containing structure");

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
    glBindBuffer(GL_ARRAY_BUFFER, buf.id);
    glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
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
    glBindBuffer(GL_ARRAY_BUFFER, buf.id);
    glBufferSubData(GL_ARRAY_BUFFER, byte_offset, byte_size, data);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
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

internal md_draw_error compile_shader_from_source(GLuint shader, const char* shader_src, const char* defines) {
    if (defines) {
        const char* src = shader_src;
        char version_str[32] = {0};
        
        if (strncmp(shader_src, "#version", 8) == 0) {
            const char* endline;
            if ((endline = strchr(shader_src, '\n')) != NULL) {
                ++endline;
                const size_t length = endline - shader_src;
                strncpy(version_str, shader_src, MIN(ARRAY_SIZE(version_str), length));
            }
            src = endline;
        }
        const char* sources[5] = {
            version_str,
            "\n",
            defines,
            "\n",
            src
        };
        glShaderSource(shader, ARRAY_SIZE(sources), sources, 0);
    } else {
        const char* c_src = shader_src;
        glShaderSource(shader, 1, &c_src, 0);
    }
    
    GLint success;
    glCompileShader(shader);
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char err_buf[256];
        glGetShaderInfoLog(shader, ARRAY_SIZE(err_buf), NULL, err_buf);
        md_printf(MD_LOG_TYPE_ERROR, "Shader compile error:\n%s\n", err_buf);
        return MD_DRAW_GL_SHADER_COMPILE_ERROR;
    }
    
    return MD_DRAW_SUCCESS;
}

internal md_draw_error compile_shader_from_file(GLuint shader, const char* filename, const char* defines) {
    md_file_o* file = md_file_open((str_t){.ptr = filename, .len = strlen(filename)}, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        char buffer[MD_SHADER_BUF_SIZE] = {0};
        uint64_t file_size = md_file_size(file);
        const uint64_t size = MIN(file_size, ARRAY_SIZE(buffer) - 1);
        md_file_read(file, buffer, size);
        md_printf(MD_LOG_TYPE_INFO, "compiling shader '%s'... ", filename);
        md_draw_error err = compile_shader_from_source(shader, buffer, defines);
        if (err == MD_DRAW_SUCCESS) {
            md_print(MD_LOG_TYPE_INFO, "OK\n");
        }
        md_file_close(file);
        return err;
    } else {
        md_printf(MD_LOG_TYPE_ERROR, "Could not open file file '%s'", filename);
        return MD_DRAW_FILE_NOT_FOUND;
    }
}

internal md_draw_error link_program(GLuint program, const GLuint shader[], uint32_t count) {
    ASSERT(program);
     
    for (uint32_t i = 0; i < count; i++) {
        GLint compile_status;
        glGetShaderiv(shader[i], GL_COMPILE_STATUS, &compile_status);
        if (glIsShader(shader[i]) && compile_status) {
            glAttachShader(program, shader[i]);
        } else {
            md_print(MD_LOG_TYPE_ERROR, "Program link error: One or more shaders are invalid\n");
            return MD_DRAW_GL_PROGRAM_LINK_ERROR;
        }
    }
    
    GLint success;
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        char err_buf[256];
        glGetProgramInfoLog(program, ARRAY_SIZE(err_buf), NULL, err_buf);
        md_printf(MD_LOG_TYPE_ERROR, "Program link error:\n%s\n", err_buf);
        return MD_DRAW_GL_PROGRAM_LINK_ERROR;
    }
    
    for (uint32_t i = 0; i < count; i++) {
        glDetachShader(program, shader[i]);
    }
    
    return MD_DRAW_SUCCESS;
}

internal md_draw_error link_program_transform_feedback(GLuint program, const GLuint shader[], uint32_t shader_count, const char* varying[],
                                            uint32_t varying_count, GLenum capture_mode) {
    ASSERT(program);
    ASSERT(capture_mode == GL_INTERLEAVED_ATTRIBS || capture_mode == GL_SEPARATE_ATTRIBS);
    
    for (uint32_t i = 0; i < shader_count; i++) {
        GLint compile_status;
        glGetShaderiv(shader[i], GL_COMPILE_STATUS, &compile_status);
        if (glIsShader(shader[i]) && compile_status) {
            glAttachShader(program, shader[i]);
        } else {
            md_print(MD_LOG_TYPE_ERROR, "One or more shaders are invalid");
            return MD_DRAW_GL_PROGRAM_LINK_ERROR;
        }
    }
    
    GLint link_status;
    glTransformFeedbackVaryings(program, varying_count, varying, capture_mode);
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &link_status);
    if (!link_status) {
        char err_buf[256];
        glGetProgramInfoLog(program, ARRAY_SIZE(err_buf), NULL, err_buf);
        md_printf(MD_LOG_TYPE_ERROR, "Program link error:\n%s\n", err_buf);
        return MD_DRAW_GL_PROGRAM_LINK_ERROR;
    }
    
    for (uint32_t i = 0; i < shader_count; i++) {
        glDetachShader(program, shader[i]);
    }
    
    return MD_DRAW_SUCCESS;
}

md_draw_error validate_context(const internal_ctx* ctx) {
    if (!ctx) {
        fprintf(stderr, "Internal context has not been initialized"); // We don't have ctx and therefore no log_i
        return MD_DRAW_CONTEXT_INVALID;
    }
    else if (ctx->version != GL_VERSION_330 && ctx->version != GL_VERSION_430) {
        md_print(MD_LOG_TYPE_ERROR, "Internal Context has not been initialized");
        return MD_DRAW_CONTEXT_INVALID;
    }
    return MD_DRAW_SUCCESS;
}

md_draw_error md_draw_molecule_set_atom_position(md_draw_molecule* ext_mol, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t byte_stride) {
    if (ext_mol && x && y && z) {
        internal_mol* mol = (internal_mol*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id) {
            return MD_DRAW_GL_INVALID_BUFFER;
        }
        if (offset + count > mol->atom_count) {
            return MD_DRAW_GL_ATTEMPTING_TO_WRITE_OUT_OF_BUFFER_RANGE;
        }
        byte_stride = MAX(sizeof(float), byte_stride);
        glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
        float* data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (data) {
            mol->flags |= MOL_FLAG_POSITION_CHANGED;
            for (uint32_t i = offset; i < count; ++i) {
                data[i * 3 + 0] = *(const float*)((const uint8_t*)x + byte_stride * i);
                data[i * 3 + 1] = *(const float*)((const uint8_t*)y + byte_stride * i);
                data[i * 3 + 2] = *(const float*)((const uint8_t*)z + byte_stride * i);
            }
            glUnmapBuffer(GL_ARRAY_BUFFER);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return MD_DRAW_SUCCESS;
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        return MD_DRAW_GL_COULD_NOT_MAP_BUFFER;
    }
    return MD_DRAW_ARGUMENT_IS_NULL;
}

md_draw_error md_draw_molecule_set_atom_radius(md_draw_molecule* ext_mol, uint32_t offset, uint32_t count, const float* radius, uint32_t byte_stride) {
    if (ext_mol && radius) {
        internal_mol* mol = (internal_mol*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS].id) {
            return MD_DRAW_GL_INVALID_BUFFER;
        }
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
                return MD_DRAW_SUCCESS;
            }
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return MD_DRAW_GL_COULD_NOT_MAP_BUFFER;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS], offset * sizeof(float), count * sizeof(float), radius);
            return MD_DRAW_SUCCESS;
        }
    }
    return MD_DRAW_ARGUMENT_IS_NULL;
}

md_draw_error md_draw_molecule_set_atom_flags(md_draw_molecule* ext_mol, uint32_t offset, uint32_t count, const uint8_t* flags, uint32_t byte_stride) {
    if (ext_mol && flags) {
        internal_mol* mol = (internal_mol*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id) {
            return MD_DRAW_GL_INVALID_BUFFER;
        }
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
                return MD_DRAW_SUCCESS;
            }
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return MD_DRAW_GL_COULD_NOT_MAP_BUFFER;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS], offset * sizeof(uint8_t), count * sizeof(uint8_t), flags);
            return MD_DRAW_SUCCESS;
        }
    }
    return MD_DRAW_ARGUMENT_IS_NULL;
}

md_draw_error md_draw_molecule_set_backbone_secondary_structure(md_draw_molecule* ext_mol, uint32_t offset, uint32_t count, const md_secondary_structure* secondary_structure, uint32_t byte_stride) {
    if (ext_mol && secondary_structure) {
        internal_mol* mol = (internal_mol*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE].id) {
            return MD_DRAW_GL_INVALID_BUFFER;
        }
        byte_stride = MAX(sizeof(md_secondary_structure), byte_stride);
        if (byte_stride > sizeof(md_secondary_structure)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE].id);
            md_secondary_structure* data = (md_secondary_structure*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (data) {
                for (uint32_t i = 0; i < mol->atom_count; ++i) {
                    data[i] = secondary_structure[byte_stride * i];
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return MD_DRAW_SUCCESS;
            }
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return MD_DRAW_GL_COULD_NOT_MAP_BUFFER;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE], offset * sizeof(md_secondary_structure), count * sizeof(md_secondary_structure), secondary_structure);
            return MD_DRAW_SUCCESS;
        }
    }
    return MD_DRAW_ARGUMENT_IS_NULL;
}

md_draw_error md_draw_molecule_update_atom_previous_position(md_draw_molecule* ext_mol) {
    if (ext_mol) {
        internal_mol* mol = (internal_mol*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id) {
            return MD_DRAW_GL_INVALID_BUFFER;
        }
        if (!mol->buffer[GL_BUFFER_MOL_ATOM_PREV_POSITION].id) {
            return MD_DRAW_GL_INVALID_BUFFER;
        }
        
        uint32_t size = gl_buffer_size(mol->buffer[GL_BUFFER_MOL_ATOM_POSITION]);
        glBindBuffer(GL_COPY_READ_BUFFER,  mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
        glBindBuffer(GL_COPY_WRITE_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_PREV_POSITION].id);
        glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, size);
        glBindBuffer(GL_COPY_READ_BUFFER, 0);
        glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
    }
    return MD_DRAW_ARGUMENT_IS_NULL;
}

/*
internal md_draw_error md_update_visible_atom_color_range(internal_rep* rep) {
    if (!rep) return MD_DRAW_REPRESENTATION_INVALID;
    if (!rep->mol) return MD_DRAW_MOLECULE_INVALID;

    glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
    void* ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);
    if (ptr) {
        uint32_t* col = (uint32_t*)ptr;
        md_range range = {0xFFFFFFFF, 0};
            
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
        return MD_DRAW_GL_COULD_NOT_MAP_BUFFER;
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    return MD_DRAW_SUCCESS;
}
*/

md_draw_error md_draw_representation_set_type_and_args(md_draw_representation* ext_rep, md_draw_representation_type type, md_draw_representation_args args) {
    internal_rep* rep = (internal_rep*)ext_rep;
    rep->type = type;
    rep->args = args;
    return MD_DRAW_SUCCESS;
}

md_draw_error md_draw_representation_set_color(md_draw_representation* ext_rep, uint32_t offset, uint32_t count, const uint32_t* color_data, uint32_t byte_stride) {
    if (ext_rep && color_data) {
        internal_rep* rep = (internal_rep*)ext_rep;
        if (rep->mol && glIsBuffer(rep->color.id)) {
            if (offset + count > rep->mol->atom_count) {
                return MD_DRAW_GL_ATTEMPTING_TO_WRITE_OUT_OF_BUFFER_RANGE;
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
                    return MD_DRAW_SUCCESS;
                }
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return MD_DRAW_GL_COULD_NOT_MAP_BUFFER;
            }
            else {
                gl_buffer_set_sub_data(rep->color, offset * sizeof(uint32_t), count * sizeof(uint32_t), color_data);
            }
            //md_update_visible_atom_color_range(rep);
        } else {
            return MD_DRAW_GL_INVALID_BUFFER;
        }
    }
    return MD_DRAW_ARGUMENT_IS_NULL;
}

md_draw_error create_permuted_program(gl_program* program_permutations, const char* vert_file, const char* geom_file, const char* frag_file) {
    GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

    for (uint32_t perm = 0; perm < MAX_SHADER_PERMUTATIONS; ++perm) {
        char define[128] = {0};
        int  at = 0;

        at += snprintf(define + at, ARRAY_SIZE(define) - at, "#define ORTHO %i\n",     (int)((perm & PERMUTATION_BIT_ORTHO) != 0));
        at += snprintf(define + at, ARRAY_SIZE(define) - at, "#define ATOM_COL %i\n",  (int)((perm & PERMUTATION_BIT_COLOR) != 0));
        at += snprintf(define + at, ARRAY_SIZE(define) - at, "#define VIEW_NORM %i\n", (int)((perm & PERMUTATION_BIT_NORMAL) != 0));
        at += snprintf(define + at, ARRAY_SIZE(define) - at, "#define ATOM_VEL %i\n",  (int)((perm & PERMUTATION_BIT_VELOCITY) != 0));
        at += snprintf(define + at, ARRAY_SIZE(define) - at, "#define ATOM_IDX %i\n",  (int)((perm & PERMUTATION_BIT_INDEX) != 0));

        md_draw_error err;
        if (vert_file && ((err = compile_shader_from_file(vert_shader, vert_file, define)) != MD_DRAW_SUCCESS)) return err;
        if (geom_file && ((err = compile_shader_from_file(geom_shader, geom_file, define)) != MD_DRAW_SUCCESS)) return err;
        if (frag_file && ((err = compile_shader_from_file(frag_shader, frag_file, define)) != MD_DRAW_SUCCESS)) return err;

        program_permutations[perm].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, geom_shader, frag_shader};
        if ((err = link_program(program_permutations[perm].id, shaders, ARRAY_SIZE(shaders))) != MD_DRAW_SUCCESS) return err;
    }

    glDeleteShader(vert_shader);
    glDeleteShader(geom_shader);
    glDeleteShader(frag_shader);

    return MD_DRAW_SUCCESS;
}

md_draw_error md_draw_context_init(md_draw_context* ext_ctx) {
    if (!ext_ctx) {
        md_print(MD_LOG_TYPE_ERROR, "Supplied draw context is NULL");
        return MD_DRAW_ARGUMENT_IS_NULL;
    }

    internal_ctx* ctx = (internal_ctx*)ext_ctx;
    memset(ctx, 0, sizeof(internal_ctx));

    if (gl3wInit() != GL3W_OK) {
        md_print(MD_LOG_TYPE_ERROR, "Could not load OpenGL extensions");
        return MD_DRAW_GL_VERSION_NOT_SUPPORTED;
    }

    GLint major, minor;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);

    if (major >= 4 && minor >= 3) {
        ctx->version = GL_VERSION_430;
    } else if (major >= 3 && minor >= 3) {
        ctx->version = GL_VERSION_330;
    } else {
        ctx->version = GL_VERSION_UNKNOWN;
        md_printf(MD_LOG_TYPE_ERROR, "OpenGL version %i.%i is not supported", major, minor);
        return MD_DRAW_GL_VERSION_NOT_SUPPORTED;
    }

    glGenVertexArrays(1, &ctx->vao);
    glGenFramebuffers(2, ctx->fbo);
    ctx->ubo = gl_buffer_create(MD_UBO_SIZE, NULL, GL_DYNAMIC_DRAW);

    for (uint32_t i = 0; i < GL_TEXTURE_COUNT; ++i) {
        glGenTextures(1, &ctx->texture[i].id);
    }

    {
        md_draw_error err;
        if ((err = create_permuted_program(ctx->permuted_program.space_fill, MD_SHADER_DIR "/space_fill.vert", MD_SHADER_DIR "/space_fill.geom", MD_SHADER_DIR "/space_fill.frag")) != MD_DRAW_SUCCESS) return err;
        if ((err = create_permuted_program(ctx->permuted_program.licorice,   MD_SHADER_DIR "/licorice.vert",   MD_SHADER_DIR "/licorice.geom",   MD_SHADER_DIR "/licorice.frag")) != MD_DRAW_SUCCESS) return err;
        if ((err = create_permuted_program(ctx->permuted_program.ribbons,    MD_SHADER_DIR "/ribbons.vert",    MD_SHADER_DIR "/ribbons.geom",    MD_SHADER_DIR "/ribbons.frag")) != MD_DRAW_SUCCESS) return err;
        if ((err = create_permuted_program(ctx->permuted_program.cartoon,    MD_SHADER_DIR "/cartoon.vert",    MD_SHADER_DIR "/cartoon.geom",    MD_SHADER_DIR "/cartoon.frag")) != MD_DRAW_SUCCESS) return err;
    }

    if (ctx->version == GL_VERSION_430) {
        GLuint comp_shader = glCreateShader(GL_COMPUTE_SHADER);

        md_draw_error err;
        if ((err = compile_shader_from_file(comp_shader, MD_SHADER_DIR "/compute_aabb.comp", NULL)) != MD_DRAW_SUCCESS) {
            return err;
        }
        ctx->program[GL_PROGRAM_COMPUTE_AABB].id = glCreateProgram();
        const GLuint shaders[] = {comp_shader};
        if ((err = link_program(ctx->program[GL_PROGRAM_COMPUTE_AABB].id, shaders, ARRAY_SIZE(shaders))) != MD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(comp_shader);
    }

    if (ctx->version == GL_VERSION_430) {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
        GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);
        GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

        md_draw_error err;
        if ((err = compile_shader_from_file(vert_shader, MD_SHADER_DIR "/draw_aabb.vert", NULL)) != MD_DRAW_SUCCESS ||
            (err = compile_shader_from_file(geom_shader, MD_SHADER_DIR "/draw_aabb.geom", NULL)) != MD_DRAW_SUCCESS ||
            (err = compile_shader_from_file(frag_shader, MD_SHADER_DIR "/cull_aabb.frag", NULL)) != MD_DRAW_SUCCESS) {
            return err;
        }
        ctx->program[GL_PROGRAM_CULL_AABB].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, geom_shader, frag_shader};
        if ((err = link_program(ctx->program[GL_PROGRAM_CULL_AABB].id, shaders, ARRAY_SIZE(shaders))) != MD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
        glDeleteShader(geom_shader);
        glDeleteShader(frag_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);

        md_draw_error err;
        if ((err = compile_shader_from_file(vert_shader, MD_SHADER_DIR "/extract_control_points.vert", NULL)) != MD_DRAW_SUCCESS) {
            return err;
        }
        ctx->program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader};
        const GLchar* varyings[] = {"out_position", "out_atom_idx", "out_velocity", "out_segment_t", "out_secondary_structure_and_flags", "out_support_and_tangent_vector"};
        if ((err = link_program_transform_feedback(ctx->program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) != MD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
        GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);

        char defines[32];
        snprintf(defines, ARRAY_SIZE(defines), "#define NUM_SUBDIVISIONS %i", MD_SPLINE_MAX_SUBDIVISION_COUNT);

        md_draw_error err;
        if ((err = compile_shader_from_file(vert_shader, MD_SHADER_DIR "/compute_spline.vert", defines)) != MD_DRAW_SUCCESS ||
            (err = compile_shader_from_file(geom_shader, MD_SHADER_DIR "/compute_spline.geom", defines)) != MD_DRAW_SUCCESS) {
            return err;
        }
        ctx->program[GL_PROGRAM_COMPUTE_SPLINE].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, geom_shader};
        const GLchar* varyings[] = {"out_position", "out_atom_idx", "out_velocity", "out_segment_t", "out_secondary_structure_and_flags", "out_support_and_tangent_vector"};
        if ((err = link_program_transform_feedback(ctx->program[GL_PROGRAM_COMPUTE_SPLINE].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) != MD_DRAW_SUCCESS) {
            return err;
        }

        glDeleteShader(vert_shader);
        glDeleteShader(geom_shader);
    }

    return MD_DRAW_SUCCESS;
}

md_draw_error md_draw_context_free(md_draw_context* ext_ctx) {
    if (ext_ctx) {
        internal_ctx* ctx = (internal_ctx*)ext_ctx;
        if (ctx->vao) glDeleteVertexArrays(1, &ctx->vao);
        if (ctx->fbo[0]) glDeleteFramebuffers(2, ctx->fbo);
        if (ctx->ubo.id) glDeleteBuffers(1, &ctx->ubo.id);
        for (uint32_t i = 0; i < GL_TEXTURE_COUNT; ++i) {
            if (ctx->texture[i].id) glDeleteTextures(1, &ctx->texture[i].id);
        }
        for (uint32_t i = 0; i < GL_PROGRAM_COUNT; ++i) {
            if (ctx->program[i].id) glDeleteProgram(ctx->program[i].id);
        }
        return MD_DRAW_SUCCESS;
    }
    return MD_DRAW_ARGUMENT_IS_NULL;
}

/*
void md_draw_molecule_desc_extract(md_draw_molecule_desc* desc, const md_molecule* md_mol) {
    desc->atom.count     = md_mol->atom.count;
    desc->atom.x         = md_mol->atom.x;
    desc->atom.y         = md_mol->atom.y;
    desc->atom.z         = md_mol->atom.z;
    desc->atom.radius    = md_mol->atom.radius;
    
    desc->residue.count               = md_mol->residue.count;
    desc->residue.atom_range          = md_mol->residue.atom_range;
    desc->residue.secondary_structure = md_mol->residue.secondary_structure;
    
    desc->bond.count     = md_mol->bond.count;
    desc->bond.atom_bond = md_mol->bond.bond;
}
*/

md_draw_error md_draw_molecule_init(md_draw_molecule* ext_mol, const md_draw_molecule_desc* desc) {
    if (ext_mol && desc) {
        md_draw_molecule_free(ext_mol);
        internal_mol* mol = (internal_mol*)ext_mol;

        mol->atom_count = desc->atom.count;
        mol->buffer[GL_BUFFER_MOL_ATOM_POSITION]      = gl_buffer_create(mol->atom_count * sizeof(float) * 3, NULL, GL_DYNAMIC_DRAW);
        mol->buffer[GL_BUFFER_MOL_ATOM_PREV_POSITION] = gl_buffer_create(mol->atom_count * sizeof(float) * 3, NULL, GL_DYNAMIC_COPY);
        mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS]        = gl_buffer_create(mol->atom_count * sizeof(float) * 1, NULL, GL_STATIC_DRAW);
        mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS]         = gl_buffer_create(mol->atom_count * sizeof(uint8_t),   NULL, GL_STATIC_DRAW);

        if (desc->atom.x && desc->atom.y && desc->atom.z) {
            md_draw_molecule_set_atom_position(ext_mol, 0, mol->atom_count, desc->atom.x, desc->atom.y, desc->atom.z, 0);
            md_draw_molecule_update_atom_previous_position(ext_mol);
        }
        if (desc->atom.radius)                              md_draw_molecule_set_atom_radius(ext_mol, 0, mol->atom_count, desc->atom.radius, 0);
        if (desc->atom.flags)                               md_draw_molecule_set_atom_flags(ext_mol, 0, mol->atom_count, desc->atom.flags, 0);

        mol->residue_count = desc->residue.count;
        mol->buffer[GL_BUFFER_MOL_RESIDUE_ATOM_RANGE]          = gl_buffer_create(mol->residue_count * sizeof(md_range),   NULL, GL_STATIC_DRAW);
        mol->buffer[GL_BUFFER_MOL_RESIDUE_AABB]                = gl_buffer_create(mol->residue_count * sizeof(float) * 6,    NULL, GL_DYNAMIC_COPY);
        mol->buffer[GL_BUFFER_MOL_RESIDUE_VISIBLE]             = gl_buffer_create(mol->residue_count * sizeof(int),          NULL, GL_DYNAMIC_COPY);

        if (desc->residue.atom_range)           gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_RESIDUE_ATOM_RANGE], 0, mol->residue_count * sizeof(uint32_t) * 2, desc->residue.atom_range);
        //if (desc->residue.backbone_atoms)       gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_RESIDUE_BACKBONE_ATOMS], 0, mol->residue_count * sizeof(uint8_t) * 4, desc->residue.backbone_atoms);

        if (desc->chain.count > 0 && desc->chain.backbone_range && desc->backbone.atoms && desc->residue.atom_range && desc->backbone.secondary_structure) {
            uint32_t backbone_residue_count = 0;
            uint32_t backbone_spline_count = 0;
            for (uint32_t i = 0; i < desc->chain.count; ++i) {
                uint32_t res_count = desc->chain.backbone_range[i].end - desc->chain.backbone_range[i].beg;
                backbone_residue_count += res_count;
                backbone_spline_count += (res_count - 1) * MD_SPLINE_MAX_SUBDIVISION_COUNT;
            }

            const uint32_t backbone_count                     = backbone_residue_count;
            const uint32_t backbone_control_point_data_count  = backbone_residue_count;
            const uint32_t backbone_control_point_index_count = backbone_residue_count + desc->chain.count * (2 + 1); // Duplicate pair first and last in each chain for adjacency + primitive restart between
            const uint32_t backbone_spline_data_count         = backbone_spline_count;
            const uint32_t backbone_spline_index_count        = backbone_spline_count + desc->chain.count * (1); // primitive restart between each chain

            mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA]                = gl_buffer_create(backbone_count                     * sizeof(gl_backbone_data), NULL, GL_STATIC_DRAW);
            mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE] = gl_buffer_create(backbone_count                     * sizeof(md_secondary_structure), NULL, GL_DYNAMIC_DRAW);
            mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_DATA]  = gl_buffer_create(backbone_control_point_data_count  * sizeof(gl_control_point), NULL, GL_DYNAMIC_COPY);
            mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_INDEX] = gl_buffer_create(backbone_control_point_index_count * sizeof(uint32_t),         NULL, GL_STATIC_DRAW);
            mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA]         = gl_buffer_create(backbone_spline_data_count         * sizeof(gl_control_point), NULL, GL_DYNAMIC_COPY);
            mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX]        = gl_buffer_create(backbone_spline_index_count        * sizeof(uint32_t),         NULL, GL_STATIC_DRAW);

            //gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE], 0, desc->backbone.count * sizeof(uint8_t) * 4, desc->backbone.secondary_structure);

            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA].id);
            gl_backbone_data* backbone_data = (gl_backbone_data*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (backbone_data) {
                uint32_t idx = 0;
                for (uint32_t i = 0; i < desc->chain.count; ++i) {
                    for (uint32_t j = desc->chain.backbone_range[i].beg; j < desc->chain.backbone_range[i].end; ++j) {
                        backbone_data[idx].residue_idx = j;
                        backbone_data[idx].segment_idx = j - desc->chain.backbone_range[i].beg;
                        backbone_data[idx].ca_idx = desc->backbone.atoms[j].ca;
                        backbone_data[idx].c_idx  = desc->backbone.atoms[j].c;
                        backbone_data[idx].o_idx  = desc->backbone.atoms[j].o;
                        ++idx;
                    }
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                return MD_DRAW_GL_COULD_NOT_MAP_BUFFER;
            }

            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE].id);
            uint32_t* secondary_structure = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (secondary_structure) {
                uint32_t idx = 0;
                for (uint32_t i = 0; i < desc->chain.count; ++i) {
                    for (uint32_t j = desc->chain.backbone_range[i].beg; j < desc->chain.backbone_range[i].end; ++j) {
                        secondary_structure[idx] = desc->backbone.secondary_structure[j];
                        ++idx;
                    }
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                return MD_DRAW_GL_COULD_NOT_MAP_BUFFER;
            }

            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_INDEX].id);
            uint32_t* control_point_index = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (control_point_index) {
                uint32_t idx = 0;
                uint32_t len = 0;
                for (uint32_t i = 0; i < desc->chain.count; ++i) {
                    control_point_index[len++] = idx;
                    for (uint32_t j = desc->chain.backbone_range[i].beg; j < desc->chain.backbone_range[i].end; ++j) {
                        control_point_index[len++] = idx++;
                    }
                    control_point_index[len++] = idx-1;
                    control_point_index[len++] = 0xFFFFFFFF;
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            }
            else {
                return MD_DRAW_GL_COULD_NOT_MAP_BUFFER;
            }

            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX].id);
            uint32_t* spline_index = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (spline_index) {
                uint32_t idx = 0;
                uint32_t len = 0;
                for (uint32_t i = 0; i < desc->chain.count; ++i) {
                    uint32_t res_count = desc->chain.backbone_range[i].end - desc->chain.backbone_range[i].beg;
                    if (res_count > 0) {
                        for (uint32_t j = 0; j < (res_count - 1) * MD_SPLINE_MAX_SUBDIVISION_COUNT; ++j) {
                            spline_index[len++] = idx++;
                        }
                    }
                    spline_index[len++] = 0xFFFFFFFF;
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                return MD_DRAW_GL_COULD_NOT_MAP_BUFFER;
            }

            glBindBuffer(GL_ARRAY_BUFFER, 0);

            mol->backbone_control_point_index_count = backbone_control_point_index_count;
            mol->backbone_spline_index_count = backbone_spline_index_count;
            mol->backbone_count = backbone_count;
        }

        mol->bond_count = desc->covalent_bond.count;
        mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES] = gl_buffer_create(mol->bond_count * sizeof(uint32_t) * 2, NULL, GL_DYNAMIC_COPY);

        if (desc->covalent_bond.atom_bond) gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES], 0, mol->bond_count * sizeof(uint32_t) * 2, desc->covalent_bond.atom_bond);
       
        return MD_DRAW_SUCCESS;
    }
    return MD_DRAW_ARGUMENT_IS_NULL;
}

md_draw_error md_draw_molecule_free(md_draw_molecule* ext_mol) {
    if (ext_mol) {
        internal_mol* mol = (internal_mol*)ext_mol;
        for (uint32_t i = 0; i < GL_BUFFER_MOL_COUNT; ++i) {
            gl_buffer_conditional_delete(&mol->buffer[i]);
        }
        uint32_t new_version = mol->version + 1;
        memset(mol, 0, sizeof(md_draw_molecule));
        mol->version = new_version;
        return MD_DRAW_SUCCESS;
    }
    return MD_DRAW_ARGUMENT_IS_NULL;
}

md_draw_error md_draw_representation_init(md_draw_representation* ext_rep, const md_draw_molecule* ext_mol) {   
    if (ext_rep && ext_mol) {
        md_draw_representation_free(ext_rep);
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
            //md_update_visible_atom_color_range(rep);
            return MD_DRAW_SUCCESS;
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        return MD_DRAW_GL_COULD_NOT_MAP_BUFFER;
    }
    return MD_DRAW_ARGUMENT_IS_NULL;
}

md_draw_error md_draw_representation_free(md_draw_representation* ext_rep) {
    if (ext_rep) {
        internal_rep* rep = (internal_rep*)ext_rep;
        gl_buffer_conditional_delete(&rep->color);
        return MD_DRAW_SUCCESS;
    }
    return MD_DRAW_ARGUMENT_IS_NULL;
}

md_draw_error gl_build_depth_mipmap(const internal_ctx* ctx, gl_texture dst, gl_texture src);
md_draw_error gl_compute_residue_aabb(const internal_ctx* ctx, const internal_mol* mol);
md_draw_error gl_cull_aabb(const internal_ctx* ctx, const internal_mol* mol);
md_draw_error gl_compute_spline(const internal_ctx* ctx, const internal_mol* mol);

md_draw_error gl_draw_space_fill(const internal_ctx* ctx, const internal_rep* rep, int program_permutation);
md_draw_error gl_draw_licorice  (const internal_ctx* ctx, const internal_rep* rep, int program_permutation);
md_draw_error gl_draw_ribbons   (const internal_ctx* ctx, const internal_rep* rep, int program_permutation);
md_draw_error gl_draw_cartoon   (const internal_ctx* ctx, const internal_rep* rep, int program_permutation);

int compare_draw_rep(const void* elem1, const void* elem2) {
    const internal_rep* r1 = (const internal_rep*)elem1;
    const internal_rep* r2 = (const internal_rep*)elem2;
    
    return (r1->type - r2->type) * 2 + (int)(r1->mol - r2->mol);
}

md_draw_error md_draw(md_draw_context* ext_ctx, const md_draw_desc* desc) {
    internal_ctx* ctx = (internal_ctx*)ext_ctx;
    if (!desc) return MD_DRAW_ARGUMENT_IS_NULL;
    md_draw_error err;
    if ((err = validate_context(ctx)) != MD_DRAW_SUCCESS) return err;

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
    memcpy(&ubo_data.view_transform.world_to_view, desc->view_transform.model_view_matrix, sizeof(mat4_t));
    memcpy(&ubo_data.view_transform.view_to_clip,  desc->view_transform.projection_matrix, sizeof(mat4_t));
    ubo_data.view_transform.world_to_clip = mat4_mul(ubo_data.view_transform.view_to_clip, ubo_data.view_transform.world_to_view);
    ubo_data.view_transform.world_to_view_normal = mat4_transpose(mat4_inverse(ubo_data.view_transform.world_to_view));
    ubo_data.view_transform.view_to_world = mat4_inverse(ubo_data.view_transform.world_to_view);
    ubo_data.view_transform.clip_to_view = mat4_inverse(ubo_data.view_transform.view_to_clip);

    if (desc->view_transform.prev_model_view_matrix && desc->view_transform.projection_matrix) {
        const mat4_t* prev_world_to_view = (const mat4_t*)desc->view_transform.prev_model_view_matrix;
        const mat4_t* prev_view_to_clip  = (const mat4_t*)desc->view_transform.prev_projection_matrix;
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
    const internal_rep* draw_rep[64] = {0};
    uint32_t            draw_rep_count = 0;
        
    const internal_mol* unique_mol_ptr[64] = {0};
    uint32_t            unique_mol_flags[64] = {0};
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
            if (rep->type == MD_DRAW_REP_RIBBONS || rep->type == MD_DRAW_REP_CARTOON) {
                unique_mol_flags[unique_mol_idx] |= MOL_FLAG_SPLINE;
            }
        }
    }
    
    qsort((void*)draw_rep, draw_rep_count, sizeof(md_draw_representation*), compare_draw_rep);
               
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
        

    // possibly resize internal depth buffer to match render_target size
    gl_texture* max_tex = &ctx->texture[GL_TEXTURE_MAX_DEPTH];
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

    int program_permutation = 0;
    if (is_ortho_proj_matrix(ubo_data.view_transform.view_to_clip)) program_permutation |= PERMUTATION_BIT_ORTHO;
        
    // If we have rendertargets -> setup and enable fbo
    if (desc->render_target) {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ctx->fbo[0]);
        GLenum draw_buffers[4] = {0};
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
            program_permutation |= PERMUTATION_BIT_COLOR;
        }
        if (desc->render_target->texture_view_normal) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, desc->render_target->texture_view_normal, 0);
            draw_buffers[count++] = GL_COLOR_ATTACHMENT1;
            program_permutation |= PERMUTATION_BIT_NORMAL;
        }
        if (desc->render_target->texture_view_velocity) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, desc->render_target->texture_view_velocity, 0);
            draw_buffers[count++] = GL_COLOR_ATTACHMENT2;
            program_permutation |= PERMUTATION_BIT_VELOCITY;
        }
        if (desc->render_target->texture_atom_index) {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_2D, desc->render_target->texture_atom_index, 0);
            draw_buffers[count++] = GL_COLOR_ATTACHMENT3;
            program_permutation |= PERMUTATION_BIT_INDEX;
        }
        glDrawBuffers(count, draw_buffers);
        glViewport(0, 0, desc->render_target->width, desc->render_target->height);
            
        GLenum fbo_status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (fbo_status != GL_FRAMEBUFFER_COMPLETE) {
            switch (fbo_status) {
            case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
                md_print(MD_LOG_TYPE_ERROR, "Incomplete fbo attachment");
                return MD_DRAW_GL_FRAMEBUFFER_ERROR;
            case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
                md_print(MD_LOG_TYPE_ERROR, "Missing fbo attachment");
                return MD_DRAW_GL_FRAMEBUFFER_ERROR;
            case GL_FRAMEBUFFER_UNSUPPORTED:
                md_print(MD_LOG_TYPE_ERROR, "Unsupported fbo attachment format");
                return MD_DRAW_GL_FRAMEBUFFER_ERROR;
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
        //glEnable(GL_DEPTH_TEST);
        //glDepthFunc(GL_LESS);
        //glColorMask(1,1,1,1);
        //glDepthMask(1);
        
        for (uint32_t i = 0; i < draw_rep_count; i++) {
            switch (draw_rep[i]->type) {
            case MD_DRAW_REP_SOLVENT_EXCLUDED_SURFACE:
                break;
            case MD_DRAW_REP_DEFAULT:
            case MD_DRAW_REP_SPACE_FILL:
                gl_draw_space_fill(ctx, draw_rep[i], program_permutation);
                break;
            case MD_DRAW_REP_RIBBONS:
                gl_draw_ribbons(ctx, draw_rep[i], program_permutation);
                break;
            case MD_DRAW_REP_CARTOON:
                gl_draw_cartoon(ctx, draw_rep[i], program_permutation);
                break;
            case MD_DRAW_REP_LICORICE:
                gl_draw_licorice(ctx, draw_rep[i], program_permutation);
                break;
            default:
                md_print(MD_LOG_TYPE_ERROR, "Representation had unexpected type");
                return MD_DRAW_UNKNOWN_ERROR;
                break;
            }
        }
    }
    POP_GPU_SECTION()

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, bound_fbo);
    glViewport(bound_viewport[0], bound_viewport[1], bound_viewport[2], bound_viewport[3]);
    glDrawBuffers(bound_draw_buffer_count, (GLenum*)bound_draw_buffer);
        
    POP_GPU_SECTION()
    
    return MD_DRAW_SUCCESS;
}

md_draw_error gl_draw_space_fill(const internal_ctx* ctx, const internal_rep* rep, int program_permutation) {
    const float radius_scale = rep->args.space_fill.radius_scale;
    gl_buffer_set_sub_data(ctx->ubo, sizeof(gl_ubo_base), sizeof(radius_scale), &radius_scale);

    glBindVertexArray(ctx->vao);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_PREV_POSITION].id);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS].id);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_BYTE, 0, 0);

    glEnableVertexAttribArray(4);
    glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
    glVertexAttribPointer(4, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
    
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glUseProgram(ctx->permuted_program.space_fill[program_permutation].id);
    glDrawArrays(GL_POINTS, 0, rep->mol->atom_count);
    glUseProgram(0);
    
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);
    
    glBindVertexArray(0);
    
    return MD_DRAW_SUCCESS;
}

md_draw_error gl_draw_licorice(const internal_ctx* ctx, const internal_rep* rep, int program_permutation) {
    const float radius = rep->args.licorice.radius;
    gl_buffer_set_sub_data(ctx->ubo, sizeof(gl_ubo_base), sizeof(radius), &radius);

    glBindVertexArray(ctx->vao);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_PREV_POSITION].id);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);
    glVertexAttribIPointer(2, 1, GL_UNSIGNED_BYTE, 0, 0);

    glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
    glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES].id);

    glUseProgram(ctx->permuted_program.licorice[program_permutation].id);
    glDrawElements(GL_LINES, rep->mol->bond_count * 2, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);

    glBindVertexArray(0);

    return MD_DRAW_SUCCESS;
}

md_draw_error gl_draw_ribbons(const internal_ctx* ctx, const internal_rep* rep, int program_permutation) {
    ASSERT(rep->mol);
    if (!rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id) {
        md_print(MD_LOG_TYPE_ERROR, "No spline present for molecule, which is required for rendering ribbons representation");
        return MD_DRAW_GL_INVALID_BUFFER;
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

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx->texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, rep->mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);

    GLuint program = ctx->permuted_program.ribbons[program_permutation].id;
    glUseProgram(program);
    glUniform1i(glGetUniformLocation(program, "u_atom_color_buffer"), 0);
    glUniform1i(glGetUniformLocation(program, "u_atom_flags_buffer"), 1);
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

    return MD_DRAW_SUCCESS;
}

md_draw_error gl_draw_cartoon(const internal_ctx* ctx, const internal_rep* rep, int program_permutation) {
    ASSERT(rep->mol);
    if (!rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id) {
        md_print(MD_LOG_TYPE_ERROR, "No spline present for molecule, which is required for rendering ribbons representation");
        return MD_DRAW_GL_INVALID_BUFFER;
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

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx->texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, rep->mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);

    GLuint program = ctx->permuted_program.cartoon[program_permutation].id;
    glUseProgram(program);
    glUniform1i(glGetUniformLocation(program, "u_atom_color_buffer"), 0);
    glUniform1i(glGetUniformLocation(program, "u_atom_flags_buffer"), 1);
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

    return MD_DRAW_SUCCESS;
}

md_draw_error gl_compute_residue_aabb(const internal_ctx* ctx, const internal_mol* mol) {   
    (void)ctx;
    (void)mol;

    return MD_DRAW_SUCCESS;
}

md_draw_error gl_cull_aabb(const internal_ctx* ctx, const internal_mol* mol) {
    (void)ctx;
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
    
    return MD_DRAW_SUCCESS;
}

md_draw_error gl_compute_spline(const internal_ctx* ctx, const internal_mol* mol) {
    if (mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA].id == 0) {
        md_print(MD_LOG_TYPE_ERROR, "Backbone data buffer is zero, which is required to compute the backbone. Is the molecule missing a backbone?");
        return MD_DRAW_GL_INVALID_BUFFER;
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
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE].id);

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
    glDrawArrays(GL_POINTS, 0, mol->backbone_count);
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

    return MD_DRAW_SUCCESS;
}