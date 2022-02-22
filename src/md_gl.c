
#include "md_gl.h"
#include "md_util.h"
#include "ext/gl3w/gl3w.h"
#include <core/md_common.h>
#include <core/md_compiler.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_file.h>
#include <core/md_allocator.h>
#include <core/md_array.inl>

#include <stdbool.h>
#include <string.h>     // memset, memcpy
#include <stdlib.h>     // qsort
#include <stdio.h>      // printf etc

#define MAGIC 0xfacb8172U
#define SPLINE_MAX_SUBDIVISION_COUNT 8

#define PUSH_GPU_SECTION(lbl)                                                                   \
{                                                                                               \
    if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
}
#define POP_GPU_SECTION()                       \
{                                           \
    if (glPopDebugGroup) glPopDebugGroup(); \
}

#define UBO_SIZE (1 << 10)
#define SHADER_BUF_SIZE KILOBYTES(14)

static const char* default_shader_output = 
"layout(location = 0) out vec4 out_color;\n"
"void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint atom_index) {\n"
"    out_color = color;\n"
"}\n";

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
    GL_PROGRAM_COMPUTE_VELOCITY,
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
    GL_BUFFER_MOL_ATOM_POSITION_PREV,              // vec3
    GL_BUFFER_MOL_ATOM_RADIUS,                     // float
    GL_BUFFER_MOL_ATOM_VELOCITY,                   // vec3
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
    PERMUTATION_BIT_ORTHO = 1
};

enum {
    MOL_FLAG_HAS_BACKBONE = 1
};

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

static inline uint32_t compute_mip_count(uint32_t width, uint32_t height) {
    uint32_t mips = 1;
    uint32_t max_dim = MAX(width, height);
    while (max_dim >>= 1) ++mips;
    return mips;
}

typedef struct gl_aabb_t {
    float min[3];
    float max[3];
} gl_aabb_t;

typedef struct gl_buffer_t {
    GLuint id;
} gl_buffer_t;

typedef struct gl_texture_t {
    GLuint id;
} gl_texture_t;

typedef struct gl_program_t {
    GLuint id;
} gl_program_t;

typedef struct gl_backbone_data_t {
    uint32_t residue_idx;               // @NOTE: the residue index which the backbone segment originates from
    uint32_t segment_idx;               // @NOTE: the segment index for the backbone, local to each chain -> 0 to N-1 where N is the number of residues per chain
    uint32_t ca_idx, c_idx, o_idx;
} gl_backbone_data_t;

typedef struct gl_control_point_t {
    float position[3];
    uint32_t atom_idx;
    float velocity[3];
    float segment_t;                    // @NOTE: stores the segment index (integer part) and the fraction within the segment (fractional part)
    uint8_t secondary_structure[3];     // @NOTE: Secondary structure as fractions within the components (0 = coil, 1 = helix, 2 = sheet) so we later can smoothly blend between them when subdividing.
    uint8_t flags;                      // @NOTE: Bitfield for setting flags, bits: [0] beg_chain, [1] end_chain, [2] beg_secondary_structure, [3] end_secondary_structure
    int16_t support_vector[3];
    int16_t tangent_vector[3];
} gl_control_point_t;

typedef struct gl_view_transform_t {
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
typedef struct gl_ubo_base_t {
    gl_view_transform_t view_transform;
    vec4_t jitter_uv;
    uint32_t atom_mask;
    uint32_t _pad[3];
} gl_ubo_base_t;

typedef uint32_t gl_version_t;

typedef struct context_t {
    GLuint vao;
    GLuint fbo;
    gl_buffer_t ubo;
    gl_texture_t texture[GL_TEXTURE_COUNT];
    gl_program_t program[GL_PROGRAM_COUNT];
    gl_version_t version;
} context_t;

static context_t ctx = {0};

typedef struct internal_shaders_t {
    gl_program_t space_fill[MAX_SHADER_PERMUTATIONS];
    gl_program_t licorice[MAX_SHADER_PERMUTATIONS];
    gl_program_t ribbons[MAX_SHADER_PERMUTATIONS];
    gl_program_t cartoon[MAX_SHADER_PERMUTATIONS];
} internal_shaders_t;

typedef struct internal_mol_t {
    gl_buffer_t buffer[GL_BUFFER_MOL_COUNT];
    uint32_t version;
    uint32_t flags;
    uint32_t magic;

    uint32_t atom_count;
    uint32_t residue_count;
    uint32_t bond_count;
    uint32_t backbone_count;
    uint32_t backbone_control_point_index_count;
    uint32_t backbone_spline_index_count;
} internal_mol_t;

typedef struct internal_rep_t {
    internal_mol_t* mol;
    md_gl_representation_args_t  args;
    md_gl_representation_type_t  type;
    uint32_t mol_version;
    uint32_t magic;
    
    //md_range_t visible_range;      // Minimal spanning visible range of atoms to draw
    gl_buffer_t color;
} internal_rep_t;

STATIC_ASSERT(sizeof(internal_shaders_t) <= sizeof(md_gl_shaders_t), "internal draw shaders does not fit into containing structure");
STATIC_ASSERT(sizeof(internal_mol_t) <= sizeof(md_gl_molecule_t), "internal draw mol does not fit into containing structure");
STATIC_ASSERT(sizeof(internal_rep_t) <= sizeof(md_gl_representation_t), "internal draw rep does not fit into containing structure");

static inline gl_buffer_t gl_buffer_create(uint32_t num_bytes, const void* data, GLenum usage_hint) {
    gl_buffer_t buf = {0};
    glGenBuffers(1, &buf.id);
    glBindBuffer(GL_ARRAY_BUFFER, buf.id);
    glBufferData(GL_ARRAY_BUFFER, num_bytes, data, usage_hint);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
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

static inline void gl_buffer_set_sub_data(gl_buffer_t buf, uint32_t byte_offset, uint32_t byte_size, const void* data) {
    glBindBuffer(GL_ARRAY_BUFFER, buf.id);
    glBufferSubData(GL_ARRAY_BUFFER, byte_offset, byte_size, data);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static inline uint32_t gl_texture_width(gl_texture_t tex, uint32_t level) {
    GLint width;
    glBindTexture(GL_TEXTURE_2D, tex.id);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, level, GL_TEXTURE_WIDTH, &width);
    glBindTexture(GL_TEXTURE_2D, 0);
    return (uint32_t)width;
}

static inline uint32_t gl_texture_height(gl_texture_t tex, uint32_t level) {
    GLint height;
    glBindTexture(GL_TEXTURE_2D, tex.id);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, level, GL_TEXTURE_HEIGHT, &height);
    glBindTexture(GL_TEXTURE_2D, 0);
    return (uint32_t)height;
}

bool extract_version_string(char version_str[], uint32_t version_cap, const char** src) {
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

static bool compile_shader_from_source(GLuint shader, char* shader_src, const char* defines, const char* extra_src) {
    const char* sources[8] = {0};
    uint32_t num_sources = 0;

    if (defines || extra_src) {
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

        if (extra_src) {
            const char search_str[] = "#pragma EXTRA_SRC";
            char* pos = strstr(src, search_str);
            if (!pos) {
                md_print(MD_LOG_TYPE_ERROR, "Missing insertion point marked as '#pragma EXTRA_SRC' in shader!");
                return false;
            }

            // Insert termination character to separate string.
            *pos = '\0';
            const char* pre_extra_src = src;
            const char* post_extra_src = pos + sizeof(search_str);

            sources[num_sources++] = pre_extra_src;
            sources[num_sources++] = extra_src;
            sources[num_sources++] = post_extra_src;
        } else {
            sources[num_sources++] = src;
        }

        glShaderSource(shader, num_sources, sources, 0);
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
        return false;
    }
    
    return true;
}

static bool compile_shader_from_file(GLuint shader, const char* filename, const char* defines, const char* extra_src) {
    md_file_o* file = md_file_open((str_t){.ptr = filename, .len = strlen(filename)}, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        char buffer[SHADER_BUF_SIZE] = {0};
        uint64_t file_size = md_file_size(file);
        const uint64_t size = MIN(file_size, ARRAY_SIZE(buffer) - 1);
        md_file_read(file, buffer, size);
        md_printf(MD_LOG_TYPE_INFO, "compiling shader '%s'... ", filename);
        bool result = compile_shader_from_source(shader, buffer, defines, extra_src);
        if (result == true) md_print(MD_LOG_TYPE_INFO, "OK");
        md_file_close(file);
        return result;
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
        char err_buf[256];
        glGetProgramInfoLog(program, ARRAY_SIZE(err_buf), NULL, err_buf);
        md_printf(MD_LOG_TYPE_ERROR, "Program link error:\n%s\n", err_buf);
        return false;
    }
    
    for (uint32_t i = 0; i < count; i++) {
        glDetachShader(program, shader[i]);
    }
    
    return true;
}

static bool link_program_transform_feedback(GLuint program, const GLuint shader[], uint32_t shader_count, const char* varying[],
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
            return false;
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
        return false;
    }
    
    for (uint32_t i = 0; i < shader_count; i++) {
        glDetachShader(program, shader[i]);
    }
    
    return true;
}

static inline bool validate_context() {
    if (ctx.version == 0) {
        md_print(MD_LOG_TYPE_ERROR, "MD GL module has not been initialized");
        return false;
    }
    return true;
}

bool md_gl_molecule_set_atom_position(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t byte_stride) {
    if (ext_mol && x && y && z) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id) {
            md_print(MD_LOG_TYPE_ERROR, "Molecule position buffer missing");
            return false;
        }
        if (!mol->buffer[GL_BUFFER_MOL_ATOM_POSITION_PREV].id) {
            md_print(MD_LOG_TYPE_ERROR, "Molecule previous position buffer missing");
            return false;
        }
        if (offset + count > mol->atom_count) {
            md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds");
            return false;
        }

        if (!mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id) {
            md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds");
            return false;
        }

        // Copy position buffer to previous position buffer
        uint32_t buffer_size = gl_buffer_size(mol->buffer[GL_BUFFER_MOL_ATOM_POSITION]);
        glBindBuffer(GL_COPY_READ_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
        glBindBuffer(GL_COPY_WRITE_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_POSITION_PREV].id);
        glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, buffer_size);
        glBindBuffer(GL_COPY_READ_BUFFER, 0);
        glBindBuffer(GL_COPY_WRITE_BUFFER, 0);

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
            return true;
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        md_print(MD_LOG_TYPE_ERROR, "Failed to map molecule position buffer");
        return false;
    }
    md_print(MD_LOG_TYPE_ERROR, "One or more arguments are missing, must pass x, y and z for position.");
    return false;
}

bool md_gl_molecule_set_atom_velocity(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t byte_stride) {
    if (ext_mol && x && y && z) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id) {
            md_print(MD_LOG_TYPE_ERROR, "Molecule position buffer missing");
            return false;
        }
        if (offset + count > mol->atom_count) {
            md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds");
            return false;
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
            return true;
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        md_print(MD_LOG_TYPE_ERROR, "Failed to map molecule position buffer");
        return false;
    }
    md_print(MD_LOG_TYPE_ERROR, "One or more arguments are missing");
    return false;
}

bool md_gl_molecule_set_atom_radius(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const float* radius, uint32_t byte_stride) {
    if (ext_mol && radius) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS].id) {
            md_print(MD_LOG_TYPE_ERROR, "Molecule radius buffer missing");
            return false;
        }
        if (offset + count > mol->atom_count) {
            md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds");
            return false;
        }
        byte_stride = MAX(sizeof(float), byte_stride);
        if (byte_stride > sizeof(float)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS].id);
            float* radius_data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (radius_data) {
                for (uint32_t i = offset; i < count; ++i) {
                    radius_data[i] = *(const float*)((const uint8_t*)radius + byte_stride * i);
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return true;
            }
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            md_print(MD_LOG_TYPE_ERROR, "Failed to map molecule radius buffer");
            return false;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS], offset * sizeof(float), count * sizeof(float), radius);
            return true;
        }
    }
    md_print(MD_LOG_TYPE_ERROR, "One or more arguments are missing");
    return false;
}

bool md_gl_molecule_set_atom_flags(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const uint8_t* flags, uint32_t byte_stride) {
    if (ext_mol && flags) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id) {
            md_print(MD_LOG_TYPE_ERROR, "Molecule flags buffer missing");
            return false;
        }
        if (offset + count > mol->atom_count) {
            md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds");
            return false;
        }
        byte_stride = MAX(sizeof(uint8_t), byte_stride);
        if (byte_stride > sizeof(uint8_t)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);
            uint8_t* data = (uint8_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (data) {
                for (uint32_t i = offset; i < count; ++i) {
                    data[i] = flags[byte_stride * i];
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return true;
            }
            md_print(MD_LOG_TYPE_ERROR, "Failed to map molecule flags buffer");
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return false;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS], offset * sizeof(uint8_t), count * sizeof(uint8_t), flags);
            return true;
        }
    }
    md_print(MD_LOG_TYPE_ERROR, "One or more arguments are missing");
    return false;
}

bool md_gl_molecule_set_covalent_bonds(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const md_bond_t* bonds, uint32_t byte_stride) {
    if (ext_mol && bonds) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES].id) {
            md_print(MD_LOG_TYPE_ERROR, "Molecule bond buffer buffer missing");
            return false;
        }
        if (offset + count > mol->bond_count) {
            md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds");
            return false;
        }
        byte_stride = MAX(sizeof(md_bond_t), byte_stride);
        if (byte_stride > sizeof(md_bond_t)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES].id);
            md_bond_t* bond_data = (md_bond_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (bond_data) {
                for (uint32_t i = offset; i < count; ++i) {
                    bond_data[i] = *(const md_bond_t*)((const uint8_t*)bond_data + byte_stride * i);
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return true;
            }
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            md_print(MD_LOG_TYPE_ERROR, "Failed to map molecule bond buffer");
            return false;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES], offset * sizeof(md_bond_t), count * sizeof(md_bond_t), bonds);
            return true;
        }
    }
    md_print(MD_LOG_TYPE_ERROR, "One or more arguments are missing");
    return false;
}

bool md_gl_molecule_set_backbone_secondary_structure(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const md_secondary_structure_t* secondary_structure, uint32_t byte_stride) {
    if (ext_mol && secondary_structure) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE].id) {
            md_print(MD_LOG_TYPE_ERROR, "Molecule secondary structure buffer missing");
            return false;
        }
        if (offset + count > mol->backbone_count) {
            md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds");
            return false;
        }
        byte_stride = MAX(sizeof(md_secondary_structure_t), byte_stride);
        if (byte_stride > sizeof(md_secondary_structure_t)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE].id);
            md_secondary_structure_t* data = (md_secondary_structure_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (data) {
                for (uint32_t i = offset; i < count; ++i) {
                    data[i] = secondary_structure[byte_stride * i];
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return true;
            }
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            md_print(MD_LOG_TYPE_ERROR, "Failed to map molecule secondary structure buffer");
            return false;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE], offset * sizeof(md_secondary_structure_t), count * sizeof(md_secondary_structure_t), secondary_structure);
            return true;
        }
    }
    md_print(MD_LOG_TYPE_ERROR, "One or more arguments are missing");
    return false;
}

bool md_gl_molecule_compute_velocity(md_gl_molecule_t* ext_mol, const float pbc_ext[3]) {
    if (!validate_context(&ctx)) {
        return false;
    }

    if (!ext_mol) {
        md_print(MD_LOG_TYPE_ERROR, "Molecule was null");
        return false;
    }
    internal_mol_t* mol = (internal_mol_t*)ext_mol;

    if (mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY].id == 0) {
        md_print(MD_LOG_TYPE_ERROR, "Velocity buffer is zero");
        return false;
    }

    if (mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY].id == 0) {
        md_print(MD_LOG_TYPE_ERROR, "Velocity buffer is zero");
        return false;
    }

    glEnable(GL_RASTERIZER_DISCARD);
    glBindVertexArray(ctx.vao);

    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_POSITION_PREV].id);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    GLuint program = ctx.program[GL_PROGRAM_COMPUTE_VELOCITY].id;

    const GLint pbc_ext_loc = glGetUniformLocation(program, "u_pbc_ext");

    glUseProgram(program);
    glUniform3fv(pbc_ext_loc, 1, pbc_ext);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY].id);
    glBeginTransformFeedback(GL_POINTS);
    glDrawArrays(GL_POINTS, 0, mol->atom_count);
    glEndTransformFeedback();
    glUseProgram(0);

    glBindVertexArray(0);
    glDisable(GL_RASTERIZER_DISCARD);

    return true;
}

bool md_gl_molecule_zero_velocity(md_gl_molecule_t* ext_mol) {
    if (ext_mol) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY].id) {
            md_print(MD_LOG_TYPE_ERROR, "Molecule position buffer missing");
            return false;
        }
        glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY].id);
        float* data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (data) {
            memset(data, 0, sizeof(float) * 3 * mol->atom_count);
            glUnmapBuffer(GL_ARRAY_BUFFER);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return true;
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        md_print(MD_LOG_TYPE_ERROR, "Failed to map molecule velocity buffer");
        return false;
    }
    md_print(MD_LOG_TYPE_ERROR, "Molecule was null");
    return false;
}

/*
static bool md_update_visible_atom_color_range(internal_rep_t* rep) {
    if (!rep) return md_gl_REPRESENTATION_INVALID;
    if (!rep->mol) return md_gl_MOLECULE_INVALID;

    glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
    void* ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);
    if (ptr) {
        uint32_t* col = (uint32_t*)ptr;
        md_range_t range = {0xFFFFFFFF, 0};
            
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
        return false;
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    return true;
}
*/

bool md_gl_representation_set_type_and_args(md_gl_representation_t* ext_rep, md_gl_representation_type_t type, md_gl_representation_args_t args) {
    internal_rep_t* rep = (internal_rep_t*)ext_rep;
    rep->type = type;
    rep->args = args;
    return true;
}

bool md_gl_representation_set_color(md_gl_representation_t* ext_rep, uint32_t offset, uint32_t count, const uint32_t* color_data, uint32_t byte_stride) {
    if (ext_rep && color_data) {
        internal_rep_t* rep = (internal_rep_t*)ext_rep;
        if (rep->mol && glIsBuffer(rep->color.id)) {
            if (offset + count > rep->mol->atom_count) {
                md_print(MD_LOG_TYPE_ERROR, "Attempting to write out of bounds");
                return false;
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
                    return true;
                }
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return false;
            }
            else {
                gl_buffer_set_sub_data(rep->color, offset * sizeof(uint32_t), count * sizeof(uint32_t), color_data);
                return true;
            }
            //md_update_visible_atom_color_range(rep);
        } else {
            md_print(MD_LOG_TYPE_ERROR, "Representation's molecule is missing");
            return false;
        }
    }
    md_print(MD_LOG_TYPE_ERROR, "One or more arguments are missing");
    return false;
}

bool create_permuted_program(gl_program_t* program_permutations, const char* vert_file, const char* geom_file, const char* frag_file, const char* frag_output_src) {
    GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

    for (uint32_t perm = 0; perm < MAX_SHADER_PERMUTATIONS; ++perm) {

        char defines[128] = {0};
        int  at = 0;

        at += snprintf(defines + at, ARRAY_SIZE(defines) - at, "#define ORTHO %i\n",     (int)((perm & PERMUTATION_BIT_ORTHO) != 0));

        if (vert_file && !compile_shader_from_file(vert_shader, vert_file, defines, NULL)) return false;
        if (geom_file && !compile_shader_from_file(geom_shader, geom_file, defines, NULL)) return false;
        if (frag_file && !compile_shader_from_file(frag_shader, frag_file, defines, frag_output_src)) return false;

        program_permutations[perm].id = glCreateProgram();
        const GLuint shaders[] = {vert_shader, geom_shader, frag_shader};
        if (!link_program(program_permutations[perm].id, shaders, ARRAY_SIZE(shaders))) return false;
    }

    glDeleteShader(vert_shader);
    glDeleteShader(geom_shader);
    glDeleteShader(frag_shader);

    return true;
}

bool md_gl_initialize() {
    if (gl3wInit() != GL3W_OK) {
        md_print(MD_LOG_TYPE_ERROR, "Could not load OpenGL extensions");
        return false;
    }

    GLint major, minor;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);

    if (major < 3 && minor < 3) {
        md_printf(MD_LOG_TYPE_ERROR, "OpenGL version %i.%i is not supported", major, minor);
        return false;
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

    if (ctx.version >= 430) {
        GLuint comp_shader = glCreateShader(GL_COMPUTE_SHADER);

        bool err;
        if ((err = compile_shader_from_file(comp_shader, MD_SHADER_DIR "/compute_aabb.comp", NULL, NULL)) != true) {
            return err;
        }
        ctx.program[GL_PROGRAM_COMPUTE_AABB].id = glCreateProgram();
        const GLuint shaders[] = { comp_shader };
        if ((err = link_program(ctx.program[GL_PROGRAM_COMPUTE_AABB].id, shaders, ARRAY_SIZE(shaders))) != true) {
            return err;
        }

        glDeleteShader(comp_shader);
    }

    if (ctx.version >= 430) {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
        GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);
        GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

        bool err;
        if ((err = compile_shader_from_file(vert_shader, MD_SHADER_DIR "/draw_aabb.vert", NULL, NULL)) != true ||
            (err = compile_shader_from_file(geom_shader, MD_SHADER_DIR "/draw_aabb.geom", NULL, NULL)) != true ||
            (err = compile_shader_from_file(frag_shader, MD_SHADER_DIR "/cull_aabb.frag", NULL, NULL)) != true) {
            return err;
        }
        ctx.program[GL_PROGRAM_CULL_AABB].id = glCreateProgram();
        const GLuint shaders[] = { vert_shader, geom_shader, frag_shader };
        if ((err = link_program(ctx.program[GL_PROGRAM_CULL_AABB].id, shaders, ARRAY_SIZE(shaders))) != true) {
            return err;
        }

        glDeleteShader(vert_shader);
        glDeleteShader(geom_shader);
        glDeleteShader(frag_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);

        if (!compile_shader_from_file(vert_shader, MD_SHADER_DIR "/extract_control_points.vert", NULL, NULL)) {
            return false;
        }
        ctx.program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id = glCreateProgram();
        const GLuint shaders[] = { vert_shader };
        const GLchar* varyings[] = { "out_position", "out_atom_idx", "out_velocity", "out_segment_t", "out_secondary_structure_and_flags", "out_support_and_tangent_vector" };
        if (!link_program_transform_feedback(ctx.program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) {
            return false;
        }

        glDeleteShader(vert_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);

        if (!compile_shader_from_file(vert_shader, MD_SHADER_DIR "/compute_velocity.vert", NULL, NULL)) {
            return false;
        }
        ctx.program[GL_PROGRAM_COMPUTE_VELOCITY].id = glCreateProgram();
        const GLuint shaders[] = { vert_shader };
        const GLchar* varyings[] = { "out_velocity" };
        if (!link_program_transform_feedback(ctx.program[GL_PROGRAM_COMPUTE_VELOCITY].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) {
            return false;
        }

        glDeleteShader(vert_shader);
    }

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
        GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);

        char defines[64];
        snprintf(defines, sizeof(defines), "#define NUM_SUBDIVISIONS %i", SPLINE_MAX_SUBDIVISION_COUNT);

        bool err;
        if ((err = compile_shader_from_file(vert_shader, MD_SHADER_DIR "/compute_spline.vert", defines, NULL)) != true ||
            (err = compile_shader_from_file(geom_shader, MD_SHADER_DIR "/compute_spline.geom", defines, NULL)) != true) {
            return err;
        }
        ctx.program[GL_PROGRAM_COMPUTE_SPLINE].id = glCreateProgram();
        const GLuint shaders[] = { vert_shader, geom_shader };
        const GLchar* varyings[] = { "out_position", "out_atom_idx", "out_velocity", "out_segment_t", "out_secondary_structure_and_flags", "out_support_and_tangent_vector" };
        if ((err = link_program_transform_feedback(ctx.program[GL_PROGRAM_COMPUTE_SPLINE].id, shaders, ARRAY_SIZE(shaders), varyings, ARRAY_SIZE(varyings), GL_INTERLEAVED_ATTRIBS)) != true) {
            return err;
        }

        glDeleteShader(vert_shader);
        glDeleteShader(geom_shader);
    }
    return true;
}

void md_gl_shutdown() {
    if (ctx.vao) glDeleteVertexArrays(1, &ctx.vao);
    if (ctx.fbo) glDeleteFramebuffers(1, &ctx.fbo);
    if (ctx.ubo.id) glDeleteBuffers(1, &ctx.ubo.id);
    for (uint32_t i = 0; i < GL_TEXTURE_COUNT; ++i) {
        if (ctx.texture[i].id) glDeleteTextures(1, &ctx.texture[i].id);
    }
    for (uint32_t i = 0; i < GL_PROGRAM_COUNT; ++i) {
        if (ctx.program[i].id) glDeleteProgram(ctx.program[i].id);
    }
}

bool md_gl_shaders_init(md_gl_shaders_t* ext_shaders, const char* custom_shader_snippet) {
    if (!ext_shaders) {
        md_print(MD_LOG_TYPE_ERROR, "Supplied shaders object is NULL");
        return false;
    }

    internal_shaders_t* shaders = (internal_shaders_t*)ext_shaders;
    memset(shaders, 0, sizeof(internal_shaders_t));

    {
        if (!custom_shader_snippet) {
            custom_shader_snippet = default_shader_output;
        }
        bool err;
        if ((err = create_permuted_program(shaders->space_fill, MD_SHADER_DIR "/space_fill.vert", MD_SHADER_DIR "/space_fill.geom", MD_SHADER_DIR "/space_fill.frag",  custom_shader_snippet)) != true) return err;
        if ((err = create_permuted_program(shaders->licorice,   MD_SHADER_DIR "/licorice.vert",   MD_SHADER_DIR "/licorice.geom",   MD_SHADER_DIR "/licorice.frag",    custom_shader_snippet)) != true) return err;
        if ((err = create_permuted_program(shaders->ribbons,    MD_SHADER_DIR "/ribbons.vert",    MD_SHADER_DIR "/ribbons.geom",    MD_SHADER_DIR "/ribbons.frag",     custom_shader_snippet)) != true) return err;
        if ((err = create_permuted_program(shaders->cartoon,    MD_SHADER_DIR "/cartoon.vert",    MD_SHADER_DIR "/cartoon.geom",    MD_SHADER_DIR "/cartoon.frag",     custom_shader_snippet)) != true) return err;
    }

    return true;
}

bool md_gl_shaders_free(md_gl_shaders_t* ext_shaders) {
    if (ext_shaders) {
        internal_shaders_t* shaders = (internal_shaders_t*)ext_shaders;
        for (int i = 0; i < MAX_SHADER_PERMUTATIONS; ++i) {
            if (glIsProgram(shaders->space_fill[i].id)) glDeleteProgram(shaders->space_fill[i].id);
            if (glIsProgram(shaders->licorice[i].id))   glDeleteProgram(shaders->licorice[i].id);
            if (glIsProgram(shaders->ribbons[i].id))    glDeleteProgram(shaders->ribbons[i].id);
            if (glIsProgram(shaders->cartoon[i].id))    glDeleteProgram(shaders->cartoon[i].id);
        }
        return true;
    }
    return false;
}

bool md_gl_molecule_init(md_gl_molecule_t* ext_mol, const md_molecule_t* mol) {
    if (ext_mol && mol) {
        if (mol->atom.count == 0) {
            md_print(MD_LOG_TYPE_ERROR, "The supplied molecule has no atoms.");
            return false;
        }
        internal_mol_t* gl_mol = (internal_mol_t*)ext_mol;

        gl_mol->atom_count = (uint32_t)mol->atom.count;
        gl_mol->buffer[GL_BUFFER_MOL_ATOM_POSITION]        = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 3, NULL, GL_DYNAMIC_DRAW);
        gl_mol->buffer[GL_BUFFER_MOL_ATOM_POSITION_PREV]   = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 3, NULL, GL_DYNAMIC_COPY);
        gl_mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY]        = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 3, NULL, GL_DYNAMIC_COPY);
        gl_mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS]          = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 1, NULL, GL_STATIC_DRAW);
        gl_mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS]           = gl_buffer_create(gl_mol->atom_count * sizeof(uint8_t),   NULL, GL_STATIC_DRAW);

        if (mol->atom.x && mol->atom.y && mol->atom.z) {
            md_gl_molecule_set_atom_position(ext_mol, 0, gl_mol->atom_count, mol->atom.x, mol->atom.y, mol->atom.z, 0);
        }
        md_gl_molecule_zero_velocity(ext_mol);

        if (mol->atom.radius) md_gl_molecule_set_atom_radius(ext_mol, 0, gl_mol->atom_count, mol->atom.radius, 0);
        if (mol->atom.flags)  md_gl_molecule_set_atom_flags(ext_mol,  0, gl_mol->atom_count, mol->atom.flags, 0);

        gl_mol->residue_count = (uint32_t)mol->residue.count;
        gl_mol->buffer[GL_BUFFER_MOL_RESIDUE_ATOM_RANGE]          = gl_buffer_create(gl_mol->residue_count * sizeof(md_range_t),   NULL, GL_STATIC_DRAW);
        gl_mol->buffer[GL_BUFFER_MOL_RESIDUE_AABB]                = gl_buffer_create(gl_mol->residue_count * sizeof(float) * 6,    NULL, GL_DYNAMIC_COPY);
        gl_mol->buffer[GL_BUFFER_MOL_RESIDUE_VISIBLE]             = gl_buffer_create(gl_mol->residue_count * sizeof(int),          NULL, GL_DYNAMIC_COPY);

        if (mol->residue.atom_range)           gl_buffer_set_sub_data(gl_mol->buffer[GL_BUFFER_MOL_RESIDUE_ATOM_RANGE], 0, gl_mol->residue_count * sizeof(uint32_t) * 2, mol->residue.atom_range);
        //if (desc->residue.backbone_atoms)       gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_RESIDUE_BACKBONE_ATOMS], 0, mol->residue_count * sizeof(uint8_t) * 4, desc->residue.backbone_atoms);

        if (mol->chain.count > 0 && mol->chain.backbone_range && mol->backbone.atoms && mol->residue.atom_range && mol->backbone.secondary_structure) {
            uint32_t backbone_residue_count = 0;
            uint32_t backbone_spline_count = 0;
            for (uint32_t i = 0; i < (uint32_t)mol->chain.count; ++i) {
                uint32_t res_count = mol->chain.backbone_range[i].end - mol->chain.backbone_range[i].beg;
                backbone_residue_count += res_count;
                backbone_spline_count += (res_count - 1) * SPLINE_MAX_SUBDIVISION_COUNT;
            }

            const uint32_t backbone_count                     = backbone_residue_count;
            const uint32_t backbone_control_point_data_count  = backbone_residue_count;
            const uint32_t backbone_control_point_index_count = backbone_residue_count + (uint32_t)mol->chain.count * (2 + 1); // Duplicate pair first and last in each chain for adjacency + primitive restart between
            const uint32_t backbone_spline_data_count         = backbone_spline_count;
            const uint32_t backbone_spline_index_count        = backbone_spline_count + (uint32_t)mol->chain.count * (1); // primitive restart between each chain

            gl_mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA]                = gl_buffer_create(backbone_count                     * sizeof(gl_backbone_data_t),         NULL, GL_STATIC_DRAW);
            gl_mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE] = gl_buffer_create(backbone_count                     * sizeof(md_secondary_structure_t),   NULL, GL_DYNAMIC_DRAW);
            gl_mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_DATA]  = gl_buffer_create(backbone_control_point_data_count  * sizeof(gl_control_point_t),         NULL, GL_DYNAMIC_COPY);
            gl_mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_INDEX] = gl_buffer_create(backbone_control_point_index_count * sizeof(uint32_t),                   NULL, GL_STATIC_DRAW);
            gl_mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA]         = gl_buffer_create(backbone_spline_data_count         * sizeof(gl_control_point_t),         NULL, GL_DYNAMIC_COPY);
            gl_mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX]        = gl_buffer_create(backbone_spline_index_count        * sizeof(uint32_t),                   NULL, GL_STATIC_DRAW);

            //gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE], 0, desc->backbone.count * sizeof(uint8_t) * 4, desc->backbone.secondary_structure);

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA].id);
            gl_backbone_data_t* backbone_data = (gl_backbone_data_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (backbone_data) {
                uint32_t idx = 0;
                for (uint32_t i = 0; i < (uint32_t)mol->chain.count; ++i) {
                    for (uint32_t j = (uint32_t)mol->chain.backbone_range[i].beg; j < (uint32_t)mol->chain.backbone_range[i].end; ++j) {
                        backbone_data[idx].residue_idx = j;
                        backbone_data[idx].segment_idx = j - mol->chain.backbone_range[i].beg;
                        backbone_data[idx].ca_idx = mol->backbone.atoms[j].ca;
                        backbone_data[idx].c_idx  = mol->backbone.atoms[j].c;
                        backbone_data[idx].o_idx  = mol->backbone.atoms[j].o;
                        ++idx;
                    }
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                return false;
            }

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE].id);
            uint32_t* secondary_structure = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (secondary_structure) {
                uint32_t idx = 0;
                for (uint32_t i = 0; i < (uint32_t)mol->chain.count; ++i) {
                    for (uint32_t j = (uint32_t)mol->chain.backbone_range[i].beg; j < (uint32_t)mol->chain.backbone_range[i].end; ++j) {
                        secondary_structure[idx] = mol->backbone.secondary_structure[j];
                        ++idx;
                    }
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                return false;
            }

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_INDEX].id);
            uint32_t* control_point_index = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (control_point_index) {
                uint32_t idx = 0;
                uint32_t len = 0;
                for (uint32_t i = 0; i < (uint32_t)mol->chain.count; ++i) {
                    control_point_index[len++] = idx;
                    for (uint32_t j = (uint32_t)mol->chain.backbone_range[i].beg; j < (uint32_t)mol->chain.backbone_range[i].end; ++j) {
                        control_point_index[len++] = idx++;
                    }
                    control_point_index[len++] = idx-1;
                    control_point_index[len++] = 0xFFFFFFFF;
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            }
            else {
                return false;
            }

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX].id);
            uint32_t* spline_index = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (spline_index) {
                uint32_t idx = 0;
                uint32_t len = 0;
                for (uint32_t i = 0; i < (uint32_t)mol->chain.count; ++i) {
                    uint32_t res_count = mol->chain.backbone_range[i].end - mol->chain.backbone_range[i].beg;
                    if (res_count > 0) {
                        for (uint32_t j = 0; j < (res_count - 1) * SPLINE_MAX_SUBDIVISION_COUNT; ++j) {
                            spline_index[len++] = idx++;
                        }
                    }
                    spline_index[len++] = 0xFFFFFFFF;
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                return false;
            }

            glBindBuffer(GL_ARRAY_BUFFER, 0);

            gl_mol->backbone_control_point_index_count = backbone_control_point_index_count;
            gl_mol->backbone_spline_index_count = backbone_spline_index_count;
            gl_mol->backbone_count = backbone_count;
            gl_mol->flags |= MOL_FLAG_HAS_BACKBONE;
        }

        gl_mol->bond_count = (uint32_t)mol->covalent_bond.count;
        gl_mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES] = gl_buffer_create(gl_mol->bond_count * sizeof(uint32_t) * 2, NULL, GL_DYNAMIC_COPY);

        if (mol->covalent_bond.bond) gl_buffer_set_sub_data(gl_mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES], 0, gl_mol->bond_count * sizeof(uint32_t) * 2, mol->covalent_bond.bond);
       
        gl_mol->magic = MAGIC;
        return true;
    }
    return false;
}

bool md_gl_molecule_free(md_gl_molecule_t* ext_mol) {
    if (ext_mol) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        ASSERT(mol->magic == MAGIC);
        for (uint32_t i = 0; i < GL_BUFFER_MOL_COUNT; ++i) {
            gl_buffer_conditional_delete(&mol->buffer[i]);
        }
        uint32_t new_version = mol->version + 1;
        memset(mol, 0, sizeof(md_gl_molecule_t));
        mol->version = new_version;
        return true;
    }
    return false;
}

bool md_gl_representation_init(md_gl_representation_t* ext_rep, const md_gl_molecule_t* ext_mol) {   
    if (ext_rep && ext_mol) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (mol->atom_count == 0) {
            md_print(MD_LOG_TYPE_ERROR, "Supplied molecule has no atoms.");
            return false;
        }
        if (mol->magic != MAGIC) {
            md_print(MD_LOG_TYPE_ERROR, "Supplied molecule is corrupt.");
            return false;
        }

        internal_rep_t* rep = (internal_rep_t*)ext_rep;
        memset(rep, 0, sizeof(internal_rep_t));
        rep->mol = (internal_mol_t*)ext_mol;
        rep->mol_version = rep->mol->version;
        rep->magic = MAGIC;
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
            return true;
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    return false;
}

bool md_gl_representation_free(md_gl_representation_t* ext_rep) {
    if (ext_rep) {
        internal_rep_t* rep = (internal_rep_t*)ext_rep;
        ASSERT(rep->magic == MAGIC);
        gl_buffer_conditional_delete(&rep->color);
        return true;
    }
    return false;
}

static bool build_depth_mipmap(gl_texture_t dst, gl_texture_t src);
static bool compute_residue_aabb(const internal_mol_t* mol);
static bool cull_aabb(const internal_mol_t* mol);
static bool compute_spline(const internal_mol_t* mol);

static bool draw_space_fill(const internal_shaders_t* shaders, const internal_rep_t* rep, int program_permutation);
static bool draw_licorice  (const internal_shaders_t* shaders, const internal_rep_t* rep, int program_permutation);
static bool draw_ribbons   (const internal_shaders_t* shaders, const internal_rep_t* rep, int program_permutation);
static bool draw_cartoon   (const internal_shaders_t* shaders, const internal_rep_t* rep, int program_permutation);

static inline void init_ubo_base_data(gl_ubo_base_t* ubo_data, const md_gl_draw_args_t* args, const mat4_t* model_matrix) {
    ASSERT(ubo_data);
    ASSERT(args);

    memcpy(&ubo_data->view_transform.world_to_view, args->view_transform.view_matrix, sizeof(mat4_t));
    if (model_matrix) {
        ubo_data->view_transform.world_to_view = mat4_mul(ubo_data->view_transform.world_to_view, *model_matrix);
    }
    memcpy(&ubo_data->view_transform.view_to_clip,  args->view_transform.projection_matrix, sizeof(mat4_t));

    ubo_data->view_transform.world_to_clip = mat4_mul(ubo_data->view_transform.view_to_clip, ubo_data->view_transform.world_to_view);
    ubo_data->view_transform.world_to_view_normal = mat4_transpose(mat4_inverse(ubo_data->view_transform.world_to_view));
    ubo_data->view_transform.view_to_world = mat4_inverse(ubo_data->view_transform.world_to_view);
    ubo_data->view_transform.clip_to_view = mat4_inverse(ubo_data->view_transform.view_to_clip);

    if (args->view_transform.prev_view_matrix && args->view_transform.projection_matrix) {
        const mat4_t* prev_world_to_view = (const mat4_t*)args->view_transform.prev_view_matrix;
        const mat4_t* prev_view_to_clip  = (const mat4_t*)args->view_transform.prev_projection_matrix;
        ubo_data->view_transform.prev_world_to_clip = mat4_mul(*prev_view_to_clip, *prev_world_to_view);
        ubo_data->view_transform.curr_view_to_prev_clip = mat4_mul(ubo_data->view_transform.prev_world_to_clip, ubo_data->view_transform.view_to_world);
        extract_jitter_uv(&ubo_data->jitter_uv.x, ubo_data->view_transform.view_to_clip);
        extract_jitter_uv(&ubo_data->jitter_uv.z, *prev_view_to_clip);
    }
    ubo_data->atom_mask = args->atom_mask;
}

typedef struct draw_mol_t {
    const internal_mol_t* mol;
    md_bitfield_t atom_mask;
    uint32_t flags;
} draw_mol_t;

bool md_gl_draw(const md_gl_draw_args_t* args) {
    if (!args) {
        md_print(MD_LOG_TYPE_ERROR, "draw args object was NULL");
        return false;
    }
    if (!validate_context()) return false;

    if (!args->shaders) {
        md_print(MD_LOG_TYPE_ERROR, "shaders object was NULL");
        return false;
    }
    internal_shaders_t* shaders = (internal_shaders_t*)args->shaders;

    PUSH_GPU_SECTION("MOLD DRAW")
            
    gl_ubo_base_t ubo_data = {0};
    init_ubo_base_data(&ubo_data, args, NULL);

    gl_buffer_set_sub_data(ctx.ubo, 0, sizeof(ubo_data), &ubo_data);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ctx.ubo.id);
        
    const uint32_t MOL_FLAG_COMPUTE_SPLINE = 1;

    md_gl_draw_op_t const** draw_ops = 0;
    draw_mol_t* draw_mols = 0;

    md_allocator_i* alloc = default_temp_allocator;
        
    for (uint32_t i = 0; i < args->draw_operations.count; i++) {
        const md_gl_draw_op_t* draw_op = &args->draw_operations.ops[i];
        const internal_rep_t* rep = (internal_rep_t*)draw_op->rep;
        if (rep && rep->mol) {
            if (rep->magic != MAGIC) {
                md_print(MD_LOG_TYPE_ERROR, "Representation is corrupt");
                continue;
            }

            if (rep->mol_version != rep->mol->version) {
                md_print(MD_LOG_TYPE_ERROR, "Representation was not created for the current molecule");
                continue;
            }

            if (rep->type == MD_GL_REP_RIBBONS || rep->type == MD_GL_REP_CARTOON) {
                if (!(rep->mol->flags & MOL_FLAG_HAS_BACKBONE)) {
                    md_print(MD_LOG_TYPE_ERROR, "Molecule is missing a backbone, therefore the chosen representation cannot be drawn");
                    continue;
                }
            }

            md_array_push(draw_ops, draw_op, default_temp_allocator);
                
            int64_t mol_idx = -1;
            for (int64_t j = 0; j < md_array_size(draw_mols); j++) {
                if (draw_mols[j].mol == rep->mol) {
                    mol_idx = j;
                    break;
                }
            }
            // Not found
            if (mol_idx == -1) {
                mol_idx = md_array_size(draw_mols);
                md_array_push(draw_mols, ((draw_mol_t){ rep->mol, 0 }), default_temp_allocator);
            }

            if (rep->type == MD_GL_REP_RIBBONS || rep->type == MD_GL_REP_CARTOON) {
                draw_mols[mol_idx].flags |= MOL_FLAG_COMPUTE_SPLINE;
            }
        }
    }
    
    //qsort((void*)draw_ent, draw_ent_count, sizeof(draw_entity_t), compare_draw_ent);
               
    // Compute residue AABBs for culling
    if (ctx.version >= 430 && (args->options & MD_GL_OPTION_RESIDUE_OCCLUSION_CULLING)) {
        PUSH_GPU_SECTION("COMPUTE RESIDUE AABB")
        for (int64_t i = 0; i < md_array_size(draw_mols); i++) {
            compute_residue_aabb(draw_mols[i].mol);
        }
        POP_GPU_SECTION()
            
        PUSH_GPU_SECTION("CULL RESIDUE AABB")
        for (uint32_t i = 0; i < md_array_size(draw_mols); i++) {
            cull_aabb(draw_mols[i].mol);
        }
        POP_GPU_SECTION()
    }
            
    PUSH_GPU_SECTION("COMPUTE SPLINE")
    for (uint32_t i = 0; i < md_array_size(draw_mols); i++) {
        if (draw_mols[i].flags & MOL_FLAG_COMPUTE_SPLINE) {
            compute_spline(draw_mols[i].mol);
        }
    }
    POP_GPU_SECTION()
        
    //bool using_internal_depth = false;

    int program_permutation = 0;
    if (is_ortho_proj_matrix(ubo_data.view_transform.view_to_clip)) program_permutation |= PERMUTATION_BIT_ORTHO;
        
    PUSH_GPU_SECTION("DRAW REPRESENTATIONS")
    for (int64_t i = 0; i < md_array_size(draw_ops); i++) {
        const internal_rep_t* rep  = (const internal_rep_t*)draw_ops[i]->rep;
        const mat4_t* model_matrix = (const mat4_t*)draw_ops[i]->model_matrix;

        if (model_matrix) {
            // If we have a model matrix, we need to recompute the entire matrix stack...
            gl_ubo_base_t ubo_tmp = {0};
            init_ubo_base_data(&ubo_tmp, args, model_matrix);
            gl_buffer_set_sub_data(ctx.ubo, 0, sizeof(gl_ubo_base_t), &ubo_tmp);
        }

        switch (rep->type) {
        case MD_GL_REP_SOLVENT_EXCLUDED_SURFACE:
            break;
        case MD_GL_REP_DEFAULT:
        case MD_GL_REP_SPACE_FILL:
            draw_space_fill(shaders, rep, program_permutation);
            break;
        case MD_GL_REP_RIBBONS:
            draw_ribbons(shaders, rep, program_permutation);
            break;
        case MD_GL_REP_CARTOON:
            draw_cartoon(shaders, rep, program_permutation);
            break;
        case MD_GL_REP_LICORICE:
            draw_licorice(shaders, rep, program_permutation);
            break;
        default:
            md_print(MD_LOG_TYPE_ERROR, "Representation had unexpected type");
            return false;
            break;
        }

        if (model_matrix) {
            // Reset to default
            gl_buffer_set_sub_data(ctx.ubo, 0, sizeof(ubo_data), &ubo_data);
        }
    }
    POP_GPU_SECTION()

    POP_GPU_SECTION()
    
    return true;
}

static bool draw_space_fill(const internal_shaders_t* shaders, const internal_rep_t* rep, int program_permutation) {
    ASSERT(shaders);
    ASSERT(shaders->space_fill[program_permutation].id);

    ASSERT(rep);
    ASSERT(rep->mol);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_ATOM_RADIUS].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);
    ASSERT(rep->color.id);

    const float radius_scale = rep->args.space_fill.radius_scale;
    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(radius_scale), &radius_scale);

    glBindVertexArray(ctx.vao);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY].id);
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
    
    glUseProgram(shaders->space_fill[program_permutation].id);
    glDrawArrays(GL_POINTS, 0, rep->mol->atom_count);
    glUseProgram(0);
    
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);
    
    glBindVertexArray(0);
    
    return true;
}

static bool draw_licorice(const internal_shaders_t* shaders, const internal_rep_t* rep, int program_permutation) {
    ASSERT(shaders);
    ASSERT(shaders->licorice[program_permutation].id);

    ASSERT(rep);
    ASSERT(rep->mol);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES].id);
    ASSERT(rep->color.id);

    const float radius = rep->args.licorice.radius;
    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(radius), &radius);

    glBindVertexArray(ctx.vao);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY].id);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);
    glVertexAttribIPointer(2, 1, GL_UNSIGNED_BYTE, 0, 0);

    glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
    glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_BOND_ATOM_INDICES].id);

    glUseProgram(shaders->licorice[program_permutation].id);
    glDrawElements(GL_LINES, rep->mol->bond_count * 2, GL_UNSIGNED_INT, 0);
    glUseProgram(0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);

    glBindVertexArray(0);

    return true;
}

bool draw_ribbons(const internal_shaders_t* shaders, const internal_rep_t* rep, int program_permutation) {
    ASSERT(shaders);
    ASSERT(shaders->ribbons[program_permutation].id);

    ASSERT(rep);
    ASSERT(rep->mol);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);
    ASSERT(rep->color.id);

    float scale[2] = {
        rep->args.ribbons.width_scale,
        rep->args.ribbons.thickness_scale * 0.1f
    };
    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(scale), &scale);

    glBindVertexArray(ctx.vao);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, position));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, atom_idx));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, velocity));

    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_BYTE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, flags));

    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 3, GL_SHORT, GL_TRUE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, support_vector));

    glEnableVertexAttribArray(5);
    glVertexAttribPointer(5, 3, GL_SHORT, GL_TRUE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, tangent_vector));

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX].id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, rep->color.id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, rep->mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);

    GLuint program = shaders->ribbons[program_permutation].id;
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

    return true;
}

static bool draw_cartoon(const internal_shaders_t* shaders, const internal_rep_t* rep, int program_permutation) {
    ASSERT(shaders);
    ASSERT(shaders->cartoon[program_permutation].id);

    ASSERT(rep);
    ASSERT(rep->mol);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);
    ASSERT(rep->color.id);

    float scale[2] = {
        rep->args.ribbons.width_scale,
        rep->args.ribbons.thickness_scale
    };
    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(scale), &scale);

    glBindVertexArray(ctx.vao);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, position));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, atom_idx));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, velocity));

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, segment_t));

    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 3, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, secondary_structure));

    glEnableVertexAttribArray(5);
    glVertexAttribIPointer(5, 1, GL_UNSIGNED_BYTE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, flags));

    glEnableVertexAttribArray(6);
    glVertexAttribPointer(6, 3, GL_SHORT, GL_TRUE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, support_vector));

    glEnableVertexAttribArray(7);
    glVertexAttribPointer(7, 3, GL_SHORT, GL_TRUE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, tangent_vector));

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_INDEX].id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, rep->color.id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, rep->mol->buffer[GL_BUFFER_MOL_ATOM_FLAGS].id);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);

    GLuint program = shaders->cartoon[program_permutation].id;
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

    return true;
}

static bool compute_residue_aabb(const internal_mol_t* mol) {   
    (void)ctx;
    (void)mol;

    return true;
}

static bool cull_aabb(const internal_mol_t* mol) {
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
    
    return true;
}

static bool compute_spline(const internal_mol_t* mol) {
    ASSERT(ctx.program[GL_PROGRAM_COMPUTE_SPLINE].id);
    ASSERT(mol);

    ASSERT(mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA].id);
    ASSERT(mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);
    ASSERT(mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY].id);
    ASSERT(mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE].id);
    ASSERT(mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_DATA].id);
    ASSERT(mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_INDEX].id);
    ASSERT(mol->buffer[GL_BUFFER_MOL_BACKBONE_SPLINE_DATA].id);

    if (mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA].id == 0) {
        md_print(MD_LOG_TYPE_ERROR, "Backbone data buffer is zero, which is required to compute the backbone. Is the molecule missing a backbone?");
        return false;
    }

    glEnable(GL_RASTERIZER_DISCARD);
    glBindVertexArray(ctx.vao);

    // Step 1: Extract control points
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_DATA].id);

    glEnableVertexAttribArray(0);
    glVertexAttribIPointer(0, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, residue_idx));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, segment_idx));

    glEnableVertexAttribArray(2);
    glVertexAttribIPointer(2, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, ca_idx));

    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, c_idx));

    glEnableVertexAttribArray(4);
    glVertexAttribIPointer(4, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, o_idx));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_MOL_ATOM_POSITION].id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_MOL_ATOM_VELOCITY].id);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_2].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, mol->buffer[GL_BUFFER_MOL_BACKBONE_SECONDARY_STRUCTURE].id);

    GLuint program = ctx.program[GL_PROGRAM_EXTRACT_CONTROL_POINTS].id;
    const GLint buf_atom_pos_loc            = glGetUniformLocation(program, "u_buf_atom_pos");
    const GLint buf_atom_vel_loc            = glGetUniformLocation(program, "u_buf_atom_vel");
    const GLint buf_secondary_structure_loc = glGetUniformLocation(program, "u_buf_secondary_structure");

    glUseProgram(program);
    glUniform1i(buf_atom_pos_loc, 0);
    glUniform1i(buf_atom_vel_loc, 1);
    glUniform1i(buf_secondary_structure_loc, 2);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_DATA].id);
    glBeginTransformFeedback(GL_POINTS);
    glDrawArrays(GL_POINTS, 0, mol->backbone_count);
    glEndTransformFeedback();
    glUseProgram(0);

    // Step 2: Compute splines using control points
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_DATA].id);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,        sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, position));

    glEnableVertexAttribArray(1);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_INT,          sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, atom_idx));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE,        sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, velocity));

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE,        sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, segment_t));

    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, secondary_structure));

    glEnableVertexAttribArray(5);
    glVertexAttribIPointer(5, 1, GL_UNSIGNED_BYTE,         sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, flags));

    glEnableVertexAttribArray(6);
    glVertexAttribPointer(6, 3, GL_SHORT, GL_TRUE,         sizeof(gl_control_point_t), (const void*)offsetof(gl_control_point_t, support_vector));
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mol->buffer[GL_BUFFER_MOL_BACKBONE_CONTROL_POINT_INDEX].id);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);

    glUseProgram(ctx.program[GL_PROGRAM_COMPUTE_SPLINE].id);
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

    return true;
}
