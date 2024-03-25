#include <md_gl.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_compiler.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_os.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_str_builder.h>
#include <md_molecule.h>

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
#ifndef MD_GL_SPLINE_SUBDIVISION_COUNT
#define MD_GL_SPLINE_SUBDIVISION_COUNT 8
#endif

#define UBO_SIZE (1 << 10)
#define SHADER_BUF_SIZE KILOBYTES(14)

#define BAKE_STR(s) {s"", sizeof(s)}

static const str_t default_shader_output = BAKE_STR(
"layout(location = 0) out vec4 out_color;\n"
"void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint atom_index) {\n"
"    out_color = color;\n"
"}\n");

#undef BAKE_STR

enum {
    GL_VERSION_UNKNOWN = 0,
    GL_VERSION_330
};

enum {
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
    GL_BUFFER_ATOM_POSITION,                   // vec3
    GL_BUFFER_ATOM_POSITION_PREV,              // vec3
    GL_BUFFER_ATOM_RADIUS,                     // float
    GL_BUFFER_ATOM_VELOCITY,                   // vec3
    GL_BUFFER_ATOM_FLAGS,                      // u8
    //GL_BUFFER_RESIDUE_ATOM_RANGE,              // u32[2]   (beg, end)
    GL_BUFFER_RESIDUE_AABB,                    // vec3[2]  (aabb_min, aabb_max)
    GL_BUFFER_RESIDUE_VISIBLE,                 // int
    GL_BUFFER_BOND_ATOM_INDICES,               // u32[2]
    GL_BUFFER_BACKBONE_DATA,                   // u32[5] residue index, segment index, CA index, C index and O Index
    GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE,    // u8[4]  (0: Unknown, 1: Coil, 2: Helix, 3: Sheet)
    GL_BUFFER_BACKBONE_CONTROL_POINT_DATA,     // Extracted control points before spline subdivision
    GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX,    // u32, LINE_STRIP_ADJACENCY Indices for spline subdivision, seperated by primitive restart index 0xFFFFFFFF
    GL_BUFFER_BACKBONE_SPLINE_DATA,            // Subdivided control points of spline
    GL_BUFFER_BACKBONE_SPLINE_INDEX,           // u32, LINE_STRIP indices for rendering, seperated by primitive restart index 0xFFFFFFFF
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
    uint32_t flags;
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
    uint32_t atom_index_base;
    uint32_t bond_index_base;
    uint32_t _pad[1];
} gl_ubo_base_t;

typedef uint32_t gl_version_t;

typedef struct context_t {
    GLuint vao;
    GLuint fbo;
    gl_buffer_t ubo;
    gl_buffer_t instance_ubo;
    gl_texture_t texture[GL_TEXTURE_COUNT];
    gl_program_t program[GL_PROGRAM_COUNT];
    gl_version_t version;
} context_t;

static context_t ctx = {0};

typedef struct internal_shaders_t {
    gl_program_t spacefill[MAX_SHADER_PERMUTATIONS];
    gl_program_t licorice[MAX_SHADER_PERMUTATIONS];
    gl_program_t ribbons[MAX_SHADER_PERMUTATIONS];
    gl_program_t cartoon[MAX_SHADER_PERMUTATIONS];
} internal_shaders_t;

typedef struct internal_mol_t {
    gl_buffer_t buffer[GL_BUFFER_COUNT];
    uint32_t version;
    uint32_t flags;
    uint32_t magic;

    uint32_t atom_index_base;
    uint32_t bond_index_base;

    uint32_t atom_count;
    uint32_t residue_count;
    uint32_t bond_count;
    uint32_t backbone_count;
    uint32_t backbone_control_point_index_count;
    uint32_t backbone_spline_index_count;
} internal_mol_t;

typedef struct internal_rep_t {
    internal_mol_t* mol;

    uint32_t mol_version;
    uint32_t magic;
    
    //md_range_t visible_range;      // Minimal spanning visible range of atoms to draw
    gl_buffer_t color;
} internal_rep_t;

typedef struct internal_bond_t {
    uint32_t atom_idx[2];
} internal_bond_t;

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

static bool compile_shader_from_source(GLuint shader, str_t src, str_t defines, str_t extra_src) {
    md_strb_t sb = md_strb_create(md_get_temp_allocator());
    
    if (!str_eq_cstr_n(src, "#version ", 9)) {
        MD_LOG_ERROR("Missing version string as first line in shader!");
        return false;
    }
    
    if (!str_empty(defines)) {
        str_t version_line;
        str_extract_line(&version_line, &src);
        md_strb_push_str(&sb, version_line);
        md_strb_push_char(&sb, '\n');
        md_strb_push_str(&sb, defines);
        md_strb_push_char(&sb, '\n');
    }

    if (!str_empty(extra_src)) {
        str_t line;
        while (str_extract_line(&line, &src)) {
            if (!str_eq_cstr_n(line, "#pragma EXTRA_SRC", 17)) {
                md_strb_push_str(&sb, line);
                md_strb_push_char(&sb, '\n');
            } else {
                md_strb_push_str(&sb, extra_src);
                md_strb_push_char(&sb, '\n');
            }
        }
    } else {
        md_strb_push_str(&sb, src);
    }

    const char* csrc = md_strb_to_cstr(sb);
    
    glShaderSource(shader, 1, &csrc, 0);
    glCompileShader(shader);
    
    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

    md_strb_free(&sb);

    if (!success) {
        char err_buf[1024];
        glGetShaderInfoLog(shader, ARRAY_SIZE(err_buf), NULL, err_buf);
        MD_LOG_ERROR("Shader compile error:\n%s\n", err_buf);
        return false;
    }
    
    return true;
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
        char err_buf[256];
        glGetProgramInfoLog(program, ARRAY_SIZE(err_buf), NULL, err_buf);
        MD_LOG_ERROR("Program link error:\n%s\n", err_buf);
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
            MD_LOG_ERROR("One or more shaders are invalid");
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
        MD_LOG_ERROR("Program link error:\n%s\n", err_buf);
        return false;
    }
    
    for (uint32_t i = 0; i < shader_count; i++) {
        glDetachShader(program, shader[i]);
    }
    
    return true;
}

static inline bool validate_context() {
    if (ctx.version == 0) {
        MD_LOG_ERROR("MD GL module has not been initialized");
        return false;
    }
    return true;
}


bool md_gl_molecule_set_index_base(md_gl_molecule_t* ext_mol, uint32_t atom_index_base, uint32_t bond_index_base) {
    if (ext_mol) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        mol->atom_index_base = atom_index_base;
        mol->bond_index_base = bond_index_base;
        return true;
    }
    return false;
}

bool md_gl_molecule_set_atom_position(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t byte_stride) {
    if (ext_mol && x && y && z) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_ATOM_POSITION].id) {
            MD_LOG_ERROR("Molecule position buffer missing");
            return false;
        }
        if (!mol->buffer[GL_BUFFER_ATOM_POSITION_PREV].id) {
            MD_LOG_ERROR("Molecule previous position buffer missing");
            return false;
        }
        if (offset + count > mol->atom_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return false;
        }

        // Copy position buffer to previous position buffer
        uint32_t buffer_size = gl_buffer_size(mol->buffer[GL_BUFFER_ATOM_POSITION]);
        glBindBuffer(GL_COPY_READ_BUFFER, mol->buffer[GL_BUFFER_ATOM_POSITION].id);
        glBindBuffer(GL_COPY_WRITE_BUFFER, mol->buffer[GL_BUFFER_ATOM_POSITION_PREV].id);
        glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, buffer_size);
        glBindBuffer(GL_COPY_READ_BUFFER, 0);
        glBindBuffer(GL_COPY_WRITE_BUFFER, 0);

        byte_stride = MAX(sizeof(float), byte_stride);
        glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_POSITION].id);
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
        MD_LOG_ERROR("Failed to map molecule position buffer");
        return false;
    }
    MD_LOG_ERROR("One or more arguments are missing, must pass x, y and z for position.");
    return false;
}

bool md_gl_molecule_set_atom_position_xyz(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const float* xyz) {
    if (ext_mol && xyz) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_ATOM_POSITION].id) {
            MD_LOG_ERROR("Molecule position buffer missing");
            return false;
        }
        if (offset + count > mol->atom_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return false;
        }
        if (!mol->buffer[GL_BUFFER_ATOM_POSITION].id) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return false;
        }

        uint32_t elem_size = sizeof(float) * 3;
        gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_ATOM_POSITION], offset * elem_size, count * elem_size, xyz);
    }
    MD_LOG_ERROR("One or more arguments are missing, must pass xyz for position.");
    return false;
}

bool md_gl_molecule_set_atom_velocity(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const float* x, const float* y, const float* z, uint32_t byte_stride) {
    if (ext_mol && x && y && z) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_ATOM_POSITION].id) {
            MD_LOG_ERROR("Molecule position buffer missing");
            return false;
        }
        if (offset + count > mol->atom_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return false;
        }
        byte_stride = MAX(sizeof(float), byte_stride);
        glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_POSITION].id);
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
        MD_LOG_ERROR("Failed to map molecule position buffer");
        return false;
    }
    MD_LOG_ERROR("One or more arguments are missing");
    return false;
}

bool md_gl_molecule_set_atom_radius(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const float* radius, uint32_t byte_stride) {
    if (ext_mol && radius) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_ATOM_RADIUS].id) {
            MD_LOG_ERROR("Molecule radius buffer missing");
            return false;
        }
        if (offset + count > mol->atom_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return false;
        }
        byte_stride = MAX(sizeof(float), byte_stride);
        if (byte_stride > sizeof(float)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_RADIUS].id);
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
            MD_LOG_ERROR("Failed to map molecule radius buffer");
            return false;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_ATOM_RADIUS], offset * sizeof(float), count * sizeof(float), radius);
            return true;
        }
    }
    MD_LOG_ERROR("One or more arguments are missing");
    return false;
}

bool md_gl_molecule_set_atom_flags(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const uint8_t* flag_data, uint32_t byte_stride) {
    if (ext_mol && flag_data) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_ATOM_FLAGS].id) {
            MD_LOG_ERROR("Molecule flags buffer missing");
            return false;
        }
        if (offset + count > mol->atom_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return false;
        }
        byte_stride = MAX(sizeof(uint8_t), byte_stride);
        if (byte_stride > sizeof(uint8_t)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
            uint8_t* data = (uint8_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (data) {
                for (uint32_t i = offset; i < count; ++i) {
                    data[i] = *(flag_data + byte_stride * i);
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                return true;
            }
            MD_LOG_ERROR("Failed to map molecule flags buffer");
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return false;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_ATOM_FLAGS], offset * sizeof(uint8_t), count * sizeof(uint8_t), flag_data);
            return true;
        }
    }
    MD_LOG_ERROR("One or more arguments are missing");
    return false;
}

bool md_gl_molecule_set_bonds(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const md_bond_pair_t* bonds, uint32_t byte_stride) {
    if (ext_mol && bonds) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_BOND_ATOM_INDICES].id) {
            MD_LOG_ERROR("Molecule bond buffer buffer missing");
            return false;
        }
        if (offset + count > mol->bond_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return false;
        }
        byte_stride = MAX(sizeof(md_bond_pair_t), byte_stride);
        glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_BOND_ATOM_INDICES].id);
        internal_bond_t* bond_data = (internal_bond_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (bond_data) {
            for (uint32_t i = offset; i < count; ++i) {
                const md_bond_pair_t* bond = (const md_bond_pair_t*)((const uint8_t*)bonds + byte_stride * i);
                bond_data[i].atom_idx[0] = bond->idx[0];
                bond_data[i].atom_idx[1] = bond->idx[1];
            }
            glUnmapBuffer(GL_ARRAY_BUFFER);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return true;
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        MD_LOG_ERROR("Failed to map molecule bond buffer");
        return false;
    }
    MD_LOG_ERROR("One or more arguments are missing");
    return false;
}

bool md_gl_molecule_set_backbone_secondary_structure(md_gl_molecule_t* ext_mol, uint32_t offset, uint32_t count, const md_secondary_structure_t* secondary_structure, uint32_t byte_stride) {
    if (ext_mol && secondary_structure) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE].id) {
            MD_LOG_ERROR("Molecule secondary structure buffer missing");
            return false;
        }
        if (offset + count > mol->backbone_count) {
            MD_LOG_ERROR("Attempting to write out of bounds");
            return false;
        }
        byte_stride = MAX(sizeof(md_secondary_structure_t), byte_stride);
        if (byte_stride > sizeof(md_secondary_structure_t)) {
            glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE].id);
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
            MD_LOG_ERROR("Failed to map molecule secondary structure buffer");
            return false;
        }
        else {
            gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE], offset * sizeof(md_secondary_structure_t), count * sizeof(md_secondary_structure_t), secondary_structure);
            return true;
        }
    }
    MD_LOG_ERROR("One or more arguments are missing");
    return false;
}

bool md_gl_molecule_compute_velocity(md_gl_molecule_t* ext_mol, const float pbc_ext[3]) {
    if (!validate_context()) {
        return false;
    }

    if (!ext_mol) {
        MD_LOG_ERROR("Molecule was null");
        return false;
    }
    internal_mol_t* mol = (internal_mol_t*)ext_mol;

    if (mol->buffer[GL_BUFFER_ATOM_VELOCITY].id == 0) {
        MD_LOG_ERROR("Velocity buffer is zero");
        return false;
    }

    if (mol->buffer[GL_BUFFER_ATOM_VELOCITY].id == 0) {
        MD_LOG_ERROR("Velocity buffer is zero");
        return false;
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

    return true;
}

bool md_gl_molecule_zero_velocity(md_gl_molecule_t* ext_mol) {
    if (ext_mol) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (!mol->buffer[GL_BUFFER_ATOM_VELOCITY].id) {
            MD_LOG_ERROR("Molecule position buffer missing");
            return false;
        }
        glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
        float* data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        if (data) {
            MEMSET(data, 0, sizeof(float) * 3 * mol->atom_count);
            glUnmapBuffer(GL_ARRAY_BUFFER);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            return true;
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        MD_LOG_ERROR("Failed to map molecule velocity buffer");
        return false;
    }
    MD_LOG_ERROR("Molecule was null");
    return false;
}

bool md_gl_representation_set_color(md_gl_representation_t* ext_rep, uint32_t offset, uint32_t count, const uint32_t* color_data, uint32_t byte_stride) {
    if (ext_rep && color_data) {
        internal_rep_t* rep = (internal_rep_t*)ext_rep;
        if (rep->mol && glIsBuffer(rep->color.id)) {
            if (offset + count > rep->mol->atom_count) {
                MD_LOG_ERROR("Attempting to write out of bounds");
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
            MD_LOG_ERROR("Representation's molecule is missing");
            return false;
        }
    }
    MD_LOG_ERROR("One or more arguments are missing");
    return false;
}

bool create_permuted_program(str_t identifier, gl_program_t* program_permutations, str_t vert_src, str_t geom_src, str_t frag_src, str_t frag_output_src) {
    GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint geom_shader = glCreateShader(GL_GEOMETRY_SHADER);
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);

    const str_t perm_str[] = {
        STR_LIT("#define ORTHO 0"),
        STR_LIT("#define ORTHO 1"),
    };
    
    ASSERT(ARRAY_SIZE(perm_str) <= MAX_SHADER_PERMUTATIONS);

    for (uint32_t perm = 0; perm < MAX_SHADER_PERMUTATIONS; ++perm) {
        const str_t defines = perm_str[perm];

        if (!str_empty(vert_src) && !compile_shader_from_source(vert_shader, vert_src, defines, (str_t){0})) {
            MD_LOG_ERROR("Error occured when compiling vertex shader for: '%.*s'", STR_ARG(identifier));
            return false;
        }
            
        if (!str_empty(geom_src) && !compile_shader_from_source(geom_shader, geom_src, defines, (str_t){0})) {
            MD_LOG_ERROR("Error occured when compiling geometry shader for: '%.*s'", STR_ARG(identifier));
            return false;
        }
        if (!str_empty(frag_src) && !compile_shader_from_source(frag_shader, frag_src, defines, frag_output_src)) {
            MD_LOG_ERROR("Error occured when compiling fragment shader for: '%.*s'", STR_ARG(identifier));
            return false;
        }

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
        MD_LOG_ERROR("Could not load OpenGL extensions");
        return false;
    }

    GLint major, minor;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);

    if (major < 3 && minor < 3) {
        MD_LOG_ERROR("OpenGL version %i.%i is not supported", major, minor);
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
    
    const str_t empty_str = {0};

    {
        GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
        if (!compile_shader_from_source(vert_shader, (str_t){(const char*)extract_control_points_vert, extract_control_points_vert_size}, empty_str, empty_str)) {
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

        if (!compile_shader_from_source(vert_shader, (str_t){(const char*)compute_velocity_vert, compute_velocity_vert_size}, empty_str, empty_str)) {
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

        char def_buf[256];
        snprintf(def_buf, ARRAY_SIZE(def_buf), "#define NUM_SUBDIVISIONS %i\n#define MAX_VERTICES %i\n", MD_GL_SPLINE_SUBDIVISION_COUNT, MD_GL_SPLINE_SUBDIVISION_COUNT+1);
        str_t defines = {def_buf, strlen(def_buf)};

        bool err;
        if ((err = compile_shader_from_source(vert_shader, (str_t){(const char*)compute_spline_vert, compute_spline_vert_size}, defines, empty_str)) != true ||
            (err = compile_shader_from_source(geom_shader, (str_t){(const char*)compute_spline_geom, compute_spline_geom_size}, defines, empty_str)) != true) {
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

bool md_gl_shaders_init(md_gl_shaders_t* ext_shaders, const char* custom_shader_snippet_ptr, int64_t custom_shader_snippet_len) {
    if (!ext_shaders) {
        MD_LOG_ERROR("Supplied shaders object is NULL");
        return false;
    }

    internal_shaders_t* shaders = (internal_shaders_t*)ext_shaders;
    MEMSET(shaders, 0, sizeof(internal_shaders_t));

    str_t custom_shader_snippet = {custom_shader_snippet_ptr, custom_shader_snippet_len};

    {
        if (str_empty(custom_shader_snippet)) {
            custom_shader_snippet = default_shader_output;
        }
        if (!create_permuted_program(STR_LIT("SpaceFill"), shaders->spacefill,  (str_t){(const char*)spacefill_vert, spacefill_vert_size}, (str_t){(const char*)spacefill_geom, spacefill_geom_size},   (str_t){(const char*)spacefill_frag, spacefill_frag_size},    custom_shader_snippet)) return false;
        if (!create_permuted_program(STR_LIT("Licorice"),  shaders->licorice,   (str_t){(const char*)licorice_vert, licorice_vert_size},   (str_t){(const char*)licorice_geom, licorice_geom_size},     (str_t){(const char*)licorice_frag, licorice_frag_size},      custom_shader_snippet)) return false;
        if (!create_permuted_program(STR_LIT("Ribbons"),   shaders->ribbons,    (str_t){(const char*)ribbons_vert, ribbons_vert_size},     (str_t){(const char*)ribbons_geom, ribbons_geom_size},       (str_t){(const char*)ribbons_frag, ribbons_frag_size},        custom_shader_snippet)) return false;
        if (!create_permuted_program(STR_LIT("Cartoon"),   shaders->cartoon,    (str_t){(const char*)cartoon_vert, cartoon_vert_size},     (str_t){(const char*)cartoon_geom, cartoon_geom_size},       (str_t){(const char*)cartoon_frag, cartoon_frag_size},        custom_shader_snippet)) return false;
    }

    return true;
}

bool md_gl_shaders_free(md_gl_shaders_t* ext_shaders) {
    if (ext_shaders) {
        internal_shaders_t* shaders = (internal_shaders_t*)ext_shaders;
        for (int i = 0; i < MAX_SHADER_PERMUTATIONS; ++i) {
            if (glIsProgram(shaders->spacefill[i].id))  glDeleteProgram(shaders->spacefill[i].id);
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
            MD_LOG_ERROR("The supplied molecule has no atoms.");
            return false;
        }
        internal_mol_t* gl_mol = (internal_mol_t*)ext_mol;

        gl_mol->bond_index_base = 0x80000000;

        gl_mol->atom_count = (uint32_t)mol->atom.count;
        gl_mol->buffer[GL_BUFFER_ATOM_POSITION]        = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 3,  NULL, GL_DYNAMIC_DRAW);
        gl_mol->buffer[GL_BUFFER_ATOM_POSITION_PREV]   = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 3,  NULL, GL_DYNAMIC_COPY);
        gl_mol->buffer[GL_BUFFER_ATOM_VELOCITY]        = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 3,  NULL, GL_DYNAMIC_COPY);
        gl_mol->buffer[GL_BUFFER_ATOM_RADIUS]          = gl_buffer_create(gl_mol->atom_count * sizeof(float) * 1,  NULL, GL_STATIC_DRAW);
        gl_mol->buffer[GL_BUFFER_ATOM_FLAGS]           = gl_buffer_create(gl_mol->atom_count * sizeof(uint8_t), NULL, GL_STATIC_DRAW);

        if (mol->atom.x && mol->atom.y && mol->atom.z) {
            md_gl_molecule_set_atom_position(ext_mol, 0, gl_mol->atom_count, mol->atom.x, mol->atom.y, mol->atom.z, 0);
        }
        md_gl_molecule_zero_velocity(ext_mol);

        if (mol->atom.radius) md_gl_molecule_set_atom_radius(ext_mol, 0, gl_mol->atom_count, mol->atom.radius, 0);
        //if (mol->atom.flags)  md_gl_molecule_set_atom_flags(ext_mol,  0, gl_mol->atom_count, mol->atom.flags, 0);

        gl_mol->residue_count = (uint32_t)mol->residue.count;
        //gl_mol->buffer[GL_BUFFER_RESIDUE_ATOM_RANGE]          = gl_buffer_create(gl_mol->residue_count * sizeof(md_range_t),   NULL, GL_STATIC_DRAW);
        //gl_mol->buffer[GL_BUFFER_RESIDUE_AABB]                = gl_buffer_create(gl_mol->residue_count * sizeof(float) * 6,    NULL, GL_DYNAMIC_COPY);
        //gl_mol->buffer[GL_BUFFER_RESIDUE_VISIBLE]             = gl_buffer_create(gl_mol->residue_count * sizeof(int),          NULL, GL_DYNAMIC_COPY);

        //if (mol->residue.atom_offset)           gl_buffer_set_sub_data(gl_mol->buffer[GL_BUFFER_RESIDUE_ATOM_RANGE], 0, gl_mol->residue_count * sizeof(uint32_t) * 2, mol->residue.atom_range);
        //if (desc->residue.backbone_atoms)       gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_RESIDUE_BACKBONE_ATOMS], 0, mol->residue_count * sizeof(uint8_t) * 4, desc->residue.backbone_atoms);

        if (mol->backbone.range.count > 0 && mol->backbone.range.offset && mol->backbone.atoms && mol->backbone.secondary_structure) {
            uint32_t backbone_residue_count = 0;
            uint32_t backbone_spline_count = 0;
            for (uint32_t i = 0; i < (uint32_t)mol->backbone.range.count; ++i) {
                uint32_t res_count = mol->backbone.range.offset[i+1] - mol->backbone.range.offset[i];
                backbone_residue_count += res_count;
                backbone_spline_count += (res_count - 1) * MD_GL_SPLINE_SUBDIVISION_COUNT + 1; // +1 For the last point
            }

            const uint32_t backbone_count                     = backbone_residue_count;
            const uint32_t backbone_control_point_data_count  = backbone_residue_count;
            const uint32_t backbone_control_point_index_count = backbone_residue_count + (uint32_t)mol->backbone.range.count * (2 + 1); // Duplicate pair first and last in each chain for adjacency + primitive restart between
            const uint32_t backbone_spline_data_count         = backbone_spline_count;
            const uint32_t backbone_spline_index_count        = backbone_spline_count + (uint32_t)mol->backbone.range.count * (1); // Primitive restart between chains

            gl_mol->buffer[GL_BUFFER_BACKBONE_DATA]                = gl_buffer_create(backbone_count                     * sizeof(gl_backbone_data_t),         NULL, GL_STATIC_DRAW);
            gl_mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE] = gl_buffer_create(backbone_count                     * sizeof(md_secondary_structure_t),   NULL, GL_DYNAMIC_DRAW);
            gl_mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA]  = gl_buffer_create(backbone_control_point_data_count  * sizeof(gl_control_point_t),         NULL, GL_DYNAMIC_COPY);
            gl_mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX] = gl_buffer_create(backbone_control_point_index_count * sizeof(uint32_t),                   NULL, GL_STATIC_DRAW);
            gl_mol->buffer[GL_BUFFER_BACKBONE_SPLINE_DATA]         = gl_buffer_create(backbone_spline_data_count         * sizeof(gl_control_point_t),         NULL, GL_DYNAMIC_COPY);
            gl_mol->buffer[GL_BUFFER_BACKBONE_SPLINE_INDEX]        = gl_buffer_create(backbone_spline_index_count        * sizeof(uint32_t),                   NULL, GL_STATIC_DRAW);

            //gl_buffer_set_sub_data(mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE], 0, desc->backbone.count * sizeof(uint8_t) * 4, desc->backbone.secondary_structure);

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_BACKBONE_DATA].id);
            gl_backbone_data_t* backbone_data = (gl_backbone_data_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (backbone_data) {
                uint32_t idx = 0;
                for (uint32_t i = 0; i < (uint32_t)mol->backbone.range.count; ++i) {
                    uint32_t beg = (uint32_t)mol->backbone.range.offset[i];
                    uint32_t end = (uint32_t)mol->backbone.range.offset[i+1];
                    for (uint32_t j = beg; j < end; ++j) {
                        const uint32_t flags = (j == beg ? 1 : 0) | (j == end - 1 ? 2 : 0);
                        backbone_data[idx].residue_idx = j;
                        backbone_data[idx].segment_idx = j - beg;
                        backbone_data[idx].ca_idx = mol->backbone.atoms[j].ca;
                        backbone_data[idx].c_idx  = mol->backbone.atoms[j].c;
                        backbone_data[idx].o_idx  = mol->backbone.atoms[j].o;
                        backbone_data[idx].flags  = flags;
                        ++idx;
                    }
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                return false;
            }

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE].id);
            uint32_t* secondary_structure = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (secondary_structure) {
                uint32_t idx = 0;
                for (uint32_t i = 0; i < (uint32_t)mol->backbone.range.count; ++i) {
                    for (uint32_t j = mol->backbone.range.offset[i]; j < mol->backbone.range.offset[i+1]; ++j) {
                        secondary_structure[idx++] = mol->backbone.secondary_structure[j];
                    }
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            } else {
                return false;
            }

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX].id);
            uint32_t* control_point_index = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (control_point_index) {
                uint32_t idx = 0;
                uint32_t len = 0;
                for (uint32_t i = 0; i < (uint32_t)mol->backbone.range.count; ++i) {
                    control_point_index[len++] = idx;
                    for (uint32_t j = mol->backbone.range.offset[i]; j < mol->backbone.range.offset[i+1]; ++j) {
                        control_point_index[len++] = idx++;
                    }
                    control_point_index[len++] = idx-1;
                    control_point_index[len++] = 0xFFFFFFFFU;
                }
                glUnmapBuffer(GL_ARRAY_BUFFER);
            }
            else {
                return false;
            }

            glBindBuffer(GL_ARRAY_BUFFER, gl_mol->buffer[GL_BUFFER_BACKBONE_SPLINE_INDEX].id);
            uint32_t* spline_index = (uint32_t*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
            if (spline_index) {
                uint32_t idx = 0;
                uint32_t len = 0;
                uint32_t range_count = (uint32_t)mol->backbone.range.count;
                for (uint32_t i = 0; i < range_count; ++i) {
                    uint32_t res_count = mol->backbone.range.offset[i+1] - mol->backbone.range.offset[i];
                    if (res_count > 0) {
                        for (uint32_t j = 0; j < (res_count - 1) * MD_GL_SPLINE_SUBDIVISION_COUNT + 1; ++j) {
                            spline_index[len++] = idx++;
                        }
                    }
                    spline_index[len++] = 0xFFFFFFFFU;
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

        gl_mol->bond_count = (uint32_t)mol->bond.count;
        gl_mol->buffer[GL_BUFFER_BOND_ATOM_INDICES] = gl_buffer_create(gl_mol->bond_count * sizeof(uint32_t) * 2, NULL, GL_DYNAMIC_COPY);

        if (mol->bond.pairs) {
            md_gl_molecule_set_bonds(ext_mol, 0, gl_mol->bond_count, mol->bond.pairs, sizeof(md_bond_pair_t));
        }
       
        gl_mol->magic = MAGIC;
        return true;
    }
    return false;
}

bool md_gl_molecule_free(md_gl_molecule_t* ext_mol) {
    if (ext_mol) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        for (uint32_t i = 0; i < GL_BUFFER_COUNT; ++i) {
            gl_buffer_conditional_delete(&mol->buffer[i]);
        }
        uint32_t new_version = mol->version + 1;
        MEMSET(mol, 0, sizeof(md_gl_molecule_t));
        mol->version = new_version;
        return true;
    }
    return false;
}

bool md_gl_representation_init(md_gl_representation_t* ext_rep, const md_gl_molecule_t* ext_mol) {   
    if (ext_rep && ext_mol) {
        internal_mol_t* mol = (internal_mol_t*)ext_mol;
        if (mol->atom_count == 0) {
            MD_LOG_ERROR("Supplied molecule has no atoms.");
            return false;
        }
        if (mol->magic != MAGIC) {
            MD_LOG_ERROR("Supplied molecule is corrupt.");
            return false;
        }

        internal_rep_t* rep = (internal_rep_t*)ext_rep;
        MEMSET(rep, 0, sizeof(internal_rep_t));
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
        if (rep->magic != MAGIC) {
            MD_LOG_ERROR("Attempting to free representation which have not been initialized.");
            return false;
        }
        gl_buffer_conditional_delete(&rep->color);
        return true;
    }
    return false;
}

static bool compute_spline(const internal_mol_t* mol);

static bool draw_space_fill(gl_program_t program, const internal_rep_t* rep, float scale);
static bool draw_licorice  (gl_program_t program, const internal_rep_t* rep, float radius, float max_length);
static bool draw_ribbons   (gl_program_t program, const internal_rep_t* rep, float width_scale, float thickness_scale);
static bool draw_cartoon   (gl_program_t program, const internal_rep_t* rep, float coil_scale, float helix_scale, float ribbon_scale);

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
}

typedef struct draw_mol_t {
    const internal_mol_t* mol;
    uint32_t flags;
} draw_mol_t;

bool md_gl_draw(const md_gl_draw_args_t* args) {
    if (!args) {
        MD_LOG_ERROR("draw args object was NULL");
        return false;
    }
    if (!validate_context()) return false;

    if (!args->shaders) {
        MD_LOG_ERROR("shaders object was NULL");
        return false;
    }
    internal_shaders_t* shaders = (internal_shaders_t*)args->shaders;

    PUSH_GPU_SECTION("MOLD DRAW")
            
    gl_ubo_base_t ubo_data = {0};
    init_ubo_base_data(&ubo_data, args, NULL);

    gl_buffer_set_sub_data(ctx.ubo, 0, sizeof(ubo_data), &ubo_data);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ctx.ubo.id);

    size_t temp_pos = md_temp_get_pos();
    md_allocator_i* alloc = md_get_temp_allocator();

    md_array(md_gl_draw_op_t const*) draw_ops = 0;
    md_array(draw_mol_t) draw_mols = 0;
        
    for (uint32_t i = 0; i < args->draw_operations.count; i++) {
        const md_gl_draw_op_t* draw_op = &args->draw_operations.ops[i];
        const internal_rep_t* rep = (internal_rep_t*)draw_op->rep;
        if (rep && rep->mol) {
            if (rep->magic != MAGIC) {
                MD_LOG_ERROR("Representation is corrupt");
                continue;
            }

            if (rep->mol_version != rep->mol->version) {
                MD_LOG_ERROR("Representation was not created for the current molecule");
                continue;
            }

            if (draw_op->type == MD_GL_REP_RIBBONS || draw_op->type == MD_GL_REP_CARTOON) {
                if (!(rep->mol->flags & MOL_FLAG_HAS_BACKBONE)) {
                    MD_LOG_ERROR("Molecule is missing a backbone, therefore the chosen representation cannot be drawn");
                    continue;
                }
            }

            md_array_push(draw_ops, draw_op, alloc);
                
            int64_t mol_idx = -1;
            for (size_t j = 0; j < md_array_size(draw_mols); j++) {
                if (draw_mols[j].mol == rep->mol) {
                    mol_idx = j;
                    break;
                }
            }
            // Not found
            if (mol_idx == -1) {
                mol_idx = md_array_size(draw_mols);
                md_array_push(draw_mols, ((draw_mol_t){ rep->mol, 0 }), alloc);
            }

            if (draw_op->type == MD_GL_REP_RIBBONS || draw_op->type == MD_GL_REP_CARTOON) {
                draw_mols[mol_idx].flags |= DRAW_FLAG_COMPUTE_BACKBONE_SPLINE;
            }
        }
    }
    
    //qsort((void*)draw_ent, draw_ent_count, sizeof(draw_entity_t), compare_draw_ent);
            
    PUSH_GPU_SECTION("COMPUTE SPLINE")
    for (uint32_t i = 0; i < md_array_size(draw_mols); i++) {
        if (draw_mols[i].flags & DRAW_FLAG_COMPUTE_BACKBONE_SPLINE) {
            compute_spline(draw_mols[i].mol);
        }
    }
    POP_GPU_SECTION()
        
    //bool using_internal_depth = false;

    int program_permutation = 0;
    if (is_ortho_proj_matrix(ubo_data.view_transform.view_to_clip)) program_permutation |= PERMUTATION_BIT_ORTHO;

    uint32_t index_base[2] = {0,0};
        
    PUSH_GPU_SECTION("DRAW REPRESENTATIONS")
    for (size_t i = 0; i < md_array_size(draw_ops); i++) {
        const md_gl_draw_op_t* draw_op = draw_ops[i];
        const internal_rep_t* rep  = (const internal_rep_t*)draw_ops[i]->rep;
        const internal_mol_t* mol  = (const internal_mol_t*)rep->mol; 
        const mat4_t* model_matrix = (const mat4_t*)draw_ops[i]->model_matrix;
        float scale = 1.0f;

        if (model_matrix) {
            // If we have a model matrix, we need to recompute the entire matrix stack...
            gl_ubo_base_t ubo_tmp = {0};
            init_ubo_base_data(&ubo_tmp, args, model_matrix);
            gl_buffer_set_sub_data(ctx.ubo, 0, sizeof(gl_view_transform_t), &ubo_tmp);
            const vec3_t model_scale = {
                vec3_length(vec3_from_vec4(ubo_tmp.view_transform.world_to_view.col[0])),
                vec3_length(vec3_from_vec4(ubo_tmp.view_transform.world_to_view.col[1])),
                vec3_length(vec3_from_vec4(ubo_tmp.view_transform.world_to_view.col[2])),
            };
            const float mean = (model_scale.x + model_scale.y + model_scale.z) / 3.0f;

            // Non uniform scale is not supported, it will cause rendering artifact for radius scaling parameters which are assumed to be uniform.
            scale = mean;
        }

        if (mol->atom_index_base != index_base[0] || mol->bond_index_base != index_base[1]) {
            index_base[0] = mol->atom_index_base;
            index_base[1] = mol->bond_index_base;
            gl_buffer_set_sub_data(ctx.ubo, offsetof(gl_ubo_base_t, atom_index_base), sizeof(index_base), &index_base);
        }

        // Maximum bond length in units (Ångström assumed)
        const float max_length = 5.0f;

        switch (draw_op->type) {
        case MD_GL_REP_SPACE_FILL:
            draw_space_fill(shaders->spacefill[program_permutation], rep, scale * draw_op->args.space_fill.radius_scale);
            break;
        case MD_GL_REP_LICORICE:
            draw_licorice(shaders->licorice[program_permutation], rep, 0.2f * scale * draw_op->args.licorice.radius, max_length);
            break;
        case MD_GL_REP_BALL_AND_STICK:
            draw_licorice(shaders->licorice[program_permutation],    rep, 0.2f * scale * draw_op->args.ball_and_stick.stick_radius, max_length);
            draw_space_fill(shaders->spacefill[program_permutation], rep, 0.2f  * scale * draw_op->args.ball_and_stick.ball_scale);
            break;
        case MD_GL_REP_RIBBONS:
            draw_ribbons(shaders->ribbons[program_permutation], rep, scale * draw_op->args.ribbons.width_scale, scale * draw_op->args.ribbons.thickness_scale);
            break;
        case MD_GL_REP_CARTOON:
            draw_cartoon(shaders->cartoon[program_permutation], rep, scale * draw_op->args.cartoon.coil_scale, draw_op->args.cartoon.helix_scale, draw_op->args.cartoon.sheet_scale);
            break;
        default:
            MD_LOG_ERROR("Representation had unexpected type");
            return false;
            break;
        }

        if (model_matrix) {
            // Reset matrix stack
            gl_buffer_set_sub_data(ctx.ubo, 0, sizeof(gl_view_transform_t), &ubo_data);
        }
    }
    POP_GPU_SECTION()

    POP_GPU_SECTION()

    md_temp_set_pos_back(temp_pos);
    
    return true;
}

static bool draw_space_fill(gl_program_t program, const internal_rep_t* rep, float scale) {
    ASSERT(rep);
    ASSERT(rep->mol);
    ASSERT(rep->mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_ATOM_RADIUS].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(rep->color.id);

    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(scale), &scale);

    glBindVertexArray(ctx.vao);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_ATOM_RADIUS].id);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    glVertexAttribIPointer(3, 1, GL_UNSIGNED_BYTE, 0, 0);

    glEnableVertexAttribArray(4);
    glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
    glVertexAttribPointer(4, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
    
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glUseProgram(program.id);
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

static bool draw_licorice(gl_program_t program, const internal_rep_t* rep, float radius, float max_length) {
    ASSERT(rep);
    ASSERT(rep->mol);
    ASSERT(rep->mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_BOND_ATOM_INDICES].id);
    ASSERT(rep->color.id);

    if (max_length == 0) {
        max_length = 1000.0f;
    }

    struct {
		float radius;
		float max_d2;
	} params = { radius, max_length * max_length };

    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(params), &params);

    glBindVertexArray(ctx.vao);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    glVertexAttribIPointer(2, 1, GL_UNSIGNED_BYTE, 0, 0);

    glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, rep->color.id);
    glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_BOND_ATOM_INDICES].id);

    glUseProgram(program.id);
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

bool draw_ribbons(gl_program_t program, const internal_rep_t* rep, float width_scale, float thickness_scale) {
    ASSERT(rep);
    ASSERT(rep->mol);
    ASSERT(rep->mol->buffer[GL_BUFFER_BACKBONE_SPLINE_DATA].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_BACKBONE_SPLINE_INDEX].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(rep->color.id);

    const float profile_scale[2] = {
        width_scale,
        thickness_scale * 0.1f,
    };
    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(profile_scale), &profile_scale);

    glBindVertexArray(ctx.vao);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_BACKBONE_SPLINE_DATA].id);

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

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_BACKBONE_SPLINE_INDEX].id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, rep->color.id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, rep->mol->buffer[GL_BUFFER_ATOM_FLAGS].id);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);

    glUseProgram(program.id);
    glUniform1i(glGetUniformLocation(program.id, "u_atom_color_buffer"), 0);
    glUniform1i(glGetUniformLocation(program.id, "u_atom_flags_buffer"), 1);
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

static bool draw_cartoon(gl_program_t program, const internal_rep_t* rep, float coil_scale, float helix_scale, float sheet_scale) {
    ASSERT(rep);
    ASSERT(rep->mol);
    ASSERT(rep->mol->buffer[GL_BUFFER_BACKBONE_SPLINE_DATA].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_BACKBONE_SPLINE_INDEX].id);
    ASSERT(rep->mol->buffer[GL_BUFFER_ATOM_FLAGS].id);
    ASSERT(rep->color.id);

    const float profile_scale[4] = {
        coil_scale,
        helix_scale,
        sheet_scale,
        0
    };
    gl_buffer_set_sub_data(ctx.ubo, sizeof(gl_ubo_base_t), sizeof(profile_scale), &profile_scale);

    glBindVertexArray(ctx.vao);
    glBindBuffer(GL_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_BACKBONE_SPLINE_DATA].id);

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

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rep->mol->buffer[GL_BUFFER_BACKBONE_SPLINE_INDEX].id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, rep->color.id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, rep->mol->buffer[GL_BUFFER_ATOM_FLAGS].id);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);

    glUseProgram(program.id);
    glUniform1i(glGetUniformLocation(program.id, "u_atom_color_buffer"), 0);
    glUniform1i(glGetUniformLocation(program.id, "u_atom_flags_buffer"), 1);
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

static bool compute_spline(const internal_mol_t* mol) {
    ASSERT(ctx.program[GL_PROGRAM_COMPUTE_SPLINE].id);
    ASSERT(mol);

    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_DATA].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_POSITION].id);
    ASSERT(mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX].id);
    ASSERT(mol->buffer[GL_BUFFER_BACKBONE_SPLINE_DATA].id);

    if (mol->buffer[GL_BUFFER_BACKBONE_DATA].id == 0) {
        MD_LOG_ERROR("Backbone data buffer is zero, which is required to compute the backbone. Is the molecule missing a backbone?");
        return false;
    }

    glEnable(GL_RASTERIZER_DISCARD);
    glBindVertexArray(ctx.vao);

    // Step 1: Extract control points
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_BACKBONE_DATA].id);

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

    glEnableVertexAttribArray(5);
    glVertexAttribIPointer(5, 1, GL_UNSIGNED_INT, sizeof(gl_backbone_data_t), (const void*)offsetof(gl_backbone_data_t, flags));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_0].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_POSITION].id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_1].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, mol->buffer[GL_BUFFER_ATOM_VELOCITY].id);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_BUFFER, ctx.texture[GL_TEXTURE_BUFFER_2].id);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA8, mol->buffer[GL_BUFFER_BACKBONE_SECONDARY_STRUCTURE].id);

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

    // Step 2: Compute splines using control points
    glBindBuffer(GL_ARRAY_BUFFER, mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_DATA].id);

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

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mol->buffer[GL_BUFFER_BACKBONE_CONTROL_POINT_INDEX].id);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);

    glUseProgram(ctx.program[GL_PROGRAM_COMPUTE_SPLINE].id);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, mol->buffer[GL_BUFFER_BACKBONE_SPLINE_DATA].id);
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
