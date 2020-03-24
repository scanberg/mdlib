#include "mold.h"

#if defined(_WIN32) && !defined(APIENTRY) && !defined(__CYGWIN__) && !defined(__SCITECH_SNAP__)
#define APIENTRY __stdcall
#endif
#ifndef APIENTRY
#define APIENTRY
#endif
#ifndef APIENTRYP
#define APIENTRYP APIENTRY *
#endif
#ifndef GLAPI
#define GLAPI extern
#endif

typedef void GLvoid;
typedef unsigned int GLenum;
typedef float GLfloat;
typedef int GLint;
typedef int GLsizei;
typedef unsigned int GLbitfield;
typedef double GLdouble;
typedef unsigned int GLuint;
typedef unsigned char GLboolean;
typedef unsigned char GLubyte;
typedef ptrdiff_t GLsizeiptr;
typedef ptrdiff_t GLintptr;

#define GL_MAJOR_VERSION                  0x821B
#define GL_MINOR_VERSION                  0x821C

#define GL_ARRAY_BUFFER                   0x8892
#define GL_ELEMENT_ARRAY_BUFFER           0x8893
#define GL_READ_ONLY                      0x88B8
#define GL_WRITE_ONLY                     0x88B9
#define GL_READ_WRITE                     0x88BA
#define GL_BUFFER_ACCESS                  0x88BB
#define GL_STREAM_DRAW                    0x88E0
#define GL_STREAM_READ                    0x88E1
#define GL_STREAM_COPY                    0x88E2
#define GL_STATIC_DRAW                    0x88E4
#define GL_STATIC_READ                    0x88E5
#define GL_STATIC_COPY                    0x88E6
#define GL_DYNAMIC_DRAW                   0x88E8
#define GL_DYNAMIC_READ                   0x88E9
#define GL_DYNAMIC_COPY                   0x88EA

// @NOTE: We rely on the fact that these functions are readily supplied
GLAPI void APIENTRY glDepthMask (GLboolean flag);
GLAPI void APIENTRY glGetInteger(GLenum pname, GLint *data);

GLAPI void APIENTRY glBindVertexArray (GLuint array);
GLAPI void APIENTRY glDeleteVertexArrays (GLsizei n, const GLuint *arrays);
GLAPI void APIENTRY glGenVertexArrays (GLsizei n, GLuint *arrays);
GLAPI void APIENTRY glDisableVertexAttribArray (GLuint index);
GLAPI void APIENTRY glEnableVertexAttribArray (GLuint index);
GLAPI GLboolean APIENTRY glIsVertexArray (GLuint array);

GLAPI void APIENTRY glBindBuffer (GLenum target, GLuint buffer);
GLAPI void APIENTRY glDeleteBuffers (GLsizei n, const GLuint *buffers);
GLAPI void APIENTRY glGenBuffers (GLsizei n, GLuint *buffers);
GLAPI GLboolean APIENTRY glIsBuffer (GLuint buffer);
GLAPI void APIENTRY glBufferData (GLenum target, GLsizeiptr size, const void *data, GLenum usage);
GLAPI void APIENTRY glBufferSubData (GLenum target, GLintptr offset, GLsizeiptr size, const void *data);
GLAPI void APIENTRY glGetBufferSubData (GLenum target, GLintptr offset, GLsizeiptr size, void *data);
GLAPI void *APIENTRY glMapBuffer (GLenum target, GLenum access);
GLAPI GLboolean APIENTRY glUnmapBuffer (GLenum target);

GLAPI void APIENTRY glDeleteShader (GLuint shader);
GLAPI GLuint APIENTRY glCreateShader (GLenum type);
GLAPI GLuint APIENTRY glCreateProgram (void);
GLAPI void APIENTRY glDeleteProgram (GLuint program);
GLAPI void APIENTRY glCompileShader (GLuint shader);
GLAPI void APIENTRY glDetachShader (GLuint program, GLuint shader);
GLAPI GLboolean APIENTRY glIsProgram (GLuint program);
GLAPI GLboolean APIENTRY glIsShader (GLuint shader);

GLAPI GLboolean APIENTRY glIsFramebuffer (GLuint framebuffer);
GLAPI void APIENTRY glBindFramebuffer (GLenum target, GLuint framebuffer);
GLAPI void APIENTRY glDeleteFramebuffers (GLsizei n, const GLuint *framebuffers);
GLAPI void APIENTRY glGenFramebuffers (GLsizei n, GLuint *framebuffers);
GLAPI GLenum APIENTRY glCheckFramebufferStatus (GLenum target);
GLAPI void APIENTRY glFramebufferTexture2D (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
GLAPI void APIENTRY glFramebufferRenderbuffer (GLenum target, GLenum attachment, GLenum renderbuffertarget, GLuint renderbuffer);

#define MOLD_MAX_REP_SLOTS 32
#define MOLD_MAX_MOL_SLOTS 32
#define MOLD_HANDLE_BASE_MOL 0xb00bf00d
#define MOLD_HANDLE_BASE_REP 0xb00bd00d

static uint32_t get_handle_base(void* handle) {
    return (uint32_t)((uintptr_t)handle & 0xFFFFFFFFU);
}

static uint32_t get_handle_index(void* handle) {
    return (uint32_t)(((uintptr_t)handle >> 32) & 0xFFFFFFFFU);
}

typedef struct mold_aabb {
    float min[3];
    float max[3];
} mold_aabb;

typedef struct mold_gl_buffer {
    GLuint  id;
} mold_gl_buffer;

typedef struct mold_gl_rep {
    mold_rep_type type;
    mold_rep_args args;

    mold_range atom_range;
    mold_gl_buffer color;
} mold_gl_rep;

typedef struct mold_gl_mol {
    struct {
        uint32_t count;
        mold_gl_buffer position; 
        mold_gl_buffer radius;
    } atom;

    struct {
        uint32_t count;
        mold_gl_buffer atom_range;
        mold_gl_buffer aabb;
    } residue;

    struct {
        uint32_t count;
        mold_gl_buffer data;
    } covalent_bond;

    struct {
        struct {
            uint32_t count;
            mold_gl_buffer atom_index;
            mold_gl_buffer secondary_structure;
        } segment;

        struct {
            uint32_t count;
            mold_gl_buffer segment_range;
        } sequence;
    } backbone;
} mold_gl_mol;

typedef uint32_t mold_gl_version;
enum mold_gl_version_ {
    MOLD_GL_VERSION_UNKNOWN = 0,
    MOLD_GL_VERSION_330,
    MOLD_GL_VERSION_430
};

typedef struct mold_gl_context {
    mold_gl_version version;

    struct {
        uint32_t count;
        mold_gl_mol data[MOLD_MAX_MOL_SLOTS];
    } molecule;

    struct {
        uint32_t count;
        mold_gl_rep data[MOLD_MAX_REP_SLOTS];
    } representation;

    GLuint vao;
    GLuint ubo;
    GLuint fbo;
} mold_gl_context;

static mold_gl_context mold_context = {0};

static mold_result mold_validate_context(mold_gl_context* ctx) {
    mold_result res = {MOLD_SUCCESS};
    if (!ctx) {
        res.error_code = MOLD_ERROR_ARGUMENT_IS_NULL;
        res.error_str = "Context is NULL";
    }
    else if (!glIsVertexArray(ctx->vao) || !glIsBuffer(ctx->ubo) || !glIsFramebuffer(ctx->fbo)) {
        res.error_code = MOLD_ERROR_CONTEXT_NOT_INITIALIZED;
        res.error_str = "Internal Context GL resources are not initialized or valid";
    }
    return res;
}

static void* mold_create_handle(uint32_t handle_base, uint32_t index) {
    return (void*)((uintptr_t)handle_base & ((uintptr_t)index << 32));
}

static mold_result mold_validate_mol_handle(mold_molecule mol) {
    mold_result res = {MOLD_SUCCESS};
    if (!mol) {
        res.error_code = MOLD_ERROR_ARGUMENT_IS_NULL;
        res.error_str = "Molecule handle is NULL";
    }
    else if (get_handle_base(mol) != MOLD_HANDLE_BASE_MOL) {
        res.error_code = MOLD_ERROR_INVALID_HANDLE;
        res.error_str = "Molecule handle is not valid";
    }
    else if (get_handle_index(mol) >= MOLD_MAX_MOL_SLOTS) {
        res.error_code = MOLD_ERROR_MOLECULE_SLOT_OUT_OF_RANGE;
        res.error_str ="Molecule handle index is out of range";
    }
    return res;
}

static mold_result mold_validate_rep_handle(mold_representation rep) {
    mold_result res = {MOLD_SUCCESS};
    if (!rep) {
        res.error_code = MOLD_ERROR_ARGUMENT_IS_NULL;
        res.error_str = "Representation handle is NULL";
    }
    else if (get_handle_base(rep) != MOLD_HANDLE_BASE_REP) {
        res.error_code = MOLD_ERROR_INVALID_HANDLE;
        res.error_str = "Representation handle is not valid";
    }
    else if (get_handle_index(rep) >= MOLD_MAX_REP_SLOTS) {
        res.error_code = MOLD_ERROR_REPRESENTATION_SLOT_OUT_OF_RANGE;
        res.error_str ="Representation handle index is out of range";
    }
    return res;
}

static mold_result mold_get_gl_mol(mold_gl_mol** gl_mol, mold_molecule mol) {
    mold_result res = mold_validate_context(&mold_context);
    if (mold_success(res)) {
        res = mold_validate_mol_handle(mol);
        if (mold_success(res)) {
            *gl_mol = mold_context.molecule.data + get_handle_index(mol);
        }
    }
    return res;
}

static mold_result mold_get_gl_rep(mold_gl_rep** gl_rep, mold_representation rep) {
    mold_result res = mold_validate_context(&mold_context);
    if (mold_success(res)) {
        res = mold_validate_rep_handle(rep);
        if (mold_success(res)) {
            *gl_rep = mold_context.representation.data + get_handle_index(rep);
        }
    }
    return res;
}

static mold_result mold_gl_init_common_resources(mold_gl_context* ctx) {
    glGenVertexArrays(1, &ctx->vao);
    glGenBuffers(1, &ctx->ubo);
    glGenFramebuffers(1, &ctx->fbo);

    mold_result res = {MOLD_SUCCESS, ""};
    return res;
}

static mold_result mold_gl_free_common_resources(mold_gl_context* ctx) {
    glDeleteVertexArrays(1, &ctx->vao);
    glDeleteBuffers(1, &ctx->ubo);
    glDeleteFramebuffers(1, &ctx->fbo);

    mold_result res = {MOLD_SUCCESS, ""};
    return res;
}

static mold_result mold_set_atom_position_data(mold_gl_mol* mol, const mold_position pos_xyz[], uint32_t count, uint32_t offset) {
    glBindBuffer(GL_ARRAY_BUFFER, mol->atom.position.id);
    glBufferSubData(GL_ARRAY_BUFFER, offset * sizeof(mold_position), count * sizeof(mold_position), pos_xyz);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static mold_result mold_set_atom_position_data_soa(mold_gl_mol* mol, const float pos_x[], const float pos_y[], const float pos_z[], uint32_t count, uint32_t offset) {
    glBindBuffer(GL_ARRAY_BUFFER, mol->atom.position.id);
    void* ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    if (ptr) {
        mold_position* pos = (mold_position*)ptr;
        for (uint32_t i = offset; i < count; i++) {
            pos[i].x = pos_x[i];
            pos[i].y = pos_y[i];
            pos[i].z = pos_z[i];
        }
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static mold_result mold_gl_set_atom_radius_data(mold_gl_mol* mol, const float radius[], uint32_t count, uint32_t offset) {
    glBindBuffer(GL_ARRAY_BUFFER, mol->atom.radius.id);
    glBufferSubData(GL_ARRAY_BUFFER, offset * sizeof(mold_radius), count * sizeof(mold_radius), radius);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static mold_result mold_gl_set_residue_atom_range_data(mold_gl_mol* mol, const mold_range atom_range[], uint32_t count, uint32_t offset) {
    glBindBuffer(GL_ARRAY_BUFFER, mol->residue.atom_range.id);
    glBufferSubData(GL_ARRAY_BUFFER, offset * sizeof(mold_range), count * sizeof(mold_range), atom_range);
    glBindBuffer(GL_ARRAY_BUFFER, 0);}

static mold_result mold_gl_set_covalent_bond_data(mold_gl_mol* mol, const mold_bond bond[], uint32_t count, uint32_t offset) {
    glBindBuffer(GL_ARRAY_BUFFER, mol->covalent_bond.data.id);
    glBufferSubData(GL_ARRAY_BUFFER, offset * sizeof(mold_bond), count * sizeof(mold_bond), bond);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static mold_result mold_gl_set_backbone_segment_index_data(mold_gl_mol* mol, const mold_segment_indices segment_indices[], uint32_t count, uint32_t offset) {
    glBindBuffer(GL_ARRAY_BUFFER, mol->backbone.segment.atom_index.id);
    glBufferSubData(GL_ARRAY_BUFFER, offset * sizeof(mold_segment_indices), count * sizeof(mold_segment_indices), segment_indices);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static mold_result mold_gl_set_backbone_segment_secondary_structure_data(mold_gl_mol* mol, const mold_secondary_structure secondary_structure[], uint32_t count, uint32_t offset) {
    glBindBuffer(GL_ARRAY_BUFFER, mol->backbone.segment.secondary_structure.id);
    glBufferSubData(GL_ARRAY_BUFFER, offset * sizeof(mold_secondary_structure), count * sizeof(mold_secondary_structure), secondary_structure);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static mold_result mold_gl_set_backbone_sequence_data(mold_gl_mol* mol, const mold_range sequence[], uint32_t count, uint32_t offset) {
    glBindBuffer(GL_ARRAY_BUFFER, mol->backbone.sequence.segment_range.id);
    glBufferSubData(GL_ARRAY_BUFFER, offset * sizeof(mold_range), count * sizeof(mold_range), sequence);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static mold_result mold_gl_init_mol_resources(mold_gl_mol* mol, const mold_molecule_desc* desc) {
    glGenBuffers(1, &mol->atom.position.id);
    glGenBuffers(1, &mol->atom.radius.id);
    mol->atom.count = desc->atom.count;
    if (mol->atom.count > 0) {
        glBindBuffer(GL_ARRAY_BUFFER, mol->atom.position.id);
        glBufferData(GL_ARRAY_BUFFER, mol->atom.count * sizeof(mold_position), NULL, GL_DYNAMIC_DRAW);
        if (desc->atom.position_data) {
            mold_gl_set_atom_position_data(mol, desc->atom.position_data, desc->atom.count, 0);
        }

        glBindBuffer(GL_ARRAY_BUFFER, mol->atom.radius.id);
        glBufferData(GL_ARRAY_BUFFER, mol->atom.count * sizeof(float), NULL, GL_STATIC_DRAW);
        if (desc->atom.radius_data) {
            mold_gl_set_atom_radius_data(mol, desc->atom.radius_data, desc->atom.count, 0);
        }
    }

    glGenBuffers(1, &mol->residue.atom_range.id);
    glGenBuffers(1, &mol->residue.aabb.id);
    mol->residue.count = desc->residue.count;
    if (mol->residue.count > 0) {
        glBindBuffer(GL_ARRAY_BUFFER, mol->residue.atom_range.id);
        glBufferData(GL_ARRAY_BUFFER, mol->residue.count * sizeof(mold_range), NULL, GL_STATIC_DRAW);
        if (desc->residue.atom_range_data) {
            mold_gl_set_residue_atom_range_data(mol, desc->residue.atom_range_data, desc->residue.count, 0);
        }

        glBindBuffer(GL_ARRAY_BUFFER, mol->residue.aabb.id);
        glBufferData(GL_ARRAY_BUFFER, mol->residue.count * sizeof(mold_aabb), NULL, GL_STATIC_DRAW);
    }

    glGenBuffers(1, &mol->covalent_bond.data.id);
    mol->covalent_bond.count = desc->covalent_bond.count;
    if (mol->covalent_bond.count > 0) {
        glBindBuffer(GL_ARRAY_BUFFER, mol->covalent_bond.data.id);
        glBufferData(GL_ARRAY_BUFFER, mol->covalent_bond.count * sizeof(mold_bond), NULL, GL_STATIC_DRAW);
        if (desc->covalent_bond.data) {
            mold_gl_set_covalent_bond_data(mol, desc->covalent_bond.data, desc->covalent_bond.count, 0);
        }
    }

    glGenBuffers(1, &mol->backbone.segment.atom_index.id);
    glGenBuffers(1, &mol->backbone.segment.secondary_structure.id);
    mol->backbone.segment.count = desc->backbone.segment.count;
    if (mol->backbone.segment.count) {
        glBindBuffer(GL_ARRAY_BUFFER, mol->backbone.segment.atom_index.id);
        glBufferData(GL_ARRAY_BUFFER, mol->backbone.segment.count * sizeof(mold_segment_indices), NULL, GL_STATIC_DRAW);
        if (desc->backbone.segment.index_data) {
            mold_gl_set_backbone_segment_index_data(mol, desc->backbone.segment.index_data, desc->backbone.segment.count, 0);
        }

        glBindBuffer(GL_ARRAY_BUFFER, mol->backbone.segment.secondary_structure.id);
        glBufferData(GL_ARRAY_BUFFER, mol->backbone.segment.count * sizeof(mold_secondary_structure), NULL, GL_STATIC_DRAW);
        if (desc->backbone.segment.index_data) {
            mold_gl_set_backbone_segment_secondary_structure_data(mol, desc->backbone.segment.secondary_structure_data, desc->backbone.segment.count, 0);
        }
    }

    glGenBuffers(1, &mol->backbone.sequence.segment_range.id);
    mol->backbone.sequence.count = desc->backbone.sequence.count;
    if (mol->backbone.sequence.count > 0) {
        glBindBuffer(GL_ARRAY_BUFFER, mol->backbone.sequence.segment_range.id);
        glBufferData(GL_ARRAY_BUFFER, mol->backbone.sequence.count * sizeof(mold_range), NULL, GL_STATIC_DRAW);
        if (desc->backbone.sequence.segment_range_data) {
            mold_gl_set_backbone_sequence_data(mol, desc->backbone.sequence.segment_range_data, desc->backbone.sequence.count, 0);
        }
    }

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    mold_result res = {MOLD_SUCCESS, ""};
    return res;
}

static mold_result mold_gl_free_mol_resources(mold_gl_mol* mol) {
    glDeleteBuffers(1, &mol->atom.position.id);
    glDeleteBuffers(1, &mol->atom.radius.id);
    glDeleteBuffers(1, &mol->residue.atom_range.id);
    glDeleteBuffers(1, &mol->residue.aabb.id);
    glDeleteBuffers(1, &mol->covalent_bond.data.id);
    glDeleteBuffers(1, &mol->backbone.segment.atom_index.id);
    glDeleteBuffers(1, &mol->backbone.segment.secondary_structure.id);
    glDeleteBuffers(1, &mol->backbone.sequence.segment_range.id);

    mold_result result = {MOLD_SUCCESS};
    return result;
}

static mold_result mold_gl_set_type_and_args(mold_gl_rep* rep, mold_rep_type type, mold_rep_args args) {
    rep->type = type;
    rep->args = args;
}

static mold_result mold_gl_set_atom_range(mold_gl_rep* rep, mold_range atom_range) {
    rep->atom_range = atom_range;
}

static mold_result mold_gl_set_color_data(mold_gl_rep* rep, const mold_color color_data[], uint32_t count, uint32_t offset) {
    rep->
}

static mold_result mold_gl_init_rep_resources(mold_gl_rep* rep, const mold_representation_desc* desc) {

}

static mold_result mold_gl_free_rep_resources(mold_gl_rep* rep) {

}

mold_result mold_initialize_context() {
    GLint major, minor;
    glGetInteger(GL_MAJOR_VERSION, &major);
    glGetInteger(GL_MINOR_VERSION, &minor);

    mold_gl_version version = MOLD_GL_VERSION_UNKNOWN;
    if (major >= 4 && minor >= 3) {
        version = MOLD_GL_VERSION_430;
    }
    else if (major >= 3 && minor >= 3) {
        version = MOLD_GL_VERSION_330;
    }
    else {
        mold_result res = {MOLD_ERROR_GL_VERSION_NOT_SUPPORTED, "OpenGL version is not supported"};
        return res;
    }

    mold_context.version = version;

    mold_result result;
    result = mold_gl_init_common_resources(&mold_context);
    if (!mold_success(result)) {
        return result;
    }

    return result;
}

mold_result mold_shutdown_context() {
    mold_result result = mold_validate_context(&mold_context);
    if (mold_success(result)) {
        for (uint32_t i = 0; i < mold_context.molecule.count; i++) {
            mold_result res = mold_gl_free_mol_resources(mold_context.molecule.data + i);
            if (!mold_success(res)) {
                result = res;
            }
        }

        for (uint32_t i = 0; i < mold_context.representation.count; i++) {
            mold_result res = mold_gl_free_mol_resources(mold_context.representation.data + i);
            if (!mold_success(res)) {
                result = res;
            }
        }

        mold_result res = mold_gl_free_common_resources(&mold_context);
        if (!mold_success(res)) {
            result = res;
        }
    }

    return result;
}

mold_result mold_create_molecule(mold_molecule* mol, const mold_molecule_desc* desc) {
    mold_result result = mold_validate_context(&mold_context);
    if (mold_success(result)) {
        if (mold_context.molecule.count < MOLD_MAX_MOL_SLOTS) {
            const uint32_t idx = mold_context.molecule.count;
            mold_gl_mol* gl_mol = &mold_context.molecule.data[idx];
            mold_result res = mold_gl_init_mol_resources(gl_mol, desc);
            if (mold_success(res)) {
                *mol = create_handle(MOLD_HANDLE_BASE_MOL, idx);
                mold_context.molecule.count++;
            } else {
                result = res;
            }
        }
        else {
            result.error_code = MOLD_ERROR_MOLECULE_SLOT_OUT_OF_RANGE;
            result.error_str = "Context is out of free molecule slots";
        }
    }
    return result;
}

mold_result mold_destroy_molecule(mold_molecule mol) {
    mold_gl_mol* gl_mol;
    mold_result result = mold_get_gl_mol(&gl_mol, mol);
    if (mold_success(result)) {
        const uint32_t idx = mold_handle_index(mol);
        mold_gl_free_mol_resources(&mold_context.molecule.data[idx]);
        mold_context.molecule.count--;
    }
    return result;
}

mold_result mold_create_representation(mold_representation* rep, const mold_representation_desc* desc) {
    mold_result result = mold_validate_context(&mold_context);
    if (mold_success(result)) {
        if (mold_context.representation.count < MOLD_MAX_REP_SLOTS) {
            const uint32_t idx = mold_context.representation.count;
            mold_gl_rep* gl_rep = &mold_context.representation.data[idx];
            mold_result res = mold_gl_init_rep_resources(gl_rep, desc);
            if (mold_success(res)) {
                *rep = create_handle(MOLD_HANDLE_BASE_REP, idx);
                mold_context.representation.count++;
            } else {
                result = res;
            }
        }
        else {
            result.error_code = MOLD_ERROR_REPRESENTATION_SLOT_OUT_OF_RANGE;
            result.error_str = "Context is out of free molecule slots";
        }
    }
    return result;
}

mold_result mold_destroy_representation(mold_representation rep) {
    mold_gl_rep* gl_rep;
    mold_result result = mold_get_gl_rep(&gl_rep, rep);
    if (mold_success(result)) {
        const uint32_t idx = mold_handle_index(rep);
        mold_gl_free_rep_resources(&mold_context.representation.data[idx]);
        mold_context.representation.count--;
    }
    return result;
}

mold_result mold_set_atom_position_data(mold_molecule mol, const mold_position pos_xyz[], uint32_t count, uint32_t offset) {
    mold_gl_mol* gl_mol;
    mold_result result = mold_get_gl_mol(&gl_mol, mol);
    if (mold_success(result)) {
        mold_set_atom_position_data(gl_mol, pos_xyz, count, offset);   
    }
    return result;
}

mold_result mold_set_atom_position_data_soa(mold_molecule mol, const float pos_x[], const float pos_y[], const float pos_z[], uint32_t count, uint32_t offset) {
    mold_gl_mol* gl_mol;
    mold_result result = mold_get_gl_mol(&gl_mol, mol);
    if (mold_success(result)) {
        mold_set_atom_position_data_soa(gl_mol, pos_x, pos_y, pos_z, count, offset);   
    }
    return result;
}

mold_result mold_set_atom_radius_data(mold_molecule mol, const float radius[], uint32_t count, uint32_t offset) {
    mold_gl_mol* gl_mol;
    mold_result result = mold_get_gl_mol(&gl_mol, mol);
    if (mold_success(result)) {
        mold_set_atom_radius_data(gl_mol, radius, count, offset);   
    }
    return result;
}

mold_result mold_set_covalent_bond_data(mold_molecule mol, const mold_bond* bond[], uint32_t count, uint32_t offset) {
    mold_gl_mol* gl_mol;
    mold_result result = mold_get_gl_mol(&gl_mol, mol);
    if (mold_success(result)) {
        mold_set_covalent_bond_data(gl_mol, bond, count, offset);
    }
    return result;
}

mold_result mold_set_backbone_segment_index_data(mold_molecule mol, const mold_segment_indices segment_indices[], uint32_t count, uint32_t offset) {
    mold_gl_mol* gl_mol;
    mold_result result = mold_get_gl_mol(&gl_mol, mol);
    if (mold_success(result)) {
        mold_set_backbone_segment_index_data(gl_mol, segment_indices, count, offset);   
    }
    return result;
}

mold_result mold_set_backbone_secondary_structure_data(mold_molecule mol, const mold_secondary_structure secondary_structure[], uint32_t count, uint32_t offset) {
    mold_gl_mol* gl_mol;
    mold_result result = mold_get_gl_mol(&gl_mol, mol);
    if (mold_success(result)) {
        mold_set_backbone_secondary_structure_data(gl_mol, secondary_structure, count, offset);
    }
    return result;
}

mold_result mold_set_backbone_sequence_data(mold_molecule mol, const mold_range segment_range[], uint32_t count, uint32_t offset) {
    mold_gl_mol* gl_mol;
    mold_result result = mold_get_gl_mol(&gl_mol, mol);
    if (mold_success(result)) {
        mold_set_backbone_sequence_dat(gl_mol, segment_range, count, offset);
    }
    return result;
}