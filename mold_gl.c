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

typedef struct mold_gl_buffer_t {
    GLuint id;
} mold_gl_buffer_t;

typedef struct mold_gl_aabb_t {
    float min[3];
    float max[3];
} mold_gl_aabb_t;

typedef struct mold_gl_mol_t {
    struct {
        uint32_t count;
        mold_gl_buffer_t position; 
        mold_gl_buffer_t radius;
    } atom;

    struct {
        uint32_t count;
        mold_gl_buffer_t atom_range;
        mold_gl_buffer_t aabb;
    } residue;

    struct {
        uint32_t count;
        mold_gl_buffer_t data;
    } covalent_bond;

    struct {
        struct {
            uint32_t count;
            mold_gl_buffer_t atom_index;
            mold_gl_buffer_t secondary_structure;
        } segment;

        struct {
            uint32_t count;
            mold_gl_buffer_t segment_range;
        } sequence;
    } backbone;
} mold_gl_mol_t;

typedef struct mold_gl_rep_t {
    mold_rep_type_t type;
    mold_rep_args_t args;

    struct {
        uint32_t count;
        GLuint data;
    } color;
} mold_gl_rep_t;

typedef uint32_t mold_gl_version_t;
typedef enum mold_gl_version {
    MOLD_UNKNOWN = 0,
    MOLD_GL_330,
    MOLD_GL_430
} mold_gl_version;

typedef struct mold_gl_context_t {
    mold_gl_version_t version;

    struct {
        uint32_t count;
        mold_gl_mol_t data[MOLD_MOL_SLOTS];
    } molecule;

    struct {
        uint32_t count;
        mold_gl_rep_t data[MOLD_REP_SLOTS];
    };

    GLuint vao;
    GLuint ubo;
    GLuint fbo;
} mold_gl_context_t;

static mold_result_t mold_gl_init_common_resources(mold_gl_context_t* ctx) {
    glGenVertexArrays(1, &ctx->vao);
    glGenBuffers(1, &ctx->ubo);
    glGenFramebuffers(1, &ctx->fbo);

    return MOLD_SUCCESS;
}

static mold_result_t mold_gl_free_common_resources(mold_gl_context_t* ctx) {
    glDeleteVertexArrays(1, &ctx->vao);
    glDeleteBuffers(1, &ctx->ubo);
    glDeleteFramebuffers(1, &ctx->fbo);

    return MOLD_SUCCESS;
}

static mold_result_t mold_gl_setup_fbo(const mold_gl_desc_t* desc) {

    return MOLD_SUCCESS;
}

static mold_result_t mold_gl_set_atom_position_data(mold_gl_mol_t* mol, const mold_position_t pos_xyz[], uint32_t count, uint32_t offset) {

}

static mold_result_t mold_gl_set_atom_position_data_soa(mold_gl_mol_t* mol, const float pos_x[], const float pos_y[], const float pos_z[], uint32_t count, uint32_t offset) {

}

static mold_result_t mold_gl_set_atom_radius_data(mold_gl_mol_t* mol, const float radius[], uint32_t count, uint32_t offset) {

}

static mold_result_t mold_gl_set_residue_atom_range_data(mold_gl_mol_t* mol, const mold_atom_range_t atom_range[], uint32_t count, uint32_t offset) {

}

static mold_result_t mold_gl_set_covalent_bond_data(mold_gl_mol_t* mol, const mold_bond_t bond[], uint32_t count, uint32_t offset) {

}

static mold_result_t mold_gl_set_backbone_segment_index_data(mold_gl_mol_t* mol, const mold_backbone_segment_indices_t indices[], uint32_t count, uint32_t offset) {

}

static mold_result_t mold_gl_set_backbone_segment_secondary_structure_data(mold_gl_mol_t* mol, const mold_secondary_structure_t secondary_structure[], uint32_t count, uint32_t offset) {

}

static mold_result_t mold_gl_set_backbone_sequence_data(mold_gl_mol_t* mol, const mold_backbone_sequence_t sequence[], uint32_t count, uint32_t offset) {

}

static mold_result_t mold_gl_init_mol_resources(mold_gl_mol_t* mol, const mold_mol_desc_t* desc) {
    glGenBuffers(1, &mol->atom.position.id);
    glGenBuffers(1, &mol->atom.radius.id);
    mol->atom.count = desc->atom.count;
    if (mol->atom.count > 0) {
        glBindBuffer(GL_ARRAY_BUFFER, mol->atom.position.id);
        glBufferData(GL_ARRAY_BUFFER, mol->atom.count * sizeof(mold_position_t), NULL, GL_DYNAMIC_DRAW);
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
        glBufferData(GL_ARRAY_BUFFER, mol->residue.count * sizeof(mold_atom_range_t), NULL, GL_STATIC_DRAW);
        if (desc->residue.atom_range_data) {
            mold_gl_set_residue_atom_range_data(mol, desc->residue.atom_range_data, desc->residue.count, 0);
        }

        glBindBuffer(GL_ARRAY_BUFFER, mol->residue.aabb.id);
        glBufferData(GL_ARRAY_BUFFER, mol->residue.count * sizeof(mold_gl_aabb_t), NULL, GL_STATIC_DRAW);
    }

    glGenBuffers(1, &mol->covalent_bond.data.id);
    mol->covalent_bond.count = desc->covalent_bond.count;
    if (mol->covalent_bond.count > 0) {
        glBindBuffer(GL_ARRAY_BUFFER, mol->covalent_bond.data.id);
        glBufferData(GL_ARRAY_BUFFER, mol->covalent_bond.count * sizeof(mold_bond_t), NULL, GL_STATIC_DRAW);
        if (desc->covalent_bond.data) {
            mold_gl_set_covalent_bond_data(mol, desc->covalent_bond.data, desc->covalent_bond.count, 0);
        }
    }

    glGenBuffers(1, &mol->backbone.segment.atom_index.id);
    glGenBuffers(1, &mol->backbone.segment.secondary_structure.id);
    mol->backbone.segment.count = desc->backbone.segment.count;
    if (mol->backbone.segment.count) {
        glBindBuffer(GL_ARRAY_BUFFER, mol->backbone.segment.atom_index.id);
        glBufferData(GL_ARRAY_BUFFER, mol->backbone.segment.count * sizeof(mold_backbone_segment_indices_t), NULL, GL_STATIC_DRAW);
        if (desc->backbone.segment.index_data) {
            mold_gl_set_backbone_segment_index_data(mol, desc->backbone.segment.index_data, desc->backbone.segment.count, 0);
        }

        glBindBuffer(GL_ARRAY_BUFFER, mol->backbone.segment.secondary_structure.id);
        glBufferData(GL_ARRAY_BUFFER, mol->backbone.segment.count * sizeof(mold_secondary_structure_t), NULL, GL_STATIC_DRAW);
        if (desc->backbone.segment.index_data) {
            mold_gl_set_backbone_segment_secondary_structure_data(mol, desc->backbone.segment.secondary_structure_data, desc->backbone.segment.count, 0);
        }
    }

    glGenBuffers(1, &mol->backbone.sequence.segment_range.id);
    mol->backbone.sequence.count = desc->backbone.sequence.count;
    if (mol->backbone.sequence.count > 0) {
        glBindBuffer(GL_ARRAY_BUFFER, mol->backbone.sequence.segment_range.id);
        glBufferData(GL_ARRAY_BUFFER, mol->backbone.sequence.count * sizeof(mold_backbone_sequence_t), NULL, GL_STATIC_DRAW);
        if (desc->backbone.sequence.data) {
            mold_gl_set_backbone_sequence_data(mol, desc->backbone.sequence.data, desc->backbone.sequence.count, 0);
        }
    }

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    return MOLD_SUCCESS;
}

static mold_result_t mold_gl_free_mol_resources(mold_gl_mol_t* mol) {
    glDeleteBuffers(1, &mol->atom.position.id);
    glDeleteBuffers(1, &mol->atom.radius.id);
    glDeleteBuffers(1, &mol->residue.atom_range.id);
    glDeleteBuffers(1, &mol->residue.aabb.id);
    glDeleteBuffers(1, &mol->covalent_bond.data.id);
    glDeleteBuffers(1, &mol->backbone.segment.atom_index.id);
    glDeleteBuffers(1, &mol->backbone.segment.secondary_structure.id);
    glDeleteBuffers(1, &mol->backbone.sequence.segment_range.id);
}

static mold_result_t mold_gl_init_rep_resources(mold_gl_rep_t* rep, const mold_rep_desc_t* desc) {

}

static mold_result_t mold_gl_free_rep_resources(mold_gl_rep_t* rep) {

}

static mold_result_t mold_context_init(const mold_gl_desc_t* desc) {
    GLint major, minor;
    glGetInteger(GL_MAJOR_VERSION, &major);
    glGetInteger(GL_MINOR_VERSION, &minor);

    mold_gl_version_t version = MOLD_UNKNOWN;
    if (major >= 4 && minor >= 3) {
        version = MOLD_GL_430;
    }
    else if (major >= 3 && minor >= 3) {
        version = MOLD_GL_330;
    }
    else {
        return MOLD_GL_VERSION_NOT_SUPPORTED;
    }
    ctx.version = version;

    mold_result_t result;
    result = mold_gl_init_common_resources();
    if (result != MOLD_SUCCESS) {
        return result;
    }

    result = mold_gl_setup_fbo(desc);
    if (result != MOLD_SUCCESS) {
        return result;
    }

    return result;
}

static mold_result_t mold_gl_free_internal(mold_gl_context_t* ctx) {
    for (uint32_t i = 0; i < ctx->molecule.count; i++) {

    }
    return mold_gl_free_common_resources(ctx);
}

static mold_result_t mold_init_molecule(mold_mol_t* mol, const mold_mol_desc_t* desc) {
    if (ctx.molecule.count == MOLD_MOL_SLOTS) {
        return MOLD_OUT_OF_MOLECULE_SLOTS;
    }

    mol->id = 0xB00BF00D + ctx->molecule.count++;

    return MOLD_SUCCESS;
}

static mold_result_t mold_gl_free_molecule(mold_gl_context_t* ctx, mold_mol_t* mol) {

}