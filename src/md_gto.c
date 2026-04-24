#include <md_gto.h>

#include <core/md_platform.h>

#include <core/md_log.h>
#include <core/md_simd.h>
#include <core/md_allocator.h>
#include <core/md_str.h>

#include <stdbool.h>
#include <float.h>

static inline void world_to_model_matrix(float out_mat[4][4], const md_grid_t* grid) {
    // There is no scaling applied in this transformation, only rotation and translation
    out_mat[0][0] = grid->orientation.elem[0][0];
    out_mat[0][1] = grid->orientation.elem[1][0];
    out_mat[0][2] = grid->orientation.elem[2][0];
    out_mat[0][3] = 0.0f;
    out_mat[1][0] = grid->orientation.elem[0][1];
    out_mat[1][1] = grid->orientation.elem[1][1];
    out_mat[1][2] = grid->orientation.elem[2][1];
    out_mat[1][3] = 0.0f;
    out_mat[2][0] = grid->orientation.elem[0][2];
    out_mat[2][1] = grid->orientation.elem[1][2];
    out_mat[2][2] = grid->orientation.elem[2][2];
    out_mat[2][3] = 0.0f;
    out_mat[3][0] = -grid->orientation.elem[0][0] * grid->origin.elem[0] - grid->orientation.elem[0][1] * grid->origin.elem[1] - grid->orientation.elem[0][2] * grid->origin.elem[2],
    out_mat[3][1] = -grid->orientation.elem[1][0] * grid->origin.elem[0] - grid->orientation.elem[1][1] * grid->origin.elem[1] - grid->orientation.elem[1][2] * grid->origin.elem[2],
    out_mat[3][2] = -grid->orientation.elem[2][0] * grid->origin.elem[0] - grid->orientation.elem[2][1] * grid->origin.elem[1] - grid->orientation.elem[2][2] * grid->origin.elem[2],
    out_mat[3][3] = 1.0f;
}

static inline void index_to_world_matrix(float out_mat[4][4], const md_grid_t* grid) {
    out_mat[0][0] = grid->orientation.elem[0][0] * grid->spacing.elem[0];
    out_mat[0][1] = grid->orientation.elem[0][1] * grid->spacing.elem[0];
    out_mat[0][2] = grid->orientation.elem[0][2] * grid->spacing.elem[0];
    out_mat[0][3] = 0.0f;
    out_mat[1][0] = grid->orientation.elem[1][0] * grid->spacing.elem[1];
    out_mat[1][1] = grid->orientation.elem[1][1] * grid->spacing.elem[1];
    out_mat[1][2] = grid->orientation.elem[1][2] * grid->spacing.elem[1];
    out_mat[1][3] = 0.0f;
    out_mat[2][0] = grid->orientation.elem[2][0] * grid->spacing.elem[2];
    out_mat[2][1] = grid->orientation.elem[2][1] * grid->spacing.elem[2];
    out_mat[2][2] = grid->orientation.elem[2][2] * grid->spacing.elem[2];
    out_mat[2][3] = 0.0f;
    out_mat[3][0] = grid->origin.elem[0];
    out_mat[3][1] = grid->origin.elem[1];
    out_mat[3][2] = grid->origin.elem[2];
    out_mat[3][3] = 1.0f;

    // Incorporate a half voxel offset to move to voxel centers
    out_mat[3][0] += 0.5f * (out_mat[0][0] + out_mat[1][0] + out_mat[2][0]);
    out_mat[3][1] += 0.5f * (out_mat[0][1] + out_mat[1][1] + out_mat[2][1]);
    out_mat[3][2] += 0.5f * (out_mat[0][2] + out_mat[1][2] + out_mat[2][2]);
}


// ---------------------------------------------------------------------------
// Spherical→Cartesian expansion tables
// Converts md_gto_basis_t (radial shells, pure radial coefficients) into the
// flat Cartesian SoA layout (md_gto_data_t / GPU buffer) used by evaluators.
// Coupling coefficients and Cartesian (l,m,n) tables ported from VeloxChem.
// ---------------------------------------------------------------------------

typedef uint8_t gto_lmn_t[3];

static const gto_lmn_t gto_S_lmn[1]  = {{0,0,0}};
static const gto_lmn_t gto_P_lmn[3]  = {{1,0,0},{0,1,0},{0,0,1}};
static const gto_lmn_t gto_D_lmn[6]  = {{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}};
static const gto_lmn_t gto_F_lmn[10] = {{3,0,0},{2,1,0},{2,0,1},{1,2,0},{1,1,1},{1,0,2},{0,3,0},{0,2,1},{0,1,2},{0,0,3}};
static const gto_lmn_t gto_G_lmn[15] = {{4,0,0},{3,1,0},{3,0,1},{2,2,0},{2,1,1},{2,0,2},{1,3,0},{1,2,1},{1,1,2},{1,0,3},{0,4,0},{0,3,1},{0,2,2},{0,1,3},{0,0,4}};

#define GTO_d3  3.464101615137754587
#define GTO_f5  1.581138830084189666
#define GTO_f15 7.745966692414833770
#define GTO_f3  1.224744871391589049
#define GTO_g35 (4.0 * 5.916079783099616042)
#define GTO_g17 (4.0 * 4.183300132670377739)
#define GTO_g5  (4.0 * 2.236067977499789696)
#define GTO_g2  (4.0 * 1.581138830084189666)

static const double  gto_S_factors[] = {1.0};
static const uint8_t gto_S_indices[] = {0};
static const uint8_t gto_S_num_fac[] = {1};

static const double  gto_P_factors[] = {1.0, 1.0, 1.0};
static const uint8_t gto_P_indices[] = {1, 2, 0};
static const uint8_t gto_P_offsets[] = {0, 1, 2};
static const uint8_t gto_P_num_fac[] = {1, 1, 1};

static const double  gto_D_factors[] = {GTO_d3, GTO_d3, -1.0, -1.0, 2.0, GTO_d3, 0.5*GTO_d3, -0.5*GTO_d3};
static const uint8_t gto_D_indices[] = {1, 4, 0, 3, 5, 2, 0, 3};
static const uint8_t gto_D_offsets[] = {0, 1, 2, 5, 6};
static const uint8_t gto_D_num_fac[] = {1, 1, 3, 1, 2};

static const double  gto_F_factors[] = {3.0*GTO_f5, -GTO_f5, GTO_f15, 4.0*GTO_f3, -GTO_f3, -GTO_f3, 2.0, -3.0, -3.0, 4.0*GTO_f3, -GTO_f3, -GTO_f3, 0.5*GTO_f15, -0.5*GTO_f15, GTO_f5, -3.0*GTO_f5};
static const uint8_t gto_F_indices[] = {1, 6, 4, 8, 1, 6, 9, 2, 7, 5, 0, 3, 2, 7, 0, 3};
static const uint8_t gto_F_offsets[] = {0, 2, 3, 6, 9, 12, 14};
static const uint8_t gto_F_num_fac[] = {2, 1, 3, 3, 3, 2, 2};

static const double gto_G_factors[] = {
    GTO_g35, -GTO_g35, 3.0*GTO_g17, -GTO_g17, 6.0*GTO_g5, -GTO_g5, -GTO_g5, 4.0*GTO_g2, -3.0*GTO_g2, -3.0*GTO_g2,
    8.0, 3.0, 3.0, 6.0, -24.0, -24.0, 4.0*GTO_g2, -3.0*GTO_g2, -3.0*GTO_g2, 3.0*GTO_g5,
    -3.0*GTO_g5, -0.5*GTO_g5, 0.5*GTO_g5, GTO_g17, -3.0*GTO_g17, 0.25*GTO_g35, 0.25*GTO_g35, -1.50*GTO_g35};
static const uint8_t gto_G_indices[] = {1, 6, 4, 11, 8, 1, 6, 13, 4, 11, 14, 0, 10, 3, 5, 12, 9, 2, 7, 5, 12, 0, 10, 2, 7, 0, 10, 3};
static const uint8_t gto_G_offsets[] = {0, 2, 4, 7, 10, 16, 19, 23, 25};
static const uint8_t gto_G_num_fac[] = {2, 2, 3, 3, 6, 3, 4, 2, 3};

#undef GTO_d3
#undef GTO_f5
#undef GTO_f15
#undef GTO_f3
#undef GTO_g35
#undef GTO_g17
#undef GTO_g5
#undef GTO_g2

static inline int gto_sph_num_factors(int l, int isph) {
    switch (l) {
    case 0: return gto_S_num_fac[isph];
    case 1: return gto_P_num_fac[isph];
    case 2: return gto_D_num_fac[isph];
    case 3: return gto_F_num_fac[isph];
    case 4: return gto_G_num_fac[isph];
    default: ASSERT(false); return 0;
    }
}

static inline const double* gto_sph_factors(int l, int isph) {
    switch (l) {
    case 0: return gto_S_factors;
    case 1: return gto_P_factors + gto_P_offsets[isph];
    case 2: return gto_D_factors + gto_D_offsets[isph];
    case 3: return gto_F_factors + gto_F_offsets[isph];
    case 4: return gto_G_factors + gto_G_offsets[isph];
    default: ASSERT(false); return NULL;
    }
}

static inline const uint8_t* gto_sph_indices(int l, int isph) {
    switch (l) {
    case 0: return gto_S_indices;
    case 1: return gto_P_indices + gto_P_offsets[isph];
    case 2: return gto_D_indices + gto_D_offsets[isph];
    case 3: return gto_F_indices + gto_F_offsets[isph];
    case 4: return gto_G_indices + gto_G_offsets[isph];
    default: ASSERT(false); return NULL;
    }
}

static inline const gto_lmn_t* gto_cart_lmn(int l) {
    switch (l) {
    case 0: return gto_S_lmn;
    case 1: return gto_P_lmn;
    case 2: return gto_D_lmn;
    case 3: return gto_F_lmn;
    case 4: return gto_G_lmn;
    default: ASSERT(false); return NULL;
    }
}

// Count the number of CGTOs and PGTOs that will result from expanding a basis.
static void gto_basis_count(uint32_t* out_num_cgtos, uint32_t* out_num_pgtos,
    const md_gto_basis_t* basis)
{
    uint32_t nc = 0, np = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        int l     = (int)basis->shells[si].l;
        int nsph  = 2 * l + 1;
        int nprims = (int)basis->shells[si].num_primitives;
        nc += (uint32_t)nsph;
        for (int isph = 0; isph < nsph; isph++) {
            np += (uint32_t)(gto_sph_num_factors(l, isph) * nprims);
        }
    }
    *out_num_cgtos = nc;
    *out_num_pgtos = np;
}

// Expand md_gto_basis_t into the flat Cartesian SoA arrays expected by the GPU shader.
// Spherical-to-Cartesian coupling factors (fcarts) are baked into pgto_coeff here.
// pgto_radius is set to FLT_MAX (no pre-culling; the shader handles distance tests).
// All output arrays must be pre-allocated to num_cgtos / num_pgtos entries respectively.
static void gto_expand_basis(
    vec4_t* cgto_xyzr, uint32_t* cgto_offset,
    float* pgto_alpha, float* pgto_coeff, float* pgto_radius, uint32_t* pgto_ijkl,
    const md_gto_basis_t* basis, const float* atom_xyz)
{
    uint32_t ci = 0, pi = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        const md_gto_shell_t* shell = &basis->shells[si];
        int l      = (int)shell->l;
        int nsph   = 2 * l + 1;
        int nprims = (int)shell->num_primitives;
        uint32_t   prim_base = shell->primitive_offset;
        const gto_lmn_t* lmn = gto_cart_lmn(l);
        float ax = atom_xyz[shell->atom_idx * 3 + 0];
        float ay = atom_xyz[shell->atom_idx * 3 + 1];
        float az = atom_xyz[shell->atom_idx * 3 + 2];

        for (int isph = 0; isph < nsph; isph++) {
            int           ncomp   = gto_sph_num_factors(l, isph);
            const double* fcarts  = gto_sph_factors(l, isph);
            const uint8_t* sidx   = gto_sph_indices(l, isph);

            int lx[8], ly[8], lz[8];
            for (int ic = 0; ic < ncomp; ic++) {
                lx[ic] = lmn[sidx[ic]][0];
                ly[ic] = lmn[sidx[ic]][1];
                lz[ic] = lmn[sidx[ic]][2];
            }

            cgto_xyzr[ci]   = (vec4_t){ax, ay, az, FLT_MAX};
            cgto_offset[ci] = pi;

            for (int ip = 0; ip < nprims; ip++) {
                float alpha = basis->alpha[prim_base + ip];
                float coef1 = basis->coeff[prim_base + ip];
                for (int ic = 0; ic < ncomp; ic++) {
                    pgto_alpha[pi]  = alpha;
                    pgto_coeff[pi]  = (float)(coef1 * fcarts[ic]);
                    pgto_radius[pi] = FLT_MAX;
                    pgto_ijkl[pi]   = md_gto_pack_ijkl(lx[ic], ly[ic], lz[ic], l);
                    pi++;
                }
            }
            ci++;
        }
    }
    cgto_offset[ci] = pi;  // sentinel
}

static size_t density_matrix_upper_tri_size(size_t n) {
    return n * (n + 1) / 2;
}

// Convert a full N×N row-major double density matrix to compact upper-triangular float.
static void density_matrix_upper_tri_extract_float(float* out, const double* dm, size_t n) {
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            out[k++] = (float)dm[i * n + j];
        }
    }
}

size_t md_gto_pgto_count(const md_gto_basis_t* basis) {
    uint32_t nc, np;
    gto_basis_count(&nc, &np, basis);
    return np;
}

size_t md_gto_expand_with_mo(md_gto_t* out, const md_gto_basis_t* basis,
    const float* atom_xyz, const double* mo_coeffs, double cutoff)
{
    ASSERT(out);
    ASSERT(basis);
    ASSERT(atom_xyz);
    ASSERT(mo_coeffs);

    size_t num_gtos = 0;
    uint32_t cgto_idx = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        const md_gto_shell_t* shell = &basis->shells[si];
        int l      = (int)shell->l;
        int nsph   = 2 * l + 1;
        int nprims = (int)shell->num_primitives;
        uint32_t prim_base = shell->primitive_offset;
        const gto_lmn_t* lmn = gto_cart_lmn(l);
        float ax = atom_xyz[shell->atom_idx * 3 + 0];
        float ay = atom_xyz[shell->atom_idx * 3 + 1];
        float az = atom_xyz[shell->atom_idx * 3 + 2];

        for (int isph = 0; isph < nsph; isph++) {
            double mo_coeff = mo_coeffs[cgto_idx++];
            int           ncomp  = gto_sph_num_factors(l, isph);
            const double* fcarts = gto_sph_factors(l, isph);
            const uint8_t* sidx  = gto_sph_indices(l, isph);

            for (int ip = 0; ip < nprims; ip++) {
                float alpha = basis->alpha[prim_base + ip];
                float coef1 = basis->coeff[prim_base + ip];
                for (int ic = 0; ic < ncomp; ic++) {
                    out[num_gtos++] = (md_gto_t){
                        .x      = ax,
                        .y      = ay,
                        .z      = az,
                        .coeff  = (float)(coef1 * fcarts[ic] * mo_coeff),
                        .alpha  = alpha,
                        .cutoff = FLT_MAX,
                        .i      = lmn[sidx[ic]][0],
                        .j      = lmn[sidx[ic]][1],
                        .k      = lmn[sidx[ic]][2],
                        .l      = (uint8_t)l,
                    };
                }
            }
        }
    }

    if (cutoff > 0.0) {
        num_gtos = md_gto_cutoff_compute_and_filter(out, num_gtos, cutoff);
    }
    return num_gtos;
}

#if !MD_PLATFORM_OSX

#include <core/md_gl_util.h>
#include <gto_shaders.inl>
#include <GL/gl3w.h>

// This should be kept in sync with the define present in segment_and_attribute_to_group.comp
#define QUANTIZATION_SCALE_FACTOR 1.0e6

static GLuint get_gto_program(void) {
    static GLuint program = 0;
    if (!program) {
        GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
        if (md_gl_shader_compile(shader, (str_t){(const char*)eval_gto_comp, eval_gto_comp_size}, 0, 0)) {
            GLuint prog = glCreateProgram();
            if (md_gl_program_attach_and_link(prog, &shader, 1)) {
                program = prog;
            }
        }
        glDeleteShader(shader);
    }
    return program;
}

static GLuint get_gto_density_program(void) {
    static GLuint program = 0;
    if (!program) {
        GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
        if (md_gl_shader_compile(shader, (str_t){(const char*)eval_gto_density_comp, eval_gto_density_comp_size}, 0, 0)) {
            GLuint prog = glCreateProgram();
            if (md_gl_program_attach_and_link(prog, &shader, 1)) {
                program = prog;
            }
        }
        glDeleteShader(shader);
    }
    return program;
}

static GLuint get_gto_density_grad_program(void) {
    static GLuint program = 0;
    if (!program) {
        GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
        if (md_gl_shader_compile(shader, (str_t) { (const char*)eval_gto_density_grad_comp, eval_gto_density_grad_comp_size }, 0, 0)) {
            GLuint prog = glCreateProgram();
            if (md_gl_program_attach_and_link(prog, &shader, 1)) {
                program = prog;
            }
        }
        glDeleteShader(shader);
    }
    return program;
}


static GLuint get_buffer(size_t size) {
    GLuint id = 0;
    glCreateBuffers(1, &id);
    glBindBuffer(GL_ARRAY_BUFFER, id);
    glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    return id;
}

static void free_buffer(GLuint id) {
    if (glIsBuffer(id)) {
        glDeleteBuffers(1, &id);
    }
}

static void gto_grid_evaluate_orb_GPU(uint32_t vol_tex, const md_grid_t* grid, const md_orbital_data_t* orb, md_gto_eval_mode_t mode, GLuint program) {
    ASSERT(grid);
    ASSERT(orb);

    if (orb->num_gtos == 0 || orb->num_orbs == 0) {
        return;
    }

    md_gl_debug_push("EVAL ORBS");

    if (!glIsTexture(vol_tex)) {
        MD_LOG_ERROR("Invalid volume texture handle");
        return;
    }

    GLenum format = 0;
    if (glGetTextureLevelParameteriv) {
        glGetTextureLevelParameteriv(vol_tex,   0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&format);
    } else {
        glBindTexture(GL_TEXTURE_3D, vol_tex);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&format);
        glBindTexture(GL_TEXTURE_3D, 0);
    }

    switch (format) {
    case GL_R16F:
    case GL_R32F:
        break;
    default:
        // Not good
        MD_LOG_ERROR("Unrecognized internal format of supplied volume texture");
        goto done;
    }

    GLintptr   ssbo_gto_offset = 0;
    GLsizeiptr ssbo_gto_size   = sizeof(md_gto_t) * orb->num_gtos;

    GLintptr   ssbo_orb_offset = ALIGN_TO(ssbo_gto_offset + ssbo_gto_size, 256);
    GLsizeiptr ssbo_orb_size   = sizeof(uint32_t) * (orb->num_orbs + 1);

    GLintptr   ssbo_scl_offset = ALIGN_TO(ssbo_orb_offset + ssbo_orb_size, 256);
    GLsizeiptr ssbo_scl_size   = sizeof(float) * (orb->num_orbs);

    size_t total_size = ALIGN_TO(ssbo_scl_offset + ssbo_scl_size, 256);
    GLuint ssbo = get_buffer(total_size);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_gto_offset, ssbo_gto_size, orb->gtos);

    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_orb_offset, ssbo_orb_size, orb->orb_offsets);
    // Fill last portion of buffer with point indices
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_scl_offset, ssbo_scl_size, orb->orb_scaling);

    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 0, ssbo, ssbo_gto_offset, ssbo_gto_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 1, ssbo, ssbo_orb_offset, ssbo_orb_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 2, ssbo, ssbo_scl_offset, ssbo_scl_size);

    glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);

    glUseProgram(program);

    float world_to_model[4][4];
    float index_to_world[4][4];

    world_to_model_matrix(world_to_model, grid);
    index_to_world_matrix(index_to_world, grid);

    glUniformMatrix4fv(0, 1, GL_FALSE, (const float*)world_to_model);
    glUniformMatrix4fv(1, 1, GL_FALSE, (const float*)index_to_world);
    glUniform3fv(2, 1, grid->spacing.elem);
    glUniform1ui(3, (GLuint)orb->num_orbs);
    glUniform1i(4, (GLint)mode);

    glBindImageTexture(0, vol_tex, 0, GL_TRUE, 0, GL_WRITE_ONLY, format);

    int num_groups[3] = {
        DIV_UP(grid->dim[0], 8),
        DIV_UP(grid->dim[1], 8),
        DIV_UP(grid->dim[2], 8),
    };

    glDispatchCompute(num_groups[0], num_groups[1], num_groups[2]);

    glUseProgram(0);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT | GL_PIXEL_BUFFER_BARRIER_BIT);

    free_buffer(ssbo);
done:
    md_gl_debug_pop();
}

void md_gto_grid_evaluate_orb_GPU(uint32_t vol_tex, const md_grid_t* vol_grid, const md_orbital_data_t* orb, md_gto_eval_mode_t mode) {
    ASSERT(vol_grid);
    ASSERT(orb);

    GLuint program = get_gto_program();
    gto_grid_evaluate_orb_GPU(vol_tex, vol_grid, orb, mode, program);
}

void md_gto_grid_evaluate_GPU(uint32_t vol_tex, const md_grid_t* vol_grid, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    uint32_t orb_offsets[2] = {0, (uint32_t)num_gtos};
    float orb_scaling[1] = {1.0f};

    md_orbital_data_t orb = {
        .num_gtos = num_gtos,
        .gtos = (md_gto_t*)gtos,
        .num_orbs = 1,
        .orb_offsets = orb_offsets,
        .orb_scaling = orb_scaling,
    };

    md_gto_grid_evaluate_orb_GPU(vol_tex, vol_grid, &orb, mode);
}

void md_gto_grid_evaluate_matrix_GPU(uint32_t vol_tex, const md_grid_t* grid, const md_gto_data_t* gto_data, const float* upper_triangular_matrix_data, size_t matrix_dim, bool include_gradients) {
    ASSERT(grid);
    ASSERT(gto_data);
    ASSERT(upper_triangular_matrix_data);

    md_gl_debug_push("EVAL DENSITY");

    if (!glIsTexture(vol_tex)) {
        MD_LOG_ERROR("Invalid volume texture handle");
        return;
    }

    GLint format = 0;
    glBindTexture(GL_TEXTURE_3D, vol_tex);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_INTERNAL_FORMAT, &format);
    glBindTexture(GL_TEXTURE_3D, 0);

    GLuint program = 0;
    if (include_gradients) {
        switch (format) {
        case GL_RGBA16F:
        case GL_RGBA32F:
            break;
        default:
            // Not good
            MD_LOG_ERROR("Unrecognized internal format of supplied volume texture");
            goto done;
        }
        program = get_gto_density_grad_program();
    }
    else {
        switch (format) {
        case GL_R16F:
        case GL_R32F:
            break;
        default:
            // Not good
            MD_LOG_ERROR("Unrecognized internal format of supplied volume texture");
            goto done;
        }
        program = get_gto_density_program();
    }

    if (!program) {
        MD_LOG_ERROR("Program not found?!");
        goto done;
    }

    size_t matrix_data_len = (matrix_dim * (matrix_dim + 1)) / 2;
    const float* matrix_data = upper_triangular_matrix_data;

    if (matrix_dim != gto_data->num_cgtos) {
        MD_LOG_ERROR("Matrix dim, cgto mismatch");
        goto done;
    }

    typedef struct {
        mat4_t world_to_model;
        mat4_t index_to_world;
        vec4_t step;
        uint32_t D_matrix_dim;
        uint32_t _pad[3];
    } uniform_block_t;

    uniform_block_t ub_data = {0};
    world_to_model_matrix(ub_data.world_to_model.elem, grid);
    index_to_world_matrix(ub_data.index_to_world.elem, grid);
    ub_data.step = vec4_from_vec3(grid->spacing, 0);
    ub_data.D_matrix_dim = (uint32_t)matrix_dim;

    GLintptr   ssbo_cgto_xyzr_base      = 0;
    GLsizeiptr ssbo_cgto_xyzr_size      = sizeof(vec4_t) * gto_data->num_cgtos;

    GLintptr   ssbo_cgto_offset_base    = ALIGN_TO(ssbo_cgto_xyzr_base + ssbo_cgto_xyzr_size, 256);
    GLsizeiptr ssbo_cgto_offset_size    = sizeof(uint32_t) * (gto_data->num_cgtos + 1);

    GLintptr   ssbo_pgto_coeff_base     = ALIGN_TO(ssbo_cgto_offset_base + ssbo_cgto_offset_size, 256);
    GLsizeiptr ssbo_pgto_coeff_size     = sizeof(float) * (gto_data->num_pgtos);

    GLintptr   ssbo_pgto_alpha_base     = ALIGN_TO(ssbo_pgto_coeff_base + ssbo_pgto_coeff_size, 256);
    GLsizeiptr ssbo_pgto_alpha_size     = sizeof(float) * (gto_data->num_pgtos); 

    GLintptr   ssbo_pgto_radius_base    = ALIGN_TO(ssbo_pgto_alpha_base + ssbo_pgto_alpha_size, 256);
    GLsizeiptr ssbo_pgto_radius_size    = sizeof(float) * (gto_data->num_pgtos); 

    GLintptr   ssbo_pgto_ijkl_base      = ALIGN_TO(ssbo_pgto_radius_base + ssbo_pgto_radius_size, 256);
    GLsizeiptr ssbo_pgto_ijkl_size      = sizeof(uint32_t) * (gto_data->num_pgtos);

    GLintptr   ssbo_matrix_base         = ALIGN_TO(ssbo_pgto_ijkl_base + ssbo_pgto_ijkl_size, 256);
    GLsizeiptr ssbo_matrix_size         = sizeof(float) * matrix_data_len;

    GLintptr   ubo_base                 = ALIGN_TO(ssbo_matrix_base + ssbo_matrix_size, 256);
    GLsizeiptr ubo_size                 = sizeof(uniform_block_t);

    size_t total_size = ALIGN_TO(ubo_base + ubo_size, 256);
    GLuint buf = get_buffer(total_size);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, buf);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_cgto_xyzr_base,      ssbo_cgto_xyzr_size,    gto_data->cgto_xyzr);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_cgto_offset_base,    ssbo_cgto_offset_size,  gto_data->cgto_offset);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_pgto_coeff_base,     ssbo_pgto_coeff_size,   gto_data->pgto_coeff);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_pgto_alpha_base,     ssbo_pgto_alpha_size,   gto_data->pgto_alpha);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_pgto_radius_base,    ssbo_pgto_radius_size,  gto_data->pgto_radius);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_pgto_ijkl_base,      ssbo_pgto_ijkl_size,    gto_data->pgto_ijkl);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_matrix_base,         ssbo_matrix_size,       matrix_data);

    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 0, buf, ssbo_cgto_xyzr_base,       ssbo_cgto_xyzr_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 1, buf, ssbo_cgto_offset_base,     ssbo_cgto_offset_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 2, buf, ssbo_pgto_coeff_base,      ssbo_pgto_coeff_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 3, buf, ssbo_pgto_alpha_base,      ssbo_pgto_alpha_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 4, buf, ssbo_pgto_radius_base,     ssbo_pgto_radius_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 5, buf, ssbo_pgto_ijkl_base,       ssbo_pgto_ijkl_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 6, buf, ssbo_matrix_base,          ssbo_matrix_size);

    glBindBuffer(GL_UNIFORM_BUFFER, buf);
    glBufferSubData(GL_UNIFORM_BUFFER, ubo_base, ubo_size, &ub_data);
    glBindBufferRange(GL_UNIFORM_BUFFER, 0, buf, ubo_base, ubo_size);

    glUseProgram(program);
    GLuint block = glGetUniformBlockIndex(program, "UniformBlock");
    glUniformBlockBinding(program, block, 0);

    glBindImageTexture(0, vol_tex, 0, GL_TRUE, 0, GL_WRITE_ONLY, format);

    GLuint query;
    glGenQueries(1, &query);

    // Start timing
    glBeginQuery(GL_TIME_ELAPSED, query);

    int num_groups[3] = {
        DIV_UP(grid->dim[0], 8),
        DIV_UP(grid->dim[1], 8),
        DIV_UP(grid->dim[2], 8),
    };
    glDispatchCompute(num_groups[0], num_groups[1], num_groups[2]);

    // End timing
    glEndQuery(GL_TIME_ELAPSED);

    // Retrieve the result (blocking until GPU finishes)
    GLuint64 elapsedTime = 0;
    glGetQueryObjectui64v(query, GL_QUERY_RESULT, &elapsedTime); // nanoseconds

	MD_LOG_DEBUG("GTO Density evaluation of [%i,%i,%i] GPU time: %.3f ms", grid->dim[0], grid->dim[1], grid->dim[2], elapsedTime / 1e6);

    glUseProgram(0);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT | GL_PIXEL_BUFFER_BARRIER_BIT);

done:
    md_gl_debug_pop();
}

void md_gto_grid_evaluate_density_GL(uint32_t vol_tex, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz,
    const double* density_matrix, bool include_gradients)
{
    ASSERT(grid);
    ASSERT(basis);
    ASSERT(atom_xyz);
    ASSERT(density_matrix);

    uint32_t num_cgtos, num_pgtos;
    gto_basis_count(&num_cgtos, &num_pgtos, basis);

    size_t temp_pos = md_temp_get_pos();
    vec4_t*   cgto_xyzr   = (vec4_t*)  md_temp_push(sizeof(vec4_t)   * num_cgtos);
    uint32_t* cgto_offset = (uint32_t*)md_temp_push(sizeof(uint32_t) * (num_cgtos + 1));
    float*    pgto_alpha  = (float*)   md_temp_push(sizeof(float)    * num_pgtos);
    float*    pgto_coeff  = (float*)   md_temp_push(sizeof(float)    * num_pgtos);
    float*    pgto_radius = (float*)   md_temp_push(sizeof(float)    * num_pgtos);
    uint32_t* pgto_ijkl   = (uint32_t*)md_temp_push(sizeof(uint32_t) * num_pgtos);
    size_t    tri_len     = ((size_t)num_cgtos * (num_cgtos + 1)) / 2;
    float*    upper_tri   = (float*)   md_temp_push(sizeof(float)    * tri_len);

    gto_expand_basis(cgto_xyzr, cgto_offset, pgto_alpha, pgto_coeff, pgto_radius, pgto_ijkl, basis, atom_xyz);
    density_matrix_to_upper_tri(upper_tri, density_matrix, num_cgtos);

    md_gto_data_t gto_data = {
        .num_cgtos   = num_cgtos,
        .cgto_xyzr   = cgto_xyzr,
        .cgto_offset = cgto_offset,
        .num_pgtos   = num_pgtos,
        .pgto_alpha  = pgto_alpha,
        .pgto_coeff  = pgto_coeff,
        .pgto_radius = pgto_radius,
        .pgto_ijkl   = pgto_ijkl,
    };

    md_gto_grid_evaluate_matrix_GPU(vol_tex, grid, &gto_data, upper_tri, num_cgtos, include_gradients);

    md_temp_set_pos_back(temp_pos);
}

#else

void md_gto_grid_evaluate_orb_GPU(uint32_t vol_tex, const md_grid_t* vol_grid, const md_orbital_data_t* orb, md_gto_eval_mode_t mode) {

}

void md_gto_segment_and_attribute_to_groups_GPU(float* out_group_values, size_t cap_groups, uint32_t vol_tex, const md_grid_t* grid, const float* point_xyzr, const uint32_t* point_group_idx, size_t num_points) {

}

void md_gto_grid_evaluate_GPU(uint32_t vol_tex, const md_grid_t* vol_grid, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {

}

void md_gto_grid_evaluate_matrix_GPU(uint32_t vol_tex, const md_grid_t* grid, const md_gto_data_t* gto_data, const float* matrix_data, size_t matrix_dim, bool output_grad_and_value) {

}

void md_gto_grid_evaluate_density_GL(uint32_t vol_tex, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz,
    const double* density_matrix, bool include_gradients) {
    (void)vol_tex; (void)grid; (void)basis; (void)atom_xyz; (void)density_matrix; (void)include_gradients;
}

#endif

#if MD_ENABLE_GPU

#include <core/md_gpu.h>
#include <gto_gpu_shaders.inl>

static md_gpu_device_t gto_cached_device = NULL;
static md_gpu_compute_pipeline_t gto_pip_density = NULL;

static md_gpu_compute_pipeline_t gto_ensure_density_pipeline(md_gpu_device_t device) {
    if (gto_cached_device != device) {
        if (gto_pip_density) { md_gpu_destroy_compute_pipeline(gto_pip_density); gto_pip_density = NULL; }
        gto_cached_device = device;
    }
    if (!gto_pip_density) {
        md_gpu_compute_pipeline_desc_t desc = {
            .shader = { .data = gto_eval_gto_density_start, .size = gto_eval_gto_density_size() },
            .threadgroup_size = { 8, 8, 8 },
        };
        gto_pip_density = md_gpu_create_compute_pipeline(device, &desc);
        if (!gto_pip_density) {
            MD_LOG_ERROR("Failed to create GTO density compute pipeline");
        }
    }
    return gto_pip_density;
}

// Sub-region layout for the packed GPU buffer.
// All offsets are 256-byte aligned to satisfy both Vulkan and Metal constraints.
typedef struct {
    size_t off_cgto_xyzr;   size_t sz_cgto_xyzr;
    size_t off_cgto_offset; size_t sz_cgto_offset;
    size_t off_pgto_coeff;  size_t sz_pgto_coeff;
    size_t off_pgto_alpha;  size_t sz_pgto_alpha;
    size_t off_pgto_radius; size_t sz_pgto_radius;
    size_t off_pgto_ijkl;   size_t sz_pgto_ijkl;
    size_t off_matrix;      size_t sz_matrix;
    size_t total_size;
} gto_buf_layout_t;

static gto_buf_layout_t gto_density_buf_compute_layout(uint32_t num_cgtos, uint32_t num_pgtos, uint32_t matrix_dim) {
    gto_buf_layout_t L;
    size_t matrix_data_len = ((size_t)matrix_dim * (matrix_dim + 1)) / 2;

    L.off_cgto_xyzr   = 0;
    L.sz_cgto_xyzr    = sizeof(float)    * 4 * num_cgtos;

    L.off_cgto_offset = ALIGN_TO(L.off_cgto_xyzr   + L.sz_cgto_xyzr,   256);
    L.sz_cgto_offset  = sizeof(uint32_t) * (num_cgtos + 1);

    L.off_pgto_coeff  = ALIGN_TO(L.off_cgto_offset  + L.sz_cgto_offset,  256);
    L.sz_pgto_coeff   = sizeof(float)    * num_pgtos;

    L.off_pgto_alpha  = ALIGN_TO(L.off_pgto_coeff   + L.sz_pgto_coeff,   256);
    L.sz_pgto_alpha   = sizeof(float)    * num_pgtos;

    L.off_pgto_radius = ALIGN_TO(L.off_pgto_alpha   + L.sz_pgto_alpha,   256);
    L.sz_pgto_radius  = sizeof(float)    * num_pgtos;

    L.off_pgto_ijkl   = ALIGN_TO(L.off_pgto_radius  + L.sz_pgto_radius,  256);
    L.sz_pgto_ijkl    = sizeof(uint32_t) * num_pgtos;

    L.off_matrix      = ALIGN_TO(L.off_pgto_ijkl    + L.sz_pgto_ijkl,    256);
    L.sz_matrix       = sizeof(float)    * matrix_data_len;

    L.total_size      = ALIGN_TO(L.off_matrix + L.sz_matrix, 256);
    return L;
}

bool md_gto_density_buf_create(md_gto_density_buf_t* buf, md_gpu_device_t device,
    const md_gto_basis_t* basis, const float* atom_xyz,
    const double* density_matrix) {
    ASSERT(buf);
    ASSERT(device);
    ASSERT(basis);
    ASSERT(atom_xyz);
    ASSERT(density_matrix);

    uint32_t num_cgtos = 0, num_pgtos = 0;
    gto_basis_count(&num_cgtos, &num_pgtos, basis);

    MEMSET(buf, 0, sizeof(*buf));
    buf->num_cgtos  = num_cgtos;
    buf->num_pgtos  = num_pgtos;
    buf->matrix_dim = num_cgtos;

    gto_buf_layout_t L = gto_density_buf_compute_layout(num_cgtos, num_pgtos, num_cgtos);

    buf->buffer = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){
        .size  = L.total_size,
        .flags = MD_GPU_BUFFER_CPU_VISIBLE,
    });
    if (!buf->buffer) {
        MD_LOG_ERROR("md_gto_density_buf_create: failed to allocate GPU buffer (%zu bytes)", L.total_size);
        return false;
    }

    size_t temp_pos = md_temp_get_pos();
    size_t tri_len  = ((size_t)num_cgtos * (num_cgtos + 1)) / 2;
    vec4_t*   cgto_xyzr   = (vec4_t*)  md_temp_push(sizeof(vec4_t)   * num_cgtos);
    uint32_t* cgto_offset = (uint32_t*)md_temp_push(sizeof(uint32_t) * (num_cgtos + 1));
    float*    pgto_alpha  = (float*)   md_temp_push(sizeof(float)    * num_pgtos);
    float*    pgto_coeff  = (float*)   md_temp_push(sizeof(float)    * num_pgtos);
    float*    pgto_radius = (float*)   md_temp_push(sizeof(float)    * num_pgtos);
    uint32_t* pgto_ijkl   = (uint32_t*)md_temp_push(sizeof(uint32_t) * num_pgtos);
    float*    matrix      = (float*)   md_temp_push(sizeof(float)    * tri_len);

    gto_expand_basis(cgto_xyzr, cgto_offset, pgto_alpha, pgto_coeff, pgto_radius, pgto_ijkl, basis, atom_xyz);
    density_matrix_to_upper_tri(matrix, density_matrix, num_cgtos);

    uint8_t* ptr = (uint8_t*)md_gpu_map_buffer(buf->buffer);
    MEMCPY(ptr + L.off_cgto_xyzr,   cgto_xyzr,   L.sz_cgto_xyzr);
    MEMCPY(ptr + L.off_cgto_offset, cgto_offset,  L.sz_cgto_offset);
    MEMCPY(ptr + L.off_pgto_coeff,  pgto_coeff,   L.sz_pgto_coeff);
    MEMCPY(ptr + L.off_pgto_alpha,  pgto_alpha,   L.sz_pgto_alpha);
    MEMCPY(ptr + L.off_pgto_radius, pgto_radius,  L.sz_pgto_radius);
    MEMCPY(ptr + L.off_pgto_ijkl,   pgto_ijkl,    L.sz_pgto_ijkl);
    MEMCPY(ptr + L.off_matrix,      matrix,        L.sz_matrix);
    md_gpu_unmap_buffer(buf->buffer);

    md_temp_set_pos_back(temp_pos);
    return true;
}

void md_gto_density_buf_update_positions(md_gto_density_buf_t* buf,
    const md_gto_basis_t* basis, const float* atom_xyz) {
    ASSERT(buf);
    ASSERT(buf->buffer);
    ASSERT(basis);
    ASSERT(atom_xyz);

    size_t temp_pos = md_temp_get_pos();
    vec4_t* cgto_xyzr = (vec4_t*)md_temp_push(sizeof(vec4_t) * buf->num_cgtos);

    uint32_t ci = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        int   nsph = 2 * (int)basis->shells[si].l + 1;
        float ax   = atom_xyz[basis->shells[si].atom_idx * 3 + 0];
        float ay   = atom_xyz[basis->shells[si].atom_idx * 3 + 1];
        float az   = atom_xyz[basis->shells[si].atom_idx * 3 + 2];
        for (int isph = 0; isph < nsph; isph++) {
            cgto_xyzr[ci++] = (vec4_t){ax, ay, az, FLT_MAX};
        }
    }

    gto_buf_layout_t L = gto_density_buf_compute_layout(buf->num_cgtos, buf->num_pgtos, buf->matrix_dim);
    uint8_t* ptr = (uint8_t*)md_gpu_map_buffer(buf->buffer);
    MEMCPY(ptr + L.off_cgto_xyzr, cgto_xyzr, L.sz_cgto_xyzr);
    md_gpu_unmap_buffer(buf->buffer);

    md_temp_set_pos_back(temp_pos);
}

void md_gto_density_buf_update_matrix(md_gto_density_buf_t* buf,
    const double* density_matrix) {
    ASSERT(buf);
    ASSERT(buf->buffer);
    ASSERT(density_matrix);

    size_t tri_len  = ((size_t)buf->matrix_dim * (buf->matrix_dim + 1)) / 2;
    size_t temp_pos = md_temp_get_pos();
    float* matrix   = (float*)md_temp_push(sizeof(float) * tri_len);
    density_matrix_to_upper_tri(matrix, density_matrix, buf->matrix_dim);

    gto_buf_layout_t L = gto_density_buf_compute_layout(buf->num_cgtos, buf->num_pgtos, buf->matrix_dim);
    uint8_t* ptr = (uint8_t*)md_gpu_map_buffer(buf->buffer);
    MEMCPY(ptr + L.off_matrix, matrix, L.sz_matrix);
    md_gpu_unmap_buffer(buf->buffer);

    md_temp_set_pos_back(temp_pos);
}

void md_gto_density_buf_destroy(md_gto_density_buf_t* buf) {
    if (!buf) return;
    if (buf->buffer) {
        md_gpu_destroy_buffer(buf->buffer);
    }
    MEMSET(buf, 0, sizeof(*buf));
}

bool md_gto_grid_evaluate_density_gpu(md_gpu_device_t device, md_gto_density_buf_t* buf, md_gpu_image_t out_image, const md_grid_t* grid, md_gpu_fence_t* out_fence) {
    if (!device || !buf || !buf->buffer || !out_image || !grid || !out_fence) {
        MD_LOG_ERROR("md_gto_grid_evaluate_density_gpu: invalid input");
        return false;
    }

    md_gpu_compute_pipeline_t pipeline = gto_ensure_density_pipeline(device);
    if (!pipeline) return false;

    typedef struct {
        float    world_to_model[4][4];
        float    index_to_world[4][4];
        float    step[4];
        uint32_t D_matrix_dim;
        uint32_t _pad[3];
    } ubo_t;

    gto_buf_layout_t L = gto_density_buf_compute_layout(buf->num_cgtos, buf->num_pgtos, buf->matrix_dim);

    ubo_t ubo = {0};
    world_to_model_matrix(ubo.world_to_model, grid);
    index_to_world_matrix(ubo.index_to_world, grid);
    ubo.step[0]      = grid->spacing.elem[0];
    ubo.step[1]      = grid->spacing.elem[1];
    ubo.step[2]      = grid->spacing.elem[2];
    ubo.step[3]      = 0.0f;
    ubo.D_matrix_dim = buf->matrix_dim;

    md_gpu_queue_t          queue = md_gpu_acquire_compute_queue(device);
    md_gpu_command_buffer_t cmd   = md_gpu_acquire_command_buffer(device);

    md_gpu_cmd_bind_compute_pipeline(cmd, pipeline);
    md_gpu_cmd_push_constants(cmd, &ubo, sizeof(ubo));
    md_gpu_cmd_bind_buffer_range(cmd, 0, buf->buffer, L.off_cgto_xyzr,   L.sz_cgto_xyzr);
    md_gpu_cmd_bind_buffer_range(cmd, 1, buf->buffer, L.off_cgto_offset,  L.sz_cgto_offset);
    md_gpu_cmd_bind_buffer_range(cmd, 2, buf->buffer, L.off_pgto_coeff,   L.sz_pgto_coeff);
    md_gpu_cmd_bind_buffer_range(cmd, 3, buf->buffer, L.off_pgto_alpha,   L.sz_pgto_alpha);
    md_gpu_cmd_bind_buffer_range(cmd, 4, buf->buffer, L.off_pgto_radius,  L.sz_pgto_radius);
    md_gpu_cmd_bind_buffer_range(cmd, 5, buf->buffer, L.off_pgto_ijkl,    L.sz_pgto_ijkl);
    md_gpu_cmd_bind_buffer_range(cmd, 6, buf->buffer, L.off_matrix,       L.sz_matrix);
    md_gpu_cmd_bind_image(cmd, 0, out_image);
    md_gpu_cmd_dispatch(cmd, grid->dim[0], grid->dim[1], grid->dim[2]);

    *out_fence = md_gpu_queue_submit(queue, cmd);
    return true;
}

#endif // MD_ENABLE_GPU

static inline float fast_powf(float base, int exp) {
    float val = 1.0f;
    switch(exp) {
    case 4: val *= base; FALLTHROUGH;
    case 3: val *= base; FALLTHROUGH;
    case 2: val *= base; FALLTHROUGH;
    case 1: val *= base; FALLTHROUGH;
    case 0: break;
    }
    return val;
}

static inline double fast_pow(double base, int exp){
    double val = 1.0;
    switch(exp) {
    case 4: val *= base; FALLTHROUGH;
    case 3: val *= base; FALLTHROUGH;
    case 2: val *= base; FALLTHROUGH;
    case 1: val *= base; FALLTHROUGH;
    case 0: break;
    }
    return val;
}


static inline md_128 md_mm_fast_pow1(md_128 base, int exp) {
    switch (exp) {
    case 1:
        return base;
    case 2:
        return md_mm_mul_ps(base, base);
    case 3:
        return md_mm_mul_ps(base, md_mm_mul_ps(base, base));
    case 4: {
        md_128 squared = md_mm_mul_ps(base, base);
        return md_mm_mul_ps(squared, squared);
    }
    case 0:
    default:
        return md_mm_set1_ps(1.0f);
    }
}

static inline md_256 md_mm256_fast_pow1(md_256 base, int exp) {
    switch (exp) {
    case 1:
        return base;
    case 2:
        return md_mm256_mul_ps(base, base);
    case 3:
        return md_mm256_mul_ps(base, md_mm256_mul_ps(base, base));
    case 4: {
        md_256 squared = md_mm256_mul_ps(base, base);
        return md_mm256_mul_ps(squared, squared);
    }
    case 0:
    default:
        return md_mm256_set1_ps(1.0f);
    }
}

static inline md_128 md_mm_fast_pow(md_128 base1, md_128i exp) {
    md_128 base2 = md_mm_mul_ps(base1, base1);
    md_128 base3 = md_mm_mul_ps(base2, base1);
    md_128 base4 = md_mm_mul_ps(base2, base2);

    md_128 mask1 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(1)));
    md_128 mask2 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(2)));
    md_128 mask3 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(3)));
    md_128 mask4 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(4)));

    md_128 res = md_mm_set1_ps(1.0f);
    res = md_mm_blendv_ps(res, base1, mask1);
    res = md_mm_blendv_ps(res, base2, mask2);
    res = md_mm_blendv_ps(res, base3, mask3);
    res = md_mm_blendv_ps(res, base4, mask4);
    return res;
}

static inline md_256 md_mm256_fast_pow(md_256 base1, md_256i exp) {
    md_256 base2 = md_mm256_mul_ps(base1, base1);
    md_256 base3 = md_mm256_mul_ps(base2, base1);
    md_256 base4 = md_mm256_mul_ps(base2, base2);

    md_256 mask1 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(1)));
    md_256 mask2 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(2)));
    md_256 mask3 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(3)));
    md_256 mask4 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(4)));

    md_256 res = md_mm256_set1_ps(1.0f);
    res = md_mm256_blendv_ps(res, base1, mask1);
    res = md_mm256_blendv_ps(res, base2, mask2);
    res = md_mm256_blendv_ps(res, base3, mask3);
    res = md_mm256_blendv_ps(res, base4, mask4);
    return res;
}

#ifdef __AVX512F__
static inline __m512 md_mm512_fast_pow(__m512 base1, __m512i exp) {
    __m512 base2 = _mm512_mul_ps(base1,  base1);
    __m512 base3 = _mm512_mul_ps(base2, base1);
    __m512 base4 = _mm512_mul_ps(base2, base2);

    __mmask16 mask1 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(1), _MM_CMPINT_EQ);
    __mmask16 mask2 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(2), _MM_CMPINT_EQ);
    __mmask16 mask3 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(3), _MM_CMPINT_EQ);
    __mmask16 mask4 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(4), _MM_CMPINT_EQ);

    __m512 res = _mm512_set1_ps(1.0f);
    res = _mm512_mask_blend_ps(mask1, res, base1);
    res = _mm512_mask_blend_ps(mask2, res, base2);
    res = _mm512_mask_blend_ps(mask3, res, base3);
    res = _mm512_mask_blend_ps(mask4, res, base4);
    return res;
}
#endif

static inline void evaluate_grid_ref(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_step_x[3], const float grid_step_y[3], const float grid_step_z[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    for (int iz = grid_idx_min[2]; iz < grid_idx_max[2]; iz++) {
        const int z_stride = iz * grid_dim[0] * grid_dim[1];
        for (int iy = grid_idx_min[1]; iy < grid_idx_max[1]; ++iy) {
            const int y_stride = iy * grid_dim[0];
            for (int ix = grid_idx_min[0]; ix < grid_idx_max[0]; ++ix) {
                const int x_stride = ix;

                float x = grid_origin[0] + ix * grid_step_x[0] + iy * grid_step_y[0] + iz * grid_step_z[0];
                float y = grid_origin[1] + ix * grid_step_x[1] + iy * grid_step_y[1] + iz * grid_step_z[1];
                float z = grid_origin[2] + ix * grid_step_x[2] + iy * grid_step_y[2] + iz * grid_step_z[2];

                double psi = 0.0;
                for (size_t i = 0; i < num_gtos; ++i) {
                    float px	= gtos[i].x;
                    float py	= gtos[i].y;
                    float pz	= gtos[i].z;
                    float alpha	= gtos[i].alpha;
                    float coeff	= gtos[i].coeff;
                    int   pi	= gtos[i].i;
                    int   pj	= gtos[i].j;
                    int   pk	= gtos[i].k;

                    float dx = x - px;
                    float dy = y - py;
                    float dz = z - pz;
                    float d2 = dx * dx + dy * dy + dz * dz;
                    float fx = powf(dx, (float)pi);
                    float fy = powf(dy, (float)pj);
                    float fz = powf(dz, (float)pk);
                    float exp_term = (alpha == 0.0f) ? 1.0f : expf(-alpha * d2);
                    float powxyz = fx * fy * fz;
                    float prod = coeff * powxyz * exp_term;
                    psi += prod;
                }

                if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                    psi *= psi;
                }

                int index = x_stride + y_stride + z_stride;
                grid_data[index] += (float)psi;
            }
        }
    }
}

#if 0
// These are old code-paths that were vectorized over the PGTOs
// It is not ideal as the occupation of vector registers is quite poor due to the large variation
// Of how many PGTOs that influence a spatial region

#ifdef __AVX512F__
static inline void evaluate_grid_512(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_stepsize[3], const md_gto_data_t* gto, md_gto_eval_mode_t mode) {
    for (int iz = grid_idx_min[2]; iz < grid_idx_max[2]; iz++) {
        float z = grid_origin[2] + iz * grid_stepsize[2];
        __m512 vz = _mm512_set1_ps(z);
        int z_stride = iz * grid_dim[0] * grid_dim[1];
        for (int iy = grid_idx_min[1]; iy < grid_idx_max[1]; ++iy) {
            float y = grid_origin[1] + iy * grid_stepsize[1];
            __m512 vy = _mm512_set1_ps(y);
            int y_stride = iy * grid_dim[0];
            for (int ix = grid_idx_min[0]; ix < grid_idx_max[0]; ++ix) {
                float x = grid_origin[0] + ix * grid_stepsize[0];
                __m512 vx = _mm512_set1_ps(x);
                int x_stride = ix;

                __m512 vpsi = _mm512_setzero_ps();
                for (size_t i = 0; i < gto->count; i += 16) {
                    __m512  px = _mm512_loadu_ps(gto->x + i);
                    __m512  py = _mm512_loadu_ps(gto->y + i);
                    __m512  pz = _mm512_loadu_ps(gto->z + i);
                    __m512  pa = _mm512_loadu_ps(gto->neg_alpha + i);
                    __m512  pc = _mm512_loadu_ps(gto->coeff + i);
                    __m512i pi = _mm512_loadu_si512(gto->i + i);
                    __m512i pj = _mm512_loadu_si512(gto->j + i);
                    __m512i pk = _mm512_loadu_si512(gto->k + i);

                    __m512 dx = _mm512_sub_ps(vx, px);
                    __m512 dy = _mm512_sub_ps(vy, py);
                    __m512 dz = _mm512_sub_ps(vz, pz);
                    __m512 d2 = _mm512_fmadd_ps(dx, dx, _mm512_fmadd_ps(dy, dy, _mm512_mul_ps(dz, dz)));
                    __m512 fx = md_mm512_fast_pow(dx, pi);
                    __m512 fy = md_mm512_fast_pow(dy, pj);
                    __m512 fz = md_mm512_fast_pow(dz, pk);
                    __m512 ex = md_mm512_exp_ps(_mm512_mul_ps(pa, d2));

                    __m512 prod = _mm512_mul_ps(_mm512_mul_ps(_mm512_mul_ps(pc, fx), _mm512_mul_ps(fy, fz)), ex);

                    if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                        prod = _mm512_mul_ps(prod, prod);
                    }

                    vpsi = _mm512_add_ps(vpsi, prod);
                }

                float psi = _mm512_reduce_add_ps(vpsi);
                int index = x_stride + y_stride + z_stride;
                grid_data[index] = psi;
            }
        }
    }
}
#endif

static inline void evaluate_grid_256(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_stepsize[3], const md_gto_data_t* gto, md_gto_eval_mode_t mode) {
    for (int iz = grid_idx_min[2]; iz < grid_idx_max[2]; iz++) {
        float z = grid_origin[2] + iz * grid_stepsize[2];
        md_256 vz = md_mm256_set1_ps(z);
        int z_stride = iz * grid_dim[0] * grid_dim[1];
        for (int iy = grid_idx_min[1]; iy < grid_idx_max[1]; ++iy) {
            float y = grid_origin[1] + iy * grid_stepsize[1];
            md_256 vy = md_mm256_set1_ps(y);
            int y_stride = iy * grid_dim[0];
            for (int ix = grid_idx_min[0]; ix < grid_idx_max[0]; ++ix) {
                float x = grid_origin[0] + ix * grid_stepsize[0];
                md_256 vx = md_mm256_set1_ps(x);
                int x_stride = ix;

                md_256 vpsi = md_mm256_setzero_ps();
                for (size_t i = 0; i < gto->count; i += 8) {
                    md_256  px = md_mm256_loadu_ps(gto->x + i);
                    md_256  py = md_mm256_loadu_ps(gto->y + i);
                    md_256  pz = md_mm256_loadu_ps(gto->z + i);
                    md_256  pa = md_mm256_loadu_ps(gto->neg_alpha + i);
                    md_256  pc = md_mm256_loadu_ps(gto->coeff + i);
                    md_256i pi = md_mm256_loadu_si256((const md_256i*)(gto->i + i));
                    md_256i pj = md_mm256_loadu_si256((const md_256i*)(gto->j + i));
                    md_256i pk = md_mm256_loadu_si256((const md_256i*)(gto->k + i));

                    md_256 dx = md_mm256_sub_ps(vx, px);
                    md_256 dy = md_mm256_sub_ps(vy, py);
                    md_256 dz = md_mm256_sub_ps(vz, pz);
                    md_256 d2 = md_mm256_fmadd_ps(dx, dx, md_mm256_fmadd_ps(dy, dy, md_mm256_mul_ps(dz, dz)));
                    md_256 fx = md_mm256_fast_pow(dx, pi);
                    md_256 fy = md_mm256_fast_pow(dy, pj);
                    md_256 fz = md_mm256_fast_pow(dz, pk);
                    md_256 ex = md_mm256_exp_ps(md_mm256_mul_ps(pa, d2));

                    md_256 prod = md_mm256_mul_ps(md_mm256_mul_ps(md_mm256_mul_ps(pc, fx), md_mm256_mul_ps(fy, fz)), ex);
                    if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                        prod = _mm256_mul_ps(prod, prod);
                    }

                    vpsi = md_mm256_add_ps(vpsi, prod);
                }

                float psi = md_mm256_reduce_add_ps(vpsi);
                int index = x_stride + y_stride + z_stride;
                grid_data[index] = psi;
            }
        }
    }
}

static inline void evaluate_grid_128(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_stepsize[3], const md_gto_data_t* gto, md_gto_eval_mode_t mode) {
    for (int iz = grid_idx_min[2]; iz < grid_idx_max[2]; iz++) {
        float z = grid_origin[2] + iz * grid_stepsize[2];
        md_128 vz = md_mm_set1_ps(z);
        int z_stride = iz * grid_dim[0] * grid_dim[1];
        for (int iy = grid_idx_min[1]; iy < grid_idx_max[1]; ++iy) {
            float y = grid_origin[1] + iy * grid_stepsize[1];
            md_128 vy = md_mm_set1_ps(y);
            int y_stride = iy * grid_dim[0];
            for (int ix = grid_idx_min[0]; ix < grid_idx_max[0]; ++ix) {
                float x = grid_origin[0] + ix * grid_stepsize[0];
                md_128 vx = md_mm_set1_ps(x);
                int x_stride = ix;

                md_128 vpsi = md_mm_setzero_ps();
                for (size_t i = 0; i < gto->count; i += 4) {
                    md_128  px = md_mm_loadu_ps(gto->x + i);
                    md_128  py = md_mm_loadu_ps(gto->y + i);
                    md_128  pz = md_mm_loadu_ps(gto->z + i);
                    md_128  pa = md_mm_loadu_ps(gto->neg_alpha + i);
                    md_128  pc = md_mm_loadu_ps(gto->coeff + i);
                    md_128i pi = md_mm_loadu_si128((const md_128i*)(gto->i + i));
                    md_128i pj = md_mm_loadu_si128((const md_128i*)(gto->j + i));
                    md_128i pk = md_mm_loadu_si128((const md_128i*)(gto->k + i));

                    md_128 dx = md_mm_sub_ps(vx, px);
                    md_128 dy = md_mm_sub_ps(vy, py);
                    md_128 dz = md_mm_sub_ps(vz, pz);
                    md_128 d2 = md_mm_fmadd_ps(dx, dx, md_mm_fmadd_ps(dy, dy, md_mm_mul_ps(dz, dz)));
                    md_128 fx = md_mm_fast_pow(dx, pi);
                    md_128 fy = md_mm_fast_pow(dy, pj);
                    md_128 fz = md_mm_fast_pow(dz, pk);
                    md_128 ex = md_mm_exp_ps(md_mm_mul_ps(pa, d2));

                    md_128 prod = md_mm_mul_ps(md_mm_mul_ps(md_mm_mul_ps(pc, fx), md_mm_mul_ps(fy, fz)), ex);
                    if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                        prod = _mm_mul_ps(prod, prod);
                    }

                    vpsi = md_mm_add_ps(vpsi, prod);
                }

                float psi = md_mm_reduce_add_ps(vpsi);
                int index = x_stride + y_stride + z_stride;
                grid_data[index] = psi;
            }
        }
    }
}
#endif

#if defined(__AVX512F__) && defined(__AVX512DQ__)

// Evaluate 8 voxels per gto
static inline void evaluate_grid_ortho_8x8x8_512(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    const md_256i vix = md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[0]), md_mm256_set_epi32(7,6,5,4,3,2,1,0));
    const md_256  vxh = md_mm256_fmadd_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step[0]), md_mm256_set1_ps(grid_origin[0]));
    const __m512  vx  = _mm512_insertf32x8(_mm512_castps256_ps512(vxh), vxh, 1);
    const int x_stride = grid_idx_min[0];

    __m512 vpsi[8][4] = {0};

    for (size_t i = 0; i < num_gtos; ++i) {
        const __m512  px = _mm512_set1_ps(gtos[i].x);
        const __m512  py = _mm512_set1_ps(gtos[i].y);
        const __m512  pz = _mm512_set1_ps(gtos[i].z);
        const __m512  pc = _mm512_set1_ps(gtos[i].coeff);
        const __m512  pa = _mm512_set1_ps(-gtos[i].alpha); // Negate alpha here
        const __m512i pi = _mm512_set1_epi32(gtos[i].i);
        const __m512i pj = _mm512_set1_epi32(gtos[i].j);
        const __m512i pk = _mm512_set1_epi32(gtos[i].k);

        for (int iz = 0; iz < 8; ++iz) {
            float z = grid_origin[2] + (grid_idx_min[2] + iz) * grid_step[2];
            __m512 vz = _mm512_set1_ps(z);
            for (int iy = 0; iy < 4; ++iy) {
                float y[2] = {
                    grid_origin[1] + (grid_idx_min[1] + iy * 2 + 0) * grid_step[1],
                    grid_origin[1] + (grid_idx_min[1] + iy * 2 + 1) * grid_step[1],
                };
                __m512 vy = _mm512_insertf32x8(_mm512_set1_ps(y[0]), _mm256_set1_ps(y[1]), 1);

                __m512 dx = _mm512_sub_ps(vx, px);
                __m512 dy = _mm512_sub_ps(vy, py);
                __m512 dz = _mm512_sub_ps(vz, pz);
                __m512 d2 = _mm512_fmadd_ps(dx, dx, _mm512_fmadd_ps(dy, dy, _mm512_mul_ps(dz, dz)));
                __m512 fx = md_mm512_fast_pow(dx, pi);
                __m512 fy = md_mm512_fast_pow(dy, pj);
                __m512 fz = md_mm512_fast_pow(dz, pk);
                __m512 ex = md_mm512_exp_ps(_mm512_mul_ps(pa, d2));
                __m512 prod = _mm512_mul_ps(_mm512_mul_ps(_mm512_mul_ps(pc, fx), _mm512_mul_ps(fy, fz)), ex);

                vpsi[iz][iy] = _mm512_add_ps(vpsi[iz][iy], prod);
            }
        }
    }

    // Write result block to memory
    for (int iz = 0; iz < 8; ++iz) {
        int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
        for (int iy = 0; iy < 4; ++iy) {
            int y_stride[2] = {
                (grid_idx_min[1] + iy * 2 + 0) * grid_dim[0],
                (grid_idx_min[1] + iy * 2 + 1) * grid_dim[0],
            };
            int index[2] = {
                x_stride + y_stride[0] + z_stride,
                x_stride + y_stride[1] + z_stride,
            };

            md_512 psi = vpsi[iz][iy];
            if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                psi = _mm512_mul_ps(psi, psi);
            }

            md_256 tpsi[2] = {
                _mm512_castps512_ps256(psi),
                _mm512_extractf32x8_ps(psi, 1),
            };

            md_mm256_storeu_ps(grid_data + index[0], tpsi[0]);
            md_mm256_storeu_ps(grid_data + index[1], tpsi[1]);
        }
    }
}

#endif

// Evaluate 8 voxels per gto
static inline void evaluate_grid_ortho_8x8x8_256(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    const md_256i vix  = md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[0]), md_mm256_set_epi32(7,6,5,4,3,2,1,0));
    const int x_stride = grid_idx_min[0];
    const md_256   vx  = md_mm256_fmadd_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step[0]), md_mm256_set1_ps(grid_origin[0]));

    // Operate on local block to avoid cache-line contention across threads
    md_256 vpsi[8][8] = {0};

    for (size_t i = 0; i < num_gtos; ++i) {
        const md_256  px = md_mm256_set1_ps(gtos[i].x);
        const md_256  py = md_mm256_set1_ps(gtos[i].y);
        const md_256  pz = md_mm256_set1_ps(gtos[i].z);
        const md_256  pc = md_mm256_set1_ps(gtos[i].coeff);
        const md_256  pa = md_mm256_set1_ps(-gtos[i].alpha); // Negate alpha here
        const md_256i pi = md_mm256_set1_epi32(gtos[i].i);
        const md_256i pj = md_mm256_set1_epi32(gtos[i].j);
        const md_256i pk = md_mm256_set1_epi32(gtos[i].k);

        for (int iz = 0; iz < 8; ++iz) {
            float z = grid_origin[2] + (grid_idx_min[2] + iz) * grid_step[2];
            md_256 vz = md_mm256_set1_ps(z);

            for (int iy = 0; iy < 8; ++iy) {
                float y = grid_origin[1] + (grid_idx_min[1] + iy) * grid_step[1];
                md_256 vy = md_mm256_set1_ps(y);

                md_256 dx = md_mm256_sub_ps(vx, px);
                md_256 dy = md_mm256_sub_ps(vy, py);
                md_256 dz = md_mm256_sub_ps(vz, pz);
                md_256 d2 = md_mm256_fmadd_ps(dx, dx, md_mm256_fmadd_ps(dy, dy, md_mm256_mul_ps(dz, dz)));
                md_256 ex = md_mm256_exp_ps(md_mm256_mul_ps(pa, d2));
                md_256 fx = md_mm256_fast_pow(dx, pi);
                md_256 fy = md_mm256_fast_pow(dy, pj);
                md_256 fz = md_mm256_fast_pow(dz, pk);

                md_256 prod_a = md_mm256_mul_ps(pc, fx);
                md_256 prod_b = md_mm256_mul_ps(fy, fz);

                vpsi[iz][iy] = md_mm256_fmadd_ps(md_mm256_mul_ps(prod_a, prod_b), ex, vpsi[iz][iy]);
            }
        }
    }

    // Write result block to memory
    for (int iz = 0; iz < 8; ++iz) {
        int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
        for (int iy = 0; iy < 8; ++iy) {
            int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
            int index = x_stride + y_stride + z_stride;
            md_256 psi = vpsi[iz][iy];

            if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                psi = md_mm256_mul_ps(psi, psi);
            }

            md_mm256_storeu_ps(grid_data + index, md_mm256_add_ps(md_mm256_loadu_ps(grid_data + index), psi));
        }
    }
}

// Evaluate 8 voxels per gto
static inline void evaluate_grid_8x8x8_256(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step_x[3], const float grid_step_y[3], const float grid_step_z[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    const md_256i vix = md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[0]), md_mm256_set_epi32(7,6,5,4,3,2,1,0));
    const int x_stride = grid_idx_min[0];

    const md_256 gsx[3] = {
        md_mm256_add_ps(md_mm256_set1_ps(grid_origin[0]), md_mm256_mul_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step_x[0]))),
        md_mm256_add_ps(md_mm256_set1_ps(grid_origin[1]), md_mm256_mul_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step_x[1]))),
        md_mm256_add_ps(md_mm256_set1_ps(grid_origin[2]), md_mm256_mul_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step_x[2]))),
    };

    const md_256 gsz[3] = {
        md_mm256_set1_ps(grid_step_z[0]),
        md_mm256_set1_ps(grid_step_z[1]),
        md_mm256_set1_ps(grid_step_z[2]),
    };

    const md_256 gsy[3] = {
        md_mm256_set1_ps(grid_step_y[0]),
        md_mm256_set1_ps(grid_step_y[1]),
        md_mm256_set1_ps(grid_step_y[2]),
    };

    // Operate on local block to avoid cache-line contention across threads
    md_256 vpsi[8][8] = {0};

    for (size_t gto_idx = 0; gto_idx < num_gtos; ++gto_idx) {
        const float px = gtos[gto_idx].x;
        const float py = gtos[gto_idx].y;
        const float pz = gtos[gto_idx].z;
        const float pc = gtos[gto_idx].coeff;
        const float pa = -gtos[gto_idx].alpha; // Negate alpha here
        const int pi = gtos[gto_idx].i;
        const int pj = gtos[gto_idx].j;
        const int pk = gtos[gto_idx].k;
    
        for (int iz = 0; iz < 8; ++iz) {
            const md_256 tz = md_mm256_cvtepi32_ps(md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[2]), md_mm256_set1_epi32(iz)));

            const md_256 xz[3] = {
                md_mm256_fmadd_ps(tz, gsz[0], gsx[0]),
                md_mm256_fmadd_ps(tz, gsz[1], gsx[1]),
                md_mm256_fmadd_ps(tz, gsz[2], gsx[2]),
            };

            for (int iy = 0; iy < 8; ++iy) {
                const md_256 ty = md_mm256_cvtepi32_ps(md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[1]), md_mm256_set1_epi32(iy)));

                md_256 vx = md_mm256_fmadd_ps(ty, gsy[0], xz[0]);
                md_256 vy = md_mm256_fmadd_ps(ty, gsy[1], xz[1]);
                md_256 vz = md_mm256_fmadd_ps(ty, gsy[2], xz[2]);

                md_256 dx = md_mm256_sub_ps(vx, md_mm256_set1_ps(px));
                md_256 dy = md_mm256_sub_ps(vy, md_mm256_set1_ps(py));
                md_256 dz = md_mm256_sub_ps(vz, md_mm256_set1_ps(pz));
                md_256 d2 = md_mm256_fmadd_ps(dx, dx, md_mm256_fmadd_ps(dy, dy, md_mm256_mul_ps(dz, dz)));
                md_256 ex = md_mm256_exp_ps(md_mm256_mul_ps(md_mm256_set1_ps(pa), d2));
                md_256 fx = md_mm256_fast_pow1(dx, pi);
                md_256 fy = md_mm256_fast_pow1(dy, pj);
                md_256 fz = md_mm256_fast_pow1(dz, pk);

                md_256 prod_a = md_mm256_mul_ps(md_mm256_set1_ps(pc), fx);
                md_256 prod_b = md_mm256_mul_ps(fy, fz);

                vpsi[iz][iy] = md_mm256_fmadd_ps(md_mm256_mul_ps(prod_a, prod_b), ex, vpsi[iz][iy]);
            }
        }
    }

    // Write result block to memory
    for (int iz = 0; iz < 8; ++iz) {
        int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
        for (int iy = 0; iy < 8; ++iy) {
            int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
            int index = x_stride + y_stride + z_stride;
            md_256 psi = vpsi[iz][iy];

            if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                psi = md_mm256_mul_ps(psi, psi);
            }

            md_mm256_storeu_ps(grid_data + index, md_mm256_add_ps(md_mm256_loadu_ps(grid_data + index), psi));
        }
    }
}

// Evaluate 8 voxels per gto
static inline void evaluate_grid_ortho_8x8x8_128(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    const md_128i vix[2] = {
        md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0] + 0), md_mm_set_epi32(3,2,1,0)),
        md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0] + 4), md_mm_set_epi32(3,2,1,0)),
    };
    const md_128   vx[2] = {
        md_mm_fmadd_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step[0]), md_mm_set1_ps(grid_origin[0])),
        md_mm_fmadd_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step[0]), md_mm_set1_ps(grid_origin[0])),
    };
    const int x_stride[2] = {
        grid_idx_min[0] + 0,
        grid_idx_min[0] + 4,
    };

    // Operate on local block to avoid cache-line contention across threads
    md_128 vpsi[8][8][2] = {0};

    for (size_t i = 0; i < num_gtos; ++i) {
        const md_128  px = md_mm_set1_ps(gtos[i].x);
        const md_128  py = md_mm_set1_ps(gtos[i].y);
        const md_128  pz = md_mm_set1_ps(gtos[i].z);
        const md_128  pc = md_mm_set1_ps(gtos[i].coeff);
        const md_128  pa = md_mm_set1_ps(-gtos[i].alpha); // Negate alpha here
        const md_128i pi = md_mm_set1_epi32(gtos[i].i);
        const md_128i pj = md_mm_set1_epi32(gtos[i].j);
        const md_128i pk = md_mm_set1_epi32(gtos[i].k);

        for (int iz = 0; iz < 8; ++iz) {
            md_128 vz = md_mm_set1_ps(grid_origin[2] + (grid_idx_min[2] + iz) * grid_step[2]);
            for (int iy = 0; iy < 8; ++iy) {
                md_128 vy = md_mm_set1_ps(grid_origin[1] + (grid_idx_min[1] + iy) * grid_step[1]);

                for (int ix = 0; ix < 2; ++ix) {
                    md_128 dx = md_mm_sub_ps(vx[ix], px);
                    md_128 dy = md_mm_sub_ps(vy,	 py);
                    md_128 dz = md_mm_sub_ps(vz,	 pz);
                    md_128 d2 = md_mm_fmadd_ps(dx, dx, md_mm_fmadd_ps(dy, dy, md_mm_mul_ps(dz, dz)));
                    md_128 fx = md_mm_fast_pow(dx, pi);
                    md_128 fy = md_mm_fast_pow(dy, pj);
                    md_128 fz = md_mm_fast_pow(dz, pk);
                    md_128 ex = md_mm_exp_ps(md_mm_mul_ps(pa, d2));
                    md_128 prod = md_mm_mul_ps(md_mm_mul_ps(md_mm_mul_ps(pc, fx), md_mm_mul_ps(fy, fz)), ex);

                    vpsi[iz][iy][ix] = md_mm_add_ps(vpsi[iz][iy][ix], prod);
                }
            }
        }
    }

    // Write result block to memory
    for (int iz = 0; iz < 8; ++iz) {
        int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
        for (int iy = 0; iy < 8; ++iy) {
            int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
            int index[2] = {
                x_stride[0] + y_stride + z_stride,
                x_stride[1] + y_stride + z_stride,
            };

            md_128 psi[2] = {
                vpsi[iz][iy][0],
                vpsi[iz][iy][1],
            };

            if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                psi[0] = md_mm_mul_ps(psi[0], psi[0]);
                psi[1] = md_mm_mul_ps(psi[1], psi[1]);
            }

            md_mm_storeu_ps(grid_data + index[0], md_mm_add_ps(md_mm_loadu_ps(grid_data + index[0]), psi[0]));
            md_mm_storeu_ps(grid_data + index[1], md_mm_add_ps(md_mm_loadu_ps(grid_data + index[1]), psi[1]));
        }
    }
}

// Evaluate 8 voxels per gto
static inline void evaluate_grid_8x8x8_128(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step_x[3], const float grid_step_y[3], const float grid_step_z[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    const int x_stride[2] = {
        grid_idx_min[0] + 0,
        grid_idx_min[0] + 4,
    };

    const md_128i vix[2] = {
        md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0])    , md_mm_set_epi32(3,2,1,0)),
        md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0] + 4), md_mm_set_epi32(3,2,1,0)),
    };

    const md_128 gsx[2][3] = {
        {
            md_mm_add_ps(md_mm_set1_ps(grid_origin[0]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step_x[0]))),
            md_mm_add_ps(md_mm_set1_ps(grid_origin[1]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step_x[1]))),
            md_mm_add_ps(md_mm_set1_ps(grid_origin[2]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step_x[2]))),
        },
        {
            md_mm_add_ps(md_mm_set1_ps(grid_origin[0]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step_x[0]))),
            md_mm_add_ps(md_mm_set1_ps(grid_origin[1]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step_x[1]))),
            md_mm_add_ps(md_mm_set1_ps(grid_origin[2]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step_x[2]))),
        },
    };

    const md_128 gsz[3] = {
        md_mm_set1_ps(grid_step_z[0]),
        md_mm_set1_ps(grid_step_z[1]),
        md_mm_set1_ps(grid_step_z[2]),
    };

    const md_128 gsy[3] = {
        md_mm_set1_ps(grid_step_y[0]),
        md_mm_set1_ps(grid_step_y[1]),
        md_mm_set1_ps(grid_step_y[2]),
    };

    // Operate on local block to avoid cache-line contention across threads
    md_128 vpsi[8][8][2] = {0};

    for (size_t i = 0; i < num_gtos; ++i) {
        const md_128  px = md_mm_set1_ps(gtos[i].x);
        const md_128  py = md_mm_set1_ps(gtos[i].y);
        const md_128  pz = md_mm_set1_ps(gtos[i].z);
        const md_128  pc = md_mm_set1_ps(gtos[i].coeff);
        const md_128  pa = md_mm_set1_ps(-gtos[i].alpha); // Negate alpha here
        const md_128i pi = md_mm_set1_epi32(gtos[i].i);
        const md_128i pj = md_mm_set1_epi32(gtos[i].j);
        const md_128i pk = md_mm_set1_epi32(gtos[i].k);

        for (int iz = 0; iz < 8; ++iz) {
            const md_128 tz = md_mm_cvtepi32_ps(md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[2]), md_mm_set1_epi32(iz)));
            const md_128 xz[2][3] = {
                {
                    md_mm_fmadd_ps(tz, gsz[0], gsx[0][0]),
                    md_mm_fmadd_ps(tz, gsz[1], gsx[0][1]),
                    md_mm_fmadd_ps(tz, gsz[2], gsx[0][2]),
                },
                {
                    md_mm_fmadd_ps(tz, gsz[0], gsx[1][0]),
                    md_mm_fmadd_ps(tz, gsz[1], gsx[1][1]),
                    md_mm_fmadd_ps(tz, gsz[2], gsx[1][2]),
                },
            };
            for (int iy = 0; iy < 8; ++iy) {
                const md_128 ty = md_mm_cvtepi32_ps(md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[1]), md_mm_set1_epi32(iy)));

                const md_128 vx[2] = {
                    md_mm_fmadd_ps(ty, gsy[0], xz[0][0]),
                    md_mm_fmadd_ps(ty, gsy[0], xz[1][0]),
                };
                const md_128 vy[2] = {
                    md_mm_fmadd_ps(ty, gsy[1], xz[0][1]),
                    md_mm_fmadd_ps(ty, gsy[1], xz[1][1]),
                };
                const md_128 vz[2] = {
                    md_mm_fmadd_ps(ty, gsy[2], xz[0][2]),
                    md_mm_fmadd_ps(ty, gsy[2], xz[1][2]),
                };

                for (int ix = 0; ix < 2; ++ix) {
                    md_128 dx = md_mm_sub_ps(vx[ix], px);
                    md_128 dy = md_mm_sub_ps(vy[ix], py);
                    md_128 dz = md_mm_sub_ps(vz[ix], pz);
                    md_128 d2 = md_mm_fmadd_ps(dx, dx, md_mm_fmadd_ps(dy, dy, md_mm_mul_ps(dz, dz)));
                    md_128 fx = md_mm_fast_pow(dx, pi);
                    md_128 fy = md_mm_fast_pow(dy, pj);
                    md_128 fz = md_mm_fast_pow(dz, pk);
                    md_128 ex = md_mm_exp_ps(md_mm_mul_ps(pa, d2));
                    md_128 prod = md_mm_mul_ps(md_mm_mul_ps(md_mm_mul_ps(pc, fx), md_mm_mul_ps(fy, fz)), ex);

                    vpsi[iz][iy][ix] = md_mm_add_ps(vpsi[iz][iy][ix], prod);
                }
            }
        }
    }

    // Write result block to memory
    for (int iz = 0; iz < 8; ++iz) {
        int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
        for (int iy = 0; iy < 8; ++iy) {
            int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
            int index[2] = {
                x_stride[0] + y_stride + z_stride,
                x_stride[1] + y_stride + z_stride,
            };

            md_128 psi[2] = {
                vpsi[iz][iy][0],
                vpsi[iz][iy][1],
            };

            if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                psi[0] = md_mm_mul_ps(psi[0], psi[0]);
                psi[1] = md_mm_mul_ps(psi[1], psi[1]);
            }

            md_mm_storeu_ps(grid_data + index[0], md_mm_add_ps(md_mm_loadu_ps(grid_data + index[0]), psi[0]));
            md_mm_storeu_ps(grid_data + index[1], md_mm_add_ps(md_mm_loadu_ps(grid_data + index[1]), psi[1]));
        }
    }
}

void md_gto_grid_evaluate_sub(float* out_values, const md_grid_t* grid, const int grid_idx_off[3], const int grid_idx_len[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    ASSERT(grid);
    ASSERT(gtos);

    const int* grid_idx_min = grid_idx_off;
    const int  grid_idx_max[3] = {
        grid_idx_off[0] + grid_idx_len[0],
        grid_idx_off[1] + grid_idx_len[1],
        grid_idx_off[2] + grid_idx_len[2],
    };

    //printf("Number of pgtos in volume region: %zu\n", gto.count);

    bool ortho =
        (grid->orientation.elem[0][1] == 0 && grid->orientation.elem[0][2] == 0) &&
        (grid->orientation.elem[1][0] == 0 && grid->orientation.elem[1][2] == 0) &&
        (grid->orientation.elem[2][0] == 0 && grid->orientation.elem[2][1] == 0);

    vec3_t step_x = vec3_mul_f(grid->orientation.col[0], grid->spacing.x);
    vec3_t step_y = vec3_mul_f(grid->orientation.col[1], grid->spacing.y);
    vec3_t step_z = vec3_mul_f(grid->orientation.col[2], grid->spacing.z);

    // There are specialized versions for evaluating 8x8x8 subgrids
    // 8x8x8 Is a good chunk size to operate on as it probably fits in L1 Cache together with the GTOs
    // Then we vectorize over the spatial domain rather than the GTOs to get better register occupation
    if (grid_idx_len[0] == 8 && grid_idx_len[1] == 8 && grid_idx_len[2] == 8) {

#if defined(__AVX512F__) && defined(__AVX512DQ__) || defined (__AVX2__) || defined(__aarch64__) || defined(_M_ARM64)
        // @TODO: Implement real AVX512 path
        if (ortho) {
            evaluate_grid_ortho_8x8x8_256(out_values, grid_idx_min, grid->dim, grid->origin.elem, grid->spacing.elem, gtos, num_gtos, mode);
        } else {
            evaluate_grid_8x8x8_256(out_values, grid_idx_min, grid->dim, grid->origin.elem, step_x.elem, step_y.elem, step_z.elem, gtos, num_gtos, mode);
        }
#elif defined(__SSE2__)
        if (ortho) {
            evaluate_grid_ortho_8x8x8_128(out_values, grid_idx_min, grid->dim, grid->origin.elem, grid->spacing.elem, gtos, num_gtos, mode);
        } else {
            evaluate_grid_8x8x8_128(out_values, grid_idx_min, grid->dim, grid->origin.elem, step_x.elem, step_y.elem, step_z.elem, gtos, num_gtos, mode);
        }
#else
        evaluate_grid_ref(out_values, grid_idx_min, grid_idx_max, grid->dim, grid->origin.elem, step_x.elem, step_y.elem, step_z.elem, gtos, num_gtos, mode);
#endif
    } else {
        // Slowpath
        evaluate_grid_ref(out_values, grid_idx_min, grid_idx_max, grid->dim, grid->origin.elem, step_x.elem, step_y.elem, step_z.elem, gtos, num_gtos, mode);
    }
}

void md_gto_grid_evaluate(float* out_values, const md_grid_t* grid, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    ASSERT(grid);
    ASSERT(gtos);

    int idx_off[3] = {0};
    int idx_len[3] = {0};

    float scl[3] = {
        (grid->orientation.elem[0][0] + grid->orientation.elem[0][1] + grid->orientation.elem[0][2]) * grid->spacing.elem[0],
        (grid->orientation.elem[1][0] + grid->orientation.elem[1][1] + grid->orientation.elem[1][2]) * grid->spacing.elem[1],
        (grid->orientation.elem[2][0] + grid->orientation.elem[2][1] + grid->orientation.elem[2][2]) * grid->spacing.elem[2]
    };

    size_t temp_pos = md_temp_get_pos();
    md_gto_t* sub_gtos = (md_gto_t*)md_temp_push(sizeof(md_gto_t) * num_gtos);

    for (idx_off[2] = 0; idx_off[2] < grid->dim[2]; idx_off[2] += 8) {
        idx_len[2] = MIN(8, grid->dim[2] - idx_off[2]);
        for (idx_off[1] = 0; idx_off[1] < grid->dim[1]; idx_off[1] += 8) {
            idx_len[1] = MIN(8, grid->dim[1] - idx_off[1]);
            for (idx_off[0] = 0; idx_off[0] < grid->dim[0]; idx_off[0] += 8) {
                idx_len[0] = MIN(8, grid->dim[0] - idx_off[0]);

                float aabb_min[3] = {
                    grid->origin.elem[0] + idx_off[0] * scl[0],
                    grid->origin.elem[1] + idx_off[1] * scl[1],
                    grid->origin.elem[2] + idx_off[2] * scl[2],
                };
                float aabb_max[3] = {
                    grid->origin.elem[0] + (idx_off[0] + idx_len[0]) * scl[0],
                    grid->origin.elem[1] + (idx_off[1] + idx_len[1]) * scl[1],
                    grid->origin.elem[2] + (idx_off[2] + idx_len[2]) * scl[2],
                };

                size_t num_sub_gtos = md_gto_aabb_test(sub_gtos, aabb_min, aabb_max, gtos, num_gtos);
                md_gto_grid_evaluate_sub(out_values, grid, idx_off, idx_len, sub_gtos, num_sub_gtos, mode);
            }
        }
    }

    md_temp_set_pos_back(temp_pos);
}

// Evaluate GTOs over a set of passed in packed XYZ coordinates with a bytestride
static void evaluate_gtos(float* out_psi, const float* in_xyz, size_t num_xyz, size_t xyz_stride, const md_gto_t* in_gto, size_t num_gtos, md_gto_eval_mode_t mode) {
    for (size_t j = 0; j < num_xyz; ++j) {
        const float* xyz = (const float*)((const char*)in_xyz + j * xyz_stride);
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];

        double psi = 0.0;
        for (size_t i = 0; i < num_gtos; ++i) {
            double cutoff	= in_gto[i].cutoff;
            double rx		= x - in_gto[i].x;
            double ry		= y - in_gto[i].y;
            double rz		= z - in_gto[i].z;
            double r2		= rx * rx + ry * ry + rz * rz;
            if (r2 > cutoff * cutoff) {
                continue;
            }

            double alpha	= in_gto[i].alpha;
            double coeff	= in_gto[i].coeff;
            int   pi		= in_gto[i].i;
            int   pj		= in_gto[i].j;
            int   pk		= in_gto[i].k;

            double fx = pow(rx, pi);
            double fy = pow(ry, pj);
            double fz = pow(rz, pk);
            double powxyz = fx * fy * fz;
            double exp_term = alpha == 0 ? 1.0 : exp(-alpha * r2);

            double prod = coeff * powxyz * exp_term;
            psi += prod;
        }

        if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
            psi = psi * psi;
        }

        out_psi[j] = (float)psi;
    }
}

void md_gto_xyz_evaluate(float* out_psi, const float* in_xyz, size_t num_xyz, size_t stride, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    if (!out_psi) {
        MD_LOG_ERROR("out_psi array is NULL!");
        return;
    }
    if (!in_xyz) {
        MD_LOG_ERROR("in_xyz base pointer is NULL!");
        return;
    }
    if (!gtos) {
        MD_LOG_ERROR("gtos is NULL!");
        return;
    }
    if (stride != 0 && stride < sizeof(float) * 3) {
        MD_LOG_ERROR("Invalid xyz stride: expected value >= 12 Bytes, got %zu", stride);
        return;
    }

    stride = (stride == 0) ? sizeof(float) * 3 : stride;
    evaluate_gtos(out_psi, in_xyz, num_xyz, stride, gtos, num_gtos, mode);
}

static inline double eval_G(double d, double C, int l, double neg_alpha) {
    return C * fast_pow(d, l) * exp(neg_alpha * d * d);
}

static inline void eval_G_and_G_prime(double* out_G, double* out_G_prime, double d, double C, int l, double neg_alpha) {
    double exp_term = exp(neg_alpha * d * d);
    *out_G		 = C * fast_pow(d, l) * exp_term;
    *out_G_prime = C * fast_pow(d, l-1) * (l + 2 * neg_alpha * d * d) * exp_term;
}

#define PRINT_RESULT 0

static double compute_distance_cutoff(double cutoff_value, int i, int j, int k, int l, double coeff, double alpha) {
    double d = 0.0;

    const double neg_alpha = -alpha;

    // Bake into single constant C
    const double C = fabs(coeff * sqrt((fast_pow(i,i) * fast_pow(j,j) * fast_pow(k,k)) / fast_pow(l,l)));

    // Compute maxima
    const double d_maxima = sqrt(l / (2.0 * fabs(neg_alpha)));

    // Check the contribution at the maxima
    const double y_max = eval_G(d_maxima, C, l, neg_alpha);
    if (y_max < cutoff_value) {
        d = 0.0;
        goto done;
    }

    // If we have an S-type orbital (l == 0) the expression collapses into an expression we can invert and evaluate
    if (l == 0) {
        double y = cutoff_value;
        double a = fabs(coeff) / y;
        double la = log(a);
        d = sqrt(fabs(la) / fabs(neg_alpha));
        goto done;
    }

    // If we end up here we need to perform a numerical search for a d value where the value G(d) < cutoff_value
    // We do not want to overestimate d which will result in a too large radius of influence for the PGTO.
    // And will directly negatively impact performance when evaluating on the grid.
    // Therefore we want to find d where G(d) < cutoff_value but within a tolerance of cutoff_value

    // Search parameters
    const double d_min = d_maxima + 0.001;
    const double d_max = d_maxima + 100.0;
    const double y_tol = cutoff_value * 0.001;
    const double d_tol = 1.0e-9;

    // Initial guess
    // This is the analytical solution for d^2/dx^2 G(d) = 0
    // Which should give us a value which corresponds to the point where we have the maximum negative slope
    d = 0.5 * sqrt(sqrt(neg_alpha*neg_alpha * (8*l + 1)) / (neg_alpha*neg_alpha) + (2*l + 1) / fabs(neg_alpha));

    // Newton-Rhapson iterative search
    for (int iter = 0; iter < 100; ++iter) {
        double y, yp;
        eval_G_and_G_prime(&y, &yp, d, C, l, neg_alpha);

        // Shift function so it intersects the x axis at the point we seek (with a bias towards values less than cutoff value)
        y = y - cutoff_value + y_tol;

        //printf ("d: %.10f, y: %.10f, yp: %.10f\n", d, y, yp);

        if (y < 0 && fabs(y) < y_tol) {
            //printf ("y tolerance met after %i iterations\n", iter);
            break;
        }

        if (fabs(yp) < DBL_EPSILON) {
            //printf ("Denominator is too small!\n");
            break;
        }

        double dn = d - y / yp;
        dn = CLAMP(dn, d_min, d_max);

        if (fabs(dn - d) < d_tol) {
            //printf ("d tolerance met after %i iterations\n", iter);
            break;
        }

        d = dn;
    }

done:
#if PRINT_RESULT
    if (d > 0.0) {
        printf("Cutoff dist and value: %15.5f, %15.12f\n", d, eval_G(d, C, l, neg_alpha));
    }
#endif
    return d;
}

double md_gto_compute_radius_of_influence(int i, int j, int k, double coeff, double alpha, double cutoff) {
    int l = i + j + k;
    return compute_distance_cutoff(cutoff, i, j, k, l, coeff, alpha);
}

size_t md_gto_cutoff_compute_and_filter(md_gto_t* gtos, size_t count, double value) {
    if (value == 0) {
        for (size_t i = 0; i < count; ++i) {
            gtos[i].cutoff = FLT_MAX;
        }
    } else {
        for (size_t i = 0; i < count;) {
            gtos[i].cutoff = (float)compute_distance_cutoff(value, gtos[i].i, gtos[i].j, gtos[i].k, gtos[i].l, gtos[i].coeff, gtos[i].alpha);
            if (gtos[i].cutoff == 0.0f) {
                gtos[i] = gtos[--count];
            } else {
                ++i;
            }
        }
    }
    return count;
}

size_t md_gto_aabb_test(md_gto_t* out_gtos, const float aabb_min[3], const float aabb_max[3], const md_gto_t* in_gtos, size_t in_num_gtos) {
    // Extract a subset of gtos that overlap with the evaluated subportion of the grid
    // @TODO: This can be vectorized, Let us pray to the compiler gods for now
    size_t num_gtos = 0;
    for (size_t i = 0; i < in_num_gtos; ++i) {
        float x  = in_gtos[i].x;
        float y  = in_gtos[i].y;
        float z  = in_gtos[i].z;
        float cutoff = in_gtos[i].cutoff;

        float cx = CLAMP(x, aabb_min[0], aabb_max[0]);
        float cy = CLAMP(y, aabb_min[1], aabb_max[1]);
        float cz = CLAMP(z, aabb_min[2], aabb_max[2]);

        float dx = x - cx;
        float dy = y - cy;
        float dz = z - cz;

        float d2 = dx * dx + dy * dy + dz * dz;

        if (d2 > cutoff * cutoff) {
            continue;
        }
        out_gtos[num_gtos++] = in_gtos[i];
    }
    return num_gtos;
}
