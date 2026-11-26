#include <md_gto_dmat.h>

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_executor.h>
#include <core/md_log.h>

#include <math.h>
#include <stdlib.h>

#define DMAT_T 16u

// The shells' sort keys hold the shell in their low bits: no state beside the arrays
static int dmat_cmp_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return x < y ? -1 : x > y;
}

static void* dmat_alloc(md_allocator_i* alloc, size_t bytes) {
    return md_alloc(alloc, MD_MAX(bytes, 1));
}

static void dmat_free_array(md_allocator_i* alloc, void* p, size_t bytes) {
    if (p) md_free(alloc, p, MD_MAX(bytes, 1));
}

void md_gto_dmat_free(md_gto_dmat_t* dm) {
    if (!dm || !dm->alloc) return;
    md_allocator_i* a = dm->alloc;
    dmat_free_array(a, dm->xyz,       sizeof(float) * 3 * dm->num_atoms);
    dmat_free_array(a, dm->shell_r,   sizeof(float) * dm->num_shells);
    dmat_free_array(a, dm->perm,      sizeof(int32_t) * DMAT_T * dm->nt);
    dmat_free_array(a, dm->tao,       sizeof(uint32_t) * dm->num_ao);
    dmat_free_array(a, dm->shell_tao, sizeof(uint32_t) * dm->num_shells);
    dmat_free_array(a, dm->row_ptr,   sizeof(uint64_t) * (dm->nt + 1));
    dmat_free_array(a, dm->cols,      sizeof(uint32_t) * dm->num_tiles);
    dmat_free_array(a, dm->tiles,     sizeof(double) * 256 * dm->num_tiles);
    md_memset(dm, 0, sizeof(*dm));
}

bool md_gto_dmat_init(md_gto_dmat_t* dm, const md_gto_dmat_desc_t* desc, md_allocator_i* alloc) {
    if (!dm) return false;
    md_memset(dm, 0, sizeof(*dm));
    if (!desc || !desc->basis || !desc->atom_xyz || !alloc || !(desc->cutoff > 0.0) || !(desc->margin >= 0.0)) {
        MD_LOG_ERROR("md_gto_dmat_init: invalid arguments (a basis, positions, an allocator, a cutoff > 0, a margin >= 0)");
        return false;
    }
    const md_gto_basis_t* basis = desc->basis;
    const uint32_t ns = basis->num_shells;
    const uint32_t num_ao = (uint32_t)md_gto_basis_num_ao(basis);
    const uint32_t num_atoms = (uint32_t)md_gto_basis_num_atoms(basis);
    if (ns == 0 || num_ao == 0) {
        MD_LOG_ERROR("md_gto_dmat_init: empty basis");
        return false;
    }
    if (ns > 0xFFFFFu || num_atoms > 0xFFFFFu) {
        MD_LOG_ERROR("md_gto_dmat_init: more than 2^20 shells or atoms");
        return false;
    }
    const uint32_t nt = (num_ao + DMAT_T - 1) / DMAT_T;
    dm->alloc = alloc;
    dm->num_ao = num_ao;
    dm->num_atoms = num_atoms;
    dm->num_shells = ns;
    dm->nt = nt;
    dm->cutoff = desc->cutoff;
    dm->margin = desc->margin;
    dm->xyz       = dmat_alloc(alloc, sizeof(float) * 3 * num_atoms);
    dm->shell_r   = dmat_alloc(alloc, sizeof(float) * ns);
    dm->perm      = dmat_alloc(alloc, sizeof(int32_t) * DMAT_T * nt);
    dm->tao       = dmat_alloc(alloc, sizeof(uint32_t) * num_ao);
    dm->shell_tao = dmat_alloc(alloc, sizeof(uint32_t) * ns);
    dm->row_ptr   = dmat_alloc(alloc, sizeof(uint64_t) * (nt + 1));
    md_allocator_i* heap = md_get_heap_allocator();
    uint64_t* key  = dmat_alloc(heap, sizeof(uint64_t) * ns);
    uint32_t* ao0  = dmat_alloc(heap, sizeof(uint32_t) * ns);       // each shell's first AO in the basis order
    float*    ra   = dmat_alloc(heap, sizeof(float) * num_atoms);   // each atom's largest radius
    uint32_t* afst = dmat_alloc(heap, sizeof(uint32_t) * (num_atoms + 1));
    uint32_t* ash  = dmat_alloc(heap, sizeof(uint32_t) * ns);
    uint8_t*  on   = dmat_alloc(heap, (size_t)nt * nt);
    bool ok = dm->xyz && dm->shell_r && dm->perm && dm->tao && dm->shell_tao && dm->row_ptr && key && ao0 && ra && afst && ash && on;
    if (ok) {
        const size_t stride = desc->atom_xyz_stride ? desc->atom_xyz_stride : sizeof(float) * 3;
        for (uint32_t a = 0; a < num_atoms; ++a) {
            const float* p = (const float*)((const uint8_t*)desc->atom_xyz + a * stride);
            dm->xyz[3 * a + 0] = p[0]; dm->xyz[3 * a + 1] = p[1]; dm->xyz[3 * a + 2] = p[2];
        }
        md_gto_basis_shell_radii(dm->shell_r, basis, desc->cutoff);

        // The tiled order: radius class (2 bohr), then atom, then shell; a shell's AOs stay together
        uint32_t acc = 0;
        for (uint32_t s = 0; s < ns; ++s) {
            ao0[s] = acc;
            acc += md_gto_num_cart_ao(basis->shells[s].l);
            key[s] = ((uint64_t)(dm->shell_r[s] * 0.5f) << 40) | ((uint64_t)basis->shells[s].atom_idx << 20) | s;
        }
        qsort(key, ns, sizeof(uint64_t), dmat_cmp_u64);
        for (uint32_t i = 0; i < DMAT_T * nt; ++i) dm->perm[i] = -1;
        acc = 0;
        for (uint32_t i = 0; i < ns; ++i) {
            const uint32_t s = (uint32_t)(key[i] & 0xFFFFFu);
            const uint32_t nc = md_gto_num_cart_ao(basis->shells[s].l);
            dm->shell_tao[s] = acc;
            for (uint32_t c = 0; c < nc; ++c) {
                dm->perm[acc + c] = (int32_t)(ao0[s] + c);
                dm->tao[ao0[s] + c] = acc + c;
            }
            acc += nc;
        }

        // The tile pairs of the shell pairs that meet (centres no further apart than their radii plus margin, both
        // radii > 0), visited atom pair by atom pair: an atom pair further apart than its largest radii plus margin
        // has no shell pair that meets
        md_memset(ra, 0, sizeof(float) * num_atoms);
        md_memset(afst, 0, sizeof(uint32_t) * (num_atoms + 1));
        for (uint32_t s = 0; s < ns; ++s) {
            const uint32_t a = basis->shells[s].atom_idx;
            ra[a] = MD_MAX(ra[a], dm->shell_r[s]);
            afst[a + 1]++;
        }
        for (uint32_t a = 0; a < num_atoms; ++a) afst[a + 1] += afst[a];
        for (uint32_t s = 0; s < ns; ++s) ash[afst[basis->shells[s].atom_idx]++] = s;   // afst[a] now the end of a
        for (uint32_t a = num_atoms; a > 0; --a) afst[a] = afst[a - 1];
        afst[0] = 0;
        md_memset(on, 0, (size_t)nt * nt);
        for (uint32_t a = 0; a < num_atoms; ++a) {
            const float* p = dm->xyz + 3 * a;
            for (uint32_t b = a; b < num_atoms; ++b) {
                const float* q = dm->xyz + 3 * b;
                const double dx = p[0] - q[0], dy = p[1] - q[1], dz = p[2] - q[2];
                const double d2 = dx * dx + dy * dy + dz * dz;
                const double am = (double)ra[a] + (double)ra[b] + desc->margin;
                if (d2 > am * am) continue;
                for (uint32_t x = afst[a]; x < afst[a + 1]; ++x) {
                    const uint32_t s = ash[x];
                    if (!(dm->shell_r[s] > 0.0f)) continue;
                    const uint32_t s0 = dm->shell_tao[s] / DMAT_T, s1 = (dm->shell_tao[s] + md_gto_num_cart_ao(basis->shells[s].l) - 1) / DMAT_T;
                    for (uint32_t y = (a == b ? x : afst[b]); y < afst[b + 1]; ++y) {
                        const uint32_t u = ash[y];
                        if (!(dm->shell_r[u] > 0.0f)) continue;
                        const double lim = (double)dm->shell_r[s] + (double)dm->shell_r[u] + desc->margin;
                        if (d2 > lim * lim) continue;
                        const uint32_t u0 = dm->shell_tao[u] / DMAT_T, u1 = (dm->shell_tao[u] + md_gto_num_cart_ao(basis->shells[u].l) - 1) / DMAT_T;
                        for (uint32_t I = s0; I <= s1; ++I) {
                            for (uint32_t J = u0; J <= u1; ++J) on[(size_t)MD_MIN(I, J) * nt + MD_MAX(I, J)] = 1;
                        }
                    }
                }
            }
        }
        dm->row_ptr[0] = 0;
        for (uint32_t I = 0; I < nt; ++I) {
            uint64_t c = 0;
            for (uint32_t J = I; J < nt; ++J) c += on[(size_t)I * nt + J];
            dm->row_ptr[I + 1] = dm->row_ptr[I] + c;
        }
        dm->num_tiles = dm->row_ptr[nt];
        dm->cols = dmat_alloc(alloc, sizeof(uint32_t) * dm->num_tiles);
        ok = dm->cols != NULL;
        if (ok) {
            uint64_t e = 0;
            for (uint32_t I = 0; I < nt; ++I) {
                for (uint32_t J = I; J < nt; ++J) if (on[(size_t)I * nt + J]) dm->cols[e++] = J;
            }
        }
    }
    dmat_free_array(heap, key,  sizeof(uint64_t) * ns);
    dmat_free_array(heap, ao0,  sizeof(uint32_t) * ns);
    dmat_free_array(heap, ra,   sizeof(float) * num_atoms);
    dmat_free_array(heap, afst, sizeof(uint32_t) * (num_atoms + 1));
    dmat_free_array(heap, ash,  sizeof(uint32_t) * ns);
    dmat_free_array(heap, on,   (size_t)nt * nt);
    if (!ok) {
        MD_LOG_ERROR("md_gto_dmat_init: out of memory");
        md_gto_dmat_free(dm);
        return false;
    }
    return true;
}

bool md_gto_dmat_pattern_valid(const md_gto_dmat_t* dm, const float* atom_xyz, size_t atom_xyz_stride) {
    if (!dm || !dm->xyz || !atom_xyz) return false;
    const size_t stride = atom_xyz_stride ? atom_xyz_stride : sizeof(float) * 3;
    const double lim = 0.25 * dm->margin * dm->margin;
    for (uint32_t a = 0; a < dm->num_atoms; ++a) {
        const float* p = (const float*)((const uint8_t*)atom_xyz + a * stride);
        const double dx = p[0] - dm->xyz[3 * a + 0], dy = p[1] - dm->xyz[3 * a + 1], dz = p[2] - dm->xyz[3 * a + 2];
        if (dx * dx + dy * dy + dz * dz > lim) return false;
    }
    return true;
}

// ---------------------------------------------------------------------------------------------------- the build

// Orbitals a chunk at a time: P[(J kc + k) 16 + c] = orbital k's coefficient of tiled AO 16 J + c, so that a tile's
// 16 x kc block of a chunk is contiguous
#define DMAT_CHUNK 128u

typedef struct dmat_build_t {
    md_gto_dmat_t* dm;
    const double*  P;
    const double*  w;       // the chunk's weights
    uint32_t       kc;      // the chunk's orbitals
} dmat_build_t;

static void dmat_build_rows(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const dmat_build_t* b = (const dmat_build_t*)arg;
    const md_gto_dmat_t* dm = b->dm;
    const uint32_t kc = b->kc;
    for (uint32_t I = beg; I < end; ++I) {
        const double* PI = b->P + (size_t)I * kc * DMAT_T;
        for (uint64_t e = dm->row_ptr[I]; e < dm->row_ptr[I + 1]; ++e) {
            const double* PJ = b->P + (size_t)dm->cols[e] * kc * DMAT_T;
            double* tile = dm->tiles + 256 * e;
            for (uint32_t a = 0; a < DMAT_T; a += 2) {
                double t0[DMAT_T], t1[DMAT_T];
                for (uint32_t c = 0; c < DMAT_T; ++c) { t0[c] = tile[a * DMAT_T + c]; t1[c] = tile[(a + 1) * DMAT_T + c]; }
                for (uint32_t k = 0; k < kc; ++k) {
                    const double q0 = b->w[k] * PI[k * DMAT_T + a], q1 = b->w[k] * PI[k * DMAT_T + a + 1];
                    const double* pj = PJ + k * DMAT_T;
                    for (uint32_t c = 0; c < DMAT_T; ++c) { t0[c] += q0 * pj[c]; t1[c] += q1 * pj[c]; }
                }
                for (uint32_t c = 0; c < DMAT_T; ++c) { tile[a * DMAT_T + c] = t0[c]; tile[(a + 1) * DMAT_T + c] = t1[c]; }
            }
        }
    }
}

bool md_gto_dmat_build(md_gto_dmat_t* dm, const double* coeff, size_t mo_stride, size_t ao_stride, const double* weight, size_t num_mo,
                       const md_executor_t* executor) {
    if (!dm || !dm->perm || !coeff || !weight) {
        MD_LOG_ERROR("md_gto_dmat_build: invalid arguments");
        return false;
    }
    if (!dm->tiles) {
        dm->tiles = dmat_alloc(dm->alloc, sizeof(double) * 256 * dm->num_tiles);
        if (!dm->tiles) {
            MD_LOG_ERROR("md_gto_dmat_build: out of memory (%llu tiles)", (unsigned long long)dm->num_tiles);
            return false;
        }
    }
    md_memset(dm->tiles, 0, sizeof(double) * 256 * dm->num_tiles);
    md_allocator_i* heap = md_get_heap_allocator();
    const size_t ntao = (size_t)DMAT_T * dm->nt;
    double*   P    = dmat_alloc(heap, sizeof(double) * ntao * DMAT_CHUNK);
    double*   w    = dmat_alloc(heap, sizeof(double) * DMAT_CHUNK);
    uint32_t* kept = dmat_alloc(heap, sizeof(uint32_t) * DMAT_CHUNK);
    bool ok = P && w && kept;
    size_t k = 0;
    while (ok && k < num_mo) {
        uint32_t kc = 0;
        for (; k < num_mo && kc < DMAT_CHUNK; ++k) {
            if (weight[k] != 0.0) { kept[kc] = (uint32_t)k; w[kc] = weight[k]; ++kc; }
        }
        if (kc == 0) break;
        for (size_t t = 0; t < ntao; ++t) {
            const int32_t ao = dm->perm[t];
            double* dst = P + ((t / DMAT_T) * kc) * DMAT_T + (t % DMAT_T);
            for (uint32_t j = 0; j < kc; ++j) dst[j * DMAT_T] = ao < 0 ? 0.0 : coeff[(size_t)kept[j] * mo_stride + (size_t)ao * ao_stride];
        }
        dmat_build_t b = { dm, P, w, kc };
        ok = md_parallel_for(executor, dm->nt, 1, dmat_build_rows, &b);
    }
    dmat_free_array(heap, P, sizeof(double) * ntao * DMAT_CHUNK);
    dmat_free_array(heap, w, sizeof(double) * DMAT_CHUNK);
    dmat_free_array(heap, kept, sizeof(uint32_t) * DMAT_CHUNK);
    if (!ok) {
        MD_LOG_ERROR("md_gto_dmat_build: out of memory or cancelled");
        return false;
    }
    // A diagonal tile exactly symmetric, as D is: its upper triangle mirrored
    for (uint32_t I = 0; I < dm->nt; ++I) {
        if (dm->row_ptr[I] == dm->row_ptr[I + 1] || dm->cols[dm->row_ptr[I]] != I) continue;
        double* tile = dm->tiles + 256 * dm->row_ptr[I];
        for (uint32_t a = 0; a < DMAT_T; ++a) for (uint32_t c = a + 1; c < DMAT_T; ++c) tile[c * DMAT_T + a] = tile[a * DMAT_T + c];
    }
    return true;
}

int64_t md_gto_dmat_find(const md_gto_dmat_t* dm, uint32_t I, uint32_t J) {
    if (!dm || !dm->row_ptr || I > J || J >= dm->nt) return -1;
    uint64_t lo = dm->row_ptr[I], hi = dm->row_ptr[I + 1];
    while (lo < hi) {
        const uint64_t mid = lo + (hi - lo) / 2;
        if (dm->cols[mid] < J) lo = mid + 1; else hi = mid;
    }
    return (lo < dm->row_ptr[I + 1] && dm->cols[lo] == J) ? (int64_t)lo : -1;
}

double md_gto_dmat_at(const md_gto_dmat_t* dm, uint32_t i, uint32_t j) {
    if (!dm || !dm->tiles || i >= dm->num_ao || j >= dm->num_ao) return 0.0;
    uint32_t ti = dm->tao[i], tj = dm->tao[j];
    if (ti / DMAT_T > tj / DMAT_T) { const uint32_t x = ti; ti = tj; tj = x; }
    const int64_t e = md_gto_dmat_find(dm, ti / DMAT_T, tj / DMAT_T);
    return e < 0 ? 0.0 : dm->tiles[256 * (size_t)e + (ti % DMAT_T) * DMAT_T + (tj % DMAT_T)];
}

bool md_gto_dmat_shell_block(double* out, const md_gto_dmat_t* dm, uint32_t sa, uint32_t na, uint32_t sb, uint32_t nb) {
    if (!out) return false;
    md_memset(out, 0, sizeof(double) * na * nb);
    if (!dm || !dm->tiles || sa >= dm->num_shells || sb >= dm->num_shells) return false;
    const uint32_t ta = dm->shell_tao[sa], tb = dm->shell_tao[sb];
    bool any = false;
    // the at most 2 x 2 tile pairs the block spans, each looked up once
    for (uint32_t I = ta / DMAT_T; I <= (ta + na - 1) / DMAT_T; ++I) {
        for (uint32_t J = tb / DMAT_T; J <= (tb + nb - 1) / DMAT_T; ++J) {
            const bool swap = I > J;
            const int64_t e = swap ? md_gto_dmat_find(dm, J, I) : md_gto_dmat_find(dm, I, J);
            if (e < 0) continue;
            any = true;
            const double* tile = dm->tiles + 256 * (size_t)e;
            for (uint32_t ca = 0; ca < na; ++ca) {
                const uint32_t ti = ta + ca;
                if (ti / DMAT_T != I) continue;
                for (uint32_t cb = 0; cb < nb; ++cb) {
                    const uint32_t tj = tb + cb;
                    if (tj / DMAT_T != J) continue;
                    out[ca * nb + cb] = swap ? tile[(tj % DMAT_T) * DMAT_T + ti % DMAT_T] : tile[(ti % DMAT_T) * DMAT_T + tj % DMAT_T];
                }
            }
        }
    }
    return any;
}

// ---------------------------------------------------------------------------------------------------- at points

typedef struct dmat_mono_t {
    uint8_t i, j, k;
    double  f;              // md_gto_cart_norm_factor
} dmat_mono_t;

typedef struct dmat_eval_t {
    const md_gto_dmat_t*  dm;
    const md_gto_basis_t* basis;
    const float*          atom_xyz;
    size_t                atom_stride;
    const float*          xyz;
    size_t                xyz_stride;
    double*               out_rho;
    double*               out_grad;
    double*               out_tau;
    bool                  deriv;
    dmat_mono_t           mono[MD_GTO_MAX_ANGULAR_MOMENTUM + 1][15];
    // per thread_idx: AO values [16 nt], gradients [3 16 nt], the tiles a point touches and their stamps
    double**              phi;
    double**              dphi;
    uint32_t**            touched;
    uint32_t**            stamp;
} dmat_eval_t;

static void dmat_eval_points(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    const dmat_eval_t* E = (const dmat_eval_t*)arg;
    const md_gto_dmat_t* dm = E->dm;
    const md_gto_basis_t* basis = E->basis;
    double* phi = E->phi[thread_idx];
    double* dphi = E->dphi[thread_idx];
    uint32_t* touched = E->touched[thread_idx];
    uint32_t* stamp = E->stamp[thread_idx];
    for (uint32_t pi = beg; pi < end; ++pi) {
        const float* pf = (const float*)((const uint8_t*)E->xyz + (size_t)pi * E->xyz_stride);
        const uint32_t tag = pi + 1;     // stamps hold the point + 1; the arrays start at 0
        uint32_t nt_touched = 0;
        for (uint32_t s = 0; s < dm->num_shells; ++s) {
            const md_gto_shell_t* sh = &basis->shells[s];
            const float* a = (const float*)((const uint8_t*)E->atom_xyz + (size_t)sh->atom_idx * E->atom_stride);
            const double d[3] = { (double)pf[0] - a[0], (double)pf[1] - a[1], (double)pf[2] - a[2] };
            const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            const double rs = dm->shell_r[s];
            if (!(r2 < rs * rs)) continue;
            double R = 0.0, dR = 0.0;      // sum c e^(-a r^2), sum -2 a c e^(-a r^2)
            for (uint32_t p = 0; p < sh->num_primitives; ++p) {
                const double al = basis->alpha[sh->primitive_offset + p];
                const double e = (double)basis->coeff[sh->primitive_offset + p] * exp(-al * r2);
                R += e;
                dR -= 2.0 * al * e;
            }
            double pw[3][MD_GTO_MAX_ANGULAR_MOMENTUM + 2];
            for (int q = 0; q < 3; ++q) {
                pw[q][0] = 1.0;
                for (uint32_t n = 1; n <= sh->l + 1; ++n) pw[q][n] = pw[q][n - 1] * d[q];
            }
            const uint32_t t0 = dm->shell_tao[s], nc = md_gto_num_cart_ao(sh->l);
            for (uint32_t c = 0; c < nc; ++c) {
                const dmat_mono_t* m = &E->mono[sh->l][c];
                const double mono = pw[0][m->i] * pw[1][m->j] * pw[2][m->k];
                phi[t0 + c] = m->f * mono * R;
                if (E->deriv) {
                    const int e3[3] = { m->i, m->j, m->k };
                    for (int q = 0; q < 3; ++q) {
                        const double o = pw[(q + 1) % 3][e3[(q + 1) % 3]] * pw[(q + 2) % 3][e3[(q + 2) % 3]];
                        const double dm_q = e3[q] > 0 ? e3[q] * pw[q][e3[q] - 1] * o : 0.0;
                        dphi[3 * (t0 + c) + q] = m->f * (dm_q * R + mono * dR * d[q]);
                    }
                }
            }
            for (uint32_t I = t0 / DMAT_T; I <= (t0 + nc - 1) / DMAT_T; ++I) {
                if (stamp[I] != tag) { stamp[I] = tag; touched[nt_touched++] = I; }
            }
        }
        double rho = 0.0, g[3] = { 0.0, 0.0, 0.0 }, tau = 0.0;
        for (uint32_t x = 0; x < nt_touched; ++x) {
            const uint32_t I = touched[x];
            const double* fI = phi + (size_t)I * DMAT_T;
            const double* gI = dphi + (size_t)I * DMAT_T * 3;
            for (uint64_t e = dm->row_ptr[I]; e < dm->row_ptr[I + 1]; ++e) {
                const uint32_t J = dm->cols[e];
                if (stamp[J] != tag) continue;
                const double s2 = J == I ? 1.0 : 2.0;
                const double* tile = dm->tiles + 256 * e;
                const double* fJ = phi + (size_t)J * DMAT_T;
                const double* gJ = dphi + (size_t)J * DMAT_T * 3;
                for (uint32_t a = 0; a < DMAT_T; ++a) {
                    const double* row = tile + a * DMAT_T;
                    double v = 0.0;
                    for (uint32_t c = 0; c < DMAT_T; ++c) v += row[c] * fJ[c];
                    rho += s2 * fI[a] * v;
                    if (!E->deriv) continue;
                    double vg[3] = { 0.0, 0.0, 0.0 };
                    for (uint32_t c = 0; c < DMAT_T; ++c) {
                        vg[0] += row[c] * gJ[3 * c + 0];
                        vg[1] += row[c] * gJ[3 * c + 1];
                        vg[2] += row[c] * gJ[3 * c + 2];
                    }
                    for (int q = 0; q < 3; ++q) g[q] += s2 * (gI[3 * a + q] * v + fI[a] * vg[q]);
                    tau += 0.5 * s2 * (gI[3 * a + 0] * vg[0] + gI[3 * a + 1] * vg[1] + gI[3 * a + 2] * vg[2]);
                }
            }
        }
        for (uint32_t x = 0; x < nt_touched; ++x) {
            const uint32_t I = touched[x];
            md_memset(phi + (size_t)I * DMAT_T, 0, sizeof(double) * DMAT_T);
            if (E->deriv) md_memset(dphi + (size_t)I * DMAT_T * 3, 0, sizeof(double) * DMAT_T * 3);
        }
        if (E->out_rho) E->out_rho[pi] = rho;
        if (E->out_grad) { E->out_grad[3 * (size_t)pi + 0] = g[0]; E->out_grad[3 * (size_t)pi + 1] = g[1]; E->out_grad[3 * (size_t)pi + 2] = g[2]; }
        if (E->out_tau) E->out_tau[pi] = tau;
    }
}

size_t md_gto_dmat_evaluate(double* out_rho, double* out_grad, double* out_tau, const md_gto_dmat_t* dm, const md_gto_basis_t* basis,
                            const float* atom_xyz, size_t atom_xyz_stride, const float* xyz, size_t count, size_t xyz_stride,
                            const md_executor_t* executor) {
    if (!dm || !dm->tiles || !basis || basis->num_shells != dm->num_shells || !atom_xyz || !xyz || count > UINT32_MAX - 1) {
        MD_LOG_ERROR("md_gto_dmat_evaluate: invalid arguments (a built dmat, its basis, positions, points)");
        return 0;
    }
    if (count == 0) return 0;
    for (uint32_t s = 0; s < basis->num_shells; ++s) {
        if (basis->shells[s].l > MD_GTO_MAX_ANGULAR_MOMENTUM) {
            MD_LOG_ERROR("md_gto_dmat_evaluate: a shell of angular momentum %u", basis->shells[s].l);
            return 0;
        }
    }
    dmat_eval_t E = {
        .dm = dm, .basis = basis, .atom_xyz = atom_xyz, .atom_stride = atom_xyz_stride ? atom_xyz_stride : sizeof(float) * 3,
        .xyz = xyz, .xyz_stride = xyz_stride ? xyz_stride : sizeof(float) * 3,
        .out_rho = out_rho, .out_grad = out_grad, .out_tau = out_tau, .deriv = out_grad || out_tau,
    };
    for (uint32_t l = 0; l <= MD_GTO_MAX_ANGULAR_MOMENTUM; ++l) {
        for (uint32_t c = 0; c < md_gto_num_cart_ao(l); ++c) {
            int i, j, k;
            md_gto_cart_ijk(&i, &j, &k, l, c);
            E.mono[l][c] = (dmat_mono_t){ (uint8_t)i, (uint8_t)j, (uint8_t)k, md_gto_cart_norm_factor(i, j, k) };
        }
    }
    const uint32_t nth = md_executor_num_threads(executor);
    const size_t ntao = (size_t)DMAT_T * dm->nt;
    md_allocator_i* heap = md_get_heap_allocator();
    double*   phi_all = dmat_alloc(heap, sizeof(double) * ntao * nth);
    double*   dph_all = E.deriv ? dmat_alloc(heap, sizeof(double) * 3 * ntao * nth) : NULL;
    uint32_t* tch_all = dmat_alloc(heap, sizeof(uint32_t) * dm->nt * nth);
    uint32_t* stm_all = dmat_alloc(heap, sizeof(uint32_t) * dm->nt * nth);
    void**    ptrs    = dmat_alloc(heap, sizeof(void*) * 4 * nth);
    size_t result = 0;
    if (phi_all && (dph_all || !E.deriv) && tch_all && stm_all && ptrs) {
        md_memset(phi_all, 0, sizeof(double) * ntao * nth);
        if (dph_all) md_memset(dph_all, 0, sizeof(double) * 3 * ntao * nth);
        md_memset(stm_all, 0, sizeof(uint32_t) * dm->nt * nth);
        E.phi = (double**)ptrs;
        E.dphi = (double**)ptrs + nth;
        E.touched = (uint32_t**)(ptrs + 2 * nth);
        E.stamp = (uint32_t**)(ptrs + 3 * nth);
        for (uint32_t t = 0; t < nth; ++t) {
            E.phi[t] = phi_all + ntao * t;
            E.dphi[t] = dph_all ? dph_all + 3 * ntao * t : NULL;
            E.touched[t] = tch_all + (size_t)dm->nt * t;
            E.stamp[t] = stm_all + (size_t)dm->nt * t;
        }
        result = md_parallel_for(executor, (uint32_t)count, 64, dmat_eval_points, &E) ? count : 0;
    } else {
        MD_LOG_ERROR("md_gto_dmat_evaluate: out of memory");
    }
    dmat_free_array(heap, phi_all, sizeof(double) * ntao * nth);
    dmat_free_array(heap, dph_all, sizeof(double) * 3 * ntao * nth);
    dmat_free_array(heap, tch_all, sizeof(uint32_t) * dm->nt * nth);
    dmat_free_array(heap, stm_all, sizeof(uint32_t) * dm->nt * nth);
    dmat_free_array(heap, ptrs, sizeof(void*) * 4 * nth);
    return result;
}
