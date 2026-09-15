#include <md_qm.h>

#include <md_system.h>
#include <md_gto.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_unit.h>
#include <core/md_parse.h>

#include <math.h>
#include <string.h>

// ---------------------------------------------------------------------------
// Units
// ---------------------------------------------------------------------------

md_unit_t md_qm_unit_wavenumber(void) { return md_unit_pow(md_unit_scl(md_unit_meter(), 1.0e-2), -1); }
md_unit_t md_qm_unit_km_per_mol(void) { return md_unit_div(md_unit_scl(md_unit_meter(), 1.0e3), md_unit_mole()); }
md_unit_t md_qm_unit_amu(void)        { return md_unit_scl(md_unit_kilogram(), 1.66053906660e-27); }

// ---------------------------------------------------------------------------
// Publishing vocabulary
// ---------------------------------------------------------------------------

md_attribute_id_t md_qm_publish(md_system_t* sys, str_t path, str_t label, md_unit_t unit,
                                md_attribute_format_t format, const void* data, size_t byte_size) {
    ASSERT(sys);
    return md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
        .path      = path,
        .format    = format,
        .unit      = unit,
        .label     = label,
        .data      = data,
        .byte_size = byte_size,
    });
}

md_attribute_id_t md_qm_publish_virtual(md_system_t* sys, str_t path, str_t label, md_unit_t unit,
                                        md_attribute_format_t format, const md_attribute_virtual_t* virt) {
    ASSERT(sys);
    return md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
        .path   = path,
        .format = format,
        .unit   = unit,
        .label  = label,
        .virt   = virt,
    });
}

md_attribute_id_t md_qm_publish_scalar(md_system_t* sys, str_t path, str_t label, md_unit_t unit, double value) {
    // The value is copied, so a local is fine.
    md_attribute_format_t format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 0 };
    return md_qm_publish(sys, path, label, unit, format, &value, sizeof(double));
}

md_attribute_id_t md_qm_publish_series(md_system_t* sys, str_t path, str_t label, md_unit_t unit, const double* values, size_t count) {
    if (!values || count == 0) {
        return MD_ATTRIBUTE_INVALID;
    }
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)count },
    };
    return md_qm_publish(sys, path, label, unit, format, values, count * sizeof(double));
}

md_attribute_id_t md_qm_publish_vec3_series(md_system_t* sys, str_t path, str_t label, md_unit_t unit, const dvec3_t* values, size_t count) {
    if (!values || count == 0) {
        return MD_ATTRIBUTE_INVALID;
    }
    // dvec3_t is three contiguous doubles, so the source array is already the interleaved layout an
    // attribute stores and this is a straight copy.
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 1, .shape = { (uint32_t)count },
    };
    return md_qm_publish(sys, path, label, unit, format, values, count * 3 * sizeof(double));
}

md_attribute_id_t md_qm_publish_matrix(md_system_t* sys, str_t path, str_t label, md_unit_t unit, const double* values, size_t rows, size_t cols) {
    if (!values || rows == 0 || cols == 0) {
        return MD_ATTRIBUTE_INVALID;
    }
    // Row major with the last index fastest, which is what every 2D quantity in the tree is.
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2, .shape = { (uint32_t)rows, (uint32_t)cols },
    };
    return md_qm_publish(sys, path, label, unit, format, values, rows * cols * sizeof(double));
}

md_attribute_id_t md_qm_publish_str(md_system_t* sys, str_t path, str_t label, str_t value) {
    if (str_empty(value)) {
        return MD_ATTRIBUTE_INVALID;
    }
    // A single string is rank 1 {1}, by the same rule that makes a single 3-vector rank 1 {1} of 3
    // components. The descriptor carries the TEXT and the table stores a handle - see STRINGS in
    // md_system.h.
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 1, .shape = { 1 },
    };
    return md_qm_publish(sys, path, label, md_unit_none(), format, &value, sizeof(str_t));
}

md_attribute_id_t md_qm_publish_strings(md_system_t* sys, str_t path, str_t label, const str_t* values, size_t count) {
    if (!values || count == 0) {
        return MD_ATTRIBUTE_INVALID;
    }
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 1, .shape = { (uint32_t)count },
    };
    return md_qm_publish(sys, path, label, md_unit_none(), format, values, count * sizeof(str_t));
}

md_attribute_id_t md_qm_publish_column(md_system_t* sys, str_t path, str_t label, md_unit_t unit, md_attribute_type_t type,
                                       const void* base, size_t stride, size_t count) {
    ASSERT(sys);
    if (!base || count == 0) {
        return MD_ATTRIBUTE_INVALID;
    }
    md_attribute_format_t format = {
        .type = type, .components = 1, .rank = 1, .shape = { (uint32_t)count },
    };
    // The source is strided and an attribute is contiguous, so the values are gathered into the
    // table's own storage rather than through a temporary which is then copied again.
    md_attribute_id_t id = md_qm_publish(sys, path, label, unit, format, NULL, 0);
    if (id == MD_ATTRIBUTE_INVALID) {
        return MD_ATTRIBUTE_INVALID;
    }
    uint8_t* dst = (uint8_t*)md_attributes_data(&sys->attributes, id, type);
    if (!dst) {
        md_attributes_remove(&sys->attributes, id);
        return MD_ATTRIBUTE_INVALID;
    }
    const size_t   elem_size = md_attribute_type_size(type);
    const uint8_t* src       = (const uint8_t*)base;
    for (size_t i = 0; i < count; ++i) {
        MEMCPY(dst + i * elem_size, src + i * stride, elem_size);
    }
    return id;
}

md_attribute_id_t md_qm_publish_origin(md_system_t* sys, str_t path, dvec3_t origin) {
    md_attribute_format_t format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 0 };
    return md_qm_publish(sys, path, (str_t){0}, md_unit_angstrom(), format, &origin, 3 * sizeof(double));
}

md_attribute_id_t md_qm_alias(md_system_t* sys, md_attribute_id_t target, str_t path) {
    ASSERT(sys);
    if (target == MD_ATTRIBUTE_INVALID) {
        return MD_ATTRIBUTE_INVALID;
    }
    const md_attribute_t* existing = md_attributes_find(&sys->attributes, path);
    if (existing) {
        md_attributes_remove(&sys->attributes, existing->id);
    }
    return md_attributes_alias(&sys->attributes, target, path, (str_t){0}, (str_t){0});
}

md_attribute_id_t md_qm_publish_or_alias(md_system_t* sys, md_attribute_id_t alpha_id, str_t path, str_t label,
                                         md_unit_t unit, const double* alpha_values, const double* beta_values, size_t count) {
    if (beta_values && beta_values == alpha_values) {
        return md_qm_alias(sys, alpha_id, path);
    }
    return md_qm_publish_series(sys, path, label, unit, beta_values, count);
}

str_t md_qm_attribute_path(char* buf, size_t cap, str_t group, str_t name) {
    int len = snprintf(buf, cap, STR_FMT "/" STR_FMT, STR_ARG(group), STR_ARG(name));
    if (len <= 0 || (size_t)len >= cap) {
        MD_LOG_ERROR("Attribute path '" STR_FMT "/" STR_FMT "' does not fit in %zu characters", STR_ARG(group), STR_ARG(name), cap - 1);
        return (str_t){0};
    }
    for (int c = (int)group.len + 1; c < len; ++c) {
        if (buf[c] == '/') buf[c] = '_';
    }
    return str_from_cstrn(buf, (size_t)len);
}

size_t md_qm_sph_to_cart_coefficients(double* dst, const double* src, size_t num_mo, const md_gto_basis_t* basis) {
    if (!dst || !src || !basis || num_mo == 0) {
        return 0;
    }
    const size_t n_sph  = md_gto_basis_num_sph_ao(basis);
    const size_t n_cart = md_gto_basis_num_ao(basis);
    if (n_sph == 0 || n_cart == 0) {
        return 0;
    }
    for (size_t mo = 0; mo < num_mo; ++mo) {
        if (md_gto_sph_to_cart_vector(dst + mo * n_cart, src + mo * n_sph, basis) != n_cart) {
            MD_LOG_ERROR("Spherical to Cartesian conversion failed for orbital %zu", mo);
            return 0;
        }
    }
    return num_mo;
}

// ---------------------------------------------------------------------------
// Normalisation
// ---------------------------------------------------------------------------

// (2n-1)!! with (-1)!! == 1. Same recurrence md_gto_cart_norm_factor is built on; kept local
// because it is three lines and importing it would mean exporting it.
static double qm_double_factorial_odd(int n) {
    double r = 1.0;
    for (int k = 2 * n - 1; k > 1; k -= 2) {
        r *= (double)k;
    }
    return r;
}

// <mono_a R | mono_b R> for two monomials on the SAME centre with a combined exponent p, divided by
// the (pi/p)^1.5 every term shares. Odd powers integrate to zero, which is what makes a Cartesian
// shell's cross terms vanish and lets the norm below be a single product.
static double qm_same_centre_monomial_overlap(const int lmn_a[3], const int lmn_b[3], double p) {
    double v = 1.0;
    for (int k = 0; k < 3; ++k) {
        const int n = lmn_a[k] + lmn_b[k];
        if (n & 1) {
            return 0.0;
        }
        v *= qm_double_factorial_odd(n / 2) / pow(2.0 * p, (double)(n / 2));
    }
    return v;
}

// <cart AO | cart AO> for one primitive pair, in md_gto's convention: the per monomial factor
// f(i,j,k) exactly cancels the (2i-1)!!(2j-1)!!(2k-1)!! the integral produces, so every Cartesian
// AO of a shell has the same norm and this closed form serves all of them.
static double qm_cart_self_overlap(uint32_t l, double p) {
    return pow(PI / p, 1.5) / pow(2.0 * p, (double)l);
}

// ||md_gto's spherical AO||^2 / ||its Cartesian AO||^2 for angular momentum l.
//
// Derived from md_gto's OWN expansion tables rather than tabulated here, by pushing a unit
// spherical coefficient through md_gto_sph_to_cart_vector and measuring what comes back, so the two
// cannot drift apart. It is a pure function of l - the radial part cancels - so it is computed once
// per l and cached.
static double qm_spherical_norm_ratio(uint32_t l) {
    static double cache[MD_GTO_MAX_ANGULAR_MOMENTUM + 1];
    static bool   cached[MD_GTO_MAX_ANGULAR_MOMENTUM + 1];

    if (l > MD_GTO_MAX_ANGULAR_MOMENTUM) {
        return 0.0;
    }
    if (cached[l]) {
        return cache[l];
    }

    const uint32_t num_cart = md_gto_num_cart_ao(l);
    const uint32_t num_sph  = md_gto_num_sph_ao(l);

    // One shell, one primitive of unit exponent: the ratio does not depend on either.
    md_gto_shell_t shell = { .atom_idx = 0, .primitive_offset = 0, .num_primitives = 1, .l = l };
    float alpha = 1.0f;
    float coeff = 1.0f;
    md_gto_basis_t basis = { .num_shells = 1, .num_primitives = 1, .shells = &shell, .alpha = &alpha, .coeff = &coeff };

    double sph[2 * MD_GTO_MAX_ANGULAR_MOMENTUM + 1] = {0};
    double cart[((MD_GTO_MAX_ANGULAR_MOMENTUM + 1) * (MD_GTO_MAX_ANGULAR_MOMENTUM + 2)) / 2] = {0};

    // The first spherical function will do: md_gto's expansion tables give all 2l+1 of them the same
    // norm, and the ratio is what is wanted, not any one of them.
    sph[0] = 1.0;
    double ratio = 0.0;
    if (md_gto_sph_to_cart_vector(cart, sph, &basis) == num_cart) {
        const double p = 2.0;  // alpha + alpha
        double num = 0.0;
        for (uint32_t a = 0; a < num_cart; ++a) {
            if (cart[a] == 0.0) continue;
            int ia, ja, ka;
            md_gto_cart_ijk(&ia, &ja, &ka, l, a);
            const int lmn_a[3] = { ia, ja, ka };
            const double fa = md_gto_cart_norm_factor(ia, ja, ka);
            for (uint32_t b = 0; b < num_cart; ++b) {
                if (cart[b] == 0.0) continue;
                int ib, jb, kb;
                md_gto_cart_ijk(&ib, &jb, &kb, l, b);
                const int lmn_b[3] = { ib, jb, kb };
                const double fb = md_gto_cart_norm_factor(ib, jb, kb);
                num += cart[a] * cart[b] * fa * fb * qm_same_centre_monomial_overlap(lmn_a, lmn_b, p) * pow(PI / p, 1.5);
            }
        }
        ratio = num / qm_cart_self_overlap(l, p);
    } else {
        MD_LOG_ERROR("Failed to expand a unit spherical coefficient for l=%u", l);
    }
    (void)num_sph;

    cache[l]  = ratio;
    cached[l] = true;
    return ratio;
}

double md_qm_primitive_norm_factor(uint32_t l, double alpha) {
    if (alpha <= 0.0) {
        return 0.0;
    }
    return 1.0 / sqrt(qm_cart_self_overlap(l, 2.0 * alpha));
}

double md_qm_shell_norm_factor(uint32_t l, const double alpha[], const double coeff[], size_t num_primitives, bool spherical) {
    if (!alpha || !coeff || num_primitives == 0) {
        return 0.0;
    }

    double sum = 0.0;
    for (size_t i = 0; i < num_primitives; ++i) {
        for (size_t j = 0; j < num_primitives; ++j) {
            sum += coeff[i] * coeff[j] * qm_cart_self_overlap(l, alpha[i] + alpha[j]);
        }
    }
    if (spherical) {
        sum *= qm_spherical_norm_ratio(l);
    }
    if (!(sum > 0.0)) {
        return 0.0;
    }
    return sqrt(sum);
}

double md_qm_cart_coeff_factor(uint32_t l, uint32_t cart_idx) {
    int i, j, k;
    if (!md_gto_cart_ijk(&i, &j, &k, l, cart_idx)) {
        return 0.0;
    }
    // f(l,0,0) / f(i,j,k), with f the per monomial factor md_gto documents.
    return sqrt(qm_double_factorial_odd(i) * qm_double_factorial_odd(j) * qm_double_factorial_odd(k) / qm_double_factorial_odd((int)l));
}

// ---------------------------------------------------------------------------
// Publishing
// ---------------------------------------------------------------------------

bool md_qm_publish_basis(md_system_t* sys, const md_gto_basis_t* basis) {
    ASSERT(sys);

    if (!basis || basis->num_shells == 0 || basis->num_primitives == 0) {
        return false;
    }
    if (!sys->attributes.alloc) {
        MD_LOG_ERROR("Attribute table allocator not set; the system has not been initialised");
        return false;
    }

    const size_t num_shells     = basis->num_shells;
    const size_t num_primitives = basis->num_primitives;
    const size_t shell_stride   = sizeof(md_gto_shell_t);

    bool ok = true;
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, STR_LIT("basis/shell/atom_index"),       STR_LIT("Atom Index"),       md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis->shells->atom_idx,         shell_stride, num_shells) && ok;
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, STR_LIT("basis/shell/primitive_offset"), STR_LIT("Primitive Offset"), md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis->shells->primitive_offset, shell_stride, num_shells) && ok;
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, STR_LIT("basis/shell/primitive_count"),  STR_LIT("Primitive Count"),  md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis->shells->num_primitives,   shell_stride, num_shells) && ok;
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, STR_LIT("basis/shell/angular_momentum"), STR_LIT("Angular Momentum"), md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis->shells->l,                shell_stride, num_shells) && ok;

    // Exponents are bohr^-2 and the contraction coefficients carry the shell's radial normalisation;
    // the per monomial factor is applied at evaluation. See the AO CONVENTION block in md_gto.h -
    // these values only mean anything against it.
    const md_unit_t inv_bohr_sq = md_unit_pow(md_unit_bohr_radius(), -2);
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, STR_LIT("basis/primitive/exponent"),    STR_LIT("Exponent"),    inv_bohr_sq,    MD_ATTRIBUTE_TYPE_F32, basis->alpha, sizeof(float), num_primitives) && ok;
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, STR_LIT("basis/primitive/coefficient"), STR_LIT("Coefficient"), md_unit_none(), MD_ATTRIBUTE_TYPE_F32, basis->coeff, sizeof(float), num_primitives) && ok;

    return ok;
}

bool md_qm_publish_atoms(md_system_t* sys, const uint8_t atomic_number[], const dvec3_t coord_angstrom[], size_t count) {
    ASSERT(sys);

    if (count == 0 || !atomic_number || !coord_angstrom) {
        return false;
    }
    if (!sys->attributes.alloc) {
        MD_LOG_ERROR("Attribute table allocator not set; the system has not been initialised");
        return false;
    }

    md_attribute_format_t z_format = {
        .type = MD_ATTRIBUTE_TYPE_U8, .components = 1, .rank = 1, .shape = { (uint32_t)count },
    };
    md_attribute_id_t z_id = md_qm_publish(sys, STR_LIT("qm/atom/atomic_number"), STR_LIT("Atomic Number"),
                                           md_unit_none(), z_format, atomic_number, count * sizeof(uint8_t));

    // Angstrom, matching the system's own coordinates rather than the bohr the evaluator works in.
    // This is the geometry the CALCULATION was run at, which need not be where the system's atoms
    // are now - a trajectory frame or an optimisation step moves them.
    md_attribute_id_t xyz_id = md_qm_publish_vec3_series(sys, STR_LIT("qm/atom/coordinate"), STR_LIT("Coordinate"),
                                                         md_unit_angstrom(), coord_angstrom, count);

    return z_id != MD_ATTRIBUTE_INVALID && xyz_id != MD_ATTRIBUTE_INVALID;
}

// ---------------------------------------------------------------------------
// The AO overlap
// ---------------------------------------------------------------------------

// The one dimensional Gaussian overlap of two monomials, integral of
// (x-Ax)^la (x-Bx)^lb exp(-p (x-Px)^2) dx, divided by the sqrt(pi/p) all three axes share.
// The expansion is written out rather than recursed: l is at most 4 here, so the double sum is 25
// terms in the worst case and needs no table.
static double qm_overlap_1d(int la, int lb, double pa, double pb, double p) {
    static const int binom[5][5] = {
        {1,0,0,0,0}, {1,1,0,0,0}, {1,2,1,0,0}, {1,3,3,1,0}, {1,4,6,4,1},
    };
    double sum = 0.0;
    for (int i = 0; i <= la; ++i) {
        for (int j = 0; j <= lb; ++j) {
            const int n = i + j;
            if (n & 1) {
                continue;
            }
            sum += binom[la][la - i] * binom[lb][lb - j]
                 * pow(pa, (double)(la - i)) * pow(pb, (double)(lb - j))
                 * qm_double_factorial_odd(n / 2) / pow(2.0 * p, (double)(n / 2));
        }
    }
    return sum;
}

// <mono_a exp(-a|r-A|^2) | mono_b exp(-b|r-B|^2)>, the whole three dimensional primitive overlap.
static double qm_primitive_overlap(double a, const double A[3], const int lmn_a[3],
                                   double b, const double B[3], const int lmn_b[3]) {
    const double p  = a + b;
    const double mu = a * b / p;
    double r2 = 0.0;
    for (int k = 0; k < 3; ++k) {
        const double d = A[k] - B[k];
        r2 += d * d;
    }
    double v = exp(-mu * r2) * pow(PI / p, 1.5);
    for (int k = 0; k < 3; ++k) {
        const double P = (a * A[k] + b * B[k]) / p;
        v *= qm_overlap_1d(lmn_a[k], lmn_b[k], P - A[k], P - B[k], p);
    }
    return v;
}

size_t md_qm_compute_overlap(double* out, const md_gto_basis_t* basis, const dvec3_t atom_coord_bohr[], size_t num_atoms) {
    if (!out || !basis || !atom_coord_bohr || basis->num_shells == 0) {
        return 0;
    }

    const size_t num_ao = md_gto_basis_num_ao(basis);
    MEMSET(out, 0, sizeof(double) * num_ao * num_ao);

    // Shell offsets into the AO axis, so the pair loop below can start anywhere.
    md_temp_scope_t temp = md_temp_begin();
    uint32_t* offset = md_temp_alloc_array(temp, uint32_t, basis->num_shells);
    if (!offset) {
        md_temp_end(temp);
        return 0;
    }
    uint32_t running = 0;
    for (uint32_t s = 0; s < basis->num_shells; ++s) {
        offset[s] = running;
        running += md_gto_num_cart_ao(basis->shells[s].l);
        if (basis->shells[s].atom_idx >= num_atoms) {
            MD_LOG_ERROR("Shell %u is centred on atom %u, and only %zu were supplied", s, basis->shells[s].atom_idx, num_atoms);
            md_temp_end(temp);
            return 0;
        }
    }

    for (uint32_t si = 0; si < basis->num_shells; ++si) {
        const md_gto_shell_t* a = &basis->shells[si];
        const dvec3_t         A = atom_coord_bohr[a->atom_idx];
        const double          Av[3] = { A.x, A.y, A.z };
        const uint32_t        na = md_gto_num_cart_ao(a->l);

        for (uint32_t sj = si; sj < basis->num_shells; ++sj) {
            const md_gto_shell_t* b = &basis->shells[sj];
            const dvec3_t         B = atom_coord_bohr[b->atom_idx];
            const double          Bv[3] = { B.x, B.y, B.z };
            const uint32_t        nb = md_gto_num_cart_ao(b->l);

            for (uint32_t ca = 0; ca < na; ++ca) {
                int ia, ja, ka;
                md_gto_cart_ijk(&ia, &ja, &ka, a->l, ca);
                const int    lmn_a[3] = { ia, ja, ka };
                const double fa = md_gto_cart_norm_factor(ia, ja, ka);

                for (uint32_t cb = 0; cb < nb; ++cb) {
                    int ib, jb, kb;
                    md_gto_cart_ijk(&ib, &jb, &kb, b->l, cb);
                    const int    lmn_b[3] = { ib, jb, kb };
                    const double fb = md_gto_cart_norm_factor(ib, jb, kb);

                    double v = 0.0;
                    for (uint32_t p = 0; p < a->num_primitives; ++p) {
                        const double alpha_a = basis->alpha[a->primitive_offset + p];
                        const double coeff_a = basis->coeff[a->primitive_offset + p];
                        for (uint32_t q = 0; q < b->num_primitives; ++q) {
                            const double alpha_b = basis->alpha[b->primitive_offset + q];
                            const double coeff_b = basis->coeff[b->primitive_offset + q];
                            v += coeff_a * coeff_b * qm_primitive_overlap(alpha_a, Av, lmn_a, alpha_b, Bv, lmn_b);
                        }
                    }
                    v *= fa * fb;

                    const size_t i = offset[si] + ca;
                    const size_t j = offset[sj] + cb;
                    out[i * num_ao + j] = v;
                    out[j * num_ao + i] = v;
                }
            }
        }
    }

    md_temp_end(temp);
    return num_ao;
}

// Reads the basis and the QM geometry back out of the table and integrates. The geometry is
// published in Angstrom, matching the system's own; the exponents are in bohr^-2, so it converts.
static size_t qm_overlap_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
    (void)slice;  // one whole {A,A} matrix, not indexed by anything a slice could fix
    md_system_t* sys = (md_system_t*)user_data;

    const md_attribute_t* coord = md_attributes_find(&sys->attributes, STR_LIT("qm/atom/coordinate"));
    if (!coord) {
        MD_LOG_ERROR("'" STR_FMT "' has no qm/atom/coordinate to integrate over", STR_ARG(attr->path));
        return 0;
    }

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    size_t written = 0;
    md_gto_basis_t basis = {0};
    const size_t num_atoms = md_attribute_value_count(&coord->format);
    double*  xyz  = md_temp_alloc_array(temp, double,  num_atoms * 3);
    dvec3_t* bohr = md_temp_alloc_array(temp, dvec3_t, num_atoms);

    if (xyz && bohr && md_gto_basis_extract_attributes(&basis, &sys->attributes, temp_alloc)
        && md_attribute_extract_f64(xyz, num_atoms * 3, coord, md_unit_none()) == num_atoms * 3) {

        const double angstrom_to_bohr = 1.0 / 0.529177210903;
        for (size_t i = 0; i < num_atoms; ++i) {
            bohr[i] = (dvec3_t){ xyz[i * 3 + 0] * angstrom_to_bohr, xyz[i * 3 + 1] * angstrom_to_bohr, xyz[i * 3 + 2] * angstrom_to_bohr };
        }
        const size_t num_ao = md_gto_basis_num_ao(&basis);
        if (num_ao * num_ao != cap) {
            MD_LOG_ERROR("'" STR_FMT "' was asked for %zu values and its basis spans %zu atomic orbitals", STR_ARG(attr->path), cap, num_ao);
        } else if (md_qm_compute_overlap((double*)dst, &basis, bohr, num_atoms) == num_ao) {
            written = cap;
        }
    }

    md_temp_end(temp);
    return written;
}

bool md_qm_publish_overlap(md_system_t* sys) {
    ASSERT(sys);

    if (!sys->attributes.alloc) {
        MD_LOG_ERROR("Attribute table allocator not set; the system has not been initialised");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    md_gto_basis_t basis = {0};
    size_t num_ao = 0;
    if (md_gto_basis_extract_attributes(&basis, &sys->attributes, md_temp_allocator(temp))) {
        num_ao = md_gto_basis_num_ao(&basis);
    }
    md_temp_end(temp);

    if (num_ao == 0 || !md_attributes_find(&sys->attributes, STR_LIT("qm/atom/coordinate"))) {
        return false;
    }

    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2, .shape = { (uint32_t)num_ao, (uint32_t)num_ao },
    };
    md_attribute_virtual_t virt = { .provider = qm_overlap_provider, .user_data = sys };
    return md_qm_publish_virtual(sys, STR_LIT("basis/overlap"), STR_LIT("AO Overlap"), md_unit_none(), format, &virt) != MD_ATTRIBUTE_INVALID;
}

// D_ij = sum_mo occ[mo] * C[mo][i] * C[mo][j]. Exact, not an approximation: an SCF density IS that
// sum, fractional occupations included, so nothing is gained by storing it separately.
static void qm_build_density(double* out, const double* coeff, const double* occ, size_t num_mo, size_t num_ao) {
    MEMSET(out, 0, sizeof(double) * num_ao * num_ao);
    for (size_t mo = 0; mo < num_mo; ++mo) {
        const double w = occ[mo];
        if (w == 0.0) continue;
        const double* c = coeff + mo * num_ao;
        for (size_t i = 0; i < num_ao; ++i) {
            const double wci = w * c[i];
            if (wci == 0.0) continue;
            for (size_t j = 0; j < num_ao; ++j) {
                out[i * num_ao + j] += wci * c[j];
            }
        }
    }
}

static size_t qm_density_provide(void* dst, size_t cap, const md_attribute_t* attr, void* user_data, bool beta) {
    md_system_t* sys = (md_system_t*)user_data;

    const str_t coeff_path = beta ? STR_LIT("orbital/beta/coefficient") : STR_LIT("orbital/alpha/coefficient");
    const str_t occ_path   = beta ? STR_LIT("orbital/beta/occupation")  : STR_LIT("orbital/alpha/occupation");

    const md_attribute_t* coeff = md_attributes_find(&sys->attributes, coeff_path);
    const md_attribute_t* occ   = md_attributes_find(&sys->attributes, occ_path);
    if (!coeff || !occ) {
        MD_LOG_ERROR("'" STR_FMT "' is missing the orbital data it reconstructs from", STR_ARG(attr->path));
        return 0;
    }

    const size_t num_mo = coeff->format.shape[0];
    const size_t num_ao = coeff->format.shape[1];
    if (md_attribute_element_count(&occ->format) != num_mo) {
        MD_LOG_ERROR("'" STR_FMT "': occupation holds %zu values, coefficients hold %zu orbitals",
                     STR_ARG(attr->path), md_attribute_element_count(&occ->format), num_mo);
        return 0;
    }
    if (num_ao * num_ao != cap) {
        MD_LOG_ERROR("'" STR_FMT "' was asked for %zu values and its coefficients span %zu atomic orbitals",
                     STR_ARG(attr->path), cap, num_ao);
        return 0;
    }

    md_temp_scope_t temp = md_temp_begin();
    double* coeff_data = md_temp_alloc_array(temp, double, num_mo * num_ao);
    double* occ_data   = md_temp_alloc_array(temp, double, num_mo);

    bool ok = coeff_data && occ_data;
    if (!ok) {
        MD_LOG_ERROR("'" STR_FMT "': failed to allocate %zu doubles of scratch", STR_ARG(attr->path), num_mo * num_ao + num_mo);
    }
    if (ok) {
        ok = md_attribute_extract_f64(coeff_data, num_mo * num_ao, coeff, md_unit_none()) == num_mo * num_ao
          && md_attribute_extract_f64(occ_data, num_mo, occ, md_unit_none()) == num_mo;
        if (ok) {
            qm_build_density((double*)dst, coeff_data, occ_data, num_mo, num_ao);
        }
    }

    md_temp_end(temp);
    return ok ? cap : 0;
}

static size_t qm_alpha_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
    (void)slice;  // one whole {A,A} matrix, not indexed by anything a slice could fix
    return qm_density_provide(dst, cap, attr, user_data, false);
}

static size_t qm_beta_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
    (void)slice;
    return qm_density_provide(dst, cap, attr, user_data, true);
}

// alpha +/- beta. Both halves are themselves virtual, so this is a derivation over derivations,
// legal because the dependency graph stays acyclic. A restricted calculation gets it for free:
// beta is an ALIAS of alpha there, so the total is twice alpha and the difference is zero with no
// special case anywhere.
static size_t qm_density_combine(void* dst, size_t cap, const md_attribute_t* attr, void* user_data, double beta_scale) {
    md_system_t* sys = (md_system_t*)user_data;

    const md_attribute_t* alpha = md_attributes_find(&sys->attributes, STR_LIT("orbital/alpha/density"));
    const md_attribute_t* beta  = md_attributes_find(&sys->attributes, STR_LIT("orbital/beta/density"));
    if (!alpha || !beta) {
        MD_LOG_ERROR("'" STR_FMT "' is missing a spin density it combines", STR_ARG(attr->path));
        return 0;
    }

    md_temp_scope_t temp = md_temp_begin();
    double* beta_data = md_temp_alloc_array(temp, double, cap);

    bool ok = beta_data != NULL;
    if (!ok) {
        MD_LOG_ERROR("'" STR_FMT "': failed to allocate %zu doubles of scratch", STR_ARG(attr->path), cap);
    }
    if (ok && md_attribute_extract_f64((double*)dst, cap, alpha, md_unit_none()) != cap) {
        MD_LOG_ERROR("'" STR_FMT "': could not read '" STR_FMT "'", STR_ARG(attr->path), STR_ARG(alpha->path));
        ok = false;
    }
    if (ok && md_attribute_extract_f64(beta_data, cap, beta, md_unit_none()) != cap) {
        MD_LOG_ERROR("'" STR_FMT "': could not read '" STR_FMT "'", STR_ARG(attr->path), STR_ARG(beta->path));
        ok = false;
    }
    if (ok) {
        double* out = (double*)dst;
        for (size_t i = 0; i < cap; ++i) {
            out[i] += beta_scale * beta_data[i];
        }
    }

    md_temp_end(temp);
    return ok ? cap : 0;
}

static size_t qm_total_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
    (void)slice;
    return qm_density_combine(dst, cap, attr, user_data, 1.0);
}

static size_t qm_difference_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
    (void)slice;
    return qm_density_combine(dst, cap, attr, user_data, -1.0);
}

bool md_qm_publish_orbital_densities(md_system_t* sys) {
    ASSERT(sys);

    if (!sys->attributes.alloc) {
        MD_LOG_ERROR("Attribute table allocator not set; the system has not been initialised");
        return false;
    }

    // EVERYTHING THIS NEEDS TO KNOW IS READ FIRST, as values. md_attributes_find returns a pointer
    // INTO the table's array, and publishing anything can move it - so a decision taken from a
    // pointer fetched before the first publish is taken from freed memory. It reads as a plausible
    // answer, which is how this came to be written that way in the first place: the restricted
    // case quietly stopped aliasing beta and reconstructed it a second time instead.
    const md_attribute_t* alpha_coeff = md_attributes_find(&sys->attributes, STR_LIT("orbital/alpha/coefficient"));
    if (!alpha_coeff || alpha_coeff->format.rank != 2) {
        return false;
    }
    const md_attribute_t* beta_coeff = md_attributes_find(&sys->attributes, STR_LIT("orbital/beta/coefficient"));
    const md_attribute_t* alpha_occ  = md_attributes_find(&sys->attributes, STR_LIT("orbital/alpha/occupation"));
    const md_attribute_t* beta_occ   = md_attributes_find(&sys->attributes, STR_LIT("orbital/beta/occupation"));

    const uint32_t num_ao   = alpha_coeff->format.shape[1];
    const bool     has_beta = beta_coeff != NULL;

    // Same coefficients AND same occupations means one density, reachable under both names. That is
    // the restricted case, and comparing the DATA rather than a spin flag is what also gets
    // restricted open shell right, where the orbitals are shared and the occupations are not.
    const bool same_density = has_beta
                           && md_attribute_same_data(alpha_coeff, beta_coeff)
                           && alpha_occ && beta_occ && md_attribute_same_data(alpha_occ, beta_occ);

    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2, .shape = { num_ao, num_ao },
    };

    md_attribute_virtual_t alpha_virt = { .provider = qm_alpha_density_provider, .user_data = sys };
    md_attribute_id_t alpha_id = md_qm_publish_virtual(sys, STR_LIT("orbital/alpha/density"), STR_LIT("Alpha Density"),
                                                       md_unit_none(), format, &alpha_virt);
    if (alpha_id == MD_ATTRIBUTE_INVALID) {
        return false;
    }

    if (!has_beta) {
        return true;
    }

    if (same_density) {
        md_qm_alias(sys, alpha_id, STR_LIT("orbital/beta/density"));
    } else {
        md_attribute_virtual_t beta_virt = { .provider = qm_beta_density_provider, .user_data = sys };
        md_qm_publish_virtual(sys, STR_LIT("orbital/beta/density"), STR_LIT("Beta Density"), md_unit_none(), format, &beta_virt);
    }

    md_attribute_virtual_t total_virt = { .provider = qm_total_density_provider,      .user_data = sys };
    md_attribute_virtual_t diff_virt  = { .provider = qm_difference_density_provider, .user_data = sys };
    md_qm_publish_virtual(sys, STR_LIT("orbital/total/density"),      STR_LIT("Total Density"),           md_unit_none(), format, &total_virt);
    md_qm_publish_virtual(sys, STR_LIT("orbital/difference/density"), STR_LIT("Spin Difference Density"), md_unit_none(), format, &diff_virt);

    return true;
}
