#include <md_qm.h>

#include <md_system.h>
#include <md_gto.h>
#include <md_gto_iso.h>
#include <md_gto_dmat.h>

#include <core/md_allocator.h>
#include <core/md_format.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_unit.h>
#include <core/md_parse.h>

#include <math.h>
#include <float.h>
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

md_attribute_id_t md_qm_publish(md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit,
                                md_attribute_format_t format, const void* data, size_t byte_size) {
    MD_ASSERT(sys);
    return md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
        .path      = path,
        .format    = format,
        .unit      = unit,
        .label     = label,
        .data      = data,
        .byte_size = byte_size,
    });
}

md_attribute_id_t md_qm_publish_virtual(md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit,
                                        md_attribute_format_t format, const md_attribute_virtual_t* virt) {
    MD_ASSERT(sys);
    return md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
        .path   = path,
        .format = format,
        .unit   = unit,
        .label  = label,
        .virt   = virt,
    });
}

md_attribute_id_t md_qm_publish_scalar(md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit, double value) {
    // The value is copied, so a local is fine.
    md_attribute_format_t format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 0 };
    return md_qm_publish(sys, path, label, unit, format, &value, sizeof(double));
}

md_attribute_id_t md_qm_publish_series(md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit, const double* values, size_t count) {
    if (!values || count == 0) {
        return MD_ATTRIBUTE_INVALID;
    }
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)count },
    };
    return md_qm_publish(sys, path, label, unit, format, values, count * sizeof(double));
}

md_attribute_id_t md_qm_publish_principal_series(md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit, const double* values, size_t count) {
    MD_ASSERT(sys);
    if (!values || count == 0) {
        return MD_ATTRIBUTE_INVALID;
    }
    return md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
        .path      = path,
        .format    = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)count }},
        .flags     = MD_ATTRIBUTE_FLAG_PRINCIPAL,
        .unit      = unit,
        .label     = label,
        .data      = values,
        .byte_size = count * sizeof(double),
    });
}

md_attribute_id_t md_qm_publish_vec3_series(md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit, const md_double3_t* values, size_t count) {
    if (!values || count == 0) {
        return MD_ATTRIBUTE_INVALID;
    }
    // md_double3_t is three contiguous doubles, so the source array is already the interleaved layout an
    // attribute stores and this is a straight copy.
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 1, .shape = { (uint32_t)count },
    };
    return md_qm_publish(sys, path, label, unit, format, values, count * 3 * sizeof(double));
}

md_attribute_id_t md_qm_publish_matrix(md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit, const double* values, size_t rows, size_t cols) {
    if (!values || rows == 0 || cols == 0) {
        return MD_ATTRIBUTE_INVALID;
    }
    // Row major with the last index fastest, which is what every 2D quantity in the tree is.
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2, .shape = { (uint32_t)rows, (uint32_t)cols },
    };
    return md_qm_publish(sys, path, label, unit, format, values, rows * cols * sizeof(double));
}

md_attribute_id_t md_qm_publish_str(md_system_t* sys, md_str_t path, md_str_t label, md_str_t value) {
    if (md_str_empty(value)) {
        return MD_ATTRIBUTE_INVALID;
    }
    // A single string is rank 1 {1}, by the same rule that makes a single 3-vector rank 1 {1} of 3
    // components. The descriptor carries the TEXT and the table stores a handle - see STRINGS in
    // md_system.h.
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 1, .shape = { 1 },
    };
    return md_qm_publish(sys, path, label, md_unit_none(), format, &value, sizeof(md_str_t));
}

md_attribute_id_t md_qm_publish_strings(md_system_t* sys, md_str_t path, md_str_t label, const md_str_t* values, size_t count) {
    if (!values || count == 0) {
        return MD_ATTRIBUTE_INVALID;
    }
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 1, .shape = { (uint32_t)count },
    };
    return md_qm_publish(sys, path, label, md_unit_none(), format, values, count * sizeof(md_str_t));
}

md_attribute_id_t md_qm_publish_column(md_system_t* sys, md_str_t path, md_str_t label, md_unit_t unit, md_attribute_type_t type,
                                       const void* base, size_t stride, size_t count) {
    MD_ASSERT(sys);
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
        md_memcpy(dst + i * elem_size, src + i * stride, elem_size);
    }
    return id;
}

md_attribute_id_t md_qm_publish_origin(md_system_t* sys, md_str_t path, md_double3_t origin) {
    md_attribute_format_t format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 0 };
    return md_qm_publish(sys, path, (md_str_t){0}, md_unit_angstrom(), format, &origin, 3 * sizeof(double));
}

md_attribute_id_t md_qm_alias(md_system_t* sys, md_attribute_id_t target, md_str_t path) {
    MD_ASSERT(sys);
    if (target == MD_ATTRIBUTE_INVALID) {
        return MD_ATTRIBUTE_INVALID;
    }
    const md_attribute_t* existing = md_attributes_find(&sys->attributes, path);
    if (existing) {
        md_attributes_remove(&sys->attributes, existing->id);
    }
    return md_attributes_alias(&sys->attributes, target, path, (md_str_t){0}, (md_str_t){0});
}

md_attribute_id_t md_qm_publish_or_alias(md_system_t* sys, md_attribute_id_t alpha_id, md_str_t path, md_str_t label,
                                         md_unit_t unit, const double* alpha_values, const double* beta_values, size_t count) {
    if (beta_values && beta_values == alpha_values) {
        return md_qm_alias(sys, alpha_id, path);
    }
    return md_qm_publish_series(sys, path, label, unit, beta_values, count);
}

md_str_t md_qm_attribute_path(char* buf, size_t cap, md_str_t group, md_str_t name) {
    int len = md_snprintf(buf, cap, MD_STR_FMT "/" MD_STR_FMT, MD_STR_ARG(group), MD_STR_ARG(name));
    if (len <= 0 || (size_t)len >= cap) {
        MD_LOG_ERROR("Attribute path '" MD_STR_FMT "/" MD_STR_FMT "' does not fit in %zu characters", MD_STR_ARG(group), MD_STR_ARG(name), cap - 1);
        return (md_str_t){0};
    }
    for (int c = (int)group.len + 1; c < len; ++c) {
        if (buf[c] == '/') buf[c] = '_';
    }
    return md_str_from_cstrn(buf, (size_t)len);
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
    if (md_gto_sph_to_cart_orbitals(dst, src, n_sph, 1, NULL, num_mo, basis, NULL) != num_mo) {
        MD_LOG_ERROR("Spherical to Cartesian conversion of %zu orbitals failed", num_mo);
        return 0;
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
    return pow(MD_PI / p, 1.5) / pow(2.0 * p, (double)l);
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
                num += cart[a] * cart[b] * fa * fb * qm_same_centre_monomial_overlap(lmn_a, lmn_b, p) * pow(MD_PI / p, 1.5);
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
    MD_ASSERT(sys);

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
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, MD_STR_LIT("basis/shell/atom_index"),       MD_STR_LIT("Atom Index"),       md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis->shells->atom_idx,         shell_stride, num_shells) && ok;
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, MD_STR_LIT("basis/shell/primitive_offset"), MD_STR_LIT("Primitive Offset"), md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis->shells->primitive_offset, shell_stride, num_shells) && ok;
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, MD_STR_LIT("basis/shell/primitive_count"),  MD_STR_LIT("Primitive Count"),  md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis->shells->num_primitives,   shell_stride, num_shells) && ok;
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, MD_STR_LIT("basis/shell/angular_momentum"), MD_STR_LIT("Angular Momentum"), md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis->shells->l,                shell_stride, num_shells) && ok;

    // Exponents are bohr^-2 and the contraction coefficients carry the shell's radial normalisation;
    // the per monomial factor is applied at evaluation. See the AO CONVENTION block in md_gto.h -
    // these values only mean anything against it.
    const md_unit_t inv_bohr_sq = md_unit_pow(md_unit_bohr_radius(), -2);
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, MD_STR_LIT("basis/primitive/exponent"),    MD_STR_LIT("Exponent"),    inv_bohr_sq,    MD_ATTRIBUTE_TYPE_F32, basis->alpha, sizeof(float), num_primitives) && ok;
    ok = MD_ATTRIBUTE_INVALID != md_qm_publish_column(sys, MD_STR_LIT("basis/primitive/coefficient"), MD_STR_LIT("Coefficient"), md_unit_none(), MD_ATTRIBUTE_TYPE_F32, basis->coeff, sizeof(float), num_primitives) && ok;

    return ok;
}

bool md_qm_publish_atoms(md_system_t* sys, const uint8_t atomic_number[], const double nuclear_charge[], const md_double3_t coord_angstrom[], size_t count) {
    MD_ASSERT(sys);

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
    md_attribute_id_t z_id = md_qm_publish(sys, MD_STR_LIT("qm/atom/atomic_number"), MD_STR_LIT("Atomic Number"),
                                           md_unit_none(), z_format, atomic_number, count * sizeof(uint8_t));

    // The charge the electrons see: the atomic number unless the reader knows better (an effective
    // core potential). Published in every case, see md_qm.h.
    md_attribute_id_t q_id = MD_ATTRIBUTE_INVALID;
    {
        md_temp_scope_t temp = md_temp_begin();
        double* charge = md_temp_alloc_array(temp, double, count);
        if (charge) {
            for (size_t i = 0; i < count; ++i) {
                charge[i] = nuclear_charge ? nuclear_charge[i] : (double)atomic_number[i];
            }
            q_id = md_qm_publish_series(sys, MD_STR_LIT("qm/atom/nuclear_charge"), MD_STR_LIT("Nuclear Charge"),
                                        md_unit_elementary_charge(), charge, count);
        }
        md_temp_end(temp);
    }

    // Angstrom, matching the system's own coordinates rather than the bohr the evaluator works in.
    // This is the geometry the CALCULATION was run at, which need not be where the system's atoms
    // are now - a trajectory frame or an optimisation step moves them.
    md_attribute_id_t xyz_id = md_qm_publish_vec3_series(sys, MD_STR_LIT("qm/atom/coordinate"), MD_STR_LIT("Coordinate"),
                                                         md_unit_angstrom(), coord_angstrom, count);

    return z_id != MD_ATTRIBUTE_INVALID && q_id != MD_ATTRIBUTE_INVALID && xyz_id != MD_ATTRIBUTE_INVALID;
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

size_t md_qm_compute_overlap(double* out, const md_gto_basis_t* basis, const md_double3_t atom_coord_bohr[], size_t num_atoms) {
    if (!out || !basis || !atom_coord_bohr || basis->num_shells == 0) {
        return 0;
    }

    const size_t num_ao = md_gto_basis_num_ao(basis);
    md_memset(out, 0, sizeof(double) * num_ao * num_ao);

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

    // Per shell pair and primitive pair, the Gaussian prefactor and the three axes' one dimensional overlaps of
    // every monomial power once (at most 5 x 5 each), and every Cartesian pair of the two shells as a product of
    // three of those: exp(-mu r^2) (pi/p)^(3/2) times the product of the three axes' qm_overlap_1d.
    // A primitive pair whose prefactor is below e^-46 (1e-20) adds nothing a double holds beside the rest.
    enum { L1 = MD_GTO_MAX_ANGULAR_MOMENTUM + 1, NC = ((MD_GTO_MAX_ANGULAR_MOMENTUM + 1) * (MD_GTO_MAX_ANGULAR_MOMENTUM + 2)) / 2 };
    int    ijk[L1][NC][3];
    double nrm[L1][NC];
    for (int l = 0; l < L1; ++l) {
        for (uint32_t c = 0; c < md_gto_num_cart_ao((uint32_t)l); ++c) {
            md_gto_cart_ijk(&ijk[l][c][0], &ijk[l][c][1], &ijk[l][c][2], (uint32_t)l, c);
            nrm[l][c] = md_gto_cart_norm_factor(ijk[l][c][0], ijk[l][c][1], ijk[l][c][2]);
        }
    }

    for (uint32_t si = 0; si < basis->num_shells; ++si) {
        const md_gto_shell_t* a = &basis->shells[si];
        const md_double3_t         A = atom_coord_bohr[a->atom_idx];
        const double          Av[3] = { A.x, A.y, A.z };
        const uint32_t        la = a->l, na = md_gto_num_cart_ao(la);
        if (la > MD_GTO_MAX_ANGULAR_MOMENTUM) {
            MD_LOG_ERROR("Shell %u has angular momentum %u, above the supported %d", si, la, MD_GTO_MAX_ANGULAR_MOMENTUM);
            md_temp_end(temp);
            return 0;
        }

        for (uint32_t sj = si; sj < basis->num_shells; ++sj) {
            const md_gto_shell_t* b = &basis->shells[sj];
            const md_double3_t         B = atom_coord_bohr[b->atom_idx];
            const double          Bv[3] = { B.x, B.y, B.z };
            const uint32_t        lb = b->l, nb = md_gto_num_cart_ao(lb);
            if (lb > MD_GTO_MAX_ANGULAR_MOMENTUM) {
                MD_LOG_ERROR("Shell %u has angular momentum %u, above the supported %d", sj, lb, MD_GTO_MAX_ANGULAR_MOMENTUM);
                md_temp_end(temp);
                return 0;
            }
            double r2 = 0.0;
            for (int k = 0; k < 3; ++k) r2 += (Av[k] - Bv[k]) * (Av[k] - Bv[k]);

            double blk[NC][NC];
            md_memset(blk, 0, sizeof(blk));
            for (uint32_t p = 0; p < a->num_primitives; ++p) {
                const double alpha_a = basis->alpha[a->primitive_offset + p];
                const double coeff_a = basis->coeff[a->primitive_offset + p];
                for (uint32_t q = 0; q < b->num_primitives; ++q) {
                    const double alpha_b = basis->alpha[b->primitive_offset + q];
                    const double coeff_b = basis->coeff[b->primitive_offset + q];
                    const double pp = alpha_a + alpha_b;
                    const double mu = alpha_a * alpha_b / pp;
                    if (mu * r2 > 46.0) continue;
                    const double pref = coeff_a * coeff_b * exp(-mu * r2) * pow(MD_PI / pp, 1.5);
                    double E[3][L1][L1];
                    for (int k = 0; k < 3; ++k) {
                        const double P = (alpha_a * Av[k] + alpha_b * Bv[k]) / pp;
                        for (uint32_t i = 0; i <= la; ++i) {
                            for (uint32_t j = 0; j <= lb; ++j) {
                                E[k][i][j] = qm_overlap_1d((int)i, (int)j, P - Av[k], P - Bv[k], pp);
                            }
                        }
                    }
                    for (uint32_t ca = 0; ca < na; ++ca) {
                        const int* ea = ijk[la][ca];
                        for (uint32_t cb = 0; cb < nb; ++cb) {
                            const int* eb = ijk[lb][cb];
                            blk[ca][cb] += pref * E[0][ea[0]][eb[0]] * E[1][ea[1]][eb[1]] * E[2][ea[2]][eb[2]];
                        }
                    }
                }
            }
            for (uint32_t ca = 0; ca < na; ++ca) {
                for (uint32_t cb = 0; cb < nb; ++cb) {
                    const double v = blk[ca][cb] * nrm[la][ca] * nrm[lb][cb];
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
static size_t qm_overlap_provider(void* dst, size_t cap, const md_attributes_t* table, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    (void)slice;  // one whole {A,A} matrix, not indexed by anything a slice could fix
    (void)user_data;

    const md_attribute_t* coord = md_attributes_find(table, MD_STR_LIT("qm/atom/coordinate"));
    if (!coord) {
        MD_LOG_ERROR("'" MD_STR_FMT "' has no qm/atom/coordinate to integrate over", MD_STR_ARG(attr->path));
        return 0;
    }

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    size_t written = 0;
    md_gto_basis_t basis = {0};
    const size_t num_atoms = md_attribute_value_count(&coord->format);
    double*  xyz  = md_temp_alloc_array(temp, double,  num_atoms * 3);
    md_double3_t* bohr = md_temp_alloc_array(temp, md_double3_t, num_atoms);

    if (xyz && bohr && md_gto_basis_extract_attributes(&basis, table, temp_alloc)
        && md_attribute_extract_f64(xyz, num_atoms * 3, table, coord, md_attribute_slice_all(), md_unit_none()) == num_atoms * 3) {

        const double angstrom_to_bohr = 1.0 / 0.529177210903;
        for (size_t i = 0; i < num_atoms; ++i) {
            bohr[i] = (md_double3_t){ xyz[i * 3 + 0] * angstrom_to_bohr, xyz[i * 3 + 1] * angstrom_to_bohr, xyz[i * 3 + 2] * angstrom_to_bohr };
        }
        const size_t num_ao = md_gto_basis_num_ao(&basis);
        if (num_ao * num_ao != cap) {
            MD_LOG_ERROR("'" MD_STR_FMT "' was asked for %zu values and its basis spans %zu atomic orbitals", MD_STR_ARG(attr->path), cap, num_ao);
        } else if (md_qm_compute_overlap((double*)dst, &basis, bohr, num_atoms) == num_ao) {
            written = cap;
        }
    }

    md_temp_end(temp);
    return written;
}

bool md_qm_publish_overlap(md_system_t* sys) {
    MD_ASSERT(sys);

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

    if (num_ao == 0 || !md_attributes_find(&sys->attributes, MD_STR_LIT("qm/atom/coordinate"))) {
        return false;
    }

    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2, .shape = { (uint32_t)num_ao, (uint32_t)num_ao },
    };
    md_attribute_virtual_t virt = { .provider = qm_overlap_provider };
    return md_qm_publish_virtual(sys, MD_STR_LIT("basis/overlap"), MD_STR_LIT("AO Overlap"), md_unit_none(), format, &virt) != MD_ATTRIBUTE_INVALID;
}

// One spin's orbitals: their coefficient and occupation attributes and their shape. Nothing is read: the
// coefficients are {M,A} over every orbital, of which a density keeps the occupied (a fifth of a protein's
// 9636 in DEF2-SV(P)), and reading them all would be the whole matrix, 742 MB there, to keep 165.
static bool qm_spin_orbital_attrs(const md_attribute_t** coeff_attr, const md_attribute_t** occ_attr, size_t* num_mo, size_t* num_ao,
                                  const md_system_t* sys, bool beta) {
    const md_attribute_t* ca = md_attributes_find(&sys->attributes, beta ? MD_STR_LIT("orbital/beta/coefficient") : MD_STR_LIT("orbital/alpha/coefficient"));
    const md_attribute_t* oa = md_attributes_find(&sys->attributes, beta ? MD_STR_LIT("orbital/beta/occupation")  : MD_STR_LIT("orbital/alpha/occupation"));
    if (!ca || !oa || ca->format.rank != 2) return false;
    const size_t nm = ca->format.shape[0], na = ca->format.shape[1];
    if (md_attribute_element_count(&oa->format) != nm || nm == 0 || na == 0) return false;
    *coeff_attr = ca;
    *occ_attr = oa;
    *num_mo = nm;
    *num_ao = na;
    return true;
}

// One spin's orbital energies, hartree, into dst (num_mo of them). False when absent, inconsistent or not finite.
static bool qm_spin_orbital_energies(double* dst, size_t num_mo, const md_system_t* sys, bool beta) {
    const md_attribute_t* ea = md_attributes_find(&sys->attributes, beta ? MD_STR_LIT("orbital/beta/energy") : MD_STR_LIT("orbital/alpha/energy"));
    if (!ea || md_attribute_element_count(&ea->format) != num_mo) return false;
    if (md_attribute_extract_f64(dst, num_mo, &sys->attributes, ea, md_attribute_slice_all(), md_unit_hartree()) != num_mo) return false;
    for (size_t k = 0; k < num_mo; ++k) if (!isfinite(dst[k])) return false;
    return true;
}

// The next at most cap rows of 'coeff' ({M,A}) whose weight is nonzero, from orbital *io_k on, in order: into rows
// (na each), their weights into w (when not NULL), a run of consecutive ones at a time (md_attribute_slice_rows: one
// read of a block where the attribute is read from a file). *out_n the rows read, *io_k where the next look starts.
// False on a failed read.
static bool qm_read_weighted_block(double* rows, double* w, size_t* out_n, size_t* io_k, const md_attributes_t* table,
                                   const md_attribute_t* coeff, const double* weight, size_t num_mo, size_t na, size_t cap) {
    size_t n = 0, k = *io_k;
    while (k < num_mo && n < cap) {
        if (weight[k] == 0.0) {
            ++k;
            continue;
        }
        size_t e = k + 1;
        while (e < num_mo && weight[e] != 0.0 && e - k < cap - n) ++e;
        const size_t len = (e - k) * na;
        if (md_attribute_extract_f64(rows + n * na, len, table, coeff, md_attribute_slice_rows((uint32_t)k, (uint32_t)(e - k)), md_unit_none()) != len) {
            return false;
        }
        if (w) md_memcpy(w + n, weight + k, sizeof(double) * (e - k));
        n += e - k;
        k = e;
    }
    *io_k = k;
    *out_n = n;
    return true;
}

// Every row of nonzero weight, as qm_read_weighted_block
static bool qm_read_weighted_rows(double* rows, const md_attributes_t* table, const md_attribute_t* coeff, const double* weight,
                                  size_t num_mo, size_t na) {
    size_t k = 0, n = 0;
    return qm_read_weighted_block(rows, NULL, &n, &k, table, coeff, weight, num_mo, na, SIZE_MAX);
}

// The orbitals of the density sa alpha + sb beta, without their rows: each channel's coefficients and each orbital's
// weight in the density - alpha's (with beta's where the spins share orbitals: restricted, restricted open shell),
// then beta's (wb NULL when shared) - with the occupations and, with energy, the orbital energies beside them.
typedef struct qm_density_channels_t {
    const md_attribute_t* ca;           // alpha's coefficients {nma, na}, NULL without alpha
    const md_attribute_t* cb;           // beta's {nmb, na}, NULL without beta or where shared
    double *wa, *wb;                    // the weights
    double *oa, *ob, *ea, *eb;          // occupations, energies (NULL without energy)
    size_t nma, nmb, na;
    bool   shared;
} qm_density_channels_t;

// From temp. canonical: every occupation 0, 1 or 2 (within 1e-3), else false. False when the orbitals are absent or
// inconsistent.
static bool qm_density_channels(qm_density_channels_t* ch, const md_system_t* sys, double sa, double sb, bool canonical, bool energy,
                                md_temp_scope_t temp) {
    md_memset(ch, 0, sizeof(*ch));
    const md_attribute_t *caa = NULL, *oaa = NULL, *cab = NULL, *oab = NULL;
    size_t nma = 0, nmb = 0, naa = 0, nab = 0;
    bool have_a = false, have_b = false, beta_is_alpha = false;
    bool ok = true;
    if (sa != 0.0) ok = have_a = qm_spin_orbital_attrs(&caa, &oaa, &nma, &naa, sys, false);
    if (ok && sb != 0.0) {
        ok = have_b = qm_spin_orbital_attrs(&cab, &oab, &nmb, &nab, sys, true);
        if (!ok && sa == 0.0) {
            // no beta channel: a closed shell's beta density is its alpha density
            ok = have_b = beta_is_alpha = qm_spin_orbital_attrs(&cab, &oab, &nmb, &nab, sys, false);
        }
    }
    if (ok && have_a && have_b && naa != nab) ok = false;
    if (!ok) return false;
    // one set of orbitals for both spins: each orbital once, with its two weights summed
    const bool shared = have_a && have_b && nma == nmb && md_attribute_same_data(caa, cab);

    // the occupations and energies (small)
    double* oa = have_a ? (double*)md_temp_alloc(temp, sizeof(double) * nma) : NULL;
    double* ob = have_b ? (double*)md_temp_alloc(temp, sizeof(double) * nmb) : NULL;
    double* ea = have_a && energy ? (double*)md_temp_alloc(temp, sizeof(double) * nma) : NULL;
    double* eb = have_b && energy ? (double*)md_temp_alloc(temp, sizeof(double) * nmb) : NULL;
    if ((have_a && (!oa || md_attribute_extract_f64(oa, nma, &sys->attributes, oaa, md_attribute_slice_all(), md_unit_none()) != nma)) ||
        (have_b && (!ob || md_attribute_extract_f64(ob, nmb, &sys->attributes, oab, md_attribute_slice_all(), md_unit_none()) != nmb)) ||
        (have_a && energy && (!ea || !qm_spin_orbital_energies(ea, nma, sys, false))) ||
        (have_b && energy && (!eb || !qm_spin_orbital_energies(eb, nmb, sys, !beta_is_alpha)))) {
        return false;
    }
    if (canonical) {
        for (size_t k = 0; have_a && k < nma; ++k) if (fabs(oa[k] - round(oa[k])) > 1.0e-3 || oa[k] < -1.0e-3 || oa[k] > 2.001) return false;
        for (size_t k = 0; have_b && k < nmb; ++k) if (fabs(ob[k] - round(ob[k])) > 1.0e-3 || ob[k] < -1.0e-3 || ob[k] > 2.001) return false;
    }
    double* wa = have_a ? (double*)md_temp_alloc(temp, sizeof(double) * nma) : NULL;
    double* wb = have_b && !shared ? (double*)md_temp_alloc(temp, sizeof(double) * nmb) : NULL;
    if ((have_a && !wa) || (have_b && !shared && !wb)) return false;
    for (size_t k = 0; have_a && k < nma; ++k) wa[k] = sa * oa[k] + (shared ? sb * ob[k] : 0.0);
    for (size_t k = 0; wb && k < nmb; ++k) wb[k] = sb * ob[k];

    ch->ca = have_a ? caa : NULL;
    ch->cb = wb ? cab : NULL;
    ch->wa = wa;
    ch->wb = wb;
    ch->oa = oa;
    ch->ob = ob;
    ch->ea = ea;
    ch->eb = eb;
    ch->nma = have_a ? nma : 0;
    ch->nmb = wb ? nmb : 0;
    ch->na = have_a ? naa : nab;
    ch->shared = shared;
    return true;
}

// The orbitals of the density sa alpha + sb beta: rows {count, A}, weights {count} (the occupations times
// sa, sb), and with out_energy_weights the same with each occupation times |its orbital's energy|. Every
// array from 'alloc'. canonical: every occupation 0, 1 or 2 (within 1e-3), else nothing.
static size_t qm_orbital_products(double** out_rows, double** out_weights, double** out_energy_weights, size_t* out_num_ao,
                                  const md_system_t* sys, double sa, double sb, bool canonical, md_allocator_i* alloc) {
    const bool energy = out_energy_weights != NULL;
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    qm_density_channels_t ch;
    if (!qm_density_channels(&ch, sys, sa, sb, canonical, energy, temp)) {
        md_temp_end(temp);
        return 0;
    }
    const size_t na = ch.na;
    size_t count = 0;
    for (size_t k = 0; k < ch.nma; ++k) count += ch.wa[k] != 0.0;
    for (size_t k = 0; k < ch.nmb; ++k) count += ch.wb[k] != 0.0;
    double* rows = count ? (double*)md_alloc(alloc, sizeof(double) * count * na) : NULL;
    double* w    = count ? (double*)md_alloc(alloc, sizeof(double) * count) : NULL;
    double* we   = count && energy ? (double*)md_alloc(alloc, sizeof(double) * count) : NULL;
    bool ok = rows && w && (!energy || we);
    size_t n = 0;
    for (size_t k = 0; ok && k < ch.nma; ++k) {
        if (ch.wa[k] == 0.0) continue;
        // a shared orbital's two spins may differ in energy (restricted open shell): each its own
        if (energy) we[n] = sa * ch.oa[k] * fabs(ch.ea[k]) + (ch.shared ? sb * ch.ob[k] * fabs(ch.eb[k]) : 0.0);
        w[n++] = ch.wa[k];
    }
    const size_t na_rows = n;
    for (size_t k = 0; ok && k < ch.nmb; ++k) {
        if (ch.wb[k] == 0.0) continue;
        if (energy) we[n] = ch.wb[k] * fabs(ch.eb[k]);
        w[n++] = ch.wb[k];
    }
    // the rows of exactly the orbitals of nonzero weight, a run of consecutive ones at a time
    ok = ok && (!ch.ca || qm_read_weighted_rows(rows, &sys->attributes, ch.ca, ch.wa, ch.nma, na)) &&
               (!ch.cb || qm_read_weighted_rows(rows + na_rows * na, &sys->attributes, ch.cb, ch.wb, ch.nmb, na));
    md_temp_end(temp);
    if (!ok) {
        if (rows) md_free(alloc, rows, sizeof(double) * count * na);
        if (w)    md_free(alloc, w, sizeof(double) * count);
        if (we)   md_free(alloc, we, sizeof(double) * count);
        return 0;
    }
    *out_rows = rows;
    *out_weights = w;
    if (energy) *out_energy_weights = we;
    *out_num_ao = na;
    return count;
}

// The weights of alpha's and beta's occupations in an SCF density, by its name; false for any other path
static bool qm_scf_density_spins(double* sa, double* sb, md_str_t density_path) {
    *sa = *sb = 0.0;
    if      (md_str_eq(density_path, MD_STR_LIT("orbital/alpha/density")))      { *sa = 1.0; }
    else if (md_str_eq(density_path, MD_STR_LIT("orbital/beta/density")))       { *sb = 1.0; }
    else if (md_str_eq(density_path, MD_STR_LIT("orbital/total/density")))      { *sa = 1.0; *sb = 1.0; }
    else if (md_str_eq(density_path, MD_STR_LIT("orbital/difference/density"))) { *sa = 1.0; *sb = -1.0; }
    else return false;
    return true;
}

size_t md_qm_density_orbitals(double** out_rows, double** out_weights, size_t* out_num_ao, const md_system_t* sys,
                              md_str_t density_path, md_allocator_i* alloc) {
    if (!out_rows || !out_weights || !out_num_ao || !sys || !alloc) return 0;
    *out_rows = NULL;
    *out_weights = NULL;
    *out_num_ao = 0;
    double sa = 0.0, sb = 0.0;
    if (!qm_scf_density_spins(&sa, &sb, density_path)) return 0;
    return qm_orbital_products(out_rows, out_weights, NULL, out_num_ao, sys, sa, sb, false, alloc);
}

size_t md_qm_density_factors(double** out_rows, double** out_weights, size_t* out_num_ao, const md_system_t* sys,
                             md_str_t density_path, const md_attribute_slice_t* slice, md_allocator_i* alloc) {
    if (!out_rows || !out_weights || !out_num_ao || !sys || !alloc) return 0;
    *out_rows = NULL;
    *out_weights = NULL;
    *out_num_ao = 0;
    if (md_str_eq(density_path, MD_STR_LIT("orbital/alpha/density")) || md_str_eq(density_path, MD_STR_LIT("orbital/beta/density")) ||
        md_str_eq(density_path, MD_STR_LIT("orbital/total/density")) || md_str_eq(density_path, MD_STR_LIT("orbital/difference/density"))) {
        return md_qm_density_orbitals(out_rows, out_weights, out_num_ao, sys, density_path, alloc);
    }
    char fp[256], wp[256];
    if (md_snprintf(fp, sizeof(fp), MD_STR_FMT "/factor", MD_STR_ARG(density_path)) >= (int)sizeof(fp) ||
        md_snprintf(wp, sizeof(wp), MD_STR_FMT "/weight", MD_STR_ARG(density_path)) >= (int)sizeof(wp)) {
        return 0;
    }
    const md_attribute_t* fa = md_attributes_find(&sys->attributes, md_str_from_cstr(fp));
    const md_attribute_t* wa = md_attributes_find(&sys->attributes, md_str_from_cstr(wp));
    if (!fa || !wa) return 0;
    const md_attribute_slice_t sl = slice ? *slice : md_attribute_slice_all();
    md_attribute_format_t ff, wf;
    if (!md_attribute_slice_format(&ff, fa, sl) || !md_attribute_slice_format(&wf, wa, sl) || ff.rank != 2 || wf.rank != 1 ||
        ff.shape[0] != wf.shape[0]) {
        MD_LOG_ERROR("md_qm_density_factors: '%s' and '%s' sliced are not {R,A} and {R}", fp, wp);
        return 0;
    }
    const size_t R = ff.shape[0], A = ff.shape[1];
    if (R == 0 || A == 0) return 0;
    md_temp_scope_t temp = md_temp_begin();
    double* f = (double*)md_temp_alloc(temp, sizeof(double) * R * A);
    double* w = (double*)md_temp_alloc(temp, sizeof(double) * R);
    size_t count = 0;
    if (f && w && md_attribute_extract_f64(f, R * A, &sys->attributes, fa, sl, md_unit_none()) == R * A &&
        md_attribute_extract_f64(w, R, &sys->attributes, wa, sl, md_unit_none()) == R) {
        // the factors that carry something: a weight, and not all zero (a rank short of R leaves zero rows)
        for (size_t r = 0; r < R; ++r) {
            bool any = false;
            for (size_t mu = 0; mu < A && w[r] != 0.0 && !any; ++mu) any = f[r * A + mu] != 0.0;
            if (any) {
                if (count != r) {
                    md_memmove(f + count * A, f + r * A, sizeof(double) * A);
                    w[count] = w[r];
                }
                ++count;
            }
        }
        if (count) {
            *out_rows = (double*)md_alloc(alloc, sizeof(double) * count * A);
            *out_weights = (double*)md_alloc(alloc, sizeof(double) * count);
            if (*out_rows && *out_weights) {
                md_memcpy(*out_rows, f, sizeof(double) * count * A);
                md_memcpy(*out_weights, w, sizeof(double) * count);
                *out_num_ao = A;
            } else {
                if (*out_rows) md_free(alloc, *out_rows, sizeof(double) * count * A);
                if (*out_weights) md_free(alloc, *out_weights, sizeof(double) * count);
                *out_rows = NULL;
                *out_weights = NULL;
                count = 0;
            }
        }
    } else {
        MD_LOG_ERROR("md_qm_density_factors: could not read the %zu factors of '" MD_STR_FMT "'", R, MD_STR_ARG(density_path));
    }
    md_temp_end(temp);
    return count;
}

// Orbital rows held at once by md_qm_density_dmat_build by default: a block of them is read, added to the tiles, and the
// next read into the same place (256 MB: 320 orbitals at 100k AOs, all of a protein's occupied ones at 9.6k)
#define QM_DMAT_BLOCK_BYTES ((size_t)256 << 20)

size_t md_qm_density_dmat_build(struct md_gto_dmat_t* dm, const md_system_t* sys, md_str_t density_path, const md_attribute_slice_t* slice,
                                size_t max_rows, const struct md_executor_t* executor) {
    if (!dm || !sys) return 0;
    md_allocator_i* heap = md_get_heap_allocator();
    double sa = 0.0, sb = 0.0;
    if (!qm_scf_density_spins(&sa, &sb, density_path)) {
        // a density published as its factors: a state's of them at once, as md_qm_density_factors reads them
        double *rows = NULL, *w = NULL;
        size_t na = 0;
        const size_t count = md_qm_density_factors(&rows, &w, &na, sys, density_path, slice, heap);
        const bool ok = count > 0 && na == dm->num_ao && md_gto_dmat_build(dm, rows, na, 1, w, count, executor);
        if (count && na != dm->num_ao) {
            MD_LOG_ERROR("md_qm_density_dmat_build: '" MD_STR_FMT "' is over %zu AOs, the tiles %u", MD_STR_ARG(density_path), na, dm->num_ao);
        }
        if (rows) md_free(heap, rows, sizeof(double) * count * na);
        if (w)    md_free(heap, w, sizeof(double) * count);
        return ok ? count : 0;
    }

    // an SCF density: its orbitals a block at a time, alpha's then beta's, in the order md_qm_density_orbitals gives them
    md_temp_scope_t temp = md_temp_begin();
    qm_density_channels_t ch;
    if (!qm_density_channels(&ch, sys, sa, sb, false, false, temp)) {
        md_temp_end(temp);
        return 0;
    }
    if (ch.na != dm->num_ao) {
        MD_LOG_ERROR("md_qm_density_dmat_build: the orbitals are over %zu AOs, the tiles %u", ch.na, dm->num_ao);
        md_temp_end(temp);
        return 0;
    }
    size_t count = 0;
    for (size_t k = 0; k < ch.nma; ++k) count += ch.wa[k] != 0.0;
    for (size_t k = 0; k < ch.nmb; ++k) count += ch.wb[k] != 0.0;
    const size_t rows_max = max_rows ? max_rows : MD_MAX(QM_DMAT_BLOCK_BYTES / (sizeof(double) * ch.na), (size_t)128);
    const size_t cap = MD_MIN(rows_max, MD_MAX(count, (size_t)1));
    double* rows = count ? (double*)md_alloc(heap, sizeof(double) * cap * ch.na) : NULL;
    double* w    = count ? (double*)md_alloc(heap, sizeof(double) * cap) : NULL;
    bool ok = rows && w && md_gto_dmat_clear(dm);
    for (int spin = 0; ok && spin < 2; ++spin) {
        const md_attribute_t* c = spin == 0 ? ch.ca : ch.cb;
        const double* weight    = spin == 0 ? ch.wa : ch.wb;
        const size_t  nm        = spin == 0 ? ch.nma : ch.nmb;
        for (size_t k = 0; ok && c && k < nm;) {
            size_t n = 0;
            ok = qm_read_weighted_block(rows, w, &n, &k, &sys->attributes, c, weight, nm, ch.na, cap);
            if (ok && n) ok = md_gto_dmat_add(dm, rows, ch.na, 1, w, n, executor);
        }
    }
    if (rows) md_free(heap, rows, sizeof(double) * cap * ch.na);
    if (w)    md_free(heap, w, sizeof(double) * cap);
    md_temp_end(temp);
    if (!ok && count) MD_LOG_ERROR("md_qm_density_dmat_build: '" MD_STR_FMT "': out of memory, a failed read or cancelled", MD_STR_ARG(density_path));
    return ok ? count : 0;
}

bool md_qm_density_factors_available(const md_system_t* sys, md_str_t density_path) {
    if (!sys) return false;
    // an SCF density: its spins' orbitals (alpha alone for a file with one channel; both for beta, the total and the
    // difference, beta perhaps a second name for alpha's)
    const bool alpha = md_str_eq(density_path, MD_STR_LIT("orbital/alpha/density"));
    const bool beta  = md_str_eq(density_path, MD_STR_LIT("orbital/beta/density"));
    if (alpha || beta || md_str_eq(density_path, MD_STR_LIT("orbital/total/density")) ||
        md_str_eq(density_path, MD_STR_LIT("orbital/difference/density"))) {
        const md_attribute_t *ca = NULL, *oa = NULL;
        size_t nm = 0, na = 0, nb = 0;
        const bool a = qm_spin_orbital_attrs(&ca, &oa, &nm, &na, sys, false);
        const bool b = qm_spin_orbital_attrs(&ca, &oa, &nm, &nb, sys, true);
        return alpha ? a : beta ? b : (a && b && na == nb);
    }
    char fp[256];
    if (md_snprintf(fp, sizeof(fp), MD_STR_FMT "/factor", MD_STR_ARG(density_path)) >= (int)sizeof(fp)) return false;
    return md_attributes_find(&sys->attributes, md_str_from_cstr(fp)) != NULL;
}

size_t md_qm_alie_orbitals(double** out_rows, double** out_weights, double** out_energy_weights, size_t* out_num_ao,
                           const md_system_t* sys, md_allocator_i* alloc) {
    if (!out_rows || !out_weights || !out_energy_weights || !out_num_ao || !sys || !alloc) return 0;
    *out_rows = NULL;
    *out_weights = NULL;
    *out_energy_weights = NULL;
    *out_num_ao = 0;
    return qm_orbital_products(out_rows, out_weights, out_energy_weights, out_num_ao, sys, 1.0, 1.0, true, alloc);
}

bool md_qm_alie_available(const md_system_t* sys) {
    if (!sys) return false;
    const md_attribute_t *ca = NULL, *oa = NULL;
    size_t nm = 0, na = 0;
    if (!qm_spin_orbital_attrs(&ca, &oa, &nm, &na, sys, false)) return false;
    const md_attribute_t* ea = md_attributes_find(&sys->attributes, MD_STR_LIT("orbital/alpha/energy"));
    if (!ea || md_attribute_element_count(&ea->format) != nm) return false;
    // the occupations of canonical orbitals: natural orbitals carry fractional ones, and no energies to weight by
    md_temp_scope_t temp = md_temp_begin();
    double* occ = (double*)md_temp_alloc(temp, sizeof(double) * nm);
    bool ok = occ && md_attribute_extract_f64(occ, nm, &sys->attributes, oa, md_attribute_slice_all(), md_unit_none()) == nm;
    for (size_t k = 0; ok && k < nm; ++k) ok = fabs(occ[k] - round(occ[k])) <= 1.0e-3 && occ[k] >= -1.0e-3 && occ[k] <= 2.001;
    md_temp_end(temp);
    return ok;
}

// A spin channel's occupations and energies: n {M}, e {M} (hartree) from 'alloc', with the coefficient attribute
// and its shape. False when anything is missing, inconsistent or not finite. beta without a beta channel is alpha.
static bool qm_spin_channel(const md_attribute_t** coeff, double** occ, double** energy, size_t* num_mo, size_t* num_ao,
                            const md_system_t* sys, bool beta, md_temp_scope_t temp) {
    const md_attribute_t* oa = NULL;
    if (!qm_spin_orbital_attrs(coeff, &oa, num_mo, num_ao, sys, beta)) {
        if (!beta || !qm_spin_orbital_attrs(coeff, &oa, num_mo, num_ao, sys, false)) return false;
        beta = false;
    }
    *occ    = (double*)md_temp_alloc(temp, sizeof(double) * *num_mo);
    *energy = (double*)md_temp_alloc(temp, sizeof(double) * *num_mo);
    return *occ && *energy &&
           md_attribute_extract_f64(*occ, *num_mo, &sys->attributes, oa, md_attribute_slice_all(), md_unit_none()) == *num_mo &&
           qm_spin_orbital_energies(*energy, *num_mo, sys, beta);
}

// A spin orbital's share of E(r): (1 - n) e when it is empty, or partly, and bound (e < 0); else nothing
static inline double qm_attachment_weight(double n, double e) {
    return (e < 0.0 && n < 1.0 - 1.0e-3) ? (1.0 - n) * e : 0.0;
}

// Per spin occupations of canonical orbitals: each 0 or 1 (the readers publish a restricted orbital's 2 as 1 per spin)
static bool qm_spin_occupations_canonical(const double* n, size_t num) {
    for (size_t k = 0; k < num; ++k) if (fabs(n[k] - round(n[k])) > 1.0e-3 || n[k] < -1.0e-3 || n[k] > 1.001) return false;
    return true;
}

size_t md_qm_attachment_orbitals(double** out_rows, double** out_weights, size_t* out_num_ao, const md_system_t* sys, md_allocator_i* alloc) {
    if (!out_rows || !out_weights || !out_num_ao || !sys || !alloc) return 0;
    *out_rows = NULL;
    *out_weights = NULL;
    *out_num_ao = 0;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    const md_attribute_t *ca = NULL, *cb = NULL;
    double *na = NULL, *nb = NULL, *ea = NULL, *eb = NULL;
    size_t ma = 0, mb = 0, aa = 0, ab = 0;
    if (!qm_spin_channel(&ca, &na, &ea, &ma, &aa, sys, false, temp) || !qm_spin_channel(&cb, &nb, &eb, &mb, &ab, sys, true, temp) ||
        aa != ab || !qm_spin_occupations_canonical(na, ma) || !qm_spin_occupations_canonical(nb, mb)) {
        md_temp_end(temp);
        return 0;
    }
    // one set of orbitals for both spins (restricted, restricted open shell): each orbital once, its two spins' shares summed
    const bool shared = ma == mb && md_attribute_same_data(ca, cb);
    size_t count = 0;
    for (size_t k = 0; k < ma; ++k) count += (qm_attachment_weight(na[k], ea[k]) + (shared ? qm_attachment_weight(nb[k], eb[k]) : 0.0)) != 0.0;
    for (size_t k = 0; !shared && k < mb; ++k) count += qm_attachment_weight(nb[k], eb[k]) != 0.0;
    if (count == 0) {
        md_temp_end(temp);
        return 0;
    }
    const size_t nao = aa;
    double* rows = (double*)md_alloc(alloc, sizeof(double) * count * nao);
    double* w    = (double*)md_alloc(alloc, sizeof(double) * count);
    double* wa   = (double*)md_temp_alloc(temp, sizeof(double) * ma);
    double* wb   = shared ? NULL : (double*)md_temp_alloc(temp, sizeof(double) * mb);
    bool ok = rows && w && wa && (shared || wb);
    size_t n = 0;
    for (size_t k = 0; ok && k < ma; ++k) {
        wa[k] = qm_attachment_weight(na[k], ea[k]) + (shared ? qm_attachment_weight(nb[k], eb[k]) : 0.0);
        if (wa[k] != 0.0) w[n++] = wa[k];
    }
    const size_t na_rows = n;
    for (size_t k = 0; ok && wb && k < mb; ++k) {
        wb[k] = qm_attachment_weight(nb[k], eb[k]);
        if (wb[k] != 0.0) w[n++] = wb[k];
    }
    // the rows of exactly the orbitals of nonzero weight, a run of consecutive ones at a time
    ok = ok && qm_read_weighted_rows(rows, &sys->attributes, ca, wa, ma, nao) &&
               (!wb || qm_read_weighted_rows(rows + na_rows * nao, &sys->attributes, cb, wb, mb, nao));
    md_temp_end(temp);
    if (!ok) {
        if (rows) md_free(alloc, rows, sizeof(double) * count * nao);
        if (w)    md_free(alloc, w, sizeof(double) * count);
        return 0;
    }
    *out_rows = rows;
    *out_weights = w;
    *out_num_ao = nao;
    return count;
}

bool md_qm_attachment_available(const md_system_t* sys) {
    if (!sys) return false;
    md_temp_scope_t temp = md_temp_begin();
    const md_attribute_t *ca = NULL, *cb = NULL;
    double *na = NULL, *nb = NULL, *ea = NULL, *eb = NULL;
    size_t ma = 0, mb = 0, aa = 0, ab = 0;
    bool ok = qm_spin_channel(&ca, &na, &ea, &ma, &aa, sys, false, temp) && qm_spin_channel(&cb, &nb, &eb, &mb, &ab, sys, true, temp) &&
              aa == ab && qm_spin_occupations_canonical(na, ma) && qm_spin_occupations_canonical(nb, mb);
    bool any = false;
    for (size_t k = 0; ok && !any && k < ma; ++k) any = qm_attachment_weight(na[k], ea[k]) != 0.0;
    for (size_t k = 0; ok && !any && k < mb; ++k) any = qm_attachment_weight(nb[k], eb[k]) != 0.0;
    md_temp_end(temp);
    return ok && any;
}

// ---------------------------------------------------------------------------
// Kinetic energy density, ELF and LOL
// ---------------------------------------------------------------------------

size_t md_qm_kinetic_orbitals(double** out_rows, double** out_weights, size_t* out_num_ao, md_gto_basis_t* out_basis,
                              const md_system_t* sys, md_str_t density_path, md_allocator_i* alloc) {
    if (!out_rows || !out_weights || !out_num_ao || !out_basis || !sys || !alloc) return 0;
    *out_rows = NULL;
    *out_weights = NULL;
    *out_num_ao = 0;
    md_memset(out_basis, 0, sizeof(*out_basis));
    if (!md_str_eq(density_path, MD_STR_LIT("orbital/alpha/density")) && !md_str_eq(density_path, MD_STR_LIT("orbital/beta/density")) &&
        !md_str_eq(density_path, MD_STR_LIT("orbital/total/density"))) {
        return 0;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* talloc = md_temp_allocator(temp);
    double *orb = NULL, *occ = NULL;
    size_t na = 0;
    const size_t num_orb = md_qm_density_orbitals(&orb, &occ, &na, sys, density_path, talloc);
    md_gto_basis_t basis = {0};
    size_t nd = 0;  // the derivative basis' AOs; stays 0 for a basis without one
    if (num_orb == 0 || !md_gto_basis_extract_attributes(&basis, &sys->attributes, talloc) || md_gto_basis_num_ao(&basis) != na ||
        (md_gto_orbital_gradients(NULL, &nd, NULL, 0, 0, &basis), nd == 0)) {
        md_temp_end(temp);
        return 0;
    }
    const size_t count = 3 * num_orb;
    double* rows = (double*)md_alloc(alloc, sizeof(double) * count * nd);
    double* w    = (double*)md_alloc(alloc, sizeof(double) * count);
    const bool ok = rows && w && md_gto_orbital_gradients(rows, &nd, orb, na, num_orb, &basis) == count &&
                    md_gto_basis_derivative(out_basis, &basis, alloc);
    md_temp_end(temp);
    if (!ok) {
        if (rows) md_free(alloc, rows, sizeof(double) * count * nd);
        if (w)    md_free(alloc, w, sizeof(double) * count);
        md_gto_basis_free(out_basis, alloc);
        return 0;
    }
    for (size_t k = 0; k < num_orb; ++k) {
        w[3 * k] = w[3 * k + 1] = w[3 * k + 2] = 0.5 * occ[k];
    }
    *out_rows = rows;
    *out_weights = w;
    *out_num_ao = nd;
    return count;
}

bool md_qm_kinetic_available(const md_system_t* sys) {
    if (!sys) return false;
    const md_attribute_t *ca = NULL, *oa = NULL;
    size_t nm = 0, na = 0;
    if (!qm_spin_orbital_attrs(&ca, &oa, &nm, &na, sys, false)) return false;
    const md_attribute_t* angl = md_attributes_find(&sys->attributes, MD_STR_LIT("basis/shell/angular_momentum"));
    if (!angl || angl->format.type != MD_ATTRIBUTE_TYPE_U32) return false;
    const size_t n = md_attribute_element_count(&angl->format);
    const uint32_t* l = (const uint32_t*)md_attribute_view(angl, MD_ATTRIBUTE_TYPE_U32, 1, 1);
    bool ok = n > 0 && l;
    for (size_t k = 0; ok && k < n; ++k) ok = l[k] < MD_GTO_MAX_ANGULAR_MOMENTUM;   // its derivative one more
    return ok;
}

// 3/10 (6 pi^2)^(2/3): the uniform electron gas' kinetic energy density of one spin, over rho_s^(5/3)
#define QM_ELF_CF_SPIN 4.557799872345596

static inline double qm_elf_pauli(double rho, double g2, double tau) {
    if (!(rho > 0.0)) return 0.0;
    const double d = tau - g2 / (8.0 * rho);
    return d > 0.0 ? d : 0.0;
}

double md_qm_elf(double rho_a, double g2_a, double tau_a, double rho_b, double g2_b, double tau_b, double rho_min, double d_min) {
    const double rho = rho_a + rho_b;
    if (!(rho > 0.0) || rho < rho_min) return 0.0;
    const double d  = qm_elf_pauli(rho_a, g2_a, tau_a) + qm_elf_pauli(rho_b, g2_b, tau_b) + (d_min > 0.0 ? d_min : 0.0);
    const double d0 = QM_ELF_CF_SPIN * ((rho_a > 0.0 ? pow(rho_a, 5.0 / 3.0) : 0.0) + (rho_b > 0.0 ? pow(rho_b, 5.0 / 3.0) : 0.0));
    if (!(d0 > 0.0)) return 0.0;
    const double chi = d / d0;
    return 1.0 / (1.0 + chi * chi);
}

double md_qm_elf_closed_shell(double rho, double g2, double tau, double rho_min, double d_min) {
    // each spin half of it: rho/2, |grad rho|^2 / 4, tau/2 (d_min is of D, both spins)
    return md_qm_elf(0.5 * rho, 0.25 * g2, 0.5 * tau, 0.5 * rho, 0.25 * g2, 0.5 * tau, rho_min, d_min);
}

double md_qm_elf_noise(double rho_a, double tau_a, double rho_b, double tau_b) {
    const double d0 = QM_ELF_CF_SPIN * ((rho_a > 0.0 ? pow(rho_a, 5.0 / 3.0) : 0.0) + (rho_b > 0.0 ? pow(rho_b, 5.0 / 3.0) : 0.0));
    const double tau = tau_a + tau_b;
    if (!(d0 > 0.0)) return INFINITY;
    return 5.9604644775390625e-8 * (tau > 0.0 ? tau : 0.0) / d0;   // float's unit roundoff, 2^-24
}

double md_qm_lol(double rho_a, double tau_a, double rho_b, double tau_b, double rho_min) {
    const double rho = rho_a + rho_b;
    if (!(rho > 0.0) || rho < rho_min) return 0.0;
    const double d0  = QM_ELF_CF_SPIN * ((rho_a > 0.0 ? pow(rho_a, 5.0 / 3.0) : 0.0) + (rho_b > 0.0 ? pow(rho_b, 5.0 / 3.0) : 0.0));
    const double tau = tau_a + tau_b;
    if (!(tau > 0.0)) return 1.0;   // t = D0 / tau without bound
    const double t = d0 / tau;
    return t / (1.0 + t);
}

double md_qm_lol_closed_shell(double rho, double tau, double rho_min) {
    return md_qm_lol(0.5 * rho, 0.5 * tau, 0.5 * rho, 0.5 * tau, rho_min);
}

double md_qm_density_scale_length(double rho, double g2, double rho_min, double length_max) {
    if (!(rho > 0.0) || rho < rho_min) return 0.0;
    const double cap = length_max > 0.0 ? MD_MIN(length_max, (double)FLT_MAX) : (double)FLT_MAX;
    const double g = sqrt(g2 > 0.0 ? g2 : 0.0);
    return g * cap > rho ? rho / g : cap;
}

bool md_qm_density_available(const md_system_t* sys) {
    if (!sys) return false;
    const md_attribute_t *ca = NULL, *oa = NULL;
    size_t nm = 0, na = 0;
    return qm_spin_orbital_attrs(&ca, &oa, &nm, &na, sys, false);
}

// One spin's (or a closed shell's total) density and kinetic energy density at the points: rho, |grad rho|^2, tau
static bool qm_elf_terms(double* rho, double* g2, double* tau, const float* xyz, size_t count, size_t stride, const md_system_t* sys,
                         md_str_t density_path, const float* atom_xyz, size_t atom_xyz_stride, md_allocator_i* alloc) {
    double *rr = NULL, *rw = NULL, *kr = NULL, *kw = NULL;
    size_t rna = 0, kna = 0;
    md_gto_basis_t basis = {0}, dbasis = {0};
    const size_t rn = md_qm_density_orbitals(&rr, &rw, &rna, sys, density_path, alloc);
    const size_t kn = rn ? md_qm_kinetic_orbitals(&kr, &kw, &kna, &dbasis, sys, density_path, alloc) : 0;
    bool ok = rn && kn && md_gto_basis_extract_attributes(&basis, &sys->attributes, alloc);
    md_gto_iso_field_t* fr = ok ? md_gto_iso_field_create_products(&basis, atom_xyz, atom_xyz_stride, rr, NULL, rw, rn, alloc) : NULL;
    md_gto_iso_field_t* fk = ok ? md_gto_iso_field_create_products(&dbasis, atom_xyz, atom_xyz_stride, kr, NULL, kw, kn, alloc) : NULL;
    ok = fr && fk;
    double* grad = ok ? (double*)md_alloc(alloc, sizeof(double) * 3 * count) : NULL;
    if (ok && grad) {
        md_gto_iso_field_evaluate(fr, xyz, count, stride, rho, grad);
        md_gto_iso_field_evaluate(fk, xyz, count, stride, tau, NULL);
        for (size_t i = 0; i < count; ++i) {
            g2[i] = grad[3 * i] * grad[3 * i] + grad[3 * i + 1] * grad[3 * i + 1] + grad[3 * i + 2] * grad[3 * i + 2];
        }
    }
    ok = ok && grad;
    if (fr) md_gto_iso_field_destroy(fr);
    if (fk) md_gto_iso_field_destroy(fk);
    return ok;
}

size_t md_qm_elf_evaluate(double* out_elf, double* out_lol, const float* xyz, size_t count, size_t stride,
                          const md_system_t* sys, const float* atom_xyz, size_t atom_xyz_stride, double rho_min, double d_min) {
    if ((!out_elf && !out_lol) || !xyz || count == 0 || !sys || !atom_xyz) return 0;
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    // A closed shell: the same orbitals and occupations in both spins, which are then one density
    double *ar = NULL, *aw = NULL, *br = NULL, *bw = NULL;
    size_t ana = 0, bna = 0;
    const size_t an = md_qm_density_orbitals(&ar, &aw, &ana, sys, MD_STR_LIT("orbital/alpha/density"), alloc);
    const size_t bn = md_qm_density_orbitals(&br, &bw, &bna, sys, MD_STR_LIT("orbital/beta/density"), alloc);
    const bool closed = an == bn && ana == bna && (an == 0 || (memcmp(ar, br, sizeof(double) * an * ana) == 0 && memcmp(aw, bw, sizeof(double) * an) == 0));

    bool ok = an > 0 || bn > 0;
    double* t = ok ? (double*)md_alloc(alloc, sizeof(double) * count * 6) : NULL;
    ok = ok && t;
    double *rho_a = t, *g2_a = t + count, *tau_a = t + 2 * count, *rho_b = t + 3 * count, *g2_b = t + 4 * count, *tau_b = t + 5 * count;
    if (ok && closed) {
        ok = qm_elf_terms(rho_a, g2_a, tau_a, xyz, count, stride, sys, MD_STR_LIT("orbital/total/density"), atom_xyz, atom_xyz_stride, alloc);
        for (size_t i = 0; ok && i < count; ++i) {
            if (out_elf) out_elf[i] = md_qm_elf_closed_shell(rho_a[i], g2_a[i], tau_a[i], rho_min, d_min);
            if (out_lol) out_lol[i] = md_qm_lol_closed_shell(rho_a[i], tau_a[i], rho_min);
        }
    } else if (ok) {
        // a spin without occupied orbitals has no density: zero, not a failure
        md_memset(t, 0, sizeof(double) * count * 6);
        if (an) ok = qm_elf_terms(rho_a, g2_a, tau_a, xyz, count, stride, sys, MD_STR_LIT("orbital/alpha/density"), atom_xyz, atom_xyz_stride, alloc);
        if (ok && bn) ok = qm_elf_terms(rho_b, g2_b, tau_b, xyz, count, stride, sys, MD_STR_LIT("orbital/beta/density"), atom_xyz, atom_xyz_stride, alloc);
        for (size_t i = 0; ok && i < count; ++i) {
            if (out_elf) out_elf[i] = md_qm_elf(rho_a[i], g2_a[i], tau_a[i], rho_b[i], g2_b[i], tau_b[i], rho_min, d_min);
            if (out_lol) out_lol[i] = md_qm_lol(rho_a[i], tau_a[i], rho_b[i], tau_b[i], rho_min);
        }
    }
    md_temp_end(temp);
    return ok ? count : 0;
}
