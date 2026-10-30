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

// The matrix 'slice' narrows 'attr' to: N, and whether it is packed. 0 when it is not one symmetric matrix.
static size_t qm_symmetric_dim(bool* out_packed, const md_attribute_t* attr, md_attribute_slice_t slice) {
    md_attribute_format_t format = {0};
    if (!attr || !md_attribute_slice_format(&format, attr, slice) || format.components != 1) {
        return 0;
    }
    const bool packed = (attr->flags & MD_ATTRIBUTE_FLAG_PACKED_SYMMETRIC) != 0;
    size_t n = 0;
    if (packed) {
        n = format.rank == 1 ? md_attribute_packed_symmetric_dim(format.shape[0]) : 0;
    } else if (format.rank == 2 && format.shape[0] == format.shape[1]) {
        n = format.shape[0];
    }
    if (n == 0) {
        MD_LOG_ERROR("'" STR_FMT "' does not slice down to one symmetric matrix", STR_ARG(attr->path));
        return 0;
    }
    if (out_packed) *out_packed = packed;
    return n;
}

size_t md_qm_extract_packed_symmetric_f32(float* dst, size_t cap, const md_attribute_t* attr, md_attribute_slice_t slice) {
    bool packed = false;
    const size_t n = qm_symmetric_dim(&packed, attr, slice);
    if (n == 0 || !dst) {
        return n;
    }
    const size_t len = n * (n + 1) / 2;
    if (cap < len) {
        MD_LOG_ERROR("'" STR_FMT "': room for %zu values, the packed matrix has %zu", STR_ARG(attr->path), cap, len);
        return 0;
    }

    if (packed) {
        return md_attribute_extract_f32(dst, len, attr, slice, md_unit_none()) == len ? n : 0;
    }

    md_temp_scope_t temp = md_temp_begin();
    double* full = md_temp_alloc_array(temp, double, n * n);
    size_t result = 0;
    if (full && md_attribute_extract_f64(full, n * n, attr, slice, md_unit_none()) == n * n) {
        size_t k = 0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                dst[k++] = (float)full[i * n + j];
            }
        }
        result = n;
    }
    md_temp_end(temp);
    return result;
}

size_t md_qm_extract_symmetric_f64(double* dst, size_t cap, const md_attribute_t* attr, md_attribute_slice_t slice) {
    bool packed = false;
    const size_t n = qm_symmetric_dim(&packed, attr, slice);
    if (n == 0 || !dst) {
        return n;
    }
    if (cap < n * n) {
        MD_LOG_ERROR("'" STR_FMT "': room for %zu values, the matrix has %zu", STR_ARG(attr->path), cap, n * n);
        return 0;
    }

    if (!packed) {
        return md_attribute_extract_f64(dst, n * n, attr, slice, md_unit_none()) == n * n ? n : 0;
    }

    const size_t len = n * (n + 1) / 2;
    md_temp_scope_t temp = md_temp_begin();
    float* tri = md_temp_alloc_array(temp, float, len);
    size_t result = 0;
    if (tri && md_attribute_extract_f32(tri, len, attr, slice, md_unit_none()) == len) {
        size_t k = 0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                dst[i * n + j] = dst[j * n + i] = (double)tri[k++];
            }
        }
        result = n;
    }
    md_temp_end(temp);
    return result;
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
    int len = md_snprintf(buf, cap, STR_FMT "/" STR_FMT, STR_ARG(group), STR_ARG(name));
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
    if (md_gto_sph_to_cart_orbitals(dst, src, n_sph, 1, NULL, num_mo, basis) != num_mo) {
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

bool md_qm_publish_atoms(md_system_t* sys, const uint8_t atomic_number[], const double nuclear_charge[], const dvec3_t coord_angstrom[], size_t count) {
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
            q_id = md_qm_publish_series(sys, STR_LIT("qm/atom/nuclear_charge"), STR_LIT("Nuclear Charge"),
                                        md_unit_elementary_charge(), charge, count);
        }
        md_temp_end(temp);
    }

    // Angstrom, matching the system's own coordinates rather than the bohr the evaluator works in.
    // This is the geometry the CALCULATION was run at, which need not be where the system's atoms
    // are now - a trajectory frame or an optimisation step moves them.
    md_attribute_id_t xyz_id = md_qm_publish_vec3_series(sys, STR_LIT("qm/atom/coordinate"), STR_LIT("Coordinate"),
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
static size_t qm_overlap_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
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
        && md_attribute_extract_f64(xyz, num_atoms * 3, coord, md_attribute_slice_all(), md_unit_none()) == num_atoms * 3) {

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
// sum, fractional occupations included, so nothing is gained by storing it separately. Blocked and
// over the cores (md_gto_density_from_orbitals): a large basis is ~num_occ num_ao^2 / 2 multiply-adds.
static bool qm_build_density(double* out, const double* coeff, const double* occ, size_t num_mo, size_t num_ao) {
    return md_gto_density_from_orbitals(out, coeff, num_ao, 1, occ, num_mo, num_ao);
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
        ok = md_attribute_extract_f64(coeff_data, num_mo * num_ao, coeff, md_attribute_slice_all(), md_unit_none()) == num_mo * num_ao
          && md_attribute_extract_f64(occ_data, num_mo, occ, md_attribute_slice_all(), md_unit_none()) == num_mo;
        if (ok) {
            ok = qm_build_density((double*)dst, coeff_data, occ_data, num_mo, num_ao);
        }
    }

    md_temp_end(temp);
    return ok ? cap : 0;
}

static size_t qm_alpha_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    (void)slice;  // one whole {A,A} matrix, not indexed by anything a slice could fix
    return qm_density_provide(dst, cap, attr, user_data, false);
}

static size_t qm_beta_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    (void)slice;
    return qm_density_provide(dst, cap, attr, user_data, true);
}

// alpha +/- beta. Both halves are themselves virtual, so this is a derivation over derivations,
// legal because the dependency graph stays acyclic. A restricted calculation gets it for free:
// beta is an ALIAS of alpha there, so the total is twice alpha and the difference is zero, and alpha
// is read once (or, for the difference, not at all) rather than once per name.
static size_t qm_density_combine(void* dst, size_t cap, const md_attribute_t* attr, void* user_data, double beta_scale) {
    md_system_t* sys = (md_system_t*)user_data;

    const md_attribute_t* alpha = md_attributes_find(&sys->attributes, STR_LIT("orbital/alpha/density"));
    const md_attribute_t* beta  = md_attributes_find(&sys->attributes, STR_LIT("orbital/beta/density"));
    if (!alpha || !beta) {
        MD_LOG_ERROR("'" STR_FMT "' is missing a spin density it combines", STR_ARG(attr->path));
        return 0;
    }

    double* out = (double*)dst;
    if (md_attribute_same_data(alpha, beta)) {
        const double scale = 1.0 + beta_scale;
        if (scale == 0.0) {
            MEMSET(out, 0, sizeof(double) * cap);
            return cap;
        }
        if (md_attribute_extract_f64(out, cap, alpha, md_attribute_slice_all(), md_unit_none()) != cap) {
            MD_LOG_ERROR("'" STR_FMT "': could not read '" STR_FMT "'", STR_ARG(attr->path), STR_ARG(alpha->path));
            return 0;
        }
        if (scale != 1.0) {
            for (size_t i = 0; i < cap; ++i) out[i] *= scale;
        }
        return cap;
    }

    md_temp_scope_t temp = md_temp_begin();
    double* beta_data = md_temp_alloc_array(temp, double, cap);

    bool ok = beta_data != NULL;
    if (!ok) {
        MD_LOG_ERROR("'" STR_FMT "': failed to allocate %zu doubles of scratch", STR_ARG(attr->path), cap);
    }
    if (ok && md_attribute_extract_f64(out, cap, alpha, md_attribute_slice_all(), md_unit_none()) != cap) {
        MD_LOG_ERROR("'" STR_FMT "': could not read '" STR_FMT "'", STR_ARG(attr->path), STR_ARG(alpha->path));
        ok = false;
    }
    if (ok && md_attribute_extract_f64(beta_data, cap, beta, md_attribute_slice_all(), md_unit_none()) != cap) {
        MD_LOG_ERROR("'" STR_FMT "': could not read '" STR_FMT "'", STR_ARG(attr->path), STR_ARG(beta->path));
        ok = false;
    }
    if (ok) {
        for (size_t i = 0; i < cap; ++i) {
            out[i] += beta_scale * beta_data[i];
        }
    }

    md_temp_end(temp);
    return ok ? cap : 0;
}

static size_t qm_total_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    (void)slice;
    return qm_density_combine(dst, cap, attr, user_data, 1.0);
}

static size_t qm_difference_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
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

    return md_qm_publish_orbital_density_combinations(sys);
}

// One spin's orbitals and occupations: C {M,A} and occ {M} into 'alloc'. False when absent or inconsistent.
// One spin's orbitals: their coefficient and occupation attributes and their shape. Nothing is read: the
// coefficients are {M,A} over every orbital, of which a density keeps the occupied (a fifth of a protein's
// 9636 in DEF2-SV(P)), and reading them all would be the whole matrix, 742 MB there, to keep 165.
static bool qm_spin_orbital_attrs(const md_attribute_t** coeff_attr, const md_attribute_t** occ_attr, size_t* num_mo, size_t* num_ao,
                                  const md_system_t* sys, bool beta) {
    const md_attribute_t* ca = md_attributes_find(&sys->attributes, beta ? STR_LIT("orbital/beta/coefficient") : STR_LIT("orbital/alpha/coefficient"));
    const md_attribute_t* oa = md_attributes_find(&sys->attributes, beta ? STR_LIT("orbital/beta/occupation")  : STR_LIT("orbital/alpha/occupation"));
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
    const md_attribute_t* ea = md_attributes_find(&sys->attributes, beta ? STR_LIT("orbital/beta/energy") : STR_LIT("orbital/alpha/energy"));
    if (!ea || md_attribute_element_count(&ea->format) != num_mo) return false;
    if (md_attribute_extract_f64(dst, num_mo, ea, md_attribute_slice_all(), md_unit_hartree()) != num_mo) return false;
    for (size_t k = 0; k < num_mo; ++k) if (!isfinite(dst[k])) return false;
    return true;
}

// The orbitals of the density sa alpha + sb beta: rows {count, A}, weights {count} (the occupations times
// sa, sb), and with out_energy_weights the same with each occupation times |its orbital's energy|. Every
// array from 'alloc'. canonical: every occupation 0, 1 or 2 (within 1e-3), else nothing.
static size_t qm_orbital_products(double** out_rows, double** out_weights, double** out_energy_weights, size_t* out_num_ao,
                                  const md_system_t* sys, double sa, double sb, bool canonical, md_allocator_i* alloc) {
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
    if (!ok) return 0;
    // one set of orbitals for both spins (restricted, restricted open shell): each orbital once, with its two
    // weights summed
    const bool shared = have_a && have_b && nma == nmb && md_attribute_same_data(caa, cab);
    const size_t na = have_a ? naa : nab;
    const bool energy = out_energy_weights != NULL;

    // the occupations and energies (small), then the rows of exactly the orbitals of nonzero weight, read one by one
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    double* oa = have_a ? (double*)md_temp_alloc(temp, sizeof(double) * nma) : NULL;
    double* ob = have_b ? (double*)md_temp_alloc(temp, sizeof(double) * nmb) : NULL;
    double* ea = have_a && energy ? (double*)md_temp_alloc(temp, sizeof(double) * nma) : NULL;
    double* eb = have_b && energy ? (double*)md_temp_alloc(temp, sizeof(double) * nmb) : NULL;
    if ((have_a && (!oa || md_attribute_extract_f64(oa, nma, oaa, md_attribute_slice_all(), md_unit_none()) != nma)) ||
        (have_b && (!ob || md_attribute_extract_f64(ob, nmb, oab, md_attribute_slice_all(), md_unit_none()) != nmb)) ||
        (have_a && energy && (!ea || !qm_spin_orbital_energies(ea, nma, sys, false))) ||
        (have_b && energy && (!eb || !qm_spin_orbital_energies(eb, nmb, sys, !beta_is_alpha)))) {
        md_temp_end(temp);
        return 0;
    }
    if (canonical) {
        for (size_t k = 0; have_a && k < nma; ++k) if (fabs(oa[k] - round(oa[k])) > 1.0e-3 || oa[k] < -1.0e-3 || oa[k] > 2.001) { md_temp_end(temp); return 0; }
        for (size_t k = 0; have_b && k < nmb; ++k) if (fabs(ob[k] - round(ob[k])) > 1.0e-3 || ob[k] < -1.0e-3 || ob[k] > 2.001) { md_temp_end(temp); return 0; }
    }
    size_t count = 0;
    for (size_t k = 0; have_a && k < nma; ++k) count += (sa * oa[k] + (shared ? sb * ob[k] : 0.0)) != 0.0;
    for (size_t k = 0; have_b && !shared && k < nmb; ++k) count += (sb * ob[k]) != 0.0;
    double* rows = count ? (double*)md_alloc(alloc, sizeof(double) * count * na) : NULL;
    double* w    = count ? (double*)md_alloc(alloc, sizeof(double) * count) : NULL;
    double* we   = count && energy ? (double*)md_alloc(alloc, sizeof(double) * count) : NULL;
    ok = rows && w && (!energy || we);
    size_t n = 0;
    for (size_t k = 0; ok && have_a && k < nma; ++k) {
        const double wk = sa * oa[k] + (shared ? sb * ob[k] : 0.0);
        if (wk == 0.0) continue;
        ok = md_attribute_extract_f64(rows + n * na, na, caa, md_attribute_slice_1((uint32_t)k), md_unit_none()) == na;
        // a shared orbital's two spins may differ in energy (restricted open shell): each its own
        if (energy) we[n] = sa * oa[k] * fabs(ea[k]) + (shared ? sb * ob[k] * fabs(eb[k]) : 0.0);
        w[n++] = wk;
    }
    for (size_t k = 0; ok && have_b && !shared && k < nmb; ++k) {
        const double wk = sb * ob[k];
        if (wk == 0.0) continue;
        ok = md_attribute_extract_f64(rows + n * na, na, cab, md_attribute_slice_1((uint32_t)k), md_unit_none()) == na;
        if (energy) we[n] = wk * fabs(eb[k]);
        w[n++] = wk;
    }
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

size_t md_qm_density_orbitals(double** out_rows, double** out_weights, size_t* out_num_ao, const md_system_t* sys,
                              str_t density_path, md_allocator_i* alloc) {
    if (!out_rows || !out_weights || !out_num_ao || !sys || !alloc) return 0;
    *out_rows = NULL;
    *out_weights = NULL;
    *out_num_ao = 0;
    // the weights of alpha's and beta's occupations in the density asked for
    double sa = 0.0, sb = 0.0;
    if      (str_eq(density_path, STR_LIT("orbital/alpha/density")))      { sa = 1.0; }
    else if (str_eq(density_path, STR_LIT("orbital/beta/density")))       { sb = 1.0; }
    else if (str_eq(density_path, STR_LIT("orbital/total/density")))      { sa = 1.0; sb = 1.0; }
    else if (str_eq(density_path, STR_LIT("orbital/difference/density"))) { sa = 1.0; sb = -1.0; }
    else return 0;
    return qm_orbital_products(out_rows, out_weights, NULL, out_num_ao, sys, sa, sb, false, alloc);
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
    const md_attribute_t* ea = md_attributes_find(&sys->attributes, STR_LIT("orbital/alpha/energy"));
    if (!ea || md_attribute_element_count(&ea->format) != nm) return false;
    // the occupations of canonical orbitals: natural orbitals carry fractional ones, and no energies to weight by
    md_temp_scope_t temp = md_temp_begin();
    double* occ = (double*)md_temp_alloc(temp, sizeof(double) * nm);
    bool ok = occ && md_attribute_extract_f64(occ, nm, oa, md_attribute_slice_all(), md_unit_none()) == nm;
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
           md_attribute_extract_f64(*occ, *num_mo, oa, md_attribute_slice_all(), md_unit_none()) == *num_mo &&
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
    bool ok = rows && w;
    size_t n = 0;
    for (size_t k = 0; ok && k < ma; ++k) {
        const double wk = qm_attachment_weight(na[k], ea[k]) + (shared ? qm_attachment_weight(nb[k], eb[k]) : 0.0);
        if (wk == 0.0) continue;
        ok = md_attribute_extract_f64(rows + n * nao, nao, ca, md_attribute_slice_1((uint32_t)k), md_unit_none()) == nao;
        w[n++] = wk;
    }
    for (size_t k = 0; ok && !shared && k < mb; ++k) {
        const double wk = qm_attachment_weight(nb[k], eb[k]);
        if (wk == 0.0) continue;
        ok = md_attribute_extract_f64(rows + n * nao, nao, cb, md_attribute_slice_1((uint32_t)k), md_unit_none()) == nao;
        w[n++] = wk;
    }
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

bool md_qm_publish_orbital_density_combinations(md_system_t* sys) {
    ASSERT(sys);
    const md_attribute_t* alpha = md_attributes_find(&sys->attributes, STR_LIT("orbital/alpha/density"));
    const md_attribute_t* beta  = md_attributes_find(&sys->attributes, STR_LIT("orbital/beta/density"));
    if (!alpha || !beta || alpha->format.rank != 2 || alpha->format.shape[0] != alpha->format.shape[1]) {
        return false;
    }
    // read as values: publishing moves the table
    const uint32_t num_ao = alpha->format.shape[0];
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2, .shape = { num_ao, num_ao },
    };
    md_attribute_virtual_t total_virt = { .provider = qm_total_density_provider,      .user_data = sys };
    md_attribute_virtual_t diff_virt  = { .provider = qm_difference_density_provider, .user_data = sys };
    md_qm_publish_virtual(sys, STR_LIT("orbital/total/density"),      STR_LIT("Total Density"),           md_unit_none(), format, &total_virt);
    md_qm_publish_virtual(sys, STR_LIT("orbital/difference/density"), STR_LIT("Spin Difference Density"), md_unit_none(), format, &diff_virt);
    return true;
}
