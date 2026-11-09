#include <core/md_linalg.h>
#include <core/md_allocator.h>
#include <core/md_common.h>

#include <float.h>
#include <math.h>

// Householder reduction to tridiagonal form and the implicit QL algorithm, after JAMA's tred2 and
// tql2 (public domain; from EISPACK). Column major so that every inner loop runs down a column, which
// is contiguous; a symmetric input reads the same either way. 'e' is n doubles of scratch. On return
// d holds the eigenvalues ascending and column j of 'a' the eigenvector of d[j].
//
// This used to exist twice, in md_vlx.c (natural transition orbitals) and md_gto_iso.c (density
// factorization); this is the md_vlx one, which had its inner loops laid out for vectorization.
static bool sym_eigen_tred2_tql2(double* a, double* d, double* e, size_t n) {
#define A_(r, c) a[(size_t)(c) * n + (size_t)(r)]
    if (n == 0) return false;

    // tred2
    for (size_t j = 0; j < n; ++j) d[j] = A_(n - 1, j);
    for (size_t i = n - 1; i > 0; --i) {
        double scale = 0.0, h = 0.0;
        for (size_t k = 0; k < i; ++k) scale += fabs(d[k]);
        if (scale == 0.0) {
            e[i] = d[i - 1];
            for (size_t j = 0; j < i; ++j) {
                d[j] = A_(i - 1, j);
                A_(i, j) = 0.0;
                A_(j, i) = 0.0;
            }
        } else {
            for (size_t k = 0; k < i; ++k) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            double f = d[i - 1];
            double g = sqrt(h);
            if (f > 0) g = -g;
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for (size_t j = 0; j < i; ++j) e[j] = 0.0;
            for (size_t j = 0; j < i; ++j) {
                f = d[j];
                A_(j, i) = f;
                const double* col = &A_(0, j);
                // four lanes for the sum, which a compiler may not reorder on its own
                double g4[4] = { 0.0, 0.0, 0.0, 0.0 };
                size_t k = j + 1;
                for (; k + 4 <= i; k += 4) {
                    for (int l = 0; l < 4; ++l) {
                        g4[l]    += col[k + l] * d[k + l];
                        e[k + l] += col[k + l] * f;
                    }
                }
                g = e[j] + A_(j, j) * f + ((g4[0] + g4[1]) + (g4[2] + g4[3]));
                for (; k < i; ++k) {
                    g    += col[k] * d[k];
                    e[k] += col[k] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (size_t j = 0; j < i; ++j) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            const double hh = f / (h + h);
            for (size_t j = 0; j < i; ++j) e[j] -= hh * d[j];
            for (size_t j = 0; j < i; ++j) {
                f = d[j];
                g = e[j];
                double* col = &A_(0, j);
                for (size_t k = j; k < i; ++k) col[k] -= (f * e[k] + g * d[k]);
                d[j] = A_(i - 1, j);
                A_(i, j) = 0.0;
            }
        }
        d[i] = h;
    }
    for (size_t i = 0; i + 1 < n; ++i) {
        A_(n - 1, i) = A_(i, i);
        A_(i, i) = 1.0;
        const double h = d[i + 1];
        double* ci1 = &A_(0, i + 1);
        if (h != 0.0) {
            for (size_t k = 0; k <= i; ++k) d[k] = ci1[k] / h;
            for (size_t j = 0; j <= i; ++j) {
                double* cj = &A_(0, j);
                double g4[4] = { 0.0, 0.0, 0.0, 0.0 };
                size_t k = 0;
                for (; k + 4 <= i + 1; k += 4) {
                    for (int l = 0; l < 4; ++l) g4[l] += ci1[k + l] * cj[k + l];
                }
                double g = (g4[0] + g4[1]) + (g4[2] + g4[3]);
                for (; k <= i; ++k) g += ci1[k] * cj[k];
                for (size_t k = 0; k <= i; ++k) cj[k] -= g * d[k];
            }
        }
        for (size_t k = 0; k <= i; ++k) ci1[k] = 0.0;
    }
    for (size_t j = 0; j < n; ++j) {
        d[j] = A_(n - 1, j);
        A_(n - 1, j) = 0.0;
    }
    A_(n - 1, n - 1) = 1.0;
    e[0] = 0.0;

    // tql2
    for (size_t i = 1; i < n; ++i) e[i - 1] = e[i];
    e[n - 1] = 0.0;
    double f = 0.0, tst1 = 0.0;
    const double eps = DBL_EPSILON;
    for (size_t l = 0; l < n; ++l) {
        tst1 = MD_MAX(tst1, fabs(d[l]) + fabs(e[l]));
        size_t m = l;
        while (m < n && fabs(e[m]) > eps * tst1) ++m;
        if (m == n) m = n - 1; // e[n-1] is zero; not reached in exact arithmetic
        if (m > l) {
            int iter = 0;
            do {
                if (++iter > 64) {
                    return false;
                }
                double g = d[l];
                double p = (d[l + 1] - g) / (2.0 * e[l]);
                double r = hypot(p, 1.0);
                if (p < 0) r = -r;
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                const double dl1 = d[l + 1];
                double h = g - d[l];
                for (size_t i = l + 2; i < n; ++i) d[i] -= h;
                f += h;

                p = d[m];
                double c = 1.0, c2 = c, c3 = c;
                const double el1 = e[l + 1];
                double s = 0.0, s2 = 0.0;
                for (size_t ii = m; ii-- > l;) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[ii];
                    h = c * p;
                    r = hypot(p, e[ii]);
                    e[ii + 1] = s * r;
                    s = e[ii] / r;
                    c = p / r;
                    p = c * d[ii] - s * g;
                    d[ii + 1] = h + s * (c * g + s * d[ii]);
                    double* x = &A_(0, ii);
                    double* y = &A_(0, ii + 1);
                    for (size_t k = 0; k < n; ++k) {
                        const double t = y[k];
                        y[k] = s * x[k] + c * t;
                        x[k] = c * x[k] - s * t;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
            } while (fabs(e[l]) > eps * tst1);
        }
        d[l] += f;
        e[l] = 0.0;
    }

    // ascending, the vectors with them
    for (size_t i = 0; i + 1 < n; ++i) {
        size_t k = i;
        double p = d[i];
        for (size_t j = i + 1; j < n; ++j) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            double* x = &A_(0, i);
            double* y = &A_(0, k);
            for (size_t j = 0; j < n; ++j) {
                const double t = x[j];
                x[j] = y[j];
                y[j] = t;
            }
        }
    }
#undef A_
    return true;
}

bool md_linalg_sym_eigen(double* a, double* eigval, size_t n) {
    if (!a || !eigval || n == 0) return false;
    md_temp_scope_t temp = md_temp_begin();
    double* e = md_temp_alloc_array(temp, double, n);
    const bool ok = e && sym_eigen_tred2_tql2(a, eigval, e, n);
    md_temp_end(temp);
    return ok;
}
