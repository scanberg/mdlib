#pragma once

#include <core/md_common.h>
#include <stddef.h>
#include <stdbool.h>

// Cell information
typedef enum {
    MD_UNITCELL_NONE          = 0,
    MD_UNITCELL_ORTHO         = 1,
    MD_UNITCELL_TRICLINIC     = 2,
    MD_UNITCELL_PBC_X         = 4,
    MD_UNITCELL_PBC_Y         = 8,
    MD_UNITCELL_PBC_Z         = 16,
    MD_UNITCELL_PBC_ALL       = 4 | 8 | 16,
} md_unitcell_flags_t;

ENUM_FLAGS(md_unitcell_flags_t)

typedef struct md_unitcell_t {
    double x, xy, xz;
    double y, yz;
    double z;
    md_unitcell_flags_t flags;
} md_unitcell_t;

#ifdef __cplusplus
extern "C" {
#endif

int md_unitcell_print(char* out_buf, size_t buf_cap, const md_unitcell_t* cell);

// --- Construction ---

static inline md_unitcell_t md_unitcell_none(void) {
    md_unitcell_t cell = {0};
    return cell;
}

md_unitcell_t md_unitcell_from_basis_parameters(double x, double y, double z, double xy, double xz, double yz);

// Assumes that all input angles (alpha, beta, gamma) are given in degrees
md_unitcell_t md_unitcell_from_extent_and_angles(double a, double b, double c, double alpha, double beta, double gamma);

static inline md_unitcell_t md_unitcell_from_extent(double x, double y, double z) {
    return md_unitcell_from_basis_parameters(x, y, z, 0, 0, 0);
}

// Construct unitcell from float matrix [3][3] (column major)
static inline md_unitcell_t md_unitcell_from_matrix_float(const float A[3][3]) {
    return md_unitcell_from_basis_parameters(A[0][0], A[1][1], A[2][2], A[1][0], A[2][0], A[2][1]);
}

// Construct unitcell from double matrix [3][3] (column major)
static inline md_unitcell_t md_unitcell_from_matrix_double(const double A[3][3]) {
    return md_unitcell_from_basis_parameters(A[0][0], A[1][1], A[2][2], A[1][0], A[2][0], A[2][1]);
}

// --- Parameter extraction ---

// Angles are returned in degrees. Any output pointer may be NULL.
void md_unitcell_extract_extent_angles(double* out_a, double* out_b, double* out_c, double* out_alpha, double* out_beta, double* out_gamma, const md_unitcell_t* cell);

// Any output pointer may be NULL.
void md_unitcell_extract_basis_parameters(double* out_x, double* out_y, double* out_z, double* out_xy, double* out_xz, double* out_yz, const md_unitcell_t* cell);

// --- Getters and helper functionality ---

static inline md_unitcell_flags_t md_unitcell_flags(const md_unitcell_t* cell) {
    if (cell) return cell->flags;
    return MD_UNITCELL_NONE;
}

static inline bool md_unitcell_is_triclinic(const md_unitcell_t* cell)    { return cell->flags & MD_UNITCELL_TRICLINIC; }
static inline bool md_unitcell_is_orthorhombic(const md_unitcell_t* cell) { return cell->flags & MD_UNITCELL_ORTHO; }

// extracts unit_cell basis matrix A in the convention A[col][row]
// (basis vectors a,b,c are columns: A[0]=a, A[1]=b, A[2]=c)
static inline void md_unitcell_A_extract_double(double out_A[3][3], const md_unitcell_t* cell) {
    if (cell) {
        out_A[0][0] = cell->x;
        out_A[0][1] = 0;
        out_A[0][2] = 0;
        out_A[1][0] = cell->xy;
        out_A[1][1] = cell->y;
        out_A[1][2] = 0;
        out_A[2][0] = cell->xz;
        out_A[2][1] = cell->yz;
        out_A[2][2] = cell->z;
    }
}

static inline void md_unitcell_A_extract_float(float out_A[3][3], const md_unitcell_t* cell) {
    if (cell) {
        out_A[0][0] = (float)cell->x;
        out_A[0][1] = 0;
        out_A[0][2] = 0;
        out_A[1][0] = (float)cell->xy;
        out_A[1][1] = (float)cell->y;
        out_A[1][2] = 0;
        out_A[2][0] = (float)cell->xz;
        out_A[2][1] = (float)cell->yz;
        out_A[2][2] = (float)cell->z;
    }
}

// extracts the unit_cell inverse basis matrix (A^-1)
void md_unitcell_I_extract_double(double out_I[3][3], const md_unitcell_t* cell);
void md_unitcell_I_extract_float (float  out_I[3][3], const md_unitcell_t* cell);

// extracts the unitcell metric tensor G=(A^T)A
void md_unitcell_G_extract_double(double out_G[3][3], const md_unitcell_t* cell);
void md_unitcell_G_extract_float (float  out_G[3][3], const md_unitcell_t* cell);

// Create a mask which represents the periodic dimensions from a unit cell.
// I.e. [0,-1,-1] -> periodic in y and z, but not x
static inline void md_unitcell_pbc_mask_extract(int out_mask[3], const md_unitcell_t* cell) {
    ASSERT(out_mask);
    ASSERT(cell);

    out_mask[0] = (cell->flags & MD_UNITCELL_PBC_X) ? -1 : 0;
    out_mask[1] = (cell->flags & MD_UNITCELL_PBC_Y) ? -1 : 0;
    out_mask[2] = (cell->flags & MD_UNITCELL_PBC_Z) ? -1 : 0;
}

static inline void md_unitcell_diag_extract_double(double out_diag[3], const md_unitcell_t* cell) {
    if (cell) {
        out_diag[0] = cell->x;
        out_diag[1] = cell->y;
        out_diag[2] = cell->z;
    }
}

static inline void md_unitcell_diag_extract_float(float out_diag[3], const md_unitcell_t* cell) {
    if (cell) {
        out_diag[0] = (float)cell->x;
        out_diag[1] = (float)cell->y;
        out_diag[2] = (float)cell->z;
    }
}

#ifdef __cplusplus
} // extern "C"
#endif

#ifdef __cplusplus

// Fuck off C++ and the standard library with your nonsense
// Making me do all this nonsense for simple bullshit.

template<typename A, typename B>
struct md_is_same { static constexpr bool value = false; };

template<typename A>
struct md_is_same<A, A> { static constexpr bool value = true; };

template<typename A, typename B>
constexpr bool md_is_same_v = md_is_same<A, B>::value;

template<typename T>
struct md_false { static constexpr bool value = false; };

// C++ template dispatch
template<typename T>
inline void md_unitcell_A_extract(T (&out_A)[3][3], const md_unitcell_t* cell) {
    if constexpr (md_is_same_v<T, double>) {
        md_unitcell_A_extract_double(out_A, cell);
    } else if constexpr (md_is_same_v<T, float>) {
        md_unitcell_A_extract_float(out_A, cell);
    } else {
        static_assert(md_false<T>::value, "Unsupported type for md_unitcell_A_extract");
    }
}

template<typename T>
inline void md_unitcell_I_extract(T (&out_I)[3][3], const md_unitcell_t* cell) {
    if constexpr (md_is_same_v<T, double>) {
        md_unitcell_I_extract_double(out_I, cell);
    } else if constexpr (md_is_same_v<T, float>) {
        md_unitcell_I_extract_float(out_I, cell);
    } else {
        static_assert(md_false<T>::value, "Unsupported type for md_unitcell_I_extract");
    }
}

template<typename T>
inline void md_unitcell_G_extract(T (&out_G)[3][3], const md_unitcell_t* cell) {
    if constexpr (md_is_same_v<T, double>) {
        md_unitcell_G_extract_double(out_G, cell);
    } else if constexpr (md_is_same_v<T, float>) {
        md_unitcell_G_extract_float(out_G, cell);
    } else {
        static_assert(md_false<T>::value , "Unsupported type for md_unitcell_G_extract");
    }
}

#else
// C11 generics

#define md_unitcell_A_extract(out_A, cell) _Generic((out_A), \
    double (*)[3]: md_unitcell_A_extract_double, \
    double (*)   : md_unitcell_A_extract_double, \
    float  (*)[3]: md_unitcell_A_extract_float,  \
    float  (*)   : md_unitcell_A_extract_float   \
)(out_A, cell)

#define md_unitcell_I_extract(out_I, cell) _Generic((out_I), \
    double (*)[3]: md_unitcell_I_extract_double, \
    float  (*)[3]: md_unitcell_I_extract_float \
)(out_I, cell)

#define md_unitcell_G_extract(out_G, cell) _Generic((out_G), \
    double (*)[3]: md_unitcell_G_extract_double, \
    float  (*)[3]: md_unitcell_G_extract_float \
)(out_G, cell)

#endif
