#include "md_unitcell.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

int md_unitcell_print(char* out_buf, size_t buf_cap, const md_unitcell_t* cell) {
    return snprintf(out_buf, buf_cap, "x: %f, y: %f, z: %f \nxy: %f, xz: %f, yz: %f \nflags: %i", cell->x, cell->y, cell->z, cell->xy, cell->xz, cell->yz, cell->flags);
}

md_unitcell_t md_unitcell_from_basis_parameters(double x, double y, double z, double xy, double xz, double yz) {
    md_unitcell_flags_t flags = MD_UNITCELL_NONE;

    if (xy == 0.0 && xz == 0.0 && yz == 0.0) {
        if (!(x == 0.0 && y == 0.0 && z == 0.0) && !(x == 1.0 && y == 1.0 && z == 1.0)) {
            flags |= MD_UNITCELL_ORTHO;
        }
    } else {
        flags |= MD_UNITCELL_TRICLINIC;
    }

    if (flags) {
        if (x != 0.0) flags |= MD_UNITCELL_PBC_X;
        if (y != 0.0) flags |= MD_UNITCELL_PBC_Y;
        if (z != 0.0) flags |= MD_UNITCELL_PBC_Z;
    }

    md_unitcell_t cell = {.x = x, .xy = xy, .xz = xz, .y = y, .yz = yz, .z = z, .flags = flags};
    return cell;
}

// From here: https://www.arianarab.com/post/crystal-structure-software
// Assumes that all input angles (alpha, beta, gamma) are given in degrees
md_unitcell_t md_unitcell_from_extent_and_angles(double a, double b, double c, double alpha, double beta, double gamma) {
    if (a == 0 && b == 0 && c == 0 && alpha == 0 && beta == 0 && gamma == 0) {
        return md_unitcell_none();
    }

    if (alpha == 90.0 && beta == 90.0 && gamma == 90.0) {
        return md_unitcell_from_basis_parameters(a, b, c, 0, 0, 0);
    }

    alpha = DEG_TO_RAD(alpha);
    beta  = DEG_TO_RAD(beta);
    gamma = DEG_TO_RAD(gamma);

    // https://docs.lammps.org/Howto_triclinic.html
    double x  = a;
    double xy = b * cos(gamma);
    double xz = c * cos(beta);
    double y  = b * sin(gamma);
    double yz = (b * c * cos(alpha) - xy * xz) / y;
    double z  = sqrt(c * c - xz * xz - yz * yz);

    return md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
}

void md_unitcell_extract_extent_angles(double* out_a, double* out_b, double* out_c, double* out_alpha, double* out_beta, double* out_gamma, const md_unitcell_t* cell) {
    if (!cell) {
        if (out_a) *out_a = 0;
        if (out_b) *out_b = 0;
        if (out_c) *out_c = 0;
        if (out_alpha) *out_alpha = 0;
        if (out_beta) *out_beta = 0;
        if (out_gamma) *out_gamma = 0;
        return;
    }

    double a = cell->x;
    double b = sqrt(cell->xy * cell->xy + cell->y * cell->y);
    double c = sqrt(cell->xz * cell->xz + cell->yz * cell->yz + cell->z * cell->z);

    double alpha = RAD_TO_DEG(acos((cell->yz * cell->xy - cell->y * cell->xz) / (b * c)));
    double beta  = RAD_TO_DEG(acos((cell->xz) / c));
    double gamma = RAD_TO_DEG(acos((cell->xy) / b));

    if (out_a) *out_a = a;
    if (out_b) *out_b = b;
    if (out_c) *out_c = c;
    if (out_alpha) *out_alpha = alpha;
    if (out_beta) *out_beta = beta;
    if (out_gamma) *out_gamma = gamma;
}

void md_unitcell_extract_basis_parameters(double* out_x, double* out_y, double* out_z, double* out_xy, double* out_xz, double* out_yz, const md_unitcell_t* cell) {
    if (!cell) {
        if (out_x) *out_x = 0;
        if (out_y) *out_y = 0;
        if (out_z) *out_z = 0;
        if (out_xy) *out_xy = 0;
        if (out_xz) *out_xz = 0;
        if (out_yz) *out_yz = 0;
        return;
    }
    if (out_x) *out_x = cell->x;
    if (out_y) *out_y = cell->y;
    if (out_z) *out_z = cell->z;
    if (out_xy) *out_xy = cell->xy;
    if (out_xz) *out_xz = cell->xz;
    if (out_yz) *out_yz = cell->yz;
}

void md_unitcell_I_extract_double(double out_I[3][3], const md_unitcell_t* cell) {
    if (cell) {
        if (!cell->flags) {
            MEMSET(out_I, 0, sizeof(double) * 9);
            return;
        }

        const double i11 = cell->x > 0.0 ? 1.0 / cell->x : 0.0;
        const double i22 = cell->y > 0.0 ? 1.0 / cell->y : 0.0;
        const double i33 = cell->z > 0.0 ? 1.0 / cell->z : 0.0;
        const double i12 = (cell->x * cell->y) > 0.0 ? -cell->xy / (cell->x * cell->y) : 0.0;
        const double i13 = (cell->x * cell->y * cell->z) > 0.0 ? (cell->xy * cell->yz - cell->xz * cell->y) / (cell->x * cell->y * cell->z) : 0.0;
        const double i23 = (cell->y * cell->z) > 0.0 ? -cell->yz / (cell->y * cell->z) : 0.0;

        out_I[0][0] = i11; out_I[0][1] = 0.0; out_I[0][2] = 0.0;
        out_I[1][0] = i12; out_I[1][1] = i22; out_I[1][2] = 0.0;
        out_I[2][0] = i13; out_I[2][1] = i23; out_I[2][2] = i33;
    }
}

void md_unitcell_I_extract_float(float out_I[3][3], const md_unitcell_t* cell) {
    if (cell) {
        double I[3][3];
        md_unitcell_I_extract_double(I, cell);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                out_I[i][j] = (float)I[i][j];
            }
        }
    }
}

void md_unitcell_G_extract_double(double out_G[3][3], const md_unitcell_t* cell) {
    if (cell) {
        const double x  = cell->x;
        const double xy = cell->xy;
        const double xz = cell->xz;
        const double y  = cell->y;
        const double yz = cell->yz;
        const double z  = cell->z;

        out_G[0][0] = x * x;
        out_G[0][1] = x * xy;
        out_G[0][2] = x * xz;

        out_G[1][0] = out_G[0][1];
        out_G[1][1] = xy * xy + y * y;
        out_G[1][2] = xy * xz + y * yz;

        out_G[2][0] = out_G[0][2];
        out_G[2][1] = out_G[1][2];
        out_G[2][2] = xz * xz + yz * yz + z * z;
    }
}

void md_unitcell_G_extract_float(float out_G[3][3], const md_unitcell_t* cell) {
    if (cell) {
        double G[3][3];
        md_unitcell_G_extract_double(G, cell);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                out_G[i][j] = (float)G[i][j];
            }
        }
    }
}
