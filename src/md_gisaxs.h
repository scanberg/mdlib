#pragma once

#include <stddef.h>
#include <stdbool.h>

// GISAXS / SAXS intensity from a particle (bead / atom) representation of a periodic slab.
//
// The scattering density is approximated as a sum of isotropic Gaussians
//     rho(r) = sum_j w_j * G_sigma_j(r - r_j)
// where w_j is the scattering weight of particle j (electrons for X-rays) and sigma_j its Gaussian width.
//
// The full 3D density is never constructed. The system is processed as a stack of XY slices (periodic in XY):
//   1. Particles are sorted along z and assigned as point masses to slices (spacing dz) and to a 2D periodic grid
//      within each slice, using cubic B-spline assignment in all three directions (as in particle mesh Ewald).
//   2. Each slice is Fourier transformed (2D real to complex, PFFFT). The in-plane B-spline window is divided out
//      and the exact in-plane Gaussian form factor exp(-q_par^2 sigma^2 / 2) is applied. Particles are grouped into
//      a few classes of equal sigma, each class is transformed separately.
//   3. The slice spectra A_m(q_par) (m = class * num_slices + slice) are reduced into a ring (|q_par|) averaged
//      cross spectral matrix
//          S_mn(ring) = < Re( A_m(q_par) conj(A_n(q_par)) ) >_{q_par in ring}
//      which corresponds to a full in-plane (azimuthal) rotational average. S is real symmetric.
//   4. Any scattering model where the amplitude is a linear combination of terms exp(-i Q z), i.e.
//          F(q_par) = sum_m c_m A_m(q_par),  c_m = sum_t coef_t * exp(-i Q_t z_k) * exp(-Q_t^2 sigma_c^2 / 2) / W(Q_t dz)
//      gives the rotationally averaged intensity as I = c^H S c. This includes the Born approximation
//      (a single term Q = q_z) and the DWBA for particles above a substrate (4 terms with Fresnel coefficients).
//      The z-direction B-spline window W and the Gaussian z form factor are applied analytically here.
//      Changing beam or substrate parameters only requires step 4, which is cheap.
//
// All lengths are in Ångström and wave vectors in 1/Ångström.
//
// Typical usage (single threaded):
//     md_gisaxs_t* ctx = md_gisaxs_create(&input, &params, alloc);
//     md_gisaxs_compute(ctx);                    // steps 1-3
//     md_gisaxs_evaluate(ctx, &model, qz, num_qz, out);   // step 4, can be called repeatedly
//     md_gisaxs_destroy(ctx);
//
// Multi threaded usage: call md_gisaxs_compute_slices over disjoint slice ranges (each worker with its own scratch),
// then md_gisaxs_compute_rings over disjoint ring ranges, and md_gisaxs_evaluate_range over disjoint q_z ranges.

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;

typedef struct md_gisaxs_t md_gisaxs_t;

typedef struct md_gisaxs_input_t {
    size_t count;
    const float* x;         // Å
    const float* y;         // Å
    const float* z;         // Å
    const float* weight;    // Scattering weight per particle (e.g. number of electrons). NULL -> 1.0
    const float* sigma;     // Gaussian width per particle (Å). NULL -> sigma_uniform
    float sigma_uniform;

    // In-plane periodic box extent (Å), the box is assumed orthorhombic with the XY plane perpendicular to z.
    double box_x;
    double box_y;
} md_gisaxs_input_t;

typedef struct md_gisaxs_params_t {
    double q_par_max;       // Maximum in-plane |q| (1/Å)
    double q_z_max;         // Maximum |q_z| used in any scattering term (1/Å), determines the slice spacing
    double oversampling;    // In-plane grid oversampling relative to q_par_max (default 2.0 if <= 0)
    double dz;              // Slice spacing (Å), 0 -> automatic: pi / (oversampling * q_z_max)
    double z_min;           // Particle z range to include (Å), if z_max <= z_min all particles are included
    double z_max;
    size_t max_slices;      // Upper limit on the number of slices, 0 -> 1024. dz is increased if required.
} md_gisaxs_params_t;

typedef struct md_gisaxs_info_t {
    int    nx, ny;          // In-plane grid dimensions
    double dx, dy;          // In-plane grid spacing (Å)
    size_t num_slices;
    double dz;              // Slice spacing (Å)
    double z0;              // z of first slice (Å)
    size_t num_rings;
    double dq_ring;         // Ring width (1/Å)
    size_t num_points;      // Number of (half plane) q_par grid points within q_par_max
    size_t num_classes;     // Number of distinct Gaussian widths (classes)
    double class_sigma[4];  // Gaussian width per class (Å)
    size_t num_particles;
    double area;            // box_x * box_y (Å^2)
    size_t spectra_bytes;   // Memory for the slice spectra (released after rings are computed)
    size_t matrix_bytes;    // Memory for the ring matrices
} md_gisaxs_info_t;

// Scattering model evaluated on top of the ring matrices.
// The ambient medium is the medium where the particles reside and where the beam enters from.
// Scattering length densities (SLD) are given in 1/Å^2. For X-rays: SLD = r_e * rho_e, with r_e = 2.8179403e-5 Å
// and rho_e the electron density (e/Å^3). Absorption is given as a positive imaginary part (SLD = re - i*im),
// related to the refractive index n = 1 - delta + i*beta through im = 2*pi*beta / lambda^2.
//
// DWBA reference medium
// - Simple (graded = false): ambient half space above a substrate at z_substrate. The particles scatter against
//   the ambient medium.
// - Graded (graded = true): the laterally averaged particle density (the slice profile) is part of the reference
//   medium, i.e. the stack is ambient / graded film (one layer per slice) / substrate. The wave fields inside the
//   film are obtained with the Parratt recursion and the particles scatter through the in-plane fluctuations
//   (q_par != 0) of the density around its lateral average. This gives refraction in the film and the film Yoneda
//   peak. The film SLD of slice k is sld_ambient + profile_sld_scale * profile[k] (absorption: profile_abs_scale *
//   profile[k]), where profile is md_gisaxs_slice_profile (weight / Å^3).
typedef struct md_gisaxs_model_t {
    double wavelength;          // Å
    double alpha_i;             // Incidence angle (radians), measured in the ambient medium

    bool   dwba;                // true: DWBA with substrate, false: Born approximation (no substrate)
    bool   graded;              // DWBA only: include the laterally averaged film in the reference medium
    double z_substrate;         // z coordinate of the substrate interface (Å)
    double sld_ambient;         // 1/Å^2 (real)
    double sld_substrate;       // 1/Å^2 (real part)
    double sld_substrate_abs;   // 1/Å^2 (absorption, imaginary part magnitude)
    double substrate_roughness; // RMS roughness of the substrate interface (Å), Nevot-Croce factor

    double profile_sld_scale;   // SLD per unit profile density (Å^-2 per weight/Å^3), e.g. r_e * contrast
    double profile_abs_scale;   // Absorption SLD per unit profile density

    double intensity_scale;     // Multiplied with the result (e.g. r_e^2 and contrast factors). 0 -> 1.0
} md_gisaxs_model_t;

md_gisaxs_t* md_gisaxs_create(const md_gisaxs_input_t* input, const md_gisaxs_params_t* params, struct md_allocator_i* alloc);
void md_gisaxs_destroy(md_gisaxs_t* ctx);

void md_gisaxs_get_info(const md_gisaxs_t* ctx, md_gisaxs_info_t* info);

// --- Steps 1-3 ---
// Single threaded convenience: computes all slices and rings, then releases the slice spectra.
bool md_gisaxs_compute(md_gisaxs_t* ctx);

// Number of bytes of scratch required per worker for md_gisaxs_compute_slices (64 byte alignment required)
size_t md_gisaxs_slice_scratch_bytes(const md_gisaxs_t* ctx);
// Compute slices in [beg, end). Thread safe for disjoint ranges given distinct scratch buffers.
void md_gisaxs_compute_slices(md_gisaxs_t* ctx, size_t beg, size_t end, void* scratch);

// Compute ring matrices in [beg, end). Requires all slices to be computed. Thread safe for disjoint ranges.
void md_gisaxs_compute_rings(md_gisaxs_t* ctx, size_t beg, size_t end);

// Releases the slice spectra (the largest allocation), after all rings have been computed.
void md_gisaxs_release_spectra(md_gisaxs_t* ctx);

// --- Results ---
size_t        md_gisaxs_num_rings(const md_gisaxs_t* ctx);
const double* md_gisaxs_ring_q(const md_gisaxs_t* ctx);           // Mean |q_par| per ring (1/Å)
const unsigned* md_gisaxs_ring_count(const md_gisaxs_t* ctx);     // Number of q_par grid points (full plane) per ring
size_t        md_gisaxs_num_slices(const md_gisaxs_t* ctx);
const double* md_gisaxs_slice_z(const md_gisaxs_t* ctx);          // z per slice (Å)
double        md_gisaxs_particle_z_min(const md_gisaxs_t* ctx);   // Extent of included particles in z (Å)
double        md_gisaxs_particle_z_max(const md_gisaxs_t* ctx);
// Laterally averaged scattering density per slice (weight/Å^3), i.e. the q_par = 0 component, B-spline smoothed
// and without the Gaussian profile applied. Useful for reflectometry.
const double* md_gisaxs_slice_profile(const md_gisaxs_t* ctx);

// --- Step 4 ---
// out: num_qz * num_rings values, row major (row = q_z index, column = ring index)
// Intensity is <|F|^2>_ring / area * intensity_scale, i.e. (weight units)^2 / Å^2
// q_z is the vacuum/ambient scattering vector component k (sin(alpha_i) + sin(alpha_f)).
// For DWBA, rows below the horizon (alpha_f < 0) are set to zero.
void md_gisaxs_evaluate(const md_gisaxs_t* ctx, const md_gisaxs_model_t* model, const double* qz, size_t num_qz, float* out);
void md_gisaxs_evaluate_range(const md_gisaxs_t* ctx, const md_gisaxs_model_t* model, const double* qz, size_t qz_beg, size_t qz_end, float* out);

// Specular reflectivity |r|^2 of the model's reference medium (DWBA stack) for specular q_z = 2 k sin(alpha).
// out: num_qz values
void md_gisaxs_reflectivity(const md_gisaxs_t* ctx, const md_gisaxs_model_t* model, const double* qz, size_t num_qz, double* out);

// Brute force reference: rotationally averaged Born intensity evaluated by explicit summation over particles,
// using the same ring/grid definition as ctx (<|F|^2>_ring / area). Only intended for testing small systems.
// out: num_qz * num_rings
void md_gisaxs_reference_born(const md_gisaxs_t* ctx, const md_gisaxs_input_t* input, const double* qz, size_t num_qz, double* out);

// Classic X-ray constants
#define MD_GISAXS_R_E 2.8179403262e-5   // Classical electron radius (Å)

#ifdef __cplusplus
}
#endif
