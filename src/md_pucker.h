#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_vec_math.h>

struct md_system_t;

// ### RING PUCKERING ###
//
// How a ring of 5 to 7 atoms departs from a plane, from the positions of its atoms in the order around it.
//
// CREMER-POPLE (md_pucker_cremer_pople). The displacements z_j of the atoms from their mean plane (Cremer & Pople,
// J. Am. Chem. Soc. 97, 1354 (1975)), as amplitudes and phases: q_m and phi_m for m = 2 .. (N - 1) / 2, and q_(N/2)
// for an even N. The total amplitude Q is the length of z. A six membered ring is (Q, theta, phi) on a sphere:
// theta = atan2(q2, q3), 0 at one chair and 180 at the other, the boats and skews around the equator.
// What the numbers mean depends on the order of the atoms: the phases on where the ring starts, theta and the sign of
// q3 on its direction. Only Q, q2 and |q3| are the same whatever the order. md_pucker_ring_order puts the rings of a
// system in the order the conventions are written for.
//
// ALTONA-SUNDARALINGAM (md_pucker_pseudorotation). The phase P and amplitude tau_m of a furanose from its five
// endocyclic torsions (Altona & Sundaralingam, J. Am. Chem. Soc. 94, 8205 (1972)), with the atoms in the order
// O4', C1', C2', C3', C4': the convention of nucleic acids, C3'-endo near P = 18 and C2'-endo near P = 162 in a
// D-sugar. P is the Cremer-Pople phi2 plus 90 degrees, to within a degree.
//
// CONFORMATIONS (md_pucker_conformer). The canonical conformation nearest the coordinates, named as IUPAC does
// (Nomenclature of Carbohydrates, 2-Carb-7): the atoms out of a reference plane, those above as superscripts before
// the letter and those below as subscripts after it, above being the side from which the ring is numbered clockwise.
// Six membered rings: two chairs (C), six boats (B), six skews (S), twelve half chairs (H) and twelve envelopes (E),
// at their (theta, phi) on the sphere. Five membered rings: ten envelopes (E) and ten twists (T) along phi2. The table
// was derived from the geometry of ideal rings at each position, and agrees with 4C1 at theta 0, 3,OB at phi 0 and
// 1S3 at phi 210. A ring less puckered than MD_PUCKER_PLANAR_AMPLITUDE is planar. Seven membered rings are not named.

#define MD_PUCKER_PLANAR_AMPLITUDE 0.1f     // Å

typedef struct md_pucker_t {
    uint32_t size;          // Atoms of the ring, 5 to 7. 0 if there is no result
    float amplitude;        // Q (Å)
    float q2;               // Å
    float phi2;             // Degrees, [0, 360)
    float q3;               // Å. Six membered rings: q_(N/2), signed. Seven: the amplitude of the third harmonic
    float phi3;             // Seven membered rings: degrees, [0, 360)
    float theta;            // Six membered rings: degrees, [0, 180]
} md_pucker_t;

typedef enum md_pucker_kind_t {
    MD_PUCKER_KIND_NONE = 0,        // Not named (seven membered rings)
    MD_PUCKER_KIND_PLANAR,
    MD_PUCKER_KIND_CHAIR,
    MD_PUCKER_KIND_BOAT,
    MD_PUCKER_KIND_SKEW,
    MD_PUCKER_KIND_HALF_CHAIR,
    MD_PUCKER_KIND_ENVELOPE,
    MD_PUCKER_KIND_TWIST,
} md_pucker_kind_t;

// A canonical conformation: its kind, and the positions in the ring (0 is the first atom of the order) of the atoms
// above and below the reference plane, -1 for none
typedef struct md_pucker_conformer_t {
    md_pucker_kind_t kind;
    int8_t above[2];
    int8_t below[2];
} md_pucker_conformer_t;

#ifdef __cplusplus
extern "C" {
#endif

// The Cremer-Pople coordinates of the ring of count (5 to 7) positions, in order around it, made consistent across
// the periodic cell by the caller. Returns false for another count or a degenerate ring.
bool md_pucker_cremer_pople(md_pucker_t* out, const vec3_t* xyz, size_t count);

// The Altona-Sundaralingam phase P and amplitude tau_m of a furanose, in degrees, from the positions of
// O4', C1', C2', C3', C4' in that order.
bool md_pucker_pseudorotation(float* out_phase, float* out_amplitude, const vec3_t xyz[5]);

// The canonical conformation nearest p (see CONFORMATIONS)
md_pucker_conformer_t md_pucker_conformer(const md_pucker_t* p);

// "chair", "boat", ...
const char* md_pucker_kind_name(md_pucker_kind_t kind);

// The IUPAC name of a conformation, its atoms written with locants[position] ("4C1" with the locants "O", "1", ...,
// "5" of a pyranose, "3,OB" for a boat). Without locants, or for a planar ring, the kind alone. Returns the length
// written (snprintf).
size_t md_pucker_conformer_name(char* buf, size_t cap, md_pucker_conformer_t c, const char* const* locants);

// The atoms of ring ring_idx of sys (md_system_t::ring) in the order the conventions are written for, written to out
// (room for MD_RING_MAX_SIZE). Returns the number of atoms.
//   - A sugar ring (MD_RING_FLAG_SUGAR) as the system has it: the ring O, the anomeric C, and on.
//   - Another ring with an atom other than C: from it (the heaviest, then the one with the most heavy atoms bonded
//     outside of the ring), towards the neighbor with the most heavy atoms outside of the ring. Proline: N, CA, CB,
//     CG, CD.
//   - A ring of C: from the C with the most heavy atoms outside of the ring, in the same way.
// out_anchored (optional) tells whether the order is fixed by the chemistry. It is not when a choice came down to the
// index of an atom (cyclohexane, benzene): the phases, theta and the names of such a ring depend on that choice, and
// only the amplitudes and the kind of conformation mean something.
size_t md_pucker_ring_order(int32_t* out, bool* out_anchored, const struct md_system_t* sys, size_t ring_idx);

#ifdef __cplusplus
}
#endif
