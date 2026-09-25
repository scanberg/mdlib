#pragma once

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_tpr_data_t;
struct md_allocator_i;

// Non-bonded pair potentials as an MD engine evaluates them: Lennard-Jones and Coulomb between two particles,
// cut off, and modified so that they go to zero at the cut-off. Evaluated this way the energy of a pair is
// exactly its share of the engine's short range energy, which makes energies of contacts comparable to (and
// summable into) what the simulation itself saw.
//
// Units are those of the force fields: nm, kJ/mol and e. c6 in kJ/mol nm^6, c12 in kJ/mol nm^12.
//
// Only pairs that interact are evaluated here. Excluded pairs (bonded neighbours within a molecule) are the
// caller's to skip. Engines add corrections for excluded pairs within the cut-off for reaction field and Ewald,
// and a self term per charge; those belong to the molecule itself, not to any contact, and are not modelled.
// Ewald (PME) is evaluated as its real space part: the reciprocal part is not a sum over pairs.
//
// The forms follow GROMACS (interaction_const.cpp and the nbnxm kernels) and are checked against its energies.

typedef enum md_nb_modifier_t {
    MD_NB_MODIFIER_NONE = 0,        // Plain cut-off: the potential jumps to zero at the cut-off
    MD_NB_MODIFIER_POT_SHIFT,       // V(r) - V(rc)
    MD_NB_MODIFIER_POT_SWITCH,      // V(r) S(r), S going smoothly from 1 at the switch distance to 0 at the cut-off
    MD_NB_MODIFIER_FORCE_SWITCH,    // The force smoothly to zero from the switch distance, the potential its integral
} md_nb_modifier_t;

typedef enum md_nb_coulomb_t {
    MD_NB_COULOMB_NONE = 0,         // No electrostatics
    MD_NB_COULOMB_CUTOFF,           // f q q / (eps_r r), modifier NONE or POT_SHIFT
    MD_NB_COULOMB_REACTION_FIELD,   // f q q / eps_r (1/r + k_rf r^2 - c_rf), zero at the cut-off
    MD_NB_COULOMB_EWALD,            // Real space Ewald: f q q / eps_r erfc(beta r) / r, modifier NONE or POT_SHIFT
} md_nb_coulomb_t;

typedef struct md_nb_desc_t {
    md_nb_modifier_t lj_modifier;
    double lj_cutoff;               // nm
    double lj_switch;               // nm, where POT_SWITCH and FORCE_SWITCH begin

    md_nb_coulomb_t coulomb;
    md_nb_modifier_t coulomb_modifier;  // NONE or POT_SHIFT, for CUTOFF and EWALD
    double coulomb_cutoff;          // nm
    double epsilon_r;               // Relative dielectric constant, 0 for infinity (no electrostatics)
    double epsilon_rf;              // Reaction field dielectric constant, 0 for infinity
    double ewald_rtol;              // Ewald: erfc(beta rc), which sets beta
} md_nb_desc_t;

// Prepared from a md_nb_desc_t. Read only, so shared between threads.
typedef struct md_nb_potential_t {
    int32_t lj_modifier;
    double lj_cutoff2;
    double lj_switch;
    // Force switch: potential of r^-p is r^-p + p (-c2/3 - c3/4 s) s^3 + cpot, s = max(r - switch, 0)
    double disp_c2, disp_c3, disp_cpot;     // p = 6
    double rep_c2, rep_c3, rep_cpot;        // p = 12
    // Potential switch: S = 1 + c3 s^3 + c4 s^4 + c5 s^5
    double sw_c3, sw_c4, sw_c5;

    int32_t coulomb;
    double coulomb_cutoff2;
    double epsfac;                  // 1 / (4 pi eps0 eps_r), kJ/mol nm / e^2
    double k_rf, c_rf;              // CUTOFF and REACTION_FIELD: 1/r + k_rf r^2 - c_rf
    double ewald_beta, ewald_shift; // EWALD: erfc(beta r)/r - shift
} md_nb_potential_t;

// 1 / (4 pi eps0) in kJ/mol nm / e^2, from the 2018 CODATA constants as GROMACS computes it
#define MD_NB_ONE_4PI_EPS0 138.93545764438198

bool md_nb_potential_init(md_nb_potential_t* pot, const md_nb_desc_t* desc);

// The potential mdrun used for a system: from the non-bonded settings of a tpr. Fails (and says why) for
// what is not a sum over pairs or not modelled: LJ-PME, tabulated potentials, a repulsion other than r^-12,
// Buckingham, and the group scheme's obsolete forms.
bool md_nb_potential_init_from_tpr(md_nb_potential_t* pot, const struct md_tpr_data_t* tpr);

// The furthest distance at which a pair interacts, nm
double md_nb_potential_cutoff(const md_nb_potential_t* pot);

// Lennard-Jones energy of a pair at squared distance r2 (nm^2), kJ/mol. Zero from the cut-off on.
double md_nb_lj_energy(const md_nb_potential_t* pot, double c6, double c12, double r2);

// Coulomb energy of a pair with charge product qq (e^2) at squared distance r2 (nm^2), kJ/mol. Zero from the cut-off on.
double md_nb_coulomb_energy(const md_nb_potential_t* pot, double qq, double r2);

// ### FORCE FIELD OF A SYSTEM ###
// The non-bonded force field of a system: what each particle is, how pairs of them interact, and which pairs do not.
// A system carries one (md_system_t::nonbonded) when it was loaded with it, today from a tpr. It is NULL otherwise,
// and then there are no energies: the parameters alone (an itp) do not say how the simulation cut them off.
typedef struct md_nb_forcefield_t {
    md_nb_potential_t potential;

    size_t    num_types;
    float*    c6;               // [num_types * num_types], kJ/mol nm^6
    float*    c12;              // [num_types * num_types], kJ/mol nm^12

    size_t    num_atoms;
    uint16_t* type;             // [num_atoms]
    float*    charge;           // [num_atoms], e

    // Exclusions, stored per molecule type: particle k is particle k - mol_beg[k] of a molecule of type mol_type[k].
    // The particles excluded from local particle i of type t are excl[excl_base[t] + excl_off[off_base[t] + i]] ..
    // up to the next offset. Types without exclusions have off_base[t] == UINT32_MAX.
    uint32_t* mol_beg;          // [num_atoms]
    uint32_t* mol_type;         // [num_atoms]
    size_t    num_mol_types;
    uint32_t* off_base;         // [num_mol_types]
    uint32_t* excl_base;        // [num_mol_types]
    uint32_t* excl_off;
    uint32_t* excl;
    size_t    excl_off_count;
    size_t    excl_count;

    struct md_allocator_i* alloc;
} md_nb_forcefield_t;

// From a tpr. Fails, quietly, when its non-bonded interactions are not supported by md_nb_potential_t.
bool md_nb_forcefield_init_from_tpr(md_nb_forcefield_t* ff, const struct md_tpr_data_t* tpr, struct md_allocator_i* alloc);
void md_nb_forcefield_free(md_nb_forcefield_t* ff);

// Whether the force field excludes the pair from non-bonded interactions (same molecule, bonded neighbours)
bool md_nb_forcefield_excluded(const md_nb_forcefield_t* ff, uint32_t a, uint32_t b);

// The Lennard-Jones and Coulomb energy (kJ/mol) of particles a and b at squared distance r2_nm (nm^2). Zero for
// excluded pairs.
void md_nb_forcefield_pair_energy(const md_nb_forcefield_t* ff, uint32_t a, uint32_t b, double r2_nm, double* e_lj, double* e_coul);

#ifdef __cplusplus
}
#endif
