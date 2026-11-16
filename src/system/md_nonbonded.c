#include <md_nonbonded.h>

#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_allocator.h>

#include <math.h>
#include <string.h>

// The forms and constants are those of GROMACS' interaction_const.cpp (force_switch_constants,
// potential_switch_constants, calc_rffac) and ewald_utils.cpp (calc_ewaldcoeff_q).

static void force_switch_constants(double p, double rsw, double rc, double* c2, double* c3, double* cpot) {
    *c2   = ((p + 1) * rsw - (p + 4) * rc) / (pow(rc, p + 2) * (rc - rsw) * (rc - rsw));
    *c3   = -((p + 1) * rsw - (p + 3) * rc) / (pow(rc, p + 2) * (rc - rsw) * (rc - rsw) * (rc - rsw));
    *cpot = -pow(rc, -p) + p * *c2 / 3 * pow(rc - rsw, 3) + p * *c3 / 4 * pow(rc - rsw, 4);
}

// The beta for which erfc(beta rc) = rtol, by bisection as GROMACS does it
static double ewald_coefficient(double rc, double rtol) {
    double beta = 5;
    int i = 0;
    do {
        i++;
        beta *= 2;
    } while (erfc(beta * rc) > rtol && i < 64);
    const int n = i + 60;
    double lo = 0, hi = beta;
    for (int k = 0; k < n; ++k) {
        beta = (lo + hi) / 2;
        if (erfc(beta * rc) > rtol) lo = beta; else hi = beta;
    }
    return beta;
}

bool md_nb_potential_init(md_nb_potential_t* pot, const md_nb_desc_t* desc) {
    if (!pot || !desc) return false;
    md_memset(pot, 0, sizeof(md_nb_potential_t));

    // Lennard-Jones
    if (!(desc->lj_cutoff > 0)) {
        MD_LOG_ERROR("Non-bonded potential: the Lennard-Jones cut-off has to be positive");
        return false;
    }
    const double rc = desc->lj_cutoff;
    const double rsw = desc->lj_switch;
    pot->lj_modifier = desc->lj_modifier;
    pot->lj_cutoff2 = rc * rc;
    pot->lj_switch = rsw;
    switch (desc->lj_modifier) {
    case MD_NB_MODIFIER_NONE:
        break;
    case MD_NB_MODIFIER_POT_SHIFT:
        pot->disp_cpot = -pow(rc, -6);
        pot->rep_cpot  = -pow(rc, -12);
        break;
    case MD_NB_MODIFIER_FORCE_SWITCH:
    case MD_NB_MODIFIER_POT_SWITCH:
        if (!(rsw >= 0 && rsw < rc)) {
            MD_LOG_ERROR("Non-bonded potential: the switch distance (%g) has to be below the cut-off (%g)", rsw, rc);
            return false;
        }
        if (desc->lj_modifier == MD_NB_MODIFIER_FORCE_SWITCH) {
            force_switch_constants(6.0,  rsw, rc, &pot->disp_c2, &pot->disp_c3, &pot->disp_cpot);
            force_switch_constants(12.0, rsw, rc, &pot->rep_c2,  &pot->rep_c3,  &pot->rep_cpot);
        } else {
            const double w = rc - rsw;
            pot->sw_c3 = -10 / (w * w * w);
            pot->sw_c4 =  15 / (w * w * w * w);
            pot->sw_c5 =  -6 / (w * w * w * w * w);
        }
        break;
    default:
        MD_LOG_ERROR("Non-bonded potential: unknown Lennard-Jones modifier %d", (int)desc->lj_modifier);
        return false;
    }

    // Coulomb
    pot->coulomb = desc->coulomb;
    if (desc->coulomb == MD_NB_COULOMB_NONE || desc->epsilon_r == 0) {
        // eps_r = 0 is an infinite dielectric: no electrostatics
        pot->coulomb = MD_NB_COULOMB_NONE;
        return true;
    }
    if (!(desc->coulomb_cutoff > 0) || !(desc->epsilon_r > 0)) {
        MD_LOG_ERROR("Non-bonded potential: the Coulomb cut-off and dielectric constant have to be positive");
        return false;
    }
    const double rcc = desc->coulomb_cutoff;
    pot->coulomb_cutoff2 = rcc * rcc;
    pot->epsfac = MD_NB_ONE_4PI_EPS0 / desc->epsilon_r;
    const bool shift = desc->coulomb_modifier == MD_NB_MODIFIER_POT_SHIFT;
    if (!shift && desc->coulomb_modifier != MD_NB_MODIFIER_NONE && desc->coulomb != MD_NB_COULOMB_REACTION_FIELD) {
        MD_LOG_ERROR("Non-bonded potential: Coulomb is shifted or not modified, switches are not supported");
        return false;
    }
    switch (desc->coulomb) {
    case MD_NB_COULOMB_CUTOFF:
        pot->k_rf = 0;
        pot->c_rf = shift ? 1 / rcc : 0;
        break;
    case MD_NB_COULOMB_REACTION_FIELD:
        if (desc->epsilon_rf == 0) {
            pot->k_rf = 1 / (2 * rcc * rcc * rcc);
        } else {
            pot->k_rf = (desc->epsilon_rf - desc->epsilon_r) / (2 * desc->epsilon_rf + desc->epsilon_r) / (rcc * rcc * rcc);
        }
        pot->c_rf = 1 / rcc + pot->k_rf * rcc * rcc;
        break;
    case MD_NB_COULOMB_EWALD:
        if (!(desc->ewald_rtol > 0 && desc->ewald_rtol < 1)) {
            MD_LOG_ERROR("Non-bonded potential: the Ewald tolerance has to be in (0, 1)");
            return false;
        }
        pot->ewald_beta = ewald_coefficient(rcc, desc->ewald_rtol);
        pot->ewald_shift = shift ? erfc(pot->ewald_beta * rcc) / rcc : 0;
        break;
    default:
        MD_LOG_ERROR("Non-bonded potential: unknown Coulomb form %d", (int)desc->coulomb);
        return false;
    }
    return true;
}

double md_nb_potential_cutoff(const md_nb_potential_t* pot) {
    if (!pot) return 0;
    return sqrt(MD_MAX(pot->lj_cutoff2, pot->coulomb != MD_NB_COULOMB_NONE ? pot->coulomb_cutoff2 : 0.0));
}

double md_nb_lj_energy(const md_nb_potential_t* pot, double c6, double c12, double r2) {
    if (!(r2 < pot->lj_cutoff2) || !(r2 > 0)) return 0;
    const double rinv6 = 1 / (r2 * r2 * r2);
    double v = c12 * (rinv6 * rinv6 + pot->rep_cpot) - c6 * (rinv6 + pot->disp_cpot);
    if (pot->lj_modifier == MD_NB_MODIFIER_FORCE_SWITCH || pot->lj_modifier == MD_NB_MODIFIER_POT_SWITCH) {
        const double s = MD_MAX(sqrt(r2) - pot->lj_switch, 0.0);
        if (pot->lj_modifier == MD_NB_MODIFIER_FORCE_SWITCH) {
            const double s3 = s * s * s;
            v += -6 * c6 * (-pot->disp_c2 / 3 - pot->disp_c3 / 4 * s) * s3
               + 12 * c12 * (-pot->rep_c2 / 3 - pot->rep_c3 / 4 * s) * s3;
        } else {
            v *= 1 + (pot->sw_c3 + (pot->sw_c4 + pot->sw_c5 * s) * s) * s * s * s;
        }
    }
    return v;
}

double md_nb_coulomb_energy(const md_nb_potential_t* pot, double qq, double r2) {
    if (pot->coulomb == MD_NB_COULOMB_NONE || !(r2 < pot->coulomb_cutoff2) || !(r2 > 0)) return 0;
    const double r = sqrt(r2);
    if (pot->coulomb == MD_NB_COULOMB_EWALD) {
        return pot->epsfac * qq * (erfc(pot->ewald_beta * r) / r - pot->ewald_shift);
    }
    return pot->epsfac * qq * (1 / r + pot->k_rf * r2 - pot->c_rf);
}

// ### FORCE FIELD OF A SYSTEM ###

void md_nb_forcefield_free(md_nb_forcefield_t* ff) {
    if (!ff || !ff->alloc) return;
    md_allocator_i* a = ff->alloc;
    const size_t T2 = ff->num_types * ff->num_types;
    if (ff->c6)        md_free(a, ff->c6,  sizeof(float) * MD_MAX(T2, 1));
    if (ff->c12)       md_free(a, ff->c12, sizeof(float) * MD_MAX(T2, 1));
    if (ff->type)      md_free(a, ff->type,     sizeof(uint16_t) * MD_MAX(ff->num_atoms, 1));
    if (ff->charge)    md_free(a, ff->charge,   sizeof(float)    * MD_MAX(ff->num_atoms, 1));
    if (ff->mol_beg)   md_free(a, ff->mol_beg,  sizeof(uint32_t) * MD_MAX(ff->num_atoms, 1));
    if (ff->mol_type)  md_free(a, ff->mol_type, sizeof(uint32_t) * MD_MAX(ff->num_atoms, 1));
    if (ff->off_base)  md_free(a, ff->off_base,  sizeof(uint32_t) * MD_MAX(ff->num_mol_types, 1));
    if (ff->excl_base) md_free(a, ff->excl_base, sizeof(uint32_t) * MD_MAX(ff->num_mol_types, 1));
    if (ff->excl_off)  md_free(a, ff->excl_off, sizeof(uint32_t) * MD_MAX(ff->excl_off_count, 1));
    if (ff->excl)      md_free(a, ff->excl,     sizeof(uint32_t) * MD_MAX(ff->excl_count, 1));
    md_memset(ff, 0, sizeof(md_nb_forcefield_t));
}

bool md_nb_forcefield_excluded(const md_nb_forcefield_t* ff, uint32_t a, uint32_t b) {
    if (a == b) return true;
    if (ff->mol_beg[a] != ff->mol_beg[b]) return false;
    const uint32_t t = ff->mol_type[a];
    if (ff->off_base[t] == UINT32_MAX) return false;
    const uint32_t la = a - ff->mol_beg[a];
    const uint32_t lb = b - ff->mol_beg[a];
    const uint32_t* off = ff->excl_off + ff->off_base[t];
    const uint32_t* ex  = ff->excl + ff->excl_base[t];
    uint32_t lo = off[la], hi = off[la + 1];
    while (lo < hi) {
        const uint32_t mid = (lo + hi) / 2;
        if (ex[mid] == lb) return true;
        if (ex[mid] < lb) lo = mid + 1; else hi = mid;
    }
    return false;
}

void md_nb_forcefield_pair_energy(const md_nb_forcefield_t* ff, uint32_t a, uint32_t b, double r2_nm, double* e_lj, double* e_coul) {
    *e_lj = 0;
    *e_coul = 0;
    if (md_nb_forcefield_excluded(ff, a, b)) return;
    const size_t idx = (size_t)ff->type[a] * ff->num_types + ff->type[b];
    *e_lj   = md_nb_lj_energy(&ff->potential, ff->c6[idx], ff->c12[idx], r2_nm);
    *e_coul = md_nb_coulomb_energy(&ff->potential, (double)ff->charge[a] * ff->charge[b], r2_nm);
}
