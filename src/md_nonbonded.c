#include <md_nonbonded.h>
#include <md_tpr.h>

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
    MEMSET(pot, 0, sizeof(md_nb_potential_t));

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
    return sqrt(MAX(pot->lj_cutoff2, pot->coulomb != MD_NB_COULOMB_NONE ? pot->coulomb_cutoff2 : 0.0));
}

double md_nb_lj_energy(const md_nb_potential_t* pot, double c6, double c12, double r2) {
    if (!(r2 < pot->lj_cutoff2) || !(r2 > 0)) return 0;
    const double rinv6 = 1 / (r2 * r2 * r2);
    double v = c12 * (rinv6 * rinv6 + pot->rep_cpot) - c6 * (rinv6 + pot->disp_cpot);
    if (pot->lj_modifier == MD_NB_MODIFIER_FORCE_SWITCH || pot->lj_modifier == MD_NB_MODIFIER_POT_SWITCH) {
        const double s = MAX(sqrt(r2) - pot->lj_switch, 0.0);
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

bool md_nb_potential_init_from_tpr(md_nb_potential_t* pot, const md_tpr_data_t* tpr) {
    if (!pot || !tpr) return false;
    MEMSET(pot, 0, sizeof(md_nb_potential_t));
    const md_tpr_nonbonded_t* nb = &tpr->nonbonded;
    if (!nb->valid) {
        MD_LOG_INFO("Non-bonded potential: the tpr has no readable simulation parameters");
        return false;
    }
    if (!tpr->nb_is_lj) {
        MD_LOG_INFO("Non-bonded potential: the force field is not Lennard-Jones");
        return false;
    }
    if (tpr->repulsion_power != 12.0) {
        MD_LOG_INFO("Non-bonded potential: a repulsion of r^-%g is not supported", tpr->repulsion_power);
        return false;
    }
    const bool verlet = nb->cutoff_scheme == MD_TPR_CUTOFF_SCHEME_VERLET;

    md_nb_desc_t desc = {
        .lj_cutoff = nb->rvdw,
        .lj_switch = nb->rvdw_switch,
        .coulomb_cutoff = nb->rcoulomb,
        .epsilon_r = nb->epsilon_r,
        .epsilon_rf = nb->epsilon_rf,
        .ewald_rtol = nb->ewald_rtol,
    };

    // The modifier as mdrun applies it. Before the Verlet scheme had modifiers, the default meant a shift for
    // Verlet and nothing for the group scheme.
    #define MAP_MODIFIER(out, in)                                                                       \
        switch (in) {                                                                                   \
        case MD_TPR_MODIFIER_POT_SHIFT_VERLET_UNSUPPORTED: out = verlet ? MD_NB_MODIFIER_POT_SHIFT : MD_NB_MODIFIER_NONE; break; \
        case MD_TPR_MODIFIER_POT_SHIFT:    out = MD_NB_MODIFIER_POT_SHIFT; break;                      \
        case MD_TPR_MODIFIER_NONE:                                                                      \
        case MD_TPR_MODIFIER_EXACT_CUTOFF: out = MD_NB_MODIFIER_NONE; break;                           \
        case MD_TPR_MODIFIER_POT_SWITCH:   out = MD_NB_MODIFIER_POT_SWITCH; break;                     \
        case MD_TPR_MODIFIER_FORCE_SWITCH: out = MD_NB_MODIFIER_FORCE_SWITCH; break;                   \
        default: MD_LOG_ERROR("Non-bonded potential: unknown modifier %d", (int)(in)); return false;    \
        }

    switch (nb->vdw_type) {
    case MD_TPR_VDW_CUT:
        MAP_MODIFIER(desc.lj_modifier, nb->vdw_modifier);
        break;
    case MD_TPR_VDW_SWITCH:     // Group scheme: what grompp turns into a potential switch for Verlet
        desc.lj_modifier = MD_NB_MODIFIER_POT_SWITCH;
        break;
    case MD_TPR_VDW_SHIFT:      // Group scheme: what grompp turns into a force switch for Verlet
        desc.lj_modifier = MD_NB_MODIFIER_FORCE_SWITCH;
        break;
    case MD_TPR_VDW_PME:
        MD_LOG_INFO("Non-bonded potential: LJ-PME is not a sum over pairs (its grid part is not modelled)");
        return false;
    case MD_TPR_VDW_USER:
        MD_LOG_INFO("Non-bonded potential: tabulated Lennard-Jones, the tables are not part of the tpr");
        return false;
    default:
        MD_LOG_INFO("Non-bonded potential: van der Waals type %d is not supported", (int)nb->vdw_type);
        return false;
    }

    md_nb_modifier_t coulomb_modifier = MD_NB_MODIFIER_NONE;
    MAP_MODIFIER(coulomb_modifier, nb->coulomb_modifier);
    #undef MAP_MODIFIER
    desc.coulomb_modifier = coulomb_modifier;

    switch (nb->coulomb_type) {
    case MD_TPR_COULOMB_CUT:
        desc.coulomb = MD_NB_COULOMB_CUTOFF;
        if (!verlet) desc.coulomb_modifier = MD_NB_MODIFIER_NONE;
        break;
    case MD_TPR_COULOMB_RF:
        desc.coulomb = MD_NB_COULOMB_REACTION_FIELD;
        break;
    case MD_TPR_COULOMB_PME:
    case MD_TPR_COULOMB_EWALD:
    case MD_TPR_COULOMB_P3M_AD:
        desc.coulomb = MD_NB_COULOMB_EWALD;
        break;
    default:
        MD_LOG_INFO("Non-bonded potential: Coulomb type %d is not supported", (int)nb->coulomb_type);
        return false;
    }
    if (desc.coulomb_modifier != MD_NB_MODIFIER_NONE && desc.coulomb_modifier != MD_NB_MODIFIER_POT_SHIFT) {
        MD_LOG_INFO("Non-bonded potential: switched Coulomb is not supported");
        return false;
    }

    return md_nb_potential_init(pot, &desc);
}

// ### FORCE FIELD OF A SYSTEM ###

void md_nb_forcefield_free(md_nb_forcefield_t* ff) {
    if (!ff || !ff->alloc) return;
    md_allocator_i* a = ff->alloc;
    const size_t T2 = ff->num_types * ff->num_types;
    if (ff->c6)        md_free(a, ff->c6,  sizeof(float) * MAX(T2, 1));
    if (ff->c12)       md_free(a, ff->c12, sizeof(float) * MAX(T2, 1));
    if (ff->type)      md_free(a, ff->type,     sizeof(uint16_t) * MAX(ff->num_atoms, 1));
    if (ff->charge)    md_free(a, ff->charge,   sizeof(float)    * MAX(ff->num_atoms, 1));
    if (ff->mol_beg)   md_free(a, ff->mol_beg,  sizeof(uint32_t) * MAX(ff->num_atoms, 1));
    if (ff->mol_type)  md_free(a, ff->mol_type, sizeof(uint32_t) * MAX(ff->num_atoms, 1));
    if (ff->off_base)  md_free(a, ff->off_base,  sizeof(uint32_t) * MAX(ff->num_mol_types, 1));
    if (ff->excl_base) md_free(a, ff->excl_base, sizeof(uint32_t) * MAX(ff->num_mol_types, 1));
    if (ff->excl_off)  md_free(a, ff->excl_off, sizeof(uint32_t) * MAX(ff->excl_off_count, 1));
    if (ff->excl)      md_free(a, ff->excl,     sizeof(uint32_t) * MAX(ff->excl_count, 1));
    MEMSET(ff, 0, sizeof(md_nb_forcefield_t));
}

bool md_nb_forcefield_init_from_tpr(md_nb_forcefield_t* ff, const md_tpr_data_t* tpr, md_allocator_i* alloc) {
    if (!ff || !tpr || !alloc) return false;
    MEMSET(ff, 0, sizeof(md_nb_forcefield_t));
    if (!md_nb_potential_init_from_tpr(&ff->potential, tpr)) return false;
    if (tpr->num_atoms > UINT32_MAX - 1) return false;

    const size_t T = tpr->num_nb_types, N = tpr->num_atoms, M = tpr->num_moltypes;
    ff->alloc = alloc;
    ff->num_types = T;
    ff->num_atoms = N;
    ff->num_mol_types = M;
    ff->c6  = md_alloc(alloc, sizeof(float) * MAX(T * T, 1));
    ff->c12 = md_alloc(alloc, sizeof(float) * MAX(T * T, 1));
    for (size_t i = 0; i < T * T; ++i) {
        ff->c6[i]  = tpr->lj[i].c6;
        ff->c12[i] = tpr->lj[i].c12;
    }

    // Exclusions of the molecule types, concatenated
    ff->off_base  = md_alloc(alloc, sizeof(uint32_t) * MAX(M, 1));
    ff->excl_base = md_alloc(alloc, sizeof(uint32_t) * MAX(M, 1));
    size_t num_off = 0, num_excl = 0;
    for (size_t t = 0; t < M; ++t) {
        const md_tpr_moltype_t* mt = &tpr->moltypes[t];
        if (mt->excl_offset) {
            num_off  += mt->num_atoms + 1;
            num_excl += mt->excl_offset[mt->num_atoms];
        }
    }
    ff->excl_off_count = num_off;
    ff->excl_count = num_excl;
    ff->excl_off = md_alloc(alloc, sizeof(uint32_t) * MAX(num_off, 1));
    ff->excl     = md_alloc(alloc, sizeof(uint32_t) * MAX(num_excl, 1));
    num_off = num_excl = 0;
    for (size_t t = 0; t < M; ++t) {
        const md_tpr_moltype_t* mt = &tpr->moltypes[t];
        if (!mt->excl_offset) {
            ff->off_base[t] = UINT32_MAX;
            ff->excl_base[t] = 0;
            continue;
        }
        ff->off_base[t] = (uint32_t)num_off;
        ff->excl_base[t] = (uint32_t)num_excl;
        MEMCPY(ff->excl_off + num_off, mt->excl_offset, sizeof(uint32_t) * (mt->num_atoms + 1));
        MEMCPY(ff->excl + num_excl, mt->excl, sizeof(uint32_t) * mt->excl_offset[mt->num_atoms]);
        num_off  += mt->num_atoms + 1;
        num_excl += mt->excl_offset[mt->num_atoms];
    }

    ff->type     = md_alloc(alloc, sizeof(uint16_t) * MAX(N, 1));
    ff->charge   = md_alloc(alloc, sizeof(float)    * MAX(N, 1));
    ff->mol_beg  = md_alloc(alloc, sizeof(uint32_t) * MAX(N, 1));
    ff->mol_type = md_alloc(alloc, sizeof(uint32_t) * MAX(N, 1));
    size_t n = 0;
    for (size_t b = 0; b < tpr->num_molblocks; ++b) {
        const uint32_t t = (uint32_t)tpr->molblocks[b].moltype;
        const md_tpr_moltype_t* mt = &tpr->moltypes[t];
        for (int32_t m = 0; m < tpr->molblocks[b].nmol; ++m) {
            const uint32_t beg = (uint32_t)n;
            for (size_t a = 0; a < mt->num_atoms && n < N; ++a, ++n) {
                ff->type[n]     = mt->atoms[a].type_idx;
                ff->charge[n]   = mt->atoms[a].charge;
                ff->mol_beg[n]  = beg;
                ff->mol_type[n] = t;
            }
        }
    }
    if (n != N) {
        md_nb_forcefield_free(ff);
        return false;
    }
    return true;
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
