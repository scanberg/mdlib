#include <md_hbond.h>

#include <md_system.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_coord_stream.h>
#include <core/md_log.h>
#include <core/md_spatial_acc.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>

// ### PARAMETERS ###

md_hbond_params_t md_hbond_params_preset(md_hbond_preset_t preset) {
    // The tools' own presets take every N and O as an acceptor, apply no exclusions and no competition: they report
    // what passes the geometry, which is what their numbers are made of.
    md_hbond_params_t p = {
        .roles             = MD_HBOND_ROLES_ALL_N_O,
        .acc_capacity_mode = MD_HBOND_CAPACITY_UNLIMITED,
    };
    switch (preset) {
    case MD_HBOND_PRESET_MDTRAJ:
        p.max_ha  = 2.5f;
        p.min_dha = 120.0f;
        break;
    case MD_HBOND_PRESET_MDANALYSIS:
        p.max_da  = 3.0f;
        p.min_dha = 150.0f;
        break;
    case MD_HBOND_PRESET_GROMACS:
        p.max_da  = 3.5f;
        p.max_hda = 30.0f;
        break;
    case MD_HBOND_PRESET_VMD:
        p.max_da  = 3.0f;
        p.min_dha = 160.0f;
        break;
    case MD_HBOND_PRESET_REALISTIC:
    default:
        p = (md_hbond_params_t){
            .roles             = MD_HBOND_ROLES_HALIDE_IONS,
            .exclude_bonds     = 3,
            .max_ha            = 2.5f,
            .min_dha           = 120.0f,
            .min_xah           = 90.0f,
            .h_capacity        = 1,
            .acc_capacity_mode = MD_HBOND_CAPACITY_LONE_PAIRS,
        };
        break;
    }
    return p;
}

const char* md_hbond_preset_name(md_hbond_preset_t preset) {
    switch (preset) {
    case MD_HBOND_PRESET_REALISTIC:  return "Realistic";
    case MD_HBOND_PRESET_MDTRAJ:     return "MDTraj (Baker-Hubbard)";
    case MD_HBOND_PRESET_MDANALYSIS: return "MDAnalysis";
    case MD_HBOND_PRESET_GROMACS:    return "GROMACS";
    case MD_HBOND_PRESET_VMD:        return "VMD";
    default:                         return "";
    }
}

static bool params_validate(const md_hbond_params_t* p) {
    if (p->max_ha < 0 || p->max_da < 0 || p->min_dha < 0 || p->max_hda < 0 || p->min_xah < 0 || p->min_strength < 0 || p->bifurcation_tol < 0) {
        MD_LOG_ERROR("Hydrogen bonds: negative parameter");
        return false;
    }
    if (p->max_ha <= 0 && p->max_da <= 0) {
        MD_LOG_ERROR("Hydrogen bonds: at least one of the distance gates (max_ha, max_da) is required");
        return false;
    }
    if (p->min_dha > 180 || p->max_hda > 180 || p->min_xah > 180) {
        MD_LOG_ERROR("Hydrogen bonds: angles are in degrees, within [0, 180]");
        return false;
    }
    if (p->min_strength > 1 || p->bifurcation_tol > 1) {
        MD_LOG_ERROR("Hydrogen bonds: min_strength and bifurcation_tol are within [0, 1]");
        return false;
    }
    if ((uint32_t)p->acc_capacity_mode > MD_HBOND_CAPACITY_UNLIMITED) {
        MD_LOG_ERROR("Hydrogen bonds: unknown acceptor capacity mode");
        return false;
    }
    return true;
}

// ### ROLES ###

static inline bool is_metal(md_atomic_number_t z) {
    return (z == 3 || z == 4) ||                // Li, Be
           (z >= 11 && z <= 13) ||              // Na, Mg, Al
           (z >= 19 && z <= 31) ||              // K .. Ga
           (z >= 37 && z <= 50) ||              // Rb .. Sn
           (z >= 55 && z <= 84) ||              // Cs .. Po
           (z >= 87);
}

static inline bool is_halogen(md_atomic_number_t z) {
    return z == MD_Z_F || z == MD_Z_Cl || z == MD_Z_Br || z == MD_Z_I;
}

// Whether the bond from an atom to a neighbour counts towards its valence: not to virtual sites (no element, such as
// the M site of TIP4P), not to metals and not coordinate bonds, which do not take a hydrogen's place.
static inline bool counts_as_neighbour(const md_system_t* sys, md_bond_iter_t* it) {
    const md_bond_flags_t f = (md_bond_flags_t)md_bond_iter_bond_flags(it);
    if (f & (MD_BOND_FLAG_COORDINATE | MD_BOND_FLAG_METAL)) return false;
    const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, md_bond_iter_atom_index(it));
    return z != 0 && !is_metal(z);
}

static void count_neighbours(int* out_h, int* out_heavy, const md_system_t* sys, size_t i) {
    int nh = 0, nx = 0;
    md_bond_iter_t it = md_bond_iter(&sys->bond, i);
    while (md_bond_iter_has_next(&it)) {
        if (counts_as_neighbour(sys, &it)) {
            if (md_atom_atomic_number(&sys->atom, md_bond_iter_atom_index(&it)) == MD_Z_H) nh += 1;
            else nx += 1;
        }
        md_bond_iter_next(&it);
    }
    *out_h = nh;
    *out_heavy = nx;
}

static inline int neighbour_degree(const md_system_t* sys, size_t i) {
    int nh, nx;
    count_neighbours(&nh, &nx, sys, i);
    return nh + nx;
}

// Sum of the three bond angles around a three connected atom, in degrees. 360 is planar, 328.4 tetrahedral.
static float angle_sum(const md_system_t* sys, const md_system_state_t* ref, size_t i) {
    vec3_t v[3];
    int n = 0;
    md_bond_iter_t it = md_bond_iter(&sys->bond, i);
    while (md_bond_iter_has_next(&it) && n < 3) {
        if (counts_as_neighbour(sys, &it)) {
            v[n++] = vec3_sub(ref->xyz[md_bond_iter_atom_index(&it)], ref->xyz[i]);
        }
        md_bond_iter_next(&it);
    }
    if (n < 3) return 360.0f;
    md_util_min_image_vec3(v, 3, &ref->unitcell);
    const float a = vec3_angle(v[0], v[1]) + vec3_angle(v[1], v[2]) + vec3_angle(v[0], v[2]);
    return (float)RAD_TO_DEG(a);
}

// Whether the lone pair of a three connected N is drawn into a pi system, from the bond graph alone: a neighbour which
// is unsaturated (a carbon with fewer than four neighbours, a two connected N), a sulfonyl or phosphoryl group, or a
// multiple or aromatic bond. Assumes explicit hydrogens on the carbons; the geometry is preferred where it is known.
static bool n_conjugated_by_graph(const md_system_t* sys, size_t i) {
    md_bond_iter_t it = md_bond_iter(&sys->bond, i);
    while (md_bond_iter_has_next(&it)) {
        if (counts_as_neighbour(sys, &it)) {
            const md_bond_flags_t f = (md_bond_flags_t)md_bond_iter_bond_flags(&it);
            const size_t x = md_bond_iter_atom_index(&it);
            const md_atomic_number_t zx = md_atom_atomic_number(&sys->atom, x);
            if (f & (MD_BOND_FLAG_DOUBLE | MD_BOND_FLAG_TRIPLE | MD_BOND_FLAG_AROMATIC)) return true;
            if (zx != MD_Z_H) {
                const int deg = neighbour_degree(sys, x);
                if (zx == MD_Z_C && deg <= 3) return true;
                if (zx == MD_Z_N && deg <= 2) return true;
                if ((zx == MD_Z_S || zx == MD_Z_P) && deg >= 3) return true;
                md_bond_iter_t it2 = md_bond_iter(&sys->bond, x);
                while (md_bond_iter_has_next(&it2)) {
                    if (md_bond_iter_bond_flags(&it2) & (MD_BOND_FLAG_DOUBLE | MD_BOND_FLAG_TRIPLE | MD_BOND_FLAG_AROMATIC)) return true;
                    md_bond_iter_next(&it2);
                }
            }
        }
        md_bond_iter_next(&it);
    }
    return false;
}

// Base of a nucleotide residue name: A, C, G, T or U, 0 if unknown. Covers A, DA, RA, DA5, DA3, ADE, ...
static char nucleobase_letter(str_t name) {
    name = str_trim(name);
    if (name.len == 0) return 0;
    const char c0 = name.ptr[0];
    if (name.len >= 3) {
        if (str_eq_cstr(name, "ADE")) return 'A';
        if (str_eq_cstr(name, "CYT")) return 'C';
        if (str_eq_cstr(name, "GUA")) return 'G';
        if (str_eq_cstr(name, "THY")) return 'T';
        if (str_eq_cstr(name, "URA") || str_eq_cstr(name, "URI")) return 'U';
    }
    const char* bases = "ACGTU";
    if (name.len == 1 || (name.len <= 3 && name.ptr[1] >= '0' && name.ptr[1] <= '9')) {
        for (const char* b = bases; *b; ++b) if (c0 == *b) return c0;
        return 0;
    }
    if ((c0 == 'D' || c0 == 'R') && name.len >= 2) {
        const char c1 = name.ptr[1];
        for (const char* b = bases; *b; ++b) if (c1 == *b) return c1;
    }
    return 0;
}

// Nitrogens of the standard residues whose lone pair is never free: the backbone N, the side chain N of Arg, Lys
// (protonated by name), Asn, Gln and Trp, and the amino and glycosidic N of the nucleobases. These are decided by name
// so that they come out right also without hydrogens, where the bond graph cannot tell an NH from a free lone pair.
static bool standard_n_without_lone_pair(str_t res, str_t atom, bool amino_acid, bool nucleotide) {
    res  = str_trim(res);
    atom = str_trim(atom);
    if (amino_acid) {
        if (str_eq_cstr(atom, "N")) return true;
        if (str_eq_cstr(res, "ARG")) return str_eq_cstr(atom, "NE") || str_eq_cstr(atom, "NH1") || str_eq_cstr(atom, "NH2");
        if (str_eq_cstr(res, "LYS")) return str_eq_cstr(atom, "NZ");
        if (str_eq_cstr(res, "ASN")) return str_eq_cstr(atom, "ND2");
        if (str_eq_cstr(res, "GLN")) return str_eq_cstr(atom, "NE2");
        if (str_eq_cstr(res, "TRP")) return str_eq_cstr(atom, "NE1");
        return false;
    }
    if (nucleotide) {
        switch (nucleobase_letter(res)) {
        case 'A': return str_eq_cstr(atom, "N9") || str_eq_cstr(atom, "N6");
        case 'G': return str_eq_cstr(atom, "N9") || str_eq_cstr(atom, "N1") || str_eq_cstr(atom, "N2");
        case 'C': return str_eq_cstr(atom, "N1") || str_eq_cstr(atom, "N4");
        case 'T':
        case 'U': return str_eq_cstr(atom, "N1") || str_eq_cstr(atom, "N3");
        default:  return false;
        }
    }
    return false;
}

static inline uint8_t element_capacity(md_atomic_number_t z) {
    switch (z) {
    case MD_Z_N: return 1;
    case MD_Z_O: return 2;
    case MD_Z_S: return 2;
    case MD_Z_F: return 3;
    default:     return 1;
    }
}

bool md_hbond_perceive_roles(uint8_t* out_role, uint8_t* out_cap, const md_system_t* sys, const md_system_state_t* ref, uint32_t role_flags) {
    if (!out_role || !sys) {
        MD_LOG_ERROR("Hydrogen bond roles: missing output or system");
        return false;
    }
    const size_t N = sys->atom.count;
    MEMSET(out_role, 0, N);
    if (out_cap) MEMSET(out_cap, 0, N);
    if (N == 0) return true;

    const bool naive    = (role_flags & MD_HBOND_ROLES_ALL_N_O) != 0;
    const bool have_ref = md_system_state_has_coords(ref) && ref->num_atoms == N;
    bool have_flags = false;
    if (!have_ref && sys->atom.flags) {
        for (size_t i = 0; i < N; ++i) {
            if (sys->atom.flags[i] & (MD_FLAG_HBOND_DONOR | MD_FLAG_HBOND_ACCEPTOR)) { have_flags = true; break; }
        }
    }

    // A three connected N with a sum of angles at least this is planar: its lone pair is in a pi system. sp3 is 328.4,
    // amides and aromatic NH sit within a few degrees of 360, anilines in between.
    const float planar_angle_sum = 345.0f;

    size_t comp = 0;
    const size_t num_comp = sys->component.atom_offset ? sys->component.count : 0;

    for (size_t i = 0; i < N; ++i) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, i);
        if (!(z == MD_Z_N || z == MD_Z_O || z == MD_Z_S || is_halogen(z))) continue;

        int nh, nx;
        count_neighbours(&nh, &nx, sys, i);
        const int deg = nh + nx;
        uint8_t role = 0;
        uint8_t cap  = 0;

        switch (z) {
        case MD_Z_O:
            if (nh > 0) role |= MD_HBOND_ROLE_DONOR;
            if (naive || deg <= 2) {
                role |= MD_HBOND_ROLE_ACCEPTOR;
                cap = 2;
            }
            break;
        case MD_Z_N: {
            if (nh > 0) role |= MD_HBOND_ROLE_DONOR;
            bool acc = false;
            if (naive) {
                acc = true;
            } else if (deg <= 2) {
                acc = true;     // Pyridine and imine type N, nitriles
            } else if (deg == 3) {
                if (have_ref)        acc = angle_sum(sys, ref, i) < planar_angle_sum;
                else if (have_flags) acc = (sys->atom.flags[i] & MD_FLAG_HBOND_ACCEPTOR) != 0;
                else                 acc = !n_conjugated_by_graph(sys, i);
            }
            // deg >= 4: ammonium, no lone pair
            if (acc && !naive) {
                while (comp < num_comp && sys->component.atom_offset[comp + 1] <= i) ++comp;
                if (comp < num_comp && sys->component.atom_offset[comp] <= i) {
                    const str_t res = md_component_name(&sys->component, comp);
                    const md_flags_t cf = md_component_flags(&sys->component, comp);
                    const bool aa  = (cf & MD_FLAG_AMINO_ACID) || md_util_resname_amino_acid(res);
                    const bool nuc = (cf & MD_FLAG_NUCLEOTIDE) || md_util_resname_nucleotide(res);
                    if ((aa || nuc) && standard_n_without_lone_pair(res, md_atom_name(&sys->atom, i), aa, nuc)) {
                        acc = false;
                    }
                }
            }
            if (acc) {
                role |= MD_HBOND_ROLE_ACCEPTOR;
                cap = 1;
            }
            break;
        }
        case MD_Z_S:
            if (role_flags & MD_HBOND_ROLES_SULFUR) {
                if (nh > 0) role |= MD_HBOND_ROLE_DONOR;
                if (deg <= 2) {
                    role |= MD_HBOND_ROLE_ACCEPTOR;
                    cap = 2;
                }
            }
            break;
        default: // Halogens
            if (deg == 0 && (role_flags & MD_HBOND_ROLES_HALIDE_IONS)) {
                role |= MD_HBOND_ROLE_ACCEPTOR;
                cap = MD_HBOND_CAPACITY_NO_LIMIT;
            } else if (z == MD_Z_F && deg == 1 && (role_flags & MD_HBOND_ROLES_FLUORINE)) {
                role |= MD_HBOND_ROLE_ACCEPTOR;
                cap = 3;
            }
            break;
        }

        out_role[i] = role;
        if (out_cap) out_cap[i] = cap;
    }
    return true;
}

void md_hbond_infer_atom_flags(md_system_t* sys, const md_system_state_t* ref) {
    if (!sys || !sys->atom.flags || sys->atom.count == 0) return;
    const size_t N = sys->atom.count;
    md_temp_scope_t temp = md_temp_begin();
    uint8_t* role = md_alloc(md_temp_allocator(temp), N);
    if (md_hbond_perceive_roles(role, NULL, sys, ref, MD_HBOND_ROLES_DEFAULT | MD_HBOND_ROLES_HALIDE_IONS)) {
        for (size_t i = 0; i < N; ++i) {
            md_flags_t f = sys->atom.flags[i] & ~(MD_FLAG_HBOND_DONOR | MD_FLAG_HBOND_ACCEPTOR);
            if (role[i] & MD_HBOND_ROLE_DONOR)    f |= MD_FLAG_HBOND_DONOR;
            if (role[i] & MD_HBOND_ROLE_ACCEPTOR) f |= MD_FLAG_HBOND_ACCEPTOR;
            sys->atom.flags[i] = f;
        }
    }
    md_temp_end(temp);
}

// ### QUERY ###

void md_hbond_query_free(md_hbond_query_t* q) {
    if (!q || !q->alloc) return;
    md_allocator_i* alloc = q->alloc;
    md_array_free(q->donor_d, alloc);
    md_array_free(q->donor_h, alloc);
    md_array_free(q->excl_off, alloc);
    md_array_free(q->excl_atom, alloc);
    md_array_free(q->acceptor, alloc);
    md_array_free(q->acceptor_cap, alloc);
    md_array_free(q->acc_nbr_off, alloc);
    md_array_free(q->acc_nbr, alloc);
    md_array_free(q->sel, alloc);
    MEMSET(q, 0, sizeof(md_hbond_query_t));
}

static bool apply_override(uint8_t* role, const md_bitfield_t* set, uint8_t bit, size_t N) {
    if (md_bitfield_end_bit(set) > N) {
        MD_LOG_ERROR("Hydrogen bonds: role override refers to atoms beyond the system (%zu atoms)", N);
        return false;
    }
    for (size_t i = 0; i < N; ++i) role[i] &= (uint8_t)~bit;
    md_bitfield_iter_t it = md_bitfield_iter_create(set);
    while (md_bitfield_iter_next(&it)) {
        role[md_bitfield_iter_idx(&it)] |= bit;
    }
    return true;
}

static bool mark_selection(uint8_t* sel, const md_bitfield_t* set, uint8_t bit, size_t N) {
    if (md_bitfield_end_bit(set) > N) {
        MD_LOG_ERROR("Hydrogen bonds: selection refers to atoms beyond the system (%zu atoms)", N);
        return false;
    }
    md_bitfield_iter_t it = md_bitfield_iter_create(set);
    while (md_bitfield_iter_next(&it)) {
        sel[md_bitfield_iter_idx(&it)] |= bit;
    }
    return true;
}

static int compare_u32(const void* a, const void* b) {
    const uint32_t x = *(const uint32_t*)a;
    const uint32_t y = *(const uint32_t*)b;
    return (x > y) - (x < y);
}

// For every donor, the acceptors within max_bonds bonds of its heavy atom, as sorted compressed rows. Donors are
// grouped by heavy atom, so consecutive donors of the same atom repeat the row.
static void build_exclusions(md_hbond_query_t* q, const md_system_t* sys, const uint8_t* is_acceptor, uint32_t max_bonds, md_allocator_i* temp_alloc) {
    const size_t N = q->num_atoms;
    md_allocator_i* alloc = q->alloc;
    uint8_t* depth = md_alloc(temp_alloc, N);
    MEMSET(depth, 0, N);
    md_array(uint32_t) queue = 0;
    md_array(uint32_t) row   = 0;
    const uint32_t depth_limit = MIN(max_bonds, 254) + 1;

    md_array_resize(q->excl_off, q->num_donors + 1, alloc);
    q->excl_off[0] = 0;
    uint32_t prev_d = UINT32_MAX;

    for (size_t k = 0; k < q->num_donors; ++k) {
        const uint32_t d = q->donor_d[k];
        if (d != prev_d) {
            md_array_shrink(queue, 0);
            md_array_shrink(row, 0);
            depth[d] = 1;
            md_array_push(queue, d, temp_alloc);
            for (size_t head = 0; head < md_array_size(queue); ++head) {
                const uint32_t cur = queue[head];
                if (depth[cur] >= depth_limit) continue;
                md_bond_iter_t it = md_bond_iter(&sys->bond, cur);
                while (md_bond_iter_has_next(&it)) {
                    const uint32_t next = (uint32_t)md_bond_iter_atom_index(&it);
                    md_bond_iter_next(&it);
                    if (next < N && depth[next] == 0) {
                        depth[next] = depth[cur] + 1;
                        md_array_push(queue, next, temp_alloc);
                        if (is_acceptor[next]) md_array_push(row, next, temp_alloc);
                    }
                }
            }
            for (size_t v = 0; v < md_array_size(queue); ++v) depth[queue[v]] = 0;
            if (md_array_size(row) > 1) qsort(row, md_array_size(row), sizeof(uint32_t), compare_u32);
            prev_d = d;
        }
        for (size_t r = 0; r < md_array_size(row); ++r) md_array_push(q->excl_atom, row[r], alloc);
        q->excl_off[k + 1] = (uint32_t)md_array_size(q->excl_atom);
    }
}

bool md_hbond_query_init(md_hbond_query_t* q, const md_hbond_desc_t* desc, const md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(q);
    MEMSET(q, 0, sizeof(md_hbond_query_t));
    if (!sys || !alloc) {
        MD_LOG_ERROR("Hydrogen bonds: missing system or allocator");
        return false;
    }
    if (sys->atom.count > UINT32_MAX - 1) {
        MD_LOG_ERROR("Hydrogen bonds: too many atoms");
        return false;
    }

    const md_hbond_params_t params = (desc && desc->params) ? *desc->params : md_hbond_params_preset(MD_HBOND_PRESET_REALISTIC);
    if (!params_validate(&params)) return false;
    if (desc && desc->set_b && !desc->set_a) {
        MD_LOG_ERROR("Hydrogen bonds: set_b requires set_a");
        return false;
    }

    const size_t N = sys->atom.count;
    q->params    = params;
    q->num_atoms = (uint32_t)N;
    q->alloc     = alloc;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool result = false;

    uint8_t* role = md_alloc(temp_alloc, MAX(N, 1));
    uint8_t* cap  = md_alloc(temp_alloc, MAX(N, 1));
    if (!md_hbond_perceive_roles(role, cap, sys, desc ? desc->reference : NULL, params.roles)) goto done;

    if (desc && desc->donors && !apply_override(role, desc->donors, MD_HBOND_ROLE_DONOR, N)) goto done;
    if (desc && desc->acceptors) {
        if (!apply_override(role, desc->acceptors, MD_HBOND_ROLE_ACCEPTOR, N)) goto done;
        for (size_t i = 0; i < N; ++i) {
            if ((role[i] & MD_HBOND_ROLE_ACCEPTOR) && cap[i] == 0) cap[i] = element_capacity(md_atom_atomic_number(&sys->atom, i));
        }
    }

    // Donors: one per hydrogen bonded to a donor atom
    for (size_t i = 0; i < N; ++i) {
        if (!(role[i] & MD_HBOND_ROLE_DONOR)) continue;
        md_bond_iter_t it = md_bond_iter(&sys->bond, i);
        while (md_bond_iter_has_next(&it)) {
            const md_atom_idx_t j = md_bond_iter_atom_index(&it);
            if (md_atom_atomic_number(&sys->atom, j) == MD_Z_H && counts_as_neighbour(sys, &it)) {
                md_array_push(q->donor_d, (uint32_t)i, alloc);
                md_array_push(q->donor_h, (uint32_t)j, alloc);
            }
            md_bond_iter_next(&it);
        }
    }
    q->num_donors = md_array_size(q->donor_d);

    // Acceptors, with the atoms bonded to them for the acceptor angle (virtual sites aside)
    md_array_push(q->acc_nbr_off, 0, alloc);
    for (size_t i = 0; i < N; ++i) {
        if (!(role[i] & MD_HBOND_ROLE_ACCEPTOR)) continue;
        md_array_push(q->acceptor, (uint32_t)i, alloc);
        md_array_push(q->acceptor_cap, cap[i], alloc);
        md_bond_iter_t it = md_bond_iter(&sys->bond, i);
        while (md_bond_iter_has_next(&it)) {
            const md_atom_idx_t j = md_bond_iter_atom_index(&it);
            if (md_atom_atomic_number(&sys->atom, j) != 0) {
                md_array_push(q->acc_nbr, (uint32_t)j, alloc);
            }
            md_bond_iter_next(&it);
        }
        md_array_push(q->acc_nbr_off, (uint32_t)md_array_size(q->acc_nbr), alloc);
    }
    q->num_acceptors = md_array_size(q->acceptor);

    if (params.exclude_bonds > 0 && q->num_donors && q->num_acceptors && sys->bond.conn.offset) {
        uint8_t* is_acc = md_alloc(temp_alloc, MAX(N, 1));
        MEMSET(is_acc, 0, N);
        for (size_t k = 0; k < q->num_acceptors; ++k) is_acc[q->acceptor[k]] = 1;
        build_exclusions(q, sys, is_acc, params.exclude_bonds, temp_alloc);
    }

    if (desc && desc->set_a) {
        md_array_resize(q->sel, N, alloc);
        MEMSET(q->sel, 0, N);
        if (!mark_selection(q->sel, desc->set_a, 1, N)) goto done;
        if (desc->set_b && !mark_selection(q->sel, desc->set_b, 2, N)) goto done;
        q->flags |= MD_HBOND_FLAG_SELECTION;
        if (desc->set_b) q->flags |= MD_HBOND_FLAG_BETWEEN;
    }

    if (q->num_donors == 0) {
        // Polar atoms but not a single hydrogen on any of them: a structure without hydrogens
        bool polar = false, polar_h = false;
        for (size_t i = 0; i < N && !polar_h; ++i) {
            const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, i);
            if (z != MD_Z_N && z != MD_Z_O) continue;
            polar = true;
            int nh, nx;
            count_neighbours(&nh, &nx, sys, i);
            polar_h = nh > 0;
        }
        if (polar && !polar_h) {
            q->flags |= MD_HBOND_FLAG_NO_HYDROGENS;
            MD_LOG_INFO("Hydrogen bonds: the system has no hydrogens on its N and O atoms, so it has no donors");
        }
    }

    result = true;
done:
    md_temp_end(temp);
    if (!result) md_hbond_query_free(q);
    return result;
}

// ### EVALUATION ###

typedef struct hb_edge_t {
    uint32_t donor;     // Index of the donor in the query
    uint32_t acc;       // Index of the acceptor in the query
    uint32_t h_atom;
    uint32_t a_atom;
    float strength;
    float r_da;
    float r_ha;
    float angle;        // D-H...A, degrees
} hb_edge_t;

typedef struct hb_eval_t {
    const md_hbond_query_t* q;
    const md_unitcell_t* cell;
    const vec3_t* xyz;

    // Active donors and acceptors, in the order of the search streams
    const uint32_t* act_don;    // -> donor index
    const vec3_t*   don_vdh;    // Minimum image H - D
    const float*    don_ldh;
    const vec3_t*   don_pt;     // Search point of the donor: H (unwrapped next to D) or D
    const uint32_t* act_acc;    // -> acceptor index
    const vec3_t*   acc_pt;

    bool  search_from_h;
    float cos_min_dha, cos_max_hda, cos_min_xah;
    float w_ha, w_da;

    md_array(hb_edge_t) edges;
    md_allocator_i* alloc;
} hb_eval_t;

static inline bool is_excluded(const md_hbond_query_t* q, uint32_t donor, uint32_t a_atom) {
    if (!q->excl_off) return false;
    size_t lo = q->excl_off[donor];
    size_t hi = q->excl_off[donor + 1];
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        const uint32_t v = q->excl_atom[mid];
        if (v == a_atom) return true;
        if (v < a_atom) lo = mid + 1; else hi = mid;
    }
    return false;
}

// The strength is 1 at these distances and below, and falls to 0 at the distance gate
#define HBOND_IDEAL_HA 1.9f
#define HBOND_IDEAL_DA 2.8f

static inline float smoothstep01(float x) {
    x = CLAMP(x, 0.0f, 1.0f);
    return x * x * (3.0f - 2.0f * x);
}

// The difference vector b - a by minimum image. dist2 is the minimum image distance squared when known (from the
// spatial structure), negative otherwise; only vectors which are longer than it pay for the general minimum image.
static inline vec3_t min_image_diff(vec3_t a, vec3_t b, float dist2, float max_len, const md_unitcell_t* cell) {
    vec3_t d = vec3_sub(b, a);
    const float l2 = vec3_dot(d, d);
    const float ref2 = dist2 >= 0.0f ? dist2 * 1.0001f + 1.0e-4f : max_len * max_len;
    if (l2 > ref2) md_util_min_image_vec3(&d, 1, cell);
    return d;
}

static float edge_strength(const hb_eval_t* e, float r_ha, float r_da, float angle_dha, float cos_hda) {
    const md_hbond_params_t* p = &e->q->params;
    float s_d = 1.0f;
    if (p->max_ha > 0) {
        s_d = smoothstep01((p->max_ha - r_ha) / e->w_ha);
    } else if (p->max_da > 0) {
        s_d = smoothstep01((p->max_da - r_da) / e->w_da);
    }
    float s_a = 1.0f;
    if (p->min_dha > 0) {
        s_a = p->min_dha < 180.0f ? smoothstep01((angle_dha - p->min_dha) / (180.0f - p->min_dha)) : 1.0f;
    } else if (p->max_hda > 0) {
        const float hda = (float)RAD_TO_DEG(acosf(CLAMP(cos_hda, -1.0f, 1.0f)));
        s_a = smoothstep01((p->max_hda - hda) / p->max_hda);
    }
    return s_d * s_a;
}

static void pair_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    hb_eval_t* e = (hb_eval_t*)user_param;
    const md_hbond_query_t* q = e->q;
    const md_hbond_params_t* p = &q->params;

    for (size_t k = 0; k < num_pairs; ++k) {
        const uint32_t i = i_idx[k];
        const uint32_t j = j_idx[k];
        const uint32_t donor = e->act_don[i];
        const uint32_t acc   = e->act_acc[j];
        const uint32_t d_atom = q->donor_d[donor];
        const uint32_t h_atom = q->donor_h[donor];
        const uint32_t a_atom = q->acceptor[acc];
        if (a_atom == d_atom || a_atom == h_atom) continue;
        if (is_excluded(q, donor, a_atom)) continue;

        const vec3_t v_dh = e->don_vdh[i];
        const float  l_dh = e->don_ldh[i];
        const vec3_t v = min_image_diff(e->don_pt[i], e->acc_pt[j], ij_dist2[k], 0, e->cell);
        vec3_t v_ha, v_da;
        if (e->search_from_h) {
            v_ha = v;
            v_da = vec3_add(v_dh, v_ha);
        } else {
            v_da = v;
            v_ha = vec3_sub(v_da, v_dh);
        }
        const float r_ha = vec3_length(v_ha);
        const float r_da = vec3_length(v_da);
        if (p->max_ha > 0 && r_ha > p->max_ha) continue;
        if (p->max_da > 0 && r_da > p->max_da) continue;
        if (r_ha <= 0.0f || r_da <= 0.0f || l_dh <= 0.0f) continue;

        // D-H...A at H: between H->D and H->A
        const float cos_dha = -vec3_dot(v_dh, v_ha) / (l_dh * r_ha);
        if (p->min_dha > 0 && cos_dha > e->cos_min_dha) continue;

        // H-D...A at D: between D->H and D->A
        const float cos_hda = vec3_dot(v_dh, v_da) / (l_dh * r_da);
        if (p->max_hda > 0 && cos_hda < e->cos_max_hda) continue;

        // X-A...H at A, for every atom X bonded to A
        if (p->min_xah > 0) {
            const vec3_t a_pos = e->acc_pt[j];
            const vec3_t v_ah = vec3_mul1(v_ha, -1.0f);
            bool ok = true;
            for (uint32_t n = q->acc_nbr_off[acc]; n < q->acc_nbr_off[acc + 1]; ++n) {
                const uint32_t x = q->acc_nbr[n];
                if (x == h_atom) continue;
                const vec3_t v_ax = min_image_diff(a_pos, e->xyz[x], -1.0f, 3.0f, e->cell);
                const float l_ax = vec3_length(v_ax);
                if (l_ax <= 0.0f) continue;
                if (vec3_dot(v_ax, v_ah) / (l_ax * r_ha) > e->cos_min_xah) { ok = false; break; }
            }
            if (!ok) continue;
        }

        const float angle = (float)RAD_TO_DEG(acosf(CLAMP(cos_dha, -1.0f, 1.0f)));
        const float s = edge_strength(e, r_ha, r_da, angle, cos_hda);
        if (p->min_strength > 0 && s < p->min_strength) continue;

        const hb_edge_t edge = {
            .donor    = donor,
            .acc      = acc,
            .h_atom   = h_atom,
            .a_atom   = a_atom,
            .strength = s,
            .r_da     = r_da,
            .r_ha     = r_ha,
            .angle    = angle,
        };
        md_array_push(e->edges, edge, e->alloc);
    }
}

// Marks every point within 'radius' of a seed point, and possibly some further away: a superset, which is all the
// competition needs (more competitors never change its outcome, fewer can). A coarse grid of planes parallel to the
// faces of the cell, spaced at least 'radius' apart, puts two points within 'radius' of each other in neighbouring
// cells, across the periodic boundary too. One cell lookup per point rather than a neighbourhood search.
static void mark_near_seeds(uint8_t* out_hit, const vec3_t* pts, size_t num_pts, const uint8_t* is_seed, float radius, const md_unitcell_t* cell, md_allocator_i* arena) {
    const bool has_cell = (md_unitcell_flags(cell) & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) != 0;
    float I[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };   // Cartesian to fractional, I[col][row]
    float width[3] = { 1, 1, 1 };                           // Perpendicular width of one fractional unit
    int   pbc[3] = { 0, 0, 0 };
    if (has_cell) {
        md_unitcell_I_extract_float(I, cell);
        md_unitcell_pbc_mask_extract(pbc, cell);
        float A[3][3];
        md_unitcell_A_extract_float(A, cell);
        const vec3_t a = { A[0][0], A[0][1], A[0][2] };
        const vec3_t b = { A[1][0], A[1][1], A[1][2] };
        const vec3_t c = { A[2][0], A[2][1], A[2][2] };
        const float vol = fabsf(vec3_dot(a, vec3_cross(b, c)));
        width[0] = vol / vec3_length(vec3_cross(b, c));
        width[1] = vol / vec3_length(vec3_cross(a, c));
        width[2] = vol / vec3_length(vec3_cross(a, b));
    }

    float* frac = md_alloc(arena, sizeof(float) * 3 * num_pts);
    float fmin[3] = { FLT_MAX, FLT_MAX, FLT_MAX };
    float fmax[3] = { -FLT_MAX, -FLT_MAX, -FLT_MAX };
    for (size_t k = 0; k < num_pts; ++k) {
        const vec3_t x = pts[k];
        for (int r = 0; r < 3; ++r) {
            float f = I[0][r] * x.x + I[1][r] * x.y + I[2][r] * x.z;
            if (pbc[r]) f -= floorf(f);
            frac[3 * k + r] = f;
            fmin[r] = MIN(fmin[r], f);
            fmax[r] = MAX(fmax[r], f);
        }
    }

    // Cells per axis: whole periods for periodic axes, the extent of the points otherwise
    uint32_t n[3];
    float origin[3], scale[3];
    for (int r = 0; r < 3; ++r) {
        const float range = pbc[r] ? 1.0f : MAX(fmax[r] - fmin[r], 1.0e-6f);
        origin[r] = pbc[r] ? 0.0f : fmin[r];
        n[r] = (uint32_t)MAX(1.0f, floorf(range * width[r] / radius));
    }
    // Keep the bit grid modest; coarser cells are still a superset
    while ((uint64_t)n[0] * n[1] * n[2] > (1ull << 24)) {
        for (int r = 0; r < 3; ++r) n[r] = MAX(1u, n[r] / 2);
    }
    for (int r = 0; r < 3; ++r) {
        const float range = pbc[r] ? 1.0f : MAX(fmax[r] - fmin[r], 1.0e-6f);
        scale[r] = (float)n[r] / range;
    }

    uint32_t* cell_of = md_alloc(arena, sizeof(uint32_t) * 3 * num_pts);
    for (size_t k = 0; k < num_pts; ++k) {
        for (int r = 0; r < 3; ++r) {
            const int64_t c = (int64_t)floorf((frac[3 * k + r] - origin[r]) * scale[r]);
            cell_of[3 * k + r] = (uint32_t)CLAMP(c, 0, (int64_t)n[r] - 1);
        }
    }

    const size_t num_cells = (size_t)n[0] * n[1] * n[2];
    uint64_t* occ = md_alloc(arena, sizeof(uint64_t) * ((num_cells + 63) / 64));
    MEMSET(occ, 0, sizeof(uint64_t) * ((num_cells + 63) / 64));
    for (size_t k = 0; k < num_pts; ++k) {
        if (!is_seed[k]) continue;
        for (int dz = -1; dz <= 1; ++dz) {
            int64_t z = (int64_t)cell_of[3 * k + 2] + dz;
            if (z < 0 || z >= (int64_t)n[2]) { if (!pbc[2]) continue; z = (z + n[2]) % n[2]; }
            for (int dy = -1; dy <= 1; ++dy) {
                int64_t y = (int64_t)cell_of[3 * k + 1] + dy;
                if (y < 0 || y >= (int64_t)n[1]) { if (!pbc[1]) continue; y = (y + n[1]) % n[1]; }
                for (int dx = -1; dx <= 1; ++dx) {
                    int64_t x = (int64_t)cell_of[3 * k + 0] + dx;
                    if (x < 0 || x >= (int64_t)n[0]) { if (!pbc[0]) continue; x = (x + n[0]) % n[0]; }
                    const size_t idx = ((size_t)z * n[1] + (size_t)y) * n[0] + (size_t)x;
                    occ[idx >> 6] |= 1ull << (idx & 63);
                }
            }
        }
    }
    for (size_t k = 0; k < num_pts; ++k) {
        const size_t idx = ((size_t)cell_of[3 * k + 2] * n[1] + cell_of[3 * k + 1]) * n[0] + cell_of[3 * k + 0];
        out_hit[k] = (occ[idx >> 6] >> (idx & 63)) & 1;
    }
}

static int compare_edge_strength(const void* a, const void* b) {
    const hb_edge_t* x = (const hb_edge_t*)a;
    const hb_edge_t* y = (const hb_edge_t*)b;
    if (x->strength != y->strength) return x->strength > y->strength ? -1 : 1;
    if (x->h_atom != y->h_atom) return x->h_atom < y->h_atom ? -1 : 1;
    return (x->a_atom > y->a_atom) - (x->a_atom < y->a_atom);
}

typedef struct key_idx_t {
    uint64_t key;
    uint64_t idx;
} key_idx_t;

static int compare_key_idx(const void* a, const void* b) {
    const uint64_t x = ((const key_idx_t*)a)->key;
    const uint64_t y = ((const key_idx_t*)b)->key;
    return (x > y) - (x < y);
}

static inline bool competition_enabled(const md_hbond_params_t* p) {
    return p->h_capacity != 0 || p->acc_capacity_mode != MD_HBOND_CAPACITY_UNLIMITED;
}

static inline uint32_t acceptor_capacity(const md_hbond_query_t* q, uint32_t acc) {
    switch (q->params.acc_capacity_mode) {
    case MD_HBOND_CAPACITY_LONE_PAIRS: return q->acceptor_cap[acc] == MD_HBOND_CAPACITY_NO_LIMIT ? 0 : q->acceptor_cap[acc];
    case MD_HBOND_CAPACITY_FIXED:      return q->params.acc_capacity_fixed;
    default:                           return 0;
    }
}

void md_hbond_set_free(md_hbond_set_t* set) {
    if (!set) return;
    if (set->alloc && set->count) {
        md_free(set->alloc, set->donor,     sizeof(uint32_t) * set->count);
        md_free(set->alloc, set->hydrogen,  sizeof(uint32_t) * set->count);
        md_free(set->alloc, set->acceptor,  sizeof(uint32_t) * set->count);
        md_free(set->alloc, set->strength,  sizeof(float) * set->count);
        md_free(set->alloc, set->dist_da,   sizeof(float) * set->count);
        md_free(set->alloc, set->dist_ha,   sizeof(float) * set->count);
        md_free(set->alloc, set->angle_dha, sizeof(float) * set->count);
    }
    MEMSET(set, 0, sizeof(md_hbond_set_t));
}

bool md_hbond_query_eval(md_hbond_set_t* out, const md_hbond_query_t* q, const md_system_state_t* state, md_allocator_i* alloc) {
    ASSERT(out);
    MEMSET(out, 0, sizeof(md_hbond_set_t));
    if (!q || !q->alloc || !state || !alloc) {
        MD_LOG_ERROR("Hydrogen bonds: missing query, state or allocator");
        return false;
    }
    if (state->num_atoms != q->num_atoms || (q->num_atoms && !state->xyz)) {
        MD_LOG_ERROR("Hydrogen bonds: the state does not match the system the query was prepared for");
        return false;
    }

    const md_hbond_params_t* p = &q->params;
    out->params = *p;
    out->flags  = q->flags;
    out->alloc  = alloc;

    if (q->num_donors == 0 || q->num_acceptors == 0) {
        return true;
    }

    // Scratch in an arena of its own, sized by the system rather than by what a temp scope has room for
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    const vec3_t* xyz = state->xyz;
    const md_unitcell_t* cell = &state->unitcell;
    const bool search_from_h = p->max_ha > 0;
    const float radius = search_from_h ? p->max_ha : p->max_da;
    const size_t nd = q->num_donors;
    const size_t na = q->num_acceptors;
    const bool compete = competition_enabled(p);

    // Every donor's hydrogen next to its heavy atom, which keeps D-H whole across the periodic boundary. A covalent
    // D-H is far shorter than 2 Å; only the few longer ones are split by the boundary and pay for the minimum image.
    vec3_t* vdh = md_alloc(arena, sizeof(vec3_t) * nd);
    for (size_t k = 0; k < nd; ++k) {
        vdh[k] = min_image_diff(xyz[q->donor_d[k]], xyz[q->donor_h[k]], -1.0f, 2.0f, cell);
    }

    // Which donors and acceptors take part: those of the selection (seeds), and with competition also those which can
    // compete with them, directly or through one other bond
    uint8_t* don_active = md_alloc(arena, nd);
    uint8_t* acc_active = md_alloc(arena, na);
    size_t num_seed = 0;
    if (q->sel) {
        const uint8_t mask = (q->flags & MD_HBOND_FLAG_BETWEEN) ? 3 : 1;
        for (size_t k = 0; k < nd; ++k) {
            don_active[k] = ((q->sel[q->donor_d[k]] | q->sel[q->donor_h[k]]) & mask) ? 1 : 0;
            num_seed += don_active[k];
        }
        for (size_t k = 0; k < na; ++k) {
            acc_active[k] = (q->sel[q->acceptor[k]] & mask) ? 1 : 0;
            num_seed += acc_active[k];
        }
    } else {
        MEMSET(don_active, 1, nd);
        MEMSET(acc_active, 1, na);
        num_seed = nd + na;
    }

    if (q->sel && compete && num_seed > 0) {
        if (2 * num_seed > nd + na) {
            // Most of the system: the shell would be all of it anyway
            MEMSET(don_active, 1, nd);
            MEMSET(acc_active, 1, na);
        } else {
            // The search points of every candidate, donors first
            vec3_t* pts = md_alloc(arena, sizeof(vec3_t) * (nd + na));
            for (size_t k = 0; k < nd; ++k) {
                const vec3_t d = xyz[q->donor_d[k]];
                pts[k] = search_from_h ? vec3_add(d, vdh[k]) : d;
            }
            for (size_t k = 0; k < na; ++k) pts[nd + k] = xyz[q->acceptor[k]];

            uint8_t* is_seed = md_alloc(arena, nd + na);
            MEMCPY(is_seed, don_active, nd);
            MEMCPY(is_seed + nd, acc_active, na);

            uint8_t* hit = md_alloc(arena, nd + na);
            mark_near_seeds(hit, pts, nd + na, is_seed, 3.0f * radius, cell, arena);
            for (size_t k = 0; k < nd; ++k) don_active[k] |= hit[k];
            for (size_t k = 0; k < na; ++k) acc_active[k] |= hit[nd + k];
        }
    }

    // Gather the active donors and acceptors
    uint32_t* act_don = md_alloc(arena, sizeof(uint32_t) * nd);
    vec3_t*   don_vdh = md_alloc(arena, sizeof(vec3_t) * nd);
    float*    don_ldh = md_alloc(arena, sizeof(float) * nd);
    vec3_t*   don_pt  = md_alloc(arena, sizeof(vec3_t) * nd);
    size_t n_don = 0;
    for (size_t k = 0; k < nd; ++k) {
        if (!don_active[k]) continue;
        const vec3_t d = xyz[q->donor_d[k]];
        act_don[n_don] = (uint32_t)k;
        don_vdh[n_don] = vdh[k];
        don_ldh[n_don] = vec3_length(vdh[k]);
        don_pt[n_don]  = search_from_h ? vec3_add(d, vdh[k]) : d;
        n_don += 1;
    }
    uint32_t* act_acc = md_alloc(arena, sizeof(uint32_t) * na);
    vec3_t*   acc_pt  = md_alloc(arena, sizeof(vec3_t) * na);
    size_t n_acc = 0;
    for (size_t k = 0; k < na; ++k) {
        if (!acc_active[k]) continue;
        act_acc[n_acc] = (uint32_t)k;
        acc_pt[n_acc]  = xyz[q->acceptor[k]];
        n_acc += 1;
    }

    hb_eval_t e = {
        .q             = q,
        .cell          = cell,
        .xyz           = xyz,
        .act_don       = act_don,
        .don_vdh       = don_vdh,
        .don_ldh       = don_ldh,
        .don_pt        = don_pt,
        .act_acc       = act_acc,
        .acc_pt        = acc_pt,
        .search_from_h = search_from_h,
        .cos_min_dha   = cosf((float)DEG_TO_RAD(p->min_dha)),
        .cos_max_hda   = cosf((float)DEG_TO_RAD(p->max_hda)),
        .cos_min_xah   = cosf((float)DEG_TO_RAD(p->min_xah)),
        .w_ha          = MAX(0.1f, p->max_ha - HBOND_IDEAL_HA),
        .w_da          = MAX(0.1f, p->max_da - HBOND_IDEAL_DA),
        .alloc         = arena,
    };

    if (n_don && n_acc) {
        md_coord_stream_t acc_stream = md_coord_stream_from_aos((const float*)acc_pt, sizeof(vec3_t), NULL, n_acc);
        md_coord_stream_t don_stream = md_coord_stream_from_aos((const float*)don_pt, sizeof(vec3_t), NULL, n_don);
        md_spatial_acc_t sa = { .alloc = arena };
        md_spatial_acc_init(&sa, &(md_spatial_acc_desc_t){ .coords = &acc_stream, .cutoff = radius, .unitcell = cell });
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&sa, &don_stream, radius, pair_callback, &e, 0);
    }

    hb_edge_t* edges = e.edges;
    size_t num_edges = md_array_size(edges);

    // Competition: strongest first, kept while the hydrogen and the acceptor have room. Ties are broken by
    // (hydrogen, acceptor), never by the order the edges were found in, which depends on what else was searched.
    //
    // Most bonds compete with nothing: their hydrogen and their acceptor have room for all of their candidates. Such
    // a bond is kept whatever the order, and takes room nobody else wants, so only the contested ones are sorted.
    if (compete && num_edges) {
        uint32_t* h_deg = md_alloc(arena, sizeof(uint32_t) * nd);
        uint32_t* a_deg = md_alloc(arena, sizeof(uint32_t) * na);
        MEMSET(h_deg, 0, sizeof(uint32_t) * nd);
        MEMSET(a_deg, 0, sizeof(uint32_t) * na);
        for (size_t k = 0; k < num_edges; ++k) {
            h_deg[edges[k].donor] += 1;
            a_deg[edges[k].acc]   += 1;
        }

        hb_edge_t* contested = md_alloc(arena, sizeof(hb_edge_t) * num_edges);
        size_t num_contested = 0;
        size_t kept = 0;
        for (size_t k = 0; k < num_edges; ++k) {
            const hb_edge_t ed = edges[k];
            const uint32_t hd = h_deg[ed.donor];
            const uint32_t acap = acceptor_capacity(q, ed.acc);
            const bool h_free = (p->h_capacity == 0 || hd <= p->h_capacity) && (hd == 1 || p->bifurcation_tol <= 0);
            const bool a_free = acap == 0 || a_deg[ed.acc] <= acap;
            if (h_free && a_free) edges[kept++] = ed;
            else contested[num_contested++] = ed;
        }

        if (num_contested) {
            qsort(contested, num_contested, sizeof(hb_edge_t), compare_edge_strength);
            // The degree arrays are reused as counts. A free bond shares its hydrogen or acceptor with a contested one
            // only where that has room for all of its bonds, so leaving the free bonds out of the count changes nothing.
            uint32_t* h_count = h_deg;
            uint32_t* a_count = a_deg;
            float*    h_best  = md_alloc(arena, sizeof(float) * nd);
            for (size_t k = 0; k < num_contested; ++k) {
                h_count[contested[k].donor] = 0;
                a_count[contested[k].acc]   = 0;
            }
            for (size_t k = 0; k < num_contested; ++k) {
                const hb_edge_t ed = contested[k];
                if (p->h_capacity && h_count[ed.donor] >= p->h_capacity) continue;
                if (h_count[ed.donor] > 0 && p->bifurcation_tol > 0 && ed.strength < p->bifurcation_tol * h_best[ed.donor]) continue;
                const uint32_t acap = acceptor_capacity(q, ed.acc);
                if (acap && a_count[ed.acc] >= acap) continue;
                if (h_count[ed.donor] == 0) h_best[ed.donor] = ed.strength;
                h_count[ed.donor] += 1;
                a_count[ed.acc]   += 1;
                edges[kept++] = ed;
            }
        }
        num_edges = kept;
    }

    // What the selection reports
    if (q->sel && num_edges) {
        const bool between = (q->flags & MD_HBOND_FLAG_BETWEEN) != 0;
        size_t kept = 0;
        for (size_t k = 0; k < num_edges; ++k) {
            const hb_edge_t ed = edges[k];
            const uint8_t sd = q->sel[q->donor_d[ed.donor]] | q->sel[ed.h_atom];
            const uint8_t sa = q->sel[ed.a_atom];
            const bool keep = between ? (((sd & 1) && (sa & 2)) || ((sd & 2) && (sa & 1))) : ((sd & 1) && (sa & 1));
            if (keep) edges[kept++] = ed;
        }
        num_edges = kept;
    }

    if (num_edges) {
        key_idx_t* order = md_alloc(arena, sizeof(key_idx_t) * num_edges);
        for (size_t k = 0; k < num_edges; ++k) {
            order[k] = (key_idx_t){ ((uint64_t)edges[k].h_atom << 32) | edges[k].a_atom, k };
        }
        qsort(order, num_edges, sizeof(key_idx_t), compare_key_idx);
        out->count     = num_edges;
        out->donor     = md_alloc(alloc, sizeof(uint32_t) * num_edges);
        out->hydrogen  = md_alloc(alloc, sizeof(uint32_t) * num_edges);
        out->acceptor  = md_alloc(alloc, sizeof(uint32_t) * num_edges);
        out->strength  = md_alloc(alloc, sizeof(float) * num_edges);
        out->dist_da   = md_alloc(alloc, sizeof(float) * num_edges);
        out->dist_ha   = md_alloc(alloc, sizeof(float) * num_edges);
        out->angle_dha = md_alloc(alloc, sizeof(float) * num_edges);
        for (size_t k = 0; k < num_edges; ++k) {
            const hb_edge_t ed = edges[order[k].idx];
            out->donor[k]     = q->donor_d[ed.donor];
            out->hydrogen[k]  = ed.h_atom;
            out->acceptor[k]  = ed.a_atom;
            out->strength[k]  = ed.strength;
            out->dist_da[k]   = ed.r_da;
            out->dist_ha[k]   = ed.r_ha;
            out->angle_dha[k] = ed.angle;
        }
    }

    md_vm_arena_destroy(arena);
    return true;
}

bool md_hbond_compute(md_hbond_set_t* out, const md_hbond_desc_t* desc, const md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    ASSERT(out);
    MEMSET(out, 0, sizeof(md_hbond_set_t));
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_hbond_query_t q;
    bool result = md_hbond_query_init(&q, desc, sys, md_temp_allocator(temp)) && md_hbond_query_eval(out, &q, state, alloc);
    md_temp_end(temp);
    return result;
}
