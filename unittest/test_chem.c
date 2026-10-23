#include "utest.h"

// Chemistry perception (md_chem): bond orders, aromaticity, delocalized groups, formal charges and hydrogen counts.
//
// The reference is the PDB Chemical Component Dictionary (unittest/ccd_chem.txt: standard residues,
// cofactors and representatives of the common functional groups, from the ideal coordinates), with the aromaticity
// RDKit perceives for each (see md_chem.h). An aromatic bond must come out aromatic, any other bond with the CCD
// order, or DELOCALIZED where the CCD has one resonance form. The charges are compared summed over delocalized groups
// and aromatic rings, where their place is a matter of resonance. Without hydrogens the protonation is a choice (pH 7
// here, the CCD has neutral forms), so h - q is compared instead, which a proton leaves unchanged.

#include <md_system.h>
#include <md_chem.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_mmcif.h>
#include <md_util.h>
#include <md_xyz.h>

#include "component_util.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>
#include <core/md_os.h>
#include <core/md_bitset.h>

#include <stdio.h>
#include <string.h>

static int find_bond(const md_system_t* sys, int a, int b) {
    for (size_t k = 0; k < sys->bond.count; ++k) {
        const md_atom_pair_t p = sys->bond.pairs[k];
        if ((p.idx[0] == a && p.idx[1] == b) || (p.idx[0] == b && p.idx[1] == a)) return (int)k;
    }
    return -1;
}

static int find_root(int* p, int x) {
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

// Number of differences from the reference, the first few printed
static int compare_component(const md_system_t* sys, const int* map, const ref_comp_t* c, bool strip_h, md_allocator_i* alloc) {
    const int na = (int)md_array_size(c->atom);
    const int nb = (int)md_array_size(c->bond);
    const int n = (int)sys->atom.count;
    int errors = 0;

    uint8_t* aromatic_atom = md_alloc(alloc, n + 1);
    MEMSET(aromatic_atom, 0, n + 1);
    for (size_t b = 0; b < sys->bond.count; ++b) {
        if (sys->bond.flags[b] & MD_BOND_FLAG_AROMATIC) aromatic_atom[sys->bond.pairs[b].idx[0]] = aromatic_atom[sys->bond.pairs[b].idx[1]] = 1;
    }
    int* group = md_alloc(alloc, sizeof(int) * (n + 1));
    for (int i = 0; i < n; ++i) group[i] = i;

    for (int r = 0; r < nb; ++r) {
        const ref_bond_t* ref = &c->bond[r];
        const int a = map[ref->i], b = map[ref->j];
        if (a < 0 || b < 0) continue;
        const int k = find_bond(sys, a, b);
        if (k < 0) {
            if (errors++ < 4) printf("  %s: bond %s-%s missing\n", c->id, c->atom[ref->i].name, c->atom[ref->j].name);
            continue;
        }
        const md_bond_flags_t f = sys->bond.flags[k];
        const int order = md_bond_order(f);
        bool ok;
        if (ref->aromatic) {
            // The CCD flags some links between aromatic rings
            ok = (f & MD_BOND_FLAG_AROMATIC) || (order == 1 && aromatic_atom[a] && aromatic_atom[b]);
        } else if (f & MD_BOND_FLAG_AROMATIC) {
            ok = false;
        } else if (f & MD_BOND_FLAG_DELOCALIZED) {
            ok = ref->order == 1 || ref->order == 2;
        } else {
            ok = order == ref->order;
        }
        if (!ok && errors++ < 4) {
            printf("  %s: bond %s-%s order %d%s, ours %d%s%s\n", c->id, c->atom[ref->i].name, c->atom[ref->j].name, ref->order, ref->aromatic ? " aromatic" : "",
                   order, (f & MD_BOND_FLAG_AROMATIC) ? " aromatic" : "", (f & MD_BOND_FLAG_DELOCALIZED) ? " delocalized" : "");
        }
        if ((f & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED)) || ref->aromatic) {
            const int x = find_root(group, a), y = find_root(group, b);
            if (x != y) group[MAX(x, y)] = MIN(x, y);
        }
    }

    int* h_ref = md_alloc(alloc, sizeof(int) * (na + 1));
    MEMSET(h_ref, 0, sizeof(int) * (na + 1));
    for (int r = 0; r < nb; ++r) {
        if (c->atom[c->bond[r].i].z == 1) h_ref[c->bond[r].j] += 1;
        if (c->atom[c->bond[r].j].z == 1) h_ref[c->bond[r].i] += 1;
    }
    int* sum_ref = md_alloc(alloc, sizeof(int) * (n + 1));
    int* sum_our = md_alloc(alloc, sizeof(int) * (n + 1));
    MEMSET(sum_ref, 0, sizeof(int) * (n + 1));
    MEMSET(sum_our, 0, sizeof(int) * (n + 1));
    for (int i = 0; i < na; ++i) {
        const int a = map[i];
        if (a < 0 || c->atom[i].z == 1) continue;
        const int g = find_root(group, a);
        const int q = md_atom_formal_charge(&sys->atom, a);
        const int h = md_atom_hydrogen_count(&sys->atom, a);
        sum_ref[g] += strip_h ? h_ref[i] - c->atom[i].q : c->atom[i].q;
        sum_our[g] += strip_h ? h - q : q;
        if (!strip_h && h != h_ref[i] && errors++ < 4) printf("  %s: atom %s has %d hydrogens, ours %d\n", c->id, c->atom[i].name, h_ref[i], h);
    }
    for (int i = 0; i < na; ++i) {
        const int a = map[i];
        if (a < 0 || c->atom[i].z == 1 || find_root(group, a) != a) continue;
        if (sum_ref[a] != sum_our[a] && errors++ < 4) {
            printf("  %s: atom %s (and its group) %s %d, ours %d\n", c->id, c->atom[i].name, strip_h ? "h - q" : "charge", sum_ref[a], sum_our[a]);
        }
    }
    return errors;
}

static int check_all(const ref_comp_t* comps, bool strip_h, bool infer, int* out_checked) {
    int failed = 0;
    *out_checked = 0;
    for (size_t i = 0; i < md_array_size(comps); ++i) {
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        md_system_state_t st = { .alloc = alloc };
        int* map = md_alloc(alloc, sizeof(int) * (md_array_size(comps[i].atom) + 1));
        build_component(&sys, &st, map, &comps[i], strip_h, infer);
        md_chem_perceive(&sys, &st, strip_h ? MD_CHEM_FLAG_PROTONATE_PH7 : MD_CHEM_FLAG_NONE);
        failed += compare_component(&sys, map, &comps[i], strip_h, alloc) != 0;
        *out_checked += 1;
        md_vm_arena_destroy(alloc);
    }
    return failed;
}

struct chem {
    md_allocator_i* alloc;
    md_array(ref_comp_t) comps;
};

UTEST_F_SETUP(chem) {
    utest_fixture->alloc = md_vm_arena_create(GIGABYTES(1));
    utest_fixture->comps = read_components(STR_LIT(MD_UNITTEST_SOURCE_DIR "/ccd_chem.txt"), utest_fixture->alloc);
    ASSERT_GT((int)md_array_size(utest_fixture->comps), 100);
}

UTEST_F_TEARDOWN(chem) {
    md_vm_arena_destroy(utest_fixture->alloc);
}

// All hydrogens given: orders, aromaticity and charges follow from the valences
UTEST_F(chem, dictionary_with_hydrogens) {
    int checked = 0;
    EXPECT_EQ(0, check_all(utest_fixture->comps, false, false, &checked));
    EXPECT_GT(checked, 100);
}

// Heavy atoms only, as in most crystal structures: hydrogens and orders from the geometry, protonation at pH 7
UTEST_F(chem, dictionary_heavy_atoms) {
    int checked = 0;
    EXPECT_EQ(0, check_all(utest_fixture->comps, true, false, &checked));
}

// As above, on bonds inferred from the coordinates
UTEST_F(chem, dictionary_inferred_bonds) {
    int checked = 0;
    EXPECT_EQ(0, check_all(utest_fixture->comps, false, true, &checked));
    EXPECT_EQ(0, check_all(utest_fixture->comps, true, true, &checked));
}

typedef struct expect_t {
    const char* id;
    int charge;         // Net formal charge
    int aromatic;       // Aromatic bonds
    int delocalized;    // Delocalized bonds
} expect_t;

// Residues at pH 7 from their heavy atoms, as a crystal structure gives them. The free amino acids of the dictionary
// are zwitterions: an NH3+ and a carboxylate (two delocalized bonds) besides the side chain.
UTEST_F(chem, residues_ph7) {
    static const expect_t expect[] = {
        { "GLY",  0,  0, 2 }, { "SER",  0,  0, 2 }, { "CYS",  0,  0, 2 },
        { "ASP", -1,  0, 4 }, { "GLU", -1,  0, 4 }, { "LYS", +1,  0, 2 }, { "ARG", +1,  0, 5 },
        { "HIS",  0,  5, 2 }, { "PHE",  0,  6, 2 }, { "TYR",  0,  6, 2 }, { "TRP",  0, 10, 2 },
        // Phosphates: 1- for each, 2- for the terminal one; the purine aromatic
        { "ATP", -4, 10, 7 },
    };
    for (size_t e = 0; e < ARRAY_SIZE(expect); ++e) {
        const ref_comp_t* c = find_component(utest_fixture->comps, expect[e].id);
        ASSERT_TRUE(c != NULL);
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        md_system_state_t st = { .alloc = alloc };
        int* map = md_alloc(alloc, sizeof(int) * (md_array_size(c->atom) + 1));
        build_component(&sys, &st, map, c, true, false);
        ASSERT_TRUE(md_chem_perceive(&sys, &st, MD_CHEM_FLAG_PROTONATE_PH7));
        int q = 0, ar = 0, dl = 0;
        for (size_t i = 0; i < sys.atom.count; ++i) q += md_atom_formal_charge(&sys.atom, i);
        for (size_t b = 0; b < sys.bond.count; ++b) {
            ar += (sys.bond.flags[b] & MD_BOND_FLAG_AROMATIC) != 0;
            dl += (sys.bond.flags[b] & MD_BOND_FLAG_DELOCALIZED) != 0;
        }
        EXPECT_EQ_MSG(expect[e].charge, q, expect[e].id);
        EXPECT_EQ_MSG(expect[e].aromatic, ar, expect[e].id);
        EXPECT_EQ_MSG(expect[e].delocalized, dl, expect[e].id);
        md_vm_arena_destroy(alloc);
    }
}

// Ions without hydrogens, built here: nitrate, phosphate, sulfate, and the monatomic ones
UTEST(chem, ions) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    static const char pdb[] =
        "HETATM    1  N   NO3 A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "HETATM    2  O1  NO3 A   1       1.250   0.000   0.000  1.00  0.00           O\n"
        "HETATM    3  O2  NO3 A   1      -0.625   1.083   0.000  1.00  0.00           O\n"
        "HETATM    4  O3  NO3 A   1      -0.625  -1.083   0.000  1.00  0.00           O\n"
        "HETATM    5  P   PO4 A   2      10.000   0.000   0.000  1.00  0.00           P\n"
        "HETATM    6  O1  PO4 A   2      10.889   0.889   0.889  1.00  0.00           O\n"
        "HETATM    7  O2  PO4 A   2       9.111  -0.889   0.889  1.00  0.00           O\n"
        "HETATM    8  O3  PO4 A   2       9.111   0.889  -0.889  1.00  0.00           O\n"
        "HETATM    9  O4  PO4 A   2      10.889  -0.889  -0.889  1.00  0.00           O\n"
        "HETATM   10  S   SO4 A   3      20.000   0.000   0.000  1.00  0.00           S\n"
        "HETATM   11  O1  SO4 A   3      20.855   0.855   0.855  1.00  0.00           O\n"
        "HETATM   12  O2  SO4 A   3      19.145  -0.855   0.855  1.00  0.00           O\n"
        "HETATM   13  O3  SO4 A   3      19.145   0.855  -0.855  1.00  0.00           O\n"
        "HETATM   14  O4  SO4 A   3      20.855  -0.855  -0.855  1.00  0.00           O\n"
        "HETATM   15 NA    NA A   4      30.000   0.000   0.000  1.00  0.00          NA\n"
        "HETATM   16 CL    CL A   5      40.000   0.000   0.000  1.00  0.00          CL\n"
        "HETATM   17 MG    MG A   6      50.000   0.000   0.000  1.00  0.00          MG\n"
        "HETATM   18  O   HOH A   7      60.000   0.000   0.000  1.00  0.00           O\n";
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, str_from_cstr(pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));

    int q[8] = { 0 };
    for (size_t ci = 0; ci < sys.residue.count; ++ci) {
        const md_urange_t r = md_residue_atom_range(&sys.residue, ci);
        for (uint32_t i = r.beg; i < r.end; ++i) q[ci] += md_atom_formal_charge(&sys.atom, i);
    }
    EXPECT_EQ(-1, q[0]);    // NO3-
    EXPECT_EQ(-3, q[1]);    // PO4 3-
    EXPECT_EQ(-2, q[2]);    // SO4 2-
    EXPECT_EQ(+1, q[3]);    // Na+
    EXPECT_EQ(-1, q[4]);    // Cl-
    EXPECT_EQ(+2, q[5]);    // Mg2+
    EXPECT_EQ(0,  q[6]);    // Water
    EXPECT_EQ(2,  md_atom_hydrogen_count(&sys.atom, 17));
    EXPECT_EQ(+1, md_atom_formal_charge(&sys.atom, 0));    // The N+ of nitrate

    // Every bond of the oxoanions is delocalized, one of each a formal double bond
    for (size_t b = 0; b < sys.bond.count; ++b) {
        EXPECT_TRUE(sys.bond.flags[b] & MD_BOND_FLAG_DELOCALIZED);
        EXPECT_TRUE(sys.bond.flags[b] & MD_BOND_FLAG_ORDER_PERCEIVED);
    }
    md_vm_arena_destroy(alloc);
}

// Orders a file gives are kept, and perception run again gives what it gave
UTEST_F(chem, given_orders_and_rerun) {
    const ref_comp_t* c = find_component(utest_fixture->comps, "PHE");
    ASSERT_TRUE(c != NULL);
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    int* map = md_alloc(alloc, sizeof(int) * (md_array_size(c->atom) + 1));
    build_component(&sys, &st, map, c, false, false);
    ASSERT_TRUE(md_chem_perceive(&sys, &st, MD_CHEM_FLAG_NONE));
    md_bond_flags_t* first = md_alloc(alloc, sizeof(md_bond_flags_t) * sys.bond.count);
    MEMCPY(first, sys.bond.flags, sizeof(md_bond_flags_t) * sys.bond.count);
    ASSERT_TRUE(md_chem_perceive(&sys, &st, MD_CHEM_FLAG_NONE));
    EXPECT_EQ(0, memcmp(first, sys.bond.flags, sizeof(md_bond_flags_t) * sys.bond.count));

    // A ring bond given as single by the file (one Kekule structure): kept, and the others follow it
    int k = -1;
    for (size_t b = 0; b < sys.bond.count && k < 0; ++b) if ((sys.bond.flags[b] & MD_BOND_FLAG_AROMATIC) && md_bond_order(sys.bond.flags[b]) == 2) k = (int)b;
    ASSERT_GE(k, 0);
    for (size_t b = 0; b < sys.bond.count; ++b) sys.bond.flags[b] = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY);
    sys.bond.flags[k] = md_bond_flags_set_order(sys.bond.flags[k], MD_BOND_ORDER_SINGLE);
    ASSERT_TRUE(md_chem_perceive(&sys, &st, MD_CHEM_FLAG_NONE));
    EXPECT_EQ(MD_BOND_ORDER_SINGLE, md_bond_order(sys.bond.flags[k]));
    EXPECT_FALSE(sys.bond.flags[k] & MD_BOND_FLAG_ORDER_PERCEIVED);
    int doubles = 0, aromatic = 0;
    for (size_t b = 0; b < sys.bond.count; ++b) {
        doubles  += md_bond_order(sys.bond.flags[b]) == 2;
        aromatic += (sys.bond.flags[b] & MD_BOND_FLAG_AROMATIC) != 0;
    }
    EXPECT_EQ(4, doubles);      // Three in the ring, the carboxyl C=O
    EXPECT_EQ(5, aromatic);     // The given bond keeps its flags (none)
    md_vm_arena_destroy(alloc);
}

// Charges the file gives (the atom/formal_charge column the loaders publish) are kept: a histidine given as
// cationic is protonated on both nitrogens, and perception run again starts from the file and not from itself
UTEST_F(chem, given_charges) {
    const ref_comp_t* c = find_component(utest_fixture->comps, "HIS");
    ASSERT_TRUE(c != NULL);
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    int* map = md_alloc(alloc, sizeof(int) * (md_array_size(c->atom) + 1));
    build_component(&sys, &st, map, c, true, false);
    int nd1 = -1, ne2 = -1;
    for (size_t i = 0; i < md_array_size(c->atom); ++i) {
        if (!strcmp(c->atom[i].name, "ND1")) nd1 = map[i];
        if (!strcmp(c->atom[i].name, "NE2")) ne2 = map[i];
    }
    ASSERT_GE(nd1, 0);
    ASSERT_GE(ne2, 0);

    float* q = md_alloc(alloc, sizeof(float) * sys.atom.count);
    for (size_t i = 0; i < sys.atom.count; ++i) q[i] = 0.0f;
    q[nd1] = 1.0f;
    sys.attributes.alloc = alloc;
    ASSERT_NE(MD_ATTRIBUTE_INVALID, md_attributes_publish_atom_column(&sys.attributes, STR_LIT("atom/formal_charge"), md_unit_none(), 1, q, sys.atom.count));

    for (int pass = 0; pass < 2; ++pass) {
        ASSERT_TRUE(md_chem_perceive(&sys, &st, MD_CHEM_FLAG_NONE));
        int total = 0, aromatic = 0;
        for (size_t i = 0; i < sys.atom.count; ++i) total += md_atom_formal_charge(&sys.atom, i);
        for (size_t b = 0; b < sys.bond.count; ++b) aromatic += (sys.bond.flags[b] & MD_BOND_FLAG_AROMATIC) != 0;
        EXPECT_EQ(+1, md_atom_formal_charge(&sys.atom, nd1));
        EXPECT_EQ(1, md_atom_hydrogen_count(&sys.atom, nd1));
        EXPECT_EQ(1, md_atom_hydrogen_count(&sys.atom, ne2));
        EXPECT_EQ(+1, total);       // The amino acid is neutral without the pH flag: NH2 and COOH
        EXPECT_EQ(5, aromatic);
    }
    md_vm_arena_destroy(alloc);
}

// A protein from a crystal structure (no hydrogens) through md_util_system_infer: charged side chains and termini,
// and the cysteines of the zinc fingers as thiolates
UTEST(chem, deposited_structure) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(2));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif")));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    ASSERT_TRUE(sys.atom.formal_charge != NULL);
    ASSERT_TRUE(sys.atom.hydrogen_count != NULL);

    int thiolate = 0, thiol = 0, lys = 0, lys_charged = 0, arg = 0, arg_charged = 0;
    for (size_t ci = 0; ci < sys.residue.count; ++ci) {
        const str_t name = md_residue_name(&sys.residue, ci);
        const md_urange_t r = md_residue_atom_range(&sys.residue, ci);
        int q = 0;
        for (uint32_t i = r.beg; i < r.end; ++i) {
            q += md_atom_formal_charge(&sys.atom, i);
            if (str_eq(name, STR_LIT("CYS")) && str_eq(md_atom_name(&sys.atom, i), STR_LIT("SG"))) {
                bool to_zinc = false;
                md_bond_iter_t it = md_bond_iter(&sys.bond, i);
                while (md_bond_iter_has_next(&it)) {
                    to_zinc |= md_atom_atomic_number(&sys.atom, md_bond_iter_atom_index(&it)) == 30;
                    md_bond_iter_next(&it);
                }
                if (to_zinc) {
                    EXPECT_EQ(-1, md_atom_formal_charge(&sys.atom, i));
                    EXPECT_EQ(0, md_atom_hydrogen_count(&sys.atom, i));
                    thiolate += 1;
                } else if (md_bond_conn_count(&sys.bond, i) == 1) {
                    EXPECT_EQ(0, md_atom_formal_charge(&sys.atom, i));
                    EXPECT_EQ(1, md_atom_hydrogen_count(&sys.atom, i));
                    thiol += 1;
                }
            }
        }
        // Side chains complete in the model; the termini add their own charge
        if (str_eq(name, STR_LIT("LYS")) && r.end - r.beg >= 9) { lys += 1; lys_charged += q >= 1; }
        if (str_eq(name, STR_LIT("ARG")) && r.end - r.beg >= 11) { arg += 1; arg_charged += q >= 1; }
    }
    EXPECT_GE(thiolate, 5);
    EXPECT_GT(thiol, 0);
    EXPECT_GT(lys, 50);
    EXPECT_EQ(lys, lys_charged);
    EXPECT_GT(arg, 20);
    EXPECT_EQ(arg, arg_charged);
    md_vm_arena_destroy(alloc);
}

// A simulation with all hydrogens: the protonation is the force field's, and every residue gets its charge
UTEST(chem, simulation_with_hydrogens) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro")));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));

    int n_res = 0, unexpected = 0;
    for (size_t ci = 0; ci < sys.residue.count; ++ci) {
        const str_t name = md_residue_name(&sys.residue, ci);
        const md_urange_t r = md_residue_atom_range(&sys.residue, ci);
        int q = 0;
        bool n_term = false, c_term = false;
        for (uint32_t i = r.beg; i < r.end; ++i) {
            q += md_atom_formal_charge(&sys.atom, i);
            const str_t an = md_atom_name(&sys.atom, i);
            n_term |= str_eq(an, STR_LIT("HT1")) || str_eq(an, STR_LIT("H1"));
            c_term |= str_eq(an, STR_LIT("OT1")) || str_eq(an, STR_LIT("OXT"));
        }
        // Amyloid beta with neutral histidines (17 atoms: HSD), and the ligand
        int expected = n_term - c_term;
        if (str_eq(name, STR_LIT("LYS")) || str_eq(name, STR_LIT("ARG"))) expected += 1;
        else if (str_eq(name, STR_LIT("ASP")) || str_eq(name, STR_LIT("GLU"))) expected -= 1;
        else if (str_eq(name, STR_LIT("PFT"))) expected = -4;
        n_res += 1;
        if (q != expected && unexpected++ < 4) printf("  %.*s %d: charge %d, expected %d\n", (int)name.len, name.ptr, (int)ci, q, expected);
    }
    EXPECT_GT(n_res, 1000);
    EXPECT_EQ(0, unexpected);
    md_vm_arena_destroy(alloc);
}

typedef struct ring_expect_t {
    const char* id;
    int count;
    struct { int size; uint32_t flags; } ring[6];    // In any order
} ring_expect_t;

static int cmp_ring(const void* a, const void* b) {
    const int* x = a;
    const int* y = b;
    return x[0] != y[0] ? x[0] - y[0] : x[1] - y[1];
}

// The rings (3 to 7 atoms, without those which are sums of smaller ones) and their flags, on dictionary components and
// on the molecules of ring_chem.txt (7 membered aromatic rings, bridged and fused ring systems, aromatic rings fused to
// rings which are not aromatic, sugars)
UTEST_F(chem, ring_flags) {
    static const ring_expect_t expect_ccd[] = {
        { "PHE", 1, { { 6, MD_RING_FLAG_AROMATIC } } },
        { "HIS", 1, { { 5, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_HETERO } } },
        { "TRP", 2, { { 5, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_HETERO | MD_RING_FLAG_FUSED | MD_RING_FLAG_FUSED_AROMATIC }, { 6, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_FUSED | MD_RING_FLAG_FUSED_AROMATIC } } },
        { "PRO", 1, { { 5, MD_RING_FLAG_HETERO } } },
        { "MES", 1, { { 6, MD_RING_FLAG_HETERO } } },                       // Morpholine: an N in the ring, no sugar
        { "NAG", 1, { { 6, MD_RING_FLAG_HETERO | MD_RING_FLAG_SUGAR } } },   // Pyranose
        // The purine, and the ribose or deoxyribose
        { "A",   3, { { 5, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_HETERO | MD_RING_FLAG_FUSED | MD_RING_FLAG_FUSED_AROMATIC }, { 6, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_HETERO | MD_RING_FLAG_FUSED | MD_RING_FLAG_FUSED_AROMATIC }, { 5, MD_RING_FLAG_HETERO | MD_RING_FLAG_SUGAR } } },
        { "DA",  3, { { 5, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_HETERO | MD_RING_FLAG_FUSED | MD_RING_FLAG_FUSED_AROMATIC }, { 6, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_HETERO | MD_RING_FLAG_FUSED | MD_RING_FLAG_FUSED_AROMATIC }, { 5, MD_RING_FLAG_HETERO | MD_RING_FLAG_SUGAR } } },
    };
    static const ring_expect_t expect_mol[] = {
        { "AZULENE", 2, { { 5, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_FUSED | MD_RING_FLAG_FUSED_AROMATIC }, { 7, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_FUSED | MD_RING_FLAG_FUSED_AROMATIC } } },
        { "TROPOLO", 1, { { 7, MD_RING_FLAG_AROMATIC } } },
        { "NORBORN", 2, { { 5, MD_RING_FLAG_FUSED }, { 5, MD_RING_FLAG_FUSED } } },               // Not the 6 around both
        { "CUBANE",  6, { { 4, MD_RING_FLAG_FUSED }, { 4, MD_RING_FLAG_FUSED }, { 4, MD_RING_FLAG_FUSED }, { 4, MD_RING_FLAG_FUSED }, { 4, MD_RING_FLAG_FUSED }, { 4, MD_RING_FLAG_FUSED } } },
        { "BICYCLO", 2, { { 3, MD_RING_FLAG_FUSED }, { 3, MD_RING_FLAG_FUSED } } },               // Bicyclo[1.1.0]butane
        { "BCO",     3, { { 6, MD_RING_FLAG_FUSED }, { 6, MD_RING_FLAG_FUSED }, { 6, MD_RING_FLAG_FUSED } } },  // All three
        { "NORCARA", 2, { { 3, MD_RING_FLAG_FUSED }, { 6, MD_RING_FLAG_FUSED } } },               // Not the 7 around both
        { "GLUCOSE", 1, { { 6, MD_RING_FLAG_HETERO | MD_RING_FLAG_SUGAR } } },
        { "DEOXYRI", 1, { { 5, MD_RING_FLAG_HETERO | MD_RING_FLAG_SUGAR } } },
        { "THF",     1, { { 5, MD_RING_FLAG_HETERO } } },
        { "MORPHOL", 1, { { 6, MD_RING_FLAG_HETERO } } },
        { "INOSITO", 1, { { 6, MD_RING_FLAG_NONE } } },
        // Aromatic rings fused only to rings which are not aromatic keep a sextet of their own: not FUSED_AROMATIC,
        // while the rings fused to them are
        { "TETRALI", 2, { { 6, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_FUSED }, { 6, MD_RING_FLAG_FUSED | MD_RING_FLAG_FUSED_AROMATIC } } },
        { "BIPHENY", 3, { { 6, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_FUSED }, { 6, MD_RING_FLAG_AROMATIC | MD_RING_FLAG_FUSED }, { 4, MD_RING_FLAG_FUSED | MD_RING_FLAG_FUSED_AROMATIC } } },
    };

    md_array(ref_comp_t) mols = read_components(STR_LIT(MD_UNITTEST_SOURCE_DIR "/ring_chem.txt"), utest_fixture->alloc);
    const struct { const ref_comp_t* comps; const ring_expect_t* expect; size_t count; bool strip_h; } sets[] = {
        { utest_fixture->comps, expect_ccd, ARRAY_SIZE(expect_ccd), false },
        { mols,                 expect_mol, ARRAY_SIZE(expect_mol), true  },
    };
    for (size_t s = 0; s < ARRAY_SIZE(sets); ++s) {
        for (size_t e = 0; e < sets[s].count; ++e) {
            const ring_expect_t* ex = &sets[s].expect[e];
            const ref_comp_t* c = find_component(sets[s].comps, ex->id);
            ASSERT_TRUE_MSG(c != NULL, ex->id);
            md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
            md_system_t sys = { .alloc = alloc };
            md_system_state_t st = { .alloc = alloc };
            int* map = md_alloc(alloc, sizeof(int) * (md_array_size(c->atom) + 1));
            build_component(&sys, &st, map, c, sets[s].strip_h, false);
            ASSERT_TRUE(md_chem_perceive(&sys, &st, sets[s].strip_h ? MD_CHEM_FLAG_PROTONATE_PH7 : MD_CHEM_FLAG_NONE));

            const int n = (int)md_index_data_num_ranges(&sys.ring);
            EXPECT_EQ_MSG(ex->count, n, ex->id);
            EXPECT_EQ_MSG((size_t)n, md_array_size(sys.ring_flags), ex->id);
            if (n == ex->count) {
                int got[6][2], want[6][2];
                for (int r = 0; r < n; ++r) {
                    got[r][0]  = (int)md_index_range_size(&sys.ring, r);
                    got[r][1]  = (int)md_system_ring_flags(&sys, r);
                    want[r][0] = ex->ring[r].size;
                    want[r][1] = (int)ex->ring[r].flags;
                }
                qsort(got, n, sizeof(got[0]), cmp_ring);
                qsort(want, n, sizeof(want[0]), cmp_ring);
                for (int r = 0; r < n; ++r) {
                    EXPECT_EQ_MSG(want[r][0], got[r][0], ex->id);
                    EXPECT_EQ_MSG(want[r][1], got[r][1], ex->id);
                }
            }
            md_vm_arena_destroy(alloc);
        }
    }
}

// The atoms of sugar rings in canonical order: the ring O, the anomeric C, and on around the ring away from the O
UTEST_F(chem, sugar_ring_order) {
    static const struct { const char* id; const char* atoms; } expect[] = {
        { "NAG", "O5 C1 C2 C3 C4 C5" },                 // Pyranose
        { "A",   "O4' C1' C2' C3' C4'" },               // Ribose
        { "DA",  "O4' C1' C2' C3' C4'" },               // 2-deoxyribose
        { "ATP", "O4' C1' C2' C3' C4'" },
        { "NAD", "O4B C1B C2B C3B C4B" },               // Both riboses, in the order of the rings
        { "NAD", "O4D C1D C2D C3D C4D" },
    };
    for (int pass = 0; pass < 2; ++pass) {          // All atoms, heavy atoms only
        const bool strip_h = pass == 1;
        for (size_t e = 0; e < ARRAY_SIZE(expect); ++e) {
            const ref_comp_t* c = find_component(utest_fixture->comps, expect[e].id);
            ASSERT_TRUE_MSG(c != NULL, expect[e].id);
            md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
            md_system_t sys = { .alloc = alloc };
            md_system_state_t st = { .alloc = alloc };
            int* map = md_alloc(alloc, sizeof(int) * (md_array_size(c->atom) + 1));
            build_component(&sys, &st, map, c, strip_h, false);
            ASSERT_TRUE(md_chem_perceive(&sys, &st, strip_h ? MD_CHEM_FLAG_PROTONATE_PH7 : MD_CHEM_FLAG_NONE));

            bool found = false;
            for (size_t r = 0; r < md_index_data_num_ranges(&sys.ring) && !found; ++r) {
                if (!(md_system_ring_flags(&sys, r) & MD_RING_FLAG_SUGAR)) continue;
                char buf[64] = {0};
                int len = 0;
                const int32_t* ring = md_index_range_beg(&sys.ring, r);
                for (size_t k = 0; k < md_index_range_size(&sys.ring, r); ++k) {
                    const str_t name = md_atom_name(&sys.atom, ring[k]);
                    len += snprintf(buf + len, sizeof(buf) - len, "%s%.*s", k ? " " : "", (int)name.len, name.ptr);
                }
                found = !strcmp(buf, expect[e].atoms);
            }
            EXPECT_TRUE_MSG(found, expect[e].atoms);
            md_vm_arena_destroy(alloc);
        }
    }
}

// ### A FRAME WITH OTHER BONDS (md_chem_perceive_delta) ###
//
// Reactions, each as its reactants (the system) and its products (a frame of it, the same atoms in the same order):
// the chemistry of the frame's region must be that of the products perceived on their own, and outside of the region
// nothing differs from the reactants. Geometries from RDKit (ETKDG, MMFF), all hydrogens explicit, the molecules of
// either side apart.

static const char delta_acid_base_before[] =
    "12\n"
    "acid_base before\n"
    "C 0.6389 -0.1087 0.0049\n"
    "H 0.4282 -0.7596 -0.8471\n"
    "H 0.0000 0.7762 -0.0651\n"
    "H 0.4218 -0.6225 0.9446\n"
    "C 2.0707 0.3130 -0.0223\n"
    "O 2.4899 1.4542 -0.1081\n"
    "O 2.9151 -0.7315 0.0607\n"
    "H 3.8049 -0.3211 0.0325\n"
    "N 9.3846 -0.0042 0.2956\n"
    "H 8.8049 -0.7391 -0.1070\n"
    "H 9.0311 0.8725 -0.0851\n"
    "H 10.3138 -0.1292 -0.1035\n"
;

static const char delta_acid_base_after[] =
    "12\n"
    "acid_base after\n"
    "C 0.4909 -0.0229 0.0769\n"
    "H 0.2134 0.4875 1.0040\n"
    "H 0.0000 0.4758 -0.7642\n"
    "H 0.1278 -1.0539 0.1276\n"
    "C 1.9980 0.0248 -0.1053\n"
    "O 2.4868 1.1875 -0.1721\n"
    "O 2.5694 -1.0989 -0.1668\n"
    "H 7.9406 0.9727 0.0055\n"
    "N 8.2692 0.0000 0.0000\n"
    "H 7.5694 -0.5951 -0.4587\n"
    "H 8.4098 -0.3192 0.9657\n"
    "H 9.1570 -0.0585 -0.5124\n"
;

static const char delta_diels_alder_before[] =
    "16\n"
    "diels_alder before\n"
    "C 0.8500 -0.0113 1.0574\n"
    "C 1.9229 0.4135 0.3803\n"
    "C 3.0805 -0.4135 0.1454\n"
    "C 4.1534 0.0113 -0.5317\n"
    "C 11.9048 0.0070 0.0359\n"
    "C 10.5710 -0.0070 -0.0359\n"
    "H 0.7872 -1.0117 1.4734\n"
    "H 0.0000 0.6463 1.2104\n"
    "H 1.9319 1.4272 -0.0150\n"
    "H 3.0715 -1.4272 0.5408\n"
    "H 5.0034 -0.6463 -0.6847\n"
    "H 4.2162 1.0117 -0.9476\n"
    "H 12.4725 -0.9167 0.0885\n"
    "H 12.4553 0.9426 0.0434\n"
    "H 10.0206 -0.9426 -0.0434\n"
    "H 10.0034 0.9167 -0.0885\n"
;

static const char delta_diels_alder_after[] =
    "16\n"
    "diels_alder after\n"
    "C 3.4709 -0.6842 0.2563\n"
    "C 3.5206 0.6899 -0.3459\n"
    "C 2.4234 1.4151 -0.6089\n"
    "C 1.0302 0.9144 -0.3607\n"
    "C 0.9723 -0.5923 -0.0956\n"
    "C 2.0979 -1.0368 0.8345\n"
    "H 4.2226 -0.7617 1.0500\n"
    "H 3.7444 -1.4108 -0.5180\n"
    "H 4.5024 1.1004 -0.5682\n"
    "H 2.5261 2.4139 -1.0257\n"
    "H 0.3992 1.1488 -1.2256\n"
    "H 0.6155 1.4576 0.4968\n"
    "H 1.0530 -1.1351 -1.0462\n"
    "H 0.0000 -0.8534 0.3374\n"
    "H 2.0379 -2.1168 1.0098\n"
    "H 1.9748 -0.5491 1.8101\n"
;

static const char delta_sn2_before[] =
    "7\n"
    "sn2 before\n"
    "C 0.3667 0.0071 0.0013\n"
    "H 0.0176 -0.4288 0.9384\n"
    "H 0.0000 -0.5643 -0.8527\n"
    "H 0.0678 1.0535 -0.0749\n"
    "Br 2.3001 -0.0675 -0.0121\n"
    "O 8.2779 0.0000 0.0000\n"
    "H 7.3001 0.0000 0.0000\n"
;

static const char delta_sn2_after[] =
    "7\n"
    "sn2 after\n"
    "C 0.7541 0.0101 -0.0209\n"
    "H 0.6911 0.8594 -0.7056\n"
    "H 0.5806 0.3422 1.0056\n"
    "H 0.0000 -0.7318 -0.2938\n"
    "Br 7.6848 0.0000 0.0000\n"
    "O 2.0371 -0.5807 -0.1166\n"
    "H 2.6848 0.1009 0.1312\n"
;

static const char delta_ring_opening_before[] =
    "10\n"
    "ring_opening before\n"
    "C 2.1816 0.8427 -0.3539\n"
    "C 2.5402 -0.4360 -0.1389\n"
    "C 1.1493 -0.8595 0.2768\n"
    "C 0.7384 0.6061 0.0303\n"
    "H 2.7166 1.7292 -0.6818\n"
    "H 3.4762 -0.9798 -0.2262\n"
    "H 0.6524 -1.5852 -0.3762\n"
    "H 1.0381 -1.1928 1.3144\n"
    "H 0.3857 1.1339 0.9231\n"
    "H 0.0000 0.7415 -0.7675\n"
;

static const char delta_ring_opening_after[] =
    "10\n"
    "ring_opening after\n"
    "C 2.0621 -0.4499 0.0968\n"
    "C 3.1721 0.4499 -0.0968\n"
    "C 4.4289 0.1388 0.2404\n"
    "C 0.8053 -0.1388 -0.2404\n"
    "H 2.2661 -1.4235 0.5377\n"
    "H 2.9681 1.4235 -0.5377\n"
    "H 5.2342 0.8479 0.0764\n"
    "H 4.6892 -0.8172 0.6832\n"
    "H 0.5450 0.8172 -0.6832\n"
    "H 0.0000 -0.8479 -0.0764\n"
;

static const char delta_pyridinium_before[] =
    "13\n"
    "pyridinium before\n"
    "N 2.0511 -1.6113 0.0130\n"
    "C 0.9374 -0.8469 -0.0006\n"
    "C 0.9536 0.5388 -0.0121\n"
    "C 2.1867 1.1824 -0.0095\n"
    "C 3.3517 0.4225 0.0043\n"
    "C 3.2336 -0.9582 0.0151\n"
    "H 0.0000 -1.3955 -0.0020\n"
    "H 0.0277 1.1022 -0.0227\n"
    "H 2.2394 2.2673 -0.0182\n"
    "H 4.3279 0.8937 0.0067\n"
    "H 4.1135 -1.5950 0.0261\n"
    "Cl 10.6935 0.0000 0.0000\n"
    "H 9.3279 0.0000 0.0000\n"
;

static const char delta_pyridinium_after[] =
    "13\n"
    "pyridinium after\n"
    "N 1.0171 -0.1170 -0.0134\n"
    "C 1.7642 -1.2332 0.0055\n"
    "C 3.1498 -1.1383 0.0195\n"
    "C 3.7416 0.1212 0.0138\n"
    "C 2.9404 1.2589 -0.0058\n"
    "C 1.5593 1.1119 -0.0193\n"
    "H 1.2346 -2.1801 0.0088\n"
    "H 3.7659 -2.0350 0.0348\n"
    "H 4.8274 0.2162 0.0247\n"
    "H 3.3917 2.2489 -0.0104\n"
    "H 0.8736 1.9525 -0.0348\n"
    "Cl 9.8274 0.0000 0.0000\n"
    "H 0.0000 -0.2060 -0.0235\n"
;

static const char delta_acetanilide_before[] =
    "21\n"
    "acetanilide before\n"
    "N 5.2194 -0.3223 -0.2951\n"
    "H 5.7711 0.5296 -0.2917\n"
    "H 5.6237 -1.0462 0.2903\n"
    "C 3.8553 -0.1087 -0.0622\n"
    "C 3.0406 -1.1601 0.3659\n"
    "C 1.6592 -0.9886 0.4891\n"
    "C 1.0747 0.2274 0.1434\n"
    "C 1.8699 1.2626 -0.3424\n"
    "C 3.2507 1.0855 -0.4635\n"
    "H 3.4762 -2.1254 0.6087\n"
    "H 1.0406 -1.8096 0.8410\n"
    "H 0.0000 0.3609 0.2326\n"
    "H 1.4160 2.2035 -0.6414\n"
    "H 3.8520 1.8915 -0.8750\n"
    "C 11.4359 -0.1305 0.0521\n"
    "H 11.1823 -0.8030 -0.7717\n"
    "H 10.7711 0.7372 0.0102\n"
    "H 11.2895 -0.6349 1.0108\n"
    "C 12.8643 0.3248 -0.0768\n"
    "O 13.1545 1.5105 -0.2060\n"
    "Cl 13.9822 -1.0041 -0.0187\n"
;

static const char delta_acetanilide_after[] =
    "21\n"
    "acetanilide after\n"
    "N 5.2257 -0.0310 0.0298\n"
    "H 5.7511 0.4758 0.7303\n"
    "H 12.9126 0.0000 0.0000\n"
    "C 3.8287 0.1614 0.0535\n"
    "C 2.9399 -0.4521 -0.8308\n"
    "C 1.5655 -0.2077 -0.7385\n"
    "C 1.0684 0.6514 0.2383\n"
    "C 1.9455 1.2680 1.1258\n"
    "C 3.3185 1.0239 1.0336\n"
    "H 3.2823 -1.1292 -1.6069\n"
    "H 0.8858 -0.6928 -1.4347\n"
    "H 0.0000 0.8383 0.3062\n"
    "H 1.5622 1.9391 1.8899\n"
    "H 3.9800 1.5187 1.7394\n"
    "C 7.4336 -0.7690 -0.5555\n"
    "H 7.6225 -0.7376 0.5214\n"
    "H 7.9126 -1.6667 -0.9570\n"
    "H 7.8644 0.1160 -1.0311\n"
    "C 5.9533 -0.8201 -0.8426\n"
    "O 5.5203 -1.4867 -1.7754\n"
    "Cl 14.2783 0.0000 0.0000\n"
;

static const char delta_biphenyl_before[] =
    "24\n"
    "biphenyl before\n"
    "Br 5.7428 -0.2589 0.0426\n"
    "C 3.8584 -0.1067 0.0175\n"
    "C 3.2664 1.1557 0.0272\n"
    "H 3.8763 2.0546 0.0489\n"
    "C 1.8752 1.2664 0.0086\n"
    "H 1.4083 2.2480 0.0160\n"
    "C 1.0837 0.1175 -0.0193\n"
    "H 0.0000 0.2051 -0.0337\n"
    "C 1.6810 -1.1436 -0.0288\n"
    "H 1.0629 -2.0376 -0.0506\n"
    "C 3.0719 -1.2577 -0.0104\n"
    "H 3.5300 -2.2429 -0.0179\n"
    "Br 16.4856 -0.2589 0.0426\n"
    "C 14.6012 -0.1067 0.0175\n"
    "C 14.0092 1.1557 0.0272\n"
    "H 14.6191 2.0546 0.0489\n"
    "C 12.6180 1.2664 0.0086\n"
    "H 12.1511 2.2480 0.0160\n"
    "C 11.8265 0.1175 -0.0193\n"
    "H 10.7428 0.2051 -0.0337\n"
    "C 12.4238 -1.1436 -0.0288\n"
    "H 11.8057 -2.0376 -0.0506\n"
    "C 13.8147 -1.2577 -0.0104\n"
    "H 14.2728 -2.2429 -0.0179\n"
;

static const char delta_biphenyl_after[] =
    "24\n"
    "biphenyl after\n"
    "Br 14.2398 0.0000 0.0000\n"
    "C 5.3525 -0.0722 -0.0098\n"
    "C 5.9546 -1.3409 -0.0504\n"
    "H 5.3373 -2.2366 -0.0799\n"
    "C 7.3453 -1.4790 -0.0696\n"
    "H 7.7909 -2.4696 -0.1047\n"
    "C 8.1584 -0.3488 -0.0474\n"
    "H 9.2398 -0.4553 -0.0619\n"
    "C 7.5824 0.9184 -0.0066\n"
    "H 8.2136 1.8028 0.0140\n"
    "C 6.1916 1.0545 0.0114\n"
    "H 5.7620 2.0533 0.0549\n"
    "Br 19.2398 0.0000 0.0000\n"
    "C 3.8873 0.0722 0.0098\n"
    "C 3.1073 -0.5892 0.9731\n"
    "H 3.5848 -1.2108 1.7280\n"
    "C 1.7165 -0.4523 0.9929\n"
    "H 1.1319 -0.9698 1.7489\n"
    "C 1.0814 0.3488 0.0475\n"
    "H 0.0000 0.4553 0.0620\n"
    "C 1.8353 1.0129 -0.9167\n"
    "H 1.3431 1.6365 -1.6582\n"
    "C 3.2261 0.8757 -0.9342\n"
    "H 3.7955 1.3941 -1.7031\n"
;

static const char delta_styrene_hbr_before[] =
    "33\n"
    "styrene_hbr before\n"
    "C 3.8096 -0.2435 -0.0073\n"
    "C 3.3144 1.0393 0.2634\n"
    "C 1.9496 1.3223 0.1431\n"
    "C 1.0598 0.3262 -0.2498\n"
    "C 1.5336 -0.9535 -0.5227\n"
    "C 2.8977 -1.2350 -0.4021\n"
    "C 5.2340 -0.6069 0.1002\n"
    "C 6.2415 0.1997 0.4595\n"
    "H 3.9764 1.8426 0.5723\n"
    "H 1.5843 2.3233 0.3573\n"
    "H 0.0000 0.5473 -0.3428\n"
    "H 0.8439 -1.7352 -0.8298\n"
    "H 5.4762 -1.6406 -0.1405\n"
    "H 6.1082 1.2431 0.7184\n"
    "H 7.2562 -0.1852 0.5025\n"
    "H 3.2402 -2.2439 -0.6217\n"
    "Br 13.7752 0.0000 0.0000\n"
    "H 12.2562 0.0000 0.0000\n"
    "C 22.6499 0.0318 -0.0386\n"
    "C 21.8830 1.2017 -0.1049\n"
    "C 20.4929 1.1433 0.0088\n"
    "C 19.8573 -0.0852 0.1791\n"
    "C 20.6108 -1.2564 0.2316\n"
    "C 22.0010 -1.1992 0.1180\n"
    "C 24.1404 0.0906 -0.1946\n"
    "H 24.6235 -0.7206 0.3603\n"
    "H 24.4096 0.0057 -1.2518\n"
    "H 24.5373 1.0329 0.1975\n"
    "H 22.3647 2.1666 -0.2448\n"
    "H 19.9060 2.0567 -0.0385\n"
    "H 18.7752 -0.1302 0.2667\n"
    "H 20.1159 -2.2155 0.3581\n"
    "H 22.5754 -2.1220 0.1533\n"
;

static const char delta_styrene_hbr_after[] =
    "33\n"
    "styrene_hbr after\n"
    "C 3.4536 0.2971 -0.0770\n"
    "C 4.4353 1.1428 0.4654\n"
    "C 5.7762 0.7541 0.4994\n"
    "C 6.1570 -0.4838 -0.0119\n"
    "C 5.1985 -1.3318 -0.5606\n"
    "C 3.8570 -0.9445 -0.5949\n"
    "C 2.0121 0.7464 -0.0981\n"
    "C 1.0131 -0.2860 0.4064\n"
    "H 4.1634 2.1176 0.8655\n"
    "H 6.5240 1.4201 0.9221\n"
    "H 7.2008 -0.7851 0.0117\n"
    "H 5.4951 -2.2949 -0.9681\n"
    "H 1.8910 1.6666 0.4846\n"
    "H 0.9792 -1.1815 -0.2237\n"
    "H 0.0000 0.1315 0.4273\n"
    "H 3.1321 -1.6208 -1.0426\n"
    "Br 1.5406 1.2555 -1.9300\n"
    "H 1.2639 -0.6034 1.4245\n"
    "C 16.0893 0.0103 -0.0376\n"
    "C 15.3831 -1.1761 -0.2716\n"
    "C 13.9876 -1.1868 -0.2371\n"
    "C 13.2871 -0.0147 0.0417\n"
    "C 13.9815 1.1677 0.2904\n"
    "C 15.3770 1.1796 0.2562\n"
    "C 17.5891 0.0149 -0.0411\n"
    "H 17.9663 -0.2088 0.9615\n"
    "H 17.9760 0.9899 -0.3558\n"
    "H 17.9804 -0.7305 -0.7413\n"
    "H 15.9165 -2.0999 -0.4832\n"
    "H 13.4474 -2.1108 -0.4251\n"
    "H 12.2008 -0.0237 0.0693\n"
    "H 13.4366 2.0809 0.5140\n"
    "H 15.9057 2.1080 0.4596\n"
;


typedef struct delta_fixture_chem_t {
    md_allocator_i* alloc;
    md_system_t before;         // The system
    md_system_state_t before_state;
    md_system_t after;          // The products on their own
    md_system_state_t after_state;
    md_system_state_t frame;    // The products' coordinates, a frame of the system
    md_bond_delta_t delta;
    md_chem_delta_t chem;
    bool ok;
} delta_fixture_chem_t;

static void delta_load(delta_fixture_chem_t* f, const char* before, const char* after) {
    MEMSET(f, 0, sizeof(*f));
    f->alloc = md_vm_arena_create(GIGABYTES(1));
    f->before = (md_system_t){ .alloc = f->alloc };
    f->after  = (md_system_t){ .alloc = f->alloc };
    f->before_state = (md_system_state_t){ .alloc = f->alloc };
    f->after_state  = (md_system_state_t){ .alloc = f->alloc };
    f->frame        = (md_system_state_t){ .alloc = f->alloc };
    if (!md_xyz_system_init_from_str(&f->before, &f->before_state, str_from_cstr(before), MD_XYZ_OPTION_NONE)) return;
    if (!md_xyz_system_init_from_str(&f->after, &f->after_state, str_from_cstr(after), MD_XYZ_OPTION_NONE)) return;
    if (!md_util_system_infer(&f->before, &f->before_state, MD_UTIL_INFER_ALL)) return;
    if (!md_util_system_infer(&f->after, &f->after_state, MD_UTIL_INFER_ALL)) return;
    const size_t n = f->before.atom.count;
    if (f->after.atom.count != n) return;
    md_system_state_init(&f->frame, n);
    for (size_t i = 0; i < n; ++i) f->frame.xyz[i] = f->after_state.xyz[i];
    const md_bond_delta_desc_t desc = { .tolerance = MD_BOND_DELTA_TOLERANCE, .breakable = MD_BOND_DELTA_BREAKABLE };
    if (!md_util_bond_delta(&f->delta, &f->before, &f->frame, &desc, f->alloc)) return;
    f->ok = md_chem_perceive_delta(&f->chem, &f->before, &f->frame, &f->delta, MD_CHEM_FLAG_PROTONATE_PH7, f->alloc);
}

static int delta_region_index(const md_chem_delta_t* chem, int atom) {
    for (size_t k = 0; k < chem->num_atoms; ++k) if (chem->atom_idx[k] == atom) return (int)k;
    return -1;
}

// The flags of the bond between atoms a and b (1 based, the map numbers of the reactions) in the frame: those of the
// region, else those of the system, -1 if there is none
static int delta_bond_flags(const delta_fixture_chem_t* f, int a, int b) {
    a -= 1; b -= 1;
    const size_t count = f->before.bond.count;
    for (size_t k = 0; k < f->chem.num_bonds; ++k) {
        const uint32_t id = f->chem.bond_idx[k];
        const md_atom_pair_t p = id < count ? f->before.bond.pairs[id] : f->delta.formed_pair[id - count];
        if ((p.idx[0] == a && p.idx[1] == b) || (p.idx[0] == b && p.idx[1] == a)) return (int)f->chem.bond_flags[k];
    }
    const int i = find_bond(&f->before, a, b);
    return (i >= 0 && !md_bitset_test(f->delta.broken, (uint32_t)i)) ? (int)f->before.bond.flags[i] : -1;
}

static int delta_order(const delta_fixture_chem_t* f, int a, int b) {
    const int fl = delta_bond_flags(f, a, b);
    return fl < 0 ? -1 : md_bond_order((md_bond_flags_t)fl);
}

static int delta_charge(const delta_fixture_chem_t* f, int a) {
    const int k = delta_region_index(&f->chem, a - 1);
    return k >= 0 ? f->chem.formal_charge[k] : md_atom_formal_charge(&f->before.atom, a - 1);
}

#define CHEM_MASK (MD_ATOM_FLAG_HYBRIDIZATION_MASK | MD_ATOM_FLAG_AROMATIC)
#define BOND_MASK (MD_BOND_FLAG_ORDER_MASK | MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED)

// The region is what the products are on their own, and the rest what the reactants are
static void delta_expect_equivalent(int* utest_result, const delta_fixture_chem_t* f, const char* name) {
    char msg[256];
    const md_system_t* before = &f->before;
    const md_system_t* after = &f->after;
    const size_t n = before->atom.count;
    for (size_t i = 0; i < n; ++i) {
        const int k = delta_region_index(&f->chem, (int)i);
        const bool lone = md_bond_conn_count(&after->bond, i) == 0;
        const md_atom_flags_t want = md_atom_flags(&after->atom, i) & CHEM_MASK;
        const md_atom_flags_t have = (k >= 0 ? f->chem.atom_flags[k] : md_atom_flags(&before->atom, i)) & CHEM_MASK;
        snprintf(msg, sizeof(msg), "%s: atom %zu (%s) chemistry", name, i + 1, k >= 0 ? "region" : "outside");
        EXPECT_EQ_MSG((int)want, (int)have, msg);
        const int hc = k >= 0 ? f->chem.hydrogen_count[k] : md_atom_hydrogen_count(&before->atom, i);
        snprintf(msg, sizeof(msg), "%s: atom %zu hydrogens", name, i + 1);
        EXPECT_EQ_MSG(md_atom_hydrogen_count(&after->atom, i), hc, msg);
        // A lone atom is a monatomic ion in the region, which md_chem_perceive gives one-atom residues alone
        if (!lone) {
            snprintf(msg, sizeof(msg), "%s: atom %zu charge", name, i + 1);
            EXPECT_EQ_MSG(md_atom_formal_charge(&after->atom, i), delta_charge(f, (int)i + 1), msg);
        }
    }
    // Every covalent bond of the products, with what it is
    size_t covalent = 0;
    for (size_t b = 0; b < after->bond.count; ++b) {
        if (md_bond_kind(after->bond.flags[b]) != MD_BOND_KIND_COVALENT) continue;
        covalent += 1;
        const int x = after->bond.pairs[b].idx[0] + 1, y = after->bond.pairs[b].idx[1] + 1;
        const int fl = delta_bond_flags(f, x, y);
        snprintf(msg, sizeof(msg), "%s: bond %d-%d", name, x, y);
        EXPECT_NE_MSG(-1, fl, msg);
        // The Kekule structure of an aromatic ring is one of several, which the bond lengths choose: only that it is
        // aromatic is compared
        md_bond_flags_t want = after->bond.flags[b] & BOND_MASK;
        md_bond_flags_t have = (md_bond_flags_t)fl & BOND_MASK;
        if ((want & MD_BOND_FLAG_AROMATIC) && (have & MD_BOND_FLAG_AROMATIC)) {
            want &= ~MD_BOND_FLAG_ORDER_MASK;
            have &= ~MD_BOND_FLAG_ORDER_MASK;
        }
        EXPECT_EQ_MSG((int)want, fl < 0 ? -1 : (int)have, msg);
    }
    // and no other
    size_t frame_covalent = 0;
    for (size_t b = 0; b < before->bond.count; ++b) {
        frame_covalent += md_bond_kind(before->bond.flags[b]) == MD_BOND_KIND_COVALENT && !md_bitset_test(f->delta.broken, (uint32_t)b);
    }
    for (size_t j = 0; j < f->delta.num_formed; ++j) frame_covalent += md_bond_kind(f->delta.formed_flags[j]) == MD_BOND_KIND_COVALENT;
    snprintf(msg, sizeof(msg), "%s: covalent bonds", name);
    EXPECT_EQ_MSG(covalent, frame_covalent, msg);

    // The rings of the region are those of the products there, with their flags
    size_t touching = 0;
    for (size_t r = 0; r < md_index_data_num_ranges(&after->ring); ++r) {
        const int32_t* ring = md_index_range_beg(&after->ring, r);
        const size_t m = md_index_range_size(&after->ring, r);
        bool in = false;
        for (size_t k = 0; k < m; ++k) in |= delta_region_index(&f->chem, ring[k]) >= 0;
        if (!in) continue;
        touching += 1;
        bool found = false;
        for (size_t q = 0; q < md_index_data_num_ranges(&f->chem.ring) && !found; ++q) {
            const int32_t* rq = md_index_range_beg(&f->chem.ring, q);
            if (md_index_range_size(&f->chem.ring, q) != m) continue;
            bool same = true;
            for (size_t k = 0; k < m; ++k) same &= rq[k] == ring[k];
            if (!same) continue;
            found = true;
            const uint32_t mask = MD_RING_FLAG_AROMATIC | MD_RING_FLAG_HETERO | MD_RING_FLAG_FUSED | MD_RING_FLAG_FUSED_AROMATIC;
            snprintf(msg, sizeof(msg), "%s: ring %zu flags", name, r);
            EXPECT_EQ_MSG((int)(after->ring_flags[r] & mask), (int)(f->chem.ring_flags[q] & mask), msg);
        }
        snprintf(msg, sizeof(msg), "%s: ring %zu of the products", name, r);
        EXPECT_TRUE_MSG(found, msg);
    }
    snprintf(msg, sizeof(msg), "%s: rings", name);
    EXPECT_EQ_MSG(touching, md_index_data_num_ranges(&f->chem.ring), msg);
}

UTEST(chem, delta_acid_base) {
    delta_fixture_chem_t f;
    delta_load(&f, delta_acid_base_before, delta_acid_base_after);
    ASSERT_TRUE(f.ok);
    delta_expect_equivalent(utest_result, &f, "acid_base");
    // Acetate and ammonium: the carboxylate delocalized, the N charged
    EXPECT_TRUE(delta_bond_flags(&f, 5, 7) & MD_BOND_FLAG_DELOCALIZED);
    EXPECT_TRUE(delta_bond_flags(&f, 5, 6) & MD_BOND_FLAG_DELOCALIZED);
    EXPECT_EQ(1, delta_charge(&f, 9));
    EXPECT_EQ(-1, delta_charge(&f, 6) + delta_charge(&f, 7));
    md_vm_arena_destroy(f.alloc);
}

UTEST(chem, delta_diels_alder) {
    delta_fixture_chem_t f;
    delta_load(&f, delta_diels_alder_before, delta_diels_alder_after);
    ASSERT_TRUE(f.ok);
    delta_expect_equivalent(utest_result, &f, "diels_alder");
    // Cyclohexene: the double bond moves to the middle of the diene, and the ring closes
    EXPECT_EQ(1, delta_order(&f, 1, 2));
    EXPECT_EQ(2, delta_order(&f, 2, 3));
    EXPECT_EQ(1, delta_order(&f, 3, 4));
    EXPECT_EQ(1, delta_order(&f, 5, 6));
    EXPECT_EQ(1, delta_order(&f, 4, 5));
    EXPECT_EQ(1, (int)md_index_data_num_ranges(&f.chem.ring));
    md_vm_arena_destroy(f.alloc);
}

UTEST(chem, delta_sn2) {
    delta_fixture_chem_t f;
    delta_load(&f, delta_sn2_before, delta_sn2_after);
    ASSERT_TRUE(f.ok);
    delta_expect_equivalent(utest_result, &f, "sn2");
    // Methanol, and the bromide gone: a monatomic ion
    EXPECT_EQ(1, delta_order(&f, 1, 6));
    EXPECT_EQ(-1, delta_order(&f, 1, 5));
    EXPECT_EQ(-1, delta_charge(&f, 5));
    EXPECT_EQ(0, delta_charge(&f, 6));
    md_vm_arena_destroy(f.alloc);
}

UTEST(chem, delta_ring_opening) {
    delta_fixture_chem_t f;
    delta_load(&f, delta_ring_opening_before, delta_ring_opening_after);
    ASSERT_TRUE(f.ok);
    delta_expect_equivalent(utest_result, &f, "ring_opening");
    // Cyclobutene to butadiene: no ring, the double bonds at the ends
    EXPECT_EQ(2, delta_order(&f, 1, 4));
    EXPECT_EQ(2, delta_order(&f, 2, 3));
    EXPECT_EQ(1, delta_order(&f, 1, 2));
    EXPECT_EQ(0, (int)md_index_data_num_ranges(&f.chem.ring));
    md_vm_arena_destroy(f.alloc);
}

UTEST(chem, delta_pyridinium) {
    delta_fixture_chem_t f;
    delta_load(&f, delta_pyridinium_before, delta_pyridinium_after);
    ASSERT_TRUE(f.ok);
    delta_expect_equivalent(utest_result, &f, "pyridinium");
    // The N+ of an aromatic pyridinium, and a chloride
    EXPECT_EQ(1, delta_charge(&f, 1));
    EXPECT_EQ(-1, delta_charge(&f, 12));
    EXPECT_TRUE(delta_bond_flags(&f, 1, 2) & MD_BOND_FLAG_AROMATIC);
    ASSERT_EQ(1, (int)md_index_data_num_ranges(&f.chem.ring));
    EXPECT_TRUE(f.chem.ring_flags[0] & MD_RING_FLAG_AROMATIC);
    md_vm_arena_destroy(f.alloc);
}

UTEST(chem, delta_acetanilide) {
    delta_fixture_chem_t f;
    delta_load(&f, delta_acetanilide_before, delta_acetanilide_after);
    ASSERT_TRUE(f.ok);
    delta_expect_equivalent(utest_result, &f, "acetanilide");
    // An amide N between the carbonyl and the ring
    EXPECT_EQ(2, delta_order(&f, 19, 20));
    EXPECT_EQ(1, delta_order(&f, 1, 19));
    EXPECT_TRUE(delta_bond_flags(&f, 4, 5) & MD_BOND_FLAG_AROMATIC);
    const int k = delta_region_index(&f.chem, 0);
    ASSERT_GE(k, 0);
    EXPECT_EQ((int)MD_HYBRIDIZATION_SP2, (int)md_atom_flags_hybridization(f.chem.atom_flags[k]));
    md_vm_arena_destroy(f.alloc);
}

UTEST(chem, delta_biphenyl) {
    delta_fixture_chem_t f;
    delta_load(&f, delta_biphenyl_before, delta_biphenyl_after);
    ASSERT_TRUE(f.ok);
    delta_expect_equivalent(utest_result, &f, "biphenyl");
    // Two aromatic rings joined by a single bond, and two bromides
    EXPECT_EQ(1, delta_order(&f, 2, 14));
    EXPECT_EQ(2, (int)md_index_data_num_ranges(&f.chem.ring));
    EXPECT_EQ(-1, delta_charge(&f, 1));
    EXPECT_EQ(-1, delta_charge(&f, 13));
    md_vm_arena_destroy(f.alloc);
}

UTEST(chem, delta_styrene_hbr) {
    delta_fixture_chem_t f;
    delta_load(&f, delta_styrene_hbr_before, delta_styrene_hbr_after);
    ASSERT_TRUE(f.ok);
    delta_expect_equivalent(utest_result, &f, "styrene_hbr");
    // The addition saturates the vinyl group; the toluene beside it is no part of the region
    EXPECT_EQ(1, delta_order(&f, 7, 8));
    EXPECT_TRUE(delta_bond_flags(&f, 1, 2) & MD_BOND_FLAG_AROMATIC);
    for (int a = 19; a <= 33; ++a) EXPECT_EQ(-1, delta_region_index(&f.chem, a - 1));
    md_vm_arena_destroy(f.alloc);
}

// Nothing changed: nothing perceived
UTEST(chem, delta_empty) {
    delta_fixture_chem_t f;
    delta_load(&f, delta_acid_base_before, delta_acid_base_before);
    ASSERT_TRUE(f.ok);
    EXPECT_EQ(0, (int)f.chem.num_atoms);
    EXPECT_EQ(0, (int)f.chem.num_bonds);
    md_vm_arena_destroy(f.alloc);
}

// The DBBA molecules touching a gold slab now and then (bonds formed to the gold in a few frames): the contact of a
// carbon with a metal changes no order, aromaticity nor hybridization, all hydrogens being there. The trajectory is
// not of the test data repository: the test runs where it has been put there.
UTEST(chem, delta_metal_contact) {
    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
    if (!md_path_is_valid(path)) UTEST_SKIP("40-40-2-ddba-dyna.xmol not in the test data");
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_xyz_data_t data = {0};
    ASSERT_TRUE(md_xyz_data_parse_file(&data, path, alloc));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_xyz_system_init_from_data(&sys, &st, &data, MD_XYZ_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    const size_t n = sys.atom.count;
    const md_bond_delta_desc_t desc = { .tolerance = MD_BOND_DELTA_TOLERANCE, .breakable = MD_BOND_DELTA_BREAKABLE };
    size_t frames = 0;
    for (size_t f = 0; f < data.num_models; ++f) {
        md_system_state_t fs = { .alloc = alloc };
        md_system_state_init(&fs, n);
        for (size_t i = 0; i < n; ++i) {
            const md_xyz_coordinate_t c = data.coordinates[data.models[f].beg_coord_index + i];
            fs.xyz[i] = vec3_set(c.x, c.y, c.z);
        }
        md_bond_delta_t d = {0};
        ASSERT_TRUE(md_util_bond_delta(&d, &sys, &fs, &desc, alloc));
        if (d.num_formed == 0 && md_bitset_empty(d.broken)) continue;
        frames += 1;
        md_chem_delta_t chem = {0};
        ASSERT_TRUE(md_chem_perceive_delta(&chem, &sys, &fs, &d, MD_CHEM_FLAG_PROTONATE_PH7, alloc));
        EXPECT_GT(chem.num_atoms, 0);
        for (size_t k = 0; k < chem.num_atoms; ++k) {
            const md_atom_flags_t mask = MD_ATOM_FLAG_HYBRIDIZATION_MASK | MD_ATOM_FLAG_AROMATIC;
            EXPECT_EQ((int)(md_atom_flags(&sys.atom, chem.atom_idx[k]) & mask), (int)(chem.atom_flags[k] & mask));
            EXPECT_EQ(md_atom_formal_charge(&sys.atom, chem.atom_idx[k]), (int)chem.formal_charge[k]);
        }
        for (size_t k = 0; k < chem.num_bonds; ++k) {
            if (chem.bond_idx[k] >= sys.bond.count) continue;
            const md_bond_flags_t mask = MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED;
            EXPECT_EQ((int)(sys.bond.flags[chem.bond_idx[k]] & mask), (int)(chem.bond_flags[k] & mask));
        }
        md_chem_delta_free(&chem, alloc);
    }
    EXPECT_GT(frames, 0);
    md_vm_arena_destroy(alloc);
}
