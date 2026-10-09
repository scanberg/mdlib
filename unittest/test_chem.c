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

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>

#include <stdio.h>
#include <string.h>

typedef struct ref_atom_t {
    char  name[8];
    int   z;
    int   q;
    vec3_t x;
} ref_atom_t;

typedef struct ref_bond_t {
    int i, j, order, aromatic;
} ref_bond_t;

typedef struct ref_comp_t {
    char id[8];
    md_array(ref_atom_t) atom;
    md_array(ref_bond_t) bond;
} ref_comp_t;

static md_array(ref_comp_t) read_components(str_t path, md_allocator_i* alloc) {
    md_array(ref_comp_t) comps = 0;
    str_t text = load_textfile(path, alloc);
    str_t line;
    char buf[256];
    while (str_extract_line(&line, &text)) {
        if (line.len == 0 || line.ptr[0] == '#' || line.len >= sizeof(buf)) continue;
        MEMCPY(buf, line.ptr, line.len);
        buf[line.len] = '\0';
        if (!strncmp(buf, "COMP ", 5)) {
            ref_comp_t c = { 0 };
            sscanf(buf + 5, "%7s", c.id);
            md_array_push(comps, c, alloc);
        } else if (buf[0] == 'A' && buf[1] == ' ' && md_array_size(comps)) {
            ref_atom_t a = { 0 };
            char el[8];
            sscanf(buf + 2, "%7s %7s %d %f %f %f", a.name, el, &a.q, &a.x.x, &a.x.y, &a.x.z);
            a.z = md_util_element_lookup((str_t){ el, strlen(el) }, true);
            md_array_push(md_array_last(comps)->atom, a, alloc);
        } else if (buf[0] == 'B' && buf[1] == ' ' && md_array_size(comps)) {
            ref_bond_t b = { 0 };
            sscanf(buf + 2, "%d %d %d %d", &b.i, &b.j, &b.order, &b.aromatic);
            md_array_push(md_array_last(comps)->bond, b, alloc);
        }
    }
    return comps;
}

static const ref_comp_t* find_component(const ref_comp_t* comps, const char* id) {
    for (size_t i = 0; i < md_array_size(comps); ++i) if (!strcmp(comps[i].id, id)) return &comps[i];
    return NULL;
}

// A system of one component; map[ref atom] = system atom or -1 (stripped hydrogens). The bonds are the CCD's, or
// inferred from the coordinates.
static void build_component(md_system_t* sys, md_system_state_t* st, int* map, const ref_comp_t* c, bool strip_h, bool infer) {
    md_allocator_i* alloc = sys->alloc;
    const int na = (int)md_array_size(c->atom);
    int n = 0;
    for (int i = 0; i < na; ++i) map[i] = (strip_h && c->atom[i].z == 1) ? -1 : n++;

    sys->atom.count = n;
    md_atom_type_add(&sys->atom.type, STR_LIT(""), (str_t){0}, 0, 0, 0, 0, 0, alloc);
    md_array_resize(sys->atom.type_idx, (size_t)n, alloc);
    md_array_resize(sys->atom.flags, (size_t)n, alloc);
    MEMSET(sys->atom.flags, 0, md_array_bytes(sys->atom.flags));
    md_system_state_init(st, n);
    for (int i = 0; i < na; ++i) {
        if (map[i] < 0) continue;
        const ref_atom_t* a = &c->atom[i];
        sys->atom.type_idx[map[i]] = md_atom_type_find_or_add(&sys->atom.type, str_from_cstr(a->name), (md_atomic_number_t)a->z, md_util_element_atomic_mass(a->z), md_util_element_vdw_radius(a->z), 0, 0, alloc);
        st->xyz[map[i]] = a->x;
    }
    st->num_atoms = n;
    sys->component.count = 1;
    md_array_push(sys->component.name, make_label(str_from_cstr(c->id)), alloc);
    md_array_push(sys->component.seq_id, 1, alloc);
    md_array_push(sys->component.flags, 0, alloc);
    md_array_push(sys->component.atom_offset, 0, alloc);
    md_array_push(sys->component.atom_offset, (uint32_t)n, alloc);

    if (infer) {
        md_util_infer_covalent_bonds(&sys->bond, st, sys, alloc);
    } else {
        for (size_t b = 0; b < md_array_size(c->bond); ++b) {
            const ref_bond_t* r = &c->bond[b];
            if (map[r->i] < 0 || map[r->j] < 0) continue;
            md_atom_pair_t p = { .idx = { map[r->i], map[r->j] } };
            md_array_push(sys->bond.pairs, p, alloc);
            md_array_push(sys->bond.flags, md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY), alloc);
            sys->bond.count += 1;
        }
        md_bond_build_connectivity(&sys->bond, n, alloc);
    }
    md_util_system_infer_rings(sys);
}

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
    for (size_t ci = 0; ci < sys.component.count; ++ci) {
        const md_urange_t r = md_component_atom_range(&sys.component, ci);
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
    for (size_t ci = 0; ci < sys.component.count; ++ci) {
        const str_t name = md_component_name(&sys.component, ci);
        const md_urange_t r = md_component_atom_range(&sys.component, ci);
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
    for (size_t ci = 0; ci < sys.component.count; ++ci) {
        const str_t name = md_component_name(&sys.component, ci);
        const md_urange_t r = md_component_atom_range(&sys.component, ci);
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
