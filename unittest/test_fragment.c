#include "utest.h"

#include <md_fragment.h>
#include <md_smiles.h>
#include <md_xyz.h>
#include <md_pdb.h>
#include <md_system.h>
#include <md_infer.h>
#include <md_element.h>
#include <core/md_allocator.h>
#include <core/md_parse.h>
#include <core/md_file.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "component_util.h"

// A system built from SMILES the way a quantum chemistry file arrives: every hydrogen an atom of its own
// (right after its heavy atom), bond orders and aromaticity known, metal - non metal bonds coordinate.
// Residues, when given, as counts of heavy atoms in SMILES order.

// What md_system_infer gives the system of its molecules and parts (MD_SYSTEM_INFER_MOLECULE_BIT, _PART_BIT)
static bool test_infer_molecules(md_system_t* sys) {
    md_molecule_data_t molecules = {0};
    const bool ok = md_molecule_infer(&molecules, sys, sys->alloc);
    md_molecule_data_free(&sys->molecule, sys->alloc);
    sys->molecule = molecules;
    return ok;
}

static bool frag_is_metal_z(md_atomic_number_t z) {
    return (z == 3 || z == 4) || (z >= 11 && z <= 13) || (z >= 19 && z <= 31) || (z >= 37 && z <= 50) || (z >= 55 && z <= 84) || z >= 87;
}

static int frag_implicit_h(const md_smiles_t* g, size_t i) {
    const md_smiles_atom_t* a = g->atoms + i;
    if (a->flags & MD_SMILES_ATOM_BRACKET) return a->h_count;
    double sum = 0.0;
    for (size_t b = 0; b < g->num_bonds; ++b) {
        const md_smiles_bond_t* bd = g->bonds + b;
        if (bd->a != i && bd->b != i) continue;
        const size_t other = bd->a == i ? bd->b : bd->a;
        if (frag_is_metal_z(g->atoms[other].z)) continue;
        const bool arom = bd->order == MD_SMILES_BOND_AROMATIC ||
            (bd->order == MD_SMILES_BOND_IMPLICIT && (g->atoms[bd->a].flags & MD_SMILES_ATOM_AROMATIC) && (g->atoms[bd->b].flags & MD_SMILES_ATOM_AROMATIC));
        sum += arom ? 1.5 : (bd->order == MD_SMILES_BOND_IMPLICIT ? 1.0 : (double)bd->order);
    }
    const int need = (int)ceil(sum - 1.0e-9);
    int valences[3] = { 0, 0, 0 };
    switch (a->z) {
    case 5:  valences[0] = 3; break;
    case 6:  valences[0] = 4; break;
    case 7:  valences[0] = 3; valences[1] = 5; break;
    case 8:  valences[0] = 2; break;
    case 15: valences[0] = 3; valences[1] = 5; break;
    case 16: valences[0] = 2; valences[1] = 4; valences[2] = 6; break;
    case 9: case 17: case 35: case 53: valences[0] = 1; break;
    default: return 0;
    }
    for (int k = 0; k < 3 && valences[k]; ++k) if (valences[k] >= need) return valences[k] - need;
    return 0;
}

static md_label_t frag_label(const char* s) {
    md_label_t l = {0};
    snprintf(l.buf, sizeof(l.buf), "%s", s);
    return l;
}

static bool frag_system_from_smiles(md_system_t* sys, const char* smiles, const char* const res_names[], const int res_heavy[], int num_res) {
    md_allocator_i* alloc = sys->alloc;
    md_smiles_t g = {0};
    md_smiles_error_t err = {0};
    if (!md_smiles_parse(&g, md_str_from_cstr(smiles), alloc, &err)) {
        printf("SMILES error in '%s' at %zu: %s\n", smiles, err.offset, err.message);
        return false;
    }

    // heavy atom i, then its hydrogens
    int32_t* map = md_alloc(alloc, sizeof(int32_t) * g.num_atoms);
    int* hcount = md_alloc(alloc, sizeof(int) * g.num_atoms);
    size_t n = 0;
    for (size_t i = 0; i < g.num_atoms; ++i) {
        hcount[i] = frag_implicit_h(&g, i);
        map[i] = (int32_t)n;
        n += 1 + (size_t)hcount[i];
    }

    sys->atom.count = n;
    md_array_resize(sys->atom.type_idx, n, alloc);
    md_array_resize(sys->atom.flags, n, alloc);
    md_array_resize(sys->atom.hydrogen_count, n, alloc);
    md_array_resize(sys->atom.formal_charge, n, alloc);
    md_memset(sys->atom.flags, 0, md_array_bytes(sys->atom.flags));
    md_memset(sys->atom.hydrogen_count, 0, md_array_bytes(sys->atom.hydrogen_count));
    md_memset(sys->atom.formal_charge, 0, md_array_bytes(sys->atom.formal_charge));

    for (size_t i = 0; i < g.num_atoms; ++i) {
        const md_atomic_number_t z = g.atoms[i].z;
        const md_str_t sym = md_atomic_number_symbol(z);
        char name[8];
        snprintf(name, sizeof(name), "%.*s%zu", (int)sym.len, sym.ptr, i + 1);
        const int32_t ai = map[i];
        sys->atom.type_idx[ai] = md_particle_type_find_or_add(&sys->atom.type, md_str_from_cstr(name), z, md_atomic_number_mass(z), md_atomic_number_vdw_radius(z), 0, 0, alloc);
        if (g.atoms[i].flags & MD_SMILES_ATOM_AROMATIC) sys->atom.flags[ai] |= MD_ATOM_FLAG_AROMATIC;
        sys->atom.hydrogen_count[ai] = (uint8_t)hcount[i];
        sys->atom.formal_charge[ai] = g.atoms[i].charge;
        for (int h = 0; h < hcount[i]; ++h) {
            char hn[8];
            snprintf(hn, sizeof(hn), "H%zu%c", i + 1, 'a' + h);
            sys->atom.type_idx[ai + 1 + h] = md_particle_type_find_or_add(&sys->atom.type, md_str_from_cstr(hn), MD_Z_H, 1.008f, 1.1f, 0, 0, alloc);
            md_atom_pair_t p = { .idx = { ai, ai + 1 + h } };
            md_array_push(sys->bond.pairs, p, alloc);
            md_bond_flags_t f = md_bond_flags_set_order(MD_BOND_FLAG_NONE, 1);
            md_array_push(sys->bond.flags, md_bond_flags_set_origin(f, MD_BOND_ORIGIN_TOPOLOGY), alloc);
            sys->bond.count += 1;
        }
    }
    for (size_t b = 0; b < g.num_bonds; ++b) {
        const md_smiles_bond_t* bd = g.bonds + b;
        const bool arom = bd->order == MD_SMILES_BOND_AROMATIC ||
            (bd->order == MD_SMILES_BOND_IMPLICIT && (g.atoms[bd->a].flags & MD_SMILES_ATOM_AROMATIC) && (g.atoms[bd->b].flags & MD_SMILES_ATOM_AROMATIC));
        const int order = arom ? 1 : (bd->order == MD_SMILES_BOND_IMPLICIT ? 1 : (int)bd->order);
        md_bond_flags_t f = md_bond_flags_set_order(MD_BOND_FLAG_NONE, order);
        if (arom) f |= MD_BOND_FLAG_AROMATIC;
        if (frag_is_metal_z(g.atoms[bd->a].z) != frag_is_metal_z(g.atoms[bd->b].z)) f = md_bond_flags_set_kind(f, MD_BOND_KIND_COORDINATE);
        md_atom_pair_t p = { .idx = { map[bd->a], map[bd->b] } };
        md_array_push(sys->bond.pairs, p, alloc);
        md_array_push(sys->bond.flags, md_bond_flags_set_origin(f, MD_BOND_ORIGIN_TOPOLOGY), alloc);
        sys->bond.count += 1;
    }

    if (num_res > 0) {
        md_array_push(sys->residue.atom_offset, 0, alloc);
        size_t heavy = 0;
        for (int r = 0; r < num_res; ++r) {
            heavy += (size_t)res_heavy[r];
            const uint32_t end = heavy < g.num_atoms ? (uint32_t)map[heavy] : (uint32_t)n;
            md_array_push(sys->residue.name, frag_label(res_names[r]), alloc);
            md_array_push(sys->residue.seq_id, r + 1, alloc);
            md_array_push(sys->residue.flags, MD_RESIDUE_FLAG_NONE, alloc);
            md_array_push(sys->residue.atom_offset, end, alloc);
        }
        sys->residue.count = (size_t)num_res;
        for (int r = 0; r < num_res; ++r) md_system_set_residue_kind(sys, (size_t)r, MD_RESIDUE_KIND_AMINO_ACID);
    }

    md_bond_build_connectivity(&sys->bond, n, alloc);
    test_infer_molecules(sys);
    {
        md_ring_data_t rings = {0};
        md_ring_infer(&rings, sys, (sys)->alloc);
        md_ring_data_free(&(sys)->ring, (sys)->alloc);
        (sys)->ring = rings;
    }
    return true;
}

typedef struct frag_case_t {
    md_allocator_i* alloc;
    md_system_t sys;
} frag_case_t;

static bool frag_case(frag_case_t* t, const char* smiles) {
    t->alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    t->sys = (md_system_t){ .alloc = t->alloc };
    return frag_system_from_smiles(&t->sys, smiles, NULL, NULL, 0) && md_util_system_infer_fragments(&t->sys);
}

static void frag_case_free(frag_case_t* t) {
    md_arena_allocator_destroy(t->alloc);
}

// The first fragment of that name, MD_FRAGMENT_NONE for none
static md_fragment_idx_t frag_find(const md_system_t* sys, const char* name) {
    for (size_t f = 0; f < md_system_fragment_count(sys); ++f) {
        if (md_str_eq_cstr(md_system_fragment_name(sys, (md_fragment_idx_t)f), name)) return (md_fragment_idx_t)f;
    }
    return MD_FRAGMENT_NONE;
}

static size_t frag_count_named(const md_system_t* sys, const char* name) {
    const md_fragment_type_idx_t t = md_fragment_type_find(sys, md_str_from_cstr(name), MD_FRAGMENT_KIND_COUNT);
    return t == MD_FRAGMENT_NONE ? 0 : md_system_fragment_find_by_type(NULL, 0, sys, t);
}

static size_t frag_heavy_count(const md_system_t* sys, md_fragment_idx_t f) {
    size_t c = 0;
    for (md_fragment_atom_iter_t it = md_fragment_atom_iter(sys, f); md_fragment_atom_iter_next(&it);) {
        c += md_atom_atomic_number(&sys->atom, (size_t)it.atom) > 1 ? 1 : 0;
    }
    return c;
}

// What every result must satisfy: the fragments' atoms are a partition of the atoms with an element, each fragment
// within one molecule, and the lookups both ways agree
static void frag_check(int* utest_result, const md_system_t* sys) {
    const size_t n = sys->atom.count;
    size_t* seen = calloc(n ? n : 1, sizeof(size_t));
    for (size_t f = 0; f < md_system_fragment_count(sys); ++f) {
        const size_t m = md_system_fragment_molecule(sys, (md_fragment_idx_t)f);
        size_t count = 0;
        for (md_fragment_atom_iter_t it = md_fragment_atom_iter(sys, (md_fragment_idx_t)f); md_fragment_atom_iter_next(&it);) {
            ASSERT_GE(it.atom, 0);
            seen[it.atom]++;
            count++;
            EXPECT_EQ(md_system_fragment_find_by_atom_idx(sys, (size_t)it.atom), (md_fragment_idx_t)f);
            EXPECT_EQ((size_t)md_system_molecule_find_by_atom_idx(sys, (size_t)it.atom), m);
        }
        EXPECT_GT(count, (size_t)0);
        EXPECT_EQ(count, md_system_fragment_atom_count(sys, (md_fragment_idx_t)f));
        EXPECT_LT((size_t)md_system_fragment_kind(sys, (md_fragment_idx_t)f), (size_t)MD_FRAGMENT_KIND_COUNT);
    }
    for (size_t a = 0; a < n; ++a) EXPECT_EQ(seen[a], md_atom_atomic_number(&sys->atom, a) ? (size_t)1 : (size_t)0);
    free(seen);
}

static void frag_print(const md_system_t* sys) {
    printf("  have:");
    for (size_t f = 0; f < md_system_fragment_count(sys); ++f) {
        const md_str_t nm = md_system_fragment_name(sys, (md_fragment_idx_t)f);
        printf(" [%.*s]", (int)nm.len, nm.ptr);
    }
    printf("\n");
}

#define EXPECT_FRAG(SYS, NAME, KIND, HEAVY) do { \
    const md_fragment_idx_t i_ = frag_find((SYS), NAME); \
    EXPECT_NE(i_, MD_FRAGMENT_NONE); \
    if (i_ != MD_FRAGMENT_NONE) { \
        EXPECT_EQ((int)md_system_fragment_kind((SYS), i_), (int)(KIND)); \
        EXPECT_EQ(frag_heavy_count((SYS), i_), (size_t)(HEAVY)); \
    } else { \
        printf("  no fragment '%s'\n", NAME); \
        frag_print(SYS); \
    } \
} while (0)

UTEST(fragment, amine_linker_phenyl) {
    // The S1 of test_data/vlx/mol.h5 moves its electron from the amine to the ring through this split
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "CNC(C)Cc1ccccc1"));
    frag_check(utest_result, &t.sys);
    ASSERT_EQ(md_system_fragment_count(&t.sys), (size_t)3);
    EXPECT_FRAG(&t.sys, "amine",     MD_FRAGMENT_GROUP,       2);   // N and its methyl
    EXPECT_FRAG(&t.sys, "C3 linker", MD_FRAGMENT_LINKER,      3);
    EXPECT_FRAG(&t.sys, "phenyl",    MD_FRAGMENT_RING_SYSTEM, 6);
    frag_case_free(&t);
}

UTEST(fragment, single_fragment_molecules) {
    // DMF: the methyls go with the amide, which is then the whole molecule
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "CN(C)C=O"));
    frag_check(utest_result, &t.sys);
    ASSERT_EQ(md_system_fragment_count(&t.sys), (size_t)1);
    EXPECT_FRAG(&t.sys, "amide", MD_FRAGMENT_GROUP, 5);
    EXPECT_EQ(md_system_fragment_atom_count(&t.sys, 0), (size_t)12);
    frag_case_free(&t);

    const struct { const char* smiles; const char* name; int kind; } cases[] = {
        { "CCOC(C)=O",        "ester",       MD_FRAGMENT_GROUP },        // the ethyl and the methyl go with it
        { "c1ccc2ccccc2c1",   "naphthyl",    MD_FRAGMENT_RING_SYSTEM },  // fused: one fragment
        { "Clc1ccccc1",       "phenyl",      MD_FRAGMENT_RING_SYSTEM },  // a single halogen goes with its ring
        { "O=C1C=CC(=O)C=C1", "benzoquinone",MD_FRAGMENT_RING_SYSTEM },  // the exocyclic O are the ring's
        { "c1ccc2c(c1)ccn2",  "indolyl",     MD_FRAGMENT_RING_SYSTEM },
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        frag_case_t u = {0};
        ASSERT_TRUE(frag_case(&u, cases[i].smiles));
        frag_check(utest_result, &u.sys);
        EXPECT_EQ(md_system_fragment_count(&u.sys), (size_t)1);
        EXPECT_TRUE(md_str_eq_cstr(md_system_fragment_name(&u.sys, 0), cases[i].name));
        EXPECT_EQ((int)md_system_fragment_kind(&u.sys, 0), cases[i].kind);
        if (!md_str_eq_cstr(md_system_fragment_name(&u.sys, 0), cases[i].name)) frag_print(&u.sys);
        frag_case_free(&u);
    }
}

UTEST(fragment, an_alkyl_is_of_a_group_but_on_a_ring_of_its_own) {
    // 4-methylpyridine: the methyl is a substituent of the ring, not a part of it
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "Cc1ccncc1"));
    frag_check(utest_result, &t.sys);
    ASSERT_EQ(md_system_fragment_count(&t.sys), (size_t)2);
    EXPECT_FRAG(&t.sys, "pyridyl", MD_FRAGMENT_RING_SYSTEM, 6);
    EXPECT_FRAG(&t.sys, "methyl",  MD_FRAGMENT_ALKYL,       1);
    frag_case_free(&t);

    // anisole: the methyl of the methoxy is of the ether
    ASSERT_TRUE(frag_case(&t, "COc1ccccc1"));
    frag_check(utest_result, &t.sys);
    ASSERT_EQ(md_system_fragment_count(&t.sys), (size_t)2);
    EXPECT_FRAG(&t.sys, "ether",  MD_FRAGMENT_GROUP,       2);
    EXPECT_FRAG(&t.sys, "phenyl", MD_FRAGMENT_RING_SYSTEM, 6);
    frag_case_free(&t);
}

UTEST(fragment, donor_bridge_acceptor_stays_apart) {
    // 4-(dimethylamino)benzonitrile, the textbook twisted charge transfer molecule
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "CN(C)c1ccc(cc1)C#N"));
    frag_check(utest_result, &t.sys);
    EXPECT_EQ(md_system_fragment_count(&t.sys), (size_t)3);
    EXPECT_FRAG(&t.sys, "amine",   MD_FRAGMENT_GROUP,       3);
    EXPECT_FRAG(&t.sys, "phenyl",  MD_FRAGMENT_RING_SYSTEM, 6);
    EXPECT_FRAG(&t.sys, "nitrile", MD_FRAGMENT_GROUP,       2);
    frag_case_free(&t);

    t = (frag_case_t){0};
    ASSERT_TRUE(frag_case(&t, "c1ccccc1[N+](=O)[O-]"));
    frag_check(utest_result, &t.sys);
    EXPECT_FRAG(&t.sys, "phenyl", MD_FRAGMENT_RING_SYSTEM, 6);
    EXPECT_FRAG(&t.sys, "nitro",  MD_FRAGMENT_GROUP,       3);
    frag_case_free(&t);

    t = (frag_case_t){0};
    ASSERT_TRUE(frag_case(&t, "FC(F)(F)c1ccccc1"));
    frag_check(utest_result, &t.sys);
    EXPECT_FRAG(&t.sys, "trifluoromethyl", MD_FRAGMENT_GROUP,       4);
    EXPECT_FRAG(&t.sys, "phenyl",          MD_FRAGMENT_RING_SYSTEM, 6);
    frag_case_free(&t);
}

UTEST(fragment, a_name_is_a_type) {
    // Two phenyls of one molecule are two fragments of one type
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "O=C(c1ccccc1)c1ccccc1"));
    frag_check(utest_result, &t.sys);
    EXPECT_FRAG(&t.sys, "ketone", MD_FRAGMENT_GROUP, 2);     // benzophenone
    EXPECT_EQ(frag_count_named(&t.sys, "phenyl"), (size_t)2);
    EXPECT_EQ(md_system_fragment_type_count(&t.sys), (size_t)2);
    frag_case_free(&t);

    t = (frag_case_t){0};
    ASSERT_TRUE(frag_case(&t, "c1ccc(cc1)-c1ccccc1"));   // biphenyl is two rings, not one system
    frag_check(utest_result, &t.sys);
    EXPECT_EQ(md_system_fragment_count(&t.sys), (size_t)2);
    EXPECT_EQ(frag_count_named(&t.sys, "phenyl"), (size_t)2);
    frag_case_free(&t);
}

UTEST(fragment, molecules_alike_share_a_template) {
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "O.O.O.CC#N.O"));
    frag_check(utest_result, &t.sys);
    EXPECT_EQ(md_system_molecule_count(&t.sys), (size_t)5);
    EXPECT_EQ(t.sys.fragment.template_count, (size_t)2);
    EXPECT_EQ(md_system_fragment_count(&t.sys), (size_t)5);
    EXPECT_EQ(frag_count_named(&t.sys, "water"), (size_t)4);
    EXPECT_FRAG(&t.sys, "nitrile", MD_FRAGMENT_GROUP, 3);   // acetonitrile: its methyl goes with it
    // every fragment of a type, ascending
    const md_fragment_type_idx_t w = md_fragment_type_find(&t.sys, MD_STR_LIT("water"), MD_FRAGMENT_GROUP);
    md_fragment_idx_t found[8];
    ASSERT_EQ(md_system_fragment_find_by_type(found, 8, &t.sys, w), (size_t)4);
    for (size_t i = 1; i < 4; ++i) EXPECT_LT(found[i - 1], found[i]);
    frag_case_free(&t);
}

UTEST(fragment, molecules_alike_share_a_template_whatever_order_their_bonds_came_in) {
    // A spatial search finds bonds in cell order, so the bonds of molecules alike come in different orders. The
    // connectivity is ordered by atom, so the molecules are walked alike and share their template.
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "CC(C)C(O)CN.CC(C)C(O)CN.CC(C)C(O)CN"));
    md_bond_data_t* bond = &t.sys.bond;
    uint32_t seed = 12345;
    for (size_t i = bond->count; i > 1; --i) {   // shuffled: each molecule's bonds in an order of their own
        seed = seed * 1664525u + 1013904223u;
        const size_t j = (seed >> 8) % i;
        const md_atom_pair_t p = bond->pairs[i - 1];
        bond->pairs[i - 1] = bond->pairs[j];
        bond->pairs[j] = p;
        const md_bond_flags_t f = bond->flags[i - 1];
        bond->flags[i - 1] = bond->flags[j];
        bond->flags[j] = f;
    }
    md_bond_build_connectivity(bond, t.sys.atom.count, t.alloc);
    ASSERT_TRUE(test_infer_molecules(&t.sys));
    ASSERT_TRUE(md_util_system_infer_fragments(&t.sys));
    frag_check(utest_result, &t.sys);
    EXPECT_EQ(md_system_molecule_count(&t.sys), (size_t)3);
    EXPECT_EQ(t.sys.fragment.template_count, (size_t)1);
    frag_case_free(&t);
}

UTEST(fragment, metal_complex) {
    // Coordination bonds join nothing: the metal and each ligand are molecules of their own
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "[Zn](n1ccccc1)n1ccccc1"));
    frag_check(utest_result, &t.sys);
    EXPECT_FRAG(&t.sys, "Zn", MD_FRAGMENT_ATOM, 1);
    EXPECT_EQ(frag_count_named(&t.sys, "pyridyl"), (size_t)2);
    EXPECT_EQ(md_system_fragment_find_by_ring(&t.sys, 0) != MD_FRAGMENT_NONE, true);
    EXPECT_TRUE(md_str_eq_cstr(md_system_fragment_name(&t.sys, md_system_fragment_find_by_ring(&t.sys, 0)), "pyridyl"));
    frag_case_free(&t);
}

UTEST(fragment, the_backbone_of_a_peptide) {
    // Gly-Gly, N, CA and C of each with the backbone role: the N-terminal amine, an alpha carbon per residue, the
    // peptide group between them (of both residues), and the C-terminal acid; all of the backbone
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    const char* names[] = { "GLY", "GLY" };
    const int heavy[] = { 4, 5 };
    ASSERT_TRUE(frag_system_from_smiles(&sys, "NCC(=O)NCC(=O)O", names, heavy, 2));
    for (size_t r = 0; r < 2; ++r) sys.residue.flags[r] |= MD_RESIDUE_FLAG_RESOLVED;
    // heavy atoms in SMILES order, each followed by its hydrogens: the backbone are the N, CA and C of each
    const md_atomic_number_t want[] = { MD_Z_N, MD_Z_C, MD_Z_C, MD_Z_N, MD_Z_C, MD_Z_C };
    md_atom_idx_t heavy_atom[9];
    size_t k = 0, heavy_seen = 0;
    for (size_t a = 0; a < sys.atom.count; ++a) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys.atom, a);
        if (z == MD_Z_H) continue;
        const size_t h = heavy_seen++;
        if (h < 9) heavy_atom[h] = (md_atom_idx_t)a;
        if (h == 0 || h == 1 || h == 2 || h == 4 || h == 5 || h == 6) {
            EXPECT_EQ(z, want[k]);
            sys.atom.flags[a] |= MD_ATOM_FLAG_BACKBONE;
            k++;
        }
    }
    ASSERT_EQ(heavy_seen, (size_t)9);
    ASSERT_TRUE(md_util_system_infer_fragments(&sys));
    frag_check(utest_result, &sys);
    EXPECT_EQ(md_system_fragment_count(&sys), (size_t)5);
    EXPECT_EQ(frag_count_named(&sys, "alpha carbon"), (size_t)2);
    EXPECT_EQ(frag_count_named(&sys, "peptide"), (size_t)1);
    for (md_fragment_idx_t f = 0; f < md_system_fragment_count(&sys); ++f) {
        EXPECT_EQ((int)md_system_fragment_kind(&sys, f), (int)MD_FRAGMENT_BACKBONE);
    }
    // the peptide group: C, O of the first residue, N of the second, and of both
    const md_fragment_idx_t pep = md_system_fragment_find_by_atom_idx(&sys, (size_t)heavy_atom[2]);
    EXPECT_TRUE(md_str_eq_cstr(md_system_fragment_name(&sys, pep), "peptide"));
    EXPECT_EQ(md_system_fragment_find_by_atom_idx(&sys, (size_t)heavy_atom[3]), pep);
    EXPECT_EQ(md_system_fragment_find_by_atom_idx(&sys, (size_t)heavy_atom[4]), pep);
    EXPECT_EQ(frag_heavy_count(&sys, pep), (size_t)3);
    md_residue_idx_t res[4];
    EXPECT_EQ(md_system_fragment_residues(res, 4, &sys, pep), (size_t)2);
    // an alpha carbon is CA alone (with its hydrogens), the termini groups of their own
    const md_fragment_idx_t ca = md_system_fragment_find_by_atom_idx(&sys, (size_t)heavy_atom[1]);
    EXPECT_TRUE(md_str_eq_cstr(md_system_fragment_name(&sys, ca), "alpha carbon"));
    EXPECT_EQ(frag_heavy_count(&sys, ca), (size_t)1);
    EXPECT_NE(md_system_fragment_find_by_atom_idx(&sys, (size_t)heavy_atom[0]), ca);
    const md_fragment_idx_t cterm = md_system_fragment_find_by_atom_idx(&sys, (size_t)heavy_atom[6]);
    EXPECT_EQ(frag_heavy_count(&sys, cterm), (size_t)3);
    md_arena_allocator_destroy(alloc);
}

UTEST(fragment, the_alpha_carbon_is_of_the_backbone_alone) {
    // Ala-Ala: the alpha carbon is CA alone, of the backbone; the CB of each a methyl of its own
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    const char* names[] = { "ALA", "ALA" };
    const int heavy[] = { 5, 6 };
    // N CA(CB) C(=O) | N CA(CB) C(=O) O
    ASSERT_TRUE(frag_system_from_smiles(&sys, "NC(C)C(=O)NC(C)C(=O)O", names, heavy, 2));
    for (size_t r = 0; r < 2; ++r) sys.residue.flags[r] |= MD_RESIDUE_FLAG_RESOLVED;
    size_t heavy_seen = 0;
    for (size_t a = 0; a < sys.atom.count; ++a) {
        if (md_atom_atomic_number(&sys.atom, a) == MD_Z_H) continue;
        const size_t h = heavy_seen++;
        if (h == 0 || h == 1 || h == 3 || h == 5 || h == 6 || h == 8) sys.atom.flags[a] |= MD_ATOM_FLAG_BACKBONE;
        if (h == 2 || h == 7) sys.atom.flags[a] |= MD_ATOM_FLAG_SIDE_CHAIN;
    }
    ASSERT_TRUE(md_util_system_infer_fragments(&sys));
    frag_check(utest_result, &sys);
    EXPECT_EQ(md_system_fragment_count(&sys), (size_t)7);
    EXPECT_EQ(frag_count_named(&sys, "alpha carbon"), (size_t)2);
    EXPECT_EQ(frag_count_named(&sys, "methyl"), (size_t)2);
    for (md_fragment_idx_t f = 0; f < md_system_fragment_count(&sys); ++f) {
        if (md_str_eq_cstr(md_system_fragment_name(&sys, f), "alpha carbon")) EXPECT_EQ(frag_heavy_count(&sys, f), (size_t)1);
        if (md_str_eq_cstr(md_system_fragment_name(&sys, f), "methyl")) {
            EXPECT_EQ((int)md_system_fragment_kind(&sys, f), (int)MD_FRAGMENT_ALKYL);
        }
    }
    md_arena_allocator_destroy(alloc);
}

UTEST(fragment, a_group_across_residues_is_of_both) {
    // Propyl ethyl ester cut into two residues through the ester: one fragment, of both residues, and in the list of
    // each residue's fragments
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    const char* names[] = { "PRO", "ETH" };
    const int heavy[] = { 4, 3 };
    ASSERT_TRUE(frag_system_from_smiles(&sys, "CCC(=O)OCC", names, heavy, 2));
    for (size_t r = 0; r < 2; ++r) md_system_set_residue_kind(&sys, r, MD_RESIDUE_KIND_OTHER);
    ASSERT_TRUE(md_util_system_infer_fragments(&sys));
    frag_check(utest_result, &sys);
    const md_fragment_idx_t ester = frag_find(&sys, "ester");
    ASSERT_NE(ester, MD_FRAGMENT_NONE);
    md_residue_idx_t res[4];
    ASSERT_EQ(md_system_fragment_residues(res, 4, &sys, ester), (size_t)2);
    EXPECT_EQ(res[0], 0);
    EXPECT_EQ(res[1], 1);
    md_fragment_idx_t frags[8];
    for (size_t r = 0; r < 2; ++r) {
        const size_t n = md_system_residue_fragments(frags, 8, &sys, r);
        bool has = false;
        for (size_t i = 0; i < n; ++i) has |= frags[i] == ester;
        EXPECT_TRUE(has);
    }
    md_arena_allocator_destroy(alloc);
}

UTEST(fragment, inferring_again_keeps_the_colours) {
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "CN(C)c1ccc(cc1)C#N"));
    const md_fragment_type_idx_t p = md_fragment_type_find(&t.sys, MD_STR_LIT("phenyl"), MD_FRAGMENT_RING_SYSTEM);
    ASSERT_NE(p, MD_FRAGMENT_NONE);
    t.sys.fragment.type.color[p] = 0xFF123456u;
    const uint64_t v = t.sys.topology_version;
    ASSERT_TRUE(md_util_system_infer_fragments(&t.sys));
    EXPECT_NE(t.sys.topology_version, v);
    const md_fragment_type_idx_t q = md_fragment_type_find(&t.sys, MD_STR_LIT("phenyl"), MD_FRAGMENT_RING_SYSTEM);
    ASSERT_NE(q, MD_FRAGMENT_NONE);
    EXPECT_EQ(md_fragment_type_color(&t.sys, q), 0xFF123456u);
    frag_case_free(&t);
}

UTEST(fragment, no_molecules_no_fragments) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    ASSERT_TRUE(md_util_system_infer_fragments(&sys));
    EXPECT_EQ(md_system_fragment_count(&sys), (size_t)0);
    EXPECT_EQ(md_system_fragment_find_by_atom_idx(&sys, 0), MD_FRAGMENT_NONE);
    md_arena_allocator_destroy(alloc);
}

// ---------------------------------------------------------------------------
// The reference molecules (unittest/frag_chem.txt): each read as a user's xyz file is, elements and coordinates
// only, and its fragments compared with those the file lists
// ---------------------------------------------------------------------------

#define FRAG_REF_LINE 96

typedef struct frag_ref_t {
    char id[8];
    md_array(char) lines;       // FRAG_REF_LINE each: "heavy kind name"
} frag_ref_t;

static const char* frag_kind_token(md_fragment_kind_t k) {
    switch (k) {
    case MD_FRAGMENT_BACKBONE:    return "backbone";
    case MD_FRAGMENT_LIPID_TAIL:  return "tail";
    case MD_FRAGMENT_RING_SYSTEM: return "ring";
    case MD_FRAGMENT_GROUP:       return "group";
    case MD_FRAGMENT_LINKER:      return "linker";
    case MD_FRAGMENT_ALKYL:       return "alkyl";
    case MD_FRAGMENT_ATOM:        return "atom";
    default:                      return "?";
    }
}

static int frag_cmp_line(const void* a, const void* b) {
    return strcmp((const char*)a, (const char*)b);
}

static md_array(frag_ref_t) frag_read_refs(md_str_t path, md_allocator_i* alloc) {
    md_array(frag_ref_t) refs = 0;
    md_str_t text = md_load_textfile(path, alloc);
    md_str_t line;
    while (md_parse_next_line(&line, &text)) {
        if (line.len >= 5 && (!strncmp(line.ptr, "COMP ", 5) || !strncmp(line.ptr, "FRAG ", 5))) {
            frag_ref_t r = { 0 };
            char buf[64] = { 0 };
            memcpy(buf, line.ptr + 5, MD_MIN(line.len - 5, sizeof(buf) - 1));
            sscanf(buf, "%7s", r.id);
            md_array_push(refs, r, alloc);
        } else if (line.len > 2 && line.ptr[0] == 'F' && line.ptr[1] == ' ' && md_array_size(refs)) {
            frag_ref_t* r = md_array_last(refs);
            const size_t base = md_array_size(r->lines);
            md_array_resize(r->lines, base + FRAG_REF_LINE, alloc);
            snprintf(r->lines + base, FRAG_REF_LINE, "%.*s", (int)(line.len - 2), line.ptr + 2);
        }
    }
    return refs;
}

// The fragments of a component read from its elements and coordinates, one line each, sorted
static size_t frag_component_lines(md_array(char)* out, const ref_comp_t* c, md_allocator_i* alloc) {
    md_array(char) xyz = 0;
    char buf[128];
    int len = snprintf(buf, sizeof(buf), "%zu\n%s\n", md_array_size(c->atom), c->id);
    md_array_push_array(xyz, buf, (size_t)len, alloc);
    for (size_t i = 0; i < md_array_size(c->atom); ++i) {
        const ref_atom_t* a = &c->atom[i];
        const md_str_t sym = md_atomic_number_symbol((md_atomic_number_t)a->z);
        len = snprintf(buf, sizeof(buf), "%.*s %.4f %.4f %.4f\n", (int)sym.len, sym.ptr, a->x.x, a->x.y, a->x.z);
        md_array_push_array(xyz, buf, (size_t)len, alloc);
    }
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    if (!md_xyz_system_init_from_str(&sys, &st, (md_str_t){ xyz, md_array_size(xyz) }, 0)) return 0;
    if (!md_system_infer(&sys, &st, MD_SYSTEM_INFER_ALL, NULL)) return 0;
    const size_t n = md_system_fragment_count(&sys);
    md_array_resize(*out, n * FRAG_REF_LINE, alloc);
    for (md_fragment_idx_t f = 0; f < n; ++f) {
        size_t heavy = 0;
        for (md_fragment_atom_iter_t it = md_fragment_atom_iter(&sys, f); md_fragment_atom_iter_next(&it);) {
            heavy += md_atom_atomic_number(&sys.atom, (size_t)it.atom) > 1;
        }
        const md_str_t name = md_system_fragment_name(&sys, f);
        snprintf(*out + (size_t)f * FRAG_REF_LINE, FRAG_REF_LINE, "%zu %s %.*s", heavy, frag_kind_token(md_system_fragment_kind(&sys, f)), (int)name.len, name.ptr);
    }
    qsort(*out, n, FRAG_REF_LINE, frag_cmp_line);
    return n;
}

UTEST(fragment, reference_molecules) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(4));
    const md_str_t frag_path = MD_STR_LIT(MD_UNITTEST_SOURCE_DIR "/frag_chem.txt");
    md_array(ref_comp_t) own = read_components(frag_path, alloc);
    md_array(ref_comp_t) ccd = read_components(MD_STR_LIT(MD_UNITTEST_SOURCE_DIR "/ccd_chem.txt"), alloc);
    md_array(frag_ref_t) refs = frag_read_refs(frag_path, alloc);
    ASSERT_GT(md_array_size(refs), (size_t)200);

    size_t failed = 0;
    for (size_t r = 0; r < md_array_size(refs); ++r) {
        const frag_ref_t* ref = &refs[r];
        const ref_comp_t* c = find_component(own, ref->id);
        if (!c) c = find_component(ccd, ref->id);
        EXPECT_TRUE(c != NULL);
        if (!c) continue;

        const size_t num_want = md_array_size(ref->lines) / FRAG_REF_LINE;
        char* want = md_alloc(alloc, (num_want ? num_want : 1) * FRAG_REF_LINE);
        if (num_want) memcpy(want, ref->lines, num_want * FRAG_REF_LINE);
        qsort(want, num_want, FRAG_REF_LINE, frag_cmp_line);

        md_allocator_i* malloc_ = md_vm_arena_create(MD_GIGABYTES(1));
        md_array(char) got = 0;
        const size_t num_got = frag_component_lines(&got, c, malloc_);
        bool same = num_got == num_want;
        for (size_t i = 0; same && i < num_got; ++i) same = strcmp(got + i * FRAG_REF_LINE, want + i * FRAG_REF_LINE) == 0;
        if (!same) {
            failed++;
            printf("  %s: want", ref->id);
            for (size_t i = 0; i < num_want; ++i) printf(" [%s]", want + i * FRAG_REF_LINE);
            printf("\n  %*s  got ", (int)strlen(ref->id), "");
            for (size_t i = 0; i < num_got; ++i) printf(" [%s]", got + i * FRAG_REF_LINE);
            printf("\n");
        }
        md_vm_arena_destroy(malloc_);
    }
    EXPECT_EQ(failed, (size_t)0);
    md_vm_arena_destroy(alloc);
}

// A peptide from a structure file: the backbone as peptide groups and alpha carbons, each residue's side chain as
// its chemistry has it
UTEST(fragment, the_residues_of_a_protein) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE, NULL));
    ASSERT_TRUE(md_system_infer(&sys, &st, MD_SYSTEM_INFER_ALL, NULL));

    // The side chain of each amino acid within a chain: its fragments of no backbone, by heavy atoms and name
    static const struct { const char* res; const char* side; } want[] = {
        { "ALA", "1 methyl" },
        { "ARG", "3 C3 linker | 4 guanidine" },
        { "ASN", "1 methylene | 3 amide" },
        { "ASP", "1 methylene | 3 carboxylate" },
        { "GLN", "2 ethylene | 3 amide" },
        { "GLU", "2 ethylene | 3 carboxylate" },
        { "GLY", "" },
        { "HIS", "1 methylene | 5 imidazolyl" },
        { "ILE", "4 sec-butyl" },
        { "LEU", "4 isobutyl" },
        { "LYS", "1 ammonium | 4 C4 linker" },
        { "MET", "2 ethylene | 2 thioether" },
        { "PHE", "1 methylene | 6 phenyl" },
        { "PRO", "3 C3 linker" },
        { "SER", "1 hydroxyl | 1 methylene" },
        { "THR", "1 hydroxyl | 2 ethylidene" },
        { "TRP", "1 methylene | 9 indolyl" },
        { "TYR", "1 hydroxyl | 1 methylene | 6 phenyl" },
        { "VAL", "3 isopropyl" },
    };
    for (size_t w = 0; w < MD_COUNTOF(want); ++w) {
        // The first residue of the name with a peptide group on either side
        bool found = false;
        for (size_t r = 0; r < sys.residue.count && !found; ++r) {
            if (!md_str_eq_cstr(md_system_residue_name(&sys, r), want[w].res)) continue;
            md_fragment_idx_t frags[16];
            const size_t num = MD_MIN(md_system_residue_fragments(frags, 16, &sys, r), (size_t)16);
            size_t num_peptide = 0, num_alpha = 0;
            char side[8][FRAG_REF_LINE];
            size_t num_side = 0;
            for (size_t k = 0; k < num; ++k) {
                const md_str_t name = md_system_fragment_name(&sys, frags[k]);
                if (md_str_eq_cstr(name, "peptide")) { num_peptide++; continue; }
                if (md_str_eq_cstr(name, "alpha carbon")) { num_alpha++; continue; }
                size_t heavy = 0;
                for (md_fragment_atom_iter_t it = md_fragment_atom_iter(&sys, frags[k]); md_fragment_atom_iter_next(&it);) {
                    heavy += md_atom_atomic_number(&sys.atom, (size_t)it.atom) > 1;
                }
                if (num_side < 8) snprintf(side[num_side++], FRAG_REF_LINE, "%zu %.*s", heavy, (int)name.len, name.ptr);
            }
            if (num_peptide != 2 || num_alpha != 1) continue;
            found = true;
            qsort(side, num_side, FRAG_REF_LINE, frag_cmp_line);
            char got[256] = "";
            for (size_t k = 0; k < num_side; ++k) {
                if (k) strcat(got, " | ");
                strcat(got, side[k]);
            }
            EXPECT_STREQ(want[w].side, got);
            if (strcmp(want[w].side, got)) printf("  %s: want '%s', got '%s'\n", want[w].res, want[w].side, got);
        }
        EXPECT_TRUE(found);
    }
    md_vm_arena_destroy(alloc);
}
