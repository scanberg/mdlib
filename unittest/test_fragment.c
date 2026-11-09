#include "utest.h"

#include <md_fragment.h>
#include <md_smiles.h>
#include <md_system.h>
#include <md_util.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

// A system built from SMILES the way a quantum chemistry file arrives: every hydrogen an atom of its own
// (right after its heavy atom), bond orders and aromaticity known, metal - non metal bonds coordinate.
// Residues, when given, as counts of heavy atoms in SMILES order.

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
        const md_str_t sym = md_util_element_symbol(z);
        char name[8];
        snprintf(name, sizeof(name), "%.*s%zu", (int)sym.len, sym.ptr, i + 1);
        const int32_t ai = map[i];
        sys->atom.type_idx[ai] = md_particle_type_find_or_add(&sys->atom.type, md_str_from_cstr(name), z, md_util_element_atomic_mass(z), md_util_element_vdw_radius(z), 0, 0, alloc);
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
            md_array_push(sys->residue.flags, md_residue_flags_set_kind(MD_RESIDUE_FLAG_NONE, MD_RESIDUE_KIND_AMINO_ACID), alloc);
            md_array_push(sys->residue.atom_offset, end, alloc);
        }
        sys->residue.count = (size_t)num_res;
    }

    md_bond_build_connectivity(&sys->bond, n, alloc);
    md_util_system_infer_rings(sys);
    return true;
}

typedef struct frag_case_t {
    md_allocator_i* alloc;
    md_system_t sys;
    md_fragments_t frags;
} frag_case_t;

static bool frag_case(frag_case_t* t, const char* smiles) {
    t->alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    t->sys = (md_system_t){ .alloc = t->alloc };
    return frag_system_from_smiles(&t->sys, smiles, NULL, NULL, 0) && md_fragments_compute(&t->frags, &t->sys, t->alloc);
}

static void frag_case_free(frag_case_t* t) {
    md_arena_allocator_destroy(t->alloc);
}

static int frag_find(const md_fragments_t* f, const char* name) {
    for (size_t i = 0; i < f->count; ++i) if (strcmp(f->name + i * MD_FRAGMENT_NAME_CAP, name) == 0) return (int)i;
    return -1;
}

static size_t frag_heavy_count(const frag_case_t* t, int frag) {
    size_t c = 0;
    for (uint32_t k = t->frags.atom_beg[frag]; k < t->frags.atom_end[frag]; ++k) {
        c += md_atom_atomic_number(&t->sys.atom, (size_t)t->frags.atom_idx[k]) > 1 ? 1 : 0;
    }
    return c;
}

// What every result must satisfy: pre-order, nested contiguous ranges, every atom in exactly one leaf
static void frag_check_tree(int* utest_result, const frag_case_t* t) {
    const md_fragments_t* f = &t->frags;
    size_t* seen = calloc(f->num_atoms, sizeof(size_t));
    for (size_t i = 0; i < f->count; ++i) {
        const uint32_t p = f->parent[i];
        if (p != MD_FRAGMENT_NONE) {
            EXPECT_LT(p, (uint32_t)i);
            EXPECT_LE(f->atom_beg[p], f->atom_beg[i]);
            EXPECT_GE(f->atom_end[p], f->atom_end[i]);
        }
        EXPECT_LT(f->atom_beg[i], f->atom_end[i]);
        if (md_fragment_is_leaf(f, i)) {
            for (uint32_t k = f->atom_beg[i]; k < f->atom_end[i]; ++k) {
                const int32_t a = f->atom_idx[k];
                seen[a]++;
                EXPECT_EQ(f->atom_fragment[a], (uint32_t)i);
            }
        }
    }
    for (size_t a = 0; a < f->num_atoms; ++a) EXPECT_EQ(seen[a], (size_t)1);
    free(seen);
}

#define EXPECT_FRAG(T, NAME, KIND, HEAVY) do { \
    const int i_ = frag_find(&(T).frags, NAME); \
    EXPECT_GE(i_, 0); \
    if (i_ >= 0) { \
        EXPECT_EQ((int)(T).frags.kind[i_], (int)(KIND)); \
        EXPECT_EQ(frag_heavy_count(&(T), i_), (size_t)(HEAVY)); \
    } else { \
        printf("  no fragment '%s'; have:", NAME); \
        for (size_t k_ = 0; k_ < (T).frags.count; ++k_) printf(" [%s]", (T).frags.name + k_ * MD_FRAGMENT_NAME_CAP); \
        printf("\n"); \
    } \
} while (0)

UTEST(fragment, amine_linker_phenyl) {
    // The S1 of test_data/vlx/mol.h5 moves its electron from the amine to the ring through this split
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "CNC(C)Cc1ccccc1"));
    frag_check_tree(utest_result, &t);
    ASSERT_EQ(t.frags.count, (size_t)4);
    EXPECT_EQ(t.frags.kind[0], (uint8_t)MD_FRAGMENT_MOLECULE);
    EXPECT_STREQ(t.frags.name, "C10H15N");
    EXPECT_FRAG(t, "amine",     MD_FRAGMENT_GROUP,       2);   // N and its methyl
    EXPECT_FRAG(t, "C3 linker", MD_FRAGMENT_LINKER,      3);
    EXPECT_FRAG(t, "phenyl",    MD_FRAGMENT_RING_SYSTEM, 6);
    // in the order of their first atoms
    EXPECT_STREQ(t.frags.name + 1 * MD_FRAGMENT_NAME_CAP, "amine");
    EXPECT_STREQ(t.frags.name + 3 * MD_FRAGMENT_NAME_CAP, "phenyl");
    frag_case_free(&t);
}

UTEST(fragment, single_unit_molecules) {
    // DMF: the methyls go with the amide, which is then the whole molecule
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "CN(C)C=O"));
    frag_check_tree(utest_result, &t);
    ASSERT_EQ(t.frags.count, (size_t)1);
    EXPECT_FRAG(t, "amide", MD_FRAGMENT_GROUP, 5);
    EXPECT_EQ(t.frags.atom_end[0] - t.frags.atom_beg[0], 12u);
    frag_case_free(&t);

    const struct { const char* smiles; const char* name; int kind; } cases[] = {
        { "CCOC(C)=O",        "ester",       MD_FRAGMENT_GROUP },        // the ethyl and the methyl go with it
        { "c1ccc2ccccc2c1",   "naphthyl",    MD_FRAGMENT_RING_SYSTEM },  // fused: one unit
        { "Clc1ccccc1",       "phenyl",      MD_FRAGMENT_RING_SYSTEM },  // a single halogen goes with its ring
        { "Cc1ccncc1",        "pyridyl",     MD_FRAGMENT_RING_SYSTEM },  // and so does a methyl
        { "O=C1C=CC(=O)C=C1", "benzoquinone",MD_FRAGMENT_RING_SYSTEM },  // the exocyclic O are the ring's
        { "c1ccc2c(c1)ccn2",  "indolyl",     MD_FRAGMENT_RING_SYSTEM },
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        frag_case_t u = {0};
        ASSERT_TRUE(frag_case(&u, cases[i].smiles));
        frag_check_tree(utest_result, &u);
        EXPECT_EQ(u.frags.count, (size_t)1);
        EXPECT_STREQ(u.frags.name, cases[i].name);
        EXPECT_EQ((int)u.frags.kind[0], cases[i].kind);
        frag_case_free(&u);
    }
}

UTEST(fragment, donor_bridge_acceptor_stays_apart) {
    // 4-(dimethylamino)benzonitrile, the textbook twisted charge transfer molecule
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "CN(C)c1ccc(cc1)C#N"));
    frag_check_tree(utest_result, &t);
    EXPECT_EQ(t.frags.count, (size_t)4);
    EXPECT_FRAG(t, "amine",   MD_FRAGMENT_GROUP,       3);
    EXPECT_FRAG(t, "phenyl",  MD_FRAGMENT_RING_SYSTEM, 6);
    EXPECT_FRAG(t, "nitrile", MD_FRAGMENT_GROUP,       2);
    frag_case_free(&t);

    t = (frag_case_t){0};
    ASSERT_TRUE(frag_case(&t, "c1ccccc1[N+](=O)[O-]"));
    frag_check_tree(utest_result, &t);
    EXPECT_FRAG(t, "phenyl", MD_FRAGMENT_RING_SYSTEM, 6);
    EXPECT_FRAG(t, "nitro",  MD_FRAGMENT_GROUP,       3);
    frag_case_free(&t);

    t = (frag_case_t){0};
    ASSERT_TRUE(frag_case(&t, "FC(F)(F)c1ccccc1"));
    frag_check_tree(utest_result, &t);
    EXPECT_FRAG(t, "trifluoromethyl", MD_FRAGMENT_GROUP,       4);
    EXPECT_FRAG(t, "phenyl",          MD_FRAGMENT_RING_SYSTEM, 6);
    frag_case_free(&t);
}

UTEST(fragment, repeated_names_are_numbered) {
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "O=C(c1ccccc1)c1ccccc1"));
    frag_check_tree(utest_result, &t);
    EXPECT_FRAG(t, "carbonyl", MD_FRAGMENT_GROUP,       2);
    EXPECT_FRAG(t, "phenyl 1", MD_FRAGMENT_RING_SYSTEM, 6);
    EXPECT_FRAG(t, "phenyl 2", MD_FRAGMENT_RING_SYSTEM, 6);
    frag_case_free(&t);

    t = (frag_case_t){0};
    ASSERT_TRUE(frag_case(&t, "c1ccc(cc1)-c1ccccc1"));   // biphenyl is two rings, not one system
    frag_check_tree(utest_result, &t);
    EXPECT_FRAG(t, "phenyl 1", MD_FRAGMENT_RING_SYSTEM, 6);
    EXPECT_FRAG(t, "phenyl 2", MD_FRAGMENT_RING_SYSTEM, 6);
    frag_case_free(&t);
}

UTEST(fragment, identical_small_molecules_form_a_set) {
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "O.O.O.CC#N"));
    frag_check_tree(utest_result, &t);
    EXPECT_FRAG(t, "water \xc3\x97" "3", MD_FRAGMENT_MOLECULE_SET, 3);
    EXPECT_FRAG(t, "water 2", MD_FRAGMENT_GROUP, 1);
    EXPECT_FRAG(t, "nitrile", MD_FRAGMENT_GROUP, 3);   // acetonitrile: its methyl goes with it
    const int w = frag_find(&t.frags, "water \xc3\x97" "3");
    if (w >= 0) EXPECT_EQ(t.frags.atom_end[w] - t.frags.atom_beg[w], 9u);
    frag_case_free(&t);
}

UTEST(fragment, metal_complex) {
    // Coordination bonds join nothing: the metal and each ligand are molecules of their own
    frag_case_t t = {0};
    ASSERT_TRUE(frag_case(&t, "[Zn](n1ccccc1)n1ccccc1"));
    frag_check_tree(utest_result, &t);
    EXPECT_FRAG(t, "Zn", MD_FRAGMENT_METAL, 1);
    EXPECT_FRAG(t, "pyridyl \xc3\x97" "2", MD_FRAGMENT_MOLECULE_SET, 12);
    EXPECT_FRAG(t, "pyridyl 1", MD_FRAGMENT_RING_SYSTEM, 6);
    frag_case_free(&t);
}

UTEST(fragment, residues) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    frag_case_t t = { .alloc = alloc, .sys = { .alloc = alloc } };
    const char* names[] = { "GLY", "GLY" };
    const int heavy[] = { 4, 5 };
    ASSERT_TRUE(frag_system_from_smiles(&t.sys, "NCC(=O)NCC(=O)O", names, heavy, 2));
    ASSERT_TRUE(md_fragments_compute(&t.frags, &t.sys, alloc));
    frag_check_tree(utest_result, &t);
    EXPECT_FRAG(t, "protein", MD_FRAGMENT_MOLECULE, 9);
    EXPECT_FRAG(t, "GLY 1", MD_FRAGMENT_RESIDUE, 4);
    EXPECT_FRAG(t, "GLY 2", MD_FRAGMENT_RESIDUE, 5);
    md_arena_allocator_destroy(alloc);
}
