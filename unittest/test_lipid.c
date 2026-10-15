#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_tpr.h>
#include <md_smiles.h>

#include <stdio.h>
#include <string.h>

// LIPIDS: the heads and tails of md_util_system_infer_lipids, on the membranes of the test data (a united atom DPPC
// bilayer whose residues are named DPP, Martini 2 and 3 membranes with cholesterol) and on molecules built from SMILES
// (heavy atoms and bonds only, as a united atom topology has them): sterols, a bolalipid, a cardiolipin, and molecules
// which are no lipids.

static size_t lipid_num_heads(const md_system_t* sys, size_t i) { return sys->lipid.head_offset[i + 1] - sys->lipid.head_offset[i]; }
static size_t lipid_num_tails(const md_system_t* sys, size_t i) { return sys->lipid.tail_offset[i + 1] - sys->lipid.tail_offset[i]; }
static size_t tail_spine_len(const md_system_t* sys, size_t t) { return md_index_range_size(&sys->lipid.tail.spine, t); }

// What holds for every lipid of every system: the heads and tails partition the heavy atoms of the lipid, the spine is a
// path of bonds which begins at an atom bonded to its anchor, the tail atoms carry the flag and nothing else does, and
// the components are lipids.
static bool lipid_invariants_hold(const md_system_t* sys) {
    const md_lipid_data_t* L = &sys->lipid;
    if (md_array_size(L->head_offset) != L->count + 1 || md_array_size(L->tail_offset) != L->count + 1) return false;
    if (md_index_data_num_ranges(&L->head.atoms) != L->head.count) return false;
    if (md_index_data_num_ranges(&L->tail.spine) != L->tail.count) return false;
    if (md_index_data_num_ranges(&L->tail.atoms) != L->tail.count) return false;

    size_t num_flagged = 0;
    for (size_t i = 0; i < sys->atom.count; ++i) num_flagged += (md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_LIPID_TAIL) != 0;
    size_t num_tail_atoms = 0;

    for (size_t l = 0; l < L->count; ++l) {
        const md_urange_t cr = L->comp_range[l];
        for (uint32_t c = cr.beg; c < cr.end; ++c) {
            if (md_component_kind(&sys->component, c) != MD_COMPONENT_KIND_LIPID) return false;
        }
        for (uint32_t t = L->tail_offset[l]; t < L->tail_offset[l + 1]; ++t) {
            const int32_t* spine = md_index_range_beg(&L->tail.spine, t);
            const size_t len = md_index_range_size(&L->tail.spine, t);
            if (len == 0) return false;
            if (md_bond_find(&sys->bond, spine[0], L->tail.anchor_beg[t]) < 0) return false;
            for (size_t k = 1; k < len; ++k) {
                if (md_bond_find(&sys->bond, spine[k - 1], spine[k]) < 0) return false;
            }
            if (L->tail.anchor_end[t] != -1 && md_bond_find(&sys->bond, spine[len - 1], L->tail.anchor_end[t]) < 0) return false;
            if (L->tail.head_beg[t] < L->head_offset[l] || L->tail.head_beg[t] >= L->head_offset[l + 1]) return false;
            const int32_t* atoms = md_index_range_beg(&L->tail.atoms, t);
            const size_t num_atoms = md_index_range_size(&L->tail.atoms, t);
            for (size_t k = 0; k < num_atoms; ++k) {
                if (!(md_atom_flags(&sys->atom, atoms[k]) & MD_ATOM_FLAG_LIPID_TAIL)) return false;
                if (k > 0 && atoms[k] <= atoms[k - 1]) return false;
            }
            num_tail_atoms += num_atoms;
        }
        for (uint32_t h = L->head_offset[l]; h < L->head_offset[l + 1]; ++h) {
            const int32_t* atoms = md_index_range_beg(&L->head.atoms, h);
            const size_t num_atoms = md_index_range_size(&L->head.atoms, h);
            if (num_atoms == 0) return false;
            for (size_t k = 0; k < num_atoms; ++k) {
                if (md_atom_flags(&sys->atom, atoms[k]) & MD_ATOM_FLAG_LIPID_TAIL) return false;
            }
        }
    }
    return num_flagged == num_tail_atoms;
}

UTEST(lipid, dppc_united_atom_pdb) {
    // 64 DPPC named DPP (the three letters of a PDB residue name), atoms numbered C1 .. C50: a lipid by its structure
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb"), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));

    EXPECT_EQ(64, (int)sys.lipid.count);
    EXPECT_EQ(64, (int)sys.lipid.head.count);
    EXPECT_EQ(128, (int)sys.lipid.tail.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    for (size_t l = 0; l < sys.lipid.count; ++l) {
        EXPECT_EQ(1, (int)lipid_num_heads(&sys, l));
        ASSERT_EQ(2, (int)lipid_num_tails(&sys, l));
        // Choline, phosphate, glycerol and the two ester groups
        EXPECT_EQ(20, (int)md_index_range_size(&sys.lipid.head.atoms, sys.lipid.head_offset[l]));
        for (uint32_t t = sys.lipid.tail_offset[l]; t < sys.lipid.tail_offset[l + 1]; ++t) {
            EXPECT_EQ(MD_LIPID_TAIL_CHAIN, sys.lipid.tail.kind[t]);
            EXPECT_EQ(-1, sys.lipid.tail.anchor_end[t]);
            // C2 .. C16 of a palmitoyl chain, hanging off its carbonyl carbon
            EXPECT_EQ(15, (int)tail_spine_len(&sys, t));
            EXPECT_EQ(MD_Z_C, md_atom_atomic_number(&sys.atom, sys.lipid.tail.anchor_beg[t]));
        }
    }
    // The first lipid: tails C17 .. C31 off C15 and C36 .. C50 off C34
    const int32_t* spine0 = md_index_range_beg(&sys.lipid.tail.spine, 0);
    EXPECT_TRUE(str_eq_cstr(md_atom_name(&sys.atom, spine0[0]), "C17"));
    EXPECT_TRUE(str_eq_cstr(md_atom_name(&sys.atom, spine0[14]), "C31"));
    EXPECT_TRUE(str_eq_cstr(md_atom_name(&sys.atom, sys.lipid.tail.anchor_beg[0]), "C15"));
    const int32_t* spine1 = md_index_range_beg(&sys.lipid.tail.spine, 1);
    EXPECT_TRUE(str_eq_cstr(md_atom_name(&sys.atom, spine1[0]), "C36"));
    EXPECT_TRUE(str_eq_cstr(md_atom_name(&sys.atom, spine1[14]), "C50"));

    // The water is no lipid
    size_t num_lipid_comps = 0;
    for (size_t c = 0; c < sys.component.count; ++c) num_lipid_comps += md_component_kind(&sys.component, c) == MD_COMPONENT_KIND_LIPID;
    EXPECT_EQ(64, (int)num_lipid_comps);

    // Inferring again gives the same: a lipid by its structure is not taken for one by name the second time
    ASSERT_TRUE(md_util_system_infer_lipids(&sys));
    EXPECT_EQ(64, (int)sys.lipid.count);
    EXPECT_EQ(128, (int)sys.lipid.tail.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    for (size_t t = 0; t < sys.lipid.tail.count; ++t) EXPECT_EQ(15, (int)tail_spine_len(&sys, t));

    md_arena_allocator_destroy(alloc);
}

// The lipids of a Martini membrane: two tailed lipids of 4 or 5 beads per tail hanging off GL1 or GL2, and cholesterols
// with the hydroxyl bead for a head and the skeleton for a tail. Prints the first thing that is not so.
static bool martini_membrane_ok(const md_system_t* sys, int num_two_tailed, int num_sterols) {
    if ((int)sys->lipid.count != num_two_tailed + num_sterols) {
        printf("lipids: %zu, expected %d\n", sys->lipid.count, num_two_tailed + num_sterols);
        return false;
    }
    if (!lipid_invariants_hold(sys)) {
        printf("invariants do not hold\n");
        return false;
    }
    int two_tailed = 0, sterols = 0;
    for (size_t l = 0; l < sys->lipid.count; ++l) {
        const str_t name = md_component_name(&sys->component, sys->lipid.comp_range[l].beg);
        if (lipid_num_heads(sys, l) != 1) {
            printf("lipid %zu (" STR_FMT "): %zu heads\n", l, STR_ARG(name), lipid_num_heads(sys, l));
            return false;
        }
        if (str_eq_cstr(name, "CHOL")) {
            const uint32_t t = sys->lipid.tail_offset[l];
            const int32_t* head = md_index_range_beg(&sys->lipid.head.atoms, sys->lipid.head_offset[l]);
            if (lipid_num_tails(sys, l) != 1 || sys->lipid.tail.kind[t] != MD_LIPID_TAIL_SKELETON ||
                md_index_range_size(&sys->lipid.head.atoms, sys->lipid.head_offset[l]) != 1 || !str_eq_cstr(md_atom_name(&sys->atom, head[0]), "ROH")) {
                printf("lipid %zu (CHOL): %zu tails, kind %d, head of %zu beads\n", l, lipid_num_tails(sys, l), (int)sys->lipid.tail.kind[t], md_index_range_size(&sys->lipid.head.atoms, sys->lipid.head_offset[l]));
                return false;
            }
            sterols += 1;
        } else {
            if (lipid_num_tails(sys, l) != 2) {
                printf("lipid %zu (" STR_FMT "): %zu tails\n", l, STR_ARG(name), lipid_num_tails(sys, l));
                return false;
            }
            for (uint32_t t = sys->lipid.tail_offset[l]; t < sys->lipid.tail_offset[l + 1]; ++t) {
                const size_t len = tail_spine_len(sys, t);
                const str_t anchor = md_atom_name(&sys->atom, sys->lipid.tail.anchor_beg[t]);
                if (sys->lipid.tail.kind[t] != MD_LIPID_TAIL_CHAIN || (len != 4 && len != 5) || !(str_eq_cstr(anchor, "GL1") || str_eq_cstr(anchor, "GL2"))) {
                    printf("lipid %zu (" STR_FMT "): tail of kind %d, %zu beads, off " STR_FMT "\n", l, STR_ARG(name), (int)sys->lipid.tail.kind[t], len, STR_ARG(anchor));
                    return false;
                }
            }
            two_tailed += 1;
        }
    }
    return two_tailed == num_two_tailed && sterols == num_sterols;
}

UTEST(lipid, martini3_membrane) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.gro")));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    // POPC 32, POPE 16, POPS 8; CHOL 16
    EXPECT_TRUE(martini_membrane_ok(&sys, 32 + 16 + 8, 16));
    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, martini3_membrane_tpr) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.tpr")));
    // The bonds are the topology's
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT));
    EXPECT_TRUE(martini_membrane_ok(&sys, 32 + 16 + 8, 16));
    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, martini2_membrane_protein) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini2_membrane_protein.gro")));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    // DPPC 50, DOPC 50; CHOL 24. The protein has no lipids.
    EXPECT_TRUE(martini_membrane_ok(&sys, 50 + 50, 24));
    md_arena_allocator_destroy(alloc);
}

// A system of one component from the heavy atoms and bonds of a SMILES string, classified and with its lipids inferred
static bool lipid_system_from_smiles(md_system_t* sys, const char* comp_name, const char* smiles) {
    md_allocator_i* alloc = sys->alloc;
    md_smiles_t g = {0};
    md_smiles_error_t err = {0};
    if (!md_smiles_parse(&g, str_from_cstr(smiles), alloc, &err)) {
        printf("SMILES error in '%s' at %zu: %s\n", smiles, err.offset, err.message);
        return false;
    }
    const size_t n = g.num_atoms;
    sys->atom.count = n;
    md_atom_type_add(&sys->atom.type, STR_LIT(""), (str_t){0}, 0, 0, 0, 0, 0, alloc);
    md_array_resize(sys->atom.type_idx, n, alloc);
    md_array_resize(sys->atom.flags, n, alloc);
    MEMSET(sys->atom.flags, 0, md_array_bytes(sys->atom.flags));
    for (size_t i = 0; i < n; ++i) {
        const md_atomic_number_t z = g.atoms[i].z;
        char name[8];
        const str_t sym = md_util_element_symbol(z);
        snprintf(name, sizeof(name), "%.*s%zu", (int)sym.len, sym.ptr, i + 1);
        sys->atom.type_idx[i] = md_atom_type_find_or_add(&sys->atom.type, str_from_cstr(name), z, md_util_element_atomic_mass(z), md_util_element_vdw_radius(z), 0, 0, alloc);
        // Aromatic as written (lowercase), as md_chem_perceive would flag it
        if (g.atoms[i].flags & MD_SMILES_ATOM_AROMATIC) sys->atom.flags[i] |= MD_ATOM_FLAG_AROMATIC;
    }
    sys->component.count = 1;
    md_array_push(sys->component.name, make_label(str_from_cstr(comp_name)), alloc);
    md_array_push(sys->component.seq_id, 1, alloc);
    md_array_push(sys->component.flags, 0, alloc);
    md_array_push(sys->component.atom_offset, 0, alloc);
    md_array_push(sys->component.atom_offset, (uint32_t)n, alloc);
    for (size_t b = 0; b < g.num_bonds; ++b) {
        md_atom_pair_t p = { .idx = { (int32_t)g.bonds[b].a, (int32_t)g.bonds[b].b } };
        md_array_push(sys->bond.pairs, p, alloc);
        md_array_push(sys->bond.flags, md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY), alloc);
        sys->bond.count += 1;
    }
    md_bond_build_connectivity(&sys->bond, n, alloc);
    md_util_system_infer_rings(sys);
    md_util_system_infer_comp_flags(sys);
    return md_util_system_infer_lipids(sys);
}

#define CHOLESTEROL     "CC(CCCC(C)C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"
#define POPC            "CCCCCCCCC=CCCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCC)COP(=O)([O-])OCC[N+](C)(C)C"
#define CARDIOLIPIN     "OC(COP(=O)([O-])OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC)COP(=O)([O-])OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
#define SPHINGOMYELIN   "CCCCCCCCCCCCCC=CC(O)C(COP(=O)([O-])OCC[N+](C)(C)C)NC(=O)CCCCCCCCCCCCCCC"
// Caldarchaeol (GDGT-0): two glycerols joined by two biphytanyl chains through ether bonds, one of them closing the ring
#define BIPHYTANYL      "CCC(C)CCCC(C)CCCC(C)CCCC(C)CCC(C)CCCC(C)CCCC(C)CCCC(C)CC"
#define CALDARCHAEOL    "OCC(O1)CO" BIPHYTANYL "OCC(CO)O" BIPHYTANYL "1"

UTEST(lipid, cholesterol) {
    // By name (CHARMM) and by its structure alone
    const char* names[] = { "CHL1", "XXX" };
    for (size_t i = 0; i < ARRAY_SIZE(names); ++i) {
        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        ASSERT_TRUE(lipid_system_from_smiles(&sys, names[i], CHOLESTEROL));
        ASSERT_EQ(1, (int)sys.lipid.count);
        EXPECT_TRUE(lipid_invariants_hold(&sys));
        EXPECT_EQ(MD_COMPONENT_KIND_LIPID, md_component_kind(&sys.component, 0));
        ASSERT_EQ(1, (int)sys.lipid.tail.count);
        EXPECT_EQ(MD_LIPID_TAIL_SKELETON, sys.lipid.tail.kind[0]);
        // The head is the hydroxyl, the anchor the C3 it is on
        ASSERT_EQ(1, (int)sys.lipid.head.count);
        const int32_t* head = md_index_range_beg(&sys.lipid.head.atoms, 0);
        EXPECT_EQ(2, (int)md_index_range_size(&sys.lipid.head.atoms, 0));
        EXPECT_EQ(MD_Z_O, md_atom_atomic_number(&sys.atom, head[1]));
        // The skeleton is the 26 carbons but C3, and its spine runs the shortest way through the rings from C4 (or C2)
        // out to the end of the side chain: C4 C5 C6 C7 C8 C14 C13 C17 C20 C22 C23 C24 C25 C26
        EXPECT_EQ(26, (int)md_index_range_size(&sys.lipid.tail.atoms, 0));
        const int32_t* spine = md_index_range_beg(&sys.lipid.tail.spine, 0);
        const size_t len = tail_spine_len(&sys, 0);
        EXPECT_EQ(14, (int)len);
        const size_t tip_degree = md_bond_conn_count(&sys.bond, spine[len - 1]);
        EXPECT_EQ(1, (int)tip_degree);
        md_arena_allocator_destroy(alloc);
    }
}

UTEST(lipid, popc_by_structure) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    ASSERT_TRUE(lipid_system_from_smiles(&sys, "XYZ", POPC));
    ASSERT_EQ(1, (int)sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    ASSERT_EQ(2, (int)sys.lipid.tail.count);
    // Oleoyl C2 .. C18 (its double bond is on the spine) and palmitoyl C2 .. C16
    EXPECT_EQ(17, (int)tail_spine_len(&sys, 0));
    EXPECT_EQ(15, (int)tail_spine_len(&sys, 1));
    EXPECT_EQ(1, (int)sys.lipid.head.count);
    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, fatty_acid_by_structure) {
    // Palmitic acid as a ligand of a crystal structure (PLM): a lipid by its structure
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    ASSERT_TRUE(lipid_system_from_smiles(&sys, "PLM", "CCCCCCCCCCCCCCCC(=O)O"));
    ASSERT_EQ(1, (int)sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    ASSERT_EQ(1, (int)sys.lipid.tail.count);
    EXPECT_EQ(15, (int)tail_spine_len(&sys, 0));
    // The head is the carboxyl group
    EXPECT_EQ(3, (int)md_index_range_size(&sys.lipid.head.atoms, 0));
    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, cardiolipin) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    ASSERT_TRUE(lipid_system_from_smiles(&sys, "CDL", CARDIOLIPIN));
    ASSERT_EQ(1, (int)sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    EXPECT_EQ(1, (int)sys.lipid.head.count);
    ASSERT_EQ(4, (int)sys.lipid.tail.count);
    for (int t = 0; t < 4; ++t) EXPECT_EQ(15, (int)tail_spine_len(&sys, t));
    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, sphingomyelin) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    ASSERT_TRUE(lipid_system_from_smiles(&sys, "SSM", SPHINGOMYELIN));
    ASSERT_EQ(1, (int)sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    // The sphingoid chain past C3 (C4 .. C18) and the N-acyl chain (C2 .. C16)
    ASSERT_EQ(2, (int)sys.lipid.tail.count);
    EXPECT_EQ(15, (int)tail_spine_len(&sys, 0));
    EXPECT_EQ(15, (int)tail_spine_len(&sys, 1));
    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, bolalipid) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    ASSERT_TRUE(lipid_system_from_smiles(&sys, "GDG", CALDARCHAEOL));
    ASSERT_EQ(1, (int)sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    // Two glycerol heads, two chains which each span them
    ASSERT_EQ(2, (int)sys.lipid.head.count);
    ASSERT_EQ(2, (int)sys.lipid.tail.count);
    for (int t = 0; t < 2; ++t) {
        EXPECT_EQ(MD_LIPID_TAIL_CHAIN, sys.lipid.tail.kind[t]);
        EXPECT_NE(-1, sys.lipid.tail.anchor_end[t]);
        EXPECT_NE(sys.lipid.tail.head_beg[t], sys.lipid.tail.head_end[t]);
        // The 32 carbons of the main chain of biphytanyl but the two bonded to the ether oxygens; the 8 methyl branches
        // are in the tail, off the spine
        EXPECT_EQ(30, (int)tail_spine_len(&sys, t));
        EXPECT_EQ(38, (int)md_index_range_size(&sys.lipid.tail.atoms, t));
    }
    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, not_lipids) {
    const char* smiles[] = {
        "CCCCCCCCCCCCCCCC",                         // Hexadecane: no head
        "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",      // Testosterone: polar at both ends of its skeleton
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",               // Ibuprofen
        "OCC1OC(O)C(O)C(O)C1O",                     // Glucose
        "CCCCCCC(=O)O",                             // Heptanoic acid: a chain too short for a lipid of no name
        "CC(C)CCCC(C)CCCC(C)CCCC(C)=CCOC(=O)CCC(=O)O[Mg]",  // A phytyl ester with a metal in its head (as chlorophyll)
        "CCCCCCCCCCCCCCCCOC(=O)CCc1ccc[nH]1",       // An aromatic head (as the porphyrin of a pheophytin)
    };
    for (size_t i = 0; i < ARRAY_SIZE(smiles); ++i) {
        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        ASSERT_TRUE(lipid_system_from_smiles(&sys, "XXX", smiles[i]));
        EXPECT_EQ(0, (int)sys.lipid.count);
        EXPECT_NE(MD_COMPONENT_KIND_LIPID, md_component_kind(&sys.component, 0));
        md_arena_allocator_destroy(alloc);
    }
}

#include <md_filter.h>
#include <core/md_bitset.h>

UTEST(lipid, script) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb"), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));

    char err[256] = {0};
    bool is_dynamic = false;
    md_bitset_t bf = {0};

    // 64 lipids of 50 united atoms: 2 x 15 in the tails, 20 in the head
    ASSERT_TRUE(md_filter(&bf, STR_LIT("lipid()"), &sys, &st, NULL, &is_dynamic, err, sizeof(err), alloc));
    EXPECT_EQ(64 * 50, (int)md_bitset_count(bf));
    ASSERT_TRUE(md_filter(&bf, STR_LIT("tail()"), &sys, &st, NULL, &is_dynamic, err, sizeof(err), alloc));
    EXPECT_EQ(64 * 30, (int)md_bitset_count(bf));
    ASSERT_TRUE(md_filter(&bf, STR_LIT("headgroup()"), &sys, &st, NULL, &is_dynamic, err, sizeof(err), alloc));
    EXPECT_EQ(64 * 20, (int)md_bitset_count(bf));
    ASSERT_TRUE(md_filter(&bf, STR_LIT("tail() and headgroup()"), &sys, &st, NULL, &is_dynamic, err, sizeof(err), alloc));
    EXPECT_EQ(0, (int)md_bitset_count(bf));
    // The phosphorus of the heads
    ASSERT_TRUE(md_filter(&bf, STR_LIT("headgroup() and element('P')"), &sys, &st, NULL, &is_dynamic, err, sizeof(err), alloc));
    EXPECT_EQ(64, (int)md_bitset_count(bf));
    // One selection per tail, 15 carbons each
    md_array(md_bitset_t) sets = 0;
    ASSERT_TRUE(md_filter_evaluate(&sets, STR_LIT("tail()"), &sys, &st, NULL, &is_dynamic, err, sizeof(err), alloc));
    ASSERT_EQ(128, (int)md_array_size(sets));
    for (size_t i = 0; i < md_array_size(sets); ++i) EXPECT_EQ(15, (int)md_bitset_count(sets[i]));
    // Written without parentheses, as the default representations write them
    ASSERT_TRUE(md_filter(&bf, STR_LIT("headgroup or tail"), &sys, &st, NULL, &is_dynamic, err, sizeof(err), alloc));
    EXPECT_EQ(64 * 50, (int)md_bitset_count(bf));
    // Within a context, the tails fully within it: those of the first lipid
    ASSERT_TRUE(md_filter(&bf, STR_LIT("tail() in residue(1)"), &sys, &st, NULL, &is_dynamic, err, sizeof(err), alloc));
    EXPECT_EQ(30, (int)md_bitset_count(bf));

    md_arena_allocator_destroy(alloc);
}
