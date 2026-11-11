#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <md_system.h>
#include <md_util.h>
#include <md_element.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_tpr.h>
#include <md_smiles.h>
#include <md_itp.h>

#include <stdio.h>
#include <string.h>

// LIPIDS: the heads and tails of md_util_system_infer_lipids, on the membranes of the test data (a united atom DPPC
// bilayer whose residues are named DPP, Martini 2 and 3 membranes with cholesterol) and on molecules built from SMILES
// (heavy atoms and bonds only, as a united atom topology has them): sterols, a bolalipid, a cardiolipin, and molecules
// which are no lipids.

static size_t lipid_num_heads(const md_system_t* sys, size_t i) { return sys->lipid.head_offset[i + 1] - sys->lipid.head_offset[i]; }
static size_t lipid_num_tails(const md_system_t* sys, size_t i) { return sys->lipid.tail_offset[i + 1] - sys->lipid.tail_offset[i]; }
static size_t tail_spine_len(const md_system_t* sys, size_t t) { return md_lipid_tail_spine_size(&sys->lipid, t); }
static size_t lipid_num_residues(const md_system_t* sys, size_t i) { return md_lipid_residue_count(&sys->lipid, i); }
static int32_t lipid_res(const md_system_t* sys, size_t i, size_t k) { return md_lipid_residues(&sys->lipid, i)[k]; }

// What holds for every lipid of every system: the heads and tails partition the heavy atoms of the lipid, the spine is a
// path of bonds which begins at an atom bonded to its anchor, the tail atoms carry the flag and nothing else does, and
// the residues are lipids.
static bool lipid_invariants_hold(const md_system_t* sys) {
    const md_lipid_data_t* L = &sys->lipid;
    if (md_array_size(L->head_offset) != L->count + 1 || md_array_size(L->tail_offset) != L->count + 1) return false;
    if (md_array_size(L->residue_offset) != L->count + 1) return false;
    if (md_array_size(L->head.offset) != L->head.count + 1) return false;
    if (md_array_size(L->tail.spine_offset) != L->tail.count + 1) return false;
    if (md_array_size(L->tail.atom_offset) != L->tail.count + 1) return false;

    size_t num_flagged = 0;
    for (size_t i = 0; i < sys->atom.count; ++i) num_flagged += (md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_LIPID_TAIL) != 0;
    size_t num_tail_atoms = 0;

    for (size_t l = 0; l < L->count; ++l) {
        // Its residues, ascending, lipids, and no other's
        const size_t num_residues = lipid_num_residues(sys, l);
        if (num_residues == 0) return false;
        for (size_t k = 0; k < num_residues; ++k) {
            const int32_t c = lipid_res(sys, l, k);
            if (k > 0 && c <= lipid_res(sys, l, k - 1)) return false;
            if (md_system_residue_kind(sys, (size_t)c) != MD_RESIDUE_KIND_LIPID) return false;
        }
        // The heads and tails are of its residues
        for (uint32_t h = L->head_offset[l]; h < L->head_offset[l + 1]; ++h) {
            const int32_t* atoms = md_lipid_head_atoms(L, h);
            for (size_t k = 0; k < md_lipid_head_atom_count(L, h); ++k) {
                const md_residue_idx_t c = md_system_residue_find_by_atom_idx(sys, (size_t)atoms[k]);
                bool of_it = false;
                for (size_t q = 0; q < num_residues; ++q) of_it |= lipid_res(sys, l, q) == c;
                if (!of_it) return false;
            }
        }
        for (uint32_t t = L->tail_offset[l]; t < L->tail_offset[l + 1]; ++t) {
            const int32_t* spine = md_lipid_tail_spine(L, t);
            const size_t len = md_lipid_tail_spine_size(L, t);
            if (len == 0) return false;
            if (md_bond_find(&sys->bond, spine[0], L->tail.anchor_beg[t]) < 0) return false;
            for (size_t k = 1; k < len; ++k) {
                if (md_bond_find(&sys->bond, spine[k - 1], spine[k]) < 0) return false;
            }
            if (L->tail.anchor_end[t] != -1 && md_bond_find(&sys->bond, spine[len - 1], L->tail.anchor_end[t]) < 0) return false;
            if (L->tail.head_beg[t] < L->head_offset[l] || L->tail.head_beg[t] >= L->head_offset[l + 1]) return false;
            const int32_t* atoms = md_lipid_tail_atoms(L, t);
            const size_t num_atoms = md_lipid_tail_atom_count(L, t);
            for (size_t k = 0; k < num_atoms; ++k) {
                if (!(md_atom_flags(&sys->atom, atoms[k]) & MD_ATOM_FLAG_LIPID_TAIL)) return false;
                if (k > 0 && atoms[k] <= atoms[k - 1]) return false;
            }
            num_tail_atoms += num_atoms;
        }
        for (uint32_t h = L->head_offset[l]; h < L->head_offset[l + 1]; ++h) {
            const int32_t* atoms = md_lipid_head_atoms(L, h);
            const size_t num_atoms = md_lipid_head_atom_count(L, h);
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb"), MD_PDB_OPTION_NONE, NULL));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL, NULL));

    EXPECT_EQ(64, (int)sys.lipid.count);
    EXPECT_EQ(64, (int)sys.lipid.head.count);
    EXPECT_EQ(128, (int)sys.lipid.tail.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    for (size_t l = 0; l < sys.lipid.count; ++l) {
        EXPECT_EQ(1, (int)lipid_num_heads(&sys, l));
        ASSERT_EQ(2, (int)lipid_num_tails(&sys, l));
        // Choline, phosphate, glycerol and the two ester groups
        EXPECT_EQ(20, (int)md_lipid_head_atom_count(&sys.lipid, sys.lipid.head_offset[l]));
        for (uint32_t t = sys.lipid.tail_offset[l]; t < sys.lipid.tail_offset[l + 1]; ++t) {
            EXPECT_EQ(MD_LIPID_TAIL_CHAIN, sys.lipid.tail.kind[t]);
            EXPECT_EQ(-1, sys.lipid.tail.anchor_end[t]);
            // C2 .. C16 of a palmitoyl chain, hanging off its carbonyl carbon
            EXPECT_EQ(15, (int)tail_spine_len(&sys, t));
            EXPECT_EQ(MD_Z_C, md_atom_atomic_number(&sys.atom, sys.lipid.tail.anchor_beg[t]));
        }
    }
    // The first lipid: tails C17 .. C31 off C15 and C36 .. C50 off C34
    const int32_t* spine0 = md_lipid_tail_spine(&sys.lipid, 0);
    EXPECT_TRUE(md_str_eq_cstr(md_atom_name(&sys.atom, spine0[0]), "C17"));
    EXPECT_TRUE(md_str_eq_cstr(md_atom_name(&sys.atom, spine0[14]), "C31"));
    EXPECT_TRUE(md_str_eq_cstr(md_atom_name(&sys.atom, sys.lipid.tail.anchor_beg[0]), "C15"));
    const int32_t* spine1 = md_lipid_tail_spine(&sys.lipid, 1);
    EXPECT_TRUE(md_str_eq_cstr(md_atom_name(&sys.atom, spine1[0]), "C36"));
    EXPECT_TRUE(md_str_eq_cstr(md_atom_name(&sys.atom, spine1[14]), "C50"));

    // The water is no lipid
    size_t num_lipid_residues = 0;
    for (size_t c = 0; c < sys.residue.count; ++c) num_lipid_residues += md_system_residue_kind(&sys, c) == MD_RESIDUE_KIND_LIPID;
    EXPECT_EQ(64, (int)num_lipid_residues);

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
        const md_str_t name = md_system_residue_name(sys, (size_t)lipid_res(sys, l, 0));
        if (lipid_num_heads(sys, l) != 1) {
            printf("lipid %zu (" MD_STR_FMT "): %zu heads\n", l, MD_STR_ARG(name), lipid_num_heads(sys, l));
            return false;
        }
        if (md_str_eq_cstr(name, "CHOL")) {
            const uint32_t t = sys->lipid.tail_offset[l];
            const int32_t* head = md_lipid_head_atoms(&sys->lipid, sys->lipid.head_offset[l]);
            if (lipid_num_tails(sys, l) != 1 || sys->lipid.tail.kind[t] != MD_LIPID_TAIL_SKELETON ||
                md_lipid_head_atom_count(&sys->lipid, sys->lipid.head_offset[l]) != 1 || !md_str_eq_cstr(md_atom_name(&sys->atom, head[0]), "ROH")) {
                printf("lipid %zu (CHOL): %zu tails, kind %d, head of %zu beads\n", l, lipid_num_tails(sys, l), (int)sys->lipid.tail.kind[t], md_lipid_head_atom_count(&sys->lipid, sys->lipid.head_offset[l]));
                return false;
            }
            sterols += 1;
        } else {
            if (lipid_num_tails(sys, l) != 2) {
                printf("lipid %zu (" MD_STR_FMT "): %zu tails\n", l, MD_STR_ARG(name), lipid_num_tails(sys, l));
                return false;
            }
            for (uint32_t t = sys->lipid.tail_offset[l]; t < sys->lipid.tail_offset[l + 1]; ++t) {
                const size_t len = tail_spine_len(sys, t);
                const md_str_t anchor = md_atom_name(&sys->atom, sys->lipid.tail.anchor_beg[t]);
                if (sys->lipid.tail.kind[t] != MD_LIPID_TAIL_CHAIN || (len != 4 && len != 5) || !(md_str_eq_cstr(anchor, "GL1") || md_str_eq_cstr(anchor, "GL2"))) {
                    printf("lipid %zu (" MD_STR_FMT "): tail of kind %d, %zu beads, off " MD_STR_FMT "\n", l, MD_STR_ARG(name), (int)sys->lipid.tail.kind[t], len, MD_STR_ARG(anchor));
                    return false;
                }
            }
            two_tailed += 1;
        }
    }
    return two_tailed == num_two_tailed && sterols == num_sterols;
}

UTEST(lipid, martini3_membrane) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.gro"), NULL));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL, NULL));
    // POPC 32, POPE 16, POPS 8; CHOL 16
    EXPECT_TRUE(martini_membrane_ok(&sys, 32 + 16 + 8, 16));
    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, martini3_membrane_tpr) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.tpr")));
    // The bonds are the topology's
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT, NULL));
    EXPECT_TRUE(martini_membrane_ok(&sys, 32 + 16 + 8, 16));
    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, martini2_membrane_protein) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini2_membrane_protein.gro"), NULL));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL, NULL));
    // DPPC 50, DOPC 50; CHOL 24. The protein has no lipids.
    EXPECT_TRUE(martini_membrane_ok(&sys, 50 + 50, 24));
    md_arena_allocator_destroy(alloc);
}

// A system from the heavy atoms of a SMILES string, in residues of the names and sizes given (one of all atoms when
// num_residues is 0), with its bonds unless told not to, classified and with its lipids inferred
static bool lipid_system_from_smiles_residues(md_system_t* sys, const char* const res_names[], const int res_sizes[], int num_residues, const char* smiles, bool bonds);

static bool lipid_system_from_smiles(md_system_t* sys, const char* res_name, const char* smiles) {
    return lipid_system_from_smiles_residues(sys, &res_name, NULL, 0, smiles, true);
}

static bool lipid_system_from_smiles_residues(md_system_t* sys, const char* const res_names[], const int res_sizes[], int num_residues, const char* smiles, bool bonds) {
    md_allocator_i* alloc = sys->alloc;
    md_smiles_t g = {0};
    md_smiles_error_t err = {0};
    if (!md_smiles_parse(&g, md_str_from_cstr(smiles), alloc, &err)) {
        printf("SMILES error in '%s' at %zu: %s\n", smiles, err.offset, err.message);
        return false;
    }
    const size_t n = g.num_atoms;
    sys->atom.count = n;
    md_particle_type_add(&sys->atom.type, MD_STR_LIT(""), (md_str_t){0}, 0, 0, 0, 0, 0, alloc);
    md_array_resize(sys->atom.type_idx, n, alloc);
    md_array_resize(sys->atom.flags, n, alloc);
    md_memset(sys->atom.flags, 0, md_array_bytes(sys->atom.flags));
    for (size_t i = 0; i < n; ++i) {
        const md_atomic_number_t z = g.atoms[i].z;
        char name[8];
        const md_str_t sym = md_util_element_symbol(z);
        snprintf(name, sizeof(name), "%.*s%zu", (int)sym.len, sym.ptr, i + 1);
        sys->atom.type_idx[i] = md_particle_type_find_or_add(&sys->atom.type, md_str_from_cstr(name), z, md_util_element_atomic_mass(z), md_util_element_vdw_radius(z), 0, 0, alloc);
        // Aromatic as written (lowercase), as md_chem_perceive would flag it
        if (g.atoms[i].flags & MD_SMILES_ATOM_AROMATIC) sys->atom.flags[i] |= MD_ATOM_FLAG_AROMATIC;
    }
    md_array_push(sys->residue.atom_offset, 0, alloc);
    if (num_residues == 0) {
        sys->residue.count = 1;
        md_array_push(sys->residue.name, md_label_make(md_str_from_cstr(res_names[0])), alloc);
        md_array_push(sys->residue.seq_id, 1, alloc);
        md_array_push(sys->residue.flags, 0, alloc);
        md_array_push(sys->residue.atom_offset, (uint32_t)n, alloc);
    } else {
        uint32_t end = 0;
        for (int c = 0; c < num_residues; ++c) {
            end += (uint32_t)res_sizes[c];
            md_array_push(sys->residue.name, md_label_make(md_str_from_cstr(res_names[c])), alloc);
            md_array_push(sys->residue.seq_id, c + 1, alloc);
            md_array_push(sys->residue.flags, 0, alloc);
            md_array_push(sys->residue.atom_offset, end, alloc);
        }
        sys->residue.count = (size_t)num_residues;
        if (end != n) {
            printf("components of %u atoms, the molecule has %zu\n", end, n);
            return false;
        }
    }
    for (size_t b = 0; b < g.num_bonds && bonds; ++b) {
        md_atom_pair_t p = { .idx = { (int32_t)g.bonds[b].a, (int32_t)g.bonds[b].b } };
        md_array_push(sys->bond.pairs, p, alloc);
        md_array_push(sys->bond.flags, md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY), alloc);
        sys->bond.count += 1;
    }
    md_bond_build_connectivity(&sys->bond, n, alloc);
    md_util_system_infer_rings(sys);
    md_util_system_infer_residue_flags(sys);
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
    for (size_t i = 0; i < MD_COUNTOF(names); ++i) {
        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        ASSERT_TRUE(lipid_system_from_smiles(&sys, names[i], CHOLESTEROL));
        ASSERT_EQ(1, (int)sys.lipid.count);
        EXPECT_TRUE(lipid_invariants_hold(&sys));
        EXPECT_EQ(MD_RESIDUE_KIND_LIPID, md_system_residue_kind(&sys, 0));
        ASSERT_EQ(1, (int)sys.lipid.tail.count);
        EXPECT_EQ(MD_LIPID_TAIL_SKELETON, sys.lipid.tail.kind[0]);
        // The head is the hydroxyl, the anchor the C3 it is on
        ASSERT_EQ(1, (int)sys.lipid.head.count);
        const int32_t* head = md_lipid_head_atoms(&sys.lipid, 0);
        EXPECT_EQ(2, (int)md_lipid_head_atom_count(&sys.lipid, 0));
        EXPECT_EQ(MD_Z_O, md_atom_atomic_number(&sys.atom, head[1]));
        // The skeleton is the 26 carbons but C3, and its spine runs the shortest way through the rings from C4 (or C2)
        // out to the end of the side chain: C4 C5 C6 C7 C8 C14 C13 C17 C20 C22 C23 C24 C25 C26
        EXPECT_EQ(26, (int)md_lipid_tail_atom_count(&sys.lipid, 0));
        const int32_t* spine = md_lipid_tail_spine(&sys.lipid, 0);
        const size_t len = tail_spine_len(&sys, 0);
        EXPECT_EQ(14, (int)len);
        const size_t tip_degree = md_bond_conn_count(&sys.bond, spine[len - 1]);
        EXPECT_EQ(1, (int)tip_degree);
        md_arena_allocator_destroy(alloc);
    }
}

UTEST(lipid, popc_by_structure) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    ASSERT_TRUE(lipid_system_from_smiles(&sys, "PLM", "CCCCCCCCCCCCCCCC(=O)O"));
    ASSERT_EQ(1, (int)sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    ASSERT_EQ(1, (int)sys.lipid.tail.count);
    EXPECT_EQ(15, (int)tail_spine_len(&sys, 0));
    // The head is the carboxyl group
    EXPECT_EQ(3, (int)md_lipid_head_atom_count(&sys.lipid, 0));
    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, cardiolipin) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
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
        EXPECT_EQ(38, (int)md_lipid_tail_atom_count(&sys.lipid, t));
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
    for (size_t i = 0; i < MD_COUNTOF(smiles); ++i) {
        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        ASSERT_TRUE(lipid_system_from_smiles(&sys, "XXX", smiles[i]));
        EXPECT_EQ(0, (int)sys.lipid.count);
        EXPECT_NE(MD_RESIDUE_KIND_LIPID, md_system_residue_kind(&sys, 0));
        md_arena_allocator_destroy(alloc);
    }
}

#include <md_filter.h>
#include <core/md_bitset.h>

UTEST(lipid, script) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb"), MD_PDB_OPTION_NONE, NULL));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL, NULL));

    const md_filter_context_t fctx = { &sys, &st, NULL };
    md_bitset_t bf = {0};

    // 64 lipids of 50 united atoms: 2 x 15 in the tails, 20 in the head
    ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("lipid()"), &fctx, alloc).success);
    EXPECT_EQ(64 * 50, (int)md_bitset_count(bf));
    ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("tail()"), &fctx, alloc).success);
    EXPECT_EQ(64 * 30, (int)md_bitset_count(bf));
    ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("head()"), &fctx, alloc).success);
    EXPECT_EQ(64 * 20, (int)md_bitset_count(bf));
    ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("tail() and head()"), &fctx, alloc).success);
    EXPECT_EQ(0, (int)md_bitset_count(bf));
    // The phosphorus of the heads
    ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("head() and element('P')"), &fctx, alloc).success);
    EXPECT_EQ(64, (int)md_bitset_count(bf));
    // One selection per lipid molecule, of all of its atoms
    md_array(md_bitset_t) lipids = 0;
    ASSERT_TRUE(md_filter_sets(&lipids, MD_STR_LIT("lipid()"), &fctx, alloc).success);
    ASSERT_EQ(64, (int)md_array_size(lipids));
    for (size_t i = 0; i < md_array_size(lipids); ++i) EXPECT_EQ(50, (int)md_bitset_count(lipids[i]));
    // Within a context, the lipids with an atom in it, whole
    md_array(md_bitset_t) near = 0;
    ASSERT_TRUE(md_filter_sets(&near, MD_STR_LIT("lipid() in atom(1:3)"), &fctx, alloc).success);
    ASSERT_EQ(1, (int)md_array_size(near));
    EXPECT_EQ(50, (int)md_bitset_count(near[0]));

    // One selection per tail, 15 carbons each
    md_array(md_bitset_t) sets = 0;
    ASSERT_TRUE(md_filter_sets(&sets, MD_STR_LIT("tail()"), &fctx, alloc).success);
    ASSERT_EQ(128, (int)md_array_size(sets));
    for (size_t i = 0; i < md_array_size(sets); ++i) EXPECT_EQ(15, (int)md_bitset_count(sets[i]));
    // Written without parentheses, as the default representations write them
    ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("head or tail"), &fctx, alloc).success);
    EXPECT_EQ(64 * 50, (int)md_bitset_count(bf));
    // Within a context, the tails fully within it: those of the first lipid
    ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("tail() in residue(1)"), &fctx, alloc).success);
    EXPECT_EQ(30, (int)md_bitset_count(bf));

    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, molecule_of_several_residues) {
    // POPC split into residues as AMBER Lipid21 splits it: the oleoyl (OL) and the palmitoyl (PA) from their carbonyls,
    // and the glycerol phosphocholine (PC), here in two runs of atoms as the SMILES orders them. One lipid molecule of
    // the four residues, whatever they are named (PA is a lipid by name, OL and PC are not).
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    const char* names[] = { "OL", "PC", "PA", "PC" };
    const int sizes[]   = { 19, 4, 17, 12 };
    ASSERT_TRUE(lipid_system_from_smiles_residues(&sys, names, sizes, 4, POPC, true));
    ASSERT_EQ(1, (int)sys.lipid.count);
    ASSERT_EQ(4, (int)lipid_num_residues(&sys, 0));
    for (int k = 0; k < 4; ++k) EXPECT_EQ(k, (int)lipid_res(&sys, 0, k));
    EXPECT_EQ(2, (int)sys.lipid.tail.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    for (size_t c = 0; c < 4; ++c) EXPECT_EQ(MD_RESIDUE_KIND_LIPID, md_system_residue_kind(&sys, c));

    md_system_state_t st = { .alloc = alloc };
    md_system_state_init(&st, sys.atom.count);
    const md_filter_context_t fctx = { &sys, &st, NULL };
    md_array(md_bitset_t) sets = 0;
    ASSERT_TRUE(md_filter_sets(&sets, MD_STR_LIT("lipid()"), &fctx, alloc).success);
    ASSERT_EQ(1, (int)md_array_size(sets));
    EXPECT_EQ((int)sys.atom.count, (int)md_bitset_count(sets[0]));
    md_arena_allocator_destroy(alloc);
}

UTEST(lipid, by_name_without_bonds) {
    // A lipid by name in a system without bonds: a lipid molecule (its residue), without heads and tails
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    const char* name = "POPC";
    ASSERT_TRUE(lipid_system_from_smiles_residues(&sys, &name, NULL, 0, POPC, false));
    ASSERT_EQ(1, (int)sys.lipid.count);
    EXPECT_EQ(0, (int)sys.lipid.head.count);
    EXPECT_EQ(0, (int)sys.lipid.tail.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));

    md_system_state_t st = { .alloc = alloc };
    md_system_state_init(&st, sys.atom.count);
    const md_filter_context_t fctx = { &sys, &st, NULL };
    md_bitset_t bf = {0};
    ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("lipid()"), &fctx, alloc).success);
    EXPECT_EQ((int)sys.atom.count, (int)md_bitset_count(bf));
    ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("head() or tail()"), &fctx, alloc).success);
    EXPECT_EQ(0, (int)md_bitset_count(bf));
    md_arena_allocator_destroy(alloc);
}

// ### MARTINI LIPIDS OF EVERY NAMING ###

// A system of Martini molecules from their bead names alone, as a GRO file has them: one residue per entry (its name,
// then its bead names), the beads 4 A apart along a zigzag
typedef struct { const char* name; const char* beads; } martini_residue_t;

static bool martini_system_from_beads(md_system_t* sys, md_system_state_t* st, const martini_residue_t* res, size_t num_res) {
    char gro[16384];
    char body[16000];
    int len = 0, num_atoms = 0;
    for (size_t r = 0; r < num_res; ++r) {
        char beads[512];
        snprintf(beads, sizeof(beads), "%s", res[r].beads);
        for (char* b = strtok(beads, " "); b; b = strtok(NULL, " ")) {
            len += snprintf(body + len, sizeof(body) - len, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", (int)r + 1, res[r].name, b, num_atoms + 1,
                            1.0 + 0.4 * num_atoms, 1.0 + 0.2 * (num_atoms & 1), 1.0 + 2.0 * r);
            num_atoms += 1;
        }
    }
    const int n = snprintf(gro, sizeof(gro), "beads\n%5d\n%s%10.5f%10.5f%10.5f\n", num_atoms, body, 60.0, 60.0, 60.0);
    return md_gro_system_init_from_str(sys, st, (md_str_t){gro, (size_t)n}, NULL) && md_util_system_infer(sys, st, MD_UTIL_INFER_ALL, NULL);
}

// The lipid whose first residue is res, -1 when none
static int lipid_of_res(const md_system_t* sys, size_t res) {
    for (size_t l = 0; l < sys->lipid.count; ++l) {
        if (lipid_res(sys, l, 0) == (int32_t)res) return (int)l;
    }
    return -1;
}

// The beads of the head of lipid l, the names in order, space separated
static void lipid_head_names(char* buf, size_t cap, const md_system_t* sys, int l) {
    buf[0] = '\0';
    size_t len = 0;
    for (uint32_t h = sys->lipid.head_offset[l]; h < sys->lipid.head_offset[l + 1]; ++h) {
        const int32_t* atoms = md_lipid_head_atoms(&sys->lipid, h);
        for (size_t k = 0; k < md_lipid_head_atom_count(&sys->lipid, h); ++k) {
            const md_str_t s = md_atom_name(&sys->atom, atoms[k]);
            len += snprintf(buf + len, cap - len, "%s" MD_STR_FMT, len ? " " : "", MD_STR_ARG(s));
        }
    }
}

// Martini molecules renamed, as the three letters of a PDB file make cholesterol CHL, PSM SPM and PCER CER: told by
// their beads, which give them their bonds, kinds and anatomy, in a Martini system (here the POPC tells)
UTEST(lipid, martini_renamed_by_beads) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    static const martini_residue_t res[] = {
        {"POPC", "NC3 PO4 GL1 GL2 C1A D2A C3A C4A C1B C2B C3B C4B"},
        {"CHL",  "ROH R1 R2 R3 R4 R5 R6 C1 C2"},
        {"SPM",  "NC3 PO4 OH1 AM2 T1A C2A C3A C4A C1B C2B C3B C4B"},
        {"CER",  "COH OH1 AM2 T1A C2A C3A C4A C1B C2B C3B C4B"},
    };
    ASSERT_TRUE(martini_system_from_beads(&sys, &st, res, MD_COUNTOF(res)));
    ASSERT_EQ(4u, sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    char head[256];

    const int chl = lipid_of_res(&sys, 1);
    ASSERT_TRUE(chl >= 0);
    EXPECT_EQ(1u, lipid_num_tails(&sys, chl));
    EXPECT_EQ(MD_LIPID_TAIL_SKELETON, sys.lipid.tail.kind[sys.lipid.tail_offset[chl]]);
    lipid_head_names(head, sizeof(head), &sys, chl);
    EXPECT_STREQ("ROH", head);
    // Its beads are those of Martini 3 cholesterol: no carbons C1 and C2, the mass of the rings on R1 and R2 and none on
    // the virtual sites of the rest of them
    for (uint32_t i = 12; i < 21; ++i) EXPECT_EQ(0, (int)md_atom_atomic_number(&sys.atom, i));
    EXPECT_NEAR(159.8f, md_atom_mass(&sys.atom, 13), 1e-3f);
    EXPECT_NEAR(0.0f, md_atom_mass(&sys.atom, 18), 1e-3f);
    EXPECT_NEAR(72.0f, md_atom_mass(&sys.atom, 20), 1e-3f);

    const int spm = lipid_of_res(&sys, 2);
    ASSERT_TRUE(spm >= 0);
    EXPECT_EQ(2u, lipid_num_tails(&sys, spm));
    lipid_head_names(head, sizeof(head), &sys, spm);
    EXPECT_STREQ("NC3 PO4 OH1 AM2", head);

    const int cer = lipid_of_res(&sys, 3);
    ASSERT_TRUE(cer >= 0);
    EXPECT_EQ(2u, lipid_num_tails(&sys, cer));
    lipid_head_names(head, sizeof(head), &sys, cer);
    EXPECT_STREQ("COH OH1 AM2", head);
    for (uint32_t t = sys.lipid.tail_offset[cer]; t < sys.lipid.tail_offset[cer + 1]; ++t) EXPECT_EQ(4u, tail_spine_len(&sys, t));

    md_vm_arena_destroy(alloc);
}

// The same through a topology: the cholesterol of martini3_membrane.tpr renamed CHL is still a sterol
UTEST(lipid, martini_renamed_tpr) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_tpr_data_t data = {0};
    ASSERT_TRUE(md_tpr_data_parse_file(&data, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.tpr"), alloc));
    size_t renamed = 0;
    for (size_t t = 0; t < data.num_moltypes; ++t) {
        for (size_t r = 0; r < data.moltypes[t].num_residues; ++r) {
            if (md_str_eq_cstr(data.moltypes[t].residues[r].name, "CHOL")) {
                data.moltypes[t].residues[r].name = MD_STR_LIT("CHL");
                renamed += 1;
            }
        }
    }
    ASSERT_GT(renamed, (size_t)0);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_tpr_system_init_from_data(&sys, &st, &data));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT, NULL));
    size_t sterols = 0;
    for (size_t l = 0; l < sys.lipid.count; ++l) {
        if (!md_str_eq_cstr(md_system_residue_name(&sys, (size_t)lipid_res(&sys, l, 0)), "CHL")) continue;
        const uint32_t t = sys.lipid.tail_offset[l];
        EXPECT_EQ(1u, lipid_num_tails(&sys, l));
        EXPECT_EQ(MD_LIPID_TAIL_SKELETON, sys.lipid.tail.kind[t]);
        sterols += 1;
    }
    EXPECT_EQ(16u, sterols);
    EXPECT_EQ(32u + 16u + 8u + 16u, sys.lipid.count);
    md_arena_allocator_destroy(alloc);
}

// The beads of residue res bonded to each other
static size_t res_num_bonds(const md_system_t* sys, size_t res) {
    const md_urange_t r = md_system_residue_atom_range(sys, res);
    size_t n = 0;
    for (uint32_t i = r.beg; i < r.end; ++i) {
        md_bond_iter_t it = md_bond_iter(&sys->bond, i);
        while (md_bond_iter_has_next(&it)) {
            const md_atom_idx_t j = md_bond_iter_atom_index(&it);
            n += (uint32_t)j > i && (uint32_t)j < r.end;
            md_bond_iter_next(&it);
        }
    }
    return n;
}

// The optimal LINCS cholesterol of Martini 2 (martini_v2.0_CHOL_02-optLINCS.itp): the mass of its rings on four dummies
// which do not interact (VOH V2 V3 V1), its ring beads virtual sites of them, so that they hold together through the
// dummies alone. Beads of the masses and bonds of its topology, the hydroxyl bead for a head and a skeleton of the
// dummies and ring beads out to the end of the side chain, as the cholesterol of 8 beads beside it has.
UTEST(lipid, martini2_cholesterol_optlincs) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    static const martini_residue_t res[] = {
        {"DPPC", "NC3 PO4 GL1 GL2 C1A C2A C3A C4A C1B C2B C3B C4B"},
        {"CHOL", "VOH V2 V3 V1 ROH R1 R2 R3 R4 R5 C1 C2"},
        {"CHOL", "ROH R1 R2 R3 R4 R5 C1 C2"},
    };
    ASSERT_TRUE(martini_system_from_beads(&sys, &st, res, MD_COUNTOF(res)));
    ASSERT_EQ(3u, sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));

    const md_urange_t r = md_system_residue_atom_range(&sys, 1);
    ASSERT_EQ(12u, r.end - r.beg);
    for (uint32_t i = r.beg; i < r.end; ++i) {
        EXPECT_EQ(0, (int)md_atom_atomic_number(&sys.atom, i));
        EXPECT_EQ(MD_PARTICLE_BEAD, md_atom_particle_kind(&sys.atom, i));
    }
    EXPECT_NEAR(125.98f, md_atom_mass(&sys.atom, r.beg + 0), 1e-3f);    // VOH
    EXPECT_NEAR(98.82f,  md_atom_mass(&sys.atom, r.beg + 3), 1e-3f);    // V1
    EXPECT_NEAR(0.0f,    md_atom_mass(&sys.atom, r.beg + 4), 1e-3f);    // ROH, a virtual site
    EXPECT_NEAR(72.0f,   md_atom_mass(&sys.atom, r.beg + 11), 1e-3f);   // C2
    // The constraints of the dummies (5), each virtual site to the first particle it is built from (7), C1-C2
    EXPECT_EQ(13u, res_num_bonds(&sys, 1));
    EXPECT_EQ(9u, res_num_bonds(&sys, 2));

    char head[256];
    for (size_t c = 1; c <= 2; ++c) {
        const int l = lipid_of_res(&sys, c);
        ASSERT_TRUE(l >= 0);
        ASSERT_EQ(1u, lipid_num_tails(&sys, l));
        const uint32_t t = sys.lipid.tail_offset[l];
        EXPECT_EQ(MD_LIPID_TAIL_SKELETON, sys.lipid.tail.kind[t]);
        lipid_head_names(head, sizeof(head), &sys, l);
        EXPECT_STREQ("ROH", head);
        EXPECT_EQ(c == 1 ? 11u : 7u, md_lipid_tail_atom_count(&sys.lipid, t));
        const int32_t* spine = md_lipid_tail_spine(&sys.lipid, t);
        EXPECT_TRUE(md_str_eq_cstr(md_atom_name(&sys.atom, spine[tail_spine_len(&sys, t) - 1]), "C2"));
    }
    md_vm_arena_destroy(alloc);
}

// A molecule of the force field whose beads are those of none of its topologies (here a cholesterol of a bead unknown),
// in a Martini system: beads of the kind of the molecule, not the elements of their names (V1 vanadium, C1 carbon), and
// without a topology no bonds, a lipid by name without heads and tails. Atoms named as an atomistic force field names them
// (a cholesterol of CHARMM, whose C1 and C2 are beads of the Martini one too) stay atoms.
UTEST(lipid, martini_unknown_variant_is_beads) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    static const martini_residue_t res[] = {
        {"DPPC", "NC3 PO4 GL1 GL2 C1A C2A C3A C4A C1B C2B C3B C4B"},
        {"CHOL", "VXH V2 V3 V1 ROH R1 R2 R3 R4 R5 C1 C2"},
        {"CHOL", "C1 C2 C3 O3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27"},
    };
    ASSERT_TRUE(martini_system_from_beads(&sys, &st, res, MD_COUNTOF(res)));
    EXPECT_TRUE(lipid_invariants_hold(&sys));

    const md_urange_t r = md_system_residue_atom_range(&sys, 1);
    for (uint32_t i = r.beg; i < r.end; ++i) {
        EXPECT_EQ(0, (int)md_atom_atomic_number(&sys.atom, i));
        EXPECT_EQ(MD_PARTICLE_BEAD, md_atom_particle_kind(&sys.atom, i));
        EXPECT_NEAR(72.0f, md_atom_mass(&sys.atom, i), 1e-3f);
    }
    EXPECT_EQ(MD_RESIDUE_KIND_LIPID, md_system_residue_kind(&sys, 1));
    EXPECT_EQ(0u, res_num_bonds(&sys, 1));
    const int l = lipid_of_res(&sys, 1);
    ASSERT_TRUE(l >= 0);
    EXPECT_EQ(0u, lipid_num_heads(&sys, l));
    EXPECT_EQ(0u, lipid_num_tails(&sys, l));

    const md_urange_t a = md_system_residue_atom_range(&sys, 2);
    EXPECT_EQ(MD_Z_C, md_atom_atomic_number(&sys.atom, a.beg));
    EXPECT_EQ(MD_Z_O, md_atom_atomic_number(&sys.atom, a.beg + 3));
    md_vm_arena_destroy(alloc);
}

// Martini 2 molecules whose chains are named without the letter of a tail: fatty acids, detergents, the lyso PC whose
// tail is one bead. The inositol of a phosphatidylinositol is named so too (C1 C2 C3) and stays in its head: a ring.
// The cardiolipins of the Martini 2 lipidome name the first glycerol bead of their second phosphatidyl GL21, as the
// second of the first is named (an erratum of the itp); CHARMM-GUI names it GL12, as the scheme GLxy (bead x of
// phosphatidyl y) has it. Of either, the beads and bonds of the topology: a head of the glycerol and the phosphatidyls,
// and the four tails each off its glycerol bead.
UTEST(lipid, martini2_cardiolipin_of_either_name) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    static const martini_residue_t res[] = {
        {"POPE", "NH3 PO4 GL1 GL2 C1A D2A C3A C4A C1B C2B C3B C4B"},
        {"CDL2", "GL0 PO41 GL11 GL21 C1A1 C2A1 D3A1 C4A1 C5A1 C1B1 C2B1 D3B1 C4B1 C5B1 PO42 GL12 GL22 C1A2 C2A2 D3A2 C4A2 C5A2 C1B2 C2B2 D3B2 C4B2 C5B2"},
        {"CDL2", "GL0 PO41 GL11 GL21 C1A1 C2A1 D3A1 C4A1 C5A1 C1B1 C2B1 D3B1 C4B1 C5B1 PO42 GL21 GL22 C1A2 C2A2 D3A2 C4A2 C5A2 C1B2 C2B2 D3B2 C4B2 C5B2"},
    };
    ASSERT_TRUE(martini_system_from_beads(&sys, &st, res, MD_COUNTOF(res)));
    ASSERT_EQ(3, (int)sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    char head[256];
    for (size_t r = 1; r <= 2; ++r) {
        const md_urange_t range = md_system_residue_atom_range(&sys, r);
        size_t num_bonds = 0;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            EXPECT_EQ(MD_PARTICLE_BEAD, md_atom_particle_kind(&sys.atom, i));
            EXPECT_NEAR(72.0f, md_atom_mass(&sys.atom, i), 1e-3f);
            md_bond_iter_t it = md_bond_iter(&sys.bond, i);
            while (md_bond_iter_has_next(&it)) {
                const md_atom_idx_t j = md_bond_iter_atom_index(&it);
                num_bonds += (uint32_t)j > i && (uint32_t)j < range.end;
                md_bond_iter_next(&it);
            }
        }
        EXPECT_EQ(26, (int)num_bonds);
        const int l = lipid_of_res(&sys, r);
        ASSERT_TRUE(l >= 0);
        ASSERT_EQ(1, (int)lipid_num_heads(&sys, l));
        lipid_head_names(head, sizeof(head), &sys, l);
        EXPECT_STREQ(r == 1 ? "GL0 PO41 GL11 GL21 PO42 GL12 GL22" : "GL0 PO41 GL11 GL21 PO42 GL21 GL22", head);
        ASSERT_EQ(4, (int)lipid_num_tails(&sys, l));
        for (uint32_t t = sys.lipid.tail_offset[l]; t < sys.lipid.tail_offset[l + 1]; ++t) {
            EXPECT_EQ(5, (int)tail_spine_len(&sys, t));
            EXPECT_EQ(MD_LIPID_TAIL_CHAIN, sys.lipid.tail.kind[t]);
        }
    }
    md_vm_arena_destroy(alloc);
}

UTEST(lipid, martini2_letterless_chains) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    static const martini_residue_t res[] = {
        {"DPPC", "NC3 PO4 GL1 GL2 C1A C2A C3A C4A C1B C2B C3B C4B"},
        {"PCA",  "COO C1 C2 C3 C4"},                            // Palmitic acid
        {"DDM",  "B1 B2 B3 B4 B5 B6 C1 C2 C3"},                 // Dodecyl maltoside
        {"SDS",  "SO3 C1 C2 C4"},
        {"CPC",  "NC3 PO4 GL1 GL2 C1B"},                        // A lyso PC of a short tail
        {"POPI", "C1 C2 C3 PO4 GL1 GL2 C1A D2A C3A C4A C1B C2B C3B C4B"},
    };
    ASSERT_TRUE(martini_system_from_beads(&sys, &st, res, MD_COUNTOF(res)));
    ASSERT_EQ(6u, sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));

    static const struct { const char* head; size_t tails; size_t spine; } expected[] = {
        {"NC3 PO4 GL1 GL2", 2, 4},
        {"COO", 1, 4},
        {"B1 B2 B3 B4 B5 B6", 1, 3},
        {"SO3", 1, 3},
        {"NC3 PO4 GL1 GL2", 1, 1},
        {"C1 C2 C3 PO4 GL1 GL2", 2, 4},
    };
    char head[256];
    for (size_t c = 0; c < MD_COUNTOF(expected); ++c) {
        const int l = lipid_of_res(&sys, c);
        ASSERT_TRUE(l >= 0);
        lipid_head_names(head, sizeof(head), &sys, l);
        EXPECT_STREQ(expected[c].head, head);
        EXPECT_EQ(expected[c].tails, lipid_num_tails(&sys, l));
        for (uint32_t t = sys.lipid.tail_offset[l]; t < sys.lipid.tail_offset[l + 1]; ++t) {
            EXPECT_EQ(MD_LIPID_TAIL_CHAIN, sys.lipid.tail.kind[t]);
            EXPECT_EQ(expected[c].spine, tail_spine_len(&sys, t));
        }
    }
    md_vm_arena_destroy(alloc);
}

// The atoms of the lipid() selections, and whether the residues are of the kind
static size_t lipid_selection_atoms(md_system_t* sys, md_allocator_i* alloc, int* num_sets) {
    md_system_state_t st = { .alloc = alloc };
    md_system_state_init(&st, sys->atom.count);
    const md_filter_context_t fctx = { sys, &st, NULL };
    md_array(md_bitset_t) sets = 0;
    if (!md_filter_sets(&sets, MD_STR_LIT("lipid()"), &fctx, alloc).success) return 0;
    *num_sets = (int)md_array_size(sets);
    size_t n = 0;
    for (size_t i = 0; i < md_array_size(sets); ++i) n += md_bitset_count(sets[i]);
    return n;
}

UTEST(lipid, glycolipid_of_several_residues) {
    // Lactosylceramide as CHARMM-GUI writes a glycolipid: the ceramide (d18:1/16:0) a residue, each sugar one of its own.
    // The lipid is the ceramide and the glucose on it, its head group; the galactose beyond is attached to the lipid but no
    // part of it, nor of its head, and is no lipid.
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    const char* names[] = { "CER", "BGLC", "BGAL" };
    const int sizes[]   = { 38, 11, 11 };
    ASSERT_TRUE(lipid_system_from_smiles_residues(&sys, names, sizes, 3,
        "CCCCCCCCCCCCCCCC(=O)NC(C(O)C=CCCCCCCCCCCCCC)CO" "C1OC(CO)C(C(O)C1O)O" "C2OC(CO)C(O)C(O)C2O", true));
    ASSERT_EQ(1, (int)sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    ASSERT_EQ(2, (int)lipid_num_residues(&sys, 0));
    EXPECT_EQ(0, lipid_res(&sys, 0, 0));
    EXPECT_EQ(1, lipid_res(&sys, 0, 1));
    EXPECT_EQ(MD_RESIDUE_KIND_LIPID, md_system_residue_kind(&sys, 0));
    EXPECT_EQ(MD_RESIDUE_KIND_LIPID, md_system_residue_kind(&sys, 1));
    EXPECT_EQ(MD_RESIDUE_KIND_OTHER, md_system_residue_kind(&sys, 2));
    // Two tails: the palmitoyl off its carbonyl, the sphingoid chain off C3. One head: the polar part of the ceramide
    // (C1 O1 C2 N C3 O3, the carbonyl and its oxygen) and the glucose, 11 of its atoms.
    EXPECT_EQ(2, (int)lipid_num_tails(&sys, 0));
    ASSERT_EQ(1, (int)lipid_num_heads(&sys, 0));
    EXPECT_EQ(8 + 11, (int)md_lipid_head_atom_count(&sys.lipid, 0));
    int num_sets = 0;
    EXPECT_EQ(38u + 11u, lipid_selection_atoms(&sys, alloc, &num_sets));
    EXPECT_EQ(1, num_sets);

    // The same again
    ASSERT_TRUE(md_util_system_infer_lipids(&sys));
    ASSERT_EQ(1, (int)sys.lipid.count);
    EXPECT_EQ(2, (int)lipid_num_residues(&sys, 0));
    EXPECT_EQ(MD_RESIDUE_KIND_OTHER, md_system_residue_kind(&sys, 2));
    md_vm_arena_destroy(alloc);
}

UTEST(lipid, molecule_of_residues_apart) {
    // A lipid whose residues are not consecutive, as a PDB entry lists a lipid A: its sugars in a chain of their own and
    // its acyl chains as ligands after the rest. Here a monopalmitoyl glycerol, its glycerol (GOL) and palmitoyl (PLM) with
    // a water between them. One lipid of the two, and the water no part of it.
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    const char* names[] = { "GOL", "HOH", "PLM" };
    const int sizes[]   = { 6, 1, 17 };
    ASSERT_TRUE(lipid_system_from_smiles_residues(&sys, names, sizes, 3, "OCC(O)CO%10.O.C%10(=O)CCCCCCCCCCCCCCC", true));
    ASSERT_EQ(1, (int)sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    ASSERT_EQ(2, (int)lipid_num_residues(&sys, 0));
    EXPECT_EQ(0, lipid_res(&sys, 0, 0));
    EXPECT_EQ(2, lipid_res(&sys, 0, 1));
    EXPECT_NE(MD_RESIDUE_KIND_LIPID, md_system_residue_kind(&sys, 1));
    ASSERT_EQ(1, (int)lipid_num_tails(&sys, 0));
    EXPECT_EQ(15, (int)tail_spine_len(&sys, 0));
    ASSERT_EQ(1, (int)lipid_num_heads(&sys, 0));
    // The glycerol and the carbonyl
    EXPECT_EQ(6 + 2, (int)md_lipid_head_atom_count(&sys.lipid, 0));
    int num_sets = 0;
    EXPECT_EQ(6u + 17u, lipid_selection_atoms(&sys, alloc, &num_sets));
    EXPECT_EQ(1, num_sets);
    md_vm_arena_destroy(alloc);
}

// The atom of the name in residue res, -1 when none
static int atom_in_res(const md_system_t* sys, size_t res, const char* name) {
    const md_urange_t r = md_system_residue_atom_range(sys, res);
    for (uint32_t i = r.beg; i < r.end; ++i) {
        if (md_str_eq_cstr(md_atom_name(&sys->atom, i), name)) return (int)i;
    }
    return -1;
}

// The lipids of AMBER (Lipid14, 17 and 21) are split into residues, each acyl chain one (OL oleoyl, PA palmitoyl, ...)
// whose C12 is bonded to the C11 or C21 of an ester of the head (PE, PC, PGR, ...), the residue before or after it: a
// lipid by its structure, a molecule of the three. Here a DOPE of a packed starting structure, whose head another one
// (PE170) overlaps: the C12 of OL142 is 1.54 A from the C11 of its own head, and 1.42 A from O22 and 1.53 A from C21 of
// the other, more bonds than a carbon takes. The bonds to the residue apart go, not the ester to the residue next to
// it, and the chain is a tail of its lipid rather than a molecule of its own.
UTEST(lipid, amber_tail_kept_on_its_head_over_a_clash) {
    static const char gro[] =
        "OL142 PE143 OL144, a DOPE of AMBER Lipid21, and PE170, whose head overlaps that of PE143\n"
        "  158\n"
        "  142OL     C12    1   3.062   6.995   2.513\n"
        "  142OL     H2R    2   3.113   6.898   2.535\n"
        "  142OL     H2S    3   3.086   7.026   2.409\n"
        "  142OL     C13    4   3.103   7.107   2.613\n"
        "  142OL     H3R    5   3.214   7.107   2.618\n"
        "  142OL     H3S    6   3.059   7.085   2.713\n"
        "  142OL     C14    7   3.050   7.248   2.593\n"
        "  142OL     H4R    8   3.077   7.279   2.490\n"
        "  142OL     H4S    9   2.939   7.243   2.596\n"
        "  142OL     C15   10   3.071   7.352   2.697\n"
        "  142OL     H5R   11   3.007   7.441   2.684\n"
        "  142OL     H5S   12   3.176   7.385   2.678\n"
        "  142OL     C16   13   3.051   7.309   2.847\n"
        "  142OL     H6R   14   3.113   7.219   2.866\n"
        "  142OL     H6S   15   2.941   7.292   2.845\n"
        "  142OL     C17   16   3.087   7.411   2.955\n"
        "  142OL     H7R   17   3.103   7.370   3.057\n"
        "  142OL     H7S   18   3.183   7.459   2.925\n"
        "  142OL     C18   19   2.976   7.515   2.969\n"
        "  142OL     H8R   20   2.886   7.453   2.987\n"
        "  142OL     H8S   21   2.955   7.569   2.874\n"
        "  142OL     C19   22   3.013   7.617   3.077\n"
        "  142OL     H9R   23   3.079   7.698   3.042\n"
        "  142OL    C110   24   2.972   7.619   3.205\n"
        "  142OL    H10R   25   3.006   7.701   3.271\n"
        "  142OL    C111   26   2.880   7.510   3.259\n"
        "  142OL    H11R   27   2.833   7.447   3.180\n"
        "  142OL    H11S   28   2.938   7.430   3.310\n"
        "  142OL    C112   29   2.788   7.565   3.376\n"
        "  142OL    H12R   30   2.693   7.508   3.388\n"
        "  142OL    H12S   31   2.750   7.668   3.361\n"
        "  142OL    C113   32   2.862   7.553   3.507\n"
        "  142OL    H13R   33   2.894   7.448   3.524\n"
        "  142OL    H13S   34   2.952   7.618   3.507\n"
        "  142OL    C114   35   2.771   7.582   3.633\n"
        "  142OL    H14R   36   2.683   7.515   3.624\n"
        "  142OL    H14S   37   2.728   7.685   3.631\n"
        "  142OL    C115   38   2.839   7.547   3.770\n"
        "  142OL    H15R   39   2.949   7.549   3.755\n"
        "  142OL    H15S   40   2.818   7.443   3.806\n"
        "  142OL    C116   41   2.791   7.648   3.883\n"
        "  142OL    H16R   42   2.837   7.744   3.850\n"
        "  142OL    H16S   43   2.829   7.615   3.982\n"
        "  142OL    C117   44   2.642   7.645   3.889\n"
        "  142OL    H17R   45   2.591   7.680   3.797\n"
        "  142OL    H17S   46   2.607   7.539   3.885\n"
        "  142OL    C118   47   2.599   7.722   4.014\n"
        "  142OL    H18R   48   2.646   7.823   4.016\n"
        "  142OL    H18S   49   2.489   7.738   4.010\n"
        "  142OL    H18T   50   2.620   7.668   4.109\n"
        "  143PE     C11   51   2.911   6.969   2.530\n"
        "  143PE     O12   52   2.843   6.985   2.628\n"
        "  143PE     O11   53   2.861   6.919   2.414\n"
        "  143PE      C1   54   2.724   6.882   2.417\n"
        "  143PE      HR   55   2.722   6.799   2.342\n"
        "  143PE      HS   56   2.680   6.850   2.514\n"
        "  143PE      C2   57   2.635   7.001   2.351\n"
        "  143PE      HX   58   2.529   6.968   2.344\n"
        "  143PE      C3   59   2.680   7.051   2.210\n"
        "  143PE      HA   60   2.782   7.088   2.233\n"
        "  143PE      HB   61   2.679   6.966   2.139\n"
        "  143PE     O31   62   2.601   7.154   2.158\n"
        "  143PE     P31   63   2.611   7.188   2.006\n"
        "  143PE     O32   64   2.758   7.247   1.998\n"
        "  143PE     C31   65   2.788   7.365   2.059\n"
        "  143PE     H1A   66   2.732   7.447   2.009\n"
        "  143PE     H1B   67   2.764   7.367   2.167\n"
        "  143PE     C32   68   2.938   7.402   2.038\n"
        "  143PE     H2A   69   2.995   7.321   2.089\n"
        "  143PE     H2B   70   2.958   7.400   1.928\n"
        "  143PE     N31   71   2.974   7.541   2.094\n"
        "  143PE    HN1A   72   2.892   7.587   2.138\n"
        "  143PE    HN1B   73   3.051   7.530   2.163\n"
        "  143PE    HN1C   74   3.007   7.605   2.019\n"
        "  143PE     O33   75   2.607   7.062   1.927\n"
        "  143PE     O34   76   2.505   7.282   1.973\n"
        "  143PE     O21   77   2.635   7.111   2.441\n"
        "  143PE     C21   78   2.559   7.101   2.550\n"
        "  143PE     O22   79   2.477   7.014   2.559\n"
        "  144OL     C12   80   2.571   7.229   2.630\n"
        "  144OL     H2R   81   2.556   7.318   2.566\n"
        "  144OL     H2S   82   2.672   7.227   2.679\n"
        "  144OL     C13   83   2.473   7.238   2.746\n"
        "  144OL     H3R   84   2.464   7.141   2.801\n"
        "  144OL     H3S   85   2.374   7.259   2.699\n"
        "  144OL     C14   86   2.522   7.326   2.865\n"
        "  144OL     H4R   87   2.608   7.287   2.924\n"
        "  144OL     H4S   88   2.438   7.346   2.935\n"
        "  144OL     C15   89   2.566   7.466   2.822\n"
        "  144OL     H5R   90   2.475   7.503   2.768\n"
        "  144OL     H5S   91   2.658   7.451   2.762\n"
        "  144OL     C16   92   2.605   7.552   2.943\n"
        "  144OL     H6R   93   2.700   7.514   2.986\n"
        "  144OL     H6S   94   2.538   7.536   3.030\n"
        "  144OL     C17   95   2.610   7.700   2.923\n"
        "  144OL     H7R   96   2.517   7.744   2.882\n"
        "  144OL     H7S   97   2.691   7.710   2.847\n"
        "  144OL     C18   98   2.638   7.781   3.048\n"
        "  144OL     H8R   99   2.675   7.711   3.125\n"
        "  144OL     H8S  100   2.545   7.821   3.095\n"
        "  144OL     C19  101   2.748   7.882   3.031\n"
        "  144OL     H9R  102   2.831   7.850   2.966\n"
        "  144OL    C110  103   2.748   8.010   3.083\n"
        "  144OL    H10R  104   2.841   8.070   3.091\n"
        "  144OL    C111  105   2.638   8.067   3.170\n"
        "  144OL    H11R  106   2.540   8.087   3.122\n"
        "  144OL    H11S  107   2.663   8.172   3.196\n"
        "  144OL    C112  108   2.626   7.996   3.312\n"
        "  144OL    H12R  109   2.632   7.886   3.301\n"
        "  144OL    H12S  110   2.530   8.035   3.352\n"
        "  144OL    C113  111   2.735   8.043   3.402\n"
        "  144OL    H13R  112   2.715   8.143   3.446\n"
        "  144OL    H13S  113   2.832   8.059   3.350\n"
        "  144OL    C114  114   2.755   7.949   3.523\n"
        "  144OL    H14R  115   2.714   7.848   3.501\n"
        "  144OL    H14S  116   2.693   7.984   3.609\n"
        "  144OL    C115  117   2.897   7.933   3.569\n"
        "  144OL    H15R  118   2.959   7.872   3.500\n"
        "  144OL    H15S  119   2.888   7.861   3.653\n"
        "  144OL    C116  120   2.964   8.059   3.617\n"
        "  144OL    H16R  121   2.898   8.115   3.687\n"
        "  144OL    H16S  122   2.990   8.129   3.535\n"
        "  144OL    C117  123   3.103   8.031   3.689\n"
        "  144OL    H17R  124   3.163   8.124   3.703\n"
        "  144OL    H17S  125   3.167   7.970   3.622\n"
        "  144OL    C118  126   3.099   7.968   3.827\n"
        "  144OL    H18R  127   3.057   8.041   3.900\n"
        "  144OL    H18S  128   3.206   7.955   3.854\n"
        "  144OL    H18T  129   3.040   7.874   3.831\n"
        "  170PE     C11  130   3.329   7.324   2.254\n"
        "  170PE     O12  131   3.341   7.247   2.346\n"
        "  170PE     O11  132   3.255   7.297   2.144\n"
        "  170PE      C1  133   3.194   7.169   2.142\n"
        "  170PE      HR  134   3.189   7.148   2.033\n"
        "  170PE      HS  135   3.243   7.086   2.197\n"
        "  170PE      C2  136   3.039   7.182   2.189\n"
        "  170PE      HX  137   2.988   7.084   2.177\n"
        "  170PE      C3  138   2.953   7.290   2.114\n"
        "  170PE      HA  139   3.008   7.382   2.141\n"
        "  170PE      HB  140   2.952   7.268   2.006\n"
        "  170PE     O31  141   2.820   7.297   2.159\n"
        "  170PE     P31  142   2.715   7.374   2.072\n"
        "  170PE     O32  143   2.768   7.522   2.084\n"
        "  170PE     C31  144   2.762   7.588   2.203\n"
        "  170PE     H1A  145   2.655   7.602   2.232\n"
        "  170PE     H1B  146   2.812   7.534   2.286\n"
        "  170PE     C32  147   2.819   7.732   2.190\n"
        "  170PE     H2A  148   2.926   7.718   2.164\n"
        "  170PE     H2B  149   2.765   7.783   2.108\n"
        "  170PE     N31  150   2.803   7.814   2.320\n"
        "  170PE    HN1A  151   2.758   7.760   2.395\n"
        "  170PE    HN1B  152   2.897   7.847   2.350\n"
        "  170PE    HN1C  153   2.743   7.898   2.305\n"
        "  170PE     O33  154   2.729   7.330   1.930\n"
        "  170PE     O34  155   2.583   7.349   2.127\n"
        "  170PE     O21  156   3.037   7.210   2.328\n"
        "  170PE     C21  157   3.063   7.109   2.411\n"
        "  170PE     O22  158   3.066   6.996   2.371\n"
        "  15.35366  15.35366   8.45366\n";
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &st, md_str_from_cstr(gro), NULL));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL, NULL));
    ASSERT_EQ(4, (int)md_system_residue_count(&sys));
    EXPECT_TRUE(lipid_invariants_hold(&sys));

    const int c12 = atom_in_res(&sys, 0, "C12");
    const int c11 = atom_in_res(&sys, 1, "C11");
    ASSERT_TRUE(c12 >= 0 && c11 >= 0);
    EXPECT_GE(md_bond_find(&sys.bond, c12, c11), 0);
    EXPECT_LT(md_bond_find(&sys.bond, c12, atom_in_res(&sys, 3, "O22")), 0);
    EXPECT_LT(md_bond_find(&sys.bond, c12, atom_in_res(&sys, 3, "C21")), 0);

    // Both chains are of the lipid of PE143 (which the overlap of the heads may join PE170 to), and the chain of OL142 is
    // a tail of it, from its C12 off the C11 of the head
    EXPECT_EQ(MD_RESIDUE_KIND_LIPID, md_system_residue_kind(&sys, 0));
    EXPECT_EQ(MD_RESIDUE_KIND_LIPID, md_system_residue_kind(&sys, 2));
    const int l = lipid_of_res(&sys, 0);
    ASSERT_TRUE(l >= 0);
    EXPECT_EQ(1, lipid_res(&sys, l, 1));
    EXPECT_EQ(2, lipid_res(&sys, l, 2));
    bool tail = false;
    for (uint32_t t = sys.lipid.tail_offset[l]; t < sys.lipid.tail_offset[l + 1]; ++t) {
        tail |= sys.lipid.tail.anchor_beg[t] == c11 && md_lipid_tail_spine(&sys.lipid, t)[0] == c12 && tail_spine_len(&sys, t) == 17;
    }
    EXPECT_TRUE(tail);
    md_vm_arena_destroy(alloc);
}

UTEST(lipid, martini_unknown_by_tail_beads) {
    // A coarse grained glycolipid the Martini database does not have (a ganglioside of a topology of its own), from a GRO
    // and its topology: its sugar beads (GM1 GM2 GM3) and amides (AM1 AM2), its tails named as those of the lipidomes. The
    // topology's masses make the particles beads (C1B was a carbon, AM1 americium); then it is a lipid by its beads.
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    static const martini_residue_t res[] = {
        {"GSL", "GM1 GM2 GM3 AM1 AM2 D1A C2A C3A C4A C1B C2B C3B C4B"},
    };
    ASSERT_TRUE(martini_system_from_beads(&sys, &st, res, MD_COUNTOF(res)));
    EXPECT_EQ(0u, sys.lipid.count);

    const char* itp =
        "[ moleculetype ]\nGSL 1\n[ atoms ]\n"
        "1 P4 1 GSL GM1 1 0 72\n2 P4 1 GSL GM2 2 0 72\n3 P4 1 GSL GM3 3 0 72\n4 P5 1 GSL AM1 4 0 72\n5 P5 1 GSL AM2 5 0 72\n"
        "6 C4 1 GSL D1A 6 0 72\n7 C1 1 GSL C2A 7 0 72\n8 C1 1 GSL C3A 8 0 72\n9 C1 1 GSL C4A 9 0 72\n"
        "10 C1 1 GSL C1B 10 0 72\n11 C1 1 GSL C2B 11 0 72\n12 C1 1 GSL C3B 12 0 72\n13 C1 1 GSL C4B 13 0 72\n"
        "[ bonds ]\n1 2 1 0.3 1e4\n1 3 1 0.3 1e4\n2 3 1 0.3 1e4\n1 4 1 0.4 1e4\n4 5 1 0.4 1e4\n4 6 1 0.47 1250\n6 7 1 0.47 1250\n"
        "7 8 1 0.47 1250\n8 9 1 0.47 1250\n5 10 1 0.47 1250\n10 11 1 0.47 1250\n11 12 1 0.47 1250\n12 13 1 0.47 1250\n";
    md_itp_data_t data = {0};
    ASSERT_TRUE(md_itp_data_parse_str(&data, md_str_from_cstr(itp), MD_STR_LIT(""), alloc));
    ASSERT_TRUE(md_itp_system_supplement(&sys, &data));

    for (size_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(MD_PARTICLE_BEAD, md_atom_particle_kind(&sys.atom, i));
        EXPECT_EQ(0, (int)md_atom_atomic_number(&sys.atom, i));
    }
    ASSERT_EQ(1u, sys.lipid.count);
    EXPECT_TRUE(lipid_invariants_hold(&sys));
    char head[256];
    lipid_head_names(head, sizeof(head), &sys, 0);
    EXPECT_STREQ("GM1 GM2 GM3 AM1 AM2", head);
    ASSERT_EQ(2u, lipid_num_tails(&sys, 0));
    for (uint32_t t = 0; t < 2; ++t) EXPECT_EQ(4u, tail_spine_len(&sys, t));

    // Inferred again, now a lipid by its kind: the same
    ASSERT_TRUE(md_util_system_infer_lipids(&sys));
    ASSERT_EQ(1u, sys.lipid.count);
    EXPECT_EQ(2u, lipid_num_tails(&sys, 0));
    lipid_head_names(head, sizeof(head), &sys, 0);
    EXPECT_STREQ("GM1 GM2 GM3 AM1 AM2", head);
    md_vm_arena_destroy(alloc);
}
