#include "utest.h"

#include <md_match.h>
#include <md_smiles.h>
#include <md_system.h>
#include <md_util.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_mmcif.h>
#include <md_filter.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_os.h>

#include <string.h>

// ### SYSTEMS ###
// Loaded once on first use and kept for the run: centered.gro alone takes a good part of a second to load and infer.

typedef enum {
    SYS_ALA = 0,        // 15 ALA, all atom (MD)
    SYS_PFTAA,          // One PFTAA, all atom
    SYS_CENTERED,       // 253 amyloid chains and 61 PFTAA, all atom (MD)
    SYS_DNA,            // DNA and protein, all atom (MD)
    SYS_1K4R,           // Crystal structure, no hydrogens
    SYS_1FEZ,           // Crystal structure, no hydrogens
    SYS_2OR2,           // Crystal structure, no hydrogens
    SYS_TUBULIN,        // Crystal structure, no hydrogens, GTP and GDP with their Mg
    SYS_TIP4P,          // A peptide in TIP4P water, whose M sites are virtual sites
    NUM_SYS,
} sys_id_t;

static md_system_t       sys_data[NUM_SYS];
static md_system_state_t sys_state[NUM_SYS];

static md_system_t* get_sys(sys_id_t id) {
    static md_allocator_i* alloc = NULL;
    static bool loaded[NUM_SYS];
    if (!alloc) alloc = md_vm_arena_create(GIGABYTES(4));
    if (!loaded[id]) {
        md_system_t* s = &sys_data[id];
        md_system_state_t* st = &sys_state[id];
        s->alloc  = alloc;
        st->alloc = alloc;
        switch (id) {
        case SYS_ALA:      md_pdb_system_init_from_file(s, st, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE); break;
        case SYS_PFTAA:    md_gro_system_init_from_file(s, st, STR_LIT(MD_UNITTEST_DATA_DIR "/pftaa.gro")); break;
        case SYS_CENTERED: md_gro_system_init_from_file(s, st, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro")); break;
        case SYS_DNA:      md_gro_system_init_from_file(s, st, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro")); break;
        case SYS_1K4R:     md_pdb_system_init_from_file(s, st, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE); break;
        case SYS_1FEZ:     md_mmcif_system_init_from_file(s, st, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif")); break;
        case SYS_2OR2:     md_mmcif_system_init_from_file(s, st, STR_LIT(MD_UNITTEST_DATA_DIR "/2or2.cif")); break;
        case SYS_TUBULIN:  md_pdb_system_init_from_file(s, st, STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE); break;
        case SYS_TIP4P:    md_gro_system_init_from_file(s, st, STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip4p.gro")); break;
        default: break;
        }
        md_util_system_infer(s, st, MD_UTIL_INFER_ALL);
        loaded[id] = true;
    }
    return &sys_data[id];
}

static const md_system_state_t* get_state(sys_id_t id) {
    get_sys(id);
    return &sys_state[id];
}

// A view of a system without the hydrogen counts and charges of md_chem_perceive, as a system has which was not
// perceived. Bond orders and aromatic flags stay.
static md_system_t without_chemistry(const md_system_t* sys) {
    md_system_t view = *sys;
    view.atom.hydrogen_count = NULL;
    view.atom.formal_charge  = NULL;
    return view;
}

static size_t count_smiles_ex(const char* smiles, const md_match_desc_t* base, const md_system_t* sys) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_match_query_t query = {0};
    md_smiles_error_t err = {0};
    if (!md_match_query_init_smiles(&query, str_from_cstr(smiles), alloc, &err)) {
        printf("SMILES error in '%s' at %zu: %s\n", smiles, err.offset, err.message);
        return SIZE_MAX;
    }
    md_match_desc_t desc = *base;
    desc.query = &query;
    size_t count = 0;
    const bool ok = md_match_for_each(&count, &desc, sys, NULL, NULL);
    md_match_query_free(&query);
    return ok ? count : SIZE_MAX;
}

static size_t count_smiles(const char* smiles, md_match_level_t level, md_match_mode_t mode, const md_system_t* sys) {
    md_match_desc_t desc = {.level = level, .mode = mode};
    return count_smiles_ex(smiles, &desc, sys);
}

// md_match_for_each as a count, SIZE_MAX for an invalid query or description
static size_t count_matches(const md_match_desc_t* desc, const md_system_t* sys, md_match_callback_t callback, void* user_param) {
    size_t count = 0;
    return md_match_for_each(&count, desc, sys, callback, user_param) ? count : SIZE_MAX;
}

static bool atoms_bonded(const md_system_t* sys, int32_t a, int32_t b) {
    md_bond_iter_t it = md_bond_iter(&sys->bond, a);
    while (md_bond_iter_has_next(&it)) {
        if (md_bond_iter_atom_index(&it) == b) return true;
        md_bond_iter_next(&it);
    }
    return false;
}

// ### SMILES ###

UTEST(match, smiles_parse) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_smiles_t g = {0};
    md_smiles_error_t err = {0};

    ASSERT_TRUE(md_smiles_parse(&g, STR_LIT("C1=CC=CC=C1"), alloc, &err));
    EXPECT_EQ(6, (int)g.num_atoms);
    EXPECT_EQ(6, (int)g.num_bonds);
    EXPECT_EQ(1, (int)g.num_components);
    EXPECT_EQ(MD_SMILES_BOND_DOUBLE,   g.bonds[0].order);
    EXPECT_EQ(MD_SMILES_BOND_IMPLICIT, g.bonds[1].order);
    // The ring closure is the last bond, between the first and the last atom
    EXPECT_EQ(0, (int)g.bonds[5].a);
    EXPECT_EQ(5, (int)g.bonds[5].b);
    EXPECT_TRUE(g.bonds[5].flags & MD_SMILES_BOND_RING);
    md_smiles_free(&g);

    // The bond of a ring closure, on either side, belongs to the ring bond and not to the next atom
    const char* ring_double[] = {"C=1CCCCC1", "C1CCCCC=1", "C=1CCCCC=1"};
    for (size_t i = 0; i < ARRAY_SIZE(ring_double); ++i) {
        ASSERT_TRUE(md_smiles_parse(&g, str_from_cstr(ring_double[i]), alloc, &err));
        ASSERT_EQ(6, (int)g.num_bonds);
        int num_double = 0;
        for (size_t b = 0; b < g.num_bonds; ++b) num_double += g.bonds[b].order == MD_SMILES_BOND_DOUBLE;
        EXPECT_EQ(1, num_double);
        EXPECT_EQ(MD_SMILES_BOND_DOUBLE, g.bonds[5].order);
        md_smiles_free(&g);
    }

    ASSERT_TRUE(md_smiles_parse(&g, STR_LIT("[13CH3:2][NH3+]"), alloc, &err));
    EXPECT_EQ(6,  g.atoms[0].z);
    EXPECT_EQ(13, g.atoms[0].isotope);
    EXPECT_EQ(3,  g.atoms[0].h_count);
    EXPECT_EQ(2,  g.atoms[0].atom_class);
    EXPECT_TRUE(g.atoms[0].flags & MD_SMILES_ATOM_BRACKET);
    EXPECT_EQ(7,  g.atoms[1].z);
    EXPECT_EQ(3,  g.atoms[1].h_count);
    EXPECT_EQ(1,  g.atoms[1].charge);
    md_smiles_free(&g);

    ASSERT_TRUE(md_smiles_parse(&g, STR_LIT("[Fe++].[O--].[Cu+2].[Cl-]"), alloc, &err));
    EXPECT_EQ(4, (int)g.num_components);
    EXPECT_EQ(0, (int)g.num_bonds);
    EXPECT_EQ(2,  g.atoms[0].charge);
    EXPECT_EQ(-2, g.atoms[1].charge);
    EXPECT_EQ(2,  g.atoms[2].charge);
    EXPECT_EQ(-1, g.atoms[3].charge);
    EXPECT_EQ(17, g.atoms[3].z);
    md_smiles_free(&g);

    // A ring bond joins the parts on either side of a '.'
    ASSERT_TRUE(md_smiles_parse(&g, STR_LIT("C1.C1"), alloc, &err));
    EXPECT_EQ(1, (int)g.num_bonds);
    EXPECT_EQ(1, (int)g.num_components);
    md_smiles_free(&g);

    ASSERT_TRUE(md_smiles_parse(&g, STR_LIT("c1cc[nH]c1"), alloc, &err));
    EXPECT_TRUE(g.atoms[0].flags & MD_SMILES_ATOM_AROMATIC);
    EXPECT_EQ(7, g.atoms[3].z);
    EXPECT_EQ(1, g.atoms[3].h_count);
    EXPECT_TRUE(g.atoms[3].flags & MD_SMILES_ATOM_AROMATIC);
    md_smiles_free(&g);

    ASSERT_TRUE(md_smiles_parse(&g, STR_LIT("[se]1cccc1"), alloc, &err));
    EXPECT_EQ(34, g.atoms[0].z);
    EXPECT_TRUE(g.atoms[0].flags & MD_SMILES_ATOM_AROMATIC);
    md_smiles_free(&g);

    ASSERT_TRUE(md_smiles_parse(&g, STR_LIT("ClCBr.[Sc].[Cs].[C@@H](F)(Cl)Br.*C"), alloc, &err));
    EXPECT_EQ(17, g.atoms[0].z);
    EXPECT_EQ(6,  g.atoms[1].z);
    EXPECT_EQ(35, g.atoms[2].z);
    EXPECT_EQ(21, g.atoms[3].z);
    EXPECT_EQ(55, g.atoms[4].z);
    EXPECT_TRUE(g.atoms[5].flags & MD_SMILES_ATOM_CHIRAL_CW);
    EXPECT_EQ(0,  g.atoms[9].z);
    md_smiles_free(&g);

    ASSERT_TRUE(md_smiles_parse(&g, STR_LIT("C%10CCCCC%10"), alloc, &err));
    EXPECT_EQ(6, (int)g.num_bonds);
    md_smiles_free(&g);

    // Branches, and offsets which count from the string as given, surrounding whitespace included
    ASSERT_TRUE(md_smiles_parse(&g, STR_LIT("  CC(=O)O  "), alloc, &err));
    ASSERT_EQ(4, (int)g.num_atoms);
    ASSERT_EQ(3, (int)g.num_bonds);
    EXPECT_EQ(2, (int)g.atoms[0].offset);
    EXPECT_EQ(6, (int)g.atoms[2].offset);
    EXPECT_EQ(1, (int)g.bonds[1].a);
    EXPECT_EQ(2, (int)g.bonds[1].b);
    EXPECT_EQ(MD_SMILES_BOND_DOUBLE, g.bonds[1].order);
    EXPECT_EQ(1, (int)g.bonds[2].a);
    EXPECT_EQ(3, (int)g.bonds[2].b);
    md_smiles_free(&g);
}

UTEST(match, smiles_errors) {
    md_allocator_i* alloc = md_get_heap_allocator();
    const struct {
        const char* str;
        size_t offset;
    } cases[] = {
        {"",            0},
        {"C(",          1},     // Unclosed '('
        {"C)",          1},     // Unbalanced ')'
        {"C1CC",        1},     // Ring bond never closed
        {"C==C",        2},     // Two bonds in a row
        {"=C",          0},     // Bond without an atom before it
        {"C=",          1},     // Bond without an atom after it
        {"C()",         2},     // Empty branch
        {"[Xx]",        1},     // Unknown element
        {"[C",          0},     // Unclosed '['
        {"Q",           0},     // Not in the organic subset
        {"CC C",        2},     // Whitespace inside
        {"C11",         2},     // Ring bond to itself
        {"C12CCCCC12",  9},     // Second ring bond duplicates the first
        {"C=1CCCCC-1",  9},     // The two ends of the ring bond disagree
        {"C.",          1},     // '.' without an atom after it
        {"C(C.)C",      3},     // Same, within a branch
        {"[CH3+",       0},     // Unclosed '['
    };
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        md_smiles_t g = {0};
        md_smiles_error_t err = {0};
        const bool ok = md_smiles_parse(&g, str_from_cstr(cases[i].str), alloc, &err);
        EXPECT_FALSE_MSG(ok, cases[i].str);
        EXPECT_EQ_MSG(cases[i].offset, err.offset, cases[i].str);
        EXPECT_NE_MSG('\0', err.message[0], cases[i].str);
        EXPECT_EQ_MSG(0, (int)g.num_atoms, cases[i].str);
    }
}

// ### QUERY ###

UTEST(match, query_from_smiles) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_match_query_t q = {0};
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("c1ccccc1C(=O)[O-]"), alloc, NULL));
    ASSERT_EQ(9, (int)q.num_atoms);
    ASSERT_EQ(9, (int)q.num_bonds);
    // Organic subset atoms leave hydrogens and charge open, bracket atoms fix them
    EXPECT_FALSE(q.atoms[0].flags & MD_MATCH_ATOM_HCOUNT);
    EXPECT_FALSE(q.atoms[0].flags & MD_MATCH_ATOM_CHARGE);
    EXPECT_TRUE (q.atoms[0].flags & MD_MATCH_ATOM_AROMATIC);
    EXPECT_FALSE(q.atoms[6].flags & MD_MATCH_ATOM_AROMATIC);
    EXPECT_TRUE (q.atoms[8].flags & MD_MATCH_ATOM_HCOUNT);
    EXPECT_TRUE (q.atoms[8].flags & MD_MATCH_ATOM_CHARGE);
    EXPECT_EQ(0,  q.atoms[8].h_max);
    EXPECT_EQ(-1, q.atoms[8].charge);
    // Implicit bonds: aromatic between aromatic atoms, single otherwise
    EXPECT_EQ(MD_MATCH_BOND_AROMATIC, q.bonds[0].order);
    md_match_query_free(&q);

    // Hydrogens written as atoms count towards the bracket atom they are bonded to
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("N[C@@]([H])(C)C=O"), alloc, NULL));
    EXPECT_EQ(1, q.atoms[1].h_min);
    EXPECT_EQ(1, q.atoms[1].h_max);
    md_match_query_free(&q);

    // In a SMILES which writes aromatic atoms lowercase, the uppercase atoms outside of rings are aliphatic. Those in
    // rings are not held to it, and neither is anything in a SMILES without lowercase atoms.
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("Cc1ccccc1C1CC1"), alloc, NULL));
    EXPECT_TRUE (q.atoms[0].flags & MD_MATCH_ATOM_ALIPHATIC);
    EXPECT_TRUE (q.atoms[1].flags & MD_MATCH_ATOM_AROMATIC);
    EXPECT_FALSE(q.atoms[1].flags & MD_MATCH_ATOM_ALIPHATIC);
    for (int i = 7; i < 10; ++i) EXPECT_FALSE(q.atoms[i].flags & (MD_MATCH_ATOM_ALIPHATIC | MD_MATCH_ATOM_AROMATIC));
    md_match_query_free(&q);
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("CC1=CC=CC=C1"), alloc, NULL));
    for (size_t i = 0; i < q.num_atoms; ++i) EXPECT_FALSE(q.atoms[i].flags & MD_MATCH_ATOM_ALIPHATIC);
    md_match_query_free(&q);
    // A ring bond across a '.' closes no ring
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("C1.C1c1ccccc1"), alloc, NULL));
    EXPECT_TRUE(q.atoms[0].flags & MD_MATCH_ATOM_ALIPHATIC);
    EXPECT_TRUE(q.atoms[1].flags & MD_MATCH_ATOM_ALIPHATIC);
    md_match_query_free(&q);

    // Atom classes are tags
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("OC(=O)[CH:1](C)[NH2:2]"), alloc, NULL));
    EXPECT_EQ(0, q.atoms[0].tag);
    EXPECT_EQ(1, q.atoms[3].tag);
    EXPECT_EQ(2, q.atoms[5].tag);
    md_match_query_free(&q);

    // A syntax error reaches the caller
    md_smiles_error_t err = {0};
    EXPECT_FALSE(md_match_query_init_smiles(&q, STR_LIT("CC(C"), alloc, &err));
    EXPECT_EQ(2, (int)err.offset);
    EXPECT_EQ(0, (int)q.num_atoms);
}

UTEST(match, query_validation) {
    const md_system_t* sys = get_sys(SYS_PFTAA);
    md_match_atom_t atoms[3] = {
        {.flags = MD_MATCH_ATOM_ELEMENT, .z = 6},
        {.flags = MD_MATCH_ATOM_ELEMENT, .z = 6},
        {.flags = MD_MATCH_ATOM_ELEMENT, .z = 6},
    };
    md_match_bond_t out_of_range[] = {{0, 3, 0}};
    md_match_bond_t self[]         = {{1, 1, 0}};
    md_match_bond_t twice[]        = {{0, 1, 0}, {2, 1, 0}, {1, 0, 0}};
    md_match_bond_t fine[]         = {{0, 1, 0}, {1, 2, 0}};

    md_match_query_t q = {.num_atoms = 3, .atoms = atoms};
    md_match_desc_t desc = {.query = &q};
    md_match_result_t res = {0};
    md_allocator_i* alloc = md_get_heap_allocator();

    q.bonds = out_of_range; q.num_bonds = ARRAY_SIZE(out_of_range);
    EXPECT_FALSE(md_match_find(&res, &desc, sys, alloc));
    q.bonds = self;         q.num_bonds = ARRAY_SIZE(self);
    EXPECT_FALSE(md_match_find(&res, &desc, sys, alloc));
    q.bonds = twice;        q.num_bonds = ARRAY_SIZE(twice);
    EXPECT_FALSE(md_match_find(&res, &desc, sys, alloc));

    q.bonds = fine;         q.num_bonds = ARRAY_SIZE(fine);
    EXPECT_TRUE(md_match_find(&res, &desc, sys, alloc));
    EXPECT_LT(0, (int)res.count);
    EXPECT_EQ(3, (int)res.width);
    md_match_result_free(&res);

    atoms[1].flags |= MD_MATCH_ATOM_HCOUNT;
    atoms[1].h_min = 2;
    atoms[1].h_max = 1;
    EXPECT_FALSE(md_match_find(&res, &desc, sys, alloc));

    md_match_query_t empty = {0};
    desc.query = &empty;
    EXPECT_FALSE(md_match_find(&res, &desc, sys, alloc));
    size_t count = 7;
    EXPECT_FALSE(md_match_for_each(&count, &desc, sys, NULL, NULL));
    EXPECT_EQ(0, (int)count);

    // A valid search which finds nothing is not a failure
    md_match_atom_t gold = {.flags = MD_MATCH_ATOM_ELEMENT, .z = 79};
    md_match_query_t none = {.num_atoms = 1, .atoms = &gold};
    desc.query = &none;
    EXPECT_TRUE(md_match_for_each(&count, &desc, sys, NULL, NULL));
    EXPECT_EQ(0, (int)count);
}

// ### MODES ###
// PFTAA: five thiophenes in a row, with an acetic acid on four of them and formyl groups at the ends

UTEST(match, modes) {
    const md_system_t* sys = get_sys(SYS_PFTAA);

    // Ground truth from the bonds
    size_t num_cc = 0;
    for (size_t i = 0; i < sys->bond.count; ++i) {
        const md_atom_pair_t p = sys->bond.pairs[i];
        num_cc += md_atom_atomic_number(&sys->atom, p.idx[0]) == MD_Z_C && md_atom_atomic_number(&sys->atom, p.idx[1]) == MD_Z_C;
    }
    ASSERT_LT(0, (int)num_cc);

    // Without bond orders every C-C bond is a match: once as a set of atoms, twice as a mapping
    md_match_desc_t desc = {.bond_orders = MD_MATCH_RESOLVE_NEVER};
    desc.mode = MD_MATCH_MODE_UNIQUE;
    EXPECT_EQ(num_cc, count_smiles_ex("CC", &desc, sys));
    desc.mode = MD_MATCH_MODE_ALL;
    EXPECT_EQ(2 * num_cc, count_smiles_ex("CC", &desc, sys));

    // A thiophene has two automorphisms (the mirror through S)
    EXPECT_EQ(5,  count_smiles("c1cccs1", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE,       sys));
    EXPECT_EQ(10, count_smiles("c1cccs1", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_ALL,          sys));
    EXPECT_EQ(5,  count_smiles("c1cccs1", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_DISJOINT,     sys));
    EXPECT_EQ(1,  count_smiles("c1cccs1", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_ONE_PER_UNIT, sys));

    // Four thiophenes in a row fit two ways along the five, and the two overlap
    EXPECT_EQ(2, count_smiles("c1ccc(s1)-c1ccc(s1)-c1ccc(s1)-c1cccs1", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, sys));
    EXPECT_EQ(1, count_smiles("c1ccc(s1)-c1ccc(s1)-c1ccc(s1)-c1cccs1", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_DISJOINT, sys));

    // Each mapping of ALL is reported once
    md_match_query_t q = {0};
    md_allocator_i* alloc = md_get_heap_allocator();
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("c1cccs1"), alloc, NULL));
    md_match_desc_t all = {.query = &q, .mode = MD_MATCH_MODE_ALL};
    md_match_result_t res = {0};
    ASSERT_TRUE(md_match_find(&res, &all, sys, alloc));
    ASSERT_EQ(10, (int)res.count);
    for (size_t i = 0; i < res.count; ++i) {
        for (size_t j = i + 1; j < res.count; ++j) {
            EXPECT_NE(0, MEMCMP(res.atom_idx + i * res.width, res.atom_idx + j * res.width, sizeof(int32_t) * res.width));
        }
    }
    md_match_result_free(&res);

    // limit, and a callback which stops the search
    all.limit = 3;
    ASSERT_TRUE(md_match_find(&res, &all, sys, alloc));
    EXPECT_EQ(3, (int)res.count);
    md_match_result_free(&res);
    md_match_query_free(&q);
}

UTEST(match, symmetry) {
    // A symmetric query matches the same atoms in several ways: UNIQUE reports them once, ALL every way
    const md_system_t* sys = get_sys(SYS_CENTERED);
    const size_t num_rings = count_smiles("c1ccccc1", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, sys);
    EXPECT_LT(0, (int)num_rings);
    EXPECT_EQ(12 * num_rings, count_smiles("c1ccccc1", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_ALL, sys));
    // Atoms which hang from one atom alike are one set: the two oxygens of the C terminal carboxylate
    const md_system_t* ala = get_sys(SYS_ALA);
    EXPECT_EQ(1, count_smiles("CC(=O)O", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
    EXPECT_EQ(2, count_smiles("CC(=O)O", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_ALL,    ala));

    // 160 residues of a chain as a reference by element: each ring which can flip over, carboxylate, guanidinium and
    // pair of methyls doubles the ways it fits its own atoms (2^40 and more), and it is found once
    const md_system_t* tub = get_sys(SYS_TUBULIN);
    md_allocator_i* alloc = md_get_heap_allocator();
    const uint32_t comp_beg = tub->instance.comp_offset[0];
    ASSERT_LE(comp_beg + 160, tub->instance.comp_offset[1]);
    const uint32_t beg = tub->component.atom_offset[comp_beg];
    const uint32_t end = tub->component.atom_offset[comp_beg + 160];
    md_array(int32_t) ref = 0;
    for (uint32_t i = beg; i < end; ++i) md_array_push(ref, (int32_t)i, alloc);
    md_match_query_t q = {0};
    ASSERT_TRUE(md_match_query_init_atoms(&q, ref, md_array_size(ref), MD_MATCH_LABEL_ELEMENT, tub, alloc));
    md_match_desc_t desc = {.query = &q, .level = MD_MATCH_LEVEL_INSTANCE};
    const md_tick_t t0 = md_tick_now();
    EXPECT_EQ(1, (int)count_matches(&desc, tub, NULL, NULL));
    printf("Reference of 160 residues (%zu atoms), found in its chain: %.2f ms\n", md_array_size(ref), md_tick_to_milliseconds(md_tick_now() - t0));
    md_match_query_free(&q);
    md_array_free(ref, alloc);
}

static bool stop_after_two(const int32_t* atom_idx, size_t width, uint32_t unit, void* user_param) {
    (void)atom_idx; (void)width; (void)unit;
    int* n = (int*)user_param;
    return ++(*n) < 2;
}

UTEST(match, callback_stops) {
    const md_system_t* sys = get_sys(SYS_PFTAA);
    md_match_query_t q = {0};
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("C"), md_get_heap_allocator(), NULL));
    md_match_desc_t desc = {.query = &q};
    int n = 0;
    EXPECT_EQ(2, (int)count_matches(&desc, sys, stop_after_two, &n));
    EXPECT_EQ(2, n);
    md_match_query_free(&q);
}

// ### CHEMISTRY ###

UTEST(match, bond_orders) {
    const md_system_t* sys = get_sys(SYS_ALA);
    md_match_desc_t orders = {0};
    md_match_desc_t never  = {.bond_orders = MD_MATCH_RESOLVE_NEVER};

    // 14 backbone carbonyls and the C terminal carboxylate, whose two C-O are delocalized
    EXPECT_EQ(16, count_smiles_ex("C=O", &orders, sys));
    EXPECT_EQ(16, count_smiles_ex("C=O", &never,  sys));
    EXPECT_EQ(2,  count_smiles_ex("C-O", &orders, sys));
    EXPECT_EQ(16, count_smiles_ex("C-O", &never,  sys));
    EXPECT_EQ(0,  count_smiles_ex("C=C", &orders, sys));
    EXPECT_EQ(1,  count_smiles_ex("CC(=O)[O-]", &orders, sys));

    // Kekule and aromatic forms find each other, and the rings are aromatic
    const md_system_t* pftaa = get_sys(SYS_PFTAA);
    EXPECT_EQ(5, count_smiles_ex("c1cccs1",    &orders, pftaa));
    EXPECT_EQ(5, count_smiles_ex("C1=CC=CS1",  &orders, pftaa));
    EXPECT_EQ(0, count_smiles_ex("c1ccccc1",   &orders, pftaa));
}

UTEST(match, hydrogens) {
    const md_system_t* ala = get_sys(SYS_ALA);
    const md_system_t* crystal = get_sys(SYS_1K4R);

    // Counted hydrogens, explicit and implicit
    EXPECT_EQ(15, count_smiles("[CH3]", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
    EXPECT_EQ(15, count_smiles("[CH]",  MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
    EXPECT_EQ(0,  count_smiles("[CH2]", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
    EXPECT_EQ(1,  count_smiles("[NH3+]", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));

    // Without counts the bonded hydrogens answer, where some atom of the element has one
    const md_system_t ala_bare = without_chemistry(ala);
    EXPECT_EQ(15, count_smiles("[CH3]", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, &ala_bare));
    EXPECT_EQ(1,  count_smiles("[NH3+]", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, &ala_bare));

    // A structure without hydrogens: counted, they are implicit and still tell methyls apart. Without counts nothing
    // can be said, every carbon is a candidate, unless asked to test regardless.
    size_t num_c = 0;
    for (size_t i = 0; i < crystal->atom.count; ++i) num_c += md_atom_atomic_number(&crystal->atom, i) == MD_Z_C;
    const size_t num_ch3 = count_smiles("[CH3]", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, crystal);
    EXPECT_LT(0, (int)num_ch3);
    EXPECT_LT(num_ch3, num_c);

    const md_system_t crystal_bare = without_chemistry(crystal);
    EXPECT_EQ(num_c, count_smiles("[CH3]", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, &crystal_bare));
    md_match_desc_t always = {.hydrogens = MD_MATCH_RESOLVE_ALWAYS};
    EXPECT_EQ(0, count_smiles_ex("[CH3]", &always, &crystal_bare));

    // Query hydrogens are mapped, not searched: one match per methyl and its own hydrogens in the row
    md_allocator_i* alloc = md_get_heap_allocator();
    md_match_query_t q = {0};
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("[H]C([H])([H])C"), alloc, NULL));
    md_match_desc_t desc = {.query = &q, .mode = MD_MATCH_MODE_ALL};
    md_match_result_t res = {0};
    ASSERT_TRUE(md_match_find(&res, &desc, ala, alloc));
    EXPECT_EQ(15, (int)res.count);
    for (size_t i = 0; i < res.count; ++i) {
        const int32_t* row = res.atom_idx + i * res.width;
        const int32_t c = row[1];
        EXPECT_EQ(MD_Z_C, md_atom_atomic_number(&ala->atom, c));
        const int32_t h[3] = {row[0], row[2], row[3]};
        for (int k = 0; k < 3; ++k) {
            EXPECT_EQ(MD_Z_H, md_atom_atomic_number(&ala->atom, h[k]));
            EXPECT_TRUE(atoms_bonded(ala, c, h[k]));
        }
        EXPECT_TRUE(h[0] < h[1] && h[1] < h[2]);
    }
    md_match_result_free(&res);

    // In a structure without hydrogen atoms they have nothing to map onto
    ASSERT_TRUE(md_match_find(&res, &desc, crystal, alloc));
    EXPECT_LT(0, (int)res.count);
    if (res.count > 0) {
        EXPECT_EQ(-1, res.atom_idx[0]);
        EXPECT_LE(0,  res.atom_idx[1]);
    }
    md_match_result_free(&res);
    md_match_query_free(&q);

    // Selection with the hydrogens of the matched atoms
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("[CH3]"), alloc, NULL));
    desc.mode = MD_MATCH_MODE_UNIQUE;
    ASSERT_TRUE(md_match_find(&res, &desc, ala, alloc));
    md_bitfield_t bf[2] = {md_bitfield_create(alloc), md_bitfield_create(alloc)};
    EXPECT_EQ(2, (int)md_match_result_select(bf, 2, &res, ala, MD_MATCH_SELECT_HYDROGENS));
    EXPECT_EQ(4, (int)md_bitfield_popcount(&bf[0]));
    EXPECT_EQ(1, (int)md_match_result_select(bf, 1, &res, ala, MD_MATCH_SELECT_NONE));
    EXPECT_EQ(1, (int)md_bitfield_popcount(&bf[0]));
    md_bitfield_free(&bf[0]);
    md_bitfield_free(&bf[1]);
    md_match_result_free(&res);
    md_match_query_free(&q);
}

UTEST(match, bracket_hydrogens) {
    const md_system_t* ala = get_sys(SYS_ALA);
    md_match_desc_t desc = {.level = MD_MATCH_LEVEL_COMPONENT, .mode = MD_MATCH_MODE_ONE_PER_UNIT};

    // The CA of every residue, with its hydrogen written either way
    EXPECT_EQ(15, count_smiles_ex("N[C@@H](C)C=O",     &desc, ala));
    EXPECT_EQ(15, count_smiles_ex("N[C@@]([H])(C)C=O", &desc, ala));
    // A CA with two hydrogens is a glycine, of which there are none
    EXPECT_EQ(0,  count_smiles_ex("N[CH2]C=O",         &desc, ala));
    EXPECT_EQ(0,  count_smiles_ex("N[CH]([H])C=O",     &desc, ala));
}

UTEST(match, aliphatic) {
    const md_system_t* sys = get_sys(SYS_DNA);
    size_t num_his = 0;
    for (size_t i = 0; i < sys->component.count; ++i) num_his += str_eq_cstr(LBL_TO_STR(sys->component.name[i]), "HIS");
    ASSERT_LT(0, (int)num_his);

    // The imidazole of a histidine with the carbon it hangs from. A purine has an imidazole too, but the carbon next
    // to it is in the other ring, and aromatic.
    md_match_desc_t desc = {.level = MD_MATCH_LEVEL_COMPONENT, .mode = MD_MATCH_MODE_ONE_PER_UNIT};
    EXPECT_EQ(num_his, count_smiles_ex("Cc1cncn1", &desc, sys));
    // Written Kekule, nothing says the carbon is not aromatic: the purines are found as well
    EXPECT_LT(num_his, count_smiles_ex("CC1=CN=CN1", &desc, sys));
}

UTEST(match, protonation) {
    // A crystal structure does not show where its protons are: the acids are found whether md_chem made them
    // carboxylates or not, and not the esters or carbonyls, whose oxygens have other bonds.
    const md_system_t* sys = get_sys(SYS_TUBULIN);
    size_t num_acid = 0;
    for (size_t i = 0; i < sys->component.count; ++i) {
        const str_t name = LBL_TO_STR(sys->component.name[i]);
        num_acid += str_eq_cstr(name, "ASP") || str_eq_cstr(name, "GLU");
    }
    ASSERT_LT(0, (int)num_acid);

    md_match_desc_t desc = {.level = MD_MATCH_LEVEL_COMPONENT, .mode = MD_MATCH_MODE_ONE_PER_UNIT};
    EXPECT_EQ(num_acid, count_smiles_ex("C(=O)[OH]", &desc, sys));
    EXPECT_EQ(num_acid, count_smiles_ex("C(=O)[O-]", &desc, sys));
    // Held to the protonation md_chem gave them (pH 7), they are carboxylates
    md_match_desc_t always = desc;
    always.hydrogens = MD_MATCH_RESOLVE_ALWAYS;
    always.charges   = MD_MATCH_RESOLVE_ALWAYS;
    EXPECT_EQ(0,        count_smiles_ex("C(=O)[OH]", &always, sys));
    EXPECT_EQ(num_acid, count_smiles_ex("C(=O)[O-]", &always, sys));

    // With hydrogen atoms the protonation is known: the only carboxylate is the C terminus
    const md_system_t* ala = get_sys(SYS_ALA);
    EXPECT_EQ(0, count_smiles_ex("C(=O)[OH]", &desc, ala));
    EXPECT_EQ(1, count_smiles_ex("C(=O)[O-]", &desc, ala));
}

UTEST(match, resonance) {
    // The C terminal carboxylate has its charge on one oxygen and shares it with the other
    const md_system_t* ala = get_sys(SYS_ALA);
    EXPECT_EQ(2, count_smiles("[O-]",    MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
    EXPECT_EQ(1, count_smiles("[O-]C=O", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
    // But it is one charge, not two, and stated for the whole group it has to be all of it
    EXPECT_EQ(0, count_smiles("[O-]C=[O-]",    MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
    EXPECT_EQ(1, count_smiles("[O-][C](=[O])C", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
    EXPECT_EQ(0, count_smiles("[O][C](=[O])C",  MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));

    // Every arginine is a guanidinium, whichever of its N carries the charge and the double bond
    const md_system_t* sys = get_sys(SYS_TUBULIN);
    md_allocator_i* alloc = md_get_heap_allocator();
    md_match_query_t q = {0};
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("NC(=[NH2+])N"), alloc, NULL));
    md_match_desc_t desc = {.query = &q, .level = MD_MATCH_LEVEL_COMPONENT, .mode = MD_MATCH_MODE_ONE_PER_UNIT};
    md_match_result_t res = {0};
    ASSERT_TRUE(md_match_find(&res, &desc, sys, alloc));
    size_t num_arg = 0, num_found = 0;
    for (size_t i = 0; i < sys->component.count; ++i) num_arg += str_eq_cstr(LBL_TO_STR(sys->component.name[i]), "ARG");
    for (size_t i = 0; i < res.count; ++i) num_found += str_eq_cstr(LBL_TO_STR(sys->component.name[res.unit[i]]), "ARG");
    EXPECT_LT(0, (int)num_arg);
    EXPECT_EQ(num_arg, num_found);
    md_match_result_free(&res);
    md_match_query_free(&q);
}

// ### THE GRAPH ###

UTEST(match, coordination) {
    // The Mg between the phosphates of the GTP is bonded to them, as a metal. That is coordination: the GTP is a
    // molecule of its own, and the whole of it.
    const md_system_t* sys = get_sys(SYS_TUBULIN);
    const char* gtp = "NC1=NC2=C(N=CN2C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)C(=O)N1";
    md_allocator_i* alloc = md_get_heap_allocator();
    md_match_query_t q = {0};
    ASSERT_TRUE(md_match_query_init_smiles(&q, str_from_cstr(gtp), alloc, NULL));
    md_match_desc_t desc = {.query = &q, .level = MD_MATCH_LEVEL_STRUCTURE, .flags = MD_MATCH_FLAG_WHOLE};
    md_match_result_t res = {0};
    ASSERT_TRUE(md_match_find(&res, &desc, sys, alloc));
    ASSERT_EQ(1, (int)res.count);

    // Its structure holds more than the GTP, the Mg at least
    const int32_t a = res.atom_idx[0];
    const uint32_t slot = (uint32_t)sys->structure.atom_slot[a];
    size_t structure_size = 0;
    for (size_t s = 0; s < sys->structure.count; ++s) {
        if (sys->structure.offset[s] <= slot && slot < sys->structure.offset[s + 1]) {
            structure_size = sys->structure.offset[s + 1] - sys->structure.offset[s];
        }
    }
    EXPECT_LT(q.num_atoms, structure_size);
    md_match_result_free(&res);
    md_match_query_free(&q);

    // The Mg is a molecule of its own too
    EXPECT_LT(0, (int)count_smiles_ex("[Mg]", &(md_match_desc_t){.flags = MD_MATCH_FLAG_WHOLE}, sys));
}

UTEST(match, molecules) {
    // Without coordination bonds and virtual sites the molecules are the structures, index for index
    const md_system_t* sys = get_sys(SYS_ALA);
    md_allocator_i* alloc = md_get_heap_allocator();
    md_match_query_t q = {0};
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("O"), alloc, NULL));
    md_match_desc_t desc = {.query = &q, .level = MD_MATCH_LEVEL_STRUCTURE};
    md_match_result_t res = {0};
    ASSERT_TRUE(md_match_find(&res, &desc, sys, alloc));
    ASSERT_LT(0, (int)res.count);
    for (size_t i = 0; i < res.count; ++i) {
        const uint32_t slot = (uint32_t)sys->structure.atom_slot[res.atom_idx[i]];
        const uint32_t u = res.unit[i];
        ASSERT_LT(u, (uint32_t)sys->structure.count);
        EXPECT_TRUE(sys->structure.offset[u] <= slot && slot < sys->structure.offset[u + 1]);
    }
    md_match_result_free(&res);
    md_match_query_free(&q);
}

UTEST(match, virtual_sites) {
    // The M sites of TIP4P water are not atoms of the graph: a water is its O and two H, whole, with or without
    // residues
    const md_system_t* sys = get_sys(SYS_TIP4P);
    size_t num_sol = 0;
    for (size_t i = 0; i < sys->component.count; ++i) num_sol += str_eq_cstr(LBL_TO_STR(sys->component.name[i]), "SOL");
    ASSERT_LT(0, (int)num_sol);
    md_match_desc_t desc = {.flags = MD_MATCH_FLAG_WHOLE};
    desc.level = MD_MATCH_LEVEL_STRUCTURE;
    EXPECT_EQ(num_sol, count_smiles_ex("[OH2]", &desc, sys));
    desc.level = MD_MATCH_LEVEL_COMPONENT;
    EXPECT_EQ(num_sol, count_smiles_ex("[OH2]", &desc, sys));
    // Nor are they any atom
    desc.flags = MD_MATCH_FLAG_NONE;
    size_t num_heavy = 0;
    for (size_t i = 0; i < sys->atom.count; ++i) num_heavy += md_atom_atomic_number(&sys->atom, i) > MD_Z_H;
    EXPECT_EQ(num_heavy, count_smiles_ex("*", &desc, sys));
}

// ### LEVELS, WHOLE AND MASK ###

UTEST(match, whole) {
    const md_system_t* sys = get_sys(SYS_CENTERED);

    // Glycine is part of every amino acid
    size_t num_aa = 0, num_gly = 0;
    for (size_t i = 0; i < sys->component.count; ++i) {
        const str_t name = LBL_TO_STR(sys->component.name[i]);
        if (md_util_resname_amino_acid(name)) num_aa++;
        if (str_eq_cstr(name, "GLY")) {
            // Without its C terminal oxygen
            size_t heavy = 0;
            const md_urange_t range = md_component_atom_range(&sys->component, i);
            for (uint32_t a = range.beg; a < range.end; ++a) heavy += md_atom_atomic_number(&sys->atom, a) != MD_Z_H;
            num_gly += heavy == 4;
        }
    }
    ASSERT_LT(0, (int)num_gly);

    md_match_desc_t desc = {.level = MD_MATCH_LEVEL_COMPONENT, .mode = MD_MATCH_MODE_ONE_PER_UNIT};
    EXPECT_EQ(num_aa, count_smiles_ex("NCC=O", &desc, sys));
    desc.flags = MD_MATCH_FLAG_WHOLE;
    EXPECT_EQ(num_gly, count_smiles_ex("NCC=O", &desc, sys));
    // With WHOLE, the modes which keep distinct atoms come to one match per unit
    desc.mode = MD_MATCH_MODE_UNIQUE;
    EXPECT_EQ(num_gly, count_smiles_ex("NCC=O", &desc, sys));
}

UTEST(match, mask) {
    const md_system_t* sys = get_sys(SYS_CENTERED);
    md_allocator_i* alloc = md_get_heap_allocator();

    // Phenyl rings: in PHE and TYR
    size_t num_rings = 0, num_rings_chain0 = 0;
    const md_urange_t chain0 = md_system_instance_atom_range(sys, 0);
    for (size_t i = 0; i < sys->component.count; ++i) {
        const str_t name = LBL_TO_STR(sys->component.name[i]);
        if (str_eq_cstr(name, "PHE") || str_eq_cstr(name, "TYR")) {
            num_rings++;
            num_rings_chain0 += md_component_atom_range(&sys->component, i).beg >= chain0.beg && md_component_atom_range(&sys->component, i).end <= chain0.end;
        }
    }

    md_match_query_t q = {0};
    ASSERT_TRUE(md_match_query_init_smiles(&q, STR_LIT("c1ccccc1"), alloc, NULL));
    md_match_desc_t desc = {.query = &q};
    EXPECT_EQ(num_rings, count_matches(&desc, sys, NULL, NULL));

    md_bitfield_t mask = md_bitfield_create(alloc);
    md_bitfield_set_range(&mask, chain0.beg, chain0.end);
    desc.mask = &mask;
    md_match_result_t res = {0};
    ASSERT_TRUE(md_match_find(&res, &desc, sys, alloc));
    EXPECT_EQ(num_rings_chain0, res.count);
    for (size_t i = 0; i < res.count * res.width; ++i) {
        EXPECT_TRUE(md_bitfield_test_bit(&mask, res.atom_idx[i]));
    }
    md_match_result_free(&res);
    md_bitfield_free(&mask);
    md_match_query_free(&q);
}

// ### REFERENCES ###

UTEST(match, reference_chain) {
    const md_system_t* sys = get_sys(SYS_CENTERED);
    md_allocator_i* alloc = md_get_heap_allocator();

    // The heavy atoms of the first chain, found in every chain
    md_structure_t structure = {0};
    md_structure_extract(&structure, &sys->structure, 0);
    md_array(int32_t) ref = 0;
    for (size_t i = 0; i < structure.count; ++i) {
        if (md_atom_atomic_number(&sys->atom, structure.atom_idx[i]) != MD_Z_H) md_array_push(ref, structure.atom_idx[i], alloc);
    }

    md_match_query_t q = {0};
    ASSERT_TRUE(md_match_query_init_atoms(&q, ref, md_array_size(ref), MD_MATCH_LABEL_ELEMENT, sys, alloc));
    md_match_desc_t desc = {.query = &q, .level = MD_MATCH_LEVEL_INSTANCE, .mode = MD_MATCH_MODE_ONE_PER_UNIT};
    const md_tick_t t0 = md_tick_now();
    const size_t count = count_matches(&desc, sys, NULL, NULL);
    const md_tick_t t1 = md_tick_now();
    printf("Reference chain of %zu atoms, %zu matches: %.3f ms\n", md_array_size(ref), count, md_tick_to_milliseconds(t1 - t0));
    EXPECT_EQ(253, (int)count);

    // The first match of the first chain is the reference itself
    md_match_result_t res = {0};
    desc.limit = 1;
    ASSERT_TRUE(md_match_find(&res, &desc, sys, alloc));
    ASSERT_EQ(1, (int)res.count);
    for (size_t i = 0; i < res.width; ++i) {
        EXPECT_EQ(md_atom_atomic_number(&sys->atom, ref[i]), md_atom_atomic_number(&sys->atom, res.atom_idx[i]));
    }
    md_match_result_free(&res);
    md_match_query_free(&q);
    md_array_free(ref, alloc);
}

UTEST(match, reference_molecule) {
    const md_system_t* sys = get_sys(SYS_CENTERED);
    md_allocator_i* alloc = md_get_heap_allocator();

    // A whole PFTAA, hydrogens included, found in every PFTAA
    md_structure_t structure = {0};
    md_structure_extract(&structure, &sys->structure, 253);
    md_match_query_t q = {0};
    ASSERT_TRUE(md_match_query_init_atoms(&q, structure.atom_idx, structure.count, MD_MATCH_LABEL_ELEMENT, sys, alloc));
    md_match_desc_t desc = {.query = &q, .level = MD_MATCH_LEVEL_COMPONENT, .mode = MD_MATCH_MODE_ONE_PER_UNIT};

    md_match_result_t res = {0};
    ASSERT_TRUE(md_match_find(&res, &desc, sys, alloc));
    EXPECT_EQ(61, (int)res.count);

    // Every hydrogen of the reference has a hydrogen of the match, bonded to the right atom, and no atom is used twice
    for (size_t i = 0; i < res.count; ++i) {
        const int32_t* row = res.atom_idx + i * res.width;
        for (size_t j = 0; j < res.width; ++j) {
            ASSERT_LE(0, row[j]);
            EXPECT_EQ(md_atom_atomic_number(&sys->atom, structure.atom_idx[j]), md_atom_atomic_number(&sys->atom, row[j]));
            for (size_t k = j + 1; k < res.width; ++k) EXPECT_NE(row[j], row[k]);
        }
        for (size_t b = 0; b < q.num_bonds; ++b) {
            EXPECT_TRUE(atoms_bonded(sys, row[q.bonds[b].a], row[q.bonds[b].b]));
        }
    }
    md_match_result_free(&res);
    md_match_query_free(&q);
}

UTEST(match, reference_identity) {
    // A reference given in the order of its atoms maps onto itself in its own unit, whatever its symmetries: the ring
    // of a TYR flipped over, the hydrogens of its CH2 and the oxygens of a carboxylate swapped are the same atoms
    md_allocator_i* alloc = md_get_heap_allocator();
    const struct { sys_id_t id; const char* name; } cases[] = {{SYS_CENTERED, "TYR"}, {SYS_CENTERED, "PHE"}, {SYS_CENTERED, "ASP"}, {SYS_CENTERED, "ARG"}, {SYS_PFTAA, "PFT"}};
    for (size_t c = 0; c < ARRAY_SIZE(cases); ++c) {
        const md_system_t* sys = get_sys(cases[c].id);
        for (size_t ci = 0; ci < sys->component.count; ++ci) {
            if (!str_eq_cstr(LBL_TO_STR(sys->component.name[ci]), cases[c].name)) continue;
            const md_urange_t range = md_component_atom_range(&sys->component, ci);
            md_array(int32_t) ref = 0;
            for (uint32_t i = range.beg; i < range.end; ++i) md_array_push(ref, (int32_t)i, alloc);
            for (int label = 0; label < 2; ++label) {
                md_match_query_t q = {0};
                ASSERT_TRUE(md_match_query_init_atoms(&q, ref, md_array_size(ref), (md_match_label_t)label, sys, alloc));
                md_bitfield_t mask = md_bitfield_create(alloc);
                md_bitfield_set_range(&mask, range.beg, range.end);
                md_match_desc_t desc = {.query = &q, .level = MD_MATCH_LEVEL_COMPONENT, .mode = MD_MATCH_MODE_UNIQUE, .flags = MD_MATCH_FLAG_WHOLE, .mask = &mask};
                md_match_result_t res = {0};
                ASSERT_TRUE(md_match_find(&res, &desc, sys, alloc));
                EXPECT_EQ(1, (int)res.count);
                for (size_t j = 0; j < res.width && res.count; ++j) EXPECT_EQ(ref[j], res.atom_idx[j]);
                md_match_result_free(&res);
                md_bitfield_free(&mask);
                md_match_query_free(&q);
            }
            md_array_free(ref, alloc);
            break;
        }
    }
}

UTEST(match, reference_ring) {
    const md_system_t* sys = get_sys(SYS_PFTAA);
    md_allocator_i* alloc = md_get_heap_allocator();

    // A thiophene ring
    {
        const int32_t ref[] = {19, 20, 21, 22, 24};
        md_match_query_t q = {0};
        ASSERT_TRUE(md_match_query_init_atoms(&q, ref, ARRAY_SIZE(ref), MD_MATCH_LABEL_ELEMENT, sys, alloc));
        md_match_desc_t desc = {.query = &q};
        EXPECT_EQ(5, (int)count_matches(&desc, sys, NULL, NULL));
        md_match_query_free(&q);
    }
    // Half of the molecule, minus the middle ring which joins the two halves
    {
        const int32_t ref[] = {0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 46, 47, 48};
        md_match_query_t q = {0};
        ASSERT_TRUE(md_match_query_init_atoms(&q, ref, ARRAY_SIZE(ref), MD_MATCH_LABEL_ELEMENT, sys, alloc));
        md_match_desc_t desc = {.query = &q};
        EXPECT_EQ(2, (int)count_matches(&desc, sys, NULL, NULL));
        md_match_query_free(&q);
    }
}

UTEST(match, reference_hydrogens) {
    const md_system_t* sys = get_sys(SYS_ALA);
    md_allocator_i* alloc = md_get_heap_allocator();

    // The fifth residue, hydrogens included, which pins the protonation of its atoms: not the N terminal residue,
    // whose N carries three hydrogens. The C terminal one has an uncharged oxygen for the carbonyl.
    const md_urange_t range = md_component_atom_range(&sys->component, 4);
    md_array(int32_t) ref = 0;
    for (uint32_t i = range.beg; i < range.end; ++i) md_array_push(ref, (int32_t)i, alloc);

    md_match_query_t q = {0};
    md_match_desc_t desc = {.query = &q, .level = MD_MATCH_LEVEL_COMPONENT, .mode = MD_MATCH_MODE_ONE_PER_UNIT};
    ASSERT_TRUE(md_match_query_init_atoms(&q, ref, md_array_size(ref), MD_MATCH_LABEL_ELEMENT, sys, alloc));
    EXPECT_EQ(14, (int)count_matches(&desc, sys, NULL, NULL));
    md_match_query_free(&q);

    // Its heavy atoms alone are a skeleton, found in every residue
    md_array_shrink(ref, 0);
    for (uint32_t i = range.beg; i < range.end; ++i) {
        if (md_atom_atomic_number(&sys->atom, i) != MD_Z_H) md_array_push(ref, (int32_t)i, alloc);
    }
    ASSERT_TRUE(md_match_query_init_atoms(&q, ref, md_array_size(ref), MD_MATCH_LABEL_ELEMENT, sys, alloc));
    EXPECT_EQ(15, (int)count_matches(&desc, sys, NULL, NULL));
    md_match_query_free(&q);
    md_array_free(ref, alloc);
}

UTEST(match, reference_names) {
    const md_system_t* sys = get_sys(SYS_ALA);
    md_allocator_i* alloc = md_get_heap_allocator();

    // A CB: by element every carbon, by name the CBs
    int32_t cb = -1;
    const md_urange_t range = md_component_atom_range(&sys->component, 4);
    for (uint32_t i = range.beg; i < range.end; ++i) {
        if (str_eq_cstr(md_atom_name(&sys->atom, i), "CB")) cb = (int32_t)i;
    }
    ASSERT_LE(0, cb);

    size_t num_c = 0;
    for (size_t i = 0; i < sys->atom.count; ++i) num_c += md_atom_atomic_number(&sys->atom, i) == MD_Z_C;

    md_match_query_t q = {0};
    md_match_desc_t desc = {.query = &q};
    ASSERT_TRUE(md_match_query_init_atoms(&q, &cb, 1, MD_MATCH_LABEL_ELEMENT, sys, alloc));
    EXPECT_EQ(num_c, count_matches(&desc, sys, NULL, NULL));
    md_match_query_free(&q);

    ASSERT_TRUE(md_match_query_init_atoms(&q, &cb, 1, MD_MATCH_LABEL_NAME, sys, alloc));
    md_match_result_t res = {0};
    ASSERT_TRUE(md_match_find(&res, &desc, sys, alloc));
    EXPECT_EQ(15, (int)res.count);
    for (size_t i = 0; i < res.count; ++i) {
        EXPECT_TRUE(str_eq_cstr(md_atom_name(&sys->atom, res.atom_idx[i]), "CB"));
    }
    md_match_result_free(&res);
    md_match_query_free(&q);
}

// ### RESIDUES ###
// Every residue told apart by its heavy atoms alone: a SMILES per residue, matched at the component level with
// MD_MATCH_FLAG_WHOLE. Bond orders are left out, so the test holds the topology and the elements and not the
// perception of the orders. The C terminal variants carry the second oxygen.

typedef struct {
    const char* name;
    const char* smiles[2];
} residue_pattern_t;

#define AA(name, side) {name, {"NC(" side ")C=O", "NC(" side ")C(=O)O"}}

static const residue_pattern_t residue_patterns[] = {
    AA("ALA", "C"),
    AA("ARG", "CCCNC(N)=N"),
    AA("ASN", "CC(N)=O"),
    AA("ASP", "CC(=O)O"),
    AA("CYS", "CS"),
    AA("GLN", "CCC(N)=O"),
    AA("GLU", "CCC(=O)O"),
    {"GLY", {"NCC=O", "NCC(=O)O"}},
    AA("HIS", "Cc1cncn1"),
    AA("ILE", "C(C)CC"),
    AA("LEU", "CC(C)C"),
    AA("LYS", "CCCCN"),
    AA("MET", "CCSC"),
    AA("PHE", "Cc1ccccc1"),
    {"PRO", {"N1CCCC1C=O", "N1CCCC1C(=O)O"}},
    AA("SER", "CO"),
    AA("THR", "C(C)O"),
    AA("TRP", "Cc1cnc2ccccc12"),
    AA("TYR", "Cc1ccc(O)cc1"),
    AA("VAL", "C(C)C"),
    // Deoxynucleotides, phosphate included; the 5' terminal ones without it
    {"DA",  {"P(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)CC1O", "OCC1OC(n2cnc3c(N)ncnc32)CC1O"}},
    {"DC",  {"P(=O)(O)OCC1OC(N2C=CC(N)=NC2=O)CC1O",  "OCC1OC(N2C=CC(N)=NC2=O)CC1O"}},
    {"DG",  {"P(=O)(O)OCC1OC(n2cnc3c2N=C(N)NC3=O)CC1O", "OCC1OC(n2cnc3c2N=C(N)NC3=O)CC1O"}},
    {"DT",  {"P(=O)(O)OCC1OC(N2C=C(C)C(=O)NC2=O)CC1O", "OCC1OC(N2C=C(C)C(=O)NC2=O)CC1O"}},
    {"HOH", {"O"}},     // Also SOL and WAT, see residue_kind
};

#undef AA

static bool residue_kind(str_t name, const char* kind) {
    if (str_eq_cstr(name, "SOL") || str_eq_cstr(name, "WAT")) name = STR_LIT("HOH");
    return str_eq_cstr(name, kind);
}

static void identify_residues(int* utest_result, sys_id_t id, bool complete) {
    const md_system_t* sys = get_sys(id);
    md_allocator_i* alloc = md_get_heap_allocator();
    const size_t num_comp = sys->component.count;
    const size_t num_pat  = ARRAY_SIZE(residue_patterns);

    // Which pattern matched each component, -1 for none, -2 for several
    int* matched = md_alloc(alloc, sizeof(int) * num_comp);
    for (size_t i = 0; i < num_comp; ++i) matched[i] = -1;

    for (size_t p = 0; p < num_pat; ++p) {
        for (size_t v = 0; v < ARRAY_SIZE(residue_patterns[p].smiles) && residue_patterns[p].smiles[v]; ++v) {
            md_match_query_t q = {0};
            md_smiles_error_t err = {0};
            if (!md_match_query_init_smiles(&q, str_from_cstr(residue_patterns[p].smiles[v]), alloc, &err)) {
                printf("%s: %s at %zu\n", residue_patterns[p].smiles[v], err.message, err.offset);
                *utest_result = UTEST_TEST_FAILURE;
                continue;
            }
            md_match_desc_t desc = {
                .query = &q,
                .level = MD_MATCH_LEVEL_COMPONENT,
                .mode  = MD_MATCH_MODE_ONE_PER_UNIT,
                .flags = MD_MATCH_FLAG_WHOLE,
                .bond_orders = MD_MATCH_RESOLVE_NEVER,
            };
            md_match_result_t res = {0};
            md_match_find(&res, &desc, sys, alloc);
            for (size_t i = 0; i < res.count; ++i) {
                const uint32_t c = res.unit[i];
                matched[c] = (matched[c] == -1 || matched[c] == (int)p) ? (int)p : -2;
            }
            md_match_result_free(&res);
            md_match_query_free(&q);
        }
    }

    size_t num_missed = 0, num_wrong = 0, num_known = 0;
    for (size_t i = 0; i < num_comp; ++i) {
        const str_t name = LBL_TO_STR(sys->component.name[i]);
        int expected = -1;
        for (size_t p = 0; p < num_pat; ++p) {
            if (residue_kind(name, residue_patterns[p].name)) expected = (int)p;
        }
        num_known += expected >= 0;
        if (matched[i] >= 0 && matched[i] != expected) {
            if (num_wrong < 8) printf("Residue %zu '" STR_FMT "' matched as %s\n", i + 1, STR_ARG(name), matched[i] >= 0 ? residue_patterns[matched[i]].name : "several");
            num_wrong++;
        } else if (matched[i] == -2) {
            num_wrong++;
        } else if (expected >= 0 && matched[i] != expected) {
            if (complete && num_missed < 8) printf("Residue %zu '" STR_FMT "' not matched\n", i + 1, STR_ARG(name));
            num_missed++;
        }
    }
    printf("%zu residues of known kinds, %zu not identified, %zu misidentified\n", num_known, num_missed, num_wrong);

    // Nothing is ever taken for something else. Residues with atoms missing (crystal structures) are not identified.
    if (num_wrong) *utest_result = UTEST_TEST_FAILURE;
    if (complete && num_missed) *utest_result = UTEST_TEST_FAILURE;
    md_free(alloc, matched, sizeof(int) * num_comp);
}

UTEST(match, residues_md) {
    identify_residues(utest_result, SYS_ALA, true);
    identify_residues(utest_result, SYS_CENTERED, true);
    identify_residues(utest_result, SYS_DNA, true);
}

UTEST(match, residues_crystal) {
    identify_residues(utest_result, SYS_1K4R, false);
    identify_residues(utest_result, SYS_1FEZ, false);
    identify_residues(utest_result, SYS_2OR2, false);
}

// ### IDENTIFICATION ###

static md_match_library_t* residue_library(md_allocator_i* alloc) {
    md_match_library_t* lib = md_match_library_create(alloc);
    for (size_t p = 0; p < ARRAY_SIZE(residue_patterns); ++p) {
        for (size_t v = 0; v < ARRAY_SIZE(residue_patterns[p].smiles) && residue_patterns[p].smiles[v]; ++v) {
            md_match_library_add_smiles(lib, str_from_cstr(residue_patterns[p].smiles[v]), str_from_cstr(residue_patterns[p].name), NULL);
        }
    }
    return lib;
}

// The first entry of the library which matches each unit, found with one search per entry
static void identify_by_search(int32_t* first, size_t num_units, const md_match_identify_desc_t* id, const md_system_t* sys, md_allocator_i* alloc) {
    for (size_t u = 0; u < num_units; ++u) first[u] = -1;
    for (size_t e = 0; e < md_match_library_count(id->library); ++e) {
        md_match_desc_t desc = {
            .query = md_match_library_query(id->library, e),
            .level = id->level,
            .mode  = MD_MATCH_MODE_ONE_PER_UNIT,
            .flags = id->flags,
            .hydrogens   = id->hydrogens,
            .bond_orders = id->bond_orders,
            .charges     = id->charges,
            .mask  = id->mask,
        };
        md_match_result_t res = {0};
        md_match_find(&res, &desc, sys, alloc);
        for (size_t i = 0; i < res.count; ++i) {
            if (first[res.unit[i]] < 0) first[res.unit[i]] = (int32_t)e;
        }
        md_match_result_free(&res);
    }
}

// names: the residues identified are those their names say. Not so in a crystal structure with side chains which are
// not all in the model: a GLN without its side chain is an ALA.
static void check_identify(int* utest_result, sys_id_t id, const md_match_identify_desc_t* desc, bool names) {
    const md_system_t* sys = get_sys(id);
    md_allocator_i* alloc = md_get_heap_allocator();
    const size_t num_units = sys->component.count;

    const md_tick_t t0 = md_tick_now();
    md_match_identify_result_t res = {0};
    const bool ok = md_match_identify(&res, desc, sys, alloc);
    const md_tick_t t1 = md_tick_now();
    EXPECT_TRUE(ok);

    int32_t* first = md_alloc(alloc, sizeof(int32_t) * num_units);
    identify_by_search(first, num_units, desc, sys, alloc);
    const md_tick_t t2 = md_tick_now();
    printf("%zu residues, %zu identified: %.2f ms in one pass, %.2f ms by one search per entry\n", num_units, res.count, md_tick_to_milliseconds(t1 - t0), md_tick_to_milliseconds(t2 - t1));

    // The same entry for the same units, and a valid mapping of the entry's atoms onto the unit's
    size_t num_first = 0;
    for (size_t u = 0; u < num_units; ++u) num_first += first[u] >= 0;
    EXPECT_EQ(num_first, res.count);
    for (size_t i = 0; i < res.count; ++i) {
        const uint32_t u = res.unit[i];
        const uint32_t e = res.entry[i];
        if (i > 0) EXPECT_LT(res.unit[i - 1], u);
        EXPECT_EQ(first[u], (int32_t)e);
        char kind[16];
        const str_t entry_name = md_match_library_name(desc->library, e);
        snprintf(kind, sizeof(kind), "%.*s", (int)entry_name.len, entry_name.ptr);
        if (names) EXPECT_TRUE(residue_kind(LBL_TO_STR(sys->component.name[u]), kind));

        const md_match_query_t* q = md_match_library_query(desc->library, e);
        ASSERT_EQ(q->num_atoms, (size_t)(res.offset[i + 1] - res.offset[i]));
        const int32_t* row = res.atom_idx + res.offset[i];
        const md_urange_t range = md_component_atom_range(&sys->component, u);
        for (size_t j = 0; j < q->num_atoms; ++j) {
            EXPECT_TRUE(row[j] >= (int32_t)range.beg && row[j] < (int32_t)range.end);
            EXPECT_EQ(q->atoms[j].z, md_atom_atomic_number(&sys->atom, row[j]));
        }
        for (size_t b = 0; b < q->num_bonds; ++b) {
            EXPECT_TRUE(atoms_bonded(sys, row[q->bonds[b].a], row[q->bonds[b].b]));
        }
    }
    md_free(alloc, first, sizeof(int32_t) * num_units);
    md_match_identify_result_free(&res);
}

UTEST(match, identify) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_match_library_t* lib = residue_library(alloc);
    EXPECT_EQ(2 * 20 + 2 * 4 + 1, (int)md_match_library_count(lib));

    md_match_identify_desc_t desc = {
        .library = lib,
        .level   = MD_MATCH_LEVEL_COMPONENT,
        .flags   = MD_MATCH_FLAG_WHOLE,
        .bond_orders = MD_MATCH_RESOLVE_NEVER,
    };
    check_identify(utest_result, SYS_CENTERED, &desc, true);
    check_identify(utest_result, SYS_DNA, &desc, true);
    check_identify(utest_result, SYS_1K4R, &desc, true);
    // With the bond orders, and in a crystal structure with ligands and truncated side chains
    desc.bond_orders = MD_MATCH_RESOLVE_AUTO;
    check_identify(utest_result, SYS_CENTERED, &desc, true);
    check_identify(utest_result, SYS_TUBULIN, &desc, false);

    // Within a part of the system
    const md_system_t* sys = get_sys(SYS_CENTERED);
    md_bitfield_t mask = md_bitfield_create(alloc);
    const md_urange_t chain0 = md_system_instance_atom_range(sys, 0);
    md_bitfield_set_range(&mask, chain0.beg, chain0.end);
    desc.mask = &mask;
    check_identify(utest_result, SYS_CENTERED, &desc, true);
    md_bitfield_free(&mask);

    md_match_library_destroy(lib);
}

UTEST(match, identify_order) {
    md_allocator_i* alloc = md_get_heap_allocator();
    const md_system_t* sys = get_sys(SYS_ALA);
    md_match_identify_result_t res = {0};

    // Without WHOLE an entry is found within a unit: the first which is, in the order added
    md_match_library_t* lib = md_match_library_create(alloc);
    EXPECT_EQ(0, md_match_library_add_smiles(lib, STR_LIT("NCC=O"), STR_LIT("backbone"), NULL));
    EXPECT_EQ(1, md_match_library_add_smiles(lib, STR_LIT("CC(N)C=O"), STR_LIT("ala"), NULL));
    md_smiles_error_t err = {0};
    EXPECT_EQ(-1, md_match_library_add_smiles(lib, STR_LIT("C1CC"), STR_LIT("broken"), &err));
    EXPECT_EQ(2, (int)md_match_library_count(lib));
    EXPECT_TRUE(str_eq_cstr(md_match_library_name(lib, 1), "ala"));

    md_match_identify_desc_t desc = {.library = lib, .level = MD_MATCH_LEVEL_COMPONENT};
    ASSERT_TRUE(md_match_identify(&res, &desc, sys, alloc));
    EXPECT_EQ(15, (int)res.count);
    for (size_t i = 0; i < res.count; ++i) EXPECT_EQ(0, (int)res.entry[i]);
    md_match_identify_result_free(&res);

    // With WHOLE neither is all of a residue (an ALA has a methyl, and the C terminal one a second oxygen) but the
    // second is all of the 14 others
    desc.flags = MD_MATCH_FLAG_WHOLE;
    ASSERT_TRUE(md_match_identify(&res, &desc, sys, alloc));
    EXPECT_EQ(14, (int)res.count);
    for (size_t i = 0; i < res.count; ++i) EXPECT_EQ(1, (int)res.entry[i]);
    md_match_identify_result_free(&res);
    md_match_library_destroy(lib);

    // An entry with a wildcard comes in its turn as well
    lib = md_match_library_create(alloc);
    md_match_library_add_smiles(lib, STR_LIT("CC(N)C=O"), STR_LIT("ala"), NULL);
    md_match_library_add_smiles(lib, STR_LIT("*C(N)C(=O)O"), STR_LIT("any C terminal"), NULL);
    md_match_library_add_smiles(lib, STR_LIT("CC(N)C(=O)O"), STR_LIT("C terminal ala"), NULL);
    desc.library = lib;
    ASSERT_TRUE(md_match_identify(&res, &desc, sys, alloc));
    EXPECT_EQ(15, (int)res.count);
    for (size_t i = 0; i < res.count; ++i) EXPECT_EQ(i + 1 < res.count ? 0 : 1, (int)res.entry[i]);
    md_match_identify_result_free(&res);
    md_match_library_destroy(lib);
}

// ### SCRIPT ###

static size_t filter_count(md_array(md_bitfield_t)* arr, const char* expr, sys_id_t id, md_allocator_i* alloc) {
    char err[256] = {0};
    bool is_dynamic = false;
    md_array_shrink(*arr, 0);
    if (!md_filter_evaluate(arr, str_from_cstr(expr), get_sys(id), get_state(id), NULL, &is_dynamic, err, sizeof(err), alloc)) {
        printf("'%s': %s\n", expr, err);
        return SIZE_MAX;
    }
    return md_array_size(*arr);
}

static bool filter_fails_with(const char* expr, const char* expected, sys_id_t id) {
    md_bitfield_t bf = md_bitfield_create(md_get_heap_allocator());
    char err[256] = {0};
    bool is_dynamic = false;
    const bool ok = md_filter(&bf, str_from_cstr(expr), get_sys(id), get_state(id), NULL, &is_dynamic, err, sizeof(err));
    md_bitfield_free(&bf);
    if (ok) {
        printf("'%s' was expected to fail\n", expr);
        return false;
    }
    if (!strstr(err, expected)) {
        printf("'%s' failed with '%s', expected '%s'\n", expr, err, expected);
        return false;
    }
    return true;
}

UTEST(match, script_smiles) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_array(md_bitfield_t) arr = 0;

    // One selection per methyl, with its hydrogens
    ASSERT_EQ(15, (int)filter_count(&arr, "smiles('[CH3]')", SYS_ALA, alloc));
    for (size_t i = 0; i < md_array_size(arr); ++i) {
        EXPECT_EQ(4, (int)md_bitfield_popcount(&arr[i]));
    }

    // Whole residues: the C terminal one has a second oxygen
    EXPECT_EQ(14, (int)filter_count(&arr, "smiles('NC(C)C=O', level='residue', mode='whole')", SYS_ALA, alloc));
    EXPECT_EQ(1,  (int)filter_count(&arr, "smiles('NC(C)C(=O)O', level='residue', mode='whole')", SYS_ALA, alloc));
    EXPECT_EQ(15, (int)filter_count(&arr, "smiles('NC(C)C=O', level='residue', mode='one_per_unit')", SYS_ALA, alloc));

    // Within a context
    EXPECT_EQ(3, (int)filter_count(&arr, "smiles('[CH3]') in residue(1:3)", SYS_ALA, alloc));

    // As one selection: the 15 carbonyl carbons and 16 oxygens, without hydrogens of their own
    md_bitfield_t bf = md_bitfield_create(alloc);
    char err[256];
    bool is_dynamic = false;
    EXPECT_TRUE(md_filter(&bf, STR_LIT("smiles('C=O')"), get_sys(SYS_ALA), get_state(SYS_ALA), NULL, &is_dynamic, err, sizeof(err)));
    EXPECT_FALSE(is_dynamic);
    EXPECT_EQ(31, (int)md_bitfield_popcount(&bf));

    EXPECT_TRUE(filter_fails_with("smiles('C1CC')", "Ring bond 1 is never closed", SYS_ALA));
    EXPECT_TRUE(filter_fails_with("smiles('CC', level='molecule')", "Unknown level", SYS_ALA));
    EXPECT_TRUE(filter_fails_with("smiles('CC', mode='first')", "Unknown mode", SYS_ALA));
    EXPECT_TRUE(filter_fails_with("smiles('CC', 'residue')", "by name", SYS_ALA));

    md_arena_allocator_destroy(alloc);
}

UTEST(match, script_match) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_array(md_bitfield_t) arr = 0;

    // The fifth residue with its hydrogens, see match.reference_hydrogens: each match has its every atom. Twice in the
    // C terminal residue, whose carboxylate shares its charge between the oxygens: either is the neutral carbonyl.
    ASSERT_EQ(15, (int)filter_count(&arr, "match(residue(5), level='residue')", SYS_ALA, alloc));
    for (size_t i = 0; i < md_array_size(arr); ++i) {
        EXPECT_EQ(10, (int)md_bitfield_popcount(&arr[i]));
    }
    // Its heavy atoms are found in every residue, and twice in the C terminal one, either oxygen as the carbonyl
    EXPECT_EQ(15, (int)filter_count(&arr, "match(residue(5) and not element('H'), level='residue', mode='one_per_unit')", SYS_ALA, alloc));
    EXPECT_EQ(16, (int)filter_count(&arr, "match(residue(5) and not element('H'))", SYS_ALA, alloc));
    EXPECT_EQ(15, (int)filter_count(&arr, "match(name('CB') in residue(5), by='name')", SYS_ALA, alloc));

    EXPECT_TRUE(filter_fails_with("match(within(3, residue(1)))", "cannot depend on the frame", SYS_ALA));
    EXPECT_TRUE(filter_fails_with("match(residue(5), by='type')", "Unknown label", SYS_ALA));

    md_arena_allocator_destroy(alloc);
}
