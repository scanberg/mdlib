#include "utest.h"

#include <md_pdb.h>
#include <md_gro.h>
#include <md_tpr.h>
#include <md_system.h>
#include <md_infer.h>
#include <core/md_arena_allocator.h>

#include <string.h>
#include <stdio.h>
#include <math.h>

#include "system_util.h"

// The coarse grained branch of covalent bond inference used to skip the first residue, test every bead
// against itself (always under the cutoff, so a self bond) and emit every other pair twice.
UTEST(martini, infer_bonds_coarse_grained) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));

    md_system_t sys = {0};
    sys.alloc = alloc;
    const char* lys[] = {"BB", "SC1", "SC2"};
    cg_add_residue(&sys, "LYS", 1, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);
    cg_add_residue(&sys, "LYS", 2, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);

    // BB beads 3.8A apart, side chain beads 3A above their BB and 3A from each other
    float x[] = {0, 0, 0,   3.8f, 3.8f, 3.8f};
    float y[] = {0, 3, 6,   0,    3,    6};
    float z[] = {0, 0, 0,   0,    0,    0};
    md_float3_t xyz[6];
    md_system_state_t state = { .num_atoms = 6, .xyz = pack_xyz(xyz, x, y, z, 6) };

    md_bond_data_t bond = {0};
    md_bond_infer(&bond, &sys, &state, alloc, NULL);

    // (0,1) (1,2) (3,4) (4,5) within residues, (0,3) between backbones. 0-2 is 6A, over the cutoff.
    EXPECT_EQ(bond.count, 5u);
    bool has[6][6] = {0};
    for (size_t i = 0; i < bond.count; ++i) {
        const int a = bond.pairs[i].idx[0];
        const int b = bond.pairs[i].idx[1];
        EXPECT_NE(a, b);
        EXPECT_FALSE(has[a][b] || has[b][a]);
        has[a][b] = true;
    }
    EXPECT_TRUE(has[0][1] && has[1][2] && has[3][4] && has[4][5] && has[0][3]);
    for (size_t i = 0; i < bond.count; ++i) {
        EXPECT_EQ((int)bond.flags[i], (int)md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_INFERRED));
    }

    md_vm_arena_destroy(alloc);
}

// martinize2 (vermouth 0.15, -ff martini3001) on residues 384-392 of 1k4r (YVGELSHQW), with the His renamed HSP as
// martinize2 names a charged one. Martini 3 has Tyr SC1-SC4 and Trp SC1-SC5. The file carries no bonds, so whether
// the chain holds together rests on the backbone beads.
static const char martini3_peptide_gro[] =
    "Martini 3 peptide (martinize2, YVGELSHQW, HIS as HSP)\n"
    "   26\n"
    "    1TYR     BB    1   1.878   1.149   1.615\n"
    "    1TYR    SC1    2   1.634   1.202   1.502\n"
    "    1TYR    SC2    3   1.450   1.367   1.498\n"
    "    1TYR    SC3    4   1.552   1.300   1.291\n"
    "    1TYR    SC4    5   1.363   1.462   1.288\n"
    "    2VAL     BB    6   1.793   1.092   1.968\n"
    "    2VAL    SC1    7   2.047   1.095   2.033\n"
    "    3GLY     BB    8   1.530   1.000   1.928\n"
    "    4GLU     BB    9   1.365   1.301   2.000\n"
    "    4GLU    SC1   10   1.000   1.167   2.124\n"
    "    5LEU     BB   11   1.625   1.416   1.973\n"
    "    5LEU    SC1   12   1.822   1.530   2.206\n"
    "    6SER     BB   13   1.746   1.709   1.799\n"
    "    6SER    SC1   14   1.552   1.747   1.679\n"
    "    7HSP     BB   15   1.960   1.818   1.504\n"
    "    7HSP    SC1   16   2.135   1.902   1.685\n"
    "    7HSP    SC2   17   2.220   2.095   1.803\n"
    "    7HSP    SC3   18   2.253   2.095   1.631\n"
    "    8GLN     BB   19   2.101   2.029   1.314\n"
    "    8GLN    SC1   20   1.781   2.141   1.166\n"
    "    9TRP     BB   21   2.304   1.959   1.000\n"
    "    9TRP    SC1   22   2.470   1.827   1.121\n"
    "    9TRP    SC2   23   2.679   1.768   1.005\n"
    "    9TRP    SC3   24   2.708   1.885   1.138\n"
    "    9TRP    SC4   25   2.709   2.018   1.298\n"
    "    9TRP    SC5   26   2.890   1.973   1.194\n"
    "   5.00000   5.00000   5.00000\n";

// The same as martinize2 wrote it, without its CONECT records
static const char martini3_peptide_pdb[] =
    "ATOM      1 BB   TYR A 384     -61.526 -17.356 229.581  1.00  0.00\n"
    "ATOM      2 SC1  TYR A 384     -63.958 -16.821 228.450  1.00  0.00\n"
    "ATOM      3 SC2  TYR A 384     -65.803 -15.175 228.406  1.00  0.00\n"
    "ATOM      4 SC3  TYR A 384     -64.786 -15.840 226.333  1.00  0.00\n"
    "ATOM      5 SC4  TYR A 384     -66.677 -14.229 226.304  1.00  0.00\n"
    "ATOM      6 BB   VAL A 385     -62.376 -17.927 233.108  1.00  0.00\n"
    "ATOM      7 SC1  VAL A 385     -59.828 -17.893 233.758  1.00  0.00\n"
    "ATOM      8 BB   GLY A 386     -65.004 -18.845 232.710  1.00  0.00\n"
    "ATOM      9 BB   GLU A 387     -66.648 -15.831 233.428  1.00  0.00\n"
    "ATOM     10 SC1  GLU A 387     -70.302 -17.177 234.666  1.00  0.00\n"
    "ATOM     11 BB   LEU A 388     -64.049 -14.684 233.156  1.00  0.00\n"
    "ATOM     12 SC1  LEU A 388     -62.078 -13.547 235.491  1.00  0.00\n"
    "ATOM     13 BB   SER A 389     -62.839 -11.758 231.415  1.00  0.00\n"
    "ATOM     14 SC1  SER A 389     -64.783 -11.373 230.216  1.00  0.00\n"
    "ATOM     15 BB   HSP A 390     -60.700 -10.669 228.464  1.00  0.00\n"
    "ATOM     16 SC1  HSP A 390     -58.952  -9.827 230.275  1.00  0.00\n"
    "ATOM     17 SC2  HSP A 390     -58.102  -7.893 231.462  1.00  0.00\n"
    "ATOM     18 SC3  HSP A 390     -57.772  -7.898 229.740  1.00  0.00\n"
    "ATOM     19 BB   GLN A 391     -59.295  -8.559 226.571  1.00  0.00\n"
    "ATOM     20 SC1  GLN A 391     -62.491  -7.433 225.082  1.00  0.00\n"
    "ATOM     21 BB   TRP A 392     -57.267  -9.252 223.427  1.00  0.00\n"
    "ATOM     22 SC1  TRP A 392     -55.603 -10.572 224.636  1.00  0.00\n"
    "ATOM     23 SC2  TRP A 392     -53.517 -11.161 223.475  1.00  0.00\n"
    "ATOM     24 SC3  TRP A 392     -53.223  -9.997 224.803  1.00  0.00\n"
    "ATOM     25 SC4  TRP A 392     -53.209  -8.667 226.404  1.00  0.00\n"
    "ATOM     26 SC5  TRP A 392     -51.401  -9.117 225.366  1.00  0.00\n"
    "TER      27      TRP A 392\n"
    "END\n";

// Every particle a bead without an element, every residue an amino acid whose BB bead is its one backbone bead and
// whose other beads are side chain, and the whole a single peptide
static void check_martini_peptide(int* utest_result, md_system_t* sys, md_system_state_t* state, size_t num_residues) {
    for (size_t t = 1; t < sys->atom.type.count; ++t) {
        EXPECT_EQ(0, sys->atom.type.z[t]);
        EXPECT_EQ(MD_PARTICLE_BEAD, md_particle_type_particle_kind(&sys->atom.type, t));
    }
    ASSERT_EQ(num_residues, md_system_residue_count(sys));
    for (size_t c = 0; c < num_residues; ++c) {
        EXPECT_EQ(MD_RESIDUE_KIND_AMINO_ACID, md_system_residue_kind(sys, c));
        const md_urange_t range = md_system_residue_atom_range(sys, c);
        size_t num_backbone = 0;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            const md_atom_flags_t flags = md_atom_flags(&sys->atom, i);
            const bool bb = md_str_eq(md_atom_name(&sys->atom, i), MD_STR_LIT("BB"));
            EXPECT_EQ((int)bb,  (int)((flags & MD_ATOM_FLAG_BACKBONE)   != 0));
            EXPECT_EQ((int)!bb, (int)((flags & MD_ATOM_FLAG_SIDE_CHAIN) != 0));
            num_backbone += (flags & MD_ATOM_FLAG_BACKBONE) != 0;
        }
        EXPECT_EQ((size_t)1, num_backbone);
    }

    ASSERT_TRUE(md_system_infer(sys, state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ((size_t)1, md_system_part_count(sys));
    EXPECT_EQ(num_residues, md_system_part_residue_range(sys, 0).end - md_system_part_residue_range(sys, 0).beg);
    EXPECT_EQ(MD_PART_KIND_PEPTIDE, md_system_part_kind(sys, 0));
}

UTEST(martini, martini3_peptide_gro) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, md_str_from_cstr(martini3_peptide_gro), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    check_martini_peptide(utest_result, &sys, &state, 9);
    md_vm_arena_destroy(alloc);
}

UTEST(martini, martini3_peptide_pdb) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(martini3_peptide_pdb), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    check_martini_peptide(utest_result, &sys, &state, 9);
    md_vm_arena_destroy(alloc);
}

// The residue names beyond the standard ones that martinize2 writes, in a chain of one each with the beads it gives
// them: the protonation states (CHARMM, Amber, GROMOS and Martini names), the lipidated residues of Martini 3 with their
// tails, Ala, which has a side chain bead in Martini 3 and none in Martini 2, and the charge dummies (SCP, SCN) of
// polarizable Martini 2. BB beads 3.5A apart, and the other beads of a residue 3A apart off its BB.
UTEST(martini, martini_residue_names) {
    static const struct { const char* name; const char* beads[11]; } res[] = {
        {"HIS", {"BB", "SC1", "SC2", "SC3"}},
        {"HSD", {"BB", "SC1", "SC2", "SC3"}},
        {"HSE", {"BB", "SC1", "SC2", "SC3"}},
        {"HSP", {"BB", "SC1", "SC2", "SC3"}},
        {"HID", {"BB", "SC1", "SC2", "SC3"}},
        {"HIE", {"BB", "SC1", "SC2", "SC3"}},
        {"HIP", {"BB", "SC1", "SC2", "SC3"}},
        {"HISE", {"BB", "SC1", "SC2", "SC3"}},
        {"HISH", {"BB", "SC1", "SC2", "SC3"}},
        {"HIH", {"BB", "SC1", "SC2", "SC3"}},
        {"LYS", {"BB", "SC1", "SC2"}},
        {"LSN", {"BB", "SC1", "SC2"}},
        {"LYN", {"BB", "SC1", "SC2"}},
        {"LYS0", {"BB", "SC1", "SC2"}},
        {"ASP", {"BB", "SC1"}},
        {"ASH", {"BB", "SC1"}},
        {"ASPP", {"BB", "SC1"}},
        {"ASP0", {"BB", "SC1"}},
        {"GLU", {"BB", "SC1"}},
        {"GLH", {"BB", "SC1"}},
        {"GLUP", {"BB", "SC1"}},
        {"GLU0", {"BB", "SC1"}},
        {"ARG", {"BB", "SC1", "SC2"}},
        {"ARG0", {"BB", "SC1", "SC2"}},
        {"CYS", {"BB", "SC1"}},
        {"CYX", {"BB", "SC1"}},
        {"CYSF", {"BB", "SC1", "C1", "C2", "C3", "C4", "C5", "C6"}},
        {"CYSG", {"BB", "SC1", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"}},
        {"CYSP", {"BB", "SC1", "C1", "C2", "C3", "C4"}},
        {"GLY", {"BB"}},
        {"GLYM", {"BB", "C1", "C2", "C3", "C4"}},
        {"ALA", {"BB", "SC1"}},
        {"PRO", {"BB", "SC1"}},
        {"ASN", {"BB", "SC1", "SCP", "SCN"}},
        {"LYS", {"BB", "SC1", "SC2", "SCP"}},
    };
    const int n = (int)MD_COUNTOF(res);

    int num_beads = 0;
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < (int)MD_COUNTOF(res[i].beads) && res[i].beads[k]; ++k) num_beads += 1;
    }
    char text[16384];
    int len = snprintf(text, sizeof(text), "Martini residue names\n%5d\n", num_beads);
    int atom_nr = 0;
    for (int i = 0; i < n; ++i) {
        const float x = 1.0f + 0.35f * (float)i;
        for (int k = 0; k < (int)MD_COUNTOF(res[i].beads) && res[i].beads[k]; ++k) {
            len += snprintf(text + len, sizeof(text) - len, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", i + 1, res[i].name, res[i].beads[k], ++atom_nr, x, 1.0f + 0.3f * (float)k, 1.0f);
        }
    }
    len += snprintf(text + len, sizeof(text) - len, "  20.00000  20.00000  20.00000\n");
    ASSERT_LT(len, (int)sizeof(text));

    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, (md_str_t){text, (size_t)len}, NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    ASSERT_EQ((size_t)n, md_system_residue_count(&sys));
    check_martini_peptide(utest_result, &sys, &state, (size_t)n);
    md_vm_arena_destroy(alloc);
}

// The bonds of the martinize2 topology of the peptide above (its CONECT records), as atom index pairs
static const int martini3_peptide_topology[][2] = {
    {0, 1}, {0, 5}, {1, 2}, {1, 3}, {2, 3}, {2, 4}, {3, 4}, {5, 6}, {5, 7}, {7, 8}, {8, 9}, {8, 10}, {10, 11}, {10, 12},
    {12, 13}, {12, 14}, {14, 15}, {14, 18}, {15, 16}, {15, 17}, {16, 17}, {18, 19}, {18, 20}, {20, 21}, {21, 22}, {21, 24},
    {22, 24}, {22, 25}, {24, 25},
};

static bool pair_in(const int (*pairs)[2], size_t num_pairs, int a, int b) {
    for (size_t i = 0; i < num_pairs; ++i) {
        if ((pairs[i][0] == a && pairs[i][1] == b) || (pairs[i][0] == b && pairs[i][1] == a)) return true;
    }
    return false;
}

// Every inferred bond is one of the topology's, and the topology's bonds it leaves out are those named
static void check_bonds_against_topology(int* utest_result, const md_system_t* sys, const int (*topology)[2], size_t num_topology, const int (*missing)[2], size_t num_missing) {
    for (size_t i = 0; i < sys->bond.count; ++i) {
        const int a = sys->bond.pairs[i].idx[0];
        const int b = sys->bond.pairs[i].idx[1];
        EXPECT_TRUE_MSG(pair_in(topology, num_topology, a, b), "a bond the topology does not have");
    }
    for (size_t i = 0; i < num_topology; ++i) {
        const bool found = md_bond_find(&sys->bond, topology[i][0], topology[i][1]) != -1;
        EXPECT_EQ((int)!pair_in(missing, num_missing, topology[i][0], topology[i][1]), (int)found);
    }
}

// Inferred from the beads alone, the bonds of the peptide are the topology's: those of each residue (the 4.7 A SC2-SC4
// of the Trp ring too) and the link from each BB to the next. Trp's SC3, a virtual site, is bonded to SC5, the first
// bead it is constructed from, as md_tpr bonds it (martinize2 writes no CONECT record for it). With the CONECT records,
// the bonds are theirs, and inference adds that one.
UTEST(martini, martini3_peptide_bonds) {
    int topology[MD_COUNTOF(martini3_peptide_topology) + 1][2];
    md_memcpy(topology, martini3_peptide_topology, sizeof(martini3_peptide_topology));
    topology[MD_COUNTOF(martini3_peptide_topology)][0] = 23;
    topology[MD_COUNTOF(martini3_peptide_topology)][1] = 25;

    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    {
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(martini3_peptide_pdb), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        check_bonds_against_topology(utest_result, &sys, (const int (*)[2])topology, MD_COUNTOF(topology), NULL, 0);
        EXPECT_EQ(1u, md_bond_conn_count(&sys.bond, 23));
    }
    {
        // The same with its CONECT records, as martinize2 wrote it
        char text[8192];
        int len = snprintf(text, sizeof(text), "%s", martini3_peptide_pdb);
        len -= (int)strlen("END\n");
        for (size_t i = 0; i < MD_COUNTOF(martini3_peptide_topology); ++i) {
            len += snprintf(text + len, sizeof(text) - len, "CONECT%5d%5d\n", martini3_peptide_topology[i][0] + 1, martini3_peptide_topology[i][1] + 1);
        }
        len += snprintf(text + len, sizeof(text) - len, "END\n");
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, (md_str_t){text, (size_t)len}, NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        check_bonds_against_topology(utest_result, &sys, (const int (*)[2])topology, MD_COUNTOF(topology), NULL, 0);
        EXPECT_EQ(MD_COUNTOF(topology), sys.bond.count);
        for (size_t i = 0; i < sys.bond.count; ++i) {
            const bool vsite = pair_in((const int (*)[2])topology + MD_COUNTOF(martini3_peptide_topology), 1, sys.bond.pairs[i].idx[0], sys.bond.pairs[i].idx[1]);
            EXPECT_EQ(vsite ? MD_BOND_ORIGIN_INFERRED : MD_BOND_ORIGIN_FILE, md_bond_origin(sys.bond.flags[i]));
        }
    }
    md_vm_arena_destroy(alloc);
}

// Polarizable Martini 2 (martinize2 -ff martini22p on villin, residues 53-57): the polar SC1 of Thr and Ser is a
// massless virtual site between two charge dummies, Arg's SC2 carries one. The bonds are the topology's as GROMACS has
// them, all of them: the virtual site is bonded to SCP, the first bead it is constructed from, and the dummies to each
// other (martinize2's CONECT records bond SCN to SC1 as well).
static const char martini22p_peptide_pdb[] =
    "ATOM      1 BB   MET A  53       9.910  -7.189  67.361  1.00  0.00\n"
    "ATOM      2 SC1  MET A  53       7.345  -7.793  70.918  1.00  0.00\n"
    "ATOM      3 BB   THR A  54       7.578  -5.503  65.552  1.00  0.00\n"
    "ATOM      4 SC1  THR A  54       6.937  -7.801  64.213  1.00  0.00\n"
    "ATOM      5 SCP  THR A  54       6.589  -7.804  63.606  1.00  0.00\n"
    "ATOM      6 SCN  THR A  54       7.285  -7.798  64.820  1.00  0.00\n"
    "ATOM      7 BB   ARG A  55       6.195  -2.964  64.028  1.00  0.00\n"
    "ATOM      8 SC1  ARG A  55       7.904  -0.987  62.207  1.00  0.00\n"
    "ATOM      9 SC2  ARG A  55       8.217  -1.441  59.751  1.00  0.00\n"
    "ATOM     10 SCP  ARG A  55       9.317  -1.441  59.751  1.00  0.00\n"
    "ATOM     11 BB   SER A  56       3.400  -3.710  64.222  1.00  0.00\n"
    "ATOM     12 SC1  SER A  56       2.364  -4.957  62.032  1.00  0.00\n"
    "ATOM     13 SCP  SER A  56       1.924  -5.331  61.636  1.00  0.00\n"
    "ATOM     14 SCN  SER A  56       2.803  -4.583  62.428  1.00  0.00\n"
    "ATOM     15 BB   ALA A  57       3.595  -4.775  67.236  1.00  0.00\n"
    "END\n";

UTEST(martini, martini22p_peptide) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(martini22p_peptide_pdb), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));

    static const int topology[][2] = {
        {0, 1}, {0, 2}, {2, 3}, {2, 6}, {3, 4}, {4, 5}, {6, 7}, {6, 10}, {7, 8}, {8, 9}, {10, 11}, {10, 14}, {11, 12}, {12, 13},
    };
    check_bonds_against_topology(utest_result, &sys, topology, MD_COUNTOF(topology), NULL, 0);

    // The masses of polarizable Martini 2
    EXPECT_EQ(72.0f, md_atom_mass(&sys.atom, 1));     // Met SC1
    EXPECT_EQ( 0.0f, md_atom_mass(&sys.atom, 3));     // Thr SC1, a virtual site
    EXPECT_EQ(36.0f, md_atom_mass(&sys.atom, 4));     // A charge dummy
    EXPECT_EQ(36.0f, md_atom_mass(&sys.atom, 8));     // Arg SC2, charged
    EXPECT_NEAR(1.0f, md_atom_radius(&sys.atom, 4), 1e-4f);
    for (size_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(MD_PARTICLE_BEAD, md_atom_particle_kind(&sys.atom, i));
    }
    check_martini_peptide(utest_result, &sys, &state, 5);
    md_vm_arena_destroy(alloc);
}

// The masses and sizes of Martini beads are those of the version of Martini the system is from: the same names stand
// for regular, small and tiny beads in Martini 3, for beads of 72 in Martini 2
UTEST(martini, martini_masses_by_version) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    {
        // Martini 3 (Tyr has SC4 and Trp SC5)
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, md_str_from_cstr(martini3_peptide_gro), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        EXPECT_EQ(72.0f, md_atom_mass(&sys.atom, 0));     // Tyr BB, regular
        EXPECT_EQ(36.0f, md_atom_mass(&sys.atom, 4));     // Tyr SC4, tiny
        EXPECT_EQ(54.0f, md_atom_mass(&sys.atom, 5));     // Val BB, small
        EXPECT_EQ(54.0f, md_atom_mass(&sys.atom, 7));     // Gly BB
        EXPECT_EQ(54.0f, md_atom_mass(&sys.atom, 11));    // Leu SC1
        EXPECT_EQ(72.0f, md_atom_mass(&sys.atom, 9));     // Glu SC1
        EXPECT_EQ(36.0f, md_atom_mass(&sys.atom, 16));    // His SC2
        EXPECT_EQ( 0.0f, md_atom_mass(&sys.atom, 23));    // Trp SC3, a virtual site
        EXPECT_EQ(36.0f, md_atom_mass(&sys.atom, 25));    // Trp SC5
        EXPECT_NEAR(2.64f, md_atom_radius(&sys.atom, 0), 0.01f);
        EXPECT_NEAR(2.30f, md_atom_radius(&sys.atom, 11), 0.01f);
        EXPECT_NEAR(1.91f, md_atom_radius(&sys.atom, 23), 0.01f);
        // The SC1 beads of Tyr (tiny) and Leu (small) are different types
        EXPECT_NE(sys.atom.type_idx[1], sys.atom.type_idx[11]);
    }
    {
        // Martini 2 (Trp has SC1-SC4, Ala no side chain bead): 72, but 45 for the small beads of the ring
        const char* text =
            "Martini 2\n"
            "    9\n"
            "    1ALA     BB    1   1.000   1.000   1.000\n"
            "    2TRP     BB    2   1.350   1.000   1.000\n"
            "    2TRP    SC1    3   1.350   1.300   1.000\n"
            "    2TRP    SC2    4   1.350   1.600   1.000\n"
            "    2TRP    SC3    5   1.650   1.600   1.000\n"
            "    2TRP    SC4    6   1.650   1.300   1.000\n"
            "    3GLY     BB    7   1.700   1.000   1.000\n"
            "    4LEU     BB    8   2.050   1.000   1.000\n"
            "    4LEU    SC1    9   2.050   1.300   1.000\n"
            "   5.00000   5.00000   5.00000\n";
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, md_str_from_cstr(text), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        for (size_t i = 0; i < sys.atom.count; ++i) {
            const bool ring = 2 <= i && i <= 5;
            EXPECT_EQ(ring ? 45.0f : 72.0f, md_atom_mass(&sys.atom, i));
        }
        EXPECT_NEAR(2.64f, md_atom_radius(&sys.atom, 6), 0.01f);   // Gly BB
        EXPECT_NEAR(2.41f, md_atom_radius(&sys.atom, 4), 0.01f);   // Trp SC3
        EXPECT_NEAR(2.64f, md_atom_radius(&sys.atom, 8), 0.01f);   // Leu SC1
    }
    md_vm_arena_destroy(alloc);
}

// Martini 2 DNA (Uusitalo et al. 2015), mapped from four nucleotides of a strand of the nucleosome in nucl-dna.gro
// (TGCA, -60 to -57): the phosphate BB1 and sugar BB2, BB3 are the backbone, the rest the base. Each bead is the centre
// of its atoms: BB1 P OP1 OP2 O5', BB2 C5' O4' C4', BB3 C3' C2' C1'; SC1 N9 C4, SC2 C2 N3 (and N2 of G), SC3 C6 N1 and
// N6 (O6 of G), SC4 C8 N7 C5 of a purine; SC1 N1 C6, SC2 N3 C2 O2, SC3 C5 C4 and N4 (O4 C7 of T) of a pyrimidine.
// Within a residue their distances are the bond and constraint lengths of the topology within 0.03 nm. The bonds are
// those of the topology (polyply's DNA_M2), BB3 to the next BB1 included.
static const char martini2_dna_pdb[] =
    "ATOM      1 BB1  DT  A -60      -3.350 101.668  23.203  1.00  0.00\n"
    "ATOM      2 BB2  DT  A -60      -2.020 101.387  20.003  1.00  0.00\n"
    "ATOM      3 BB3  DT  A -60      -1.460  99.483  19.917  1.00  0.00\n"
    "ATOM      4 SC1  DT  A -60       0.850  99.875  21.345  1.00  0.00\n"
    "ATOM      5 SC2  DT  A -60       2.520  99.690  19.910  1.00  0.00\n"
    "ATOM      6 SC3  DT  A -60       2.758  99.835  23.080  1.00  0.00\n"
    "ATOM      7 BB1  DG  A -59      -3.655  97.578  17.750  1.00  0.00\n"
    "ATOM      8 BB2  DG  A -59      -1.057  97.917  15.617  1.00  0.00\n"
    "ATOM      9 BB3  DG  A -59      -0.357  96.060  15.630  1.00  0.00\n"
    "ATOM     10 SC1  DG  A -59       1.545  96.355  18.015  1.00  0.00\n"
    "ATOM     11 SC2  DG  A -59       4.450  96.427  17.973  1.00  0.00\n"
    "ATOM     12 SC3  DG  A -59       3.530  96.343  20.680  1.00  0.00\n"
    "ATOM     13 SC4  DG  A -59       0.887  96.270  19.463  1.00  0.00\n"
    "ATOM     14 BB1  DC  A -58      -0.695  94.367  12.558  1.00  0.00\n"
    "ATOM     15 BB2  DC  A -58       2.600  93.820  13.640  1.00  0.00\n"
    "ATOM     16 BB3  DC  A -58       2.597  91.947  14.360  1.00  0.00\n"
    "ATOM     17 SC1  DC  A -58       2.015  92.455  17.020  1.00  0.00\n"
    "ATOM     18 SC2  DC  A -58       3.943  92.957  18.023  1.00  0.00\n"
    "ATOM     19 SC3  DC  A -58       1.660  92.757  19.563  1.00  0.00\n"
    "ATOM     20 BB1  DA  A -57       3.948  90.553  10.992  1.00  0.00\n"
    "ATOM     21 BB2  DA  A -57       6.963  91.600  12.027  1.00  0.00\n"
    "ATOM     22 BB3  DA  A -57       7.717  89.830  12.417  1.00  0.00\n"
    "ATOM     23 SC1  DA  A -57       6.215  89.725  15.085  1.00  0.00\n"
    "ATOM     24 SC2  DA  A -57       6.710  90.325  17.270  1.00  0.00\n"
    "ATOM     25 SC3  DA  A -57       4.267  89.427  17.783  1.00  0.00\n"
    "ATOM     26 SC4  DA  A -57       4.760  89.077  14.940  1.00  0.00\n"
    "END\n";

UTEST(martini, martini2_dna) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(martini2_dna_pdb), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    ASSERT_EQ(4u, md_system_residue_count(&sys));

    for (size_t c = 0; c < 4; ++c) {
        EXPECT_EQ(MD_RESIDUE_KIND_NUCLEOTIDE, md_system_residue_kind(&sys, c));
        const md_urange_t range = md_system_residue_atom_range(&sys, c);
        for (uint32_t i = range.beg; i < range.end; ++i) {
            const md_str_t name = md_atom_name(&sys.atom, i);
            const md_atom_flags_t flags = md_atom_flags(&sys.atom, i);
            EXPECT_EQ(MD_PARTICLE_BEAD, md_atom_particle_kind(&sys.atom, i));
            EXPECT_EQ(0, md_atom_atomic_number(&sys.atom, i));
            if (md_str_eq(name, MD_STR_LIT("BB1"))) {
                EXPECT_EQ((int)MD_ATOM_FLAG_BACKBONE, (int)flags);
                EXPECT_EQ(72.0f, md_atom_mass(&sys.atom, i));
            } else if (md_str_eq(name, MD_STR_LIT("BB2")) || md_str_eq(name, MD_STR_LIT("BB3"))) {
                EXPECT_EQ((int)(MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_NUCLEOSIDE), (int)flags);
                EXPECT_EQ(45.0f, md_atom_mass(&sys.atom, i));
            } else {
                EXPECT_EQ((int)(MD_ATOM_FLAG_NUCLEOBASE | MD_ATOM_FLAG_NUCLEOSIDE), (int)flags);
                EXPECT_EQ(45.0f, md_atom_mass(&sys.atom, i));
            }
        }
    }

    static const int topology[][2] = {
        {0, 1}, {1, 2}, {2, 3}, {2, 6}, {3, 4}, {3, 5}, {4, 5}, {6, 7}, {7, 8}, {8, 9}, {8, 13}, {9, 10}, {9, 12}, {10, 11},
        {10, 12}, {11, 12}, {13, 14}, {14, 15}, {15, 16}, {15, 19}, {16, 17}, {16, 18}, {17, 18}, {19, 20}, {20, 21},
        {21, 22}, {22, 23}, {22, 25}, {23, 24}, {23, 25}, {24, 25},
    };
    check_bonds_against_topology(utest_result, &sys, topology, MD_COUNTOF(topology), NULL, 0);

    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ(1u, md_system_part_count(&sys));
    EXPECT_EQ(MD_PART_KIND_DNA, md_system_part_kind(&sys, 0));
    md_vm_arena_destroy(alloc);
}

// Martini 2 DNA likewise: a second strand under the chain id of the first, numbered on, begins without its phosphate
// bead. The chain ends there: two molecules of one chain, each strand whole.
static const char martini2_dna_second_strand_pdb[] =
    "ATOM     27 BB2  DT  A -56      37.980 101.387  20.003  1.00  0.00\n"
    "ATOM     28 BB3  DT  A -56      38.540  99.483  19.917  1.00  0.00\n"
    "ATOM     29 SC1  DT  A -56      40.850  99.875  21.345  1.00  0.00\n"
    "ATOM     30 SC2  DT  A -56      42.520  99.690  19.910  1.00  0.00\n"
    "ATOM     31 SC3  DT  A -56      42.758  99.835  23.080  1.00  0.00\n"
    "ATOM     32 BB1  DG  A -55      36.345  97.578  17.750  1.00  0.00\n"
    "ATOM     33 BB2  DG  A -55      38.943  97.917  15.617  1.00  0.00\n"
    "ATOM     34 BB3  DG  A -55      39.643  96.060  15.630  1.00  0.00\n"
    "ATOM     35 SC1  DG  A -55      41.545  96.355  18.015  1.00  0.00\n"
    "ATOM     36 SC2  DG  A -55      44.450  96.427  17.973  1.00  0.00\n"
    "ATOM     37 SC3  DG  A -55      43.530  96.343  20.680  1.00  0.00\n"
    "ATOM     38 SC4  DG  A -55      40.887  96.270  19.463  1.00  0.00\n"
    "ATOM     39 BB1  DC  A -54      39.305  94.367  12.558  1.00  0.00\n"
    "ATOM     40 BB2  DC  A -54      42.600  93.820  13.640  1.00  0.00\n"
    "ATOM     41 BB3  DC  A -54      42.597  91.947  14.360  1.00  0.00\n"
    "ATOM     42 SC1  DC  A -54      42.015  92.455  17.020  1.00  0.00\n"
    "ATOM     43 SC2  DC  A -54      43.943  92.957  18.023  1.00  0.00\n"
    "ATOM     44 SC3  DC  A -54      41.660  92.757  19.563  1.00  0.00\n"
    "ATOM     45 BB1  DA  A -53      43.948  90.553  10.992  1.00  0.00\n"
    "ATOM     46 BB2  DA  A -53      46.963  91.600  12.027  1.00  0.00\n"
    "ATOM     47 BB3  DA  A -53      47.717  89.830  12.417  1.00  0.00\n"
    "ATOM     48 SC1  DA  A -53      46.215  89.725  15.085  1.00  0.00\n"
    "ATOM     49 SC2  DA  A -53      46.710  90.325  17.270  1.00  0.00\n"
    "ATOM     50 SC3  DA  A -53      44.267  89.427  17.783  1.00  0.00\n"
    "ATOM     51 SC4  DA  A -53      44.760  89.077  14.940  1.00  0.00\n"
    "END\n";

UTEST(martini, martini2_dna_strands_under_one_chain_id) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    char text[sizeof(martini2_dna_pdb) + sizeof(martini2_dna_second_strand_pdb)];
    const size_t first_len = sizeof(martini2_dna_pdb) - 1 - 4;     // without its END line
    memcpy(text, martini2_dna_pdb, first_len);
    memcpy(text + first_len, martini2_dna_second_strand_pdb, sizeof(martini2_dna_second_strand_pdb));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(text), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ(8u, md_system_residue_count(&sys));
    ASSERT_EQ(51u, (unsigned)sys.atom.count);
    EXPECT_TRUE(md_system_is_coarse_grained(&sys));
    EXPECT_EQ(2u, md_system_part_count(&sys));
    EXPECT_EQ(1u, md_system_chain_count(&sys));
    EXPECT_EQ(51u, md_system_chain_atom_count(&sys, 0));
    EXPECT_EQ(2u, md_system_molecule_count(&sys));
    for (size_t i = 0; i < md_array_size(sys.molecule.atom_idx); ++i) {
        EXPECT_TRUE(0 <= sys.molecule.atom_idx[i] && sys.molecule.atom_idx[i] < 51);
    }
    check_molecule_invariants(utest_result, &sys);
    md_vm_arena_destroy(alloc);
}

// Between residues: the numbering may jump forward over a bond (an engineered insert) but starting over is a new chain,
// however close its first residue is. Two cysteines of different chains 2.4 A apart are bridged.
UTEST(martini, martini_chain_links) {
    const char* text =
        "Chains\n"
        "   10\n"
        "    1ALA     BB    1   1.000   1.000   1.000\n"
        "    1ALA    SC1    2   1.000   0.750   1.000\n"
        "    2CYS     BB    3   1.350   1.000   1.000\n"
        "    2CYS    SC1    4   1.350   1.100   1.200\n"
        "  921ALA     BB    5   1.700   1.000   1.000\n"
        "  921ALA    SC1    6   1.700   0.750   1.000\n"
        "    1ALA     BB    7   1.700   1.350   1.000\n"
        "    1ALA    SC1    8   1.950   1.350   1.000\n"
        "    2CYS     BB    9   1.350   1.350   1.000\n"
        "    2CYS    SC1   10   1.350   1.250   1.380\n"
        "   5.00000   5.00000   5.00000\n";
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, md_str_from_cstr(text), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    static const int expected[][2] = {
        {0, 1}, {2, 3}, {4, 5}, {6, 7}, {8, 9},     // BB-SC1
        {0, 2}, {2, 4},                             // 1-2, and 2-921 over the jump
        {6, 8},                                     // the second chain
        {3, 9},                                     // the disulfide
    };
    check_bonds_against_topology(utest_result, &sys, expected, MD_COUNTOF(expected), NULL, 0);
    EXPECT_EQ(-1, md_bond_find(&sys.bond, 4, 6));  // 3.5 A, but the numbering starts over
    md_vm_arena_destroy(alloc);
}

// A residue whose beads are not those of a topology (here Trp of Martini 3 written without its virtual site SC3) is bonded
// by its geometry: the tree of its beads and the ring bonds besides, but for the 4.7 A SC2-SC4 of the Trp ring. Between
// residues numbered one after the other a link may be stretched to 5.5 A, as in a simulation; over a jump in the
// numbering it is one up to 4.3 A.
UTEST(martini, martini_bonds_by_geometry) {
    const char* text =
        "Geometry\n"
        "   11\n"
        "    9TRP     BB    1   2.304   1.959   1.000\n"
        "    9TRP    SC1    2   2.470   1.827   1.121\n"
        "    9TRP    SC2    3   2.679   1.768   1.005\n"
        "    9TRP    SC4    4   2.709   2.018   1.298\n"
        "    9TRP    SC5    5   2.890   1.973   1.194\n"
        "   10ALA     BB    6   2.304   1.959   1.500\n"
        "   10ALA    SC1    7   2.104   1.959   1.500\n"
        "   11ALA     BB    8   2.304   1.959   1.950\n"
        "   11ALA    SC1    9   2.104   1.959   1.950\n"
        "   20ALA     BB   10   2.304   1.959   2.410\n"
        "   20ALA    SC1   11   2.104   1.959   2.410\n"
        "   5.00000   5.00000   5.00000\n";
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, md_str_from_cstr(text), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    static const int expected[][2] = {
        {0, 1}, {1, 2}, {2, 4}, {3, 4}, {1, 3},     // Trp, without SC2-SC4 (2-3)
        {5, 6}, {7, 8}, {9, 10},                    // Ala BB-SC1
        {0, 5},                                     // 5.0 A, numbered next
        {5, 7},                                     // 4.5 A, numbered next
    };
    check_bonds_against_topology(utest_result, &sys, expected, MD_COUNTOF(expected), NULL, 0);
    EXPECT_EQ(-1, md_bond_find(&sys.bond, 7, 9));  // 4.6 A over a jump in the numbering
    md_vm_arena_destroy(alloc);
}

// Martini systems simulated with GROMACS (test_data/martini, see its README): a frame of each as .gro, and its topology
// as .tpr. Read from the frame alone, every bead gets the mass, size, particle kind and role the topology gives it,
// every residue its kind, and the bonds are those of the topology, no more and no fewer. The frames are of a running
// simulation: bonds stretched, molecules broken across the box.
UTEST(martini, martini_reference_systems) {
    static const char* names[] = {
        "martini3_membrane",            // POPC, POPE, POPS, CHOL (virtual sites), W, NA, CL
        "martini3_lipids",              // 24 lipids of the Martini 3 lipidome: phospholipids, sphingolipids, glycerides, ...
        "martini3_peptides",            // Charged termini (Pro, Ala), neutral termini (Gly), a disulfide
        "martini3_rna",                 // An RNA duplex (reForge), virtual sites in the bases
        "martini3_small_molecules",     // Small molecules (rings of virtual sites), solvents, sugars, ions, SW, TW, W
        "martini3_polymer",             // PEO with hydroxyl ends and polystyrene (polyply): residues linked by the force field
        "martini2_membrane_protein",    // Ubiquitin on DPPC, DOPC, CHOL; W, WF, NA+, CL-
        "martini2_lipids",              // 24 lipids of the Martini 2 lipidome, surfactants, sterols
        "martini2p_protein",            // Polarizable Martini 2: villin, PW, charge dummies, virtual sites
        "martini2_dna",                 // A DNA duplex (martinize-dna), its elastic network left out
        "martini2_rna",                 // An RNA duplex (martinize-nucleotide)
        "martini2_small_molecules",     // The solvents, sugars and ions of Martini 2
    };
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    for (size_t n = 0; n < MD_COUNTOF(names); ++n) {
        char tpr_path[256], gro_path[256];
        snprintf(tpr_path, sizeof(tpr_path), MD_UNITTEST_DATA_DIR "/martini/%s.tpr", names[n]);
        snprintf(gro_path, sizeof(gro_path), MD_UNITTEST_DATA_DIR "/martini/%s.gro", names[n]);
        md_system_t top = { .alloc = alloc }, gro = { .alloc = alloc };
        md_system_state_t top_state = { .alloc = alloc }, gro_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&top, &top_state, md_str_from_cstr(tpr_path), NULL));
        ASSERT_TRUE(md_gro_system_init_from_file(&gro, &gro_state, md_str_from_cstr(gro_path), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&gro, &gro_state));
        ASSERT_EQ(top.atom.count, gro.atom.count);
        ASSERT_EQ(top.residue.count, gro.residue.count);

        size_t bad_atoms = 0;
        for (size_t i = 0; i < top.atom.count; ++i) {
            const bool same = fabsf(md_atom_mass(&top.atom, i) - md_atom_mass(&gro.atom, i)) < 1e-3f &&
                fabsf(md_atom_radius(&top.atom, i) - md_atom_radius(&gro.atom, i)) < 0.02f &&
                md_atom_particle_kind(&top.atom, i) == md_atom_particle_kind(&gro.atom, i) &&
                md_atom_flags(&top.atom, i) == md_atom_flags(&gro.atom, i) &&
                md_atom_atomic_number(&gro.atom, i) == 0;
            bad_atoms += !same;
        }
        EXPECT_EQ_MSG((size_t)0, bad_atoms, names[n]);

        size_t bad_residues = 0;
        for (size_t c = 0; c < top.residue.count; ++c) {
            bad_residues += md_system_residue_kind(&top, c) != md_system_residue_kind(&gro, c);
        }
        EXPECT_EQ_MSG((size_t)0, bad_residues, names[n]);

        size_t missing = 0, extra = 0;
        for (size_t b = 0; b < top.bond.count; ++b) {
            missing += md_bond_find(&gro.bond, top.bond.pairs[b].idx[0], top.bond.pairs[b].idx[1]) == -1;
        }
        for (size_t b = 0; b < gro.bond.count; ++b) {
            extra += md_bond_find(&top.bond, gro.bond.pairs[b].idx[0], gro.bond.pairs[b].idx[1]) == -1;
        }
        EXPECT_EQ_MSG((size_t)0, missing, names[n]);
        EXPECT_EQ_MSG((size_t)0, extra, names[n]);
        EXPECT_GT(top.bond.count, (size_t)0);
    }
    md_vm_arena_destroy(alloc);
}

// The nucleic acids of Martini (test_data/martini): each strand a coarse grained nucleic backbone of all its
// nucleotides, read from the frame alone or with its topology. A segment holds the phosphate bead BB1 in p (none at the
// 5' end, which starts a strand) and the sugar beads BB2 and BB3 in c4 and c1, the C4' and C1' they hold; no breaks.
UTEST(martini, martini_nucleic_backbones) {
    static const struct {
        const char* name;
        size_t      strand;     // Nucleotides of each of its two strands
        md_part_kind_t kind;
    } systems[] = {
        { "martini2_dna", 12, MD_PART_KIND_DNA },
        { "martini2_rna", 14, MD_PART_KIND_RNA },
        { "martini3_rna", 22, MD_PART_KIND_RNA },
    };
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    for (size_t n = 0; n < MD_COUNTOF(systems); ++n) {
        for (int with_tpr = 0; with_tpr < 2; ++with_tpr) {
            char path[256];
            snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/martini/%s.%s", systems[n].name, with_tpr ? "tpr" : "gro");
            md_system_t sys = { .alloc = alloc };
            md_system_state_t state = { .alloc = alloc };
            ASSERT_TRUE(with_tpr ? md_tpr_system_init_from_file(&sys, &state, md_str_from_cstr(path), NULL) : md_gro_system_init_from_file(&sys, &state, md_str_from_cstr(path), NULL));
            ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));

            const md_nucleic_backbone_data_t* bb = &sys.nucleic_backbone;
            ASSERT_EQ_MSG((size_t)2, bb->range.count, path);
            ASSERT_EQ_MSG(2 * systems[n].strand, bb->segment.count, path);
            for (size_t r = 0; r < bb->range.count; ++r) {
                EXPECT_TRUE(md_nucleic_backbone_range_is_coarse_grained(bb, r));
                EXPECT_EQ(systems[n].strand, (size_t)(bb->range.offset[r + 1] - bb->range.offset[r]));
                EXPECT_TRUE(md_part_kind_is_nucleic_acid(md_system_part_kind(&sys, bb->range.part_idx[r])));
                EXPECT_EQ((int)systems[n].kind, (int)md_system_part_kind(&sys, bb->range.part_idx[r]));
                for (uint32_t s = bb->range.offset[r]; s < bb->range.offset[r + 1]; ++s) {
                    const md_nucleic_acid_atoms_t a = bb->segment.atoms[s];
                    const md_residue_idx_t res = bb->segment.res_idx[s];
                    EXPECT_EQ(-1, a.o5);
                    EXPECT_EQ(-1, a.o3);
                    EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, a.c4), MD_STR_LIT("BB2")));
                    EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, a.c1), MD_STR_LIT("BB3")));
                    EXPECT_EQ(res, (md_residue_idx_t)md_system_residue_find_by_atom_idx(&sys, (size_t)a.c4));
                    if (s == bb->range.offset[r]) {
                        EXPECT_EQ(-1, a.p);
                    } else {
                        ASSERT_GE(a.p, 0);
                        EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, a.p), MD_STR_LIT("BB1")));
                    }
                    EXPECT_EQ(0, (int)(md_system_residue_flags(&sys, (size_t)res) & (MD_RESIDUE_FLAG_BREAK_BEFORE | MD_RESIDUE_FLAG_BREAK_AFTER)));
                }
            }
        }
    }
    md_vm_arena_destroy(alloc);
}

// A strand of Martini 2 DNA (one chain, the first seven nucleotides of martini/martini2_dna but the fourth) is broken
// where the fourth is missing: no bead of the third is bonded to the BB1 of the fifth, 8.8 A away. Two ranges, the
// residues on either side of the break marked as such, and the link over it from the BB3 to the BB1. A nucleotide with
// its BB1 does not start a strand.
static const char martini2_dna_gap_pdb[] =
    "ATOM      1 BB2   DC A   1      33.620  44.890  45.110  1.00  0.00\n"
    "ATOM      2 BB3   DC A   1      33.350  43.570  43.660  1.00  0.00\n"
    "ATOM      3 SC1   DC A   1      31.610  41.910  45.080  1.00  0.00\n"
    "ATOM      4 SC2   DC A   1      32.610  40.180  45.990  1.00  0.00\n"
    "ATOM      5 SC3   DC A   1      29.800  40.440  46.390  1.00  0.00\n"
    "ATOM      6 BB1   DG A   2      36.200  43.640  41.430  1.00  0.00\n"
    "ATOM      7 BB2   DG A   2      36.750  40.590  42.920  1.00  0.00\n"
    "ATOM      8 BB3   DG A   2      36.400  39.240  41.580  1.00  0.00\n"
    "ATOM      9 SC1   DG A   2      33.990  38.330  42.990  1.00  0.00\n"
    "ATOM     10 SC2   DG A   2      33.420  35.860  44.500  1.00  0.00\n"
    "ATOM     11 SC3   DG A   2      30.930  37.410  44.160  1.00  0.00\n"
    "ATOM     12 SC4   DG A   2      32.640  39.140  42.680  1.00  0.00\n"
    "ATOM     13 BB1   DC A   3      39.410  37.770  40.920  1.00  0.00\n"
    "ATOM     14 BB2   DC A   3      37.700  34.940  42.220  1.00  0.00\n"
    "ATOM     15 BB3   DC A   3      36.840  33.840  40.970  1.00  0.00\n"
    "ATOM     16 SC1   DC A   3      34.700  35.320  40.920  1.00  0.00\n"
    "ATOM     17 SC2   DC A   3      33.060  34.200  41.880  1.00  0.00\n"
    "ATOM     18 SC3   DC A   3      32.370  36.510  40.360  1.00  0.00\n"
    "ATOM     19 BB1   DA A   5      35.490  25.450  38.650  1.00  0.00\n"
    "ATOM     20 BB2   DA A   5      32.050  25.760  39.770  1.00  0.00\n"
    "ATOM     21 BB3   DA A   5      30.960  25.540  38.160  1.00  0.00\n"
    "ATOM     22 SC1   DA A   5      30.910  28.390  37.490  1.00  0.00\n"
    "ATOM     23 SC2   DA A   5      29.070  29.740  37.340  1.00  0.00\n"
    "ATOM     24 SC3   DA A   5      30.630  31.410  35.990  1.00  0.00\n"
    "ATOM     25 SC4   DA A   5      32.100  29.000  36.570  1.00  0.00\n"
    "ATOM     26 BB1   DA A   6      29.490  22.800  38.560  1.00  0.00\n"
    "ATOM     27 BB2   DA A   6      26.770  24.910  38.020  1.00  0.00\n"
    "ATOM     28 BB3   DA A   6      26.250  24.770  36.130  1.00  0.00\n"
    "ATOM     29 SC1   DA A   6      27.410  27.350  35.400  1.00  0.00\n"
    "ATOM     30 SC2   DA A   6      26.480  29.330  34.710  1.00  0.00\n"
    "ATOM     31 SC3   DA A   6      28.770  29.740  33.420  1.00  0.00\n"
    "ATOM     32 SC4   DA A   6      28.780  27.100  34.580  1.00  0.00\n"
    "ATOM     33 BB1   DT A   7      23.110  23.480  35.870  1.00  0.00\n"
    "ATOM     34 BB2   DT A   7      22.770  26.750  34.350  1.00  0.00\n"
    "ATOM     35 BB3   DT A   7      22.660  26.780  32.400  1.00  0.00\n"
    "ATOM     36 SC1   DT A   7      25.320  26.770  32.060  1.00  0.00\n"
    "ATOM     37 SC2   DT A   7      25.960  28.690  31.290  1.00  0.00\n"
    "ATOM     38 SC3   DT A   7      27.780  26.040  31.410  1.00  0.00\n"
    "END\n";

UTEST(martini, martini_nucleic_break) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(martini2_dna_gap_pdb), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ((size_t)6, md_system_residue_count(&sys));
    EXPECT_EQ(-1, md_bond_find(&sys.bond, 14, 18));    // BB3 of 3 and BB1 of 5
    EXPECT_NE(-1, md_bond_find(&sys.bond, 7, 12));     // BB3 of 2 and BB1 of 3

    const md_nucleic_backbone_data_t* bb = &sys.nucleic_backbone;
    ASSERT_EQ((size_t)2, bb->range.count);
    EXPECT_EQ((size_t)6, bb->segment.count);
    EXPECT_TRUE(md_nucleic_backbone_range_is_coarse_grained(bb, 0));
    EXPECT_TRUE(md_nucleic_backbone_range_is_coarse_grained(bb, 1));
    EXPECT_EQ(3u, bb->range.offset[1]);
    EXPECT_EQ(-1, bb->segment.atoms[0].p);
    EXPECT_EQ(18, bb->segment.atoms[3].p);
    EXPECT_TRUE(md_system_residue_flags(&sys, 2) & MD_RESIDUE_FLAG_BREAK_AFTER);
    EXPECT_TRUE(md_system_residue_flags(&sys, 3) & MD_RESIDUE_FLAG_BREAK_BEFORE);
    md_atom_idx_t link[2];
    ASSERT_TRUE(md_system_chain_break_atoms(link, &sys, 2));
    EXPECT_EQ(14, link[0]);
    EXPECT_EQ(18, link[1]);
    EXPECT_FALSE(md_system_chain_break_atoms(link, &sys, 1));
    md_vm_arena_destroy(alloc);
}
