#include "utest.h"

#include <md_io.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_xyz.h>
#include <md_mmcif.h>
#include <md_tpr.h>
#include <md_lammps.h>
#include <md_system.h>
#include <md_infer.h>
#include <md_element.h>
#include <core/md_geom.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_path.h>
#include <core/md_bitset.h>

#include <string.h>

#include "system_invariants.h"
#include "system_util.h"

// What md_system_infer gives the system of its parts (MD_SYSTEM_INFER_PART_BIT)
static bool test_infer_parts(md_system_t* sys) {
    md_part_data_t parts = {0};
    const bool ok = md_part_infer(&parts, sys, sys->alloc);
    md_part_data_free(&sys->part, sys->alloc);
    sys->part = parts;
    return ok;
}

UTEST(hierarchy, part) {
    reference_systems_t refsys = reference_systems_load();
    EXPECT_EQ(1,        refsys.mol_ala.part.count);
	EXPECT_EQ(1,        refsys.mol_pftaa.part.count);
	EXPECT_EQ(2,        refsys.mol_nucleotides.part.count);
	EXPECT_EQ(253 + 61, refsys.mol_centered.part.count);

    const md_system_t* sys = &refsys.mol_centered;
    MD_ASSERT(sys->part.count > 253);

    const md_urange_t ref = md_system_part_atom_range(sys, 0);
    for (size_t i = 1; i < 253; ++i) {
        const md_urange_t r = md_system_part_atom_range(sys, i);
        EXPECT_EQ(ref.end - ref.beg, r.end - r.beg);
        EXPECT_EQ(md_system_part_type_idx(sys, 0), md_system_part_type_idx(sys, i));
    }
    // The 253 chains are made alike: one molecule type
    EXPECT_EQ(md_system_molecule_type_idx(sys, 0), md_system_molecule_type_idx(sys, 252));
    reference_systems_free(&refsys);
}

UTEST(hierarchy, molecule_count) {
    reference_systems_t refsys = reference_systems_load();
    size_t num_molecules_pftaa = md_system_molecule_count(&refsys.mol_pftaa);
    EXPECT_EQ(1, num_molecules_pftaa);
    size_t num_molecules_nucleotides = md_system_molecule_count(&refsys.mol_nucleotides);
    EXPECT_EQ(2, num_molecules_nucleotides);
    size_t num_molecules_ala = md_system_molecule_count(&refsys.mol_ala);
	EXPECT_EQ(1, num_molecules_ala);
	size_t num_molecules_centered = md_system_molecule_count(&refsys.mol_centered);
	EXPECT_EQ(253+61, num_molecules_centered);
    reference_systems_free(&refsys);
}

UTEST(hierarchy, molecules) {
    reference_systems_t refsys = reference_systems_load();
    size_t num_molecules = 0;

    num_molecules = md_system_molecule_count(&refsys.mol_nucleotides);
    EXPECT_EQ(num_molecules, 2);

    num_molecules = md_system_molecule_count(&refsys.mol_ala);
	EXPECT_EQ(num_molecules, 1);

	num_molecules = md_system_molecule_count(&refsys.mol_pftaa);
	EXPECT_EQ(num_molecules, 1);

	num_molecules = md_system_molecule_count(&refsys.mol_centered);
    EXPECT_EQ(num_molecules, 253 + 61); // Chains + PFTAA
    reference_systems_free(&refsys);
}

// The molecule hierarchy must be rooted at the graph center and every atom must appear after the
// atom it was reached from. md_geom_unwrap depends on
// that ordering, and atom_slot must be a consistent inverse of atom_idx.
UTEST(hierarchy, molecule_hierarchy) {
    md_allocator_i* alloc = md_get_heap_allocator();

    {
        // A chain of nine has a single center, the middle atom
        md_system_t sys = {0};
        make_chain_system(&sys, alloc, 1, 9);
        ASSERT_TRUE(test_infer_molecules(&sys));

        ASSERT_EQ(md_system_molecule_count(&sys), 1);
        EXPECT_EQ(sys.molecule.atom_idx[0], 4);
        EXPECT_EQ(sys.molecule.parent_idx[0], sys.molecule.atom_idx[0]); // the root is its own parent

        for (size_t slot = 0; slot < 9; ++slot) {
            const int32_t atom   = sys.molecule.atom_idx[slot];
            const int32_t parent = sys.molecule.parent_idx[slot];
            EXPECT_EQ(md_system_molecule_atom_slot(&sys, atom), (int32_t)slot);
            EXPECT_EQ(md_system_molecule_atom_parent(&sys, atom), parent);
            if (parent != atom) {
                // slot(parent) < slot(child), the invariant a single forward pass relies on
                EXPECT_LT(md_system_molecule_atom_slot(&sys, parent), (int32_t)slot);
            }
        }

        md_system_free(&sys);
    }

    {
        // Separate residues are rooted independently and stored back to back
        md_system_t sys = {0};
        make_chain_system(&sys, alloc, 2, 5);
        ASSERT_TRUE(test_infer_molecules(&sys));

        ASSERT_EQ(md_system_molecule_count(&sys), 2);
        EXPECT_EQ(sys.molecule.offset[0], 0u);
        EXPECT_EQ(sys.molecule.offset[1], 5u);
        EXPECT_EQ(sys.molecule.offset[2], 10u);
        EXPECT_EQ(sys.molecule.atom_idx[0], 2);
        EXPECT_EQ(sys.molecule.atom_idx[5], 7);

        md_system_free(&sys);
    }

    {
        // A lone atom and a bonded pair must not trip the double sweep used to locate the center
        md_system_t sys = {0};
        sys.alloc = alloc;
        sys.atom.count = 3;
        md_atom_pair_t pair = {1, 2};
        md_array_push(sys.bond.pairs, pair, alloc);
        sys.bond.count = 1;
        md_bond_build_connectivity(&sys.bond, 3, alloc);
        ASSERT_TRUE(test_infer_molecules(&sys));

        ASSERT_EQ(md_system_molecule_count(&sys), 2);
        for (size_t slot = 0; slot < 3; ++slot) {
            EXPECT_EQ(md_system_molecule_atom_slot(&sys, sys.molecule.atom_idx[slot]), (int32_t)slot);
        }

        md_system_free(&sys);
    }
}

// Coarse grained systems frequently carry no bonds. The molecule hierarchy must then come from the
// residue hierarchy: beads hang off their residue's anchor, consecutive polymer residues are
// chained, non polymer residues stay separate, and nothing is written to sys->bond.
UTEST(hierarchy, molecule_hierarchy_coarse_grained_without_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));

    md_system_t sys = {0};
    sys.alloc = alloc;

    const char* ala[] = {"BB"};
    const char* lys[] = {"BB", "SC1", "SC2"};
    const char* phe[] = {"BB", "SC1", "SC2", "SC3"};
    const char* popc[] = {"NC3", "PO4", "GL1", "GL2", "C1A", "C1B"};
    const char* w[] = {"W"};

    cg_add_residue(&sys, "LYS", 1, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc); // atoms 0-2
    cg_add_residue(&sys, "ALA", 2, MD_RESIDUE_KIND_AMINO_ACID, ala, 1, alloc); // atom      3
    cg_add_residue(&sys, "PHE", 3, MD_RESIDUE_KIND_AMINO_ACID, phe, 4, alloc); // atoms 4-7
    cg_add_residue(&sys, "ALA", 7, MD_RESIDUE_KIND_AMINO_ACID, ala, 1, alloc); // atom      8, sequence gap: new chain
    cg_add_residue(&sys, "ALA", 8, MD_RESIDUE_KIND_AMINO_ACID, ala, 1, alloc); // atom      9
    cg_add_residue(&sys, "POPC", 9, 0, popc, 6, alloc);                  // atoms 10-15
    cg_add_residue(&sys, "POPC", 10, 0, popc, 6, alloc);                 // atoms 16-21, must not join the previous lipid
    cg_add_residue(&sys, "W", 11, MD_RESIDUE_KIND_WATER, w, 1, alloc);             // atom  22

    ASSERT_TRUE(test_infer_molecules(&sys));
    EXPECT_EQ(sys.bond.count, 0u);

    ASSERT_EQ(md_system_molecule_count(&sys), 5u);
    check_molecule_invariants(utest_result, &sys);

    // The first chain is a path of three anchors, so the middle one is the root
    EXPECT_EQ(sys.molecule.offset[1] - sys.molecule.offset[0], 8u);
    EXPECT_EQ(sys.molecule.atom_idx[0], 3);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 0), 3);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 1), 0);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 2), 0);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 4), 3);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 7), 4);

    // Second chain, broken off by the sequence gap
    EXPECT_EQ(sys.molecule.offset[2] - sys.molecule.offset[1], 2u);

    // Each lipid on its own, rooted at its first bead, and the water bead alone
    EXPECT_EQ(sys.molecule.offset[3] - sys.molecule.offset[2], 6u);
    EXPECT_EQ(sys.molecule.offset[4] - sys.molecule.offset[3], 6u);
    EXPECT_EQ(sys.molecule.offset[5] - sys.molecule.offset[4], 1u);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 15), 10);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 21), 16);

    md_vm_arena_destroy(alloc);
}

// Where bonds exist they are used, and links only fill in what they leave disconnected. A system that is
// not coarse grained gets no hierarchy: only what makes a molecule of a residue (synthesize_molecule_links).
UTEST(hierarchy, molecule_hierarchy_coarse_grained_partial_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));

    {
        md_system_t sys = {0};
        sys.alloc = alloc;
        const char* lys[] = {"BB", "SC1", "SC2"};
        cg_add_residue(&sys, "LYS", 1, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);
        cg_add_residue(&sys, "LYS", 2, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);

        // SC1-SC2 bonded within each residue, nothing else
        md_atom_pair_t p0 = {{1, 2}};
        md_atom_pair_t p1 = {{4, 5}};
        md_array_push(sys.bond.pairs, p0, alloc);
        md_array_push(sys.bond.pairs, p1, alloc);
        sys.bond.count = 2;
        md_bond_build_connectivity(&sys.bond, sys.atom.count, alloc);

        ASSERT_TRUE(test_infer_molecules(&sys));
        EXPECT_EQ(sys.bond.count, 2u);
        ASSERT_EQ(md_system_molecule_count(&sys), 1u);
        check_molecule_invariants(utest_result, &sys);

        // SC2 is reached through its bond to SC1, not through a redundant link to BB
        EXPECT_EQ(md_system_molecule_atom_parent(&sys, 2), 1);
        EXPECT_EQ(md_system_molecule_atom_parent(&sys, 5), 4);
    }

    {
        md_system_t sys = {0};
        sys.alloc = alloc;
        const char* lys[] = {"BB", "SC1", "SC2"};
        cg_add_residue(&sys, "LYS", 1, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);
        cg_add_residue(&sys, "LYS", 2, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);
        for (size_t i = 0; i < sys.atom.type.count; ++i) {
            sys.atom.type.flags[i] = md_particle_type_flags_set_particle_kind(sys.atom.type.flags[i], MD_PARTICLE_ATOM);
        }

        // Atoms without bonds: each residue a molecule, its atoms joined to its first. Without instances there are no
        // chains to join the residues.
        ASSERT_TRUE(test_infer_molecules(&sys));
        EXPECT_EQ(md_system_molecule_count(&sys), 2u);
        check_molecule_invariants(utest_result, &sys);
    }

    md_vm_arena_destroy(alloc);
}

// A metal coordinated by residues is a molecule of its own, and so are the waters of its shell: the coordination bonds
// between residues hold no molecule together. Unless a topology gives them, as bonds within one of its molecules.
static size_t molecule_size_of_atom(const md_system_t* sys, int32_t atom);

UTEST(hierarchy, coordinated_metals_are_molecules_of_their_own) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    {
        // The Mg of tubulin, between the phosphates of its GTP and four waters
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        size_t num_mg = 0;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            if (md_atom_atomic_number(&sys.atom, i) != MD_Z_Mg) continue;
            num_mg += 1;
            EXPECT_EQ(1u, molecule_size_of_atom(&sys, (int32_t)i));
        }
        EXPECT_GT(num_mg, 0u);
        for (size_t c = 0; c < sys.residue.count; ++c) {
            if (md_system_residue_kind(&sys, c) != MD_RESIDUE_KIND_WATER) continue;
            const md_urange_t r = md_system_residue_atom_range(&sys, c);
            EXPECT_EQ((size_t)(r.end - r.beg), molecule_size_of_atom(&sys, (int32_t)r.beg));
        }
        check_molecule_invariants(utest_result, &sys);
    }
    {
        // The Mg of 1fez, coordinated by its protein; and then the same bonds as a topology would give them
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        int32_t mg = -1;
        for (size_t i = 0; i < sys.atom.count && mg < 0; ++i) if (md_atom_atomic_number(&sys.atom, i) == MD_Z_Mg) mg = (int32_t)i;
        ASSERT_GE(mg, 0);
        EXPECT_EQ(1u, molecule_size_of_atom(&sys, mg));
        size_t num_coord = 0;
        for (size_t b = 0; b < sys.bond.count; ++b) {
            if (md_bond_kind(sys.bond.flags[b]) != MD_BOND_KIND_COORDINATE) continue;
            sys.bond.flags[b] = md_bond_flags_set_origin(sys.bond.flags[b], MD_BOND_ORIGIN_TOPOLOGY);
            num_coord += 1;
        }
        EXPECT_GT(num_coord, 0u);
        ASSERT_TRUE(test_infer_molecules(&sys));
        EXPECT_LT(1000u, molecule_size_of_atom(&sys, mg));
        check_molecule_invariants(utest_result, &sys);
    }
    md_vm_arena_destroy(alloc);
}

static size_t molecule_size_of_atom(const md_system_t* sys, int32_t atom) {
    return md_system_molecule_atom_count(sys, (size_t)md_system_molecule_find_by_atom_idx(sys, (size_t)atom));
}

// A molecule holds its residues and its polymer chains whole, whatever its bonds say: a chain with residues missing
// from a crystal structure is one molecule, and the waters of a TIP4P model are one each with their virtual sites
UTEST(hierarchy, molecules_hold_residues_and_chains) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));

    // The molecule of an atom, by the slots of the structures
    #define MOLECULE_OF(sys, atom) md_system_molecule_find_by_atom_idx(&(sys), (size_t)(atom))

    {
        // 8g7u: chains with residues missing, each one molecule over its gaps
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        size_t breaks = 0;
        for (size_t i = 0; i < sys.residue.count; ++i) breaks += (md_system_residue_flags(&sys, i) & MD_RESIDUE_FLAG_BREAK_AFTER) != 0;
        EXPECT_GT(breaks, 0u);
        size_t polymers = 0;
        for (size_t i = 0; i < md_system_part_count(&sys); ++i) {
            if (!md_part_kind_is_polymer(md_system_part_kind(&sys, i))) continue;
            polymers += 1;
            const md_urange_t r = md_system_part_atom_range(&sys, i);
            const size_t mol = MOLECULE_OF(sys, r.beg);
            for (uint32_t a = r.beg; a < r.end; ++a) EXPECT_EQ(mol, MOLECULE_OF(sys, a));
        }
        EXPECT_GT(polymers, 1u);
        EXPECT_LE(md_system_molecule_count(&sys), md_system_part_count(&sys));
        check_molecule_invariants(utest_result, &sys);
    }

    {
        // Two DNA strands under one chain id, numbered on: the 5' end of the second has no phosphate, so the chain ends
        // there rather than going on over a gap. Two molecules, one chain.
        static const char pdb[] =
            "ATOM      1  O5'  DA A   1      24.130 139.750  27.200  1.00  0.00           O\n"
            "ATOM      2  C5'  DA A   1      23.020 138.940  27.630  1.00  0.00           C\n"
            "ATOM      3  C4'  DA A   1      21.700 139.500  27.160  1.00  0.00           C\n"
            "ATOM      4  O4'  DA A   1      21.830 139.890  25.780  1.00  0.00           O\n"
            "ATOM      5  C1'  DA A   1      21.020 139.070  24.950  1.00  0.00           C\n"
            "ATOM      6  C3'  DA A   1      20.540 138.510  27.210  1.00  0.00           C\n"
            "ATOM      7  C2'  DA A   1      20.540 137.920  25.820  1.00  0.00           C\n"
            "ATOM      8  O3'  DA A   1      19.300 139.190  27.440  1.00  0.00           O\n"
            "ATOM      9  P    DT A   2      17.930 138.350  27.500  1.00  0.00           P\n"
            "ATOM     10  O5'  DT A   2      17.220 138.660  26.110  1.00  0.00           O\n"
            "ATOM     11  C5'  DT A   2      16.740 139.970  25.830  1.00  0.00           C\n"
            "ATOM     12  C4'  DT A   2      15.860 139.970  24.600  1.00  0.00           C\n"
            "ATOM     13  O4'  DT A   2      16.660 139.700  23.420  1.00  0.00           O\n"
            "ATOM     14  C1'  DT A   2      16.150 138.560  22.730  1.00  0.00           C\n"
            "ATOM     15  C3'  DT A   2      14.730 138.940  24.610  1.00  0.00           C\n"
            "ATOM     16  C2'  DT A   2      15.240 137.830  23.710  1.00  0.00           C\n"
            "ATOM     17  O3'  DT A   2      13.540 139.550  24.110  1.00  0.00           O\n"
            "ATOM     18  O5'  DA A   3      76.380 128.690 -24.150  1.00  0.00           O\n"
            "ATOM     19  C5'  DA A   3      76.750 129.990 -24.620  1.00  0.00           C\n"
            "ATOM     20  C4'  DA A   3      78.030 130.440 -23.960  1.00  0.00           C\n"
            "ATOM     21  O4'  DA A   3      77.720 131.090 -22.700  1.00  0.00           O\n"
            "ATOM     22  C1'  DA A   3      78.420 130.480 -21.630  1.00  0.00           C\n"
            "ATOM     23  C3'  DA A   3      78.960 129.270 -23.630  1.00  0.00           C\n"
            "ATOM     24  C2'  DA A   3      78.820 129.100 -22.130  1.00  0.00           C\n"
            "ATOM     25  O3'  DA A   3      80.320 129.520 -23.980  1.00  0.00           O\n"
            "ATOM     26  P    DT A   4      81.440 128.440 -23.560  1.00  0.00           P\n"
            "ATOM     27  O5'  DT A   4      81.940 128.990 -22.150  1.00  0.00           O\n"
            "ATOM     28  C5'  DT A   4      82.550 130.270 -22.050  1.00  0.00           C\n"
            "ATOM     29  C4'  DT A   4      83.510 130.280 -20.890  1.00  0.00           C\n"
            "ATOM     30  O4'  DT A   4      82.760 130.310 -19.650  1.00  0.00           O\n"
            "ATOM     31  C1'  DT A   4      83.000 129.130 -18.890  1.00  0.00           C\n"
            "ATOM     32  C3'  DT A   4      84.360 129.010 -20.840  1.00  0.00           C\n"
            "ATOM     33  C2'  DT A   4      83.640 128.130 -19.840  1.00  0.00           C\n"
            "ATOM     34  O3'  DT A   4      85.690 129.290 -20.420  1.00  0.00           O\n"
            "END\n";
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(pdb), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        ASSERT_EQ(sys.atom.count, 34u);
        EXPECT_EQ(md_system_instance_count(&sys), 1u);
        EXPECT_EQ(md_system_part_count(&sys), 2u);
        EXPECT_EQ(md_system_chain_count(&sys), 1u);
        EXPECT_EQ(md_system_chain_atom_count(&sys, 0), 34u);
        EXPECT_EQ(md_system_molecule_count(&sys), 2u);
        EXPECT_NE(MOLECULE_OF(sys, 0), MOLECULE_OF(sys, 17));
        EXPECT_EQ(sys.molecule.offset[sys.molecule.count], 34u);
        for (size_t s = 0; s < md_array_size(sys.molecule.atom_idx); ++s) {
            EXPECT_TRUE(0 <= sys.molecule.atom_idx[s] && sys.molecule.atom_idx[s] < 34);
        }
        check_molecule_invariants(utest_result, &sys);
    }

    {
        // Ions and a water given without any bonds: a molecule each, as their residues are
        static const char pdb[] =
            "HETATM    1 NA    NA A   1       0.000   0.000   0.000  1.00  0.00          NA\n"
            "HETATM    2 CL    CL A   2       5.000   0.000   0.000  1.00  0.00          CL\n"
            "HETATM    3  O   HOH A   3      10.000   0.000   0.000  1.00  0.00           O\n"
            "END\n";
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(pdb), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        EXPECT_EQ(md_system_bond_count(&sys), 0u);
        EXPECT_EQ(md_system_molecule_count(&sys), 3u);
        check_molecule_invariants(utest_result, &sys);
    }

    #undef MOLECULE_OF
    md_vm_arena_destroy(alloc);
}

// A branched molecule of no polymer (a glycan) is one instance: its third residue is bonded to the first, not to the
// second before it. A second one, bonded to nothing of the first, is an instance of its own.
UTEST(hierarchy, branched_molecule_is_one_instance) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    static const char pdb[] =
        "HETATM    1  C1  SGR A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  SGR A   1       1.400   0.000   0.000  1.00  0.00           O\n"
        "HETATM    3  C1  SGR A   2       2.800   0.000   0.000  1.00  0.00           C\n"
        "HETATM    4  O2  SGR A   2       4.200   0.000   0.000  1.00  0.00           O\n"
        "HETATM    5  C1  SGR A   3       0.700   1.212   0.000  1.00  0.00           C\n"
        "HETATM    6  O3  SGR A   3       0.700   2.612   0.000  1.00  0.00           O\n"
        "HETATM    7  C1  SGR A   4      20.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    8  O1  SGR A   4      21.400   0.000   0.000  1.00  0.00           O\n"
        "HETATM    9  C1  SGR A   5      22.800   0.000   0.000  1.00  0.00           C\n"
        "HETATM   10  O2  SGR A   5      24.200   0.000   0.000  1.00  0.00           O\n"
        "END\n";
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(pdb), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ(5u, md_system_residue_count(&sys));
    // The branch: the third bonded to the first
    EXPECT_TRUE(md_bond_find(&sys.bond, 0, 4) >= 0 || md_bond_find(&sys.bond, 1, 4) >= 0);
    EXPECT_LT(md_bond_find(&sys.bond, 2, 4), 0);
    ASSERT_EQ(2u, md_system_part_count(&sys));
    EXPECT_EQ(3u, md_system_part_residue_range(&sys, 0).end - md_system_part_residue_range(&sys, 0).beg);
    EXPECT_EQ(2u, md_system_part_residue_range(&sys, 1).end - md_system_part_residue_range(&sys, 1).beg);
    EXPECT_EQ(2u, md_system_molecule_count(&sys));
    md_vm_arena_destroy(alloc);
}

static inline bool init_system(md_system_t* sys, md_system_state_t* sys_state, md_str_t path) {
    md_str_t ext = md_path_ext(path);
    if (md_str_empty(ext)) {
        return false;
    }

    if (md_str_eq_ignore_case(ext, MD_STR_LIT("pdb"))) {
        return md_pdb_system_init_from_file(sys, sys_state, path, NULL) && test_infer_as_loaded(sys, sys_state);
    } else
    if (md_str_eq_ignore_case(ext, MD_STR_LIT("gro"))) {
        return md_gro_system_init_from_file(sys, sys_state, path, NULL) && test_infer_as_loaded(sys, sys_state);
    } else
    if (md_str_eq_ignore_case(ext, MD_STR_LIT("cif"))) {
        return md_mmcif_system_init_from_file(sys, sys_state, path, NULL) && test_infer_as_loaded(sys, sys_state);
    }

    return false;
}

static md_part_type_flags_t inferred_kind(md_part_kind_t kind) {
    return md_part_type_flags_set_kind(MD_PART_TYPE_FLAG_INFERRED, kind);
}

// The chain of the first residue of a part, "" for one in no chain
static md_str_t part_chain_id(const md_system_t* sys, size_t part_idx) {
    const md_chain_idx_t c = md_system_chain_find_by_residue_idx(sys, md_system_part_residue_range(sys, part_idx).beg);
    return c >= 0 ? md_system_chain_id(sys, (size_t)c) : MD_STR_LIT("");
}

// Units only where the data gives them (the chains of a PDB file), parts and part types always inferred from them
UTEST(hierarchy, part_and_part_type) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    MD_ASSERT(alloc);

    {
        // No chain ids: no units, and the parts are found over the whole system
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);
        EXPECT_EQ(md_system_instance_count(&sys), 0);
        EXPECT_EQ(md_system_chain_count(&sys), 0);
        ASSERT_TRUE(test_infer_parts(&sys));
        ASSERT_EQ(md_system_part_type_count(&sys), 1);
        EXPECT_EQ(md_system_part_type_flags(&sys, 0), inferred_kind(MD_PART_KIND_PEPTIDE));
        EXPECT_TRUE(md_str_eq(md_system_part_type_description(&sys, 0), MD_STR_LIT("peptide")));
        ASSERT_EQ(md_system_part_count(&sys), 1);
        EXPECT_EQ(md_system_part_type_idx(&sys, 0), 0);
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);

        // The units are the chains of the file, named by them and by no made up label
        ASSERT_EQ(md_system_instance_count(&sys), 3);
        static const char* ids[] = {"A", "B", "C"};
        for (size_t i = 0; i < 3; ++i) {
            EXPECT_TRUE(md_str_empty(md_system_instance_id(&sys, i)));
            EXPECT_TRUE(md_str_eq(md_system_instance_auth_id(&sys, i), md_str_from_cstr(ids[i])));
        }
        ASSERT_TRUE(test_infer_parts(&sys));
        ASSERT_EQ(md_system_part_type_count(&sys), 1);
        EXPECT_EQ(md_system_part_type_flags(&sys, 0), inferred_kind(MD_PART_KIND_PEPTIDE));
        ASSERT_EQ(md_system_part_count(&sys), 3);
        for (size_t p = 0; p < 3; ++p) {
            EXPECT_EQ(md_system_part_type_idx(&sys, p), 0);
            EXPECT_TRUE(md_str_eq(part_chain_id(&sys, p), md_str_from_cstr(ids[p])));
        }
    }

    {
        // One chain, the arginine bound to it and its 106 waters, each water a part of its own. The arginine is numbered
        // on after the chain under its id, but is a free amino acid with a carboxyl group of its own: the chain ends
        // before it (CHAIN BREAKS in md_backbone.h)
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1LAF.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);
        ASSERT_TRUE(test_infer_parts(&sys));
        ASSERT_EQ(md_system_part_type_count(&sys), 3);
        EXPECT_EQ(md_system_part_type_flags(&sys, 0), inferred_kind(MD_PART_KIND_PEPTIDE));
        EXPECT_EQ(md_system_part_type_flags(&sys, 2), inferred_kind(MD_PART_KIND_WATER));
        EXPECT_TRUE(md_str_eq(md_system_part_type_description(&sys, 1), MD_STR_LIT("ARG")));

        ASSERT_EQ(md_system_part_count(&sys), 1 + 1 + 106);
        EXPECT_EQ(md_system_part_residue_range(&sys, 0).end - md_system_part_residue_range(&sys, 0).beg, 238);
        EXPECT_EQ(md_system_part_residue_range(&sys, 1).end - md_system_part_residue_range(&sys, 1).beg, 1);
        for (size_t p = 0; p < md_system_part_count(&sys); ++p) {
            EXPECT_TRUE(md_str_eq(part_chain_id(&sys, p), MD_STR_LIT("E")));
            if (p < 2) continue;
            EXPECT_EQ(md_system_part_type_idx(&sys, p), 2);
            EXPECT_EQ(md_system_part_residue_range(&sys, p).end - md_system_part_residue_range(&sys, p).beg, 1);
        }
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);

        // GTP and GDP are named like nucleotides atom by atom, but are not linked into a chain: ligands
        size_t num_nucleotides = 0;
        for (size_t i = 0; i < md_system_residue_count(&sys); ++i) {
            num_nucleotides += md_system_residue_kind(&sys, i) == MD_RESIDUE_KIND_NUCLEOTIDE;
        }
        EXPECT_EQ(num_nucleotides, 0);

        ASSERT_TRUE(test_infer_parts(&sys));
        static const struct { md_part_kind_t kind; const char* desc; } types[] = {
            {MD_PART_KIND_PEPTIDE, "peptide"}, {MD_PART_KIND_PEPTIDE, "peptide"}, {MD_PART_KIND_NON_POLYMER, "GTP"}, {MD_PART_KIND_NON_POLYMER, "MG"},
            {MD_PART_KIND_NON_POLYMER, "SO4"}, {MD_PART_KIND_WATER, "water"}, {MD_PART_KIND_NON_POLYMER, "GDP"}, {MD_PART_KIND_NON_POLYMER, "VLB"},
        };
        ASSERT_EQ(md_system_part_type_count(&sys), MD_COUNTOF(types));
        for (size_t i = 0; i < MD_COUNTOF(types); ++i) {
            EXPECT_EQ(md_system_part_type_flags(&sys, i), inferred_kind(types[i].kind));
            EXPECT_TRUE(md_str_eq(md_system_part_type_description(&sys, i), md_str_from_cstr(types[i].desc)));
        }

        // The strands, then the parts of chain A, of chain B and of chain C: the two SO4 and the five waters of chain A
        // are parts of their own
        static const struct { const char* chain; int type; } parts[] = {
            {"A", 0}, {"B", 1},
            {"A", 2}, {"A", 3}, {"A", 4}, {"A", 4}, {"A", 5}, {"A", 5}, {"A", 5}, {"A", 5}, {"A", 5},
            {"B", 6}, {"B", 4}, {"B", 5}, {"B", 5}, {"B", 5},
            {"C", 7},
        };
        ASSERT_EQ(md_system_part_count(&sys), MD_COUNTOF(parts));
        for (size_t i = 0; i < MD_COUNTOF(parts); ++i) {
            EXPECT_TRUE(md_str_eq(part_chain_id(&sys, i), md_str_from_cstr(parts[i].chain)));
            EXPECT_EQ(md_system_part_type_idx(&sys, i), parts[i].type);
        }
    }

    {
        // 64 lipids, the first 26 with chain ids of their own, and the waters: a unit for each id, and one for the rest,
        // which the file names nothing
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);
        EXPECT_EQ(md_system_instance_count(&sys), 27);
        EXPECT_EQ(md_system_chain_count(&sys), 26);
        EXPECT_TRUE(md_str_empty(md_system_instance_auth_id(&sys, 26)));

        ASSERT_TRUE(test_infer_parts(&sys));
        ASSERT_EQ(md_system_part_type_count(&sys), 2);
        EXPECT_EQ(md_system_part_type_flags(&sys, 0), inferred_kind(MD_PART_KIND_NON_POLYMER));
        EXPECT_EQ(md_system_part_type_flags(&sys, 1), inferred_kind(MD_PART_KIND_WATER));
        EXPECT_TRUE(md_str_eq(md_system_part_type_description(&sys, 0), MD_STR_LIT("DPP")));

        ASSERT_EQ(md_system_part_count(&sys), md_system_residue_count(&sys));
        for (size_t i = 0; i < md_system_part_count(&sys); ++i) {
            EXPECT_EQ(md_system_part_type_idx(&sys, i), i < 64 ? 0 : 1);
        }
    }

    md_temp_end(temp_scope);
}

// Part types the file does not define are inferred and flagged as such, file defined ones are not
UTEST(hierarchy, part_type_derived) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);

    {
        // Nucleosome: two DNA strands and eight histones, a .gro carries no units and no types
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &sys_state));
        md_system_infer(&sys, &sys_state, MD_SYSTEM_INFER_ALL, NULL);
        EXPECT_EQ(md_system_instance_count(&sys), 0);

        ASSERT_EQ(md_system_part_type_count(&sys), 6);
        // The DNA strands are DNA (by their residue names), which the nucleic backbone extraction relies upon
        EXPECT_EQ(md_system_part_type_flags(&sys, 0), inferred_kind(MD_PART_KIND_DNA));
        EXPECT_EQ(md_system_part_type_flags(&sys, 1), inferred_kind(MD_PART_KIND_DNA));
        EXPECT_TRUE(md_str_eq(md_system_part_type_description(&sys, 0), MD_STR_LIT("DNA")));
        for (size_t i = 2; i < 6; ++i) {
            EXPECT_EQ(md_system_part_type_flags(&sys, i), inferred_kind(MD_PART_KIND_PEPTIDE));
        }
        EXPECT_EQ(md_system_part_count(&sys), 10);
        EXPECT_EQ(sys.nucleic_backbone.range.count, 2);
        EXPECT_EQ(sys.protein_backbone.range.count, 8);

        // Deriving again starts over rather than appending
        ASSERT_TRUE(test_infer_parts(&sys));
        EXPECT_EQ(md_system_part_type_count(&sys), 6);
        EXPECT_EQ(md_system_part_count(&sys), 10);
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &sys_state));
        md_system_infer(&sys, &sys_state, MD_SYSTEM_INFER_ALL, NULL);

        // Each asym is one part, of the entity of the file
        ASSERT_GT(md_system_part_type_count(&sys), 0);
        for (size_t i = 0; i < md_system_part_type_count(&sys); ++i) {
            EXPECT_FALSE(md_system_part_type_flags(&sys, i) & MD_PART_TYPE_FLAG_INFERRED);
        }
        ASSERT_EQ(md_system_part_count(&sys), md_system_instance_count(&sys));
        for (size_t p = 0; p < md_system_part_count(&sys); ++p) {
            EXPECT_EQ(md_system_part_type_idx(&sys, p), md_system_instance_part_type_idx(&sys, p));
        }
    }

    md_temp_end(temp_scope);
}

// mmCIF: the part types are the file's entities, and the molecules of a non-polymer asym (the waters of a chain) are
// instances of their own which share the asym's id
UTEST(hierarchy, part_type_instance_mmcif) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &sys_state));

        static const struct { md_part_kind_t kind; const char* desc; } types[] = {
            {MD_PART_KIND_PEPTIDE, "PHOSPHONOACETALDEHYDE HYDROLASE"}, {MD_PART_KIND_NON_POLYMER, "MAGNESIUM ION"},
            {MD_PART_KIND_NON_POLYMER, "TUNGSTATE(VI)ION"}, {MD_PART_KIND_WATER, "water"},
        };
        ASSERT_EQ(md_system_part_type_count(&sys), MD_COUNTOF(types));
        for (size_t i = 0; i < MD_COUNTOF(types); ++i) {
            EXPECT_EQ(md_system_part_type_kind(&sys, i), types[i].kind);
            EXPECT_TRUE(md_str_eq(md_system_part_type_description(&sys, i), md_str_from_cstr(types[i].desc)));
        }

        static const struct { const char* id; const char* auth; int type; } instances[] = {
            {"A", "A", 0}, {"B", "B", 0}, {"E", "A", 1}, {"F", "B", 2}, {"G", "B", 1},
            {"K", "A", 3}, {"K", "A", 3}, {"K", "A", 3}, {"K", "A", 3},
            {"L", "B", 3}, {"L", "B", 3}, {"L", "B", 3}, {"L", "B", 3},
        };
        ASSERT_EQ(md_system_instance_count(&sys), MD_COUNTOF(instances));
        for (size_t i = 0; i < MD_COUNTOF(instances); ++i) {
            EXPECT_TRUE(md_str_eq(md_system_instance_id(&sys, i),      md_str_from_cstr(instances[i].id)));
            EXPECT_TRUE(md_str_eq(md_system_instance_auth_id(&sys, i), md_str_from_cstr(instances[i].auth)));
            EXPECT_EQ(md_system_instance_part_type_idx(&sys, i),    instances[i].type);
            if (i >= 2) EXPECT_EQ(md_system_instance_residue_count(&sys, i), 1);
        }
        // MG is an ion as a residue, whatever its entity is
        EXPECT_EQ(md_system_residue_kind(&sys, md_system_instance_residue_range(&sys, 2).beg), MD_RESIDUE_KIND_ION);
    }

    {
        // RNA strands which begin with a GTP (the 5' triphosphate): linked into the chain, so a nucleotide
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &sys_state));
        md_system_infer(&sys, &sys_state, MD_SYSTEM_INFER_ALL, NULL);
        ASSERT_GE(md_system_part_type_count(&sys), 5);
        EXPECT_EQ(md_system_part_type_kind(&sys, 2), MD_PART_KIND_RNA);
        EXPECT_EQ(md_system_part_type_kind(&sys, 3), MD_PART_KIND_RNA);
        size_t resolved = 0;
        for (size_t p = 0; p < md_system_part_count(&sys); ++p) {
            if (!md_part_kind_is_nucleic_acid(md_system_part_kind(&sys, p))) continue;
            const md_urange_t range = md_system_part_residue_range(&sys, p);
            EXPECT_TRUE(md_str_eq(md_system_residue_name(&sys, range.beg), MD_STR_LIT("GTP")) || md_str_eq(md_system_residue_name(&sys, range.beg), MD_STR_LIT("UTP")));
            for (uint32_t c = range.beg; c < range.end; ++c) {
                EXPECT_EQ(md_system_residue_kind(&sys, c), MD_RESIDUE_KIND_NUCLEOTIDE);
                resolved += (md_system_residue_flags(&sys, c) & MD_RESIDUE_FLAG_RESOLVED) != 0;
            }
        }
        EXPECT_EQ(resolved, 48);
        EXPECT_EQ(sys.nucleic_backbone.range.count, 2);
    }

    md_vm_arena_destroy(alloc);
}

// A hypermodified nucleotide of a tRNA: wybutosine (YYG) at 37 of the anticodon loop, between A36 and A38, its heavy
// atoms only (39 of them, more than any standard nucleotide has with its hydrogens). Its name is no standard one, but it
// is linked into the chain by its phosphodiester bonds: a nucleotide, on the backbone with its neighbours, its base the
// base. From the tRNA of rendering_all_atom.gro, a membrane complex of a deposited structure.
UTEST(hierarchy, hypermodified_nucleotide) {
    const char* text =
        "tRNA A36 YYG37 A38\n"
        "   83\n"
        "   36A        P 4250  15.017  13.816   1.786\n"
        "   36A      OP1 4251  14.910  13.913   1.824\n"
        "   36A      OP2 4252  15.149  13.867   1.738\n"
        "   36A      O5' 4253  15.043  13.715   1.907\n"
        "   36A      C5' 4254  14.936  13.633   1.952\n"
        "   36A      C4' 4255  14.986  13.521   2.041\n"
        "   36A      O4' 4256  15.062  13.421   1.969\n"
        "   36A      C3' 4257  15.077  13.551   2.159\n"
        "   36A      O3' 4258  15.009  13.616   2.268\n"
        "   36A      C2' 4259  15.122  13.410   2.199\n"
        "   36A      O2' 4260  15.033  13.344   2.287\n"
        "   36A      C1' 4261  15.120  13.334   2.065\n"
        "   36A       N9 4262  15.254  13.293   2.026\n"
        "   36A       C8 4263  15.349  13.351   1.945\n"
        "   36A       N7 4264  15.464  13.287   1.944\n"
        "   36A       C5 4265  15.442  13.178   2.028\n"
        "   36A       C6 4266  15.525  13.070   2.070\n"
        "   36A       N6 4267  15.652  13.053   2.030\n"
        "   36A       N1 4268  15.471  12.979   2.155\n"
        "   36A       C2 4269  15.344  12.994   2.192\n"
        "   36A       N3 4270  15.257  13.089   2.159\n"
        "   36A       C4 4271  15.313  13.179   2.076\n"
        "   37YYG      P 4272  15.093  13.697   2.379\n"
        "   37YYG    OP1 4273  14.999  13.793   2.445\n"
        "   37YYG    OP2 4274  15.218  13.749   2.316\n"
        "   37YYG    O5' 4275  15.137  13.584   2.483\n"
        "   37YYG    C5' 4276  15.038  13.502   2.550\n"
        "   37YYG    C4' 4277  15.103  13.379   2.608\n"
        "   37YYG    O4' 4278  15.148  13.291   2.502\n"
        "   37YYG    C3' 4279  15.229  13.406   2.688\n"
        "   37YYG    O3' 4280  15.207  13.449   2.822\n"
        "   37YYG    C2' 4281  15.304  13.273   2.681\n"
        "   37YYG    O2' 4282  15.260  13.178   2.775\n"
        "   37YYG    C1' 4283  15.270  13.227   2.539\n"
        "   37YYG     N9 4284  15.375  13.266   2.444\n"
        "   37YYG     C8 4285  15.370  13.366   2.351\n"
        "   37YYG     N7 4286  15.479  13.380   2.282\n"
        "   37YYG     C5 4287  15.563  13.281   2.332\n"
        "   37YYG     C6 4288  15.697  13.247   2.295\n"
        "   37YYG     O6 4289  15.767  13.299   2.207\n"
        "   37YYG     N1 4290  15.746  13.142   2.373\n"
        "   37YYG     C2 4291  15.678  13.080   2.473\n"
        "   37YYG     N2 4292  15.750  12.987   2.538\n"
        "   37YYG     N3 4293  15.553  13.107   2.506\n"
        "   37YYG     C3 4294  15.481  13.032   2.611\n"
        "   37YYG     C4 4295  15.501  13.209   2.432\n"
        "   37YYG    C10 4296  15.995  12.939   2.545\n"
        "   37YYG    C11 4297  15.877  13.004   2.490\n"
        "   37YYG    C12 4298  15.876  13.098   2.377\n"
        "   37YYG    C13 4299  15.995  13.119   2.280\n"
        "   37YYG    C14 4300  15.988  13.001   2.179\n"
        "   37YYG    C15 4301  16.023  13.027   2.031\n"
        "   37YYG    C16 4302  16.005  12.889   1.971\n"
        "   37YYG    O17 4303  15.889  12.843   1.960\n"
        "   37YYG    O18 4304  16.103  12.820   1.932\n"
        "   37YYG    C19 4305  16.182  12.868   1.821\n"
        "   37YYG    N20 4306  15.929  13.130   1.971\n"
        "   37YYG    C21 4307  15.923  13.166   1.837\n"
        "   37YYG    O22 4308  15.817  13.160   1.771\n"
        "   37YYG    O23 4309  16.028  13.206   1.781\n"
        "   37YYG    C24 4310  16.092  13.322   1.834\n"
        "   38A        P 4311  15.313  13.550   2.889\n"
        "   38A      OP1 4312  15.266  13.584   3.026\n"
        "   38A      OP2 4313  15.340  13.658   2.792\n"
        "   38A      O5' 4314  15.445  13.462   2.903\n"
        "   38A      C5' 4315  15.447  13.344   2.984\n"
        "   38A      C4' 4316  15.584  13.282   2.978\n"
        "   38A      O4' 4317  15.607  13.230   2.845\n"
        "   38A      C3' 4318  15.693  13.385   2.998\n"
        "   38A      O3' 4319  15.716  13.413   3.135\n"
        "   38A      C2' 4320  15.810  13.328   2.919\n"
        "   38A      O2' 4321  15.885  13.229   2.988\n"
        "   38A      C1' 4322  15.736  13.269   2.799\n"
        "   38A       N9 4323  15.716  13.370   2.695\n"
        "   38A       C8 4324  15.597  13.432   2.664\n"
        "   38A       N7 4325  15.606  13.519   2.566\n"
        "   38A       C5 4326  15.740  13.513   2.529\n"
        "   38A       C6 4327  15.814  13.580   2.429\n"
        "   38A       N6 4328  15.759  13.669   2.345\n"
        "   38A       N1 4329  15.946  13.553   2.419\n"
        "   38A       C2 4330  15.999  13.462   2.503\n"
        "   38A       N3 4331  15.940  13.391   2.600\n"
        "   38A       C4 4332  15.809  13.422   2.608\n"
        "122.00000 52.00000 36.00000\n"
        ;
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, md_str_from_cstr(text), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ((size_t)3, md_system_residue_count(&sys));
    EXPECT_TRUE(md_str_eq(md_system_residue_name(&sys, 1), MD_STR_LIT("YYG")));
    for (size_t c = 0; c < 3; ++c) {
        EXPECT_EQ(MD_RESIDUE_KIND_NUCLEOTIDE, md_system_residue_kind(&sys, c));
        EXPECT_TRUE(md_system_residue_flags(&sys, c) & MD_RESIDUE_FLAG_RESOLVED);
    }
    ASSERT_EQ((size_t)1, md_system_part_count(&sys));
    EXPECT_EQ(MD_PART_KIND_RNA, md_system_part_kind(&sys, 0));
    ASSERT_EQ((size_t)1, sys.nucleic_backbone.range.count);
    EXPECT_EQ((size_t)3, sys.nucleic_backbone.segment.count);

    // The backbone of the YYG is its P O5' C5' C4' C3' O3', its base all of it but the phosphate and the sugar: the
    // tricyclic base and its side chain
    const md_urange_t range = md_system_residue_atom_range(&sys, 1);
    size_t base = 0, backbone = 0;
    for (uint32_t i = range.beg; i < range.end; ++i) {
        base     += (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_NUCLEOBASE) != 0;
        backbone += (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_BACKBONE) != 0;
    }
    EXPECT_EQ((size_t)27, base);
    EXPECT_EQ((size_t)6, backbone);
    md_vm_arena_destroy(alloc);
}

// The chains are the chain ids the data gives (the author's), whatever carries them, and nothing where it gives none
static size_t chain_atoms(const md_system_t* sys, size_t c) {
    return md_system_chain_atom_count(sys, c);
}

UTEST(hierarchy, chains) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));

    {
        // One chain of a polymer, its arginine and its 106 waters: they carry its id, they are of it
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1LAF.pdb")));
        ASSERT_EQ(1, (int)md_system_chain_count(&sys));
        EXPECT_TRUE(md_str_eq(md_system_chain_id(&sys, 0), MD_STR_LIT("E")));
        EXPECT_EQ(md_system_atom_count(&sys), chain_atoms(&sys, 0));
        // The protein, the arginine and the waters are consecutive residues: one run
        EXPECT_EQ(1, (int)(md_system_chain_run_end(&sys, 0) - md_system_chain_run_beg(&sys, 0)));
        for (size_t c = 0; c < md_system_residue_count(&sys); ++c) EXPECT_EQ(0, md_system_chain_find_by_residue_idx(&sys, c));
        EXPECT_EQ(0, md_system_chain_find_by_atom_idx(&sys, md_system_atom_count(&sys) - 1));
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb")));
        ASSERT_EQ(3, (int)md_system_chain_count(&sys));
        static const char* ids[] = {"A", "B", "C"};
        size_t total = 0;
        for (size_t c = 0; c < 3; ++c) {
            EXPECT_TRUE(md_str_eq(md_system_chain_id(&sys, c), md_str_from_cstr(ids[c])));
            total += chain_atoms(&sys, c);
        }
        EXPECT_EQ(md_system_atom_count(&sys), total);
    }

    {
        // mmCIF: the chains are the author's ids, not the asyms. Asyms E, F, G, K and L hold ligands, ions and waters
        // the author filed under the chains A and B.
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif")));
        ASSERT_EQ(2, (int)md_system_chain_count(&sys));
        EXPECT_TRUE(md_str_eq(md_system_chain_id(&sys, 0), MD_STR_LIT("A")));
        EXPECT_TRUE(md_str_eq(md_system_chain_id(&sys, 1), MD_STR_LIT("B")));
        EXPECT_EQ(md_system_atom_count(&sys), chain_atoms(&sys, 0) + chain_atoms(&sys, 1));
        for (size_t i = 0; i < md_system_instance_count(&sys); ++i) {
            const md_urange_t residues = md_system_instance_residue_range(&sys, i);
            for (uint32_t k = residues.beg; k < residues.end; ++k) {
                const md_chain_idx_t c = md_system_chain_find_by_residue_idx(&sys, k);
                ASSERT_GE(c, 0);
                EXPECT_TRUE(md_str_eq(md_system_chain_id(&sys, (size_t)c), md_system_instance_auth_id(&sys, i)));
            }
        }
        // Chain A is its protein; after the protein of B, its Mg; after the WO4 and Mg of B, its waters: three runs
        EXPECT_EQ(3, (int)(md_system_chain_run_end(&sys, 0) - md_system_chain_run_beg(&sys, 0)));
    }

    {
        // No ids in the data: no chains, and the strands of the polymers are parts
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro")));
        ASSERT_TRUE(md_system_infer(&sys, &sys_state, MD_SYSTEM_INFER_ALL, NULL));
        EXPECT_EQ(0, (int)md_system_chain_count(&sys));
        size_t polymers = 0;
        for (size_t p = 0; p < md_system_part_count(&sys); ++p) {
            polymers += md_part_kind_is_polymer(md_system_part_kind(&sys, p));
            EXPECT_EQ(-1, md_system_chain_find_by_residue_idx(&sys, md_system_part_residue_range(&sys, p).beg));
        }
        EXPECT_EQ(253, (int)polymers);
    }

    {
        // A PDB file without chain ids (as simulations write them) is the same
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb")));
        EXPECT_EQ(0, (int)md_system_chain_count(&sys));
    }

    {
        // A topology names no chains either: the peptide is a polymer part, the waters and ions parts of their own
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip4p.tpr"), NULL));
        ASSERT_TRUE(md_system_infer(&sys, &sys_state, MD_SYSTEM_INFER_ALL & ~MD_SYSTEM_INFER_BOND_BIT, NULL));
        EXPECT_EQ(0, (int)md_system_chain_count(&sys));
        EXPECT_EQ(MD_PART_KIND_PEPTIDE, md_system_part_kind(&sys, 0));
        EXPECT_EQ(MD_PART_KIND_WATER, md_system_part_kind(&sys, 1));
    }

    {
        // martinize-dna writes the two strands of a duplex as one molecule type: one molecule, as GROMACS has it, of
        // two parts, a strand each, as the .gro of the same system gives them
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini2_dna.tpr"), NULL));
        ASSERT_TRUE(md_system_infer(&sys, &sys_state, MD_SYSTEM_INFER_ALL & ~MD_SYSTEM_INFER_BOND_BIT, NULL));
        EXPECT_EQ(0, (int)md_system_chain_count(&sys));
        const md_urange_t duplex = md_system_instance_residue_range(&sys, 0);
        EXPECT_EQ(24u, duplex.end - duplex.beg);
        EXPECT_EQ(0, md_system_part_find_by_residue_idx(&sys, duplex.beg));
        EXPECT_EQ(1, md_system_part_find_by_residue_idx(&sys, duplex.end - 1));
        EXPECT_EQ(12u, md_system_part_residue_range(&sys, 0).end - md_system_part_residue_range(&sys, 0).beg);
        const md_urange_t atoms = md_system_instance_atom_range(&sys, 0);
        const md_molecule_idx_t mol = md_system_molecule_find_by_atom_idx(&sys, atoms.beg);
        EXPECT_EQ(atoms.end - atoms.beg, md_system_molecule_atom_count(&sys, (size_t)mol));
        EXPECT_EQ(mol, md_system_molecule_find_by_atom_idx(&sys, atoms.end - 1));

        md_system_t gro = {.alloc = alloc};
        md_system_state_t gro_state = { .alloc = alloc };
        ASSERT_TRUE(md_gro_system_init_from_file(&gro, &gro_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini2_dna.gro"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&gro, &gro_state));
        ASSERT_TRUE(md_system_infer(&gro, &gro_state, MD_SYSTEM_INFER_ALL, NULL));
        for (size_t p = 0; p < 2; ++p) {
            const md_urange_t a = md_system_part_atom_range(&sys, p), b = md_system_part_atom_range(&gro, p);
            EXPECT_EQ(a.end - a.beg, b.end - b.beg);
            EXPECT_TRUE(md_part_kind_is_nucleic_acid(md_system_part_kind(&gro, p)));
        }
    }

    md_vm_arena_destroy(alloc);
}

// A molecule of a topology is one molecule whatever its bonds, as GROMACS has it: here two waters made one molecule of
// the topology by hand, which nothing else would join
UTEST(hierarchy, topology_molecules_are_whole) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip3p.tpr"), NULL));
    ASSERT_TRUE(md_system_infer(&sys, &sys_state, MD_SYSTEM_INFER_ALL & ~MD_SYSTEM_INFER_BOND_BIT, NULL));
    const size_t num_inst = md_system_instance_count(&sys);
    EXPECT_EQ(num_inst, md_system_molecule_count(&sys));
    const md_molecule_type_idx_t water = md_system_instance_molecule_type_idx(&sys, 1);
    ASSERT_GE(water, 0);
    ASSERT_EQ(water, md_system_instance_molecule_type_idx(&sys, 2));

    // Instances 1 and 2 as one molecule of the topology
    for (size_t i = 2; i < num_inst; ++i) {
        sys.instance.res_offset[i] = sys.instance.res_offset[i + 1];
        sys.instance.id[i] = sys.instance.id[i + 1 < num_inst ? i + 1 : i];
        sys.instance.auth_id[i] = sys.instance.auth_id[i + 1 < num_inst ? i + 1 : i];
        sys.instance.part_type_idx[i] = sys.instance.part_type_idx[i + 1 < num_inst ? i + 1 : i];
        sys.instance.mol_type_idx[i] = sys.instance.mol_type_idx[i + 1 < num_inst ? i + 1 : i];
    }
    sys.instance.count -= 1;
    ASSERT_TRUE(test_infer_molecules(&sys));
    EXPECT_EQ(num_inst - 1, md_system_molecule_count(&sys));
    check_molecule_invariants(utest_result, &sys);

    // Not a topology's: two molecules
    sys.instance.mol_type_idx[1] = -1;
    ASSERT_TRUE(test_infer_molecules(&sys));
    EXPECT_EQ(num_inst, md_system_molecule_count(&sys));
    md_vm_arena_destroy(alloc);
}

// Topologies: each molecule type is a molecule type of the system, each molecule an instance of it, and what the parts
// of the molecules are is inferred from the residues
UTEST(hierarchy, molecule_type_topology) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip4p.tpr"), NULL));

        static const struct { md_part_kind_t kind; const char* name; size_t count; } types[] = {
            {MD_PART_KIND_PEPTIDE, "Protein_chain_U", 1}, {MD_PART_KIND_WATER, "SOL", 506}, {MD_PART_KIND_NON_POLYMER, "NA", 2}, {MD_PART_KIND_NON_POLYMER, "CL", 2},
        };
        ASSERT_EQ(md_system_molecule_type_count(&sys), MD_COUNTOF(types));
        size_t inst = 0;
        for (size_t t = 0; t < MD_COUNTOF(types); ++t) {
            EXPECT_EQ(md_system_molecule_type_flags(&sys, t), MD_MOLECULE_TYPE_FLAG_NONE);
            EXPECT_TRUE(md_str_eq(md_system_molecule_type_name(&sys, t), md_str_from_cstr(types[t].name)));
            for (size_t k = 0; k < types[t].count; ++k, ++inst) {
                EXPECT_EQ(md_system_instance_molecule_type_idx(&sys, inst), (int)t);
                EXPECT_TRUE(md_str_empty(md_system_instance_id(&sys, inst)));
            }
        }
        ASSERT_EQ(md_system_instance_count(&sys), inst);
        EXPECT_EQ(md_system_instance_residue_count(&sys, 0), 2);

        // Inference leaves the topology's own alone: every molecule is of its type, and the parts are classified
        ASSERT_TRUE(md_system_infer(&sys, &sys_state, MD_SYSTEM_INFER_ALL & ~MD_SYSTEM_INFER_BOND_BIT, NULL));
        EXPECT_EQ(md_system_molecule_type_count(&sys), MD_COUNTOF(types));
        EXPECT_EQ(md_system_instance_count(&sys), inst);
        ASSERT_EQ(md_system_molecule_count(&sys), inst);
        inst = 0;
        for (size_t t = 0; t < MD_COUNTOF(types); ++t) {
            for (size_t k = 0; k < types[t].count; ++k, ++inst) {
                const md_urange_t atoms = md_system_instance_atom_range(&sys, inst);
                const md_molecule_idx_t mol = md_system_molecule_find_by_atom_idx(&sys, atoms.beg);
                EXPECT_EQ(md_system_molecule_type_idx(&sys, (size_t)mol), (int)t);
                const md_part_idx_t part = md_system_part_find_by_atom_idx(&sys, atoms.beg);
                EXPECT_EQ(md_system_part_kind(&sys, (size_t)part), types[t].kind);
            }
        }
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/martini3.tpr"), NULL));
        // The peptide, two lipids, four ions and 583 waters, as the structures are
        ASSERT_EQ(md_system_molecule_type_count(&sys), 5);
        EXPECT_EQ(md_system_instance_count(&sys), 1 + 2 + 4 + 583);
        ASSERT_TRUE(md_system_infer(&sys, &sys_state, MD_SYSTEM_INFER_ALL & ~MD_SYSTEM_INFER_BOND_BIT, NULL));
        EXPECT_EQ(md_system_part_kind(&sys, (size_t)md_system_part_find_by_residue_idx(&sys, md_system_instance_residue_range(&sys, 0).beg)), MD_PART_KIND_PEPTIDE);
        EXPECT_EQ(md_system_part_kind(&sys, (size_t)md_system_part_find_by_residue_idx(&sys, md_system_instance_residue_range(&sys, 1).beg)), MD_PART_KIND_NON_POLYMER);
        EXPECT_EQ(md_system_part_kind(&sys, md_system_part_count(&sys) - 1), MD_PART_KIND_WATER);
        EXPECT_EQ(md_system_molecule_type_count(&sys), 5);
    }

    {
        // LAMMPS knows its molecules but not what they are: a molecule type per sequence of particle types, named by
        // its formula, and water told by its atoms
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        const char* atom_format = md_lammps_atom_format_strings()[MD_LAMMPS_ATOM_FORMAT_FULL];
        ASSERT_TRUE(md_lammps_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"), &(md_load_desc_t){ .format = md_str_from_cstr(atom_format) }));
        ASSERT_EQ(md_system_molecule_type_count(&sys), 2);
        EXPECT_TRUE(md_str_eq(md_system_molecule_type_name(&sys, 0), MD_STR_LIT("C2H6")));
        EXPECT_TRUE(md_str_eq(md_system_molecule_type_name(&sys, 1), MD_STR_LIT("H2O")));
        ASSERT_EQ(md_system_instance_count(&sys), md_system_residue_count(&sys));
        EXPECT_EQ(md_system_instance_molecule_type_idx(&sys, 0), 0);
        EXPECT_EQ(md_system_instance_molecule_type_idx(&sys, md_system_instance_count(&sys) - 1), 1);
        ASSERT_TRUE(test_infer_parts(&sys));
        EXPECT_EQ(md_system_part_kind(&sys, 0), MD_PART_KIND_NON_POLYMER);
        EXPECT_EQ(md_system_part_kind(&sys, md_system_part_count(&sys) - 1), MD_PART_KIND_WATER);
    }

    md_vm_arena_destroy(alloc);
}

// The molecules of the same make share a type, inferred where no topology gives one: named by their parts, or by their
// formula where they have none. Inferring them again gives the same types, not more of them.
UTEST(hierarchy, molecule_type_inferred) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif"), NULL));
        ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        EXPECT_VALID_SYSTEM(&sys);
        // The two chains of the hydrolase, the Mg, the WO4 and the waters
        const size_t num_types = md_system_molecule_type_count(&sys);
        ASSERT_EQ(num_types, 4u);
        for (size_t t = 0; t < num_types; ++t) {
            EXPECT_TRUE(md_system_molecule_type_flags(&sys, t) & MD_MOLECULE_TYPE_FLAG_INFERRED);
        }
        const md_molecule_type_idx_t a = md_system_molecule_type_idx(&sys, (size_t)md_system_molecule_find_by_atom_idx(&sys, md_system_instance_atom_range(&sys, 0).beg));
        const md_molecule_type_idx_t b = md_system_molecule_type_idx(&sys, (size_t)md_system_molecule_find_by_atom_idx(&sys, md_system_instance_atom_range(&sys, 1).beg));
        EXPECT_EQ(a, b);
        EXPECT_TRUE(md_str_eq(md_system_molecule_type_name(&sys, (size_t)a), MD_STR_LIT("PHOSPHONOACETALDEHYDE HYDROLASE")));
        const md_molecule_type_idx_t w = md_system_molecule_type_idx(&sys, md_system_molecule_count(&sys) - 1);
        EXPECT_TRUE(md_str_eq(md_system_molecule_type_name(&sys, (size_t)w), MD_STR_LIT("water")));

        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        EXPECT_EQ(md_system_molecule_type_count(&sys), num_types);
        EXPECT_VALID_SYSTEM(&sys);
    }

    {
        // No residues, no parts: the type is named by the formula
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_xyz_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/C60-Ih.xyz"), NULL));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        ASSERT_EQ(md_system_molecule_count(&sys), 1u);
        ASSERT_EQ(md_system_molecule_type_count(&sys), 1u);
        EXPECT_TRUE(md_str_eq(md_system_molecule_type_name(&sys, 0), MD_STR_LIT("C60")));
    }

    md_vm_arena_destroy(alloc);
}

// Every change of the topology gives the system a version no consumer has seen
UTEST(hierarchy, topology_version) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    EXPECT_EQ(sys.topology_version, 0u);

    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &st));
    const uint64_t loaded = sys.topology_version;
    EXPECT_NE(loaded, 0u);

    ASSERT_TRUE(md_system_infer(&sys, &st, MD_SYSTEM_INFER_ALL, NULL));
    const uint64_t inferred = sys.topology_version;
    EXPECT_NE(inferred, loaded);

    // Inferring again derives the backbones from scratch rather than appending to them
    const size_t num_protein = sys.protein_backbone.segment.count, num_protein_ranges = sys.protein_backbone.range.count;
    const size_t num_nucleic = sys.nucleic_backbone.segment.count, num_nucleic_ranges = sys.nucleic_backbone.range.count;
    ASSERT_TRUE(md_system_infer(&sys, &st, MD_SYSTEM_INFER_ALL, NULL));
    EXPECT_NE(sys.topology_version, inferred);
    EXPECT_EQ(sys.protein_backbone.segment.count, num_protein);
    EXPECT_EQ(sys.protein_backbone.range.count, num_protein_ranges);
    EXPECT_EQ(sys.nucleic_backbone.segment.count, num_nucleic);
    EXPECT_EQ(sys.nucleic_backbone.range.count, num_nucleic_ranges);
    EXPECT_EQ(md_array_size(sys.protein_backbone.range.offset), num_protein_ranges + 1);

    // Bonds added and removed by hand
    const uint64_t before_insert = sys.topology_version;
    md_system_bond_insert(&sys, 0, 100, md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_USER));
    EXPECT_NE(sys.topology_version, before_insert);
    const uint64_t before_remove = sys.topology_version;
    md_system_bond_remove(&sys, (md_bond_idx_t)(sys.bond.count - 1));
    EXPECT_NE(sys.topology_version, before_remove);

    // Read only operations leave it alone
    const uint64_t read_only = sys.topology_version;
    md_bitmask_t mask = md_bitmask_create(0, (uint32_t)sys.atom.count, alloc);
    md_bitmask_set_bit(&mask, 0);
    md_bond_mask_grow(&mask, &sys.bond, sys.atom.count, 2, NULL);
    EXPECT_EQ(sys.topology_version, read_only);

    // A system loaded anew, into the same struct, has a version none of the above had
    const uint64_t last = sys.topology_version;
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &st));
    EXPECT_GT(sys.topology_version, last);

    md_system_free(&sys);
    EXPECT_EQ(sys.topology_version, 0u);
    md_vm_arena_destroy(alloc);
}

// The binary searches from an atom to its residue and from a residue to its instance
UTEST(hierarchy, find_by_index) {
    const uint32_t offset[] = {0, 3, 3, 7, 10};   // [0,3) [3,3) [3,7) [7,10)
    EXPECT_EQ(md_offset_range_find(offset, 4, 0), 0);
    EXPECT_EQ(md_offset_range_find(offset, 4, 2), 0);
    EXPECT_EQ(md_offset_range_find(offset, 4, 3), 2);   // The empty range holds nothing
    EXPECT_EQ(md_offset_range_find(offset, 4, 6), 2);
    EXPECT_EQ(md_offset_range_find(offset, 4, 7), 3);
    EXPECT_EQ(md_offset_range_find(offset, 4, 9), 3);
    EXPECT_EQ(md_offset_range_find(offset, 4, 10), -1);
    EXPECT_EQ(md_offset_range_find(offset, 0, 0), -1);
    EXPECT_EQ(md_offset_range_find(NULL, 4, 0), -1);

    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &sys_state));
    md_system_infer(&sys, &sys_state, MD_SYSTEM_INFER_ALL, NULL);
    for (size_t c = 0; c < sys.residue.count; ++c) {
        const md_urange_t range = md_system_residue_atom_range(&sys, c);
        EXPECT_EQ(md_system_residue_find_by_atom_idx(&sys, range.beg), (int)c);
        EXPECT_EQ(md_system_residue_find_by_atom_idx(&sys, range.end - 1), (int)c);
    }
    for (size_t i = 0; i < sys.instance.count; ++i) {
        const md_urange_t range = md_system_instance_atom_range(&sys, i);
        EXPECT_EQ(md_system_instance_find_by_atom_idx(&sys, range.beg), (int)i);
        EXPECT_EQ(md_system_instance_find_by_atom_idx(&sys, range.end - 1), (int)i);
    }
    EXPECT_EQ(md_system_residue_find_by_atom_idx(&sys, sys.atom.count), -1);
    md_vm_arena_destroy(alloc);
}

// The walk keeps a per atom depth in temp memory, where zero means 'not reached'. Temp memory is not cleared on
// allocation, so a walk which followed some other use of temp memory used to see stale depths and stop early.
UTEST(hierarchy, mask_grow_by_bonds_dirty_temp) {
    reference_systems_t refsys = reference_systems_load();
    const md_system_t* sys = &refsys.mol_ala;
    md_bitmask_t mask = md_bitmask_create(0, (uint32_t)sys->atom.count, refsys.alloc);
    md_bitmask_t ref  = md_bitmask_create(0, (uint32_t)sys->atom.count, refsys.alloc);

    // Stale depths of 1 make every neighbour look as if it had already been reached at a lower depth
    const size_t junk_size = MD_MEGABYTES(1);
    {
        md_temp_scope_t temp = md_temp_begin();
        memset(md_temp_alloc(temp, junk_size), 1, junk_size);
        md_temp_end(temp);
    }
    md_bitmask_set_bit(&mask, 0);
    md_bond_mask_grow(&mask, &sys->bond, sys->atom.count, 3, NULL);

    {
        md_temp_scope_t temp = md_temp_begin();
        memset(md_temp_alloc(temp, junk_size), 0, junk_size);
        md_temp_end(temp);
    }
    md_bitmask_set_bit(&ref, 0);
    md_bond_mask_grow(&ref, &sys->bond, sys->atom.count, 3, NULL);

    EXPECT_GT(md_bitmask_count(&ref), (size_t)1);
    EXPECT_EQ(md_bitmask_count(&ref), md_bitmask_count(&mask));
    reference_systems_free(&refsys);
}
