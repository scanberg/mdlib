#include "utest.h"

#include <md_pdb.h>
#include <md_gro.h>
#include <md_xyz.h>
#include <md_mmcif.h>
#include <md_system.h>
#include <md_infer.h>
#include <md_dssp.h>
#include <md_itp.h>
#include <md_element.h>
#include <core/md_writer.h>
#include <core/md_arena_allocator.h>

#include <string.h>
#include <stdio.h>

#include "system_invariants.h"
#include "system_util.h"

UTEST(backbone, backbone) {
    reference_systems_t refsys = reference_systems_load();
	EXPECT_EQ(0,        refsys.mol_pftaa.protein_backbone.range.count);
	EXPECT_EQ(0,        refsys.mol_nucleotides.protein_backbone.range.count);
    EXPECT_EQ(1,        refsys.mol_ala.protein_backbone.range.count);
    EXPECT_EQ(15,       refsys.mol_ala.protein_backbone.segment.count);
	EXPECT_EQ(253,      refsys.mol_centered.protein_backbone.range.count);
    EXPECT_EQ(10626,    refsys.mol_centered.protein_backbone.segment.count); // Should be equal to the total count of residues in chains
    reference_systems_free(&refsys);
}

// The caps of a peptide (ACE, NME) are of its strand: bonded into it by its backbone. The alanine dipeptide of tleap.
UTEST(backbone, caps_are_of_their_strand) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    static const char pdb[] =
        "ATOM      1 HH31 ACE     1       2.000   1.000  -0.000  1.00  0.00\n"
        "ATOM      2  CH3 ACE     1       2.000   2.090   0.000  1.00  0.00\n"
        "ATOM      3 HH32 ACE     1       1.486   2.454   0.890  1.00  0.00\n"
        "ATOM      4 HH33 ACE     1       1.486   2.454  -0.890  1.00  0.00\n"
        "ATOM      5  C   ACE     1       3.427   2.641  -0.000  1.00  0.00\n"
        "ATOM      6  O   ACE     1       4.391   1.877  -0.000  1.00  0.00\n"
        "ATOM      7  N   ALA     2       3.555   3.970  -0.000  1.00  0.00\n"
        "ATOM      8  H   ALA     2       2.733   4.556  -0.000  1.00  0.00\n"
        "ATOM      9  CA  ALA     2       4.853   4.614  -0.000  1.00  0.00\n"
        "ATOM     10  HA  ALA     2       5.408   4.316   0.890  1.00  0.00\n"
        "ATOM     11  CB  ALA     2       5.661   4.221  -1.232  1.00  0.00\n"
        "ATOM     12  HB1 ALA     2       5.123   4.521  -2.131  1.00  0.00\n"
        "ATOM     13  HB2 ALA     2       6.630   4.719  -1.206  1.00  0.00\n"
        "ATOM     14  HB3 ALA     2       5.809   3.141  -1.241  1.00  0.00\n"
        "ATOM     15  C   ALA     2       4.713   6.129   0.000  1.00  0.00\n"
        "ATOM     16  O   ALA     2       3.601   6.653   0.000  1.00  0.00\n"
        "ATOM     17  N   NME     3       5.846   6.835   0.000  1.00  0.00\n"
        "ATOM     18  H   NME     3       6.737   6.359  -0.000  1.00  0.00\n"
        "ATOM     19  CH3 NME     3       5.846   8.284   0.000  1.00  0.00\n"
        "ATOM     20 HH31 NME     3       4.819   8.648   0.000  1.00  0.00\n"
        "ATOM     21 HH32 NME     3       6.360   8.648   0.890  1.00  0.00\n"
        "ATOM     22 HH33 NME     3       6.360   8.648  -0.890  1.00  0.00\n"
        "END\n";
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(pdb), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ(3u, md_system_residue_count(&sys));
    EXPECT_EQ((int)MD_RESIDUE_KIND_AMINO_ACID, (int)md_system_residue_kind(&sys, 1));
    EXPECT_EQ(1u, md_system_part_count(&sys));
    EXPECT_EQ(1u, md_system_molecule_count(&sys));
    md_vm_arena_destroy(alloc);
}

// A coarse grained strand numbered over a gap, without bonds, in a chain the file names: one part, and so one molecule
// (the hierarchy links follow the part, not the numbering), as every part is within one molecule
UTEST(backbone, coarse_grained_strand_over_a_gap) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_str_writer_t sb = md_str_writer_create(alloc);
    static const int resids[] = {1, 2, 3, 4, 5, 10, 11, 12, 13, 14};
    for (int i = 0; i < 10; ++i) {
        const float x = 3.8f * i + (i >= 5 ? 20.0f : 0.0f);
        md_str_writer_fmt(&sb, "ATOM  %5d  BB  ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00            \n", i + 1, resids[i], x, 0.0, 0.0);
    }
    md_str_writer_push_cstr(&sb, "END\n");
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_writer_to_str(&sb), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    EXPECT_TRUE(md_system_is_coarse_grained(&sys));
    EXPECT_VALID_SYSTEM(&sys);
    ASSERT_EQ(1u, md_system_part_count(&sys));
    EXPECT_EQ(MD_PART_KIND_PEPTIDE, md_system_part_kind(&sys, 0));
    EXPECT_EQ(1u, md_system_molecule_count(&sys));
    md_vm_arena_destroy(alloc);
}

// A topology supplementing a coarse grained GRO splits what the numbering ran together into its molecules: the parts and
// the backbones follow, each backbone range of one part
UTEST(backbone, topology_supplement_reinfers_backbones) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_str_writer_t sb = md_str_writer_create(alloc);
    md_str_writer_push_cstr(&sb, "cg\n12\n");
    for (int i = 0; i < 12; ++i) {
        const float x = 0.70f * (i % 6) + (i >= 6 ? 3.0f : 0.0f);
        const float y = i >= 6 ? 3.0f : 0.0f;
        md_str_writer_fmt(&sb, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", i + 1, "ALA", "BB", i + 1, x + 1.0f, y + 1.0f, 1.0f);
    }
    md_str_writer_push_cstr(&sb, "   8.00000   8.00000   8.00000\n");
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, md_str_writer_to_str(&sb), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &state));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    EXPECT_EQ(1u, md_system_part_count(&sys));

    md_str_writer_t top = md_str_writer_create(alloc);
    md_str_writer_push_cstr(&top, "[ defaults ]\n1 1\n[ atomtypes ]\nP5 72.0 0.0 A 0.47 5.0\n[ moleculetype ]\nPROT 1\n[ atoms ]\n");
    for (int i = 0; i < 6; ++i) md_str_writer_fmt(&top, "%d P5 %d ALA BB %d 0.0 72.0\n", i + 1, i + 1, i + 1);
    md_str_writer_push_cstr(&top, "[ bonds ]\n");
    for (int i = 0; i < 5; ++i) md_str_writer_fmt(&top, "%d %d 1 0.35 1250\n", i + 1, i + 2);
    md_str_writer_push_cstr(&top, "[ system ]\nx\n[ molecules ]\nPROT 2\n");
    md_itp_data_t data = {0};
    ASSERT_TRUE(md_itp_data_parse_str(&data, md_str_writer_to_str(&top), MD_STR_LIT(""), alloc));
    ASSERT_TRUE(md_itp_system_supplement(&sys, &data));

    EXPECT_VALID_SYSTEM(&sys);
    ASSERT_EQ(2u, md_system_instance_count(&sys));
    EXPECT_EQ(2u, md_system_part_count(&sys));
    EXPECT_EQ(2u, md_system_molecule_count(&sys));
    EXPECT_EQ(1u, md_system_molecule_type_count(&sys));
    const md_protein_backbone_data_t* bb = &sys.protein_backbone;
    ASSERT_EQ(2u, bb->range.count);
    for (size_t r = 0; r < bb->range.count; ++r) {
        const md_part_idx_t part = bb->range.part_idx[r];
        ASSERT_GE(part, 0);
        const md_urange_t residues = md_system_part_residue_range(&sys, (size_t)part);
        for (uint32_t k = bb->range.offset[r]; k < bb->range.offset[r + 1]; ++k) {
            const md_residue_idx_t res = bb->segment.res_idx[k];
            EXPECT_TRUE(residues.beg <= (uint32_t)res && (uint32_t)res < residues.end);
        }
    }
    md_vm_arena_destroy(alloc);
}

// The backbone angles and secondary structure of a frame are the state's, in its attributes
UTEST(backbone, state_backbone) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&sys, &st));
    ASSERT_TRUE(md_system_infer(&sys, &st, MD_SYSTEM_INFER_ALL, NULL));
    const size_t n = sys.protein_backbone.segment.count;
    ASSERT_GT(n, 0);

    // Nothing until computed
    EXPECT_TRUE(md_system_state_backbone_angles(&st, &sys) == NULL);
    EXPECT_TRUE(md_system_state_secondary_structure(&st, &sys) == NULL);

    ASSERT_TRUE(md_system_state_backbone_compute(&st, &sys));
    const md_backbone_angles_t*     angle = md_system_state_backbone_angles(&st, &sys);
    const md_secondary_structure_t* ss    = md_system_state_secondary_structure(&st, &sys);
    ASSERT_TRUE(angle != NULL);
    ASSERT_TRUE(ss != NULL);

    // The same as computing them by hand
    md_backbone_angles_t*     ref_angle = md_alloc(alloc, n * sizeof(md_backbone_angles_t));
    md_secondary_structure_t* ref_ss    = md_alloc(alloc, n * sizeof(md_secondary_structure_t));
    md_backbone_angles_compute(ref_angle, n, st.xyz, &st.unitcell, &sys.protein_backbone);
    md_dssp_compute(ref_ss, n, st.xyz, &st.unitcell, &sys.protein_backbone);
    size_t structured = 0;
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(angle[i].phi, ref_angle[i].phi);
        EXPECT_EQ(angle[i].psi, ref_angle[i].psi);
        EXPECT_EQ(ss[i], ref_ss[i]);
        structured += ss[i] == MD_SECONDARY_STRUCTURE_HELIX_ALPHA || ss[i] == MD_SECONDARY_STRUCTURE_BETA_SHEET;
    }
    EXPECT_GT(structured, n / 4);

    // A producer writing every frame reuses the storage, and every write is a new version
    const md_attribute_t* attr = md_attributes_find(&st.attributes, MD_STR_LIT(MD_BACKBONE_SECONDARY_STRUCTURE_PATH));
    ASSERT_TRUE(attr != NULL);
    const uint64_t version = attr->version;
    md_secondary_structure_t* dst = md_system_state_secondary_structure_write(&st, &sys);
    EXPECT_TRUE(dst == ss);
    EXPECT_GT(md_attributes_version(&st.attributes, attr->id), version);

    // A view has no table to write into
    md_system_state_t view = { .num_atoms = st.num_atoms, .xyz = st.xyz, .unitcell = st.unitcell };
    EXPECT_TRUE(md_system_state_secondary_structure_write(&view, &sys) == NULL);
    EXPECT_FALSE(md_system_state_backbone_compute(&view, &sys));

    // What a state carries for another backbone does not fit this one
    md_system_t other = {.alloc = alloc};
    md_system_state_t other_st = {.alloc = alloc};
    ASSERT_TRUE(md_gro_system_init_from_file(&other, &other_st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"), NULL));
    ASSERT_TRUE(test_infer_as_loaded(&other, &other_st));
    ASSERT_TRUE(md_system_infer(&other, &other_st, MD_SYSTEM_INFER_ALL, NULL));
    EXPECT_TRUE(md_system_state_secondary_structure(&st, &other) == NULL);

    md_vm_arena_destroy(alloc);
}

UTEST(backbone, ramachandran_classify) {
    // Two chains, both starting with a proline, the second right after a residue which is not one:
    //   chain 0: PRO GLY PRO ALA PRO PRO SER
    //   chain 1: PRO ALA
    const char* names[9] = { "PRO", "GLY", "PRO", "ALA", "PRO", "PRO", "SER", "PRO", "ALA" };
    md_label_t labels[9];
    md_residue_idx_t res_idx[9];
    for (int i = 0; i < 9; ++i) {
        labels[i]   = md_label_make(md_str_from_cstr(names[i]));
        res_idx[i] = i;
    }
    uint32_t range_offset[3] = { 0, 7, 9 };

    md_system_t sys = {0};
    sys.residue.count = 9;
    sys.residue.name    = labels;
    sys.protein_backbone.range.count      = 2;
    sys.protein_backbone.range.offset     = range_offset;
    sys.protein_backbone.segment.count    = 9;
    sys.protein_backbone.segment.res_idx = res_idx;

    // MolProbity's types: a glycine or a proline before a proline keeps its own type, and the last residue of a chain
    // is not pre-proline because of the first residue of the next
    enum {
        G = MD_RAMACHANDRAN_TYPE_GENERAL,
        Y = MD_RAMACHANDRAN_TYPE_GLYCINE,
        P = MD_RAMACHANDRAN_TYPE_PROLINE,
        R = MD_RAMACHANDRAN_TYPE_PREPROL,
    };
    const md_ramachandran_type_t expected[9] = { P, Y, P, R, P, P, G, P, G };

    // A guard on each side of the buffer catches a write outside it
    md_ramachandran_type_t buf[1 + 9 + 1];
    md_memset(buf, 0xAB, sizeof(buf));
    EXPECT_TRUE(md_backbone_ramachandran_classify(buf + 1, 9, &sys));
    EXPECT_EQ(0xAB, buf[0]);
    EXPECT_EQ(0xAB, buf[10]);
    for (int i = 0; i < 9; ++i) EXPECT_EQ(expected[i], buf[1 + i]);

    // A capacity short of the segments: the last type written is still decided by the segment after it
    md_memset(buf, 0xAB, sizeof(buf));
    EXPECT_TRUE(md_backbone_ramachandran_classify(buf + 1, 4, &sys));
    EXPECT_EQ(0xAB, buf[0]);
    for (int i = 0; i < 4; ++i) EXPECT_EQ(expected[i], buf[1 + i]);
    EXPECT_EQ(0xAB, buf[5]);
}

// BACKBONES FROM THE BONDS (md_util.c). A file read with its residue and atom names, and the same file as an XYZ file
// would give it, the element and the position of each atom and nothing else, have the same backbones: the same atoms
// in each segment, consecutive along the same chains, and so the same secondary structure and ramachandran types.

// The XYZ text of a system, its atoms in their order, or by element (as VASP and ase.build.sort write them). Atom k of
// the text is atom xyz_to_sys[k] of the system; atoms without an element are left out. Returns the number of atoms.
static size_t xyz_text_of(md_str_t* out, int32_t* xyz_to_sys, const md_system_t* sys, const md_system_state_t* state, bool by_element, md_allocator_i* alloc) {
    size_t num = 0;
    for (int z = 1; z < 256; ++z) {
        for (size_t i = 0; i < sys->atom.count; ++i) {
            const md_atomic_number_t zi = md_atom_atomic_number(&sys->atom, i);
            if (zi != 0 && (!by_element || zi == z)) xyz_to_sys[num++] = (int32_t)i;
        }
        if (!by_element) break;
    }
    char* buf = md_alloc(alloc, 64 * (num + 1));
    size_t len = (size_t)snprintf(buf, 64, "%zu\nxyz\n", num);
    for (size_t k = 0; k < num; ++k) {
        const md_str_t sym = md_atomic_number_symbol(md_atom_atomic_number(&sys->atom, (size_t)xyz_to_sys[k]));
        const md_float3_t p = md_state_coord(state, (size_t)xyz_to_sys[k]);
        len += (size_t)snprintf(buf + len, 64, "%.*s %.4f %.4f %.4f\n", (int)sym.len, sym.ptr, p.x, p.y, p.z);
    }
    *out = (md_str_t){ buf, len };
    return num;
}

// The segment of each atom which is the CA (C4') of one, -1 for none
static int32_t* segment_of_atom(size_t num_atoms, md_allocator_i* alloc) {
    int32_t* seg = md_alloc(alloc, sizeof(int32_t) * num_atoms);
    md_memset(seg, -1, sizeof(int32_t) * num_atoms);
    return seg;
}

static void expect_backbones_from_bonds(int* utest_result, const char* path, bool by_element, bool compare_secondary_structure) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t ref = { .alloc = alloc };
    md_system_state_t ref_state = { .alloc = alloc };
    const md_str_t file = md_str_from_cstr(path);
    bool loaded = false;
    if (md_str_ends_with(file, MD_STR_LIT(".pdb"))) loaded = md_pdb_system_init_from_file(&ref, &ref_state, file, NULL);
    if (md_str_ends_with(file, MD_STR_LIT(".gro"))) loaded = md_gro_system_init_from_file(&ref, &ref_state, file, NULL);
    if (md_str_ends_with(file, MD_STR_LIT(".cif"))) loaded = md_mmcif_system_init_from_file(&ref, &ref_state, file, NULL);
    ASSERT_TRUE(loaded);
    ASSERT_TRUE(md_system_infer(&ref, &ref_state, MD_SYSTEM_INFER_ALL, NULL));

    int32_t* xyz_to_ref = md_alloc(alloc, sizeof(int32_t) * ref.atom.count);
    md_str_t text;
    const size_t num_xyz = xyz_text_of(&text, xyz_to_ref, &ref, &ref_state, by_element, alloc);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_xyz_system_init_from_str(&sys, &state, text, NULL));
    ASSERT_EQ(num_xyz, sys.atom.count);
    ASSERT_EQ(0, sys.residue.count);
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));

    // An atom of the reference as an atom of the XYZ
    int32_t* to_xyz = md_alloc(alloc, sizeof(int32_t) * ref.atom.count);
    md_memset(to_xyz, -1, sizeof(int32_t) * ref.atom.count);
    for (size_t k = 0; k < num_xyz; ++k) to_xyz[xyz_to_ref[k]] = (int32_t)k;
#define X(a) ((a) < 0 ? -1 : to_xyz[a])
#define SEG(a) (X(a) < 0 ? -1 : seg[X(a)])

    {
        const md_protein_backbone_data_t* rb = &ref.protein_backbone;
        const md_protein_backbone_data_t* bb = &sys.protein_backbone;
        EXPECT_EQ(rb->range.count, bb->range.count);
        EXPECT_EQ(rb->segment.count, bb->segment.count);
        for (size_t s = 0; s < bb->segment.count; ++s) EXPECT_EQ(-1, bb->segment.res_idx[s]);
        for (size_t r = 0; r < bb->range.count; ++r) EXPECT_EQ(-1, bb->range.part_idx[r]);

        int32_t* seg = segment_of_atom(sys.atom.count, alloc);
        for (size_t s = 0; s < bb->segment.count; ++s) seg[bb->segment.atoms[s].ca] = (int32_t)s;
        size_t same = 0;
        for (size_t r = 0; r < rb->range.count; ++r) {
            for (uint32_t j = rb->range.offset[r]; j < rb->range.offset[r + 1]; ++j) {
                const md_amino_acid_atoms_t a = rb->segment.atoms[j];
                const int32_t s = SEG(a.ca);
                if (s < 0) continue;
                const md_amino_acid_atoms_t b = bb->segment.atoms[s];
                // Either oxygen of a C terminal carboxylate
                const bool o = b.o == X(a.o) || (X(a.o) >= 0 && md_bond_find(&sys.bond, b.c, b.o) >= 0 && md_atom_atomic_number(&sys.atom, (size_t)b.o) == MD_Z_O);
                // Consecutive along the chain
                const bool chain = j == rb->range.offset[r] || (s > 0 && SEG(rb->segment.atoms[j - 1].ca) == s - 1);
                same += b.n == X(a.n) && b.c == X(a.c) && o && b.cb == X(a.cb) && chain;
            }
        }
        EXPECT_EQ(rb->segment.count, same);

        if (rb->segment.count && rb->segment.count == bb->segment.count) {
            md_ramachandran_type_t* rama_ref = md_alloc(alloc, sizeof(md_ramachandran_type_t) * rb->segment.count);
            md_ramachandran_type_t* rama = md_alloc(alloc, sizeof(md_ramachandran_type_t) * bb->segment.count);
            EXPECT_TRUE(md_backbone_ramachandran_classify(rama_ref, rb->segment.count, &ref));
            EXPECT_TRUE(md_backbone_ramachandran_classify(rama, bb->segment.count, &sys));
            md_secondary_structure_t* ss_ref = md_alloc(alloc, sizeof(md_secondary_structure_t) * rb->segment.count);
            md_secondary_structure_t* ss = md_alloc(alloc, sizeof(md_secondary_structure_t) * bb->segment.count);
            if (compare_secondary_structure) {
                EXPECT_TRUE(md_dssp_compute(ss_ref, rb->segment.count, ref_state.xyz, &ref_state.unitcell, rb));
                EXPECT_TRUE(md_dssp_compute(ss, bb->segment.count, state.xyz, &state.unitcell, bb));
            }
            size_t same_rama = 0, same_ss = 0;
            for (size_t j = 0; j < rb->segment.count; ++j) {
                const int32_t s = SEG(rb->segment.atoms[j].ca);
                if (s < 0) continue;
                // A residue whose side chain the crystal did not resolve, not even its CB, is a glycine by its bonds
                // (LYS 156 of 9gfu)
                const bool truncated = rb->segment.atoms[j].cb < 0 && rama_ref[j] != MD_RAMACHANDRAN_TYPE_GLYCINE;
                same_rama += rama[s] == rama_ref[j] || (truncated && rama[s] == MD_RAMACHANDRAN_TYPE_GLYCINE);
                same_ss += !compare_secondary_structure || md_memcmp(&ss[s], &ss_ref[j], sizeof(md_secondary_structure_t)) == 0;
            }
            EXPECT_EQ(rb->segment.count, same_rama);
            EXPECT_EQ(rb->segment.count, same_ss);
        }
    }

    {
        const md_nucleic_backbone_data_t* rb = &ref.nucleic_backbone;
        const md_nucleic_backbone_data_t* bb = &sys.nucleic_backbone;
        EXPECT_EQ(rb->range.count, bb->range.count);
        EXPECT_EQ(rb->segment.count, bb->segment.count);
        for (size_t s = 0; s < bb->segment.count; ++s) EXPECT_EQ(-1, bb->segment.res_idx[s]);

        int32_t* seg = segment_of_atom(sys.atom.count, alloc);
        for (size_t s = 0; s < bb->segment.count; ++s) seg[bb->segment.atoms[s].c4] = (int32_t)s;
        size_t same = 0;
        for (size_t r = 0; r < rb->range.count; ++r) {
            for (uint32_t j = rb->range.offset[r]; j < rb->range.offset[r + 1]; ++j) {
                const md_nucleic_acid_atoms_t a = rb->segment.atoms[j];
                const int32_t s = SEG(a.c4);
                if (s < 0) continue;
                const md_nucleic_acid_atoms_t b = bb->segment.atoms[s];
                const bool chain = j == rb->range.offset[r] || (s > 0 && SEG(rb->segment.atoms[j - 1].c4) == s - 1);
                // The P on the O5' of a 5' triphosphate, which the names do not tell (PA of the GTP which begins a
                // chain of 8g7u)
                const bool p = b.p == X(a.p) || (a.p < 0 && b.p >= 0 && md_bond_find(&sys.bond, b.o5, b.p) >= 0);
                same += p && b.o5 == X(a.o5) && b.c5 == X(a.c5) && b.c3 == X(a.c3) && b.o3 == X(a.o3) &&
                        b.c1 == X(a.c1) && b.c2 == X(a.c2) && b.o4 == X(a.o4) && chain;
            }
        }
        EXPECT_EQ(rb->segment.count, same);
    }
#undef SEG
#undef X

    md_vm_arena_destroy(alloc);
}

// Crystal structures, without hydrogens. In 1k4r the backbone O of HIS 86 of chain C is 1.4 Å from the ND1 of HIS 86 of
// chain A, close enough to be taken for a bond: the residue is still found, and its chain not broken there.
UTEST(backbone, backbones_from_bonds_protein) {
    const char* files[] = {
        MD_UNITTEST_DATA_DIR "/1k4r.pdb",
        MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb",
        MD_UNITTEST_DATA_DIR "/1fez.cif",
        MD_UNITTEST_DATA_DIR "/2or2.cif",
        MD_UNITTEST_DATA_DIR "/9gfu.cif",
    };
    for (size_t i = 0; i < MD_COUNTOF(files); ++i) {
        expect_backbones_from_bonds(utest_result, files[i], false, true);
    }
}

// Proteins with DNA (a nucleosome, with hydrogens) and with RNA. With hydrogens the HN of a residue is its hydrogen
// rather than the one the names leave out, which changes the secondary structure DSSP finds now and then.
UTEST(backbone, backbones_from_bonds_nucleic) {
    expect_backbones_from_bonds(utest_result, MD_UNITTEST_DATA_DIR "/nucl-dna.gro", false, false);
    expect_backbones_from_bonds(utest_result, MD_UNITTEST_DATA_DIR "/8g7u.cif", false, true);
    expect_backbones_from_bonds(utest_result, MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb", false, false);
}

// The order of the atoms in the file plays no part
UTEST(backbone, backbones_from_bonds_by_element) {
    expect_backbones_from_bonds(utest_result, MD_UNITTEST_DATA_DIR "/nucl-dna.gro", true, false);
    expect_backbones_from_bonds(utest_result, MD_UNITTEST_DATA_DIR "/1k4r.pdb", true, true);
}

// Residues whose atoms are not named as the standard ones: the atoms of a LAMMPS data file are named by their types, its
// residues "MOL". The backbones are found from the bonds all the same, each segment of the residue of its CA.
UTEST(backbone, backbones_from_bonds_with_residues) {
    const char* files[] = {
        MD_UNITTEST_DATA_DIR "/1k4r.pdb",
        MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb",
        MD_UNITTEST_DATA_DIR "/nucl-dna.gro",
    };
    for (size_t i = 0; i < MD_COUNTOF(files); ++i) {
        for (int rename_residues = 0; rename_residues < 2; ++rename_residues) {
            md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
            md_system_t ref = { .alloc = alloc }, sys = { .alloc = alloc };
            md_system_state_t ref_state = { .alloc = alloc }, state = { .alloc = alloc };
            const md_str_t file = md_str_from_cstr(files[i]);
            const bool gro = md_str_ends_with(file, MD_STR_LIT(".gro"));
            ASSERT_TRUE(gro ? md_gro_system_init_from_file(&ref, &ref_state, file, NULL) : md_pdb_system_init_from_file(&ref, &ref_state, file, NULL));
            ASSERT_TRUE(gro ? md_gro_system_init_from_file(&sys, &state, file, NULL) : md_pdb_system_init_from_file(&sys, &state, file, NULL));
            for (size_t t = 0; t < sys.atom.type.count; ++t) {
                char buf[16];
                snprintf(buf, sizeof(buf), "%u", (unsigned)(t + 1));
                sys.atom.type.name[t] = md_label_make(md_str_from_cstr(buf));
            }
            if (rename_residues) {
                for (size_t r = 0; r < sys.residue.count; ++r) sys.residue.name[r] = md_label_make(MD_STR_LIT("MOL"));
            }
            ASSERT_TRUE(md_system_infer(&ref, &ref_state, MD_SYSTEM_INFER_ALL, NULL));
            ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));

            const md_protein_backbone_data_t* ra = &ref.protein_backbone;
            const md_protein_backbone_data_t* pa = &sys.protein_backbone;
            ASSERT_TRUE(ra->segment.count > 0);
            ASSERT_EQ(ra->range.count, pa->range.count);
            ASSERT_EQ(ra->segment.count, pa->segment.count);
            size_t same = 0;
            for (size_t s = 0; s < pa->segment.count; ++s) {
                const md_amino_acid_atoms_t a = ra->segment.atoms[s], b = pa->segment.atoms[s];
                same += a.n == b.n && a.ca == b.ca && a.c == b.c && a.o == b.o && a.cb == b.cb && ra->segment.res_idx[s] == pa->segment.res_idx[s];
            }
            EXPECT_EQ(ra->segment.count, same);
            EXPECT_EQ(ref.nucleic_backbone.segment.count, sys.nucleic_backbone.segment.count);
            md_vm_arena_destroy(alloc);
        }
    }
}

// The atoms a segment stands for without residues, by the roles of its atoms: those of its residue as the file gives it,
// the atoms in the order of the file or sorted by element. The segments are matched by the atom of their CA.
static void expect_segment_atoms_by_roles(int* utest_result, const char* path, bool by_element, size_t num_different) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t ref = { .alloc = alloc };
    md_system_state_t ref_state = { .alloc = alloc };
    const md_str_t file = md_str_from_cstr(path);
    const bool gro = md_str_ends_with(file, MD_STR_LIT(".gro"));
    ASSERT_TRUE(gro ? md_gro_system_init_from_file(&ref, &ref_state, file, NULL) : md_pdb_system_init_from_file(&ref, &ref_state, file, NULL));
    ASSERT_TRUE(md_system_infer(&ref, &ref_state, MD_SYSTEM_INFER_ALL, NULL));

    int32_t* xyz_to_ref = md_alloc(alloc, sizeof(int32_t) * ref.atom.count);
    md_str_t text;
    xyz_text_of(&text, xyz_to_ref, &ref, &ref_state, by_element, alloc);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_xyz_system_init_from_str(&sys, &state, text, NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ(ref.protein_backbone.segment.count, sys.protein_backbone.segment.count);

    int32_t* ref_seg = segment_of_atom(ref.atom.count, alloc);
    for (size_t s = 0; s < ref.protein_backbone.segment.count; ++s) ref_seg[ref.protein_backbone.segment.atoms[s].ca] = (int32_t)s;
    uint8_t* in_ref = md_alloc(alloc, ref.atom.count);
    size_t same = 0;
    for (size_t s = 0; s < sys.protein_backbone.segment.count; ++s) {
        md_atom_idx_t atoms[512], ref_atoms[512];
        const size_t n = md_system_protein_segment_atoms(atoms, 512, &sys, s);
        const int32_t rs = ref_seg[xyz_to_ref[sys.protein_backbone.segment.atoms[s].ca]];
        if (rs < 0) continue;
        const size_t m = md_system_protein_segment_atoms(ref_atoms, 512, &ref, (size_t)rs);
        md_memset(in_ref, 0, ref.atom.count);
        for (size_t k = 0; k < m; ++k) in_ref[ref_atoms[k]] = 1;
        bool equal = n == m;
        for (size_t k = 0; k < n && equal; ++k) equal = in_ref[xyz_to_ref[atoms[k]]];
        same += equal;
    }
    EXPECT_EQ(sys.protein_backbone.segment.count - num_different, same);
    md_vm_arena_destroy(alloc);
}

UTEST(backbone, protein_segment_atoms_by_roles) {
    // Without hydrogens and with 18 disulfides (1k4r), and with hydrogens and termini (the nucleosome, whose histones are
    // free at both ends). In 1k4r
    // the backbone O of HIS 86 of chain C is 1.4 Å from the ND1 of HIS 86 of chain A, taken for a bond: without
    // residues the two amino acids reach each other over it (the O is of the one, the side chain of the other)
    expect_segment_atoms_by_roles(utest_result, MD_UNITTEST_DATA_DIR "/1k4r.pdb", false, 2);
    expect_segment_atoms_by_roles(utest_result, MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb", true, 0);
    expect_segment_atoms_by_roles(utest_result, MD_UNITTEST_DATA_DIR "/nucl-dna.gro", false, 0);
}

// A reference for the backbones found from the bonds, independent of them: the names the files give the atoms. Every
// atom of every segment is the one its name says it is, and there are as many segments as the atom names found before
// the backbones came from the bonds (when the names were what found them).
static bool name_is(md_str_t name, const char* const* names) {
    for (; *names; ++names) {
        if (md_str_eq_cstr(name, *names)) return true;
    }
    return false;
}

// The atom of the name, or primed as of a sugar (C4', C4*), -1 counting as named (a field the segment does not have)
static bool atom_named(const md_system_t* sys, md_atom_idx_t atom, const char* base, bool primed) {
    if (atom < 0) return true;
    const md_str_t name = md_atom_name(&sys->atom, (size_t)atom);
    const size_t n = strlen(base);
    if (!md_str_begins_with(name, (md_str_t){ base, n })) return false;
    if (!primed) return name.len == n;
    return name.len == n + 1 && (name.ptr[n] == '\'' || name.ptr[n] == '*');
}

UTEST(backbone, backbones_named_as_their_atoms) {
    const struct { const char* file; size_t protein; size_t nucleic; } refs[] = {
        { MD_UNITTEST_DATA_DIR "/1k4r.pdb",              1185,   0 },
        { MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb",        865,   0 },
        { MD_UNITTEST_DATA_DIR "/1LAF.pdb",               238,   0 },
        { MD_UNITTEST_DATA_DIR "/1a64.pdb",               188,   0 },
        { MD_UNITTEST_DATA_DIR "/1fez.cif",               512,   0 },
        { MD_UNITTEST_DATA_DIR "/2or2.cif",               592,   0 },
        { MD_UNITTEST_DATA_DIR "/9gfu.cif",               520,   0 },
        { MD_UNITTEST_DATA_DIR "/8g7u.cif",              1634,  48 },
        { MD_UNITTEST_DATA_DIR "/nucl-dna.gro",           974, 294 },
        { MD_UNITTEST_DATA_DIR "/centered.gro",         10626,   0 },
        { MD_UNITTEST_DATA_DIR "/inside-md-pullout.gro",  119,   0 },
    };
    static const char* const O[]  = { "O", "OXT", "O1", "O2", "OT1", "OT2", "OC1", "OC2", NULL };
    static const char* const HN[] = { "H", "HN", NULL };
    for (size_t i = 0; i < MD_COUNTOF(refs); ++i) {
        md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        const md_str_t file = md_str_from_cstr(refs[i].file);
        bool loaded = false;
        if (md_str_ends_with(file, MD_STR_LIT(".pdb"))) loaded = md_pdb_system_init_from_file(&sys, &state, file, NULL);
        if (md_str_ends_with(file, MD_STR_LIT(".gro"))) loaded = md_gro_system_init_from_file(&sys, &state, file, NULL);
        if (md_str_ends_with(file, MD_STR_LIT(".cif"))) loaded = md_mmcif_system_init_from_file(&sys, &state, file, NULL);
        ASSERT_TRUE(loaded);
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));

        const md_protein_backbone_data_t* pb = &sys.protein_backbone;
        EXPECT_EQ(refs[i].protein, pb->segment.count);
        size_t named = 0;
        for (size_t s = 0; s < pb->segment.count; ++s) {
            const md_amino_acid_atoms_t a = pb->segment.atoms[s];
            named += atom_named(&sys, a.n, "N", false) && atom_named(&sys, a.ca, "CA", false) && atom_named(&sys, a.c, "C", false) &&
                     name_is(md_atom_name(&sys.atom, (size_t)a.o), O) && atom_named(&sys, a.cb, "CB", false) &&
                     (a.hn < 0 || name_is(md_atom_name(&sys.atom, (size_t)a.hn), HN));
        }
        EXPECT_EQ(pb->segment.count, named);

        const md_nucleic_backbone_data_t* nb = &sys.nucleic_backbone;
        EXPECT_EQ(refs[i].nucleic, nb->segment.count);
        named = 0;
        for (size_t s = 0; s < nb->segment.count; ++s) {
            const md_nucleic_acid_atoms_t a = nb->segment.atoms[s];
            // The P on the O5' of a 5' triphosphate is its PA
            const bool p = a.p < 0 || atom_named(&sys, a.p, "P", false) || atom_named(&sys, a.p, "PA", false);
            named += p && atom_named(&sys, a.o5, "O5", true) && atom_named(&sys, a.c5, "C5", true) && atom_named(&sys, a.c4, "C4", true) &&
                     atom_named(&sys, a.c3, "C3", true) && atom_named(&sys, a.o3, "O3", true) && atom_named(&sys, a.c1, "C1", true) &&
                     atom_named(&sys, a.c2, "C2", true) && atom_named(&sys, a.o4, "O4", true);
        }
        EXPECT_EQ(nb->segment.count, named);
        md_vm_arena_destroy(alloc);
    }
}

// Nothing which is no polymer: fullerenes, and a tryptophan alone (one residue is no chain)
UTEST(backbone, backbones_from_bonds_none) {
    const char* files[] = {
        MD_UNITTEST_DATA_DIR "/c720.xyz",
        MD_UNITTEST_DATA_DIR "/tryptophan.xyz",
    };
    for (size_t i = 0; i < MD_COUNTOF(files); ++i) {
        md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_xyz_system_init_from_file(&sys, &state, md_str_from_cstr(files[i]), NULL));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        EXPECT_TRUE(sys.bond.count > 0);
        EXPECT_EQ(0, sys.protein_backbone.range.count);
        EXPECT_EQ(0, sys.protein_backbone.segment.count);
        EXPECT_EQ(0, sys.nucleic_backbone.range.count);
        md_vm_arena_destroy(alloc);
    }
}
