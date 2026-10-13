#include "utest.h"
#include <string.h>

#include <md_pdb.h>
#include <md_system.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

// Secondary structure, atomistic (DSSP) and of coarse grained traces (CG-DSSP), see md_util.c

// The letters of md_secondary_structure_t: unknown, coil, turn, bend, 3-10 helix, alpha helix, pi helix, sheet, bridge
static const char ss_letters[] = "?CTSGHIEB";

// Helix, sheet (or bridge) and coil: what a cartoon tells apart
static char ss_three_state(md_secondary_structure_t ss) {
    switch (ss) {
    case MD_SECONDARY_STRUCTURE_HELIX_310:
    case MD_SECONDARY_STRUCTURE_HELIX_ALPHA:
    case MD_SECONDARY_STRUCTURE_HELIX_PI:    return 'H';
    case MD_SECONDARY_STRUCTURE_BETA_SHEET:
    case MD_SECONDARY_STRUCTURE_BETA_BRIDGE: return 'E';
    default:                                 return 'C';
    }
}

// Ubiquitin (1UBQ) in Martini 3, as martinize2 (vermouth 0.15) maps it: a BB bead per residue at the center of its
// N, CA, C and O, and the side chain beads. From the integration tests of vermouth (tier-1/1UBQ).
static const char martini3_ubiquitin_pdb[] =
    "ATOM      1 BB   MET A   1      27.168  25.742   3.357  1.00  0.00\n"
    "ATOM      2 SC1  MET A   1      24.481  24.285   5.673  1.00  0.00\n"
    "ATOM      3 BB   GLN A   2      25.962  28.749   4.446  1.00  0.00\n"
    "ATOM      4 SC1  GLN A   2      26.776  32.120   2.397  1.00  0.00\n"
    "ATOM      5 BB   ILE A   3      27.033  30.748   7.177  1.00  0.00\n"
    "ATOM      6 SC1  ILE A   3      26.903  28.414   9.087  1.00  0.00\n"
    "ATOM      7 BB   PHE A   4      26.587  32.877  10.061  1.00  0.00\n"
    "ATOM      8 SC1  PHE A   4      25.404  34.579   7.877  1.00  0.00\n"
    "ATOM      9 SC2  PHE A   4      24.102  34.183   5.906  1.00  0.00\n"
    "ATOM     10 SC3  PHE A   4      25.636  35.236   5.594  1.00  0.00\n"
    "ATOM     11 BB   VAL A   5      28.795  34.925  12.172  1.00  0.00\n"
    "ATOM     12 SC1  VAL A   5      30.044  32.850  13.157  1.00  0.00\n"
    "ATOM     13 BB   LYS A   6      28.034  36.809  14.994  1.00  0.00\n"
    "ATOM     14 SC1  LYS A   6      25.484  38.706  14.521  1.00  0.00\n"
    "ATOM     15 SC2  LYS A   6      23.063  40.721  14.834  1.00  0.00\n"
    "ATOM     16 BB   THR A   7      29.486  39.423  16.577  1.00  0.00\n"
    "ATOM     17 SC1  THR A   7      31.896  39.157  15.771  1.00  0.00\n"
    "ATOM     18 BB   LEU A   8      29.841  41.902  19.151  1.00  0.00\n"
    "ATOM     19 SC1  LEU A   8      29.027  40.021  21.711  1.00  0.00\n"
    "ATOM     20 BB   THR A   9      31.098  44.070  16.796  1.00  0.00\n"
    "ATOM     21 SC1  THR A   9      33.260  43.398  17.341  1.00  0.00\n"
    "ATOM     22 BB   GLY A  10      29.381  43.791  14.001  1.00  0.00\n"
    "ATOM     23 BB   LYS A  11      30.588  41.272  12.826  1.00  0.00\n"
    "ATOM     24 SC1  LYS A  11      33.571  41.201  11.734  1.00  0.00\n"
    "ATOM     25 SC2  LYS A  11      35.337  40.430   9.627  1.00  0.00\n"
    "ATOM     26 BB   THR A  12      30.302  39.143   9.931  1.00  0.00\n"
    "ATOM     27 SC1  THR A  12      27.622  38.894  10.317  1.00  0.00\n"
    "ATOM     28 BB   ILE A  13      30.808  35.743   9.215  1.00  0.00\n"
    "ATOM     29 SC1  ILE A  13      33.486  35.168  10.482  1.00  0.00\n"
    "ATOM     30 BB   THR A  14      31.524  33.585   6.560  1.00  0.00\n"
    "ATOM     31 SC1  THR A  14      30.429  34.850   4.854  1.00  0.00\n"
    "ATOM     32 BB   LEU A  15      30.772  30.109   6.105  1.00  0.00\n"
    "ATOM     33 SC1  LEU A  15      32.488  29.742   9.143  1.00  0.00\n"
    "ATOM     34 BB   GLU A  16      31.830  26.882   5.026  1.00  0.00\n"
    "ATOM     35 SC1  GLU A  16      31.729  28.107   0.924  1.00  0.00\n"
    "ATOM     36 BB   VAL A  17      29.533  24.119   5.426  1.00  0.00\n"
    "ATOM     37 SC1  VAL A  17      29.747  25.074   7.920  1.00  0.00\n"
    "ATOM     38 BB   GLU A  18      28.555  21.077   7.047  1.00  0.00\n"
    "ATOM     39 SC1  GLU A  18      28.860  20.012   3.364  1.00  0.00\n"
    "ATOM     40 BB   PRO A  19      26.326  19.303   8.931  1.00  0.00\n"
    "ATOM     41 SC1  PRO A  19      24.827  19.858   6.910  1.00  0.00\n"
    "ATOM     42 BB   SER A  20      28.698  17.128   9.627  1.00  0.00\n"
    "ATOM     43 SC1  SER A  20      28.765  16.177   7.488  1.00  0.00\n"
    "ATOM     44 BB   ASP A  21      29.994  19.105  11.258  1.00  0.00\n"
    "ATOM     45 SC1  ASP A  21      31.897  20.425   8.722  1.00  0.00\n"
    "ATOM     46 BB   THR A  22      31.696  20.027  13.870  1.00  0.00\n"
    "ATOM     47 SC1  THR A  22      32.964  17.964  14.694  1.00  0.00\n"
    "ATOM     48 BB   ILE A  23      32.150  22.357  16.298  1.00  0.00\n"
    "ATOM     49 SC1  ILE A  23      29.691  22.748  17.483  1.00  0.00\n"
    "ATOM     50 BB   GLU A  24      35.221  22.167  16.349  1.00  0.00\n"
    "ATOM     51 SC1  GLU A  24      37.212  19.122  18.014  1.00  0.00\n"
    "ATOM     52 BB   ASN A  25      35.536  22.822  13.283  1.00  0.00\n"
    "ATOM     53 SC1  ASN A  25      35.563  21.422  10.874  1.00  0.00\n"
    "ATOM     54 BB   VAL A  26      34.258  25.560  13.239  1.00  0.00\n"
    "ATOM     55 SC1  VAL A  26      31.632  25.476  13.061  1.00  0.00\n"
    "ATOM     56 BB   LYS A  27      36.267  26.867  15.221  1.00  0.00\n"
    "ATOM     57 SC1  LYS A  27      34.856  26.242  18.143  1.00  0.00\n"
    "ATOM     58 SC2  LYS A  27      33.340  26.095  20.889  1.00  0.00\n"
    "ATOM     59 BB   ALA A  28      38.713  26.410  13.302  1.00  0.00\n"
    "ATOM     60 SC1  ALA A  28      39.285  24.336  13.566  1.00  0.00\n"
    "ATOM     61 BB   LYS A  29      37.669  28.208  10.979  1.00  0.00\n"
    "ATOM     62 SC1  LYS A  29      35.701  25.986   9.407  1.00  0.00\n"
    "ATOM     63 SC2  LYS A  29      34.858  24.630   7.053  1.00  0.00\n"
    "ATOM     64 BB   ILE A  30      37.604  30.767  12.684  1.00  0.00\n"
    "ATOM     65 SC1  ILE A  30      34.773  31.012  14.047  1.00  0.00\n"
    "ATOM     66 BB   GLN A  31      40.554  30.857  13.305  1.00  0.00\n"
    "ATOM     67 SC1  GLN A  31      42.686  28.923  15.862  1.00  0.00\n"
    "ATOM     68 BB   ASP A  32      41.622  30.874  10.453  1.00  0.00\n"
    "ATOM     69 SC1  ASP A  32      42.604  28.593   8.949  1.00  0.00\n"
    "ATOM     70 BB   LYS A  33      40.000  33.645   9.631  1.00  0.00\n"
    "ATOM     71 SC1  LYS A  33      37.395  33.079   8.328  1.00  0.00\n"
    "ATOM     72 SC2  LYS A  33      36.024  33.483   5.945  1.00  0.00\n"
    "ATOM     73 BB   GLU A  34      40.230  35.519  12.664  1.00  0.00\n"
    "ATOM     74 SC1  GLU A  34      37.378  37.071  11.301  1.00  0.00\n"
    "ATOM     75 BB   GLY A  35      42.134  34.179  15.169  1.00  0.00\n"
    "ATOM     76 BB   ILE A  36      40.297  32.743  17.186  1.00  0.00\n"
    "ATOM     77 SC1  ILE A  36      38.067  34.832  17.587  1.00  0.00\n"
    "ATOM     78 BB   PRO A  37      40.452  30.909  19.795  1.00  0.00\n"
    "ATOM     79 SC1  PRO A  37      42.307  32.248  20.426  1.00  0.00\n"
    "ATOM     80 BB   PRO A  38      38.576  28.222  20.676  1.00  0.00\n"
    "ATOM     81 SC1  PRO A  38      40.069  27.147  18.801  1.00  0.00\n"
    "ATOM     82 BB   ASP A  39      38.580  28.746  23.873  1.00  0.00\n"
    "ATOM     83 SC1  ASP A  39      40.554  28.602  25.687  1.00  0.00\n"
    "ATOM     84 BB   GLN A  40      36.961  31.623  23.409  1.00  0.00\n"
    "ATOM     85 SC1  GLN A  40      40.364  33.805  23.361  1.00  0.00\n"
    "ATOM     86 BB   GLN A  41      34.467  30.166  22.116  1.00  0.00\n"
    "ATOM     87 SC1  GLN A  41      35.843  30.189  18.265  1.00  0.00\n"
    "ATOM     88 BB   ARG A  42      30.949  30.935  21.949  1.00  0.00\n"
    "ATOM     89 SC1  ARG A  42      29.843  31.014  24.911  1.00  0.00\n"
    "ATOM     90 SC2  ARG A  42      27.111  32.424  26.552  1.00  0.00\n"
    "ATOM     91 BB   LEU A  43      28.193  29.117  20.658  1.00  0.00\n"
    "ATOM     92 SC1  LEU A  43      30.128  28.857  17.747  1.00  0.00\n"
    "ATOM     93 BB   ILE A  44      24.942  30.412  19.404  1.00  0.00\n"
    "ATOM     94 SC1  ILE A  44      24.811  31.210  22.222  1.00  0.00\n"
    "ATOM     95 BB   PHE A  45      21.716  28.605  19.082  1.00  0.00\n"
    "ATOM     96 SC1  PHE A  45      21.723  27.916  16.308  1.00  0.00\n"
    "ATOM     97 SC2  PHE A  45      20.484  28.426  14.315  1.00  0.00\n"
    "ATOM     98 SC3  PHE A  45      19.953  26.813  15.122  1.00  0.00\n"
    "ATOM     99 BB   ALA A  46      18.627  28.819  19.893  1.00  0.00\n"
    "ATOM    100 SC1  ALA A  46      17.864  27.977  18.346  1.00  0.00\n"
    "ATOM    101 BB   GLY A  47      19.668  29.225  23.002  1.00  0.00\n"
    "ATOM    102 BB   LYS A  48      22.231  27.401  22.566  1.00  0.00\n"
    "ATOM    103 SC1  LYS A  48      20.213  24.827  22.126  1.00  0.00\n"
    "ATOM    104 SC2  LYS A  48      18.064  24.455  20.590  1.00  0.00\n"
    "ATOM    105 BB   GLN A  49      25.049  25.843  23.458  1.00  0.00\n"
    "ATOM    106 SC1  GLN A  49      27.856  27.471  25.815  1.00  0.00\n"
    "ATOM    107 BB   LEU A  50      27.608  24.637  21.687  1.00  0.00\n"
    "ATOM    108 SC1  LEU A  50      25.939  25.527  18.770  1.00  0.00\n"
    "ATOM    109 BB   GLU A  51      29.078  21.449  21.268  1.00  0.00\n"
    "ATOM    110 SC1  GLU A  51      28.624  19.354  24.704  1.00  0.00\n"
    "ATOM    111 BB   ASP A  52      32.050  20.088  20.010  1.00  0.00\n"
    "ATOM    112 SC1  ASP A  52      34.187  22.155  21.368  1.00  0.00\n"
    "ATOM    113 BB   GLY A  53      30.918  16.866  19.371  1.00  0.00\n"
    "ATOM    114 BB   ARG A  54      28.848  18.003  17.705  1.00  0.00\n"
    "ATOM    115 SC1  ARG A  54      26.519  17.453  19.841  1.00  0.00\n"
    "ATOM    116 SC2  ARG A  54      27.646  14.411  20.672  1.00  0.00\n"
    "ATOM    117 BB   THR A  55      26.733  18.658  15.109  1.00  0.00\n"
    "ATOM    118 SC1  THR A  55      27.031  16.520  13.692  1.00  0.00\n"
    "ATOM    119 BB   LEU A  56      24.860  20.641  13.180  1.00  0.00\n"
    "ATOM    120 SC1  LEU A  56      27.152  22.891  11.817  1.00  0.00\n"
    "ATOM    121 BB   SER A  57      22.530  18.578  12.777  1.00  0.00\n"
    "ATOM    122 SC1  SER A  57      23.548  16.794  11.497  1.00  0.00\n"
    "ATOM    123 BB   ASP A  58      21.786  18.228  15.811  1.00  0.00\n"
    "ATOM    124 SC1  ASP A  58      24.174  16.003  16.613  1.00  0.00\n"
    "ATOM    125 BB   TYR A  59      20.601  21.380  15.641  1.00  0.00\n"
    "ATOM    126 SC1  TYR A  59      22.515  22.102  17.183  1.00  0.00\n"
    "ATOM    127 SC2  TYR A  59      24.662  21.414  18.210  1.00  0.00\n"
    "ATOM    128 SC3  TYR A  59      22.820  22.032  19.629  1.00  0.00\n"
    "ATOM    129 SC4  TYR A  59      24.935  21.374  20.693  1.00  0.00\n"
    "ATOM    130 BB   ASN A  60      19.243  22.209  12.952  1.00  0.00\n"
    "ATOM    131 SC1  ASN A  60      16.554  21.255  12.303  1.00  0.00\n"
    "ATOM    132 BB   ILE A  61      21.241  23.577  10.907  1.00  0.00\n"
    "ATOM    133 SC1  ILE A  61      23.203  24.990  12.869  1.00  0.00\n"
    "ATOM    134 BB   GLN A  62      20.668  25.662   8.252  1.00  0.00\n"
    "ATOM    135 SC1  GLN A  62      16.876  23.760   7.212  1.00  0.00\n"
    "ATOM    136 BB   LYS A  63      20.846  27.381   5.691  1.00  0.00\n"
    "ATOM    137 SC1  LYS A  63      20.703  26.551   2.854  1.00  0.00\n"
    "ATOM    138 SC2  LYS A  63      18.492  26.227   0.990  1.00  0.00\n"
    "ATOM    139 BB   GLU A  64      21.610  30.746   6.606  1.00  0.00\n"
    "ATOM    140 SC1  GLU A  64      21.159  30.803   2.737  1.00  0.00\n"
    "ATOM    141 BB   SER A  65      22.366  30.729   9.345  1.00  0.00\n"
    "ATOM    142 SC1  SER A  65      21.012  28.422  10.185  1.00  0.00\n"
    "ATOM    143 BB   THR A  66      22.828  32.169  12.668  1.00  0.00\n"
    "ATOM    144 SC1  THR A  66      22.926  34.695  11.648  1.00  0.00\n"
    "ATOM    145 BB   LEU A  67      25.477  32.626  14.755  1.00  0.00\n"
    "ATOM    146 SC1  LEU A  67      26.234  29.331  14.277  1.00  0.00\n"
    "ATOM    147 BB   HIS A  68      26.716  33.258  18.001  1.00  0.00\n"
    "ATOM    148 SC1  HIS A  68      24.596  34.843  18.453  1.00  0.00\n"
    "ATOM    149 SC2  HIS A  68      22.405  35.021  17.454  1.00  0.00\n"
    "ATOM    150 SC3  HIS A  68      23.309  36.522  17.513  1.00  0.00\n"
    "ATOM    151 BB   LEU A  69      29.616  35.063  19.208  1.00  0.00\n"
    "ATOM    152 SC1  LEU A  69      32.296  33.940  18.057  1.00  0.00\n"
    "ATOM    153 BB   VAL A  70      31.210  34.555  22.223  1.00  0.00\n"
    "ATOM    154 SC1  VAL A  70      29.114  35.531  23.778  1.00  0.00\n"
    "ATOM    155 BB   LEU A  71      33.678  35.753  24.378  1.00  0.00\n"
    "ATOM    156 SC1  LEU A  71      35.644  36.178  21.620  1.00  0.00\n"
    "ATOM    157 BB   ARG A  72      36.072  33.977  26.534  1.00  0.00\n"
    "ATOM    158 SC1  ARG A  72      35.360  31.797  27.830  1.00  0.00\n"
    "ATOM    159 SC2  ARG A  72      34.377  29.113  26.651  1.00  0.00\n"
    "ATOM    160 BB   LEU A  73      38.267  35.183  28.623  1.00  0.00\n"
    "ATOM    161 SC1  LEU A  73      39.151  37.446  26.048  1.00  0.00\n"
    "ATOM    162 BB   ARG A  74      41.535  34.542  30.191  1.00  0.00\n"
    "ATOM    163 SC1  ARG A  74      41.718  31.351  30.456  1.00  0.00\n"
    "ATOM    164 SC2  ARG A  74      44.179  29.920  29.559  1.00  0.00\n"
    "ATOM    165 BB   GLY A  75      41.317  37.180  32.109  1.00  0.00\n"
    "ATOM    166 BB   GLY A  76      40.212  39.684  35.087  1.00  0.00\n"
    "TER     167      GLY A  76\n"
    "END\n";

// The DSSP of the atomistic structure the beads were mapped from (mdlib's, in 3 states): strands 1-7, 10-17, 40-45,
// 48-49 and 66-71, the alpha helix 23-34 and the 3-10 helices 38-40 and 57-59.
static const char ubiquitin_atomistic_ss[] = "EEEEEEECCCCEEEEEECCCCEHHHHHHHHHHHHCCCHHHEEEEECCEECCCCCECHHHCCCCCCEEEEEECCCCC";

static bool load_martini3_ubiquitin(md_system_t* sys, md_system_state_t* state) {
    if (!md_pdb_system_init_from_str(sys, state, str_from_cstr(martini3_ubiquitin_pdb), MD_PDB_OPTION_NONE)) return false;
    return md_util_system_infer(sys, state, MD_UTIL_INFER_ALL);
}

// The BB beads of a Martini protein make a trace: one segment per residue, the bead in ca and no atoms
UTEST(secondary_structure, martini_trace_backbone) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_martini3_ubiquitin(&sys, &state));

    const md_protein_backbone_data_t* bb = &sys.protein_backbone;
    ASSERT_EQ((size_t)1, bb->range.count);
    ASSERT_EQ((size_t)76, bb->segment.count);
    EXPECT_TRUE(md_protein_backbone_range_is_trace(bb, 0));
    for (size_t i = 0; i < bb->segment.count; ++i) {
        const md_amino_acid_atoms_t a = bb->segment.atoms[i];
        EXPECT_EQ((md_component_idx_t)i, bb->segment.comp_idx[i]);
        EXPECT_TRUE(str_eq(md_atom_name(&sys.atom, a.ca), STR_LIT("BB")));
        EXPECT_EQ(-1, a.n);
        EXPECT_EQ(-1, a.c);
        EXPECT_EQ(-1, a.o);
        EXPECT_EQ(-1, a.hn);
        EXPECT_EQ(-1, a.cb);
    }

    // A trace has no phi and psi, and so no Ramachandran type
    md_backbone_angles_t angles[76];
    EXPECT_TRUE(md_util_backbone_angles_compute(angles, 76, state.xyz, &state.unitcell, bb));
    md_ramachandran_type_t types[76];
    EXPECT_TRUE(md_util_backbone_ramachandran_classify(types, 76, &sys));
    for (size_t i = 0; i < 76; ++i) {
        EXPECT_EQ(0.0f, angles[i].phi);
        EXPECT_EQ(0.0f, angles[i].psi);
        EXPECT_EQ(MD_RAMACHANDRAN_TYPE_UNKNOWN, types[i]);
    }

    md_arena_allocator_destroy(alloc);
}

// CG-DSSP on the beads against the DSSP of the atoms they were mapped from. On this structure 68 of the 76 residues agree
// (89 %): the helices, the five strands, and the coil between them; a residue more or less at the ends of strands and
// helices is all that differs.
UTEST(secondary_structure, martini_cg_dssp) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_martini3_ubiquitin(&sys, &state));
    ASSERT_EQ((size_t)76, sys.protein_backbone.segment.count);

    md_secondary_structure_t ss[76];
    ASSERT_TRUE(md_util_backbone_secondary_structure_infer(ss, 76, state.xyz, &state.unitcell, &sys.protein_backbone));

    char three[77] = {0};
    size_t agree = 0;
    for (size_t i = 0; i < 76; ++i) {
        three[i] = ss_three_state(ss[i]);
        agree += three[i] == ubiquitin_atomistic_ss[i];
    }
    printf("CG-DSSP  %s\nDSSP     %s\n", three, ubiquitin_atomistic_ss);
    EXPECT_GE(agree, (size_t)66);

    // The alpha helix and the first hairpin, whole
    for (size_t i = 23; i < 33; ++i) EXPECT_EQ('H', three[i]);
    for (size_t i = 2; i < 6; ++i) EXPECT_EQ('E', three[i]);
    for (size_t i = 11; i < 15; ++i) EXPECT_EQ('E', three[i]);

    md_arena_allocator_destroy(alloc);
}

// A trace split by the periodic boundary gets the secondary structure it has whole
UTEST(secondary_structure, martini_cg_dssp_periodic) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_martini3_ubiquitin(&sys, &state));
    const size_t n = sys.protein_backbone.segment.count;
    ASSERT_EQ((size_t)76, n);

    md_secondary_structure_t whole[76];
    const md_unitcell_t no_cell = {0};
    ASSERT_TRUE(md_util_backbone_secondary_structure_infer(whole, n, state.xyz, &no_cell, &sys.protein_backbone));

    // Centered on a corner of a 60 A box, and wrapped into it: every axis cuts the protein in two
    vec3_t com = {0};
    for (size_t i = 0; i < state.num_atoms; ++i) com = vec3_add(com, state.xyz[i]);
    com = vec3_mul1(com, 1.0f / (float)state.num_atoms);
    for (size_t i = 0; i < state.num_atoms; ++i) state.xyz[i] = vec3_sub(state.xyz[i], com);
    const md_unitcell_t cell = md_unitcell_from_extent(60.0, 60.0, 60.0);
    ASSERT_TRUE(md_util_pbc(state.xyz, NULL, state.num_atoms, &cell));

    md_secondary_structure_t split[76];
    ASSERT_TRUE(md_util_backbone_secondary_structure_infer(split, n, state.xyz, &cell, &sys.protein_backbone));
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(whole[i], split[i]);
    }

    md_arena_allocator_destroy(alloc);
}

// The atomistic DSSP of 1LAF as it was before traces were added beside it: atomistic ranges keep their bonds and patterns
UTEST(secondary_structure, atomistic_dssp_unchanged) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1LAF.pdb"), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));

    static const char expected[] = "CCCCEEEEEECCCBTTTBBBSTTTSBBSHHHHHHHHHHHHHTTEEEEEESSGGGHHHHHHTTTCCBBCCCCBSSHHHHHHTBBCCCCBCCCEEEEEETTTCCCSSHHHHTTTEEEEETTTHHHHHHIIIIIGGTTEEEEECCHHHHHHHHHTTTCCEEEEEHHHHIIIIITTGGGGGEEECCCCBCSHHHHTCCBCCBBSTTTHHHHHHHHHHHHHHHHHTHHHHHHHGTTSCCTTTC";
    const size_t n = sys.protein_backbone.segment.count;
    ASSERT_EQ(strlen(expected), n);
    ASSERT_EQ((size_t)1, sys.protein_backbone.range.count);
    EXPECT_FALSE(md_protein_backbone_range_is_trace(&sys.protein_backbone, 0));

    md_secondary_structure_t* ss = md_alloc(alloc, n * sizeof(md_secondary_structure_t));
    ASSERT_TRUE(md_util_backbone_secondary_structure_infer(ss, n, state.xyz, &state.unitcell, &sys.protein_backbone));
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(expected[i], ss_letters[ss[i]]);
    }

    md_arena_allocator_destroy(alloc);
}

// A protein split by the periodic boundary, residues and all (every atom wrapped into the cell on its own), gets the
// secondary structure it has whole: the bonds between residues in different images are those of the nearest image
UTEST(secondary_structure, atomistic_dssp_periodic) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1LAF.pdb"), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    const size_t n = sys.protein_backbone.segment.count;
    ASSERT_GT(n, (size_t)200);

    md_secondary_structure_t* whole = md_alloc(alloc, n * sizeof(md_secondary_structure_t));
    md_secondary_structure_t* split = md_alloc(alloc, n * sizeof(md_secondary_structure_t));
    const md_unitcell_t no_cell = {0};
    ASSERT_TRUE(md_util_backbone_secondary_structure_infer(whole, n, state.xyz, &no_cell, &sys.protein_backbone));

    // Centered on a corner of a cell wider than the protein and the reach of the bonds: every axis cuts it in two
    vec3_t com = {0};
    for (size_t i = 0; i < state.num_atoms; ++i) com = vec3_add(com, state.xyz[i]);
    com = vec3_mul1(com, 1.0f / (float)state.num_atoms);
    for (size_t i = 0; i < state.num_atoms; ++i) state.xyz[i] = vec3_sub(state.xyz[i], com);
    const md_unitcell_t cell = md_unitcell_from_extent(120.0, 120.0, 120.0);
    ASSERT_TRUE(md_util_pbc(state.xyz, NULL, state.num_atoms, &cell));

    ASSERT_TRUE(md_util_backbone_secondary_structure_infer(split, n, state.xyz, &cell, &sys.protein_backbone));
    size_t differ = 0;
    for (size_t i = 0; i < n; ++i) differ += whole[i] != split[i];
    EXPECT_EQ((size_t)0, differ);

    md_arena_allocator_destroy(alloc);
}

// A trace and an atomistic protein in one backbone, overlapping in space: each gets the secondary structure it has on
// its own, as the bonds of one kind are never looked for between residues of the other
UTEST(secondary_structure, mixed_trace_and_atomistic) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t cg = { .alloc = alloc };
    md_system_state_t cg_state = { .alloc = alloc };
    ASSERT_TRUE(load_martini3_ubiquitin(&cg, &cg_state));
    md_system_t aa = { .alloc = alloc };
    md_system_state_t aa_state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&aa, &aa_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1LAF.pdb"), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&aa, &aa_state, MD_UTIL_INFER_ALL));

    const md_protein_backbone_data_t* cg_bb = &cg.protein_backbone;
    const md_protein_backbone_data_t* aa_bb = &aa.protein_backbone;
    const size_t n_cg = cg_bb->segment.count;
    const size_t n_aa = aa_bb->segment.count;
    const size_t n = n_cg + n_aa;
    md_unitcell_t no_cell = {0};

    md_secondary_structure_t* ref = md_alloc(alloc, n * sizeof(md_secondary_structure_t));
    ASSERT_TRUE(md_util_backbone_secondary_structure_infer(ref, n_cg, cg_state.xyz, &no_cell, cg_bb));
    ASSERT_TRUE(md_util_backbone_secondary_structure_infer(ref + n_cg, n_aa, aa_state.xyz, &no_cell, aa_bb));

    // The beads first, then the atoms, both centered on the origin
    const size_t num_atoms = cg_state.num_atoms + aa_state.num_atoms;
    vec3_t* xyz = md_alloc(alloc, num_atoms * sizeof(vec3_t));
    vec3_t c_cg = {0}, c_aa = {0};
    for (size_t i = 0; i < cg_state.num_atoms; ++i) c_cg = vec3_add(c_cg, cg_state.xyz[i]);
    for (size_t i = 0; i < aa_state.num_atoms; ++i) c_aa = vec3_add(c_aa, aa_state.xyz[i]);
    c_cg = vec3_mul1(c_cg, 1.0f / (float)cg_state.num_atoms);
    c_aa = vec3_mul1(c_aa, 1.0f / (float)aa_state.num_atoms);
    for (size_t i = 0; i < cg_state.num_atoms; ++i) xyz[i] = vec3_sub(cg_state.xyz[i], c_cg);
    for (size_t i = 0; i < aa_state.num_atoms; ++i) xyz[cg_state.num_atoms + i] = vec3_sub(aa_state.xyz[i], c_aa);

    const size_t num_ranges = cg_bb->range.count + aa_bb->range.count;
    uint32_t* offset = md_alloc(alloc, (num_ranges + 1) * sizeof(uint32_t));
    md_amino_acid_atoms_t* atoms = md_alloc(alloc, n * sizeof(md_amino_acid_atoms_t));
    for (size_t r = 0; r < cg_bb->range.count; ++r) offset[r] = cg_bb->range.offset[r];
    for (size_t r = 0; r <= aa_bb->range.count; ++r) offset[cg_bb->range.count + r] = (uint32_t)n_cg + aa_bb->range.offset[r];
    for (size_t i = 0; i < n_cg; ++i) atoms[i] = cg_bb->segment.atoms[i];
    const md_atom_idx_t shift = (md_atom_idx_t)cg_state.num_atoms;
    for (size_t i = 0; i < n_aa; ++i) {
        md_amino_acid_atoms_t a = aa_bb->segment.atoms[i];
        a.n  = a.n  < 0 ? a.n  : a.n  + shift;
        a.ca = a.ca < 0 ? a.ca : a.ca + shift;
        a.c  = a.c  < 0 ? a.c  : a.c  + shift;
        a.o  = a.o  < 0 ? a.o  : a.o  + shift;
        a.cb = a.cb < 0 ? a.cb : a.cb + shift;
        a.hn = a.hn < 0 ? a.hn : a.hn + shift;
        atoms[n_cg + i] = a;
    }
    md_protein_backbone_data_t mixed = {0};
    mixed.range.count = num_ranges;
    mixed.range.offset = offset;
    mixed.segment.count = n;
    mixed.segment.atoms = atoms;
    EXPECT_TRUE(md_protein_backbone_range_is_trace(&mixed, 0));
    EXPECT_FALSE(md_protein_backbone_range_is_trace(&mixed, num_ranges - 1));

    md_secondary_structure_t* ss = md_alloc(alloc, n * sizeof(md_secondary_structure_t));
    ASSERT_TRUE(md_util_backbone_secondary_structure_infer(ss, n, xyz, &no_cell, &mixed));
    size_t differ = 0;
    for (size_t i = 0; i < n; ++i) differ += ss[i] != ref[i];
    EXPECT_EQ((size_t)0, differ);

    md_arena_allocator_destroy(alloc);
}
