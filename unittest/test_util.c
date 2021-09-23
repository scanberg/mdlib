#include "utest.h"
#include <string.h>

#include <md_pdb.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <md_util.h>
#include <core/md_allocator.h>
#include <core/md_file.h>

#include "rmsd.h"

UTEST(util, rmsd) {
    md_allocator_i* alloc = default_allocator;
    const str_t path = make_cstr(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {};
    ASSERT_TRUE(md_pdb_data_parse_file(path, &pdb_data, alloc));

    md_molecule_t mol = {};
    EXPECT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));

    md_trajectory_i traj = {};
    EXPECT_TRUE(md_pdb_trajectory_open(&traj, path, alloc));

    const int64_t stride = mol.atom.count;
    const int64_t mem_size = stride * 6 * sizeof(float) + stride * 6 * sizeof(double);
    float* mem = md_alloc(alloc, mem_size);
    float* x0 = mem + stride * 0;
    float* y0 = mem + stride * 1;
    float* z0 = mem + stride * 2;
    float* x1 = mem + stride * 3;
    float* y1 = mem + stride * 4;
    float* z1 = mem + stride * 5;

    double* xyz0 = (double*)(mem + stride * 6);
    double* xyz1 = (double*)(mem + stride * 12);

    md_trajectory_load_frame(&traj, 0, NULL, x0, y0, z0, mol.atom.count);
    md_trajectory_load_frame(&traj, 1, NULL, x1, y1, z1, mol.atom.count);

    for (int64_t i = 0; i < mol.atom.count; ++i) {
        xyz0[i * 3 + 0] = x0[i];
        xyz0[i * 3 + 1] = y0[i];
        xyz0[i * 3 + 2] = z0[i];

        xyz1[i * 3 + 0] = x1[i];
        xyz1[i * 3 + 1] = y1[i];
        xyz1[i * 3 + 2] = z1[i];
    }

    // Reference
    double ref_rmsd;
    fast_rmsd((double(*)[3])xyz0, (double(*)[3])xyz1, (int)mol.atom.count, &ref_rmsd);

    // Our implementation
    double rmsd = md_util_compute_rmsd(x0, y0, z0, x1, y1, z1, mol.atom.mass, mol.atom.count);

    md_pdb_molecule_free(&mol, alloc);
    md_pdb_data_free(&pdb_data, alloc);
}