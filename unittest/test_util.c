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
    const str_t path = MAKE_STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, path, alloc));

    md_molecule_t mol = {0};
    EXPECT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));

    md_trajectory_i* traj = md_pdb_trajectory_create(path, alloc);
    EXPECT_TRUE(traj);

    const int64_t stride = mol.atom.count;
    const int64_t mem_size = stride * 6 * sizeof(float) + stride * 6 * sizeof(double);
    float* mem = md_alloc(alloc, mem_size);
    md_vec3_soa_t coord[2] = {
        {
            mem + stride * 0,
            mem + stride * 1,
            mem + stride * 2,
        },
        {
            mem + stride * 3,
            mem + stride * 4,
            mem + stride * 5,
        },
    };

    double* xyz0 = (double*)(mem + stride * 6);
    double* xyz1 = (double*)(mem + stride * 6) + stride * 3;

    md_trajectory_load_frame(traj, 0, NULL, coord[0].x, coord[0].y, coord[0].z);
    md_trajectory_load_frame(traj, 1, NULL, coord[1].x, coord[1].y, coord[1].z);

    for (int64_t i = 0; i < mol.atom.count; ++i) {
        xyz0[i * 3 + 0] = coord[0].x[i];
        xyz0[i * 3 + 1] = coord[0].y[i];
        xyz0[i * 3 + 2] = coord[0].z[i];

        xyz1[i * 3 + 0] = coord[1].x[i];
        xyz1[i * 3 + 1] = coord[1].y[i];
        xyz1[i * 3 + 2] = coord[1].z[i];
    }

    // Reference
    double ref_rmsd;
    fast_rmsd((double(*)[3])xyz0, (double(*)[3])xyz1, (int)mol.atom.count, &ref_rmsd);

    // Our implementation
    vec3_t com[2] = {
        md_util_compute_com(coord[0].x, coord[0].y, coord[0].z, mol.atom.mass, mol.atom.count),
        md_util_compute_com(coord[1].x, coord[1].y, coord[1].z, mol.atom.mass, mol.atom.count),
    };
    double rmsd = md_util_compute_rmsd(coord, com, mol.atom.mass, mol.atom.count);
    
    EXPECT_LE(fabs(ref_rmsd - rmsd), 0.1);

    md_molecule_free(&mol, alloc);
    md_pdb_trajectory_free(traj);
    md_pdb_data_free(&pdb_data, alloc);
}