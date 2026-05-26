#include "utest.h"

#include <md_gto.h>
#include <core/md_gpu.h>
#include <core/md_vec_math.h>

#include <math.h>
#include <stdio.h>

UTEST(gpu_pass, gto_mo_root) {
    md_gpu_device_t device = md_gpu_device_create();
    if (!device) {
        UTEST_SKIP("No GPU device available");
    }

    md_gto_gpu_initialize(device);

    md_gto_shell_t shells[1] = {{
        .atom_idx = 0,
        .primitive_offset = 0,
        .num_primitives = 1,
        .l = 0,
    }};
    float alpha[1] = {1.0f};
    float coeff[1] = {1.0f};
    md_gto_basis_t basis = {
        .num_shells = 1,
        .num_primitives = 1,
        .shells = shells,
        .alpha = alpha,
        .coeff = coeff,
    };

    md_gto_gpu_basis_t gpu_basis = md_gto_gpu_basis_create(device, &(md_gto_gpu_basis_desc_t){
        .basis = &basis,
        .cutoff = 0.0,
    });
    ASSERT_TRUE(gpu_basis != NULL);
    ASSERT_EQ(1u, md_gto_gpu_basis_num_cgtos(gpu_basis));

    md_gpu_buffer_t atom_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size = md_gto_gpu_atom_buffer_size(1),
        .flags = MD_GPU_BUFFER_CPU_VISIBLE,
    });
    md_gpu_buffer_t coeff_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size = md_gto_gpu_coeff_size_mo(1, 1),
        .flags = MD_GPU_BUFFER_CPU_VISIBLE,
    });
    ASSERT_TRUE(atom_buf != NULL);
    ASSERT_TRUE(coeff_buf != NULL);

    float atom_xyz[3] = {0.0f, 0.0f, 0.0f};
    md_gto_gpu_atom_pack((float*)md_gpu_buffer_cpu_ptr(atom_buf), atom_xyz, 0, 1);
    ((float*)md_gpu_buffer_cpu_ptr(coeff_buf))[0] = 1.0f;

    md_grid_t grid = {
        .orientation = mat3_ident(),
        .origin = {-1.0f, -1.0f, -1.0f},
        .spacing = {0.5f, 0.5f, 0.5f},
        .dim = {4, 4, 4},
    };
    const size_t voxel_count = (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2];
    const size_t readback_size = sizeof(float) * voxel_count;

    md_gpu_image_t out_image = md_gpu_image_create(device, &(md_gpu_image_desc_t){
        .width = (uint32_t)grid.dim[0],
        .height = (uint32_t)grid.dim[1],
        .depth = (uint32_t)grid.dim[2],
        .format = MD_GPU_IMAGE_FORMAT_R32_FLOAT,
        .flags = MD_GPU_IMAGE_STORAGE,
    });
    md_gpu_buffer_t readback_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size = readback_size,
        .flags = MD_GPU_BUFFER_CPU_VISIBLE,
    });
    ASSERT_TRUE(out_image != NULL);
    ASSERT_TRUE(readback_buf != NULL);

    md_gpu_pass_buffer_usage_t buffers[4] = {
        { .buffer = md_gto_gpu_basis_buffer(gpu_basis), .usage = MD_GPU_PASS_RESOURCE_READ },
        { .buffer = atom_buf, .usage = MD_GPU_PASS_RESOURCE_READ },
        { .buffer = coeff_buf, .usage = MD_GPU_PASS_RESOURCE_READ },
        { .buffer = readback_buf, .usage = MD_GPU_PASS_RESOURCE_TRANSFER_DST },
    };
    md_gpu_pass_image_usage_t images[1] = {{
        .image = out_image,
        .usage = MD_GPU_PASS_RESOURCE_WRITE | MD_GPU_PASS_RESOURCE_TRANSFER_SRC,
    }};
    md_gpu_pass_desc_t pass_desc = {
        .label = "gto mo root pass test",
        .buffers = buffers,
        .buffer_count = 4,
        .images = images,
        .image_count = 1,
    };

    md_gpu_pass_t pass = md_gpu_pass_begin(device, &pass_desc);
    ASSERT_TRUE(pass != NULL);
    md_gto_gpu_mo_root_pass_record(pass, gpu_basis, atom_buf, coeff_buf, 1, out_image, &grid, MD_GTO_EVAL_MODE_PSI, MD_GTO_OP_SET);
    md_gpu_pass_barrier_image(pass, out_image, MD_GPU_BARRIER_STAGE_COMPUTE, MD_GPU_BARRIER_STAGE_TRANSFER);
    md_gpu_pass_copy_image_region_to_buffer(pass, out_image, (md_gpu_image_region_t){0}, readback_buf, 0);
    md_gpu_pass_id_t id = md_gpu_pass_end(pass);
    ASSERT_TRUE(id.value != 0);
    md_gpu_pass_wait(device, id);

    const float* values = (const float*)md_gpu_buffer_cpu_ptr(readback_buf);
    float max_delta = 0.0f;
    for (int z = 0; z < grid.dim[2]; ++z) {
        for (int y = 0; y < grid.dim[1]; ++y) {
            for (int x = 0; x < grid.dim[0]; ++x) {
                size_t idx = (size_t)x + (size_t)grid.dim[0] * ((size_t)y + (size_t)grid.dim[1] * (size_t)z);
                float px = grid.origin.x + ((float)x + 0.5f) * grid.spacing.x;
                float py = grid.origin.y + ((float)y + 0.5f) * grid.spacing.y;
                float pz = grid.origin.z + ((float)z + 0.5f) * grid.spacing.z;
                float expected = expf(-(px * px + py * py + pz * pz));
                float delta = fabsf(values[idx] - expected);
                if (delta > max_delta) max_delta = delta;
            }
        }
    }

    if (max_delta >= 5.0e-4f) {
        printf("GTO root pass max delta: %g\n", (double)max_delta);
    }
    EXPECT_TRUE(max_delta < 5.0e-4f);

    md_gpu_buffer_destroy(readback_buf);
    md_gpu_image_destroy(out_image);
    md_gpu_buffer_destroy(coeff_buf);
    md_gpu_buffer_destroy(atom_buf);
    md_gto_gpu_basis_destroy(gpu_basis);
    md_gto_gpu_shutdown();
    md_gpu_device_destroy(device);
}
