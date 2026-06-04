#include "utest.h"

#include <md_gto.h>
#include <core/md_gpu.h>
#include <core/md_vec_math.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

static md_gpu_event_t gpu_test_upload_buffer(md_gpu_device_t device, md_gpu_buffer_t dst_buffer, size_t dst_offset, const void* src_data, size_t size) {
    md_gpu_event_t event = {0};
    if (!device || !dst_buffer || !src_data || size == 0) return event;

    md_gpu_queue_t queue = md_gpu_queue_transfer(device);
    if (!queue) return event;

    md_gpu_cmd_t cmd = md_gpu_cmd_begin(queue, "gpu test upload buffer");
    if (!cmd) return event;

    md_gpu_transient_t upload = md_gpu_cmd_temp_alloc(cmd, size);
    if (!upload.buffer || !upload.cpu_ptr || upload.size < size) {
        md_gpu_cmd_discard(cmd);
        return event;
    }

    memcpy(upload.cpu_ptr, src_data, size);
    if (!md_gpu_cmd_copy_buffer(cmd, upload.buffer, dst_buffer, size, upload.offset, dst_offset) || !md_gpu_cmd_end(cmd)) {
        md_gpu_cmd_discard(cmd);
        return event;
    }

    return md_gpu_queue_submit_one(queue, cmd);
}

static md_gpu_event_t gpu_test_upload_image(md_gpu_device_t device, md_gpu_image_t dst_image, const void* src_data, size_t size) {
    md_gpu_event_t event = {0};
    if (!device || !dst_image || !src_data || size == 0) return event;

    md_gpu_queue_t queue = md_gpu_queue_transfer(device);
    if (!queue) return event;

    md_gpu_cmd_t cmd = md_gpu_cmd_begin(queue, "gpu test upload image");
    if (!cmd) return event;

    md_gpu_transient_t upload = md_gpu_cmd_temp_alloc(cmd, size);
    if (!upload.buffer || !upload.cpu_ptr || upload.size < size) {
        md_gpu_cmd_discard(cmd);
        return event;
    }

    memcpy(upload.cpu_ptr, src_data, size);
    if (!md_gpu_cmd_copy_buffer_to_image(cmd, upload.buffer, dst_image) || !md_gpu_cmd_end(cmd)) {
        md_gpu_cmd_discard(cmd);
        return event;
    }

    return md_gpu_queue_submit_one(queue, cmd);
}

/* =========================================================
 * gpu.device
 * Basic device lifecycle: create and destroy without crashing.
 * ========================================================= */

UTEST(gpu, device_create_destroy) {
    md_gpu_device_t device = md_gpu_device_create();
    if (!device) {
        UTEST_SKIP("No GPU device available");
    }
    md_gpu_device_destroy(device);
}

/* =========================================================
 * gpu.buffer_upload_download
 * Round-trip: upload a known pattern into a GPU-only buffer via an
 * explicit transfer command, then read it back with md_gpu_readback_buffer
 * and verify every byte matches.
 * ========================================================= */

UTEST(gpu, buffer_upload_download) {
    md_gpu_device_t device = md_gpu_device_create();
    if (!device) {
        UTEST_SKIP("No GPU device available");
    }

    const size_t N = 1024;
    const size_t buf_size = N * sizeof(uint32_t);

    /* Source pattern. */
    uint32_t src[1024];
    for (size_t i = 0; i < N; ++i) {
        src[i] = (uint32_t)(i * 2654435761u); /* Knuth multiplicative hash */
    }

    /* GPU buffer (device-local; no CPU_VISIBLE flag so the transfer path is exercised). */
    md_gpu_buffer_t gpu_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size  = buf_size,
        .flags = MD_GPU_BUFFER_NONE,
    });
    ASSERT_TRUE(gpu_buf != NULL);

    /* Upload. */
    md_gpu_event_t up_event = gpu_test_upload_buffer(device, gpu_buf, 0, src, buf_size);
    ASSERT_TRUE(up_event.value != 0);

    /* Read back through the new async readback handle API. */
    md_gpu_readback_t readback = md_gpu_readback_buffer(gpu_buf, 0, buf_size, up_event);
    ASSERT_TRUE(readback != NULL);
    md_gpu_readback_wait(readback);
    const uint32_t* dst = (const uint32_t*)md_gpu_readback_cpu_ptr(readback);
    ASSERT_TRUE(dst != NULL);

    /* Verify. */
    for (size_t i = 0; i < N; ++i) {
        EXPECT_EQ(src[i], dst[i]);
    }

    md_gpu_readback_destroy(readback);
    md_gpu_buffer_destroy(gpu_buf);
    md_gpu_device_destroy(device);
}

/* =========================================================
 * gpu.buffer_upload_download_with_offset
 * Upload into the middle of a buffer and download back with a
 * matching offset to verify the offset arithmetic is correct.
 * ========================================================= */

UTEST(gpu, buffer_upload_download_with_offset) {
    md_gpu_device_t device = md_gpu_device_create();
    if (!device) {
        UTEST_SKIP("No GPU device available");
    }

    const size_t ELEM = 64;
    const size_t OFFSET_ELEM = 16;
    const size_t buf_size = ELEM * sizeof(uint32_t);
    const size_t offset_bytes = OFFSET_ELEM * sizeof(uint32_t);
    const size_t payload_elems = ELEM - OFFSET_ELEM;
    const size_t payload_bytes = payload_elems * sizeof(uint32_t);

    uint32_t src[64 - 16]; /* payload only */
    for (size_t i = 0; i < payload_elems; ++i) {
        src[i] = (uint32_t)(i + 1) * 0xDEADBEEFu;
    }

    md_gpu_buffer_t gpu_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size  = buf_size,
        .flags = MD_GPU_BUFFER_NONE,
    });
    ASSERT_TRUE(gpu_buf != NULL);

    md_gpu_event_t up_event = gpu_test_upload_buffer(device, gpu_buf, offset_bytes, src, payload_bytes);
    ASSERT_TRUE(up_event.value != 0);

    md_gpu_readback_t readback = md_gpu_readback_buffer(gpu_buf, offset_bytes, payload_bytes, up_event);
    ASSERT_TRUE(readback != NULL);
    md_gpu_readback_wait(readback);
    const uint32_t* dst = (const uint32_t*)md_gpu_readback_cpu_ptr(readback);
    ASSERT_TRUE(dst != NULL);

    for (size_t i = 0; i < payload_elems; ++i) {
        EXPECT_EQ(src[i], dst[i]);
    }

    md_gpu_readback_destroy(readback);
    md_gpu_buffer_destroy(gpu_buf);
    md_gpu_device_destroy(device);
}

/* =========================================================
 * gpu.image_upload_download
 * Upload a known float pattern into an R32_FLOAT storage image via an
 * explicit transfer command, then read it back via md_gpu_readback_image
 * and verify all values match.
 * ========================================================= */

UTEST(gpu, image_upload_download) {
    md_gpu_device_t device = md_gpu_device_create();
    if (!device) {
        UTEST_SKIP("No GPU device available");
    }

    const uint32_t W = 8, H = 8, D = 8;
    const size_t voxel_count = W * H * D;
    const size_t img_size = voxel_count * sizeof(float);

    float src[8 * 8 * 8];
    for (size_t i = 0; i < voxel_count; ++i) {
        src[i] = (float)i * 0.001f;
    }

    md_gpu_image_t img = md_gpu_image_create(device, &(md_gpu_image_desc_t){
        .width  = W,
        .height = H,
        .depth  = D,
        .format = MD_GPU_IMAGE_FORMAT_R32_FLOAT,
        .flags  = MD_GPU_IMAGE_STORAGE,
    });
    ASSERT_TRUE(img != NULL);

    /* Upload. */
    md_gpu_event_t up_event = gpu_test_upload_image(device, img, src, img_size);
    ASSERT_TRUE(up_event.value != 0);

    /* Download. */
    md_gpu_readback_t readback = md_gpu_readback_image(img, (md_gpu_image_region_t){0}, up_event);
    ASSERT_TRUE(readback != NULL);
    md_gpu_readback_wait(readback);
    const float* dst = (const float*)md_gpu_readback_cpu_ptr(readback);
    ASSERT_TRUE(dst != NULL);

    /* Verify — floating-point exact equality is fine here as no GPU math is involved. */
    float max_delta = 0.0f;
    for (size_t i = 0; i < voxel_count; ++i) {
        float d = fabsf(dst[i] - src[i]);
        if (d > max_delta) max_delta = d;
    }
    EXPECT_TRUE(max_delta == 0.0f);

    md_gpu_readback_destroy(readback);
    md_gpu_image_destroy(img);
    md_gpu_device_destroy(device);
}

UTEST(gpu, transient_scratch_buffer_copy) {
    md_gpu_device_t device = md_gpu_device_create();
    if (!device) {
        UTEST_SKIP("No GPU device available");
    }

    const size_t N = 256;
    const size_t buf_size = N * sizeof(uint32_t);
    uint32_t src[256];
    uint32_t dst[256];
    for (size_t i = 0; i < N; ++i) {
        src[i] = (uint32_t)(i ^ 0xA5A5A5A5u);
    }
    memset(dst, 0, sizeof(dst));

    md_gpu_buffer_t gpu_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size = buf_size,
        .flags = MD_GPU_BUFFER_NONE,
    });
    ASSERT_TRUE(gpu_buf != NULL);

    md_gpu_queue_t queue = md_gpu_queue_transfer(device);
    ASSERT_TRUE(queue != NULL);

    md_gpu_cmd_t upload_cmd = md_gpu_cmd_begin(queue, "transient scratch upload");
    ASSERT_TRUE(upload_cmd != NULL);
    md_gpu_transient_t up = md_gpu_cmd_temp_alloc(upload_cmd, buf_size);
    ASSERT_TRUE(up.buffer != NULL);
    ASSERT_TRUE(up.cpu_ptr != NULL);
    ASSERT_EQ(buf_size, up.size);
    memcpy(up.cpu_ptr, src, buf_size);
    ASSERT_TRUE(md_gpu_cmd_copy_buffer(upload_cmd, up.buffer, gpu_buf, up.size, up.offset, 0));
    ASSERT_TRUE(md_gpu_cmd_end(upload_cmd));
    md_gpu_event_t upload_event = md_gpu_queue_submit_one(queue, upload_cmd);
    ASSERT_TRUE(upload_event.value != 0);

    md_gpu_cmd_t readback_cmd = md_gpu_cmd_begin(queue, "transient scratch readback");
    ASSERT_TRUE(readback_cmd != NULL);
    md_gpu_transient_t rb = md_gpu_cmd_temp_alloc(readback_cmd, buf_size);
    ASSERT_TRUE(rb.buffer != NULL);
    ASSERT_TRUE(rb.cpu_ptr != NULL);
    ASSERT_TRUE(md_gpu_cmd_copy_buffer(readback_cmd, gpu_buf, rb.buffer, rb.size, 0, rb.offset));
    ASSERT_TRUE(md_gpu_cmd_end(readback_cmd));
    md_gpu_event_t readback_event = md_gpu_queue_submit_one_after(queue, readback_cmd, upload_event, MD_GPU_BARRIER_STAGE_TRANSFER);
    ASSERT_TRUE(readback_event.value != 0);
    md_gpu_event_wait(readback_event);

    memcpy(dst, rb.cpu_ptr, buf_size);
    for (size_t i = 0; i < N; ++i) {
        EXPECT_EQ(src[i], dst[i]);
    }

    md_gpu_buffer_destroy(gpu_buf);
    md_gpu_device_destroy(device);
}

UTEST(gpu, transient_scratch_large_allocation) {
    md_gpu_device_t device = md_gpu_device_create();
    if (!device) {
        UTEST_SKIP("No GPU device available");
    }

    const size_t size = (size_t)(3 * 1024 * 1024);
    md_gpu_queue_t queue = md_gpu_queue_transfer(device);
    ASSERT_TRUE(queue != NULL);

    md_gpu_cmd_t cmd = md_gpu_cmd_begin(queue, "transient scratch large allocation");
    ASSERT_TRUE(cmd != NULL);

    md_gpu_transient_t a = md_gpu_cmd_temp_alloc(cmd, size);
    md_gpu_transient_t b = md_gpu_cmd_temp_alloc(cmd, size);
    ASSERT_TRUE(a.buffer != NULL);
    ASSERT_TRUE(b.buffer != NULL);
    ASSERT_TRUE(a.cpu_ptr != NULL);
    ASSERT_TRUE(b.cpu_ptr != NULL);
    EXPECT_TRUE(a.buffer != b.buffer);

    memset(a.cpu_ptr, 0x11, size);
    memset(b.cpu_ptr, 0x22, size);

    ASSERT_TRUE(md_gpu_cmd_end(cmd));
    md_gpu_event_t event = md_gpu_queue_submit_one(queue, cmd);
    ASSERT_TRUE(event.value != 0);
    md_gpu_event_wait(event);

    md_gpu_device_destroy(device);
}

/* =========================================================
 * gpu.gto_mo_root
 * Evaluate a single s-type GTO (alpha=1, coeff=1) centred at
 * the origin on a 4^3 grid and verify the GPU output matches
 * the expected analytical values within a tight tolerance.
 * (Migrated from the former gpu_pass test suite.)
 * ========================================================= */

UTEST(gpu, gto_mo_root) {
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
        .basis  = &basis,
        .cutoff = 0.0,
    });
    ASSERT_TRUE(gpu_basis != NULL);
    ASSERT_EQ(1u, md_gto_gpu_basis_num_cgtos(gpu_basis));

    md_gpu_buffer_t atom_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size  = md_gto_gpu_atom_buffer_size(1),
        .flags = MD_GPU_BUFFER_CPU_VISIBLE,
    });
    md_gpu_buffer_t coeff_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size  = md_gto_gpu_coeff_size_mo(1, 1),
        .flags = MD_GPU_BUFFER_CPU_VISIBLE,
    });
    ASSERT_TRUE(atom_buf != NULL);
    ASSERT_TRUE(coeff_buf != NULL);

    float atom_xyz[3] = {0.0f, 0.0f, 0.0f};
    md_gto_gpu_atom_pack((float*)md_gpu_buffer_cpu_ptr(atom_buf), atom_xyz, 0, 1);
    ((float*)md_gpu_buffer_cpu_ptr(coeff_buf))[0] = 1.0f;

    md_grid_t grid = {
        .orientation = mat3_ident(),
        .origin      = {-1.0f, -1.0f, -1.0f},
        .spacing     = {0.5f, 0.5f, 0.5f},
        .dim         = {4, 4, 4},
    };
    const size_t voxel_count = (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2];
    const size_t readback_size = sizeof(float) * voxel_count;

    md_gpu_image_t out_image = md_gpu_image_create(device, &(md_gpu_image_desc_t){
        .width  = (uint32_t)grid.dim[0],
        .height = (uint32_t)grid.dim[1],
        .depth  = (uint32_t)grid.dim[2],
        .format = MD_GPU_IMAGE_FORMAT_R32_FLOAT,
        .flags  = MD_GPU_IMAGE_STORAGE,
    });
    md_gpu_buffer_t readback_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size  = readback_size,
        .flags = MD_GPU_BUFFER_CPU_VISIBLE,
    });
    ASSERT_TRUE(out_image != NULL);
    ASSERT_TRUE(readback_buf != NULL);

    md_gpu_queue_t queue = md_gpu_queue_compute(device);
    ASSERT_TRUE(queue != NULL);

    md_gpu_cmd_t eval_cmd = md_gpu_cmd_begin(queue, "gto mo root eval");
    ASSERT_TRUE(eval_cmd != NULL);
    md_gto_gpu_mo_record(eval_cmd, gpu_basis, atom_buf, coeff_buf, 1, out_image, &grid, MD_GTO_EVAL_MODE_PSI, MD_GTO_OP_SET);
    ASSERT_TRUE(md_gpu_cmd_end(eval_cmd));
    md_gpu_event_t eval_event = md_gpu_queue_submit_one(queue, eval_cmd);
    ASSERT_TRUE(eval_event.value != 0);

    md_gpu_cmd_t copy_cmd = md_gpu_cmd_begin(queue, "gto mo root readback");
    ASSERT_TRUE(copy_cmd != NULL);
    md_gpu_cmd_copy_image_region_to_buffer(copy_cmd, out_image, (md_gpu_image_region_t){0}, readback_buf, 0);
    ASSERT_TRUE(md_gpu_cmd_end(copy_cmd));
    md_gpu_event_t copy_event = md_gpu_queue_submit_one_after(queue, copy_cmd, eval_event, MD_GPU_BARRIER_STAGE_TRANSFER);
    ASSERT_TRUE(copy_event.value != 0);
    md_gpu_event_wait(copy_event);

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
        printf("GTO root eval max delta: %g\n", (double)max_delta);
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
