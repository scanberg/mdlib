/*
md_gpu_util.h -- removed.

This held md_gpu_bump_t, a staging bump allocator for the old md_gpu command
buffer API. md_gpu.h now owns staging directly:

    md_gpu_upload_begin / md_gpu_upload_end   zero-copy upload into a device
                                              pointer, staged only when the
                                              destination is busy
    md_gpu_memcpy_async                       host<->device, direction inferred

Nothing includes this header any more. The file remains only because the tool
that performed the migration could not delete it; it is safe to remove.
*/
#pragma once
