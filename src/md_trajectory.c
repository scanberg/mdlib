#include "md_trajectory.h"

#include <stdlib.h> // malloc
#include <stdarg.h> // VA_ARGS
#include <string.h>

#ifndef ASSERT
#include <assert.h>
#define ASSERT assert
#endif

#include <core/file.h>
#include <ext/xtc/xdrfile_xtc.h>
#include <ext/dcd/readdcd.h>

#ifdef __cplusplus
extern "C" {
#endif

#define internal static

#define MD_TRAJECTORY_MAGIC 0xbadb00b5

/*
* 
typedef uint32_t md_trajectory_type;

enum {
    TRAJECTORY_TYPE_UNDEFINED,
    TRAJECTORY_TYPE_PDB,
    TRAJECTORY_TYPE_XTC,
    TRAJECTORY_TYPE_DCD,
};

struct md_trajectory {
    uint32_t magic;
    md_trajectory_type type;
    uint32_t num_atoms;
    uint32_t num_frames;
    float start_time;           // Probably almost always zero?
    float frame_delta_time;
    md_file* file;
    uint64_t* frame_offsets;
};

bool str_equal_ignore_case(const char *a, const char *b) {
    int ca, cb;
    do {
        ca = (unsigned char) *a++;
        cb = (unsigned char) *b++;
        ca = tolower(ca);
        cb = tolower(cb);
    } while (ca == cb && ca != '\0');
    return ca == cb;
}

internal void log_error(md_log_i* err_log, const char* format, ...) {
    char buf[1024];

    va_list ap;
    va_start(ap, format);

    if (err_log) {
        int len = vsnprintf(buf, ARRAY_SIZE(buf), format, ap);
        err_log->print(MD_LOG_SEVERITY_ERROR, buf, len);
    } else {
        vfprintf(stderr, format, ap);
    }

    va_end(ap);
}

md_error md_trajectory_init(md_trajectory* traj, const char* filename, uint64_t cache_size, md_context* ctx) {
    // Determine trajectory type based on file ending
    int len = strlen(filename);
    const char* ext = strrchr(filename, '.');
    const md_log_i* err_log = ctx ? &ctx->print : NULL;
    if (!ext) {
        log_error(err_log, "Could not resolve extension of supplied filename: '%s'", filename);
        return MD_TRAJECTORY_INVALID_FILE;
    }

    md_trajectory_type type = TRAJECTORY_TYPE_UNDEFINED;
    if (str_equal_ignore_case(ext, "pdb")) {
        type = TRAJECTORY_TYPE_PDB;
    }
    else if (str_equal_ignore_case(ext, "xtc")) {
        type = TRAJECTORY_TYPE_XTC;
    }
    else if (str_equal_ignore_case(ext, "dcd")) {
        type = TRAJECTORY_TYPE_DCD;
    }

    if (type == TRAJECTORY_TYPE_UNDEFINED) {
        log_error(err_log, "Unsuported file extension for trajectory: '%s'", ext);
        return MD_TRAJECTORY_INVALID_FILE;
    }

    // Read or Generate frame index cache
    char cache_file[512];
    snprintf(cache_file, sizeof(cache_file), "%s.cache");
}

md_error md_trajectory_free(md_trajectory* traj) {

}

md_error md_trajectory_read_header(const md_trajectory* traj, md_trajectory_header* header) {
    ASSERT(traj);
    ASSERT(header);
    header->num_atoms = traj->num_atoms;
    header->num_frames = traj->num_frames;
    header->start_time = traj->start_time;
    header->frame_delta_time = traj->frame_delta_time;
    return MD_TRAJECTORY_SUCCESS;
}

md_error md_trajectory_load_frame(md_trajectory* traj, uint32_t frame_idx, md_trajectory_target_storage* target) {

}

*/

#ifdef __cplusplus
}
#endif