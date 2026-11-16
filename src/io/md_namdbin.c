#include <md_namdbin.h>

#include <md_system.h>

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_file.h>

#include <math.h>
#include <string.h>

// The atom count, and whether the file is of the other byte order; false when the size does not fit either
static bool namdbin_header(md_file_t file, size_t* out_atoms, bool* out_swap) {
    const uint64_t size = (uint64_t)md_file_size(file);
    uint32_t n = 0;
    if (size < 4 + 24 || md_file_read_at(file, 0, &n, 4) != 4) return false;
    const uint32_t swapped = md_bswap32(n);
    if (n > 0 && (uint64_t)n * 24 + 4 == size) {
        *out_atoms = n;
        *out_swap = false;
        return true;
    }
    if (swapped > 0 && (uint64_t)swapped * 24 + 4 == size) {
        *out_atoms = swapped;
        *out_swap = true;
        return true;
    }
    return false;
}

bool md_namdbin_file_is_coor(md_str_t filename) {
    return md_namdbin_file_atom_count(filename) > 0;
}

size_t md_namdbin_file_atom_count(md_str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return 0;
    size_t n = 0;
    bool swap = false;
    const bool ok = namdbin_header(file, &n, &swap);
    md_file_close(&file);
    return ok ? n : 0;
}

bool md_namdbin_state_init_from_file(md_system_state_t* state, md_str_t filename, size_t num_atoms) {
    MD_ASSERT(state);
    if (!state->alloc) {
        MD_LOG_ERROR("NAMD coordinates: the state has no allocator");
        return false;
    }
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("NAMD coordinates: cannot open '" MD_STR_FMT "'", MD_STR_ARG(filename));
        return false;
    }
    bool result = false;
    size_t N = 0;
    bool swap = false;
    double* xyz = NULL;
    if (!namdbin_header(file, &N, &swap)) {
        MD_LOG_ERROR("NAMD coordinates: '" MD_STR_FMT "' is not a NAMD binary coordinate file (its size is not that of its atom count)", MD_STR_ARG(filename));
        goto done;
    }
    if (num_atoms && N != num_atoms) {
        MD_LOG_ERROR("NAMD coordinates: '" MD_STR_FMT "' holds %zu atoms, the system %zu", MD_STR_ARG(filename), N, num_atoms);
        goto done;
    }
    xyz = md_alloc(md_get_heap_allocator(), N * 24);
    if (md_file_read_at(file, 4, xyz, N * 24) != N * 24) {
        MD_LOG_ERROR("NAMD coordinates: failed to read '" MD_STR_FMT "'", MD_STR_ARG(filename));
        goto done;
    }
    if (!md_system_state_init(state, N)) goto done;
    float* out = (float*)state->xyz;
    for (size_t i = 0; i < 3 * N; ++i) {
        uint64_t u;
        md_memcpy(&u, &xyz[i], 8);
        if (swap) u = md_bswap64(u);
        double v;
        md_memcpy(&v, &u, 8);
        out[i] = (float)v;
    }
    state->num_atoms = N;
    state->unitcell = md_unitcell_none();
    result = true;
done:
    if (xyz) md_free(md_get_heap_allocator(), xyz, N * 24);
    md_file_close(&file);
    return result;
}
