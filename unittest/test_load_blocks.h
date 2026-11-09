#pragma once

// For the tests of the loaders which read their lines in blocks and pieces (md_reader.h): what two loads of the
// same text are to agree on, and files large enough for several blocks

#include <md_system.h>
#include <core/md_os.h>
#include <core/md_array.h>
#include <string.h>
#include <stdio.h>

#define LOAD_BLOCKS_CHECK(cond, what) do { if (!(cond)) { printf("  load_blocks_diff: %s\n", what); d += 1; } } while (0)

// The number of things two loads disagree on, each told: 0 when they gave the same system
static size_t load_blocks_diff(const md_system_t* a, const md_system_state_t* as, const md_system_t* b, const md_system_state_t* bs) {
    size_t d = 0;
    if (a->atom.count != b->atom.count || a->residue.count != b->residue.count || a->atom.type.count != b->atom.type.count) {
        printf("  load_blocks_diff: %zu vs %zu atoms, %zu vs %zu residues, %zu vs %zu types\n", a->atom.count, b->atom.count,
               a->residue.count, b->residue.count, a->atom.type.count, b->atom.type.count);
        return 1;
    }
    const size_t n = a->atom.count, r = a->residue.count;
    LOAD_BLOCKS_CHECK(memcmp(as->xyz, bs->xyz, n * sizeof(vec3_t)) == 0, "xyz");
    LOAD_BLOCKS_CHECK(memcmp(a->atom.type_idx, b->atom.type_idx, n * sizeof(md_particle_type_idx_t)) == 0, "atom.type_idx");
    LOAD_BLOCKS_CHECK(memcmp(a->atom.flags, b->atom.flags, n * sizeof(md_atom_flags_t)) == 0, "atom.flags");
    for (size_t t = 0; t < a->atom.type.count; ++t) {
        LOAD_BLOCKS_CHECK(md_str_eq(md_particle_type_name(&a->atom.type, t), md_particle_type_name(&b->atom.type, t)), "type name");
        LOAD_BLOCKS_CHECK(a->atom.type.z[t] == b->atom.type.z[t], "type element");
    }
    LOAD_BLOCKS_CHECK(memcmp(a->residue.atom_offset, b->residue.atom_offset, (r + 1) * sizeof(uint32_t)) == 0, "residue.atom_offset");
    LOAD_BLOCKS_CHECK(memcmp(a->residue.name, b->residue.name, r * sizeof(md_label_t)) == 0, "residue.name");
    LOAD_BLOCKS_CHECK(memcmp(a->residue.seq_id, b->residue.seq_id, r * sizeof(md_sequence_id_t)) == 0, "residue.seq_id");
    LOAD_BLOCKS_CHECK(memcmp(a->residue.flags, b->residue.flags, r * sizeof(md_residue_flags_t)) == 0, "residue.flags");
    LOAD_BLOCKS_CHECK(a->bond.count == b->bond.count, "bond count");
    if (a->bond.count == b->bond.count) d += memcmp(a->bond.pairs, b->bond.pairs, a->bond.count * sizeof(md_atom_pair_t)) != 0;
    LOAD_BLOCKS_CHECK(a->instance.count == b->instance.count, "instance count");
    if (a->instance.count == b->instance.count && a->instance.count) {
        LOAD_BLOCKS_CHECK(memcmp(a->instance.res_offset, b->instance.res_offset, (a->instance.count + 1) * sizeof(uint32_t)) == 0, "instance.res_offset");
    }
    LOAD_BLOCKS_CHECK(as->unitcell.x == bs->unitcell.x && as->unitcell.y == bs->unitcell.y && as->unitcell.z == bs->unitcell.z &&
                      as->unitcell.xy == bs->unitcell.xy && as->unitcell.xz == bs->unitcell.xz && as->unitcell.yz == bs->unitcell.yz &&
                      as->unitcell.flags == bs->unitcell.flags, "unitcell");
    return d;
}

static bool load_blocks_write(const char* path, md_str_t text) {
    md_file_t f = {0};
    if (!md_file_open(&f, md_str_from_cstr(path), MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
    const bool ok = md_file_write(f, text.ptr, text.len) == text.len;
    md_file_close(&f);
    return ok;
}
