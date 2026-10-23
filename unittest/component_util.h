#pragma once

// Molecules as the chemistry tests read them (unittest/ccd_chem.txt, unittest/ring_chem.txt), and systems of one of
// them. Shared by test_chem.c and test_ring.c.
// Format: COMP id / A name element charge x y z / B i j order aromatic, '#' starts a comment line.

#include <md_system.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>

#include <stdio.h>
#include <string.h>

typedef struct ref_atom_t {
    char  name[8];
    int   z;
    int   q;
    vec3_t x;
} ref_atom_t;

typedef struct ref_bond_t {
    int i, j, order, aromatic;
} ref_bond_t;

typedef struct ref_comp_t {
    char id[8];
    md_array(ref_atom_t) atom;
    md_array(ref_bond_t) bond;
} ref_comp_t;

static inline md_array(ref_comp_t) read_components(str_t path, md_allocator_i* alloc) {
    md_array(ref_comp_t) comps = 0;
    str_t text = load_textfile(path, alloc);
    str_t line;
    char buf[256];
    while (str_extract_line(&line, &text)) {
        if (line.len == 0 || line.ptr[0] == '#' || line.len >= sizeof(buf)) continue;
        MEMCPY(buf, line.ptr, line.len);
        buf[line.len] = '\0';
        if (!strncmp(buf, "COMP ", 5)) {
            ref_comp_t c = { 0 };
            sscanf(buf + 5, "%7s", c.id);
            md_array_push(comps, c, alloc);
        } else if (buf[0] == 'A' && buf[1] == ' ' && md_array_size(comps)) {
            ref_atom_t a = { 0 };
            char el[8];
            sscanf(buf + 2, "%7s %7s %d %f %f %f", a.name, el, &a.q, &a.x.x, &a.x.y, &a.x.z);
            a.z = md_util_element_lookup((str_t){ el, strlen(el) }, true);
            md_array_push(md_array_last(comps)->atom, a, alloc);
        } else if (buf[0] == 'B' && buf[1] == ' ' && md_array_size(comps)) {
            ref_bond_t b = { 0 };
            sscanf(buf + 2, "%d %d %d %d", &b.i, &b.j, &b.order, &b.aromatic);
            md_array_push(md_array_last(comps)->bond, b, alloc);
        }
    }
    return comps;
}

static inline const ref_comp_t* find_component(const ref_comp_t* comps, const char* id) {
    for (size_t i = 0; i < md_array_size(comps); ++i) if (!strcmp(comps[i].id, id)) return &comps[i];
    return NULL;
}

// A system of one component; map[ref atom] = system atom or -1 (stripped hydrogens). The bonds are the CCD's, or
// inferred from the coordinates.
static inline void build_component(md_system_t* sys, md_system_state_t* st, int* map, const ref_comp_t* c, bool strip_h, bool infer) {
    md_allocator_i* alloc = sys->alloc;
    const int na = (int)md_array_size(c->atom);
    int n = 0;
    for (int i = 0; i < na; ++i) map[i] = (strip_h && c->atom[i].z == 1) ? -1 : n++;

    sys->atom.count = n;
    md_particle_type_add(&sys->atom.type, STR_LIT(""), (str_t){0}, 0, 0, 0, 0, 0, alloc);
    md_array_resize(sys->atom.type_idx, (size_t)n, alloc);
    md_array_resize(sys->atom.flags, (size_t)n, alloc);
    MEMSET(sys->atom.flags, 0, md_array_bytes(sys->atom.flags));
    md_system_state_init(st, n);
    for (int i = 0; i < na; ++i) {
        if (map[i] < 0) continue;
        const ref_atom_t* a = &c->atom[i];
        sys->atom.type_idx[map[i]] = md_particle_type_find_or_add(&sys->atom.type, str_from_cstr(a->name), (md_atomic_number_t)a->z, md_util_element_atomic_mass(a->z), md_util_element_vdw_radius(a->z), 0, 0, alloc);
        st->xyz[map[i]] = a->x;
    }
    st->num_atoms = n;
    sys->residue.count = 1;
    md_array_push(sys->residue.name, make_label(str_from_cstr(c->id)), alloc);
    md_array_push(sys->residue.seq_id, 1, alloc);
    md_array_push(sys->residue.flags, 0, alloc);
    md_array_push(sys->residue.atom_offset, 0, alloc);
    md_array_push(sys->residue.atom_offset, (uint32_t)n, alloc);

    if (infer) {
        md_util_infer_covalent_bonds(&sys->bond, st, sys, alloc);
    } else {
        for (size_t b = 0; b < md_array_size(c->bond); ++b) {
            const ref_bond_t* r = &c->bond[b];
            if (map[r->i] < 0 || map[r->j] < 0) continue;
            md_atom_pair_t p = { .idx = { map[r->i], map[r->j] } };
            md_array_push(sys->bond.pairs, p, alloc);
            md_array_push(sys->bond.flags, md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY), alloc);
            sys->bond.count += 1;
        }
        md_bond_build_connectivity(&sys->bond, n, alloc);
    }
    md_util_system_infer_rings(sys);
}

