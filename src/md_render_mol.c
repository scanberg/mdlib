#include "md_render_mol.h"

#include <md_system.h>
#include <md_bond.h>
#include <md_system_state.h>
#include <core/md_log.h>

#include <string.h>

static uint32_t bond_class(md_bond_flags_t flags) {
    switch (md_bond_kind(flags)) {
    case MD_BOND_KIND_COORDINATE:  return MD_RENDER_MOL_BOND_COORDINATE;
    case MD_BOND_KIND_ELASTIC:     return MD_RENDER_MOL_BOND_ELASTIC;
    case MD_BOND_KIND_METALLIC:    return MD_RENDER_MOL_BOND_METALLIC;
    case MD_BOND_KIND_INTERACTION: return MD_RENDER_MOL_BOND_INTERACTION;
    case MD_BOND_KIND_HYDROGEN:    return MD_RENDER_MOL_BOND_HYDROGEN;
    default:                       return MD_RENDER_MOL_BOND_COVALENT;
    }
}

static uint64_t align256(uint64_t x) { return (x + 255u) & ~(uint64_t)255u; }

bool md_render_mol_plan(md_render_mol_plan_t* out, const md_system_t* sys, md_temp_scope_t temp) {
    if (!out || !sys) return false;
    memset(out, 0, sizeof(*out));
    md_render_mol_topo_t* L = &out->layout;
    L->topology_version = sys->topology_version;
    L->atom_count       = (uint32_t)sys->atom.count;

    const size_t n = sys->bond.count + sys->bond.noncovalent_count;
    uint32_t* pairs = NULL;
    uint32_t* ids   = NULL;
    if (n) {
        pairs = md_temp_alloc(temp, n * 2 * sizeof(uint32_t));
        ids   = md_temp_alloc(temp, n * sizeof(uint32_t));
        if (!pairs || !ids) return false;
        /* Counting sort by class, stable: within a class the system's order */
        uint32_t count[MD_RENDER_MOL_BOND_CLASS_COUNT] = {0};
        for (size_t i = 0; i < n; ++i) count[bond_class(sys->bond.flags ? sys->bond.flags[i] : 0)]++;
        uint32_t off = 0;
        for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) {
            L->bond_class_offset[c] = off;
            off += count[c];
        }
        L->bond_class_offset[MD_RENDER_MOL_BOND_CLASS_COUNT] = off;
        uint32_t fill[MD_RENDER_MOL_BOND_CLASS_COUNT];
        memcpy(fill, L->bond_class_offset, sizeof(fill));
        for (size_t i = 0; i < n; ++i) {
            uint32_t c = bond_class(sys->bond.flags ? sys->bond.flags[i] : 0);
            uint32_t j = fill[c]++;
            pairs[2 * j + 0] = (uint32_t)sys->bond.pairs[i].idx[0];
            pairs[2 * j + 1] = (uint32_t)sys->bond.pairs[i].idx[1];
            ids[j]           = (uint32_t)i;
        }
    }
    L->bond_count = (uint32_t)n;
    L->bond_pairs = 0;
    L->bond_ids   = align256((uint64_t)n * 8u);
    out->size     = align256(L->bond_ids + (uint64_t)n * 4u);
    if (out->size == 0) out->size = 256;
    L->size       = out->size;
    out->bond_pairs = pairs;
    out->bond_ids   = ids;
    return true;
}

bool md_render_mol_write(md_gpu_stream_t s, md_render_mol_topo_t* out, const md_render_mol_plan_t* plan,
                         md_gpu_addr_t base, uint64_t content[MD_RENDER_CONTENT_COUNT]) {
    if (!s || !out || !plan || !base) return false;
    *out = plan->layout;
    out->base       = base;
    out->bond_pairs = base + plan->layout.bond_pairs;
    out->bond_ids   = base + plan->layout.bond_ids;
    if (plan->layout.bond_count) {
        if (!md_gpu_upload(s, out->bond_pairs, plan->bond_pairs, (size_t)plan->layout.bond_count * 8u)) return false;
        if (!md_gpu_upload(s, out->bond_ids,   plan->bond_ids,   (size_t)plan->layout.bond_count * 4u)) return false;
    }
    if (content) content[MD_RENDER_CONTENT_GEOMETRY]++;
    return true;
}

bool md_render_mol_upload_state(md_gpu_stream_t s, md_gpu_addr_t dst, const md_system_state_t* state,
                                uint64_t content[MD_RENDER_CONTENT_COUNT]) {
    if (!s || !dst || !state) return false;
    _Static_assert(sizeof(md_float3_t) == 12, "packed float3");
    if (state->num_atoms && !md_gpu_upload(s, dst, state->xyz, state->num_atoms * 12u)) return false;
    if (content) content[MD_RENDER_CONTENT_GEOMETRY]++;
    return true;
}

/* ---- Emit ------------------------------------------------------------------------ */

static md_render_points_t atom_points(const md_render_mol_topo_t* topo, const md_render_mol_atoms_t* atoms, const md_render_mol_rep_t* rep) {
    md_render_points_t p;
    memset(&p, 0, sizeof(p));
    p.count           = topo->atom_count;
    p.position        = atoms->position;
    p.stride          = atoms->stride;
    p.velocity        = atoms->velocity;
    p.velocity_stride = atoms->velocity_stride;
    p.radius          = atoms->radius;
    p.color           = rep->atom_color;
    p.palette_index   = rep->atom_palette_index;
    p.visible[0]      = rep->visible;
    p.visible[1]      = atoms->visible;
    p.state[0]        = atoms->selection;
    p.state[1]        = atoms->highlight;
    return p;
}

static void put(md_render_batch_t* out, uint32_t cap, uint32_t* n, const md_render_batch_t* b) {
    if (*n < cap) out[*n] = *b;
    (*n)++;
}

uint32_t md_render_mol_emit(const md_render_mol_topo_t* topo, const md_render_mol_atoms_t* atoms,
                            const md_render_mol_rep_t* rep, md_render_batch_t* out, uint32_t cap) {
    if (!topo || !atoms || !rep || (cap && !out)) return 0;
    if (rep->key == 0 || rep->key >= (1ull << 56)) {
        MD_LOG_ERROR("md_render_mol: representation key must be in [1, 2^56)");
        return 0;
    }
    if (!atoms->position) {
        MD_LOG_ERROR("md_render_mol: no atom positions");
        return 0;
    }

    md_render_batch_t base;
    memset(&base, 0, sizeof(base));
    base.points          = atom_points(topo, atoms, rep);
    base.color           = rep->color;
    base.transform       = rep->transform;
    base.transform_count = rep->transform_count;

    bool  spheres = false, sticks = false;
    float sphere_scale = 0.0f, sphere_add = 0.0f, stick_radius = 0.0f;
    uint32_t classes = rep->bond_classes;
    const uint32_t default_classes = MD_RENDER_MOL_BONDS(MD_RENDER_MOL_BOND_COVALENT) | MD_RENDER_MOL_BONDS(MD_RENDER_MOL_BOND_COORDINATE);

    switch (rep->type) {
    case MD_RENDER_MOL_REP_SPACEFILL:
        spheres = true;
        sphere_scale = rep->params.spacefill.radius_scale;
        break;
    case MD_RENDER_MOL_REP_LICORICE:
        spheres = sticks = true;
        sphere_add = stick_radius = rep->params.licorice.radius;
        if (!classes) classes = default_classes;
        break;
    case MD_RENDER_MOL_REP_BALL_AND_STICK:
        spheres = sticks = true;
        sphere_scale = rep->params.ball_and_stick.ball_scale;
        stick_radius = rep->params.ball_and_stick.stick_radius;
        if (!classes) classes = default_classes;
        break;
    case MD_RENDER_MOL_REP_BONDS:
        sticks = true;
        stick_radius = rep->params.bonds.radius;
        if (!classes) classes = MD_RENDER_MOL_BONDS_ALL;
        break;
    default:
        MD_LOG_ERROR("md_render_mol: representation type %d is not implemented yet", (int)rep->type);
        return 0;
    }

    /* Sticks first, atoms last: where a stick's round end and an atom's sphere
       are the same surface (licorice joints) the later draw wins the depth
       tie, and a joint should pick as its atom. */
    uint32_t n = 0;
    if (sticks) {
        for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) {
            if (!(classes & MD_RENDER_MOL_BONDS(c))) continue;
            uint32_t beg = topo->bond_class_offset[c], end = topo->bond_class_offset[c + 1];
            if (end <= beg) continue;
            md_render_batch_t b = base;
            b.key   = (rep->key << 8) | (MD_RENDER_MOL_PART_BONDS + c);
            b.kind  = MD_RENDER_CAPSULES;
            b.count = end - beg;
            b.index = topo->bond_pairs + (md_gpu_addr_t)beg * 8u;
            b.id    = topo->bond_ids   + (md_gpu_addr_t)beg * 4u;
            b.style.capsules.radius = stick_radius;
            put(out, cap, &n, &b);
        }
    }
    if (spheres && topo->atom_count) {
        if (sphere_scale > 0.0f && !atoms->radius) {
            MD_LOG_ERROR("md_render_mol: radius-scaled spheres need atom radii");
            return 0;
        }
        md_render_batch_t b = base;
        b.key   = (rep->key << 8) | MD_RENDER_MOL_PART_ATOMS;
        b.kind  = MD_RENDER_SPHERES;
        b.count = topo->atom_count;
        b.style.spheres.radius_scale = sphere_scale;
        b.style.spheres.radius_add   = sphere_add;
        put(out, cap, &n, &b);
    }
    return n;
}

bool md_render_mol_decode_pick(const md_render_pick_t* pick, md_render_mol_hit_t* hit) {
    if (!pick || !hit) return false;
    memset(hit, 0, sizeof(*hit));
    hit->bond = UINT32_MAX;
    if (pick->key == 0) return false;
    uint32_t part = (uint32_t)(pick->key & 0xFFu);
    hit->rep_key  = pick->key >> 8;
    hit->instance = pick->instance;
    memcpy(hit->position, pick->position, sizeof(hit->position));
    if (part == MD_RENDER_MOL_PART_ATOMS && pick->kind == MD_RENDER_SPHERES) {
        hit->kind = MD_RENDER_MOL_HIT_ATOM;
        hit->atom = pick->attr;
        return true;
    }
    if (part >= MD_RENDER_MOL_PART_BONDS && part < MD_RENDER_MOL_PART_BONDS + MD_RENDER_MOL_BOND_CLASS_COUNT && pick->kind == MD_RENDER_CAPSULES) {
        hit->kind = MD_RENDER_MOL_HIT_BOND;
        hit->bond = pick->id;
        hit->atom = pick->attr;
        return true;
    }
    return false;
}
