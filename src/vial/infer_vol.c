// The volumes, sdf() and density() (volume.h): their arguments checked, and their grids worked out when compiling, so
// a volume's type has its extents and what is too large is said before anything is evaluated.

#include "infer_internal.h"
#include "volume.h"
#include "builtins.h"

#include <md_system.h>
#include <md_mcs.h>
#include <core/md_unit.h>
#include <core/md_common.h>

#include <math.h>
#include <string.h>

typedef vial_type_t T;

static T err(void) { return vial_type_error(); }

static const T* type_of(const infer_t* x, vial_node_id n) { return &x->out->types[n - x->out->base]; }

static md_unit_t angstrom(void) {
    md_unit_t u;
    vial_unit_name(STR_LIT("angstrom"), &u);
    return u;
}

// A number known when compiling, as a length in Å: a bare number, or one with a unit of length
static bool known_length(infer_t* x, vial_node_id n, double* out) {
    if (!n) return false;
    const vial_value_t v = vx_value_of(x, n);
    double f;
    if (v.kind == VIAL_VAL_INT) f = (double)v.i;
    else if (v.kind == VIAL_VAL_FLOAT) f = v.f;
    else return false;
    const T* t = type_of(x, n);
    if (t->unit_kind == VIAL_UNIT_KNOWN && !md_unit_is_none(t->unit)) {
        double c;
        if (!md_unit_conversion_factor(&c, t->unit, angstrom())) return false;
        f *= c;
    }
    *out = f;
    return true;
}

static bool known_number(infer_t* x, vial_node_id n, double* out) {
    if (!n) return false;
    const vial_value_t v = vx_value_of(x, n);
    if (v.kind == VIAL_VAL_INT) { *out = (double)v.i; return true; }
    if (v.kind == VIAL_VAL_FLOAT) { *out = v.f; return true; }
    return false;
}

static bool is_none(infer_t* x, vial_node_id n) {
    return !n || vx_value_of(x, n).kind == VIAL_VAL_NULL;
}

// Three lengths: one for all three, or a tuple of three
static bool three_lengths(infer_t* x, vial_node_id n, double out[3]) {
    if (vx_kind(x, n) == VIAL_N_TUPLE) {
        if (vx_num_children(x, n) != 3) return false;
        for (uint32_t a = 0; a < 3; ++a) if (!known_length(x, vx_child(x, n, a), &out[a])) return false;
        return true;
    }
    if (!known_length(x, n, &out[0])) return false;
    out[1] = out[2] = out[0];
    return true;
}

// Points, float[n][3] in a unit of length, or one set: what a volume adds up
static bool target_ok(infer_t* x, vial_node_id n, const char* fname, const char* pname) {
    const T* t = type_of(x, n);
    if (vial_type_is_error(t)) return false;
    if (t->dtype == VIAL_DT_SET) {
        if (t->ndim == 0) return true;
        vx_diag_fix(x, "E0313", n, NULL, "%s= of `%s` is one set or points; this is %s, an array of sets", pname, fname, vx_tstr(t));
        return false;
    }
    const bool fl = t->dtype == VIAL_DT_FLOAT || t->dtype == VIAL_DT_FLOAT64;
    const bool pts = fl && t->ndim == 2 && t->dims[1].kind == VIAL_EXT_STATIC && t->dims[1].n == 3;
    const bool len = t->unit_kind != VIAL_UNIT_KNOWN || md_unit_base_equal(t->unit, angstrom());
    if (pts && len) return true;
    vx_diag_fix(x, "E0313", n, NULL, "%s= of `%s` is one set or points, float[n][3] Å; this is %s", pname, fname, vx_tstr(t));
    return false;
}

// What one point weighs, as weights= is everywhere: 'mass', None, or one number per point, float[n]; the unit of
// the volume is the weights' over Å³
static bool weights_ok(infer_t* x, vial_node_id n, const char* fname, md_unit_t* unit) {
    *unit = md_unit_count();
    if (is_none(x, n)) return true;
    if (!vx_weights_ok(x, n, fname, NULL)) return false;
    const T* t = type_of(x, n);
    const vial_value_t w = vx_value_of(x, n);
    if (w.kind == VIAL_VAL_STR) { vial_unit_name(STR_LIT("u"), unit); return true; }
    if (t->unit_kind == VIAL_UNIT_KNOWN && !md_unit_is_none(t->unit)) *unit = t->unit;
    return true;
}

static bool samples_of(infer_t* x, vial_node_id n, double* out) {
    *out = VIAL_VOLUME_SAMPLES;
    if (!n) return true;
    double s;
    if (!known_number(x, n, &s) || !(s > 0.0)) {
        vx_diag_fix(x, "E0313", n, NULL, "samples= is a number of samples per Å, more than 0, known when compiling");
        return false;
    }
    *out = s;
    return true;
}

static bool sigma_ok(infer_t* x, vial_node_id n) {
    double s;
    if (!n || is_none(x, n)) return true;
    if (known_length(x, n, &s) && s >= 0.0) return true;
    vx_diag_fix(x, "E0313", n, NULL, "sigma= is a length, 0 or more, known when compiling");
    return false;
}

// The grid of lengths len at samples per Å, W0407 when it is sampled coarser to fit
static void grid(infer_t* x, vial_node_id call, const double len[3], double samples, uint32_t dims[3]) {
    bool coarser = false;
    const double used = vial_volume_grid(len, samples, dims, &coarser);
    if (!coarser) return;
    double want = 1.0;
    for (int a = 0; a < 3; ++a) want *= ceil(len[a] * samples);
    vx_diag_fix(x, "W0407", call, NULL, "%.0f voxels at %.3g samples per Å are more than 512³; it is sampled at %.3g per Å, %u×%u×%u",
        want, samples, used, dims[0], dims[1], dims[2]);
}

// The structures of an sdf: alike, particle by particle, as v1 has them: as many particles, of the same elements, or
// of the same names where an element is not known
static bool structures_alike(infer_t* x, vial_node_id n, const vial_sets_t* s) {
    const md_system_t* sys = x->c->sys;
    const uint32_t n0 = vial_sets_size(s, 0);
    const uint32_t* a0 = vial_sets_items(s, 0);
    for (uint32_t i = 1; i < s->count; ++i) {
        if (vial_sets_size(s, i) != n0) {
            vx_diag_fix(x, "E0320", n, NULL, "the structures of `sdf` are fitted particle by particle; structure %u has %u particles and the first %u",
                i, vial_sets_size(s, i), n0);
            return false;
        }
        if (!sys) continue;
        const uint32_t* ai = vial_sets_items(s, i);
        for (uint32_t k = 0; k < n0; ++k) {
            const uint32_t p = ai[k], q = a0[k];
            if (p >= sys->atom.count || q >= sys->atom.count) continue;
            const md_atomic_number_t zp = md_atom_atomic_number(&sys->atom, p), zq = md_atom_atomic_number(&sys->atom, q);
            const bool same = (zp && zq) ? zp == zq : str_eq(md_atom_name(&sys->atom, p), md_atom_name(&sys->atom, q));
            if (same) continue;
            const str_t np = md_atom_name(&sys->atom, p), nq = md_atom_name(&sys->atom, q);
            vx_diag_fix(x, "E0320", n, NULL, "the structures of `sdf` are fitted particle by particle; particle %u of structure %u is %.*s, and of the first %.*s",
                k, i, (int)np.len, np.ptr, (int)nq.len, nq.ptr);
            return false;
        }
    }
    return true;
}

// Two particles alike as structures pair them: of one element, or of one name where an element is not known
static bool alike(const md_system_t* sys, uint32_t p, uint32_t q) {
    if (!sys || p >= sys->atom.count || q >= sys->atom.count) return true;
    const md_atomic_number_t zp = md_atom_atomic_number(&sys->atom, p), zq = md_atom_atomic_number(&sys->atom, q);
    return (zp && zq) ? zp == zq : str_eq(md_atom_name(&sys->atom, p), md_atom_name(&sys->atom, q));
}

// One kind of an sdf's structures: a set, an array of sets alike, or set[k][w] whose columns are alike, static; its
// first structure's particles as a set into *first, when they are known
static bool kind_ok(infer_t* x, vial_node_id at, const T* t, const vial_value_t* v, const vial_sets_t** first) {
    *first = NULL;
    if (vial_type_is_error(t)) return false;
    const bool ordered = t->dtype == VIAL_DT_SET && t->ndim == 2 && (t->flags & VIAL_TYPE_ONE_PARTICLE);
    if (t->dtype != VIAL_DT_SET || (t->ndim > 1 && !ordered)) {
        vx_diag_fix(x, "E0313", at, NULL, "the structures of `sdf` are a set, an array of sets, set[k][w] of one particle each, or a tuple of kinds of them; this is %s", vx_tstr(t));
        return false;
    }
    if (t->stage != VIAL_STAGE_STATIC) {
        vx_diag_fix(x, "E0313", at, NULL, "the structures of `sdf` are fitted particle by particle, so they are static; these change per frame");
        return false;
    }
    if (v->kind != VIAL_VAL_SETS || !v->sets || !v->sets->count) return true;
    const vial_sets_t* s = v->sets;
    if (!ordered) {
        if (!structures_alike(x, at, s)) return false;
        *first = vial_sets_one(vial_sets_items(s, 0), vial_sets_size(s, 0), x->alloc);
        return true;
    }
    const uint32_t w = t->dims[1].kind == VIAL_EXT_STATIC ? t->dims[1].n : 0;
    if (!w || s->count % w) return true;
    const md_system_t* sys = x->c->sys;
    for (uint32_t i = 1; i < s->count / w; ++i) {
        for (uint32_t j = 0; j < w; ++j) {
            const uint32_t p = vial_sets_items(s, i * w + j)[0], q = vial_sets_items(s, j)[0];
            if (alike(sys, p, q)) continue;
            const str_t np = md_atom_name(&sys->atom, p), nq = md_atom_name(&sys->atom, q);
            vx_diag_fix(x, "E0320", at, NULL, "the structures of `sdf` are fitted column by column; column %u of structure %u is %.*s, and of the first %.*s",
                j, i, (int)np.len, np.ptr, (int)nq.len, nq.ptr);
            return false;
        }
    }
    uint32_t* row = md_alloc(x->alloc, (w + 1) * sizeof(uint32_t));
    for (uint32_t j = 0; j < w; ++j) row[j] = vial_sets_items(s, j)[0];
    for (uint32_t i = 1; i < w; ++i) { const uint32_t r = row[i]; uint32_t j = i; while (j > 0 && row[j - 1] > r) { row[j] = row[j - 1]; --j; } row[j] = r; }
    *first = vial_sets_one(row, w, x->alloc);
    return true;
}

static bool ref_ok(infer_t* x, vial_node_id at) {
    if (!at || is_none(x, at)) return true;
    const T* rt = type_of(x, at);
    if (vial_type_is_error(rt)) return false;
    if (rt->stage != VIAL_STAGE_STATIC) {
        vx_diag_fix(x, "E0313", at, NULL, "ref= is the pose the structures are fitted onto, so it is static: at(..., frame=k)");
        return false;
    }
    if (!((rt->dtype == VIAL_DT_FLOAT || rt->dtype == VIAL_DT_FLOAT64) && rt->ndim == 2 && rt->dims[1].kind == VIAL_EXT_STATIC && rt->dims[1].n == 3)) {
        vx_diag_fix(x, "E0313", at, NULL, "ref= is the positions of one structure, float[n][3]; this is %s", vx_tstr(rt));
        return false;
    }
    return true;
}

// A pose of a kind: None, or (R, t), float[3][3] and float[3], static
static bool pose_ok(infer_t* x, vial_node_id at) {
    if (!at || is_none(x, at)) return true;
    const T* t = type_of(x, at);
    if (vial_type_is_error(t)) return false;
    const bool ok = t->dtype == VIAL_DT_TUPLE && t->num_elems == 2 && t->stage == VIAL_STAGE_STATIC &&
                    t->elems[0].ndim == 2 && t->elems[0].dims[0].kind == VIAL_EXT_STATIC && t->elems[0].dims[0].n == 3 &&
                    t->elems[0].dims[1].kind == VIAL_EXT_STATIC && t->elems[0].dims[1].n == 3 &&
                    t->elems[1].ndim == 1 && t->elems[1].dims[0].kind == VIAL_EXT_STATIC && t->elems[1].dims[0].n == 3;
    if (!ok) vx_diag_fix(x, "E0313", at, NULL, "a pose is (R, t), float[3][3] and float[3] Å, static, that places a kind's reference on the first's: x = R @ y + t; this is %s", vx_tstr(t));
    return ok;
}

static T volume_type(const uint32_t dims[3], bool known, vial_node_id call, md_unit_t unit, const T* from) {
    T t = vial_type_scalar(VIAL_DT_FLOAT);
    for (int a = 0; a < 3; ++a) {
        vial_type_push_dim(&t, known ? vial_extent_static(dims[a]) : (vial_extent_t){.kind = VIAL_EXT_UNKNOWN, .id = (uint32_t)call + (uint32_t)a});
    }
    t.unit_kind = VIAL_UNIT_KNOWN;
    t.unit = md_unit_div(unit, md_unit_pow(angstrom(), 3));
    t.stage = VIAL_STAGE_TRAJECTORY;
    t.pass = vx_pass(from);
    return t;
}

static T call_sdf(infer_t* x, vial_node_id call, const char* fname) {
    static const vx_param_t p[] = {{"structures", VX_KIND_WHOLE, true}, {"target", VX_KIND_WHOLE, true}, {"extent", VX_KIND_WHOLE, true},
                                   {"samples", VX_KIND_WHOLE, false}, {"ref", VX_KIND_WHOLE, false}, {"weights", VX_KIND_WHOLE, false},
                                   {"sigma", VX_KIND_WHOLE, false}, {"poses", VX_KIND_WHOLE, false}};
    vial_node_id a[8] = {0};
    if (!vx_bind(x, call, fname, p, 8, a)) return err();
    bool ok = true;
    // The structures: static, so the particles fitted are the same in every frame; one kind, or a tuple of kinds
    const T* st = type_of(x, a[0]);
    enum { MAX_KINDS = 64 };
    const vial_sets_t* firsts[MAX_KINDS] = {0};
    uint32_t nk = 1;
    const bool literal = vx_kind(x, a[0]) == VIAL_N_TUPLE;
    if (vial_type_is_error(st)) {
        ok = false;
    } else if (st->dtype == VIAL_DT_TUPLE) {
        nk = st->num_elems;
        if (!nk || nk > MAX_KINDS) {
            vx_diag_fix(x, "E0313", a[0], NULL, "the kinds of structure of `sdf` are 1 to %d", MAX_KINDS);
            ok = false;
        }
        for (uint32_t i = 0; ok && i < nk; ++i) {
            const vial_node_id c = literal ? vx_child(x, a[0], i) : 0;
            const vial_value_t v = c ? vx_value_of(x, c) : (vial_value_t){0};
            if (!kind_ok(x, c ? c : a[0], &st->elems[i], &v, &firsts[i])) ok = false;
        }
    } else {
        const vial_value_t v = vx_value_of(x, a[0]);
        if (!kind_ok(x, a[0], st, &v, &firsts[0])) ok = false;
    }
    // A reference per kind, and a pose per kind but the first
    if (a[4] && vx_kind(x, a[4]) == VIAL_N_TUPLE && type_of(x, a[4])->dtype == VIAL_DT_TUPLE) {
        if (vx_num_children(x, a[4]) != nk) { vx_diag_fix(x, "E0313", a[4], NULL, "ref= gives one reference per kind, %u", nk); ok = false; }
        for (uint32_t i = 0; ok && i < nk; ++i) if (!ref_ok(x, vx_child(x, a[4], i))) ok = false;
    } else if (!ref_ok(x, a[4])) {
        ok = false;
    }
    bool posed[MAX_KINDS] = {0};
    if (a[7] && !is_none(x, a[7])) {
        if (vx_kind(x, a[7]) != VIAL_N_TUPLE || vx_num_children(x, a[7]) != nk) {
            vx_diag_fix(x, "E0313", a[7], NULL, "poses= is a tuple of a pose per kind, %u, the first None", nk);
            ok = false;
        } else {
            for (uint32_t i = 0; i < nk; ++i) {
                const vial_node_id c = vx_child(x, a[7], i);
                if (i == 0 && !is_none(x, c)) { vx_diag_fix(x, "E0313", c, NULL, "the first kind is where the others are placed: its pose is None"); ok = false; }
                else if (!pose_ok(x, c)) ok = false;
                posed[i] = !is_none(x, c);
            }
        }
    }
    // What each kind has in common with the first: enough to place it, or a pose given
    for (uint32_t i = 1; ok && i < nk; ++i) {
        if (posed[i] || !firsts[0] || !firsts[i] || !x->c->sys) continue;
        vial_mcs_t m;
        const vial_node_id at = literal ? vx_child(x, a[0], i) : a[0];
        if (!vial_mcs(&m, x->c->sys, firsts[0], firsts[i], false, NULL, x->alloc)) {
            vx_diag_fix(x, "E0321", at, NULL, "kind %u of `sdf` is too large to find what it has in common with the first; give poses=", i);
            ok = false;
            continue;
        }
        const size_t common = m.a->count, smaller = m.a_heavy < m.b_heavy ? m.a_heavy : m.b_heavy;
        if (common < 3 || 2 * common < smaller) {
            vx_diag_fix(x, "E0321", at, NULL, "kind %u of `sdf` has %zu heavy atoms in common with the first, of %zu and %zu: too little to place it; give poses=",
                i, common, m.b_heavy, m.a_heavy);
            ok = false;
        }
    }
    if (!target_ok(x, a[1], fname, "target")) ok = false;
    double ext[3];
    if (!three_lengths(x, a[2], ext) || !(ext[0] > 0.0 && ext[1] > 0.0 && ext[2] > 0.0)) {
        vx_diag_fix(x, "E0313", a[2], NULL, "extent= is how far the volume reaches from the structure's centre: a length more than 0, or three, known when compiling");
        ok = false;
    }
    double samples;
    if (!samples_of(x, a[3], &samples)) ok = false;
    md_unit_t unit;
    if (!weights_ok(x, a[5], fname, &unit)) ok = false;
    if (!sigma_ok(x, a[6])) ok = false;
    if (!ok) return err();
    uint32_t dims[3];
    const double len[3] = {2.0 * ext[0], 2.0 * ext[1], 2.0 * ext[2]};
    grid(x, call, len, samples, dims);
    return volume_type(dims, true, call, unit, type_of(x, a[1]));
}

static T call_density(infer_t* x, vial_node_id call, const char* fname) {
    static const vx_param_t p[] = {{"points", VX_KIND_WHOLE, true}, {"samples", VX_KIND_WHOLE, false}, {"box", VX_KIND_WHOLE, false},
                                   {"weights", VX_KIND_WHOLE, false}, {"sigma", VX_KIND_WHOLE, false}};
    vial_node_id a[5] = {0};
    if (!vx_bind(x, call, fname, p, 5, a)) return err();
    bool ok = target_ok(x, a[0], fname, "points");
    double samples;
    if (!samples_of(x, a[1], &samples)) ok = false;
    double lo[3], hi[3];
    const bool box = !is_none(x, a[2]);
    if (box) {
        const vial_node_id b = a[2];
        bool bok = vx_kind(x, b) == VIAL_N_TUPLE && vx_num_children(x, b) == 2;
        if (bok) {
            const vial_node_id l = vx_child(x, b, 0), h = vx_child(x, b, 1);
            bok = vx_kind(x, l) == VIAL_N_TUPLE && vx_kind(x, h) == VIAL_N_TUPLE && three_lengths(x, l, lo) && three_lengths(x, h, hi);
        }
        for (int q = 0; bok && q < 3; ++q) bok = hi[q] > lo[q];
        if (!bok) {
            vx_diag_fix(x, "E0313", b, NULL, "box= is ((x0, y0, z0), (x1, y1, z1)), lengths known when compiling, each hi above its lo");
            ok = false;
        }
    }
    md_unit_t unit;
    if (!weights_ok(x, a[3], fname, &unit)) ok = false;
    if (!sigma_ok(x, a[4])) ok = false;
    if (!ok) return err();
    uint32_t dims[3] = {0};
    if (box) {
        const double len[3] = {hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
        grid(x, call, len, samples, dims);
    }
    // Over the cell: the grid of frame 0's cell, known when evaluating
    return volume_type(dims, box, call, unit, type_of(x, a[0]));
}

// mcs(a, b, by='element', least=0.5): known when compiling, from the topology; evaluation chooses among symmetric
// mappings by frame 0's positions
static T call_mcs(infer_t* x, vial_node_id call, const char* fname) {
    static const vx_param_t p[] = {{"a", VX_KIND_WHOLE, true}, {"b", VX_KIND_WHOLE, true}, {"match", VX_KIND_WHOLE, false}, {"least", VX_KIND_WHOLE, false}};
    vial_node_id a[4] = {0};
    if (!vx_bind(x, call, fname, p, 4, a)) return err();
    bool ok = true;
    const vial_sets_t* s[2] = {0};
    for (int i = 0; i < 2; ++i) {
        const T* t = type_of(x, a[i]);
        if (vial_type_is_error(t)) { ok = false; continue; }
        const vial_value_t v = vx_value_of(x, a[i]);
        if (t->dtype != VIAL_DT_SET || t->ndim != 0 || t->stage != VIAL_STAGE_STATIC || v.kind != VIAL_VAL_SETS || !v.sets || v.sets->count != 1) {
            vx_diag_fix(x, "E0313", a[i], NULL, "`mcs` takes two sets known when compiling, a structure each; this is %s", vx_tstr(t));
            ok = false;
            continue;
        }
        s[i] = v.sets;
    }
    bool by_name = false;
    if (a[2]) {
        const vial_value_t v = vx_value_of(x, a[2]);
        if (v.kind != VIAL_VAL_STR || !(str_eq(v.s, STR_LIT("element")) || str_eq(v.s, STR_LIT("name")))) {
            vx_diag_fix(x, "E0313", a[2], NULL, "match= takes 'element' or 'name'");
            ok = false;
        }
        by_name = v.kind == VIAL_VAL_STR && str_eq(v.s, STR_LIT("name"));
    }
    double least = 0.5;
    if (a[3] && (!known_number(x, a[3], &least) || !(least >= 0.0 && least <= 1.0))) {
        vx_diag_fix(x, "E0313", a[3], NULL, "`least` is the share of the smaller structure's heavy atoms that must be common, from 0 to 1");
        ok = false;
    }
    if (!ok || !x->c->sys) return err();
    vial_mcs_t m;
    if (!vial_mcs(&m, x->c->sys, s[0], s[1], by_name, NULL, x->alloc)) {
        vx_diag_fix(x, "E0321", call, NULL, "a structure of more than %d heavy atoms is matched residue by residue, and these are not polymers", MD_MCS_GRAPH_MAX);
        return err();
    }
    const size_t common = m.a->count, smaller = m.a_heavy < m.b_heavy ? m.a_heavy : m.b_heavy;
    if (common < 3 || (double)common < least * (double)smaller) {
        vx_diag_fix(x, "E0321", call, NULL, "the two have %zu heavy atoms in common, of %zu and %zu; least=%g asks for %.0f", common, m.a_heavy, m.b_heavy, least,
            ceil(least * (double)smaller) < 3 ? 3 : ceil(least * (double)smaller));
        return err();
    }
    if (m.incomplete) vx_diag_fix(x, "W0408", call, NULL, "the search gave up before it was sure: %zu heavy atoms in common is the most it found", common);
    T* elems = md_alloc(x->alloc, 2 * sizeof(T));
    elems[0] = vx_sets_type(x, m.a, call);
    elems[1] = vx_sets_type(x, m.b, call);
    T t = vial_type_scalar(VIAL_DT_TUPLE);
    t.elems = elems;
    t.num_elems = 2;
    return t;
}

bool vx_call_volume(infer_t* x, vial_node_id call, const vial_builtin_t* b, T* out) {
    if (!strcmp(b->name, "mcs")) { *out = call_mcs(x, call, b->name); return true; }
    if (!strcmp(b->name, "sdf")) { *out = call_sdf(x, call, b->name); return true; }
    if (!strcmp(b->name, "density")) { *out = call_density(x, call, b->name); return true; }
    return false;
}
