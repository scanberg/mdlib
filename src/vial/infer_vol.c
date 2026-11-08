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
static bool structures_alike(infer_t* x, vial_node_id n, const vial_sets_t* s, const char* fname) {
    const md_system_t* sys = x->c->sys;
    const uint32_t n0 = vial_sets_size(s, 0);
    const uint32_t* a0 = vial_sets_items(s, 0);
    const char* how = fname[0] == 's' ? "fitted" : "matched";
    for (uint32_t i = 1; i < s->count; ++i) {
        if (vial_sets_size(s, i) != n0) {
            vx_diag_fix(x, "E0320", n, NULL, "the structures of `%s` are %s particle by particle; structure %u has %u particles and the first %u",
                fname, how, i, vial_sets_size(s, i), n0);
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
            vx_diag_fix(x, "E0320", n, NULL, "the structures of `%s` are %s particle by particle; particle %u of structure %u is %.*s, and of the first %.*s",
                fname, how, k, i, (int)np.len, np.ptr, (int)nq.len, nq.ptr);
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

// An sdf's structures: a set, an array of sets alike, or set[k][w] whose columns are alike, static
static bool structures_ok(infer_t* x, vial_node_id at, const T* t, const vial_value_t* v) {
    if (vial_type_is_error(t)) return false;
    const bool ordered = t->dtype == VIAL_DT_SET && t->ndim == 2 && (t->flags & VIAL_TYPE_ONE_PARTICLE);
    if (t->dtype != VIAL_DT_SET || (t->ndim > 1 && !ordered)) {
        vx_diag_fix(x, "E0313", at, NULL, "the structures of `sdf` are a set, an array of alike sets, or set[k][w] of one particle each, as mcs and "
            "smiles_atoms give them; this is %s", vx_tstr(t));
        return false;
    }
    if (t->stage != VIAL_STAGE_STATIC) {
        vx_diag_fix(x, "E0313", at, NULL, "the structures of `sdf` are fitted particle by particle, so they are static; these change per frame");
        return false;
    }
    if (v->kind != VIAL_VAL_SETS || !v->sets || !v->sets->count) return true;
    const vial_sets_t* s = v->sets;
    if (!ordered) return structures_alike(x, at, s, "sdf");
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
                                   {"sigma", VX_KIND_WHOLE, false}};
    vial_node_id a[7] = {0};
    if (!vx_bind(x, call, fname, p, 7, a)) return err();
    bool ok = true;
    // The structures: static, so the particles fitted are the same in every frame
    const vial_value_t sv = vx_value_of(x, a[0]);
    if (!structures_ok(x, a[0], type_of(x, a[0]), &sv)) ok = false;
    if (!ref_ok(x, a[4])) ok = false;
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
    // mcs(*structures, match='element', least=0.5): every positional argument is a set or an array of sets, static
    static const vx_param_t sig[] = {{"structures", VX_KIND_WHOLE, true}, {"match", VX_KIND_WHOLE, false}, {"least", VX_KIND_WHOLE, false}};
    vx_record_sig(x, call, sig, 3);
    bool ok = true;
    const vial_sets_t* parts[64];
    uint32_t np = 0;
    vial_node_id only = 0, mn = 0, ln = 0;
    for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
        const vial_node_id a = vx_child(x, call, i);
        if (!a) continue;
        if (vx_kind(x, a) == VIAL_N_KWARG) {
            const vial_node_id kn = vx_child(x, a, 0), vn = vx_child(x, a, 1);
            if (vn) vx_infer(x, vn);
            const str_t kw = kn ? vx_text(x, kn) : STR_LIT("");
            if (str_eq(kw, STR_LIT("match"))) mn = vn;
            else if (str_eq(kw, STR_LIT("least"))) ln = vn;
            else {
                char fix[160] = "";
                if (str_eq(kw, STR_LIT("by")) && vn) {
                    const str_t vt = vx_text(x, vn);
                    if (vt.len < 140) snprintf(fix, sizeof(fix), "match=%.*s", (int)vt.len, vt.ptr);
                }
                vx_diag_fix(x, "E0202", a, fix[0] ? fix : NULL, "`mcs` takes no keyword `%.*s`; it takes structures, match and least", (int)kw.len, kw.ptr);
                ok = false;
            }
            continue;
        }
        const T t = vx_infer(x, a);
        if (vial_type_is_error(&t)) { ok = false; continue; }
        const vial_value_t v = vx_value_of(x, a);
        if (t.dtype != VIAL_DT_SET || t.ndim > 1 || t.stage != VIAL_STAGE_STATIC || v.kind != VIAL_VAL_SETS || !v.sets || !v.sets->count) {
            vx_diag_fix(x, "E0313", a, NULL, "`mcs` takes structures known when compiling, a set or an array of sets each; this is %s", vx_tstr(&t));
            ok = false;
            continue;
        }
        if (np < 64) parts[np++] = v.sets;
        only = np == 1 ? a : 0;
    }
    if (!np) {
        if (ok) vx_diag_fix(x, "E0204", call, NULL, "`mcs` needs a structure");
        return err();
    }
    bool by_name = false;
    if (mn) {
        const vial_value_t v = vx_value_of(x, mn);
        if (v.kind != VIAL_VAL_STR || !(str_eq(v.s, STR_LIT("element")) || str_eq(v.s, STR_LIT("name")))) {
            vx_diag_fix(x, "E0313", mn, NULL, "match= takes 'element' or 'name'");
            ok = false;
        }
        by_name = v.kind == VIAL_VAL_STR && str_eq(v.s, STR_LIT("name"));
    }
    double least = 0.5;
    if (ln && (!known_number(x, ln, &least) || !(least >= 0.0 && least <= 1.0))) {
        vx_diag_fix(x, "E0313", ln, NULL, "`least` is the share of the smaller structure's heavy atoms that must be common, from 0 to 1");
        ok = false;
    }
    if (!ok || !x->c->sys) return err();
    const vial_sets_t* st = np == 1 ? parts[0] : vial_sets_concat(parts, np, x->alloc);
    vial_mcs_t m;
    if (!vial_mcs(&m, x->c->sys, st, by_name, NULL, x->alloc)) {
        vx_diag_fix(x, "E0321", call, NULL, "a structure of more than %d heavy atoms is matched residue by residue, and these are not polymers", MD_MCS_GRAPH_MAX);
        return err();
    }
    if (m.least_at) {
        // The structure that shares the least with the first must still share enough
        const size_t need = ceil(least * (double)m.least_heavy) < 3 ? 3 : (size_t)ceil(least * (double)m.least_heavy);
        if (m.least_common < need) {
            vx_diag_fix(x, "E0321", call, NULL, "structure %u has %zu heavy atoms in common with the first, of %zu and %zu; least=%g asks for %zu",
                m.least_at, m.least_common, m.first_heavy, m.least_heavy, least, need);
            return err();
        }
    }
    if (m.width < 3) {
        vx_diag_fix(x, "E0321", call, NULL, "the structures have %u heavy atoms in common, all of them; fewer than 3 is nothing to place", m.width);
        return err();
    }
    if (m.incomplete) vx_diag_fix(x, "W0408", call, NULL, "the search gave up before it was sure: the atoms in common are the most it found");
    // set[k][w]: the structures, then the w particles matched, one each. One array's own dimension keeps its labels
    T t = vial_type_scalar(VIAL_DT_SET);
    t.flags |= VIAL_TYPE_ONE_PARTICLE;
    const T* ot = only ? type_of(x, only) : NULL;
    if (ot && ot->ndim == 1) vial_type_push_dim(&t, ot->dims[0]);
    else {
        vial_extent_t e = vial_extent_static(m.k);
        e.desc = vx_desc(x, VIAL_LEVEL_NONE, st);
        vial_type_push_dim(&t, e);
    }
    vial_type_push_dim(&t, vial_extent_static(m.width));
    return t;
}

bool vx_call_volume(infer_t* x, vial_node_id call, const vial_builtin_t* b, T* out) {
    if (!strcmp(b->name, "mcs")) { *out = call_mcs(x, call, b->name); return true; }
    if (!strcmp(b->name, "sdf")) { *out = call_sdf(x, call, b->name); return true; }
    if (!strcmp(b->name, "density")) { *out = call_density(x, call, b->name); return true; }
    return false;
}
