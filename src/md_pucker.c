#include <md_pucker.h>

#include <md_system.h>
#include <md_util.h>
#include <md_types.h>

#include <core/md_common.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct pucker_canon_t {
    float theta;
    float phi;
    char  kind;
    int8_t above[2];
    int8_t below[2];
} pucker_canon_t;

// Generated from the geometry of ideal rings at each position (see md_pucker.h, CONFORMATIONS)

// Pyranose conformations (IUPAC), at their Cremer-Pople (theta, phi): ring positions 0 (the O) to 5, in the order
// of the ring from the O through the anomeric C. Above: the side from which that order runs clockwise.
static const pucker_canon_t pucker_canon6[] = {
    {   0.0f,    0.0f, 'C', {4, -1}, {1, -1}},
    { 180.0f,    0.0f, 'C', {1, -1}, {4, -1}},
    {  90.0f,    0.0f, 'B', {0, 3}, {-1, -1}},
    {  90.0f,   30.0f, 'S', {3, -1}, {1, -1}},
    {  90.0f,   60.0f, 'B', {-1, -1}, {1, 4}},
    {  90.0f,   90.0f, 'S', {5, -1}, {1, -1}},
    {  90.0f,  120.0f, 'B', {2, 5}, {-1, -1}},
    {  90.0f,  150.0f, 'S', {2, -1}, {0, -1}},
    {  90.0f,  180.0f, 'B', {-1, -1}, {0, 3}},
    {  90.0f,  210.0f, 'S', {1, -1}, {3, -1}},
    {  90.0f,  240.0f, 'B', {1, 4}, {-1, -1}},
    {  90.0f,  270.0f, 'S', {1, -1}, {5, -1}},
    {  90.0f,  300.0f, 'B', {-1, -1}, {2, 5}},
    {  90.0f,  330.0f, 'S', {0, -1}, {2, -1}},
    {  54.7f,    0.0f, 'E', {0, -1}, {-1, -1}},
    {  50.8f,   30.0f, 'H', {0, -1}, {1, -1}},
    {  54.7f,   60.0f, 'E', {-1, -1}, {1, -1}},
    {  50.8f,   90.0f, 'H', {2, -1}, {1, -1}},
    {  54.7f,  120.0f, 'E', {2, -1}, {-1, -1}},
    {  50.8f,  150.0f, 'H', {2, -1}, {3, -1}},
    {  54.7f,  180.0f, 'E', {-1, -1}, {3, -1}},
    {  50.8f,  210.0f, 'H', {4, -1}, {3, -1}},
    {  54.7f,  240.0f, 'E', {4, -1}, {-1, -1}},
    {  50.8f,  270.0f, 'H', {4, -1}, {5, -1}},
    {  54.7f,  300.0f, 'E', {-1, -1}, {5, -1}},
    {  50.8f,  330.0f, 'H', {0, -1}, {5, -1}},
    { 125.3f,    0.0f, 'E', {3, -1}, {-1, -1}},
    { 129.2f,   30.0f, 'H', {3, -1}, {4, -1}},
    { 125.3f,   60.0f, 'E', {-1, -1}, {4, -1}},
    { 129.2f,   90.0f, 'H', {5, -1}, {4, -1}},
    { 125.3f,  120.0f, 'E', {5, -1}, {-1, -1}},
    { 129.2f,  150.0f, 'H', {5, -1}, {0, -1}},
    { 125.3f,  180.0f, 'E', {-1, -1}, {0, -1}},
    { 129.2f,  210.0f, 'H', {1, -1}, {0, -1}},
    { 125.3f,  240.0f, 'E', {1, -1}, {-1, -1}},
    { 129.2f,  270.0f, 'H', {1, -1}, {2, -1}},
    { 125.3f,  300.0f, 'E', {-1, -1}, {2, -1}},
    { 129.2f,  330.0f, 'H', {3, -1}, {2, -1}},
};

// Furanose conformations (IUPAC), at their Cremer-Pople phi2: positions 0 (the O) to 4, as above
static const pucker_canon_t pucker_canon5[] = {
    {   0.0f,    0.0f, 'E', {0, -1}, {-1, -1}},
    {   0.0f,   18.0f, 'T', {0, -1}, {1, -1}},
    {   0.0f,   36.0f, 'E', {-1, -1}, {1, -1}},
    {   0.0f,   54.0f, 'T', {2, -1}, {1, -1}},
    {   0.0f,   72.0f, 'E', {2, -1}, {-1, -1}},
    {   0.0f,   90.0f, 'T', {2, -1}, {3, -1}},
    {   0.0f,  108.0f, 'E', {-1, -1}, {3, -1}},
    {   0.0f,  126.0f, 'T', {4, -1}, {3, -1}},
    {   0.0f,  144.0f, 'E', {4, -1}, {-1, -1}},
    {   0.0f,  162.0f, 'T', {4, -1}, {0, -1}},
    {   0.0f,  180.0f, 'E', {-1, -1}, {0, -1}},
    {   0.0f,  198.0f, 'T', {1, -1}, {0, -1}},
    {   0.0f,  216.0f, 'E', {1, -1}, {-1, -1}},
    {   0.0f,  234.0f, 'T', {1, -1}, {2, -1}},
    {   0.0f,  252.0f, 'E', {-1, -1}, {2, -1}},
    {   0.0f,  270.0f, 'T', {3, -1}, {2, -1}},
    {   0.0f,  288.0f, 'E', {3, -1}, {-1, -1}},
    {   0.0f,  306.0f, 'T', {3, -1}, {4, -1}},
    {   0.0f,  324.0f, 'E', {-1, -1}, {4, -1}},
    {   0.0f,  342.0f, 'T', {0, -1}, {4, -1}},
};

static inline double deg(double rad) { return rad * (180.0 / 3.14159265358979323846); }
static inline double rad(double deg) { return deg * (3.14159265358979323846 / 180.0); }

static inline double wrap360(double a) {
    a = fmod(a, 360.0);
    return a < 0.0 ? a + 360.0 : a;
}

bool md_pucker_cremer_pople(md_pucker_t* out, const vec3_t* xyz, size_t count) {
    ASSERT(out);
    MEMSET(out, 0, sizeof(md_pucker_t));
    if (!xyz || count < 5 || count > 7) return false;
    const int N = (int)count;

    double c[3] = {0};
    for (int j = 0; j < N; ++j) {
        c[0] += xyz[j].x; c[1] += xyz[j].y; c[2] += xyz[j].z;
    }
    for (int k = 0; k < 3; ++k) c[k] /= N;

    double r[7][3];
    double rp[3] = {0}, rpp[3] = {0};
    for (int j = 0; j < N; ++j) {
        r[j][0] = xyz[j].x - c[0];
        r[j][1] = xyz[j].y - c[1];
        r[j][2] = xyz[j].z - c[2];
        const double s = sin(2.0 * PI * j / N), co = cos(2.0 * PI * j / N);
        for (int k = 0; k < 3; ++k) {
            rp[k]  += r[j][k] * s;
            rpp[k] += r[j][k] * co;
        }
    }
    double n[3] = {
        rp[1] * rpp[2] - rp[2] * rpp[1],
        rp[2] * rpp[0] - rp[0] * rpp[2],
        rp[0] * rpp[1] - rp[1] * rpp[0],
    };
    const double len = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (!(len > 1.0e-12)) return false;
    for (int k = 0; k < 3; ++k) n[k] /= len;

    double z[7];
    double Q2 = 0.0;
    for (int j = 0; j < N; ++j) {
        z[j] = r[j][0] * n[0] + r[j][1] * n[1] + r[j][2] * n[2];
        Q2 += z[j] * z[j];
    }

    // The harmonics m = 2 .. (N - 1) / 2
    double qm[2] = {0}, phim[2] = {0};
    for (int m = 2; m <= (N - 1) / 2; ++m) {
        double a = 0.0, b = 0.0;
        for (int j = 0; j < N; ++j) {
            a += z[j] * cos(2.0 * PI * m * j / N);
            b -= z[j] * sin(2.0 * PI * m * j / N);
        }
        a *= sqrt(2.0 / N);
        b *= sqrt(2.0 / N);
        qm[m - 2]   = sqrt(a * a + b * b);
        phim[m - 2] = wrap360(deg(atan2(b, a)));
    }

    out->size      = (uint32_t)N;
    out->amplitude = (float)sqrt(Q2);
    out->q2        = (float)qm[0];
    out->phi2      = (float)phim[0];
    if (N == 6) {
        double q3 = 0.0;
        for (int j = 0; j < N; ++j) q3 += (j & 1) ? -z[j] : z[j];
        q3 *= sqrt(1.0 / N);
        out->q3    = (float)q3;
        out->theta = (float)deg(atan2(qm[0], q3));
    } else if (N == 7) {
        out->q3   = (float)qm[1];
        out->phi3 = (float)phim[1];
    }
    return true;
}

static double torsion(const vec3_t a, const vec3_t b, const vec3_t c, const vec3_t d) {
    const double b0[3] = {a.x - b.x, a.y - b.y, a.z - b.z};
    const double b1[3] = {c.x - b.x, c.y - b.y, c.z - b.z};
    const double b2[3] = {d.x - c.x, d.y - c.y, d.z - c.z};
    const double l1 = sqrt(b1[0] * b1[0] + b1[1] * b1[1] + b1[2] * b1[2]);
    const double u[3] = {b1[0] / l1, b1[1] / l1, b1[2] / l1};
    const double d0 = b0[0] * u[0] + b0[1] * u[1] + b0[2] * u[2];
    const double d2 = b2[0] * u[0] + b2[1] * u[1] + b2[2] * u[2];
    const double v[3] = {b0[0] - d0 * u[0], b0[1] - d0 * u[1], b0[2] - d0 * u[2]};
    const double w[3] = {b2[0] - d2 * u[0], b2[1] - d2 * u[1], b2[2] - d2 * u[2]};
    const double x = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
    const double uv[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
    const double y = uv[0] * w[0] + uv[1] * w[1] + uv[2] * w[2];
    return deg(atan2(y, x));
}

bool md_pucker_pseudorotation(float* out_phase, float* out_amplitude, const vec3_t xyz[5]) {
    if (!xyz) return false;
    const vec3_t O = xyz[0], C1 = xyz[1], C2 = xyz[2], C3 = xyz[3], C4 = xyz[4];
    const double nu0 = torsion(C4, O, C1, C2);
    const double nu1 = torsion(O, C1, C2, C3);
    const double nu2 = torsion(C1, C2, C3, C4);
    const double nu3 = torsion(C2, C3, C4, O);
    const double nu4 = torsion(C3, C4, O, C1);
    const double num = (nu4 + nu1) - (nu3 + nu0);
    const double den = 2.0 * nu2 * (sin(rad(36.0)) + sin(rad(72.0)));
    if (num == 0.0 && den == 0.0) return false;
    const double P = wrap360(deg(atan2(num, den)));
    const double cp = cos(rad(P));
    if (out_phase) *out_phase = (float)P;
    // tau_m from nu2 when cos P is not small, otherwise from the amplitude of all five
    double tm = 0.0;
    if (fabs(cp) > 0.3) {
        tm = nu2 / cp;
    } else {
        const double nu[5] = {nu0, nu1, nu2, nu3, nu4};
        double s = 0.0;
        for (int j = 0; j < 5; ++j) {
            const double c = cos(rad(P + 144.0 * (j - 2)));
            s += c * c;
            tm += nu[j] * c;
        }
        tm = s > 0.0 ? tm / s : 0.0;
    }
    if (out_amplitude) *out_amplitude = (float)tm;
    return true;
}

static md_pucker_kind_t kind_of(char k) {
    switch (k) {
    case 'C': return MD_PUCKER_KIND_CHAIR;
    case 'B': return MD_PUCKER_KIND_BOAT;
    case 'S': return MD_PUCKER_KIND_SKEW;
    case 'H': return MD_PUCKER_KIND_HALF_CHAIR;
    case 'E': return MD_PUCKER_KIND_ENVELOPE;
    case 'T': return MD_PUCKER_KIND_TWIST;
    default:  return MD_PUCKER_KIND_NONE;
    }
}

static md_pucker_conformer_t conformer_of(const pucker_canon_t* c) {
    md_pucker_conformer_t r = {kind_of(c->kind), {c->above[0], c->above[1]}, {c->below[0], c->below[1]}};
    return r;
}

md_pucker_conformer_t md_pucker_conformer(const md_pucker_t* p) {
    md_pucker_conformer_t none = {MD_PUCKER_KIND_NONE, {-1, -1}, {-1, -1}};
    if (!p || p->size < 5 || p->size > 6) return none;
    if (p->amplitude < MD_PUCKER_PLANAR_AMPLITUDE) {
        none.kind = MD_PUCKER_KIND_PLANAR;
        return none;
    }
    if (p->size == 6) {
        // Nearest on the sphere
        const double t = rad(p->theta), f = rad(p->phi2);
        double best = -2.0;
        size_t bi = 0;
        for (size_t i = 0; i < ARRAY_SIZE(pucker_canon6); ++i) {
            const double tc = rad(pucker_canon6[i].theta), fc = rad(pucker_canon6[i].phi);
            const double c = cos(t) * cos(tc) + sin(t) * sin(tc) * cos(f - fc);
            if (c > best) { best = c; bi = i; }
        }
        return conformer_of(&pucker_canon6[bi]);
    }
    // Nearest along phi2
    double best = 1.0e9;
    size_t bi = 0;
    for (size_t i = 0; i < ARRAY_SIZE(pucker_canon5); ++i) {
        double d = fabs(wrap360(p->phi2 - pucker_canon5[i].phi + 180.0) - 180.0);
        if (d < best) { best = d; bi = i; }
    }
    return conformer_of(&pucker_canon5[bi]);
}

const char* md_pucker_kind_name(md_pucker_kind_t kind) {
    switch (kind) {
    case MD_PUCKER_KIND_PLANAR:     return "planar";
    case MD_PUCKER_KIND_CHAIR:      return "chair";
    case MD_PUCKER_KIND_BOAT:       return "boat";
    case MD_PUCKER_KIND_SKEW:       return "skew";
    case MD_PUCKER_KIND_HALF_CHAIR: return "half chair";
    case MD_PUCKER_KIND_ENVELOPE:   return "envelope";
    case MD_PUCKER_KIND_TWIST:      return "twist";
    default:                        return "";
    }
}

static char kind_letter(md_pucker_kind_t kind) {
    switch (kind) {
    case MD_PUCKER_KIND_CHAIR:      return 'C';
    case MD_PUCKER_KIND_BOAT:       return 'B';
    case MD_PUCKER_KIND_SKEW:       return 'S';
    case MD_PUCKER_KIND_HALF_CHAIR: return 'H';
    case MD_PUCKER_KIND_ENVELOPE:   return 'E';
    case MD_PUCKER_KIND_TWIST:      return 'T';
    default:                        return '?';
    }
}

// Numbers first, ascending, then the others (the O of "3,O")
static bool locant_before(const char* a, const char* b) {
    const bool da = a[0] >= '0' && a[0] <= '9';
    const bool db = b[0] >= '0' && b[0] <= '9';
    if (da != db) return da;
    if (da) return atoi(a) < atoi(b);
    return strcmp(a, b) < 0;
}

static size_t write_locants(char* buf, size_t cap, const int8_t pos[2], const char* const* locants) {
    const char* l[2] = {0};
    int n = 0;
    for (int i = 0; i < 2; ++i) if (pos[i] >= 0) l[n++] = locants[pos[i]];
    if (n == 2 && !locant_before(l[0], l[1])) { const char* t = l[0]; l[0] = l[1]; l[1] = t; }
    size_t len = 0;
    for (int i = 0; i < n; ++i) {
        len += (size_t)md_snprintf(buf + MIN(len, cap), cap - MIN(len, cap), "%s%s", i ? "," : "", l[i]);
    }
    return len;
}

size_t md_pucker_canonical(md_pucker_canonical_t* out, size_t cap, size_t size) {
    const pucker_canon_t* table = size == 6 ? pucker_canon6 : size == 5 ? pucker_canon5 : NULL;
    const size_t count = size == 6 ? ARRAY_SIZE(pucker_canon6) : size == 5 ? ARRAY_SIZE(pucker_canon5) : 0;
    for (size_t i = 0; out && i < MIN(cap, count); ++i) {
        out[i].theta = table[i].theta;
        out[i].phi2  = table[i].phi;
        out[i].conf  = conformer_of(&table[i]);
    }
    return count;
}

size_t md_pucker_locants(char* buf, size_t cap, const int8_t pos[2], const char* const* locants) {
    if (!buf || cap == 0) return 0;
    buf[0] = '\0';
    if (!pos || !locants) return 0;
    return write_locants(buf, cap, pos, locants);
}

size_t md_pucker_conformer_name(char* buf, size_t cap, md_pucker_conformer_t c, const char* const* locants) {
    if (!buf || cap == 0) return 0;
    buf[0] = '\0';
    if (c.kind == MD_PUCKER_KIND_NONE) return 0;
    if (!locants || c.kind == MD_PUCKER_KIND_PLANAR) {
        return (size_t)md_snprintf(buf, cap, "%s", md_pucker_kind_name(c.kind));
    }
    size_t len = write_locants(buf, cap, c.above, locants);
    len += (size_t)md_snprintf(buf + MIN(len, cap), cap - MIN(len, cap), "%c", kind_letter(c.kind));
    len += write_locants(buf + MIN(len, cap), cap - MIN(len, cap), c.below, locants);
    return len;
}

// Heavy atoms bonded to atom a which are not in the ring
static int exocyclic_heavy(const md_system_t* sys, int32_t a, const int32_t* ring, size_t n) {
    int count = 0;
    md_bond_iter_t it = md_bond_iter(&sys->bond, a);
    while (md_bond_iter_has_next(&it)) {
        const md_atom_idx_t x = md_bond_iter_atom_index(&it);
        md_bond_iter_next(&it);
        bool in_ring = false;
        for (size_t k = 0; k < n; ++k) in_ring |= ring[k] == x;
        if (!in_ring && md_atom_atomic_number(&sys->atom, x) != MD_Z_H) count += 1;
    }
    return count;
}

size_t md_pucker_ring_order(int32_t* out, bool* out_anchored, const md_system_t* sys, size_t ring_idx) {
    ASSERT(out);
    ASSERT(sys);
    if (out_anchored) *out_anchored = false;
    if (ring_idx >= md_index_data_num_ranges(&sys->ring)) return 0;
    const int32_t* ring = md_index_range_beg(&sys->ring, ring_idx);
    const size_t n = md_index_range_size(&sys->ring, ring_idx);
    if (n == 0 || n > MD_RING_MAX_SIZE) return 0;

    if (md_system_ring_flags(sys, ring_idx) & MD_RING_FLAG_SUGAR) {
        MEMCPY(out, ring, n * sizeof(int32_t));
        if (out_anchored) *out_anchored = true;
        return n;
    }

    int exo[MD_RING_MAX_SIZE];
    int z[MD_RING_MAX_SIZE];
    bool hetero = false;
    for (size_t k = 0; k < n; ++k) {
        exo[k] = exocyclic_heavy(sys, ring[k], ring, n);
        z[k] = md_atom_atomic_number(&sys->atom, ring[k]);
        hetero |= z[k] != MD_Z_C;
    }

    // The first atom: the heaviest atom other than C (else any C), then the most heavy atoms outside of the ring
    size_t start = 0;
    int ties = 0;
    for (size_t k = 0; k < n; ++k) {
        if (hetero && z[k] == MD_Z_C) continue;
        if (ties == 0) { start = k; ties = 1; continue; }
        const int zk = hetero ? z[k] : 0, zs = hetero ? z[start] : 0;
        if (zk > zs || (zk == zs && exo[k] > exo[start])) { start = k; ties = 1; }
        else if (zk == zs && exo[k] == exo[start]) {
            ties += 1;
            if (ring[k] < ring[start]) start = k;
        }
    }
    bool anchored = ties == 1;

    // The direction: towards the neighbor with the most heavy atoms outside of the ring, then the next ones
    int dir = 1;
    bool decided = false;
    for (size_t step = 1; step <= n / 2 && !decided; ++step) {
        const int fwd = exo[(start + step) % n];
        const int bwd = exo[(start + n - step) % n];
        if (fwd != bwd) {
            dir = fwd > bwd ? 1 : -1;
            decided = true;
        }
    }
    anchored = anchored && decided;

    for (size_t k = 0; k < n; ++k) {
        int idx = ((int)start + dir * (int)k) % (int)n;
        if (idx < 0) idx += (int)n;
        out[k] = ring[idx];
    }
    if (out_anchored) *out_anchored = anchored;
    return n;
}

// The locant of the atom at position k of a sugar ring: "O" for the ring O, the number in the atom's name otherwise
static void sugar_locant(char out[8], const md_system_t* sys, int32_t atom, size_t k) {
    if (k == 0) { out[0] = 'O'; out[1] = '\0'; return; }
    const str_t name = md_atom_name(&sys->atom, atom);
    size_t n = 0;
    for (size_t i = 0; i < name.len && n < 7; ++i) {
        if (name.ptr[i] >= '0' && name.ptr[i] <= '9') out[n++] = name.ptr[i];
    }
    if (n == 0) n = (size_t)md_snprintf(out, 8, "%d", (int)k);
    out[n] = '\0';
}

static void atom_label(char* out, size_t cap, const md_system_t* sys, int32_t atom) {
    const str_t s = md_atom_name(&sys->atom, atom);
    md_snprintf(out, cap, "%.*s", (int)MIN(s.len, cap - 1), s.ptr);
}

static void describe(md_pucker_ring_t* r, const md_system_t* sys, const md_system_state_t* state, const vec3_t* xyz, vec3_t c) {
    const md_pucker_conformer_t conf = md_pucker_conformer(&r->cp);
    char* buf = r->name;
    const size_t cap = sizeof(r->name);
    const size_t n = r->size;
    buf[0] = '\0';
    if (conf.kind == MD_PUCKER_KIND_NONE) return;
    if (conf.kind == MD_PUCKER_KIND_PLANAR || (!r->anchored && !r->sugar)) {
        md_snprintf(buf, cap, "%s", md_pucker_kind_name(conf.kind));
        return;
    }
    char name[7][8];
    const char* locants[7];
    for (size_t k = 0; k < n; ++k) {
        if (r->sugar) sugar_locant(name[k], sys, r->atoms[k], k);
        else          atom_label(name[k], sizeof(name[k]), sys, r->atoms[k]);
        locants[k] = name[k];
    }
    if (r->sugar) {
        size_t len = md_pucker_conformer_name(buf, cap, conf, locants);
        len += (size_t)md_snprintf(buf + MIN(len, cap), cap - MIN(len, cap), " %s", md_pucker_kind_name(conf.kind));
        if (n != 5 || (conf.kind != MD_PUCKER_KIND_ENVELOPE && conf.kind != MD_PUCKER_KIND_TWIST)) return;
        // Endo: the side of the C outside of the ring on the last C (C5' of a nucleotide)
        int32_t ref = -1;
        md_bond_iter_t it = md_bond_iter(&sys->bond, r->atoms[4]);
        while (md_bond_iter_has_next(&it)) {
            const md_atom_idx_t x = md_bond_iter_atom_index(&it);
            md_bond_iter_next(&it);
            bool in_ring = false;
            for (size_t k = 0; k < n; ++k) in_ring |= r->atoms[k] == x;
            if (!in_ring && md_atom_atomic_number(&sys->atom, x) == MD_Z_C) ref = x;
        }
        if (ref < 0) return;
        // Above: the side from which the ring is numbered clockwise, against the normal of the polygon
        vec3_t nw = {0};
        for (size_t k = 0; k < n; ++k) nw = vec3_add(nw, vec3_cross(vec3_sub(xyz[k], c), vec3_sub(xyz[(k + 1) % n], c)));
        vec3_t d = vec3_sub(state->xyz[ref], xyz[4]);
        md_util_min_image_vec3(&d, 1, &state->unitcell);
        const bool ref_above = vec3_dot(vec3_add(vec3_sub(xyz[4], c), d), nw) < 0.0f;
        char an[16] = "", bn[16] = "";
        if (conf.above[0] >= 0) atom_label(an, sizeof(an), sys, r->atoms[conf.above[0]]);
        if (conf.below[0] >= 0) atom_label(bn, sizeof(bn), sys, r->atoms[conf.below[0]]);
        const char* up = ref_above ? "endo" : "exo";
        const char* dn = ref_above ? "exo" : "endo";
        if (an[0] && bn[0]) len += (size_t)md_snprintf(buf + MIN(len, cap), cap - MIN(len, cap), ", %s-%s %s-%s", an, up, bn, dn);
        else if (an[0])     len += (size_t)md_snprintf(buf + MIN(len, cap), cap - MIN(len, cap), ", %s-%s", an, up);
        else if (bn[0])     len += (size_t)md_snprintf(buf + MIN(len, cap), cap - MIN(len, cap), ", %s-%s", bn, dn);
        return;
    }
    // Another ring: the atoms out of the plane by name
    size_t len = (size_t)md_snprintf(buf, cap, "%s (", md_pucker_kind_name(conf.kind));
    bool first = true;
    for (int i = 0; i < 2; ++i) if (conf.above[i] >= 0) {
        len += (size_t)md_snprintf(buf + MIN(len, cap), cap - MIN(len, cap), "%s%s above", first ? "" : ", ", locants[conf.above[i]]);
        first = false;
    }
    for (int i = 0; i < 2; ++i) if (conf.below[i] >= 0) {
        len += (size_t)md_snprintf(buf + MIN(len, cap), cap - MIN(len, cap), "%s%s below", first ? "" : ", ", locants[conf.below[i]]);
        first = false;
    }
    md_snprintf(buf + MIN(len, cap), cap - MIN(len, cap), ")");
}

static bool compute_ring(md_pucker_ring_t* r, const md_system_t* sys, const md_system_state_t* state) {
    const size_t n = r->size;
    vec4_t xyzw[7];
    for (size_t k = 0; k < n; ++k) {
        const vec3_t x = state->xyz[r->atoms[k]];
        xyzw[k] = vec4_set(x.x, x.y, x.z, 1.0f);
    }
    vec3_t com = {0};
    md_util_deperiodize_self_vec4(xyzw, n, &state->unitcell, &com);
    vec3_t xyz[7];
    vec3_t c = {0};
    for (size_t k = 0; k < n; ++k) {
        xyz[k] = vec3_from_vec4(xyzw[k]);
        c = vec3_add(c, xyz[k]);
    }
    c = vec3_mul1(c, 1.0f / (float)n);
    if (!md_pucker_cremer_pople(&r->cp, xyz, n)) {
        r->cp.size = (uint32_t)n;
        return false;
    }
    if (r->has_pseudo && n == 5) {
        md_pucker_pseudorotation(&r->pseudo_phase, &r->pseudo_amplitude, xyz);
    }
    describe(r, sys, state, xyz, c);
    return true;
}

bool md_pucker_ring(md_pucker_ring_t* out, const md_system_t* sys, const md_system_state_t* state, size_t ring_idx) {
    ASSERT(out);
    MEMSET(out, 0, sizeof(md_pucker_ring_t));
    if (!sys || !state || !state->xyz || ring_idx >= md_index_data_num_ranges(&sys->ring)) return false;
    const size_t n = md_index_range_size(&sys->ring, ring_idx);
    if (n < 5 || n > 7) return false;
    bool anchored = false;
    md_pucker_ring_order(out->atoms, &anchored, sys, ring_idx);
    out->size = (uint8_t)n;
    out->anchored = anchored;
    out->sugar = (md_system_ring_flags(sys, ring_idx) & MD_RING_FLAG_SUGAR) != 0;
    out->has_pseudo = out->sugar && n == 5;
    return compute_ring(out, sys, state);
}

bool md_pucker_ring_atoms(md_pucker_ring_t* out, const md_system_t* sys, const md_system_state_t* state, const int32_t* atoms, size_t count) {
    ASSERT(out);
    MEMSET(out, 0, sizeof(md_pucker_ring_t));
    if (!sys || !state || !state->xyz || !atoms || count < 5 || count > 7) return false;
    MEMCPY(out->atoms, atoms, count * sizeof(int32_t));
    out->size = (uint8_t)count;
    out->anchored = true;
    out->has_pseudo = count == 5;
    return compute_ring(out, sys, state);
}
