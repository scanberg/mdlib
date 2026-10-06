/* md_gpu_select.inl — adapter selection shared by the md_gpu backends.
   Included once by md_gpu_vulkan.c and md_gpu_metal.m (internal, not API). */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

/* The selector of a device description: desc->adapter, else $MD_GPU_DEVICE,
   with surrounding whitespace removed. Returns false when there is none. */
static bool md_gpu_sel_get(const md_gpu_device_desc_t* desc, char* out, size_t cap) {
    const char* s = (desc && desc->adapter && desc->adapter[0]) ? desc->adapter : getenv("MD_GPU_DEVICE");
    out[0] = 0;
    if (!s) return false;
    while (*s && isspace((unsigned char)*s)) ++s;
    size_t n = strlen(s);
    while (n && isspace((unsigned char)s[n - 1])) --n;
    if (n == 0) return false;
    if (n >= cap) n = cap - 1;
    memcpy(out, s, n);
    out[n] = 0;
    return true;
}

/* All digits: an index into the adapter list. */
static bool md_gpu_sel_index(const char* sel, uint32_t* out_index) {
    if (!sel[0]) return false;
    uint32_t v = 0;
    for (const char* p = sel; *p; ++p) {
        if (!isdigit((unsigned char)*p)) return false;
        v = v * 10u + (uint32_t)(*p - '0');
    }
    *out_index = v;
    return true;
}

static bool md_gpu_sel_matches(const char* sel, uint32_t index, const char* name) {
    uint32_t want;
    if (md_gpu_sel_index(sel, &want)) return want == index;
    const size_t n = strlen(sel);
    for (const char* h = name; *h; ++h) {
        size_t k = 0;
        while (k < n && h[k] && tolower((unsigned char)h[k]) == tolower((unsigned char)sel[k])) ++k;
        if (k == n) return true;
    }
    return false;
}

/* Higher is better. */
static int md_gpu_sel_rank(md_gpu_device_type_t type, md_gpu_device_preference_t pref) {
    const bool low_power = pref == MD_GPU_DEVICE_PREFER_LOW_POWER;
    switch (type) {
    case MD_GPU_DEVICE_TYPE_DISCRETE:   return low_power ? 50 : 100;
    case MD_GPU_DEVICE_TYPE_INTEGRATED: return low_power ? 100 : 50;
    case MD_GPU_DEVICE_TYPE_VIRTUAL:    return 30;
    case MD_GPU_DEVICE_TYPE_OTHER:      return 20;
    case MD_GPU_DEVICE_TYPE_CPU:        return 10;
    }
    return 0;
}

static const char* md_gpu_sel_type_str(md_gpu_device_type_t type) {
    switch (type) {
    case MD_GPU_DEVICE_TYPE_DISCRETE:   return "discrete";
    case MD_GPU_DEVICE_TYPE_INTEGRATED: return "integrated";
    case MD_GPU_DEVICE_TYPE_VIRTUAL:    return "virtual";
    case MD_GPU_DEVICE_TYPE_CPU:        return "cpu";
    case MD_GPU_DEVICE_TYPE_OTHER:      return "other";
    }
    return "?";
}

/* Picks an adapter from `list`. Returns its index, or -1 with a description of
   why in `why` (which names the adapters, so the caller's error is useful). */
static int md_gpu_sel_pick(const md_gpu_adapter_info_t* list, uint32_t count, const md_gpu_device_desc_t* desc,
                           char* why, size_t why_cap) {
    char sel[128];
    const bool have_sel = md_gpu_sel_get(desc, sel, sizeof(sel));
    const md_gpu_device_preference_t pref = desc ? desc->preference : MD_GPU_DEVICE_PREFER_DEFAULT;

    int best = -1, best_rank = -1;
    int matched_unusable = -1;
    for (uint32_t i = 0; i < count; ++i) {
        if (have_sel && !md_gpu_sel_matches(sel, i, list[i].name)) continue;
        if (!list[i].usable) { if (matched_unusable < 0) matched_unusable = (int)i; continue; }
        const int r = md_gpu_sel_rank(list[i].type, pref);
        if (r > best_rank) { best_rank = r; best = (int)i; }
    }
    if (best >= 0) return best;

    /* Nothing: explain, naming every adapter. */
    size_t len = 0;
    if (have_sel && matched_unusable >= 0) {
        len += (size_t)snprintf(why + len, why_cap - len, "adapter '%s' (selected by '%s') lacks %s",
                                list[matched_unusable].name, sel, list[matched_unusable].missing);
    } else if (have_sel) {
        len += (size_t)snprintf(why + len, why_cap - len, "no adapter matches '%s'", sel);
    } else {
        len += (size_t)snprintf(why + len, why_cap - len, "no usable adapter");
    }
    if (len < why_cap) len += (size_t)snprintf(why + len, why_cap - len, "; adapters:");
    for (uint32_t i = 0; i < count && len < why_cap; ++i) {
        len += (size_t)snprintf(why + len, why_cap - len, " [%u] %s (%s%s%s)", i, list[i].name,
                                md_gpu_sel_type_str(list[i].type), list[i].usable ? "" : ", lacks ",
                                list[i].usable ? "" : list[i].missing);
    }
    if (count == 0 && len < why_cap) snprintf(why + len, why_cap - len, " none");
    return -1;
}
