// What is added up over the frames: the volumes and rdf, behind one interface for the interpreter (accum.h)

#include "accum.h"

#include <core/md_common.h>

#include <string.h>

bool vial_accum_builtin(const char* name) {
    return vial_volume_builtin(name) || !strcmp(name, "rdf");
}

bool vial_accum_begin(vial_kctx_t* k, vial_accum_t* acc, const char* name, const vial_val_t* const* args, const md_str_t* names, uint32_t n,
                      const vial_frame_t* first, md_allocator_i* alloc) {
    md_memset(acc, 0, sizeof(*acc));
    if (!strcmp(name, "rdf")) {
        acc->kind = VIAL_ACCUM_RDF;
        const vial_args_t a = {args, names, n};
        return vial_rdf_begin(k, &acc->rdf, &a);
    }
    acc->kind = VIAL_ACCUM_VOLUME;
    return vial_volume_begin(k, &acc->vol, name, args, names, n, first, true, alloc);
}

bool vial_accum_add(vial_kctx_t* k, vial_accum_t* acc, const vial_val_t* const* args, const md_str_t* names, uint32_t n) {
    switch (acc->kind) {
    case VIAL_ACCUM_VOLUME: return vial_volume_add(k, &acc->vol, args, names, n);
    case VIAL_ACCUM_RDF: {
        const vial_args_t a = {args, names, n};
        return vial_rdf_add(k, &acc->rdf, &a);
    }
    default: return false;
    }
}

bool vial_accum_end(vial_kctx_t* k, vial_accum_t* acc, vial_val_t* out) {
    switch (acc->kind) {
    case VIAL_ACCUM_VOLUME: return vial_volume_end(k, &acc->vol, out);
    case VIAL_ACCUM_RDF:    return vial_rdf_end(k, &acc->rdf, out);
    default:                return false;
    }
}

bool vial_accum_like(vial_accum_t* out, const vial_accum_t* acc, md_allocator_i* alloc) {
    md_memset(out, 0, sizeof(*out));
    out->kind = acc->kind;
    switch (acc->kind) {
    case VIAL_ACCUM_VOLUME: return vial_volume_like(&out->vol, &acc->vol, alloc);
    case VIAL_ACCUM_RDF:    return vial_rdf_like(&out->rdf, &acc->rdf);
    default:                return false;
    }
}

void vial_accum_free(vial_accum_t* acc) {
    if (!acc) return;
    switch (acc->kind) {
    case VIAL_ACCUM_VOLUME: vial_volume_free(&acc->vol); break;
    case VIAL_ACCUM_RDF:    vial_rdf_free(&acc->rdf); break;
    default: break;
    }
    acc->kind = 0;
}

size_t vial_accum_count(const vial_accum_t* acc) {
    if (!acc) return 0;
    switch (acc->kind) {
    case VIAL_ACCUM_VOLUME: return (size_t)acc->vol.dims[0] * acc->vol.dims[1] * acc->vol.dims[2];
    case VIAL_ACCUM_RDF:    return acc->rdf.bins;
    default:                return 0;
    }
}
