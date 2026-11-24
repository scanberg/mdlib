#pragma once

#include <md_system.h>
#include <md_element.h>
#include <core/md_array.h>
#include <core/md_log.h>

// For the readers (md_io.h, LOAD). Not part of the interface.

// A system of its own builds onto a system without atoms: onto one with atoms it is MD_LOAD_FAILED, and nothing changes
static inline bool md_io_onto_no_atoms(const md_system_t* sys, md_str_t path) {
    if (sys->atom.count == 0) return true;
    MD_LOG_ERROR("'" MD_STR_FMT "' is a system of its own: it builds onto a system without atoms, and joins one of its atoms through md_system_merge", MD_STR_ARG(path));
    return false;
}

// Data of a system's atoms loads onto a system with atoms: onto one without it is MD_LOAD_FAILED
static inline bool md_io_onto_atoms(const md_system_t* sys, md_str_t path) {
    if (sys->atom.count != 0) return true;
    MD_LOG_ERROR("'" MD_STR_FMT "' is data of a system's atoms: load a system first", MD_STR_ARG(path));
    return false;
}

// A system known by its atom count alone (md_io.h, LOAD: coordinates without names, as a restart has them): N atoms
// without residues, of which state holds the coordinates. Each of one type, of no element: what reads the atoms through
// their types finds one (type 0 is the unknown, as every loader has it, of no size; this one is drawn at the size of an
// atom of no element). It names them no more than their count does: a system of named atoms they join takes the names.
static inline void md_io_points_system(md_system_t* sys, size_t num_atoms) {
    md_system_reset(sys);
    sys->atom.count = num_atoms;
    md_particle_type_find_or_add(&sys->atom.type, MD_STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);
    const md_particle_type_idx_t point = md_particle_type_add(&sys->atom.type, MD_STR_LIT("X"), (md_str_t){0}, 0, 0.0f,
        md_atomic_number_vdw_radius(0), md_atomic_number_cpk_color(0), 0, sys->alloc);
    md_array_resize(sys->atom.type_idx, num_atoms, sys->alloc);
    for (size_t i = 0; i < num_atoms; ++i) sys->atom.type_idx[i] = point;
}

// Whether the system has the run, its frame axis "<run>/time": what a series is along
static inline bool md_io_has_run(const md_system_t* sys, md_str_t run) {
    if (md_str_empty(run)) return false;
    const md_attribute_t* axis = md_attributes_find_in(&sys->attributes, run, MD_STR_LIT("time"));
    return axis && md_attributes_axis(&sys->attributes, axis) == axis;
}

// A series which does not cover the run: of another run
static inline void md_io_log_not_along(const char* format, md_str_t path, md_str_t run) {
    MD_LOG_ERROR("%s: '" MD_STR_FMT "' is not along '" MD_STR_FMT "': the run has frames at times it holds no values of", format, MD_STR_ARG(path), MD_STR_ARG(run));
}
