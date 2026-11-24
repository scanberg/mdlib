#pragma once

#include <md_system.h>
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
// without types or residues, of which state holds the coordinates
static inline void md_io_points_system(md_system_t* sys, size_t num_atoms) {
    md_system_reset(sys);
    sys->atom.count = num_atoms;
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
