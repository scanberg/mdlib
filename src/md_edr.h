#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_unit.h>

struct md_allocator_i;
struct md_system_t;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_edr_energy_t {
	str_t name;
	str_t unit_str;
    md_unit_t unit;
    // Has length num_frames. Double because the file may be: a total energy of -6e5 kJ/mol kept in
    // float moves in steps of 0.06, which is the size of the drift anyone plotting it is looking for.
    double* values;
} md_edr_energy_t;

typedef struct md_edr_energies_t {
	size_t			 num_frames;
	double*			 frame_time;	// timestamp of frames

	size_t			 num_energies;
	md_edr_energy_t* energy;
	
	struct md_allocator_i* alloc;
} md_edr_energies_t;

bool md_edr_energies_parse_file(md_edr_energies_t* energies, str_t filename, struct md_allocator_i* alloc);
void md_edr_energies_free(md_edr_energies_t* energies);

// Publishes the energies into the system's attribute table under "<run>/edr", replacing whatever
// an earlier energy file put there. The file keeps its OWN frame axis rather than being resampled
// onto the trajectory's, so nothing it holds is dropped when it was written more often than the
// coordinates; a consumer relates the two by time (md_attribute_axis_map).
//
//     <run>/edr/time          F64 {R}       ps, the file's frames
//     <run>/edr/<term>        F64 {R}       one per energy term, in the term's unit
//     <run>/edr/<tensor>      F64 {R,3,3}   nine terms <name>-XX .. <name>-ZZ, e.g. "vir", "pres"
//     <run>/edr/<vector>      F64 {R}  c3   three terms <name>-X/-Y/-Z, e.g. "box", or the diagonal
//                                           <name>-XX/-YY/-ZZ when that is all the file holds
//
// <term> is the GROMACS name lowercased with every run of other characters folded to '_'
// ("Kinetic En." -> "kinetic_en", "dVcoul/dl" -> "dvcoul_dl"), and the original is the label.
// Terms are only grouped when every member is present once and all share a unit.
//
// When "<run>/time" exists every one of its frames has to find a row in the file, and nothing is
// published otherwise: an energy file that does not cover the trajectory belongs to some other
// run. Without it the energies simply form the run. Frame times must not decrease.
bool md_edr_system_supplement(struct md_system_t* sys, const md_edr_energies_t* energies, str_t run);

// The same from a file, which additionally publishes "<run>/edr/source": the path it was loaded
// from, so that whoever saves the session can find what to load again. It lives and dies with the
// energies it describes, so it can never name a file whose data is gone.
bool md_edr_system_supplement_from_file(struct md_system_t* sys, str_t filename, str_t run);

#if 0
// Represents the string block inside the edr file which gives us the names and units (encoded as text)
// of all energy fields present in the file.
typedef struct md_edr_str_block_t {
	int64_t count;
	str_t* names;
	str_t* units;
} md_edr_str_block_t;

bool md_edr_str_block_read_file(md_edr_str_block_t* str_block, str_t filename, struct md_allocator_i* alloc);
void md_edr_str_block_free(md_edr_str_block_t* str_block, struct md_allocator_i* alloc);
#endif


#ifdef __cplusplus
}
#endif
