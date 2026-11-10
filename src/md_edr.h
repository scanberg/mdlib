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
	md_str_t name;
	md_str_t unit_str;
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

bool md_edr_energies_parse_file(md_edr_energies_t* energies, md_str_t filename, struct md_allocator_i* alloc);
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
bool md_edr_system_supplement(struct md_system_t* sys, const md_edr_energies_t* energies, md_str_t run);

// The same from a file, which additionally publishes "<run>/edr/source": the path it was loaded
// from, so that whoever saves the session can find what to load again. It lives and dies with the
// energies it describes, so it can never name a file whose data is gone.
bool md_edr_system_supplement_from_file(struct md_system_t* sys, md_str_t filename, md_str_t run);

// THE GROUPS OF A TERM. GROMACS names the terms it writes for groups by the groups' names:
//
//     Coul-SR:Protein-SOL     a non-bonded term between two energy groups (energygrps), Coul-SR, LJ-SR, Buck-SR,
//                             Coul-14 and LJ-14, for every pair in the order of energygrps, a group with itself too
//     T-Protein, Lamb-Protein the temperature and the scaling factor of a temperature coupling group (tc-grps)
//
// A group's name may hold '-' (Protein-H, non-Water), so the two of a pair are found among the groups of the system
// (md_group_data_t): the energy groups of its role, in order, where it has them (a tpr), else any groups of those names
// (an index file). Of a pair whose names could be split in more than one way, the split of two energy groups in their
// order is taken. The groups of the T-coupling role likewise, else by name.
//
// term is the name as the file has it, which the published attribute keeps as its label. False when it is not a term of
// groups, or a group it names is not in the system ('rest' of a role the system does not have, say).
typedef struct md_edr_term_groups_t {
    uint32_t role;          // md_group_role_t: MD_GROUP_ROLE_ENERGY for a pair, MD_GROUP_ROLE_TC for T- and Lamb-
    uint32_t count;         // 2 for a pair, 1 for the term of one group
    uint32_t group[2];      // Indices into sys->group, the pair in the term's order; a pair of a group with itself has it twice
} md_edr_term_groups_t;

bool md_edr_term_groups(md_edr_term_groups_t* out, const struct md_system_t* sys, md_str_t term);

// The groups of several terms, each once, in the order the terms first name them: what drawing an array of terms shows,
// such as every Coul-SR term, where a group is in many pairs and is one group all the same. Terms of no groups add none.
// Writes at most cap group indices and returns how many there are.
size_t md_edr_terms_groups(uint32_t out_groups[], size_t cap, const struct md_system_t* sys, const md_str_t terms[], size_t count);

// THE ROLES OF THE GROUPS, FROM THE TERMS. A tpr says which groups were the energy groups and which the T-coupling groups
// (md_group_data_t); an index file only what groups there are. Of a system without the role, an energy file tells it by
// its terms: every pair of energy groups has its terms, each group with itself too (Coul-SR:SOL-SOL), in the order of
// energygrps, so those of a group with itself give the names in order, whatever '-' they hold. The T-coupling groups
// likewise, from T-<name>. Each name is found among the groups of the system (the first of the name, as grompp takes
// it, else the first of the name in another case), and 'rest', the role's own, is the atoms in none of the others.
//
// A role is set only where every name is found and the groups found do not overlap: a role is a partition, which grompp
// makes sure of, so overlapping groups are of another index file than the one grompp was given, and nothing is set
// rather than something wrong. A role the system has (from a tpr) is left as it is. Returns the number of roles set.
// term_names are the names of the terms as the file has them (md_edr_energy_t::name, or an attribute's label).
size_t md_edr_system_infer_group_roles(struct md_system_t* sys, const md_str_t* term_names, size_t count);

// The same from the terms published for a run ("<run>/edr", their labels): for a system an index file was added to
// after its energy file. md_edr_system_supplement does it itself.
size_t md_edr_system_infer_group_roles_of_run(struct md_system_t* sys, md_str_t run);

#if 0
// Represents the string block inside the edr file which gives us the names and units (encoded as text)
// of all energy fields present in the file.
typedef struct md_edr_str_block_t {
	int64_t count;
	md_str_t* names;
	md_str_t* units;
} md_edr_str_block_t;

bool md_edr_str_block_read_file(md_edr_str_block_t* str_block, md_str_t filename, struct md_allocator_i* alloc);
void md_edr_str_block_free(md_edr_str_block_t* str_block, struct md_allocator_i* alloc);
#endif


#ifdef __cplusplus
}
#endif
