#pragma once

#include <stdbool.h>

#include <core/md_str.h>

struct md_system_t;
struct md_load_desc_t;

#ifdef __cplusplus
extern "C" {
#endif

// AMBER ENERGIES: the mdout log of sander and pmemd (written every ntpr steps), and the mden file (every ntwe steps)
//
// Published into the system's attribute table under "<run>/amber_energy", replacing what an earlier energy file put
// there. The file keeps its OWN frame axis, as an EDR file does (md_edr.h):
//
//     <run>/amber_energy/time      F64 {R}   ps, the records
//     <run>/amber_energy/<term>    F64 {R}   one per term present in every record
//     <run>/amber_energy/source    STR       the file it was loaded from
//
// <term> is the mdout name lowercased with every run of other characters folded to '_' ("EPtot" -> "eptot",
// "1-4 EEL" -> "1_4_eel", "TEMP(K)" -> "temp_k"), and the original is the label. An mden file publishes its terms
// under the mdout names where mdout has them (E_pot -> eptot, E_vdw -> vdwaals), so that a script reads either.
// Energies are in kcal/mol, temperatures in K, pressures in bar, volume in Angstrom^3, density in g/cm^3.
//
// Of an mdout file only the records of the run are read (section 4. RESULTS), not its averages and fluctuations.
// Minimization output has no time and is refused. When "<run>/time" exists every frame of it has to find a record,
// and nothing is published otherwise: the file belongs to some other run.

// An mdout file: '4.  RESULTS' or an 'NSTEP =' record near its top
bool md_amber_energy_file_is_mdout(md_str_t filename);

// An mden file: it starts with "L0  Nsteps"
bool md_amber_energy_file_is_mden(md_str_t filename);

bool md_amber_energy_system_supplement_from_file(struct md_system_t* sys, md_str_t filename, const struct md_load_desc_t* desc);

#ifdef __cplusplus
}
#endif
