#include <md_types.h>
#include <md_util.h>
#include <md_element.h>

#include <core/md_hash.h>
#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_format.h>
#include <core/md_log.h>

#include <string.h>
#include <ctype.h>

typedef struct {
    md_str_t str;
    md_atomic_number_t z;
} Mapping;


// Predefined mappings (applies to a combination from res_name:atom_name)
static const Mapping predefined_mappings[] = {
    // Common ions
    {MD_STR_INIT("NA:NA"), MD_Z_Na},
    {MD_STR_INIT("SOD:SOD"), MD_Z_Na},
    {MD_STR_INIT("SOD:NA"),MD_Z_Na},
    {MD_STR_INIT("K:K"),   MD_Z_K},
    {MD_STR_INIT("POT:K"), MD_Z_K},
    {MD_STR_INIT("POT:POT"), MD_Z_K},  // CHARMM names the ion as its residue
    {MD_STR_INIT("CL:CL"), MD_Z_Cl},
    {MD_STR_INIT("CLA:CL"),MD_Z_Cl},
    {MD_STR_INIT("CLA:CLA"),MD_Z_Cl},
    {MD_STR_INIT("CA:CA"), MD_Z_Ca},
    {MD_STR_INIT("CAL:CA"),MD_Z_Ca},
    {MD_STR_INIT("MG:MG"), MD_Z_Mg},
    {MD_STR_INIT("MAG:MG"),MD_Z_Mg},
    {MD_STR_INIT("ZN:ZN"), MD_Z_Zn},
    {MD_STR_INIT("FE:FE"), MD_Z_Fe},
    {MD_STR_INIT("CU:CU"), MD_Z_Cu},
    {MD_STR_INIT("MN:MN"), MD_Z_Mn},
    {MD_STR_INIT("F: F"),  MD_Z_F},
    {MD_STR_INIT("BR:BR"), MD_Z_Br},
    {MD_STR_INIT("I: I"),  MD_Z_I},
    {MD_STR_INIT("IOD:I"), MD_Z_I},

    // Other
    {MD_STR_INIT("POPC:HS"), MD_Z_H},
    {MD_STR_INIT("HOH:O"), MD_Z_O},
    {MD_STR_INIT("WAT:O"), MD_Z_O},
    {MD_STR_INIT("H2O:O"), MD_Z_O},
    {MD_STR_INIT("DOD:O"), MD_Z_O}, // Heavy water
    {MD_STR_INIT("CO3:C"), MD_Z_C}, // Carbonate
    {MD_STR_INIT("CO3:O"), MD_Z_O},
    {MD_STR_INIT("CO2:C"), MD_Z_C}, // Carbon dioxide
    {MD_STR_INIT("CO2:O"), MD_Z_O},
    {MD_STR_INIT("NH4:N"), MD_Z_N}, // Ammonium
    {MD_STR_INIT("NH4:H"), MD_Z_H},
    {MD_STR_INIT("NO3:N"), MD_Z_N}, // Nitrate
    {MD_STR_INIT("NO3:O"), MD_Z_O},
    {MD_STR_INIT("SO4:S"), MD_Z_S}, // Sulfate
    {MD_STR_INIT("SO4:O"), MD_Z_O},

    // Hem
    {MD_STR_INIT("HEM:FE"), MD_Z_Fe},
};

static const md_str_t resname_fallback_to_single_character[] = {
    MD_STR_INIT("ALA"), MD_STR_INIT("ARG"),
    MD_STR_INIT("ASN"), MD_STR_INIT("ASP"), MD_STR_INIT("CYS"), MD_STR_INIT("GLN"),
    MD_STR_INIT("GLU"), MD_STR_INIT("GLY"), MD_STR_INIT("HIS"), MD_STR_INIT("ILE"), MD_STR_INIT("LEU"), MD_STR_INIT("LYS"), MD_STR_INIT("MET"), MD_STR_INIT("PHE"), MD_STR_INIT("PRO"), MD_STR_INIT("SER"), MD_STR_INIT("THR"), MD_STR_INIT("TRP"), MD_STR_INIT("TYR"), MD_STR_INIT("VAL"),
    MD_STR_INIT("DA"),  MD_STR_INIT("DC"),  MD_STR_INIT("DG"),  MD_STR_INIT("DT"),  MD_STR_INIT("A"),   MD_STR_INIT("C"),   MD_STR_INIT("G"),   MD_STR_INIT("U"),
    MD_STR_INIT("HOH"), MD_STR_INIT("H2O"), MD_STR_INIT("DOD"), MD_STR_INIT("WAT"),
    
    MD_STR_INIT("HEM"), MD_STR_INIT("HEME"),
    MD_STR_INIT("ADP"), MD_STR_INIT("ATP"), MD_STR_INIT("GDP"), MD_STR_INIT("GTP"),

    MD_STR_INIT("POPC"), MD_STR_INIT("DPPC"), MD_STR_INIT("DMPC"), MD_STR_INIT("DOPC"), MD_STR_INIT("POPE"), MD_STR_INIT("DOPE"),

    MD_STR_INIT("CHOL"),  // Cholesterol
};

static inline md_atomic_number_t find_in_mappings(md_str_t key, const Mapping* mappings, size_t num_mappings) {
    for (size_t i = 0; i < num_mappings; ++i) {
        if (md_str_eq_ignore_case(key, mappings[i].str)) {
            return mappings[i].z;
        }
    }
    return MD_Z_X;
}

static inline bool resname_in_fallback_list(md_str_t res_name) {
    for (size_t i = 0; i < MD_COUNTOF(resname_fallback_to_single_character); ++i) {
        if (md_str_eq_ignore_case(res_name, resname_fallback_to_single_character[i])) {
            return true;
        }
    }
    return false;
}

// Inference functions
md_atomic_number_t md_atomic_number_infer_from_label(md_str_t atom_name, md_str_t res_name, size_t res_size) {    

    // Special case: if atom name is empty but residue name is an element (ion case)
    if (atom_name.len == 0 && res_name.len > 0) {
        md_atomic_number_t res_element = md_atomic_number_from_symbol(res_name, true);
        if (res_element != MD_Z_X) {
            return res_element;
        }
        return MD_Z_X;
    }
    
    if (atom_name.len == 0) return MD_Z_X;

    if (md_ascii_is_digit(atom_name.ptr[0])) {
        atom_name = md_str_substr(atom_name, 1, SIZE_MAX); // Skip leading digits
    }
    
    // First try residue+atom combination
    if (res_name.len > 0) {
        // Special case: if water, amino acid or nucleotide: Match against first character only
        if (md_util_resname_water(res_name) || md_util_resname_amino_acid(res_name) || md_util_resname_nucleotide(res_name)) {
            md_str_t first_char = md_str_substr(atom_name, 0, 1);
            md_atomic_number_t z = md_atomic_number_from_symbol(first_char, true);
            if (z != MD_Z_X) return z;
        }

		char combined_key_buf[16];
		int len = md_snprintf(combined_key_buf, sizeof(combined_key_buf), MD_STR_FMT ":" MD_STR_FMT, MD_STR_ARG(res_name), MD_STR_ARG(atom_name));
		md_str_t combined_key = { combined_key_buf, (size_t)len };

		// Check predefined resname : atomlabel mappings first
        {
        md_atomic_number_t z = find_in_mappings(combined_key, predefined_mappings, MD_COUNTOF(predefined_mappings));
            if (z != MD_Z_X) return z;
        }

        // Check if we should fallback to single-character atom name
        if (resname_in_fallback_list(res_name)) {
            md_str_t first_char = md_str_substr(atom_name, 0, 1);
            md_atomic_number_t z = md_atomic_number_from_symbol(first_char, true);
            if (z != MD_Z_X) return z;
        }
        
        // If residue name itself is an element (ions for example)
        md_atomic_number_t res_element = md_atomic_number_from_symbol(res_name, true);
        if (res_element != MD_Z_X) {
            // If atom name is empty or equals residue, return that element
            if (md_str_eq_ignore_case(atom_name, res_name)) {
                return res_element;
            }
        }
    }

    // User has provided a residue size and its likely not an ion
    if (res_size > 2) {
        md_str_t two_letter_str = md_str_substr(atom_name, 0, 2);
        md_atomic_number_t two_z = md_atomic_number_from_symbol(two_letter_str, false);
        if (two_z != MD_Z_X) {
            return two_z;
        }

        // Try single-letter element heuristic (e.g., C1 => C, N2 => N)
        md_str_t first_letter_str = md_str_substr(atom_name, 0, 1);
        md_atomic_number_t first_z = md_atomic_number_from_symbol(first_letter_str, true);
        if (first_z != MD_Z_X) {
            return first_z;
        }
    }
    
    // Try two-letter element heuristic (e.g., CL12 => Cl, BR1 => Br)
    if (md_str_len(atom_name) >= 2) {
        md_str_t two_letter_str = md_str_substr(atom_name, 0, 2);
        md_atomic_number_t two_z = md_atomic_number_from_symbol(two_letter_str, true);
        if (two_z != MD_Z_X) {
            return two_z;
        }
    }
    
    // Final fallback: first-letter element mapping
    md_str_t first_letter_str = md_str_substr(atom_name, 0, 1);
    return md_atomic_number_from_symbol(first_letter_str, true);
}

md_atomic_number_t md_atomic_number_infer_from_mass(float mass) {
    md_atomic_number_t res = MD_Z_X;
    md_atomic_number_infer_from_mass_batch(&res, &mass, 1);
    return res;
}

