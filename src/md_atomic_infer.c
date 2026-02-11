#include <md_types.h>
#include <md_util.h>

#include <core/md_hash.h>
#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_log.h>

#include <string.h>
#include <ctype.h>

typedef struct {
    str_t str;
    md_atomic_number_t z;
} Mapping;

#define STR(str) {str, sizeof(str)-1}

// Predefined mappings (applies to a combination from res_name:atom_name)
static const Mapping predefined_mappings[] = {
    // Common ions
    {STR("NA:NA"), MD_Z_Na},
    {STR("SOD:SOD"), MD_Z_Na},
    {STR("SOD:NA"),MD_Z_Na},
    {STR("K:K"),   MD_Z_K},
    {STR("POT:K"), MD_Z_K},
    {STR("CL:CL"), MD_Z_Cl},
    {STR("CLA:CL"),MD_Z_Cl},
    {STR("CLA:CLA"),MD_Z_Cl},
    {STR("CA:CA"), MD_Z_Ca},
    {STR("CAL:CA"),MD_Z_Ca},
    {STR("MG:MG"), MD_Z_Mg},
    {STR("MAG:MG"),MD_Z_Mg},
    {STR("ZN:ZN"), MD_Z_Zn},
    {STR("FE:FE"), MD_Z_Fe},
    {STR("CU:CU"), MD_Z_Cu},
    {STR("MN:MN"), MD_Z_Mn},
    {STR("F: F"),  MD_Z_F},
    {STR("BR:BR"), MD_Z_Br},
    {STR("I: I"),  MD_Z_I},
    {STR("IOD:I"), MD_Z_I},

    // Other
    {STR("POPC:HS"), MD_Z_H},
    {STR("HOH:O"), MD_Z_O},
    {STR("WAT:O"), MD_Z_O},
    {STR("H2O:O"), MD_Z_O},
    {STR("DOD:O"), MD_Z_O}, // Heavy water
    {STR("CO3:C"), MD_Z_C}, // Carbonate
    {STR("CO3:O"), MD_Z_O},
    {STR("CO2:C"), MD_Z_C}, // Carbon dioxide
    {STR("CO2:O"), MD_Z_O},
    {STR("NH4:N"), MD_Z_N}, // Ammonium
    {STR("NH4:H"), MD_Z_H},
    {STR("NO3:N"), MD_Z_N}, // Nitrate
    {STR("NO3:O"), MD_Z_O},
    {STR("SO4:S"), MD_Z_S}, // Sulfate
    {STR("SO4:O"), MD_Z_O},

    // Hem
    {STR("HEM:FE"), MD_Z_Fe},
};

static const str_t resname_fallback_to_single_character[] = {
    STR("ALA"), STR("ARG"),
    STR("ASN"), STR("ASP"), STR("CYS"), STR("GLN"),
    STR("GLU"), STR("GLY"), STR("HIS"), STR("ILE"), STR("LEU"), STR("LYS"), STR("MET"), STR("PHE"), STR("PRO"), STR("SER"), STR("THR"), STR("TRP"), STR("TYR"), STR("VAL"),
    STR("DA"),  STR("DC"),  STR("DG"),  STR("DT"),  STR("A"),   STR("C"),   STR("G"),   STR("U"),
    STR("HOH"), STR("H2O"), STR("DOD"), STR("WAT"),
    
    STR("HEM"), STR("HEME"),
    STR("ADP"), STR("ATP"), STR("GDP"), STR("GTP"),

    STR("POPC"), STR("DPPC"), STR("DMPC"), STR("DOPC"), STR("POPE"), STR("DOPE"),

    STR("CHOL"),  // Cholesterol
};

// Core atomic number functions using existing md_util tables
md_atomic_number_t md_atomic_number_from_symbol(str_t sym, bool ignore_case) {
    return md_util_element_lookup(sym, ignore_case);
}

str_t md_atomic_number_symbol(md_atomic_number_t z) {
    return md_util_element_symbol(z);
}

str_t md_atomic_number_name(md_atomic_number_t z) {
    return md_util_element_name(z);
}

float md_atomic_number_mass(md_atomic_number_t z) {
    return md_util_element_atomic_mass(z);
}

float md_atomic_number_vdw_radius(md_atomic_number_t z) {
    return md_util_element_vdw_radius(z);
}

float md_atomic_number_covalent_radius(md_atomic_number_t z) {
    return md_util_element_covalent_radius(z);
}

int md_atomic_number_max_valence(md_atomic_number_t z) {
    return md_util_element_max_valence(z);
}

uint32_t md_atomic_number_cpk_color(md_atomic_number_t z) {
    return md_util_element_cpk_color(z);
}

static inline md_atomic_number_t find_in_mappings(str_t key, const Mapping* mappings, size_t num_mappings) {
    for (size_t i = 0; i < num_mappings; ++i) {
        if (str_eq_ignore_case(key, mappings[i].str)) {
            return mappings[i].z;
        }
    }
    return MD_Z_X;
}

static inline bool resname_in_fallback_list(str_t res_name) {
    for (size_t i = 0; i < ARRAY_SIZE(resname_fallback_to_single_character); ++i) {
        if (str_eq_ignore_case(res_name, resname_fallback_to_single_character[i])) {
            return true;
        }
    }
    return false;
}

// Inference functions
md_atomic_number_t md_atomic_number_infer_from_label(str_t atom_name, str_t res_name, size_t res_size) {    

    // Special case: if atom name is empty but residue name is an element (ion case)
    if (atom_name.len == 0 && res_name.len > 0) {
        md_atomic_number_t res_element = md_atomic_number_from_symbol(res_name, true);
        if (res_element != MD_Z_X) {
            return res_element;
        }
        return MD_Z_X;
    }
    
    if (atom_name.len == 0) return MD_Z_X;

    if (is_digit(atom_name.ptr[0])) {
        atom_name = str_substr(atom_name, 1, SIZE_MAX); // Skip leading digits
    }
    
    // First try residue+atom combination
    if (res_name.len > 0) {
        // Special case: if water, amino acid or nucleotide: Match against first character only
        if (md_util_resname_water(res_name) || md_util_resname_amino_acid(res_name) || md_util_resname_nucleotide(res_name)) {
            str_t first_char = str_substr(atom_name, 0, 1);
            md_atomic_number_t z = md_atomic_number_from_symbol(first_char, true);
            if (z != MD_Z_X) return z;
        }

		char combined_key_buf[16];
		int len = snprintf(combined_key_buf, sizeof(combined_key_buf), STR_FMT ":" STR_FMT, STR_ARG(res_name), STR_ARG(atom_name));
		str_t combined_key = { combined_key_buf, (size_t)len };

		// Check predefined resname : atomlabel mappings first
        {
        md_atomic_number_t z = find_in_mappings(combined_key, predefined_mappings, ARRAY_SIZE(predefined_mappings));
            if (z != MD_Z_X) return z;
        }

        // Check if we should fallback to single-character atom name
        if (resname_in_fallback_list(res_name)) {
            str_t first_char = str_substr(atom_name, 0, 1);
            md_atomic_number_t z = md_atomic_number_from_symbol(first_char, true);
            if (z != MD_Z_X) return z;
        }
        
        // If residue name itself is an element (ions for example)
        md_atomic_number_t res_element = md_atomic_number_from_symbol(res_name, true);
        if (res_element != MD_Z_X) {
            // If atom name is empty or equals residue, return that element
            if (str_eq_ignore_case(atom_name, res_name)) {
                return res_element;
            }
        }
    }

    // User has provided a component size and its likely not an ion
    if (res_size > 2) {
        // Try single-letter element heuristic (e.g., C1 => C, N2 => N)
        str_t first_letter_str = str_substr(atom_name, 0, 1);
        md_atomic_number_t first_z = md_atomic_number_from_symbol(first_letter_str, true);
        if (first_z != MD_Z_X) {
            return first_z;
        }
    }
    
    // Try two-letter element heuristic (e.g., CL12 => Cl, BR1 => Br)
    if (str_len(atom_name) >= 2) {
        str_t two_letter_str = str_substr(atom_name, 0, 2);
        md_atomic_number_t two_z = md_atomic_number_from_symbol(two_letter_str, true);
        if (two_z != MD_Z_X) {
            return two_z;
        }
    }
    
    // Final fallback: first-letter element mapping
    str_t first_letter_str = str_substr(atom_name, 0, 1);
    return md_atomic_number_from_symbol(first_letter_str, true);
}

md_atomic_number_t md_atomic_number_infer_from_mass(float mass) {
    md_atomic_number_t res = MD_Z_X;
    md_util_element_from_mass(&res, &mass, 1);
    return res;
}

size_t md_atomic_number_infer_from_mass_batch(md_atomic_number_t out_z[], const float masses[], size_t count) {
    if (!out_z) {
        MD_LOG_ERROR("Output array is null");
        return 0;
    }
    if (!masses) {
        MD_LOG_ERROR("Input masses array is null");
        return 0;
    }
    return md_util_element_from_mass(out_z, masses, count);
}
