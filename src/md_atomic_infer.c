#include <md_types.h>
#include <md_util.h>

#include <core/md_hash.h>
#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_log.h>

#include <string.h>
#include <ctype.h>

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

// Inference functions
md_atomic_number_t md_atomic_number_infer_from_label(str_t atom_name, str_t res_name) {    

    // Special case: if atom name is empty but residue name is an element (ion case)
    if (atom_name.len == 0 && res_name.len > 0) {
        md_atomic_number_t res_element = md_atomic_number_from_symbol(res_name, true);
        if (res_element != MD_Z_X) {
            return res_element;
        }
        return MD_Z_X;
    }
    
    if (atom_name.len == 0) return MD_Z_X;
    
    // First try residue+atom combination
    if (res_name.len > 0) {       
        // Special case: if water, amino acid or nucleotide: Match against first character only
        if (md_util_resname_water(res_name) || md_util_resname_amino_acid(res_name) || md_util_resname_nucleotide(res_name)) {
            str_t first_char = str_substr(atom_name, 0, 1);
            md_atomic_number_t z = md_atomic_number_from_symbol(first_char, true);
            if (z != MD_Z_X) return z;
        }
        
        // If residue name itself is an element (ion case)
        md_atomic_number_t res_element = md_atomic_number_from_symbol(res_name, true);
        if (res_element != MD_Z_X) {
            // If atom name is empty or equals residue, return that element
            if (atom_name.len == 0 || str_eq_ignore_case(atom_name, res_name)) {
                return res_element;
            }
        }
    }
    
    // Try two-letter element heuristic (e.g., CL12 => Cl, BR1 => Br)
    if (atom_name.len >= 2) {
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

size_t md_atomic_number_infer_from_label_batch(md_atomic_number_t out_z[], const str_t atom_names[], const str_t atom_resnames[], size_t count) {
    if (!out_z) {
        MD_LOG_ERROR("Output array is null");
        return 0;
    }
    if (!atom_names) {
        MD_LOG_ERROR("Input atom names array is null");
        return 0;
    }

    size_t temp_pos = md_temp_get_pos();
    md_hashmap32_t cache = {.allocator = md_get_temp_allocator() };
    md_hashmap_reserve(&cache, 256);
    
    size_t num_success = 0;
    for (size_t i = 0; i < count; ++i) {
        str_t atom_name = atom_names[i];
        str_t res_name = atom_resnames ? atom_resnames[i] : STR_LIT("");
        uint64_t key = md_hash64_str(atom_name, md_hash64_str(res_name, 0));

        uint32_t* cached_z = md_hashmap_get(&cache, key);
        if (cached_z) {
            out_z[i] = (md_atomic_number_t)*cached_z;
            num_success += (size_t)(*cached_z != MD_Z_X);
            continue;
        }

        out_z[i] = md_atomic_number_infer_from_label(atom_name, res_name);
        md_hashmap_add(&cache, key, out_z[i]);
        num_success += (size_t)(out_z[i] != MD_Z_X);
    }
    
    md_temp_set_pos_back(temp_pos);

    return num_success;
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
