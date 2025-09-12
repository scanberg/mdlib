#include <core/md_atomic.h>
#include <md_molecule.h>
#include <md_util.h>

#include <core/md_hash.h>
#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_str.h>

#include <string.h>
#include <ctype.h>

// Core atomic number functions using existing md_util tables
md_atomic_number_t md_atomic_number_from_symbol(str_t sym) {
    return md_util_element_lookup(sym);
}

md_atomic_number_t md_atomic_number_from_symbol_icase(str_t sym) {
    return md_util_element_lookup_ignore_case(sym);
}

str_t md_symbol_from_atomic_number(md_atomic_number_t z) {
    return md_util_element_symbol(z);
}

str_t md_name_from_atomic_number(md_atomic_number_t z) {
    return md_util_element_name(z);
}

float md_atomic_mass(md_atomic_number_t z) {
    return md_util_element_atomic_mass(z);
}

float md_vdw_radius(md_atomic_number_t z) {
    return md_util_element_vdw_radius(z);
}

float md_covalent_radius(md_atomic_number_t z) {
    return md_util_element_covalent_radius(z);
}

int md_max_valence(md_atomic_number_t z) {
    return md_util_element_max_valence(z);
}

uint32_t md_cpk_color(md_atomic_number_t z) {
    return md_util_element_cpk_color(z);
}

// Inference functions
md_atomic_number_t md_atom_infer_atomic_number(str_t atom_name, str_t res_name) {    

    // Special case: if atom name is empty but residue name is an element (ion case)
    if (atom_name.len == 0 && res_name.len > 0) {
        md_atomic_number_t res_element = md_atomic_number_from_symbol_icase(res_name);
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
            md_atomic_number_t z = md_atomic_number_from_symbol_icase(first_char);
            if (z != MD_Z_X) return z;
        }
        
        // If residue name itself is an element (ion case)
        md_atomic_number_t res_element = md_atomic_number_from_symbol_icase(res_name);
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
        md_atomic_number_t two_z = md_atomic_number_from_symbol_icase(two_letter_str);
        if (two_z != MD_Z_X) {
            return two_z;
        }
    }
    
    // Final fallback: first-letter element mapping
    str_t first_letter_str = str_substr(atom_name, 0, 1);
    return md_atomic_number_from_symbol_icase(first_letter_str);
}

bool md_atoms_infer_atomic_numbers(md_atomic_number_t out[], size_t n, const struct md_molecule_t* mol) {
    if (!out || !mol || n == 0) return false;
    
    size_t count = MIN(n, mol->atom.count);
    for (size_t i = 0; i < count; ++i) {
        str_t atom_name = LBL_TO_STR(mol->atom.type[i]);
        str_t res_name = mol->atom.resname ? LBL_TO_STR(mol->atom.resname[i]) : (str_t){0};
        out[i] = md_atom_infer_atomic_number(atom_name, res_name);
    }
    
    return true;
}