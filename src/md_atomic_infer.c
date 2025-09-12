#include <core/md_atomic.h>
#include <md_molecule.h>
#include <md_util.h>

#include <core/md_hash.h>
#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_str.h>

#include <string.h>
#include <ctype.h>

// Static hashmaps for element inference
static md_hashmap32_t residue_atom_map = {0};
static md_hashmap32_t atom_only_map = {0};
static bool maps_initialized = false;

// Helper functions
static inline char to_upper_c(char c) {
    return (c >= 'a' && c <= 'z') ? (c - 'a' + 'A') : c;
}

static inline bool is_alpha_c(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

static inline bool is_digit_c(char c) {
    return c >= '0' && c <= '9';
}

// Strip digits from end of string
static str_t strip_digits(str_t str) {
    while (str.len > 0 && is_digit_c(str.ptr[str.len - 1])) {
        str.len--;
    }
    return str;
}

// Create normalized uppercase string for hashing
static void normalize_to_upper(char* dst, str_t src, bool strip_digits_flag) {
    str_t s = strip_digits_flag ? strip_digits(src) : src;
    for (size_t i = 0; i < s.len; ++i) {
        dst[i] = to_upper_c(s.ptr[i]);
    }
    dst[s.len] = '\0';
}

// Initialize lookup tables
static void init_lookup_tables(void) {
    if (maps_initialized) return;
    
    md_allocator_i* alloc = md_get_heap_allocator();
    residue_atom_map.allocator = alloc;
    atom_only_map.allocator = alloc;
    
    md_hashmap_reserve(&residue_atom_map, 512);
    md_hashmap_reserve(&atom_only_map, 256);
    
    // Water variants
    struct { const char* res; const char* atom; md_atomic_number_t z; } water_entries[] = {
        {"HOH", "O", MD_Z_O}, {"HOH", "OW", MD_Z_O}, {"HOH", "OH2", MD_Z_O},
        {"HOH", "H", MD_Z_H}, {"HOH", "H1", MD_Z_H}, {"HOH", "H2", MD_Z_H}, {"HOH", "HW1", MD_Z_H}, {"HOH", "HW2", MD_Z_H}, {"HOH", "HW", MD_Z_H},
        {"WAT", "O", MD_Z_O}, {"WAT", "OW", MD_Z_O}, {"WAT", "OH2", MD_Z_O},
        {"WAT", "H", MD_Z_H}, {"WAT", "H1", MD_Z_H}, {"WAT", "H2", MD_Z_H}, {"WAT", "HW1", MD_Z_H}, {"WAT", "HW2", MD_Z_H}, {"WAT", "HW", MD_Z_H},
        {"TIP3", "O", MD_Z_O}, {"TIP3", "OW", MD_Z_O}, {"TIP3", "OH2", MD_Z_O},
        {"TIP3", "H", MD_Z_H}, {"TIP3", "H1", MD_Z_H}, {"TIP3", "H2", MD_Z_H}, {"TIP3", "HW1", MD_Z_H}, {"TIP3", "HW2", MD_Z_H}, {"TIP3", "HW", MD_Z_H},
        {"TIP4", "O", MD_Z_O}, {"TIP4", "OW", MD_Z_O}, {"TIP4", "OH2", MD_Z_O},
        {"TIP4", "H", MD_Z_H}, {"TIP4", "H1", MD_Z_H}, {"TIP4", "H2", MD_Z_H}, {"TIP4", "HW1", MD_Z_H}, {"TIP4", "HW2", MD_Z_H}, {"TIP4", "HW", MD_Z_H},
        {"TIP5", "O", MD_Z_O}, {"TIP5", "OW", MD_Z_O}, {"TIP5", "OH2", MD_Z_O},
        {"TIP5", "H", MD_Z_H}, {"TIP5", "H1", MD_Z_H}, {"TIP5", "H2", MD_Z_H}, {"TIP5", "HW1", MD_Z_H}, {"TIP5", "HW2", MD_Z_H}, {"TIP5", "HW", MD_Z_H},
        {"SPC", "O", MD_Z_O}, {"SPC", "OW", MD_Z_O}, {"SPC", "OH2", MD_Z_O},
        {"SPC", "H", MD_Z_H}, {"SPC", "H1", MD_Z_H}, {"SPC", "H2", MD_Z_H}, {"SPC", "HW1", MD_Z_H}, {"SPC", "HW2", MD_Z_H}, {"SPC", "HW", MD_Z_H},
        {"SOL", "O", MD_Z_O}, {"SOL", "OW", MD_Z_O}, {"SOL", "OH2", MD_Z_O},
        {"SOL", "H", MD_Z_H}, {"SOL", "H1", MD_Z_H}, {"SOL", "H2", MD_Z_H}, {"SOL", "HW1", MD_Z_H}, {"SOL", "HW2", MD_Z_H}, {"SOL", "HW", MD_Z_H},
        {"H2O", "O", MD_Z_O}, {"H2O", "OW", MD_Z_O}, {"H2O", "OH2", MD_Z_O},
        {"H2O", "H", MD_Z_H}, {"H2O", "H1", MD_Z_H}, {"H2O", "H2", MD_Z_H}, {"H2O", "HW1", MD_Z_H}, {"H2O", "HW2", MD_Z_H}, {"H2O", "HW", MD_Z_H},
    };
    
    // Add water entries to residue+atom map
    for (size_t i = 0; i < ARRAY_SIZE(water_entries); ++i) {
        char key_str[64];
        snprintf(key_str, sizeof(key_str), "%s\t%s", water_entries[i].res, water_entries[i].atom);
        uint64_t key = md_hash64(key_str, strlen(key_str), 0);
        md_hashmap_add(&residue_atom_map, key, water_entries[i].z);
    }
    
    // Amino acid backbone and sidechain entries
    struct { const char* res; const char* atom; md_atomic_number_t z; } amino_entries[] = {
        // Common backbone atoms for all amino acids
        {"ALA", "CA", MD_Z_C}, {"ALA", "N", MD_Z_N}, {"ALA", "C", MD_Z_C}, {"ALA", "O", MD_Z_O}, {"ALA", "OXT", MD_Z_O},
        {"ARG", "CA", MD_Z_C}, {"ARG", "N", MD_Z_N}, {"ARG", "C", MD_Z_C}, {"ARG", "O", MD_Z_O}, {"ARG", "OXT", MD_Z_O},
        {"ASN", "CA", MD_Z_C}, {"ASN", "N", MD_Z_N}, {"ASN", "C", MD_Z_C}, {"ASN", "O", MD_Z_O}, {"ASN", "OXT", MD_Z_O},
        {"ASP", "CA", MD_Z_C}, {"ASP", "N", MD_Z_N}, {"ASP", "C", MD_Z_C}, {"ASP", "O", MD_Z_O}, {"ASP", "OXT", MD_Z_O},
        {"CYS", "CA", MD_Z_C}, {"CYS", "N", MD_Z_N}, {"CYS", "C", MD_Z_C}, {"CYS", "O", MD_Z_O}, {"CYS", "OXT", MD_Z_O},
        {"GLN", "CA", MD_Z_C}, {"GLN", "N", MD_Z_N}, {"GLN", "C", MD_Z_C}, {"GLN", "O", MD_Z_O}, {"GLN", "OXT", MD_Z_O},
        {"GLU", "CA", MD_Z_C}, {"GLU", "N", MD_Z_N}, {"GLU", "C", MD_Z_C}, {"GLU", "O", MD_Z_O}, {"GLU", "OXT", MD_Z_O},
        {"GLY", "CA", MD_Z_C}, {"GLY", "N", MD_Z_N}, {"GLY", "C", MD_Z_C}, {"GLY", "O", MD_Z_O}, {"GLY", "OXT", MD_Z_O},
        {"HIS", "CA", MD_Z_C}, {"HIS", "N", MD_Z_N}, {"HIS", "C", MD_Z_C}, {"HIS", "O", MD_Z_O}, {"HIS", "OXT", MD_Z_O},
        {"ILE", "CA", MD_Z_C}, {"ILE", "N", MD_Z_N}, {"ILE", "C", MD_Z_C}, {"ILE", "O", MD_Z_O}, {"ILE", "OXT", MD_Z_O},
        {"LEU", "CA", MD_Z_C}, {"LEU", "N", MD_Z_N}, {"LEU", "C", MD_Z_C}, {"LEU", "O", MD_Z_O}, {"LEU", "OXT", MD_Z_O},
        {"LYS", "CA", MD_Z_C}, {"LYS", "N", MD_Z_N}, {"LYS", "C", MD_Z_C}, {"LYS", "O", MD_Z_O}, {"LYS", "OXT", MD_Z_O},
        {"MET", "CA", MD_Z_C}, {"MET", "N", MD_Z_N}, {"MET", "C", MD_Z_C}, {"MET", "O", MD_Z_O}, {"MET", "OXT", MD_Z_O},
        {"PHE", "CA", MD_Z_C}, {"PHE", "N", MD_Z_N}, {"PHE", "C", MD_Z_C}, {"PHE", "O", MD_Z_O}, {"PHE", "OXT", MD_Z_O},
        {"PRO", "CA", MD_Z_C}, {"PRO", "N", MD_Z_N}, {"PRO", "C", MD_Z_C}, {"PRO", "O", MD_Z_O}, {"PRO", "OXT", MD_Z_O},
        {"SER", "CA", MD_Z_C}, {"SER", "N", MD_Z_N}, {"SER", "C", MD_Z_C}, {"SER", "O", MD_Z_O}, {"SER", "OXT", MD_Z_O},
        {"THR", "CA", MD_Z_C}, {"THR", "N", MD_Z_N}, {"THR", "C", MD_Z_C}, {"THR", "O", MD_Z_O}, {"THR", "OXT", MD_Z_O},
        {"TRP", "CA", MD_Z_C}, {"TRP", "N", MD_Z_N}, {"TRP", "C", MD_Z_C}, {"TRP", "O", MD_Z_O}, {"TRP", "OXT", MD_Z_O},
        {"TYR", "CA", MD_Z_C}, {"TYR", "N", MD_Z_N}, {"TYR", "C", MD_Z_C}, {"TYR", "O", MD_Z_O}, {"TYR", "OXT", MD_Z_O},
        {"VAL", "CA", MD_Z_C}, {"VAL", "N", MD_Z_N}, {"VAL", "C", MD_Z_C}, {"VAL", "O", MD_Z_O}, {"VAL", "OXT", MD_Z_O},
        // Sidechain specific atoms
        {"SER", "OG", MD_Z_O}, {"THR", "OG1", MD_Z_O}, {"TYR", "OH", MD_Z_O},
        {"CYS", "SG", MD_Z_S}, {"MET", "SD", MD_Z_S},
    };
    
    // Add amino acid entries to residue+atom map
    for (size_t i = 0; i < ARRAY_SIZE(amino_entries); ++i) {
        char key_str[64];
        snprintf(key_str, sizeof(key_str), "%s\t%s", amino_entries[i].res, amino_entries[i].atom);
        uint64_t key = md_hash64(key_str, strlen(key_str), 0);
        md_hashmap_add(&residue_atom_map, key, amino_entries[i].z);
    }
    
    // Nucleic acid entries
    struct { const char* res; const char* atom; md_atomic_number_t z; } nucleic_entries[] = {
        // DNA
        {"DA", "P", MD_Z_P}, {"DA", "OP1", MD_Z_O}, {"DA", "OP2", MD_Z_O}, {"DA", "O1P", MD_Z_O}, {"DA", "O2P", MD_Z_O},
        {"DC", "P", MD_Z_P}, {"DC", "OP1", MD_Z_O}, {"DC", "OP2", MD_Z_O}, {"DC", "O1P", MD_Z_O}, {"DC", "O2P", MD_Z_O},
        {"DG", "P", MD_Z_P}, {"DG", "OP1", MD_Z_O}, {"DG", "OP2", MD_Z_O}, {"DG", "O1P", MD_Z_O}, {"DG", "O2P", MD_Z_O},
        {"DT", "P", MD_Z_P}, {"DT", "OP1", MD_Z_O}, {"DT", "OP2", MD_Z_O}, {"DT", "O1P", MD_Z_O}, {"DT", "O2P", MD_Z_O},
        // RNA
        {"A", "P", MD_Z_P}, {"A", "OP1", MD_Z_O}, {"A", "OP2", MD_Z_O}, {"A", "O1P", MD_Z_O}, {"A", "O2P", MD_Z_O},
        {"C", "P", MD_Z_P}, {"C", "OP1", MD_Z_O}, {"C", "OP2", MD_Z_O}, {"C", "O1P", MD_Z_O}, {"C", "O2P", MD_Z_O},
        {"G", "P", MD_Z_P}, {"G", "OP1", MD_Z_O}, {"G", "OP2", MD_Z_O}, {"G", "O1P", MD_Z_O}, {"G", "O2P", MD_Z_O},
        {"U", "P", MD_Z_P}, {"U", "OP1", MD_Z_O}, {"U", "OP2", MD_Z_O}, {"U", "O1P", MD_Z_O}, {"U", "O2P", MD_Z_O},
    };
    
    // Add nucleic acid entries to residue+atom map
    for (size_t i = 0; i < ARRAY_SIZE(nucleic_entries); ++i) {
        char key_str[64];
        snprintf(key_str, sizeof(key_str), "%s\t%s", nucleic_entries[i].res, nucleic_entries[i].atom);
        uint64_t key = md_hash64(key_str, strlen(key_str), 0);
        md_hashmap_add(&residue_atom_map, key, nucleic_entries[i].z);
    }
    
    // Selenomethionine
    char mse_key[64];
    snprintf(mse_key, sizeof(mse_key), "MSE\tSE");
    uint64_t mse_hash = md_hash64(mse_key, strlen(mse_key), 0);
    md_hashmap_add(&residue_atom_map, mse_hash, MD_Z_Se);
    
    // Atom-only fallbacks
    struct { const char* atom; md_atomic_number_t z; } atom_entries[] = {
        {"H", MD_Z_H}, {"C", MD_Z_C}, {"N", MD_Z_N}, {"O", MD_Z_O}, {"S", MD_Z_S}, {"P", MD_Z_P},
        {"F", MD_Z_F}, {"CL", MD_Z_Cl}, {"BR", MD_Z_Br}, {"I", MD_Z_I},
        {"OW", MD_Z_O}, {"OH", MD_Z_O}, {"HW", MD_Z_H},
        {"CA", MD_Z_C}, {"CB", MD_Z_C}, {"CG", MD_Z_C}, {"CD", MD_Z_C}, {"CE", MD_Z_C}, {"CZ", MD_Z_C},
        {"OXT", MD_Z_O},
        // Common ions
        {"NA", MD_Z_Na}, {"K", MD_Z_K}, {"MG", MD_Z_Mg}, {"ZN", MD_Z_Zn}, {"FE", MD_Z_Fe},
        {"MN", MD_Z_Mn}, {"CU", MD_Z_Cu}, {"CO", MD_Z_Co}, {"NI", MD_Z_Ni}, {"CD", MD_Z_Cd},
        {"SR", MD_Z_Sr}, {"BA", MD_Z_Ba}, {"LI", MD_Z_Li}, {"CS", MD_Z_Cs}, {"RB", MD_Z_Rb},
        {"AL", MD_Z_Al}, {"TI", MD_Z_Ti}, {"CR", MD_Z_Cr}, {"HG", MD_Z_Hg}, {"PB", MD_Z_Pb},
        {"AG", MD_Z_Ag}, {"AU", MD_Z_Au}, {"PT", MD_Z_Pt},
    };
    
    // Add atom-only entries
    for (size_t i = 0; i < ARRAY_SIZE(atom_entries); ++i) {
        char norm_atom[16];
        str_t atom_str = {atom_entries[i].atom, strlen(atom_entries[i].atom)};
        normalize_to_upper(norm_atom, atom_str, true);
        uint64_t key = md_hash64(norm_atom, strlen(norm_atom), 0);
        md_hashmap_add(&atom_only_map, key, atom_entries[i].z);
    }
    
    maps_initialized = true;
}

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
    init_lookup_tables();
    
    // Special case: if atom name is empty but residue name is an element (ion case)
    if (atom_name.len == 0 && res_name.len > 0) {
        md_atomic_number_t res_element = md_atomic_number_from_symbol_icase(res_name);
        if (res_element != MD_Z_X) {
            return res_element;
        }
        return MD_Z_X;
    }
    
    if (atom_name.len == 0) return MD_Z_X;
    
    // Normalize inputs
    char norm_res[16] = {0};
    char norm_atom[16] = {0};
    char norm_atom_stripped[16] = {0};
    
    if (res_name.len > 0) {
        normalize_to_upper(norm_res, res_name, false);
    }
    normalize_to_upper(norm_atom, atom_name, false);
    normalize_to_upper(norm_atom_stripped, atom_name, true);
    
    // First try residue+atom combination
    if (res_name.len > 0) {
        char res_atom_key[64];
        snprintf(res_atom_key, sizeof(res_atom_key), "%s\t%s", norm_res, norm_atom);
        uint64_t key = md_hash64(res_atom_key, strlen(res_atom_key), 0);
        uint32_t* result = md_hashmap_get(&residue_atom_map, key);
        if (result) {
            return (md_atomic_number_t)*result;
        }
        
        // Special case: if residue is water
        if (md_util_resname_water(res_name)) {
            if (norm_atom[0] == 'O') return MD_Z_O;
            if (norm_atom[0] == 'H') return MD_Z_H;
        }
        
        // Special case: if residue is amino acid and atom is CA, return carbon
        if (md_util_resname_amino_acid(res_name) && strcmp(norm_atom, "CA") == 0) {
            return MD_Z_C;
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
    
    // Try atom-only map with digits stripped
    uint64_t atom_key = md_hash64(norm_atom_stripped, strlen(norm_atom_stripped), 0);
    uint32_t* atom_result = md_hashmap_get(&atom_only_map, atom_key);
    if (atom_result) {
        md_atomic_number_t z = (md_atomic_number_t)*atom_result;
        // Override CA (Calcium) to Carbon if in amino acid context
        if (z == MD_Z_Ca && res_name.len > 0 && md_util_resname_amino_acid(res_name)) {
            return MD_Z_C;
        }
        return z;
    }
    
    // Try two-letter element heuristic (e.g., CL12 => Cl, BR1 => Br)
    if (strlen(norm_atom_stripped) >= 2) {
        char two_letter[3] = {norm_atom_stripped[0], norm_atom_stripped[1], '\0'};
        str_t two_letter_str = {two_letter, 2};
        md_atomic_number_t two_z = md_atomic_number_from_symbol_icase(two_letter_str);
        if (two_z != MD_Z_X) {
            // Override CA (Calcium) to Carbon if in amino acid context
            if (two_z == MD_Z_Ca && res_name.len > 0 && md_util_resname_amino_acid(res_name)) {
                return MD_Z_C;
            }
            return two_z;
        }
    }
    
    // Final fallback: first-letter element mapping
    char first_letter[2] = {norm_atom_stripped[0], '\0'};
    str_t first_letter_str = {first_letter, 1};
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