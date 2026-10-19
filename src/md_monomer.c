#include <md_types.h>

#include <core/md_common.h>
#include <core/md_str.h>

// The predefined table of standard monomers (md_monomer_t), and the names residues of these go by in files.
// What md_util_element_* is to the atom types, this is to the component types.

#define MONOMER_RGB(r,g,b) ((0xFFu << 24) | (((uint32_t)(b) & 255) << 16) | (((uint32_t)(g) & 255) << 8) | ((uint32_t)(r) & 255))

typedef struct monomer_info_t {
    str_t id;
    str_t name;
    char  letter;
    md_component_kind_t kind;
    uint32_t color;
} monomer_info_t;

#define AA MD_COMPONENT_KIND_AMINO_ACID
#define NT MD_COMPONENT_KIND_NUCLEOTIDE

// Default colors:
//  - Amino acids: a palette made for viamd in OKLCH. The hue tells the physicochemical family (acidic red, basic blue,
//    hydroxyl orange, amide teal, aliphatic green, sulfur yellow, aromatic violet, conformationally special neutral),
//    the lightness tells the members of a family apart (the aliphatics darken with their size), and the chroma stays
//    moderate. Tuned so that the smallest difference between any two (in OKLab) is about 0.08, and 0.06 as the Linear
//    color input of viamd shows them or with deuteranopia. The classic schemes (shapely, amino, Taylor) are of fully
//    saturated RGB: shapely colors PHE and PRO nearly alike (0.02). SEC and PYL take shades of the families of the
//    amino acids they are analogues of (CYS, LYS), UNK a neutral grey, as it is the color of every modified residue and
//    cap of a chain.
//  - Nucleotides: the colors the sequence view of viamd used, one hue per base, the same for DNA and RNA. Shapely
//    colors G and U nearly alike (FF7070, FF8080). Inosine takes the shapely color, N the neutral grey of UNK.
static const monomer_info_t monomer_table[MD_MONOMER_COUNT] = {
    [MD_MONOMER_NONE] = { STR_INIT(""),    STR_INIT(""),                      0,   MD_COMPONENT_KIND_OTHER, 0 },

    [MD_MONOMER_ALA] = { STR_INIT("ALA"), STR_INIT("Alanine"),               'A', AA, MONOMER_RGB(186, 231, 165) },
    [MD_MONOMER_ARG] = { STR_INIT("ARG"), STR_INIT("Arginine"),              'R', AA, MONOMER_RGB( 46,  82, 136) },
    [MD_MONOMER_ASN] = { STR_INIT("ASN"), STR_INIT("Asparagine"),            'N', AA, MONOMER_RGB( 94, 218, 217) },
    [MD_MONOMER_ASP] = { STR_INIT("ASP"), STR_INIT("Aspartic acid"),         'D', AA, MONOMER_RGB(191,  86,  56) },
    [MD_MONOMER_CYS] = { STR_INIT("CYS"), STR_INIT("Cysteine"),              'C', AA, MONOMER_RGB(242, 222,  78) },
    [MD_MONOMER_GLN] = { STR_INIT("GLN"), STR_INIT("Glutamine"),             'Q', AA, MONOMER_RGB( 56, 175, 167) },
    [MD_MONOMER_GLU] = { STR_INIT("GLU"), STR_INIT("Glutamic acid"),         'E', AA, MONOMER_RGB(162,  48,  50) },
    [MD_MONOMER_GLY] = { STR_INIT("GLY"), STR_INIT("Glycine"),               'G', AA, MONOMER_RGB(245, 245, 245) },
    [MD_MONOMER_HIS] = { STR_INIT("HIS"), STR_INIT("Histidine"),             'H', AA, MONOMER_RGB(114, 177, 246) },
    [MD_MONOMER_ILE] = { STR_INIT("ILE"), STR_INIT("Isoleucine"),            'I', AA, MONOMER_RGB( 30, 135,  86) },
    [MD_MONOMER_LEU] = { STR_INIT("LEU"), STR_INIT("Leucine"),               'L', AA, MONOMER_RGB( 94, 161,  93) },
    [MD_MONOMER_LYS] = { STR_INIT("LYS"), STR_INIT("Lysine"),                'K', AA, MONOMER_RGB( 90, 137, 216) },
    [MD_MONOMER_MET] = { STR_INIT("MET"), STR_INIT("Methionine"),            'M', AA, MONOMER_RGB(189, 173,  62) },
    [MD_MONOMER_PHE] = { STR_INIT("PHE"), STR_INIT("Phenylalanine"),         'F', AA, MONOMER_RGB(150, 110, 175) },
    [MD_MONOMER_PRO] = { STR_INIT("PRO"), STR_INIT("Proline"),               'P', AA, MONOMER_RGB(138, 133, 128) },
    [MD_MONOMER_SER] = { STR_INIT("SER"), STR_INIT("Serine"),                'S', AA, MONOMER_RGB(251, 184, 113) },
    [MD_MONOMER_THR] = { STR_INIT("THR"), STR_INIT("Threonine"),             'T', AA, MONOMER_RGB(204, 127,  51) },
    [MD_MONOMER_TRP] = { STR_INIT("TRP"), STR_INIT("Tryptophan"),            'W', AA, MONOMER_RGB(104,  87, 167) },
    [MD_MONOMER_TYR] = { STR_INIT("TYR"), STR_INIT("Tyrosine"),              'Y', AA, MONOMER_RGB(198, 152, 208) },
    [MD_MONOMER_VAL] = { STR_INIT("VAL"), STR_INIT("Valine"),                'V', AA, MONOMER_RGB(145, 183, 135) },
    [MD_MONOMER_SEC] = { STR_INIT("SEC"), STR_INIT("Selenocysteine"),        'U', AA, MONOMER_RGB(239, 188,  58) },
    [MD_MONOMER_PYL] = { STR_INIT("PYL"), STR_INIT("Pyrrolysine"),           'O', AA, MONOMER_RGB( 89, 105, 174) },
    [MD_MONOMER_UNK] = { STR_INIT("UNK"), STR_INIT("Unknown amino acid"),    'X', AA, MONOMER_RGB(187, 183, 175) },

    [MD_MONOMER_DA]  = { STR_INIT("DA"),  STR_INIT("Deoxyadenosine"),        'A', NT, MONOMER_RGB(252, 105, 122) },
    [MD_MONOMER_DC]  = { STR_INIT("DC"),  STR_INIT("Deoxycytidine"),         'C', NT, MONOMER_RGB(248, 238,  92) },
    [MD_MONOMER_DG]  = { STR_INIT("DG"),  STR_INIT("Deoxyguanosine"),        'G', NT, MONOMER_RGB(213, 179, 239) },
    [MD_MONOMER_DT]  = { STR_INIT("DT"),  STR_INIT("Thymidine"),             'T', NT, MONOMER_RGB(159, 243, 160) },
    [MD_MONOMER_DU]  = { STR_INIT("DU"),  STR_INIT("Deoxyuridine"),          'U', NT, MONOMER_RGB(254, 157,  45) },
    [MD_MONOMER_DI]  = { STR_INIT("DI"),  STR_INIT("Deoxyinosine"),          'I', NT, MONOMER_RGB(128, 255, 255) },
    [MD_MONOMER_A]   = { STR_INIT("A"),   STR_INIT("Adenosine"),             'A', NT, MONOMER_RGB(252, 105, 122) },
    [MD_MONOMER_C]   = { STR_INIT("C"),   STR_INIT("Cytidine"),              'C', NT, MONOMER_RGB(248, 238,  92) },
    [MD_MONOMER_G]   = { STR_INIT("G"),   STR_INIT("Guanosine"),             'G', NT, MONOMER_RGB(213, 179, 239) },
    [MD_MONOMER_U]   = { STR_INIT("U"),   STR_INIT("Uridine"),               'U', NT, MONOMER_RGB(254, 157,  45) },
    [MD_MONOMER_I]   = { STR_INIT("I"),   STR_INIT("Inosine"),               'I', NT, MONOMER_RGB(128, 255, 255) },
    [MD_MONOMER_N]   = { STR_INIT("N"),   STR_INIT("Unknown nucleotide"),    'N', NT, MONOMER_RGB(187, 183, 175) },
};

#undef AA
#undef NT

typedef struct monomer_alias_t {
    const char*  name;
    md_monomer_t monomer;
} monomer_alias_t;

// Names of forms of the amino acids, besides their ids
static const monomer_alias_t amino_acid_aliases[] = {
    // Histidine: Amber (HID, HIE, HIP), CHARMM (HSD, HSE, HSP), GROMACS (HISD, HISE, HISH, HIS1, HIS2), GROMOS (HISA,
    // HISB), Martini (HISH, HIH)
    {"HID",  MD_MONOMER_HIS}, {"HIE",  MD_MONOMER_HIS}, {"HIP",  MD_MONOMER_HIS},
    {"HSD",  MD_MONOMER_HIS}, {"HSE",  MD_MONOMER_HIS}, {"HSP",  MD_MONOMER_HIS},
    {"HISD", MD_MONOMER_HIS}, {"HISE", MD_MONOMER_HIS}, {"HISH", MD_MONOMER_HIS}, {"HISA", MD_MONOMER_HIS},
    {"HISB", MD_MONOMER_HIS}, {"HIS1", MD_MONOMER_HIS}, {"HIS2", MD_MONOMER_HIS}, {"HIH",  MD_MONOMER_HIS},
    // Cysteine: disulfide bonded and deprotonated forms, the lipidated ones of Martini 3, and modified residues
    {"CYX",  MD_MONOMER_CYS}, {"CYM",  MD_MONOMER_CYS}, {"CYS1", MD_MONOMER_CYS}, {"CYS2", MD_MONOMER_CYS},
    {"CYSH", MD_MONOMER_CYS}, {"CYSF", MD_MONOMER_CYS}, {"CYSG", MD_MONOMER_CYS}, {"CYSP", MD_MONOMER_CYS},
    {"CSO",  MD_MONOMER_CYS}, {"CSD",  MD_MONOMER_CYS}, {"CME",  MD_MONOMER_CYS}, {"OCS",  MD_MONOMER_CYS},
    // Protonated aspartic and glutamic acid
    {"ASH",  MD_MONOMER_ASP}, {"ASPP", MD_MONOMER_ASP}, {"ASPH", MD_MONOMER_ASP}, {"ASP0", MD_MONOMER_ASP},
    {"GLH",  MD_MONOMER_GLU}, {"GLUP", MD_MONOMER_GLU}, {"GLUH", MD_MONOMER_GLU}, {"GLU0", MD_MONOMER_GLU},
    // Neutral lysine and arginine, and modified lysines
    {"LYN",  MD_MONOMER_LYS}, {"LSN",  MD_MONOMER_LYS}, {"LYSN", MD_MONOMER_LYS}, {"LYSH", MD_MONOMER_LYS},
    {"LYS0", MD_MONOMER_LYS}, {"KCX",  MD_MONOMER_LYS}, {"MLY",  MD_MONOMER_LYS}, {"M3L",  MD_MONOMER_LYS},
    {"ALY",  MD_MONOMER_LYS}, {"LLP",  MD_MONOMER_LYS},
    {"ARG0", MD_MONOMER_ARG}, {"ARGN", MD_MONOMER_ARG},
    // Modified residues of the PDB
    {"MSE",  MD_MONOMER_MET}, {"SEP",  MD_MONOMER_SER}, {"TPO",  MD_MONOMER_THR}, {"PTR",  MD_MONOMER_TYR},
    {"HYP",  MD_MONOMER_PRO},
    // Myristoylated glycine of Martini 3
    {"GLYM", MD_MONOMER_GLY},
};

// The caps of peptide chains, which are not amino acids but are part of their chains
static const char* peptide_caps[] = { "ACE", "NME", "NHE", "NH2", "FOR" };

// Names of forms of the nucleotides, besides their ids. The 5' and 3' termini of Amber (DA5, RA3, DAN, ...) are not
// listed: they are found by their suffix.
static const monomer_alias_t nucleotide_aliases[] = {
    // CHARMM, which uses the same names for DNA and RNA
    {"ADE", MD_MONOMER_A}, {"CYT", MD_MONOMER_C}, {"GUA", MD_MONOMER_G}, {"URA", MD_MONOMER_U}, {"THY", MD_MONOMER_DT},
    // Amber RNA
    {"RA",  MD_MONOMER_A}, {"RC",  MD_MONOMER_C}, {"RG",  MD_MONOMER_G}, {"RU",  MD_MONOMER_U},
    // Thymidine of older PDB files, and the unknown deoxynucleotide
    {"T",   MD_MONOMER_DT}, {"DN",  MD_MONOMER_N},
    // Modified residues of the PDB
    {"PSU", MD_MONOMER_U}, {"H2U", MD_MONOMER_U}, {"5MU", MD_MONOMER_U}, {"4SU", MD_MONOMER_U}, {"OMU", MD_MONOMER_U},
    {"5BU", MD_MONOMER_U}, {"5MC", MD_MONOMER_C}, {"OMC", MD_MONOMER_C}, {"5CM", MD_MONOMER_DC},
    {"1MA", MD_MONOMER_A}, {"1MG", MD_MONOMER_G}, {"2MG", MD_MONOMER_G}, {"7MG", MD_MONOMER_G}, {"M2G", MD_MONOMER_G},
    {"OMG", MD_MONOMER_G}, {"YG",  MD_MONOMER_G},
};

static inline bool valid_monomer(md_monomer_t monomer) {
    return MD_MONOMER_NONE < monomer && monomer < MD_MONOMER_COUNT;
}

str_t md_monomer_id(md_monomer_t monomer) {
    return valid_monomer(monomer) ? monomer_table[monomer].id : STR_LIT("");
}

str_t md_monomer_name(md_monomer_t monomer) {
    return valid_monomer(monomer) ? monomer_table[monomer].name : STR_LIT("");
}

char md_monomer_letter(md_monomer_t monomer) {
    return valid_monomer(monomer) ? monomer_table[monomer].letter : 0;
}

md_component_kind_t md_monomer_kind(md_monomer_t monomer) {
    return valid_monomer(monomer) ? monomer_table[monomer].kind : MD_COMPONENT_KIND_OTHER;
}

uint32_t md_monomer_color(md_monomer_t monomer) {
    return valid_monomer(monomer) ? monomer_table[monomer].color : 0;
}

md_monomer_t md_monomer_from_id(str_t id) {
    id = str_trim(id);
    if (str_empty(id)) return MD_MONOMER_NONE;
    for (int m = MD_MONOMER_NONE + 1; m < MD_MONOMER_COUNT; ++m) {
        if (str_eq_ignore_case(id, monomer_table[m].id)) return (md_monomer_t)m;
    }
    return MD_MONOMER_NONE;
}

// The monomer of the given kind whose id or alias the name is
static md_monomer_t match_name(str_t name, md_component_kind_t kind) {
    if (str_empty(name)) return MD_MONOMER_NONE;
    for (int m = MD_MONOMER_NONE + 1; m < MD_MONOMER_COUNT; ++m) {
        if (monomer_table[m].kind == kind && str_eq_ignore_case(name, monomer_table[m].id)) return (md_monomer_t)m;
    }
    const monomer_alias_t* aliases = kind == MD_COMPONENT_KIND_AMINO_ACID ? amino_acid_aliases : nucleotide_aliases;
    const size_t num_aliases = kind == MD_COMPONENT_KIND_AMINO_ACID ? ARRAY_SIZE(amino_acid_aliases) : ARRAY_SIZE(nucleotide_aliases);
    for (size_t i = 0; i < num_aliases; ++i) {
        if (str_eq_cstr_ignore_case(name, aliases[i].name)) return aliases[i].monomer;
    }
    return MD_MONOMER_NONE;
}

static inline bool is_peptide_cap(str_t name) {
    for (size_t i = 0; i < ARRAY_SIZE(peptide_caps); ++i) {
        if (str_eq_cstr_ignore_case(name, peptide_caps[i])) return true;
    }
    return false;
}

static md_monomer_t match_amino_acid(str_t name) {
    md_monomer_t m = match_name(name, MD_COMPONENT_KIND_AMINO_ACID);
    // The N- and C-terminal residues of Amber: NALA, CALA, NHIE, ...
    if (!m && name.len == 4) {
        const int c = to_upper(name.ptr[0]);
        if (c == 'N' || c == 'C') {
            m = match_name(str_substr(name, 1, SIZE_MAX), MD_COMPONENT_KIND_AMINO_ACID);
        }
    }
    return m;
}

static md_monomer_t match_nucleotide(str_t name) {
    md_monomer_t m = match_name(name, MD_COMPONENT_KIND_NUCLEOTIDE);
    // The 5' and 3' terminal residues and the free nucleosides of Amber: DA5, DA3, DAN, RA5, A3, ...
    if (!m && name.len >= 2) {
        const int c = to_upper(name.ptr[name.len - 1]);
        if (c == '5' || c == '3' || c == 'N') {
            m = match_name(str_substr(name, 0, name.len - 1), MD_COMPONENT_KIND_NUCLEOTIDE);
        }
    }
    return m;
}

md_monomer_t md_monomer_from_comp_name(str_t comp_name, md_component_kind_t kind) {
    const str_t name = str_trim(comp_name);
    switch (kind) {
    case MD_COMPONENT_KIND_AMINO_ACID: {
        const md_monomer_t m = match_amino_acid(name);
        return m ? m : MD_MONOMER_UNK;
    }
    case MD_COMPONENT_KIND_NUCLEOTIDE: {
        const md_monomer_t m = match_nucleotide(name);
        return m ? m : MD_MONOMER_N;
    }
    case MD_COMPONENT_KIND_OTHER:
        return is_peptide_cap(name) ? MD_MONOMER_UNK : MD_MONOMER_NONE;
    default:
        return MD_MONOMER_NONE;
    }
}
