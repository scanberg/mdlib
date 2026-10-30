// The names VIAL gives a script before it binds any: the built-in functions, the selectors, the predefined sets, the
// properties that exist only after a dot, and the kinds a parameter may be annotated with. Their signatures and type
// rules arrive with the steps of the checker that use them; this table is what name resolution needs.

#include "builtins.h"

#include <core/md_common.h>

#include <string.h>

#define F   VIAL_BUILTIN_FUNC
#define SEL (VIAL_BUILTIN_FUNC | VIAL_BUILTIN_SELECTOR | VIAL_BUILTIN_VALUE)
#define SET VIAL_BUILTIN_VALUE
#define P   VIAL_BUILTIN_PROPERTY
#define K   VIAL_BUILTIN_KIND

static const vial_builtin_t builtins[] = {
    // Selectors, usable bare for all their units
    {"particles", SEL}, {"atoms", SEL}, {"beads", SEL}, {"residues", SEL}, {"parts", SEL}, {"molecules", SEL},
    {"chains", SEL},
    // Predefined sets
    {"protein", SET}, {"nucleic", SET}, {"water", SET}, {"ion", SET}, {"lipid", SET}, {"backbone", SET},
    {"sidechain", SET}, {"system", SET | F},
    // Geometry and structure
    {"positions", F}, {"mass", F}, {"charge", F}, {"element", F}, {"name", F}, {"bfactor", F}, {"com", F | P},
    {"displacement", F}, {"closest_pair", F}, {"farthest_pair", F}, {"angle", F}, {"dihedral", F}, {"within", F},
    {"fit", F},
    // Reductions
    {"sum", F}, {"mean", F}, {"min", F}, {"max", F}, {"std", F}, {"var", F}, {"any", F}, {"all", F}, {"argmin", F},
    {"argmax", F}, {"len", F}, {"count", F},
    // Math
    {"abs", F}, {"floor", F}, {"ceil", F}, {"round", F}, {"clip", F}, {"sqrt", F}, {"exp", F}, {"log", F},
    {"log10", F}, {"sin", F}, {"cos", F}, {"tan", F}, {"asin", F}, {"acos", F}, {"atan", F}, {"atan2", F},
    {"dot", F}, {"cross", F}, {"norm", F}, {"normalize", F}, {"where", F}, {"unitless", F}, {"float", F | K},
    {"int", F | K},
    // Arrays
    {"diff", F}, {"zip", F}, {"range", F}, {"common", F}, {"labels", F}, {"spread", F}, {"missing_where", F},
    {"valid", F}, {"filled", F}, {"zeros", F}, {"ones", F}, {"full", F}, {"zeros_like", F}, {"ones_like", F},
    {"shape", F}, {"linspace", F}, {"arange", F}, {"eye", F}, {"reshape", F}, {"flatten", F}, {"concatenate", F},
    {"covariance", F}, {"eigh", F}, {"trace", F}, {"det", F}, {"inv", F}, {"solve", F},
    // Sets, patterns, bonds, attributes
    {"smiles", F}, {"smiles_atoms", F}, {"bonds", F}, {"formed", F}, {"broken", F}, {"attr", F},
    // Frames
    {"histogram", F}, {"smooth", F}, {"derivative", F}, {"autocorrelation", F}, {"fold", F}, {"at", F},
    {"time", F}, {"frame", F},
    // Drawing
    {"shown", F}, {"marker", F}, {"line", F}, {"arrow", F}, {"sphere", F}, {"arc", F}, {"disk", F}, {"label", F},
    {"highlight", F},
    // Properties that exist only after a dot
    {"T", P}, {"xyz", P}, {"x", P}, {"y", P}, {"z", P},
    // Kinds, for annotations only
    {"set", K}, {"group", K}, {"point", K}, {"points", K}, {"bool", K},
};

#undef F
#undef SEL
#undef SET
#undef P
#undef K

const vial_builtin_t* vial_builtin_find(str_t name) {
    for (size_t i = 0; i < ARRAY_SIZE(builtins); ++i) {
        if (str_eq(name, str_from_cstr(builtins[i].name))) return &builtins[i];
    }
    return NULL;
}

size_t vial_builtin_count(void) { return ARRAY_SIZE(builtins); }
const vial_builtin_t* vial_builtin_at(size_t i) { return i < ARRAY_SIZE(builtins) ? &builtins[i] : NULL; }

// v1's names, with what to write in VIAL
static const struct { const char* v1; const char* message; const char* fix; } v1_names[] = {
    {"residue",       "`residue` is v1; selectors are plural, and positions are subscripts", "residues"},
    {"chain",         "`chain` is v1; selectors are plural", "chains"},
    {"atom",          "`atom` is v1; write particles", "particles"},
    {"molecule",      "`molecule` is v1; selectors are plural", "molecules"},
    {"distance_min",  "`distance_min` is v1; write distance(a, b, between='closest')", "distance"},
    {"distance_max",  "`distance_max` is v1; write distance(a, b, between='farthest')", "distance"},
    {"distance_pair", "`distance_pair` is v1; every pair is distance(a.xyz[:, None], b.xyz)", "distance"},
    {"coord",         "`coord` is v1; write s.xyz", "positions"},
    {"resname",       "`resname` is v1; it is a keyword, particles(resname=...)", "particles"},
    {"resid",         "`resid` is v1; it is a keyword, particles(resid=...)", "particles"},
};

bool vial_builtin_v1(str_t name, const char** message, const char** fix) {
    for (size_t i = 0; i < ARRAY_SIZE(v1_names); ++i) {
        if (str_eq(name, str_from_cstr(v1_names[i].v1))) {
            if (message) *message = v1_names[i].message;
            if (fix) *fix = v1_names[i].fix;
            return true;
        }
    }
    return false;
}

// Words a script may reach for that VIAL spells differently
static const struct { const char* word; const char* name; } synonyms[] = {
    {"centre", "com"}, {"center", "com"}, {"center_of_mass", "com"}, {"centroid_of", "centroid"},
    {"rg", "gyration"}, {"radius_of_gyration", "gyration"}, {"coords", "positions"}, {"coordinates", "positions"},
    {"pos", "positions"}, {"union", "any"}, {"intersection", "all"}, {"stack", "concatenate"},
};

const char* vial_builtin_synonym(str_t word) {
    for (size_t i = 0; i < ARRAY_SIZE(synonyms); ++i) {
        if (str_eq(word, str_from_cstr(synonyms[i].word))) return synonyms[i].name;
    }
    return NULL;
}
