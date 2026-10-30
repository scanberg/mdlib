// The checker, step by step

#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>

#include <stdio.h>
#include <string.h>

#include <md_vial.h>

struct md_system_t;

// The diagnostics of a script as "CODE line: message | fix", one per line; compiled against sys, which may be NULL
static const char* diags_on(const char* src, const struct md_system_t* sys, md_allocator_i* alloc) {
    const md_vial_compile_desc_t desc = {.sys = sys};
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &desc, alloc);
    const md_vial_diagnostic_t* d = NULL;
    const size_t n = md_vial_diagnostics(prog, &d);
    char* out = NULL;
    for (size_t i = 0; i < n; ++i) {
        uint32_t line = 1;
        for (uint32_t k = 0; k < d[i].span.beg; ++k) line += src[k] == '\n';
        char buf[512];
        int len = snprintf(buf, sizeof(buf), "%s %u: %.*s%s%.*s\n", d[i].code, line, (int)d[i].message.len, d[i].message.ptr,
            d[i].fix.len ? " | " : "", (int)d[i].fix.len, d[i].fix.ptr);
        md_array_push_array(out, buf, (size_t)len, alloc);
    }
    md_array_push(out, '\0', alloc);
    md_vial_program_free(prog);
    return out;
}
static const char* diags_of(const char* src, md_allocator_i* alloc) { return diags_on(src, NULL, alloc); }

struct vial_check { md_allocator_i* alloc; };
UTEST_F_SETUP(vial_check) { utest_fixture->alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4)); }
UTEST_F_TEARDOWN(vial_check) { md_arena_allocator_destroy(utest_fixture->alloc); }

#define EXPECT_DIAGS(src, want) do { \
    const char* got_ = diags_of(src, utest_fixture->alloc); \
    if (strcmp(got_, want) != 0) printf("for:\n%s\ngot:\n%s\nwanted:\n%s\n", src, got_, want); \
    EXPECT_STREQ(want, got_); \
} while (0)

UTEST_F(vial_check, top_level_bindings) {
    EXPECT_DIAGS("n = 1\nm = n + 1\n", "");
    EXPECT_DIAGS("n = 1\nn = 2\n", "E0305 2: `n` is already bound on line 1\n");
    EXPECT_DIAGS("n = 1\nn += 1\n", "E0305 2: `n` is already bound on line 1; at top level a name is bound once\n");
    EXPECT_DIAGS("a, b = 1, 2\nb = 3\n", "E0305 2: `b` is already bound on line 1\n");
    // The branches of a static if each bind the name once
    EXPECT_DIAGS("c = True\nif c:\n    x = 1\nelif not c:\n    x = 2\nelse:\n    x = 3\ny = x\n", "");
    EXPECT_DIAGS("c = True\nif c:\n    x = 1\n    x = 2\n", "E0305 4: `x` is already bound on line 3\n");
    EXPECT_DIAGS("x = 0\nc = True\nif c:\n    x = 1\n", "E0305 4: `x` is already bound on line 1\n");
    // A top-level loop's variable belongs to the loop
    EXPECT_DIAGS("for r in range(3):\n    q = r\nfor r in range(4):\n    q = r\n", "");
    EXPECT_DIAGS("__all__ = ['x']\n", "E0801 1: `__all__` is reserved: a later version may list the published names with it\n");
}

UTEST_F(vial_check, visibility) {
    // A top-level name is seen from the statement after its binding
    EXPECT_DIAGS("y = x\nx = 1\n", "E0201 1: `x` is used before it is bound, on line 2\n");
    // A def body sees every top-level name: it runs when called
    EXPECT_DIAGS("def f():\n    return k\nk = 2\nv = f()\n", "");
    // Hiding a built-in: the lines above still see it
    EXPECT_DIAGS("a = residues\nresidues = 5\nb = residues\n", "W0201 2: `residues` now names a variable; the built-in is hidden from here on\n");
    EXPECT_DIAGS("def distance(a, b):\n    return 0\n", "W0201 1: `distance` now names a def; the built-in is hidden from here on\n");
    // Locals, parameters, nested defs and comprehensions
    EXPECT_DIAGS("def f(a, b=1):\n    c = a + b\n    def g(x):\n        return x + c\n    return g(c)\nv = f(1)\n", "");
    EXPECT_DIAGS("v = [x * y for x in range(3) for y in range(x)]\nw = x\n", "E0201 2: `x` is not defined\n");
    EXPECT_DIAGS("def f():\n    return z\n", "E0201 2: `z` is not defined\n");
}

UTEST_F(vial_check, unknown_names) {
    EXPECT_DIAGS("c = centre(particles)\n", "E0201 1: `centre` is not defined | com\n");
    EXPECT_DIAGS("c = residuez\n", "E0201 1: `residuez` is not defined | residues\n");
    EXPECT_DIAGS("d = 1.2 * nm\n", "E0201 1: `nm` is a unit, not a name: write it after a number, as 10 nm | 1.2 nm\n");
    EXPECT_DIAGS("c = com(residues)\nxy = c.xy\n", "E0201 2: `.xy` is not defined: VIAL has no swizzles; write c[..., 0:2] | c[..., 0:2]\n");
    EXPECT_DIAGS("c = com(residues)\nxz = c.xz\n", "E0201 2: `.xz` is not defined: VIAL has no swizzles; write c[..., [0, 2]] | c[..., [0, 2]]\n");
    EXPECT_DIAGS("c = residues.cmo\n", "E0201 1: `.cmo` is not a function or property | com\n");
    EXPECT_DIAGS("d = distance_min(residues, residues)\n", "E0701 1: `distance_min` is v1; write distance(a, b, between='closest') | distance\n");
    EXPECT_DIAGS("print(1)\n", "E0803 1: `print` is not part of VIAL; hover shows any value\n");
    // A method is any function, built in or a def
    EXPECT_DIAGS("def twice(x):\n    return 2 * x\nv = (3).twice()\nc = residues.com\n", "");
}

UTEST_F(vial_check, annotations) {
    EXPECT_DIAGS("def f(a: set, b: points, c: float = 1.0):\n    return a\n", "");
    EXPECT_DIAGS("def f(a: str):\n    return a\n",
        "E0312 1: `str` is not a parameter kind; the kinds are set, group, point, points, float, int and bool\n");
    EXPECT_DIAGS("def f(a) -> float:\n    return a\n", "E0312 1: a def's result is not annotated; it is inferred at each call\n");
}

// A message quotes at most the first line of source text, cut on a character, so it stays one line
UTEST_F(vial_check, quoted_source) {
    EXPECT_DIAGS("x = 1 '''kJ\nmol'''\n", "E0502 1: `'''kJ…` is not a unit\n");
    EXPECT_DIAGS("x = 2 'abcdefghijklmnopqrstuvwxyz\xC3\xA5\xC3\xA4\xC3\xB6" "abc'\n",
        "E0502 1: `'abcdefghijklmnopqrstuvwxyz\xC3\xA5\xC3\xA4…` is not a unit\n");
    EXPECT_DIAGS("def f(a: [1,\n         2]):\n    return a\n",
        "E0312 1: `[1,…` is not a parameter kind; the kinds are set, group, point, points, float, int and bool\n");
}

UTEST_F(vial_check, recursion) {
    EXPECT_DIAGS("def walk(n):\n    return walk(n - 1)\n", "E0802 2: `walk` calls itself; VIAL has no recursion\n");
    EXPECT_DIAGS("def a(n):\n    return b(n)\ndef b(n):\n    return a(n)\n",
        "E0802 4: `b` calls `a`, which calls `b` again; VIAL has no recursion\n");
    EXPECT_DIAGS("def a(n):\n    return n\ndef b(n):\n    return a(n) + a(1)\nv = b(2)\n", "");
}

#include <vial/vial_internal.h>

// The type the checker gives the statement on a line, as "type, stage = value"; compiled against sys, which may be NULL
static const char* type_on(const char* src, uint32_t line, const struct md_system_t* sys, md_allocator_i* alloc) {
    const md_vial_compile_desc_t desc = {.sys = sys};
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &desc, alloc);
    md_vial__result_t r;
    static char out[512];
    if (!md_vial__line_result(prog, line, &r)) { snprintf(out, sizeof(out), "(none)"); }
    else {
        int n = 0;
        for (int k = 0; k < r.num_types; ++k) n += snprintf(out + n, sizeof(out) - n, "%s%s", k ? " and " : "", r.types[k]);
        n += snprintf(out + n, sizeof(out) - n, ", %s", r.stage);
        if (r.has_value) snprintf(out + n, sizeof(out) - n, " = %s", r.value);
    }
    md_vial_program_free(prog);
    return out;
}
static const char* type_at(const char* src, uint32_t line, md_allocator_i* alloc) { return type_on(src, line, NULL, alloc); }

#define EXPECT_TYPE(src, line, want) do { \
    const char* got_ = type_at(src, line, utest_fixture->alloc); \
    if (strcmp(got_, want) != 0) printf("for:\n%s\nline %d got:  %s\nwanted:      %s\n", src, line, got_, want); \
    EXPECT_STREQ(want, got_); \
} while (0)

UTEST_F(vial_check, literals_and_units) {
    EXPECT_TYPE("big = 1_000_000\n", 1, "int, static = 1000000");
    EXPECT_TYPE("x = 0x1F\n", 1, "int, static = 31");
    EXPECT_TYPE("f = 2.5e-3\n", 1, "float, static = 0.0025");
    EXPECT_TYPE("t = 10 ns\n", 1, "int ns, static = 10");
    EXPECT_TYPE("d = 0.5 nm + 2.0\n", 1, "float nm, static = 2.5");
    EXPECT_TYPE(u8"d = 0.5 nm + 3.0 Å\n", 1, "float nm, static = 0.8");
    EXPECT_TYPE(u8"c = 0.5 nm > 3.0 Å\n", 1, "bool, static = True");
    EXPECT_TYPE(u8"r = 3.0 Å / 2 ps\n", 1, u8"float Å/ps, static = 1.5");
    EXPECT_TYPE(u8"a = (3.0 Å)**2\n", 1, u8"float Å², static = 9");
    EXPECT_TYPE("e = 1.5 'kJ/mol'\n", 1, "float kJ/mol, static = 1.5");
    EXPECT_TYPE("h = 7 // 2\n", 1, "int, static = 3");
    EXPECT_TYPE("q = 7 / 2\n", 1, "float, static = 3.5");
    EXPECT_TYPE("m = -7 % 3\n", 1, "int, static = 2");
    EXPECT_TYPE("s = 'CA' 'B'\n", 1, "str, static = 'CAB'");
    EXPECT_TYPE("b = not True\n", 1, "bool, static = False");
    EXPECT_TYPE("n = 442\nok = 2 < n <= 500\n", 2, "bool, static = True");
    EXPECT_TYPE("k = 5 in 1..10\n", 1, "bool, static");
    EXPECT_TYPE("n = None\n", 1, "None, static = None");
    EXPECT_DIAGS(u8"x = 1.0 Å + 90 deg\n", u8"W0501 1: adds ° to Å; the result's unit becomes unknown | unitless(...)\n");
    EXPECT_TYPE(u8"x = 1.0 Å + 90 deg\ny = x * 2.0\n", 2, "float ?, static");
    EXPECT_DIAGS("x = 'CA' + 1\n", "E0319 1: `+` takes numbers; this is str and int\n");
    EXPECT_DIAGS("x = 2.0 ** 1 nm\n", "W0501 1: a power of nm; an exponent has no unit | unitless(...)\n");
}

UTEST_F(vial_check, shapes) {
    EXPECT_TYPE("ids = [45, 48, 112]\n", 1, "int[3], static");
    EXPECT_TYPE(u8"c = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]\nd = c * 2 Å\n", 2, u8"float[2][3] Å, static");
    EXPECT_TYPE("c = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]\nt = c.T\n", 2, "float[3][2], static");
    EXPECT_TYPE("c = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]\nz = c.z\n", 2, "float[2], static");
    EXPECT_TYPE("c = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]\nr = c[:, None] - c\n", 2, "float[2][2][3], static");
    EXPECT_TYPE("c = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]\nr = c[1:]\n", 2, "float[1][3], static");
    EXPECT_TYPE("c = [1.0, 2.0, 3.0, 4.0, 5.0]\nr = c[::2]\n", 2, "float[3], static");
    EXPECT_TYPE("c = [1.0, 2.0, 3.0, 4.0, 5.0]\nr = c[-3:]\n", 2, "float[3], static");
    EXPECT_TYPE("c = [1.0, 2.0, 3.0, 4.0, 5.0]\nr = c[[3, 0, 1]]\n", 2, "float[3], static");
    EXPECT_TYPE("c = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]\nr = c[..., 2]\n", 2, "float[2], static");
    EXPECT_TYPE("m = [[1.0, 0.0], [0.0, 1.0]]\nv = [1.0, 2.0]\nr = m @ v\n", 3, "float[2], static");
    EXPECT_TYPE("m = [[1.0, 0.0], [0.0, 1.0]]\nr = m @ m\n", 2, "float[2][2], static");
    EXPECT_TYPE("v = [1.0, 2.0]\nr = v @ v\n", 2, "float, static");
    EXPECT_DIAGS("d = [1.0, 2.0]\nc = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]\nbad = d + c\n",
        "E0301 3: `float[2]` and `float[2][3]` do not broadcast | d[:, None]\n");
    EXPECT_DIAGS("c = [1.0, 2.0, 3.0]\nout = c[3]\n", "E0304 2: position 3 is outside `[3]`\n");
    EXPECT_DIAGS("a = [1.0, [2.0]]\n", "E0301 1: the items of a list have one shape; these are `float` and `float[1]`\n");
    EXPECT_DIAGS("c = [1.0, 2.0]\nz = c.z\n", "E0319 2: `.z` takes coordinates, a last dimension of 3, or a set; this is `float[2]`\n");
    EXPECT_DIAGS("a = [[[[[[[1.0]]]]]]]\n", "E0308 1: the result would have 7 dimensions; the limit is 6\n");
}

UTEST_F(vial_check, static_conditions) {
    // The branch a static condition does not take is not type-checked
    EXPECT_TYPE("v = 1.0 if True else 'x' + 1\n", 1, "float, static = 1");
    EXPECT_DIAGS("v = 1.0 if True else 'x' + 1\n", "");
    EXPECT_DIAGS("v = 1.0 if [True, False] else 2.0\n", "E0303 1: a conditional needs one bool; this is `bool[2]` | where(...)\n");
    EXPECT_TYPE("c = 3\nv = 1.0 if c > 2 else 2\n", 2, "float, static = 1");
}

// ### STEP 4: THE SYSTEM AND THE SELECTORS ###

#include "vial_test_systems.h"

#define EXPECT_TYPE_ON(sys, src, line, want) do { \
    const char* got_ = type_on(src, line, sys, utest_fixture->alloc); \
    if (strcmp(got_, want) != 0) printf("for:\n%s\nline %d got:  %s\nwanted:      %s\n", src, line, got_, want); \
    EXPECT_STREQ(want, got_); \
} while (0)

#define EXPECT_DIAGS_ON(sys, src, want) do { \
    const char* got_ = diags_on(src, sys, utest_fixture->alloc); \
    if (strcmp(got_, want) != 0) printf("for:\n%s\ngot:\n%s\nwanted:\n%s\n", src, got_, want); \
    EXPECT_STREQ(want, got_); \
} while (0)

UTEST_F(vial_check, selectors) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE_ON(T, "n = len(residues)\n", 1, "int, static = 880");
    EXPECT_TYPE_ON(T, "r = residues(chain='A')\n", 1, "set[442], static");
    EXPECT_TYPE_ON(T, "r = residues(1..4, chain='A')\n", 1, "set[4], static");
    EXPECT_TYPE_ON(T, "n = len(residues([1, 3, 5], chain='B'))\n", 1, "int, static = 3");
    EXPECT_TYPE_ON(T, "n = count(particles(name='C*', resname='GTP'))\n", 1, "int, static = 10");
    EXPECT_TYPE_ON(T, "n = count(particles(resname=['GTP', 'GDP']))\n", 1, "int, static = 60");
    EXPECT_TYPE_ON(T, "n = count(particles(element='MG'))\n", 1, "int, static = 1");    // elements ignore case
    EXPECT_TYPE_ON(T, "n = count(particles(element=['Mg', 's']))\n", 1, "int, static = 52");     // 51 S, 1 Mg
    EXPECT_TYPE_ON(T, "m = particles(name='MG')\n", 1, "set, static");
    EXPECT_TYPE_ON(T, "n = len(residues(scope=particles(chain='B')))\n", 1, "int, static = 437");
    EXPECT_TYPE_ON(T, "n = len(parts(kind='water'))\n", 1, "int, static = 8");
    EXPECT_TYPE_ON(T, "n = len(parts(kind='ion'))\n", 1, "int, static = 1");
    EXPECT_TYPE_ON(T, "n = len(molecules('SO*'))\n", 1, "int, static = 3");
    EXPECT_TYPE_ON(T, "n = len(residues(particles(resname='GTP')))\n", 1, "int, static = 1");
    EXPECT_TYPE_ON(T, "n = count(particles(residues(1..2, chain='A')))\n", 1, "int, static = 19");     // MET 8, ARG 11
    EXPECT_TYPE_ON(T, "n = count(atoms(resname='VLB'))\n", 1, "int, static = 59");
    EXPECT_TYPE_ON(T, "n = count(beads)\n", 1, "int, static = 0");
}

UTEST_F(vial_check, predefined_sets) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE_ON(T, "n = count(protein)\n", 1, "int, static = 6784");
    EXPECT_TYPE_ON(T, "n = count(water)\n", 1, "int, static = 8");
    EXPECT_TYPE_ON(T, "n = count(ion)\n", 1, "int, static = 1");
    EXPECT_TYPE_ON(T, "n = count(system)\n", 1, "int, static = 6927");
    EXPECT_TYPE_ON(T, "n = count(protein and not backbone and not sidechain)\n", 1, "int, static = 865");
    EXPECT_TYPE_ON(T, "n = count(protein or water)\n", 1, "int, static = 6792");
    EXPECT_TYPE_ON(T, "n = count(not protein)\n", 1, "int, static = 143");
    EXPECT_TYPE_ON(T, "n = count(backbone and particles(chain='A'))\n", 1, "int, static = 1299");
}

UTEST_F(vial_check, set_indexing_and_labels) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE_ON(T, "s = particles[0]\n", 1, "set, static");
    EXPECT_TYPE_ON(T, "n = count(particles(resname='GTP')[0:3])\n", 1, "int, static = 3");
    EXPECT_TYPE_ON(T, "n = count(particles(resname='GTP')[[0, 5, 5]])\n", 1, "int, static = 2");
    EXPECT_TYPE_ON(T, "r = residues[-3:]\n", 1, "set[3], static");
    EXPECT_TYPE_ON(T, "n = count(residues[0])\n", 1, "int, static = 8");
    EXPECT_TYPE_ON(T, "c = count(residues(1..3, chain='A'))\n", 1, "int[3], static");
    EXPECT_TYPE_ON(T, "l = labels(chains)[1]\n", 1, "str, static = 'B'");
    EXPECT_TYPE_ON(T, "l = labels(parts)[2]\n", 1, "str, static = 'GTP501/A'");
    EXPECT_TYPE_ON(T, "l = labels(parts)[0]\n", 1, "str, static = 'peptide/A'");
    EXPECT_TYPE_ON(T, "l = labels(residues(chain='B'))[-1]\n", 1, "str, static = 'HOH603/B'");
    EXPECT_DIAGS_ON(T, "l = labels(particles(name='MG'))\n", "E0316 1: `labels` takes a dimension of set, which has 0\n");
    EXPECT_DIAGS_ON(T, "s = particles(resname='GTP')[32]\n", "E0304 1: position 32 is outside a set of 32 particles\n");
    EXPECT_DIAGS_ON(T, "s = residues[880]\n", "E0304 1: position 880 is outside `[880]`\n");
}

UTEST_F(vial_check, selector_diagnostics) {
    const md_system_t* T = vial_test_system('T');
    const md_system_t* P = vial_test_system('P');
    ASSERT_TRUE(T != NULL && P != NULL);
    EXPECT_DIAGS_ON(T, "b = chains(name='A')\n", "E0202 1: `chains` takes no keyword `name`; it takes id and scope | id='A'\n");
    EXPECT_DIAGS_ON(T, "b = residues(resname='ALA')\n", "E0202 1: `residues` takes no keyword `resname`; it takes id, name, kind, chain and scope | name='ALA'\n");
    EXPECT_DIAGS_ON(T, "t = residues(45, 'LYS')\n", "E0203 1: `residues` takes one positional argument; `'LYS'` is the second | name='LYS'\n");
    EXPECT_DIAGS_ON(T, "n = chains(0)\n", "E0205 1: `chains` takes an id string, not int; positions are subscripts | chains[0]\n");
    EXPECT_DIAGS_ON(T, "p = parts('GTP')\n", "E0205 1: `parts` takes a set; parts have no name or id, not str\n");
    EXPECT_DIAGS_ON(T, "r = residues(kind='amino')\n",
        "E0206 1: 'amino' is not a residue kind; the kinds are amino_acid, nucleotide, water, ion, lipid and other | 'amino_acid'\n");
    EXPECT_DIAGS_ON(T, "r = residues(kind='')\n",
        "E0206 1: '' is not a residue kind; the kinds are amino_acid, nucleotide, water, ion, lipid and other\n");
    EXPECT_DIAGS_ON(T, "r = residues(id='12')\n", "E0206 1: `id` takes an int, a range a..b, or a list of ints; this is str\n");
    EXPECT_DIAGS_ON(T, "r = residues(39..45, chain='A')\n", "W0601 1: `residues(39..45, chain='A')` matches nothing in its scope\n");
    EXPECT_DIAGS_ON(T, "p = particles(resid=501)\n",
        "W0602 1: `resid=501` matches residues in chains A and B; the merged set hides which | resid=501, chain='A'\n");
    EXPECT_DIAGS_ON(T, "p = particles(resid=501, chain='B')\nq = particles(resid=1..2)\n",
        "W0602 2: `resid=1..2` matches residues in chains A and B; the merged set hides which | resid=1..2, chain='A'\n");
    EXPECT_DIAGS_ON(P, "c = chains('A')\n",
        "E0207 1: the data gives no chain ids, and chains are never inferred; a strand is found as a part | parts(kind='polymer')\n");
    EXPECT_DIAGS_ON(P, "c = particles(chain='A')\n",
        "E0207 1: the data gives no chain ids, and chains are never inferred; a strand is found as a part | parts(kind='polymer')\n");
    EXPECT_DIAGS_ON(P, "s = parts(kind='peptide')\nn = len(molecules)\n", "");
    EXPECT_TYPE_ON(P, "l = labels(residues)[0]\n", 1, "str, static = 'ALA1'");
    EXPECT_TYPE_ON(P, "l = labels(parts)[0]\n", 1, "str, static = 'peptide1'");
}

UTEST_F(vial_check, selectors_without_a_system) {
    EXPECT_TYPE("r = residues(1..4)\n", 1, "set[n], static");
    EXPECT_TYPE("p = particles(name='CA')\n", 1, "set, static");
    EXPECT_TYPE("n = len(residues)\n", 1, "int, static");
    EXPECT_TYPE("n = count(protein)\n", 1, "int, static");
    EXPECT_TYPE("w = within(5.0, of=protein)\n", 1, "set, frame");
    EXPECT_TYPE("r = residues(within(5.0, of=protein))\n", 1, "set[?], frame");
    EXPECT_DIAGS("w = within(5.0, protein)\n", "E0203 1: `within` takes the radius as its one positional argument; of= is always written, of=protein | of=protein\n");
    EXPECT_DIAGS("w = within(5.0 ns, of=protein)\n", "E0313 1: `within` takes a radius, one length; this is float ns\n");
}

// ### STEP 5: KINDS AND GEOMETRY ###

UTEST_F(vial_check, positions) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE_ON(T, "p = particles(resname='GTP').xyz\n", 1, u8"float[32][3] Å, frame");
    EXPECT_TYPE_ON(T, "z = protein.z\n", 1, u8"float[6784] Å, frame");
    EXPECT_TYPE_ON(T, "w = residues(kind='water').xyz\n", 1, u8"float[8][3] Å, frame");
    EXPECT_TYPE_ON(T, "p = positions(particles(name='MG'))\n", 1, u8"float[1][3] Å, frame");
    EXPECT_TYPE_ON(T, "n = residues(within(5.0, of=protein)).xyz\n", 1, "error, static");
    EXPECT_TYPE_ON(T, "n = within(5.0, of=protein).xyz\n", 1, u8"float[?][3] Å, frame");
    EXPECT_DIAGS_ON(T, "r = residues(chain='A', kind='amino_acid').xyz\n",
        u8"E0318 1: `residues(chain='A', kind='amino_…` holds 433 sets of 4 to 14 particles; their positions do not form one array"
        " | particles(residues(chain='A', kind='amino_acid')).xyz\n");
    EXPECT_TYPE_ON(T, "m = mass(particles(resname='GTP'))\n", 1, "float[32] u, static");
    EXPECT_TYPE_ON(T, "e = element(particles(resname='GTP'))\n", 1, "str[32], static");
    EXPECT_TYPE_ON(T, "q = charge(protein)\n", 1, "float[6784] e, static");
}

UTEST_F(vial_check, centres_and_points) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE_ON(T, "c = particles(resname='GTP').com\n", 1, u8"float[3] Å, frame");
    EXPECT_TYPE_ON(T, "c = com(residues(chain='A', kind='amino_acid'))\n", 1, u8"float[433][3] Å, frame");
    EXPECT_TYPE_ON(T, "c = com(residues(chain='A'), weights=None)\n", 1, u8"float[442][3] Å, frame");
    EXPECT_TYPE_ON(T, "c = [1.0, 2.0, 3.0].com\n", 1, u8"float[3] Å, static");
    EXPECT_DIAGS_ON(T, "c = com(protein, weights='charge')\n", "E0313 1: `weights` takes 'mass' or None\n");
    EXPECT_TYPE_ON(T, "d = displacement(particles(name='MG'), particles(resname='GTP').com)\n", 1, u8"float[3] Å, frame");
    EXPECT_TYPE_ON(T, "d = displacement(residues(kind='water'), particles(name='MG'))\n", 1, u8"float[8][3] Å, frame");
    EXPECT_DIAGS_ON(T, "d = displacement(particles(resname='GTP'), particles(name='MG'))\n",
        "E0302 1: `displacement` needs a point for a; `particles(resname='GTP')` may hold several particles, so write .com for its centre"
        " or between='closest' | particles(resname='GTP').com\n");
    EXPECT_DIAGS_ON(T, "d = displacement(12, particles(name='MG'))\n",
        "E0302 1: `displacement` needs a point for a; an int is not one | particles(serial=12)\n");
    EXPECT_DIAGS_ON(T, "d = displacement(residues(1..4, chain='A').com, residues(1..3, chain='B').com)\n",
        "E0301 1: the dimensions `displacement` maps over do not broadcast: `int[4]` and `int[3]`\n");
}

UTEST_F(vial_check, angles_and_pairs) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE_ON(T, "a = angle(particles(resname='GTP').com, particles(resname='GDP').com, particles(name='MG'))\n", 1, u8"float °, frame");
    EXPECT_TYPE_ON(T, "a = angle([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0])\n", 1, u8"float °, static");
    EXPECT_TYPE_ON(T, "c = com(residues(1..8, chain='A'))\nt = dihedral(c[0:4])\n", 2, u8"float °, frame");
    EXPECT_TYPE_ON(T, "c = com(residues(1..8, chain='A'))\nt = dihedral(c[0], c[1], c[2], c[3:5])\n", 2, u8"float[2] °, frame");
    EXPECT_DIAGS_ON(T, "a = angle(particles(name='MG'), particles(name='MG'))\n", "E0204 1: `angle` takes 3 points; it got 2\n");
    EXPECT_TYPE_ON(T, "d, i, j = closest_pair(particles(resname='GTP'), particles(resname='VLB'))\n", 1, u8"float Å and int and int, frame");
    EXPECT_TYPE_ON(T, "c = com(residues(chain='A', kind='amino_acid'))\nd, i, j = closest_pair(c[:, None], c, distinct=True)\n", 2,
        u8"float[433] Å and int[433] and int[433], frame");
    EXPECT_DIAGS_ON(T, "c = com(residues(chain='A', kind='amino_acid'))\nd, i, j = closest_pair(c[:, None], c)\n",
        "W0608 2: `closest_pair` pairs every point with itself, so the shortest distance is zero; distinct=True leaves those pairs out"
        " | closest_pair(c[:, None], c, distinct=True)\n");
    EXPECT_DIAGS_ON(T, "d, i, j = farthest_pair(protein, protein)\n", "");
    EXPECT_TYPE_ON(T, "f = fit(particles(resname='GTP'), particles(resname='GTP').xyz)\n", 1, u8"float[32][3] Å, frame");
}
