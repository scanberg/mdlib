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
    EXPECT_TYPE(u8"a = (3.0 Å)**2\n", 1, u8"float Å², static = 9.0");
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

// Folding never overflows: a result beyond 32 bits, or a float beyond an int, is typed and not folded
UTEST_F(vial_check, folding_limits) {
    EXPECT_TYPE("x = 2147483647 + 1\n", 1, "int, static");
    EXPECT_TYPE("x = 99999999999999999999 * 3\n", 1, "int, static");
    EXPECT_TYPE("x = 2 ** 40\n", 1, "int, static");
    EXPECT_TYPE("x = 2 ** 30\n", 1, "int, static = 1073741824");
    EXPECT_TYPE("x = int(1e30)\n", 1, "int, static");
    EXPECT_TYPE("x = 1e300 * 1e300\n", 1, "float, static");
    EXPECT_TYPE("x = -2147483647 - 1\n", 1, "int, static = -2147483648");
    EXPECT_TYPE("x = 1 << 40\n", 1, "int, static");
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
    EXPECT_TYPE("v = 1.0 if True else 'x' + 1\n", 1, "float, static = 1.0");
    EXPECT_DIAGS("v = 1.0 if True else 'x' + 1\n", "");
    EXPECT_DIAGS("v = 1.0 if [True, False] else 2.0\n", "E0303 1: a conditional needs one bool; this is `bool[2]` | where(...)\n");
    EXPECT_TYPE("c = 3\nv = 1.0 if c > 2 else 2\n", 2, "float, static = 1.0");
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

// ### STEP 6: REDUCTIONS, MATH AND ARRAYS ###

UTEST_F(vial_check, reductions) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE_ON(T, "c = com(residues(chain='A'))\ns = sum(c)\n", 2, u8"float Å, frame");
    EXPECT_TYPE_ON(T, "c = com(residues(chain='A'))\ns = mean(c, axis=0)\n", 2, u8"float[3] Å, frame");
    EXPECT_TYPE_ON(T, "c = com(residues(chain='A'))\ns = var(c, axis=-1)\n", 2, u8"float[442] Å², frame");
    EXPECT_TYPE_ON(T, "c = com(residues(chain='A'))\ns = mean(c, axis='frame')\n", 2, u8"float[442][3] Å, trajectory");
    EXPECT_TYPE_ON(T, "c = com(residues(chain='A'))\ns = std(c, axis=(0, 'frame'))\n", 2, u8"float[3] Å, trajectory");
    EXPECT_TYPE_ON(T, "c = com(residues(chain='A'))\ni = argmin(c, axis=0)\n", 2, "int[3], frame");
    EXPECT_TYPE_ON(T, "z = protein.z\nn = sum(z > 10.0)\n", 2, "int, frame");
    EXPECT_TYPE_ON(T, "z = protein.z\nm = min(z, 5.0)\n", 2, u8"float[6784] Å, frame");
    EXPECT_TYPE_ON(T, "b = any(residues(501))\n", 1, "set, static");
    EXPECT_TYPE_ON(T, "n = count(any(residues(501)))\n", 1, "int, static = 60");
    EXPECT_TYPE_ON(T, "n = count(all(residues(chain='A')[0:2]))\n", 1, "int, static = 0");
    EXPECT_DIAGS_ON(T, "c = com(residues(chain='A'))\ns = sum(c, axis=2)\n",
        u8"E0316 2: axis=2 on `float[442][3] Å`, which has 2 dimensions | axis=1\n");
    EXPECT_DIAGS_ON(T, "s = sum(5.0, axis='frame')\n", "E0316 1: axis='frame' on a static value, which has no frames\n");
    EXPECT_DIAGS_ON(T, "z = protein.z\nm = mean(z, by=residues)\n",
        "W0405 2: the groups with no values are missing\nW0603 2: 15 of 880 groups have no values in `z`\n");
    EXPECT_TYPE_ON(T, "z = protein.z\nm = mean(z, by=residues)\n", 2, u8"float[880] Å | missing, frame");
    EXPECT_DIAGS_ON(T, "m = mean(zeros((10,)), by=residues)\n", "E0311 1: by= groups particles, and `float[10]` has no dimension over particles\n");
    EXPECT_TYPE_ON(T, "n = within(5.0, of=protein).z\nm = max(n)\n", 2, u8"float Å | missing, frame");
}

UTEST_F(vial_check, math) {
    EXPECT_TYPE("d = sqrt(4.0 'Å^2')\n", 1, u8"float Å, static");
    EXPECT_TYPE("a = asin(0.5)\n", 1, u8"float °, static");
    EXPECT_TYPE("s = sin(30 deg)\n", 1, "float, static");
    EXPECT_DIAGS("e = exp(2.0 nm)\n", "W0501 1: `exp` takes a value without a unit; this is float nm | unitless(...)\n");
    EXPECT_TYPE("d = dot([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])\n", 1, "float, static");
    EXPECT_TYPE("c = cross([1.0, 0.0, 0.0] * 1 nm, [0.0, 1.0, 0.0])\n", 1, "float[3] nm, static");
    EXPECT_TYPE("n = norm([[3.0, 4.0, 0.0]] * 1 nm)\n", 1, "float[1] nm, static");
    EXPECT_TYPE("u = normalize([3.0, 4.0, 0.0] * 1 nm)\n", 1, "float[3], static");
    EXPECT_TYPE("w = where([True, False], [1.0, 2.0], 0.0)\n", 1, "float[2], static");
    EXPECT_TYPE("u = unitless(3.0 nm)\n", 1, "float, static = 3.0");
    EXPECT_TYPE("i = int(3.7)\n", 1, "int, static = 3");
    EXPECT_TYPE("a = atan2(1.0 nm, 2.0 nm)\n", 1, u8"float °, static");
}

UTEST_F(vial_check, building_and_joining) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE_ON(T, "w = zeros(shape(residues(chain='A')))\n", 1, "float[442], static");
    EXPECT_TYPE_ON(T, "n = shape(residues(chain='A'))[0]\n", 1, "int, static = 442");
    EXPECT_TYPE("g = zeros((10, 3))\n", 1, "float[10][3], static");
    EXPECT_TYPE("g = full(4, 2.0 nm)\n", 1, "float[4] nm, static");
    EXPECT_TYPE("g = ones(3, unit='kJ/mol')\n", 1, "float[3] kJ/mol, static");
    EXPECT_TYPE("r = linspace(0.0, 20.0, 81)\n", 1, "float[81], static");
    EXPECT_TYPE("r = arange(0, 10, 3)\n", 1, "float[4], static");
    EXPECT_TYPE("i = eye(3)\n", 1, "float[3][3], static");
    EXPECT_TYPE("r = reshape(zeros((4, 3)), (2, -1))\n", 1, "float[2][6], static");
    EXPECT_DIAGS("r = reshape(zeros((433, 3)), (400, 3))\n", "E0307 1: 1299 values do not fill 1200\n");
    EXPECT_TYPE("f = flatten(zeros((4, 3)))\n", 1, "float[12], static");
    EXPECT_TYPE("d = diff(zeros((4, 3)), axis=0)\n", 1, "float[3][3], static");
    EXPECT_TYPE("c = concatenate([zeros((4, 3)), zeros((2, 3))])\n", 1, "float[6][3], static");
    EXPECT_DIAGS("c = concatenate([zeros((4, 3)), zeros((4, 2))])\n", "E0309 1: `float[4][3]` and `float[4][2]` differ in dimension 1\n");
    EXPECT_TYPE_ON(T, "ab = concatenate([residues(1..4, chain='A'), residues(1..4, chain='B')])\n", 1, "set[8], static");
    EXPECT_TYPE_ON(T, "a, b = common(residues(chain='A', kind='amino_acid'), residues(chain='B', kind='amino_acid'), by='id')\n", 1,
        "set[424] and set[424], static");
    EXPECT_TYPE("w, v = eigh(eye(3) * 1 nm)\n", 1, "float[3] nm and float[3][3], static");
    EXPECT_DIAGS("d = det(zeros((3, 2)))\n", "E0310 1: `det` needs a square matrix; this is `float[3][2]`\n");
    EXPECT_TYPE("m = missing_where([True, False], [1.0, 2.0])\n", 1, "float[2] | missing, static");
    EXPECT_TYPE("m = filled(missing_where([True, False], [1.0, 2.0]), 0.0)\n", 1, "float[2], static");
}

// ### STEP 7: DEFS ###

UTEST_F(vial_check, defs) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE("def twice(v):\n    return 2 * v\nd = twice(3.0 nm)\n", 3, "float nm, static = 6.0");
    EXPECT_TYPE("def f(a, b=2.0):\n    return a + b\nd = f(1.0 nm)\n", 3, "float nm, static = 3.0");
    EXPECT_TYPE("def f(a, b=2.0):\n    return a + b\nd = f(b=1 ns, a=2 ns)\ne = f(5 ns, 1 ns)\n", 3, "int ns, static = 3");
    EXPECT_TYPE("def f(a, b=2.0):\n    return a + b\nd = f(b=1 ns, a=2 ns)\ne = f(5 ns, 1 ns)\n", 4, "int ns, static = 6");
    EXPECT_TYPE("def f(a, op='sum'):\n    if op == 'sum':\n        return a\n    return a > 0\nd = f(1.0)\ne = f(1.0, op='x')\n", 5, "float, static = 1.0");
    EXPECT_TYPE("def f(a, op='sum'):\n    if op == 'sum':\n        return a\n    return a > 0\nd = f(1.0)\ne = f(1.0, op='x')\n", 6, "bool, static = True");
    EXPECT_DIAGS("def f(op):\n    assert op in ('a', 'b'), \"op= takes 'a' or 'b'\"\n    return 1\nx = f('c')\n",
        "E0314 2: op= takes 'a' or 'b'\n");
    EXPECT_DIAGS("def f(a):\n    return a\nx = f()\ny = f(1, 2)\nz = f(b=1)\n",
        "E0204 3: `f` needs a\nE0204 4: `f` takes 1 argument\nE0202 5: `f` takes no keyword `b`\n");
    EXPECT_TYPE_ON(T, "def size(r: set):\n    return count(r)\nn = size(residues(chain='A'))\n", 3, "int[442], static");
    EXPECT_TYPE_ON(T, "def mid(p: points):\n    return mean(p, axis=0)\nm = mid(residues(chain='A', kind='amino_acid'))\n", 3,
        u8"float[433][3] Å, frame");
    EXPECT_TYPE_ON(T, "def f(a: point, b: point):\n    return norm(b - a)\nd = f(residues(kind='water'), particles(name='MG'))\n", 3,
        u8"float[8] Å, frame");
}

UTEST_F(vial_check, prelude_defs) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE_ON(T, "d = distance(particles(name='MG'), particles(resname='GTP').com)\n", 1, u8"float Å, frame");
    EXPECT_TYPE_ON(T, "d = distance(particles(resname='GTP'), particles(resname='VLB'), between='closest')\n", 1, u8"float Å, frame");
    EXPECT_TYPE_ON(T, "d = distance(residues(chain='A', kind='amino_acid'), particles(resname='VLB'), between='farthest')\n", 1,
        u8"float[433] Å, frame");
    EXPECT_TYPE_ON(T, "g = gyration(residues(chain='A', kind='amino_acid'))\n", 1, u8"float[433] Å, frame");
    EXPECT_TYPE_ON(T, "c = centroid(particles(resname='GTP'))\n", 1, u8"float[3] Å, frame");
    EXPECT_DIAGS_ON(T, "g = particles(resname='GTP')\nd = distance(g, particles(resname='GDP'))\n",
        "E0302 2: `displacement` needs a point for a; `g` may hold several particles, so write .com for its centre or between='closest' | g.com\n"
        "E0302 2: `displacement` needs a point for b; `particles(resname='GDP')` may hold several particles, so write .com for its centre"
        " or between='closest' | particles(resname='GDP').com\n");
    EXPECT_DIAGS_ON(T, "g = particles(resname='GTP')\nd = distance(g, g, between='closest')\n",
        "W0608 2: `closest_pair` pairs every point with itself, so the shortest distance is zero; distinct=True leaves those pairs out\n");
    EXPECT_DIAGS_ON(T, "d = distance(particles(name='MG'), particles(name='MG'), between='nearest')\n",
        "E0314 1: between= takes 'closest' or 'farthest'\n");
    EXPECT_DIAGS_ON(T, "g = gyration(protein, weights='charge')\n", "E0314 1: weights= takes 'mass' or None\n");
}

UTEST_F(vial_check, comprehensions) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE_ON(T, "cas = [particles(r, name='CA') for r in residues(chain='A', kind='amino_acid')]\n", 1, "set[433], static");
    EXPECT_TYPE_ON(T, "cas = [particles(r, name='CA') for r in residues(chain='A', kind='amino_acid')]\np = cas.xyz\n", 2,
        u8"float[433][3] Å, frame");
    EXPECT_TYPE_ON(T, "cas = [particles(r, name='CA') for r in residues(chain='A', kind='amino_acid')]\nl = labels(cas)[2]\n", 2,
        "str, static = 'GLU3/A'");
    EXPECT_TYPE_ON(T, "m = [count(r) for r in residues(chain='A')]\n", 1, "int[442], static");
    EXPECT_TYPE_ON(T, "m = [count(r) for r in residues(chain='A') if count(r) > 5]\n", 1, "int[n], static");
    EXPECT_TYPE_ON(T, "s = [i * 2 for i in range(10)]\n", 1, "int[10], static");
    EXPECT_TYPE_ON(T, "s = [i for i in range(2, 11, 3)]\n", 1, "int[3], static");
    EXPECT_TYPE_ON(T, "p = [distance(a.com, b.com) for a, b in zip(residues[0:7], residues[8:15])]\n", 1, u8"float[7] Å, frame");
    EXPECT_DIAGS_ON(T, "p = [a for a, b in zip(residues[0:7], residues[8:14])]\n", "E0301 1: `zip` pairs elements in order, and these have 7 and 6\n");
    EXPECT_TYPE_ON(T, "o = [p for p in particles(resname='GTP')]\n", 1, "set[32], static");
    EXPECT_DIAGS_ON(T, "o = [p for p in 5]\n", "E0313 1: a comprehension iterates over an array, a set, range() or zip(); this is int\n");
}

UTEST_F(vial_check, folds_and_loops) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    const char* run =
        "e = distance(particles(name='MG'), particles(resname='GTP').com)\n"
        "def step(s, v):\n"
        "    return where(v > 30.0, s + 1, 0)\n"
        "run = fold(step, e, init=0)\n"
        "top = max(run, axis='frame')\n";
    EXPECT_TYPE_ON(T, run, 4, "int, trajectory, per frame");
    EXPECT_TYPE_ON(T, run, 5, "int, trajectory");
    EXPECT_DIAGS_ON(T, "e = particles(resname='GTP').com\ndef step(s, v):\n    return v\nr = fold(step, e, init=0.0)\n",
        u8"E0315 4: the step returns float[3] Å; init is float, and a fold's state keeps its type\n");
    EXPECT_TYPE_ON(T, "e = distance(particles(name='MG'), particles(resname='GTP').com)\np = previous(e)\n", 2,
        u8"float Å | missing, trajectory, per frame");
    EXPECT_DIAGS_ON(T, "e = distance(particles(name='MG'), particles(resname='GTP').com)\np = previous(e)\nq = p * 2\n",
        "W0405 2: this may be missing in some frames; the series shows gaps there\n");
    EXPECT_TYPE_ON(T, "e = distance(particles(name='MG'), particles(resname='GTP').com)\nc = cumulative(e, op='max')\n", 2,
        u8"float Å, trajectory, per frame");
    EXPECT_TYPE_ON(T, "f = frame()\n", 1, "int, frame");
    EXPECT_TYPE_ON(T, "t = time()\n", 1, "float64 ps, frame");
    EXPECT_TYPE("def total(n):\n    s = 0\n    for i in range(n):\n        s = s + i\n    return s\nx = total(4)\n", 6, "int, static");
}

UTEST_F(vial_check, w0607) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    const char* src =
        "res_a = residues(chain='A', kind='amino_acid')\n"
        "cut = full(shape(res_a), 4.0)\n"
        "def near(r: set, cutoff=4.0):\n"
        "    return count(r) > cutoff\n"
        "c = near(res_a, cut)\n"
        "def near2(r: set, cutoff: float = 4.0):\n"
        "    return count(r) > cutoff\n"
        "d = near2(res_a, cut)\n";
    EXPECT_TYPE_ON(T, src, 5, "bool[433][433], static");
    EXPECT_TYPE_ON(T, src, 8, "bool[433], static");
    EXPECT_DIAGS_ON(T, src, "W0607 5: `near` maps over int[433], and cutoff runs over the same units but is passed whole, so the result has both;"
        " annotate it, cutoff: float, to pair them\n");
}

// ### STEP 8: STAGES, PASSES AND OUTPUTS ###

static const char* diags_frames(const char* src, const md_system_t* sys, uint64_t frames, md_allocator_i* alloc) {
    const md_vial_compile_desc_t desc = {.sys = sys, .num_frames = frames};
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &desc, alloc);
    const md_vial_diagnostic_t* d = NULL;
    const size_t n = md_vial_diagnostics(prog, &d);
    char* out = NULL;
    for (size_t i = 0; i < n; ++i) {
        uint32_t line = 1;
        for (uint32_t k = 0; k < d[i].span.beg; ++k) line += src[k] == '\n';
        char buf[512];
        int len = snprintf(buf, sizeof(buf), "%s %u\n", d[i].code, line);
        md_array_push_array(out, buf, (size_t)len, alloc);
    }
    md_array_push(out, '\0', alloc);
    md_vial_program_free(prog);
    return out;
}

UTEST_F(vial_check, stages) {
    const md_system_t* P = vial_test_system('P');
    ASSERT_TRUE(P != NULL);
    const char* base = "ca = particles(name='CA')\nends = distance(ca[0], ca[14])\n";
    char src[1024];
#define WITH(tail) (snprintf(src, sizeof(src), "%s%s", base, tail), src)
    EXPECT_TYPE_ON(P, WITH("m = mean(ends, axis='frame')\n"), 3, u8"float Å, trajectory");
    EXPECT_TYPE_ON(P, WITH("m = ends - mean(ends, axis='frame')\n"), 3, u8"float Å, trajectory, per frame");
    EXPECT_TYPE_ON(P, WITH("s = smooth(ends, 5)\n"), 3, u8"float Å, trajectory, per frame");
    EXPECT_TYPE_ON(P, WITH("v = derivative(ends)\n"), 3, u8"float Å/ps, trajectory, per frame");
    EXPECT_TYPE_ON(P, WITH("h = histogram(ends, axis='frame', bins=32, range=(0.0, 50.0))\n"), 3, "float[32], trajectory");
    EXPECT_TYPE_ON(P, WITH("a = autocorrelation(ends, max_lag=10)\n"), 3, "float[11], trajectory");
    EXPECT_TYPE_ON(P, WITH("f = at(ends, frame=0)\n"), 3, u8"float Å, static");
    EXPECT_STREQ("E0405 3\n", diags_frames(WITH("f = at(ends, frame=38)\n"), P, 38, utest_fixture->alloc));
    EXPECT_STREQ("", diags_frames(WITH("f = at(ends, frame=37)\n"), P, 38, utest_fixture->alloc));
    EXPECT_STREQ("E0404 3\n", diags_frames(WITH("f = at(mean(ends, axis='frame'), frame=0)\n"), P, 0, utest_fixture->alloc));
    EXPECT_STREQ("E0404 3\n", diags_frames(WITH("f = at(at(ends, frame=1) + ends, frame=0)\n"), P, 0, utest_fixture->alloc));
    EXPECT_STREQ("E0313 3\n", diags_frames(WITH("s = smooth(ends, 4)\n"), P, 0, utest_fixture->alloc));
    EXPECT_STREQ("E0403 3\n", diags_frames(WITH("if ends < 5.0:\n    pass\n"), P, 0, utest_fixture->alloc));
    EXPECT_STREQ("E0317 3\n", diags_frames(WITH("b = ends is None\n"), P, 0, utest_fixture->alloc));
    EXPECT_STREQ("E0401 3\n", diags_frames(WITH("d = distance(within(3.0, of=ca[0]).xyz, ca[0])\n"), P, 0, utest_fixture->alloc));
    EXPECT_STREQ("", diags_frames(WITH("_d = distance(within(3.0, of=ca[0]).xyz, ca[0])\n"), P, 0, utest_fixture->alloc));
    EXPECT_STREQ("W0404 3\n", diags_frames(WITH("dm = distance(particles.xyz[:, None], particles.xyz)\n"), P, 1000, utest_fixture->alloc));
    EXPECT_STREQ("W0406 5\nW0406 6\n", diags_frames(WITH("a = ends - mean(ends, axis='frame')\nb = a - mean(a, axis='frame')\n"
        "c = b - mean(b, axis='frame')\nd = c - mean(c, axis='frame')\n"), P, 0, utest_fixture->alloc));
#undef WITH
}

// ### STEP 9: DRAWINGS ###

UTEST_F(vial_check, drawings) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_TYPE_ON(T, "l = line(particles(resname='GTP').com, particles(name='MG'))\n", 1, "drawing, frame");
    EXPECT_TYPE_ON(T, "l = line(residues(chain='A').com, particles(name='MG'), label=1.0)\n", 1, "drawing[442], frame");
    EXPECT_TYPE_ON(T, "s = sphere(particles(name='MG'), 3.0) + marker([0.0, 0.0, 0.0])\n", 1, "drawing, frame");
    EXPECT_TYPE_ON(T, "a = arrow([0.0, 0.0, 0.0], direction=[0.0, 0.0, 1.0], length=5.0)\n", 1, "drawing, static");
    EXPECT_TYPE_ON(T, "h = highlight(residues(chain='A'))\n", 1, "drawing[442], static");
    EXPECT_TYPE_ON(T, "t = label(particles(name='MG'), 'MG')\n", 1, "drawing, frame");
    EXPECT_TYPE_ON(T, "d = shown(2.0 nm, sphere([0.0, 0.0, 0.0], 1.0))\n", 1, "float nm, static = 2.0");
    EXPECT_DIAGS_ON(T, "s = sphere(particles(name='MG'), 3.0 ps)\n", "E0313 1: `sphere` takes a length for r; this is float ps\n");
    EXPECT_DIAGS_ON(T, "a = arrow([0.0, 0.0, 0.0])\n", "E0204 1: `arrow` runs from a to b, or from a along direction= for length=\n");
    EXPECT_DIAGS_ON(T, "d = shown(2.0, 3.0)\n", "E0313 1: `shown` attaches a drawing, made by a primitive such as line(); this is float\n");
    EXPECT_DIAGS_ON(T, "b = arrow(particles(resname='GTP'), particles(name='MG'))\n",
        "E0302 1: `arrow` needs a point for a; `particles(resname='GTP')` may hold several particles, so write .com for its centre or"
        " between='closest' | particles(resname='GTP').com\n");
}

UTEST_F(vial_check, bonds_spread_and_set_indices) {
    const md_system_t* T = vial_test_system('T');
    const md_system_t* P = vial_test_system('P');
    ASSERT_TRUE(T != NULL && P != NULL);
    EXPECT_TYPE_ON(T, "b = bonds(particles(resname='GTP'))\nn = len(b)\n", 2, "int, static = 34");
    EXPECT_TYPE_ON(T, "b = bonds(particles(resname='GTP'))\nl = distance(b[:, 0], b[:, 1])\n", 2, u8"float[34] Å, frame");
    EXPECT_TYPE_ON(P, "p = particles.xyz\nz = p[particles(name='CA')]\n", 2, u8"float[15][3] Å, frame");
    EXPECT_DIAGS_ON(P, "z = zeros(4)[particles(name='CA')]\n", "E0311 1: a set picks values by particle, and `float[4]` has no dimension over particles\n");
    EXPECT_TYPE_ON(T, "s = spread(count(residues(chain='A')), residues(chain='A'))\n", 1, "int[3441], static");   // 3393 protein, 32 GTP, 1 MG, 10 SO4, 5 water
    EXPECT_DIAGS_ON(T, "r = residues(chain='A')[0:2]\ns = spread(count(concatenate([r, r])), concatenate([r, r]))\n",
        "E0306 2: the groups share particles, so each of those would get several values; combine= says how | combine='mean'\n");
}

// attr() on a run whose energy terms are attributes, as an energy file gives them
#include <md_attributes.h>

static const char* type_run(const char* src, uint32_t line, const md_system_t* sys, md_allocator_i* alloc, const char** diags) {
    const md_vial_compile_desc_t desc = {.sys = sys, .run = STR_LIT("run/md")};
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &desc, alloc);
    static char out[256];
    md_vial__result_t r;
    if (!md_vial__line_result(prog, line, &r)) snprintf(out, sizeof(out), "(none)");
    else snprintf(out, sizeof(out), "%s, %s", r.types[0], r.stage);
    const md_vial_diagnostic_t* d = NULL;
    const size_t n = md_vial_diagnostics(prog, &d);
    static char dbuf[512];
    dbuf[0] = '\0';
    for (size_t i = 0; i < n; ++i) snprintf(dbuf + strlen(dbuf), sizeof(dbuf) - strlen(dbuf), "%s ", d[i].code);
    *diags = dbuf;
    md_vial_program_free(prog);
    return out;
}

UTEST_F(vial_check, attributes) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    md_system_t sys = *T;
    sys.attributes = (md_attributes_t){.alloc = utest_fixture->alloc};
    md_unit_t kj_mol;
    vial_unit_string(STR_LIT("kJ/mol"), &kj_mol);
    const md_attribute_desc_t time = {.path = STR_LIT("run/md/time"), .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {38}}};
    ASSERT_TRUE(md_attributes_create(&sys.attributes, &time) != MD_ATTRIBUTE_INVALID);
    const char* terms[] = {"run/md/edr/Potential", "run/md/edr/coul_sr_a", "run/md/edr/coul_sr_b"};
    for (int i = 0; i < 3; ++i) {
        const md_attribute_desc_t d = {.path = str_from_cstr(terms[i]), .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = kj_mol,
            .format = {.type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 1, .shape = {38}}};
        ASSERT_TRUE(md_attributes_create(&sys.attributes, &d) != MD_ATTRIBUTE_INVALID);
    }
    const char* diags;
    EXPECT_STREQ("float kJ/mol, frame", type_run("pot = attr('edr/Potential')\n", 1, &sys, utest_fixture->alloc, &diags));
    EXPECT_STREQ("", diags);
    EXPECT_STREQ("float[2] kJ/mol, frame", type_run("c = attr('edr/coul_sr_*')\n", 1, &sys, utest_fixture->alloc, &diags));
    EXPECT_STREQ("float kJ/mol, frame", type_run("c = sum(attr('edr/coul_sr_*'))\n", 1, &sys, utest_fixture->alloc, &diags));
    type_run("p = attr('edr/Kinetic')\n", 1, &sys, utest_fixture->alloc, &diags);
    EXPECT_STREQ("E0206 ", diags);
    // The run's frame axis gives the number of frames
    type_run("e = attr('edr/Potential')\nf = at(e, frame=38)\n", 2, &sys, utest_fixture->alloc, &diags);
    EXPECT_STREQ("E0405 ", diags);
}

// ### STEP 10: PAIRING, OUTPUTS AND HOVER ###

UTEST_F(vial_check, w0605) {
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    EXPECT_DIAGS_ON(T, "a = residues(1..4, chain='A').com\nb = residues(1..4, chain='B').com\nd = a - b\n",
        "W0605 3: dimensions of 4 that run over different units are paired by position: 0 is MET1/A against MET1/B;"
        " common(a, b, by='id') pairs them by id\n");
    EXPECT_DIAGS_ON(T, "c = residues(chain='A').com\nd = c[1:5] - c[0:4]\n", "");
    EXPECT_DIAGS_ON(T, "r = residues(1..4, chain='A')\ns = residues(1..4, chain='A')\nd = r.com - s.com\n", "");
    EXPECT_DIAGS_ON(T, "a, b = common(residues(1..4, chain='A'), residues(1..4, chain='B'), by='id')\nd = distance(a.com, b.com)\n", "");
    EXPECT_DIAGS_ON(T, "g = particles(resname='GTP')\nd = distance(g[0:3].xyz, g[3:6].xyz)\n", "");
    EXPECT_DIAGS_ON(T, "d = distance(particles(name='CA', resid=1..4, chain='A').xyz, particles(name='CB', resid=1..4, chain='A').xyz)\n",
        "W0605 1: dimensions of 4 that run over different units are paired by position: 0 is CA/MET1/A against CB/MET1/A;"
        " common(a, b, by='id') pairs them by id\n");
}

UTEST_F(vial_check, outputs_and_hover) {
    const md_system_t* P = vial_test_system('P');
    ASSERT_TRUE(P != NULL);
    const char* src =
        "ca    = particles(name='CA')\n"
        "_p    = particles.xyz\n"
        "rmsf  = sqrt(sum(var(_p, axis='frame'), axis=1))\n"
        "ends  = distance(ca[0], ca[14])\n"
        "h     = histogram(ends, axis='frame', bins=32, range=(0.0, 50.0))\n"
        "near  = residues(within(4.0, of=ca[0]))\n"
        "name0 = labels(residues)[0]\n"
        "dev   = ends - mean(ends, axis='frame')\n"
        "big   = zeros((2, 2, 2))\n";
    const md_vial_compile_desc_t desc = {.sys = P};
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &desc, utest_fixture->alloc);
    const md_vial_output_t* o = NULL;
    const size_t n = md_vial_outputs(prog, &o);
    ASSERT_EQ(6, (int)n);
    EXPECT_TRUE(str_eq(o[0].name, STR_LIT("ca")));   EXPECT_TRUE(o[0].selection);  EXPECT_EQ(MD_VIAL_STAGE_STATIC, o[0].stage);
    EXPECT_TRUE(str_eq(o[1].name, STR_LIT("rmsf"))); EXPECT_STREQ(u8"float[153] Å", o[1].type); EXPECT_EQ(MD_VIAL_STAGE_TRAJECTORY, o[1].stage);
    EXPECT_TRUE(str_eq(o[2].name, STR_LIT("ends"))); EXPECT_EQ(MD_VIAL_STAGE_FRAME, o[2].stage); EXPECT_EQ(1, (int)o[2].pass);
    EXPECT_TRUE(str_eq(o[3].name, STR_LIT("h")));    EXPECT_STREQ("float[32]", o[3].type);
    EXPECT_TRUE(str_eq(o[4].name, STR_LIT("near"))); EXPECT_TRUE(o[4].selection); EXPECT_EQ(MD_VIAL_STAGE_FRAME, o[4].stage);
    EXPECT_TRUE(str_eq(o[5].name, STR_LIT("dev")));  EXPECT_EQ(MD_VIAL_STAGE_TRAJECTORY_PER_FRAME, o[5].stage); EXPECT_EQ(2, (int)o[5].pass);

    // Hover: the innermost expression at an offset
    md_vial_hover_t hv;
    const uint32_t at_ca0 = (uint32_t)(strstr(src, "ca[0], ca[14]") - src) + 1;    // inside `ca` of ca[0]
    ASSERT_TRUE(md_vial_hover(prog, at_ca0, &hv));
    EXPECT_STREQ("set", hv.type);
    const uint32_t at_idx = (uint32_t)(strstr(src, "ca[0], ca[14]") - src) + 3;    // the 0
    ASSERT_TRUE(md_vial_hover(prog, at_idx, &hv));
    EXPECT_STREQ("int", hv.type);
    EXPECT_STREQ("0", hv.value);
    const uint32_t at_name0 = (uint32_t)(strstr(src, "labels") - src);
    ASSERT_TRUE(md_vial_hover(prog, at_name0, &hv));
    EXPECT_STREQ("str[15]", hv.type);
    const uint32_t at_dist = (uint32_t)(strstr(src, "distance") - src);
    ASSERT_TRUE(md_vial_hover(prog, at_dist, &hv));
    EXPECT_STREQ(u8"float Å", hv.type);
    EXPECT_EQ(MD_VIAL_STAGE_FRAME, hv.stage);
    EXPECT_FALSE(md_vial_hover(prog, 0xFFFFFF, &hv));
    md_vial_program_free(prog);
}
