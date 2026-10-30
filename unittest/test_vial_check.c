// The checker, step by step

#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>

#include <stdio.h>
#include <string.h>

#include <md_vial.h>

// The diagnostics of a script as "CODE line: message | fix", one per line
static const char* diags_of(const char* src, md_allocator_i* alloc) {
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), NULL, alloc);
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

UTEST_F(vial_check, recursion) {
    EXPECT_DIAGS("def walk(n):\n    return walk(n - 1)\n", "E0802 2: `walk` calls itself; VIAL has no recursion\n");
    EXPECT_DIAGS("def a(n):\n    return b(n)\ndef b(n):\n    return a(n)\n",
        "E0802 4: `b` calls `a`, which calls `b` again; VIAL has no recursion\n");
    EXPECT_DIAGS("def a(n):\n    return n\ndef b(n):\n    return a(n) + a(1)\nv = b(2)\n", "");
}

#include <vial/vial_internal.h>

// The type the checker gives the statement on a line, as "type, stage = value"
static const char* type_at(const char* src, uint32_t line, md_allocator_i* alloc) {
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), NULL, alloc);
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
