#include "ir.h"
#include "ast.h"
#include "builtins.h"

#include <core/md_allocator.h>
#include <core/md_array.h>

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

void vial_ir_init(vial_ir_t* ir, md_allocator_i* alloc) {
    memset(ir, 0, sizeof(*ir));
    ir->alloc = alloc;
    md_array_push(ir->values, (vial_ir_value_t){0}, alloc);
    md_array_push(ir->regions, (vial_ir_region_t){0}, alloc);
    const uint32_t f = vial_ir_func(ir, STR_LIT(""), (md_vial_span_t){0}, NULL, 0);
    ir->insert = ir->funcs[f].body;
}

uint32_t vial_ir_region(vial_ir_t* ir, uint32_t owner, uint32_t func, const vial_type_t* params, uint32_t num_params) {
    const uint32_t r = (uint32_t)md_array_size(ir->regions);
    vial_ir_region_t reg = {.param_beg = (vial_ssa_t)md_array_size(ir->values), .num_params = num_params, .owner = owner, .func = func};
    md_array_push(ir->regions, reg, ir->alloc);
    for (uint32_t i = 0; i < num_params; ++i) {
        md_array_push(ir->values, ((vial_ir_value_t){.type = params[i], .op = VIAL_IR_NO_OP, .region = r, .index = i}), ir->alloc);
    }
    return r;
}

void vial_ir_attach(vial_ir_t* ir, uint32_t op, uint32_t k, uint32_t region) {
    ir->ops[op].region[k] = region;
    ir->regions[region].owner = op;
}

uint32_t vial_ir_func(vial_ir_t* ir, str_t name, md_vial_span_t span, const vial_type_t* params, uint32_t num_params) {
    const uint32_t f = (uint32_t)md_array_size(ir->funcs);
    md_array_push(ir->funcs, ((vial_ir_func_t){.name = name, .span = span}), ir->alloc);
    ir->funcs[f].body = vial_ir_region(ir, VIAL_IR_NO_OP, f, params, num_params);
    return f;
}

uint32_t vial_ir_emit(vial_ir_t* ir, const vial_ir_emit_t* e) {
    const uint32_t op = (uint32_t)md_array_size(ir->ops);
    vial_ir_op_t o = {
        .kind = e->kind, .sub = e->sub, .aux = e->aux,
        .arg_beg = (uint32_t)md_array_size(ir->args), .num_args = e->num_args,
        .result = (vial_ssa_t)md_array_size(ir->values), .num_results = e->num_results,
        .parent = ir->insert, .span = e->span, .attr = (uint32_t)md_array_size(ir->ints),
    };
    for (uint32_t i = 0; i < e->num_attr; ++i) md_array_push(ir->ints, e->attr[i], ir->alloc);
    for (uint32_t i = 0; i < e->num_args; ++i) {
        md_array_push(ir->args, e->args[i], ir->alloc);
        md_array_push(ir->names, e->names ? e->names[i] : (str_t){0}, ir->alloc);
    }
    for (uint32_t i = 0; i < e->num_results; ++i) {
        md_array_push(ir->values, ((vial_ir_value_t){.type = e->types[i], .op = op, .region = ir->insert, .index = i}), ir->alloc);
    }
    if (!e->num_results) o.result = 0;
    md_array_push(ir->ops, o, ir->alloc);
    md_array_push(ir->regions[ir->insert].ops, op, ir->alloc);
    return op;
}

vial_ssa_t vial_ir_const(vial_ir_t* ir, const vial_val_t* v, const vial_type_t* t, str_t text, md_vial_span_t span) {
    const uint32_t c = (uint32_t)md_array_size(ir->consts);
    md_array_push(ir->consts, *v, ir->alloc);
    md_array_push(ir->const_text, text, ir->alloc);
    const uint32_t op = vial_ir_emit(ir, &(vial_ir_emit_t){.kind = VIAL_IR_CONST, .aux = c, .types = t, .num_results = 1, .span = span});
    return ir->ops[op].result;
}

vial_ssa_t vial_ir_op1(vial_ir_t* ir, uint16_t kind, uint16_t sub, uint32_t aux, const vial_ssa_t* args, uint32_t n,
                       const vial_type_t* t, md_vial_span_t span) {
    const uint32_t op = vial_ir_emit(ir, &(vial_ir_emit_t){.kind = kind, .sub = sub, .aux = aux, .args = args, .num_args = n,
                                                           .types = t, .num_results = 1, .span = span});
    return ir->ops[op].result;
}

void vial_ir_yield(vial_ir_t* ir, const vial_ssa_t* args, uint32_t n, md_vial_span_t span) {
    vial_ir_emit(ir, &(vial_ir_emit_t){.kind = VIAL_IR_YIELD, .args = args, .num_args = n, .span = span});
}

const char* vial_ir_kind_str(vial_ir_kind_t kind) {
    static const char* s[] = {"none", "const", "unary", "not", "binary", "fma", "compare", "logic", "setop", "convert",
                              "index", "list", "tuple", "elem", "builtin", "call", "apply", "if", "for", "map", "kind", "fold", "pairs", "yield"};
    return kind < VIAL_IR_KIND_COUNT ? s[kind] : "?";
}

// ### VERIFIER ###

typedef struct verifier_t {
    const vial_ir_t* ir;
    uint8_t*         defined;   // by ssa
    char*            err;
    size_t           cap;
    bool             ok;
} verifier_t;

static void fail(verifier_t* v, const char* fmt, uint32_t a, uint32_t b) {
    if (!v->ok) return;
    v->ok = false;
    if (v->cap) snprintf(v->err, v->cap, fmt, a, b);
}

static uint32_t yields_of(const vial_ir_t* ir, uint32_t region, bool* has) {
    const vial_ir_region_t* r = &ir->regions[region];
    const size_t n = md_array_size(r->ops);
    *has = n && ir->ops[r->ops[n - 1]].kind == VIAL_IR_YIELD;
    return *has ? ir->ops[r->ops[n - 1]].num_args : 0;
}

static void verify_region(verifier_t* v, uint32_t region, int expect_yields);

static void verify_op(verifier_t* v, uint32_t opi, bool last) {
    const vial_ir_t* ir = v->ir;
    const vial_ir_op_t* op = &ir->ops[opi];
    for (uint32_t i = 0; i < op->num_args; ++i) {
        const vial_ssa_t a = ir->args[op->arg_beg + i];
        if (!a) {
            if (op->kind != VIAL_IR_BUILTIN && !(op->kind == VIAL_IR_PAIRS && i == 2)) fail(v, "op %u: argument %u is none", opi, i);
            continue;
        }
        if (a >= md_array_size(ir->values) || !v->defined[a]) fail(v, "op %u uses %%%u before it is defined", opi, a);
    }
    if (op->kind == VIAL_IR_YIELD && !last) fail(v, "op %u: a yield before the end of region %u", opi, op->parent);
    switch (op->kind) {
    case VIAL_IR_UNARY: case VIAL_IR_NOT: case VIAL_IR_CONVERT: case VIAL_IR_ELEM: case VIAL_IR_KIND:
        if (op->num_args != 1) fail(v, "op %u takes 1 argument, has %u", opi, op->num_args);
        break;
    case VIAL_IR_BINARY: case VIAL_IR_COMPARE: case VIAL_IR_LOGIC: case VIAL_IR_SETOP:
        if (op->num_args != 2) fail(v, "op %u takes 2 arguments, has %u", opi, op->num_args);
        break;
    case VIAL_IR_FMA:
        if (op->num_args != 3) fail(v, "op %u takes 3 arguments, has %u", opi, op->num_args);
        break;
    case VIAL_IR_PAIRS:
        if (op->num_args != 3 || !op->num_results) fail(v, "op %u takes 3 arguments and gives a value, has %u", opi, op->num_args);
        break;
    case VIAL_IR_CONST:
        if (op->aux >= md_array_size(ir->consts)) fail(v, "op %u: constant %u does not exist", opi, op->aux);
        break;
    case VIAL_IR_CALL:
        if (op->aux == 0 || op->aux >= md_array_size(ir->funcs)) fail(v, "op %u calls function %u, which does not exist", opi, op->aux);
        else if (ir->regions[ir->funcs[op->aux].body].num_params != op->num_args) fail(v, "op %u: %u arguments for a function of other arity", opi, op->num_args);
        break;
    case VIAL_IR_IF:
        if (op->num_args != 1 || !op->region[0] || !op->region[1]) fail(v, "op %u: an if needs a condition and two regions", opi, 0);
        break;
    case VIAL_IR_FOR: case VIAL_IR_MAP: case VIAL_IR_FOLD: case VIAL_IR_APPLY:
        if (!op->region[0]) fail(v, "op %u needs a region", opi, 0);
        break;
    default: break;
    }
    for (int k = 0; k < 2; ++k) {
        if (!op->region[k]) continue;
        if (ir->regions[op->region[k]].owner != opi) fail(v, "op %u: region %u has another owner", opi, op->region[k]);
        const int expect = op->kind == VIAL_IR_FOR ? (int)op->num_args - 1 : op->kind == VIAL_IR_FOLD ? 1 : (int)op->num_results;
        verify_region(v, op->region[k], op->kind == VIAL_IR_MAP ? -1 : expect);
    }
    for (uint32_t i = 0; i < op->num_results; ++i) v->defined[op->result + i] = 1;
}

static void verify_region(verifier_t* v, uint32_t region, int expect_yields) {
    const vial_ir_t* ir = v->ir;
    const vial_ir_region_t* r = &ir->regions[region];
    for (uint32_t i = 0; i < r->num_params; ++i) v->defined[r->param_beg + i] = 1;
    const size_t n = md_array_size(r->ops);
    for (size_t i = 0; i < n; ++i) {
        if (ir->ops[r->ops[i]].parent != region) fail(v, "op %u is listed in region %u but belongs to another", r->ops[i], region);
        verify_op(v, r->ops[i], i + 1 == n);
    }
    bool has;
    const uint32_t ny = yields_of(ir, region, &has);
    if (!has) fail(v, "region %u does not end in a yield", region, 0);
    else if (expect_yields >= 0 && ny != (uint32_t)expect_yields) fail(v, "region %u yields %u values", region, ny);
    // Its values go out of scope
    for (uint32_t i = 0; i < r->num_params; ++i) v->defined[r->param_beg + i] = 0;
    for (size_t i = 0; i < n; ++i) {
        const vial_ir_op_t* op = &ir->ops[r->ops[i]];
        for (uint32_t k = 0; k < op->num_results; ++k) v->defined[op->result + k] = 0;
    }
}

bool vial_ir_verify(const vial_ir_t* ir, char* err, size_t cap) {
    verifier_t v = {.ir = ir, .err = err, .cap = cap, .ok = true};
    if (cap) err[0] = 0;
    const size_t nv = md_array_size(ir->values);
    v.defined = md_alloc(ir->alloc, nv ? nv : 1);
    memset(v.defined, 0, nv ? nv : 1);
    // A def's body may read the script's top-level names, which are bound once each before any call that reads them
    // (the check sees to that): the top level's values are defined in every function
    const vial_ir_region_t* top = &ir->regions[ir->funcs[0].body];
    for (size_t f = 1; f < md_array_size(ir->funcs); ++f) {
        for (size_t i = 0; i < md_array_size(top->ops); ++i) {
            const vial_ir_op_t* op = &ir->ops[top->ops[i]];
            for (uint32_t k = 0; k < op->num_results; ++k) v.defined[op->result + k] = 2;
        }
        const uint32_t body = ir->funcs[f].body;
        bool has;
        const uint32_t ny = yields_of(ir, body, &has);
        verify_region(&v, body, (int)(has ? ny : 1));
        for (size_t i = 0; i < md_array_size(top->ops); ++i) {
            const vial_ir_op_t* op = &ir->ops[top->ops[i]];
            for (uint32_t k = 0; k < op->num_results; ++k) v.defined[op->result + k] = 0;
        }
    }
    verify_region(&v, ir->funcs[0].body, 0);
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) {
        const vial_ssa_t o = ir->outputs[i].value;
        if (!o || o >= nv || ir->values[o].region != ir->funcs[0].body) fail(&v, "output %u is not a value of the top level", (uint32_t)i, 0);
    }
    md_free(ir->alloc, v.defined, nv ? nv : 1);
    return v.ok;
}

// ### DUMP ###

typedef struct dumper_t {
    const vial_ir_t* ir;
    str_t            src;
    char*            buf;   // md_array
    md_allocator_i*  alloc;
} dumper_t;

static void sb_fmt(dumper_t* d, const char* fmt, ...) {
    char tmp[512];
    va_list args;
    va_start(args, fmt);
    int len = vsnprintf(tmp, sizeof(tmp), fmt, args);
    va_end(args);
    if (len < 0) return;
    if ((size_t)len >= sizeof(tmp)) len = (int)sizeof(tmp) - 1;
    md_array_push_array(d->buf, tmp, (size_t)len, d->alloc);
}
static void sb_cstr(dumper_t* d, const char* s) { md_array_push_array(d->buf, s, strlen(s), d->alloc); }
static void sb_char(dumper_t* d, char c) { md_array_push(d->buf, c, d->alloc); }

static void line_col(str_t src, uint32_t off, uint32_t* line, uint32_t* col) {
    uint32_t l = 1, c = 1;
    for (uint32_t i = 0; i < off && i < src.len; ++i) {
        if (src.ptr[i] == '\n') { ++l; c = 1; } else if (((uint8_t)src.ptr[i] & 0xC0) != 0x80) ++c;
    }
    *line = l; *col = c;
}

static void dump_ssa(dumper_t* d, vial_ssa_t v) {
    if (v) sb_fmt(d, "%%%u", v);
    else sb_cstr(d, "_");
}

static const char* op_text(uint16_t kind, uint16_t sub) {
    switch (kind) {
    case VIAL_IR_UNARY:   return vial_op_str((vial_op_t)sub);
    case VIAL_IR_BINARY:  return vial_op_str((vial_op_t)sub);
    case VIAL_IR_COMPARE: return vial_cmp_str((vial_cmp_t)sub);
    case VIAL_IR_LOGIC:   return sub == VIAL_KW_AND ? "and" : "or";
    case VIAL_IR_SETOP:   return sub == VIAL_KW_AND ? "and" : sub == VIAL_KW_OR ? "or" : sub == VIAL_KW_NOT ? "not" : "and not";
    default:              return vial_ir_kind_str((vial_ir_kind_t)kind);
    }
}

static void dump_region(dumper_t* d, uint32_t region, int indent);

static void dump_op(dumper_t* d, uint32_t opi, int indent) {
    const vial_ir_t* ir = d->ir;
    const vial_ir_op_t* op = &ir->ops[opi];
    const size_t line_start = md_array_size(d->buf);
    sb_fmt(d, "%*s", indent * 2, "");
    if (op->kind == VIAL_IR_YIELD) {
        sb_cstr(d, "yield");
        for (uint32_t i = 0; i < op->num_args; ++i) { sb_cstr(d, i ? ", " : " "); dump_ssa(d, ir->args[op->arg_beg + i]); }
        sb_char(d, '\n');
        return;
    }
    for (uint32_t i = 0; i < op->num_results; ++i) {
        char tbuf[160];
        vial_type_print(tbuf, sizeof(tbuf), &ir->values[op->result + i].type);
        sb_fmt(d, "%s%%%u: %s", i ? ", " : "", op->result + i, tbuf);
    }
    // Pad to a column, so the right-hand sides line up
    {
        const size_t len = md_array_size(d->buf) - line_start;
        const size_t col = 32 + indent * 2;
        if (len < col) sb_fmt(d, "%*s", (int)(col - len), "");
    }
    sb_cstr(d, " = ");
    const vial_ssa_t* args = op->num_args ? ir->args + op->arg_beg : NULL;
    switch (op->kind) {
    case VIAL_IR_CONST: {
        char vbuf[96];
        const str_t text = ir->const_text[op->aux];
        if (text.len) sb_fmt(d, "const %.*s", (int)(text.len < 60 ? text.len : 60), text.ptr);
        else { vial_val_print(vbuf, sizeof(vbuf), &ir->consts[op->aux], 8); sb_fmt(d, "const %s", vbuf); }
        break;
    }
    case VIAL_IR_UNARY: case VIAL_IR_NOT:
        sb_fmt(d, "%s ", op_text(op->kind, op->sub));
        dump_ssa(d, args[0]);
        break;
    case VIAL_IR_BINARY: case VIAL_IR_COMPARE: case VIAL_IR_LOGIC:
        dump_ssa(d, args[0]);
        sb_fmt(d, " %s ", op_text(op->kind, op->sub));
        dump_ssa(d, args[1]);
        break;
    case VIAL_IR_SETOP:
        if (op->sub == VIAL_KW_NOT) { sb_cstr(d, "not "); dump_ssa(d, args[1]); sb_cstr(d, " in "); dump_ssa(d, args[0]); }
        else { dump_ssa(d, args[0]); sb_fmt(d, " %s ", op_text(op->kind, op->sub)); dump_ssa(d, args[1]); }
        break;
    case VIAL_IR_INDEX: {
        dump_ssa(d, args[0]);
        sb_char(d, '[');
        for (uint32_t i = 0; i < op->sub; ++i) {
            const vial_index_item_t* it = &ir->index_items[op->aux + i];
            if (i) sb_cstr(d, ", ");
            switch (it->kind) {
            case VIAL_INDEX_SCALAR: case VIAL_INDEX_ARRAY: dump_ssa(d, it->a); break;
            case VIAL_INDEX_SET: sb_cstr(d, "set "); dump_ssa(d, it->a); sb_cstr(d, " in "); dump_ssa(d, it->b); break;
            case VIAL_INDEX_NEWAXIS: sb_cstr(d, "None"); break;
            case VIAL_INDEX_ELLIPSIS: sb_cstr(d, "..."); break;
            case VIAL_INDEX_SLICE:
                if (it->a) dump_ssa(d, it->a);
                sb_char(d, ':');
                if (it->b) dump_ssa(d, it->b);
                if (it->c) { sb_char(d, ':'); dump_ssa(d, it->c); }
                break;
            }
        }
        sb_char(d, ']');
        break;
    }
    case VIAL_IR_ELEM:
        dump_ssa(d, args[0]);
        sb_fmt(d, ".%u", op->aux);
        break;
    case VIAL_IR_KIND: {
        static const char* kinds[] = {"whole", "set", "group", "point", "points", "float", "int", "bool"};
        sb_fmt(d, "as %s ", op->sub == VIAL_KIND_EACH ? "each" : op->sub < 8 ? kinds[op->sub] : "?");
        dump_ssa(d, args[0]);
        break;
    }
    default: {
        const char* name = vial_ir_kind_str((vial_ir_kind_t)op->kind);
        if (op->kind == VIAL_IR_BUILTIN) {
            const vial_builtin_t* b = vial_builtin_at(op->aux);
            name = b ? b->name : "?";
        }
        if (op->kind == VIAL_IR_CALL) sb_fmt(d, "call f%u", op->aux);
        else sb_cstr(d, name);
        sb_char(d, '(');
        bool first = true;
        for (uint32_t i = 0; i < op->num_args; ++i) {
            if (!args[i] && op->kind == VIAL_IR_BUILTIN) continue;
            if (!first) sb_cstr(d, ", ");
            first = false;
            const str_t nm = ir->names[op->arg_beg + i];
            if (nm.len && op->kind == VIAL_IR_BUILTIN && i) sb_fmt(d, "%.*s=", (int)nm.len, nm.ptr);
            dump_ssa(d, args[i]);
        }
        sb_char(d, ')');
        if (op->kind == VIAL_IR_MAP) {
            sb_cstr(d, " mapping");
            for (uint32_t i = 0; i < op->num_args; ++i) sb_fmt(d, " %u", ir->ints[op->attr + i]);
        }
        if (op->kind == VIAL_IR_PAIRS) {
            static const char* reds[] = {"?", "min", "max", "argmin", "argmax", "count", "any"};
            const uint32_t axis = ir->ints[op->attr], test = ir->ints[op->attr + 1];
            sb_fmt(d, " %s %s%s", op->sub < 7 ? reds[op->sub] : "?", axis == VIAL_PAIRS_ALL ? "all" : axis ? "rows" : "columns",
                test == VIAL_CMP_LT ? " d < r" : test == VIAL_CMP_LE ? " d <= r" : "");
        }
    }
    }
    if (op->span.end) {
        uint32_t l0, c0, l1, c1;
        line_col(d->src, op->span.beg, &l0, &c0);
        line_col(d->src, op->span.end, &l1, &c1);
        sb_fmt(d, "  @ %u:%u-%u:%u", l0, c0, l1, c1);
    }
    if (op->num_results) sb_fmt(d, "  %s", vial_stage_str(&ir->values[op->result].type));
    sb_char(d, '\n');
    for (int k = 0; k < 2; ++k) {
        if (!op->region[k]) continue;
        sb_fmt(d, "%*s%s", indent * 2, "", op->kind == VIAL_IR_IF ? (k ? "else" : "then") : "do");
        const vial_ir_region_t* r = &ir->regions[op->region[k]];
        for (uint32_t i = 0; i < r->num_params; ++i) sb_fmt(d, "%s%%%u", i ? ", " : " ", r->param_beg + i);
        sb_cstr(d, ":\n");
        dump_region(d, op->region[k], indent + 1);
    }
}

static void dump_region(dumper_t* d, uint32_t region, int indent) {
    const vial_ir_region_t* r = &d->ir->regions[region];
    for (size_t i = 0; i < md_array_size(r->ops); ++i) dump_op(d, r->ops[i], indent);
}

str_t vial_ir_dump(const vial_ir_t* ir, str_t src, md_allocator_i* alloc) {
    dumper_t d = {.ir = ir, .src = src, .alloc = alloc};
    for (size_t f = 1; f < md_array_size(ir->funcs); ++f) {
        const vial_ir_func_t* fn = &ir->funcs[f];
        const vial_ir_region_t* r = &ir->regions[fn->body];
        sb_fmt(&d, "def f%u %.*s(", (uint32_t)f, (int)fn->name.len, fn->name.ptr);
        for (uint32_t i = 0; i < r->num_params; ++i) {
            char tbuf[160];
            vial_type_print(tbuf, sizeof(tbuf), &ir->values[r->param_beg + i].type);
            sb_fmt(&d, "%s%%%u: %s", i ? ", " : "", r->param_beg + i, tbuf);
        }
        sb_cstr(&d, "):\n");
        dump_region(&d, fn->body, 1);
    }
    dump_region(&d, ir->funcs[0].body, 0);
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) {
        const vial_ir_output_t* o = &ir->outputs[i];
        sb_fmt(&d, "out %.*s = %%%u\n", (int)o->name.len, o->name.ptr, o->value);
    }
    return (str_t){d.buf, md_array_size(d.buf)};
}
