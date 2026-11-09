#!/usr/bin/env python3
"""The md_common.h and md_str.h namespacing rename, re-runnable on any checkout of mdlib or VIAMD.

    tools/rename_md_prefix.py [--dry-run] [--cpp-headers] <dir or file> ...

mdlib:  tools/rename_md_prefix.py src unittest benchmark examples
VIAMD:  ext/mdlib/tools/rename_md_prefix.py --cpp-headers src test

Renames the macros, types and functions of md_common.h and md_str.h to their MD_ / md_ names in C, C++ and Objective-C
sources (.c .h .inl .cpp .hpp .m .mm). It tokenizes enough of the language to never touch a string,
character or raw string literal, and to tell code from comments:

  code      md_common.h's functions, object-like macros and str_t: every token. Function-like macros and
            md_str.h's functions: where the name is followed by '(' (so 'static constexpr int MAX',
            'PanelMeshes::MAX', an allocator parameter named str_alloc or a local str_len stay), and every
            occurrence on a preprocessor directive line (#ifdef MIN, #define MY_MAX MAX).
  comments  only names that cannot be prose: the function names, and the macro names that are
            unmistakably identifiers (ARRAY_SIZE, MEMCPY, str_t, ...). MIN, MAX, PI, DEBUG, RELEASE, CLAMP,
            CONCAT, STRINGIFY, and the str_ functions that double as variable names (str_len, str_beg,
            str_end, str_ptr, str_alloc) only where followed by '('. 'defer' never (it is English).
  'defer'   C++ only (.cpp .hpp .mm, and .h under a C++ tree passed with --cpp-headers); in C it is an
            ordinary identifier.

Directories named 'shaders' are skipped: GLSL/Slang, whose shared headers define their own names.
md_common.h itself is renamed like any file; its definitions are then edited by hand.
"""
import os, re, sys

FUNCS = ['popcnt32', 'popcnt64', 'clz32', 'clz64', 'ctz32', 'ctz64', 'bsf32', 'bsf64', 'bsr32', 'bsr64',
         'find_first_zero_nibble16', 'find_first_zero_nibble32', 'find_first_zero_nibble64',
         'find_first_zero_byte32', 'find_first_zero_byte64', 'next_power_of_two32', 'next_power_of_two64']
OBJ_MACROS = ['STATIC_ASSERT', 'THREAD_LOCAL', 'FALLTHROUGH', 'RESTRICT', 'FORCE_INLINE', 'NOINLINE', 'DEBUG',
              'RELEASE', 'PI', 'TWO_PI', 'MEMCPY', 'MEMSET', 'MEMMOVE', 'MEMCMP', 'BSWAP16', 'BSWAP32', 'BSWAP64']
FN_MACROS = ['ALIGNAS', 'TYPEOF', 'ARRAY_SIZE', 'TYPE_COMPATIBLE', 'IS_POW2', 'ALIGN_TO', 'IS_ALIGNED',
             'NEXT_ALIGNED_ADDRESS', 'MIN', 'MAX', 'CLAMP', 'DIV_UP', 'ROUND_UP', 'ROUND_DOWN', 'DEG_TO_RAD',
             'RAD_TO_DEG', 'KILOBYTES', 'MEGABYTES', 'GIGABYTES', 'CONCAT', 'CONCAT_INTERNAL', 'STRINGIFY',
             'STRINGIFY_VAL', 'ASSERT', 'DEBUG_BREAK', 'ENUM_FLAGS']
SPECIAL = {'ARRAY_SIZE': 'MD_COUNTOF'}      # MD_ARRAY_SIZE is md_array.h's, for its dynamic arrays
# Aliases of a function or intrinsic are named like functions: a call behaves exactly as a function call (each
# argument evaluated once, nothing compiled away), unlike MD_MIN or MD_ASSERT, whose capitals say they do not.
# The MD_ spellings are those of the first version of this rename, mapped too so that it runs on either base.
ALIASES = {'MEMCPY': 'md_memcpy', 'MEMSET': 'md_memset', 'MEMMOVE': 'md_memmove', 'MEMCMP': 'md_memcmp',
           'BSWAP16': 'md_bswap16', 'BSWAP32': 'md_bswap32', 'BSWAP64': 'md_bswap64', 'DEBUG_BREAK': 'md_debug_break'}
SPECIAL.update(ALIASES)
INTERMEDIATE = {'MD_' + k: v for k, v in ALIASES.items()}
CPP_ONLY = {'defer': 'MD_DEFER', 'ExitScopeHelp': 'md_defer_helper', 'array_size_impl': 'md_array_size_impl'}

# md_str.h: the type, the macros, and every function it declares, str_ or not
STR_TYPES = ['str_t']
STR_OBJ_MACROS = ['STR_FMT']
STR_FN_MACROS = ['STR_INIT', 'STR_LIT', 'STR_ARG', 'DEF_VAL']
STR_FUNCS = ['str_empty', 'str_beg', 'str_end', 'str_len', 'str_ptr', 'str_swap', 'str_trim_beg', 'str_trim_end',
             'str_trim', 'str_eq', 'str_eq_ignore_case', 'str_eq_n', 'str_eq_n_ignore_case', 'str_eq_cstr',
             'str_eq_cstr_ignore_case', 'str_eq_cstr_n', 'str_eq_cstr_n_ignore_case', 'str_cmp_lex', 'str_begins_with',
             'str_ends_with', 'str_count_equal_chars', 'str_count_occur_char', 'str_substr', 'str_from_cstr',
             'str_from_cstrn', 'str_join', 'str_extract_line', 'str_peek_line', 'str_skip_line', 'str_extract_i32',
             'str_extract_i64', 'str_extract_f32', 'str_extract_f64', 'str_find_char', 'str_rfind_char', 'str_find_str',
             'str_alloc', 'str_free', 'str_copy', 'str_copy_cstr', 'str_copy_cstrn', 'str_printf',
             'str_copy_to_char_buf', 'str_edit_distance', 'load_textfile', 'extract_ext', 'extract_file',
             'extract_file_path_without_ext', 'extract_folder_path', 'replace_char', 'convert_to_lower',
             'convert_to_upper']

PROSE = {'MIN', 'MAX', 'PI', 'TWO_PI', 'DEBUG', 'RELEASE', 'CLAMP', 'CONCAT', 'STRINGIFY', 'RESTRICT',
         'str_len', 'str_beg', 'str_end', 'str_ptr', 'str_alloc'}

NEW = {}
for n in FUNCS + STR_TYPES + STR_FUNCS: NEW[n] = 'md_' + n
for n in OBJ_MACROS + FN_MACROS + STR_OBJ_MACROS + STR_FN_MACROS: NEW[n] = SPECIAL.get(n, 'MD_' + n)
NEW.update(INTERMEDIATE)
FN = set(FN_MACROS) | set(STR_FN_MACROS) | set(STR_FUNCS) | {'MD_DEBUG_BREAK'}   # renamed where called or on a directive line

EXT = ('.c', '.h', '.inl', '.cpp', '.hpp', '.m', '.mm')
# Vendored, self-contained sources that do not include md_common.h: their names are their own
VENDORED = {'rmsd.c', 'rmsd.h'}
CPP_EXT = ('.cpp', '.hpp', '.mm')

IDENT = re.compile(r'[A-Za-z_]\w*')


def segments(s):
    """Split a source into (kind, text) with kind in code, comment, literal."""
    out, i, n, start = [], 0, len(s), 0
    def flush(kind_end, kind):
        nonlocal start
        if kind_end > start: out.append((kind, s[start:kind_end]))
        start = kind_end
    while i < n:
        c = s[i]
        if c == '/' and i + 1 < n and s[i + 1] == '/':
            flush(i, 'code'); j = i
            while j < n and s[j] != '\n':
                j += 2 if (s[j] == '\\' and j + 1 < n and s[j + 1] == '\n') else 1
            flush(j, 'comment'); i = j
        elif c == '/' and i + 1 < n and s[i + 1] == '*':
            flush(i, 'code'); j = s.find('*/', i + 2); j = n if j < 0 else j + 2
            flush(j, 'comment'); i = j
        elif c == 'R' and i + 1 < n and s[i + 1] == '"' and re.search(r'(?:^|[^\w])(?:u8|u|U|L)?$', s[max(0, i - 3):i]):
            m = re.match(r'R"([^()\\\s]{0,16})\(', s[i:])
            if not m:
                i += 1; continue
            flush(i, 'code'); end = s.find(')' + m.group(1) + '"', i + m.end())
            j = n if end < 0 else end + len(m.group(1)) + 2
            flush(j, 'literal'); i = j
        elif c == '"' or (c == "'" and not (i > 0 and s[i - 1].isdigit())):
            # a ' after a digit is a C++14 digit separator, not a character literal
            flush(i, 'code'); j = i + 1
            while j < n and s[j] != c and s[j] != '\n':
                j += 2 if s[j] == '\\' else 1
            j = min(j + 1, n)
            flush(j, 'literal'); i = j
        else:
            i += 1
    flush(n, 'code')
    return out


def directive_mask(s):
    """For each character, whether it is on a preprocessor directive line (continuations included)."""
    mask = bytearray(len(s)); i = 0; in_dir = False
    for line in s.splitlines(keepends=True):
        stripped = line.lstrip()
        if not in_dir and stripped.startswith('#'):
            in_dir = True
        if in_dir:
            mask[i:i + len(line)] = b'\x01' * len(line)
        in_dir = in_dir and line.rstrip('\r\n').endswith('\\')
        i += len(line)
    return mask


def rename(s, cpp):
    names = dict(NEW)
    if cpp: names.update(CPP_ONLY)
    mask = directive_mask(s)
    out, pos = [], 0
    for kind, text in segments(s):
        if kind == 'literal':
            out.append(text); pos += len(text); continue
        res, last = [], 0
        for m in IDENT.finditer(text):
            w = m.group(0)
            if w not in names: continue
            after = text[m.end():].lstrip(' \t')
            called = after.startswith('(')
            if kind == 'code':
                ok = (w not in FN) or called or mask[pos + m.start()]
            else:   # comment
                ok = w != 'defer' and (called or w not in PROSE)
            if ok:
                res.append(text[last:m.start()]); res.append(names[w]); last = m.end()
        res.append(text[last:]); out.append(''.join(res)); pos += len(text)
    return ''.join(out)


def main(argv):
    dry = '--dry-run' in argv
    cpp_headers = '--cpp-headers' in argv
    paths = [a for a in argv if not a.startswith('--')]
    files = []
    for p in paths:
        if os.path.isfile(p): files.append(p); continue
        for d, dirs, fs in os.walk(p):
            # Shader sources are GLSL/Slang, with their own DIV_UP etc.; a .h there may be shared with C (it then
            # defines its own names under #ifndef guards), so the whole directory is left alone
            dirs[:] = [x for x in dirs if x != 'shaders']
            files += [os.path.join(d, f) for f in fs if f.endswith(EXT) and f not in VENDORED]
    changed = 0
    for f in sorted(files):
        raw = open(f, 'rb').read()
        bom = raw.startswith(b'\xef\xbb\xbf')
        s = raw[3:].decode('utf-8') if bom else raw.decode('utf-8', errors='surrogateescape')
        cpp = f.endswith(CPP_EXT) or (cpp_headers and f.endswith(('.h', '.inl')))
        t = rename(s, cpp)
        if t != s:
            changed += 1
            if not dry:
                open(f, 'wb').write((b'\xef\xbb\xbf' if bom else b'') + t.encode('utf-8', errors='surrogateescape'))
    print(f'{changed} of {len(files)} files changed' + (' (dry run)' if dry else ''))


if __name__ == '__main__':
    main(sys.argv[1:])
