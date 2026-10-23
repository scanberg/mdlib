#!/usr/bin/env python3
"""Renames the vocabulary of a system in C and C++ sources, as mdlib renamed it.

    component      -> residue          (md_component_* -> md_residue_*, sys->component -> sys->residue, comp_idx -> res_idx)
    structure      -> molecule         (md_structure_* -> md_molecule_*, sys->structure -> sys->molecule)
    atom type      -> particle type    (md_atom_type_* -> md_particle_type_*, MD_ATOM_TYPE_* -> MD_PARTICLE_TYPE_*)

Run it over a branch made before the rename, then merge or rebase it onto one made after:

    python3 tools/rename_vocabulary.py <file or directory>...           rewrite the files in place
    python3 tools/rename_vocabulary.py --dry-run <file or directory>...  list what would change
    python3 tools/rename_vocabulary.py --prose <file or directory>...    the words of the comments as well

mdlib and VIAMD were renamed with --prose over mdlib's src, unittest, benchmark and examples and VIAMD's src and test.

What it renames:
  - Identifiers, in code and in comments: every identifier with 'component' or 'comp' as one of its parts ('_'
    separated, or a capitalised part of a CamelCase name), and every one with the parts 'atom' 'type'. The names the
    file formats give their own things stay: mmCIF's comp_id and chem_comp, LAMMPS's and GROMACS's atom types
    (md_lammps_atom_type_t, md_itp_atomtype_t), and words which only begin with comp (compute, compare, ...).
  - The members of md_system_t, read through '.' or '->': component, structure.
  - The molecule API by name (md_structure_*, md_util_system_infer_structures, MD_MATCH_LEVEL_STRUCTURE, ...). Other
    structures (secondary structure, a GPU structure, an electronic structure) are not molecules and stay.
  - With --prose, the words in comments of the files given: component(s) -> residue(s), atom type(s) -> particle
    type(s); a line which speaks of connected components, the Chemical Component Dictionary or the components of a
    vector is left alone.
Strings are never touched: the names a file or a workspace stores keep.

A local name is renamed in its scope (a top level declaration, or a function with its body: from a line which begins
with a name to the next): where its short form is taken there (res_idx), it is spelled out (residue_idx), and where both
are it is left as it is and reported, to be renamed by hand. A member read through '.' or '->' is always renamed. The
columns of a renamed line (runs of two or more spaces between text) are kept where one space allows.
"""

import os
import re
import sys

SOURCE_EXT = ('.c', '.h', '.cpp', '.hpp', '.inl', '.cc')

# Identifiers which have 'comp' or 'component' as a part but are not residues
KEEP = {
    # mmCIF's names for the residue name and the Chemical Component Dictionary
    'comp_id', 'label_comp_id', 'auth_comp_id', 'ATOM_SITE_LABEL_COMP_ID', 'ATOM_SITE_AUTH_COMP_ID', 'chem_comp',
    # The components of a vector, of a transition density, of an attribute's shape; compressed files
    'ncomp', 'n_comp', 'w_comp', 'h5_comp', 'comp_len',
    'ElectronicStructureTransitionDensityComponent', 'ElectronicStructureNtoComponent', 'transition_density_component',
    'nto_component', 'electronic_structure_transition_density_component_str', 'electronic_structure_nto_component_str',
    'no_components', 'slice_keeps_components', 'components_must_be_stated',
    # OpenGL and Vulkan
    'GL_DEPTH_COMPONENT', 'GL_DEPTH_COMPONENT24', 'GL_DEPTH_COMPONENT32F', 'DEPTH_COMPONENT32F', 'VkColorComponentFlags',
    # The components of a value of an attribute (md_attribute_format_t), a float[3] has three; a connected component of
    # the bond graph. The member of md_system_t is renamed as a member (MEMBER).
    'components', 'component', 'component_size',
    # The fragments of a SMILES string, connected components of a graph
    'num_components', 'count_components', 'collect_component', 'match_component',
    # Chemical components (CCD definitions) of the tests of the chemistry
    'ref_comp_t', 'compare_component', 'build_component', 'find_component', 'read_components', 'component_util',
    # The format's own atom types
    'md_lammps_atom_type_t',
}

# Files (by name) whose 'comp' is not a residue: shader programs, compressed files, a composition
KEEP_FILES = {
    'md_gto.c', 'md_gto.h', 'md_gto_iso.c', 'md_gto_iso.h', 'md_gto_int.c', 'md_topo.c', 'md_topo.h', 'md_flow.c',
    'md_gfx.c', 'md_gfx.h', 'md_smiles.c', 'md_smiles.h', 'md_chem.c', 'md_chem.h', 'md_attributes.c', 'md_attributes.h',
    'md_attributes_internal.h', 'md_pot.c', 'md_pot.h', 'test_vlx.c', 'test_chem.c', 'test_ring.c', 'component_util.h',
    'test_attributes.c', 'quantum_chemistry.cpp', 'md_qm.c', 'md_qm.h',
    # VIAMD: the completions of the filter input
    'filter_input.cpp',
}
# Identifiers kept in one file only: md_match's composition of a query (z, count pairs)
KEEP_IN_FILE = {
    'md_match.c': {'comp', 'comp_off', 'comp_len'},
    # The component of a transition density or a natural transition orbital
    'representation.cpp': {'component'},
}
# The files whose comments speak of residues when they say component (--prose); elsewhere a component is of a vector,
# of a graph, of a value of an attribute
PROSE_FILES = {
    # mdlib
    'md_system.h', 'md_system.c', 'md_types.h', 'md_util.h', 'md_util.c', 'md_match.h', 'md_match.c',
    'md_script_functions.inl', 'md_script_complete.inl', 'md_pdb.c', 'md_pdb.h', 'md_gro.c', 'md_gro.h', 'md_mmcif.c',
    'md_mmcif.h', 'md_tpr.c', 'md_tpr.h', 'md_lammps.c', 'md_lammps.h', 'md_itp.c', 'md_itp.h', 'md_gl.c', 'md_gl.h',
    'md_dssp.c', 'md_dssp.h', 'md_dssp_cg.inl', 'md_hbond.c', 'md_hbond.h', 'md_monomer.c', 'md_martini_db.inl',
    'md_atomic_infer.c', 'md_xyz.c', 'md_ase_traj.c', 'md_molden.c', 'md_trexio.c', 'md_nonbonded.c', 'md_contact.c',
    'md_pucker.c', 'md_water.c',
    'test_util.c', 'test_match.c', 'test_pdb.c', 'test_gro.c', 'test_mmcif.c', 'test_tpr.c', 'test_lipid.c',
    'test_script.c', 'test_filter.c', 'test_monomer.c', 'system_invariants.h', 'test_system_invariants.c',
    'test_secondary_structure.c', 'test_itp.c', 'test_lammps.c', 'test_h5md.c', 'test_hbond.c', 'test_water.c',
    'test_bonds.c', 'test_xyz.c', 'test_infer_atomic_number.c', 'test_contact.c', 'test_nonbonded.c',
    # VIAMD
    'viamd.h', 'viamd.cpp', 'script_suggestions.cpp', 'script_suggestions.h', 'color_utils.cpp', 'color_utils.h',
    'cartoon_glyphs.cpp', 'cartoon_glyphs.h', 'dataset.cpp', 'contacts.cpp', 'representation.cpp', 'pucker.cpp',
    'ramachandran.cpp', 'test_script_suggestions.cpp', 'test_cartoon_glyphs.cpp', 'test_color_utils.cpp',
}
# Directories never entered
SKIP_DIRS = {'core', 'shaders', 'ext', '.git', 'build', 'deprecated', 'app', 'gfx'}

# Explicit renames, ahead of the rules by parts
EXPLICIT = {
    # The type of an atom (an accessor of the atoms) rather than of a type
    'md_atom_type_idx': 'md_atom_particle_type_idx',
    'md_monomer_from_comp_name': 'md_monomer_from_residue_name',
    # Molecules
    'md_structure_data_t': 'md_molecule_data_t',
    'md_structure_t': 'md_molecule_t',
    'md_structure_count': 'md_molecule_count',
    'md_structure_extract': 'md_molecule_extract',
    'md_structure_atom_slot': 'md_molecule_atom_slot',
    'md_structure_atom_parent': 'md_molecule_atom_parent',
    'md_util_system_infer_structures': 'md_util_system_infer_molecules',
    'MD_UTIL_INFER_STRUCTURE_BIT': 'MD_UTIL_INFER_MOLECULE_BIT',
    'MD_UTIL_INFER_UNWRAP_STRUCTURE_BIT': 'MD_UTIL_INFER_UNWRAP_MOLECULE_BIT',
    'md_util_unwrap_structure': 'md_util_unwrap_molecule',
    'MD_MATCH_LEVEL_STRUCTURE': 'MD_MATCH_LEVEL_MOLECULE',
    'check_structure_invariants': 'check_molecule_invariants',
    'structure_size_of_atom': 'molecule_size_of_atom',
}

IDENT = re.compile(r'[A-Za-z_][A-Za-z0-9_]*')


def split_parts(ident):
    """'md_component_idx_t' -> ['md', '_', 'component', '_', 'idx', '_', 't']; CamelCase parts split too"""
    out = []
    for i, chunk in enumerate(ident.split('_')):
        if i:
            out.append('_')
        if chunk:
            out.extend(re.findall(r'[A-Z]+(?![a-z])|[A-Z]?[a-z0-9]+|[0-9]+', chunk) or [chunk])
    return out


PART = {
    'component': 'residue', 'components': 'residues', 'Component': 'Residue', 'Components': 'Residues',
    'COMPONENT': 'RESIDUE', 'COMPONENTS': 'RESIDUES',
    'comp': 'res', 'comps': 'residues', 'Comp': 'Res', 'COMP': 'RES', 'COMPS': 'RESIDUES',
}


PART_API = {'comp': 'residue', 'Comp': 'Residue', 'COMP': 'RESIDUE'}


def rename_ident(ident, spelled_out=False):
    if ident in KEEP:
        return ident
    if ident in EXPLICIT:
        return EXPLICIT[ident]
    parts = split_parts(ident)
    if ''.join(parts) != ident:
        return ident
    api = spelled_out or ident.startswith(('md_', 'MD_'))
    changed = False
    out = []
    i = 0
    while i < len(parts):
        p = parts[i]
        # atom type -> particle type, in any case
        if p.lower() == 'atom' and i + 2 < len(parts) + 1:
            nxt = [q for q in parts[i + 1:i + 3]]
            if len(nxt) == 2 and nxt[0] == '_' and nxt[1].lower() in ('type', 'types') or \
               len(nxt) >= 1 and nxt[0] in ('Type', 'Types'):
                out.append({'atom': 'particle', 'Atom': 'Particle', 'ATOM': 'PARTICLE'}.get(p, p))
                changed = changed or p in ('atom', 'Atom', 'ATOM')
                i += 1
                continue
        if p in PART:
            # The API spells the word out: md_system_chain_find_by_residue_idx; a local abbreviates it: res_idx
            out.append(PART_API.get(p, PART[p]) if api else PART[p])
            changed = True
        else:
            out.append(p)
        i += 1
    return ''.join(out) if changed else ident


# The members of md_system_t, by the receiver they are read through; a GPU operation has a structure of its own
MEMBER = re.compile(r'(?<![A-Za-z0-9_])((?:[A-Za-z_][A-Za-z0-9_]*(?:\[[^\]]*\])?|\))\s*(?:\.|->)\s*)(component|structure)\b')
MEMBER_NEW = {'component': 'residue', 'structure': 'molecule'}
# The declaration of the member in md_system_t
SYSTEM_MEMBER_DECL = re.compile(r'(md_(?:structure|molecule|component|residue)_data_t\s+)(structure|component)(\s*;)')
MEMBER_KEEP_RECEIVER = {'op', 'in_op', 'cmd', 'draw_op'}
DESIGNATED = re.compile(r'(?<![A-Za-z0-9_\]\)])(\.)(component|structure)(\s*=)')

PROSE_SKIP = re.compile(r'\bevent\b|EventType|this build|compiled in|whichever component|Initialize component|no component and no event|\blanes?\b|Cartesian|triclinic basis|% 3|three components|per component scalars|already reduced|current component|across components|over the component|to the component rather|Sweep \d|component is smallest|Components of two atoms|component\(|connected component|Chemical Component|Component Dictionary|COMPONENT WISE|component wise|component-wise|[Cc]omponents? of (a|the|one) (vector|float|value|velocit|position|coordinate|shape)|COMPONENTS of|[xyzw][ -]component|vector component|UI component|component wise|Component wise|COMPONENT WISE|components/')
PROSE = [
    (re.compile(r'\bcomponent \(residue\)'), 'residue'), (re.compile(r'\bComponent \(residue\)', re.I), 'Residue'),
    (re.compile(r'\bcomponents \(residues\)'), 'residues'),
    (re.compile(r'\ban atom type'), 'a particle type'), (re.compile(r'\bAn atom type'), 'A particle type'),
    (re.compile(r'\bcomponent types\b'), 'residue types'), (re.compile(r'\bcomponent type\b'), 'residue type'),
    (re.compile(r'\bComponent types\b'), 'Residue types'), (re.compile(r'\bComponent type\b'), 'Residue type'),
    (re.compile(r'\bCOMPONENT TYPES\b'), 'RESIDUE TYPES'), (re.compile(r'\bCOMPONENT TYPE\b'), 'RESIDUE TYPE'),
    (re.compile(r'\bcomponents\b'), 'residues'), (re.compile(r'\bcomponent\b'), 'residue'),
    (re.compile(r'\bComponents\b'), 'Residues'), (re.compile(r'\bComponent\b'), 'Residue'),
    (re.compile(r'\bCOMPONENTS\b'), 'RESIDUES'), (re.compile(r'\bCOMPONENT\b'), 'RESIDUE'),
    (re.compile(r'\batom types\b'), 'particle types'), (re.compile(r'\batom type\b'), 'particle type'),
    (re.compile(r'\bAtom types\b'), 'Particle types'), (re.compile(r'\bAtom type\b'), 'Particle type'),
    (re.compile(r'\bATOM TYPES\b'), 'PARTICLE TYPES'), (re.compile(r'\bATOM TYPE\b'), 'PARTICLE TYPE'),
]


def lex(text):
    """Splits C/C++ into (kind, text) segments: 'code', 'comment', 'string'"""
    segs = []
    i, n = 0, len(text)
    start = 0
    kind = 'code'

    def push(k, a, b):
        if b > a:
            if segs and segs[-1][0] == k:
                segs[-1] = (k, segs[-1][1] + text[a:b])
            else:
                segs.append((k, text[a:b]))

    while i < n:
        c = text[i]
        if text.startswith('//', i):
            push('code', start, i)
            j = text.find('\n', i)
            j = n if j < 0 else j
            push('comment', i, j)
            i = start = j
        elif text.startswith('/*', i):
            push('code', start, i)
            j = text.find('*/', i + 2)
            j = n if j < 0 else j + 2
            push('comment', i, j)
            i = start = j
        elif c == 'R' and text.startswith('R"', i) and (i == 0 or not (text[i - 1].isalnum() or text[i - 1] == '_')):
            m = re.match(r'R"([^()\\\s]{0,16})\(', text[i:])
            if not m:
                i += 1
                continue
            push('code', start, i)
            end = text.find(')' + m.group(1) + '"', i + m.end())
            j = n if end < 0 else end + len(m.group(1)) + 2
            push('string', i, j)
            i = start = j
        elif c == '"' or (c == "'" and not (i > 0 and (text[i - 1].isalnum()) and text[i - 1] not in 'uUL')):
            # A char literal, but not the digit separator of C++14 (1'000)
            push('code', start, i)
            j = i + 1
            while j < n and text[j] != c and text[j] != '\n':
                j += 2 if text[j] == '\\' else 1
            j = min(j + 1, n)
            push('string', i, j)
            i = start = j
        else:
            i += 1
    push('code', start, n)
    return segs


def rewrite(path, text, prose, report):
    name = os.path.basename(path)
    keep_here = KEEP_IN_FILE.get(name, set())
    segs = lex(text)
    # The scopes a local name lives in: a top level declaration, a function with its body, each from a line which
    # begins with a name at its first column to the next. By the layout rather than the braces, which the branches
    # of a preprocessor conditional leave uneven.
    unit_of = []            # [segment] -> (offset, unit) of each identifier of a code segment, in order
    unit_names = [set()]
    at_line_start = True
    for kind, seg in segs:
        if kind != 'code':
            unit_of.append(None)
            at_line_start = seg.endswith('\n') if kind == 'comment' else False
            continue
        marks = []
        for m in re.finditer(r'\n|[A-Za-z_][A-Za-z0-9_]*', seg):
            t = m.group(0)
            if t == '\n':
                at_line_start = True
                continue
            if at_line_start and m.start() > 0 and seg[m.start() - 1] == '\n' or at_line_start and m.start() == 0:
                unit_names.append(set())
            at_line_start = False
            unit_names[-1].add(t)
            marks.append((m.start(), len(unit_names) - 1))
        if seg and seg[-1] != '\n':
            at_line_start = False
        unit_of.append(marks)
    skip_tokens = name in KEEP_FILES
    collisions = set()

    cur = {'marks': [], 'k': 0}
    seg_now = {'text': ''}

    def unit_at(offset):
        marks = cur['marks']
        k = cur['k']
        while k + 1 < len(marks) and marks[k + 1][0] <= offset:
            k += 1
        cur['k'] = k
        return marks[k][1] if marks else 0

    def ident_sub(m):
        tok = m.group(0)
        if tok in keep_here:
            return tok
        if skip_tokens and tok not in EXPLICIT and not tok.startswith(('md_', 'MD_')):
            return tok
        new = rename_ident(tok)
        # A member is read through its struct, which no local name can collide with
        before = seg_now['text'][:m.start()].rstrip()
        if before.endswith(('.', '->')):
            return new
        if new != tok and tok not in EXPLICIT and not tok.startswith(('md_', 'MD_')):
            names = unit_names[unit_at(m.start())]
            if new in names and rename_ident(new) == new:
                # The word spelled out: residue_idx where res_idx is taken
                alt = rename_ident(tok, spelled_out=True)
                if alt != new and alt not in names:
                    collisions.add((tok, alt, 'renamed'))
                    return alt
                collisions.add((tok, new, 'kept'))
                return tok
        return new

    def member_sub(m):
        recv = m.group(1)
        base = re.match(r'[A-Za-z_][A-Za-z0-9_]*', recv)
        if base and base.group(0) in MEMBER_KEEP_RECEIVER:
            return m.group(0)
        if name in ('md_gfx.c', 'md_gfx.h', 'main.cpp') and m.group(2) == 'structure' and base and base.group(0) not in ('sys', 'mol', 'mold'):
            return m.group(0)
        return recv + MEMBER_NEW[m.group(2)]

    # The names the code of the file renames, which its comments speak of by name: 'the residues comps'
    renamed_locals = {}
    for (kind, seg), marks in zip(segs, unit_of):
        if kind == 'code':
            for t in set(IDENT.findall(seg)):
                if t in keep_here or t in KEEP or (skip_tokens and not t.startswith(('md_', 'MD_'))):
                    continue
                n = rename_ident(t)
                if n != t and re.search(r'(^|_)(comp|comps)($|_)', t):
                    renamed_locals[t] = n

    out = []
    for (kind, seg), marks in zip(segs, unit_of):
        if kind == 'string':
            out.append(seg)
            continue
        if kind == 'code':
            # Identifiers first, at the offsets the scopes were found at; then the members of the system
            cur['marks'] = marks
            cur['k'] = 0
            seg_now['text'] = seg
            seg = IDENT.sub(ident_sub, seg)
            seg = MEMBER.sub(member_sub, seg)
            seg = DESIGNATED.sub(lambda m: m.group(1) + MEMBER_NEW[m.group(2)] + m.group(3), seg)
            seg = SYSTEM_MEMBER_DECL.sub(lambda m: m.group(1) + MEMBER_NEW[m.group(2)] + m.group(3), seg)
        else:
            # Comments: identifiers which are of the API (md_..., MD_..., the explicit renames), and the words
            def comment_ident(m):
                tok = m.group(0)
                if tok in renamed_locals:
                    return renamed_locals[tok]
                if tok in EXPLICIT or tok.startswith(('md_', 'MD_')) or '_' in tok:
                    if skip_tokens and not tok.startswith(('md_', 'MD_')) and tok not in EXPLICIT:
                        return tok
                    if tok in keep_here:
                        return tok
                    return rename_ident(tok)
                return tok
            seg = re.sub(r'(?<![A-Za-z0-9_])(sys|mol)(->|\.)(component|structure)\b', lambda m: m.group(1) + m.group(2) + MEMBER_NEW[m.group(3)], seg)
            seg = IDENT.sub(comment_ident, seg)
            if prose and name in PROSE_FILES:
                lines = seg.split('\n')
                for k, line in enumerate(lines):
                    if PROSE_SKIP.search(line):
                        continue
                    for rx, rep in PROSE:
                        line = rx.sub(rep, line)
                    lines[k] = line
                seg = '\n'.join(lines)
        out.append(seg)
    for tok, new, what in sorted(collisions):
        if what == 'renamed':
            report.append(f'{path}: {tok} -> {new}, as its short form is taken where it is')
        else:
            report.append(f'{path}: {tok} collides with a name where it is, left as it is')
    return ''.join(out)


GAP = re.compile(r'(?<=\S)( {2,})(?=\S)')


def keep_columns(old, new):
    """Where a renamed line had columns (runs of two or more spaces between text), the text after each run is put back
    in the column it was in, as far as one space allows"""
    old_lines = old.split('\n')
    new_lines = new.split('\n')
    if len(old_lines) != len(new_lines):
        return new
    out = []
    for a, b in zip(old_lines, new_lines):
        if a == b:
            out.append(b)
            continue
        ga = [m for m in GAP.finditer(a)]
        gb = [m for m in GAP.finditer(b)]
        if not ga or len(ga) != len(gb):
            out.append(b)
            continue
        line = ''
        pos = 0
        for ma, mb in zip(ga, gb):
            line += b[pos:mb.start()]
            # The column the text after the gap began at in the old line
            want = ma.end() - len(line)
            line += ' ' * max(1, want)
            pos = mb.end()
        line += b[pos:]
        out.append(line)
    return '\n'.join(out)


def files_of(paths):
    for p in paths:
        if os.path.isfile(p):
            yield p
            continue
        for root, dirs, files in os.walk(p):
            dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
            for f in sorted(files):
                if f.endswith(SOURCE_EXT):
                    yield os.path.join(root, f)


def main(argv):
    dry = '--dry-run' in argv
    prose = '--prose' in argv
    paths = [a for a in argv if not a.startswith('--')]
    if not paths:
        print(__doc__)
        return 1
    report = []
    changed = 0
    for path in files_of(paths):
        with open(path, encoding='utf-8', errors='surrogateescape', newline='') as f:
            text = f.read()
        new = keep_columns(text, rewrite(path, text, prose, report))
        if new != text:
            changed += 1
            if dry:
                print(path)
            else:
                with open(path, 'w', encoding='utf-8', errors='surrogateescape', newline='') as f:
                    f.write(new)
    for line in report:
        print(line, file=sys.stderr)
    print(f'{changed} files {"would change" if dry else "changed"}', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
