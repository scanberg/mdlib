#!/usr/bin/env python3
"""Ports code to the mdlib folder split: src/ in folders of layers (core, system, analysis, qm, io, gfx, script).

    python3 port_folder_split.py [--move] <dir or file> ...

What it does to the C, C++ and Objective-C sources under the paths given:
  - includes of the headers moved into core, which is included with its folder: <md_geom.h> -> <core/md_geom.h>,
    the same for md_unitcell.h and md_sort.h
  - includes through a relative path into src/: "../src/md_attributes_internal.h" -> "../src/system/..."
  - the functions renamed on the way: md_geom_unwrap_system -> md_system_unwrap, md_geom_unwrap_molecule(state, mol)
    -> md_geom_unwrap on the coordinates of the state (a call it cannot rewrite is reported),
    md_nb_potential_init_from_tpr -> md_tpr_nb_potential_init, md_nb_forcefield_init_from_tpr -> md_tpr_nb_forcefield_init
With --move, run in an mdlib checkout, it also moves (git mv) the files still directly in src/ into their folders: for
a branch that added files to src/ before the split. A file it does not know is listed, to be placed by hand (and added
to the list of its layer in CMakeLists.txt).

Includes of the other headers, <md_pdb.h>, stay as they are: every folder is on the include path.
"""
import os
import re
import subprocess
import sys

LAYERS = {
    'core':     'unitcell geom sort',
    'system':   'types system system_fields system_state atom bond hierarchy backbone ring lipid martini martini_db infer '
                'infer_internal attributes attributes_internal run group assembly topology_build gmx_restypes fragment '
                'chem dssp dssp_cg nonbonded element monomer smiles smiles_write match',
    'analysis': 'contact hbond pucker mcs',
    'qm':       'qm qm_ct gto gto_int gto_iso gto_iso_gpu topo flow',
    'io':       'io amber_energy amber_netcdf amber_netcdf4 amber_prmtop amber_rst7 ase_traj csv cube dcd edr gro h5md '
                'hdf5 hdf5_direct itp lammps mmcif molden namdbin ndx pdb pot psf tpr trexio trr vlx xdr xtc xvg xyz',
    'gfx':      'gl vector_graphics',
    'script':   'script script_functions script_complete filter vial',
}
LAYER_OF = {'md_' + n: layer for layer, names in LAYERS.items() for n in names.split()}
INTO_CORE = ['md_' + n for n in LAYERS['core'].split()]

SOURCE_EXT = ('.c', '.h', '.inl', '.cpp', '.hpp', '.cc', '.cxx', '.m', '.mm')

RENAMES = [
    (r'\bmd_geom_unwrap_system\b', 'md_system_unwrap'),
    (r'\bmd_nb_potential_init_from_tpr\b', 'md_tpr_nb_potential_init'),
    (r'\bmd_nb_forcefield_init_from_tpr\b', 'md_tpr_nb_forcefield_init'),
]


def layer_of(name):
    return LAYER_OF.get(re.sub(r'\.(c|h|inl|m)$', '', name))


def port_text(text):
    # <md_geom.h> or "md_geom.h" -> <core/md_geom.h>
    for stem in INTO_CORE:
        text = re.sub(r'(#\s*include\s*)[<"]' + stem + r'\.h[>"]', r'\1<core/' + stem + '.h>', text)

    # "../src/md_x.h" -> "../src/<layer>/md_x.h"
    def relative(m):
        name = m.group(3)
        layer = layer_of(name)
        return m.group(0) if layer is None else f'{m.group(1)}{m.group(2)}src/{layer}/{name}'
    text = re.sub(r'(#\s*include\s*["<])((?:\.\./)*(?:[\w.-]+/)*?)src/(md_[\w]+\.(?:h|c|inl))', relative, text)

    for old, new in RENAMES:
        text = re.sub(old, new, text)

    # md_geom_unwrap_molecule(&state, &mol) and (state_ptr, mol_ptr): on the coordinates of the state, along the tree
    text = re.sub(r'\bmd_geom_unwrap_molecule\(\s*&(\w+)\s*,\s*&(\w+)\s*\)',
                  r'md_geom_unwrap((float*)\1.xyz, sizeof(md_float3_t), \2.atom_idx, \2.parent_idx, \2.count, &\1.unitcell)', text)
    text = re.sub(r'\bmd_geom_unwrap_molecule\(\s*(\w+)\s*,\s*(\w+)\s*\)',
                  r'md_geom_unwrap((float*)\1->xyz, sizeof(md_float3_t), \2->atom_idx, \2->parent_idx, \2->count, &\1->unitcell)', text)
    return text


def port_file(path):
    with open(path, 'rb') as f:
        raw = f.read()
    try:
        text = raw.decode('utf-8')
    except UnicodeDecodeError:
        return False
    new = port_text(text)
    if new != text:
        with open(path, 'wb') as f:
            f.write(new.encode('utf-8'))
    for n, line in enumerate(new.splitlines(), 1):
        if 'md_geom_unwrap_molecule' in line:
            print(f'{path}:{n}: md_geom_unwrap_molecule(state, mol) is md_geom_unwrap((float*)state->xyz, '
                  f'sizeof(md_float3_t), mol->atom_idx, mol->parent_idx, mol->count, &state->unitcell)')
    return new != text


def move_files(root):
    files = subprocess.run(['git', 'ls-files', 'src'], cwd=root, capture_output=True, text=True, check=True).stdout.split()
    unknown = []
    for f in files:
        d, name = os.path.split(f)
        if d != 'src':
            continue
        layer = layer_of(name)
        if layer is None:
            unknown.append(f)
            continue
        os.makedirs(os.path.join(root, 'src', layer), exist_ok=True)
        subprocess.run(['git', 'mv', f, f'src/{layer}/{name}'], cwd=root, check=True)
        print(f'moved {f} -> src/{layer}/{name}')
    for f in unknown:
        print(f'{f}: not of the split, to be placed in a folder by hand (and in CMakeLists.txt)')


def main(argv):
    move = '--move' in argv
    paths = [a for a in argv if a != '--move']
    if not paths:
        print(__doc__)
        return 1
    if move:
        move_files(paths[0] if os.path.isdir(paths[0]) else '.')
    changed = 0
    for p in paths:
        if os.path.isfile(p):
            changed += port_file(p)
            continue
        for dirpath, dirnames, filenames in os.walk(p):
            dirnames[:] = [d for d in dirnames if d not in ('.git', 'build', 'ext')]
            for name in filenames:
                if name.endswith(SOURCE_EXT):
                    changed += port_file(os.path.join(dirpath, name))
    print(f'{changed} files changed')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
