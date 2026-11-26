"""Inputs for md_basis_census_bench: a copy of a VeloxChem file's root datasets (the geometry, the electron counts and
the rest of its core block, nothing of its SCF) under another basis set's name, one per basis set.

    pip install h5py
    python basis_census_inputs.py [--library dir] source.h5 out_dir CC-PVTZ AUG-CC-PVTZ DEF2-QZVP ...

The files are out_dir/<stem>.<basis>.h5, <stem> the source's name up to its first dot. mdlib reads the basis set itself
from VeloxChem's basis library beside the executable (bin/basis), so a name must be one of those files. With --library
(that folder), each basis file is also copied into out_dir with only the molecule's elements, as <NAME>-TRIM, and the
files name that copy: mdlib's parser refuses a whole file that has h functions anywhere, which VeloxChem's CC-PVQZ and
AUG-CC-PVQZ have for the transition metals.
"""

import os
import sys

import h5py
import numpy as np

SYMBOLS = ('X H HE LI BE B C N O F NE NA MG AL SI P S CL AR K CA SC TI V CR MN FE CO NI CU ZN GA GE AS SE BR KR '
           'RB SR Y ZR NB MO TC RU RH PD AG CD IN SN SB TE I XE').split()


def trim(src, dst, ident, symbols):
    out, keep = [], True
    with open(src) as f:
        for line in f:
            tok = line.split()
            if tok[:1] == ['@BASIS_SET']:
                line = '@BASIS_SET %s\n' % ident
            elif tok[:1] == ['@ATOMBASIS']:
                keep = tok[1].upper() in symbols
            if keep:
                out.append(line)
            if tok[:1] == ['@END']:
                keep = True
    with open(dst, 'w') as f:
        f.writelines(out)


def main():
    args = sys.argv[1:]
    library = None
    if len(args) >= 2 and args[0] == '--library':
        library, args = args[1], args[2:]
    if len(args) < 3:
        print(__doc__)
        return 1
    src, out_dir, names = args[0], args[1], args[2:]
    os.makedirs(out_dir, exist_ok=True)
    with h5py.File(src, 'r') as f:
        root = {k: np.array(f[k]) for k in f.keys() if isinstance(f[k], h5py.Dataset)}
    symbols = {SYMBOLS[int(round(z))] for z in root['nuclear_charges']}
    stem = os.path.basename(src).split('.')[0]
    for name in names:
        ident = name.upper()
        if library:
            ident += '-TRIM'
            trim(os.path.join(library, name.upper()), os.path.join(out_dir, ident), ident, symbols)
        path = os.path.join(out_dir, '%s.%s.h5' % (stem, name.lower()))
        with h5py.File(path, 'w') as f:
            for k, v in root.items():
                f[k] = np.array([ident.encode()]) if k == 'basis_set' else v
        print(path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
