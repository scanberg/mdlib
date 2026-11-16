#!/usr/bin/env python3
# Inputs for md_topo_gto_bench beyond test_data: RHF densities (pyscf) written as molden files, each molecule
# axis-aligned and randomly rotated, so that results which depend on the cube lattice's alignment show up.
#
#   water        cc-pVDZ   3 attractors, 2 BCPs
#   dimer        cc-pVDZ   water dimer: 6 attractors, 5 BCPs (one of them the hydrogen bond, small eigenvalues)
#   benzene      6-31G*    12 attractors, 12 BCPs, 1 RCP
#   cubane       6-31G*    16 attractors, 20 BCPs, 6 RCPs, 1 CCP
#
# These are the inputs the simplification of the certified search was checked on (claude/topo-gto-simplify):
# with every test that decided nothing on its own removed, each run gave the same critical points and edges,
# bit for bit (the result hashes md_topo_gto_bench prints), on the CPU and on the GPU sweep.
#
#   python3 topo_gto_inputs.py <out dir>                 (pyscf; written with 2.6.2)
#   md_topo_gto_bench <out dir>/*.molden
#   md_topo_gto_bench --dump a <out dir>/cubane.molden   ... then on another build: --compare a.cpu.txt

import sys, os
import numpy as np
from pyscf import gto, scf
from pyscf.tools import molden

def ring(n_atoms, r_c, r_h, z=0.0):
    out = []
    for k in range(n_atoms):
        t = 2.0 * np.pi / n_atoms * k
        out.append(('C', r_c * np.cos(t), r_c * np.sin(t), z))
        out.append(('H', r_h * np.cos(t), r_h * np.sin(t), z))
    return out

def cube(r_c, r_h):
    out = []
    for sx in (-1, 1):
        for sy in (-1, 1):
            for sz in (-1, 1):
                out.append(('C', r_c * sx, r_c * sy, r_c * sz))
                out.append(('H', r_h * sx, r_h * sy, r_h * sz))
    return out

# Angstrom
MOLECULES = {
    'water':   ([('O', 0.0, 0.0, 0.1173), ('H', 0.0, 0.7572, -0.4692), ('H', 0.0, -0.7572, -0.4692)], 'cc-pvdz'),
    'dimer':   ([('O', -1.551007, -0.114520, 0.0), ('H', -1.934259, 0.762503, 0.0), ('H', -0.599677, 0.040712, 0.0),
                 ('O', 1.350625, 0.111469, 0.0), ('H', 1.680398, -0.373741, -0.758561), ('H', 1.680398, -0.373741, 0.758561)], 'cc-pvdz'),
    'benzene': (ring(6, 1.397, 2.481), '6-31g*'),
    'cubane':  (cube(0.78, 1.409), '6-31g*'),           # C-C 1.56, C-H 1.09 along the body diagonals
}

def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(out_dir, exist_ok=True)
    q, _ = np.linalg.qr(np.random.default_rng(7).normal(size=(3, 3)))
    R = q * np.sign(np.linalg.det(q))                    # a fixed proper rotation
    for name, (atoms, basis) in MOLECULES.items():
        for rotated in (False, True):
            xyz = [(a, *(np.array([x, y, z]) @ R.T if rotated else (x, y, z))) for a, x, y, z in atoms]
            mol = gto.M(atom=[(a, (x, y, z)) for a, x, y, z in xyz], basis=basis, unit='Angstrom', verbose=0)
            mf = scf.RHF(mol).run()
            path = os.path.join(out_dir, name + ('_rot' if rotated else '') + '.molden')
            molden.from_scf(mf, path)
            print(f'{path}: {mol.natm} atoms, {mol.nao} AOs, E = {mf.e_tot:.10f}')

if __name__ == '__main__':
    main()
