# Reference ELF / LOL at points along lines, from PySCF's AOs (double) and Multiwfn's formulas (function.f90,
# ELF_LOL: Kohout-Savin spin form, Fc_pol = 4.557799872, D += 1e-5 when ELF_addminimal = 1): an implementation
# independent of mdlib's, for unittest/elf_ref_data (gto_localization.matches_independent).
#
#   pip install pyscf
#   python3 tools/elf_reference.py test_data/molden/<name>.molden unittest/elf_ref_data/<name>.ref.txt
#
# Columns: x y z (bohr), rho, ELF, ELF with D + 1e-5, LOL; 5 lines from the molecule's centre out to 12 bohr,
# 49 points each, 0.25 bohr apart.
import sys, numpy as np
from pyscf.tools import molden
from pyscf import dft
FC_POL = 4.557799872
def terms(mol, C, occ, pts):
    ao = mol.eval_gto('GTOval_sph_deriv1' if not mol.cart else 'GTOval_cart_deriv1', pts)   # (4, n, nao)
    phi = np.einsum('xpa,ai->xpi', ao, C)            # (4, n, nmo)
    rho = np.einsum('pi,i->p', phi[0]**2, occ)
    grad = 2.0 * np.einsum('pi,xpi,i->xp', phi[0], phi[1:], occ)
    tau = 0.5 * np.einsum('xpi,i->p', phi[1:]**2, occ)
    return rho, (grad**2).sum(0), tau
def elf_lol(ra, ga, ta, rb, gb, tb, addmin):
    D = ta + tb
    D = D - np.where(ra > 0, ga / np.where(ra > 0, ra, 1) / 8, 0) - np.where(rb > 0, gb / np.where(rb > 0, rb, 1) / 8, 0)
    Dh = FC_POL * (ra**(5/3) + rb**(5/3))
    if addmin: D = D + 1e-5
    elf = 1.0 / (1.0 + (D / Dh)**2)
    t = Dh / (ta + tb)
    lol = t / (1 + t)
    return elf, lol
f = sys.argv[1]
mol, e, C, occ, irrep, spins = molden.load(f)
xyz = mol.atom_coords()
c = xyz.mean(0)
dirs = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,1],[1,-2,0.5]], float)
dirs /= np.linalg.norm(dirs, axis=1)[:, None]
r = np.arange(0, 12.0001, 0.25)
pts = np.concatenate([c + np.outer(r, d) for d in dirs])
if isinstance(C, tuple) or (np.ndim(C) == 3):
    Ca, Cb = C; oa, ob = occ
else:
    spins = np.asarray(spins) if spins is not None else None
    if spins is not None and len(spins) == C.shape[1] and (spins == 1).any():
        Ca, oa = C[:, spins == 0], occ[spins == 0]; Cb, ob = C[:, spins == 1], occ[spins == 1]
    else:
        Ca, oa = C, occ / 2; Cb, ob = C, occ / 2
ra, ga, ta = terms(mol, Ca, oa, pts)
rb, gb, tb = terms(mol, Cb, ob, pts)
out = []
for addmin in (0, 1):
    out.append(elf_lol(ra, ga, ta, rb, gb, tb, addmin))
np.savetxt(sys.argv[2], np.column_stack([pts, ra + rb, out[0][0], out[1][0], out[0][1]]),
           header='x y z (bohr) rho elf elf_addminimal lol', fmt='%.10e')
print(f, 'nelec', mol.nelectron, 'alpha/beta occ', oa.sum(), ob.sum(), 'points', len(pts), 'centre', c)
