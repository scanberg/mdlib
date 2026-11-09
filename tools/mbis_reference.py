# Reference MBIS charges and pro-atoms, an implementation independent of md_mbis: PySCF's density (from the Molden
# file) on PySCF's own molecular grid (Treutler radial, Lebedev angular, Becke partition), and the fixed point of
# Verstraelen et al., JCTC 12, 3894 (2016) in numpy, from the same initial guess. Gives the values unittest/test_mbis.c
# asserts (water, OH, He), converged in the grid: levels 3, 5 and 7 agree to 1e-6.
#
#   pip install pyscf
#   python3 tools/mbis_reference.py test_data/molden/h2o_ccpvdz.molden
import sys, numpy as np
from pyscf.tools import molden
from pyscf import dft
from pyscf.dft import numint

def nshell(z):
    return int(np.searchsorted(np.array([2, 10, 18, 36, 54, 86, 118]), z)) + 1

def mbis(path, level=5, thresh=1e-10, maxiter=5000):
    mol, e, C, occ, _, _ = molden.load(path)
    if isinstance(C, tuple) or (isinstance(C, np.ndarray) and C.ndim == 3):
        dm = (C[0] * occ[0]) @ C[0].T + (C[1] * occ[1]) @ C[1].T
    else:
        dm = (C * occ) @ C.T
    g = dft.gen_grid.Grids(mol); g.level = level; g.prune = None; g.build()
    ao = numint.eval_ao(mol, g.coords)
    rho = numint.eval_rho(mol, ao, dm)
    w = g.weights
    Z = mol.atom_charges(); R = mol.atom_coords()
    dist = np.linalg.norm(g.coords[None, :, :] - R[:, None, :], axis=2)   # [A, P]
    pars = []
    for A, z in enumerate(Z):
        n = nshell(z); S0 = 2.0 * z
        alpha = (2.0 / S0) ** (1.0 / (n - 1)) if n > 1 else 1.0
        pars.append([[z / n, S0 * alpha ** i] for i in range(n)])
    for it in range(maxiter):
        terms = [[N * S**3 / (8 * np.pi) * np.exp(-S * dist[A]) for N, S in pars[A]] for A in range(len(Z))]
        rho0 = sum(t for ts in terms for t in ts)
        rho0 = np.maximum(rho0, 1e-300)
        new = []
        change = 0.0
        for A in range(len(Z)):
            pa = []
            for i, (N, S) in enumerate(pars[A]):
                t = w * rho * terms[A][i] / rho0
                Nn = t.sum(); rbar = (t * dist[A]).sum() / Nn
                pa.append([Nn, 3.0 / rbar])
                change = max(change, abs(Nn - N), abs(3.0 / rbar - S) / S)
            new.append(pa)
        pars = new
        if change < thresh: break
    pop = np.array([sum(N for N, S in p) for p in pars])
    return mol, Z - pop, pars, it, (w * rho).sum()

if __name__ == "__main__":
    for path in sys.argv[1:]:
        for level in (3, 5, 7):
            mol, q, pars, it, ne = mbis(path, level)
            print(path.split('/')[-1], "level", level, "iters", it, "N", "%.8f" % ne, "q", " ".join("%.6f" % x for x in q))
        for A, p in enumerate(pars):
            print("  atom", A, " ".join("N=%.6f sigma=%.6f" % (N, 1 / S) for N, S in p))
