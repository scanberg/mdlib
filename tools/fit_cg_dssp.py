#!/usr/bin/env python3
"""
Fits the virtual carbonyls of CG-DSSP (md_dssp.c) and writes them to src/md_dssp_cg.inl.

CG-DSSP assigns secondary structure to a trace of coarse grained backbone beads (the BB beads of Martini) from bonds
between its peptide groups, and a peptide group's bond direction is that of its C=O. The beads do not give it: the
peptide plane turns about the bond between two beads with the conformation. This fits where the C and O of each residue
sit, in a frame of the beads around it, as linear functions of trig features of the trace there, against atomistic
structures and their coarse grained mappings.

The frame of residue i, from its bead and its neighbors i-1, i+1:
    t = unit(x[i+1] - x[i-1])
    n = unit(x[i-1] + x[i+1] - 2 x[i], without its component along t)
    b = t x n
The features (radians): tau, the bond angle at i; theta_a, theta_b the dihedrals (i-2, i-1, i, i+1) and (i-1, i, i+1,
i+2); tau_p, tau_n the bond angles at i-1 and i+1. Residues with both dihedrals (kind 1), only theta_b (kind 2, the
second residue of a chain) or only theta_a (kind 3, the second to last) have a table each; tau_p (tau_n) is tau where
there is no i-2 (i+2), theta_a (theta_b) is 0 there. features() below is the reference for the order, which
cg_dssp_dipole in md_dssp.c follows.

Each input is a folder with aa.pdb, an atomistic structure, and cg.pdb, its mapping to Martini beads as martinize2 writes
it (-x). Their residues are matched chain by chain in order (martinize2 may write the chains in another order). The fit
adds copies of each trace with Gaussian noise on the beads, which makes the model expect a trace as it is in dynamics
rather than as it is mapped from a crystal structure.

Run from the repository root:  python3 tools/fit_cg_dssp.py <folder> [<folder> ...] [--noise 0.35] [--copies 4]
Requires numpy.
"""
import argparse
import math
import os
import sys

import numpy as np

CHAIN_BREAK = 4.6   # A, a longer bond between consecutive BB beads ends a chain
ATOMS = ("N", "CA", "C", "O")


def read_aa(path):
    """Backbone atoms of every residue with a CA, by chain in file order (first model, first altloc)"""
    res, order = {}, []
    for line in open(path):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[16] not in (" ", "A"):
            continue
        name = line[12:16].strip()
        key = (line[21], int(line[22:26]), line[26])
        if key not in res:
            res[key] = {}
            order.append(key)
        p = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        if name in ATOMS and name not in res[key]:
            res[key][name] = p
        if name in ("OT1", "OC1", "O1") and "O" not in res[key]:
            res[key]["O"] = p
    chains = {}
    for key in order:
        if "CA" in res[key]:
            chains.setdefault(key[0], []).append(res[key])
    return chains


def read_cg(path):
    """BB beads in file order: chain id and position"""
    out = []
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "BB":
            out.append((line[21], np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return out


def load(folder):
    """Beads (n, 3), atoms (n, 4, 3) with NaN where missing, and the chains as [beg, end) ranges"""
    aa = read_aa(os.path.join(folder, "aa.pdb"))
    cg = read_cg(os.path.join(folder, "cg.pdb"))
    x = np.array([p for _, p in cg])
    atoms = np.full((len(cg), 4, 3), np.nan)
    used = {}
    for i, (chain, _) in enumerate(cg):
        r = used.get(chain, 0)
        used[chain] = r + 1
        if r >= len(aa.get(chain, [])):
            raise ValueError(f"{folder}: chain {chain} of cg.pdb has more residues than aa.pdb")
        for a, name in enumerate(ATOMS):
            if name in aa[chain][r]:
                atoms[i, a] = aa[chain][r][name]
    chains, beg = [], 0
    for i in range(1, len(cg) + 1):
        if i == len(cg) or cg[i][0] != cg[i - 1][0] or np.linalg.norm(x[i] - x[i - 1]) > CHAIN_BREAK:
            chains.append((beg, i))
            beg = i
    return x, atoms, chains


def unit(v):
    return v / np.linalg.norm(v)


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, unit(p2 - p1), p3 - p2
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w))


def angle(a, b, c):
    u, v = a - b, c - b
    return math.acos(max(-1.0, min(1.0, np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))))


def features(tau, ta, tb, tp, tn, kind):
    c, s = math.cos, math.sin
    f = [1.0, c(tau), s(tau), c(tau) ** 2,
         c(tp), c(tn), c(tp) * c(tb), c(tn) * c(tb), c(tp) * c(ta), c(tn) * c(ta),
         c(tp) * s(tb), c(tn) * s(tb), c(tp) * s(ta), c(tn) * s(ta)]
    da = [c(ta), s(ta), c(2 * ta), s(2 * ta)]
    db = [c(tb), s(tb), c(2 * tb), s(2 * tb)]
    if kind == 1:
        f += da + db + [c(ta) * c(tb), s(ta) * s(tb), s(ta) * c(tb), c(ta) * s(tb)] + [c(tau) * v for v in da + db]
    elif kind == 2:
        f += db + [c(tau) * v for v in db]
    else:
        f += da + [c(tau) * v for v in da]
    return f


def residues(x, chains):
    """Per residue with a frame: (index, kind, frame rows t n b, features)"""
    for (beg, end) in chains:
        for i in range(beg + 1, end - 1):
            has_a, has_b = i - 2 >= beg, i + 2 < end
            if not (has_a or has_b):
                continue
            t = unit(x[i + 1] - x[i - 1])
            m = x[i - 1] + x[i + 1] - 2 * x[i]
            n = unit(m - np.dot(m, t) * t)
            frame = np.stack([t, n, np.cross(t, n)])
            tau = angle(x[i - 1], x[i], x[i + 1])
            ta = dihedral(x[i - 2], x[i - 1], x[i], x[i + 1]) if has_a else 0.0
            tb = dihedral(x[i - 1], x[i], x[i + 1], x[i + 2]) if has_b else 0.0
            tp = angle(x[i - 2], x[i - 1], x[i]) if has_a else tau
            tn = angle(x[i], x[i + 1], x[i + 2]) if has_b else tau
            kind = 1 if (has_a and has_b) else 2 if has_b else 3
            yield i, kind, frame, features(tau, ta, tb, tp, tn, kind)


def fit(systems, noise, copies, ridge=1e-3, seed=7):
    rng = np.random.default_rng(seed)
    X, Y = {1: [], 2: [], 3: []}, {1: [], 2: [], 3: []}
    for (x, atoms, chains) in systems:
        for c in range(1 + copies):
            xc = x + rng.normal(0.0, noise, x.shape) if c > 0 else x
            for i, kind, frame, f in residues(xc, chains):
                local = frame @ (atoms[i, 2:4] - xc[i]).T      # C, O in (t, n, b), as columns
                if np.isnan(local).any():
                    continue
                X[kind].append(f)
                Y[kind].append(local.T.reshape(6))
    W = {}
    for kind in (1, 2, 3):
        A, B = np.array(X[kind]), np.array(Y[kind])
        W[kind] = np.linalg.solve(A.T @ A + ridge * np.eye(A.shape[1]), A.T @ B)
    return W, {k: len(v) for k, v in X.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folders", nargs="+")
    ap.add_argument("--noise", type=float, default=0.35, help="Gaussian noise on the beads of the copies, A per coordinate")
    ap.add_argument("--copies", type=int, default=4, help="noisy copies of each structure")
    ap.add_argument("--out", default=os.path.join("src", "md_dssp_cg.inl"))
    args = ap.parse_args()

    systems = [load(f) for f in args.folders]
    W, counts = fit(systems, args.noise, args.copies)
    num_residues = sum(len(s[0]) for s in systems)
    names = ", ".join(os.path.basename(os.path.normpath(f)) for f in args.folders)

    lines = [
        "// Generated by tools/fit_cg_dssp.py, do not edit by hand.",
        "//",
        "// The virtual carbonyls of CG-DSSP (md_dssp.c): where the C and O of a residue sit in the frame (t, n, b) of its",
        "// bead, {C.t, C.n, C.b, O.t, O.n, O.b} = sum_f feature[f] * W[f], for Martini BB beads (the center of N, CA, C",
        "// and O as martinize2 maps them). See the script for the frame and the order of the features.",
        f"// Fit on {len(systems)} structures, {num_residues} residues ({names}), with {args.copies} copies of each",
        f"// with {args.noise} A of noise on the beads.",
        "",
    ]
    for kind, label in ((1, "both dihedrals"), (2, "no dihedral before: the second residue of a chain"), (3, "no dihedral after: the second to last")):
        rows = ",\n".join("    {" + ", ".join(f"{v: .8e}f" for v in r) + "}" for r in W[kind])
        lines.append(f"// Kind {kind}, {label} ({counts[kind]} samples)")
        lines.append(f"static const float cg_dssp_carbonyl_w{kind}[{W[kind].shape[0]}][6] = {{\n{rows}\n}};")
        lines.append("")
    with open(args.out, "w") as f:
        f.write("\n".join(lines))
    print(f"wrote {args.out}: {num_residues} residues, samples per kind {counts}")


if __name__ == "__main__":
    sys.exit(main())
