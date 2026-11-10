#!/usr/bin/env python3
# Prints a COMP block for unittest/frag_chem.txt, to append to it with the fragments the molecule should have (F
# lines) written by hand. The test reads the elements and coordinates only, so a structure of any source will do:
#
#   frag_chem_add.py ID structure.xyz                   an xyz file (its first frame), no dependencies
#   frag_chem_add.py ID --smiles "CN(C)c1ccc(C#N)cc1"   a conformer by RDKit (ETKDG, MMFF), with its bonds
#
# Appended without F lines, the molecule fails fragment.reference_molecules, which prints the fragments it has
# ("got"): check each against the chemistry before taking it as the F lines, and fix md_fragment where it is wrong.
import argparse, os, sys

def from_xyz(path):
    with open(path) as f:
        n = int(f.readline().split()[0])
        f.readline()
        atoms = []
        for _ in range(n):
            el, x, y, z = f.readline().split()[:4]
            atoms.append((el.capitalize(), 0, float(x), float(y), float(z)))
    return atoms, []

def from_smiles(smi):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 7
    if AllChem.EmbedMolecule(m, ps) != 0:
        sys.exit("embedding failed")
    if AllChem.MMFFHasAllMoleculeParams(m):
        AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    else:
        AllChem.UFFOptimizeMolecule(m, maxIters=2000)
    conf = m.GetConformer()
    atoms = [(a.GetSymbol(), a.GetFormalCharge(), *conf.GetAtomPosition(a.GetIdx())) for a in m.GetAtoms()]
    Chem.Kekulize(m, clearAromaticFlags=False)
    order = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), order.get(b.GetBondType(), 1), int(b.GetIsAromatic())) for b in m.GetBonds()]
    return atoms, bonds

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("id", help="at most 7 characters, unique in frag_chem.txt and ccd_chem.txt")
    ap.add_argument("xyz", nargs="?")
    ap.add_argument("--smiles")
    ap.add_argument("--note", default="", help="what the molecule covers, for the comment line")
    args = ap.parse_args()
    if len(args.id) > 7:
        sys.exit("the id is at most 7 characters")
    if args.smiles:
        atoms, bonds = from_smiles(args.smiles)
        source = args.smiles
    elif args.xyz:
        atoms, bonds = from_xyz(args.xyz)
        source = os.path.basename(args.xyz)
    else:
        sys.exit("give an xyz file or --smiles")
    print(f"COMP {args.id} {len(atoms)} {len(bonds)}")
    print(f"# {args.note + ': ' if args.note else ''}{source}")
    count = {}
    for el, q, x, y, z in atoms:
        count[el] = count.get(el, 0) + 1
        print(f"A {el}{count[el]} {el} {q} {x:.3f} {y:.3f} {z:.3f}")
    for i, j, o, ar in bonds:
        print(f"B {i} {j} {o} {ar}")
    print("# F heavy-atoms kind name, one per fragment, from the chemistry")
    print("END")

if __name__ == "__main__":
    main()
