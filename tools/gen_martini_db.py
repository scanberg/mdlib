#!/usr/bin/env python3
"""Generates src/md_martini_db.inl: the residues of the Martini force fields, with the mass, size and particle kind of
each bead and the bonds within each residue, and the links (bonds) between consecutive residues of the polymers and of
the molecules of several residues, for md_util_system_infer_atom_types and the coarse grained bond inference.

Every value is read from the force field files themselves:
- masses from the [ atoms ] of a molecule where it gives them, else from the [ atomtypes ] of the main itp
- the radius of a bead from the sigma of the Lennard-Jones interaction of its type with itself, as md_tpr gives it:
  half the distance of the minimum of the potential, 2^(1/6) sigma / 2
- a virtual site is a bead without mass; one without Lennard-Jones is a charge site (MD_PARTICLE_VIRTUAL_SITE)
- bonds are the chemical ones, as GROMACS flags them (IF_CHEMBOND): [ bonds ] of type 1, 2, 3, 4, 5, 7, 8 and
  [ constraints ] of type 1, and a virtual site to the first atom it is constructed from, as md_tpr bonds it.
  Harmonic potentials (type 6: elastic networks) and restraints are not bonds.

Sources (the directory of each is an argument, see --help):
  Martini 3   martini_v3.0.0.itp, ions, solvents and sterols v1 and the v2 lipidome of
              github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters (ITPs/);
              martini_v3.0.0_phospholipids_v1, small_molecules_v1, sugars_v1 of martini_v300.zip,
              martini_v3.0.0_sugars_v2.itp and small_molecules_v1 (2021-11-01, 90 molecules) from cgmartini.nl;
              the small (SW) and tiny (TW) water, of a bead of their atom types; proteins from vermouth
              (martini3001); RNA from github.com/DanYev/reForge (the itp its martinize_rna.py writes, --nucleic);
              polymers from polyply (data/martini3)
  Martini 2   martini_v2.2.itp, martini_v2.2P.itp, martini_v2.0_ions, solvents and sugars and the lipidome from
              cgmartini.nl; proteins from vermouth (martini22, elnedyn22, martini22p); DNA and RNA from the
              topologies martinize-dna.py and martinize-nucleotide.py (cgmartini.nl) write (--nucleic); polymers
              from polyply (data/martini2)

A residue name may stand for different molecules in different versions (POPC of Martini 2 and 3 have the same beads
of other sizes); the versions of a residue say which. Where two sources of the same version disagree on a residue,
the first wins and the other is reported.
"""

import argparse, collections, glob, math, os, re, sys

M2, M2P, M3 = 1, 2, 4
KINDS = {'OTHER': 'MD_COMPONENT_KIND_OTHER', 'AMINO_ACID': 'MD_COMPONENT_KIND_AMINO_ACID', 'NUCLEOTIDE': 'MD_COMPONENT_KIND_NUCLEOTIDE',
         'WATER': 'MD_COMPONENT_KIND_WATER', 'ION': 'MD_COMPONENT_KIND_ION', 'LIPID': 'MD_COMPONENT_KIND_LIPID'}
WATER_NAMES = {'W', 'PW', 'WF', 'SW', 'TW', 'BMW'}
CHEM_BONDS = {1, 2, 3, 4, 5, 7, 8}
K = 0.5 * 2 ** (1 / 6) * 10.0   # sigma (nm) to the radius (A)


def strip(line):
    return line.split(';')[0].strip()


MACROS = {}   # #define NAME VALUE of the itps read so far (the named bond types of the Martini 3 lipidome)


def read_macros(path):
    for raw in open(path, errors='ignore'):
        s = strip(raw)
        if s.startswith('#define'):
            w = s.split(None, 2)
            if len(w) == 3:
                MACROS[w[1]] = w[2]


def sections(path, defines=()):
    """Yields (section, fields) of an itp, following #ifdef/#ifndef/#else/#endif"""
    skip = []
    sec = None
    read_macros(path)
    for raw in open(path, errors='ignore'):
        s = strip(raw)
        if not s:
            continue
        if s.startswith('#'):
            w = s.split()
            if w[0] == '#ifdef':
                skip.append(w[1] not in defines)
            elif w[0] == '#ifndef':
                skip.append(w[1] in defines)
            elif w[0] == '#else' and skip:
                skip[-1] = not skip[-1]
            elif w[0] == '#endif' and skip:
                skip.pop()
            continue
        if any(skip):
            continue
        m = re.match(r'\[\s*(\w+)\s*\]', s)
        if m:
            sec = m.group(1).lower()
            continue
        yield sec, s.split()


class ForceField:
    """The atom types of a main itp: their masses and the sigma of each with itself"""
    def __init__(self, *paths):
        self.mass = {}
        self.self_lj = {}
        self.own_lj = {}
        comb = 1
        for path in paths:
            for sec, f in sections(path):
                if sec == 'defaults':
                    comb = int(f[1])
                elif sec == 'atomtypes':
                    # name [at.num] mass charge ptype V W
                    try:
                        if len(f) >= 7 and f[4] in ('A', 'S', 'V', 'D'):
                            mass, v, w = float(f[2]), float(f[5]), float(f[6])
                        else:
                            mass, v, w = float(f[1]), float(f[4]), float(f[5])
                    except (ValueError, IndexError):
                        continue
                    self.mass[f[0]] = mass
                    self.own_lj[f[0]] = (comb, v, w)
                elif sec == 'nonbond_params' and len(f) >= 5 and f[0] == f[1]:
                    self.self_lj[f[0]] = (comb, float(f[3]), float(f[4]))

    def sigma(self, t):
        comb, v, w = self.self_lj.get(t) or self.own_lj.get(t, (1, 0.0, 0.0))
        if comb == 1:
            return (w / v) ** (1 / 6) if v > 0 and w > 0 else 0.0
        return v if w > 0 else 0.0

    def bead(self, t, mass, vsite):
        if t not in self.mass:
            raise KeyError(t)
        m = self.mass[t] if mass is None else mass
        s = self.sigma(t)
        if s > 0:
            return (round(m, 3), round(K * s, 3), 'BEAD')
        if vsite and m == 0:
            return (0.0, 0.0, 'VSITE')     # A charge site
        return (round(m, 3), 1.0, 'BEAD')  # No Lennard-Jones: the radius of element 0, as md_tpr gives it


Residue = collections.namedtuple('Residue', 'name kind beads bonds')


def molecules(path, defines=()):
    """The molecule types of an itp: residues of (name, type, mass, vsite) and bonds"""
    mols = []
    mol = None
    for sec, f in sections(path, defines):
        if sec == 'moleculetype':
            mol = {'name': f[0], 'atoms': [], 'bonds': [], 'vsites': set()}
            mols.append(mol)
        elif mol is None:
            continue
        elif sec == 'atoms' and len(f) >= 5:
            mass = None
            if len(f) >= 8:
                try:
                    mass = float(f[7])
                except ValueError:
                    pass
            mol['atoms'].append({'type': f[1], 'resnr': f[2], 'resname': f[3], 'name': f[4], 'mass': mass})
        elif sec in ('bonds', 'constraints') and len(f) >= 3:
            ft = f[2] if f[2] not in MACROS else MACROS[f[2]].split()[0]
            try:
                a, b, ft = int(f[0]), int(f[1]), int(ft)
            except ValueError:
                continue
            if (sec == 'bonds' and ft in CHEM_BONDS) or (sec == 'constraints' and ft == 1):
                mol['bonds'].append((a - 1, b - 1))
        elif sec.startswith('virtual_sites') or sec.startswith('dummies'):
            # A site, and the first atom it is constructed from (the third field of an N-site, after the function),
            # which md_tpr bonds it to
            try:
                site, first = int(f[0]) - 1, int(f[2] if sec == 'virtual_sitesn' else f[1]) - 1
            except (ValueError, IndexError):
                continue
            mol['vsites'].add(site)
            mol['bonds'].append((site, first))
    return mols


def residue_runs(mol):
    """The residue of each atom of a molecule type: runs of atoms of one residue number and name, as GRO, PDB and TPR
    readers split them (a residue whose atoms are not together is several)"""
    run, prev = [], None
    for a in mol['atoms']:
        key = (a['resnr'], a['resname'])
        run.append(run[-1] + (key != prev) if run else 0)
        prev = key
    return run


def residue_links(mol):
    """The bonds of a molecule type between consecutive residues of it (a disaccharide, a lipid of several residues), as
    links: (resnames, bead, bead of the next residue, atom resnames)"""
    out = []
    run = residue_runs(mol)
    for i, j in mol['bonds']:
        if not (0 <= i < len(mol['atoms']) and 0 <= j < len(mol['atoms'])):
            continue
        a, b = mol['atoms'][i], mol['atoms'][j]
        ka, kb = run[i], run[j]
        if abs(ka - kb) != 1:
            continue
        if ka > kb:
            a, b = b, a
        out.append((frozenset((a['resname'], b['resname'])), a['name'], b['name'], {a['resname']}, {b['resname']}))
    return out


def residues_of(mol, ff, kind, report):
    """Splits a molecule type into its residues"""
    out = []
    groups = collections.OrderedDict()
    for i, r in enumerate(residue_runs(mol)):
        groups.setdefault(r, []).append(i)
    for idx in groups.values():
        resname = mol['atoms'][idx[0]]['resname']
        local = {i: k for k, i in enumerate(idx)}
        beads = []
        try:
            for i in idx:
                a = mol['atoms'][i]
                beads.append((a['name'],) + ff.bead(a['type'], a['mass'], i in mol['vsites']))
        except KeyError as e:
            report.append('%s %s: unknown type %s' % (mol['name'], resname, e))
            continue
        bonds = sorted({tuple(sorted((local[a], local[b]))) for a, b in mol['bonds'] if a in local and b in local})
        k = kind
        if k == 'SOLVENT':
            k = 'WATER' if resname in WATER_NAMES else 'OTHER'
        if k == 'ION' and len(idx) > 1:
            k = 'OTHER'     # A molecular ion of several beads (an ionic liquid, phosphate)
        out.append(Residue(resname, k, tuple(beads), tuple(bonds)))
    return out


def protein_residues(ffdir, ff, report):
    """The residues of a vermouth force field (aminoacids.ff): bead names and types, masses where it gives them, and
    the bonds within each (of its [ bonds ] and [ constraints ], under any #meta ifdef, and of a virtual site to the
    first bead it is constructed from). Atoms are referred to by name or by number."""
    out = []
    macros = {}
    mol = None
    sec = None
    want_name = False
    for raw in open(os.path.join(ffdir, 'aminoacids.ff'), errors='ignore'):
        s = strip(raw)
        if not s or s.startswith('#'):
            continue
        m = re.match(r'\[\s*(\S+?)\s*\]', s)
        if m:
            sec = m.group(1)
            want_name = sec == 'moleculetype'
            if sec in ('link', 'modification'):
                mol = None
            continue
        f = s.split()
        if sec == 'macros':
            macros[f[0]] = f[1].strip('"')
        elif want_name:
            mol = {'name': f[0], 'beads': [], 'pairs': []}
            out.append(mol)
            want_name = False
        elif mol is None:
            continue
        elif sec == 'atoms' and len(f) >= 5:
            t = f[1][1:] if f[1].startswith('$') else f[1]
            t = macros.get(t, f[1])
            mass = None
            if len(f) >= 8:
                try:
                    mass = float(f[7])
                except ValueError:
                    pass
            mol['beads'].append((f[4], t, mass))
        elif sec in ('bonds', 'constraints') and len(f) >= 3:
            try:
                ft = int(f[2])
            except ValueError:
                continue
            if (sec == 'bonds' and ft in CHEM_BONDS) or (sec == 'constraints' and ft == 1):
                mol['pairs'].append((f[0], f[1]))
        elif sec.startswith('virtual_sites') and len(f) >= 3:
            # vermouth writes an N-site as 'site atoms... -- funct', and the others as GROMACS does: 'site atoms... funct'
            mol['pairs'].append((f[0], f[1]))
    res = []
    for mol in out:
        try:
            beads = tuple((n,) + ff.bead(t, m, False) for n, t, m in mol['beads'])
        except KeyError as e:
            report.append('%s: unknown type %s' % (mol['name'], e))
            continue
        names = [b[0] for b in mol['beads']]
        def index(ref):
            return int(ref) - 1 if ref.isdigit() else names.index(ref) if ref in names else -1
        bonds = set()
        for a, b in mol['pairs']:
            i, j = index(a), index(b)
            if i < 0 or j < 0 or i >= len(names) or j >= len(names) or i == j:
                report.append('%s: bond %s-%s of no beads of it' % (mol['name'], a, b))
                continue
            bonds.add((min(i, j), max(i, j)))
        res.append(Residue(mol['name'], 'AMINO_ACID', beads, tuple(sorted(bonds))))
    return res


def ref_offset(ref):
    """An atom of a block or a link: its name, and the residue it is in relative to the first ('+', '>' next, '-', '<'
    previous)"""
    n = len(ref) - len(ref.lstrip('+>'))
    m = len(ref) - len(ref.lstrip('-<'))
    return ref[max(n, m):], n - m


def polymer_library(path, ff, report):
    """The blocks (residues) and links of a polyply force field: a .ff file (vermouth's format, blocks and [ link ]s of
    atoms by name) or an .itp of one block, whose atom numbers past its last stand for the atoms of the residues after
    it. Returns the residues and the links, (resnames, bead, bead of the next residue, atom resnames) for each bond
    between consecutive residues."""
    blocks, links = [], []
    cur = None
    sec = None
    want_name = False
    for raw in open(path, errors='ignore'):
        s = strip(raw)
        if not s or s.startswith('#'):
            continue
        m = re.match(r'\[\s*(\S+?)\s*\]', s)
        if m:
            sec = m.group(1)
            if sec == 'moleculetype':
                want_name = True
            elif sec == 'link':
                cur = {'link': True, 'resnames': None, 'pairs': [], 'atoms': {}}
                links.append(cur)
            elif sec == 'modification':
                cur = None
            continue
        f = s.split()
        if want_name:
            cur = {'link': False, 'name': f[0], 'beads': [], 'pairs': []}
            blocks.append(cur)
            want_name = False
        elif cur is None:
            continue
        elif cur['link'] and f[0] == 'resname' and len(f) >= 2:
            cur['resnames'] = set(f[1].strip('"').split('|'))
        elif sec == 'atoms' and cur['link']:
            m = re.search(r'"resname"\s*:\s*"([^"]+)"', s)
            if m:
                cur['atoms'][f[0]] = set(m.group(1).split('|'))
        elif sec == 'atoms' and len(f) >= 5:
            mass = None
            if len(f) >= 8:
                try:
                    mass = float(f[7])
                except ValueError:
                    pass
            cur['beads'].append((f[4], f[1], mass))
        elif sec in ('bonds', 'constraints') and len(f) >= 3:
            try:
                ft = int(f[2])
            except ValueError:
                continue
            if (sec == 'bonds' and ft in CHEM_BONDS) or (sec == 'constraints' and ft == 1):
                cur['pairs'].append((f[0], f[1]))
        elif sec.startswith('virtual_sites') and not cur['link'] and len(f) >= 3:
            cur['pairs'].append((f[0], f[1]))

    residues = []
    for b in blocks:
        names = [x[0] for x in b['beads']]
        try:
            beads = tuple((n,) + ff.bead(t, m, False) for n, t, m in b['beads'])
        except KeyError as e:
            report.append('%s (%s): unknown type %s' % (b['name'], os.path.basename(path), e))
            continue
        bonds = set()
        for x, y in b['pairs']:
            if x.isdigit() and y.isdigit():
                # An .itp block: the atoms of the residues after it continue its numbering
                i, j = int(x) - 1, int(y) - 1
                ri, rj = i // len(names), j // len(names)
                if ri == rj == 0:
                    bonds.add((min(i, j), max(i, j)))
                elif abs(ri - rj) == 1:
                    a, c = (i, j) if ri < rj else (j, i)
                    links.append({'link': True, 'resnames': {b['name']}, 'pairs': [(names[a % len(names)], '+' + names[c % len(names)])], 'atoms': {}})
                continue
            (nx, ox), (ny, oy) = ref_offset(x), ref_offset(y)
            if ox == oy == 0 and nx in names and ny in names and nx != ny:
                i, j = names.index(nx), names.index(ny)
                bonds.add((min(i, j), max(i, j)))
        residues.append(Residue(b['name'], 'OTHER', beads, tuple(sorted(bonds))))

    out = []
    for l in links:
        for x, y in l['pairs']:
            (nx, ox), (ny, oy) = ref_offset(x), ref_offset(y)
            if abs(ox - oy) != 1:
                continue
            if ox > oy:
                nx, ny = ny, nx
            out.append((frozenset(l['resnames'] or ()), nx, ny, l['atoms'].get(nx), l['atoms'].get(ny)))
    return residues, out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vermouth', required=True, help='vermouth/data/force_fields of vermouth 0.15')
    ap.add_argument('--m3-lipids', required=True, help='ITPs/ of M3-Lipid-Parameters')
    ap.add_argument('--m3', required=True, help='martini_v300/ of martini_v300.zip')
    ap.add_argument('--m3-sugars-v2', required=True, help='martini_v3.0.0_sugars_v2.itp')
    ap.add_argument('--m3-small-molecules', required=True, help='martini_v3.0.0_small_molecules_v1.itp of 2021-11-01 (90 molecules)')
    ap.add_argument('--m3-rna', required=True, help='reforge/martini of reForge')
    ap.add_argument('--m2', required=True, help='directory of martini_v2.2.itp, martini_v2.2P.itp, martini_v2.0_{ions,solvents,sugars}.itp')
    ap.add_argument('--m2-lipids', required=True, help='the Martini 2 lipidome, an itp per lipid in any subdirectory')
    ap.add_argument('--polyply', required=True, help='polyply/data of polyply 1.8 (its martini3 and martini2 polymer libraries)')
    ap.add_argument('--nucleic', nargs='*', default=[], help='version:itp of nucleic acids, e.g. M2:dna.itp M3:rna.itp')
    ap.add_argument('-o', '--output', required=True)
    args = ap.parse_args()

    report = []
    for p in glob.glob(os.path.join(args.m3_lipids, '*ffbonded*.itp')):
        read_macros(p)
    m3_main = os.path.join(args.m3_lipids, 'martini_v3.0.0.itp')
    ff3 = ForceField(m3_main)
    ff3_rna = ForceField(m3_main, os.path.join(args.m3_rna, 'martinize_rna_v3.0.0_itps', 'martini_v3.0.0_rna.itp'))
    ff2 = ForceField(os.path.join(args.m2, 'martini_v2.2.itp'))
    ff2p = ForceField(os.path.join(args.m2, 'martini_v2.2P.itp'))
    dna_itp = os.path.join(args.m2, 'martini_v2.1-dna.itp')
    ff2_dna = ForceField(dna_itp) if os.path.exists(dna_itp) else ff2

    # (versions, forcefield, kind, itp) in priority order: where two of one version disagree, the first wins
    sources = []
    lip3 = sorted(glob.glob(os.path.join(args.m3_lipids, '*.itp')))
    for p in lip3:
        b = os.path.basename(p)
        if b in ('martini_v3.0.0.itp', 'martini_v3.0.0_ffbonded_v2.itp'):
            continue
        kind = 'ION' if 'ions' in b else 'SOLVENT' if 'solvents' in b else 'OTHER' if 'hydrocarbons' in b else 'LIPID'
        sources.append((M3, ff3, kind, p))
    sources.append((M3, ff3, 'LIPID', os.path.join(args.m3, 'martini_v3.0.0_phospholipids_v1.itp')))
    sources.append((M3, ff3, 'OTHER', args.m3_small_molecules))
    sources.append((M3, ff3, 'OTHER', os.path.join(args.m3, 'martini_v3.0.0_small_molecules_v1.itp')))
    sources.append((M3, ff3, 'OTHER', args.m3_sugars_v2))
    sources.append((M3, ff3, 'OTHER', os.path.join(args.m3, 'martini_v3.0.0_sugars_v1.itp')))
    sources.append((M3, ff3, 'OTHER', os.path.join(args.m3, 'martini_v3.0.0_nucleobases_v1.itp')))
    sources.append((M2 | M2P, ff2, 'ION', os.path.join(args.m2, 'martini_v2.0_ions.itp')))
    sources.append((M2 | M2P, ff2, 'SOLVENT', os.path.join(args.m2, 'martini_v2.0_solvents.itp')))
    sources.append((M2 | M2P, ff2, 'OTHER', os.path.join(args.m2, 'martini_v2.0_sugars.itp')))
    sources.append((M2, ff2, 'SOLVENT', os.path.join(args.m2, 'martini_v2.2.itp')))      # W and WF
    sources.append((M2P, ff2p, 'SOLVENT', os.path.join(args.m2, 'martini_v2.2P.itp')))   # PW
    for p in sorted(glob.glob(os.path.join(args.m2_lipids, '**', '*.itp'), recursive=True)):
        sources.append((M2 | M2P, ff2, 'LIPID', p))

    residues = []   # (versions, Residue, source)
    links = []      # (versions, resnames, bead, bead of the next, atom resnames)
    for versions, ff, kind, path in sources:
        for mol in molecules(path):
            links.extend((versions,) + l for l in residue_links(mol))
            for r in residues_of(mol, ff, kind, report):
                residues.append((versions, r, os.path.basename(path)))
                # An ion as gmx genion names it: the residue and the atom after the molecule type (NA NA, NA+ NA+)
                if r.kind == 'ION' and len(r.beads) == 1 and (r.name, r.beads[0][0]) != (mol['name'], mol['name']):
                    alias = Residue(mol['name'], 'ION', ((mol['name'],) + r.beads[0][1:],), ())
                    residues.append((versions, alias, os.path.basename(path)))

    # The small (SW) and tiny (TW) water of Martini 3: atom types of martini_v3.0.0.itp, of no molecule type there
    for name in ('SW', 'TW'):
        residues.append((M3, Residue(name, 'WATER', ((name,) + ff3.bead(name, None, False),), ()), 'martini_v3.0.0.itp'))

    for name, versions in (('martini3001', M3), ('martini22', M2), ('elnedyn22', M2), ('martini22p', M2P)):
        ff = ff3 if versions == M3 else ff2p if versions == M2P else ff2
        for r in protein_residues(os.path.join(args.vermouth, name), ff, report):
            residues.append((versions, r, name))

    # The polymers of polyply: their monomers, and the links between them
    for sub, versions, ff in (('martini3', M3, ff3), ('martini2', M2 | M2P, ff2)):
        for path in sorted(glob.glob(os.path.join(args.polyply, sub, '*.ff')) + glob.glob(os.path.join(args.polyply, sub, '*.itp'))):
            if os.path.basename(path) in ('aminoacids.ff', 'modifications.ff', 'DNA_M2.ff'):
                continue
            res, lnk = polymer_library(path, ff, report)
            for r in res:
                residues.append((versions, r, os.path.basename(path)))
            links.extend((versions,) + l for l in lnk)

    for spec in args.nucleic:
        v, path = spec.split(':', 1)
        versions = {'M2': M2 | M2P, 'M3': M3}[v]
        ff = ff3_rna if versions == M3 else ff2_dna
        for mol in molecules(path):
            for r in residues_of(mol, ff, 'NUCLEOTIDE', report):
                residues.append((versions, r, os.path.basename(path)))

    # Merge: the same residue in several versions is one entry; the first of a version wins over a later one
    merged = collections.OrderedDict()   # Residue -> versions
    owner = {}                           # (name, bead names, version bit) -> Residue
    for versions, r, src in residues:
        key_names = (r.name, tuple(sorted(b[0] for b in r.beads)))
        keep = 0
        for bit in (M2, M2P, M3):
            if not versions & bit:
                continue
            k = key_names + (bit,)
            if k in owner:
                if owner[k] != r:
                    report.append('%s (%s) in %s: kept the earlier definition' % (r.name, ' '.join(b[0] for b in r.beads), src))
                continue
            owner[k] = r
            keep |= bit
        if keep:
            merged[r] = merged.get(r, 0) | keep

    # Sorted by name, as martini_db_find searches them
    entries = sorted(((r, v) for r, v in merged.items()), key=lambda rv: (rv[0].name, -rv[1]))

    # The links, between the residues (of the version) that have their beads
    has_bead = collections.defaultdict(set)   # (name, bead) -> versions
    for r, v in entries:
        for b in r.beads:
            has_bead[(r.name, b[0])].add(v)
    def of_version(name, bead, versions):
        return any(v & versions for v in has_bead.get((name, bead), ()))
    names = sorted({r.name for r, _ in entries})
    link_rows = {}
    for versions, resnames, a, b, ra, rb in links:
        cand = sorted(resnames) if resnames else names
        for x in cand:
            if (ra and x not in ra) or not of_version(x, a, versions):
                continue
            for y in cand:
                if (rb and y not in rb) or not of_version(y, b, versions):
                    continue
                key = (x, a, y, b)
                link_rows[key] = link_rows.get(key, 0) | versions
    link_rows = sorted(link_rows.items())

    bead_defs = collections.OrderedDict()
    def bead_index(b):
        if b not in bead_defs:
            bead_defs[b] = len(bead_defs)
        return bead_defs[b]

    res_beads, res_bonds, rows = [], [], []
    for r, v in entries:
        bead_off, bond_off = len(res_beads), len(res_bonds)
        res_beads.extend(bead_index(b) for b in r.beads)
        res_bonds.extend(r.bonds)
        rows.append((r.name, v, r.kind, bead_off, len(r.beads), bond_off, len(r.bonds)))
    assert len(res_beads) < 65536 and len(res_bonds) < 65536 and len(bead_defs) < 65536
    assert all(n < 256 for *_, n, _, _ in [(0, 0, 0, 0, row[4], 0, 0) for row in rows])

    vname = lambda v: ' | '.join(n for bit, n in ((M2, 'MARTINI_2'), (M2P, 'MARTINI_2P'), (M3, 'MARTINI_3')) if v & bit)
    with open(args.output, 'w') as o:
        o.write('// Generated by tools/gen_martini_db.py from the Martini force field files (see there). Do not edit.\n')
        o.write('// %d residues, %d beads of %d kinds, %d bonds, %d links.\n\n' % (len(rows), len(res_beads), len(bead_defs), len(res_bonds), len(link_rows)))
        o.write('static const martini_db_bead_t martini_db_beads[] = {\n')
        for (name, mass, radius, kind), i in bead_defs.items():
            o.write('    {"%s", %.3ff, %.3ff, %s},\n' % (name, mass, radius, 'MD_PARTICLE_VIRTUAL_SITE' if kind == 'VSITE' else 'MD_PARTICLE_BEAD'))
        o.write('};\n\n')
        o.write('static const uint16_t martini_db_residue_beads[] = {')
        for i, b in enumerate(res_beads):
            o.write(('\n    ' if i % 24 == 0 else ' ') + '%d,' % b)
        o.write('\n};\n\n')
        o.write('static const uint8_t martini_db_residue_bonds[][2] = {')
        for i, (a, b) in enumerate(res_bonds):
            o.write(('\n    ' if i % 16 == 0 else ' ') + '{%d,%d},' % (a, b))
        o.write('\n};\n\n')
        o.write('// Sorted by name\n')
        o.write('static const martini_db_residue_t martini_db_residues[] = {\n')
        for name, v, kind, bo, bn, do, dn in rows:
            o.write('    {"%s", %s, %s, %d, %d, %d, %d},\n' % (name, vname(v), KINDS[kind], bo, bn, do, dn))
        o.write('};\n\n')
        o.write('// The bonds between consecutive residues: of the polymers of polyply, and of molecules of several residues\n')
        o.write('static const martini_db_link_t martini_db_links[] = {\n')
        for (x, a, y, b), v in link_rows:
            o.write('    {"%s", "%s", "%s", "%s", %s},\n' % (x, a, y, b, vname(v)))
        o.write('};\n')

    for line in report:
        print('note:', line, file=sys.stderr)
    print('%d residues, %d bead entries, %d bead kinds, %d bonds, %d links' % (len(rows), len(res_beads), len(bead_defs), len(res_bonds), len(link_rows)), file=sys.stderr)


if __name__ == '__main__':
    main()
