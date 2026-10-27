# mdlib Script Reference

The mdlib scripting language is a small, declarative language for describing **what to measure** in a molecular
system: selections of atoms, geometric quantities, distributions and volumes. You write a handful of assignments,
mdlib compiles them against a loaded system and evaluates them over the trajectory.

This document is the reference for every procedure the language provides. It lives next to the implementation
(`src/md_script.c`, `src/md_script_functions.inl`) and is meant to be the single source of truth for both the
in-application help and the wiki.

- [Quick start](#quick-start)
- [Language](#language): [statements](#statements-and-comments), [literals](#literals), [ranges](#ranges),
  [arrays](#arrays-and-subscripts), [operators](#operators), [named arguments](#named-arguments)
- [Types and units](#types-and-units)
- [Contexts: `in` and `out`](#contexts-in-and-out)
- [What becomes a property](#what-becomes-a-property)
- [Procedure index](#procedure-index)
- Reference: [Selectors](#selectors), [Bond sets](#bond-sets), [Properties](#properties), [Geometry](#geometry),
  [Math](#math), [Linear algebra and constructors](#linear-algebra-and-constructors),
  [Combining and reshaping](#combining-and-reshaping), [Data from outside the structure](#data-from-outside-the-structure)
- [Implicit conversions](#implicit-conversions)
- [Worked examples](#worked-examples)
- [Known issues](#known-issues)

> **Conventions used in this file.** Every procedure has exactly one `###` heading whose text is the procedure
> name (so the anchor is `#name`), an HTML comment `<!-- proc ... -->` with machine-readable metadata, a
> `text` block with one signature per overload, and at least one ```` ```mdscript ```` example. Signatures read
> `name(arg: type, ...) -> type [unit]`. All indices are **1-based** and all ranges are **inclusive**.
> Lengths are in Ångström (Å), angles in radians, unless stated otherwise.

---

## Quick start

```mdscript
# A property is an assignment whose value changes over the trajectory.
d  = distance(1, 5);                                  # Å between atom 1 and atom 5
phi = dihedral(1, 2, 3, 4);                           # radians

# The same measurement on every residue: 'in' evaluates the left side once per residue.
res_d = distance(1, 2) in residue(:);                 # float[num_residues]

# Selections are values too; keep them in variables and reuse them.
backbone_c = element("C") and backbone();
n_bb = count(backbone_c);

# Radial distribution function, cutoff 10 Å
g = rdf(element("C"), element("O"), 10.0);
```

Each assignment ends with a semicolon. Text after `#` is a comment. Results that vary over the trajectory
(or are distributions/volumes) are published by name, e.g. `res_d` above is available as `script/res_d`.

---

## Language

### Statements and comments

A script is a sequence of **assignments** terminated by `;`. There is no control flow — the language is
declarative and every expression is evaluated for every frame.

```mdscript
# comments run to the end of the line
a = distance(1, 2);          # statement
b = a * 2.0 +                # statements may span lines
    1.0;
```

Identifiers start with a letter or `_` and continue with letters, digits or `_`. The words `in`, `of`, `out`,
`and`, `or`, `xor`, `not` are keywords and cannot be used as identifiers or as parameter names. An identifier
must be defined before it is used.

**Multiple assignment** unpacks an array or vector of known length into several identifiers:

```mdscript
{a, b} = {distance(1, 2), distance(1, 3)};      # two properties
{lin, plan, iso} = shape_weights(all());        # the three components of a float[3]
{first, second} = residue(1:2);                 # two selections, one residue each
```

The number of identifiers on the left must match the length on the right.

### Literals

| Kind | Examples | Type |
|---|---|---|
| Integer | `1`, `42` | `int` |
| Float | `2.5`, `1e-3`, `0.5` | `float` |
| String | `"CA"`, `'CA'` | `string` (either quote, but they must match; no line breaks) |
| Range | `1:5`, `2.0:8.0`, `:10`, `4:` | `irange` / `frange` |
| Array | `{1, 2, 3}`, `{"CA", "N"}` | `int[3]`, `string[2]` |
| Constant | `PI`, `TAU`, `E` | `float` (`PI` and `TAU` carry the unit *rad*) |

### Ranges

A range is written `beg:end` and is **inclusive** at both ends. Either side may be left out to mean "unbounded":
`:10` is everything up to 10, `4:` is everything from 4 on, `:` is everything.

An optional stride goes in the middle, `beg:step:end`:

```mdscript
a = count(atom(3:10));        # atoms 3,4,...,10   -> 8
b = count(atom(1:2:10));      # atoms 1,3,5,7,9    -> 5
c = count(atom(3:));          # atom 3 to the last atom
```

Integer ranges convert implicitly to float ranges, so `within(2.0:6.0, ...)` and `atom(1:5)` both work.

Indices are never counted from the end. A negative bound has to be parenthesised (`(-3):5`; the bare `-3:5` is read
as a subtraction and rejected), and a range or index outside the valid range of its context — `atom(1:500)` in a
153-atom system, `atom(0:5)` — is a compile error that names the valid range.

### Arrays and subscripts

`{ }` builds an array. Elements must share a type after implicit conversion (`{1, 2, 3.0}` becomes a float
array, `{1, 2, 4, 7:9}` an integer-range array).

`[ ]` subscripts an array. Indices are 1-based, and a range (with optional stride) selects several elements:

```mdscript
first_res = residue(:)[1];              # first residue
some_res  = residue(:)[2:4];            # residues 2, 3, 4
odd_res   = residue(:)[1:2:9];          # residues 1,3,5,7,9
d = distance(1, 2) in residue(1:3);
d2 = d[2];                              # second element of d
```

### Operators

From tightest to loosest binding:

| Level | Operators | Notes |
|---|---|---|
| 1 | calls, `[ ]`, `( )` | |
| 2 | unary `-`, `not` | |
| 3 | `*`, `/`, `//` | `//` is integer division (integers only) |
| 4 | `+`, `-` | |
| 5 | `<`, `>`, `<=`, `>=` | comparison of `float`s, giving `bool` |
| 6 | `==` | equality of `float`s, giving `bool` |
| 7 | `and` | |
| 8 | `xor` | |
| 9 | `or` | |
| 10 | `out` | see [Contexts](#contexts-in-and-out) |
| 11 | `in` | |
| 12 | `=` | |

Arithmetic (`+ - * /`) works on `int`, `float`, and their arrays, element-wise; a scalar is broadcast over an
array. `volume` values can be added, subtracted, multiplied and divided with each other or with a scalar. A
`distribution` can only be scaled by a number (`2.0 * g`, `g / 2.0`) and combined with another by
[`pool`](#pool); see [Distributions](#distributions) for why. Comparisons work on `float`s and arrays of them (a scalar is compared with every element) and give a
`bool` or an array of `bool`. `and`, `or`, `xor`, `not` combine **atomsets** (atom selections) and, separately,
`bool`s and arrays of `bool` element by element. A `bool` value cannot be published as a property, so comparisons are
building blocks for other expressions rather than results in their own right.

```mdscript
d = distance(1, 2) in residue(1:3);
close = d < 1.01;                       # bool[3]
both = (d > 1.005) and (d < 1.02);      # element-wise logic on bool arrays
```

`in` binds looser than every other operator (only `=` is looser), so `a + b in c` means `(a + b) in c`. Use
parentheses whenever the grouping is not obvious:

```mdscript
notO   = (not element("O")) in resname("ALA");
mixed  = (residue(1) or residue(2)) and element("O");    # `and` binds tighter than `or`
```

### Named arguments

Procedures listed with a *Parameters* table can be called with named arguments. Positional arguments come first,
then named ones, in any order. Skipping an optional parameter is allowed.

```mdscript
a = distance(a=1, b=2) in residue(1:3);
sel = within(radius=3.0, around=atom(1));
n = count(sel=protein(), unit="residue");
```

---

## Types and units

| Type | Description |
|---|---|
| `int`, `float`, `bool` | Scalars. `int` converts implicitly to `float`. |
| `string` | Text, used to match names. |
| `irange`, `frange` | Integer / float ranges; `irange` converts to `frange`. |
| `atomset` | A set of atoms — the value of every selector. |
| `T[n]` | An array of `n` elements of type `T`. `atomset[n]` is *n* selections, typically *n* residues or chains. |
| `float[2]`, `float[3]`, `float[4]`, `float[4][4]` | Vectors and matrices (`vec2`, `vec3`, `vec4`, 4×4 matrix). |
| `position` | Not a real type: the argument of a geometric procedure. It accepts atom indices or ranges (`1:3`), atomsets (their centre of mass is used), or `float[3]` vectors. See [Coordinate arguments](#coordinate-arguments). |
| `distribution` | A 1024-bin histogram (for example `rdf`, `density_x`), with a weight per bin; see [Distributions](#distributions). `distribution[n]` is *n* of them, e.g. `rdf(...) in residue(1:3)`. |
| `volume` | `float[128][128][128]` grid (for example `sdf`). |
| `contact` | The group pairs in contact in a frame, from [`contacts`](#contacts); read with [`count`](#count) and [`degree`](#degree). |
| `plane` | The plane fitted to a group of atoms or to positions, from [`plane`](#plane): its centroid, oriented normal, radius and how far the atoms are from it; read with [`centroid`](#centroid), [`normal`](#normal), [`radius`](#radius), [`deviation`](#deviation) and [`angle`](#angle). `plane[n]` is *n* of them. |
| `bondset` | A set of bonds, from [`bond`](#bond): the persistent bonds of the system it picks, and bonds computed for a frame (those [`formed`](#formed) in it). Not a set of atoms, which would bring in every bond among its atoms; [`atom`](#atom) gives the atoms of the bonds, and a filter of atoms (a representation, a selection) which picks bonds picks their atoms. Within an expression the two do not mix: `protein() and bond()` is an error. `bondset[n]` is *n* sets, for example the bonds within each residue. |
| `bond_delta` | The bonds formed and broken in a frame, from [`bond_delta`](#bond_delta); read with [`formed`](#formed), [`broken`](#broken) and [`bond`](#bond). |
| `pucker` | How a ring of 5 to 7 atoms is puckered, from [`pucker`](#pucker): its Cremer-Pople coordinates and the name of its conformation; read with [`amplitude`](#amplitude), [`phase`](#phase), [`theta`](#theta) and [`pseudorotation`](#pseudorotation). |

Units propagate through arithmetic: `distance(1,2) - distance(1,3)` is in Å, `distance(1,2) * distance(1,3)` in Å²,
`distance(1,2) / distance(1,3)` is dimensionless, and scaling by a plain number (`2.0 * distance(1,2)`) keeps the
unit, and so does unary minus. Other operations produce a plain number without a unit tag: adding a bare number to a
quantity (`distance(1,2) + 1.0`) and the math functions (`abs`, `sqrt`, `min`, ...). Reported units: `distance`
in Å, `angle` and `dihedral` in rad, `density_*` in kg/m³, `sdf` in count/Å³; `PI` and `TAU` carry the unit rad.
`count`, `rdf` and `porosity` are dimensionless. Treat a missing unit tag as "unknown", not as a guarantee that the
value is dimensionless.

Periodic values propagate the same way. A value is periodic when it is the same as itself plus a whole period, as an
angle about an axis is: `dihedral` (in `[-π, π)`), `phase` and `pseudorotation` (in `[0, 2π)`) and `atan2`. Moved
along or scaled by a constant (`dihedral(1,2,3,4) * 180 / PI` has a period of 360), negated, added to or subtracted
from another of the same period (`dihedral(1,2,3,4) - dihedral(2,3,4,5)`), gathered into an array, indexed, or the
`min`, `max` or `sum` of an array of them, it stays periodic. What no longer lies on a circle (`abs`, `cos`, a
product of two angles) does not. A periodic property is plotted cut where its values are least likely to be: a
dihedral at 0, where torsions are eclipsed, a pseudorotation at 270 degrees (O4'-exo). Its mean, spread and extent
over a population are circular: the mean of 179 and -179 degrees is 180, not 0.

### Coordinate arguments

Every procedure documented with `position` arguments (`distance`, `angle`, `dihedral`, `com`, `plane`, `coord`,
`rdf`, `within(radius, around)`, ...) accepts any of:

- **atom indices** such as `1`, `1:3` or `{1, 5, 9}`. Inside a context these are relative to that context;
- **an atomset** such as `residue(1)` or `element("C")`;
- **a `float[3]` vector**, e.g. `vec3(0, 0, 0)`, or the result of `com(...)`;
- **an array of atomsets**, e.g. `residue(1:3)` or `ring()`.

How a selection is turned into positions depends on the procedure, but the rule of thumb is:

| Argument | `distance`, `angle`, `dihedral`, `com` | `distance_pair`, `coord*` | `rdf`, `within`, `plane`, `shape_weights` |
|---|---|---|---|
| single atomset / index range | one point: the mass-weighted centre of mass | every atom, individually | every atom, individually (see each entry) |
| array of atomsets (`residue(1:3)`) | see below | one centre of mass **per element** | see each entry |

`distance`, `com` and `within` *merge* an array of atomsets into one selection first, so `com(residue(1:3))` is the
centre of mass of all three residues together and `distance(residue(1:3), vec3(0,0,0))` is a single value. To get
one value per residue, use a [context](#contexts-in-and-out): `distance(1, vec3(0,0,0)) in residue(1:3)`, or use a
procedure that keeps the elements apart such as `distance_pair(residue(1:3), vec3(0,0,0))` (three values) or
`coord_x(residue(1:3))` (three values: the x coordinate of each residue's centre of mass).

Procedures that describe a structure rather than reduce it to a point, [`plane`](#plane), [`rmsd`](#rmsd) and
[`pucker`](#pucker), keep the atomsets of an array apart and give one result per atomset. Take their
[`union`](#union) for one result over all of them: `plane(union(residue(1:3)))`.

`distance_min` and `distance_max` follow a rule of their own: they measure atoms, not centres of mass, and give one
value per element of their first argument (`distance_min(residue(1:3), residue(5))` is three values), while their
second argument is taken as a whole.

---

## Contexts: `in` and `out`

`expr in selection` evaluates `expr` **once for each structure in `selection`**, and inside `expr` all atom indices are
**relative to that structure**. The result gains a leading dimension the size of the context.

```mdscript
# Atom 1 and 2 of every residue: one distance per residue -> float[num_residues]
d = distance(1, 2) in residue(:);

# Restrict to a subset of residues
d2 = distance(1, 2) in residue(1:3);

# Local atom indices are relative to the context: atoms 1..4 of every residue
sel = atom(1:4) in residue(:);
```

If the context is a single structure the result is a single value. An index that does not exist in the context
(e.g. atom 30 in a 12-atom residue) is a compile error.

`out expr` **lifts** the enclosed expression out of the surrounding context so that it is evaluated globally:

```mdscript
# Distance from atom 1 of each residue to the centre of mass of *global* atom 1:
d = distance(1, out com(atom(1))) in residue(1:3);

# The second operand is atom 4 of the whole system, not of a residue:
e = distance(com(atom(1:3)), out atom(4));
```

The keyword `of` is reserved for future use.

---

## What becomes a property

Assignments compile to *values*. Only some of them are published as **properties** and can be plotted or exported:

- a result that varies with the frame (**temporal**): `distance`, `angle`, `count(within(...))`, `com`, ...
- a `distribution` (`rdf`, `density_x/y/z`), accumulated over all frames, or an array of them,
- a `volume` (`sdf`), accumulated over all frames.

Constants, plain selections (an `atomset`), strings and ranges are *not* properties; they exist to be used by other
statements. The engine publishes each property under the attribute `script/<name>`. A script with no property at
all (`a = count(protein());`, which is constant) cannot be evaluated and reports "no properties".

Some values have no property form yet: `float[2]`/`float[3]` arrays (`coord`, `coord_xy`, `shape_weights` without
destructuring) cannot be plotted directly — use them as inputs to other procedures or
[destructure](#statements-and-comments) them. An array of distributions (`density(...)`, `rdf(...) in residue(1:3)`)
is published, as `[n][1024]`, but not yet plotted.

### Distributions

A bin of a distribution is a value and a weight, and what it shows is the value divided by the weight: for `rdf` the
pair count of a shell over what a uniform distribution would put there. Both are sums, over the frames and over bins
merged for display, so the estimate is the ratio of the sums. That is also why a distribution only scales and
pools: `(a + b)` or `(a - b)` per frame, then summed over the frames, is not `a + b` of the estimates. Scaling
(the values, not the weights) and pooling (both, see [`pool`](#pool)) come out the same either way.

---

## Procedure index

| Category | Procedures |
|---|---|
| [Selectors: atoms](#selectors-atom-level) | [`all`](#all), [`atom`](#atom), [`element`](#element), [`name`](#name) (`label`, `type`), [`backbone`](#backbone), [`side`](#side) (`sidechain`), [`ion`](#ion), [`nucleoside`](#nucleoside), [`nucleobase`](#nucleobase), [`head`](#head), [`tail`](#tail), [`ring`](#ring) |
| [Selectors: by kind](#selectors-by-kind) | [`polymer`](#polymer), [`protein`](#protein), [`nucleic`](#nucleic), [`water`](#water), [`lipid`](#lipid), [`glycan`](#glycan) |
| [Selectors: residues](#selectors-residue-level) | [`nucleotide`](#nucleotide), [`qm`](#qm), [`environment`](#environment), [`resname`](#resname) (`residue`, `component`), [`resid`](#resid), [`residue`](#residue), [`component`](#component) |
| [Selectors: chains, molecules and parts](#selectors-chains-molecules-and-parts) | [`chain`](#chain) (`chain_id`, `auth_id`), [`molecule`](#molecule), [`part`](#part) |
| [Selectors: structure matching](#selectors-structure-matching) | [`smiles`](#smiles), [`match`](#match) |
| [Selectors: spatial](#selectors-spatial) | [`within`](#within), [`within_x`](#within_x), [`within_y`](#within_y), [`within_z`](#within_z), [`within_xyz`](#within_xyz) |
| [Properties](#properties) | [`distance`](#distance), [`distance_min`](#distance_min), [`distance_max`](#distance_max), [`distance_pair`](#distance_pair), [`angle`](#angle), [`dihedral`](#dihedral), [`rmsd`](#rmsd), [`rdf`](#rdf), [`density_x`](#density_x) / [`density_y`](#density_y) / [`density_z`](#density_z), [`density`](#density), [`pool`](#pool), [`sdf`](#sdf), [`count`](#count), [`contact_count`](#contact_count), [`contacts`](#contacts), [`degree`](#degree), [`porosity`](#porosity) |
| [Geometry](#geometry) | [`com`](#com), [`plane`](#plane), [`centroid`](#centroid), [`normal`](#normal), [`radius`](#radius), [`deviation`](#deviation), [`pucker`](#pucker), [`amplitude`](#amplitude), [`phase`](#phase), [`theta`](#theta), [`pseudorotation`](#pseudorotation), [`shape_weights`](#shape_weights), [`coord`](#coord), [`coord_x`](#coord_x) / [`coord_y`](#coord_y) / [`coord_z`](#coord_z), [`coord_xy`](#coord_xy) / [`coord_xz`](#coord_xz) / [`coord_yz`](#coord_yz) |
| [Math](#math) | [`sqrt`](#sqrt), [`cbrt`](#cbrt), [`abs`](#abs), [`floor`](#floor), [`ceil`](#ceil), [`sin`](#sin), [`cos`](#cos), [`asin`](#asin), [`acos`](#acos), [`atan`](#atan), [`atan2`](#atan2), [`log`](#log), [`log2`](#log2), [`log10`](#log10), [`exp`](#exp), [`exp2`](#exp2), [`pow`](#pow), [`min`](#min), [`max`](#max) |
| [Linear algebra](#linear-algebra-and-constructors) | [`vec2`](#vec2), [`vec3`](#vec3), [`vec4`](#vec4), [`dot`](#dot), [`cross`](#cross), [`length`](#length), [`normalize`](#normalize), [`mul`](#mul) |
| [Bond sets](#bond-sets) | [`bond`](#bond), [`atom`](#atom), [`bond_delta`](#bond_delta), [`formed`](#formed), [`broken`](#broken), and [`count`](#count), [`length`](#length) and `and` `or` `xor` `not` on bond sets |
| [Combining](#combining-and-reshaping) | [`union`](#union), [`flatten`](#flatten), [`split`](#split), [`chunks`](#chunks), [`transpose`](#transpose), and the operators `and` `or` `xor` `not` on selections |
| [External data](#data-from-outside-the-structure) | [`attr`](#attr) |

---

## Selectors

A **selector** returns an `atomset` (a set of atoms) or a `atomset[]` (several sets — one per residue, chain, ring,
...). Selectors take no positions and do not change during a run, except the [spatial selectors](#selectors-spatial),
which are re-evaluated for every frame.

A name, residue name or identifier that matches nothing in the loaded system is a **compile error** ("The string
'X' did not match any ..."), which catches typos early but means a script written for a protein will not compile
against a system without one.

Selections combine with `and`, `or`, `xor`, `not`:

```mdscript
ca_of_ala  = name("CA") and resname("ALA");
not_water  = not water();
polar      = element({"N", "O"}) or name("S*");
```

### Selectors: atom level

These return a single `atomset`.

### all

<!-- proc name=all category=selector.atom -->

```text
all() -> atomset
```

Every atom in the system (or, inside a context, every atom of that context).

```mdscript
n_atoms = count(all());
per_residue = count(all()) in residue(:);      # atoms per residue
```

### atom

<!-- proc name=atom category=selector.atom -->

```text
atom(indices: irange[]) -> atomset
atom(indices: int[])    -> atomset
atom(b: bondset[])      -> atomset[]      # the atoms of bonds
```

Select atoms by 1-based index. Indices are relative to the current context, so inside `in residue(:)` `atom(1)` is
the first atom of each residue. An index outside the valid range is a compile error. Note that `atom` takes one
argument: for several atoms write an array, `atom({1, 5, 9})`, or a range, `atom(1:9)`.

Given [bond sets](#bond-sets), the atoms of their bonds, one set of atoms per set of bonds, as
[`residue`](#residue) grows a selection to its residues. This is how bonds become atoms within an expression: a set of
bonds and a set of atoms do not combine as they are. The complement of a bond set is the other bonds, so
`not atom(bond(kind="hydrogen"))` is every atom outside of a hydrogen bond, while `atom(not bond(kind="hydrogen"))` is
the atoms of every other bond, which a hydrogen bond's atoms have too.

```mdscript
first_ten = atom(1:10);
every_other = atom(1:2:20);                    # 1,3,5,...,19
n_terminal = atom(1:4) in residue(1);          # first four atoms of residue 1
metal_ligands = atom(bond(kind="coordinate")); # the metals and their ligand atoms
```

### element

<!-- proc name=element category=selector.atom -->

```text
element(symbol: string[]) -> atomset
element(z: irange[])      -> atomset
```

Select by chemical element, either by symbol (`"C"`) or atomic number (`6`). Several values may be given as an
array.

```mdscript
carbons = element("C");
heavy_polar = element({"N", "O"});
also_carbons = element(6);
```

### name

<!-- proc name=name aliases=label,type category=selector.atom -->
<a id="label"></a><a id="type"></a>

```text
name(pattern: string[])  -> atomset
label(pattern: string[]) -> atomset
type(pattern: string[])  -> atomset
```

Select atoms by atom name (for example `"CA"`, `"OW"`). `label` and `type` are aliases. Patterns may contain the
wildcard `*`, which stands for any run of characters, anywhere and as often as needed: `"C*"` begins with C, `"*1"`
ends with 1, `"*H*"` contains H, `"C*A"` begins with C and ends with A. Names are matched case and all.

```mdscript
alpha_carbons = name("CA");
any_carbon = name("C*");
backbone_atoms = name({"N", "CA", "C", "O"});
```

### backbone

<!-- proc name=backbone category=selector.atom -->

```text
backbone() -> atomset
```

The atoms flagged as backbone in protein and nucleic-acid residues.

```mdscript
bb = backbone();
n_bb = count(backbone());
```

### side

<!-- proc name=side aliases=sidechain category=selector.atom -->
<a id="sidechain"></a>

```text
side()      -> atomset
sidechain() -> atomset
```

`sidechain()` selects the side-chain atoms of protein residues. `side()` selects the same plus the nucleoside
(sugar) atoms of nucleic acids, so it is the complement of `backbone()` for both kinds of polymer. Despite the
similar names, `sidechain` is **not** an alias of `side` in nucleic acids.

```mdscript
sc = sidechain();
sc_of_first_residues = side() and residue(1:5);
```

### ion

<!-- proc name=ion category=selector.atom -->

```text
ion() -> atomset
```

Atoms of residues that were recognised as ions (for example Na⁺, Cl⁻).

```mdscript
ions = ion();
solvent_and_ions = water() or ion();
```

### nucleoside

<!-- proc name=nucleoside category=selector.atom -->

```text
nucleoside() -> atomset
```

Atoms flagged as part of the nucleoside (sugar) portion of nucleic-acid residues.

```mdscript
sugars = nucleoside();
```

### nucleobase

<!-- proc name=nucleobase category=selector.atom -->

```text
nucleobase() -> atomset
```

Atoms flagged as part of the nucleobase of nucleic-acid residues.

```mdscript
bases = nucleobase();
```

### head

<!-- proc name=head category=selector.atom -->

```text
head() -> atomset
```

The atoms of the heads of the lipids: all of a lipid but its tails (see [`tail`](#tail)), hydrogens included. For a
phospholipid the head group, the phosphate, the glycerol and the ester or ether linkages; for a sterol its hydroxyl;
for a bolalipid both of its heads; for a Martini lipid the beads of its head and linker (NC3, PO4, GL1, GL2, ROH, ...).
`head() or tail()` is every lipid whose heads and tails are known: a lipid by name in a system without bonds has none.

```mdscript
heads = head();
p_z = coord_z(head() and element("P"));
```

### tail

<!-- proc name=tail category=selector.atom -->

```text
tail() -> atomset[]
```

One selection per tail of the lipids, with its branches and hydrogens: each acyl, alkyl or sphingoid chain, and the ring
skeleton of a sterol with its side chain. A tail is a part of the lipid made of carbon bonded to carbon and hydrogen only
which hangs off the rest of it (found from the bonds, whatever the atoms are named), or for a Martini lipid the beads of
a tail (C1A, D2A, ..., R1, ...). The chains of a bolalipid span its two heads. The tails come in the order of the lipids
and, within a lipid, of their atoms. Within a context only the tails fully within it count, as for [`ring`](#ring).
Visualized as the line of each tail from its head to its tip.

```mdscript
tails = tail();
tail_z = coord_z(tail());                   # the z of the centre of each tail
n_tail_atoms = count(tail() and element("C"));
```

### ring

<!-- proc name=ring category=selector.atom -->

```text
ring(pattern: string[], kind: string[], size: irange) -> atomset[]
ring(index: irange[], kind: string[], size: irange)   -> atomset[]      # by position
```

**Parameters:** `pattern` or `index` (optional), `kind` and `size` (optional, given by name).

One selection per ring of the molecular graph (rings of 3 to 7 atoms, see below), useful to measure something per
ring. Inside a context only rings that lie completely within the context are returned. Each parameter narrows the
rings down, and a ring has to pass all of those given; a list given to one of them is any of its entries.

- `pattern`: [SMILES](#smiles), read as there, along the ring. A chain is consecutive atoms of the ring in either
  direction and from any atom: `"CCNCC"` is in pyridine as in piperidine. A chain closed back onto its first atom is
  the whole ring: `"c1ccncc1"`. Parts separated by `.` are anywhere in the ring, on atoms of their own: `"N.O"` is a
  ring with an N and an O (oxazole, morpholine). The atoms of a ring have two neighbours in it, so a pattern which
  branches is a compile error: substituents are for [`smiles`](#smiles). `~` is a bond of any order: `"C~N"` takes
  a C=N in a ring which is not aromatic, which `"CN"` (a single bond) does not.
- `kind`: `"aromatic"` (every bond of the ring aromatic), `"aliphatic"` (not aromatic), `"hetero"` (an atom other than
  C), `"heteroaromatic"` (both), `"sugar"` (the ring of a furanose or pyranose).
- `size`: the number of atoms, an int or a range.
- `index`: rings by their **1-based position** among the rings the other parameters let through, in the order of
  the system, as [`residue(1:5)`](#residue) counts residues: `ring(3)` is the third ring of the system, which is the
  index the circle of an aromatic ring shows when hovered. Inside a context the positions count the rings within
  it: `ring(1) in residue(:)` is the first ring of every residue that has one. A position past the last ring is a
  compile error without a context, and selects nothing within one; `ring(:)` is every ring.

The rings are the smallest ones, without those made of smaller rings: naphthalene has two rings and not the ten
membered one around both, norbornane two five membered rings and not the six membered one around them. Larger cycles
are loops and not rings. The atoms of a sugar ring are in order from the ring O and the anomeric C.

```mdscript
n_ring_atoms = count(ring());                           # the rings are merged into one selection
d = distance_pair(ring(), vec3(0, 0, 0));               # centre of each ring to the origin, one value per ring
aromatic = ring(kind="aromatic");
pyridines = ring("c1ccncc1");
azoles = ring({"c[nH]c", "cnc"}, size=5);
sugars = ring(kind="sugar") in residue(1:20);
third = ring(3);                                        # the third ring of the system
first_aromatic = ring(1, kind="aromatic") in residue(:); # the first aromatic ring of each residue that has one
```

### Selectors: by kind

These return a `atomset[]` with one selection per [part](#part) of the kind: a strand of a polymer (with its caps and
the modified residues bonded into it), a water. [`lipid`](#lipid) is one per lipid molecule, [`glycan`](#glycan) one per
glycan tree. Within a context, the parts with an atom in it, whole. By position (`protein(2)`) the parts of the kind are
counted as [`part`](#part) counts them, within a context those with an atom in it, each cut to it. Use them as the
right-hand side of `in`, or pass them where an atomset is accepted (an array of atomsets is then merged into one
selection; [`rmsd`](#rmsd), [`plane`](#plane) and [`pucker`](#pucker) are the exceptions and give one value per
selection). One selection per residue is [`residue`](#residue) of them: `residue(protein())`.

### polymer

<!-- proc name=polymer category=selector.kind -->

```text
polymer()                 -> atomset[]
polymer(index: irange[])  -> atomset[]
```

Every strand of a polymer: a peptide, a nucleic acid or any other polymer, one selection per strand. A strand goes on
over a gap (residues missing from a crystal structure) and ends where the chemistry says the chain does. This is how a
strand is selected where the data gives no chain ids (a GRO file, a topology, a PDB file written by a simulation): its
position among the strands, `polymer(2)`. Where it does, `polymer() in chain("A")` is the strand of chain A without the
ligands and waters filed under it.

```mdscript
strands = polymer();                                # one selection per strand
second = polymer(2);                                # the second strand of the system
strand_a = polymer() in chain("A");                 # the strand of chain A alone
```

### protein

<!-- proc name=protein category=selector.kind -->

```text
protein()                 -> atomset[]
protein(index: irange[])  -> atomset[]
```

Every peptide strand, one selection per strand: its amino acids (those whose backbone was found and verified by its
bonds, and those only named like an amino acid: missing atoms, a coarse grained model), with the caps (ACE, NME) and
the modified residues bonded into it. A free amino acid (a single residue bound as a ligand) is no strand.

```mdscript
n_protein_atoms = count(protein());                 # the strands are merged into one selection
first = protein(1);                                 # the first peptide strand
d = distance(1, 2) in residue(protein());           # one value per protein residue
```

### nucleic

<!-- proc name=nucleic category=selector.kind -->

```text
nucleic()                 -> atomset[]
nucleic(index: irange[])  -> atomset[]
```

Every strand of a nucleic acid (DNA, RNA or a hybrid), one selection per strand. A free nucleotide ligand (GTP, ATP)
is no strand. For the nucleotides one by one, see [`nucleotide`](#nucleotide).

```mdscript
dna = nucleic();
strands = count(nucleic(), "part");                 # the strands of a duplex
```

### water

<!-- proc name=water category=selector.kind -->

```text
water() -> atomset[]
```

Every water molecule, one selection each.

```mdscript
solvent = water();
n_waters = count(water(), "residue");
```

### lipid

<!-- proc name=lipid category=selector.kind -->

```text
lipid() -> atomset[]
```

Every lipid, sterols included: the lipids of the Martini lipidomes, coarse grained or atomistic by the same names
(POPC, DPPC, CHOL, ...), the cholesterol of CHARMM (CHL1), and the molecules which are lipids by their structure
whatever they are named: a polar head with a hydrocarbon chain of at least 8 carbons (a fatty acid, a phospholipid
whose residue is named DPP in a PDB file), or a sterol. One selection per lipid MOLECULE, as for
[`molecule`](#molecule): a lipid split into several residues (head and tails, as the AMBER lipids are) is one selection
of all of them. Within a context, the lipids with an atom in it, whole. See [`head`](#head) and [`tail`](#tail) for its
parts.

```mdscript
membrane = lipid();
n_lipids = count(lipid(), "molecule");                # lipid molecules
```

### glycan

<!-- proc name=glycan category=selector.kind -->

```text
glycan() -> atomset[]
```

Every glycan, one selection per glycan [part](#part): a tree of saccharides attached to a protein (N- or O-linked), an
oligosaccharide on its own, a single sugar. Its residues are those with the ring of a monosaccharide, a furanose or a
pyranose, whatever they are named (NAG, BGLCNA, MAN, ...). A glycan of a glycoprotein is of the protein's molecule,
so `glycan()` is how it is selected alone. Within a context, the glycans with an atom in it, whole.

```mdscript
glycans = glycan();                                   # one selection per glycan tree
glycans_a = glycan() in chain("A");                   # the glycans of chain A
n = count(glycan(), "part");
```

### Selectors: residue level

These return a `atomset[]` with one selection per matching residue.

### nucleotide

<!-- proc name=nucleotide category=selector.residue -->

```text
nucleotide()  -> atomset[]
```

Every nucleotide residue: those whose backbone was found and verified by its bonds, and those which are only named
like a nucleotide. A free nucleotide ligand (GTP, ATP) is not one unless it is named like a residue of a nucleic acid.

```mdscript
bases = nucleotide();
n = count(nucleotide());                              # the nucleotides, merged into one selection
```

### qm

<!-- proc name=qm category=selector.residue -->

```text
qm() -> atomset[]
```

The QM region: the atoms of the quantum chemistry calculation the system was loaded from (VeloxChem, Molden,
TREXIO), one selection per residue with atoms in it. A system from a QM calculation alone is QM throughout; in a
polarizable embedding or QM/MM system it is the region the embedding surrounds. A residue need not lie wholly on one
side: a QM/MM boundary can cut through it, and its selection then holds the QM atoms only.

Only a loader says which atoms are QM, nothing infers it: in a system without a QM region (a structure or
trajectory file), `qm` is a **compile error** ("The system has no QM region") rather than an empty selection.

```mdscript
d = distance_min(qm(), water());     # the closest water atom to each residue of the QM region
near_qm = within(5, qm());
```

### environment

<!-- proc name=environment category=selector.residue -->

```text
environment() -> atomset[]
```

The environment of the QM region: every atom which is not in it - the sites of a polarizable embedding, the MM
atoms of a QM/MM system - one selection per residue (an embedding fragment is one). A **compile error** where the
system has no QM region, or is QM throughout ("The system is QM throughout: it has no environment").

```mdscript
first_shell = within(3.5, qm()) and environment();
n_fragments = count(environment(), "residue");
```

### resname

<!-- proc name=resname aliases=residue,component category=selector.residue -->

```text
resname(pattern: string[])   -> atomset[]
residue(pattern: string[])   -> atomset[]
component(pattern: string[]) -> atomset[]
```

Select residues by name. `residue` and `component` are aliases (see [`residue`](#residue) for the other overloads).
The wildcard `*` works as in [`name`](#name).

```mdscript
alanines = resname("ALA");
hydrophobic = resname({"ALA", "VAL", "LEU", "ILE"});
any_l = resname("L*");
```

### resid

<!-- proc name=resid category=selector.residue -->

```text
resid(ids: irange[]) -> atomset[]
```

Select residues by their **sequence number** as stored in the file (the number printed in a PDB or GRO file).
Compare with [`residue(ids)`](#residue), which uses the 1-based *position* in the system.

```mdscript
first_two = resid(1:2);
loop = resid(10:20);
```

### residue

<!-- proc name=residue category=selector.residue -->

```text
residue(pattern: string[])   -> atomset[]      # same as resname
residue(index: irange[])     -> atomset[]      # by position
residue(sel: atomset[])     -> atomset[]      # grow to whole residues
```

- `residue("ALA")` — select residues by name.
- `residue(1:5)` — select residues by their **1-based position** in the system, `residue(:)` selects all. Inside a
  context the positions count the residues with an atom in the context, each cut to it: `residue(2) in chain(1)` is
  the second residue of the chain. A position past the last is a compile error without a context, and selects nothing
  within one (a chain need not have a tenth residue). [`chain`](#chain), [`molecule`](#molecule) and
  [`ring`](#ring) count their positions the same way.
- `residue(sel)` — every residue that the selection touches, expanded to its whole residue.

```mdscript
first = residue(1);
first_five = residue(1:5);
all_res = residue(:);
second_of_chain = residue(2) in chain(1);
ca_residues = residue(name("CA"));           # whole residues that contain a CA atom
near = residue(within(4.0, residue(1)));     # whole residues within 4 Å of residue 1
```

### component

<!-- proc name=component category=selector.residue -->

```text
component(pattern: string[]) -> atomset[]
component(index: irange[])   -> atomset[]
```

Alias of [`residue`](#residue) for the `string` and `irange` overloads. Coarse-grained and non-protein systems
often call their building blocks *components*.

```mdscript
first_component = component(1);
named = component("PFT");
```

### Selectors: chains, molecules and parts

A **chain** is what the data files under one chain identifier: the chain id of a PDB file, the author's chain id
(`auth_asym_id`) of an mmCIF file. Whatever carries the id is of the chain, a protein and the ligands, ions and waters
filed with it alike, and `chain("A")` is all of them, as **one** selection however the atoms are spread through the file.
A chain id is a label, as an atom name is.

Where the data gives no chain ids (a GRO file, a topology, a PDB file written by a simulation), none are made up: the
system has no chains, and `chain(...)` is a compile error, in each of its forms. A strand of a polymer is then found as what it is, a part of
a polymer: [`polymer`](#polymer), `polymer(1)` by its position.

A **molecule** is atoms joined by bonds, and what holds together where the bonds do not say so: the atoms of a residue,
and the residues of a polymer chain over a gap (residues missing from a crystal structure). A residue and a polymer
chain are always within one molecule; a molecule may hold several chains (joined by a disulfide), and a ligand bonded
to a protein is of its molecule. A metal coordinated by residues (the Mg of a GTP, the Zn of a zinc finger) is a molecule
of its own, as are the waters around it, unless a topology bonds it into a molecule. An ion or a water is a molecule of
its own. Where the system comes with a topology (a GROMACS `.tpr`, or a `.top` beside a `.gro`), its molecules are the
molecules of the topology, as GROMACS has them, each of its molecule type: the two strands of a DNA duplex which
martinize-dna writes as one molecule type are one molecule of two strands.

A **part** is one polymer strand, one glycan tree, one ligand, one ion or one water: the parts the molecules are built
of. A glycoprotein is one molecule of many parts, its strand and each of its glycans; a chain holds its strand and the
ligands and waters filed under its id, each a part. A strand goes on over a gap and ends where the chemistry says the
chain does; residues of no polymer are one part while each is bonded to one of it (the branches of a glycan). A bond
from a side chain (a glycosylation, a disulfide) and a coordination bond join no parts.

### chain

<!-- proc name=chain aliases=chain_id,auth_id category=selector.chain -->
<a id="chain_id"></a><a id="auth_id"></a>

```text
chain(id: string[])     -> atomset
chain(index: irange[])  -> atomset[]
chain(sel: atomset[])   -> atomset[]     # grow to whole chains
chain_id(id: string[])  -> atomset
auth_id(id: string[])   -> atomset
```

- `chain("A")` — every atom of the chain id, as one selection. Several ids, or a pattern (`"*"` as in
  [`name`](#name)), select all of their atoms together. An id no chain has is a compile error which lists the ids
  there are, and so is any id in a system whose data gives none (see [`polymer`](#polymer)). Inside a context, the atoms
  of the context which carry the id. `chain_id` and `auth_id` are other names of it.
- `chain(1:2)` — chains by their 1-based position, one selection per chain, counted as [`residue`](#residue) counts
  its positions. The chains are in the order of their first atoms.
- `chain(sel)` — every chain that the selection touches, whole, one selection per chain. Atoms in no chain add none.

```mdscript
chain_a = chain("A");                          # the protein of chain A, and its ligands and waters
protein_a = protein() in chain("A");           # the protein alone
first_two = chain(1:2);                        # two selections
its_chain = chain(atom(1));                    # the chain that contains atom 1
d = distance(com(chain("A")), com(chain("B")));
```

### molecule

<!-- proc name=molecule category=selector.chain -->

```text
molecule(index: irange[])  -> atomset[]
molecule(sel: atomset[])   -> atomset[]     # grow to whole molecules
```

- `molecule(1:3)` — molecules by their 1-based position, one selection per molecule, counted as
  [`residue`](#residue) counts its positions. The molecules are in the order of their first atoms, and their atoms need
  not be consecutive.
- `molecule(sel)` — every molecule that the selection touches, whole, one selection per molecule.

```mdscript
its_molecule = molecule(atom(1));             # the molecule that contains atom 1
near = molecule(within(4.0, resname("LIG"))); # the molecules within 4 Å of the ligand, whole
n = count(water(), "molecule");
```

### part

<!-- proc name=part category=selector.chain -->

```text
part(index: irange[])  -> atomset[]
part(sel: atomset[])   -> atomset[]     # grow to whole parts
```

- `part(1:3)` — parts by their 1-based position, one selection per part, counted as [`residue`](#residue) counts its
  positions. The parts are in the order of their residues.
- `part(sel)` — every part that the selection touches, whole, one selection per part: the glycan an atom is of,
  without the protein it is attached to.

```mdscript
its_glycan = part(atom(19818));                # the glycan tree that atom 19818 is of
strands = part(nucleic());                     # each strand of a duplex
n = count(all(), "part");
```

### Selectors: structure matching

These find where a structure occurs in the bonds of the system, written as a SMILES pattern or given as a selection
of the system to use as a reference. They return a `atomset[]` with one selection per occurrence. The occurrences
depend on the bonds only, so they are found once, when the script is compiled. Inside a context only the atoms of
the context are searched. The full semantics are in `src/md_match.h`.

Both take two named arguments:

- `level`: where an occurrence lies. `"molecule"` (the default): within one [molecule](#molecule). Patterns follow
  covalent bonds only: a metal coordinated by another residue is a molecule of its own here even where a topology
  bonds it into one. `"residue"` (or `"component"`): within one residue, whose bonds to its
  neighbours are not part of the graph, so a residue pattern ends at the peptide bond. `"part"`: within one
  [part](#part) (a strand of a polymer, a glycan tree, a ligand), the bonds leaving it (a glycosylation, a disulfide)
  not followed. `"chain"`: within one [chain](#chain), all that the data gives its id, the bonds leaving it (a disulfide
  to another chain) not followed; atoms the data puts in no chain are not searched, and a system whose data gives no
  chain ids is a compile error.
- `mode`: which occurrences there are. `"unique"` (the default): one per set of atoms; a symmetric pattern fits the
  same atoms in several ways (a benzene ring in 12) and counts once. `"all"`: every way it fits. `"one_per_unit"`: the
  first in each unit of `level`. `"disjoint"`: occurrences which share no atoms, the first found winning.
  `"whole"`: occurrences which are all of their unit, its every atom (hydrogens aside) and no more bonds between them
  than the pattern has. `smiles("NCC=O", level="residue", mode="whole")` is the residues which are glycine; without
  `"whole"` every amino acid has those atoms.

What a system cannot tell is not compared. Hydrogen counts, bond orders, aromaticity and formal charges come from
the file or from the chemistry perceived when the system is loaded; a hydrogen count counts implicit hydrogens too,
so a methyl is a methyl in a structure without hydrogen atoms as well. Such a structure does not show where its
protons are, though, and a proton more or less is not held against an atom there: `C(=O)[OH]` and `C(=O)[O-]` both
find the carboxyl groups of a crystal structure. A charge in a resonance group (a carboxylate, a guanidinium)
belongs to the group, whichever of its atoms the charge was put on: `[O-]` finds both oxygens of a carboxylate.

### smiles

<!-- proc name=smiles category=selector.structure -->

```text
smiles(pattern: string, level: string, mode: string) -> atomset[]
```

**Parameters:** `pattern` (required), `level` (default `"molecule"`) and `mode` (default `"unique"`), given by name.

One selection per occurrence of a SMILES pattern, each with the hydrogens of its atoms. The pattern is read as a
query, the way SMARTS reads the same string:

- Atoms of the organic subset (`C`, `N`, `O`, `c`, `n`, ...) leave their hydrogens and charge open; atoms in brackets
  fix both. `CC(=O)O` finds acids, carboxylates and esters, `C(=O)[OH]` only the acids, `[NH3+]` an ammonium.
- Bonds have the order written, single where none is: `C=O` is a double bond and `CO` a single one. Aromatic bonds
  and resonance (a carboxylate, a guanidinium) match either, so a Kekule ring such as `C1=CC=CC=C1` finds an
  aromatic ring as well.
- Lowercase atoms have to be aromatic. In a pattern with lowercase atoms, the uppercase ones outside of rings have
  to be aliphatic: `Cc1cncn1` is the side chain of a histidine, and not a part of a purine.
- `~` is a bond of any order (from SMARTS): `C~N` finds C-N, C=N, C#N and aromatic C-N bonds alike.
- Isotopes, chirality and atom classes are read but not compared. `*` is any atom but a hydrogen. Hydrogens written
  as atoms count towards a bracket atom: `N[C@@]([H])(C)C=O` is `N[C@@H](C)C=O`.

A pattern which is not valid SMILES is a compile error, which points at the character.

```mdscript
methyls = smiles("[CH3]");
acids   = smiles("C(=O)[OH]");
gly     = smiles("NCC=O", level="residue", mode="whole");     # the residues which are glycine
rings   = smiles("c1ccccc1");                                  # phenyl rings, of PHE, TYR and any ligand
d = distance_pair(smiles("c1ccccc1"), vec3(0, 0, 0));          # centre of each ring to the origin
```

### match

<!-- proc name=match category=selector.structure -->

```text
match(ref: atomset, by: string, level: string, mode: string) -> atomset[]
```

**Parameters:** `ref` (required), `by` (default `"element"`), `level` and `mode` as for [`smiles`](#smiles), given by
name.

One selection per occurrence of `ref`: its atoms, with the bonds between them. Each selection holds the atoms which
take the place of those of `ref`, and no others. `by="name"` compares the atom names as well as the elements, which
only means something between residues which name their atoms the same way.

The bonds keep their order where the system knows it. Hydrogens in `ref` pin the protonation: an atom whose
hydrogens are all in `ref` has to have as many, and the same charge. Without any hydrogens `ref` is a skeleton, and
its atoms may have any. `ref` has to be the same for every frame.

```mdscript
like_5 = match(residue(5), level="residue", mode="one_per_unit");   # residues built like residue 5
cbs    = match(name("CB") in residue(5), by="name");
skel   = match(residue(5) and not element("H"));
```

### Selectors: spatial

Spatial selectors depend on atom positions, so they are re-evaluated for **every frame**; their size is not known
in advance. Consequently their results can be counted or measured, but not used as a context whose length must be
known up front.

### within

<!-- proc name=within category=selector.spatial -->

```text
within(radius: float,  around: position[]) -> atomset
within(radius: frange, around: position[]) -> atomset
within(radius: float)                       -> atomset      # centre from the context
within(radius: frange)                      -> atomset
```

**Parameters:** `radius` (required), `around` (optional).

All atoms within `radius` Å of the atoms (or points) in `around`. A range gives a shell: `within(2.0:4.0, x)`.
When `around` is a selection, its own atoms are **not** part of the result; a point (`vec3(...)`, `com(...)`) has no
atoms to exclude. Without `around`, the centre is the context: `within(3.0) in residue(:)` finds the neighbours of
each residue. Periodic boundary conditions are respected.

```mdscript
# Water molecules around the protein
hydration = residue(within(3.5, protein()) and element("O"));

# Neighbours of atom 1 in a 4 Å shell
shell = within(4.0, atom(1));

# Between 2 and 4 Å of each residue's centre
ring_shell = within(2.0:4.0, com(residue(:)));

# Neighbours of every residue, measured per residue
n_neighbours = count(within(3.0) in residue(1:5));
```

### within_x

<!-- proc name=within_x category=selector.spatial -->

```text
within_x(range: frange) -> atomset
```

Atoms whose x coordinate lies in `range` (Å, in the coordinate system of the trajectory).

```mdscript
slab = within_x(20:25);
left_half = within_x(:30.0);
```

### within_y

<!-- proc name=within_y category=selector.spatial -->

```text
within_y(range: frange) -> atomset
```

Atoms whose y coordinate lies in `range`.

```mdscript
slab = within_y(10:15);
```

### within_z

<!-- proc name=within_z category=selector.spatial -->

```text
within_z(range: frange) -> atomset
```

Atoms whose z coordinate lies in `range`.

```mdscript
slab = within_z(20:25);
n_in_slab = count(within_z(20:25));
```

### within_xyz

<!-- proc name=within_xyz category=selector.spatial -->

```text
within_xyz(x: frange, y: frange, z: frange) -> atomset
```

**Parameters:** `x`, `y`, `z` (all required).

Atoms inside an axis-aligned box.

```mdscript
box = within_xyz(20:25, 20:25, 20:25);
n_in_box = count(within_xyz(20:25, 20:25, 20:25));
```

---

## Bond sets

Bonds are values of their own, a `bondset`: a set of atoms cannot say which bonds are meant, as the atoms of a bond
bring in every other bond among them. A bond is one of the system's persistent bonds, which the bond graph holds
(covalent and coordination bonds, from which the rings, structures and chemistry are derived) together with the
non-covalent bonds a file gives or a user marks (hydrogen bonds, interactions), or a bond computed for a frame (one
[`formed`](#formed) in it). Bond sets combine with `and`, `or`, `xor` and `not` like selections, and give atoms with
[`atom`](#atom), numbers with [`count`](#count) and [`length`](#length).

The kinds, and what they select:

| Kind | Bonds |
|---|---|
| `"covalent"` | Of the bond graph, not coordination |
| `"coordinate"` | Metal-ligand bonds |
| `"elastic"` | Springs of the force field which hold a coarse grained molecule in shape (an elastic network), or restraints: kept with the bonds, but no bonds of the molecule |
| `"metallic"` | Between the atoms of a metal in bulk, at a surface or in a particle, to every neighbour of the first coordination shell: kept with the bonds, but no bonds of a molecule. What touches the metal has no bonds to it but those formed in a frame (`bond_delta`) |
| `"noncovalent"` | Every non-covalent bond: `"interaction"` and `"hydrogen"` |
| `"interaction"` | Non-covalent bonds of no particular kind: one a user marks between two atoms, a salt bridge, a halogen bond |
| `"hydrogen"` | Hydrogen bonds |

A bond is of one kind, so a charge-assisted hydrogen bond is a `"hydrogen"` bond.

### bond

<!-- proc name=bond category=bonds -->

```text
bond(a: atomset[],   b: atomset[], kind: string[], order: irange, resonance: string[], origin: string[]) -> bondset[]
bond(d: bond_delta,  b: atomset[], kind: string[], order: irange, resonance: string[], origin: string[]) -> bondset
```

**Parameters:** `a` or `d`, and `b` (optional), `kind`, `order`, `resonance` and `origin` (optional, by name).

The persistent bonds of the system which every filter given lets through, within the context (both atoms in it).

- Without `a` and `b`: every bond.
- With `a`: the bonds with both atoms in `a`, one set per element of `a`, so `bond(residue(:))` is the bonds within
  each residue.
- With `a` and `b`: the bonds with one atom in the element of `a` and the other in `b`.
- With `b` alone: the bonds with an atom in `b`.

An atom index stands for that atom, as for [`distance`](#distance): `bond(1, 4)` is the bond between the first and the
fourth atom, if they have one.

With a [`bond_delta`](#bond_delta) `d` in place of `a`: the bonds of the frame within its scope, the persistent bonds
with an atom in it which are not broken and those formed, with an atom in `b` when it is given. A bond formed is of the
`"inferred"` origin, covalent or coordinate, and of no known order.

`kind` is any of the kinds above. `order` is the order of a bond of the graph, a number or a range: 1 single, 2 double,
3 triple, 4 quadruple, 0 unknown; non-covalent bonds have none. `resonance` is any of `"aromatic"`, `"delocalized"`
(carboxylate, guanidinium, nitro, ...) and `"localized"` (neither), of a bond of the graph. `origin` is any of `"file"`,
`"topology"`, `"user"` and `"inferred"`. An unknown name is a compile error which lists the valid ones.

```mdscript
metal_sites = bond(kind="coordinate");
disulfides  = bond(element("S"), element("S"), kind="covalent");
adducts     = bond(protein(), not protein(), kind="covalent");   # covalent ligands, crosslinks
drawn       = bond(kind="hydrogen", origin="user");               # the hydrogen bonds a user marked
kekule      = bond(resonance="aromatic", order=2);
n_per_res   = count(bond(residue(1:10)));
one         = bond(1, 4);                                         # the bond between atoms 1 and 4
d           = bond_delta(resid(10:20));
lengths     = length(bond(d));                                    # the bonds of each frame, within resid 10 to 20
metal_now   = bond(d, kind="coordinate");                         # those to a metal, formed ones included
```

### bond_delta

<!-- proc name=bond_delta category=bonds -->

```text
bond_delta(scope: atomset[], tolerance: float, breakable: string[]) -> bond_delta
```

**Parameters:** `scope` (optional), `tolerance` (default 0.1) and `breakable` (default `{"file", "inferred"}`), by name.

The bonds formed and broken in each frame, against the bond graph of the system, which is the reference: only the
bonds with an atom in `scope` are looked at (all of them without it), so a bond to an atom just outside of it forms
too. A bond of the graph breaks once it is longer than its cutoff by `tolerance`, and a pair forms a bond once it is
shorter than its cutoff by `tolerance` (the cutoffs of the bond inference, from the covalent radii). The band keeps a
bond which vibrates near its cutoff from flickering; it is measured against the graph rather than an earlier frame,
so a frame gives the same bonds whatever was evaluated before it. Only the bonds of the `breakable` origins break: a
bond a user added stays, and so by default does one of a force field topology, as in a classical simulation. A bond
forming as another breaks (a proton moving, the transition state of a substitution) shows both for the frames it
takes. A pair with a non-covalent bond in the system forms no bond. In a coarse grained system no bonds form or
break: it is a compile error.

The result is read with [`formed`](#formed), [`broken`](#broken) and [`bond`](#bond).

```mdscript
d = bond_delta(within(6.0, resname("LIG")));
n_formed = count(formed(d));
n_broken = count(broken(d));
```

### formed

<!-- proc name=formed category=bonds -->

```text
formed(d: bond_delta) -> bondset
```

**Parameters:** `d` (required).

The bonds formed in the frame, which the system does not have: computed bonds, of the inferred origin.

```mdscript
new_bonds = formed(bond_delta(resid(10:20)));
```

### broken

<!-- proc name=broken category=bonds -->

```text
broken(d: bond_delta) -> bondset
```

**Parameters:** `d` (required).

The bonds of the bond graph broken in the frame.

```mdscript
lost = count(broken(bond_delta()));
```

---

## Properties

Procedures that measure something. Most of them are **dynamic**: they are evaluated for every frame and give a
temporal property. `rdf`, `density_*` and `sdf` accumulate over all evaluated frames into a distribution or volume.

### distance

<!-- proc name=distance category=property -->

```text
distance(a: position[], b: position[]) -> float [Å]
```

**Parameters:** `a`, `b` (both required).

The distance between two points. A selection is reduced to its mass-weighted centre of mass, and an array of
selections is merged first (see [Coordinate arguments](#coordinate-arguments)). Periodic boundary conditions are
respected when the system has a unit cell. Use a [context](#contexts-in-and-out) to measure something per residue.

```mdscript
d_atoms = distance(1, 5);                                       # atom 1 to atom 5
d_com   = distance(com(residue(1)), com(residue(2)));           # between two residue centres
d_point = distance(residue(1), vec3(0, 0, 0));                  # residue 1 to the origin
per_res = distance(1, 2) in residue(:);                         # float[num_residues]
```

### distance_min

<!-- proc name=distance_min category=property -->

```text
distance_min(a: position[N], b: position[]) -> float[N] [Å]
```

**Parameters:** `a`, `b`.

One value per element of `a`: the shortest distance from any atom of that element to any atom of `b` (the closest
approach of each element to `b`). The elements of `a` are the selections of an array of selections (`residue(1:3)`
gives three values, `protein()` one per protein residue), the atoms of a list of indices (`{1, 5, 9}` gives three
values) or the points of a list of `float[3]`; a single selection or index range is one element and gives a single
value. `b` is taken as a whole: an array of selections is merged into one set of atoms. Atoms are measured, not
centres of mass, and periodic boundary conditions are respected when the system has a unit cell.

The number of values does not change when the selections change size from frame to frame (`within(...)`), and an
element without atoms, or an empty `b`, gives 0. If the *number* of elements changes from frame to frame
(`residue(within(...))`) the result is not published as a property.

```mdscript
closest = distance_min(residue(1), residue(2));                 # one value
per_res = distance_min(residue(1:10), residue(20));             # float[10]: residues 1 to 10, each to residue 20
to_ligand = distance_min(protein(), resname("LIG"));            # one value per protein residue
gap = distance_min(within_x(:20), within_x(30:));
```

### distance_max

<!-- proc name=distance_max category=property -->

```text
distance_max(a: position[N], b: position[]) -> float[N] [Å]
```

**Parameters:** `a`, `b`.

One value per element of `a`: the largest distance from any atom of that element to any atom of `b` (the extent of
each element and `b` together). The elements of `a` and the handling of `b` are those of
[`distance_min`](#distance_min).

```mdscript
extent = distance_max(residue(1), residue(2));                  # one value
per_res = distance_max(residue(1:10), residue(20));             # float[10]
```

### distance_pair

<!-- proc name=distance_pair category=property -->

```text
distance_pair(a: position[], b: position[]) -> float[N*M] [Å]
```

**Parameters:** `a`, `b`.

Every pairwise distance between the `N` positions of `a` and the `M` positions of `b`, laid out row by row (`a`
varies slowest). A selection contributes each of its atoms; an array of selections contributes one centre of mass per
element. The number of positions must be known when the script is compiled, so dynamic selections such as
`within(...)` are not allowed here.

```mdscript
all_pairs = distance_pair(atom(1:3), atom(4:6));               # 9 values
to_origin = distance_pair(residue(1:3), vec3(0, 0, 0));        # 3 values, one per residue
```

### angle

<!-- proc name=angle category=property -->

```text
angle(a: position[], b: position[], c: position[]) -> float [rad]
angle(a: plane[], b: plane[]) -> float[] [rad]
```

**Parameters:** `a`, `b`, `c` (three positions), or `a`, `b` (two sets of planes, from [`plane`](#plane)).

The angle at `b` between the directions to `a` and `c`, in `[0, π]`. Of two planes, the angle between them, in
`[0, π/2]`: one value per pair, a single plane against each of the other.

```mdscript
bond_angle = angle(1, 2, 3);
per_res = angle(1, 2, 3) in residue(:);
in_degrees = angle(1, 2, 3) * 180.0 / PI;
```

### dihedral

<!-- proc name=dihedral category=property -->

```text
dihedral(a: position[], b: position[], c: position[], d: position[]) -> float [rad]
```

**Parameters:** `a`, `b`, `c`, `d` (all required).

The signed torsion angle around the axis `b`–`c`, in `(-π, π]`.

```mdscript
phi = dihedral(1, 2, 3, 4);
side_chain_torsion = dihedral(2, 3, 6, 10) in resname("ALA");
```

### rmsd

<!-- proc name=rmsd category=property -->

```text
rmsd(sel: atomset[]) -> float[] [Å]
```

The mass-weighted root-mean-square deviation of the selected atoms from the system's reference structure (the first
frame), after removing translation and finding the optimal rotation. Each atomset of the argument is a structure of
its own, fitted and reported separately: the result has one value per atomset. Unlike most procedures, which merge
an array of selections into one, `rmsd` keeps them apart; use [`union`](#union) to fit them together as one
structure. A single selection (`backbone()`, `within(...)`, a combination with `and`/`or`) gives a single value.

The selection may change from frame to frame, and an empty one gives 0. If the *number* of selections changes from
frame to frame (`residue(within(...))`) the result is not published as a property; take its `union`. Inside a context
only the atoms of the context are used. Periodic boundaries are handled as long as a structure spans less than half
the cell.

```mdscript
backbone_rmsd = rmsd(backbone());               # one value
chain_rmsd = rmsd(chain(:));                    # one value per chain, each fitted on its own
residue_rmsd = rmsd(protein());                 # one value per protein residue
protein_rmsd = rmsd(union(protein()));          # the whole protein as one structure
ligand_rmsd = rmsd(resname("LIG"));             # one value per LIG residue
```

### rdf

<!-- proc name=rdf category=property -->

```text
rdf(a: position[], b: position[], cutoff: float)  -> distribution
rdf(a: position[], b: position[], cutoff: frange) -> distribution
```

**Parameters:** `a`, `b`, `cutoff` (all required).

The radial distribution function of `b` around `a`: a histogram with 1024 bins between `0` and `cutoff` (or the
two ends of a range), normalised per shell so that a uniform distribution gives a constant value, and averaged over
all frames. `cutoff` must be a constant. Zero distances are ignored, so an atom is not counted against itself.
If `a` is an array of selections, the reference points are the centre of mass of each one.

```mdscript
# Carbon–oxygen pair distribution out to 10 Å
g_co = rdf(element("C"), element("O"), 10.0);

# Restrict to a shell
g_shell = rdf(atom(1), all(), 2.0:10.0);

# Around the centre of mass of each residue
g_res = rdf(com(all()) in residue(:), element("O"), 20.0);
```

### density_x

<!-- proc name=density_x category=property -->

```text
density_x(sel: atomset[]) -> distribution [kg/m³]
```

Mass density of the selection along the x axis of the simulation cell, as a 1024-bin histogram accumulated over
all frames. The selection is merged, so an array of selections is treated as one. Requires a unit cell.

```mdscript
rho_x = density_x(protein());
```

### density_y

<!-- proc name=density_y category=property -->

```text
density_y(sel: atomset[]) -> distribution [kg/m³]
```

As [`density_x`](#density_x), along the y axis.

```mdscript
rho_y = density_y(water());
```

### density_z

<!-- proc name=density_z category=property -->

```text
density_z(sel: atomset[]) -> distribution [kg/m³]
```

As [`density_x`](#density_x), along the z axis. The classic use is the density profile across a membrane or slab.

```mdscript
rho_z = density_z(all());
rho_water = density_z(water());
```

### density

<!-- proc name=density category=property -->

```text
density(sel: atomset[]) -> distribution[3] [kg/m³]
```

The three axis profiles (x, y, z) of the selection in a single value. It is published as an array of distributions,
which is not plotted yet; [`density_x`](#density_x), [`density_y`](#density_y) and [`density_z`](#density_z) give
one axis each.

```mdscript
rho = density(all());
```

### pool

<!-- proc name=pool category=property -->

```text
pool(a: distribution[], b: distribution[]) -> distribution[]
```

**Parameters:** `a`, `b` (both required, the same length and bins).

The two distributions as one: values and weights added bin by bin, as if their samples had been taken together. Of
two `rdf` over different reference groups, it is the `rdf` over both. It is not the sum of the two curves, which a
distribution cannot give (see [Distributions](#distributions)).

```mdscript
g_ab = pool(rdf(element("C"), element("O"), 10.0), rdf(element("N"), element("O"), 10.0));
```

### sdf

<!-- proc name=sdf category=property -->

```text
sdf(structures: atomset[], target: atomset, extent: float) -> volume [count/Å³]
```

**Parameters:** `structures`, `target`, `extent` (all required).

The spatial distribution function: every element of `structures` is superimposed onto the first, and the density of
`target` atoms within `extent` Å of it is accumulated on a 128×128×128 grid over all frames. The elements of
`structures` must be equivalent (same atoms in the same order), which is checked when the script is compiled.
The result is a `volume` property that can be rendered as an isosurface.

```mdscript
# Distribution of oxygen atoms around every PFT molecule
sdf1 = sdf(resname("PFT"), element("O"), 20.0);

# Around each chain, considering only residues named "PFT"
sdf2 = sdf(chain(:), resname("PFT"), 30.0);
```

### count

<!-- proc name=count category=property -->

```text
count(sel: atomset[])                  -> float
count(sel: atomset[], unit: string)    -> float
count(c: contact)                       -> float
count(c: contact, unit: string)         -> float
count(b: bondset[])                     -> float
```

**Parameters:** `sel` or `c` (required), `unit` (optional).

The number of atoms in the selection (an array is merged first). With `unit` the number of larger structures the
selection touches instead: `"atom"`, `"residue"`, `"part"`, `"chain"` or `"molecule"`. An unknown unit is a compile error that
lists the valid ones. A `count` of a dynamic selection changes over time and is a temporal property.

Given a bond set from [`bond`](#bond), `count` is the number of bonds (an array is merged first).

Given a contact set from [`contacts`](#contacts), `count` is the number of group pairs in contact. With `unit`
`"group"` it is the same, and with `"atom"` the number of particle pairs behind them.

```mdscript
n_atoms = count(within(4.0, atom(1)));                    # temporal: neighbours of atom 1
n_res = count(protein(), "residue");                      # number of residues
n_hydrated = count(within(3.5, protein()) and element("O"));
```

### contact_count

<!-- proc name=contact_count category=property -->

```text
contact_count(a: atomset[], b: atomset[], cutoff: float) -> float[N]
```

**Parameters:** `a`, `b`, `cutoff` (all required).

For each element of `a`, the number of atom pairs `(atom in that element, atom in b)` that are closer than `cutoff`
Å. Atoms in `a` itself and atoms up to four bonds away from it are never counted as contacts, so the covalent
neighbourhood does not contribute. `b` is merged into one selection. The number of elements of `a` must be known
when the script is compiled.

```mdscript
contacts = contact_count(residue(1), residue(2:3), 4.0);
```

### contacts

<!-- proc name=contacts category=property -->

```text
contacts(a: atomset[], b: atomset[], cutoff: float, exclude_bonds: int, min_separation: int,
         parent: atomset[], exclude_within: atomset[]) -> contact
```

**Parameters:** `a` and `cutoff` (required), `b`, `exclude_bonds` (default 3), `min_separation` (default 0),
`parent` and `exclude_within` (optional). Everything after `a` and `b` is given by name.

The groups in contact at each frame: pairs of groups with at least one pair of particles closer than `cutoff` Å.
A group is any selection - a residue, a chain, a slice of a fibril - and groups may overlap.

- Without `b`: the pairs of distinct groups of `a`, each pair once. Particle pairs within one group never count,
  so `a` has to be several groups (`residue(:)`, not `protein()`).
- With `b`: every group of `a` against every group of `b`.

`exclude_bonds` skips particle pairs joined by at most that many bonds (3 skips 1-2, 1-3 and 1-4 pairs, as a force
field's exclusions do; 0 skips nothing). Within one set, `min_separation` skips pairs of groups whose indices differ
by less than it (3 ignores a residue's nearest neighbours along its chain); `parent` says which groups share a
sequence, and groups of different parents are never neighbours. `exclude_within` skips particle pairs that lie in
the same element of it, for example the same chain or fibril.

The groups and every parameter have to be the same for every frame, so the query is prepared once when the script
is compiled. The result is a contact set, which is reduced with [`count`](#count) or [`degree`](#degree).

```mdscript
c = contacts(residue(:), cutoff=4.5);
n_pairs = count(c);                                     # residue pairs in contact
n_atom_pairs = count(c, "atom");                        # the particle pairs behind them

# Between two sets, without excluding bonded neighbours
c2 = contacts(residue(1:5), residue(6:15), cutoff=4.0, exclude_bonds=0);

# Slices of five residues, ignoring their near neighbours and contacts within a chain
c3 = contacts(chunks(residue(:), 5), cutoff=5, min_separation=3, parent=residue(:), exclude_within=chain(:));
```

### degree

<!-- proc name=degree category=property -->

```text
degree(c: contact) -> float[N]
```

**Parameters:** `c` (required).

For each group of `a` in a contact set from [`contacts`](#contacts), the number of groups it is in contact with.
The length is the number of groups of `a`, which is known when the script is compiled.

```mdscript
c = contacts(residue(:), cutoff=4.5);
d = degree(c);                                          # contacts per residue, per frame
```

### porosity

<!-- proc name=porosity category=property -->

```text
porosity(sel: atomset) -> float
```

The void ratio of the selection inside its bounding box: the fraction of the box not occupied by atoms, in `[0, 1]`.

```mdscript
void_fraction = porosity(all());
```

---

## Geometry

Procedures that turn atoms into points, vectors and planes. Their results are values you feed into other
procedures (`com` into `distance`, `plane` into `angle`) rather than properties to plot; those that result in a plain
`float` or `float[3]` are published as properties when they change over the trajectory.

### com

<!-- proc name=com category=geometry -->

```text
com(a: position[]) -> float[3] [Å]
```

The mass-weighted centre of mass of the positions (an array of selections is merged first). Periodic boundary
conditions are taken into account so that molecules that span the cell edge get a correct centre.

```mdscript
centre = com(residue(1));
d = distance(com(residue(1)), com(residue(2)));
centre_of_selection = com(within(4.0, atom(1)));
```

### plane

<!-- proc name=plane category=geometry -->

```text
plane(groups: atomset[]) -> plane[]
plane(points: position[]) -> plane
```

The least squares plane of a group of atoms, not weighted by mass: one plane per atomset of `groups`, or one through
`points` (atom indices, or `float[3]` positions such as [`coord`](#coord) gives). Each needs at least three atoms
or positions. A plane record holds what [`centroid`](#centroid), [`normal`](#normal),
[`radius`](#radius), [`deviation`](#deviation) and [`angle`](#angle) read of it, and hovering it shows them.

Like [`rmsd`](#rmsd), `plane` keeps the groups apart: `plane(residue(1:10))` is ten planes. For one plane through all
of their atoms take their [`union`](#union) first, `plane(union(residue(1:10)))`; for one through their centres pass
the centres, `plane(coord(residue(1:10)))`. The groups cannot depend on the frame; positions can
(`plane(coord(within(...)))`).

The normal is oriented, so that it keeps to the atoms from frame to frame and a time series of it does not flip:

- A group which is a ring of the system ([`ring`](#ring)): the normal is on the side from which the ring goes around
  anticlockwise, in its order. A sugar ring goes from the ring O to the anomeric C.
- Any other group: by two of its atoms chosen in the reference structure, the one farthest from the centroid and the
  one which spans the largest area with it.
- Positions: on the side from which they go around anticlockwise in the order given, so give them in order around
  the plane (`plane({1, 2, 5, 4})`).

```mdscript
p = plane(ring(kind="aromatic"));                  # one per aromatic ring
tilt = angle(plane(ring("c1ccccc1") in resid(10)), plane(ring("c1ccccc1") in resid(20)));
flat = deviation(plane(residue(1:10)));            # how far each residue is from flat, in Å
whole = plane(union(resname("PFT")));              # one plane through all PFT atoms
per_res = plane(1:4) in resname("PFT");            # atoms 1 to 4 of each PFT residue, in that order
sheet = plane(coord(residue(1:10)));               # one plane through the centres of the residues
```

### centroid

<!-- proc name=centroid category=geometry -->

```text
centroid(p: plane[]) -> float[3][]
```

The centre of each plane: the mean position of its atoms, in Å. A `position`, so `distance_pair(centroid(p), ...)`
measures from it.

### normal

<!-- proc name=normal category=geometry -->

```text
normal(p: plane[]) -> float[3][]
```

The unit normal of each plane, oriented as [`plane`](#plane) describes.

### radius

<!-- proc name=radius category=geometry -->

```text
radius(p: plane[]) -> float[] [Å]
```

The mean distance of each plane's atoms from its centroid: the radius of a ring.

### deviation

<!-- proc name=deviation category=geometry -->

```text
deviation(p: plane[]) -> float[] [Å]
```

The root mean square distance of each plane's atoms from it: 0 for a flat group, a few tenths of an Å for a puckered
ring.

### pucker

<!-- proc name=pucker category=geometry -->

```text
pucker(rings: atomset[]) -> pucker[]
pucker(atoms: int[]) -> pucker
```

How each ring of 5 to 7 atoms is puckered: its Cremer-Pople coordinates (Cremer & Pople, 1975), read with
[`amplitude`](#amplitude), [`phase`](#phase) and [`theta`](#theta), the pseudorotation of a furanose
(Altona & Sundaralingam, 1972), read with [`pseudorotation`](#pseudorotation), and the name of its conformation, which
hovering it shows (`4C1 chair`, `3E envelope, C3'-endo`). The rings are rings of the system ([`ring`](#ring)), and
cannot depend on the frame; `pucker({i, j, ...})` takes the atoms of one ring in the order given instead.

What the phase and theta mean depends on where the ring starts and which way it goes, so the atoms are taken in the
order the conventions are written for:

- A sugar ring from its ring O and the anomeric C (O5, C1, ..., C5 of a pyranose; O4', C1', ..., C4' of the furanose of
  a nucleotide). A 4C1 chair is at theta 0, a 1C4 chair at theta 180.
- Another ring from the heaviest atom other than C (the N of proline, then CA), or from the C with the most atoms
  bonded outside of the ring, towards the neighbour with the most of those.
- A ring alike all around (benzene, cyclohexane) has no such atom: its amplitude and the kind of its conformation (a
  chair, a boat) mean something, its phase and theta depend on which atom happened to be taken first. `phase` and
  `theta` warn of it.

The conformations are those of IUPAC (2-Carb-7): the 38 of a six membered ring (two chairs, six boats, six skews,
twelve half chairs, twelve envelopes) by the nearest on the Cremer-Pople sphere, the 20 of a five membered ring
(envelopes and twists) by phi2. A ring with an amplitude under 0.1 Å is planar. Seven membered rings get their
coordinates but no name.

```mdscript
p = pucker(ring(kind="sugar"));
P = pseudorotation(pucker(ring(kind="sugar", size=5)));     # nucleotides: C3'-endo near 18, C2'-endo near 162 degrees
t = theta(pucker(ring(kind="sugar", size=6)));              # pyranoses: 4C1 near 0, 1C4 near 180 degrees
q = amplitude(pucker(ring("CCNCC", size=5)));               # pyrrolidines, prolines
```

### amplitude

<!-- proc name=amplitude category=geometry -->

```text
amplitude(p: pucker[]) -> float[] [Å]
```

The total Cremer-Pople puckering amplitude Q of each ring: the root of the sum of the squared displacements of its
atoms from their mean plane. The same whatever the order of the atoms.

### phase

<!-- proc name=phase category=geometry -->

```text
phase(p: pucker[]) -> float[] [rad]
```

The Cremer-Pople phase phi2 of each ring, in `[0, 2π)`: where along the pseudorotation of a five membered ring it is,
or the longitude of a six membered ring on the sphere.

### theta

<!-- proc name=theta category=geometry -->

```text
theta(p: pucker[]) -> float[] [rad]
```

The Cremer-Pople theta of each six membered ring, in `[0, π]`: 0 and π the two chairs, π/2 the boats and skews. An
error for rings of another size.

### pseudorotation

<!-- proc name=pseudorotation category=geometry -->

```text
pseudorotation(p: pucker[]) -> float[] [rad]
```

The Altona-Sundaralingam phase P of each furanose, in `[0, 2π)`, the convention of nucleic acids: about 18 degrees
for C3'-endo (A-form), 162 for C2'-endo (B-form). P is phi2 shifted by a quarter turn. An error for a ring which is not a
five membered sugar ring, unless its five atoms were given in order (O4', C1', C2', C3', C4').

### shape_weights

<!-- proc name=shape_weights category=geometry -->

```text
shape_weights(a: position[]) -> float[3]
```

How linear, planar or isotropic the point cloud is, from the eigenvalues of its gyration tensor. The three weights
sum to 1: `(linear, planar, isotropic)`. Destructure the result to get them as separate properties.

```mdscript
{lin, plan, iso} = shape_weights(all());
```

### coord

<!-- proc name=coord category=geometry -->

```text
coord(a: position[]) -> float[3][N] [Å]
```

The coordinates of each position: every atom of a selection, or one centre of mass per element of an array of
selections. Not published as a property on its own; use it to build vectors.

```mdscript
xyz = coord(atom(1:2));
```

### coord_x

<!-- proc name=coord_x category=geometry -->

```text
coord_x(a: position[]) -> float[N] [Å]
```

The x coordinate of each position (see [`coord`](#coord)).

```mdscript
x = coord_x(residue(1));                 # x of every atom in residue 1
x_com = coord_x(residue(1:3));           # x of the centre of mass of each of three residues
```

### coord_y

<!-- proc name=coord_y category=geometry -->

```text
coord_y(a: position[]) -> float[N] [Å]
```

The y coordinate of each position.

```mdscript
y = coord_y(atom(1:5));
```

### coord_z

<!-- proc name=coord_z category=geometry -->

```text
coord_z(a: position[]) -> float[N] [Å]
```

The z coordinate of each position.

```mdscript
z = coord_z(atom(1:5));
height = max(coord_z(protein()));
```

### coord_xy

<!-- proc name=coord_xy category=geometry -->

```text
coord_xy(a: position[]) -> float[2][N] [Å]
```

The x and y coordinates of each position. Not published as a property.

```mdscript
xy = coord_xy(atom(1:3));
```

### coord_xz

<!-- proc name=coord_xz category=geometry -->

```text
coord_xz(a: position[]) -> float[2][N] [Å]
```

The x and z coordinates of each position. Not published as a property.

```mdscript
xz = coord_xz(atom(1:3));
```

### coord_yz

<!-- proc name=coord_yz category=geometry -->

```text
coord_yz(a: position[]) -> float[2][N] [Å]
```

The y and z coordinates of each position. Not published as a property.

```mdscript
yz = coord_yz(atom(1:3));
```

---

## Math

Scalar functions take and return `float` (an `int` is converted implicitly). A few also work on whole arrays.
Note that these functions return a plain number without a unit tag (see [Types and units](#types-and-units)).

```mdscript
# Normalised distance and a Gaussian weight
d = distance(1, 5);
w = exp(-pow(d - 1.5, 2.0));
```

### sqrt

<!-- proc name=sqrt category=math -->

```text
sqrt(x: float) -> float
```

Square root.

```mdscript
v = sqrt(distance(1, 2));
```

### cbrt

<!-- proc name=cbrt category=math -->

```text
cbrt(x: float) -> float
```

Cube root.

```mdscript
v = cbrt(count(all()) * 1.0);
```

### abs

<!-- proc name=abs category=math -->

```text
abs(x: float) -> float
abs(x: float[]) -> float[]
abs(x: volume) -> volume
```

Absolute value. Also works element-wise on arrays and on volumes.

```mdscript
v = abs(dihedral(1, 2, 3, 4));
```

### floor

<!-- proc name=floor category=math -->

```text
floor(x: float) -> float
floor(x: float[]) -> float[]
```

Largest integer not greater than the argument. Also element-wise on arrays.

```mdscript
v = floor(distance(1, 2));
```

### ceil

<!-- proc name=ceil category=math -->

```text
ceil(x: float) -> float
ceil(x: float[]) -> float[]
```

Smallest integer not less than the argument. Also element-wise on arrays.

```mdscript
v = ceil(distance(1, 2));
```

### sin

<!-- proc name=sin category=math -->

```text
sin(x: float) -> float
```

Sine of an angle in radians.

```mdscript
v = sin(angle(1, 2, 3));
```

### cos

<!-- proc name=cos category=math -->

```text
cos(x: float) -> float
```

Cosine of an angle in radians.

```mdscript
v = cos(angle(1, 2, 3));
```

### asin

<!-- proc name=asin category=math -->

```text
asin(x: float) -> float
```

Inverse sine, result in radians.

```mdscript
v = asin(0.5);
```

### acos

<!-- proc name=acos category=math -->

```text
acos(x: float) -> float
```

Inverse cosine, result in radians.

```mdscript
v = acos(cos(angle(1, 2, 3)));
```

### atan

<!-- proc name=atan category=math -->

```text
atan(x: float) -> float
atan(y: float, x: float) -> float
```

Inverse tangent, result in radians. With two arguments it is the same as `atan2`.

```mdscript
v = atan(1.0, 2.0);
```

### log

<!-- proc name=log category=math -->

```text
log(x: float) -> float
```

Natural logarithm.

```mdscript
v = log(distance(1, 2));
```

### log2

<!-- proc name=log2 category=math -->

```text
log2(x: float) -> float
```

Base-2 logarithm.

```mdscript
v = log2(distance(1, 2));
```

### log10

<!-- proc name=log10 category=math -->

```text
log10(x: float) -> float
```

Base-10 logarithm.

```mdscript
v = log10(distance(1, 2));
```

### exp

<!-- proc name=exp category=math -->

```text
exp(x: float) -> float
```

`e` raised to the argument.

```mdscript
v = exp(-distance(1, 2));
```

### exp2

<!-- proc name=exp2 category=math -->

```text
exp2(x: float) -> float
```

2 raised to the argument.

```mdscript
v = exp2(distance(1, 2));
```

### atan2

<!-- proc name=atan2 category=math -->

```text
atan2(y: float, x: float) -> float
```

Two-argument inverse tangent in `(-π, π]`, using the signs of both arguments to pick the quadrant.

```mdscript
heading = atan2(coord_y(atom(1)) - coord_y(atom(2)), coord_x(atom(1)) - coord_x(atom(2)));
```

### pow

<!-- proc name=pow category=math -->

```text
pow(x: float, y: float) -> float
```

`x` raised to the power `y`.

```mdscript
squared = pow(distance(1, 2), 2.0);
```

### min

<!-- proc name=min category=math -->

```text
min(a: float, b: float)   -> float
min(x: float[])           -> float
min(a: volume, b: float)  -> volume
min(a: volume, b: volume) -> volume
```

The smaller of two numbers, or the smallest element of an array (with one argument), in the unit of the array.

```mdscript
closest = min(distance(1, 2), distance(1, 3));
smallest = min(distance(1, 2) in residue(:));
clipped = min(distance(1, 2), 5.0);
```

### max

<!-- proc name=max category=math -->

```text
max(a: float, b: float)   -> float
max(x: float[])           -> float
max(a: volume, b: float)  -> volume
max(a: volume, b: volume) -> volume
```

The larger of two numbers, or the largest element of an array (with one argument), in the unit of the array. The volume overloads are
useful to clamp an `sdf`.

```mdscript
farthest = max(distance(1, 2), distance(1, 3));
largest = max(distance(1, 2) in residue(:));
floor_at_zero = max(dihedral(1, 2, 3, 4), 0.0);
```

### sum

<!-- proc name=sum category=math -->

```text
sum(x: float[]) -> float
```

The sum of the elements of an array, 0 for an empty one, in the unit of the array. With an `attr` pattern it adds
up every attribute the pattern matches, at every frame.

```mdscript
total = sum(distance(1, 2) in residue(:));
# coulomb = sum(attr("edr/coul_sr_*"));    # requires the terms in the loaded system
```

---

## Linear algebra and constructors

Vectors are `float[2]`, `float[3]` and `float[4]`; matrices are `float[4][4]`.

### vec2

<!-- proc name=vec2 category=linalg -->

```text
vec2(x: float, y: float) -> float[2]
```

```mdscript
p = vec2(1.0, 2.0);
```

### vec3

<!-- proc name=vec3 category=linalg -->

```text
vec3(x: float, y: float, z: float) -> float[3]
```

A point or direction; also usable wherever a `position` is expected.

```mdscript
origin = vec3(0, 0, 0);
d_origin = distance(com(residue(1)), vec3(0, 0, 0));
```

### vec4

<!-- proc name=vec4 category=linalg -->

```text
vec4(x: float, y: float, z: float, w: float) -> float[4]
```

```mdscript
q = vec4(0, 1, 13, 5);
```

### dot

<!-- proc name=dot category=linalg -->

```text
dot(a: float[], b: float[]) -> float
```

The dot product of two vectors of equal length.

```mdscript
along_y = dot(vec3(1, 2, 3), vec3(0, 1, 0));
```

### cross

<!-- proc name=cross category=linalg -->

```text
cross(a: float[3], b: float[3]) -> float[3]
```

The cross product of two 3-vectors.

```mdscript
n = cross(vec3(1, 0, 0), vec3(0, 1, 0));
area = length(cross(vec3(1, 0, 0), vec3(0, 2, 0)));
```

### length

<!-- proc name=length category=linalg -->

```text
length(v: float[])   -> float
length(b: bondset[]) -> float[] [Å]
```

The Euclidean norm of a vector. Given a bond set, the length of each of its bonds in the frame (the shortest across
the periodic cell), in the order of the set: the persistent bonds by index, then the bonds computed for the frame. An
array of sets is merged first.

```mdscript
norm = length({3.0, 4.0});                     # 5
```

### normalize

<!-- proc name=normalize category=linalg -->

```text
normalize(v: float[]) -> float[]
```

The vector scaled to unit length.

```mdscript
unit = normalize(vec3(3, 0, 4));
```

### mul

<!-- proc name=mul category=linalg -->

```text
mul(a: float[4][4], b: float[4][4]) -> float[4][4]
mul(m: float[4][4], v: float[4])    -> float[4]
```

Matrix–matrix and matrix–vector product. There is no matrix constructor, so matrices come from data such as
[`attr`](#attr).

```mdscript
# Apply a 4x4 transform that is stored per frame in the trajectory (requires such an attribute)
# moved = mul(flatten(attr("run/a/transform")), vec4(1, 0, 0, 1));
```

---

## Combining and reshaping

### union

<!-- proc name=union aliases=join category=combine -->
<a id="join"></a>

```text
union(sel: atomset[]) -> atomset
```

The atoms of all of the selections, each once: one selection made of an array of them, `or` over all of its
elements. Where a procedure keeps the selections of an array apart ([`rmsd`](#rmsd), [`plane`](#plane),
[`pucker`](#pucker)), `union` makes them one structure. `join` is its old name, which still compiles with a warning.

```mdscript
all_protein = union(protein());
protein_rmsd = rmsd(union(protein()));          # the whole protein as one structure
```

### flatten

<!-- proc name=flatten category=combine -->

```text
flatten(x: T[..]) -> T[]          # any array, one-dimensional result
```

The elements of a multi-dimensional array in one dimension: the per-residue results of a context, or a tensor read
with [`attr`](#attr). Of selections it gives their [`union`](#union), with a warning: that is not a change of shape
(an atom in two of them is in the result once), and `union` says what it is.

```mdscript
per_res = flatten(distance(1, 2) in residue(:));
```

### split

<!-- proc name=split category=combine -->

```text
split(sel: atomset, parts: int) -> atomset[]
```

**Parameters:** `sel`, `parts` (both required).

Cut a selection into `parts` pieces of (nearly) equal size by atom index. `parts` must be a positive constant.

```mdscript
halves = split(residue(1), 2);
d = distance(1, 2) in split(chain(1), 4);              # measure each quarter of a chain
```

### chunks

<!-- proc name=chunks category=combine -->

```text
chunks(sel: atomset[], size: int) -> atomset[]
```

**Parameters:** `sel`, `size` (both required).

Cut each element of `sel` into consecutive runs of `size` atoms, in atom index order; a last run shorter than `size`
is kept. Where [`split`](#split) makes a given number of parts, `chunks` makes parts of a given size: groups the
topology does not name, such as the slices of a fibril that is a single residue. `size` must be a positive constant.

```mdscript
slices = chunks(residue(:), 5);
c = contacts(chunks(residue(:), 5), cutoff=5.0);
```

### transpose

<!-- proc name=transpose category=combine -->

```text
transpose(x: T[n][m]) -> T[m][n]
```

Swap the two dimensions of a two-dimensional array (for example the rows and columns of a tensor read with
[`attr`](#attr)). A one-dimensional argument is returned unchanged. The argument must have a fixed size.

```mdscript
per_res = distance(1, 2) in residue(:);
same = transpose(per_res);
```

### Logical operators on selections

<!-- proc name=and aliases=or,xor,not category=combine -->
<a id="and"></a><a id="or"></a><a id="xor"></a><a id="not"></a>

```text
a and b -> atomset      # atoms in both
a or b  -> atomset      # atoms in either
a xor b -> atomset      # atoms in exactly one
not a   -> atomset      # every atom not in a
a and b -> bool | bool[] # element-wise logic on bool values (for example comparisons)
```

Set operations on selections. Arrays of selections are merged before combining.

```mdscript
hydrophobic_c = element("C") and not backbone();
either = residue(1) or residue(2);
outside = not protein();
```

---

## Data from outside the structure

### attr

<!-- proc name=attr category=external -->

```text
attr(path: string) -> float | float[n] | float[n][m]
```

Read a *temporal attribute* of the system as a value that follows the trajectory: the trajectory's own
quantities, and anything loaded next to it from outside. `path` must be a constant string and is either absolute
(`"run/a/edr/potential"`) or relative to a run (`"edr/potential"`), in which case it has to exist in exactly one run.
Attributes that do not vary with time or do not hold numbers are refused. A tensor per frame comes back as
`float[3][3]`; wrap it in `flatten` to get a flat array.

A path with a `*` in it is a *pattern*, and reads every attribute it matches as one array, in the order of their
paths: `attr("edr/coul_sr_*")` is `float[n]` for n energy terms, also when n is 1. A `*` stands for any run of
characters within one segment of the path, so `"edr/*"` is what `edr` holds and not what the groups below it hold.
What it matches has to be in one run, agree on what one value is (all scalars, or all vectors) and on its unit; a frame
axis (`time`) is never one of the matches. The members keep their names: a property that holds the array, or
arithmetic on it that keeps its members (`e * 2`, `e - e0`), is shown member by member under the names their file
gives them (`Coul-SR:Protein-SOL`), or else under the part of the path the `*` stood in. Hovering the pattern in the
editor lists what it matched. Loading another file that matches changes the array, and with it every index into it.

An attribute sampled at other times than the frames - energies written every step, velocities written every fifth
frame - is read at the row whose time matches the frame being evaluated.

Files are not read by the script. They are loaded next to the trajectory (drop them on the window, or load them
from the File menu), and what they hold becomes attributes of the run:

| File | Where its values land | Read with |
|---|---|---|
| Energy file (`.edr`) | `<run>/edr/<term>`, on the file's own time axis | `attr("edr/potential")` |
| xmgrace columns (`.xvg`) | `<run>/xvg/<file>/<legend>`; its own time axis when the x axis is time, one row per frame otherwise | `attr("xvg/energy/coul_sr_prot_sol")` |
| Comma separated (`.csv`) | `<run>/csv/<file>/<column>`; its own time axis when the first column is named time, one row per frame otherwise | `attr("csv/distances/d1_nm")` |

Names are folded to lower case letters, digits and `_`. Units in parentheses in a label or column name, like
`Distance (nm)`, are kept. A file with a time axis must cover every frame of the run it is loaded into.

```mdscript
# Requires an attribute in the loaded system (path is illustrative)
# potential = attr("edr/potential");
# stress = flatten(attr("run/a/edr/pressure_tensor"));
# coulomb = attr("xvg/energy/coul_sr_protein_sol");
# coulomb_terms = attr("edr/coul_sr_*");     # one per energy group pair
# coulomb_lowest = min(attr("edr/coul_sr_*"));
```

`import()` has been removed: the files it read are loaded as above and read with `attr`.

---

## Implicit conversions

The compiler inserts these conversions where an argument of another type is expected:

| From | To |
|---|---|
| `int` | `float`, `irange` |
| `irange` | `frange` |
| `int[]` | `float[]`, `irange[]` |
| `irange[]` | `frange[]` |
| `int[]`, `irange[]` | `atomset` (atom indices) |
| `atomset[]` | `atomset` (merged) |

That is why `atom(1:5)`, `within(2.0:4.0, ...)` and `count(protein())` work without any explicit conversion.

---

## Worked examples

### Flexibility of a peptide

```mdscript
# Distance between the first and last residue centres, and how far each residue drifts from the reference
end_to_end = distance(com(residue(1)), com(residue(15)));
drift = rmsd(backbone());
per_res_size = count(all()) in residue(:);
```

### A local environment

```mdscript
# How many atoms are within 4 Å of residue 3, and the closest approach of residue 1
neighbours = count(within(4.0, residue(3)));
closest = distance_min(residue(1), residue(3));

# The same neighbourhood measured for every residue at once
n_per_res = count(within(4.0) in residue(:));
```

### Planarity of a molecule (per residue)

Four dihedrals are measured in every residue called `PFT` and combined into one number that is zero for a perfectly
flat molecule. Each dihedral is `float[N]` with one value per residue, and arithmetic works element-wise.

```mdscript
d1 = dihedral(22, 20,  1,  2) in resname("PFT");
d2 = dihedral( 2,  3,  6, 10) in resname("PFT");
d3 = dihedral(10,  9, 27, 29) in resname("PFT");
d4 = dihedral(29, 31, 33, 35) in resname("PFT");

# Each term is 0 when the dihedral is ±90°
planarity = abs(abs(d1) - PI / 2) / (PI / 2) +
            abs(abs(d2) - PI / 2) / (PI / 2) +
            abs(abs(d3) - PI / 2) / (PI / 2) +
            abs(abs(d4) - PI / 2) / (PI / 2);
```

### Shape of a molecule

```mdscript
{linear, planar, isotropic} = shape_weights(all());
axis = normal(plane(all()));
```

### Structure of the neighbourhood

```mdscript
# Pair distribution of carbon around oxygen out to 10 Å
g_co = rdf(element("C"), element("O"), 10.0);

# Spatial distribution of oxygen around the residues named PFT
cloud = sdf(resname("PFT"), element("O"), 20.0);

# Density profile of everything along z
profile = density_z(all());
```

### Building on other properties

A property can be used in later statements; the result has the same shape and follows the same rules.

```mdscript
d_a = distance(1, 2) in residue(:);
d_b = distance(1, 3) in residue(:);
ratio = d_a / d_b;
mean_ratio_first_five = (d_a[1] + d_a[2] + d_a[3] + d_a[4] + d_a[5]) / 5.0;
```

---

## Known issues

Defects and gaps in the current implementation. They are listed so that the documentation does not promise more than
the code delivers; remove an entry when it is fixed.

1. **Units are lost in several places.** Unary minus, adding a bare number to a quantity, and the math functions
   return values without a unit tag.
2. **Some values cannot be plotted.** `density(...)`, `coord`, `coord_xy/xz/yz`, `bool` values and arrays of `float[3]`
   are not published as properties (see [What becomes a property](#what-becomes-a-property)).
3. **Constants are not properties.** A script whose statements are all constant (`a = count(protein());`) reports
   "no properties" and cannot be evaluated; add a per-frame quantity or use the value inside another property.
4. **`!=` is not available.** The comparison operators are `<`, `>`, `<=`, `>=` and `==`.

Reports of anything that behaves differently from this reference are very welcome — the fastest fix is a test in
`unittest/test_script.c` and a correction here.

---

## Appendix: reserved words and limits

- Keywords: `in`, `of`, `out`, `and`, `or`, `xor`, `not`.
- Constants: `PI`, `TAU`, `E`.
- Distributions have 1024 bins; volumes are 128 × 128 × 128.
- A script can define as many identifiers as it likes, but each identifier must be defined before it is used, and
  defining the same identifier twice is an error.
- A wildcard `*` stands for any run of characters, in name patterns and in `attr` patterns, where it stays within one
  path segment.
