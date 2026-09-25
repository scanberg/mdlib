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
- Reference: [Selectors](#selectors), [Properties](#properties), [Geometry](#geometry),
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
array. `distribution` and `volume` values can be added, subtracted, multiplied and divided with each other or
with a scalar. Comparisons work on `float`s and arrays of them (a scalar is compared with every element) and give a
`bool` or an array of `bool`. `and`, `or`, `xor`, `not` combine **bitfields** (atom selections) and, separately,
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
| `bitfield` | A set of atoms — the value of every selector. |
| `T[n]` | An array of `n` elements of type `T`. `bitfield[n]` is *n* selections, typically *n* residues or chains. |
| `float[2]`, `float[3]`, `float[4]`, `float[4][4]` | Vectors and matrices (`vec2`, `vec3`, `vec4`, 4×4 matrix). |
| `position` | Not a real type: the argument of a geometric procedure. It accepts atom indices or ranges (`1:3`), bitfields (their centre of mass is used), or `float[3]` vectors. See [Coordinate arguments](#coordinate-arguments). |
| `distribution` | `float[1024]` histogram (for example `rdf`, `density_x`). |
| `volume` | `float[128][128][128]` grid (for example `sdf`). |

Units propagate through arithmetic: `distance(1,2) - distance(1,3)` is in Å, `distance(1,2) * distance(1,3)` in Å²,
`distance(1,2) / distance(1,3)` is dimensionless, and scaling by a plain number (`2.0 * distance(1,2)`) keeps the
unit. Other operations produce a plain number without a unit tag: adding a bare number to a quantity
(`distance(1,2) + 1.0`), unary minus, and the math functions (`abs`, `sqrt`, `min`, ...). Reported units: `distance`
in Å, `angle` and `dihedral` in rad, `density_*` in kg/m³, `sdf` in count/Å³; `PI` and `TAU` carry the unit rad.
`count`, `rdf` and `porosity` are dimensionless. Treat a missing unit tag as "unknown", not as a guarantee that the
value is dimensionless.

### Coordinate arguments

Every procedure documented with `position` arguments (`distance`, `angle`, `dihedral`, `com`, `plane`, `coord`,
`rdf`, `within(radius, around)`, ...) accepts any of:

- **atom indices** such as `1`, `1:3` or `{1, 5, 9}`. Inside a context these are relative to that context;
- **a bitfield** such as `residue(1)` or `element("C")`;
- **a `float[3]` vector**, e.g. `vec3(0, 0, 0)`, or the result of `com(...)`;
- **an array of bitfields**, e.g. `residue(1:3)` or `ring()`.

How a selection is turned into positions depends on the procedure, but the rule of thumb is:

| Argument | `distance`, `angle`, `dihedral`, `com`, `distance_min/max` | `distance_pair`, `coord*` | `rdf`, `within`, `plane`, `shape_weights` |
|---|---|---|---|
| single bitfield / index range | one point: the mass-weighted centre of mass | every atom, individually | every atom, individually (see each entry) |
| array of bitfields (`residue(1:3)`) | see below | one centre of mass **per element** | see each entry |

`distance`, `com` and `within` *merge* an array of bitfields into one selection first, so `com(residue(1:3))` is the
centre of mass of all three residues together and `distance(residue(1:3), vec3(0,0,0))` is a single value. To get
one value per residue, use a [context](#contexts-in-and-out): `distance(1, vec3(0,0,0)) in residue(1:3)`, or use a
procedure that keeps the elements apart such as `distance_pair(residue(1:3), vec3(0,0,0))` (three values) or
`coord_x(residue(1:3))` (three values: the x coordinate of each residue's centre of mass).

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
- a `distribution` (`rdf`, `density_x/y/z`), accumulated over all frames,
- a `volume` (`sdf`), accumulated over all frames.

Constants, plain selections (a `bitfield`), strings and ranges are *not* properties; they exist to be used by other
statements. The engine publishes each property under the attribute `script/<name>`. A script with no property at
all (`a = count(protein());`, which is constant) cannot be evaluated and reports "no properties".

Some values have no property form yet: `float[2]`/`float[3]` arrays (`coord`, `coord_xy`, `shape_weights` without
destructuring) and the combined `density(...)` cannot be plotted directly — use them as inputs to other procedures
or [destructure](#statements-and-comments) them.

---

## Procedure index

| Category | Procedures |
|---|---|
| [Selectors: atoms](#selectors-atom-level) | [`all`](#all), [`atom`](#atom), [`element`](#element), [`name`](#name) (`label`, `type`), [`backbone`](#backbone), [`side`](#side) (`sidechain`), [`ion`](#ion), [`nucleoside`](#nucleoside), [`nucleobase`](#nucleobase), [`ring`](#ring) |
| [Selectors: residues](#selectors-residue-level) | [`protein`](#protein), [`nucleic`](#nucleic) (`nucleotide`), [`water`](#water), [`resname`](#resname) (`residue`, `component`), [`resid`](#resid), [`residue`](#residue), [`component`](#component) |
| [Selectors: instances](#selectors-instance-level) | [`instance`](#instance), [`chain`](#chain), [`chain_id`](#chain_id), [`auth_id`](#auth_id) |
| [Selectors: spatial](#selectors-spatial) | [`within`](#within), [`within_x`](#within_x), [`within_y`](#within_y), [`within_z`](#within_z), [`within_xyz`](#within_xyz) |
| [Properties](#properties) | [`distance`](#distance), [`distance_min`](#distance_min), [`distance_max`](#distance_max), [`distance_pair`](#distance_pair), [`angle`](#angle), [`dihedral`](#dihedral), [`rmsd`](#rmsd), [`rdf`](#rdf), [`density_x`](#density_x) / [`density_y`](#density_y) / [`density_z`](#density_z), [`density`](#density), [`sdf`](#sdf), [`count`](#count), [`contact_count`](#contact_count), [`contacts`](#contacts), [`degree`](#degree), [`porosity`](#porosity) |
| [Geometry](#geometry) | [`com`](#com), [`plane`](#plane), [`shape_weights`](#shape_weights), [`coord`](#coord), [`coord_x`](#coord_x) / [`coord_y`](#coord_y) / [`coord_z`](#coord_z), [`coord_xy`](#coord_xy) / [`coord_xz`](#coord_xz) / [`coord_yz`](#coord_yz) |
| [Math](#math) | [`sqrt`](#sqrt), [`cbrt`](#cbrt), [`abs`](#abs), [`floor`](#floor), [`ceil`](#ceil), [`sin`](#sin), [`cos`](#cos), [`asin`](#asin), [`acos`](#acos), [`atan`](#atan), [`atan2`](#atan2), [`log`](#log), [`log2`](#log2), [`log10`](#log10), [`exp`](#exp), [`exp2`](#exp2), [`pow`](#pow), [`min`](#min), [`max`](#max) |
| [Linear algebra](#linear-algebra-and-constructors) | [`vec2`](#vec2), [`vec3`](#vec3), [`vec4`](#vec4), [`dot`](#dot), [`cross`](#cross), [`length`](#length), [`normalize`](#normalize), [`mul`](#mul) |
| [Combining](#combining-and-reshaping) | [`join`](#join) (`flatten`), [`split`](#split), [`chunks`](#chunks), [`transpose`](#transpose), and the operators `and` `or` `xor` `not` on selections |
| [External data](#data-from-outside-the-structure) | [`attr`](#attr) |

---

## Selectors

A **selector** returns a `bitfield` (a set of atoms) or a `bitfield[]` (several sets — one per residue, chain, ring,
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

These return a single `bitfield`.

### all

<!-- proc name=all category=selector.atom -->

```text
all() -> bitfield
```

Every atom in the system (or, inside a context, every atom of that context).

```mdscript
n_atoms = count(all());
per_residue = count(all()) in residue(:);      # atoms per residue
```

### atom

<!-- proc name=atom category=selector.atom -->

```text
atom(indices: irange[]) -> bitfield
atom(indices: int[])    -> bitfield
```

Select atoms by 1-based index. Indices are relative to the current context, so inside `in residue(:)` `atom(1)` is
the first atom of each residue. An index outside the valid range is a compile error. Note that `atom` takes one
argument: for several atoms write an array, `atom({1, 5, 9})`, or a range, `atom(1:9)`.

```mdscript
first_ten = atom(1:10);
every_other = atom(1:2:20);                    # 1,3,5,...,19
n_terminal = atom(1:4) in residue(1);          # first four atoms of residue 1
```

### element

<!-- proc name=element category=selector.atom -->

```text
element(symbol: string[]) -> bitfield
element(z: irange[])      -> bitfield
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
name(pattern: string[])  -> bitfield
label(pattern: string[]) -> bitfield
type(pattern: string[])  -> bitfield
```

Select atoms by atom name (for example `"CA"`, `"OW"`). `label` and `type` are aliases. Patterns may contain the
wildcard `*` at the beginning, the end, or both ends (`"C*"`, `"*1"`, `"*H*"`); at most two `*` are allowed.

```mdscript
alpha_carbons = name("CA");
any_carbon = name("C*");
backbone_atoms = name({"N", "CA", "C", "O"});
```

### backbone

<!-- proc name=backbone category=selector.atom -->

```text
backbone() -> bitfield
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
side()      -> bitfield
sidechain() -> bitfield
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
ion() -> bitfield
```

Atoms of residues that were recognised as ions (for example Na⁺, Cl⁻).

```mdscript
ions = ion();
solvent_and_ions = water() or ion();
```

### nucleoside

<!-- proc name=nucleoside category=selector.atom -->

```text
nucleoside() -> bitfield
```

Atoms flagged as part of the nucleoside (sugar) portion of nucleic-acid residues.

```mdscript
sugars = nucleoside();
```

### nucleobase

<!-- proc name=nucleobase category=selector.atom -->

```text
nucleobase() -> bitfield
```

Atoms flagged as part of the nucleobase of nucleic-acid residues.

```mdscript
bases = nucleobase();
```

### ring

<!-- proc name=ring category=selector.atom -->

```text
ring() -> bitfield[]
```

One selection per ring found in the molecular graph, useful to measure something per ring. Inside a context only
rings that lie completely within the context are returned.

```mdscript
n_ring_atoms = count(ring());                   # the rings are merged into one selection
d = distance_pair(ring(), vec3(0, 0, 0));       # centre of each ring to the origin, one value per ring
```

### Selectors: residue level

These return a `bitfield[]` with one selection per matching residue. Use them as the right-hand side of `in`, or
pass them where a bitfield is accepted (an array of bitfields is then flattened into one selection; [`rmsd`](#rmsd)
is the exception and gives one value per selection).

### protein

<!-- proc name=protein category=selector.residue -->

```text
protein() -> bitfield[]
```

Every residue that is part of a protein.

```mdscript
n_protein_atoms = count(protein());                 # the residues are flattened into one selection
d = distance(1, 2) in protein();                    # one value per protein residue
```

### nucleic

<!-- proc name=nucleic aliases=nucleotide category=selector.residue -->
<a id="nucleotide"></a>

```text
nucleic()     -> bitfield[]
nucleotide()  -> bitfield[]
```

Every residue that is part of a nucleic acid. `nucleotide` is an alias.

```mdscript
dna = nucleic();
n_dna_atoms = count(nucleic());
```

### water

<!-- proc name=water category=selector.residue -->

```text
water() -> bitfield[]
```

Every water molecule.

```mdscript
solvent = water();
n_waters = count(water(), "residue");
```

### resname

<!-- proc name=resname aliases=residue,component category=selector.residue -->

```text
resname(pattern: string[])   -> bitfield[]
residue(pattern: string[])   -> bitfield[]
component(pattern: string[]) -> bitfield[]
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
resid(ids: irange[]) -> bitfield[]
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
residue(pattern: string[])   -> bitfield[]      # same as resname
residue(index: irange[])     -> bitfield[]      # by position
residue(sel: bitfield[])     -> bitfield[]      # grow to whole residues
```

- `residue("ALA")` — select residues by name.
- `residue(1:5)` — select residues by their **1-based position** in the system, `residue(:)` selects all.
- `residue(sel)` — every residue that the selection touches, expanded to its whole residue.

```mdscript
first = residue(1);
first_five = residue(1:5);
all_res = residue(:);
ca_residues = residue(name("CA"));           # whole residues that contain a CA atom
near = residue(within(4.0, residue(1)));     # whole residues within 4 Å of residue 1
```

### component

<!-- proc name=component category=selector.residue -->

```text
component(pattern: string[]) -> bitfield[]
component(index: irange[])   -> bitfield[]
```

Alias of [`residue`](#residue) for the `string` and `irange` overloads. Coarse-grained and non-protein systems
often call their building blocks *components*.

```mdscript
first_component = component(1);
named = component("PFT");
```

### Selectors: instance level

Instances are the largest building blocks in the file: chains in a PDB file, molecules in a GRO file.

### instance

<!-- proc name=instance category=selector.instance -->

```text
instance(id: string[])     -> bitfield[]
instance(index: irange[])  -> bitfield[]
```

Select instances by identifier or 1-based position.

```mdscript
first_instance = instance(1);
by_id = instance("A");
```

### chain

<!-- proc name=chain category=selector.instance -->

```text
chain(id: string[])     -> bitfield[]
chain(index: irange[])  -> bitfield[]
chain(sel: bitfield[])  -> bitfield[]      # grow to whole chains
```

Select chains by identifier (`"A"`), by 1-based position (`chain(1:2)`), or expand a selection to every chain it
touches.

```mdscript
chain_a = chain("A");
first_two = chain(1:2);
its_chain = chain(atom(1));                    # the chain that contains atom 1
```

### chain_id

<!-- proc name=chain_id category=selector.instance -->

```text
chain_id(id: string[]) -> bitfield[]
```

Select chains by their identifier (the label in the file). Equivalent to `chain("...")`.

```mdscript
a = chain_id("A");
```

### auth_id

<!-- proc name=auth_id category=selector.instance -->

```text
auth_id(id: string[]) -> bitfield[]
```

Select instances by their author-assigned identifier (as found in mmCIF files), which can differ from the label
identifier used by `instance` and `chain_id`. A string that matches no instance is a compile error.

```mdscript
by_author = auth_id("A");
```

### Selectors: spatial

Spatial selectors depend on atom positions, so they are re-evaluated for **every frame**; their size is not known
in advance. Consequently their results can be counted or measured, but not used as a context whose length must be
known up front.

### within

<!-- proc name=within category=selector.spatial -->

```text
within(radius: float,  around: position[]) -> bitfield
within(radius: frange, around: position[]) -> bitfield
within(radius: float)                       -> bitfield      # centre from the context
within(radius: frange)                      -> bitfield
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
within_x(range: frange) -> bitfield
```

Atoms whose x coordinate lies in `range` (Å, in the coordinate system of the trajectory).

```mdscript
slab = within_x(20:25);
left_half = within_x(:30.0);
```

### within_y

<!-- proc name=within_y category=selector.spatial -->

```text
within_y(range: frange) -> bitfield
```

Atoms whose y coordinate lies in `range`.

```mdscript
slab = within_y(10:15);
```

### within_z

<!-- proc name=within_z category=selector.spatial -->

```text
within_z(range: frange) -> bitfield
```

Atoms whose z coordinate lies in `range`.

```mdscript
slab = within_z(20:25);
n_in_slab = count(within_z(20:25));
```

### within_xyz

<!-- proc name=within_xyz category=selector.spatial -->

```text
within_xyz(x: frange, y: frange, z: frange) -> bitfield
```

**Parameters:** `x`, `y`, `z` (all required).

Atoms inside an axis-aligned box.

```mdscript
box = within_xyz(20:25, 20:25, 20:25);
n_in_box = count(within_xyz(20:25, 20:25, 20:25));
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
distance_min(a: position[], b: position[]) -> float [Å]
```

**Parameters:** `a`, `b`.

The shortest distance between any atom of `a` and any atom of `b` (closest approach between two groups). Unlike
`distance`, the result is a single value even when the selections change size from frame to frame.

```mdscript
closest = distance_min(residue(1), residue(2));
gap = distance_min(within_x(:20), within_x(30:));
```

### distance_max

<!-- proc name=distance_max category=property -->

```text
distance_max(a: position[], b: position[]) -> float [Å]
```

**Parameters:** `a`, `b`.

The largest distance between any atom of `a` and any atom of `b` (the extent between two groups).

```mdscript
extent = distance_max(residue(1), residue(2));
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
```

**Parameters:** `a`, `b`, `c` (all required).

The angle at `b` between the directions to `a` and `c`, in `[0, π]`.

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
rmsd(sel: bitfield[]) -> float[] [Å]
```

The mass-weighted root-mean-square deviation of the selected atoms from the system's reference structure (the first
frame), after removing translation and finding the optimal rotation. Each bitfield of the argument is a structure of
its own, fitted and reported separately: the result has one value per bitfield. Unlike most procedures, which merge
an array of selections into one, `rmsd` keeps them apart; use [`flatten`](#flatten) to fit them together as one
structure. A single selection (`backbone()`, `within(...)`, a combination with `and`/`or`) gives a single value.

The selection may change from frame to frame, and an empty one gives 0. If the *number* of selections changes from
frame to frame (`residue(within(...))`) the result is not published as a property; `flatten` it. Inside a context
only the atoms of the context are used. Periodic boundaries are handled as long as a structure spans less than half
the cell.

```mdscript
backbone_rmsd = rmsd(backbone());               # one value
chain_rmsd = rmsd(chain(:));                    # one value per chain, each fitted on its own
residue_rmsd = rmsd(protein());                 # one value per protein residue
protein_rmsd = rmsd(flatten(protein()));        # the whole protein as one structure
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
density_x(sel: bitfield[]) -> distribution [kg/m³]
```

Mass density of the selection along the x axis of the simulation cell, as a 1024-bin histogram accumulated over
all frames. The selection is merged, so an array of selections is treated as one. Requires a unit cell.

```mdscript
rho_x = density_x(protein());
```

### density_y

<!-- proc name=density_y category=property -->

```text
density_y(sel: bitfield[]) -> distribution [kg/m³]
```

As [`density_x`](#density_x), along the y axis.

```mdscript
rho_y = density_y(water());
```

### density_z

<!-- proc name=density_z category=property -->

```text
density_z(sel: bitfield[]) -> distribution [kg/m³]
```

As [`density_x`](#density_x), along the z axis. The classic use is the density profile across a membrane or slab.

```mdscript
rho_z = density_z(all());
rho_water = density_z(water());
```

### density

<!-- proc name=density category=property -->

```text
density(sel: bitfield[]) -> float[3][2][1024]
```

The three axis profiles of the selection in a single value. It is not published as a property; use
[`density_x`](#density_x), [`density_y`](#density_y) and [`density_z`](#density_z) instead.

```mdscript
rho = density(all());
```

### sdf

<!-- proc name=sdf category=property -->

```text
sdf(structures: bitfield[], target: bitfield, extent: float) -> volume [count/Å³]
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
count(sel: bitfield[])                  -> float
count(sel: bitfield[], unit: string)    -> float
count(c: contact)                       -> float
count(c: contact, unit: string)         -> float
```

**Parameters:** `sel` or `c` (required), `unit` (optional).

The number of atoms in the selection (an array is merged first). With `unit` the number of larger structures the
selection touches instead: `"atom"`, `"residue"`, `"chain"` or `"structure"`. An unknown unit is a compile error that
lists the valid ones. A `count` of a dynamic selection changes over time and is a temporal property.

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
contact_count(a: bitfield[], b: bitfield[], cutoff: float) -> float[N]
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
contacts(a: bitfield[], b: bitfield[], cutoff: float, exclude_bonds: int, min_separation: int,
         parent: bitfield[], exclude_within: bitfield[]) -> contact
```

**Parameters:** `a` and `cutoff` (required), `b`, `exclude_bonds` (default 3), `min_separation` (default 0),
`parent` and `exclude_within` (optional). Everything after `a` and `b` is given by name.

The groups in contact at each frame: pairs of groups with at least one pair of particles closer than `cutoff` �.
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
porosity(sel: bitfield) -> float
```

The void ratio of the selection inside its bounding box: the fraction of the box not occupied by atoms, in `[0, 1]`.

```mdscript
void_fraction = porosity(all());
```

---

## Geometry

Procedures that turn atoms into points, vectors and planes. Their results are values you feed into other
procedures (`com` into `distance`, `plane` into `dot`) rather than properties to plot; those that result in a plain
`float`, `float[3]` or `float[4]` are published as properties when they change over the trajectory.

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
plane(a: position[]) -> float[4]
```

The best-fit plane through at least three positions, returned as `(nx, ny, nz, d)` where `n` is the unit normal (the
direction of least variance) and `d = n · c` with `c` the centre of the points, so points on the plane satisfy
`n · x = d`.

```mdscript
p = plane(resname("PFT"));                          # one plane through all PFT atoms
planes = plane(1:4) in resname("PFT");              # one plane per residue
```

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

The smaller of two numbers, or the smallest element of an array (with one argument).

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

The larger of two numbers, or the largest element of an array (with one argument). The volume overloads are
useful to clamp an `sdf`.

```mdscript
farthest = max(distance(1, 2), distance(1, 3));
largest = max(distance(1, 2) in residue(:));
floor_at_zero = max(dihedral(1, 2, 3, 4), 0.0);
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
length(v: float[]) -> float
```

The Euclidean norm of a vector.

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

### join

<!-- proc name=join aliases=flatten category=combine -->
<a id="flatten"></a>

```text
join(sel: bitfield[])    -> bitfield
flatten(sel: bitfield[]) -> bitfield
flatten(x: T[..])        -> T[]          # any array, one-dimensional result
```

Merge an array of selections into one. `flatten` is an alias for selections and additionally flattens any
multi-dimensional array of numbers into one dimension, for example the per-residue results of a context.

```mdscript
all_protein = join(protein());
merged = flatten(residue(1:3));
per_res = flatten(distance(1, 2) in residue(:));
```

### split

<!-- proc name=split category=combine -->

```text
split(sel: bitfield, parts: int) -> bitfield[]
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
chunks(sel: bitfield[], size: int) -> bitfield[]
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
a and b -> bitfield      # atoms in both
a or b  -> bitfield      # atoms in either
a xor b -> bitfield      # atoms in exactly one
not a   -> bitfield      # every atom not in a
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
| `int[]`, `irange[]` | `bitfield` (atom indices) |
| `bitfield[]` | `bitfield` (merged) |

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
axis = plane(all());
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
- Wildcards (`*`) in name queries are limited to one leading and/or one trailing `*`.
