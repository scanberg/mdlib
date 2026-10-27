# VIAL language reference

Draft for phase 0. Every example ends with its expected result, and the phase 1 test suite runs them all.

## How to read this reference

Examples run on two structures from mdlib's test data, and each example line ends with the result the test suite expects.

| Name | File | What it holds |
| --- | --- | --- |
| T | `test_data/tubulin-A-B.pdb` (PDB 4EB6) | Tubulin with vinblastine, X-ray, no hydrogens. Chain A: 433 amino acids numbered 1–440 with 39–45 missing; GTP 501, MG 502, SO4 503 and 504, waters 601–605. Chain B: 432 amino acids numbered 1–442 with 45–46 and 361–368 missing; GDP 501, SO4 502, waters 601–603. Chain C: vinblastine, VLB 503. Some residues have alternate locations. 17 parts: two peptide strands, six ligands, one ion and eight waters, each part its own molecule |
| P | `test_data/1ALA-560ns.pdb` | A 15-residue alanine peptide with hydrogens: 153 particles, 38 frames. The file gives no chain ids, so P has no chains; its strand is one part and one molecule |

A comment above each example names its structure: `# on T` or `# on P`.

| Comment | Means |
| --- | --- |
| `# -> float[433][3] Å` | the type: dtype, shape, unit |
| `# -> 442` | the value |
| `# -> float Å, frame` | type and stage; where no stage is written, it is not checked |
| `# -> float Å, trajectory, per frame` | a trajectory value that keeps the frame |
| `# -> float[n] Å` | a letter as an extent is not checked |
| `# error E0301` | compiling fails with this code |
| `# warning W0501 -> float ?` | it compiles, with this warning and this type |
| a warning the comment does not list | fails the test, so every warning an example gives is written down |
| `a, b = … # -> set[424] and set[424]` | tuple unpacking: one type per name, in order |
| pass 2 | the pass, checked where it is written |
| `# warning W0201` with no type | only the warning is checked |
| ; R0101 in frames where … | a runtime code, checked on the frames the comment names |
| a code block | one script, unless its first comment says `with the names above`, which continues the block before |
| text after a colon, or after a comma that follows the type | commentary, not checked |
| a line with no expectation | an illustration, not checked |
| `# on T`, `# on P` | the structure from here on; `# on a run with an energy file` has no test file yet |

The counts come from the files, as mdlib loads them: it keeps only the first alternate location, so T's protein has 6,784 particles rather than 6,799 records. They must be confirmed against mdlib's classification when the tests are written. T has a CRYST1 cell of 64.7 × 129.7 × 252.8 Å, so, like any system with a unit cell, it is periodic and distances use the minimum image. The examples on T stay well under half the box; a test of long distances turns periodicity off for the system, which a user can do too.

## Lexical structure and statements

A script is a Python 3 module: UTF-8 text, `#` comments, indentation for blocks, and Python's rules for names and literals. The only addition is the inclusive range `a..b`.

### Literals and operators

| Form | Examples | Notes |
| --- | --- | --- |
| Integers | `42`, `1_000_000` | 32-bit |
| Floats | `0.5`, `1e-3` | Stored in 32 bits |
| Strings | `'CA'`, `"GTP"` | Names and labels; no f-strings |
| Booleans | `True`, `False` |  |
| None | `None` | The absent default, as in `between=None` |
| Identity | `x is None`, `x is not None` | Compile-time only: `is` compares with `None` or a static value |
| Arrays | `[45, 48, 112]`, `[5.0, 0.0, 2.3]` | Every item has one type |
| Tuples | `(0.0, 20.0)` | For pair arguments and several return values |
| Inclusive ranges | `10..20` | Filter values and `in` tests. Binds looser than arithmetic and tighter than comparisons, so `1..n+1` is `1..(n+1)` |
| Arithmetic | `+ - * / // % **` | Elementwise on arrays |
| Matrix product | `a @ b` | Binds as `*` does, as in Python |
| Comparisons | `< <= > >= == !=`, chained | Elementwise on arrays |
| Membership | `x in 10..20`, `x not in names` |  |
| Logic | `and`, `or`, `not` | Scalars as in Python, stopping early. Elementwise on bool arrays; intersection, union and complement on sets |
| Bits | `& \| ^ ~ << >>` | Integers only |
| Other | `a if c else b`, `x[i]`, `f(x, k=v)`, `x.f(...)`, `x.p` | Methods and properties are another spelling of a call |

- A number followed by `..` ends there, so `1..4` is a range, not `1.` and `.4`. Ranges of floats, such as `0.0..5.0`, are allowed in `in` tests.
- `int` with `float` gives `float`. `/` always gives `float`; `//` keeps `int`.
- `sum` of a bool array is an `int`.

### Statements

| Statement | Available | Notes |
| --- | --- | --- |
| Assignment, tuple unpacking | First version | A top-level name is bound once (E0305) |
| `def`, `return`, defaults, keyword arguments, docstrings | First version | Specialized per call. No recursion (E0802). A def may be nested in another and read its names; it cannot be returned or stored, but may be passed as an argument, to fold or to another def, and each def passed gives its own specialization. Parameters may be annotated with a kind, which maps a call over arrays; see Mapping a def over arrays |
| `if`, `elif`, `else` | First version | At top level the condition must be static (E0403); inside functions it may change per frame. Code that a static condition does not reach, the untaken branch of an if or a ternary and what follows a return taken on a static condition, is parsed and its names checked, but not type-checked. A condition on an array is E0303, which points to `where()`. A ternary `a if c else b` on a frame-stage scalar is allowed anywhere and lowers to `where()` |
| `for`, `break`, `continue` | First version | Over a range, an array or a set; over a set, each particle as a one-particle set. A top-level loop variable belongs to the loop: it is not published, and it does not count for E0305 |
| List comprehensions | First version | They build an array whose first dimension is what they iterate over; a nested comprehension adds a dimension per level. An `if` clause gives a dynamic extent |
| `pass` | First version |  |
| `assert c, 'message'` | First version | c must be static; when it is false, compiling fails with E0314 and the message, at the user's call. A condition that changes per frame is E0403. The prelude uses it to refuse unknown options |
| Augmented assignment `x += 1` | First version, inside functions | At top level it rebinds a name: E0305 |
| `while` | Later | Parsed from the start; until then E0801 |
| Item assignment `a[i] = x` | Later, inside functions | Parsed from the start; until then E0801 |
| Classes, `import`, `lambda`, exceptions, generators, `global`, decorators, `with`, `async`, dicts, set literals, `*args`, \*\*kwargs, `print` | Never | E0803 names the construct |

### Examples

```python
# on T
n    = len(residues(chain='A'))                 # -> 442: 433 amino acids, GTP, MG, two SO4, five waters
a    = residues(36..48, chain='A')              # -> set[6]: 39 to 45 are missing
ok   = 2 < n <= 500                             # -> True
big  = 1_000_000                                # -> 1000000
ids  = [45, 48, 112]                            # -> int[3]
n += 1                                          # error E0305

def centre_distance(a, b):
    """Distance between the centres of mass of a and b."""
    return distance(com(a), com(b))

d    = centre_distance(particles(resname='GTP'), particles(resname='GDP'))   # -> float Å

def halve(x):
    while x > 1.0:                              # error E0801: until while is supported
        x = x / 2.0
    return x
```

## Types and how they print

A type is a dtype, a shape and a unit, printed C-style: `float[433][3] Å`. Shapes, indexing, broadcasting and reductions follow NumPy. The frame is not a dimension: a frame value is an array that exists once per frame.

| Printed | Means |
| --- | --- |
| `float` | a scalar with no unit |
| `float Å` | a scalar length |
| `float[433] Å` | 433 lengths |
| `float[433][3] Å` | 433 positions; the last dimension runs fastest |
| `float[?][3]` | the first extent is known only per frame |
| `set`, `set[433]`, `set[14][4]` | one set; 433 sets; 14 matches of four particles each |
| `float ?` | a unit that is unknown, or was dropped after a mismatch |
| `float Å \| missing` | lengths that may be missing in some frames or entries; they carry a validity bitmask |
| `float[4..14][3] Å` | inside a mapped call: an extent that changes from call to call, shown with its range |
| `float64 ps` | `time()`, the one 64-bit value |

### Dtypes

| Dtype | Holds |
| --- | --- |
| `bool` | true or false |
| `int` | 32-bit integers |
| `float` | 32-bit values |
| `str` | names and labels; never published |
| `set` | a set of particles |
| `range` | an inclusive range of values |
| `drawing` | lines, shapes and labels for hover; never published |

`int` with `float` gives `float`. A number literal takes the dtype of the array it meets when both are integers or both are floats, as NumPy 2 does, so `0.5 * a` keeps a float array's dtype, and 0.5 times an int array is float. Reductions, including those over frames, accumulate in 64 bits internally and round once at the end, and an integer result that does not fit in 32 bits is R0102. `time()` is the one 64-bit value, since a 32-bit float resolves a 10 µs run to about 1 ps; what is computed from it alone stays 64-bit and rounds to 32 bits where it meets other data.

### Dimensions

- Broadcasting is NumPy's: shapes align from the right, a dimension of 1 stretches, and equal extents pair in order. Anything else is E0301.
- `axis=` takes an int, negative from the end, or a tuple. `.T` reverses the dimensions, and `a[:, None]` adds one.
- Reductions reduce every dimension by default, as in NumPy. The frame is the exception: it is only reduced when named, `axis='frame'`.
- Positions and other 3-vectors are a last dimension of 3. `c.x`, `c.y` and `c.z` are `c[..., 0]`, `c[..., 1]` and `c[..., 2]`.
- Every dimension can carry a description of what it runs over, such as the residues of chain A, a pattern's matches, or x, y and z. Descriptions give labels in plots and exports, show in hover, and are the identity that `rmsf[r]`, `by=`, `spread`, `fit(on=)` and `common` need; those refuse an array without one (E0311). Descriptions never refuse arithmetic. Two dimensions of equal extent that describe different particles or units, paired by position, give W0605 and run, unless both derive from slices of one array or set, through elementwise operations, positions(), fit() and calls, such as c\[1:\] - c\[:-1\] or positions(ca\[0:7\]) against positions(ca\[8:15\]); or both come from one common() call; or a zip in a comprehension pairs them.
- A result's labels come from the inputs that have each dimension; an extent of 1 or a missing dimension does not count. One such input gives its labels, inputs that describe the same units keep them, and inputs that describe different units give pair labels such as `MET1/A–GLU3/A`. Hover says what was paired with what, for example “dim 0: residues 1–4 of chain A with residues 1–4 of chain B”.
- A value has at most six dimensions (E0308).
- A dynamic extent `?` has an identity, the selection that produced it. Elementwise operations, `.com`, `positions()` indexing by a set, particles(x), per-particle properties such as mass(), and shape() keep it. A slice by fixed offsets, such as \[1:\] or \[:-1\], gives a related identity, ?−1, which pairs with any slice of the same ? and length, so c\[1:\] - c\[:-1\] and diff() work on a dynamic extent. Other slicing, a comprehension with an `if`, and a new selector give a new one. Two different `?` never broadcast (E0301).

```text
c : float[433][3] Å                                    frame
    dim 0   433   residues of chain A: MET1/A … VAL440/A   (ids 1–440; 39–45 missing)
    dim 1     3   x, y, z
    float32, row-major, 5.2 KB per frame; stored as float[1][433][3], the frame outermost
```

### Units

- `+`, `-` and comparisons convert compatible units, such as nm and Å. Incompatible units give W0501, and the result's unit becomes unknown, so later steps do not warn again.
- A bare number, and an array built by zeros, ones or full without unit=, takes the unit the operation expects. That is the other operand's unit in `+`, `-`, comparisons, `min`, `max`, `where` and `clip`, and the parameter's declared unit elsewhere, such as Å for a point, the radius of `within` and the `range` of a histogram. In `*` and `/` it has no unit.
- Trigonometric functions take angles in any angle unit. `exp`, `log` and the like take values without a unit; anything else gives W0501.
- `unitless(x)` drops the unit on purpose, without a warning.

```python
# on T
res_a = residues(chain='A', kind='amino_acid')           # -> set[433]
vlb   = particles(resname='VLB')                         # -> set
c     = com(res_a)                                       # -> float[433][3] Å
cv    = vlb.com                                          # -> float[3] Å
d     = distance(c, cv)                                  # -> float[433] Å
o     = distance(c, [0.0, 0.0, 0.0])                     # -> float[433] Å
z     = c.z                                              # -> float[433] Å
mid   = mean(c, axis=0)                                  # -> float[3] Å
all_  = mean(c)                                          # -> float Å: every dimension, as in NumPy
dm    = distance(c[:, None], c[None, :])                 # -> float[433][433] Å: a distance matrix by broadcasting
sym   = dm - dm.T                                        # -> float[433][433] Å
step  = c[1:] - c[:-1]                                   # -> float[432][3] Å
bad   = d + c                                            # error E0301: float[433] and float[433][3] do not broadcast
sc    = d[:, None] * c                                   # -> float[433][3] Å²
g     = distance(residues(501, chain='A').com, cv)       # -> float[1] Å
nope  = distance(res_a, vlb)                             # error E0302: residues of several particles are not points

gtp   = particles(resname='GTP').com                     # -> float[3] Å
gdp   = particles(resname='GDP').com                     # -> float[3] Å
dg    = distance(gtp, gdp)                               # -> float Å
near  = dg < 12.0                                        # -> bool: 12.0 is read as Å
scaled = dg / 10.0                                       # -> float Å
mixed = dg + angle(gtp, gdp, cv)                         # warning W0501 -> float ?
sw    = 1 / (1 + exp((dg - 4.0) / 0.5))                  # warning W0501 -> float ?: exp of a length
sw2   = 1 / (1 + exp(unitless(dg - 4.0) / 0.5))          # -> float
```

## Stages

Every value has one of three stages, inferred from what it reads. Nobody declares a stage; hover and a mark in the gutter show it.

| Stage | Reads | Computed | Typical sources |
| --- | --- | --- | --- |
| static | the topology and constants | once, when compiling | selectors, literals, `smiles()`, `at(x, frame=k)`, values pinned from another system |
| frame | one frame | once per frame | `positions()`, `attr()`, `within()`, anything that uses them |
| trajectory | many frames | over the frames it needs | reductions with `axis='frame'`, `histogram` over frames, `smooth`, `autocorrelation` |

- A value's stage is the latest stage among its inputs and its operation. A trajectory value that keeps the frame axis, as `smooth()` gives, is marked per frame.
- A frame value carries the frame axis, and every operation maps over it unless the operation names it.
- A value's pass is the read over the trajectory in which it is complete. Static values are pass 0, and frame values are pass 1 at the earliest. A frame value, a fold, and a reduction that streams over frames, such as `mean` or `histogram` with `axis='frame'`, are complete in the same pass as their latest input. A bounded window, as in `smooth` and `derivative`, also stays in its input's pass: what reads it waits until the window has arrived, as long as frames arrive in order. A value that reads a finished whole-trajectory result of pass k, such as `ends - mean_e`, is pass k + 1, since that result is known only when pass k ends. An operation on whole-trajectory results with no frame axis, such as `eigh(C)`, keeps pass k. So `ends` and `mean_e` are pass 1; `dev = ends - mean_e` and `mean(dev, axis='frame')` are pass 2; `ends - previous(ends)` and `ends - smooth(ends, 5)` are pass 1. Hover and the IR dump show the pass, and more than three passes is W0406.
- A published frame value must have the same shape in every frame (E0401). Selections are the exception: a frame-stage set or array of sets is published as a recipe and evaluated for the frame being shown.
- A value from another system must be pinned to one of its frames with `at()` (E0402).
- `at(x, frame=k)` is evaluated when compiling by reading frame k, so compiling waits for that frame. `x` may depend on static values and on frame k alone: a trajectory value or another `at()` inside it is E0404, and a k past the run's last frame is E0405.
- A value can be missing in a frame: an empty set where a point is needed, a group with no values, a position outside a dynamic extent (R0101), an integer out of range (R0102), or a set empty in that frame under `between=`. Evaluation goes on. A value that may be missing has an optional type, shown in hover as `float[433] | missing`, and carries a validity bitmask whatever its dtype, float included; a value the compiler knows is always present carries none. `missing_where(c, x)` marks values missing, `valid(x)` gives the mask as a bool array, and `filled(x, v)` replaces missing values and drops the mask. Elementwise operations are missing where an input is; reductions skip missing entries and frames and count them, published series show gaps, and W0405 is given at compile time, once per statement, where a value may be missing; values computed from it do not warn again. By default a fold skips a frame whose input is missing: `step` is not called and the state carries over. NaN from arithmetic is a value, not missing, and reductions do not skip it, so a bug shows as NaN.
- `fold(step, x, init)` keeps one state from frame to frame, in order. With skip=False, step is also called for a frame whose x is missing, and the state may hold missing values. An init that is a frame value is taken from the first frame where it is present. Folds run once over the whole run and their per-frame results are stored, so evaluating a subset of frames reads them rather than running the fold again. Its result is a trajectory value per frame, and it is the only operation that needs frames in order.
- Everything evaluated within a frame-stage scope is frame stage.
- Only published and hovered values are computed.

```python
# on P
ca     = particles(name='CA')                          # -> set, static
ends   = distance(ca[0], ca[14])                       # -> float Å, frame: single particles are points
rg     = gyration(particles)                           # -> float Å, frame
mean_e = mean(ends, axis='frame')                      # -> float Å, trajectory
sm     = smooth(ends, 5)                               # -> float Å, trajectory, per frame
hist   = histogram(ends, axis='frame', bins=32, range=(0.0, 50.0))   # -> float[32], trajectory
first  = at(ends, frame=0)                             # -> float Å, static: read from frame 0 when compiling
drift  = ends - first                                  # -> float Å, frame
dev    = ends - mean_e                                 # -> float Å, trajectory, per frame, pass 2
near   = residues(within(4.0, residues[0]))            # -> set[?], frame: published as a recipe, evaluated per frame
n_near = len(near)                                     # -> int, frame
d_near = distance(near.com, ca[0])                     # error E0401: a value, not a selection, and its length changes
_d     = distance(near.com, ca[0])                     # -> float[?] Å, frame
```

## Selectors, keywords and scope

A system is built of particles, residues, parts and molecules, each nested in the next. A part is one polymer strand, one glycan tree, one ligand, one lipid, one ion or one water; a strand goes on over a gap and ends where the chemistry says it does. A molecule is what the bonds hold together, and with a topology it is the topology's molecule, so a glycoprotein is one molecule of several parts. A chain is not a level of it: it is the author's chain id, a label that the data files particles of any level under, and it is never inferred. Where the data gives no chain ids, as in a GRO file or a PDB file written by a simulation, the system has no chains, and a strand is found as a part.

`particles()` returns one set; the grouping selectors return one set per unit. A selector written without parentheses means all its units, and every selector searches its scope.

| Selector | Returns | One unit is |
| --- | --- | --- |
| `particles(...)` | `set` | an atom, a bead or a virtual site; given sets or arrays of sets, it merges them |
| `atoms(...)`, `beads(...)` | `set` | `particles(kind='atom')`, `particles(kind='bead')` |
| `residues(...)` | `set[n]`, one per residue | a residue |
| `parts(...)` | `set[n]`, one per part | a polymer strand, glycan tree, ligand, lipid, ion or water |
| `molecules(...)` | `set[n]`, one per molecule | what the bonds hold together, or a topology's molecule |
| `chains(...)` | `set[n]`, one per chain | everything filed under one author chain id; E0207 where the data gives none |

### Keywords

| Selector | Keyword | Matches | Value |
| --- | --- | --- | --- |
| `particles` | `name` | the atom or bead name | string with `*` and `?`, or a list |
|  | `serial` | the serial number in the file | int, `a..b` or a list |
|  | `element` | the element | string or list |
|  | `kind` | atom, bead or virtual | string or list |
|  | `resid` | the number of its residue in the file | int, `a..b` or a list |
|  | `resname` | the name of its residue | string with `*` and `?`, or a list |
|  | `chain` | the author id of its chain | string or list |
| `residues` | `id` | the sequence number in the file | int, `a..b` or a list |
|  | `name` | the residue name | string with `*` and `?`, or a list |
|  | `kind` | amino\_acid, nucleotide, water, ion, lipid, other | string or list |
|  | `chain` | the author id of its chain | string or list |
| `parts` | `kind` | polymer, peptide, nucleic, dna, rna, glycan, ligand, lipid, ion, water; see below | string or list |
|  | `chain` | the author id of its chain | string or list |
| `molecules` | `name` | the name of its molecule type: the topology's, else what its parts are, such as peptide, GTP or water | string with `*` and `?`, or a list |
|  | `chain` | the author id of its chain | string or list |
| `chains` | `id` | the author chain id | string or list |
| all of them | `scope` | the set to search instead of the default scope | a set, or `system` |

A list means any of its values, and separate keywords must all hold.

The kinds of a part: `polymer` is any polymer strand, `peptide`, `nucleic`, `dna` and `rna` the strands of those kinds (`nucleic` is DNA, RNA or a hybrid), and `glycan` a tree of saccharides or one sugar. A part that is no polymer is a `lipid` when its residues are lipids, an `ion` when it is one ion, a `water`, or else a `ligand`. A molecule has no kind of its own, since a glycoprotein holds a strand and its glycans: select the parts, and `molecules(parts(kind='peptide'))` gives the molecules they are in.

### Positional arguments

1. At most one, written first (E0203).
2. A set or an array of sets, bonds included: `particles()` merges their particles; a grouping selector keeps the units that contain any of them, whole within the scope.
3. A string is the name, as a glob. For `chains` it is the author id.
4. An int, a range or a list of ints is the id for `residues` and the serial for `particles`.
5. A positional value the selector does not take, such as a string for `parts` or an int for `chains`, is E0205.

### Predefined names

`protein`, `nucleic`, `water`, `ion`, `lipid`, `backbone` and `sidechain` are sets: `protein` is the particles of the peptide parts, with their caps and the modified residues bonded into them, `nucleic` those of the nucleic acid parts, and `water`, `ion` and `lipid` those of the parts of that kind. `system` is every particle, whatever the scope; `system('name')` names another loaded system once viamd holds several.

### Scope

- In a script the scope is the system. In a filter it is the set viamd evaluates it within: a representation group, or one reference structure of a density volume.
- Positions count from the scope's start; ids are kept. The scope is all there is, so a unit partly inside it appears with only its particles inside.
- `scope=system` reaches outside. Spatial selectors search the scope, but their points may lie outside it.
- Scopes nest. A scope that changes per frame makes everything evaluated within it frame stage.

```python
# on T
n_all  = len(residues)                                 # -> 880: 442 in A, 437 in B, 1 in C
ch     = chains                                        # -> set[3]
a      = chains('A')                                   # -> set[1]
pa     = particles(chain='A')                          # -> set
n_gtp  = count(particles(resname='GTP'))               # -> 32
both   = residues(501)                                 # -> set[2]: GTP in A, GDP in B; the shape shows it
p501   = particles(resid=501)                          # warning W0602: the merge would hide two chains
wat    = residues(kind='water')                        # -> set[8]
so4    = residues('SO4')                               # -> set[3]
n_part = len(parts)                                    # -> 17: two strands, six ligands, one ion, eight waters
str_   = parts(kind='polymer')                         # -> set[2]: chain B's strand goes on over its gap
lig    = parts(kind='ligand')                          # -> set[6]: GTP, three SO4, GDP and VLB
mol    = molecules(parts(kind='peptide'))              # -> set[2]: the molecules the strands are in
n_mol  = len(molecules)                                # -> 17
first5 = residues[0:5]                                 # -> set[5]
firstb = residues(scope=particles(chain='B'))[0]       # -> set: the first residue of chain B
miss   = particles(name='CA', resid=39, chain='A')     # warning W0601: matches nothing
bad    = chains(name='A')                              # error E0202: chains have no name
two    = residues(45, 'LYS')                           # error E0203
num    = chains(0)                                     # error E0205: chains take an id string; positions are subscripts
_nearv = residues(within(5.0, particles(resname='VLB')))   # -> set[?], frame
```

```python
# on P
strand = parts(kind='peptide')                         # -> set[1]: no chain ids, but the strand is a part
n_mol  = len(molecules)                                # -> 1
c      = chains('A')                                   # error E0207: P gives no chain ids; parts finds its strand
```

```python
# on T, as filters in a representation group whose own filter is
# residues(within(8.0, particles(resname='VLB'))): 22 residues and 241 particles,
# twenty amino acids of chain B, GDP and VLB itself
particles(name='CA')                                   # the CAs of the group's twenty amino acids
residues[0]                                            # the group's first residue: residue 11 of chain B
residues(177)                                          # chain B's residue 177; chain A's lies outside the group
residues(250)                                          # warning W0601: not among this group's 241 particles
within(4.0, particles(resid=181, chain='B', scope=system))   # the group's particles near residue 181 of chain B, which lies outside it
```

## Indexing: positions, ids, at(), sets as indices, by=

Positions count from 0 and are written as subscripts; ids are the file's numbers and names, written as filter values. A per-particle array can also be indexed by a set.

| Form | Means | Result |
| --- | --- | --- |
| `a[i]`, `a[i:j]`, `a[::2]`, `a[-1]` | positions along the first dimension, half-open, as in NumPy | one element, or a shorter dimension |
| `a[i, j]`, `a[:, 0]`, `a[..., 2]` | positions along several dimensions |  |
| `a[[3, 0, 7]]` | the positions listed, in that order |  |
| `a[:, None]` | a new dimension of 1, for broadcasting |  |
| `a.T` | the dimensions reversed |  |
| `c.x`, `c.y`, `c.z` | `c[..., 0]`, `c[..., 1]`, `c[..., 2]` |  |
| `s[i]`, `s[i:j]` on a set | its particles by position | a set |
| `at(x, frame=k)` | frame k of a frame value, read when compiling | a static value |
| `rmsf[s]` with a set | the values of the particles in set s; an intersection | a shorter dimension over those particles |
| `a[mask]` with a bool array matching a's leading dimensions | the entries where the mask is true | a shorter dimension; per frame if the mask is |
| `mean(a, by=g)` | for each set in array g, the reduction of `a[set]` | one value per set in g |
| `labels(a, axis=0)` | the labels of a dimension | `str` array |

- An out-of-range position on a static extent is E0304. On a dynamic extent it is R0101 in that frame, and the value is missing there; the compiler says so beforehand with W0405.
- `by=` works with every reduction and is defined as `[f(a[s]) for s in g]`. A group with no values gives a missing value and W0603, which counts them. `a[s]` and `by=` act on the first dimension of `a` that describes particles; an array with no such dimension is E0311.
- `var` and `std` over frames run in one pass, with Welford's algorithm.

```python
# on T
res_a = residues(chain='A', kind='amino_acid')         # -> set[433]
r0    = res_a[0]                                       # -> set: MET1/A
name0 = labels(res_a)[0]                               # -> 'MET1/A'
last3 = res_a[-3:]                                     # -> set[3]
every = res_a[::2]                                     # -> set[217]
c     = com(res_a)                                     # -> float[433][3] Å
x     = c[:, 0]                                        # -> float[433] Å
x2    = c.x                                            # -> float[433] Å: the same
first = particles[0:10]                                # -> set: the first ten particles
out   = res_a[433]                                     # error E0304
_z    = positions(particles(protein)).z                # -> float[6784] Å
z_res = mean(_z, by=residues)                          # warning W0603, W0405 -> float[880] Å | missing: 15 residues have no values
```

```python
# on P
_p     = positions(particles)                          # -> float[153][3] Å, frame
rmsf   = sqrt(sum(var(_p, axis='frame'), axis=1))      # -> float[153] Å, trajectory
ca     = rmsf[particles(name='CA')]                    # -> float[15] Å
per_r  = mean(rmsf, by=residues)                       # -> float[15] Å
same   = [mean(rmsf[r]) for r in residues]             # -> float[15] Å
end0   = at(positions(particles(name='CA')), frame=0)  # -> float[15][3] Å, static
_close = residues[distance(residues.com, residues[0].com) < 6.0]   # -> set[?], frame
```

## Mapping a def over arrays

A parameter can be annotated with a kind. The kind says how much of an argument one call takes. When an argument has more, the call maps over the extra dimensions on its left and the result gains them. The mapped dimensions of different parameters broadcast as for any function and keep their descriptions, so the result is labelled. This is NumPy's generalized ufunc, and it is how the built-ins already behave: `com` maps over an array of sets, and a `point` parameter broadcasts its leading dimensions.

| Kind | One call takes | Given `set[433]` | Given `float[32][3]` |
| --- | --- | --- | --- |
| `set` | one set | maps over the 433 | E0313 |
| `group` | one set, or one point as a group of one location | maps over the 433 | maps over the 32 points |
| `point` | `float[3]`, or a set known to hold one particle | E0302 | maps over the 32 points |
| `points` | `float[n][3]`, or one set as its particles' positions | maps over the 433 | one call |
| `float`, `int`, `bool` | one element | E0313 | maps over every dimension |

- An unannotated parameter takes its argument whole, as in Python. When an unannotated argument has a dimension that describes the same units as a mapped one, W0607 suggests the annotation that pairs them.
- Inside the body the mapped dimensions are not visible, so a reduction there reduces only what one call took. An extent that changes from call to call, such as the number of particles in each residue, behaves like `?`: it broadcasts only with itself or 1, and hover shows its range, `float[4..14][3]`.
- A group that is empty in some frame gives a missing entry, which reductions skip.
- Any other annotation, and a return annotation, is E0312. An int where a point or group is expected is E0302, with the particles(serial=) fix; any other argument that does not fit its kind, such as a float for int, is E0313. particles(p), per-particle properties and shape() keep the identity of p's extents, so the weights in gyration pair with its positions.

```python
# on T
res_a = residues(chain='A', kind='amino_acid')       # -> set[433]
vlb   = particles(resname='VLB')                     # -> set

def contacts(r: set, lig: set, cutoff: float = 4.0):
    """Particle pairs closer than cutoff."""
    return sum(distance(positions(r)[:, None], positions(lig)) < cutoff)

n1  = contacts(res_a, vlb)                           # -> int[433], frame: one count per residue
n2  = contacts(res_a[0], vlb)                        # -> int, frame: one residue, no mapping
cut = full(shape(res_a), 4.0)                        # -> float[433]: one cutoff per residue
n3  = contacts(res_a, vlb, cut)                      # -> int[433], frame: cutoff: float pairs it per residue
rg  = gyration(res_a)                                # -> float[433] Å, frame: points maps over the residues

def in_contact(r: set, lig: set, cutoff=4.0):
    return distance(r, lig, between='closest') < cutoff

c   = in_contact(res_a, vlb, cut)                    # warning W0607 -> bool[433][433]: cutoff is passed whole
```

## The prelude

Functions marked *prelude* in the tables below are written in VIAL and compiled before every script. They are ordinary `def`s: specialized per call, fused like a user's code, shown in hover with their source. This section is generated from `vial/prelude.vial`, and the examples that use each function elsewhere in this reference are its tests. An error or warning raised inside a prelude function is reported at the user's call and names the function, as “distance: residues(1..4, chain='A') holds 4 sets of several particles”.

```python
def distance(a, b, between=None):
    if between is None:
        r = norm(displacement(a, b))
        return shown(r, line(a, b, label=r))
    return _between(a, b, between)

def _between(a: group, b: group, op):
    assert op in ('closest', 'farthest'), "between= takes 'closest' or 'farthest'"
    pa = positions(a)                            # one group each; maps over arrays of groups
    pb = positions(b)
    d, i, j = closest_pair(pa, pb) if op == 'closest' else farthest_pair(pa, pb)   # i and j only when drawn
    return shown(d, line(pa[i], pb[j], label=d))

def centroid(s):
    return com(s, weights=None)

def gyration(p: points, weights='mass'):
    assert weights in ('mass', None), "weights= takes 'mass' or None"
    w = mass(particles(p)) if weights == 'mass' else ones(shape(p)[0])
    c = com(p, weights)
    r = norm(displacement(c, p))                 # minimum image from a periodic-aware centre
    rg = sqrt(sum(w * r**2) / sum(w))
    return shown(rg, sphere(c, rg))

def normal(p: points):
    c = com(p, weights=None)
    d = displacement(c, p)                       # float[n][3]
    w, v = eigh(d.T @ d)                         # eigenvalues ascending
    n = v[:, 0]                                  # the direction of least spread
    return shown(n, arrow(c, direction=n, length=5.0))

def rmsd(p: points, ref: points, on=None, weights=None, superpose=True):
    assert weights in ('mass', None), "weights= takes 'mass' or None"
    q = fit(p, ref, on=on, weights=weights) if superpose else p
    w = mass(particles(q)) if weights == 'mass' else ones(shape(q)[0])
    return sqrt(sum(w * norm(q - ref)**2) / sum(w))

def previous(x, n=1):
    def step(s, v):
        return concatenate((s[1:], v[None]))
    h = fold(step, x, init=[x for _ in range(n + 1)], skip=False)
    return missing_where(frame() < n, h[0])

def cumulative(x, op='sum'):
    assert op in ('sum', 'min', 'max'), "op= takes 'sum', 'min' or 'max'"
    def step(s, v):
        if op == 'max':
            return max(s, v)
        if op == 'min':
            return min(s, v)
        return s + v
    return fold(step, x, init=zeros_like(x) if op == 'sum' else x)

def ema(x, alpha):
    def step(s, v):
        return s + alpha * (v - s)
    return fold(step, x, init=x)
```

The prelude uses nothing a script cannot, drawing included: `shown` attaches what hover draws, as any def may (see Drawing). `_between` is VIAL; only the search for the closest or farthest pair, `closest_pair` and `farthest_pair`, is built in, since a spatial search is not an equation over arrays. `_between`, `gyration`, `normal` and `rmsd` annotate their parameters with kinds, so they map over arrays of sets: `_between` only ever sees one group per argument. `_between` calls `norm(displacement())` rather than `distance`, since a cycle of calls is recursion (E0802). `gyration` and `normal` measure from a periodic-aware centre under the minimum image, so they are right for any group smaller than half the box, whether or not the trajectory was made whole; positions themselves are as the trajectory gives them. A fold's step is a nested `def`, which may read the names of the function around it, as `cumulative`'s reads `op` and `ema`'s reads `alpha`; there is no `lambda`. `atoms()` and `beads()` are built-in shorthands for `particles(kind=...)`, since forwarding keywords would need `**kwargs`.

## Built-in functions: geometry and structure

Geometry takes points, and a parameter declared `point` or `points` is the only place where a selection turns into coordinates. A set of several particles never quietly becomes one point: write `s.com` or `centroid(s)` for its centre, or `between='closest'` for its closest particles. MDAnalysis, MDTraj and GROMACS also ask for centres explicitly.

| Parameter kind | Accepts | Becomes |
| --- | --- | --- |
| `point` | a float array whose last dimension is 3, such as `[5.0, 0.0, 2.3]` or `s.com` | itself; leading dimensions broadcast |
| `point` | a set the compiler knows holds one particle, such as `ca[0]` or `particles(name='MG')` | that particle's position |
| `point` | an array of such sets, such as `smiles_atoms(...)[:, 0]` | one position per set |
| `point` | a set that may hold several particles | E0302, with fixes that write `.com` or `between='closest'` |
| `point` | an int | E0302, with a fix that writes `particles(serial=12)` |
| `points` | `float[n][3]` | itself |
| `points` | a set | the positions of its particles; an array of sets maps, one result per set |

### Geometry

| Function | Result | Notes |
| --- | --- | --- |
| `positions(s)` | `float[n][3]` Å | frame stage; on a point, itself as float\[1\]\[3\]; n is `count(s)` |
| `mass(s)`, `charge(s)`, `element(s)`, `name(s)`, `bfactor(s)` | `float[n]` u, e, or `str[n]` | per-particle properties of a set; static, except `bfactor` where the file gives one per frame |
| `particles(x)` on an array | `set` | the particles an array's first particle dimension describes, so `mass(particles(fit(positions(s), ref)))` works after a fit; an array with no such dimension is E0311 |
| `com(s, weights='mass')` | `float[3]` Å | periodic-aware; takes a set or points, and maps over arrays of sets; also written `s.com` |
| `centroid(s)` | `float[3]` Å | prelude: `com(s, weights=None)` |
| `distance(a: point, b: point)` | `float` Å | prelude: `norm(displacement(a, b))`; leading dimensions broadcast. Drawn as a line labelled with the value |
| `distance(a, b, between='closest')` | `float` Å | prelude: `_between(a: group, b: group)`, which is `min(distance(positions(a)[:, None], positions(b)))` for one group each and maps over arrays of groups; `max` for `'farthest'`. A group is a set, or a point on its own. Replaces v1's `distance_min`, `distance_max` and `distance_pair`. The pair comes from `closest_pair`, a spatial search, so no n×m matrix is computed or stored. Drawn as a line between the pair |
| `displacement(a: point, b: point)` | `float[3]` Å | b minus a, minimum image |
| `closest_pair(p: points, q: points)` | `(float Å, int, int)` | the shortest distance from a point of p to a point of q under the minimum image, and the two points' indices; built in, with a spatial search, so its time grows with n + m for compact groups. The indices are computed only where they are used, as when hover draws the pair. Missing, with W0601, when p or q is empty |
| `farthest_pair(p: points, q: points)` | `(float Å, int, int)` | the longest distance, and its two indices |
| `angle(a, b, c: point)` | `float` ° | also takes one `float[3][3]`, three points as rows; the one-argument form maps over leading dimensions, such as `set[k][3]` from `smiles_atoms` |
| `dihedral(a, b, c, d: point)` | `float` °, in (−180, 180\] | also takes one `float[4][3]`, or four one-particle sets; the one-argument form maps over leading dimensions, such as `set[k][4]` from `smiles_atoms` |
| `gyration(p: points, weights='mass')` | `float` Å | prelude: `sqrt(sum(w * norm(displacement(com(p, weights), p))**2) / sum(w))`. Drawn as a sphere of that radius at the centre |
| `normal(p: points)` | `float[3]` | prelude: the direction of least spread, `v[:, 0]` of `eigh(d.T @ d)` with `d = displacement(com(p, weights=None), p)`. Drawn as an arrow from the centre |
| `within(r, x: points)` | `set`, frame | the scope's particles within r of any particle or point of x |

A matrix of distances needs no special function; broadcasting gives it: `distance(c[:, None], c[None, :])`.

### Every form of distance

`distance` is a prelude function with two forms. `distance(a, b)` is `norm(displacement(a, b))` between points. `distance(a, b, between='closest')` calls `_between`, whose parameters are `group`s: one call takes one set, or one point on its own, and computes what `min(distance(positions(a)[:, None], positions(b)))` would, through `closest_pair`; `'farthest'` is `max`, through `farthest_pair`. An array of groups maps, adding its dimensions, so residues of different sizes never meet in one array, and no matrix of every pair is ever made. Everything else follows from the point rule, the kinds and broadcasting: each argument's kind decides what it becomes and which dimension it adds, and the dimensions broadcast as for any function. Every form gives Å under the minimum image and is frame stage because positions are.

The examples use these names:

```python
# on T
p      = particles(resname='GTP').com               # -> float[3] Å
ps     = positions(particles(resname='GTP'))        # -> float[32][3] Å
mg     = particles(name='MG')                       # -> set: one particle
gtp    = particles(resname='GTP')                   # -> set: 32 particles
vlb    = particles(resname='VLB')                   # -> set: 59 particles
res_a  = residues(chain='A', kind='amino_acid')     # -> set[433]
cas    = [particles(r, name='CA') for r in res_a]   # -> set[433]: one particle each
near   = within(5.0, vlb)                           # -> set, frame: its size changes per frame
near_r = residues(within(5.0, vlb))                 # -> set[?], frame
```

| Kind | On T | Without `between=` | With `between=` | Adds |
| --- | --- | --- | --- | --- |
| a position, `float[3]` | `p`, `[0.0, 0.0, 0.0]` | itself | one location | — |
| positions, `float[n][3]` | `ps`, `res_a.com` | n points | n locations, each on its own | `[n]` |
| a set of one particle | `mg`, `cas[0]`, `gtp[0]` | its position | its particle | — |
| a set of several particles | `gtp`, `vlb` | E0302; `gtp.com` is a point | its particles | — |
| an array of one-particle sets | `cas` | one position per set | one particle per set | `[k]` |
| an array of sets | `res_a` | E0302; `res_a.com` is `float[433][3]` | each set's particles | `[k]` |
| a set that changes per frame | `near` | E0302; `near.com` is a point | its particles in each frame; missing in a frame where it is empty | — |
| an array of sets that changes per frame | `near_r` | E0302; `near_r.com` is `float[?][3]` | each set's particles | `[?]` |
| an int | `12` | E0302 | E0302 | — |
| a set known to be empty | `particles(resname='ATP')` | missing, with W0601 | missing, with W0601 | — |

Under `between=` a set is a group of its particles and a position is a location on its own, so `ps` against `vlb` is 32 closest approaches, `float[32]`, never one.

So a point against a point is `float`, a point against `cas` or `res_a.com` is `float[433]`, two arrays of the same extent pair element by element, and two different `?` or two different extents are E0301. A dimension of unknown length, `[?]`, broadcasts only with 1 or with the same `?`. A result with `[?]` cannot be published (E0401) until it is reduced. The forms beyond one dimension, and the edge cases:

```python
# on T, with the names above
cm    = distance(res_a[:, None], res_a, between='closest')    # -> float[433][433] Å: a residue contact map
cc    = distance(res_a.com[:, None], res_a.com)               # -> float[433][433] Å: centre to centre
far   = distance(gtp, vlb, between='farthest')                # -> float Å
one   = distance(res_a[0:1], res_a.com)                       # error E0302: a set[1] of several particles
st    = distance(res_a[0:1].com, res_a.com)                   # -> float[433] Å: the [1] stretches
lab   = distance(cas[0:2], cas[2:4])                          # -> float[2] Å, labelled MET1/A–GLU3/A, ARG2/A–CYS4/A: slices of one array, so no W0605
zero  = distance(cas, cas)                                    # -> float[433] Å, all zero
self_ = distance(gtp, gtp, between='closest')                 # -> float Å, zero
n3    = sum(distance(near_r, mg, between='closest') < 3.0)     # -> int, frame: reduced, so it can be published
bad   = distance(near_r, cas, between='closest')              # error E0301: [?] and [433] do not broadcast
empty = distance(particles(resname='ATP'), mg)                # warning W0601, W0405 -> float Å | missing: empty on T, so missing in every frame

# v1's distance_pair: every particle of a against every particle of b, minimum image
pairs = distance(positions(gtp)[:, None], positions(vlb))     # -> float[32][59] Å: row i is GTP particle i
flat  = flatten(pairs)                                        # -> float[1888] Å: v1's layout, b varying fastest
four  = distance(positions(particles(res_a[0:4]))[:, None], positions(vlb))   # -> float[34][59] Å: v1 merged an array of sets; v2 says so with particles()
```

### Structure

| Function | Result | Notes |
| --- | --- | --- |
| `fit(p, ref, on=None, weights=None)` | `float[n][3]` Å | superposes p onto ref using the particles in `on`, applied to all of p |
| `rmsd(p, ref, on=None, weights=None, superpose=True)` | `float` Å | prelude: the weighted root mean square of `q - ref`, with `q = fit(p, ref)`, or `p` itself when `superpose=False`; weights are equal, or by mass with `weights='mass'`, in the fit and the mean alike |

Positions from the same set pair by particle. Positions from different sets pair in order: different counts do not broadcast (E0301), and different particles paired by position give W0605, which names the first pair that differs; slices of one set do not warn. To pair by residue id and particle name, pair the sets first with `common(sa, sb, by='label')`; it keeps only the particles present in both.

```python
# on T
gtp  = particles(resname='GTP')
gdp  = particles(resname='GDP')
vlb  = particles(resname='VLB')
d1   = distance(gtp.com, gdp.com)                      # -> float Å: centre to centre
d2   = distance(gtp, vlb, between='closest')           # -> float Å: the closest pair of particles
d3   = distance(gtp, gdp)                              # error E0302: GTP has 32 particles; write .com or between=
pg   = particles(resname='GTP', name='PG')             # -> set: one particle
d4   = distance(pg, particles(name='MG'))              # -> float Å: single particles are points
a    = angle(gtp.com, gdp.com, vlb.com)                # -> float °
pos  = positions(gtp)                                  # -> float[32][3] Å, frame
cg   = com(gtp)                                        # -> float[3] Å
near = within(5.0, vlb)                                # -> set, frame
bad  = distance(12, gtp.com)                           # error E0302
```

```python
# on P
ca   = particles(name='CA')
tw   = dihedral(ca[0], ca[1], ca[2], ca[3])            # -> float °, frame
ref0 = at(positions(ca), frame=0)                      # -> float[15][3] Å, static
r    = rmsd(positions(ca), ref0)                       # -> float Å, frame
raw  = rmsd(positions(ca), ref0, superpose=False)      # -> float Å, frame
_f   = fit(positions(particles), at(positions(particles), frame=0), on=ca)   # -> float[153][3] Å, frame
half = rmsd(positions(ca[0:7]), positions(ca[8:15]))   # -> float Å, frame: slices of one set pair in order without W0605
odd  = rmsd(positions(ca[0:7]), positions(ca[7:15]))   # error E0301: [7][3] and [8][3] do not broadcast
```

Later, ported from v1 after the first version: `rdf`, `sdf`, `density`, contact maps, hydrogen bonds, secondary structure and ring pucker.

## Built-in functions: reductions, math and arrays

Reductions remove dimensions, math works elementwise, and the array functions say how values pair or where they go.

### Reductions

| Function | Result | Unit |
| --- | --- | --- |
| `sum`, `mean`, `min`, `max` | the dimensions not reduced | the argument's |
| `std`, `var` | the dimensions not reduced | the argument's; squared for `var` |
| `any`, `all` | `bool` | none |
| `argmin`, `argmax` | `int`: a position along the reduced dimension | none |
| `len(a)` | `int`: the extent of the first dimension | none |

`count(s)` is for sets only: the number of particles in each set, mapping over an array of sets like `com`. How many values of a bool array are true is `sum(mask)`, an `int`.

- `axis=` takes an int, negative from the end, a tuple, or `'frame'`. Without it every dimension is reduced, as in NumPy, but never the frame.
- `by=` reduces in groups; see the indexing section. A reduction over an extent that is empty in some frame gives its identity where it has one, 0 for sum and len, false for any and true for all, and otherwise a missing value, for min, max, mean, std, var, argmin and argmax, with W0405.
- `len(a)` is the extent of the first dimension, so `len(residues)` is the number of residues. On a single set it is the number of particles, the same as `count(s)`.
- `min(a, b)` and `max(a, b)` with two arguments work elementwise, as Python's do; with one argument they reduce.
- Sums accumulate in 64 bits in one fixed order, the same in every backend. An integer sum that does not fit in 32 bits is R0102 rather than wrapping around.

### Math

| Function | Takes | Gives |
| --- | --- | --- |
| `abs`, `floor`, `ceil`, `round`, `clip(x, lo, hi)` | any unit | the same unit |
| `sqrt` | any unit | the unit's square root |
| `exp`, `log`, `log10` | no unit | no unit |
| `sin`, `cos`, `tan` | an angle | no unit |
| `asin`, `acos`, `atan`, `atan2(y, x)` | no unit | degrees |
| `dot(u, v)` | arrays whose last dimensions match | the sum of products over the last dimension, row by row; NumPy calls it `vecdot` |
| `a @ b` | matrices and vectors | the matrix product, as NumPy's `matmul` |
| `cross(u, v)` | last dimension 3 | the cross product, row by row |
| `norm(v, axis=-1)`, `normalize(v)` | any float array | the length over a dimension; v divided by it |
| `where(c, a, b)` | a bool array and two values | elementwise choice |
| `unitless(x)` | any unit | no unit, without a warning |
| `float(x)`, `int(x)` | any | the same values as float or int, keeping the unit; int truncates toward zero |

### Arrays

| Function | Does |
| --- | --- |
| `diff(a, axis=0)` | differences of neighbours along a dimension, as NumPy's |
| `zip(a, b)` | Python's: pairs elements in order, for comprehensions |
| `range(n)` | Python's: positions 0 to n−1 |
| `common(a, b, by='id')` | the elements of a and b that match, as two arrays in a's order. `by='id'` matches the residue id, and for particles also the particle name; `by='label'` adds the residue name. Neither looks at the chain, so chain A pairs with chain B |
| `labels(a, axis=0)` | the labels of a dimension |
| `spread(values, groups, combine=None)` | gives each particle of a group that group's value, the reverse of `by=`; the result runs over `particles(groups)`. `combine='mean'`, `'max'` or `'first'` is required where groups overlap (E0306) |
| `missing_where(c, x)` | x with the entries where c is true marked missing; the type becomes optional |
| `valid(x)` | a bool array: true where x is present |
| `filled(x, v)` | x with its missing entries replaced by v; the type is no longer optional |

```python
# on P
ca    = particles(name='CA')
c1    = distance(ca.com, ca[0])                # -> float Å: from the centre of the CAs to the first CA
dd    = distance(positions(ca), ca[0])         # -> float[15] Å, frame: one distance per CA
total = sum(dd)                                # -> float Å, frame
far   = argmax(dd)                             # -> int, frame
n5    = sum(dd > 5.0)                          # -> int, frame
lo    = min(dd, 5.0)                           # -> float[15] Å, frame
mn    = min(dd)                                # -> float Å, frame
w     = where(dd > 5.0, dd, 0.0)               # -> float[15] Å, frame
bv    = diff(positions(ca), axis=0)            # -> float[14][3] Å, frame: CA to CA vectors
bl    = norm(bv)                               # -> float[14] Å, frame
cosb  = dot(normalize(bv[:-1]), normalize(bv[1:]))   # -> float[13], frame
t_m   = mean(dd, axis='frame')                 # -> float[15] Å, trajectory
all_m = mean(dd, axis=(0, 'frame'))            # -> float Å, trajectory
cs    = cos(angle(ca[0], ca[1], ca[2]))        # -> float, frame
pairs = [distance(x.com, y.com) for x, y in zip(residues[0:7], residues[8:15])]   # -> float[7] Å, frame
```

```python
# on T
res_a = residues(chain='A', kind='amino_acid')
vlb   = particles(resname='VLB')
d_res = distance(res_a, vlb, between='closest')     # -> float[433] Å: each residue's closest approach to VLB
shade = spread(d_res, res_a)                        # -> float[n] Å over the particles of res_a; the groups do not overlap
rb    = residues(chain='B', kind='amino_acid')      # -> set[432]
a, b  = common(res_a, rb, by='id')                  # -> set[424] and set[424]: the ids both chains have
dab   = distance(a.com, b.com)                      # -> float[424] Å, no warning: common paired them; labelled MET1/A–MET1/B, ARG2/A–ARG2/B, …, CYS4/A–ILE4/B, …
```

## Arrays: building, reshaping, joining and linear algebra

Arrays can be built, reshaped, joined and used in linear algebra. Everything follows NumPy unless a row says otherwise.

### Building

| Function | Result |
| --- | --- |
| `zeros(shape)`, `ones(shape)`, `full(shape, value)` | a float array of that shape, with no unit until it meets one, as a bare number; unit= gives one. zeros\_like(x) and ones\_like(x) keep x's dtype, shape, unit and descriptions |
| `shape(a)` | a's extents; they also remember what each dimension runs over |
| `linspace(lo, hi, n)`, `arange(lo, hi, step)` | a one-dimensional float array |
| `eye(n)` | the n by n identity matrix |

A shape is an int, a tuple of ints such as `(10, 3)`, or `shape(a)`. `zeros(shape(res_a))` is `float[433]` and carries the residue labels; `zeros((10, 3))` carries no descriptions.

### Reshaping

| Function | Result |
| --- | --- |
| `reshape(a, shape)` | the same values in a new shape with as many elements; `-1` infers one extent (E0307) |
| `flatten(a)` | one dimension; labels combine, such as `CA/ALA12/A:x` |

Values are never changed in place, so NumPy's distinction between views and copies never shows.

### Joining

| Function | Result |
| --- | --- |
| `concatenate(list, axis=0)` | the arrays of a tuple, or of a list literal written in the call, joined along an existing dimension rather than stacked; the other extents must match (E0309) |
| `[a, b, c]` | the list literal stacks: a new first dimension, one entry per item of the same shape. Items that are top-level names give their names as labels, so `[d_gtp, d_gdp]` is labelled `d_gtp`, `d_gdp`. There is no `stack()` |

### Linear algebra

| Function | Result |
| --- | --- |
| `a @ b` | the matrix product, as NumPy's `matmul` |
| `dot(u, v)` | the sum of products over the last dimension, row by row |
| `covariance(x, axis='frame')` | `float[n][n]` for x of shape `[n]`; over frames it streams, in one pass |
| `eigh(m)` | eigenvalues in ascending order and eigenvectors in the columns, as NumPy's |
| `trace(m)`, `det(m)`, `inv(m)`, `solve(m, b)` | over square matrices (E0310) |

An outer product is broadcasting: `u[:, None] * v[None, :]`.

### Example: principal components

```python
# on P
ca    = particles(name='CA')
ref0  = at(positions(ca), frame=0)
_p    = fit(positions(ca), ref0)                         # -> float[15][3] Å, frame
_q    = flatten(_p)                                      # -> float[45] Å, frame
C     = covariance(_q, axis='frame')                     # -> float[45][45] Å², trajectory
w, v  = eigh(C)                                          # -> float[45] Å² and float[45][45], trajectory
frac  = w[::-1] / sum(w)                                 # -> float[45], trajectory: the share of variance, largest first
pcs   = (_q - mean(_q, axis='frame')) @ v[:, [-1, -2]]   # -> float[2] Å, trajectory, per frame, pass 2: PC1 and PC2
mode1 = reshape(v[:, -1], shape(_p))                     # -> float[15][3]: the largest mode as motion per particle
```

A covariance over n values takes n² numbers: 3,000 coordinates, the CAs of a 1,000-residue protein, take 36 MB.

```python
# on T
res_a = residues(chain='A', kind='amino_acid')
w0    = zeros(shape(res_a))                              # -> float[433]: with the residue labels
grid  = zeros((10, 3))                                   # -> float[10][3]
r     = linspace(0.0, 20.0, 81)                          # -> float[81]
I     = eye(3)                                           # -> float[3][3]
gtp   = particles(resname='GTP')
gdp   = particles(resname='GDP')
vlb   = particles(resname='VLB')
d_gtp = distance(gtp.com, vlb.com)
d_gdp = distance(gdp.com, vlb.com)
dd    = [d_gtp, d_gdp]                                   # -> float[2] Å, labelled d_gtp, d_gdp
ab    = concatenate([residues(1..4, chain='A'), residues(1..4, chain='B')])     # -> set[8]
bad   = reshape(com(res_a), (400, 3))                    # error E0307: 1,299 values do not fill 1,200
```

## Built-in functions: sets, patterns, bonds, attributes and frames

### Sets

| Form | Gives |
| --- | --- |
| `a and b`, `a or b`, `not a` | intersection, union, and complement within the scope |
| `any(arr, axis=)`, `all(arr, axis=)` | the union, or the intersection, of an array of sets along a dimension; the reductions of `or` and `and`. There is no `union()` |
| `particles(x)` | the particles of a set or an array of sets, as one set |
| `count(s)` | `int`: the number of particles |

```python
# on T
lig  = particles(resname=['GTP', 'GDP', 'VLB'])        # -> set
sa   = particles(chain='A') and protein                # -> set
side = protein and not backbone                        # -> set
both = any(residues(501))                              # -> set: GTP and GDP merged
none = count(particles(chain='C') and protein)         # -> 0
```

### Patterns

| Function | Result |
| --- | --- |
| `smiles(pattern, scope=None, level='molecule', mode='unique')` | `set[k]`: one set per match |
| `smiles_atoms(pattern, ...)` | `set[k][w]`: one particle per pattern atom, in pattern order; only the atoms tagged with a class such as `[C:1]` when the pattern has classes, in class order |

`level` is where a match may lie: within a `'molecule'`, a `'residue'` or a `'part'`. `mode` is `'unique'`, one match per distinct set of atoms, or of tagged atoms when the pattern has classes; `'all'`, every mapping; `'one_per_unit'`; or `'disjoint'`, which keeps the first of overlapping matches. Matches may overlap, and the compiler knows where. A hydrogen in a pattern is mapped onto one hydrogen of its atom, never searched.

Hydrogens and charges follow md\_match, with one exception. An organic-subset atom such as `C` leaves its hydrogen count and charge open, and a bracket atom such as `[CH2]` fixes them. A bracket atom that only adds a class, such as `[C:1]`, leaves them open too, so tagging an atom never changes what it matches; `[CH1:1]` fixes the count. md\_match needs this exception added. Where the file has no bond orders, md\_chem perceives them and the pattern's orders are tested on those, as md\_match does, so a carboxylate's two oxygens both match =O.

```python
# on T
po4  = smiles("O=P(O)(O)O")                            # -> set[5]: three phosphates of GTP, two of GDP
# on P
phi  = dihedral(smiles_atoms("[C:1](=O)[N:2][C:3][C:4]=O"))   # -> float[14] °, frame: φ of residues 2 to 15; residue 15's O and OXT give one match, since matches are unique over the tagged atoms
```

### Bonds

| Function | Result |
| --- | --- |
| `bonds(s=None, kind=None)` | `set[n][2]`: per bond, its two particles as one-particle sets, for the bonds with both ends in s, filtered by kind. There is no `bonds` dtype |
| `formed(b, ref)`, `broken(b, ref)` | `set[?][2]`: the bonds present now and absent in ref, or the reverse; a set difference by identity |

Where bonds change between frames, as in a reactive run, these are frame stage.

```python
# on T
cov  = bonds(particles(resname='VLB'))                 # -> set[n][2]: per bond, its two particles
nb   = len(cov)                                        # -> int: the number of bonds
bl   = distance(cov[:, 0], cov[:, 1])                  # -> float[n] Å, frame: every bond length
long = cov[bl > 1.6]                                   # -> set[?][2], frame
```

### Attributes

`attr(path)` reads an attribute of the run, such as an energy term, with its unit. A path with `*` matches several and gives one member each, labelled by name.

```python
# on a run with an energy file
pot  = attr("edr/Potential")                           # -> float kJ/mol, frame
coul = attr("edr/coul_sr_*")                           # -> float[k] kJ/mol, frame: labelled by the part * matched
lj   = attr("edr/lj_sr_*")                             # -> float[k] kJ/mol, frame
tot  = sum(coul)                                       # -> float kJ/mol, frame: summed within each frame
inter = sum(coul) + sum(lj)                            # -> float kJ/mol, frame: coul + lj would pair members by position
```

### Frames

| Function | Result | Stage |
| --- | --- | --- |
| `mean`, `var`, `std`, `min`, `max`, `sum` with `axis='frame'` | the frame removed | trajectory; streams while frames load |
| `histogram(x, axis='frame', bins=, range=)` | `float[bins]`, with bin edges: counts, or a density with density=True; values outside range are counted apart, below and above | trajectory; streams |
| `smooth(x, n)` | keeps the frame | trajectory, per frame; a centred window of n frames, n odd, that shrinks at the ends |
| `derivative(x)` | keeps the frame; unit per ps | trajectory, per frame; central differences over time(), one-sided at the ends |
| `autocorrelation(x, max_lag=)` | `float[max_lag+1]` | trajectory; the mean removed and normalized to 1 at lag 0, so it has no unit; complete when its pass ends |
| `fold(step, x, init=)` | the state after each frame, with the type and shape of `init`; `step(state, x)` is a `def` | trajectory, per frame; frames in order, streams as they arrive in order; skip=False also steps over missing frames |
| `previous(x, n=1)`, `cumulative(x, op='sum')`, `ema(x, alpha)` | x from n frames back, missing in the first n and where that frame was missing; a running sum, min or max; an exponential moving average | trajectory, per frame; prelude, as folds |
| `at(x, frame=k)` | frame k | static |
| `time()`, `frame()` | `float64` ps, the one 64-bit value, rounded to 32 bits where it meets other data; `int` | frame |

```python
# on P
ca   = particles(name='CA')
ends = distance(ca[0], ca[14])                         # -> float Å, frame
h    = histogram(ends, axis='frame', bins=32, range=(0.0, 50.0))   # -> float[32], trajectory
sm   = smooth(ends, 5)                                 # -> float Å, trajectory, per frame
ac   = autocorrelation(ends, max_lag=10)               # -> float[11], trajectory
f    = frame()                                         # -> int, frame
rate = ends - previous(ends)                           # warning W0405 -> float Å, trajectory, per frame: missing in frame 0

def step(state, x):
    return where(x > 30.0, state + 1, 0)
run  = fold(step, ends, init=0)                        # -> int, trajectory, per frame: frames in a row above 30 Å
long = max(run, axis='frame')                          # -> int, trajectory
```

## Outputs: what is published and how

Every top-level name bound to a number or a set is published under `script/<name>`, and a name starting with `_` is not.

| Value | Published as | Alongside |
| --- | --- | --- |
| Static number | a constant | unit |
| Static array | one constant per label | labels, unit |
| Frame scalar, or a trajectory scalar per frame such as `smooth()` gives | a series over frames | unit |
| Frame value with one dimension | one series per label | labels, unit |
| Frame value with two dimensions | a matrix per frame, for heat maps | labels for each dimension |
| Trajectory scalar | one value | the number of frames it covers |
| Trajectory matrix, such as a covariance | a matrix, up to 64 MB; larger stays in hover (W0404) | labels for each dimension |
| Trajectory value over particles or residues | one value per label; viamd can colour by it | labels, unit |
| Trajectory value per frame with one dimension, such as `pcs` | one series per label | labels, unit |
| Trajectory array not over particles, such as eigenvalues | one constant per label | labels, unit |
| Static matrix, such as `eye(3)` | a matrix | labels for each dimension |
| bool | as an int, 0 or 1, in any of the forms above | — |
| Histogram or autocorrelation | a curve | bin edges, or lags |
| `set`, static or frame, and a frame-stage array of sets, `set[?]` | a selection under its name; a frame-stage one is a recipe, evaluated for the frame being shown |  |
| Static array of sets | one selection under its name, holding the array: filters use the name for the union, or name\[i\] | labels |
| Three dimensions or more, `str`, `range`, tuples, drawings | not published; shown in hover |  |

- A published frame value must have the same shape in every frame (E0401), since it is stored as a series. Selections are the exception: a set or array of sets whose length changes per frame, such as `residues(within(5.0, particles(resname='GTP')))`, is published as a recipe. It holds no values; representations and filters evaluate it for the frame they show. Other values whose length changes are kept private, or reduced with `len` or another reduction.
- A published per-frame value that would take more than 64 MB over the whole run, 16 million values, gives W0404 with its size. It is usually an intermediate missing its `_`.
- Each published value has a fingerprint of what it depends on, so an edit re-evaluates only what it feeds.
- viamd's stored selections are static sets in the script, read-only. Filters see published names and stored selections, never private names.

```python
# on P
ca    = particles(name='CA')                           # published: a selection
_p    = positions(particles)                           # private
rmsf  = sqrt(sum(var(_p, axis='frame'), axis=1))       # published: 153 values, one per particle
ends  = distance(ca[0], ca[14])                        # published: a series of 38 values
h     = histogram(ends, axis='frame', bins=32, range=(0.0, 50.0))   # published: a curve with its bin edges
_c    = positions(ca)
dm    = distance(_c[:, None], _c[None, :])             # -> float[15][15] Å, frame: a matrix per frame
```

## Drawing: what hover shows

Hovering an expression in the editor draws its value in the 3D view, at the frame being shown. A value is drawn by its type, unless a built-in or a `def` attached a drawing with `shown`.

| Value | Drawn as |
| --- | --- |
| `set` | its particles highlighted |
| An array of sets | each set highlighted, the hovered element emphasized |
| A point, or points: `float[3]` Å or `float[n][3]` Å | dots |
| A number per particle or per residue | the particles or residues coloured by value |
| `bonds(...)` | the bonds |
| `angle`, `dihedral` | an arc at the vertex; for a dihedral also the two planes |
| `within(r, x)` | a sphere of radius r around each point of x, and the particles inside highlighted |
| anything else | nothing; hover shows the value as text |

`shown(x, d)` is `x`, with `d` drawn when it is hovered. Its type, unit, stage and pass are `x`'s, and published values are computed exactly as without it: `d` is evaluated only for hover, at the frame being shown, so a drawing costs evaluation nothing. A name bound to `shown(...)`, and a call of a `def` that returns it, are drawn with `d`; passing the value through names and returns keeps its drawing, and any operation on it gives a new value, drawn by its type. Hover draws the hovered expression's own drawing only; its inputs show theirs when their spans are hovered.

| Primitive | Draws |
| --- | --- |
| `dot(p: point)` | a point |
| `line(a: point, b: point, label=None)` | a segment, labelled at its middle when `label=` is given |
| `arrow(a: point, b: point)`, `arrow(a: point, direction=v, length=r)` | a segment with a head; the second form starts at a and runs r along v |
| `sphere(c: point, r)` | a sphere of radius r |
| `arc(c: point, a: point, b: point, label=None)` | the angle at c from a to b |
| `disk(c: point, n, r)` | a disk of radius r through c, with normal n: a plane |
| `label(p: point, x)` | text at p: a number with its unit, which viamd formats as it formats values, or a `str` |
| `highlight(s: set)` | a set's particles |

- A primitive's value is a `drawing`. Primitives broadcast like elementwise functions, so an array of points gives an array of drawings, `drawing[433]`, and `+` combines drawings elementwise. A `def` with kinds maps, so `distance(res_a, vlb, between='closest')` draws 433 lines.
- Lengths are in Å: a bare number for `r` or `length=` is read as Å, as a radius is for `within`. A segment from a to b runs to b's nearest image, as `displacement` does, so no line spans the box.
- A drawing has the stage of its inputs and may be later than the value it is attached to; it is drawn for the frame being shown.
- Drawings carry no colours or styles: viamd draws them in its own, and emphasizes the hovered element of an array.
- Hover draws at most 10,000 primitives. Above that it draws the hovered element only, and says how many were left out.
- A drawing is never published; a top-level name bound to one is shown in hover only.

```python
# on T
gtp   = particles(resname='GTP')                                 # -> set
vlb   = particles(resname='VLB')                                 # -> set
mg    = particles(resname='MG')                                  # -> set
res_a = residues(chain='A', kind='amino_acid')                   # -> set[433]
d     = distance(gtp.com, vlb.com)                               # -> float Å, frame: drawn as a labelled line
dc, i, j = closest_pair(positions(gtp), positions(vlb))          # -> float Å and int and int
l     = line(gtp.com, vlb.com, label=d)                          # -> drawing, frame
both  = l + sphere(mg, 3.0)                                      # -> drawing, frame: 3.0 is read as Å
dd    = shown(d, both)                                           # -> float Å, frame: hover draws the line and the sphere
spokes = line(res_a.com, vlb.com)                                # -> drawing[433], frame: one line per residue
tag   = label(vlb.com, 'VLB')                                    # -> drawing, frame
nv    = normal(vlb)                                              # -> float[3]: hover draws an arrow from VLB's centre
bad   = arrow(gtp, vlb.com)                                      # error E0302: GTP has 32 particles; write .com
```

A script draws exactly as the prelude does:

```python
# on P
def end_to_end(s: set):
    """The distance from the first to the last particle, drawn with both ends marked."""
    a = s[0]
    b = s[-1]
    r = distance(a, b)
    return shown(r, line(a, b, label=r) + dot(a) + dot(b))

ca = particles(name='CA')                                        # -> set
e  = end_to_end(ca)                                              # -> float Å, frame: published as a series; hover draws it
```

## Diagnostics: every code, its message and its fix

Codes are stable: E for errors, W for warnings, R for errors found while evaluating a frame. The hundreds group them: 01 syntax, 02 names and calls, 03 types and shapes, 04 stages and outputs, 05 units, 06 selections and pairing, 07 v1 spellings, 08 unsupported constructs. Where a call is both an unknown name and a v1 spelling, such as `residue(1:4)`, E0701 is reported, since it says more. The parser recognizes v1 call forms such as residue(1:4), whose colon is not Python, so they reach E0701 rather than E0101.

| Code | When | Message, as an example | Fix offered |
| --- | --- | --- | --- |
| E0101 | Unexpected token | unexpected `)` |  |
| E0102 | Unclosed bracket or string | `(` opened on line 3 is never closed |  |
| E0103 | Inconsistent indentation | this line is indented by 3 spaces; the block uses 4 |  |
| E0104 | Nesting too deep: more than 100 levels of brackets, blocks and prefix operators, or an expression more than 1000 operations deep | this nests more than 100 levels deep |  |
| E0201 | Unknown name | `centre` is not defined | `com` |
| W0201 | A name that shadows a built-in | `residues` now names a variable | another name |
| E0202 | A keyword the function does not take | chains have no name; they take `id` | `chains(id='A')` |
| E0203 | A second positional argument | `residues` takes one positional argument; `'LYS'` is the second | `name='LYS'` |
| E0204 | Wrong number of arguments | `angle` takes 3 points; it got 2 |  |
| E0205 | A positional value the selector does not take | `chains` takes an id string, not a position | `chains[0]` |
| E0206 | A keyword value the selector does not know | 'amino' is not a residue kind | kind='amino\_acid' |
| E0207 | A chain where the data gives no chain ids | P gives no chain ids; chains are labels of the data and never inferred | `parts(kind='polymer')` |
| E0301 | Shapes that do not broadcast | `float[433]` and `float[433][3]` do not broadcast | `d[:, None]` |
| E0302 | Not a point: an int, or a set that may hold several particles | `distance` needs a point; `residues(1..4)` holds 4 sets of several particles | `.com` or `between='closest'`; `particles(serial=12)` for an int |
| E0303 | A condition on an array | `if` needs one bool; this is `bool[442]` | `where(...)` |
| E0304 | A position outside a static extent | position 433 is outside `[433]` |  |
| E0305 | A top-level name bound twice | `n` is already bound on line 2 |  |
| E0306 | `spread` on overlapping groups without `combine=` | 12 groups share particles | `combine='mean'` |
| E0307 | `reshape` to a different number of elements | 1,299 values do not fill 1,200 |  |
| E0308 | A seventh dimension | the result would have 7 dimensions; the limit is 6 |  |
| E0309 | `concatenate` of arrays whose other extents differ | `[4][3]` and `[4][2]` differ in dimension 1 |  |
| E0310 | A matrix function on a matrix that is not square | `eigh` needs a square matrix; this is `float[433][3]` |  |
| E0311 | A set index, by=, spread, fit(on=), common or particles(x) on an array with no particle description | zeros((10,)) has no dimension over particles | zeros(shape(res\_a)) |
| E0312 | An annotation that is not a kind, or a return annotation | str is not a parameter kind; the kinds are set, group, point, points, float, int and bool | remove it |
| E0313 | An argument that does not fit its parameter's kind | r: set got float\[433\]\[3\], which is not a set |  |
| E0314 | A static assert that fails | cumulative: op= takes 'sum', 'min' or 'max'; got 'avg' | the assert's message |
| E0315 | A fold step whose result differs from init in dtype, shape or unit | step returns float\[3\] Å; init is float Å |  |
| E0316 | An axis the array does not have | axis=2 on float\[433\]\[3\] | axis=1 |
| E0317 | is on a value that is not static | is compares static values; d changes per frame | == |
| E0401 | A published value, other than a selection, whose shape changes between frames | `d_near` has a different length in each frame | a reduction, or `_d_near` |
| E0402 | A value from another system not pinned to a frame | values from `crystal` need a frame | `at(..., frame=0)` |
| E0403 | A top-level `if` on a value that changes per frame | the condition changes per frame; outputs cannot | `where(...)` |
| W0404 | A published per-frame value over 64 MB for the run | `dm` takes 2.1 GB over 22,000 frames | `_dm` |
| E0404 | at() of a value that depends on other frames | mean\_e is a trajectory value; at() reads one frame |  |
| E0405 | at() past the last frame | frame 500 of a run with 38 frames |  |
| W0405 | A value that may be missing in some frames | near\[4\] is missing in frames with fewer than five residues nearby | none; the series shows gaps |
| W0406 | More than three passes over the trajectory | pcs2 is pass 4; each pass reads every frame again |  |
| W0501 | Units that do not match | adds Å to degrees; the result's unit becomes unknown | `unitless(...)` |
| W0601 | A selector that matches nothing in its scope; set algebra such as a and b never warns | `residues(250)` is not among this group's 241 particles |  |
| W0602 | `particles(resid=...)` with an id in several chains | `resid=501` matches chains A and B; the merge would hide it | `chain='A'` |
| W0603 | Groups with no values in `by=` | 15 residues have no values in `_z` | a narrower group |
| W0605 | Equal extents that describe different units, paired by position | pairs the 433 residues of chain A with the 433 of chain B by position | common(a, b, by='id') |
| W0607 | An unannotated argument over the same units as a dimension the call maps over | cutoff is float\[433\] over the residues of chain A, and the call maps over them; it is passed whole, so the result is bool\[433\]\[433\] | cutoff: float |
| E0701 | A v1 spelling | `residue(1:4)` is v1; positions are subscripts | `residues[0:4]` |
| E0801 | A construct not supported yet | `while` is not supported yet |  |
| E0802 | Recursion | `walk` calls itself |  |
| E0803 | A construct never supported | classes are not part of VIAL |  |
| R0101 | A position outside a dynamic extent | frame 12: position 5 is outside `[4]` |  |
| R0102 | An integer result outside 32 bits | sum over frames is 3,118,402,117, more than int holds | sum(float(x)) |

```python
# on T
residue(1:4)                                            # error E0701
residues(1..4, 'ALA')                                   # error E0203
angle(particles(resname='GTP'), particles(resname='GDP'))   # error E0204
centre(particles(resname='GTP'))                        # error E0201
distance(residues(1..4, chain='A'), particles(name='MG'))   # error E0302: residues of several particles are not points
if distance(residues(chain='A').com, particles(resname='VLB').com) < 5.0:   # error E0303
    pass
if distance(particles(resname='GTP').com, particles(resname='GDP').com) < 5.0:   # error E0403
    pass
class Foo:                                              # error E0803
    pass

def walk(n):
    return walk(n - 1)                                  # error E0802

m = smiles("O=P(O)(O)O")
x = spread(count(m), m)                                 # error E0306: the phosphates of GTP share bridging oxygens
n = rmsd(positions(particles(name='CA', resid=1..4, chain='A')), positions(particles(name='CB', resid=1..4, chain='A')))   # warning W0605: particle 0: CA against CB
residues = 5                                            # warning W0201: last, so the lines above still see the selector
```

```python
# on P
ca = particles(name='CA')
_near = residues(within(3.0, ca[0]))                    # -> set[?], frame
r = distance(_near[4].com, ca[0])                       # warning W0405; R0101 and a missing value in frames where fewer than five residues are that close
```
