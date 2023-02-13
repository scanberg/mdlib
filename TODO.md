### CORE ###
    [X] Allocator interface
    [X] Default allocator (system allocator)
    [X] Default temp allocator (Thread-local, Ringbuffer)
        [ ] Expose raw interface...
        [ ] Examine potential performance benefits by using direct calls...
    [X] Pool allocator
        [ ] Expose raw interface to avoid indirect function calls
    [X] Arena allocator
        [ ] Expose raw interface to avoid indirect function calls
    [ ] Stack allocator (Could be very useful for evaluating trees, where you can pop once a branch is evaluated)
    [ ] String builder
    [X] Tracking allocator for tracking allocations making sure things are freed properly

### BITFIELD / BITOP ###
    [ ] Core API changes
        [ ] Change interface from offset + length to beg + end
        [ ] Fix bit-scan
    [X] Implement and mitigate to an opaque 'semisparse-bitfield' type

### MOLECULE + TRAJECTORY (1 Week) ###
    [X] Revise trajectory interface (implement in mdlib)
        [X] Implement XTC
        [X] Implement PDB
    [X] Molecule format and conversions from native data
        [X] PDB
        [X] GRO
    [ ] Support negative speed for playback
    [ ] Implement Frame offset caching for XTC (and PDB)

### Basic Script (--- Weeks) ### 
    [X] Tokenizer
    [X] Function calls
    [X] Numerical Constants (E, PI, TAU)
    [X] Parse expressions sequentially -> AST for each expression.
    [X] Fix bug where tokenizer never stops if script ends with whitespace.
    [X] Support array types in arithmetic operators
    [X] Resolve ambiguities within the syntax regarding expressions within a local scope vs. global
    [X] Fix the line bug in tokenizer
    [X] Implement FLAG_SYMMETRIC_ARGS to mark procedures where arguments are symmetric
    [X] Implement FLAG_RET_TYPE_EQUAL_LENGTH to mark procedures that return a length equivalent to that of the input
    [X] Implement FLAG_ARG_TYPES_EQUAL_LENGTH to mark procedures where the input arguments should match in length
    [X] Modify compatible_type() to support matching of float[4][1] to float[4].
    [X] Finalize static check for context nodes (determine the type and size)
    [X] Evaluate syntax tree
    [X] Implement selections
        [X] all
        [X] type/name/label
        [X] element(str/irng/int)
        [X] resname
        [X] resid
        [X] residue
        [X] chain(str/irng/int)
        [X] x, y, z
        [ ] within
            [ ] Implement spatial hashing as acceleration structure
    [X] Implement implicit conversion to float3 via position functions
        [X] Int (Atom index)
        [X] Irange ()
        [X] Bitfield -> atoms
    [X] Implement fully the FLAG_QUERYABLE for length and validation
    [ ] Implement full subset support in array [] operator.
    [/] Implement Geometric operations
        [X] com
        [ ] plane

    [/] Implement Compute properties
        [X] distance
        [X] distance_min
        [X] distance_max
        [X] distance_pair
        [X] angle
        [X] dihedral
        [ ] rmsd
        [X] rdf
    [X] Expose simplified function for selection queries.
    [/] IMPLEMENTS TESTSSS!!!
    [X] Implement proper context within evaluation to provide meaningful errors at compile time


#### Further improvements (much later on) ####
    [ ] Compile time evaluate AST tree, evaluate EVERYTHING that can be evaluated at compile time.
    [ ] Parallelize expression parsing.
    [ ] Resolve identifier references and complete type check for unresolved expressions. (YIELD)
    [ ] Construct expression dependency tree based on references to other identifiers. This is important to see what parts can be evaluated in parallel.

### Scripting Editor (1 Week) ###
    [X] Error messages in script
    [/] Syntax highlighting
    [ ] Context (Molecule structure) aware suggestions
        [ ] Map all unique residue names
        [ ] Map all unique atom labels == types
        [ ] Map all unique chains

### Timeline touch up (1 Week) ###
    [ ] Always show current time in the back (vertical marker)
    [ ] Show small time stamps
    [ ] Maintain a fixed point to pixel ratio for line plot. (Performance)

    [ ] Add timestamp with annotations

### Distributions (1 Week) ###
    [ ] Properly display periodic histograms (control which range is shown) (Possibly sync this offset among all distributions which are the same type)
    [ ] Apply adjustable kernel filter when rendering histogram (Do render as histogram to be transparent to the user in what has actually been computed)
    [ ] Implement fixed increments on x-axis based on type (15 degree for angles)

### Volume Rendering (Spatial Distribution) (1 Week) ###
    [ ] Create a dedicated control window for volume rendering
    [ ] Resolution control
    [ ] Spatial-extent control (We may only be interested in a subregion)
    [ ] Expose a Separate window where one can view Spatial Distribution Functions in their native reference frames
        [ ] Option: Show average structure as reference
        [ ] Option: Show all structures superimposed (from the current timestep) as reference (Like we did in the article).

### Property Export (1 Week) ###
    [ ] GUI

### Add DCD support (1 Week) ###

### Graphics (1 Week) ###
    [ ] Add coordinate system widget, user needs to be aware of the orientation of the coordinate system
    [ ] Possibly exchange TAA for SMAA ??? (Will simplify alot of things)

### BUG FIXING AND TWEAKS (3+ Weeks) ###

    [X] Fix bug with unbounded chain and residue as context
    [X] Fix bug where evaluation of subexpressions within large contexts (100+) runs out of temp memory and starts eating its own tail causing corruption.
        Examle: atom(1:10) and element('O') in residue(:)
        atom and element are both bitfields and will be allocated (in full!) for each residue during evaluation
    [X] Fix bug with mold_draw init secondary structures writing or reading out of range. 
    [ ] Why are functions which only take single values matched with arrays? example distance(residue(1), residue(2)) should never have been matched as a valid function since both arguments types are arrays of vec3...

    [ ] Resolve the issue of some some node that are unable to get their types fully resolved since they are the result of some function call which produces a variable length.
        - Propagate flags
        - Yield the type checking if the flag is found on a node.
        - Assert on evaluation that if the length is unknown (-1) then make sure that the flag is set (otherwise a bug).


### Hydrogen Bonds ###

Notes:


Compute Hydrogen bond order

Molecular Geometry

Coordination number
   (1)- Terminal
    2 - Linear, Bent
    3 - Trigonal Planar, Trogonal Pyramidal, T-shaped
    4 - Tetrahedral, Square planar, Seesaw
    5 - Trigonal bipyramidal, Square pyramidal, Pentagonal Planar
    6 - Octahedral, Trigonal Prismatic, Pentagonal Pyramidal
    7 - Pentagonal Bipyramidal, Capped Octahedral, Capped Trigonal Prismatic
    8 - Square antiprismatic, Dodecahedral, Bicapped Trigonal Prismatic,
    9 - Tricapped Trigonal Prismatic, Capped Square Antiprismatic

Tetrahedral - 4 connected components (substituents) to a central one
    When all connected components are of the same type, the bond angles are acos(-1/3) = 109.4712206

