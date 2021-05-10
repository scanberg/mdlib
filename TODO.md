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

### BITFIELD ###
    [ ] Change interface from offset + length to beg + end
    [ ] Fix bit-scan
    [ ] Implement and mitigate to an opaque 'sparse-bitfield' type

### MOLECULE + TRAJECTORY (1 Week) ###
    [X] Revise trajectory interface (implement in mdlib)
        [X] Implement XTC
        [X] Implement PDB
    [/] Molecule format and conversions from native data
        [X] PDB
        [ ] GRO

### Basic Script (2 Weeks) ### 
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
    [ ] Implement fully the FLAG_QUERYABLE for length and validation
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
        [ ] dihedral
        [ ] rmsd
    [ ] Expose simplified function for selection queries.
    [/] IMPLEMENTS TESTSSS!!!


#### Further improvements (much later on) ####
    [ ] Compile time evaluate AST tree, evaluate EVERYTHING that can be evaluated at compile time.
    [ ] Parallelize expression parsing.
    [ ] Resolve identifier references and complete type check for unresolved expressions. (YIELD)
    [ ] Construct expression dependency tree based on references to other identifiers. This is important to see what parts can be evaluated in parallel.

### Scripting Editor (1 Week) ###
    [ ] Error messages in script
    [ ] Syntax highlighting
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
    [ ] Set histogram resolution
    [ ] 

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