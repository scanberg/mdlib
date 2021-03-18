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

### MOLECULE + TRAJECTORY ###
    [ ] Revise trajectory interface (implement in mdlib)

### Basic Script (2 Weeks) ### 
    [X] Tokenizer
    [X] Function calls
    [X] Numerical Constants (E, PI), are these really needed? The user can just trivially define these whenever they are needed.
    [X] Parse expressions sequentially -> AST for each expression.
    [X] Fix bug where tokenizer never stops if script ends with whitespace.
    [X] Support array types in arithmetic operators
    [/] Resolve ambiguities within the syntax regarding expressions within a local scope vs. global, usage of keyword in
        Partly solved by using keyword in to supply an explicit context for the operation.
    [X] Fix the line bug in tokenizer
    [ ] Compile time evaluate AST tree, evaluate EVERYTHING that can be evaluated at compile time. Propagate evaluated values back up.
    [ ] Evaluate syntax tree
    [ ] Compute selections
    [ ] Compute properties

#### Further improvements (much later on) ####
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