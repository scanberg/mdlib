#!/usr/bin/env python3
"""
spirv_alias.py -- decorate the variables that hold buffer pointers AliasedPointer.

Why this exists
---------------
SPV_KHR_physical_storage_buffer: a variable whose pointee is a pointer, or an
array of pointers, in the PhysicalStorageBuffer storage class must be decorated
with exactly one of AliasedPointer or RestrictPointer. Slang (2026.8.1) leaves
the Function variables it makes for such locals undecorated: a kernel's
`uint* o = a.out;`, and temporaries of its own. The Vulkan validation layer
reports every one (spirv-val of SPIRV-Tools 2025.1 does not check it yet).
Drivers accept the modules, but they are not valid SPIR-V.

AliasedPointer claims nothing: the pointer may alias, as any pointer may.
RestrictPointer would be a promise nobody checked.

The decorations are added at the end of the annotation section, one per
variable that has neither; nothing else in the module changes.

    spirv_alias.py in.spv out.spv      (in and out may be the same file)
"""
import struct
import sys

OP_TYPE_ARRAY, OP_TYPE_RUNTIME_ARRAY, OP_TYPE_POINTER, OP_VARIABLE, OP_DECORATE = 28, 29, 32, 59, 71
STORAGE_FUNCTION, STORAGE_PHYSICAL_STORAGE_BUFFER = 7, 5349
DECORATION_RESTRICT_POINTER, DECORATION_ALIASED_POINTER = 5355, 5356

# What follows the annotation section: types, constants, global variables, OpLine and OpNoLine, the module-scope
# OpExtInst of debug information, and the functions. An annotation never comes after the first of these.
AFTER_ANNOTATIONS = set(range(19, 40)) | {1, 8, 12, 41, 42, 43, 44, 46, 48, 49, 50, 51, 52, 54, 59, 317, 322,
                                         327, 4433, 4456, 5341, 5358}


def instructions(words):
    """(index, opcode, operands) of each instruction after the header."""
    i = 5
    while i < len(words):
        op, n = words[i] & 0xFFFF, words[i] >> 16
        if n == 0 or i + n > len(words):
            raise ValueError("spirv_alias: malformed instruction at word %d" % i)
        yield i, op, list(words[i + 1:i + n])
        i += n


def undecorated(words):
    """The ids of the variables that hold buffer pointers and carry neither decoration, and the word index at which
    the annotation section ends."""
    psb_pointers = set()     # pointer types in the PhysicalStorageBuffer storage class
    elements = {}            # array type -> element type
    pointees = {}            # pointer type -> pointee type (other storage classes)
    variables = []           # (result type, id) of Function variables
    decorated = set()
    end = None
    for i, op, a in instructions(words):
        if end is None and op in AFTER_ANNOTATIONS:
            end = i
        if op == OP_TYPE_POINTER:
            if a[1] == STORAGE_PHYSICAL_STORAGE_BUFFER:
                psb_pointers.add(a[0])
            else:
                pointees[a[0]] = a[2]
        elif op in (OP_TYPE_ARRAY, OP_TYPE_RUNTIME_ARRAY):
            elements[a[0]] = a[1]
        elif op == OP_VARIABLE and a[2] == STORAGE_FUNCTION:
            variables.append((a[0], a[1]))
        elif op == OP_DECORATE and a[1] in (DECORATION_RESTRICT_POINTER, DECORATION_ALIASED_POINTER):
            decorated.add(a[0])

    def holds_pointer(t):
        while t in elements:
            t = elements[t]
        return t in psb_pointers

    ids = [v for t, v in variables if v not in decorated and holds_pointer(pointees.get(t))]
    return ids, end if end is not None else len(words)


def decorate(words):
    """The module with an OpDecorate AliasedPointer for each such variable, and how many were added."""
    ids, end = undecorated(words)
    added = []
    for v in ids:
        added += [(3 << 16) | OP_DECORATE, v, DECORATION_ALIASED_POINTER]
    out = list(words[:end]) + added + list(words[end:])
    left, _ = undecorated(out)
    if left:
        raise ValueError("spirv_alias: variables %s still undecorated" % left)
    return out, len(ids)


def main():
    if len(sys.argv) != 3:
        sys.stderr.write(__doc__)
        return 2
    data = open(sys.argv[1], "rb").read()
    out, n = decorate(struct.unpack("<%dI" % (len(data) // 4), data))
    open(sys.argv[2], "wb").write(struct.pack("<%dI" % len(out), *out))
    print("spirv_alias: %d variables decorated" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
