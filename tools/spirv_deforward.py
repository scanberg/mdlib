#!/usr/bin/env python3
"""
spirv_deforward.py -- drop OpTypeForwardPointer where no pointer type is recursive.

Why this exists
---------------
Slang declares every physical-storage-buffer pointer to a struct with
OpTypeForwardPointer, also when the pointee is not recursive. That is valid
SPIR-V, but the Intel Windows driver for Gen9 (HD 530, driver 31.0.101.2111,
Vulkan 1.3.215) cannot compile a shader in which such a forward-declared
pointer type is the result of an OpAccessChain or OpPtrAccessChain:
vkCreateComputePipelines returns VK_SUCCESS and a VK_NULL_HANDLE pipeline.
In the GTO kernels that is every read of a float4x4 argument (Slang wraps it
in a _MatrixStorage struct) and every args.pgto[i]. The same modules compile
once the declarations are in dependency order and the forward pointers are
gone, which is what this script produces.

The type/constant/global section is re-emitted in dependency order (stable:
unrelated declarations keep their relative order), so a pointer type is
defined before the struct that holds it. A forward declaration is kept only
for a pointer that closes a genuine cycle (a recursive struct). OpLine and
OpNoLine inside that section are dropped; nothing else in the module changes.

    spirv_deforward.py in.spv out.spv      (in and out may be the same file)
"""
import struct
import sys

# Declarations whose result id is the first operand (types).
TYPE_RESULT_FIRST = {19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
                     35, 36, 37, 38, 322, 327, 4456, 5341, 5358}
# Everything that may appear in the type/constant/global-variable section.
SECTION = TYPE_RESULT_FIRST | {1, 39, 41, 42, 43, 44, 46, 48, 49, 50, 51, 52, 59, 8, 317}
OP_LINE, OP_NOLINE, OP_FUNCTION, OP_TYPE_POINTER, OP_TYPE_FORWARD_POINTER = 8, 317, 54, 32, 39


def id_operands(op, a):
    """Ids a declaration refers to (its own result id excluded)."""
    if op in (19, 20, 21, 22, 26, 31, 34, 35, 36, 37, 38): return []
    if op in (23, 24, 29): return [a[1]]                 # Vector, Matrix, RuntimeArray
    if op in (25, 27): return [a[1]]                     # Image (sampled type), SampledImage
    if op == 28: return [a[1], a[2]]                     # Array: element, length
    if op in (30, 33): return a[1:]                      # Struct members; Function return + params
    if op == 32: return [a[2]]                           # Pointer: pointee
    if op == 4456: return a[1:]                          # CooperativeMatrixKHR: type, scope/rows/cols/use ids
    if op in (1, 41, 42, 43, 46, 48, 49, 50): return [a[0]]
    if op in (44, 51): return [a[0]] + a[2:]             # (Spec)ConstantComposite
    if op == 52: return [a[0]] + a[3:]                   # SpecConstantOp
    if op == 59: return [a[0]] + a[3:4]                  # Variable (+ initializer)
    raise ValueError("spirv_deforward: unhandled opcode %d in the type section" % op)


def result_id(op, a):
    return a[0] if op in TYPE_RESULT_FIRST else a[1]


def deforward(words):
    """Returns (new words, forward pointers before, forward pointers after)."""
    if len(words) < 5 or words[0] != 0x07230203:
        raise ValueError("spirv_deforward: not a SPIR-V module")
    hdr, ins, i = list(words[:5]), [], 5
    while i < len(words):
        n = words[i] >> 16
        if n == 0 or i + n > len(words):
            raise ValueError("spirv_deforward: malformed instruction stream")
        ins.append((words[i] & 0xffff, list(words[i + 1:i + n])))
        i += n

    first = next((k for k, (op, _) in enumerate(ins) if op in SECTION and op not in (OP_LINE, OP_NOLINE)), None)
    if first is None:
        return list(words), 0, 0
    last = next((k for k, (op, _) in enumerate(ins) if op == OP_FUNCTION), len(ins))
    sec = [(op, a) for op, a in ins[first:last] if op not in (OP_LINE, OP_NOLINE)]
    num_fwd = sum(1 for op, _ in sec if op == OP_TYPE_FORWARD_POINTER)
    decls = [(op, a) for op, a in sec if op != OP_TYPE_FORWARD_POINTER]
    by_id = {result_id(op, a): k for k, (op, a) in enumerate(decls)}

    state, order, need_fwd = {}, [], []
    def visit(k):
        stack = [(k, iter(id_operands(*decls[k])))]
        state[k] = 1
        while stack:
            node, deps = stack[-1]
            for dep in deps:
                j = by_id.get(dep)
                if j is None:
                    continue
                if state.get(j) == 1:                    # back edge: only a pointer may close a cycle
                    if decls[j][0] != OP_TYPE_POINTER:
                        raise ValueError("spirv_deforward: non-pointer type cycle")
                    if dep not in need_fwd:
                        need_fwd.append(dep)
                    continue
                if j not in state:
                    state[j] = 1
                    stack.append((j, iter(id_operands(*decls[j]))))
                    break
            else:
                stack.pop()
                state[node] = 2
                order.append(node)
    for k in range(len(decls)):
        if k not in state:
            visit(k)

    new_sec = [(OP_TYPE_FORWARD_POINTER, [p, decls[by_id[p]][1][1]]) for p in need_fwd]
    new_sec += [decls[k] for k in order]
    out = hdr
    for op, a in ins[:first] + new_sec + ins[last:]:
        out += [((len(a) + 1) << 16) | op] + a
    return out, num_fwd, len(need_fwd)


def main(argv):
    if len(argv) != 3:
        sys.stderr.write("usage: spirv_deforward.py in.spv out.spv\n")
        return 2
    data = open(argv[1], "rb").read()
    words = struct.unpack("<%dI" % (len(data) // 4), data)
    out, _, _ = deforward(words)
    with open(argv[2], "wb") as f:
        f.write(struct.pack("<%dI" % len(out), *out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
