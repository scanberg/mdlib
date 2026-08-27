#!/usr/bin/env python3
"""
check_gpu_arg_layout.py -- portability check for md_gpu argument structs.

Why this exists
---------------
An md_gpu kernel takes one pointer to a caller-defined argument struct, and the
host mirrors that struct in C. A single C mirror has to be correct for *both*
backends, so the struct must lay out identically in SPIR-V and in Metal Shading
Language. Slang emits SPIR-V with scalar layout -- a vector's alignment is its
component's, and uint3 occupies 12 bytes -- while MSL gives uint2 8/8, uint3
16/16 and uint4 16/16. So:

    struct { uint n; uint4 v; ... }
        SPIR-V: n@0  v@4
        MSL:    n@0  v@16

and the one C struct the host fills is wrong on one backend -- silently, with no
diagnostic from either toolchain.

md_gpu.h answers this with an ABI rule rather than with per-backend types. The
vector *sizes* already agree, so the two layouts agree exactly when every vector
member sits at an offset satisfying the stricter (Metal) alignment:

    vec2, texture/sampler handle    offset must be a multiple of 8
    vec4, float4x4                  offset must be a multiple of 16

(Scalars and 8-byte device pointers place themselves.)

Three-component vectors are not allowed at all. Their sizes differ -- 12 bytes
against 16 -- and no padding reconciles that, because MSL's slack sits *inside*
the vector while SPIR-V's sits after it, so every following member shifts by 4
on one target only. Landing one on a 16-byte boundary fixes where it starts and
nothing else. Use a 4-vector with an unused .w, or scalars.

This script is what makes that a rule rather than a convention. For each entry
point it compiles the shader for both targets and rejects an argument struct if

  1. it contains a 3-component vector, or
  2. the SPIR-V and MSL layouts disagree at any member.

Whatever passes can be mirrored by one C struct built from md_gpu.h's
md_gpu_uint4 / md_gpu_float4x4 / ... -- which carry the required alignment and
no #ifdef, so the header does not need to know which backend it was built for.

  * SPIR-V offsets come from OpMemberDecorate Offset in the compiled module.
  * Metal offsets are computed from Slang's generated MSL.

The measured SPIR-V is also compared against the model, so a future Slang that
changes its layout rules is caught here rather than in a wrong result several
launches downstream.

Usage:
    python3 tools/check_gpu_arg_layout.py --slangc <path> <file.slang> <entry>...
"""

import argparse
import re
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

STRUCT_RE = re.compile(r'^struct\s+(\w+)\s*\{(.*?)^\};', re.M | re.S)

# Metal Shading Language: scalars natural; vec2 8/8; vec3 and vec4 16/16;
# device/constant pointers 8/8; texture and sampler handles 8/8.
MSL_TYPES = {
    'bool': (1, 1), 'char': (1, 1), 'uchar': (1, 1),
    'short': (2, 2), 'ushort': (2, 2), 'half': (2, 2),
    'int': (4, 4), 'uint': (4, 4), 'float': (4, 4),
    'long': (8, 8), 'ulong': (8, 8),
    'int2': (8, 8), 'uint2': (8, 8), 'float2': (8, 8),
    'int3': (16, 16), 'uint3': (16, 16), 'float3': (16, 16),
    'int4': (16, 16), 'uint4': (16, 16), 'float4': (16, 16),
    'half2': (4, 4), 'half3': (8, 8), 'half4': (8, 8),
    # Slang lowers a float4x4 to a wrapper holding array<float4,4>: four
    # 16-byte columns. Measured, not assumed.
    '_MatrixStorage_float4x4': (64, 16),
    # A device pointer and a bindless texture/sampler handle are both plain
    # 8-byte, 8-aligned values here.
    'ptr': (8, 8), 'handle': (8, 8),
}


def run_slang(slangc, source, entry, target, extra, suffix):
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / ("o" + suffix)
        cmd = [slangc, str(source), "-target", target, "-entry", entry,
               "-Wno-40100", "-o", str(out)] + extra
        p = subprocess.run(cmd, capture_output=True, text=True)
        if not out.exists():
            raise SystemExit(f"slangc failed for {source}:{entry} [{target}]\n{p.stderr}")
        return out.read_bytes()


def spirv_struct_offsets(blob):
    """{struct name: [(field, offset), ...]} from OpMemberDecorate Offset."""
    words = struct.unpack('<%dI' % (len(blob) // 4), blob)
    i, ins = 5, []
    while i < len(words):
        op, count = words[i] & 0xFFFF, words[i] >> 16
        if count == 0:
            break
        ins.append((op, list(words[i:i + count])))
        i += count

    def text(ws):
        return b''.join(struct.pack('<I', v) for v in ws).split(b'\0')[0].decode('utf8', 'replace')

    names, members, offsets = {}, {}, {}
    for op, x in ins:
        if op == 5:                                  # OpName
            names[x[1]] = text(x[2:])
        elif op == 6:                                # OpMemberName
            members.setdefault(x[1], {})[x[2]] = text(x[3:])
        elif op == 72 and x[3] == 35:                # OpMemberDecorate Offset
            offsets.setdefault(x[1], {})[x[2]] = x[4]

    out = {}
    for sid, offs in offsets.items():
        name = names.get(sid)
        if not name:
            continue
        out[base_name(name)] = [(members.get(sid, {}).get(k, f'm{k}'), offs[k]) for k in sorted(offs)]
    return out


def msl_structs(text):
    return dict(STRUCT_RE.findall(text))


def msl_arg_struct_names(structs):
    """Pointee types of the Root struct's members -- the argument structs."""
    names = set()
    for _, body in structs.items():
        for m in re.finditer(r'^\s*(\w+)\s+(?:device|constant)\s*\*\s*\w+\s*;', body, re.M):
            names.add(m.group(1))
    return names


def base_name(name):
    """Normalise a struct name to compare across targets.

    Slang decorates the SPIR-V type with the layout it chose ('RootArgs_natural',
    'FillRoot_std430') and the MSL type with a disambiguating index
    ('RootArgs_natural_0'), so the suffixes can stack. Strip until stable."""
    prev = None
    while prev != name:
        prev = name
        name = re.sub(r'_(natural|std430|std140|scalar|packed)$', '', name)
        name = re.sub(r'_\d+$', '', name)
    return name


def normalise_member_type(base):
    """Map Slang's generated matrix-storage wrapper onto a modelled type."""
    if base.startswith('_MatrixStorage_float4x4'):
        return '_MatrixStorage_float4x4'
    return base


# SPIR-V as Slang emits it for a pointer-accessed struct: scalar layout, so a
# vector's alignment is its component's and a 3-vector is not rounded up. The
# 3-vector rows exist so a violation can be reported precisely, not because a
# 3-vector is ever allowed to pass.
SCALAR_TYPES = {
    'bool': (1, 1), 'char': (1, 1), 'uchar': (1, 1),
    'short': (2, 2), 'ushort': (2, 2), 'half': (2, 2),
    'int': (4, 4), 'uint': (4, 4), 'float': (4, 4),
    'long': (8, 8), 'ulong': (8, 8),
    'int2': (8, 4), 'uint2': (8, 4), 'float2': (8, 4),
    'int3': (12, 4), 'uint3': (12, 4), 'float3': (12, 4),
    'int4': (16, 4), 'uint4': (16, 4), 'float4': (16, 4),
    'half2': (4, 2), 'half3': (6, 2), 'half4': (8, 2),
    # Same 64 bytes in SPIR-V, but scalar layout aligns it to its component.
    # So a float4x4 agrees across targets exactly when it sits on a 16-byte
    # offset -- being the first member is the usual way that happens.
    '_MatrixStorage_float4x4': (64, 4),
    # A device pointer is a 64-bit PhysicalStorageBuffer address: 8/8. A
    # bindless handle is NOT -- Slang lowers DescriptorHandle<T> to a pair of
    # 32-bit words, so scalar layout gives it 8/4 where MSL gives 8/8. Like a
    # vec2, it agrees across the targets only on an 8-byte boundary.
    'ptr': (8, 8), 'handle': (8, 4),
}


# Types whose *size* differs between the targets (12 bytes in SPIR-V, 16 in
# MSL). Nothing the surrounding struct can do reconciles that, so they are
# rejected outright rather than modelled.
VEC3_TYPES = {'int3', 'uint3', 'float3', 'half3'}


def msl_fields(body):
    """[(field, msl_type_name), ...] or None if a member is not understood."""
    out = []
    for line in body.splitlines():
        line = line.strip()
        if not line or line.startswith('//'):
            continue
        m = re.match(r'(?:const\s+)?(\w+)(?:<[^>]*>)?\s+(?:device|constant|thread)\s*\*\s*(\w+)\s*;', line)
        if m:
            out.append((m.group(2), '*'))
            continue
        m = re.match(r'(\w+)(<[^>]*>)?\s+(\w+)\s*;', line)
        if not m:
            continue
        base, field = normalise_member_type(m.group(1)), m.group(3)
        if base in MSL_TYPES or base.startswith('texture') or base == 'sampler':
            out.append((field, base))
        else:
            return None
    return out


def plain(name):
    """Strip the disambiguating suffix Slang appends to MSL field names."""
    return re.sub(r'_\d+$', '', name)


def kind(ty):
    """Collapse a member's spelling onto a row in the layout tables."""
    if ty == '*':
        return 'ptr'
    if ty.startswith('texture') or ty == 'sampler':
        return 'handle'
    return ty


def model_offsets(fields, table):
    off, out = 0, []
    for field, ty in fields:
        size, align = table[kind(ty)]
        off = (off + align - 1) // align * align
        out.append((field, off))
        off += size
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slangc", required=True)
    ap.add_argument("--bindless-space", default="0")
    ap.add_argument("source")
    ap.add_argument("entries", nargs="+")
    args = ap.parse_args()

    src = Path(args.source)
    problems, checked = [], 0

    for entry in args.entries:
        spv = run_slang(args.slangc, src, entry, "spirv",
                        ["-emit-spirv-directly", "-profile", "glsl_450",
                         "-bindless-space-index", args.bindless_space], ".spv")
        msl = run_slang(args.slangc, src, entry, "metal", [], ".metal").decode('utf8', 'replace')

        spv_layout = spirv_struct_offsets(spv)
        structs    = msl_structs(msl)
        targets    = msl_arg_struct_names(structs)
        if not targets:
            problems.append(f"{src.name}:{entry}: no argument struct found in the generated MSL")
            continue

        for mname in sorted(targets):
            body = structs.get(mname)
            if body is None:
                continue
            base   = base_name(mname)
            fields = msl_fields(body)
            s_off  = spv_layout.get(base)
            if fields is None:
                problems.append(f"{src.name}:{entry}: {base} has a member type this check does "
                                f"not model; verify its layout by hand")
                continue
            if s_off is None:
                problems.append(f"{src.name}:{entry}: {base} has no SPIR-V member offsets "
                                f"(is it used by this entry point?)")
                continue
            checked += 1

            # Rule 1: no 3-component vectors. Report every offender rather
            # than the first, since the fix is the same for all of them and the
            # author would otherwise rebuild once per member.
            bad3 = [f for f, ty in fields if ty in VEC3_TYPES]
            if bad3:
                for f in bad3:
                    problems.append(
                        f"{src.name}:{entry}: {base}.{plain(f)} is a 3-component vector, which is "
                        f"12 bytes in SPIR-V and 16 in MSL. No padding fixes this -- use a "
                        f"4-vector with an unused .w, or separate scalars.")
                continue

            # SPIR-V is measured from the module; Metal is modelled from the MSL
            # layout rules, because there is no Metal compiler here. Check the
            # model against the measurement first: if the two disagree, Slang has
            # changed its layout and every conclusion below is worthless.
            spv_model = model_offsets(fields, SCALAR_TYPES)
            msl_model = model_offsets(fields, MSL_TYPES)

            if len(s_off) != len(spv_model):
                problems.append(f"{src.name}:{entry}: {base}: SPIR-V reports {len(s_off)} members, "
                                f"the model expects {len(spv_model)}")
                continue
            stale = False
            for (n1, o1), (_, o2) in zip(s_off, spv_model):
                if o1 != o2:
                    problems.append(
                        f"{src.name}:{entry}: {base}.{n1} is at offset {o1} in SPIR-V but this "
                        f"script's model of scalar layout says {o2}. Slang's layout rules have "
                        f"changed; SCALAR_TYPES here (and the ABI rule in md_gpu.h) need "
                        f"revisiting before anything else in this report is trustworthy.")
                    stale = True
                    break
            if stale:
                continue

            # Rule 2: the two targets must agree member for member. Under the
            # ABI rule they can, and where they do not it is always the same
            # cause -- a vector sitting at an offset that does not satisfy its
            # Metal alignment -- so say what to pad and by how much.
            for (n1, o1), (_, o2) in zip(spv_model, msl_model):
                if o1 != o2:
                    need = MSL_TYPES.get(kind(dict(fields).get(n1, '')), (0, 0))[1]
                    hint = (f" Pad the members before it so {plain(n1)} starts on a "
                            f"{need}-byte boundary." if need else "")
                    problems.append(
                        f"{src.name}:{entry}: {base}.{plain(n1)} lands at offset {o1} in "
                        f"SPIR-V but {o2} in MSL, so one C mirror cannot be right for "
                        f"both.{hint} (offsets in {base}: "
                        f"SPIR-V {[o for _, o in spv_model]}, MSL {[o for _, o in msl_model]})")
                    break

    if problems:
        print("md_gpu argument structs are not portable across backends:\n", file=sys.stderr)
        for p in sorted(set(problems)):
            print("  " + p, file=sys.stderr)
        return 1
    if checked == 0:
        print(f"{src.name}: no argument structs found to check", file=sys.stderr)
        return 1

    print(f"md_gpu arg layout OK: {src.name} "
          f"({len(args.entries)} entry points, {checked} argument structs compared field by field)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
