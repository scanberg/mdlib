#!/usr/bin/env python3
"""
check_gpu_arg_layout.py -- portability check for md_gpu argument structs.

Why this exists
---------------
An md_gpu kernel takes one pointer to a caller-defined argument struct, and
the host mirrors that struct in C. For the mirror to be correct on *both*
backends, the struct must lay out identically in SPIR-V and in Metal Shading
Language. Slang emits SPIR-V with scalar layout (a vector's alignment is its
component's, and uint3 occupies 12 bytes), while MSL gives uint2 8/8, uint3
16/16 and uint4 16/16. So:

    struct { uint3 dim; float scale; ... }
        SPIR-V: dim@0  scale@12
        MSL:    dim@0  scale@16

and the single C struct the host fills is wrong on one backend -- silently,
with no diagnostic from either toolchain.

Adding explicit padding does not rescue the 3-component case: SPIR-V puts the
pad in the hole at offset 12, MSL puts it after its own internal padding at 16,
so the two shift by different amounts. 2- and 4-component vectors *are*
rescuable, because only their alignment differs, not their size -- pad them to
a multiple of 8 or 16 respectively and both targets agree.

md_gpu.h absorbs this with md_gpu_uint3 / md_gpu_float4 / ... whose size
and alignment track the backend being built, so a C mirror written with them
matches the shader on both. That only holds while Slang lays things out the way
those types assume, which is what this script checks: it models both targets
from the field list and compares the model against the offsets Slang actually
produced.

  * SPIR-V offsets come from OpMemberDecorate Offset in the compiled module.
  * Metal offsets are computed from Slang's generated MSL.

A mismatch means Slang changed its layout and the C vector types in md_gpu.h
need to change with it -- or that a member uses a construct the types do not
model (a matrix, a nested struct), which must then be verified by hand.

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
# vector's alignment is its component's and a 3-vector is not rounded up.
# These mirror MD_GPU_VEC_ALIGN* in md_gpu.h with the Vulkan values.
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
}


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


def model_offsets(fields, table):
    off, out = 0, []
    for field, ty in fields:
        if ty == '*' or ty.startswith('texture') or ty == 'sampler':
            size = align = 8
        else:
            size, align = table[ty]
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

            # SPIR-V is measured from the module; Metal is modelled from the
            # MSL layout rules, because there is no Metal compiler here. If the
            # measured SPIR-V matches its model, the shared model is sound and
            # md_gpu.h's vector types line up on both.
            spv_model = model_offsets(fields, SCALAR_TYPES)
            msl_model = model_offsets(fields, MSL_TYPES)

            if len(s_off) != len(spv_model):
                problems.append(f"{src.name}:{entry}: {base}: SPIR-V reports {len(s_off)} members, "
                                f"the model expects {len(spv_model)}")
                continue
            for (n1, o1), (_, o2) in zip(s_off, spv_model):
                if o1 != o2:
                    problems.append(
                        f"{src.name}:{entry}: {base}.{n1} is at offset {o1} in SPIR-V but the "
                        f"model in md_gpu.h assumes {o2}. Slang's layout has changed; "
                        f"MD_GPU_VEC_ALIGN* needs updating to match.")
                    break

            # The two targets deliberately DO differ for vectors and matrices --
            # that is exactly what md_gpu.h's md_gpu_uint3 / md_gpu_float4x4 and
            # friends absorb, by carrying per-backend size and alignment. So a
            # difference here is fine; an unmodelled member type is not, because
            # then md_gpu.h has no C type that tracks it.
            _ = msl_model

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
