#!/usr/bin/env python3
"""
md_gpu_build_shader.py -- build one entry point of an md_gpu shader.

compile_gpu_shaders() (cmake/CompileGpuShaders.cmake) runs this once per entry
point, and the one run does everything the entry point needs:

  1. slangc compiles the entry point for SPIR-V and MSL in a single process.
     The front end -- parsing, checking, the md_gpu.slang prelude -- runs once,
     and slangc's fixed start-up cost (~0.4 s) is paid once rather than once
     per target. The output is byte-identical to two separate slangc runs.
  2. check_gpu_arg_layout.check_entry() compares the argument-struct layouts of
     the two targets; a mismatch is a build error (see that script for why).
     Every descriptor the SPIR-V declares must sit in the bindless set, where
     md_gpu's heap is (check_descriptor_sets).
  3. On Vulkan, spirv_deforward reorders the SPIR-V type section for drivers
     that reject forward-declared buffer pointers (see that script).
  4. The kernel or shader descriptor -- group size and argument-struct size --
     is written to --desc.
  5. With --embed, the backend's binary is written as a C array: SPIR-V on
     Vulkan, the MSL itself on Metal without an offline compiler. (With one,
     compile_gpu_shaders() builds a metallib from --msl and embeds that.)

Why a single step: every entry point becomes one build command that no other
command in the target waits on, so any generator can run them all at once --
including MSBuild, whose parallel custom build does not order items by their
inputs and outputs.

Nothing is written unless every step succeeds. A failed check leaves the
previous outputs older than their inputs, so the step reruns on the next build.

Usage:
    md_gpu_build_shader.py --slangc <path> --backend VULKAN|METAL
        --stage compute|vertex|fragment --namespace <ns>
        --spv <out.spv> --msl <out.metal> --desc <out.inl> [--embed <out.c>]
        [--bindless-space <n>] [--precise-math] [--debug-info] <file.slang> <entry>
"""

import argparse
import os
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

# The modules below live in the source tree; do not litter it with __pycache__.
sys.dont_write_bytecode = True
import bin2c                                    # noqa: E402
import check_gpu_arg_layout as layout           # noqa: E402
import spirv_deforward                          # noqa: E402


def compile_both(slangc, source, entry, bindless_space, debug_info=False):
    """(spirv bytes, msl bytes) from one slangc run, or exit on failure. With
    debug_info, slangc -g: source files, lines and names for shader profilers
    (SPIR-V: NonSemantic.Shader.DebugInfo.100)."""
    with tempfile.TemporaryDirectory() as td:
        spv, msl = Path(td) / "o.spv", Path(td) / "o.metal"
        cmd = [slangc, str(source), "-entry", entry,
               # Metal reserves 'main', so Slang renames entry points.
               "-Wno-40100",
               # The bindless space is a SPIR-V notion; the Metal target in the
               # same run warns that it cannot use it. Its MSL is unaffected.
               "-Wno-39012",
               "-target", "spirv", "-emit-spirv-directly", "-profile", "glsl_450",
               # md_gpu's subgroup baseline. Declared so Slang warns (E41012) only when a
               # kernel needs more than this.
               "-capability", "subgroup_basic+subgroup_ballot+subgroup_arithmetic",
               "-bindless-space-index", bindless_space, "-o", str(spv),
               "-target", "metal", "-o", str(msl)]
        if debug_info:
            cmd.insert(2, "-g")
        p = subprocess.run(cmd, capture_output=True, encoding="utf-8", errors="replace")
        # Pass slangc's own diagnostics through, warnings included.
        sys.stdout.write(p.stdout)
        sys.stderr.write(p.stderr)
        if p.returncode != 0 or not spv.exists() or not msl.exists():
            sys.stderr.write(f"slangc failed for {source}:{entry}\n")
            sys.exit(1)
        return spv.read_bytes(), msl.read_bytes()


def check_descriptor_sets(spv, bindless_space):
    """The DescriptorSet of every descriptor in the SPIR-V, if one is not the
    bindless space -- a problem md_gpu would otherwise only meet at pipeline
    creation (the heap missing from the layout). Slang moves its heap off the
    requested space when the file declares [vk::constant_id] constants (it
    warns E39012, which the build silences for the Metal target's sake), so
    that is what the message points at."""
    words = struct.unpack("<%dI" % (len(spv) // 4), spv)
    sets = set()
    i = 5
    while i < len(words):
        op, n = words[i] & 0xFFFF, words[i] >> 16
        if n == 0:
            break
        if op == 71 and n >= 4 and words[i + 2] == 34:      # OpDecorate <id> DescriptorSet <set>
            sets.add(words[i + 3])
        i += n
    wrong = sorted(x for x in sets if x != int(bindless_space))
    if not wrong:
        return None
    return (f"descriptors in set {', '.join(map(str, wrong))}, not the bindless set {bindless_space}; "
            "declare specialization constants with [SpecializationConstant], not [vk::constant_id] "
            "(which moves Slang's bindless heap)")


def write_atomic(path, data):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser(description="Build one entry point of an md_gpu shader.")
    ap.add_argument("--slangc", required=True)
    ap.add_argument("--bindless-space", default="0")
    ap.add_argument("--backend", required=True, choices=["VULKAN", "METAL"])
    ap.add_argument("--stage", required=True, choices=["compute", "vertex", "fragment"])
    ap.add_argument("--namespace", required=True)
    ap.add_argument("--spv", required=True, help="SPIR-V output (deforwarded on Vulkan)")
    ap.add_argument("--msl", required=True, help="MSL output")
    ap.add_argument("--desc", required=True, help="descriptor header output")
    ap.add_argument("--embed", help="write the backend binary as a C array to this file")
    ap.add_argument("--precise-math", action="store_true",
                    help="the descriptor asks for IEEE floating point (Metal: no fast math)")
    ap.add_argument("--debug-info", action="store_true",
                    help="compile with debug information (slangc -g), for shader profilers")
    ap.add_argument("source")
    ap.add_argument("entry")
    args = ap.parse_args()

    src, entry = Path(args.source), args.entry
    spv, msl = compile_both(args.slangc, src, entry, args.bindless_space, args.debug_info)

    problems, group_size, args_size = layout.check_entry(
        src.name, entry, args.stage != "compute", spv, msl.decode("utf8", "replace"))
    if problems:
        layout.report(problems)
        return 1
    problem = check_descriptor_sets(spv, args.bindless_space)
    if problem:
        sys.stderr.write(f"{src.name}:{entry}: {problem}\n")
        return 1

    if args.backend == "VULKAN":
        try:
            words = struct.unpack("<%dI" % (len(spv) // 4), spv)
            out, _, _ = spirv_deforward.deforward(words)
        except ValueError as e:
            sys.stderr.write(f"{src.name}:{entry}: {e}\n")
            return 1
        spv = struct.pack("<%dI" % len(out), *out)

    # Must match the symbol compile_gpu_shaders() declares in its header.
    symbol = f"{args.namespace}_{src.stem}_{entry}"
    outputs = [
        (args.spv, spv),
        (args.msl, msl),
        (args.desc, layout.descriptor_source(
            args.namespace, src.stem, entry, group_size, args_size,
            f"tools/md_gpu_build_shader.py from {src.name} [{entry}]",
            args.precise_math).encode()),
    ]
    if args.embed:
        blob = spv if args.backend == "VULKAN" else msl
        outputs.append((args.embed, bin2c.c_source(blob, symbol).encode()))
    for path, data in outputs:
        write_atomic(path, data)
    return 0


if __name__ == "__main__":
    sys.exit(main())
