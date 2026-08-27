#!/usr/bin/env python3
"""
gen_md_gpu_builtin.py -- regenerate src/core/md_gpu_builtin_spv.inl.

The Vulkan backend needs the built-in md_gpu_make_grid kernel available before
any user shader is compiled, and before CMake's shader step has run: it is
created during device creation. So its SPIR-V is checked in as a C array rather
than built from source, and this script is what produces that array.

Run it after editing src/shaders/gpu/md_gpu_make_grid.slang (or the md_gpu.slang
prelude, which it includes):

    python3 tools/gen_md_gpu_builtin.py --slangc <path-to-slangc>

The flags below must stay in step with cmake/CompileGpuShaders.cmake -- the
bindless space index in particular, or the descriptors silently stop matching
the set layout. Metal has no equivalent step: its built-in is the hand-written
MSL string in src/core/md_gpu_builtin_msl.inl, which has to be edited by hand
to match whenever the argument struct changes.
"""

import argparse
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT    = Path(__file__).resolve().parent.parent
SOURCE  = ROOT / "src/shaders/gpu/md_gpu_make_grid.slang"
OUTPUT  = ROOT / "src/core/md_gpu_builtin_spv.inl"
ENTRY   = "main"
SYMBOL  = "md_gpu_make_grid_spv"
BINDLESS_SPACE = "0"          # == MD_GPU_BINDLESS_SPACE in CompileGpuShaders.cmake
WORDS_PER_LINE = 8


def compile_spirv(slangc, source, entry):
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "out.spv"
        cmd = [slangc, str(source),
               "-Wno-40100",
               "-target", "spirv",
               "-emit-spirv-directly",
               "-profile", "glsl_450",
               "-bindless-space-index", BINDLESS_SPACE,
               "-entry", entry,
               "-o", str(out)]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if not out.exists():
            raise SystemExit(f"slangc failed for {source}:{entry}\n{p.stderr}")
        return out.read_bytes()


def slangc_version(slangc):
    p = subprocess.run([slangc, "-v"], capture_output=True, text=True)
    return (p.stdout + p.stderr).strip().splitlines()[0] if (p.stdout or p.stderr) else "unknown"


def emit(blob, version):
    if len(blob) % 4:
        raise SystemExit("SPIR-V module is not a whole number of words")
    words = struct.unpack("<%dI" % (len(blob) // 4), blob)
    lines = []
    for i in range(0, len(words), WORDS_PER_LINE):
        lines.append("    " + " ".join("0x%08xu," % w for w in words[i:i + WORDS_PER_LINE]))
    rel = SOURCE.relative_to(ROOT)
    return (f"/* Generated from {rel} by slangc {version}.\n"
            f"   Regenerate with tools/gen_md_gpu_builtin.py -- do not edit by hand. */\n"
            f"\n"
            f"static const uint32_t {SYMBOL}[] = {{\n"
            + "\n".join(lines) + "\n};\n\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slangc", required=True, help="path to the slangc executable")
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit non-zero if the checked-in file is stale")
    a = ap.parse_args()

    text = emit(compile_spirv(a.slangc, SOURCE, ENTRY), slangc_version(a.slangc))

    if a.check:
        current = OUTPUT.read_text() if OUTPUT.exists() else ""
        if current != text:
            print(f"{OUTPUT.relative_to(ROOT)} is out of date; "
                  f"rerun tools/gen_md_gpu_builtin.py", file=sys.stderr)
            return 1
        print(f"{OUTPUT.relative_to(ROOT)} is up to date")
        return 0

    OUTPUT.write_text(text)
    print(f"wrote {OUTPUT.relative_to(ROOT)} ({len(text.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
