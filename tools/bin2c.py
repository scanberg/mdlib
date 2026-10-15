#!/usr/bin/env python3
"""
bin2c.py -- embed a binary file in a C source file.

    bin2c.py <input> <output.c> <symbol>

The output defines

    const uint8_t <symbol>_start[];
    const size_t  <symbol>_byte_size;

which compile_gpu_shaders() (cmake/CompileGpuShaders.cmake) declares in the
header it generates. md_gpu_build_shader.py imports c_source() and embeds
straight from memory; the command line is for binaries produced by another
tool, such as a .metallib.

This replaces cmake/BinToC.cmake, which walked the file one byte per CMake loop
iteration and copied the whole hex string on every iteration -- quadratic, 11 s
for a 77 KB SPIR-V module against a few ms here.
"""
import sys


def c_source(data, symbol):
    """The C translation unit defining <symbol>_start and <symbol>_byte_size."""
    if not data:
        raise ValueError(f"bin2c: {symbol}: nothing to embed (empty input)")
    rows = ("    " + ",".join("0x%02x" % b for b in data[i:i + 16]) + ","
            for i in range(0, len(data), 16))
    return ("#include <stdint.h>\n"
            "#include <stddef.h>\n"
            "\n"
            f"const uint8_t {symbol}_start[] = {{\n"
            + "\n".join(rows) +
            "\n};\n"
            "\n"
            f"const size_t {symbol}_byte_size = sizeof({symbol}_start);\n")


def main(argv):
    if len(argv) != 4:
        sys.stderr.write("usage: bin2c.py <input> <output.c> <symbol>\n")
        return 2
    with open(argv[1], "rb") as f:
        data = f.read()
    try:
        text = c_source(data, argv[3])
    except ValueError as e:
        sys.stderr.write(f"{e}\n")
        return 1
    with open(argv[2], "w", newline="\n") as f:
        f.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
