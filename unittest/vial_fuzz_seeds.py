#!/usr/bin/env python3
"""Writes a seed corpus for fuzz_vial: every code block of docs/vial_reference.md, the prelude and the golden inputs.

    python3 unittest/vial_fuzz_seeds.py <out_dir>
"""
import pathlib, re, sys

root = pathlib.Path(__file__).resolve().parent.parent
out = pathlib.Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
ref = (root / 'docs' / 'vial_reference.md').read_text(encoding='utf-8')
blocks = re.findall(r'^```(?:python|vial)[^\n]*\n(.*?)^```', ref, flags=re.S | re.M)
files = [root / 'src' / 'vial' / 'prelude.vial', *sorted((root / 'unittest' / 'vial_golden').glob('*.vial'))]
for i, b in enumerate(blocks):
    (out / f'ref_{i:03}.vial').write_text(b, encoding='utf-8')
for f in files:
    (out / f.name).write_bytes(f.read_bytes())
print(f'{len(blocks)} reference blocks and {len(files)} files in {out}')
