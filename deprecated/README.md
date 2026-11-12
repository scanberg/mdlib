# Deprecated

Code kept for reference, which is not compiled and not kept up to date with the rest of mdlib.

- `gfx/`: `md_gfx`, an experimental OpenGL renderer of spacefill with GPU culling of instances and clusters
  (viamd's `EXPERIMENTAL_GFX_API`, removed), and its shaders under `gfx/shaders/gfx/`, which is where `md_gfx.c`
  expects them (`<shaders/gfx/common.h>` with `gfx/` on the include path, and `<shader base dir>/gfx/...` at run time).
