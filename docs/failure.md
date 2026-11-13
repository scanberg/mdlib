# Failure in mdlib

A procedure returns `bool` when its caller can do something about a failure, and `void` when it cannot.

**`bool`**: what happens depends on what the procedure is given at runtime, and the caller has to know.

- Reading and writing files, and parsing: a file can be missing, of another format, or cut short.
- Input that can be wrong at runtime: a selection expression, a file a user picked, a basis that does not fit the
  orbitals.
- Acquiring something: memory, a GPU buffer, an OS object (`md_mutex_init`).
- A question whose answer is yes or no (`md_mutex_try_lock`, `md_unitcell_is_triclinic`).

On `false` the output is in a defined state, empty or as it was, and the reason is in the log unless it is the
caller's to tell.

**`void`**: everything else.

- What fails only when the caller is wrong (NULL, an index out of range, a mutex unlocked twice) asserts.
  `md_mutex_lock`, `md_semaphore_release`.
- A computation that has nothing to do does nothing: no cell, nothing to wrap (`md_geom_pbc`). One that cannot do all
  of its work does what it can and logs the rest (`md_gl_draw` of a representation it does not know).

A `bool` that nobody can act on is noise: every caller has to decide whether to check it, and most do not, so the ones
that matter are lost among them. `md_file_close` stays `bool`: closing a written file is where some file systems report
that the writes did not make it.
