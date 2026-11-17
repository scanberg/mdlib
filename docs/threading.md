# Threads in mdlib

mdlib owns no threads and never waits on any. It provides the building blocks; the application composes them, and
decides how they run in parallel. There is no thread pool in mdlib, and no procedure reaches for a global one.

Where a procedure's work runs in parallel is decided by its shape, in this order of preference:

1. **The caller parallelizes around the call.** Cheap, independent work: per frame, per structure, per residue. The
   procedure is serial and knows nothing of threads. This is most of mdlib.

2. **The procedure is split into prepare / run a range / finish**, when the caller must decide *when* pieces run:
   spread over several frames, shown as they come, or interleaved with other work. `md_gto_iso_cache_info` with
   `md_gto_iso_cache_evaluate(cache, beg, count)` is the model.

3. **The procedure takes an `md_executor_t`** (`core/md_executor.h`), when its parallelism is inside the algorithm:
   phases with serial merges between them, where exposing the phases would leak its internals into every caller.
   `md_topo_compute_extremum_graph_gto` and `md_vial_eval` work this way. A NULL executor runs the loops on the
   calling thread, and is always correct, so the procedure stays usable (and testable) with no threads at all.

A procedure never creates threads of its own (`md_thread_create`) to get its work done. One that does is a sign the
right one of the three above is missing.

## Executors through generic interfaces

A procedure reached through a generic interface gets its executor from the caller of that interface:

- **Virtual attributes**: the extraction context (`md_attribute_io_t`) carries one. A caller that wants a computed
  attribute built in parallel extracts it through `md_attribute_io_create(executor, alloc)` and the `_io` extract
  functions; a provider asks `md_attribute_io_executor(io)`, and hands `io` on to whatever it extracts in turn. A
  plain extract has no context: the provider then runs on the calling thread.
- **Loaders** take it in their `md_load_desc_t` (`md_io.h`: `md_gro_system_init_from_file(sys, state, path, &(md_load_desc_t){ .executor = ex })`). The text loaders (`md_gro`,
  `md_pdb`) read a file a block at a time and each block in pieces at line boundaries, one piece per loop index
  (`core/md_reader.h`: blocks and pieces): a piece counts its lines, which places its atoms, then reads them where they go.

## The executor contract

`md_parallel_for(ex, count, grain, body, arg)` calls `body(arg, beg, end, thread_idx)` over `[0, count)`:

- in half open pieces of at least `grain` (but the last), in any order, on up to `num_threads` threads at once;
- `thread_idx` is in `[0, md_executor_num_threads(ex))`, and no two pieces running at the same time share one, so a
  procedure may keep scratch and accumulators per `thread_idx`, allocated before the loop and merged after it;
- it returns false if the executor stopped handing out pieces (the caller's work was cancelled). Some of the range
  was then never called, and the procedure must treat its results as incomplete.

What a procedure must keep in return:

- **One level deep.** A body does not issue parallel work of its own (no nested `md_parallel_for`) and does not wait
  for any. A thread which waits may run other pieces meanwhile, of the same loop too, and they would share the
  waiting piece's `thread_idx`. A procedure's phases run one after the other.
- **Reproducible where it matters.** Results do not depend on which piece ran where: outputs are addressed by index,
  and per thread accumulators are merged in a fixed order, or are exact (integers).

## Wide sections (experiment)

A procedure of the third kind may also be written *wide*: one function which every thread of a group runs, the
multi-core-by-default (SPMD) model, rather than a sequence of loops handed to the executor one at a time.
`md_run_lanes(ex, max_lanes, fn, arg)` runs `fn(lane, arg)` on up to `max_lanes` lanes; inside it, work is split with
`MD_LANE_FOR(lane, count, grain, r)` (pieces claimed from a shared counter, a sync at the end), serial parts are narrowed
to a lane (`md_lane_idx(lane) == 0`), and lanes meet at `md_lane_sync` and exchange values with `md_lane_bcast_u64`.
Locals live on each lane's stack, and the procedure's control flow reads top to bottom: no phase is a function of its
own, and no context struct is passed between them.

The lanes are the pieces of one `md_parallel_for` on the executor, so any executor runs them. Membership is elastic:
a piece joins when it starts, until the group's first sync closes it, and one which starts later returns at once. A
lane therefore only waits for lanes which have started, and an executor which runs its pieces one after another runs
the group with a single lane. What follows from it, for the one writing a wide function:

- **Open with a loop.** Lanes which start late join only before the first sync. Do setup before going wide (on the
  calling thread), not on lane 0 behind a broadcast at the top, which would close the group before the executor's
  threads have woken.
- **Every lane reaches the same syncs.** A sync inside a branch not every lane takes waits forever; every sync records
  its line, and a group which meets at different ones stops with an assertion instead. A body of `MD_LANE_FOR` is left
  by its end only.
- **Results by piece, not by lane**, where they are to be reproducible: how many lanes a group gets depends on timing.
  `md_lane_count` is known after the first sync, and the even split `md_lane_split` with it, for uniform work only.
- **A lane waits for nothing but its group**: no `md_parallel_for` inside, no other wait that could run pieces.
- **Cancellation is decided at syncs**, the same for every lane.

`md_spatial_acc_init` (the sort of its points into cells) and `md_linalg_sym_eigen_ex` (its tridiagonalization) have a
wide version next to their fork-join one, chosen by `md_executor_t.fork_join` while the two are compared:
`md_lanes_bench` times both on the same pool.

## Tests

`unittest/test_executor.h` has two executors for tests: real threads (`test_executor_threads`), and one thread
running the pieces out of order and handing out `thread_idx` in turn (`test_executor_shuffled`). A procedure whose
result depends on the order of its pieces, or on which `thread_idx` a piece had, differs under the second every
time, not once in a while.
