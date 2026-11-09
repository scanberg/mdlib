// md_load_bench: how a structure file loads over threads.
//
//   md_load_bench [--reps n] [--threads a,b,c] [--str | --mmap | --mmap-prefetch] [--reuse-memory] [--phases] file ...
//
//   --reps <n>       timed loads per thread count, the best and the median are reported (default 5)
//   --threads <list> thread counts to load with (default 1, 2, 4, ... and every logical core)
//   --str            the file read into memory first and loaded from the string: the parse alone, not the reads
//   --reuse-memory   every load into the same memory, reset between loads, so that only the first has pages to fault
//                    in: against a run without it, what the page faults of a load cost
//   --phases         the loader's own account of where the time of the best load of each count went (its debug line)
//   --mmap           the file mapped into memory and loaded from the mapping as a string: no copy out of the OS's cache,
//                    a page fault for every page instead (the mapping and unmapping are timed with the load)
//   --mmap-prefetch  the same, the OS asked to bring the pages in first (PrefetchVirtualMemory on Windows, MAP_POPULATE
//                    elsewhere)
//
// .pdb files load with md_pdb_system_init_from_file, .gro with md_gro_system_init_from_file. Each line is a thread
// count: the best and median wall time of a load, the speedup over one thread, the bytes of the file per second and the
// page faults of the best load. The system each count loads is summarized by a hash of its atoms, residues, bonds and
// types, which every count must agree on ('!' marks one that does not). The header names the benchmark's revision (bumped
// with every change of it) and the mdlib commit it was built from.

#include <md_pdb.h>
#include <md_gro.h>
#include <md_system.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>

#include <core/md_executor.h>
#include <core/md_common.h>
#include "bench_rev.h"                    // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if MD_PLATFORM_WINDOWS
#   define WIN32_LEAN_AND_MEAN
#   include <windows.h>
#   include <psapi.h>
#else
#   include <sys/resource.h>
#   include <sys/mman.h>
#   include <sys/stat.h>
#   include <fcntl.h>
#   include <unistd.h>
#endif

// The benchmark's own revision: bumped with every change of it, so that a result names the build it came from
#define MD_LOAD_BENCH_REV "R4"

#define MAX_REPS 64
#define MAX_THREADS 64

// The benchmark's own threads (mdlib owns none): a pool which keeps them, as an application's task system does. Workers
// sleep on a semaphore of their own between loops; a loop wakes them, takes pieces on the calling thread as well, and
// waits for them. Threads made and joined for every loop cost Windows about half a millisecond a loop at 8 threads, and a
// load runs hundreds of loops: that would be measured, and is no load's.
typedef struct bench_pool_t {
    md_executor_t     ex;
    uint32_t          num_workers;          // threads besides the calling one
    md_thread_t*      thread[MAX_THREADS];
    md_semaphore_t    start[MAX_THREADS];
    md_semaphore_t    done;
    volatile int32_t  quit;
    // The loop
    uint32_t          count, grain;
    md_range_fn       body;
    void*             arg;
    volatile uint64_t next;
} bench_pool_t;

typedef struct bench_worker_t {
    bench_pool_t* pool;
    uint32_t      idx;
} bench_worker_t;

static bench_worker_t bench_workers[MAX_THREADS];

static void bench_pool_run(bench_pool_t* pool, uint32_t thread_idx) {
    const uint64_t num_pieces = ((uint64_t)pool->count + pool->grain - 1) / pool->grain;
    for (;;) {
        const uint64_t p = md_atomic_fetch_add64(&pool->next, 1);
        if (p >= num_pieces) break;
        const uint32_t beg = (uint32_t)(p * pool->grain);
        const uint32_t end = (uint32_t)MIN((p + 1) * pool->grain, (uint64_t)pool->count);
        pool->body(pool->arg, beg, end, thread_idx);
    }
}

static void bench_worker_main(void* data) {
    bench_worker_t* w = (bench_worker_t*)data;
    bench_pool_t* pool = w->pool;
    for (;;) {
        md_semaphore_aquire(&pool->start[w->idx]);
        if (pool->quit) break;
        bench_pool_run(pool, w->idx + 1);
        md_semaphore_release(&pool->done);
    }
}

static bool bench_pool_for(void* user, uint32_t count, uint32_t grain, md_range_fn body, void* arg) {
    bench_pool_t* pool = (bench_pool_t*)user;
    if (count == 0) return true;
    pool->count = count;
    pool->grain = grain ? grain : 1;
    pool->body  = body;
    pool->arg   = arg;
    pool->next  = 0;
    // As many workers as there are pieces beyond the calling thread's
    const uint64_t num_pieces = ((uint64_t)count + pool->grain - 1) / pool->grain;
    const uint32_t woken = (uint32_t)MIN((uint64_t)pool->num_workers, num_pieces - 1);
    for (uint32_t i = 0; i < woken; ++i) md_semaphore_release(&pool->start[i]);
    bench_pool_run(pool, 0);
    for (uint32_t i = 0; i < woken; ++i) md_semaphore_aquire(&pool->done);
    return true;
}

static void bench_pool_init(bench_pool_t* pool, uint32_t num_threads) {
    MEMSET(pool, 0, sizeof(*pool));
    pool->num_workers = num_threads > 1 ? MIN(num_threads, MAX_THREADS) - 1 : 0;
    pool->ex = (md_executor_t){ .parallel_for = bench_pool_for, .user = pool, .num_threads = pool->num_workers + 1 };
    md_semaphore_init(&pool->done, 0);
    for (uint32_t i = 0; i < pool->num_workers; ++i) {
        md_semaphore_init(&pool->start[i], 0);
        bench_workers[i] = (bench_worker_t){ pool, i };
        pool->thread[i] = md_thread_create(bench_worker_main, &bench_workers[i]);
    }
}

static void bench_pool_free(bench_pool_t* pool) {
    pool->quit = 1;
    for (uint32_t i = 0; i < pool->num_workers; ++i) md_semaphore_release(&pool->start[i]);
    for (uint32_t i = 0; i < pool->num_workers; ++i) {
        md_thread_join(pool->thread[i]);
        md_semaphore_destroy(&pool->start[i]);
    }
    md_semaphore_destroy(&pool->done);
}

// The page faults of the process so far: soft and hard (Windows counts them as one)
static uint64_t page_faults(void) {
#if MD_PLATFORM_WINDOWS
    PROCESS_MEMORY_COUNTERS pmc = { .cb = sizeof(pmc) };
    return GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)) ? (uint64_t)pmc.PageFaultCount : 0;
#else
    struct rusage ru;
    return getrusage(RUSAGE_SELF, &ru) == 0 ? (uint64_t)(ru.ru_minflt + ru.ru_majflt) : 0;
#endif
}
#define MAX_COUNTS 32

static uint64_t fnv(uint64_t h, const void* p, size_t n) {
    const uint8_t* b = (const uint8_t*)p;
    for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ull; }
    return h;
}

static uint64_t system_hash(const md_system_t* sys, const md_system_state_t* st) {
    uint64_t h = 1469598103934665603ull;
    const size_t na = sys->atom.count, nr = sys->residue.count;
    h = fnv(h, &na, sizeof(na));
    h = fnv(h, st->xyz, na * sizeof(vec3_t));
    h = fnv(h, sys->atom.type_idx, na * sizeof(md_particle_type_idx_t));
    for (size_t t = 0; t < sys->atom.type.count; ++t) {
        const str_t name = md_particle_type_name(&sys->atom.type, t);
        h = fnv(h, name.ptr, name.len);
        h = fnv(h, &sys->atom.type.z[t], sizeof(sys->atom.type.z[t]));
    }
    if (nr) h = fnv(h, sys->residue.atom_offset, (nr + 1) * sizeof(uint32_t));
    h = fnv(h, sys->residue.seq_id, nr * sizeof(md_sequence_id_t));
    h = fnv(h, sys->bond.pairs, sys->bond.count * sizeof(md_atom_pair_t));
    if (sys->instance.count) h = fnv(h, sys->instance.res_offset, (sys->instance.count + 1) * sizeof(uint32_t));
    return h;
}

// A file mapped read only, as a string
typedef struct mapped_file_t {
    str_t text;
#if MD_PLATFORM_WINDOWS
    HANDLE file, mapping;
#else
    int fd;
#endif
} mapped_file_t;

#if MD_PLATFORM_WINDOWS
typedef struct bench_memory_range_t { void* addr; size_t bytes; } bench_memory_range_t;
typedef BOOL (WINAPI *prefetch_fn)(HANDLE, ULONG_PTR, bench_memory_range_t*, ULONG);
#endif

static bool map_file(mapped_file_t* m, const char* path, bool prefetch) {
    MEMSET(m, 0, sizeof(*m));
#if MD_PLATFORM_WINDOWS
    wchar_t wpath[1024];
    if (!MultiByteToWideChar(CP_UTF8, 0, path, -1, wpath, (int)ARRAY_SIZE(wpath))) return false;
    m->file = CreateFileW(wpath, GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (m->file == INVALID_HANDLE_VALUE) return false;
    LARGE_INTEGER size;
    if (!GetFileSizeEx(m->file, &size) || size.QuadPart == 0) { CloseHandle(m->file); return false; }
    m->mapping = CreateFileMappingW(m->file, NULL, PAGE_READONLY, 0, 0, NULL);
    if (!m->mapping) { CloseHandle(m->file); return false; }
    void* view = MapViewOfFile(m->mapping, FILE_MAP_READ, 0, 0, 0);
    if (!view) { CloseHandle(m->mapping); CloseHandle(m->file); return false; }
    m->text = (str_t){ (const char*)view, (size_t)size.QuadPart };
    if (prefetch) {
        // Windows 8 and later; looked up so that the benchmark builds against any SDK
        static prefetch_fn fn = NULL;
        if (!fn) fn = (prefetch_fn)(void*)GetProcAddress(GetModuleHandleW(L"kernel32.dll"), "PrefetchVirtualMemory");
        bench_memory_range_t range = { view, (size_t)size.QuadPart };
        if (fn) fn(GetCurrentProcess(), 1, &range, 0);
    }
    return true;
#else
    m->fd = open(path, O_RDONLY);
    if (m->fd < 0) return false;
    struct stat st;
    if (fstat(m->fd, &st) != 0 || st.st_size == 0) { close(m->fd); return false; }
    int flags = MAP_PRIVATE;
#   ifdef MAP_POPULATE
    if (prefetch) flags |= MAP_POPULATE;
#   endif
    void* view = mmap(NULL, (size_t)st.st_size, PROT_READ, flags, m->fd, 0);
    if (view == MAP_FAILED) { close(m->fd); return false; }
#   ifndef MAP_POPULATE
    if (prefetch) madvise(view, (size_t)st.st_size, MADV_WILLNEED);
#   endif
    m->text = (str_t){ (const char*)view, (size_t)st.st_size };
    return true;
#endif
}

static void unmap_file(mapped_file_t* m) {
#if MD_PLATFORM_WINDOWS
    UnmapViewOfFile(m->text.ptr);
    CloseHandle(m->mapping);
    CloseHandle(m->file);
#else
    munmap((void*)m->text.ptr, m->text.len);
    close(m->fd);
#endif
    MEMSET(m, 0, sizeof(*m));
}

static bool load(md_system_t* sys, md_system_state_t* st, str_t path, str_t text, bool gro, const md_executor_t* ex) {
    if (text.len) {
        return gro ? md_gro_system_init_from_str(sys, st, text, ex) : md_pdb_system_init_from_str(sys, st, text, MD_PDB_OPTION_NONE, ex);
    }
    return gro ? md_gro_system_init_from_file(sys, st, path, ex) : md_pdb_system_init_from_file(sys, st, path, MD_PDB_OPTION_NONE, ex);
}

static int cmp_double(const void* a, const void* b) {
    const double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

// The loaders' debug lines of the load being kept (--phases), errors always
static bool keep_debug = false;
static char kept_line[1024];

static void quiet(struct md_logger_o* inst, enum md_log_type_t type, const char* msg) {
    (void)inst;
    if (type == MD_LOG_TYPE_ERROR) fprintf(stderr, "%s\n", msg);
    if (keep_debug && type == MD_LOG_TYPE_DEBUG && (strstr(msg, "PDB: ") == msg || strstr(msg, "GRO: ") == msg)) {
        snprintf(kept_line, sizeof(kept_line), "%s", msg);
    }
}

int main(int argc, char** argv) {
    int reps = 5;
    bool from_str = false;
    int  mmap_mode = 0;     // 1: --mmap, 2: --mmap-prefetch
    bool reuse = false;
    bool phases = false;
    int counts[MAX_COUNTS];
    int num_counts = 0;
    const char* files[64];
    int num_files = 0;

    md_os_sys_info_t si = {0};
    md_os_sys_info_query(&si);
    const int cores = si.num_virtual_cores > 0 ? si.num_virtual_cores : 1;

    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (!strcmp(a, "--reps") && i + 1 < argc) {
            reps = atoi(argv[++i]);
        } else if (!strcmp(a, "--threads") && i + 1 < argc) {
            for (const char* c = argv[++i]; *c && num_counts < MAX_COUNTS; ) {
                const int n = atoi(c);
                if (n > 0) counts[num_counts++] = MIN(n, MAX_THREADS);
                while (*c && *c != ',') ++c;
                if (*c == ',') ++c;
            }
        } else if (!strcmp(a, "--str")) {
            from_str = true;
        } else if (!strcmp(a, "--mmap")) {
            mmap_mode = 1;
        } else if (!strcmp(a, "--mmap-prefetch")) {
            mmap_mode = 2;
        } else if (!strcmp(a, "--reuse-memory")) {
            reuse = true;
        } else if (!strcmp(a, "--phases")) {
            phases = true;
        } else if (a[0] == '-') {
            printf("usage: %s [--reps n] [--threads a,b,c] [--str | --mmap | --mmap-prefetch] [--reuse-memory] [--phases] file ...\n", argv[0]);
            return a[1] == 'h' ? 0 : 1;
        } else if (num_files < 64) {
            files[num_files++] = a;
        }
    }
    reps = CLAMP(reps, 1, MAX_REPS);
    if (num_counts == 0) {
        for (int n = 1; n < cores && num_counts < MAX_COUNTS - 1; n *= 2) counts[num_counts++] = n;
        counts[num_counts++] = MIN(cores, MAX_THREADS);
    }
    if (num_files == 0) {
        printf("usage: %s [--reps n] [--threads a,b,c] [--str | --mmap | --mmap-prefetch] [--reuse-memory] [--phases] file ...\n", argv[0]);
        return 1;
    }

    // Errors only: the loaders' notes on each of many loads are not the point
    static md_logger_i logger = { .log = quiet };
    md_log_unregister(default_logger);
    md_log_register(&logger);

    printf("md_load_bench %s: structure files loaded over threads\n", MD_LOAD_BENCH_REV);
    printf("mdlib %s\n", MD_BENCH_SOURCE_REV);
    printf("CPU %d logical cores | %d timed load%s per count | %s | %s\n\n", cores, reps, reps == 1 ? "" : "s",
           mmap_mode == 2 ? "from a mapping of the file, prefetched" : mmap_mode == 1 ? "from a mapping of the file" :
           from_str ? "from a string (the parse, not the reads)" : "from the file",
           reuse ? "into the same memory (faulted in once)" : "into fresh memory each load");

    for (int f = 0; f < num_files; ++f) {
        const str_t path = str_from_cstr(files[f]);
        const bool gro = str_ends_with(path, STR_LIT(".gro"));
        md_file_info_t info = {0};
        if (!md_file_info_extract_from_path(path, &info)) {
            fprintf(stderr, "%s: cannot be read\n", files[f]);
            continue;
        }
        str_t text = {0};
        if (from_str) text = load_textfile(path, md_get_heap_allocator());
        const double mb = (double)info.size / (1024.0 * 1024.0);

        printf("%s (%.1f MB)\n", files[f], mb);
        printf("  threads      best    median   speedup      MB/s    faults   system\n");
        double base = 0;
        uint64_t first_hash = 0;
        for (int c = 0; c < num_counts; ++c) {
            static bench_pool_t pool;
            bench_pool_init(&pool, (uint32_t)counts[c]);
            const md_executor_t* ex = counts[c] > 1 ? &pool.ex : NULL;
            double ms[MAX_REPS];
            uint64_t hash = 0;
            size_t num_atoms = 0;
            uint64_t best_faults = 0;
            double best_ms = 1e300;
            char best_line[1024] = "";
            bool ok = true;
            md_allocator_i* shared = reuse ? md_vm_arena_create(GIGABYTES(64)) : NULL;
            // One untimed load first: the page cache and the allocators as a load after another finds them
            for (int r = -1; r < reps && ok; ++r) {
                md_allocator_i* alloc = shared ? shared : md_vm_arena_create(GIGABYTES(64));
                if (shared) md_vm_arena_reset(shared);
                md_system_t sys = { .alloc = alloc };
                md_system_state_t st = { .alloc = alloc };
                keep_debug = phases;
                kept_line[0] = '\0';
                const uint64_t f0 = page_faults();
                const md_tick_t t0 = md_tick_now();
                if (mmap_mode) {
                    // The mapping is part of what is timed: made, read from and unmade by every load
                    mapped_file_t m;
                    ok = map_file(&m, files[f], mmap_mode == 2);
                    if (ok) {
                        ok = load(&sys, &st, path, m.text, gro, ex);
                        unmap_file(&m);
                    }
                } else {
                    ok = load(&sys, &st, path, text, gro, ex);
                }
                const md_tick_t t1 = md_tick_now();
                const uint64_t f1 = page_faults();
                keep_debug = false;
                if (ok && r >= 0) {
                    ms[r] = md_tick_to_milliseconds(t1 - t0);
                    if (ms[r] < best_ms) {
                        best_ms = ms[r];
                        best_faults = f1 - f0;
                        snprintf(best_line, sizeof(best_line), "%s", kept_line);
                    }
                }
                if (ok && r == -1) {
                    hash = system_hash(&sys, &st);
                    num_atoms = sys.atom.count;
                }
                if (!shared) md_vm_arena_destroy(alloc);
            }
            if (shared) md_vm_arena_destroy(shared);
            bench_pool_free(&pool);
            if (!ok) {
                printf("  %7d   failed to load\n", counts[c]);
                continue;
            }
            qsort(ms, (size_t)reps, sizeof(double), cmp_double);
            const double best = ms[0], median = ms[reps / 2];
            if (c == 0) {
                base = best;
                first_hash = hash;
            }
            printf("  %7d %9.1f %9.1f %8.2fx %9.0f %9llu   %zu atoms %016llx%s\n", counts[c], best, median, base / best, mb / (best * 1e-3),
                   (unsigned long long)best_faults, num_atoms, (unsigned long long)hash, hash == first_hash ? "" : " !");
            if (phases && best_line[0]) {
                // Without the source location the log appends
                char* at = strstr(best_line, " [");
                if (at) *at = '\0';
                printf("            %s\n", best_line);
            }
        }
        printf("\n");
        if (text.ptr) str_free(text, md_get_heap_allocator());
    }
    return 0;
}
