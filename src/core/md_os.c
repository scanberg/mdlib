#include "md_os.h"

#include "md_common.h"
#include "md_platform.h"
#include "md_log.h"
#include "md_allocator.h"

#include <stdio.h>
#include <stdarg.h>

#define MD_MAX_PATH 4096

#include <time.h>

#if MD_PLATFORM_WINDOWS

#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN 1
#endif

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN 1
#endif

#include <Windows.h>
#include <Shlwapi.h>
#include <direct.h>

// If we need to explicitly link against some lib
#pragma comment(lib, "User32.lib")
#pragma comment(lib, "Shlwapi.lib")

typedef LONG NTSTATUS;

typedef NTSTATUS (NTAPI *_NtQuerySemaphore)(
    HANDLE SemaphoreHandle, 
    DWORD SemaphoreInformationClass, /* Would be SEMAPHORE_INFORMATION_CLASS */
    PVOID SemaphoreInformation,      /* but this is to much to dump here     */
    ULONG SemaphoreInformationLength, 
    PULONG ReturnLength OPTIONAL
    ); 

typedef struct _SEMAPHORE_BASIC_INFORMATION {   
    ULONG CurrentCount; 
    ULONG MaximumCount;
} SEMAPHORE_BASIC_INFORMATION;

STATIC_ASSERT(sizeof(CRITICAL_SECTION) <= sizeof(md_mutex_t), "Win32 CRITICAL_SECTION does not fit into md_mutex_t!");

#elif MD_PLATFORM_UNIX

#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <pthread.h>

#endif

#if MD_PLATFORM_LINUX
#include <semaphore.h>
STATIC_ASSERT(sizeof(sem_t) <= sizeof(md_semaphore_t), "Linux sem_t does not fit into md_semaphore_t!");
#endif

#if MD_PLATFORM_OSX
#include <sys/param.h>
#include <mach/mach_init.h>
#include <mach/task.h>
#include <mach/semaphore.h>

STATIC_ASSERT(sizeof(semaphore_t) <= sizeof(md_semaphore_t), "MacOS semaphore_t does not fit into md_semaphore_t!");
#endif

#include <string.h>

#define SZBUF_FROM_STR(buf, str) strncpy(buf, str.ptr, MIN(sizeof(buf)-1, str.len))

static char CWD[4096];

#if MD_PLATFORM_WINDOWS
// https://docs.microsoft.com/en-us/windows/win32/debug/retrieving-the-last-error-code
void print_windows_error() {
    LPVOID msg_buf;
    DWORD err_code = GetLastError();

    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER |
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        err_code,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR)&msg_buf,
        0, NULL);

    MD_LOG_ERROR("Error code (%d): %s", err_code, msg_buf);
    LocalFree(msg_buf);
}
#endif

str_t md_path_cwd() {
#if MD_PLATFORM_WINDOWS
    char* val = _getcwd(CWD, sizeof(CWD));
    ASSERT(val != 0);
#elif MD_PLATFORM_UNIX
    getcwd(CWD, sizeof(CWD));
#else
    ASSERT(false);
#endif
    str_t res = { .ptr = CWD, .len = (int64_t)strnlen(CWD, sizeof(CWD)) };
    return res;
}

static inline str_t internal_fullpath(str_t path, md_allocator_i* alloc) {
    path = str_copy(path, default_temp_allocator); // Zero terminate

    char sz_buf[MD_MAX_PATH] = "";
#if MD_PLATFORM_WINDOWS
    int64_t len = (int64_t)GetFullPathName(path.ptr, sizeof(sz_buf), sz_buf, NULL);
    if (len == 0) {
        print_windows_error();
    }
#elif MD_PLATFORM_UNIX
    int64_t len = 0;
    if (realpath(path.ptr, sz_buf) != NULL) {
        len = (int64_t)strnlen(sz_buf, sizeof(sz_buf));
    }
    else {
        switch (errno) {
        case EACCES:        MD_LOG_ERROR("Read or search permission was denied for a component of the path prefix."); break;
        case EINVAL:        MD_LOG_ERROR("path is NULL or resolved_path is NULL."); break;
        case EIO:           MD_LOG_ERROR("An I/O error occurred while reading from the filesystem."); break;
        case ELOOP:         MD_LOG_ERROR("Too many symbolic links were encountered in translating the pathname."); break;
        case ENAMETOOLONG:  MD_LOG_ERROR("A component of a pathname exceeded NAME_MAX characters, or an entire pathname exceeded PATH_MAX characters."); break;
        case ENOENT:        MD_LOG_ERROR("The named file does not exist."); break;
        case ENOMEM:        MD_LOG_ERROR("Out of memory."); break;
        case ENOTDIR:       MD_LOG_ERROR("A component of the path prefix is not a directory."); break;
        default:            MD_LOG_ERROR("Undefined error"); break;
        }
    }
#else
    ASSERT(false);
#endif
    str_t result = { 0 };
    if (len > 0) {
        result = str_copy((str_t) { sz_buf, len }, alloc);
    }
    return result;
}

str_t md_path_make_canonical(str_t path, struct md_allocator_i* alloc) {
    ASSERT(alloc);

    path = internal_fullpath(path, alloc);

    if (path.len > 0) {
#if MD_PLATFORM_WINDOWS
        convert_backslashes((char*)path.ptr, path.len);
#endif
    }
    else {
        MD_LOG_ERROR("Failed to create canonical path!");
    }

    return path;
}

str_t md_path_make_relative(str_t from, str_t to, struct md_allocator_i* alloc) {
    char sz_buf[MD_MAX_PATH] = "";

    str_t relative_str = { 0 };
    bool result = false;

    // Make 2 canonical paths
    from = internal_fullpath(from, default_temp_allocator);
    to = internal_fullpath(to, default_temp_allocator);

#if MD_PLATFORM_WINDOWS
    result = PathRelativePathTo(sz_buf, from.ptr, FILE_ATTRIBUTE_NORMAL, to.ptr, FILE_ATTRIBUTE_NORMAL);
    convert_backslashes(sz_buf, sizeof(sz_buf));
#elif MD_PLATFORM_UNIX

    // Find the common base
    int64_t count = str_count_equal_chars(from, to);
    result = count > 0;

    if (result) {
        // Count number of folders as N in from and add N times '../'
        str_t rfrom = str_substr(from, count, -1);
        str_t rto = str_substr(to, count, -1);

        const int64_t folder_count = str_count_occur_char(rfrom, '/');

        str_t folder_up = STR("../");
        int len = 0;
        for (int64_t i = 0; i < folder_count; ++i) {
            len += snprintf(sz_buf + len, sizeof(sz_buf) - len, "%.*s", (int)folder_up.len, folder_up.ptr);
        }
        len += snprintf(sz_buf + len, sizeof(sz_buf) - len, "%.*s", (int)rto.len, rto.ptr);
    }
#else
    ASSERT(false);
#endif
    if (result) {
        relative_str = str_copy((str_t) { .ptr = sz_buf, .len = (int64_t)strnlen(sz_buf, sizeof(sz_buf)) }, alloc);
    }
    else {
        MD_LOG_ERROR("Failed to extract relative path.");
    }
    return relative_str;
}

bool md_path_is_valid(str_t path) {
#if MD_PLATFORM_WINDOWS
    path = str_copy(path, default_temp_allocator);
    bool result = PathFileExists(path.ptr);
#elif MD_PLATFORM_UNIX
    bool result = (access(path.ptr, F_OK) == 0);
#else
    ASSERT(false);
#endif
    return result;
}

bool md_path_is_directory(str_t path) {
#if MD_PLATFORM_WINDOWS
    path = str_copy(path, default_temp_allocator);
    bool result = PathIsDirectory(path.ptr);
#elif MD_PLATFORM_UNIX
    bool result = false;
    struct stat s;
    if (stat(path.ptr, &s) == 0) {
        if (s.st_mode & S_IFDIR) {
            result = true;
        }
    }
#else
    ASSERT(false);
#endif
    return result;
}


// ### FILE ###
struct md_file_o {
    FILE handle;
};

md_file_o* md_file_open(str_t filename, uint32_t flags) {
    if (!filename.ptr || !filename.len) return NULL;

#if MD_PLATFORM_WINDOWS
    wchar_t w_file[MD_MAX_PATH] = {0};
    const int w_file_len = MultiByteToWideChar(CP_UTF8, 0, filename.ptr, (int)filename.len, w_file, ARRAY_SIZE(w_file));
    if (w_file_len >= ARRAY_SIZE(w_file)) {
        MD_LOG_ERROR("File path exceeds limit!");
        return NULL;
    }
    w_file[w_file_len] = L'\0';

    const wchar_t* w_mode = 0;
    switch (flags) {
    case MD_FILE_READ:
        w_mode = L"r";
        break;
    case (MD_FILE_READ | MD_FILE_BINARY):
        w_mode = L"rb";
        break;
    case MD_FILE_WRITE:
        w_mode = L"w";
        break;
    case (MD_FILE_WRITE | MD_FILE_BINARY):
        w_mode = L"wb";
        break;
    case MD_FILE_APPEND:
        w_mode = L"a";
        break;
    case (MD_FILE_APPEND | MD_FILE_BINARY):
        w_mode = L"ab";
        break;
    default:
        MD_LOG_ERROR("Invalid combination of file access flags!");
        return NULL;
    }

    return (md_file_o*) _wfopen(w_file, w_mode);
#else
    const char* mode = 0;
    switch (flags) {
    case MD_FILE_READ:
        mode = "r";
        break;
    case (MD_FILE_READ | MD_FILE_BINARY):
        mode = "rb";
        break;
    case MD_FILE_WRITE:
        mode = "w";
        break;
    case (MD_FILE_WRITE | MD_FILE_BINARY):
        mode = "wb";
        break;
    case MD_FILE_APPEND:
        mode = "a";
        break;
    case (MD_FILE_APPEND | MD_FILE_BINARY):
        mode = "ab";
        break;
    default:
        MD_LOG_ERROR("Invalid combination of file access flags!");
        return NULL;
    }

    char file[MD_MAX_PATH] = "";
    strncpy(file, filename.ptr, ARRAY_SIZE(file)-1);
    return (md_file_o*)fopen(file, mode);
#endif
}

void md_file_close(md_file_o* file) {
    fclose((FILE*)file);
}

int64_t md_file_tell(md_file_o* file) {
#if MD_PLATFORM_WINDOWS
    return _ftelli64((FILE*)file);
#else
    return ftello((FILE*)file);
#endif
}

bool md_file_seek(md_file_o* file, int64_t offset, md_file_seek_origin_t origin) {
    ASSERT(file);

    int o = 0;
    switch(origin) {
    case MD_FILE_BEG:
        o = SEEK_SET;
        break;
    case MD_FILE_CUR:
        o = SEEK_CUR;
        break;
    case MD_FILE_END:
        o = SEEK_END;
        break;
    default:
        MD_LOG_ERROR("Invalid seek origin value!");
        return false;
    }

#if MD_PLATFORM_WINDOWS
    return _fseeki64((FILE*)file, offset, o) == 0;
#else
    return fseeko((FILE*)file, offset, o) == 0;
#endif
}

int64_t md_file_size(md_file_o* file) {
    ASSERT(file);
    int64_t cur = md_file_tell(file);
    md_file_seek(file, 0, SEEK_END);
    int64_t end = md_file_tell(file);
    md_file_seek(file, cur, SEEK_SET);
    return end;
}

// Returns the number of successfully written/read bytes
int64_t md_file_read(md_file_o* file, void* ptr, int64_t num_bytes) {
    ASSERT(file);
    return fread(ptr, 1, num_bytes, (FILE*)file);
}

int64_t md_file_read_line(md_file_o* file, char* buf, int64_t cap) {
    int64_t pos = md_file_tell(file);
    char* res = fgets(buf, (int)cap, (FILE*)file);
    int64_t len = md_file_tell(file) - pos;
    return res ? len : 0;
}

int64_t md_file_read_lines(md_file_o* file, char* buf, int64_t cap, int64_t line_count) {
    if (!file || !buf || cap < 1 || line_count < -1) return 0;
    line_count = (line_count == -1) ? INT64_MAX : line_count;
    int64_t tot_size = 0;
    for (int64_t i = 0; i < line_count && cap > 0; ++i) {
        int64_t read_size = md_file_read_line(file, buf, cap);
        if (read_size == 0) break;
        buf += read_size;
        cap -= read_size;
        tot_size += read_size;
    }
    return tot_size;
}

int64_t md_file_write(md_file_o* file, const void* ptr, int64_t num_bytes) {
    ASSERT(file);
    return fwrite(ptr, 1, num_bytes, (FILE*)file);
}

int64_t md_file_printf(md_file_o* file, const char* format, ...) {
    ASSERT(file);
    va_list args;
    va_start (args, format);
    int res = vfprintf((FILE*)file, format, args);
    va_end (args);
    return res;
}


// ### TIME ###
md_timestamp_t md_time_current() {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER t;
    QueryPerformanceCounter(&t);
    return t.QuadPart;
#elif MD_PLATFORM_UNIX
    struct timespec t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
    return t.tv_sec * 1000000000 + t.tv_nsec;
#else
    ASSERT(false);
#endif
}

double md_time_as_nanoseconds(md_timestamp_t t) {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    return (t * 1E6) / (double)frequency.QuadPart;
#elif MD_PLATFORM_UNIX
    return t * 1.0e-3;
#else
    ASSERT(false);
#endif
}

double md_time_as_milliseconds(md_timestamp_t t) {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    return (t * 1E3) / (double)frequency.QuadPart;
#elif MD_PLATFORM_UNIX
    return t * 1.0e-6;
#else
    ASSERT(false);
#endif
}

double md_time_as_seconds(md_timestamp_t t) {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    return (double)(t) / (double)frequency.QuadPart;
#elif MD_PLATFORM_UNIX
    return t * 1.0e-9;
#else
    ASSERT(false);
#endif
}

// MEM
// Linux equivalent of virtual memory allocation is partly taken from here
// https://forums.pcsx2.net/Thread-blog-VirtualAlloc-on-Linux

uint64_t md_physical_ram() {
#if MD_PLATFORM_WINDOWS
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys;
#elif MD_PLATFORM_UNIX
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
#else
    ASSERT(false);
    return 0;
#endif
}

uint64_t md_vm_page_size(void) {
#if MD_PLATFORM_WINDOWS
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    return info.dwPageSize;
#elif MD_PLATFORM_UNIX
    return sysconf(_SC_PAGE_SIZE);
#else
    ASSERT(false);
#endif
}

void* md_vm_reserve(uint64_t size) {
    uint64_t gb_snapped_size = ALIGN_TO(size, GIGABYTES(1));
#if MD_PLATFORM_WINDOWS
    return VirtualAlloc(0, gb_snapped_size, MEM_RESERVE, PAGE_NOACCESS);
#elif MD_PLATFORM_UNIX
    void* result = mmap(0, gb_snapped_size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    ASSERT(result != MAP_FAILED);
    return result;
#else
    ASSERT(false);
#endif
}

void md_vm_release(void* ptr, uint64_t size) {
#if MD_PLATFORM_WINDOWS
    (void)size;
    VirtualFree(ptr, 0, MEM_RELEASE);
#elif MD_PLATFORM_UNIX
    uint64_t gb_snapped_size = ALIGN_TO(size, GIGABYTES(1));
    int result;
    //result = madvise(ptr, -1, MADV_DONTNEED);
    //ASSERT(result == 0);
    result = munmap(ptr, gb_snapped_size);
    ASSERT(result == 0);
#else
    ASSERT(false);
#endif
}

void md_vm_commit(void* ptr, uint64_t size) {
    uint64_t page_snapped_size = ALIGN_TO(size, md_vm_page_size());
#if MD_PLATFORM_WINDOWS
    void* result = VirtualAlloc(ptr, page_snapped_size, MEM_COMMIT, PAGE_READWRITE);
    ASSERT(result != NULL);
#elif MD_PLATFORM_UNIX
    int result = mprotect(ptr, page_snapped_size, PROT_READ | PROT_WRITE);
    ASSERT(result == 0);
#else
    ASSERT(false);
#endif
}

void md_vm_decommit(void* ptr, uint64_t size) {
#if MD_PLATFORM_WINDOWS
    VirtualFree(ptr, size, MEM_DECOMMIT);
#elif MD_PLATFORM_UNIX
    int res;
    uint64_t page_snapped_size = ALIGN_TO(size, md_vm_page_size());
    res = mprotect(ptr, page_snapped_size, PROT_NONE);
    ASSERT(res == 0);
    res = madvise(ptr, page_snapped_size, MADV_DONTNEED);
    ASSERT(res == 0);
#else
    ASSERT(false);
#endif
}



// ### THREAD ###

// Thread
md_thread_t* md_thread_create(md_thread_entry func, void* user_data) {
#if MD_PLATFORM_WINDOWS
    DWORD unused;
    HANDLE id = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)func, user_data, 0, &unused);
    return (md_thread_t*)id;
#elif MD_PLATFORM_UNIX
#if MD_COMPILER_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-function-type"
#endif
    pthread_t thread;
    pthread_create(&thread, NULL, (void* (*)(void*))func, user_data);
    return (md_thread_t*)thread;
#if MD_COMPILER_GCC
#pragma GCC diagnostic pop
#endif
#else
    ASSERT(false);
#endif
}

void md_thread_detach(md_thread_t* thread) {
#if MD_PLATFORM_WINDOWS
    CloseHandle((HANDLE)thread);
#elif MD_PLATFORM_UNIX
    pthread_detach((pthread_t)thread);
#else
    ASSERT(false);
#endif
}

bool md_thread_join(md_thread_t* thread) {
#if MD_PLATFORM_WINDOWS
    WaitForSingleObject((HANDLE)thread, INFINITE);
    return CloseHandle((HANDLE)thread);
#elif MD_PLATFORM_UNIX
    return pthread_join((pthread_t)thread, NULL) == 0;
#else
    ASSERT(false);
#endif
}

md_thread_id_t md_thread_get_id(md_thread_t* thread) {
#if MD_PLATFORM_WINDOWS
    return GetThreadId((HANDLE)thread);
#elif MD_PLATFORM_UNIX
    return (md_thread_id_t)thread;
#else
    ASSERT(false);
#endif
}

md_thread_id_t md_thread_id(void) {
#if MD_PLATFORM_WINDOWS
    return GetCurrentThreadId();
#elif MD_PLATFORM_UNIX
    return (md_thread_id_t)pthread_self();
#else
    ASSERT(false);
#endif
}

bool md_thread_on_exit(md_thread_exit callback) {
    if (!callback) {
        MD_LOG_ERROR("OS: callback was NULL");
        return false;
    }

#if MD_PLATFORM_WINDOWS
    DWORD idx = FlsAlloc(callback);
    return idx != FLS_OUT_OF_INDEXES;
#elif MD_PLATFORM_UNIX
    pthread_key_t key;
    int ret = pthread_key_create(&key, callback);
    return ret == 0;
#else
    ASSERT(false);
    return false;
#endif
}

void md_thread_sleep(uint64_t milliseconds) {
#if MD_PLATFORM_WINDOWS
    Sleep((DWORD)milliseconds);
#elif MD_PLATFORM_UNIX
    usleep(milliseconds * 1000);
#else
    ASSERT(false);
#endif
}

// Mutex

bool md_mutex_init(md_mutex_t* mutex) {
#if MD_PLATFORM_WINDOWS
    return InitializeCriticalSectionAndSpinCount((CRITICAL_SECTION*)mutex, 2000);
#elif MD_PLATFORM_UNIX
    return pthread_mutex_init((pthread_mutex_t*)mutex, NULL) == 0;
#else
    ASSERT(false);
#endif
}

md_mutex_t md_mutex_create() {
#if MD_PLATFORM_WINDOWS
    md_mutex_t mutex;
    md_mutex_init(&mutex);
    return mutex;
#elif MD_PLATFORM_UNIX
    md_mutex_t mutex;
    md_mutex_init(&mutex);
    return mutex;
#else
    ASSERT(false);
#endif
}

bool md_mutex_destroy(md_mutex_t* mutex) {
#if MD_PLATFORM_WINDOWS
    DeleteCriticalSection((CRITICAL_SECTION*)mutex);
    return true;
#elif MD_PLATFORM_UNIX
    return pthread_mutex_destroy((pthread_mutex_t*)mutex) == 0;
#else
    ASSERT(false);
#endif
}

bool md_mutex_lock(md_mutex_t* mutex) {
#if MD_PLATFORM_WINDOWS
    EnterCriticalSection((CRITICAL_SECTION*)mutex);
    return true;
#elif MD_PLATFORM_UNIX
    return pthread_mutex_lock((pthread_mutex_t*)mutex) == 0;
#else
    ASSERT(false);
#endif
}

bool md_mutex_try_lock(md_mutex_t* mutex) {
#if MD_PLATFORM_WINDOWS
    return TryEnterCriticalSection((CRITICAL_SECTION*)mutex);
#elif MD_PLATFORM_UNIX
    return pthread_mutex_trylock((pthread_mutex_t*)mutex) == 0;
#else
    ASSERT(false);
#endif
}

bool md_mutex_unlock(md_mutex_t* mutex) {
#if MD_PLATFORM_WINDOWS
    LeaveCriticalSection((CRITICAL_SECTION*)mutex);
    return true;
#elif MD_PLATFORM_UNIX
    return pthread_mutex_unlock((pthread_mutex_t*)mutex) == 0;
#else
    ASSERT(false);
#endif
}



// Semaphore

bool md_semaphore_init(md_semaphore_t* semaphore, int32_t initial_count) {
#if MD_PLATFORM_WINDOWS
    semaphore->_data[0] = CreateSemaphoreA(NULL, (LONG)initial_count, MAXLONG, NULL);
    return semaphore != NULL;
#elif MD_PLATFORM_LINUX
    return sem_init((sem_t*)semaphore, 0, (uint32_t)initial_count) == 0;
#elif MD_PLATFORM_OSX
    semaphore_t sema;
    mach_port_t self = mach_task_self();
    kern_return_t ret = semaphore_create(self, &sema, SYNC_POLICY_FIFO, initial_count);
    if (ret != KERN_SUCCESS) {
        MD_LOG_ERROR("Failed to initialize semaphore");
    }
    semaphore->_data[0] = (void*)(uint64_t)sema;
    return ret == KERN_SUCCESS;
#endif
}

#if MD_PLATFORM_WINDOWS
static inline bool win32_semaphore_wait(md_semaphore_t* semaphore, DWORD milliseconds) {
    return WaitForSingleObjectEx((HANDLE)semaphore->_data[0], milliseconds, FALSE) == WAIT_OBJECT_0;
}
#endif

bool md_semaphore_destroy(md_semaphore_t* semaphore) {
#if MD_PLATFORM_WINDOWS
    return CloseHandle((HANDLE)semaphore->_data[0]);
#elif MD_PLATFORM_LINUX
    return sem_destroy((sem_t*)semaphore) == 0;
#elif MD_PLATFORM_OSX
    mach_port_t self = mach_task_self();
    return semaphore_destroy(self, (semaphore_t)semaphore->_data[0]) == KERN_SUCCESS;
#endif
}

bool md_semaphore_aquire(md_semaphore_t* semaphore) {
#if MD_PLATFORM_WINDOWS
    return win32_semaphore_wait(semaphore, INFINITE);
#elif MD_PLATFORM_LINUX
    return sem_wait((sem_t*)semaphore) == 0;
#elif MD_PLATFORM_OSX
    return semaphore_wait((semaphore_t)semaphore->_data[0]) == KERN_SUCCESS;
#endif
}

bool md_semaphore_try_aquire(md_semaphore_t* semaphore) {
#if MD_PLATFORM_WINDOWS
    return win32_semaphore_wait(semaphore, 0);
#elif MD_PLATFORM_LINUX
    return sem_trywait((sem_t*)semaphore) == 0;
#elif MD_PLATFORM_OSX
    mach_timespec_t mts;
    mts.tv_sec = 0;
    mts.tv_nsec = 0;
    return semaphore_timedwait((semaphore_t)semaphore->_data[0], mts) == KERN_SUCCESS;
#endif
}

bool md_semaphore_query_count(md_semaphore_t* semaphore, int32_t* count) {
    ASSERT(semaphore);
    ASSERT(count);
#if MD_PLATFORM_WINDOWS

    _NtQuerySemaphore NtQuerySemaphore = (_NtQuerySemaphore)GetProcAddress(GetModuleHandle("ntdll.dll"), "NtQuerySemaphore");

    if (NtQuerySemaphore) {
        SEMAPHORE_BASIC_INFORMATION basic_info = {0};
        NTSTATUS status = NtQuerySemaphore((HANDLE)semaphore->_data[0], 0 /*SemaphoreBasicInformation*/,  &basic_info, sizeof (SEMAPHORE_BASIC_INFORMATION), NULL);
        if (status == ERROR_SUCCESS) {
            *count = basic_info.CurrentCount;
            return true;
        }
    }
    return false;
#elif MD_PLATFORM_LINUX
    return sem_getvalue((sem_t*)semaphore, count) == 0;
#elif MD_PLATFORM_OSX
    (void)semaphore;
    (void)count;
    // Not supported
    return false;
#endif
}

bool md_semaphore_release(md_semaphore_t* semaphore) {
#if MD_PLATFORM_WINDOWS
    return ReleaseSemaphore((HANDLE)semaphore->_data[0], 1, NULL);
#elif MD_PLATFORM_LINUX
    return sem_post((sem_t*)semaphore) == 0;
#elif MD_PLATFORM_OSX
    return semaphore_signal((semaphore_t)semaphore->_data[0]) == KERN_SUCCESS;
#endif
}

// Common high level operations

md_semaphore_t md_semaphore_create(int32_t initial_count) {
    md_semaphore_t semaphore;
    md_semaphore_init(&semaphore, initial_count);
    return semaphore;
}

bool md_semaphore_try_aquire_n(md_semaphore_t* semaphore, int32_t count) {
    ASSERT(count > 0);
    int32_t aquired_count = 0;
    for (int32_t i = 0; i < count; ++i) {
        aquired_count += md_semaphore_try_aquire(semaphore) == true ? 1 : 0;
    }
    if (aquired_count == count) {
        return true;
    }
    md_semaphore_release_n(semaphore, aquired_count);
    return false;
}

bool md_semaphore_release_n(md_semaphore_t* semaphore, int32_t count) {
    bool result = true;
    for (int32_t i = 0; i < count; ++i) {
        result |= md_semaphore_release(semaphore);
    }
    return result;
}
