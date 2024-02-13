#include <core/md_os.h>

#include <core/md_common.h>
#include <core/md_platform.h>
#include <core/md_log.h>
#include <core/md_allocator.h>

#include <stdio.h>
#include <stdarg.h>
#include <time.h>

#define MD_MAX_PATH 4096

#if MD_PLATFORM_WINDOWS

#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN 1
#endif

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN 1
#endif

#include <Windows.h>
#include <Shlwapi.h>
#include <ShlObj.h>
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

static size_t win32_allocation_granularity(void) {
    static size_t allocation_granularity = 0;
    if (!allocation_granularity) {
        SYSTEM_INFO info;
        GetSystemInfo(&info);
        allocation_granularity = info.dwAllocationGranularity;
    }
    return allocation_granularity;
}

#elif MD_PLATFORM_UNIX

#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <pthread.h>
#include <pwd.h>
#include <fcntl.h>

#endif

#if MD_PLATFORM_LINUX
#include <semaphore.h>
STATIC_ASSERT(sizeof(sem_t) <= sizeof(md_semaphore_t), "Linux sem_t does not fit into md_semaphore_t!");
#endif

#if MD_PLATFORM_OSX
#include <limits.h>
#include <sys/param.h>
#include <mach/mach_init.h>
#include <mach/task.h>
#include <mach/semaphore.h>
#include <mach-o/dyld.h>

STATIC_ASSERT(sizeof(semaphore_t) <= sizeof(md_semaphore_t), "MacOS semaphore_t does not fit into md_semaphore_t!");
#endif

#include <string.h>

#define SZBUF_FROM_STR(buf, str) strncpy(buf, str.ptr, MIN(sizeof(buf)-1, str.len))

#if MD_PLATFORM_WINDOWS
// https://docs.microsoft.com/en-us/windows/win32/debug/retrieving-the-last-error-code
static void print_windows_error() {
    LPVOID msg_buf = 0;
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

static size_t fullpath(char* buf, size_t cap, str_t path) {
    str_t zpath = str_copy(path, md_temp_allocator); // Zero terminate
    if (zpath.len == 0) return 0;
    
#if MD_PLATFORM_WINDOWS
    size_t len = GetFullPathName(zpath.ptr, (DWORD)cap, buf, NULL);
    if (len == 0) {
        print_windows_error();
        return 0;
    }
    return len;

#elif MD_PLATFORM_UNIX
    size_t len = 0;
    if (realpath(zpath.ptr, buf) != NULL) {
        // realpath will not append a trailing '/' if the path is a directory. We want this to
        // be able to resolve relative paths more easily
        len = strnlen(buf, cap);
        if (len > 0 && md_path_is_directory((str_t){buf, len}) && buf[len-1] != '/') {
            if (len < cap) {
                buf[len++] = '/';
                buf[len]   = '\0';
            }
        }
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
    return len;
#else
#error "Platform not supported"
#endif
}

size_t md_path_write_cwd(char* buf, size_t cap) {
    char* val;
#if MD_PLATFORM_WINDOWS
    val = _getcwd(buf, (int)cap);
#elif MD_PLATFORM_UNIX
    val = getcwd(buf, (size_t)cap);
#else
    ASSERT(false);
#endif
    if (!val) {
    	MD_LOG_ERROR("Failed to get current working directory");
		return 0;
    }
    return strnlen(buf, cap);
}

size_t md_path_write_exe(char* buf, size_t buf_cap) {
#if MD_PLATFORM_WINDOWS
    DWORD res = GetModuleFileName(NULL, buf, (DWORD)buf_cap);
    if (res != 0) {
        return (size_t)res;
    }
#elif MD_PLATFORM_LINUX
    ssize_t res = readlink("/proc/self/exe", buf, buf_cap);
    if (res != -1) {
    	return (size_t)res;
    }
#elif MD_PLATFORM_OSX
	uint32_t buf_cap32 = (uint32_t)buf_cap;
	int res = _NSGetExecutablePath(buf, &buf_cap32);
	if (res == 0) {
    	return (size_t)buf_cap32;
    }
#else
    ASSERT(false);
#endif
    return 0;
}

// https://stackoverflow.com/questions/2899013/how-do-i-get-the-application-data-path-in-windows-using-c
// https://stackoverflow.com/questions/2910377/get-home-directory-in-linux
size_t md_path_write_user_dir(char* buf, size_t buf_cap) {
#if MD_PLATFORM_WINDOWS
    WCHAR* homedir = NULL;
    HRESULT res = SHGetKnownFolderPath(&FOLDERID_Profile, 0, NULL, &homedir);
    if (res != S_OK) {
		MD_LOG_ERROR("Failed to get home directory");
        return 0;
	}
	res = WideCharToMultiByte(CP_UTF8, 0, homedir, -1, buf, (int)buf_cap, NULL, NULL);
	CoTaskMemFree(homedir);
	if (res == 0) {
		MD_LOG_ERROR("Failed to convert home directory to UTF8");
		return 0;
	}
	return strnlen(buf, buf_cap);
#elif MD_PLATFORM_UNIX
    const char *homedir = NULL;
    if ((homedir = getenv("HOME")) == NULL) {
        struct passwd* pwd = getpwuid(getuid());
        if (pwd) {
			homedir = pwd->pw_dir;
		}
    }
    if (!homedir) {
        MD_LOG_ERROR("Failed to get home directory");
		return 0;
    }
    return (size_t)snprintf(buf, buf_cap, "%s", homedir);
#else
    ASSERT(false);
#endif
}

bool md_path_set_cwd(str_t path) {
    // Ensure zero termination of path
    char buf[MD_MAX_PATH];
    str_copy_to_char_buf(buf, sizeof(buf), path);
#if MD_PLATFORM_WINDOWS
    if (SetCurrentDirectory(buf)) {
        return true;
    }
    print_windows_error();
#elif MD_PLATFORM_UNIX
    if (chdir(buf) == 0) {
        return true;
    }
    switch (errno) {
    case EACCES:        MD_LOG_ERROR("CWD: Permission denied for some component of the path '%.*s'", (int)path.len, path.ptr); break;
    case EFAULT:        MD_LOG_ERROR("CWD: Path points outside of accessible address space '%.*s'", (int)path.len, path.ptr); break;
    case EIO:           MD_LOG_ERROR("CWD: I/O Error '%.*s'", (int)path.len, path.ptr); break;
    case ELOOP:         MD_LOG_ERROR("CWD: Too many symbolic links encountered in resolving path '%.*s'", (int)path.len, path.ptr); break;
    case ENAMETOOLONG:  MD_LOG_ERROR("CWD: Path is too long '%.*s'", (int)path.len, path.ptr); break;
    case ENOENT:        MD_LOG_ERROR("CWD: Path does not exist '%.*s'", (int)path.len, path.ptr); break;
    case ENOMEM:        MD_LOG_ERROR("CWD: Insufficient kernel memory '%.*s'", (int)path.len, path.ptr); break;
    case ENOTDIR:       MD_LOG_ERROR("CWD: Path is not a directory '%.*s'", (int)path.len, path.ptr); break;
    default:            MD_LOG_ERROR("CWD: Undefined error"); break;
    }
#else
    ASSERT(false);
#endif
    return false;
}

size_t md_path_write_canonical(char* buf, size_t cap, str_t path) {
    size_t len = fullpath(buf, cap, path);
#if MD_PLATFORM_WINDOWS
    replace_char(buf, len, '\\', '/');
#endif
    return len;
}

str_t md_path_make_canonical(str_t path, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    char buf[MD_MAX_PATH];
    const size_t len = fullpath(buf, sizeof(buf), path);
    str_t result = {0};

    if (path.len > 0 && len == 0) {
        MD_LOG_ERROR("Failed to create canonical path from '%.*s'", (int)path.len, path.ptr);
        return result;
    }

#if MD_PLATFORM_WINDOWS
    replace_char(buf, len, '\\', '/');
#endif
    
    return str_copy_cstrn(buf, len, alloc);
}

size_t md_path_write_relative(char* out_buf, size_t out_cap, str_t from, str_t to) {
    char from_buf[MD_MAX_PATH];
    char   to_buf[MD_MAX_PATH];

    bool success = false;
    size_t len = 0;

    // Make 2 canonical paths
    const size_t from_len = fullpath(from_buf, sizeof(from_buf), from);
    const size_t to_len   = fullpath(to_buf, sizeof(to_buf), to);

#if MD_PLATFORM_WINDOWS
    (void)from_len;
    (void)to_len;
    //MD_LOG_DEBUG("rel_from: '%s'", from_buf);
    //MD_LOG_DEBUG("rel_to:   '%s'", to_buf);
    
    success = PathRelativePathTo(out_buf, from_buf, FILE_ATTRIBUTE_NORMAL, to_buf, FILE_ATTRIBUTE_NORMAL);
    len = strnlen(out_buf, out_cap);
    replace_char(out_buf, len, '\\', '/');
#elif MD_PLATFORM_UNIX

    // Find the common base
    str_t can_from = {0};
    str_t can_to   = {to_buf, to_len};

    if (!extract_folder_path(&can_from, (str_t){from_buf, from_len})) {
    	MD_LOG_ERROR("Failed to extract folder path from '" STR_FMT "'", STR_ARG(from));
		return 0;
    }

    //MD_LOG_DEBUG("rel_from: '%.*s'", (int)can_from.len, can_from.ptr);
    //MD_LOG_DEBUG("rel_to:   '%.*s'", (int)can_to.len, can_to.ptr);
    
    size_t count = str_count_equal_chars(can_from, can_to);
    success = count > 0;

    if (success) {
        str_t rel_from = str_substr(can_from, count, SIZE_MAX);
        str_t rel_to   = str_substr(can_to,   count, SIZE_MAX);

        // Count number of folders as N in from and add N times '../'
        size_t folder_count = str_count_occur_char(rel_from, '/');
        if (folder_count) {
            while (folder_count-- > 0) {
                len += snprintf(out_buf + len, out_cap - len, "../");
            }
        } else {
            len += snprintf(out_buf + len, out_cap - len, "./");
        }

        len += snprintf(out_buf + len, out_cap - len, "%.*s", (int)rel_to.len, rel_to.ptr);
    }
#else
    ASSERT(false);
#endif
    if (!success) {
        MD_LOG_ERROR("Failed to extract relative path.");
    }
    return len;
}
    

str_t md_path_make_relative(str_t from, str_t to, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    char  rel_buf[MD_MAX_PATH];
    size_t len = md_path_write_relative(rel_buf, sizeof(rel_buf), from, to);
    if (!len) {
        return (str_t){0};
    }
    return str_copy_cstrn(rel_buf, len, alloc);
}

bool md_path_is_valid(str_t path) {
#if MD_PLATFORM_WINDOWS
    path = str_copy(path, md_temp_allocator);
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
    path = str_copy(path, md_temp_allocator);
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

md_file_o* md_file_open(str_t filename, uint32_t in_flags) {
    if (!filename.ptr || !filename.len) return NULL;

    if (in_flags == 0) {
        MD_LOG_ERROR("No flags passed to file_open, expected at least MD_FILE_READ or MD_FILE_WRITE");
        return NULL;
    }

#if MD_PLATFORM_WINDOWS
    wchar_t w_file[MD_MAX_PATH] = {0};
    const int w_file_len = MultiByteToWideChar(CP_UTF8, 0, filename.ptr, (int)filename.len, w_file, ARRAY_SIZE(w_file));
    if (w_file_len >= ARRAY_SIZE(w_file)) {
        MD_LOG_ERROR("File path exceeds limit!");
        return NULL;
    }
    w_file[w_file_len] = L'\0';

    DWORD desired_access = 0;
    DWORD share_mode = 0;
    LPSECURITY_ATTRIBUTES security_attributes = NULL;
    DWORD creation_disposition = 0;
    DWORD flags_and_attributes = FILE_ATTRIBUTE_NORMAL;
    HANDLE template_file = NULL;

    if (in_flags & MD_FILE_READ) {
		desired_access |= GENERIC_READ;
        share_mode = FILE_SHARE_READ;
		creation_disposition = OPEN_EXISTING;
    }
    if (in_flags & MD_FILE_WRITE) {
        desired_access |= GENERIC_WRITE;
        if (in_flags & MD_FILE_APPEND) {
            creation_disposition = CREATE_NEW;
		} else {
			creation_disposition = CREATE_ALWAYS;
        }
    }

    HANDLE file = CreateFileW(w_file, desired_access, share_mode, security_attributes, creation_disposition, flags_and_attributes, template_file);
    if (file == INVALID_HANDLE_VALUE) {
		MD_LOG_ERROR("Failed to open file");
		return NULL;
	}

    return (md_file_o*)file;

#elif MD_PLATFORM_UNIX
    int flags = 0;

    if (in_flags == MD_FILE_READ) {
        flags |= O_RDONLY;
	}
    if (in_flags & MD_FILE_WRITE) {
		if (in_flags & MD_FILE_APPEND) {
			flags |= O_APPEND;
		} else {
			flags |= O_CREAT;
		}
	}

    mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;

    // Ensure zero termination of path
    str_t path = str_copy(filename, md_temp_allocator);
    int result = open(path.ptr, flags, mode);

    if (result == -1) {
		MD_LOG_ERROR("Failed to open file");
		return NULL;
	}

    uint64_t fd = (uint64_t)result;
	return (md_file_o*)fd;

#else
    const char* mode = 0;

    if (in_flags & MD_FILE_APPEND) {
        mode = "ab";
    } else if (in_flags & MD_FILE_WRITE) {
        mode = "wb";
    } else if (in_flags & MD_FILE_READ) {
		mode = "rb";
	} else {
        MD_LOG_ERROR("Invalid combination of file access flags!");
        return NULL;
    }

    char file[MD_MAX_PATH] = "";
    strncpy(file, filename.ptr, ARRAY_SIZE(file)-1);
    return (md_file_o*)fopen(file, mode);
#endif
}

bool md_file_close(md_file_o* file) {
#if MD_PLATFORM_WINDOWS
    return CloseHandle((HANDLE*)file) == 0;
#elif MD_PLATFORM_UNIX
    uint64_t fd = (uint64_t)file;
    if (close((int)fd) == 0) {
        return true;
    }
    MD_LOG_ERROR("Failed to close file");
    return false;
#else
    return fclose((FILE*)file) == 0;
#endif
}

int64_t md_file_tell(md_file_o* file) {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER pos  = {0};
    LARGE_INTEGER dist = {0};
    if (SetFilePointerEx((HANDLE*)file, dist, &pos, FILE_CURRENT)) {
        return pos.QuadPart;
    }
    MD_LOG_ERROR("Failed to get file position");
    print_windows_error();

    return 0;
#elif MD_PLATFORM_UNIX
    uint64_t fd = (uint64_t)file;
    return lseek((int)fd, 0, SEEK_CUR);
#else
    return ftello((FILE*)file);
#endif
}

bool md_file_seek(md_file_o* file, int64_t offset, md_file_seek_origin_t origin) {
#if MD_PLATFORM_WINDOWS
    DWORD move_method = 0;
    switch (origin) {
        case MD_FILE_BEG:
        case MD_FILE_CUR:
        case MD_FILE_END:
			move_method = (DWORD)origin;
        break;
        default:
            MD_LOG_ERROR("Invalid seek origin value!");
			return false;
    }

    LARGE_INTEGER dist = {.QuadPart = offset};
    if (SetFilePointerEx((HANDLE*)file, dist, NULL, move_method)) {
        return true;
    }
    MD_LOG_ERROR("Failed to seek file");
    print_windows_error();

    return false;
#elif MD_PLATFORM_UNIX
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
    uint64_t fd = (uint64_t)file;
	return lseek((int)fd, offset, o) != -1;
#else
    if (!file) {
        MD_LOG_ERROR("File handle was NULL");
        return false;
    }

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
    return fseeko((FILE*)file, offset, o) == 0;
#endif
}

// Creates a mapping of a file in memory
md_file_mapping_o* md_file_mem_map(md_file_o* file) {
#if MD_PLATFORM_WINDOWS
    HANDLE handle = CreateFileMapping((HANDLE)file, NULL, PAGE_READONLY, 0, 0, NULL);
    return handle;
#elif MD_PLATFORM_UNIX
    uint64_t fd = (uint64_t)file;
    return (md_file_mapping_o*)fd;
#endif
}

bool md_file_mem_unmap(md_file_mapping_o* mapping) {
#if MD_PLATFORM_WINDOWS
    if (!CloseHandle((HANDLE)mapping)) {
        MD_LOG_ERROR("Failed to unmap file");
		print_windows_error();
		return false;
    }
    return true;
#elif MD_PLATFORM_UNIX
    return mapping != 0;
#else
	ASSERT(false);
#endif
}

// Creates a view of a file in memory
const char* md_file_mem_map_view(md_file_mapping_o* mapping, size_t offset, size_t length) {
#if MD_PLATFORM_WINDOWS
    size_t alloc_gran = win32_allocation_granularity();
    const int64_t aligned_offset = ROUND_DOWN(offset, alloc_gran);
    const int64_t length_to_map = (int64_t)offset - aligned_offset + (int64_t)length;
    const char* ptr = MapViewOfFile((HANDLE)mapping, FILE_MAP_READ, (DWORD)(aligned_offset >> 32), (DWORD)(aligned_offset & 0xffffffff), length_to_map);
    if (!ptr) {
        MD_LOG_ERROR("Failed to map view of file");
		print_windows_error();
        return NULL;
    }
    return ptr + (int64_t)offset - aligned_offset;
#elif MD_PLATFORM_UNIX
    size_t page_size = md_vm_page_size();
    const int64_t aligned_offset = ROUND_DOWN(offset, page_size);
    const int64_t length_to_map = (int64_t)offset - aligned_offset + (int64_t)length;
    uint64_t fd = (uint64_t)mapping;
	const char* ptr = mmap(0, length_to_map, PROT_READ, MAP_SHARED, (int)fd, aligned_offset);
	if (ptr == MAP_FAILED) {
		MD_LOG_ERROR("Failed to map view of file");
		return NULL;
	}
    return ptr + (int64_t)offset - aligned_offset;
#else
    ASSERT(false);
#endif
}

bool md_file_mem_unmap_view(const char* addr, size_t size) {
#if MD_PLATFORM_WINDOWS
    (void)size;
    size_t alloc_gran = win32_allocation_granularity();
    const uint64_t aligned_beg = ROUND_DOWN((uint64_t)addr, alloc_gran);
    void* beg = (void*)aligned_beg;
    if (!UnmapViewOfFile(beg)) {
		MD_LOG_ERROR("Failed to unmap view of file");
        print_windows_error();
        return false;
    }
    return true;
#elif MD_PLATFORM_UNIX
    size_t page_size = md_vm_page_size();
	const uint64_t aligned_beg = ROUND_DOWN((uint64_t)addr, page_size);
    const uint64_t aligned_end = ROUND_UP((uint64_t)addr + size, page_size);
    void* beg = (void*)aligned_beg;
    const uint64_t len = aligned_end - aligned_beg;
    return munmap(beg, len) == 0;
#else
    ASSERT(false);
#endif
}

size_t md_file_size(md_file_o* file) {
    ASSERT(file);
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER size = {0};
    if (GetFileSizeEx((HANDLE*)file, &size)) {
    	return (size_t)size.QuadPart;
    }
    MD_LOG_ERROR("Failed to get file size");
    print_windows_error();
    return 0;
#elif MD_PLATFORM_UNIX
    struct stat st;
    uint64_t fd = (uint64_t)file;
    if (fstat((int)fd, &st) == 0) {
    	return (size_t)st.st_size;
    }
    MD_LOG_ERROR("Failed to get file size");
    return 0;
#else
    if (!file) {
        MD_LOG_ERROR("File handle was NULL");
        return 0;
    }
    int64_t cur = md_file_tell(file);
    md_file_seek(file, 0, SEEK_END);
    int64_t end = md_file_tell(file);
    md_file_seek(file, cur, SEEK_SET);
    return (size_t)end;
#endif
}

// Returns the number of successfully written/read bytes
size_t md_file_read(md_file_o* file, void* in_ptr, size_t num_bytes) {
    ASSERT(file);
#if MD_PLATFORM_WINDOWS
    char* ptr = (char*)in_ptr;
    while (num_bytes > 0) {
        DWORD bytes_to_read = (DWORD)MIN(num_bytes, UINT32_MAX);
        DWORD bytes_read = 0;
        if (!ReadFile((HANDLE*)file, ptr, bytes_to_read, &bytes_read, NULL)) {
			MD_LOG_ERROR("Failed to read from file");
            print_windows_error();
			return 0;
		}
        num_bytes -= bytes_read;
        ptr += bytes_read;
        if (bytes_read < bytes_to_read) {
        	break;
        }
    }
    ASSERT(ptr >= (char*)in_ptr);
    return (size_t)(ptr - (char*)in_ptr);
#elif MD_PLATFORM_UNIX
    uint64_t fd = (uint64_t)file;
    char* ptr = (char*)in_ptr;
    while (num_bytes > 0) {
        size_t bytes_to_read = MIN(num_bytes, 0x7ffff000);
    	ssize_t bytes_read = read((int)fd, ptr, bytes_to_read);
    	if (bytes_read == -1) {
            MD_LOG_ERROR("Failed to read from file");
			return 0;
		}
        num_bytes -= bytes_read;
		ptr += bytes_read;
        if ((size_t)bytes_read < bytes_to_read) {
            break;
        }
	}
    ASSERT(ptr >= (char*)in_ptr);
    return (size_t)(ptr - (char*)in_ptr);
#else
    return fread(ptr, 1, num_bytes, (FILE*)file);
#endif
}

size_t md_file_write(md_file_o* file, const void* in_ptr, size_t num_bytes) {
    ASSERT(file);
#if MD_PLATFORM_WINDOWS
    char* ptr = (char*)in_ptr;
    while (num_bytes > 0) {
    	DWORD bytes_to_write = (DWORD)MIN(num_bytes, UINT32_MAX);
    	DWORD bytes_written = 0;
    	if (!WriteFile((HANDLE*)file, ptr, bytes_to_write, &bytes_written, NULL)) {
            MD_LOG_ERROR("Failed to write to file");
            print_windows_error();
            goto done;
        }
        num_bytes -= bytes_written;
        ptr += bytes_written;
        if (bytes_written < bytes_to_write) {
            break;
        }
    }
done:
    ASSERT(ptr >= (char*)in_ptr);
    return (size_t)(ptr - (char*)in_ptr);
#elif MD_PLATFORM_UNIX
    uint64_t fd = (uint64_t)file;
    char* ptr = (char*)in_ptr;
    while (num_bytes > 0) {
        size_t bytes_to_write = MIN(num_bytes, 0x7ffff000);
        ssize_t bytes_written = write((int)fd, ptr, bytes_to_write);
        if (bytes_written < 0) {
			MD_LOG_ERROR("Failed to write to file");
            goto done;
        }
        MD_LOG_DEBUG("bytes_written: %zu", (size_t)bytes_written);
        num_bytes -= bytes_written;
        ptr += bytes_written;
        if ((size_t)bytes_written < bytes_to_write) {
            break;
        }
    }
done:
    ASSERT(ptr >= (char*)in_ptr);
    return (size_t)(ptr - (char*)in_ptr);
#else 
    return fwrite(ptr, 1, num_bytes, (FILE*)file);
#endif
}

size_t md_file_printf(md_file_o* file, const char* format, ...) {
    ASSERT(file);
    char buf[4096];

    va_list args;
    va_start (args, format);
    int res = vsnprintf(buf, sizeof(buf), format, args);
    va_end (args);
    if (res < 0) {
		MD_LOG_ERROR("Failed to write to file");
		return 0;
	}
    size_t len = (size_t)res;
    return md_file_write(file, buf, len);
}


// ### TIME ###

md_timestamp_t md_time_current() {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER t;
    QueryPerformanceCounter(&t);
    return t.QuadPart;
#elif MD_PLATFORM_UNIX
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
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
    return (double)t;
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

size_t md_os_physical_ram(void) {
#if MD_PLATFORM_WINDOWS
    MEMORYSTATUSEX status = {
        .dwLength = sizeof(MEMORYSTATUSEX),
    };
    GlobalMemoryStatusEx(&status);
    return (size_t)status.ullTotalPhys;
#elif MD_PLATFORM_UNIX
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return (size_t)(pages * page_size);
#else
    ASSERT(false);
    return 0;
#endif
}

size_t md_os_num_processors(void) {
#if MD_PLATFORM_WINDOWS
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    int num_cpu = sysinfo.dwNumberOfProcessors;
    return (size_t)num_cpu;
#elif MD_PLATFORM_UNIX
    int num_cpu = sysconf(_SC_NPROCESSORS_ONLN);
    return (size_t)num_cpu;
#else
    ASSERT(false);
    return 0;
#endif
}

static size_t page_size = 0;
static size_t allocation_granularity = 0;

size_t md_vm_page_size(void) {
    if (!page_size) {
#if MD_PLATFORM_WINDOWS
        SYSTEM_INFO info;
        GetSystemInfo(&info);
        page_size = info.dwPageSize;
        allocation_granularity = info.dwAllocationGranularity;
#elif MD_PLATFORM_UNIX
        page_size = sysconf(_SC_PAGE_SIZE);
#else
        ASSERT(false);
#endif
    }
    return page_size;
}

void* md_vm_reserve(size_t size) {
    const size_t gb_snapped_size = ALIGN_TO(size, GIGABYTES(1));
#if MD_PLATFORM_WINDOWS
    return VirtualAlloc(0, gb_snapped_size, MEM_RESERVE, PAGE_NOACCESS);
#elif MD_PLATFORM_UNIX
    void* result = mmap(0, gb_snapped_size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    ASSERT(result != MAP_FAILED);
    return result;
#else
    ASSERT(false);
#endif
}

void md_vm_release(void* ptr, size_t size) {
#if MD_PLATFORM_WINDOWS
    (void)size;
    VirtualFree(ptr, 0, MEM_RELEASE);
#elif MD_PLATFORM_UNIX
    size_t gb_snapped_size = ALIGN_TO(size, GIGABYTES(1));
    int result;
    //result = madvise(ptr, -1, MADV_DONTNEED);
    //ASSERT(result == 0);
    result = munmap(ptr, gb_snapped_size);
    ASSERT(result == 0);
    (void)result;
#else
    ASSERT(false);
#endif
}

void md_vm_commit(void* ptr, size_t size) {
    const size_t page_snapped_size = ALIGN_TO(size, md_vm_page_size());
#if MD_PLATFORM_WINDOWS
    void* result = VirtualAlloc(ptr, page_snapped_size, MEM_COMMIT, PAGE_READWRITE);
    (void)result;
    ASSERT(result != NULL);
#elif MD_PLATFORM_UNIX
    int result = mprotect(ptr, page_snapped_size, PROT_READ | PROT_WRITE);
    ASSERT(result == 0);
    (void)result;
#else
    ASSERT(false);
#endif
}

void md_vm_decommit(void* ptr, size_t size) {
#if MD_PLATFORM_WINDOWS
#if MD_COMPILER_MSVC
#   pragma warning(suppress : 6250)
#endif
    VirtualFree(ptr, size, MEM_DECOMMIT);
#elif MD_PLATFORM_UNIX
    int result;
    size_t page_snapped_size = ALIGN_TO(size, md_vm_page_size());
    result = mprotect(ptr, page_snapped_size, PROT_NONE);
    ASSERT(result == 0);
    result = madvise(ptr, page_snapped_size, MADV_DONTNEED);
    ASSERT(result == 0);
    (void)result;
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

void md_thread_sleep(size_t milliseconds) {
#if MD_PLATFORM_WINDOWS
    Sleep((DWORD)milliseconds);
#elif MD_PLATFORM_UNIX
    usleep(milliseconds * 1000);
#else
    ASSERT(false);
#endif
}

// ### MUTEX ###

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



// ### SEMAPHORE ###

bool md_semaphore_init(md_semaphore_t* semaphore, size_t initial_count) {
    ASSERT(initial_count < INT32_MAX);
#if MD_PLATFORM_WINDOWS
    semaphore->_data[0] = CreateSemaphoreA(NULL, (LONG)initial_count, MAXLONG, NULL);
    return semaphore != NULL;
#elif MD_PLATFORM_LINUX
    return sem_init((sem_t*)semaphore, 0, (uint32_t)initial_count) == 0;
#elif MD_PLATFORM_OSX
    semaphore_t sema;
    mach_port_t self = mach_task_self();
    kern_return_t ret = semaphore_create(self, &sema, SYNC_POLICY_FIFO, (int)initial_count);
    if (ret != KERN_SUCCESS) {
        MD_LOG_ERROR("Failed to initialize semaphore");
    }
    MEMCPY(semaphore, &sema, sizeof(semaphore_t));
    return ret == KERN_SUCCESS;
#endif
}

bool md_semaphore_destroy(md_semaphore_t* semaphore) {
#if MD_PLATFORM_WINDOWS
    return CloseHandle((HANDLE)semaphore->_data[0]);
#elif MD_PLATFORM_LINUX
    return sem_destroy((sem_t*)semaphore) == 0;
#elif MD_PLATFORM_OSX
    semaphore_t sema;
    MEMCPY(&sema, semaphore, sizeof(semaphore_t));
    mach_port_t self = mach_task_self();
    return semaphore_destroy(self, sema) == KERN_SUCCESS;
#endif
}

bool md_semaphore_aquire(md_semaphore_t* semaphore) {
#if MD_PLATFORM_WINDOWS
    return WaitForSingleObjectEx((HANDLE)semaphore->_data[0], INFINITE, FALSE) == WAIT_OBJECT_0;
#elif MD_PLATFORM_LINUX
    return sem_wait((sem_t*)semaphore) == 0;
#elif MD_PLATFORM_OSX
    semaphore_t sema;
    MEMCPY(&sema, semaphore, sizeof(semaphore_t));
    return semaphore_wait(sema) == KERN_SUCCESS;
#endif
}

bool md_semaphore_try_aquire(md_semaphore_t* semaphore) {
#if MD_PLATFORM_WINDOWS
    return WaitForSingleObjectEx((HANDLE)semaphore->_data[0], 0, FALSE) == WAIT_OBJECT_0;
#elif MD_PLATFORM_LINUX
    return sem_trywait((sem_t*)semaphore) == 0;
#elif MD_PLATFORM_OSX
    semaphore_t sema;
    MEMCPY(&sema, semaphore, sizeof(semaphore_t));
    mach_timespec_t mts;
    mts.tv_sec = 0;
    mts.tv_nsec = 0;
    return semaphore_timedwait(sema, mts) == KERN_SUCCESS;
#endif
}

bool md_semaphore_query_count(md_semaphore_t* semaphore, size_t* count) {
    ASSERT(semaphore);
    ASSERT(count);
#if MD_PLATFORM_WINDOWS
    HMODULE hntdll = GetModuleHandle("ntdll.dll");
    if (hntdll) {
        _NtQuerySemaphore NtQuerySemaphore = (_NtQuerySemaphore)GetProcAddress(hntdll, "NtQuerySemaphore");
        if (NtQuerySemaphore) {
            SEMAPHORE_BASIC_INFORMATION basic_info = {0};
            NTSTATUS status = NtQuerySemaphore((HANDLE)semaphore->_data[0], 0 /*SemaphoreBasicInformation*/,  &basic_info, sizeof (SEMAPHORE_BASIC_INFORMATION), NULL);
            if (status == ERROR_SUCCESS) {
                *count = basic_info.CurrentCount;
                return true;
            }
        }
    }
    return false;
#elif MD_PLATFORM_LINUX
    int val;
    if (sem_getvalue((sem_t*)semaphore, &val) != 0) {
        MD_LOG_ERROR("Failed to query semaphore value");
        return false;
	}
    *count = (size_t)val;
	return true;
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
    semaphore_t sema;
    MEMCPY(&sema, semaphore, sizeof(semaphore_t));
    return semaphore_signal(sema) == KERN_SUCCESS;
#endif
}

// Common high level Semaphore operations

md_semaphore_t md_semaphore_create(size_t initial_count) {
    md_semaphore_t semaphore;
    md_semaphore_init(&semaphore, initial_count);
    return semaphore;
}

bool md_semaphore_try_aquire_n(md_semaphore_t* semaphore, size_t count) {
    ASSERT(count > 0);
    size_t aquired_count = 0;
    for (size_t i = 0; i < count; ++i) {
        if (!md_semaphore_try_aquire(semaphore)) {
            break;
        }
        aquired_count += 1;
    }
    if (aquired_count == count) {
        return true;
    }
    // If we did not manage to aquire the sufficient count, we release it again.
    md_semaphore_release_n(semaphore, aquired_count);
    return false;
}

bool md_semaphore_release_n(md_semaphore_t* semaphore, size_t count) {
#if MD_PLATFORM_WINDOWS
    return ReleaseSemaphore((HANDLE)semaphore->_data[0], (LONG)count, NULL);
#else
    for (size_t i = 0; i < count; ++i) {
        if (!md_semaphore_release(semaphore)) {
            return false;
        }
    }
    return true;
#endif
}
