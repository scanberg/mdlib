#include <stdint.h>
#include <stdbool.h>
#include <core/md_str.h>

typedef struct md_file_o md_file_o;

typedef enum md_file_flags_t {
    MD_FILE_READ   = 1,
    MD_FILE_WRITE  = 2,
    MD_FILE_APPEND = 4,
    MD_FILE_BINARY = 8
} md_file_flags_t;

typedef enum md_file_seek_origin_t {
    MD_FILE_BEG,
    MD_FILE_CUR,
    MD_FILE_END
} md_file_seek_origin_t;

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

md_file_o* md_file_open(str_t filename, uint32_t flags);
void md_file_close(md_file_o* file);

int64_t md_file_tell(md_file_o* file);
bool md_file_seek(md_file_o* file, int64_t offset, md_file_seek_origin_t origin);
int64_t md_file_size(md_file_o* file);

// Reads a line from the file by searching for the first occurrence of '\n'.
// Returns the number of characters read (up to cap-1)
int64_t md_file_read_line(md_file_o* file, char* buf, int64_t cap);
int64_t md_file_read_lines(md_file_o* file, char* buf, int64_t cap, int64_t line_count);

int64_t md_file_read(md_file_o* file, void* ptr, int64_t num_bytes);
int64_t md_file_write(md_file_o* file, const void* ptr, int64_t num_bytes);

int64_t md_file_printf(md_file_o* file, const char* format, ...);

#ifdef __cplusplus
}
#endif // __cplusplus
