#ifndef _MD_FILE_
#define _MD_FILE_

#include <stdint.h>
#include <stdbool.h>

// The only reason that this file exists is to overcome the flaws of the C standard library's file functions:
// - Support utf-8 filenames on windows
// - Support non zero terminated strings through ptr + len
// - Support 64-bit offsets for seek and tell of large files

#ifndef SEEK_SET
#define SEEK_SET 0
#endif
#ifndef SEEK_CUR
#define SEEK_CUR 1
#endif
#ifndef SEEK_END
#define SEEK_END 2
#endif

// This is just same as FILE* which is just a void*
typedef struct md_file md_file;
struct md_file;

// Similar interface as fopen, but explicit length and supports utf-8 encoding
md_file* md_file_open(const char* filename, uint32_t filename_len, const char* mode);
void       md_file_close(md_file* file);

// supports 64-bit offsets
int64_t md_file_tell(md_file* file);
bool    md_file_seek(md_file* file, int64_t offset, int origin);

uint64_t md_file_size(md_file* file);

// Returns the number of successfully written/read bytes
// [NOTE] The order of the arguments compared to fread/fwrite
uint64_t md_file_read(md_file* file, void* ptr, uint64_t num_bytes);
uint64_t md_file_write(md_file* file, const void* ptr, uint64_t num_bytes);

#endif