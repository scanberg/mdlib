#ifndef _MOLD_FILE_
#define _MOLD_FILE_

#include <stdint.h>
#include <stdbool.h>

// The only reason that this file exists is to overcome the flaws of the C standard libraries file functions:
// - Support utf-8 filenames on windows
// - Support 64-bit offsets for seek and tell of large files

// This is just same as FILE*
typedef struct mold_file mold_file;
struct mold_file;

// Same interface as fopen, but supports utf-8 encoding
mold_file* mold_file_open(const char* filename, const char* mode);
void       mold_file_close(mold_file* file);

// supports 64-bit offsets
int64_t mold_file_tell(mold_file* file);
bool    mold_file_seek(mold_file* file, int64_t offset, int origin);

// Returns the number of successfully written/read bytes
// Note the change in order of arguments compared to fread/fwrite
uint64_t mold_file_read(mold_file* file, void* ptr, uint64_t num_bytes);
uint64_t mold_file_write(mold_file* file, const void* ptr, uint64_t num_bytes);

#endif