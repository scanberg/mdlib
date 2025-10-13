#pragma once

struct md_system_loader_i;

// Utils for reading PDBX/mmCIF files
// https://mmcif.wwpdb.org/
// This is only meant to extract subset of data from the file, not the entire file and is thus
// not a complete implementation of the format.

#ifdef __cplusplus
extern "C" {
#endif

struct md_system_loader_i* md_mmcif_system_loader(void);

#ifdef __cplusplus
}
#endif
