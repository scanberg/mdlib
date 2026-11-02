#pragma once

// Amber NetCDF trajectories and restarts in the NetCDF-4 encoding, which is HDF5 underneath: what md_amber_netcdf.c
// hands a file to when it starts with the HDF5 signature. Built with MD_HDF5 only. Internal to the Amber reader; the
// API is md_amber_netcdf.h.

#include <md_amber_netcdf.h>

#ifdef __cplusplus
extern "C" {
#endif

// The kind of the file from its Conventions, and its atom count. MD_AMBER_NETCDF_HDF5 for an HDF5 file which is not
// Amber NetCDF.
md_amber_netcdf_kind_t md_amber_nc4_file_kind(str_t filename, size_t* out_num_atoms);

bool md_amber_nc4_state_init_from_file(struct md_system_state_t* state, str_t filename, size_t num_atoms);

bool md_amber_nc4_system_publish_run(struct md_system_t* sys, str_t filename, str_t run, uint32_t flags);

#ifdef __cplusplus
}
#endif
