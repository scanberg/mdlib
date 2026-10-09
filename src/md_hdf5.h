#pragma once

// The one lock every HDF5 call in mdlib is made under.
//
// HDF5 is not built thread safe by default - not by vcpkg, not by most distributions - and its
// global state is shared by every reader in this library, so a lock private to one reader does not
// make it safe: the VeloxChem, TREXIO and H5MD readers, and the attribute providers they leave
// behind, can each be running on a different thread. Holding this one lock around every HDF5 call
// is what serialises them.
//
// While it is held, HDF5's own error printing is silenced: it prints its stack to stderr on every
// failed probe, and probing for optional objects is most of what reading these files is. The
// previous handler is restored on unlock.
//
// Not recursive. Take it at the outermost function that touches HDF5 and nowhere beneath it.

#include <hdf5.h>

typedef struct md_hdf5_lock_t {
    H5E_auto2_t func;
    void*       client_data;
} md_hdf5_lock_t;

#ifdef __cplusplus
extern "C" {
#endif

md_hdf5_lock_t md_hdf5_lock(void);
void           md_hdf5_unlock(md_hdf5_lock_t lock);

#ifdef __cplusplus
}
#endif
