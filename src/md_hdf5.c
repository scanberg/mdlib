#include "md_hdf5.h"

#include <core/md_platform.h>

// Statically initialised rather than created on first use: two readers taking the lock for the first
// time on two threads would otherwise race to initialise it.
#if MD_PLATFORM_WINDOWS
#   ifndef WIN32_LEAN_AND_MEAN
#       define WIN32_LEAN_AND_MEAN
#   endif
#   include <windows.h>
static SRWLOCK hdf5_lock = SRWLOCK_INIT;
static void lock_acquire(void) { AcquireSRWLockExclusive(&hdf5_lock); }
static void lock_release(void) { ReleaseSRWLockExclusive(&hdf5_lock); }
#elif MD_PLATFORM_UNIX
#   include <pthread.h>
static pthread_mutex_t hdf5_lock = PTHREAD_MUTEX_INITIALIZER;
static void lock_acquire(void) { pthread_mutex_lock(&hdf5_lock); }
static void lock_release(void) { pthread_mutex_unlock(&hdf5_lock); }
#else
#   error "md_hdf5: no lock for this platform"
#endif

md_hdf5_lock_t md_hdf5_lock(void) {
    lock_acquire();
    md_hdf5_lock_t lock = {0};
    H5Eget_auto2(H5E_DEFAULT, &lock.func, &lock.client_data);
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
    return lock;
}

void md_hdf5_unlock(md_hdf5_lock_t lock) {
    H5Eset_auto2(H5E_DEFAULT, lock.func, lock.client_data);
    lock_release();
}
