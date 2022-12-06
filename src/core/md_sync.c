#include "md_sync.h"
#include "md_compiler.h"
#include "md_platform.h"
#include "md_log.h"
#include "md_common.h"

#include <string.h>

#if MD_PLATFORM_WINDOWS

//#pragma comment(lib, "ntdll.lib")

#ifndef NOMINMAX
#define NOMINMAX
#endif

#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN
#endif

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include <Windows.h>

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


// Mutex
bool md_mutex_init(md_mutex_t* mutex) {
	return InitializeCriticalSectionAndSpinCount((CRITICAL_SECTION*)mutex, 2000);
}

md_mutex_t md_mutex_create() {
	md_mutex_t mutex;
	md_mutex_init(&mutex);
	return mutex;
}

bool md_mutex_destroy(md_mutex_t* mutex) {
	DeleteCriticalSection((CRITICAL_SECTION*)mutex);
	return true;
}

bool md_mutex_lock(md_mutex_t* mutex) {
	EnterCriticalSection((CRITICAL_SECTION*)mutex);
	return true;
}

bool md_mutex_try_lock(md_mutex_t* mutex) {
	return TryEnterCriticalSection((CRITICAL_SECTION*)mutex);
}

bool md_mutex_unlock(md_mutex_t* mutex) {
	LeaveCriticalSection((CRITICAL_SECTION*)mutex);
	return true;
}

// Semaphore
bool md_semaphore_init(md_semaphore_t* semaphore, int32_t initial_count) {
	semaphore->_data[0] = CreateSemaphoreA(NULL, (LONG)initial_count, MAXLONG, NULL);
	return semaphore != NULL;
}

static inline bool semaphore_wait(md_semaphore_t* semaphore, DWORD milliseconds) {
	return WaitForSingleObjectEx((HANDLE)semaphore->_data[0], milliseconds, FALSE) == WAIT_OBJECT_0;
}

bool md_semaphore_destroy(md_semaphore_t* semaphore) {
	return CloseHandle((HANDLE)semaphore->_data[0]);
}

bool md_semaphore_aquire(md_semaphore_t* semaphore) {
	return semaphore_wait(semaphore, INFINITE);
}

bool md_semaphore_try_aquire(md_semaphore_t* semaphore) {
	return semaphore_wait(semaphore, 0);
}

bool md_semaphore_query_count(md_semaphore_t* semaphore, int32_t* count) {
	ASSERT(semaphore);
	ASSERT(count);

	_NtQuerySemaphore NtQuerySemaphore;
	SEMAPHORE_BASIC_INFORMATION BasicInfo;
	NTSTATUS Status;

	NtQuerySemaphore = (_NtQuerySemaphore)GetProcAddress(GetModuleHandle("ntdll.dll"), "NtQuerySemaphore");

	if (NtQuerySemaphore) {
		Status = NtQuerySemaphore ((HANDLE)semaphore->_data[0], 0 /*SemaphoreBasicInformation*/, 
			&BasicInfo, sizeof (SEMAPHORE_BASIC_INFORMATION), NULL);
		if (Status == ERROR_SUCCESS) {
			*count = BasicInfo.CurrentCount;
			return true;
		}
	}

	return false;
}

bool md_semaphore_release(md_semaphore_t* semaphore) {
	return ReleaseSemaphore((HANDLE)semaphore->_data[0], 1, NULL);
}

#elif MD_PLATFORM_UNIX

#include <pthread.h>

STATIC_ASSERT(sizeof(pthread_mutex_t) <= sizeof(md_mutex_t), "pthread_mutex_t does not fit into md_mutex_t!");

bool md_mutex_init(md_mutex_t* mutex) {
	return pthread_mutex_init((pthread_mutex_t*)mutex, NULL) == 0;
}

md_mutex_t md_mutex_create() {
	md_mutex_t mutex;
	md_mutex_init(&mutex);
	return mutex;
}

bool md_mutex_destroy(md_mutex_t* mutex) {
	return pthread_mutex_destroy((pthread_mutex_t*)mutex) == 0;
}

bool md_mutex_lock(md_mutex_t* mutex) {
	return pthread_mutex_lock((pthread_mutex_t*)mutex) == 0;
}

bool md_mutex_try_lock(md_mutex_t* mutex) {
	return pthread_mutex_trylock((pthread_mutex_t*)mutex) == 0;
}

bool md_mutex_unlock(md_mutex_t* mutex) {
	return pthread_mutex_unlock((pthread_mutex_t*)mutex) == 0;
}

#if MD_PLATFORM_LINUX
#include <semaphore.h>

STATIC_ASSERT(sizeof(sem_t) <= sizeof(md_semaphore_t), "Linux sem_t does not fit into md_semaphore_t!");

// Semaphore
bool md_semaphore_init(md_semaphore_t* semaphore, int32_t initial_count) {
	ASSERT(semaphore);
	return sem_init((sem_t*)semaphore, 0, (uint32_t)initial_count) == 0;
}

bool md_semaphore_destroy(md_semaphore_t* semaphore) {
	ASSERT(semaphore);
	return sem_destroy((sem_t*)semaphore) == 0;
}

bool md_semaphore_aquire(md_semaphore_t* semaphore) {
	ASSERT(semaphore);
	return sem_wait((sem_t*)semaphore) == 0;
}

bool md_semaphore_try_aquire(md_semaphore_t* semaphore) {
	ASSERT(semaphore);
	return sem_trywait((sem_t*)semaphore) == 0;
}

bool md_semaphore_query_count(md_semaphore_t* semaphore, int32_t* count) {
	ASSERT(semaphore);
	ASSERT(count);
	return sem_getvalue((sem_t*)semaphore, count) == 0;
}

bool md_semaphore_release(md_semaphore_t* semaphore) {
	ASSERT(semaphore);
	return sem_post((sem_t*)semaphore) == 0;
}
#endif

#endif

#if MD_PLATFORM_OSX
// MacOS deprecated pthreads semaphores
#include <mach/mach_init.h>
#include <mach/task.h>
#include <mach/semaphore.h>

STATIC_ASSERT(sizeof(semaphore_t) <= sizeof(md_semaphore_t), "MacOS semaphore_t does not fit into md_semaphore_t!");

// Semaphore
bool md_semaphore_init(md_semaphore_t* semaphore, int32_t initial_count) {
	ASSERT(semaphore);
	semaphore_t sema;
	mach_port_t self = mach_task_self();
	kern_return_t ret = semaphore_create(self, &sema, SYNC_POLICY_FIFO, initial_count);
	if (ret != KERN_SUCCESS) {
		md_print(MD_LOG_TYPE_ERROR, "Failed to initialize semaphore");
	}
	semaphore->_data[0] = (void*)(uint64_t)sema;
	return ret == KERN_SUCCESS;
}

bool md_semaphore_destroy(md_semaphore_t* semaphore) {
	mach_port_t self = mach_task_self();
	return semaphore_destroy(self, (semaphore_t)semaphore->_data[0]) == KERN_SUCCESS;
}

bool md_semaphore_aquire(md_semaphore_t* semaphore) {
	return semaphore_wait((semaphore_t)semaphore->_data[0]) == KERN_SUCCESS;
}

bool md_semaphore_try_aquire(md_semaphore_t* semaphore) {
	mach_timespec_t mts;
	mts.tv_sec = 0;
	mts.tv_nsec = 0;
	return semaphore_timedwait((semaphore_t)semaphore->_data[0], mts) == KERN_SUCCESS;
}

bool md_semaphore_release(md_semaphore_t* semaphore) {
	return semaphore_signal((semaphore_t)semaphore->_data[0]) == KERN_SUCCESS;
}

bool md_semaphore_query_count(md_semaphore_t* semaphore, int32_t* count) {
	(void)semaphore;
	(void)count;
	// Not supported
	return false;
}

#endif

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

#if MD_COMPILER_GCC
#pragma GCC diagnostic pop
#endif
