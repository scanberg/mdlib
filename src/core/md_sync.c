#include "md_sync.h"
#include "md_compiler.h"
#include "md_platform.h"
#include "md_log.h"
#include "md_common.h"

#if MD_PLATFORM_WINDOWS

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

STATIC_ASSERT(sizeof(CRITICAL_SECTION) <= sizeof(md_mutex_t), "Win32 CRITICAL_SECTION does not fit into md_mutex_t!");

// Thread
md_thread_t* md_thread_create(md_thread_func func, void* user_data) {
	DWORD unused;
	HANDLE id = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)func, user_data, 0, &unused);
	return (md_thread_t*)id;
}

void md_thread_detach(md_thread_t* thread) {
	CloseHandle((HANDLE)thread);
}

int md_thread_join(md_thread_t* thread) {
	WaitForSingleObject((HANDLE)thread, INFINITE);
	CloseHandle((HANDLE)thread);
	return 1;
}

md_thread_id_t md_thread_get_id(md_thread_t* thread) {
	return GetThreadId((HANDLE)thread);
}

md_thread_id_t md_thread_id(void) {
	return GetCurrentThreadId();
}

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
	semaphore->_id = CreateSemaphoreA(NULL, (LONG)initial_count, MAXLONG, NULL);
	return semaphore->_id != NULL;
}

md_semaphore_t md_semaphore_create(int32_t initial_count) {
	md_semaphore_t semaphore;
	md_semaphore_init(&semaphore, initial_count);
	return semaphore;
}

static inline bool semaphore_wait(md_semaphore_t* semaphore, DWORD milliseconds) {
	return WaitForSingleObjectEx(semaphore->_id, milliseconds, FALSE) == WAIT_OBJECT_0;
}

bool md_semaphore_destroy(md_semaphore_t* semaphore) {
	CloseHandle((HANDLE)semaphore->_id);
}

bool md_semaphore_aquire(md_semaphore_t* semaphore) {
	return semaphore_wait(semaphore, INFINITE);
}

bool md_semaphore_try_aquire(md_semaphore_t* semaphore) {
	return semaphore_wait(semaphore, 0);
}

bool md_semaphore_release(md_semaphore_t* semaphore) {
	return ReleaseSemaphore(semaphore->_id, 1, NULL);
}

#elif MD_PLATFORM_UNIX

#include <pthread.h>

md_thread_t* md_thread_create(md_thread_func fn, void* user_data) {
	pthread_t thread;
	pthread_create(&thread, NULL, (void* (*)(void*))fn, user_data);
	return (md_thread_t*)thread;
}

void md_thread_detach(md_thread_t* thread) {
	pthread_detach((pthread_t)thread);
}

int md_thread_join(md_thread_t* thread) {
	pthread_join((pthread_t)thread, NULL);
	return 1;
}

md_thread_id_t md_thread_get_id(md_thread_t* thread) {
	return (md_thread_id_t)thread;
}

md_thread_id_t md_thread_id(void) {
	return (md_thread_id_t)pthread_self();
}

md_mutex_t md_mutex_create() {
	STATIC_ASSERT(sizeof(pthread_mutex_t) <= sizeof(md_mutex_t), "pthread_mutex_t does not fit into md_mutex_t!");
	md_mutex_t mutex;
	pthread_mutex_init((pthread_mutex_t*)&mutex, NULL);
	return mutex;
}

bool md_mutex_destroy(md_mutex_t* mutex) {
	pthread_mutex_destroy((pthread_mutex_t*)mutex);
	return true;
}

bool md_mutex_lock(md_mutex_t* mutex) {
	pthread_mutex_lock((pthread_mutex_t*)mutex);
	return true;
}

bool md_mutex_try_lock(md_mutex_t* mutex) {
	return pthread_mutex_trylock((pthread_mutex_t*)mutex);
}

bool md_mutex_unlock(md_mutex_t* mutex) {
	pthread_mutex_unlock((pthread_mutex_t*)mutex);
	return true;
}

#endif

#if MD_PLATFORM_LINUX

#elif MD_PLATFORM_OSX
// MacOS deprecated pthreads semaphores
#include <mach/mach_init.h>
#include <mach/task.h>
#include <mach/semaphore.h>

// Semaphore
bool md_semaphore_init(md_semaphore_t* semaphore, int32_t initial_count) {
	ASSERT(semaphore);
	semaphore_t sema;
	mach_port_t self = mach_task_self();
	kern_return_t ret = semaphore_create(self, &sema, SYNC_POLICY_FIFO, initial_count);
	if (ret != KERN_SUCCESS) {
		md_print(MD_LOG_TYPE_ERROR, "Failed to initialize semaphore");
	}
	semaphore->_id = (void*)sema;
	return ret == KERN_SUCCESS;
}

md_semaphore_t md_semaphore_create(int32_t initial_count) {
	md_semaphore_t semaphore;
	md_semaphore_init(&semaphore, initial_count);
	return semaphore;
}

bool md_semaphore_destroy(md_semaphore_t* semaphore) {
	mach_port_t self = mach_task_self();
	return semaphore_destroy(self, (semaphore_t)semaphore->_id) == KERN_SUCCESS;
}

bool md_semaphore_aquire(md_semaphore_t* semaphore) {
	return semaphore_wait((semaphore_t)semaphore->_id) == KERN_SUCCESS;
}

bool md_semaphore_try_aquire(md_semaphore_t* semaphore) {
	mach_timespec_t mts;
	mts.tv_sec = 0;
	mts.tv_nsec = 0;
	return semaphore_timedwait((semaphore_t)semaphore->_id, mts) == KERN_SUCCESS;
}

bool md_semaphore_release(md_semaphore_t* semaphore) {
	return semaphore_signal((semaphore_t)semaphore->_id) == KERN_SUCCESS;
}

#endif