#include "md_sync.h"
#include "md_compiler.h"
#include "md_platform.h"

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
md_thread_t* md_thread_create(md_thread_func func, const char* name, void* user_data) {
	(void)name;
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

#else

#include <pthread.h>

md_thread_t* md_thread_create(md_thread_func func, const char* name, void* user_data) {
	pthread_t thread;
	pthread_create(&thread, NULL, (void* (*)(void*))fn, udata);
#if MD_PLATFORM_OSX
	if (name) pthread_setname_np(thread, name);
#else
	(void)name;
#endif
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