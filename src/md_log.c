#define _CRT_SECURE_NO_WARNINGS

#include "md_log.h"
#include "core/common.h"
#include "core/compiler.h"
#if MD_COMPILER_MSVC
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
#include <debugapi.h>
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>

#define MAX_LOGGERS 64

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


static void _log(struct md_logger_o* inst, enum md_log_type log_type, const char* msg) {
    (void)inst;

    time_t now = time(0);
    struct tm tstruct;
    tstruct = *localtime(&now);

    char time_buf[64];
    uint64_t count = strftime(time_buf, ARRAY_SIZE(time_buf), "[%T]", &tstruct);
    ASSERT(count <= ARRAY_SIZE(time_buf)); // We should never have to trunkate the time string

    fprintf(stderr, "%.*s", (uint32_t)ARRAY_SIZE(time_buf), time_buf);

    switch (log_type) {
    case MD_LOG_TYPE_INFO:
        fprintf(stderr, ANSI_COLOR_GREEN "[info]" ANSI_COLOR_RESET);
        break;
    case MD_LOG_TYPE_DEBUG:
        fprintf(stderr, ANSI_COLOR_MAGENTA "[debug]" ANSI_COLOR_RESET);
        break;
    case MD_LOG_TYPE_ERROR:
        fprintf(stderr, ANSI_COLOR_RED "[error]" ANSI_COLOR_RESET);
        break;
    default:
        break;
    }
    fprintf(stderr, " %s\n",  msg);

#if MD_COMPILER_MSVC
    OutputDebugString(msg);
    OutputDebugString("\n");
#endif
}

static struct md_logger_i _default_logger = {
    NULL,
    _log
};

static int num_loggers = 1;
static md_logger_i* loggers[MAX_LOGGERS] = {&_default_logger};
md_logger_i* default_logger = &_default_logger;

void md_add_logger(const md_logger_i* logger) {
    if (num_loggers < MAX_LOGGERS) {
        loggers[num_loggers++] = logger; // push back
        return;
    }
    md_printf(MD_LOG_TYPE_ERROR, "Failed to add logger, maximum capacity reached. (%d)", MAX_LOGGERS);
}

void md_remove_logger(const md_logger_i* logger) {
    for (int i = 0; i < num_loggers; ++i) {
        if (loggers[i] == logger) {
            loggers[i] = loggers[--num_loggers]; // swap and pop back
            return;
        }
    }
    md_print(MD_LOG_TYPE_ERROR, "Failed to remove logger, logger was not registered.");
}

int md_print(enum md_log_type log_type, const char* msg) {
    for (int i = 0; i < num_loggers; ++i) {
        loggers[i]->log(loggers[i]->inst, log_type, msg);
    }
    return (int)strlen(msg);
}

int md_printf(enum md_log_type log_type, const char* format, ...) {
    char msg[512] = {0};
    va_list args;
    va_start(args, format);
    int res = vsnprintf(msg, ARRAY_SIZE(msg), format, args);
    va_end(args);

    for (int i = 0; i < num_loggers; ++i) {
        loggers[i]->log(loggers[i]->inst, log_type, msg);
    }

    return res;
}
