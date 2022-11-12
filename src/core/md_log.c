#include "md_log.h"
#include "core/md_common.h"
#include "core/md_compiler.h"
#include "core/md_platform.h"

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
#include <debugapi.h>
#endif

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#define MAX_LOGGERS 64

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

typedef enum console_color_t {
  COLOR_DEFAULT = 0,
  COLOR_RED,
  COLOR_GREEN,
  COLOR_BLUE,
  COLOR_YELLOW,
  COLOR_MAGENTA,
  COLOR_CYAN,
} console_color_t;

// Taken from Ghoul (thanks Alex!) (https://github.com/OpenSpace/Ghoul/blob/master/src/logging/consolelog.cpp#L151)
static void set_console_text_color(console_color_t color, bool intense) {
#ifdef WIN32
    WORD fg_color = 0;

    switch (color) {
    case COLOR_RED:
        fg_color = FOREGROUND_RED;
        break;
    case COLOR_GREEN:
        fg_color = FOREGROUND_GREEN;
        break;
    case COLOR_BLUE:
        fg_color = FOREGROUND_BLUE;
        break;
    case COLOR_YELLOW:
        fg_color = FOREGROUND_RED | FOREGROUND_GREEN;
        break;
    case COLOR_MAGENTA:
        fg_color = FOREGROUND_RED | FOREGROUND_BLUE;
        break;
    case COLOR_CYAN:
        fg_color = FOREGROUND_GREEN | FOREGROUND_BLUE;
        break;
    case COLOR_DEFAULT:
        fg_color = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE;
        break;
    default:
        ASSERT(false);
    }

    fg_color |= intense ? FOREGROUND_INTENSITY : 0;

    HANDLE console = GetStdHandle(STD_OUTPUT_HANDLE);
    // Get the old color information
    CONSOLE_SCREEN_BUFFER_INFO info = { 0 };
    GetConsoleScreenBufferInfo(console, &info);
    // Or-ing the new foreground color with the old values for the background
    const WORD bg_color = BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED | BACKGROUND_INTENSITY;
    SetConsoleTextAttribute(console, fg_color | (info.wAttributes & bg_color));
#else
(void)intense;
    switch (color) {
    case COLOR_RED:
        fprintf(stderr, ANSI_COLOR_RED);
        break;
    case COLOR_GREEN:
        fprintf(stderr, ANSI_COLOR_GREEN);
        break;
    case COLOR_BLUE:
        fprintf(stderr, ANSI_COLOR_BLUE);
        break;
    case COLOR_YELLOW:
        fprintf(stderr, ANSI_COLOR_YELLOW);
        break;
    case COLOR_MAGENTA:
        fprintf(stderr, ANSI_COLOR_MAGENTA);
        break;
    case COLOR_CYAN:
        fprintf(stderr, ANSI_COLOR_CYAN);
        break;
    case COLOR_DEFAULT:
        fprintf(stderr, ANSI_COLOR_RESET);
        break;
    default:
        ASSERT(false);
    }
    
#endif
}

static void _log(md_logger_o* inst, md_log_type_t log_type, const char* msg) {
    (void)inst;

    time_t now = time(0);
    struct tm tstruct;
    tstruct = *localtime(&now);

    char time_buf[64];
    const uint64_t time_len = strftime(time_buf, sizeof(time_buf), "[%T]", &tstruct);

    fprintf(stderr, "%.*s", (int)time_len, time_buf);

    switch (log_type) {
    case MD_LOG_TYPE_INFO:
        set_console_text_color(COLOR_GREEN, true);
        fprintf(stderr, "[info]:  ");
        break;
    case MD_LOG_TYPE_DEBUG:
        set_console_text_color(COLOR_MAGENTA, true);
        fprintf(stderr, "[debug]: ");
        break;
    case MD_LOG_TYPE_ERROR:
        set_console_text_color(COLOR_RED, true);
        fprintf(stderr, "[error]: ");
        break;
    default:
        ASSERT(false);
        break;
    }
    set_console_text_color(COLOR_DEFAULT, false);

    fprintf(stderr, "%s\n",  msg);

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

void md_logger_add(const md_logger_i* logger) {
    if (num_loggers < MAX_LOGGERS) {
        loggers[num_loggers++] = (md_logger_i*)logger; // push back
        return;
    }
    md_printf(MD_LOG_TYPE_ERROR, "Failed to add logger, maximum capacity reached. (%d)", MAX_LOGGERS);
}

void md_logger_remove(const md_logger_i* logger) {
    for (int i = 0; i < num_loggers; ++i) {
        if (loggers[i] == logger) {
            loggers[i] = loggers[--num_loggers]; // swap and pop back
            return;
        }
    }
    md_print(MD_LOG_TYPE_ERROR, "Failed to remove logger, logger was not registered.");
}

int md_print(md_log_type_t log_type, const char* msg) {
    for (int i = 0; i < num_loggers; ++i) {
        loggers[i]->log(loggers[i]->inst, log_type, msg);
    }
    return (int)strlen(msg);
}

int md_printf(md_log_type_t log_type, const char* format, ...) {
    char buf[4096];

    va_list args;
    va_start(args, format);
    int res = vsnprintf(buf, sizeof(buf), format, args);
    va_end(args);

    for (int i = 0; i < num_loggers; ++i) {
        loggers[i]->log(loggers[i]->inst, log_type, buf);
    }

    return res;
}
