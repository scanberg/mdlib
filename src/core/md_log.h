#pragma once

#include <core/md_common.h>

// @TODO: Use a selection macro which picks the correct version to use based on the number of arguments...

#define MD_LOG_INFO(fmt, ...)  md_logf(MD_LOG_TYPE_INFO,  fmt           , ##__VA_ARGS__)
#define MD_LOG_DEBUG(fmt, ...) md_logf(MD_LOG_TYPE_DEBUG, fmt " [%s:%d]", ##__VA_ARGS__, __FILE__, __LINE__)
#define MD_LOG_ERROR(fmt, ...) md_logf(MD_LOG_TYPE_ERROR, fmt " [%s]"   , ##__VA_ARGS__, __func__)

typedef struct md_logger_o md_logger_o;

typedef enum md_log_type_t {
    MD_LOG_TYPE_INFO,
    MD_LOG_TYPE_DEBUG,
    MD_LOG_TYPE_ERROR
} md_log_type_t;

typedef struct md_logger_i {
    struct md_logger_o* inst;
    void (*log)(struct md_logger_o* inst, enum md_log_type_t log_type, const char* msg);
} md_logger_i;

#ifdef __cplusplus
extern "C" {
#endif

void md_logger_add(const md_logger_i* logger);
void md_logger_remove(const md_logger_i* logger);

int  md_log (md_log_type_t log_type, const char* msg);
int  md_logf(md_log_type_t log_type, const char* fmt, ...);

// This is added from the beginning,
// If you do not want to have this, you have to remove it explicitly
extern struct md_logger_i* default_logger;

#ifdef __cplusplus
}
#endif
