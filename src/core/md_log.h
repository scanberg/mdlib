#ifndef _MD_LOGGER_H_
#define _MD_LOGGER_H_

#ifdef __cplusplus
extern "C" {
#endif

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

void md_add_logger(const md_logger_i* logger);
void md_remove_logger(const md_logger_i* logger);
int  md_print(enum md_log_type_t log_type, const char* msg);
int  md_printf(enum md_log_type_t log_type, const char* format, ...);

extern struct md_logger_i* default_logger;

#ifdef __cplusplus
}
#endif

#endif
