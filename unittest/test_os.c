#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_os.h>

UTEST(os, path) {
    {
        char buf[4096];
        str_t cwd = md_os_current_working_directory();
        snprintf(buf, sizeof(buf), "%.*s/../../", (int)cwd.len, cwd.ptr);
        str_t path = {buf, strnlen(buf, sizeof(buf))};
        str_t result = md_os_path_make_canonical(path, default_temp_allocator);
        printf("result: '%.*s'", (int)result.len, result.ptr);
    }
    {
        /*
        str_t path = make_cst
        r("/cool/fool/../bool/../file.txt");
        str_t result = md_os_path_make_canonical(path, default_temp_allocator);
        str_t ref = make_cstr("/cool/file.txt");
        EXPECT_TRUE(compare_str(result, ref));
        */
    }
}
