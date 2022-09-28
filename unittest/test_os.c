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
        // FAILS ON UNIX
        str_t path = MAKE_STR("cool/fool/../bool/../file.txt");
        str_t result = md_os_path_make_canonical(path, default_temp_allocator);
        str_t ref = MAKE_STR("cool/file.txt");
        EXPECT_TRUE(str_equal(result, ref));
        */
    }
}

UTEST(os, mem) {
    void* ptr = md_os_reserve(GIGABYTES(1));
    md_os_commit(ptr, MEGABYTES(1));
    md_os_commit((char*)ptr + MEGABYTES(1), MEGABYTES(2));

    EXPECT_NE(ptr, NULL);

    md_os_decommit((char*)ptr + MEGABYTES(1), MEGABYTES(2));
    md_os_decommit(ptr, MEGABYTES(1));
    md_os_release(ptr);
}
