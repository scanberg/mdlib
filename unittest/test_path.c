#include "utest.h"

#include <core/md_path.h>

#define STR(lit) MD_STR_LIT(lit)
#define EXPECT_STR(str, cstr) EXPECT_TRUE(md_str_eq_cstr((str), (cstr)))

UTEST(path, ext) {
    EXPECT_STR(md_path_ext(STR("/data/run.1/traj.nc")), "nc");
    EXPECT_STR(md_path_ext(STR("C:\\runs\\v1.2\\sys.parm7")), "parm7");
    EXPECT_STR(md_path_ext(STR("file.tar.gz")), "gz");

    // Only in the file name: a dot in a folder is not an extension of a file without one (Amber's mdout, prmtop)
    EXPECT_STR(md_path_ext(STR("/data/run.1/mdout")), "");
    EXPECT_STR(md_path_ext(STR("C:\\runs\\v1.2\\prmtop")), "");

    // A leading '.' marks a hidden file, not an extension
    EXPECT_STR(md_path_ext(STR(".bashrc")), "");
    EXPECT_STR(md_path_ext(STR("/home/user/.bashrc")), "");
    EXPECT_STR(md_path_ext(STR(".config.json")), "json");
    EXPECT_STR(md_path_ext(STR("..x")), "x");

    EXPECT_STR(md_path_ext(STR("file.")), "");
    EXPECT_STR(md_path_ext(STR("Makefile")), "");
    EXPECT_STR(md_path_ext((md_str_t){0}), "");
}

UTEST(path, filename) {
    EXPECT_STR(md_path_filename(STR("c:/folder/file.ext")), "file.ext");
    EXPECT_STR(md_path_filename(STR("c:\\folder\\file.ext")), "file.ext");

    // The last separator of either kind ends the folder
    EXPECT_STR(md_path_filename(STR("C:/data\\run\\x.pdb")), "x.pdb");
    EXPECT_STR(md_path_filename(STR("C:\\data/run/x.pdb")), "x.pdb");

    // No folder: the whole path is the file name
    EXPECT_STR(md_path_filename(STR("x.pdb")), "x.pdb");

    EXPECT_STR(md_path_filename(STR("c:/folder/")), "");
    EXPECT_STR(md_path_filename((md_str_t){0}), "");
}

UTEST(path, without_ext) {
    EXPECT_STR(md_path_without_ext(STR("c:/folder/file.ext")), "c:/folder/file");
    EXPECT_STR(md_path_without_ext(STR("file.tar.gz")), "file.tar");
    EXPECT_STR(md_path_without_ext(STR("/data/run.1/traj.nc")), "/data/run.1/traj");

    // No extension: all of it, the folder never cut
    EXPECT_STR(md_path_without_ext(STR("/data/run.1/mdout")), "/data/run.1/mdout");
    EXPECT_STR(md_path_without_ext(STR("C:\\runs\\v1.2\\prmtop")), "C:\\runs\\v1.2\\prmtop");
    EXPECT_STR(md_path_without_ext(STR("./file")), "./file");
    EXPECT_STR(md_path_without_ext(STR("../file")), "../file");
    EXPECT_STR(md_path_without_ext(STR(".bashrc")), ".bashrc");
    EXPECT_STR(md_path_without_ext(STR("file.")), "file");

    // path is without_ext + '.' + ext, whenever there is an extension
    const md_str_t paths[] = { STR("a.b"), STR("dir.x/a.b.c"), STR("C:\\v1.2\\sys.parm7"), STR(".config.json") };
    for (size_t i = 0; i < sizeof(paths) / sizeof(paths[0]); ++i) {
        const md_str_t ext  = md_path_ext(paths[i]);
        const md_str_t stem = md_path_without_ext(paths[i]);
        EXPECT_GT(ext.len, (size_t)0);
        EXPECT_EQ(stem.len + 1 + ext.len, paths[i].len);
        EXPECT_EQ(stem.ptr, paths[i].ptr);
    }
}

UTEST(path, stem) {
    EXPECT_STR(md_path_stem(STR("c:/folder/file.tar.gz")), "file.tar");
    EXPECT_STR(md_path_stem(STR("/data/run.1/mdout")), "mdout");
    EXPECT_STR(md_path_stem(STR("traj.xtc")), "traj");
    EXPECT_STR(md_path_stem(STR("/home/user/.bashrc")), ".bashrc");
    EXPECT_STR(md_path_stem(STR("file.")), "file");
    EXPECT_STR(md_path_stem(STR("c:/folder/")), "");
}

UTEST(path, folder) {
    EXPECT_STR(md_path_folder(STR("c:/folder/file.ext")), "c:/folder/");
    EXPECT_STR(md_path_folder(STR("c:\\folder\\file.ext")), "c:\\folder\\");

    // The last separator of either kind ends the folder
    EXPECT_STR(md_path_folder(STR("C:/data\\run\\x.pdb")), "C:/data\\run\\");
    EXPECT_STR(md_path_folder(STR("C:\\data/run/x.pdb")), "C:\\data/run/");

    EXPECT_STR(md_path_folder(STR("/")), "/");

    // A bare file name has no folder
    EXPECT_STR(md_path_folder(STR("x.pdb")), "");
    EXPECT_STR(md_path_folder((md_str_t){0}), "");
}

// path is folder + filename
UTEST(path, folder_and_filename_partition_the_path) {
    const md_str_t paths[] = { STR("a/b/c.d"), STR("a\\b/c"), STR("c"), STR("/"), STR("a/"), STR("C:\\x\\y.z") };
    for (size_t i = 0; i < sizeof(paths) / sizeof(paths[0]); ++i) {
        const md_str_t folder = md_path_folder(paths[i]);
        const md_str_t name   = md_path_filename(paths[i]);
        EXPECT_EQ(folder.len + name.len, paths[i].len);
        EXPECT_EQ(name.ptr, paths[i].ptr + folder.len);
    }
}

UTEST(path, is_absolute) {
    // Recognized the same way on every platform
    EXPECT_TRUE(md_path_is_absolute(STR("/home/user/file.txt")));
    EXPECT_TRUE(md_path_is_absolute(STR("C:/Users/user/file.txt")));
    EXPECT_TRUE(md_path_is_absolute(STR("C:\\Users\\user\\file.txt")));
    EXPECT_TRUE(md_path_is_absolute(STR("//server/share/file.txt")));
    EXPECT_TRUE(md_path_is_absolute(STR("\\\\server\\share\\file.txt")));
    EXPECT_FALSE(md_path_is_absolute(STR("./file.txt")));
    EXPECT_FALSE(md_path_is_absolute(STR("../dir/file.txt")));
    EXPECT_FALSE(md_path_is_absolute(STR("dir/file.txt")));
    EXPECT_FALSE(md_path_is_absolute(STR("C:file.txt")));
    EXPECT_FALSE(md_path_is_absolute((md_str_t){0}));
}
