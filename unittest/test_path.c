#include "utest.h"

#include <core/md_path.h>

#define STR(lit) MD_STR_LIT(lit)

UTEST(path, ext) {
    md_str_t ext = {0};
    EXPECT_TRUE(md_path_ext(&ext, STR("/data/run.1/traj.nc")));
    EXPECT_TRUE(md_str_eq_cstr(ext, "nc"));
    // A dot in a folder is not an extension of a file without one (Amber's mdout, prmtop)
    EXPECT_FALSE(md_path_ext(&ext, STR("/data/run.1/mdout")));
    EXPECT_FALSE(md_path_ext(&ext, STR("C:\\runs\\v1.2\\prmtop")));
    EXPECT_TRUE(md_path_ext(&ext, STR("C:\\runs\\v1.2\\sys.parm7")));
    EXPECT_TRUE(md_str_eq_cstr(ext, "parm7"));

    EXPECT_TRUE(md_path_ext(&ext, STR("file.tar.gz")));
    EXPECT_TRUE(md_str_eq_cstr(ext, "gz"));
    EXPECT_TRUE(md_path_ext(&ext, STR("file.")));
    EXPECT_TRUE(md_str_eq_cstr(ext, ""));
    EXPECT_TRUE(md_path_ext(NULL, STR("file.pdb")));
    EXPECT_FALSE(md_path_ext(&ext, (md_str_t){0}));
}

UTEST(path, filename) {
    md_str_t name = {0};
    EXPECT_TRUE(md_path_filename(&name, STR("c:/folder/file.ext")));
    EXPECT_TRUE(md_str_eq_cstr(name, "file.ext"));
    EXPECT_TRUE(md_path_filename(&name, STR("c:\\folder\\file.ext")));
    EXPECT_TRUE(md_str_eq_cstr(name, "file.ext"));

    // The last separator of either kind ends the folder
    EXPECT_TRUE(md_path_filename(&name, STR("C:/data\\run\\x.pdb")));
    EXPECT_TRUE(md_str_eq_cstr(name, "x.pdb"));
    EXPECT_TRUE(md_path_filename(&name, STR("C:\\data/run/x.pdb")));
    EXPECT_TRUE(md_str_eq_cstr(name, "x.pdb"));

    // No folder: the whole path is the file name
    EXPECT_TRUE(md_path_filename(&name, STR("x.pdb")));
    EXPECT_TRUE(md_str_eq_cstr(name, "x.pdb"));

    EXPECT_TRUE(md_path_filename(&name, STR("c:/folder/")));
    EXPECT_TRUE(md_str_eq_cstr(name, ""));

    EXPECT_FALSE(md_path_filename(&name, (md_str_t){0}));
}

UTEST(path, without_ext) {
    md_str_t out = {0};
    EXPECT_TRUE(md_path_without_ext(&out, STR("c:/folder/file.ext")));
    EXPECT_TRUE(md_str_eq_cstr(out, "c:/folder/file"));
    EXPECT_TRUE(md_path_without_ext(&out, STR("file.tar.gz")));
    EXPECT_TRUE(md_str_eq_cstr(out, "file.tar"));

    // A dot in a folder is not an extension: the path must not be cut inside the folder
    out = STR("untouched");
    EXPECT_FALSE(md_path_without_ext(&out, STR("/data/run.1/mdout")));
    EXPECT_FALSE(md_path_without_ext(&out, STR("C:\\runs\\v1.2\\prmtop")));
    EXPECT_FALSE(md_path_without_ext(&out, STR("./file")));
    EXPECT_FALSE(md_path_without_ext(&out, STR("../file")));
    EXPECT_TRUE(md_str_eq_cstr(out, "untouched"));

    EXPECT_TRUE(md_path_without_ext(&out, STR("/data/run.1/traj.nc")));
    EXPECT_TRUE(md_str_eq_cstr(out, "/data/run.1/traj"));

    // Agrees with md_path_ext: path == without_ext + '.' + ext
    const md_str_t paths[] = { STR("a.b"), STR("dir.x/a.b.c"), STR("C:\\v1.2\\sys.parm7"), STR(".hidden") };
    for (size_t i = 0; i < sizeof(paths) / sizeof(paths[0]); ++i) {
        md_str_t ext, stem;
        EXPECT_TRUE(md_path_ext(&ext, paths[i]));
        EXPECT_TRUE(md_path_without_ext(&stem, paths[i]));
        EXPECT_EQ(stem.len + 1 + ext.len, paths[i].len);
    }
}

UTEST(path, folder) {
    md_str_t folder = {0};
    EXPECT_TRUE(md_path_folder(&folder, STR("c:/folder/file.ext")));
    EXPECT_TRUE(md_str_eq_cstr(folder, "c:/folder/"));
    EXPECT_TRUE(md_path_folder(&folder, STR("c:\\folder\\file.ext")));
    EXPECT_TRUE(md_str_eq_cstr(folder, "c:\\folder\\"));

    // The last separator of either kind ends the folder
    EXPECT_TRUE(md_path_folder(&folder, STR("C:/data\\run\\x.pdb")));
    EXPECT_TRUE(md_str_eq_cstr(folder, "C:/data\\run\\"));
    EXPECT_TRUE(md_path_folder(&folder, STR("C:\\data/run/x.pdb")));
    EXPECT_TRUE(md_str_eq_cstr(folder, "C:\\data/run/"));

    EXPECT_TRUE(md_path_folder(&folder, STR("/")));
    EXPECT_TRUE(md_str_eq_cstr(folder, "/"));

    // A bare file name has no folder, and the output is left as it was
    folder = STR("untouched");
    EXPECT_FALSE(md_path_folder(&folder, STR("x.pdb")));
    EXPECT_FALSE(md_path_folder(&folder, (md_str_t){0}));
    EXPECT_TRUE(md_str_eq_cstr(folder, "untouched"));
}

UTEST(path, folder_and_filename_partition_the_path) {
    const md_str_t paths[] = { STR("a/b/c.d"), STR("a\\b/c"), STR("c"), STR("/"), STR("a/"), STR("C:\\x\\y.z") };
    for (size_t i = 0; i < sizeof(paths) / sizeof(paths[0]); ++i) {
        md_str_t folder = {0}, name = {0};
        md_path_folder(&folder, paths[i]);
        EXPECT_TRUE(md_path_filename(&name, paths[i]));
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
