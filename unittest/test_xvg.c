#include "utest.h"
#include <string.h>

#include <md_xvg.h>
#include <md_csv.h>
#include <md_system.h>
#include <md_script.h>
#include <core/md_arena_allocator.h>
#include <stdio.h>
#include <core/md_allocator.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_array.h>

UTEST(xvg, rdf) {
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/rdf.xvg");
    md_xvg_t xvg = {0};
    bool result = md_xvg_parse_file(&xvg, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    
    EXPECT_EQ(2,    xvg.num_fields);
    EXPECT_EQ(1024, xvg.num_values);

    EXPECT_NEAR(0.010, xvg.fields[0][10], 1.0e-6f);
    EXPECT_NEAR(0.000, xvg.fields[1][10], 1.0e-6f);
    
    EXPECT_NEAR(0.049, xvg.fields[0][50], 1.0e-6f);
    EXPECT_NEAR(0.000, xvg.fields[1][50], 1.0e-6f);

    ASSERT_EQ(1, xvg.header_info.num_legends);
    EXPECT_TRUE(str_eq(xvg.header_info.legends[0], STR_LIT("OW")));

    EXPECT_TRUE(str_eq(xvg.header_info.title, STR_LIT("Radial distribution")));
    EXPECT_TRUE(str_eq(xvg.header_info.xaxis_label, STR_LIT("r (nm)")));
    EXPECT_TRUE(str_eq(xvg.header_info.yaxis_label, STR_LIT("g(r)")));

    md_xvg_free(&xvg, md_get_heap_allocator());
}

UTEST(xvg, energy) {
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/energy.xvg");
    md_xvg_t xvg = {0};
    bool result = md_xvg_parse_file(&xvg, path, md_get_heap_allocator());
    ASSERT_TRUE(result);

    EXPECT_EQ(5,   xvg.num_fields);
    EXPECT_EQ(356, xvg.num_values);

    //    18.000000  -427.499512  -170.860535  -3487.963135  -2470.888916
    EXPECT_NEAR(   18.000000, xvg.fields[0][9], 1.0e-6f);
    EXPECT_NEAR( -427.499512, xvg.fields[1][9], 1.0e-6f);
    EXPECT_NEAR( -170.860535, xvg.fields[2][9], 1.0e-6f);
    EXPECT_NEAR(-3487.963135, xvg.fields[3][9], 1.0e-6f);
    EXPECT_NEAR(-2470.888916, xvg.fields[4][9], 1.0e-6f);

    //   710.000000  -414.190063  -272.062927  -3277.823975  -2500.806641
    EXPECT_NEAR(  710.000000, xvg.fields[0][355], 1.0e-6f);
    EXPECT_NEAR( -414.190063, xvg.fields[1][355], 1.0e-6f);
    EXPECT_NEAR( -272.062927, xvg.fields[2][355], 1.0e-6f);
    EXPECT_NEAR(-3277.823975, xvg.fields[3][355], 1.0e-6f);
    EXPECT_NEAR(-2500.806641, xvg.fields[4][355], 1.0e-6f);

    ASSERT_EQ(4, xvg.header_info.num_legends);
    EXPECT_TRUE(str_eq(xvg.header_info.legends[0], STR_LIT("Coul-SR:2S29-2S29")));
    EXPECT_TRUE(str_eq(xvg.header_info.legends[1], STR_LIT("LJ-SR:2S29-2S29")));
    EXPECT_TRUE(str_eq(xvg.header_info.legends[2], STR_LIT("Coul-SR:2S29-SOL")));
    EXPECT_TRUE(str_eq(xvg.header_info.legends[3], STR_LIT("LJ-SR:2S29-SOL")));

    EXPECT_TRUE(str_eq(xvg.header_info.title, STR_LIT("GROMACS Energies")));
    EXPECT_TRUE(str_eq(xvg.header_info.xaxis_label, STR_LIT("Time (ps)")));
    EXPECT_TRUE(str_eq(xvg.header_info.yaxis_label, STR_LIT("(kJ/mol)")));

    md_xvg_free(&xvg, md_get_heap_allocator());
}

UTEST(xvg, lj_sr_lig_protein) {
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/LJ-SR_LIG-Protein.xvg");
    md_xvg_t xvg = {0};
    bool result = md_xvg_parse_file(&xvg, path, md_get_heap_allocator());
    ASSERT_TRUE(result);

    EXPECT_EQ(2,   xvg.num_fields);
    EXPECT_EQ(1001, xvg.num_values);

    // Test first data point: 0.000000  -129.218613
    EXPECT_NEAR(0.000000, xvg.fields[0][0], 1.0e-6f);
    EXPECT_NEAR(-129.218613, xvg.fields[1][0], 1.0e-6f);

    // Test last data point: 100000.000000  -143.968918
    EXPECT_NEAR(100000.000000, xvg.fields[0][1000], 1.0e-6f);
    EXPECT_NEAR(-143.968918, xvg.fields[1][1000], 1.0e-6f);

    // Test some middle data point
    EXPECT_NEAR(500.000000, xvg.fields[0][5], 1.0e-6f);
    EXPECT_NEAR(-159.923813, xvg.fields[1][5], 1.0e-6f);

    ASSERT_EQ(1, xvg.header_info.num_legends);
    EXPECT_TRUE(str_eq(xvg.header_info.legends[0], STR_LIT("LJ-SR:LIG-Protein")));

    EXPECT_TRUE(str_eq(xvg.header_info.title, STR_LIT("GROMACS Energies")));
    EXPECT_TRUE(str_eq(xvg.header_info.xaxis_label, STR_LIT("Time (ps)")));
    EXPECT_TRUE(str_eq(xvg.header_info.yaxis_label, STR_LIT("(kJ/mol)")));

    md_xvg_free(&xvg, md_get_heap_allocator());
}

// ### SERIES ALONG A RUN ###

// A run of the given times, and nothing else: all a series needs to be published along.
static void series_run(md_system_t* sys, str_t run, const double* times, size_t n, md_unit_t unit) {
    char buf[256];
    md_attributes_create(&sys->attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, STR_LIT("time")),
        .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {(uint32_t)n}},
        .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = unit, .data = times, .byte_size = n * sizeof(double)});
}

// energy.xvg: time in ps every 2 ps, four energy columns. Along a run with a frame every 4 ps each
// frame finds its row by time; along a run offset from it, none does and nothing is published.
UTEST(xvg, series_along_a_run) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    sys.attributes.alloc = arena;

    enum { F = 100 };
    double every4[F], offset[F];
    for (int i = 0; i < F; ++i) { every4[i] = 4.0 * i; offset[i] = 4.0 * i + 1.0; }
    series_run(&sys, STR_LIT("run/r"), every4, F, md_unit_picosecond());
    series_run(&sys, STR_LIT("run/o"), offset, F, md_unit_picosecond());

    const str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/energy.xvg");
    ASSERT_TRUE(md_xvg_system_supplement_from_file(&sys, path, STR_LIT("run/r")));
    EXPECT_FALSE(md_xvg_system_supplement_from_file(&sys, path, STR_LIT("run/o")));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/o/xvg/energy/time")) == NULL);

    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/r/xvg/energy/time"));
    const md_attribute_t* coul = md_attributes_find(&sys.attributes, STR_LIT("run/r/xvg/energy/coul_sr_2s29_2s29"));
    const md_attribute_t* lj   = md_attributes_find(&sys.attributes, STR_LIT("run/r/xvg/energy/lj_sr_2s29_sol"));
    const md_attribute_t* src  = md_attributes_find(&sys.attributes, STR_LIT("run/r/xvg/energy/source"));
    ASSERT_TRUE(time && coul && lj && src);
    EXPECT_EQ(356u, time->format.shape[0]);
    EXPECT_EQ(time, md_attributes_axis(&sys.attributes, coul));
    EXPECT_TRUE(md_unit_equal(time->unit, md_unit_picosecond()));
    EXPECT_TRUE(str_eq(coul->label, STR_LIT("Coul-SR:2S29-2S29")));
    md_unit_t kj_per_mol;
    ASSERT_TRUE(md_unit_parse(&kj_per_mol, STR_LIT("kJ/mol")));
    EXPECT_TRUE(md_unit_equal(coul->unit, kj_per_mol));

    // Frame 1 of the run is at 4 ps, row 2 of the file.
    const str_t paths[] = { STR_INIT("xvg/energy/coul_sr_2s29_2s29") };
    md_system_extract_t* ex = md_system_extract_begin(&sys, STR_LIT("run/r"), paths, 1, md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);
    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, 0);
    ASSERT_TRUE(md_system_extract_frame(ex, 1, &st));
    const md_attribute_t* v = md_attributes_find(&st.attributes, STR_LIT("xvg/energy/coul_sr_2s29_2s29"));
    ASSERT_TRUE(v != NULL);
    EXPECT_NEAR(-443.881683f, ((const float*)v->data)[0], 1.0e-3f);
    md_system_extract_end(ex);

    // And through the script, which is what the file is loaded for
    md_script_ir_t* ir = md_script_ir_create(arena);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("e = attr(\"run/r/xvg/energy/lj_sr_2s29_sol\");"), &sys, NULL));
    EXPECT_TRUE(md_script_ir_valid(ir));
    md_script_ir_free(ir);

    // Loading it again replaces it
    ASSERT_TRUE(md_xvg_system_supplement_from_file(&sys, path, STR_LIT("run/r")));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/r/xvg/energy/coul_sr_2s29_2s29")) != NULL);

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

// A CSV with a time column goes on its own axis; one without has a row per frame, or is refused.
UTEST(xvg, csv_series_along_a_run) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    sys.attributes.alloc = arena;
    const double times[4] = { 0.0, 1.0, 2.0, 3.0 };
    series_run(&sys, STR_LIT("run/r"), times, 4, md_unit_picosecond());

    const str_t with_time = STR_INIT("md_unittest_series_time.csv");
    const str_t per_frame = STR_INIT("md_unittest_series_frames.csv");
    const str_t too_short = STR_INIT("md_unittest_series_short.csv");
    md_file_t f = {0};
    ASSERT_TRUE(md_file_open(&f, with_time, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    md_file_printf(f, "Time (ps),Distance (nm),Angle\n");
    for (int i = 0; i <= 6; ++i) md_file_printf(f, "%g,%g,%g\n", 0.5 * i, 1.0 + i, 10.0 * i);
    md_file_close(&f);
    ASSERT_TRUE(md_file_open(&f, per_frame, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    md_file_printf(f, "Count\n1\n2\n3\n4\n");
    md_file_close(&f);
    ASSERT_TRUE(md_file_open(&f, too_short, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    md_file_printf(f, "Count\n1\n2\n3\n");
    md_file_close(&f);

    ASSERT_TRUE(md_csv_system_supplement_from_file(&sys, with_time, STR_LIT("run/r")));
    const md_attribute_t* dist = md_attributes_find(&sys.attributes, STR_LIT("run/r/csv/md_unittest_series_time/distance_nm"));
    ASSERT_TRUE(dist != NULL);
    EXPECT_TRUE(md_unit_equal(dist->unit, md_unit_nanometer()));
    // Frame 3 is at 3 ps, row 6
    EXPECT_EQ(7.0f, ((const float*)dist->data)[6]);

    ASSERT_TRUE(md_csv_system_supplement_from_file(&sys, per_frame, STR_LIT("run/r")));
    const md_attribute_t* count = md_attributes_find(&sys.attributes, STR_LIT("run/r/csv/md_unittest_series_frames/count"));
    ASSERT_TRUE(count != NULL);
    EXPECT_EQ(md_attributes_find(&sys.attributes, STR_LIT("run/r/time")), md_attributes_axis(&sys.attributes, count));

    EXPECT_FALSE(md_csv_system_supplement_from_file(&sys, too_short, STR_LIT("run/r")));
    EXPECT_FALSE(md_csv_system_supplement_from_file(&sys, per_frame, STR_LIT("run/none")));

    remove(with_time.ptr);
    remove(per_frame.ptr);
    remove(too_short.ptr);
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}
