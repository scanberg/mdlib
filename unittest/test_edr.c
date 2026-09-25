#include "utest.h"
#include <math.h>

#include <md_edr.h>
#include <core/md_allocator.h>

#define MAX_VALIDATION_SAMPLES 100

UTEST(edr, pullout) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR_LIT(MD_UNITTEST_DATA_DIR "/inside-md-pullout.edr"), md_get_heap_allocator());
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}

UTEST(edr, orires) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR_LIT(MD_UNITTEST_DATA_DIR "/orires.edr"), md_get_heap_allocator());
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}

UTEST(edr, ener) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR_LIT(MD_UNITTEST_DATA_DIR "/ener.edr"), md_get_heap_allocator());
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}

UTEST(edr, dhdl) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR_LIT(MD_UNITTEST_DATA_DIR "/dhdl.edr"), md_get_heap_allocator());
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}

UTEST(edr, comprehensive_validation) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Test all available EDR files systematically
    str_t edr_files[] = {
        STR_INIT(MD_UNITTEST_DATA_DIR "/ener.edr"),
        STR_INIT(MD_UNITTEST_DATA_DIR "/dhdl.edr"),
        STR_INIT(MD_UNITTEST_DATA_DIR "/inside-md-pullout.edr"),
        STR_INIT(MD_UNITTEST_DATA_DIR "/orires.edr")
    };
    
    for (int f = 0; f < 4; ++f) {
        md_edr_energies_t energies = {0};
        bool result = md_edr_energies_parse_file(&energies, edr_files[f], alloc);
        ASSERT_TRUE(result);
        
        EXPECT_GT(energies.num_energies, 0);
        EXPECT_GT(energies.num_frames, 0);
        
        // Validate that frame times are reasonable
        for (int64_t i = 0; i < MIN(MAX_VALIDATION_SAMPLES, energies.num_frames); ++i) {
            EXPECT_FALSE(isnan(energies.frame_time[i]));
            EXPECT_FALSE(isinf(energies.frame_time[i]));
            EXPECT_GE(energies.frame_time[i], 0.0);
        }
        
        // Validate that energy data is accessible
        for (int64_t j = 0; j < energies.num_energies; ++j) {
            if (energies.energy[j].values) {
                for (int64_t i = 0; i < MIN(10, energies.num_frames); ++i) {
                    EXPECT_FALSE(isnan(energies.energy[j].values[i]));
                    EXPECT_FALSE(isinf(energies.energy[j].values[i]));
                }
            }
        }
        
        md_edr_energies_free(&energies);
    }
}

UTEST(edr, nonexistent_file) {
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/nonexistent.edr");
    
    md_edr_energies_t energies = {0};
    bool result = md_edr_energies_parse_file(&energies, path, alloc);
    EXPECT_FALSE(result);
    
    md_edr_energies_free(&energies);
}

// ### ATTRIBUTES ###

#include <md_system.h>

static const md_attribute_t* edr_attr(const md_system_t* sys, const char* path) {
    return md_attributes_find(&sys->attributes, str_from_cstr(path));
}

static const md_edr_energy_t* edr_term(const md_edr_energies_t* e, const char* name) {
    for (size_t i = 0; i < e->num_energies; ++i) {
        if (str_eq_cstr(e->energy[i].name, name)) return &e->energy[i];
    }
    return NULL;
}

// Every term lands under <run>/edr along the file's own axis, tensors and vectors as one value, and
// the numbers are the file's own - in double, as the file may hold them.
UTEST(edr, supplement_publishes_terms_along_their_own_axis) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_edr_energies_t e = {0};
    ASSERT_TRUE(md_edr_energies_parse_file(&e, STR_LIT(MD_UNITTEST_DATA_DIR "/ener.edr"), alloc));

    md_system_t sys = {.alloc = alloc};
    ASSERT_TRUE(md_edr_system_supplement(&sys, &e, STR_LIT("run/test")));

    const md_attribute_t* axis = edr_attr(&sys, "run/test/edr/time");
    ASSERT_TRUE(axis != NULL);
    EXPECT_EQ((uint32_t)e.num_frames, axis->format.shape[0]);
    EXPECT_TRUE(md_unit_equal(axis->unit, md_unit_picosecond()));

    // 40 terms: Vir-* and Pres-* fold into two tensors, Box-X/Y/Z and Box-Vel-XX/YY/ZZ into two
    // vectors, which leaves 20 attributes beside the axis.
    EXPECT_EQ(md_attributes_query(NULL, 0, &sys.attributes, STR_LIT("run/test/edr")), 21u);

    const md_attribute_t* pot = edr_attr(&sys, "run/test/edr/potential");
    ASSERT_TRUE(pot != NULL);
    EXPECT_EQ(axis, md_attributes_axis(&sys.attributes, pot));
    EXPECT_TRUE(str_eq_cstr(pot->label, "Potential"));
    EXPECT_TRUE(md_unit_equal(pot->unit, edr_term(&e, "Potential")->unit));
    double v[9];
    md_attribute_slice_t row = md_attribute_slice_1(7);
    ASSERT_EQ(md_attribute_extract_slice_f64(v, 9, pot, &row, md_unit_none()), 1u);
    EXPECT_EQ(v[0], edr_term(&e, "Potential")->values[7]);

    // Folded names: the label keeps what GROMACS called it.
    const md_attribute_t* lj = edr_attr(&sys, "run/test/edr/lj_sr");
    ASSERT_TRUE(lj != NULL);
    EXPECT_TRUE(str_eq_cstr(lj->label, "LJ (SR)"));
    EXPECT_TRUE(edr_attr(&sys, "run/test/edr/kinetic_en") != NULL);
    EXPECT_TRUE(edr_attr(&sys, "run/test/edr/surf_surften") != NULL);

    // The virial is one 3x3 value per frame, row major: [0][1] is Vir-XY.
    const md_attribute_t* vir = edr_attr(&sys, "run/test/edr/vir");
    ASSERT_TRUE(vir != NULL);
    EXPECT_EQ(3u, vir->format.rank);
    EXPECT_EQ(3u, vir->format.shape[1]);
    EXPECT_EQ(3u, vir->format.shape[2]);
    EXPECT_EQ(1u, vir->format.components);
    ASSERT_EQ(md_attribute_extract_slice_f64(v, 9, vir, &row, md_unit_none()), 9u);
    EXPECT_EQ(v[1], edr_term(&e, "Vir-XY")->values[7]);
    EXPECT_EQ(v[5], edr_term(&e, "Vir-YZ")->values[7]);
    EXPECT_TRUE(edr_attr(&sys, "run/test/edr/vir_xx") == NULL);

    // Box-X/Y/Z is one 3-vector, and so is the diagonal Box-Vel-XX/YY/ZZ.
    const md_attribute_t* box = edr_attr(&sys, "run/test/edr/box");
    ASSERT_TRUE(box != NULL);
    EXPECT_EQ(1u, box->format.rank);
    EXPECT_EQ(3u, box->format.components);
    ASSERT_EQ(md_attribute_extract_slice_f64(v, 9, box, &row, md_unit_none()), 3u);
    EXPECT_EQ(v[2], edr_term(&e, "Box-Z")->values[7]);
    const md_attribute_t* box_vel = edr_attr(&sys, "run/test/edr/box_vel");
    ASSERT_TRUE(box_vel != NULL);
    EXPECT_EQ(3u, box_vel->format.components);

    md_attributes_free(&sys.attributes);
    md_edr_energies_free(&e);
}

static md_attribute_id_t edr_publish_run_axis(md_system_t* sys, const double* times, uint32_t n) {
    return md_attributes_create(&sys->attributes, &(md_attribute_desc_t){
        .path = STR_INIT("run/test/time"),
        .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {n}},
        .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_picosecond(),
        .data = times, .byte_size = n * sizeof(double)});
}

// A run that already has frames only takes an energy file that covers every one of them, and a
// refusal leaves the table as it was.
UTEST(edr, supplement_requires_every_frame_of_the_run) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_edr_energies_t e = {0};
    ASSERT_TRUE(md_edr_energies_parse_file(&e, STR_LIT(MD_UNITTEST_DATA_DIR "/ener.edr"), alloc));

    // ener.edr is written every 0.2 ps from 0 to 10 ps.
    {
        md_system_t sys = {.alloc = alloc};
        sys.attributes.alloc = alloc;
        const double frames[3] = {0.0, 1.0, 2.0};
        ASSERT_NE(edr_publish_run_axis(&sys, frames, 3), MD_ATTRIBUTE_INVALID);
        EXPECT_TRUE(md_edr_system_supplement(&sys, &e, STR_LIT("run/test")));

        // Frame 2 of the run is 2 ps, which is row 10 of the file.
        size_t idx = 0;
        EXPECT_TRUE(md_attribute_axis_map(&idx, edr_attr(&sys, "run/test/time"), 2, edr_attr(&sys, "run/test/edr/time")));
        EXPECT_EQ(idx, 10u);
        md_attributes_free(&sys.attributes);
    }
    {
        md_system_t sys = {.alloc = alloc};
        sys.attributes.alloc = alloc;
        const double frames[2] = {0.0, 1.1};
        ASSERT_NE(edr_publish_run_axis(&sys, frames, 2), MD_ATTRIBUTE_INVALID);
        EXPECT_FALSE(md_edr_system_supplement(&sys, &e, STR_LIT("run/test")));
        EXPECT_EQ(md_attributes_query(NULL, 0, &sys.attributes, STR_LIT("run/test/edr")), 0u);
        md_attributes_free(&sys.attributes);
    }

    md_edr_energies_free(&e);
}

// A second energy file replaces the first rather than mixing its terms into an axis that only
// describes one of them.
UTEST(edr, supplement_replaces_the_previous_file) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_system_t sys = {.alloc = alloc};

    ASSERT_TRUE(md_edr_system_supplement_from_file(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/ener.edr"), STR_LIT("run/test")));
    EXPECT_TRUE(edr_attr(&sys, "run/test/edr/box") != NULL);

    const md_attribute_t* src = edr_attr(&sys, "run/test/edr/source");
    ASSERT_TRUE(src != NULL);
    EXPECT_TRUE(str_eq_cstr(md_attribute_str(&sys.attributes, src, 0), MD_UNITTEST_DATA_DIR "/ener.edr"));

    ASSERT_TRUE(md_edr_system_supplement_from_file(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/dhdl.edr"), STR_LIT("run/test")));
    EXPECT_TRUE(edr_attr(&sys, "run/test/edr/box") == NULL);
    EXPECT_TRUE(edr_attr(&sys, "run/test/edr/dvcoul_dl") != NULL);
    EXPECT_EQ(101u, edr_attr(&sys, "run/test/edr/time")->format.shape[0]);
    src = edr_attr(&sys, "run/test/edr/source");
    ASSERT_TRUE(src != NULL);
    EXPECT_TRUE(str_eq_cstr(md_attribute_str(&sys.attributes, src, 0), MD_UNITTEST_DATA_DIR "/dhdl.edr"));

    md_attributes_free(&sys.attributes);
}

// Times that go backwards cannot be searched, and a restarted run written with overlap is the usual
// way to get them.
UTEST(edr, supplement_refuses_decreasing_times) {
    md_allocator_i* alloc = md_get_heap_allocator();
    double times[3]  = {0.0, 2.0, 1.0};
    double values[3] = {1.0, 2.0, 3.0};
    md_edr_energy_t term = {.name = STR_INIT("Potential"), .unit_str = STR_INIT("kJ/mol"), .values = values};
    md_edr_energies_t e = {.num_frames = 3, .frame_time = times, .num_energies = 1, .energy = &term};

    md_system_t sys = {.alloc = alloc};
    EXPECT_FALSE(md_edr_system_supplement(&sys, &e, STR_LIT("run/test")));
    times[2] = 4.0;
    EXPECT_TRUE(md_edr_system_supplement(&sys, &e, STR_LIT("run/test")));
    EXPECT_TRUE(edr_attr(&sys, "run/test/edr/potential") != NULL);

    md_attributes_free(&sys.attributes);
}
