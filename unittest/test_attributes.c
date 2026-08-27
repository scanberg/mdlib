#include "utest.h"

#include <md_system.h>
#include <core/md_unit.h>
#include <md_pdb.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>

static md_attribute_format_t fmt1(md_attribute_type_t type, uint32_t n) {
    md_attribute_format_t f = {0};
    f.type = type;
    f.rank = 1;
    f.shape[0] = n;
    return f;
}

static md_attribute_format_t fmt2(md_attribute_type_t type, uint32_t n, uint32_t c) {
    md_attribute_format_t f = {0};
    f.type = type;
    f.rank = 2;
    f.shape[0] = n;
    f.shape[1] = c;
    return f;
}

UTEST(attributes, format_math) {
    md_attribute_format_t vec = fmt2(MD_ATTRIBUTE_TYPE_F32, 10, 3);
    EXPECT_EQ(md_attribute_components(&vec),    3u);
    EXPECT_EQ(md_attribute_value_count(&vec),   10u);
    EXPECT_EQ(md_attribute_element_count(&vec), 30u);
    EXPECT_EQ(md_attribute_byte_size(&vec),     120u);

    md_attribute_format_t scalars = fmt1(MD_ATTRIBUTE_TYPE_F64, 10);
    EXPECT_EQ(md_attribute_components(&scalars),    1u);
    EXPECT_EQ(md_attribute_value_count(&scalars),   10u);
    EXPECT_EQ(md_attribute_element_count(&scalars), 10u);
    EXPECT_EQ(md_attribute_byte_size(&scalars),     80u);

    // rank 0 is a single value, not an empty one
    md_attribute_format_t single = {.type = MD_ATTRIBUTE_TYPE_F64, .rank = 0};
    EXPECT_EQ(md_attribute_components(&single),    1u);
    EXPECT_EQ(md_attribute_value_count(&single),   1u);
    EXPECT_EQ(md_attribute_element_count(&single), 1u);
    EXPECT_EQ(md_attribute_byte_size(&single),     8u);

    // {M,N,3}: the trailing axis is components, everything before it indexes
    md_attribute_format_t modes = {.type = MD_ATTRIBUTE_TYPE_F32, .rank = 3, .shape = {42, 10, 3}};
    EXPECT_EQ(md_attribute_components(&modes),    3u);
    EXPECT_EQ(md_attribute_value_count(&modes),   420u);
    EXPECT_EQ(md_attribute_element_count(&modes), 1260u);

    EXPECT_EQ(md_attribute_type_size(MD_ATTRIBUTE_TYPE_U8),   1u);
    EXPECT_EQ(md_attribute_type_size(MD_ATTRIBUTE_TYPE_I64),  8u);
    EXPECT_EQ(md_attribute_type_size(MD_ATTRIBUTE_TYPE_NONE), 0u);
}

UTEST(attributes, create_and_reject) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_id_t id = md_attributes_create(&t, STR_LIT("atom/velocity"), fmt2(MD_ATTRIBUTE_TYPE_F32, 10, 3), md_unit_none(), NULL, 0);
    EXPECT_NE(id, MD_ATTRIBUTE_INVALID);
    EXPECT_EQ(md_attributes_count(&t), 1u);

    // create zeroes the storage
    const float* data = (const float*)md_attributes_get(&t, id)->data;
    for (int i = 0; i < 30; ++i) {
        EXPECT_EQ(data[i], 0.0f);
    }

    EXPECT_EQ(md_attributes_create(&t, STR_LIT("atom/velocity"), fmt1(MD_ATTRIBUTE_TYPE_F32, 10), md_unit_none(), NULL, 0), MD_ATTRIBUTE_INVALID); // duplicate
    EXPECT_EQ(md_attributes_create(&t, STR_LIT("/atom/x"), fmt1(MD_ATTRIBUTE_TYPE_F32, 10), md_unit_none(), NULL, 0), MD_ATTRIBUTE_INVALID); // leading separator
    EXPECT_EQ(md_attributes_create(&t, STR_LIT("atom/x/"), fmt1(MD_ATTRIBUTE_TYPE_F32, 10), md_unit_none(), NULL, 0), MD_ATTRIBUTE_INVALID); // trailing separator
    EXPECT_EQ(md_attributes_create(&t, STR_LIT("atom//x"), fmt1(MD_ATTRIBUTE_TYPE_F32, 10), md_unit_none(), NULL, 0), MD_ATTRIBUTE_INVALID); // empty segment
    EXPECT_EQ(md_attributes_create(&t, STR_LIT(""), fmt1(MD_ATTRIBUTE_TYPE_F32, 10), md_unit_none(), NULL, 0), MD_ATTRIBUTE_INVALID); // empty path
    EXPECT_EQ(md_attributes_create(&t, STR_LIT("a/b"), fmt1(MD_ATTRIBUTE_TYPE_NONE, 10), md_unit_none(), NULL, 0), MD_ATTRIBUTE_INVALID); // no type
    EXPECT_EQ(md_attributes_create(&t, STR_LIT("a/b"), fmt1(MD_ATTRIBUTE_TYPE_F32, 0), md_unit_none(), NULL, 0),  MD_ATTRIBUTE_INVALID); // zero extent
    EXPECT_EQ(md_attributes_count(&t), 1u);

    // an unset allocator is a rejection, not a crash
    md_attributes_t no_alloc = {0};
    EXPECT_EQ(md_attributes_create(&no_alloc, STR_LIT("a/b"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 0), MD_ATTRIBUTE_INVALID);

    md_attributes_free(&t);
    EXPECT_EQ(md_attributes_count(&t), 0u);
}

UTEST(attributes, create_with_data) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const float charges[4] = {-0.5f, 0.25f, 0.0f, 1.5f};
    md_attribute_id_t id = md_attributes_create(&t, STR_LIT("atom/charge"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), charges, sizeof(charges));
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    const md_attribute_t* a = md_attributes_get(&t, id);
    ASSERT_TRUE(a != NULL);
    for (int i = 0; i < 4; ++i) {
        EXPECT_EQ(((const float*)a->data)[i], charges[i]);
    }

    // the copy is the table's own, not a view onto the caller's buffer
    EXPECT_TRUE(a->data != (const void*)charges);

    // NULL reserves and zeroes
    md_attribute_id_t vel = md_attributes_create(&t, STR_LIT("atom/velocity"), fmt2(MD_ATTRIBUTE_TYPE_F32, 4, 3), md_unit_none(), NULL, 0);
    ASSERT_NE(vel, MD_ATTRIBUTE_INVALID);
    const float* v = (const float*)md_attributes_get(&t, vel)->data;
    for (int i = 0; i < 12; ++i) {
        EXPECT_EQ(v[i], 0.0f);
    }

    md_attributes_free(&t);
}

UTEST(attributes, create_rejects_a_size_that_disagrees_with_the_format) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const double wrong_type[4] = {1.0, 2.0, 3.0, 4.0};
    const float  values[4]     = {1.0f, 2.0f, 3.0f, 4.0f};

    // the f64-buffer-into-an-f32-format mistake, which is the whole point of taking a size
    EXPECT_EQ(md_attributes_create(&t, STR_LIT("a/b"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), wrong_type, sizeof(wrong_type)), MD_ATTRIBUTE_INVALID);
    // too few elements
    EXPECT_EQ(md_attributes_create(&t, STR_LIT("a/b"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), values, 2 * sizeof(float)), MD_ATTRIBUTE_INVALID);
    // a size with no pointer is a caller who meant to pass one
    EXPECT_EQ(md_attributes_create(&t, STR_LIT("a/b"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 4 * sizeof(float)), MD_ATTRIBUTE_INVALID);
    // a pointer with no size likewise
    EXPECT_EQ(md_attributes_create(&t, STR_LIT("a/b"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), values, 0), MD_ATTRIBUTE_INVALID);

    // nothing was registered by any of them
    EXPECT_EQ(md_attributes_count(&t), 0u);

    // and the correct call still works
    EXPECT_NE(md_attributes_create(&t, STR_LIT("a/b"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), values, sizeof(values)), MD_ATTRIBUTE_INVALID);
    EXPECT_EQ(md_attributes_count(&t), 1u);

    md_attributes_free(&t);
}

UTEST(attributes, path_and_lookup) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_id_t mul = md_attributes_create(&t, STR_LIT("atom/charge/mulliken"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 0);
    md_attribute_id_t scf = md_attributes_create(&t, STR_LIT("scf_energy"), fmt1(MD_ATTRIBUTE_TYPE_F64, 1), md_unit_none(), NULL, 0);

    const md_attribute_t* a = md_attributes_find(&t, STR_LIT("atom/charge/mulliken"));
    ASSERT_TRUE(a != NULL);
    EXPECT_EQ(a->id, mul);
    EXPECT_TRUE(str_eq(md_attribute_category(a), STR_LIT("atom/charge")));
    EXPECT_TRUE(str_eq(md_attribute_field(a),    STR_LIT("mulliken")));

    // a path with no separator is all field, no category
    const md_attribute_t* e = md_attributes_find(&t, STR_LIT("scf_energy"));
    ASSERT_TRUE(e != NULL);
    EXPECT_EQ(e->id, scf);
    EXPECT_TRUE(str_empty(md_attribute_category(e)));
    EXPECT_TRUE(str_eq(md_attribute_field(e), STR_LIT("scf_energy")));

    // an interior path is a prefix, not an attribute
    EXPECT_TRUE(md_attributes_find(&t, STR_LIT("atom/charge")) == NULL);
    EXPECT_TRUE(md_attributes_get(&t, MD_ATTRIBUTE_INVALID) == NULL);
    EXPECT_TRUE(md_attributes_get(&t, 0xDEADBEEF) == NULL);

    md_attributes_free(&t);
}

UTEST(attributes, unit_round_trips) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const float b[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    md_unit_t angstrom_sq = md_unit_pow(md_unit_angstrom(), 2);

    md_attribute_id_t bf  = md_attributes_create(&t, STR_LIT("atom/b_factor"),  fmt1(MD_ATTRIBUTE_TYPE_F32, 4), angstrom_sq,     b, sizeof(b));
    md_attribute_id_t occ = md_attributes_create(&t, STR_LIT("atom/occupancy"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(),  b, sizeof(b));

    EXPECT_TRUE(md_unit_equal(md_attributes_get(&t, bf)->unit, angstrom_sq));
    EXPECT_TRUE(md_unit_is_none(md_attributes_get(&t, occ)->unit));

    // a zeroed unit is dimensionless, not an invalid state to guard against
    md_attribute_id_t z = md_attributes_create(&t, STR_LIT("atom/anything"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), (md_unit_t){0}, NULL, 0);
    EXPECT_TRUE(md_unit_is_none(md_attributes_get(&t, z)->unit));

    md_attributes_free(&t);
}

// A dipole has no index space to anchor it, so its origin is a sibling of the same shape.
// The two disagree on units, which is why the unit sits on the attribute and not the group.
UTEST(attributes, anchored_vector_group) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const double vector[3] = {0.0, 0.0, 1.85};
    const double origin[3] = {1.0, 2.0, 3.0};
    md_attribute_format_t vec3 = fmt2(MD_ATTRIBUTE_TYPE_F64, 1, 3);

    md_unit_t au = md_unit_scl(md_unit_mul(md_unit_mul(md_unit_ampere(), md_unit_second()), md_unit_bohr_radius()), 1.602176634e-19);

    md_attributes_create(&t, STR_LIT("dipole/ground_state/vector"), vec3, au,                vector, sizeof(vector));
    md_attributes_create(&t, STR_LIT("dipole/ground_state/origin"), vec3, md_unit_angstrom(), origin, sizeof(origin));

    // one 3-vector, not three scalars: this is what rank 2 buys over rank 1 {3}
    const md_attribute_t* v = md_attributes_find(&t, STR_LIT("dipole/ground_state/vector"));
    const md_attribute_t* o = md_attributes_find(&t, STR_LIT("dipole/ground_state/origin"));
    ASSERT_TRUE(v != NULL);
    ASSERT_TRUE(o != NULL);
    EXPECT_EQ(md_attribute_components(&v->format),  3u);
    EXPECT_EQ(md_attribute_value_count(&v->format), 1u);
    EXPECT_EQ(md_attribute_byte_size(&v->format),   24u);

    // origin mirrors the vector's shape, so a consumer never has to branch on it
    EXPECT_EQ(md_attribute_value_count(&o->format), md_attribute_value_count(&v->format));
    EXPECT_EQ(md_attribute_components(&o->format),  md_attribute_components(&v->format));

    // and they legitimately differ in unit
    EXPECT_FALSE(md_unit_base_equal(v->unit, o->unit));
    EXPECT_TRUE(md_unit_equal(o->unit, md_unit_angstrom()));

    EXPECT_EQ(((const double*)v->data)[2], 1.85);
    EXPECT_EQ(((const double*)o->data)[0], 1.0);

    // the group is discoverable without knowing either name
    str_t children[4];
    EXPECT_EQ(md_attributes_query_children(children, ARRAY_SIZE(children), &t, STR_LIT("dipole")), 1u);
    EXPECT_TRUE(str_eq(children[0], STR_LIT("ground_state")));
    EXPECT_EQ(md_attributes_query_children(children, ARRAY_SIZE(children), &t, STR_LIT("dipole/ground_state")), 2u);
    EXPECT_TRUE(str_eq(children[0], STR_LIT("origin")));
    EXPECT_TRUE(str_eq(children[1], STR_LIT("vector")));

    md_attribute_id_t ids[4];
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("dipole")), 2u);

    md_attributes_free(&t);
}

UTEST(attributes, extract_converts_type) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const double src[4] = {1.5, -2.25, 0.0, 1.0e-8};
    md_attribute_id_t f64 = md_attributes_create(&t, STR_LIT("a/f64"), fmt1(MD_ATTRIBUTE_TYPE_F64, 4), md_unit_none(), src, sizeof(src));

    const int16_t ints[4] = {-3, 0, 7, 32767};
    md_attribute_id_t i16 = md_attributes_create(&t, STR_LIT("a/i16"), fmt1(MD_ATTRIBUTE_TYPE_I16, 4), md_unit_none(), ints, sizeof(ints));

    float dst[8];

    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, f64), md_unit_none()), 4u);
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(dst[i], (float)src[i], 1.0e-12f);
    }

    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, i16), md_unit_none()), 4u);
    for (int i = 0; i < 4; ++i) {
        EXPECT_EQ(dst[i], (float)ints[i]);
    }

    // the whole attribute, components included: {N,3} yields N*3 floats
    const float vecs[6] = {1,2,3, 4,5,6};
    md_attribute_id_t v = md_attributes_create(&t, STR_LIT("a/vec"), fmt2(MD_ATTRIBUTE_TYPE_F32, 2, 3), md_unit_none(), vecs, sizeof(vecs));
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, v), md_unit_none()), 6u);
    EXPECT_EQ(dst[5], 6.0f);

    // a destination which cannot hold it writes nothing
    EXPECT_EQ(md_attribute_extract_f32(dst, 3, md_attributes_get(&t, v), md_unit_none()), 0u);

    md_attributes_free(&t);
}

UTEST(attributes, extract_converts_unit) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const float lengths[2] = {10.0f, 2.5f};
    md_attribute_id_t len = md_attributes_create(&t, STR_LIT("a/length"), fmt1(MD_ATTRIBUTE_TYPE_F32, 2), md_unit_angstrom(), lengths, sizeof(lengths));

    float dst[4];

    // Angstrom into nanometre is a factor of ten
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, len), md_unit_nanometer()), 2u);
    EXPECT_NEAR(dst[0], 1.0f,  1.0e-6f);
    EXPECT_NEAR(dst[1], 0.25f, 1.0e-6f);

    // same unit is the identity
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, len), md_unit_angstrom()), 2u);
    EXPECT_NEAR(dst[0], 10.0f, 1.0e-6f);

    // none as the target means as stored, whatever the attribute carries
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, len), md_unit_none()), 2u);
    EXPECT_NEAR(dst[0], 10.0f, 1.0e-6f);

    md_attributes_free(&t);
}

// The refusal is the reason unit conversion lives in the extractor at all: a dipole rendered at
// 2.5x the right length looks entirely plausible.
UTEST(attributes, extract_refuses_incompatible_units) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const double dipole[3] = {0.0, 0.0, 0.7273};   // e a0
    md_attribute_id_t d = md_attributes_create(&t, STR_LIT("dipole/ground_state/vector"), fmt2(MD_ATTRIBUTE_TYPE_F64, 1, 3), md_unit_elementary_charge_bohr(), dipole, sizeof(dipole));

    const float occ[2] = {0.5f, 1.0f};
    md_attribute_id_t o = md_attributes_create(&t, STR_LIT("atom/occupancy"), fmt1(MD_ATTRIBUTE_TYPE_F32, 2), md_unit_none(), occ, sizeof(occ));

    float dst[4] = {-1.0f, -1.0f, -1.0f, -1.0f};

    // a dipole is not a length
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, d), md_unit_angstrom()), 0u);
    EXPECT_EQ(dst[0], -1.0f);   // nothing written

    // and a dimensionless quantity cannot be given one
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, o), md_unit_angstrom()), 0u);

    // but the dipole does convert to Debye
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, d), md_unit_debye()), 3u);
    EXPECT_NEAR(dst[2], 0.7273f * 2.541746473f, 1.0e-4f);

    md_attributes_free(&t);
}

// The b factor is the round trip a colormap will actually do.
UTEST(attributes, extract_from_pdb) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"), MD_PDB_OPTION_NONE));

    const md_attribute_t* bfc = md_attributes_find(&sys.attributes, STR_LIT("atom/b_factor"));
    ASSERT_TRUE(bfc != NULL);

    float* dst = (float*)md_alloc(alloc, sys.atom.count * sizeof(float));
    md_unit_t angstrom_sq = md_unit_pow(md_unit_angstrom(), 2);

    EXPECT_EQ(md_attribute_extract_f32(dst, sys.atom.count, bfc, angstrom_sq), sys.atom.count);
    EXPECT_NEAR(dst[0],   201.24f, 1.0e-3f);
    EXPECT_NEAR(dst[980], 150.42f, 1.0e-3f);

    // into nm^2, a factor of 1/100
    EXPECT_EQ(md_attribute_extract_f32(dst, sys.atom.count, bfc, md_unit_pow(md_unit_nanometer(), 2)), sys.atom.count);
    EXPECT_NEAR(dst[0], 2.0124f, 1.0e-4f);

    // and the occupancy alongside it is dimensionless, so Angstrom squared is a refusal
    const md_attribute_t* occ = md_attributes_find(&sys.attributes, STR_LIT("atom/occupancy"));
    ASSERT_TRUE(occ != NULL);
    EXPECT_EQ(md_attribute_extract_f32(dst, sys.atom.count, occ, angstrom_sq), 0u);
    EXPECT_EQ(md_attribute_extract_f32(dst, sys.atom.count, occ, md_unit_none()), sys.atom.count);
    EXPECT_NEAR(dst[980], 0.50f, 1.0e-5f);

    md_system_free(&sys);
    md_arena_allocator_destroy(alloc);
}

UTEST(attributes, query_matches_at_segment_boundary) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_id_t vel  = md_attributes_create(&t, STR_LIT("atom/velocity"), fmt2(MD_ATTRIBUTE_TYPE_F32, 4, 3), md_unit_none(), NULL, 0);
    md_attribute_id_t mul  = md_attributes_create(&t, STR_LIT("atom/charge/mulliken"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 0);
    md_attribute_id_t low  = md_attributes_create(&t, STR_LIT("atom/charge/lowdin"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 0);
    md_attribute_id_t zed  = md_attributes_create(&t, STR_LIT("atomic_number/z"), fmt1(MD_ATTRIBUTE_TYPE_U8,  4), md_unit_none(), NULL, 0);
    md_attribute_id_t dip  = md_attributes_create(&t, STR_LIT("dipole/magnetic"), fmt1(MD_ATTRIBUTE_TYPE_F64, 3), md_unit_none(), NULL, 0);
    // sorts BETWEEN "atom" and "atom/..." because '-' is below '/', so a naive walk stops early
    md_attribute_id_t dash = md_attributes_create(&t, STR_LIT("atom-x/weird"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 0);
    EXPECT_EQ(md_attributes_count(&t), 6u);

    md_attribute_id_t ids[16];

    size_t n = md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("atom"));
    EXPECT_EQ(n, 3u);
    for (size_t i = 0; i < n; ++i) {
        EXPECT_NE(ids[i], zed);
        EXPECT_NE(ids[i], dash);
        EXPECT_NE(ids[i], dip);
    }

    // a trailing separator on the prefix is ignored
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("atom/")),   3u);
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("atom///")), 3u);

    // results come back in name order
    EXPECT_EQ(md_attributes_get(&t, ids[0])->id, low);
    EXPECT_EQ(md_attributes_get(&t, ids[1])->id, mul);
    EXPECT_EQ(md_attributes_get(&t, ids[2])->id, vel);

    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("atom/charge")),   2u);
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("atomic_number")), 1u);
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("atom-x")),        1u);
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("")),              6u);
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("nope")),          0u);
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("ato")),           0u); // partial segment

    // the return is the total, the cap only bounds the write
    EXPECT_EQ(md_attributes_query(NULL, 0, &t, STR_LIT("atom")), 3u);
    ids[1] = MD_ATTRIBUTE_INVALID;
    EXPECT_EQ(md_attributes_query(ids, 1, &t, STR_LIT("atom")), 3u);
    EXPECT_EQ(ids[1], MD_ATTRIBUTE_INVALID);

    md_attributes_free(&t);
}

UTEST(attributes, query_children) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attributes_create(&t, STR_LIT("atom/velocity"), fmt2(MD_ATTRIBUTE_TYPE_F32, 4, 3), md_unit_none(), NULL, 0);
    md_attributes_create(&t, STR_LIT("atom/charge/mulliken"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 0);
    md_attributes_create(&t, STR_LIT("atom/charge/lowdin"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 0);
    md_attributes_create(&t, STR_LIT("dipole/magnetic"), fmt1(MD_ATTRIBUTE_TYPE_F64, 3), md_unit_none(), NULL, 0);

    str_t children[8];

    size_t n = md_attributes_query_children(children, ARRAY_SIZE(children), &t, STR_LIT("atom"));
    EXPECT_EQ(n, 2u);
    EXPECT_TRUE(str_eq(children[0], STR_LIT("charge")));
    EXPECT_TRUE(str_eq(children[1], STR_LIT("velocity")));

    n = md_attributes_query_children(children, ARRAY_SIZE(children), &t, STR_LIT("atom/charge"));
    EXPECT_EQ(n, 2u);
    EXPECT_TRUE(str_eq(children[0], STR_LIT("lowdin")));
    EXPECT_TRUE(str_eq(children[1], STR_LIT("mulliken")));

    n = md_attributes_query_children(children, ARRAY_SIZE(children), &t, STR_LIT(""));
    EXPECT_EQ(n, 2u);
    EXPECT_TRUE(str_eq(children[0], STR_LIT("atom")));
    EXPECT_TRUE(str_eq(children[1], STR_LIT("dipole")));

    // a leaf has no children
    EXPECT_EQ(md_attributes_query_children(children, ARRAY_SIZE(children), &t, STR_LIT("atom/velocity")), 0u);

    md_attributes_free(&t);
}

UTEST(attributes, data_is_type_checked) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_id_t id = md_attributes_create(&t, STR_LIT("atom/charge"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 0);

    float* p = (float*)md_attributes_data(&t, id, MD_ATTRIBUTE_TYPE_F32);
    ASSERT_TRUE(p != NULL);
    p[0] = 1.0f;
    p[3] = -2.5f;

    const md_attribute_t* a = md_attributes_get(&t, id);
    EXPECT_EQ(((const float*)a->data)[0],  1.0f);
    EXPECT_EQ(((const float*)a->data)[3], -2.5f);

    EXPECT_TRUE(md_attributes_data(&t, id, MD_ATTRIBUTE_TYPE_F64) == NULL);
    EXPECT_TRUE(md_attributes_data(&t, id, MD_ATTRIBUTE_TYPE_I32) == NULL);
    EXPECT_TRUE(md_attributes_data(&t, 0xBADC0DE, MD_ATTRIBUTE_TYPE_F32) == NULL);

    md_attributes_free(&t);
}

UTEST(attributes, remove_keeps_the_table_sorted) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_id_t vel = md_attributes_create(&t, STR_LIT("atom/velocity"), fmt2(MD_ATTRIBUTE_TYPE_F32, 4, 3), md_unit_none(), NULL, 0);
    md_attribute_id_t mul = md_attributes_create(&t, STR_LIT("atom/charge/mulliken"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 0);
    md_attribute_id_t low = md_attributes_create(&t, STR_LIT("atom/charge/lowdin"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 0);

    EXPECT_TRUE(md_attributes_remove(&t, mul));
    EXPECT_FALSE(md_attributes_remove(&t, mul));
    EXPECT_EQ(md_attributes_count(&t), 2u);
    EXPECT_TRUE(md_attributes_find(&t, STR_LIT("atom/charge/mulliken")) == NULL);

    md_attribute_id_t ids[8];
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("atom")), 2u);
    EXPECT_EQ(md_attributes_get(&t, ids[0])->id, low);
    EXPECT_EQ(md_attributes_get(&t, ids[1])->id, vel);

    // the id is the path hash, so recreating the same path yields the same id
    md_attribute_id_t again = md_attributes_create(&t, STR_LIT("atom/charge/mulliken"), fmt1(MD_ATTRIBUTE_TYPE_F32, 4), md_unit_none(), NULL, 0);
    EXPECT_EQ(again, mul);
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("atom")), 3u);

    md_attributes_free(&t);
}

// The attributes a PDB actually carries. tubulin-A-B.pdb is the interesting one: it has both
// partial occupancies and altLoc entries, which the loader skips, so the columns must line up
// with the atoms that were kept rather than with the coordinate records.
UTEST(attributes, pdb_occupancy_and_b_factor) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"), MD_PDB_OPTION_NONE));

    // 6942 coordinate records, 15 of them the 'B' alternate of a pair
    EXPECT_EQ(sys.atom.count, 6927u);

    md_attribute_id_t ids[8];
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &sys.attributes, STR_LIT("atom")), 2u);

    const md_attribute_t* occ = md_attributes_find(&sys.attributes, STR_LIT("atom/occupancy"));
    const md_attribute_t* bfc = md_attributes_find(&sys.attributes, STR_LIT("atom/b_factor"));
    ASSERT_TRUE(occ != NULL);
    ASSERT_TRUE(bfc != NULL);

    EXPECT_EQ(occ->format.type, MD_ATTRIBUTE_TYPE_F32);
    EXPECT_EQ(occ->format.rank, 1u);
    EXPECT_EQ(occ->format.shape[0], (uint32_t)sys.atom.count);
    EXPECT_EQ(md_attribute_components(&occ->format), 1u);
    EXPECT_TRUE(str_eq(md_attribute_category(occ), STR_LIT("atom")));
    EXPECT_TRUE(str_eq(md_attribute_field(occ),    STR_LIT("occupancy")));

    const float* o = (const float*)occ->data;
    const float* b = (const float*)bfc->data;

    EXPECT_NEAR(o[0], 1.00f, 1.0e-5f);
    EXPECT_NEAR(b[0], 201.24f, 1.0e-3f);

    // Atom 980 is the 'A' alternate of GLN 133 CA. Atom 981 must be the next ACCEPTED atom,
    // the C of the same residue, not the discarded 'B' alternate whose b factor is 150.83.
    EXPECT_NEAR(o[980], 0.50f, 1.0e-5f);
    EXPECT_NEAR(b[980], 150.42f, 1.0e-3f);
    EXPECT_NEAR(o[981], 1.00f, 1.0e-5f);
    EXPECT_NEAR(b[981], 152.47f, 1.0e-3f);

    // the tail is still aligned after 15 skipped records
    EXPECT_NEAR(b[sys.atom.count - 1], 154.54f, 1.0e-3f);

    size_t partial = 0;
    for (size_t i = 0; i < sys.atom.count; ++i) {
        if (o[i] < 1.0f) partial += 1;
    }
    EXPECT_EQ(partial, 27u); // 4 at 0.25, 15 at 0.50, 4 at 0.60, 4 at 0.75

    md_system_free(&sys);
    md_arena_allocator_destroy(alloc);
}

// A uniform column carries nothing worth listing, so it is not published. 1k4r has real
// b factors but occupancy 1.00 throughout.
UTEST(attributes, pdb_uniform_column_is_not_published) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE));

    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("atom/occupancy")) == NULL);

    const md_attribute_t* bfc = md_attributes_find(&sys.attributes, STR_LIT("atom/b_factor"));
    ASSERT_TRUE(bfc != NULL);
    EXPECT_EQ(bfc->format.shape[0], (uint32_t)sys.atom.count);
    EXPECT_NEAR(((const float*)bfc->data)[0], 22.56f, 1.0e-3f);

    md_system_free(&sys);
    md_arena_allocator_destroy(alloc);
}

// c60.pdb has occupancy 1.00 and b factor 0.00 everywhere, so it publishes nothing at all.
UTEST(attributes, pdb_without_useful_columns_publishes_nothing) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/c60.pdb"), MD_PDB_OPTION_NONE));

    EXPECT_EQ(md_attributes_count(&sys.attributes), 0u);
    EXPECT_EQ(md_attributes_query(NULL, 0, &sys.attributes, STR_LIT("atom")), 0u);

    md_system_free(&sys);
    md_arena_allocator_destroy(alloc);
}
