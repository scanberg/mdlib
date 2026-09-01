#include "utest.h"

#include <md_system.h>
#include <core/md_unit.h>
#include <md_pdb.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>

// N scalars: one index axis of extent n, values one component wide.
static md_attribute_format_t fmt_scalars(md_attribute_type_t type, uint32_t n) {
    md_attribute_format_t f = {0};
    f.type = type;
    f.components = 1;
    f.rank = 1;
    f.shape[0] = n;
    return f;
}

// N values of c components each. Still ONE index axis: the component count is not an axis.
static md_attribute_format_t fmt_vec(md_attribute_type_t type, uint32_t n, uint32_t c) {
    md_attribute_format_t f = {0};
    f.type = type;
    f.components = c;
    f.rank = 1;
    f.shape[0] = n;
    return f;
}

UTEST(attributes, format_math) {
    md_attribute_format_t vec = fmt_vec(MD_ATTRIBUTE_TYPE_F32, 10, 3);
    EXPECT_EQ(md_attribute_components(&vec),    3u);
    EXPECT_EQ(md_attribute_value_count(&vec),   10u);
    EXPECT_EQ(md_attribute_element_count(&vec), 30u);
    EXPECT_EQ(md_attribute_byte_size(&vec),     120u);

    md_attribute_format_t scalars = fmt_scalars(MD_ATTRIBUTE_TYPE_F64, 10);
    EXPECT_EQ(md_attribute_components(&scalars),    1u);
    EXPECT_EQ(md_attribute_value_count(&scalars),   10u);
    EXPECT_EQ(md_attribute_element_count(&scalars), 10u);
    EXPECT_EQ(md_attribute_byte_size(&scalars),     80u);

    // rank 0 is a single value, not an empty one
    md_attribute_format_t single = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 0};
    EXPECT_EQ(md_attribute_components(&single),    1u);
    EXPECT_EQ(md_attribute_value_count(&single),   1u);
    EXPECT_EQ(md_attribute_element_count(&single), 1u);
    EXPECT_EQ(md_attribute_byte_size(&single),     8u);

    // {M,N} of 3 components: two index axes, and the 3 is not one of them
    md_attribute_format_t modes = {.type = MD_ATTRIBUTE_TYPE_F32, .components = 3, .rank = 2, .shape = {42, 10}};
    EXPECT_EQ(md_attribute_components(&modes),    3u);
    EXPECT_EQ(md_attribute_value_count(&modes),   420u);
    EXPECT_EQ(md_attribute_element_count(&modes), 1260u);

    EXPECT_EQ(md_attribute_type_size(MD_ATTRIBUTE_TYPE_U8),   1u);
    EXPECT_EQ(md_attribute_type_size(MD_ATTRIBUTE_TYPE_I64),  8u);
    EXPECT_EQ(md_attribute_type_size(MD_ATTRIBUTE_TYPE_NONE), 0u);
}

UTEST(attributes, create_and_reject) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/velocity"), .format = fmt_vec(MD_ATTRIBUTE_TYPE_F32, 10, 3), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    EXPECT_NE(id, MD_ATTRIBUTE_INVALID);
    EXPECT_EQ(md_attributes_count(&t), 1u);

    // create zeroes the storage
    const float* data = (const float*)md_attributes_get(&t, id)->data;
    for (int i = 0; i < 30; ++i) {
        EXPECT_EQ(data[i], 0.0f);
    }

    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/velocity"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 10), .unit = md_unit_none(), .data = NULL, .byte_size = 0}), MD_ATTRIBUTE_INVALID); // duplicate
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("/atom/x"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 10), .unit = md_unit_none(), .data = NULL, .byte_size = 0}), MD_ATTRIBUTE_INVALID); // leading separator
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/x/"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 10), .unit = md_unit_none(), .data = NULL, .byte_size = 0}), MD_ATTRIBUTE_INVALID); // trailing separator
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom//x"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 10), .unit = md_unit_none(), .data = NULL, .byte_size = 0}), MD_ATTRIBUTE_INVALID); // empty segment
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT(""), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 10), .unit = md_unit_none(), .data = NULL, .byte_size = 0}), MD_ATTRIBUTE_INVALID); // empty path
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/b"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_NONE, 10), .unit = md_unit_none(), .data = NULL, .byte_size = 0}), MD_ATTRIBUTE_INVALID); // no type
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/b"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 0), .unit = md_unit_none(), .data = NULL, .byte_size = 0}),  MD_ATTRIBUTE_INVALID); // zero extent
    EXPECT_EQ(md_attributes_count(&t), 1u);

    // an unset allocator is a rejection, not a crash
    md_attributes_t no_alloc = {0};
    EXPECT_EQ(md_attributes_create(&no_alloc, &(md_attribute_desc_t){.path = STR_LIT("a/b"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 0}), MD_ATTRIBUTE_INVALID);

    md_attributes_free(&t);
    EXPECT_EQ(md_attributes_count(&t), 0u);
}

UTEST(attributes, create_with_data) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const float charges[4] = {-0.5f, 0.25f, 0.0f, 1.5f};
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/charge"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = charges, .byte_size = sizeof(charges)});
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    const md_attribute_t* a = md_attributes_get(&t, id);
    ASSERT_TRUE(a != NULL);
    for (int i = 0; i < 4; ++i) {
        EXPECT_EQ(((const float*)a->data)[i], charges[i]);
    }

    // the copy is the table's own, not a view onto the caller's buffer
    EXPECT_TRUE(a->data != (const void*)charges);

    // NULL reserves and zeroes
    md_attribute_id_t vel = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/velocity"), .format = fmt_vec(MD_ATTRIBUTE_TYPE_F32, 4, 3), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
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
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/b"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = wrong_type, .byte_size = sizeof(wrong_type)}), MD_ATTRIBUTE_INVALID);
    // too few elements
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/b"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = values, .byte_size = 2 * sizeof(float)}), MD_ATTRIBUTE_INVALID);
    // a size with no pointer is a caller who meant to pass one
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/b"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 4 * sizeof(float)}), MD_ATTRIBUTE_INVALID);
    // a pointer with no size likewise
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/b"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = values, .byte_size = 0}), MD_ATTRIBUTE_INVALID);

    // nothing was registered by any of them
    EXPECT_EQ(md_attributes_count(&t), 0u);

    // and the correct call still works
    EXPECT_NE(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/b"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = values, .byte_size = sizeof(values)}), MD_ATTRIBUTE_INVALID);
    EXPECT_EQ(md_attributes_count(&t), 1u);

    md_attributes_free(&t);
}

UTEST(attributes, path_and_lookup) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_id_t mul = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/charge/mulliken"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    md_attribute_id_t scf = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("scf_energy"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F64, 1), .unit = md_unit_none(), .data = NULL, .byte_size = 0});

    const md_attribute_t* a = md_attributes_find(&t, STR_LIT("atom/charge/mulliken"));
    ASSERT_TRUE(a != NULL);
    EXPECT_EQ(a->id, mul);
    EXPECT_TRUE(str_eq(md_attribute_group(a), STR_LIT("atom/charge")));
    EXPECT_TRUE(str_eq(md_attribute_leaf(a),    STR_LIT("mulliken")));

    // a path with no separator is all field, no category
    const md_attribute_t* e = md_attributes_find(&t, STR_LIT("scf_energy"));
    ASSERT_TRUE(e != NULL);
    EXPECT_EQ(e->id, scf);
    EXPECT_TRUE(str_empty(md_attribute_group(e)));
    EXPECT_TRUE(str_eq(md_attribute_leaf(e), STR_LIT("scf_energy")));

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

    md_attribute_id_t bf  = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/b_factor"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = angstrom_sq, .data = b, .byte_size = sizeof(b)});
    md_attribute_id_t occ = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/occupancy"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = b, .byte_size = sizeof(b)});

    EXPECT_TRUE(md_unit_equal(md_attributes_get(&t, bf)->unit, angstrom_sq));
    EXPECT_TRUE(md_unit_is_none(md_attributes_get(&t, occ)->unit));

    // a zeroed unit is dimensionless, not an invalid state to guard against
    md_attribute_id_t z = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/anything"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = (md_unit_t){0}, .data = NULL, .byte_size = 0});
    EXPECT_TRUE(md_unit_is_none(md_attributes_get(&t, z)->unit));

    md_attributes_free(&t);
}

// A dipole has no index space to anchor it, so its origin is a sibling of the same shape.
// The two disagree on units, which is why the unit sits on the attribute and not the group.
UTEST(attributes, anchored_vector_group) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const double vector[3] = {0.0, 0.0, 1.85};
    const double origin[3] = {1.0, 2.0, 3.0};
    md_attribute_format_t vec3 = fmt_vec(MD_ATTRIBUTE_TYPE_F64, 1, 3);

    md_unit_t au = md_unit_scl(md_unit_mul(md_unit_mul(md_unit_ampere(), md_unit_second()), md_unit_bohr_radius()), 1.602176634e-19);

    md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("dipole/ground_state/vector"), .format = vec3, .unit = au, .data = vector, .byte_size = sizeof(vector)});
    md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("dipole/ground_state/origin"), .format = vec3, .unit = md_unit_angstrom(), .data = origin, .byte_size = sizeof(origin)});

    // one 3-vector, not three scalars: this is what components buys over an extra axis
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
    md_attribute_id_t f64 = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/f64"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F64, 4), .unit = md_unit_none(), .data = src, .byte_size = sizeof(src)});

    const int16_t ints[4] = {-3, 0, 7, 32767};
    md_attribute_id_t i16 = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/i16"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_I16, 4), .unit = md_unit_none(), .data = ints, .byte_size = sizeof(ints)});

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
    md_attribute_id_t v = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/vec"), .format = fmt_vec(MD_ATTRIBUTE_TYPE_F32, 2, 3), .unit = md_unit_none(), .data = vecs, .byte_size = sizeof(vecs)});
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, v), md_unit_none()), 6u);
    EXPECT_EQ(dst[5], 6.0f);

    // a destination which cannot hold it writes nothing
    EXPECT_EQ(md_attribute_extract_f32(dst, 3, md_attributes_get(&t, v), md_unit_none()), 0u);

    md_attributes_free(&t);
}

UTEST(attributes, extract_converts_unit) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const float lengths[2] = {10.0f, 2.5f};
    md_attribute_id_t len = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/length"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 2), .unit = md_unit_angstrom(), .data = lengths, .byte_size = sizeof(lengths)});

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
    md_attribute_id_t d = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("dipole/ground_state/vector"), .format = fmt_vec(MD_ATTRIBUTE_TYPE_F64, 1, 3), .unit = md_unit_elementary_charge_bohr(), .data = dipole, .byte_size = sizeof(dipole)});

    const float occ[2] = {0.5f, 1.0f};
    md_attribute_id_t o = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/occupancy"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 2), .unit = md_unit_none(), .data = occ, .byte_size = sizeof(occ)});

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

    md_attribute_id_t vel  = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/velocity"), .format = fmt_vec(MD_ATTRIBUTE_TYPE_F32, 4, 3), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    md_attribute_id_t mul  = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/charge/mulliken"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    md_attribute_id_t low  = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/charge/lowdin"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    md_attribute_id_t zed  = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atomic_number/z"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_U8,  4), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    md_attribute_id_t dip  = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("dipole/magnetic"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F64, 3), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    // sorts BETWEEN "atom" and "atom/..." because '-' is below '/', so a naive walk stops early
    md_attribute_id_t dash = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom-x/weird"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
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

    md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/velocity"), .format = fmt_vec(MD_ATTRIBUTE_TYPE_F32, 4, 3), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/charge/mulliken"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/charge/lowdin"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("dipole/magnetic"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F64, 3), .unit = md_unit_none(), .data = NULL, .byte_size = 0});

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

    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/charge"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 0});

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

    md_attribute_id_t vel = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/velocity"), .format = fmt_vec(MD_ATTRIBUTE_TYPE_F32, 4, 3), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    md_attribute_id_t mul = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/charge/mulliken"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
    md_attribute_id_t low = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/charge/lowdin"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 0});

    EXPECT_TRUE(md_attributes_remove(&t, mul));
    EXPECT_FALSE(md_attributes_remove(&t, mul));
    EXPECT_EQ(md_attributes_count(&t), 2u);
    EXPECT_TRUE(md_attributes_find(&t, STR_LIT("atom/charge/mulliken")) == NULL);

    md_attribute_id_t ids[8];
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("atom")), 2u);
    EXPECT_EQ(md_attributes_get(&t, ids[0])->id, low);
    EXPECT_EQ(md_attributes_get(&t, ids[1])->id, vel);

    // the id is the path hash, so recreating the same path yields the same id
    md_attribute_id_t again = md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("atom/charge/mulliken"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .data = NULL, .byte_size = 0});
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
    EXPECT_TRUE(str_eq(md_attribute_group(occ), STR_LIT("atom")));
    EXPECT_TRUE(str_eq(md_attribute_leaf(occ),    STR_LIT("occupancy")));

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

// The one mistake the old encoding could not catch: extents declared, value width forgotten.
UTEST(attributes, components_must_be_stated) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_format_t no_components = {.type = MD_ATTRIBUTE_TYPE_F32, .rank = 1, .shape = {10}};
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/b"), .format = no_components, .unit = md_unit_none()}), MD_ATTRIBUTE_INVALID);

    // An extent past the declared rank is the old {N,3} spelling surviving a port.
    md_attribute_format_t stale = {.type = MD_ATTRIBUTE_TYPE_F32, .components = 3, .rank = 1, .shape = {10, 3}};
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){.path = STR_LIT("a/c"), .format = stale, .unit = md_unit_none()}), MD_ATTRIBUTE_INVALID);

    EXPECT_EQ(md_attributes_count(&t), 0u);
    md_attributes_free(&t);
}

// {S,N} scalars is a shape the old encoding could not spell at all: two index axes and a scalar
// value. A slice of it is one state's per atom values.
UTEST(attributes, multi_axis_slice) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    // 3 states x 4 atoms, one scalar each
    const float charges[12] = {
        0.0f, 0.1f, 0.2f, 0.3f,
        1.0f, 1.1f, 1.2f, 1.3f,
        2.0f, 2.1f, 2.2f, 2.3f,
    };
    md_attribute_format_t per_state = {.type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 2, .shape = {3, 4}};
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("atom/charge/mulliken_per_state"), .format = per_state, .unit = md_unit_none(),
        .data = charges, .byte_size = sizeof(charges)});
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    const md_attribute_t* a = md_attributes_get(&t, id);
    EXPECT_EQ(md_attribute_components(&a->format),    1u);
    EXPECT_EQ(md_attribute_value_count(&a->format),   12u);
    EXPECT_EQ(md_attribute_element_count(&a->format), 12u);

    // Ask how big the slice is before allocating for it. The answer comes from the format, so it
    // does not depend on the data being there.
    const md_attribute_slice_t state1 = md_attribute_slice_1(1);
    EXPECT_EQ(md_attribute_slice_count(a, &state1), 4u);

    // and the slice has a format of its own: the fixed axis is gone, the value is untouched
    md_attribute_format_t slice_fmt = {0};
    ASSERT_TRUE(md_attribute_slice_format(&slice_fmt, a, &state1));
    EXPECT_EQ(slice_fmt.rank, 1u);
    EXPECT_EQ(slice_fmt.shape[0], 4u);
    EXPECT_EQ(slice_fmt.components, 1u);

    float dst[16] = {0};
    EXPECT_EQ(md_attribute_extract_slice_f32(dst, ARRAY_SIZE(dst), a, &state1, md_unit_none()), 4u);
    EXPECT_NEAR(dst[0], 1.0f, 1.0e-6f);
    EXPECT_NEAR(dst[3], 1.3f, 1.0e-6f);

    // fixing nothing is the whole attribute, which is exactly md_attribute_extract_f32
    const md_attribute_slice_t all = md_attribute_slice_all();
    EXPECT_EQ(md_attribute_slice_count(a, &all), 12u);
    EXPECT_EQ(md_attribute_extract_slice_f32(dst, ARRAY_SIZE(dst), a, &all, md_unit_none()), 12u);
    EXPECT_NEAR(dst[11], 2.3f, 1.0e-6f);

    // out of range is a refusal, not a clamp; so is asking for more axes than there are
    const md_attribute_slice_t state9 = md_attribute_slice_1(9);
    EXPECT_EQ(md_attribute_slice_count(a, &state9), 0u);
    EXPECT_EQ(md_attribute_extract_slice_f32(dst, ARRAY_SIZE(dst), a, &state9, md_unit_none()), 0u);

    md_attribute_slice_t three = {0};
    three.num_idx = 3;
    EXPECT_EQ(md_attribute_extract_slice_f32(dst, ARRAY_SIZE(dst), a, &three, md_unit_none()), 0u);

    md_attributes_free(&t);
}

// Slicing never splits a value: the components stay together whatever is fixed.
UTEST(attributes, slice_keeps_components) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    // 2 modes x 3 atoms, values 3 components wide
    float disp[18];
    for (int i = 0; i < 18; ++i) disp[i] = (float)i;

    md_attribute_format_t modes = {.type = MD_ATTRIBUTE_TYPE_F32, .components = 3, .rank = 2, .shape = {2, 3}};
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("atom/normal_mode"), .format = modes, .unit = md_unit_none(),
        .data = disp, .byte_size = sizeof(disp)});
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    const md_attribute_t* a = md_attributes_get(&t, id);
    EXPECT_EQ(md_attribute_value_count(&a->format),   6u);
    EXPECT_EQ(md_attribute_element_count(&a->format), 18u);

    float dst[32] = {0};
    const md_attribute_slice_t mode1 = md_attribute_slice_1(1);
    // one mode is 3 atoms x 3 components = 9 floats, starting at flat 9
    EXPECT_EQ(md_attribute_extract_slice_f32(dst, ARRAY_SIZE(dst), a, &mode1, md_unit_none()), 9u);
    EXPECT_NEAR(dst[0],  9.0f, 1.0e-6f);
    EXPECT_NEAR(dst[8], 17.0f, 1.0e-6f);

    // fixing both index axes yields exactly one value, components intact
    const md_attribute_slice_t mode1_atom2 = md_attribute_slice_2(1, 2);
    md_attribute_format_t one_value = {0};
    ASSERT_TRUE(md_attribute_slice_format(&one_value, a, &mode1_atom2));
    EXPECT_EQ(one_value.rank, 0u);
    EXPECT_EQ(one_value.components, 3u);
    EXPECT_EQ(md_attribute_extract_slice_f32(dst, ARRAY_SIZE(dst), a, &mode1_atom2, md_unit_none()), 3u);
    EXPECT_NEAR(dst[0], 15.0f, 1.0e-6f);
    EXPECT_NEAR(dst[2], 17.0f, 1.0e-6f);

    md_attributes_free(&t);
}

// Presentation only: stored, returned, and never used to find anything.
UTEST(attributes, label_and_description) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_id_t a = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("atom/charge/mulliken"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4),
        .unit = md_unit_none(), .label = STR_LIT("Mulliken"),
        .description = STR_LIT("Partial charge by Mulliken population analysis")});
    ASSERT_NE(a, MD_ATTRIBUTE_INVALID);

    // duplicate labels are legal and not a warning: nothing keys on them
    md_attribute_id_t b = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("atom/charge/lowdin"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4),
        .unit = md_unit_none(), .label = STR_LIT("Mulliken")});
    EXPECT_NE(b, MD_ATTRIBUTE_INVALID);

    EXPECT_TRUE(str_eq(md_attributes_get(&t, a)->label, STR_LIT("Mulliken")));
    EXPECT_TRUE(str_eq(md_attributes_get(&t, b)->label, STR_LIT("Mulliken")));
    EXPECT_TRUE(str_eq(md_attributes_get(&t, a)->description, STR_LIT("Partial charge by Mulliken population analysis")));

    // an omitted label is empty: a valid state, and the consumer's cue to prettify the leaf
    md_attribute_id_t c = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("atom/charge/hirshfeld"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none()});
    ASSERT_NE(c, MD_ATTRIBUTE_INVALID);
    EXPECT_TRUE(str_empty(md_attributes_get(&t, c)->label));
    EXPECT_TRUE(str_empty(md_attributes_get(&t, c)->description));

    // removal frees the strings along with the data
    EXPECT_TRUE(md_attributes_remove(&t, a));
    EXPECT_EQ(md_attributes_count(&t), 2u);

    md_attributes_free(&t);
}

// Doubles survive the round trip. A total energy needs more significant digits than a float carries,
// which is the whole reason the f64 path exists.
UTEST(attributes, extract_f64_keeps_precision) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    // Distinguishable in double, indistinguishable once narrowed to float.
    const double src[2] = { -76.0266327408, -76.0266327409 };
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("scf/energy"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F64, 2),
        .unit = md_unit_none(), .data = src, .byte_size = sizeof(src)});
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    const md_attribute_t* a = md_attributes_get(&t, id);

    double d[2] = {0};
    EXPECT_EQ(md_attribute_extract_f64(d, ARRAY_SIZE(d), a, md_unit_none()), 2u);
    EXPECT_EQ(d[0], src[0]);
    EXPECT_EQ(d[1], src[1]);
    EXPECT_TRUE(d[0] != d[1]);

    float f[2] = {0};
    EXPECT_EQ(md_attribute_extract_f32(f, ARRAY_SIZE(f), a, md_unit_none()), 2u);
    EXPECT_TRUE(f[0] == f[1]);  // the loss the f64 path exists to avoid

    md_attributes_free(&t);
}

// A slice of a {M,A} matrix is one row, which is how an MO coefficient vector comes out.
UTEST(attributes, extract_slice_f64_row) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    // 3 orbitals x 4 AOs
    const double coeff[12] = {
        0.0, 0.1, 0.2, 0.3,
        1.0, 1.1, 1.2, 1.3,
        2.0, 2.1, 2.2, 2.3,
    };
    md_attribute_format_t format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2, .shape = {3, 4}};
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("orbital/alpha/coefficient"), .format = format, .unit = md_unit_none(),
        .data = coeff, .byte_size = sizeof(coeff)});
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    const md_attribute_t* a = md_attributes_get(&t, id);

    double row[8] = {0};
    const md_attribute_slice_t mo = md_attribute_slice_1(2);
    EXPECT_EQ(md_attribute_slice_count(a, &mo), 4u);
    EXPECT_EQ(md_attribute_extract_slice_f64(row, ARRAY_SIZE(row), a, &mo, md_unit_none()), 4u);
    EXPECT_EQ(row[0], 2.0);
    EXPECT_EQ(row[3], 2.3);

    // the same slice against a sibling of the same shape: nothing about it is bound to one attribute
    const md_attribute_slice_t bad = md_attribute_slice_1(3);
    EXPECT_EQ(md_attribute_slice_count(a, &bad), 0u);
    EXPECT_EQ(md_attribute_extract_slice_f64(row, ARRAY_SIZE(row), a, &bad, md_unit_none()), 0u);

    md_attributes_free(&t);
}

// --- VIRTUAL ATTRIBUTES ---------------------------------------------------------------------

// 3 states x 4 atoms, value(state, atom) = state + atom / 10. Computed rather than stored, so the
// same formula must come back whether asked for whole or one state's row - that agreement is the
// point of the test, not the formula itself.
static size_t provider_state_ramp_f32(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
    (void)user_data;
    float* out = (float*)dst;
    size_t num_atoms = attr->format.shape[1];
    if (slice && slice->num_idx > 0) {
        // one fixed state, cap == num_atoms
        for (size_t atom = 0; atom < cap; ++atom) {
            out[atom] = (float)slice->idx[0] + (float)atom * 0.1f;
        }
    } else {
        for (size_t i = 0; i < cap; ++i) {
            out[i] = (float)(i / num_atoms) + (float)(i % num_atoms) * 0.1f;
        }
    }
    return cap;
}

// Reads a single constant out of user_data, so the provider itself carries no state.
static size_t provider_constant_f32(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
    (void)attr;
    (void)slice;
    float value = *(const float*)user_data;
    float* out = (float*)dst;
    for (size_t i = 0; i < cap; ++i) {
        out[i] = value;
    }
    return cap;
}

// A provider that never honours cap, to exercise the "provider lied" error path.
static size_t provider_wrong_count(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
    (void)dst;
    (void)attr;
    (void)slice;
    (void)user_data;
    return cap > 0 ? cap - 1 : 0;
}

UTEST(attributes, virtual_create_and_reject) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_format_t per_state = {.type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 2, .shape = {3, 4}};

    // no provider is a rejection, not a virtual attribute that always fails
    md_attribute_virtual_t no_provider = {0};
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("a/b"), .format = per_state, .unit = md_unit_none(), .virt = &no_provider}), MD_ATTRIBUTE_INVALID);

    // virt and resident data at the same time is ambiguous about which one is authoritative
    const float data[12] = {0};
    md_attribute_virtual_t virt = {.provider = provider_state_ramp_f32};
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("a/b"), .format = per_state, .unit = md_unit_none(), .virt = &virt, .data = data, .byte_size = sizeof(data)}), MD_ATTRIBUTE_INVALID);
    EXPECT_EQ(md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("a/b"), .format = per_state, .unit = md_unit_none(), .virt = &virt, .byte_size = sizeof(data)}), MD_ATTRIBUTE_INVALID);

    EXPECT_EQ(md_attributes_count(&t), 0u);

    // the correct call publishes a virtual attribute with no resident storage
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("atom/charge/mulliken_per_state"), .format = per_state, .unit = md_unit_none(), .virt = &virt});
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    const md_attribute_t* a = md_attributes_get(&t, id);
    ASSERT_TRUE(a != NULL);
    EXPECT_EQ(a->storage, MD_ATTRIBUTE_STORAGE_VIRTUAL);
    EXPECT_TRUE(a->data == NULL);

    md_attributes_free(&t);
}

// The provider must be reachable through both the whole-attribute and the sliced extraction path,
// and the two must agree, since nothing else here would catch them disagreeing.
UTEST(attributes, virtual_extract_whole_and_slice_agree) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_format_t per_state = {.type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 2, .shape = {3, 4}};
    md_attribute_virtual_t virt = {.provider = provider_state_ramp_f32};
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("atom/charge/mulliken_per_state"), .format = per_state, .unit = md_unit_none(), .virt = &virt});
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    const md_attribute_t* a = md_attributes_get(&t, id);

    float whole[12] = {0};
    EXPECT_EQ(md_attribute_extract_f32(whole, ARRAY_SIZE(whole), a, md_unit_none()), 12u);
    EXPECT_NEAR(whole[0],  0.0f, 1.0e-6f);
    EXPECT_NEAR(whole[11], 2.3f, 1.0e-6f);

    float row[4] = {0};
    const md_attribute_slice_t state1 = md_attribute_slice_1(1);
    EXPECT_EQ(md_attribute_extract_slice_f32(row, ARRAY_SIZE(row), a, &state1, md_unit_none()), 4u);
    EXPECT_NEAR(row[0], whole[4], 1.0e-6f);
    EXPECT_NEAR(row[3], whole[7], 1.0e-6f);

    md_attributes_free(&t);
}

// Unit conversion happens centrally regardless of where the values came from, so a virtual
// attribute refuses and rescales exactly like a resident one.
UTEST(attributes, virtual_extract_converts_unit) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    float length_angstrom = 10.0f;
    md_attribute_virtual_t virt = {.provider = provider_constant_f32, .user_data = &length_angstrom};
    md_attribute_format_t single = {.type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 1, .shape = {2}};
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("a/length"), .format = single, .unit = md_unit_angstrom(), .virt = &virt});
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    const md_attribute_t* a = md_attributes_get(&t, id);

    float dst[2] = {0};
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), a, md_unit_nanometer()), 2u);
    EXPECT_NEAR(dst[0], 1.0f, 1.0e-6f);

    // a dipole is not a length, virtual or otherwise
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), a, md_unit_debye()), 0u);

    md_attributes_free(&t);
}

// A provider that does not honour cap must fail the whole extraction, not return a partial buffer.
UTEST(attributes, virtual_provider_lying_about_count_fails) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_virtual_t virt = {.provider = provider_wrong_count};
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("a/b"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4), .unit = md_unit_none(), .virt = &virt});
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    float dst[4] = {-1.0f, -1.0f, -1.0f, -1.0f};
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, id), md_unit_none()), 0u);
    EXPECT_EQ(dst[0], -1.0f);   // nothing written

    md_attributes_free(&t);
}

// There is no buffer to hand out for a computed attribute, so this is a refusal, not a NULL that
// happens to also mean "not found".
UTEST(attributes, virtual_attributes_data_returns_null) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    md_attribute_virtual_t virt = {.provider = provider_state_ramp_f32};
    md_attribute_format_t per_state = {.type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 2, .shape = {3, 4}};
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("a/b"), .format = per_state, .unit = md_unit_none(), .virt = &virt});
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    EXPECT_TRUE(md_attributes_data(&t, id, MD_ATTRIBUTE_TYPE_F32) == NULL);

    md_attributes_free(&t);
}

// user_data allocated through the table is owned by it and must survive to every provider call;
// a borrowed pointer must be left alone by removal and by freeing the whole table.
UTEST(attributes, virtual_user_data_lifecycle) {
    md_attributes_t no_alloc = {0};
    EXPECT_TRUE(md_attributes_alloc_user_data(&no_alloc, sizeof(float)) == NULL);

    md_attributes_t t = {.alloc = md_get_heap_allocator()};
    EXPECT_TRUE(md_attributes_alloc_user_data(&t, 0) == NULL);

    float* owned = (float*)md_attributes_alloc_user_data(&t, sizeof(float));
    ASSERT_TRUE(owned != NULL);
    *owned = 42.0f;

    md_attribute_virtual_t owned_virt = {.provider = provider_constant_f32, .user_data = owned, .user_data_size = sizeof(float)};
    md_attribute_id_t owned_id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("a/owned"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 1), .unit = md_unit_none(), .virt = &owned_virt});
    ASSERT_NE(owned_id, MD_ATTRIBUTE_INVALID);

    float dst[1] = {0};
    EXPECT_EQ(md_attribute_extract_f32(dst, 1, md_attributes_get(&t, owned_id), md_unit_none()), 1u);
    EXPECT_EQ(dst[0], 42.0f);

    // a borrowed pointer (user_data_size 0) must not be freed by remove or by md_attributes_free
    float borrowed = 7.0f;
    md_attribute_virtual_t borrowed_virt = {.provider = provider_constant_f32, .user_data = &borrowed, .user_data_size = 0};
    md_attribute_id_t borrowed_id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("a/borrowed"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 1), .unit = md_unit_none(), .virt = &borrowed_virt});
    ASSERT_NE(borrowed_id, MD_ATTRIBUTE_INVALID);

    EXPECT_TRUE(md_attributes_remove(&t, borrowed_id));
    EXPECT_EQ(borrowed, 7.0f);  // still alive after removal

    md_attributes_free(&t);
    EXPECT_EQ(borrowed, 7.0f);  // still alive after the table is gone
}

// A second name for one datum: same storage, same numbers, its own label - and NOT a copy.
UTEST(attributes, alias_shares_storage) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const float charges[4] = {0.5f, -0.25f, 0.125f, 1.0f};
    md_attribute_id_t src = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("vlx/atom/mulliken"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 4),
        .unit = md_unit_none(), .label = STR_LIT("Mulliken (VeloxChem)"),
        .data = charges, .byte_size = sizeof(charges)});
    ASSERT_NE(src, MD_ATTRIBUTE_INVALID);

    md_attribute_id_t alias = md_attributes_alias(&t, src, STR_LIT("atom/charge/mulliken"), STR_LIT("Mulliken"), (str_t){0});
    ASSERT_NE(alias, MD_ATTRIBUTE_INVALID);
    EXPECT_EQ(md_attributes_count(&t), 2u);

    const md_attribute_t* a = md_attributes_get(&t, src);
    const md_attribute_t* b = md_attributes_get(&t, alias);
    ASSERT_TRUE(a && b);

    // Same datum, different names. The id differs - that is exactly why comparing ids is the wrong
    // question and md_attribute_same_data is the right one.
    EXPECT_NE(a->id, b->id);
    EXPECT_TRUE(md_attribute_same_data(a, b));
    EXPECT_TRUE(md_attribute_same_data(a, a));

    // Shared, not copied: the resident fast path works through either name and finds one buffer.
    EXPECT_TRUE(a->data == b->data);
    EXPECT_EQ(b->storage, MD_ATTRIBUTE_STORAGE_ALIAS);

    // Format and unit come from the target; the label is the alias's own.
    EXPECT_EQ(md_attribute_element_count(&b->format), 4u);
    EXPECT_TRUE(str_eq(b->label, STR_LIT("Mulliken")));
    EXPECT_TRUE(str_eq(a->label, STR_LIT("Mulliken (VeloxChem)")));

    // and it reads as data, through the name a consumer of the neutral path would use
    float dst[4] = {0};
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), b, md_unit_none()), 4u);
    EXPECT_NEAR(dst[1], -0.25f, 1.0e-6f);

    // Both names are in the table and the prefix queries see each under its own group.
    md_attribute_id_t ids[8];
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("atom")), 1u);
    EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t, STR_LIT("vlx")), 1u);

    md_attributes_free(&t);
}

// An alias cannot outlive the storage it reads, so removing the target takes it with it.
UTEST(attributes, alias_cannot_outlive_its_target) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const float v[2] = {1.0f, 2.0f};
    md_attribute_id_t src = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("vlx/a"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 2),
        .unit = md_unit_none(), .data = v, .byte_size = sizeof(v)});
    ASSERT_NE(src, MD_ATTRIBUTE_INVALID);

    md_attribute_id_t one = md_attributes_alias(&t, src, STR_LIT("a/one"), (str_t){0}, (str_t){0});
    md_attribute_id_t two = md_attributes_alias(&t, src, STR_LIT("a/two"), (str_t){0}, (str_t){0});
    ASSERT_NE(one, MD_ATTRIBUTE_INVALID);
    ASSERT_NE(two, MD_ATTRIBUTE_INVALID);

    // Aliasing an alias flattens: three names, one owner.
    md_attribute_id_t three = md_attributes_alias(&t, two, STR_LIT("a/three"), (str_t){0}, (str_t){0});
    ASSERT_NE(three, MD_ATTRIBUTE_INVALID);
    EXPECT_TRUE(md_attribute_same_data(md_attributes_get(&t, three), md_attributes_get(&t, src)));
    EXPECT_EQ(md_attributes_get(&t, three)->root, src);

    EXPECT_EQ(md_attributes_count(&t), 4u);
    EXPECT_TRUE(md_attributes_remove(&t, src));
    EXPECT_EQ(md_attributes_count(&t), 0u);   // every name went with the datum

    md_attributes_free(&t);
}

// Removing an alias leaves the datum, and every other name for it, alone.
UTEST(attributes, removing_an_alias_leaves_the_target) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const float v[2] = {3.0f, 4.0f};
    md_attribute_id_t src = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("vlx/b"), .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F32, 2),
        .unit = md_unit_none(), .data = v, .byte_size = sizeof(v)});
    md_attribute_id_t alias = md_attributes_alias(&t, src, STR_LIT("b/neutral"), (str_t){0}, (str_t){0});
    ASSERT_NE(alias, MD_ATTRIBUTE_INVALID);

    EXPECT_TRUE(md_attributes_remove(&t, alias));
    EXPECT_EQ(md_attributes_count(&t), 1u);

    float dst[2] = {0};
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, src), md_unit_none()), 2u);
    EXPECT_NEAR(dst[0], 3.0f, 1.0e-6f);

    // A duplicate path is refused whichever way it is published, and so is a missing target.
    EXPECT_EQ(md_attributes_alias(&t, src, STR_LIT("vlx/b"), (str_t){0}, (str_t){0}), MD_ATTRIBUTE_INVALID);
    EXPECT_EQ(md_attributes_alias(&t, alias, STR_LIT("b/gone"), (str_t){0}, (str_t){0}), MD_ATTRIBUTE_INVALID);

    // Writing goes through the owner, never through a second name for it.
    md_attribute_id_t again = md_attributes_alias(&t, src, STR_LIT("b/neutral"), (str_t){0}, (str_t){0});
    ASSERT_NE(again, MD_ATTRIBUTE_INVALID);
    EXPECT_TRUE(md_attributes_data(&t, again, MD_ATTRIBUTE_TYPE_F32) == NULL);
    EXPECT_TRUE(md_attributes_data(&t, src,   MD_ATTRIBUTE_TYPE_F32) != NULL);

    md_attributes_free(&t);
}

// An id is a function of the PATH, which is what makes it the right thing for a consumer to store.
// It resolves before anything is published, it is the id create hands back, and it survives the
// attribute being removed and published again - none of which an index into a list can do.
UTEST(attributes, id_from_path_resolves_without_the_table) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    const str_t path = STR_LIT("vlx/density_property/relaxed");
    const md_attribute_id_t expected = md_attributes_id_from_path(path);

    EXPECT_NE(expected, MD_ATTRIBUTE_INVALID);
    EXPECT_TRUE(md_attributes_get(&t, expected) == NULL);   // names a path, not a thing that exists

    const double v[4] = {1.0, 2.0, 3.0, 4.0};
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = path, .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F64, 4),
        .unit = md_unit_none(), .data = v, .byte_size = sizeof(v)});
    EXPECT_EQ(id, expected);

    // And again after a reload, which is the case the whole thing exists for.
    EXPECT_TRUE(md_attributes_remove(&t, id));
    id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = path, .format = fmt_scalars(MD_ATTRIBUTE_TYPE_F64, 4),
        .unit = md_unit_none(), .data = v, .byte_size = sizeof(v)});
    EXPECT_EQ(id, expected);

    // Different paths are different ids; the same path is the same id however it is spelled in the
    // caller's own storage.
    EXPECT_NE(md_attributes_id_from_path(STR_LIT("vlx/density_property/relaxed2")), expected);
    EXPECT_EQ(md_attributes_id_from_path(str_from_cstr("vlx/density_property/relaxed")), expected);

    md_attributes_free(&t);
}

// A rank 3 attribute sliced by its leading axis is a matrix, and the format says so without the
// caller redoing the arithmetic. This is what the density and NTO paths rely on.
UTEST(attributes, slice_format_narrows_rank_3_to_a_matrix) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    double v[2 * 3 * 3] = {0};
    for (int i = 0; i < 2 * 3 * 3; ++i) v[i] = (double)i;

    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 3, .shape = {2, 3, 3},
    };
    md_attribute_id_t id = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("rsp/transition_density"), .format = format,
        .unit = md_unit_none(), .data = v, .byte_size = sizeof(v)});
    ASSERT_NE(id, MD_ATTRIBUTE_INVALID);

    const md_attribute_t* attr = md_attributes_get(&t, id);
    md_attribute_slice_t slice = md_attribute_slice_1(1);

    md_attribute_format_t sliced = {0};
    ASSERT_TRUE(md_attribute_slice_format(&sliced, attr, &slice));
    EXPECT_EQ(sliced.rank, 2u);
    EXPECT_EQ(sliced.shape[0], 3u);
    EXPECT_EQ(sliced.shape[1], 3u);
    EXPECT_EQ(md_attribute_slice_count(attr, &slice), 9u);

    double dst[9] = {0};
    EXPECT_EQ(md_attribute_extract_slice_f64(dst, ARRAY_SIZE(dst), attr, &slice, md_unit_none()), 9u);
    EXPECT_NEAR(dst[0], 9.0, 1.0e-12);   // plane 1 starts at 1 * 3 * 3
    EXPECT_NEAR(dst[8], 17.0, 1.0e-12);

    // No slice at all is the whole thing, still rank 3.
    md_attribute_format_t whole = {0};
    ASSERT_TRUE(md_attribute_slice_format(&whole, attr, NULL));
    EXPECT_EQ(whole.rank, 3u);
    EXPECT_EQ(md_attribute_slice_count(attr, NULL), 18u);

    // An index past the axis selects nothing rather than reading past the end.
    md_attribute_slice_t bad = md_attribute_slice_1(2);
    EXPECT_EQ(md_attribute_slice_count(attr, &bad), 0u);
    EXPECT_FALSE(md_attribute_slice_format(&whole, attr, &bad));

    md_attributes_free(&t);
}

// A provider that reads two OTHER attributes out of the table and combines them. This is the shape
// the SCF total/difference densities have: a derivation over derivations, which is legal as long as
// the graph stays acyclic.
typedef struct combine_ctx_t {
    md_attributes_t* table;
    str_t            lhs;
    str_t            rhs;
    double           rhs_scale;
} combine_ctx_t;

static size_t provider_combine_f32(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
    (void)attr;
    (void)slice;
    combine_ctx_t* ctx = (combine_ctx_t*)user_data;

    const md_attribute_t* lhs = md_attributes_find(ctx->table, ctx->lhs);
    const md_attribute_t* rhs = md_attributes_find(ctx->table, ctx->rhs);
    if (!lhs || !rhs) return 0;

    md_temp_scope_t temp = md_temp_begin();
    float* rhs_data = md_temp_alloc_array(temp, float, cap);
    bool ok = rhs_data
           && md_attribute_extract_f32((float*)dst, cap, lhs, md_unit_none()) == cap
           && md_attribute_extract_f32(rhs_data,    cap, rhs, md_unit_none()) == cap;
    if (ok) {
        float* out = (float*)dst;
        for (size_t i = 0; i < cap; ++i) out[i] += (float)ctx->rhs_scale * rhs_data[i];
    }
    md_temp_end(temp);
    return ok ? cap : 0;
}

// An ALIAS is a second NAME, not a different kind of storage - so aliasing a VIRTUAL attribute has
// to read through the same provider. It did not: the extract branched on the storage tag, sent the
// alias down the resident path, found no data and returned 0 with no diagnostic. Every consumer of
// the alias saw an attribute that existed and held nothing.
UTEST(attributes, alias_of_a_virtual_attribute_reads_through_its_provider) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    float constant = 2.5f;
    md_attribute_virtual_t virt = {.provider = provider_constant_f32, .user_data = &constant};
    md_attribute_format_t fmt = {.type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 1, .shape = {4}};

    md_attribute_id_t alpha = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("orbital/alpha/density"), .format = fmt, .unit = md_unit_none(), .virt = &virt});
    ASSERT_NE(alpha, MD_ATTRIBUTE_INVALID);

    md_attribute_id_t beta = md_attributes_alias(&t, alpha, STR_LIT("orbital/beta/density"), (str_t){0}, (str_t){0});
    ASSERT_NE(beta, MD_ATTRIBUTE_INVALID);

    float dst[4] = {0};
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, beta), md_unit_none()), 4u);
    EXPECT_NEAR(dst[0], 2.5f, 1.0e-6f);
    EXPECT_NEAR(dst[3], 2.5f, 1.0e-6f);

    // And through a slice, which takes the same path.
    md_attribute_slice_t all = md_attribute_slice_all();
    MEMSET(dst, 0, sizeof(dst));
    EXPECT_EQ(md_attribute_extract_slice_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, beta), &all, md_unit_none()), 4u);
    EXPECT_NEAR(dst[2], 2.5f, 1.0e-6f);

    // Writing still goes through the owner, and neither name has resident storage to hand out.
    EXPECT_TRUE(md_attributes_data(&t, beta,  MD_ATTRIBUTE_TYPE_F32) == NULL);
    EXPECT_TRUE(md_attributes_data(&t, alpha, MD_ATTRIBUTE_TYPE_F32) == NULL);

    md_attributes_free(&t);
}

// The case the bug actually showed up in: a virtual attribute combining a virtual one with an ALIAS
// of it. In a restricted SCF calculation beta IS alpha, so the total comes out as twice alpha and
// the difference as exactly zero - with no special case anywhere.
UTEST(attributes, virtual_over_an_alias_of_a_virtual) {
    md_attributes_t t = {.alloc = md_get_heap_allocator()};

    float constant = 3.0f;
    md_attribute_virtual_t spin_virt = {.provider = provider_constant_f32, .user_data = &constant};
    md_attribute_format_t fmt = {.type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 1, .shape = {4}};

    md_attribute_id_t alpha = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("orbital/alpha/density"), .format = fmt, .unit = md_unit_none(), .virt = &spin_virt});
    ASSERT_NE(md_attributes_alias(&t, alpha, STR_LIT("orbital/beta/density"), (str_t){0}, (str_t){0}), MD_ATTRIBUTE_INVALID);

    combine_ctx_t total_ctx = {&t, STR_LIT("orbital/alpha/density"), STR_LIT("orbital/beta/density"),  1.0};
    combine_ctx_t diff_ctx  = {&t, STR_LIT("orbital/alpha/density"), STR_LIT("orbital/beta/density"), -1.0};

    md_attribute_virtual_t total_virt = {.provider = provider_combine_f32, .user_data = &total_ctx};
    md_attribute_virtual_t diff_virt  = {.provider = provider_combine_f32, .user_data = &diff_ctx};

    md_attribute_id_t total = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("orbital/total/density"), .format = fmt, .unit = md_unit_none(), .virt = &total_virt});
    md_attribute_id_t diff = md_attributes_create(&t, &(md_attribute_desc_t){
        .path = STR_LIT("orbital/difference/density"), .format = fmt, .unit = md_unit_none(), .virt = &diff_virt});
    ASSERT_NE(total, MD_ATTRIBUTE_INVALID);
    ASSERT_NE(diff,  MD_ATTRIBUTE_INVALID);

    float dst[4] = {0};
    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, total), md_unit_none()), 4u);
    EXPECT_NEAR(dst[0], 6.0f, 1.0e-6f);
    EXPECT_NEAR(dst[3], 6.0f, 1.0e-6f);

    EXPECT_EQ(md_attribute_extract_f32(dst, ARRAY_SIZE(dst), md_attributes_get(&t, diff), md_unit_none()), 4u);
    EXPECT_NEAR(dst[0], 0.0f, 1.0e-6f);
    EXPECT_NEAR(dst[3], 0.0f, 1.0e-6f);

    md_attributes_free(&t);
}
