#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_string_builder.h>

UTEST(string_builder, all) {
	md_string_builder_t sb = {0};
	md_string_builder_init(&sb, default_allocator);
	str_t str;

	str = md_string_builder_to_string(&sb);

	EXPECT_EQ(0, str.ptr);
	EXPECT_EQ(0, str.len);

	md_string_builder_print(&sb, "Hej");
	str = md_string_builder_to_string(&sb);

	EXPECT_STREQ("Hej", str.ptr);

	md_string_builder_printf(&sb, " %d Cool", 2);
	str = md_string_builder_to_string(&sb);

	EXPECT_STREQ("Hej 2 Cool", str.ptr);

	md_string_builder_append_str(&sb, MAKE_STR(" 4 School."));
	str = md_string_builder_to_string(&sb);

	EXPECT_STREQ("Hej 2 Cool 4 School.", str.ptr);

	md_string_builder_reset(&sb);
	str = md_string_builder_to_string(&sb);

	EXPECT_EQ(NULL, str.ptr);
	EXPECT_EQ(0, str.len);

	md_string_builder_free(&sb);
}
