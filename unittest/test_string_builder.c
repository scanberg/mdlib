#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_str_builder.h>

UTEST(string_builder, all) {
	md_str_builder_t sb = {0};
	md_str_builder_init(&sb, default_allocator);
	str_t str;

	str = md_str_builder_to_str(&sb);

	EXPECT_EQ(NULL, str.ptr);
	EXPECT_EQ(0, str.len);

	md_str_builder_append_cstr(&sb, "Hej");
	str = md_str_builder_to_str(&sb);

	EXPECT_STREQ("Hej", str.ptr);

	md_str_builder_printf(&sb, " %d Cool", 2);
	str = md_str_builder_to_str(&sb);

	EXPECT_STREQ("Hej 2 Cool", str.ptr);

	md_str_builder_append_str(&sb, STR(" 4 School."));
	str = md_str_builder_to_str(&sb);

	EXPECT_STREQ("Hej 2 Cool 4 School.", str.ptr);

	md_str_builder_append_char(&sb, '!');
	str = md_str_builder_to_str(&sb);

	EXPECT_STREQ("Hej 2 Cool 4 School.!", str.ptr);

	md_str_builder_reset(&sb);
	str = md_str_builder_to_str(&sb);

	EXPECT_EQ(NULL, str.ptr);
	EXPECT_EQ(0, str.len);

	md_str_builder_free(&sb);
}
