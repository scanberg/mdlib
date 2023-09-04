#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_str_builder.h>

UTEST(string_builder, all) {
	md_strb_t sb = {0};
	md_strb_init(&sb, md_heap_allocator);
	str_t str;

	str = md_strb_to_str(&sb);

	EXPECT_EQ(NULL, str.ptr);
	EXPECT_EQ(0, str.len);

	md_strb_push_cstr(&sb, "Hej");
	str = md_strb_to_str(&sb);

	EXPECT_STREQ("Hej", str.ptr);

	md_strb_fmt(&sb, " %d Cool", 2);
	str = md_strb_to_str(&sb);

	EXPECT_STREQ("Hej 2 Cool", str.ptr);

	md_strb_push_str(&sb, STR(" 4 School."));
	str = md_strb_to_str(&sb);

	EXPECT_STREQ("Hej 2 Cool 4 School.", str.ptr);

	md_strb_push_char(&sb, '!');
	str = md_strb_to_str(&sb);

	EXPECT_STREQ("Hej 2 Cool 4 School.!", str.ptr);

	md_strb_reset(&sb);
	str = md_strb_to_str(&sb);

	EXPECT_EQ(NULL, str.ptr);
	EXPECT_EQ(0, str.len);

	md_strb_free(&sb);
}
