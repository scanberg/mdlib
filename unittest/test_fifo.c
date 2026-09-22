#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_fifo.h>

UTEST(fifo, create_empty) {
    md_fifo_t fifo = md_fifo_create(8, md_get_heap_allocator());

    EXPECT_TRUE(md_fifo_empty(&fifo));
    EXPECT_FALSE(md_fifo_full(&fifo));
    EXPECT_GE(fifo.capacity, (uint32_t)8);

    md_fifo_free(&fifo);
}

UTEST(fifo, push_pop_order) {
    md_fifo_t fifo = md_fifo_create(8, md_get_heap_allocator());

    // md_fifo_create enforces a minimum capacity, so fill to whatever capacity was actually granted.
    const int count = (int)fifo.capacity;
    for (int i = 0; i < count; ++i) {
        md_fifo_push(&fifo, i);
    }

    EXPECT_TRUE(md_fifo_full(&fifo));

    for (int i = 0; i < count; ++i) {
        ASSERT_FALSE(md_fifo_empty(&fifo));
        EXPECT_EQ(md_fifo_pop(&fifo), i);
    }

    EXPECT_TRUE(md_fifo_empty(&fifo));

    md_fifo_free(&fifo);
}

UTEST(fifo, clear_resets_state) {
    md_fifo_t fifo = md_fifo_create(8, md_get_heap_allocator());

    md_fifo_push(&fifo, 1);
    md_fifo_push(&fifo, 2);
    md_fifo_clear(&fifo);

    EXPECT_TRUE(md_fifo_empty(&fifo));

    md_fifo_push(&fifo, 42);
    EXPECT_EQ(md_fifo_pop(&fifo), 42);

    md_fifo_free(&fifo);
}

// Regression test for a bug where growing the ring buffer while its contents
// were wrapped (tail > head) corrupted subsequent pop order, since the old
// implementation only reallocated the backing array without re-linearizing
// the wrapped entries relative to the new capacity.
UTEST(fifo, grow_preserves_order_when_wrapped) {
    md_fifo_t fifo = md_fifo_create(8, md_get_heap_allocator());
    const int capacity = (int)fifo.capacity;

    // Fill to capacity, then drain most of it so head sits near the end of the
    // buffer and tail sits a few slots behind it (tail < head, small remaining size).
    for (int i = 0; i < capacity; ++i) md_fifo_push(&fifo, i);
    for (int i = 0; i < capacity - 2; ++i) EXPECT_EQ(md_fifo_pop(&fifo), i);

    // head==capacity-0 wrapped to 0 only if pushes wrapped; with no further pushes yet
    // head == capacity & (capacity-1) == 0, tail == capacity-2. Pushing now wraps head
    // around the end of the buffer before the fifo becomes full and grows.
    const int values[] = {100, 101, 102, 103, 104, 105, 106, 107, 108, 109};
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i) {
        md_fifo_push(&fifo, values[i]);
    }

    // Remaining original entries must come out before the newly pushed ones.
    EXPECT_EQ(md_fifo_pop(&fifo), capacity - 2);
    EXPECT_EQ(md_fifo_pop(&fifo), capacity - 1);

    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i) {
        ASSERT_FALSE(md_fifo_empty(&fifo));
        EXPECT_EQ(md_fifo_pop(&fifo), values[i]);
    }

    EXPECT_TRUE(md_fifo_empty(&fifo));

    md_fifo_free(&fifo);
}

UTEST(fifo, many_push_pop_cycles_stress) {
    md_fifo_t fifo = md_fifo_create(16, md_get_heap_allocator());

    int next_push = 0;
    int next_pop = 0;

    // Repeatedly push more than we pop to force several grow events, then
    // drain fully and verify strict FIFO ordering held throughout.
    for (int round = 0; round < 200; ++round) {
        const int push_count = 1 + (round % 5);
        for (int i = 0; i < push_count; ++i) {
            md_fifo_push(&fifo, next_push++);
        }

        const int pop_count = round % 3;
        for (int i = 0; i < pop_count; ++i) {
            if (md_fifo_empty(&fifo)) break;
            EXPECT_EQ(md_fifo_pop(&fifo), next_pop++);
        }
    }

    while (!md_fifo_empty(&fifo)) {
        EXPECT_EQ(md_fifo_pop(&fifo), next_pop++);
    }

    EXPECT_EQ(next_pop, next_push);

    md_fifo_free(&fifo);
}
