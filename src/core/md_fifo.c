#include <core/md_fifo.h>

#include <core/md_common.h>
#include <core/md_intrinsics.h>

bool md_fifo_empty(const md_fifo_t* fifo) {
    ASSERT(fifo);
    return MD_FIFO_RAW_EMPTY(fifo->size);
}

bool md_fifo_full(const md_fifo_t* fifo) {
    ASSERT(fifo);
    return MD_FIFO_RAW_FULL(fifo->size, fifo->capacity);
}

md_fifo_t md_fifo_create(size_t capacity, md_allocator_i* alloc) {
    ASSERT(alloc);

    md_fifo_t fifo = {0};
    MD_FIFO_RAW_INIT_WITH_ALLOC(fifo.capacity, fifo.size, fifo.head, fifo.tail, fifo.data, fifo.alloc, alloc);
    MD_FIFO_RAW_ENSURE(fifo.capacity, fifo.size, fifo.head, fifo.tail, fifo.data, MAX((size_t)16, capacity), fifo.alloc);

#if DEBUG
    MEMSET(fifo.data, 0, sizeof(int) * fifo.capacity);
#endif

    return fifo;
}

void md_fifo_free(md_fifo_t* fifo) {
    ASSERT(fifo);
    ASSERT(fifo->alloc);

    MD_FIFO_RAW_FREE(fifo->capacity, fifo->size, fifo->head, fifo->tail, fifo->data, fifo->alloc);
    MEMSET(fifo, 0, sizeof(*fifo));
}

void md_fifo_clear(md_fifo_t* fifo) {
    ASSERT(fifo);

    MD_FIFO_RAW_CLEAR(fifo->size, fifo->head, fifo->tail);

#if DEBUG
    if (fifo->data) {
        MEMSET(fifo->data, 0, sizeof(int) * fifo->capacity);
    }
#endif
}

void md_fifo_push(md_fifo_t* fifo, int value) {
    ASSERT(fifo);
    ASSERT(fifo->data);

    MD_FIFO_RAW_PUSH(fifo->capacity, fifo->size, fifo->head, fifo->tail, fifo->data, value, fifo->alloc);
}

int md_fifo_pop(md_fifo_t* fifo) {
    ASSERT(fifo);
    ASSERT(!md_fifo_empty(fifo));

    const int value = MD_FIFO_RAW_FRONT(fifo->data, fifo->tail);
    MD_FIFO_RAW_POP(fifo->capacity, fifo->size, fifo->head, fifo->tail);

    return value;
}
