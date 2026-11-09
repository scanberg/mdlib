#include <core/md_fifo.h>

#include <core/md_common.h>

bool md_fifo_empty(const md_fifo_t* fifo) {
    MD_ASSERT(fifo);
    return MD_FIFO_RAW_EMPTY(fifo->size);
}

bool md_fifo_full(const md_fifo_t* fifo) {
    MD_ASSERT(fifo);
    return MD_FIFO_RAW_FULL(fifo->size, fifo->capacity);
}

md_fifo_t md_fifo_create(size_t capacity, md_allocator_i* alloc) {
    MD_ASSERT(alloc);

    md_fifo_t fifo = {0};
    MD_FIFO_RAW_INIT_WITH_ALLOC(fifo.capacity, fifo.size, fifo.head, fifo.tail, fifo.data, fifo.alloc, alloc);
    MD_FIFO_RAW_ENSURE(fifo.capacity, fifo.size, fifo.head, fifo.tail, fifo.data, MD_MAX((size_t)16, capacity), fifo.alloc);

#if MD_DEBUG
    md_memset(fifo.data, 0, sizeof(int) * fifo.capacity);
#endif

    return fifo;
}

void md_fifo_free(md_fifo_t* fifo) {
    MD_ASSERT(fifo);
    MD_ASSERT(fifo->alloc);

    MD_FIFO_RAW_FREE(fifo->capacity, fifo->size, fifo->head, fifo->tail, fifo->data, fifo->alloc);
    md_memset(fifo, 0, sizeof(*fifo));
}

void md_fifo_clear(md_fifo_t* fifo) {
    MD_ASSERT(fifo);

    MD_FIFO_RAW_CLEAR(fifo->size, fifo->head, fifo->tail);

#if MD_DEBUG
    if (fifo->data) {
        md_memset(fifo->data, 0, sizeof(int) * fifo->capacity);
    }
#endif
}

void md_fifo_push(md_fifo_t* fifo, int value) {
    MD_ASSERT(fifo);
    MD_ASSERT(fifo->data);

    MD_FIFO_RAW_PUSH(fifo->capacity, fifo->size, fifo->head, fifo->tail, fifo->data, value, fifo->alloc);
}

int md_fifo_pop(md_fifo_t* fifo) {
    MD_ASSERT(fifo);
    MD_ASSERT(!md_fifo_empty(fifo));

    const int value = MD_FIFO_RAW_FRONT(fifo->data, fifo->tail);
    MD_FIFO_RAW_POP(fifo->capacity, fifo->size, fifo->head, fifo->tail);

    return value;
}
