#pragma once
#include <core/std.h>

struct Block {
    void* cursor;
    void* end;
    Block* next;
};

// Each thread gets one large block of memory for allocations
class Arena {
   public:
    static void* add(void* a, size_t b) { return reinterpret_cast<void*>(reinterpret_cast<size_t>(a) + b); }

    static void* sub(void* a, size_t b) { return reinterpret_cast<void*>(reinterpret_cast<size_t>(a) - b); }

    static size_t sub(void* a, void* b) { return reinterpret_cast<size_t>(a) - reinterpret_cast<size_t>(b); }

    // per thread capacity in bytes
    Arena(size_t capacity) : m_capacity(capacity) {}

    // each thread must call initThread() before allocate()
    void initThread() {
        Block*& block = getBlock();
        if (!block) {
            block = new Block;
            block->cursor = malloc(m_capacity);
            block->end = add(block->cursor, m_capacity);

            block->next = m_blocks.load();
            while (!m_blocks.compare_exchange_weak(block->next, block)) {
            }
        }
    }

    // not thread safe with resetGlobal()
    void* allocate(size_t size) {
        Block* block = getBlock();

        // round up to ensure alignment
        constexpr size_t ALIGN = 4;
        size = (size + ALIGN - 1) & ~(ALIGN - 1);
        void* alloc = block->cursor;
        assert(add(alloc, size) <= block->end);
        block->cursor = add(alloc, size);
        return alloc;
    }

    // rewind cursor in the current thread
    void resetThread() {
        Block* block = getBlock();
        block->cursor = sub(block->end, m_capacity);
    }

    // rewinds cursors in all threads
    // not thread safe with allocate()
    void resetGlobal() {
        for (Block* block = m_blocks.load(); block; block = block->next) {
            block->cursor = sub(block->end, m_capacity);
        }
    }

    size_t numBlocks() {
        size_t count = 0;
        for (Block* block = m_blocks.load(); block; block = block->next) {
            count += 1;
        }
        return count;
    }

    size_t totalSize() {
        size_t size = 0;
        for (Block* block = m_blocks.load(); block; block = block->next) {
            size += sub(block->cursor, sub(block->end, m_capacity));
        }
        return size;
    }

    size_t smallestFreeSpace() {
        size_t size = std::numeric_limits<size_t>::max();
        for (Block* block = m_blocks.load(); block; block = block->next) {
            minimize(size, sub(block->end, block->cursor));
        }
        return size;
    }

   private:
    const size_t m_capacity;
    std::atomic<Block*> m_blocks;

    Block*& getBlock() {
        thread_local Block* block = nullptr;
        return block;
    }
};

Arena g_arena(size_t(1) << 30);

template <class T>
class ArenaAllocator {
   public:
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T value_type;

    ArenaAllocator() {}
    ArenaAllocator(const ArenaAllocator&) {}

    pointer allocate(size_type n, const void* = 0) { return (T*)g_arena.allocate(n * sizeof(T)); }

    void deallocate(void* p, size_type) {}

    pointer address(reference x) const { return &x; }
    const_pointer address(const_reference x) const { return &x; }
    ArenaAllocator<T>& operator=(const ArenaAllocator&) { return *this; }
    void construct(pointer p, const T& val) { new ((T*)p) T(val); }
    void destroy(pointer p) { p->~T(); }

    size_type max_size() const { return size_t(-1); }

    template <class U>
    struct rebind {
        typedef ArenaAllocator<U> other;
    };

    template <class U>
    ArenaAllocator(const ArenaAllocator<U>&) {}

    template <class U>
    ArenaAllocator& operator=(const ArenaAllocator<U>&) {
        return *this;
    }
};
