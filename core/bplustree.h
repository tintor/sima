#pragma once
#include <core/exception.h>
#include <core/format.h>
#include <core/range.h>
#include <core/std.h>

template <typename T>
class FreeList {
   public:
    ~FreeList() { clear(); }

    T* get() {
        if (!_data) return new T;
        auto p = reinterpret_cast<T**>(_data);
        T* e = _data;
        _data = *p;
        return e;
    }

    void put(T* e) {
        auto p = reinterpret_cast<T**>(e);
        *p = _data;
        _data = e;
    }

    void clear() {
        while (_data) delete get();
    }

   private:
    static_assert(sizeof(T) >= sizeof(T*));
    T* _data = nullptr;
};

// TODO bool uniqueness
// TODO bool reverse_iterable
// TODO remove operation and node merging
// TODO binary search in leaves and branches
// TODO interleave keys and pointers or not? (same for leaves)
template <typename K, typename V, uint NL, uint NB = NL>
class bplustree_map {
   private:
    static_assert(NL >= 2);
    static_assert(NB >= 2);
    constexpr static uint MinNL = std::min(2u, NL - (NL / 2));
    constexpr static uint MinNB = std::min(2u, NB - (NB / 2));
    struct leaf;
    struct branch;
    using KV = pair<const K, V>;

   public:
    class iterator {
       public:
        iterator(leaf* leaf, uint index) : _leaf(leaf), _index(index) {}
        pair<const K, V>& operator*() { return _leaf->pairs[_index]; }
        pair<const K, V>* operator->() { return &_leaf->pairs[_index]; }
        bool operator==(iterator e) const { return _leaf == e._leaf && _index == e._index; }
        bool operator!=(iterator e) const { return _leaf != e._leaf || _index != e._index; }
        void operator++() {
            if (++_index == _leaf->size) {
                _leaf = _leaf->next;
                _index = 0;
            }
        }

       private:
        leaf* _leaf;
        uint _index;
    };

    class const_iterator {
       public:
        const_iterator(const leaf* leaf, uint index) : _leaf(leaf), _index(index) {}
        const pair<const K, V>& operator*() const { return _leaf->pairs[_index]; }
        bool operator==(const_iterator e) const { return _leaf == e._leaf && _index == e._index; }
        bool operator!=(const_iterator e) const { return _leaf != e._leaf || _index != e._index; }
        void operator++() {
            if (++_index == _leaf->size) {
                _leaf = _leaf->next;
                _index = 0;
            }
        }

       private:
        const leaf* _leaf;
        uint _index;
    };

    bool contains(K key) const { return find(key) != end(); }

    /*	void remove(K key) {
                    auto it = find(key);
                    if (it != end())
                            remove(it);
            }

            iterator remove(iterator it) {
                    // TODO scan from root to find branch pointers and positions within branch
                    leaf* l = it->_leaf;
                    uint p = it->_index;
                    // TODO remove element from leaf node
                    if (l->size >= N / 2)
                            return (p == l->size) ? iterator(l->next, 0) : iterator(l, p);
                    // TODO look at prev and next nodes
                    THROW(not_implemented);
            }*/

    const V& at(K key) const { return find(key)->second; }
    V& operator[](K key) { return find_or_insert(key)->second; }
    const V& operator[](K key) const { return find(key)->second; }

    iterator find(K key) {
        if (_root == nullptr) return end();
        branch* b = _root;
        for (uint i : range(_branch_levels)) b = b->find(key);
        leaf* l = reinterpret_cast<leaf*>(b);
        for (uint i : range(l->size))
            if (key == l->pairs[i].first) return iterator(l, i);
        return end();
    }

    const_iterator find(K key) const {
        if (_root == nullptr) return end();
        branch* b = _root;
        for (uint i : range(_branch_levels)) b = b->find(key);
        leaf* l = reinterpret_cast<leaf*>(b);
        for (uint i : range(l->size))
            if (key == l->pairs[i].first) return const_iterator(l, i);
        return end();
    }

    // 590ms @ 24
    iterator find_or_insert(K key) {
        if (_root == nullptr) {
            leaf* l = _min_leaf = _unused_leaves.get();
            _root = reinterpret_cast<branch*>(l);
            l->size = 1;
            l->init(0, key);
            l->next = nullptr;
            _size = 1;
            return iterator(l, 0);
        }
        branch* b = _root;
        array<pair<branch*, uint>, 30> stack;
        for (uint i : range(_branch_levels)) {
            uint pos = b->findi(key);
            stack[i] = {b, pos};
            b = b->pointers[pos];
        }
        leaf* l = reinterpret_cast<leaf*>(b);
        uint index = l->size;
        for (uint i : range(l->size)) {
            if (key < l->pairs[i].first) {
                index = i;
                break;
            }
            if (key == l->pairs[i].first) return iterator(l, i);
        }
        return insert_into_leaf(l, index, key, stack);
    }

    uint size() const { return _size; }

    iterator begin() { return iterator(_min_leaf, 0); }
    iterator end() { return iterator(nullptr, 0); }
    const_iterator begin() const { return const_iterator(_min_leaf, 0); }
    const_iterator end() const { return const_iterator(nullptr, 0); }

    void assert_invariant() {
        if (_size == 0 || _root == nullptr || _min_leaf == nullptr) {
            ASSERT_ALWAYS(_size == 0);
            ASSERT_ALWAYS(_min_leaf == nullptr);
            ASSERT_ALWAYS(_root == nullptr);
            ASSERT_ALWAYS(_branch_levels == 0);
        }

        // _size and no loops
        uint sum = 0;
        for (leaf* a = _min_leaf; a; a = a->next) {
            ASSERT_ALWAYS(a->size > 0);
            sum += a->size;
            ASSERT_ALWAYS(sum <= _size);
        }
        ASSERT_ALWAYS(_size == sum);

        if (_root) ASSERT_ALWAYS(assert_branch(_root, 0, nullopt, nullopt, _min_leaf) == nullptr);
    }

    void debug() { debug(_root, 0); }

   private:
    void debug(branch* b, uint level) {
        for (auto i : range(level)) print("    ");
        print("0x%p ", b);
        if (level == _branch_levels) {
            auto l = reinterpret_cast<leaf*>(b);
            for (auto i : range(l->size)) print("%s:%s ", l->pairs[i].first, l->pairs[i].second);
            print("\n");
            return;
        }
        for (auto i : range(b->size - 1)) print("%s ", b->keys[i]);
        print("\n");
        for (auto i : range(b->size)) debug(b->pointers[i], level + 1);
    }

    void insert_into_branch(branch* b, uint level, K key, const array<pair<branch*, uint>, 30>& stack) {
        if (level == 0) {
            // tree grows one level
            _branch_levels += 1;
            auto a = _root;
            _root = _unused_branches.get();
            _root->size = 2;
            _root->keys[0] = key;
            _root->pointers[0] = a;
            _root->pointers[1] = b;
            return;
        }

        branch* p = stack[level - 1].first;
        uint pos = stack[level - 1].second + 1;
        ASSERT_ALWAYS(0 < pos && pos <= p->size);

        if (p->size < NB) {
            move(p->pointers + pos + 1, p->pointers + pos, p->pointers + p->size);
            p->pointers[pos] = b;
            move(p->keys + pos, p->keys + pos - 1, p->keys + p->size - 1);
            p->keys[pos - 1] = key;
            p->size += 1;
            return;
        }

        // TODO make it faster by avoiding this extra copy (or using a smaller index array)
        uint w = 0;
        array<K, NB + 1> keys;
        array<branch*, NB + 1> pointers;
        keys[w] = K();
        pointers[w++] = p->pointers[0];
        for (auto i : range(1u, p->size)) {
            if (i == pos) {
                keys[w] = key;
                pointers[w++] = b;
            }
            keys[w] = p->keys[i - 1];
            pointers[w++] = p->pointers[i];
        }
        if (p->size == pos) {
            keys[w] = key;
            pointers[w++] = b;
        }

        // need to split the branch node
        branch* q = _unused_branches.get();
        p->size = NB + 1 - MinNB;
        q->size = MinNB;

        copy(p->pointers, pointers.data(), pointers.data() + p->size);
        copy(q->pointers, pointers.data() + p->size, pointers.data() + NB + 1);
        copy(p->keys, keys.data() + 1, keys.data() + p->size);
        copy(q->keys, keys.data() + 1 + p->size, keys.data() + NB + 1);
        K extra_key = keys[p->size];

        /*if (pos < MinN) {
                // new pointer goes to P
                copy(q->pointers, p->pointers + MinN - 1, p->pointers + N);
                move(p->pointers + pos + 1, p->pointers + pos, p->pointers + MinN - 1);
                p->pointers[pos] = b;

                copy(q->keys, p->keys + MinN, p->keys + N);
                extra_key = p->keys[MinN - 1];
                p->keys[pos - 1] = key;
        } else {
                // new pointer goes to Q
                copy(q->pointers + pos - MinN + 1, p->pointers + pos + 1, p->pointers + N);
                q->pointers[pos - MinN] = b;
                copy(q->pointers, p->pointers + MinN, p->pointers + pos + 1);

                q->keys[pos - MinN - 1] = key;
                extra_key = ...
        }*/

        insert_into_branch(q, level - 1, extra_key, stack);
    }

    template <typename T>
    static void move(T* dest, const T* src_begin, const T* src_end) {
        memmove(dest, src_begin, (src_end - src_begin) * sizeof(T));
    }
    template <typename T>
    static void copy(T* dest, const T* src_begin, const T* src_end) {
        memcpy(dest, src_begin, (src_end - src_begin) * sizeof(T));
    }

    iterator insert_into_leaf(leaf* a, uint pos, K key, const array<pair<branch*, uint>, 30>& stack) {
        _size += 1;
        if (a->size < NL) {
            move(a->pairs + pos + 1, a->pairs + pos, a->pairs + a->size);
            a->init(pos, key);
            a->size += 1;
            return iterator(a, pos);
        }

        leaf* b = _unused_leaves.get();
        b->next = a->next;
        a->next = b;
        a->size = NL + 1 - MinNL;
        b->size = MinNL;
        if (pos < a->size) {
            // new key goes to A
            copy(b->pairs, a->pairs + NL - MinNL, a->pairs + NL);
            move(a->pairs + pos + 1, a->pairs + pos, a->pairs + NL - MinNL);
            a->init(pos, key);
            insert_into_branch(reinterpret_cast<branch*>(b), _branch_levels, b->pairs[0].first, stack);
            return iterator(a, pos);
        } else {
            // new key goes to B
            copy(b->pairs + pos + MinNL - NL, a->pairs + pos, a->pairs + NL);
            b->init(pos + MinNL - NL - 1, key);
            copy(b->pairs, a->pairs + NL + 1 - MinNL, a->pairs + pos);
            insert_into_branch(reinterpret_cast<branch*>(b), _branch_levels, b->pairs[0].first, stack);
            return iterator(b, pos - a->size);
        }
    }

    leaf* assert_branch(branch* b, uint level, optional<K> lower, optional<K> upper, leaf* next_leaf) {
        ASSERT_ALWAYS(b != nullptr);
        if (level == _branch_levels) {
            auto l = reinterpret_cast<leaf*>(b);
            ASSERT_ALWAYS(l == next_leaf, "l=%p next=%p", l, next_leaf);

            // verify leaf size
            const uint min_size = (b == _root) ? 1 : MinNL;
            ASSERT_ALWAYS(min_size <= l->size && l->size <= NL);

            // keys in leaves are in proper range
            if (lower) ASSERT_ALWAYS(*lower <= l->pairs[0].first);
            if (upper) ASSERT_ALWAYS(l->pairs[l->size - 1].first < *upper);

            // keys inside leaves are in order
            for (auto i : range(1u, l->size)) ASSERT_ALWAYS(l->pairs[i - 1].first < l->pairs[i].first);
            // keys between leaves are in order
            if (l->next) ASSERT_ALWAYS(l->pairs[l->size - 1].first < l->next->pairs[0].first);

            return l->next;
        }

        // verify node size
        const uint min_size = (b == _root) ? 2 : MinNB;
        ASSERT_ALWAYS(min_size <= b->size && b->size <= NB, "branch size %s not in [%s %s]", b->size, min_size, NB);

        // verify keys are in order
        for (auto i : range(1u, b->size - 1)) ASSERT_ALWAYS(b->keys[i - 1] < b->keys[i]);
        if (lower) ASSERT_ALWAYS(*lower <= b->keys[0]);
        if (upper) ASSERT_ALWAYS(b->keys[b->size - 2] < *upper);

        // recursive descent
        next_leaf = assert_branch(b->pointers[0], level + 1, lower, b->keys[0], next_leaf);
        for (auto i : range(1u, b->size - 1))
            next_leaf = assert_branch(b->pointers[i], level + 1, b->keys[i - 1], b->keys[i], next_leaf);
        return assert_branch(b->pointers[b->size - 1], level + 1, b->keys[b->size - 2], upper, next_leaf);
    }

   private:
    struct branch {
        uint size;  // number of pointers
        K keys[NB - 1];
        branch* pointers[NB];

        branch* find(K key) const {
            for (auto i : range(0u, size - 1))
                if (key < keys[i]) return pointers[i];
            return pointers[size - 1];
        }
        uint findi(K key) const {
            for (auto i : range(0u, size - 1))
                if (key < keys[i]) return i;
            return size - 1;
        }
    };

    struct leaf {
        uint size;
        KV pairs[NL];
        leaf* next;

        void init(uint i, K key) {
            reinterpret_cast<pair<K, V>&>(pairs[i]).first = key;
            pairs[i].second = V();
        }
    };

    uint _size = 0;
    uint _branch_levels = 0;
    branch* _root = nullptr;
    leaf* _min_leaf = nullptr;

    FreeList<branch> _unused_branches;
    FreeList<leaf> _unused_leaves;
};

template <typename T>
T quick_select(span<T> m, size_t kth) {
    ASSERT_ALWAYS(kth < m.size());
    T pivot = median(m[0], m[m.size() / 2], m[m.size() - 1]);
    auto it = std::partition(m.begin(), m.end(), [pivot](const T& e) { return e < pivot; });
    size_t index = std::distance(m.begin(), it);
    if (kth < index) return quick_select(m.subspan(0, index), kth);
    if (kth > index) return quick_select(m.subspan(index + 1), kth - index - 1);
    return *it;
}

template <typename T>
T percentile(span<T> m, double p) {
    ASSERT_ALWAYS(0 <= p && p <= 1);
    size_t rank = std::max(0, floor(p * m.size()));
    return quick_select(m, std::min(m.size() - 1, rank));
}

template <typename K, uint NL, uint NB = NL>
class ranked_set {
   private:
    static_assert(NL >= 2);
    static_assert(NB >= 2);
    constexpr static uint MinNL = std::min(2u, NL - (NL / 2));
    constexpr static uint MinNB = std::min(2u, NB - (NB / 2));
    struct leaf;
    struct branch;

   public:
    uint rank(K key) {
        branch* b = _root;
        uint count = 0;
        for (uint i : range(_branch_levels)) {
            uint pos = b->find(key);
            if (pos > 0)
                for (uint j : range(pos - 1)) count += b->sizes[j];
            b = b->pointers[pos];
        }
        leaf* l = reinterpret_cast<leaf*>(b);
        for (uint i : range(l->size))
            if (key <= l->keys[i]) return count + i;
        return count + l->size - 1;
    }

    K findByRank(uint count) {
        ASSERT_ALWAYS(count < _size);
        branch* b = _root;
        for (uint i : range(_branch_levels)) {
            for (uint j : range(b->size)) {
                if (count < b->sizes[j]) {
                    b = b->pointers[j];
                    break;
                }
                count -= b->sizes[j];
            }
        }
        return reinterpret_cast<leaf*>(b)->keys[count];
    }

    K percentile(double percentile) {
        uint rank = std::min(_size - 1, uint(floor(percentile * _size)));
        return findByRank(rank);
    }

    void insert(K key) {
        if (_root == nullptr) {
            leaf* l = new leaf;
            _root = reinterpret_cast<branch*>(l);
            l->size = 1;
            l->keys[0] = key;
            _size = 1;

            debug();
            verify(_root, _size, 0);
            return;
        }
        branch* b = _root;
        array<pair<branch*, uint>, 30> stack;
        for (uint i : range(_branch_levels)) {
            uint pos = b->find(key);
            stack[i] = {b, pos};
            b->sizes[pos] += 1;
            b = b->pointers[pos];
        }
        leaf* l = reinterpret_cast<leaf*>(b);
        uint index = l->size;
        for (uint i : range(l->size))
            if (key <= l->keys[i]) {
                index = i;
                break;
            }
        insert_into_leaf(l, index, key, stack);
        debug();
        verify(_root, _size, 0);
    }

    void debug() { debug(_root, 0); }

   private:
    void debug(branch* b, uint level) {
        for (auto i : range(level)) print("    ");
        print("0x%p ", b);
        if (level == _branch_levels) {
            auto l = reinterpret_cast<leaf*>(b);
            for (auto i : range(l->size)) print("%s ", l->keys[i]);
            print("\n");
            return;
        }
        print("*[%s] ", b->sizes[0]);
        for (auto i : range(b->size - 1)) print("%s[%s] ", b->keys[i], b->sizes[i + 1]);
        print("\n");
        for (auto i : range(b->size)) debug(b->pointers[i], level + 1);
    }

    void verify(branch* b, uint size, uint level) {
        if (level == _branch_levels) {
            ASSERT_ALWAYS(size == reinterpret_cast<leaf*>(b)->size, "%s == %s", size, reinterpret_cast<leaf*>(b)->size);
            return;
        }
        uint sum = 0;
        for (auto i : range(b->size)) {
            verify(b->pointers[i], b->sizes[i], level + 1);
            sum += b->sizes[i];
        }
        ASSERT_ALWAYS(size == sum);
    }

    void insert_into_branch(uint a_size, uint b_size, branch* b, uint level, K key,
                            const array<pair<branch*, uint>, 30>& stack) {
        if (level == 0) {
            // tree grows one level
            _branch_levels += 1;
            auto a = _root;
            _root = new branch;
            _root->size = 2;
            _root->keys[0] = key;
            _root->pointers[0] = a;
            _root->pointers[1] = b;
            _root->sizes[0] = a_size;
            _root->sizes[1] = b_size;
            return;
        }

        branch* p = stack[level - 1].first;
        uint pos = stack[level - 1].second + 1;
        ASSERT_ALWAYS(0 < pos && pos <= p->size);

        // print("a_size %s, b_size %s, pos %s\n", a_size, b_size, pos);
        if (p->size < NB) {
            move(p->pointers + pos + 1, p->pointers + pos, p->pointers + p->size);
            move(p->sizes + pos + 1, p->sizes + pos, p->sizes + p->size);
            p->pointers[pos] = b;
            move(p->keys + pos, p->keys + pos - 1, p->keys + p->size - 1);
            p->keys[pos - 1] = key;
            p->sizes[pos - 1] = a_size;
            p->sizes[pos] = b_size;
            p->size += 1;
            return;
        }

        // TODO make it faster by avoiding this extra copy (or using a smaller index array)
        uint w = 0;
        array<K, NB + 1> keys;
        array<branch*, NB + 1> pointers;
        array<uint, NB + 1> sizes;
        keys[w] = K();
        pointers[w] = p->pointers[0];
        sizes[w++] = p->sizes[0];
        for (auto i : range(1u, p->size)) {
            if (i == pos) {
                keys[w] = key;
                pointers[w] = b;
                sizes[w++] = b_size;
            }
            keys[w] = p->keys[i - 1];
            pointers[w] = p->pointers[i];
            sizes[w++] = (i == pos - 1) ? a_size : p->sizes[i];
        }
        if (p->size == pos) {
            keys[w] = key;
            pointers[w] = b;
            sizes[w++] = b_size;
        }

        // need to split the branch node
        branch* q = new branch;
        p->size = NB + 1 - MinNB;
        q->size = MinNB;

        copy(p->pointers, pointers.data(), pointers.data() + p->size);
        copy(q->pointers, pointers.data() + p->size, pointers.data() + NB + 1);
        copy(p->keys, keys.data() + 1, keys.data() + p->size);
        copy(q->keys, keys.data() + 1 + p->size, keys.data() + NB + 1);
        copy(p->sizes, sizes.data(), sizes.data() + p->size);
        copy(q->sizes, sizes.data() + p->size, sizes.data() + NB + 1);
        K extra_key = keys[p->size];

        uint p_size = 0;
        for (uint i : range(p->size)) p_size += p->sizes[i];
        uint q_size = 0;
        for (uint i : range(q->size)) q_size += q->sizes[i];

        insert_into_branch(p_size, q_size, q, level - 1, extra_key, stack);
    }

    template <typename T>
    static void move(T* dest, const T* src_begin, const T* src_end) {
        memmove(dest, src_begin, (src_end - src_begin) * sizeof(T));
    }
    template <typename T>
    static void copy(T* dest, const T* src_begin, const T* src_end) {
        memcpy(dest, src_begin, (src_end - src_begin) * sizeof(T));
    }

    void insert_into_leaf(leaf* a, uint pos, K key, const array<pair<branch*, uint>, 30>& stack) {
        _size += 1;
        if (a->size < NL) {
            move(a->keys + pos + 1, a->keys + pos, a->keys + a->size);
            a->keys[pos] = key;
            a->size += 1;
            return;
        }

        leaf* b = new leaf;
        a->size = NL + 1 - MinNL;
        b->size = MinNL;
        if (pos < a->size) {
            // new key goes to A
            copy(b->keys, a->keys + NL - MinNL, a->keys + NL);
            move(a->keys + pos + 1, a->keys + pos, a->keys + NL - MinNL);
            a->keys[pos] = key;
            insert_into_branch(a->size, b->size, reinterpret_cast<branch*>(b), _branch_levels, b->keys[0], stack);
        } else {
            // new key goes to B
            copy(b->keys + pos + MinNL - NL, a->keys + pos, a->keys + NL);
            b->keys[pos + MinNL - NL - 1] = key;
            copy(b->keys, a->keys + NL + 1 - MinNL, a->keys + pos);
            insert_into_branch(a->size, b->size, reinterpret_cast<branch*>(b), _branch_levels, b->keys[0], stack);
        }
    }

   private:
    struct branch {
        uint size;  // number of pointers, TODO get rid of this as it is tored in parent
        K keys[NB - 1];
        branch* pointers[NB];
        uint sizes[NB];  // number of leaf keys in every child node

        uint find(K key) const {
            for (auto i : range(0u, size - 1))
                if (key < keys[i]) return i;
            return size - 1;
        }
    };

    struct leaf {
        uint size;  // TODO get rid of this as it is stored in parent
        K keys[NL];
    };

    uint _size = 0;
    uint _branch_levels = 0;
    branch* _root = nullptr;
};
