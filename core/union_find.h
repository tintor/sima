#pragma once

class UnionFind {
   public:
    UnionFind() : _parent(this), _rank(0) {}

    void merge(UnionFind& b) {
        UnionFind* pa = find();
        UnionFind* pb = b.find();

        if (pa->_rank < pb->_rank) {
            pa->_parent = pb;
        } else if (pa->_rank > pb->_rank) {
            pb->_parent = pa;
        } else {
            pa->_rank += 1;
            pb->_parent = pa;
        }
    }

    UnionFind* find() {
        // Path halving faster than path compression (from Wikipedia)
        UnionFind* e = this;
        while (e->_parent != e) {
            e->_parent = e->_parent->_parent;
            e = e->_parent;
        }
        return e;
    }

   private:
    UnionFind* _parent;
    uint _rank;
};
