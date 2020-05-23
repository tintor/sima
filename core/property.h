#pragma once

template<typename P, unsigned long tag>
class PropertyT {
protected:
    [[no_unique_address]]
    struct {
        operator P*() { return reinterpret_cast<P*>(this); }
        operator const P*() const { return reinterpret_cast<const P*>(this); }
        P* operator->() { return operator P*(); }
        const P* operator->() const { return operator const P*(); }
    } parent;
};

#define TProperty(N, M) [[no_unique_address]] struct N : public PropertyT<M, __LINE__>
#define Property(M) [[no_unique_address]] struct : public PropertyT<M, __LINE__>
