#pragma once
// Define vector with sizeof(vector) == 16 instead of 24
//#include <boost/container/vector.hpp>
//template <typename T, typename Allocator = std::allocator<T>>
//using vector =
//    boost::container::vector<T, Allocator,
//                             typename boost::container::vector_options<boost::container::stored_size<uint32_t>>::type>;
//static_assert(sizeof(vector<int>) == 16);

// vector implemented in array
#include <boost/container/static_vector.hpp>
using boost::container::static_vector;

// static storage for small number of elements
#include <boost/container/small_vector.hpp>
using boost::container::small_vector;
