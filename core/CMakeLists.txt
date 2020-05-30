add_library(core
	callstack.cc
	diff.cc
	file.cc
    model.cc
    optimizer.cc
	util.cc
	timestamp.cc)

add_executable(core_test
	../test.cc
	test_array_deque.cc
	test_array_set.cc
 	test_auto.cc
	test_bits.cc
	test_diff.cc
	test_dynamic_array.cc
	test_format.cc
	test_map.cc
    test_property.cc
	test_range.cc
	test_sparse_matrix.cc
	test_tensor.cc
	test_util.cc
	test_string_util.cc
	test_zip.cc
    test_ipp.cc
	test_bplustree.cc)

target_link_libraries(core_test core ${LIBRARIES})
