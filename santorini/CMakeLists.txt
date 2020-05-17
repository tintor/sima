project(santorini)

add_executable(test_santorini
	../test.cc
	test_santorini.cc)

add_executable(santorini santorini.cc)
add_executable(old old.cc)

target_link_libraries(test_santorini ${LIBRARIES} core)
target_link_libraries(old       ${LIBRARIES} core)
target_link_libraries(santorini ${LIBRARIES} core view)
