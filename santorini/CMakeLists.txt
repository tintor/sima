project(santorini)

add_executable(santorini santorini.cc)
add_executable(santorini2 santorini2.cc)

target_link_libraries(santorini   ${LIBRARIES} core)
target_link_libraries(santorini2  ${LIBRARIES} core view)
