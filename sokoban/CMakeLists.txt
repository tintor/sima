project(sokoban)

add_executable(soko corrals.cc level.cc util.cc main.cc)
target_link_libraries(soko core ${LIBRARIES})
