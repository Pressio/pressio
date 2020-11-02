
# to include google test fixtures
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/fixtures)

add_subdirectory(vector)
add_subdirectory(dense_matrix)
add_subdirectory(multi_vector)
add_subdirectory(expressions)
add_subdirectory(various)
