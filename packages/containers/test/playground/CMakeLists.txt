
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/build)


# # This is needed so that I can access the native cmake inside pybind11
# add_subdirectory(${PYBIND11_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/build)
# # this is what I can use to add module
# pybind11_add_module(main main.cc)

configure_file(test.py test.py)
