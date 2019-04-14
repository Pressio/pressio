
# #include the lumped headers to access classes of solvers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)

#include the utilities in utils
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../test/utils)

ADD_SUBDIRECTORY(steady_linear_adv_diff_2d)
