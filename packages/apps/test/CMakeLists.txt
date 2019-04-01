
#------------------------------------------------------------------------
# to include the lumped headers to have access to package classes
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)

#include the utilities in utils
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/utils)

ADD_SUBDIRECTORIES(burgers1d)
ADD_SUBDIRECTORIES(steady_linear_adv_diff_1d)
ADD_SUBDIRECTORIES(steady_linear_adv_diff_2d)
