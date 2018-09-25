
#include the lumped headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)

ADD_SUBDIRECTORIES(
   linear_direct
   linear_iterative
   # leastsquares
)
