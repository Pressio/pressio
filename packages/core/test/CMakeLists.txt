

#include the lumped headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lumped_include_headers)

ADD_SUBDIRECTORIES(
  kokkos_works
  #vector
  #eigen1
)
