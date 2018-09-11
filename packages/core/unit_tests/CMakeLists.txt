
# this function creates and adds an individual test executable
# this is needed for cases where we need to specify the number
# of MPI processes, for example to do testing for distributed objects
#
#  if more than 2 arguments, the third one HAS to be the num of ranks
#
#
function(add_unit_test test_name mainfile )

  if(ARGC EQUAL 2)
    SET(NUMMPI 1)
  else(ARGC GREATER 2)
    SET(NUMMPI ${ARGV2})
  endif()
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    utest_${test_name}
    SOURCES ${test_name}.cc ${mainfile}.cc
    COMM serial mpi
    NUM_MPI_PROCS ${NUMMPI}
    PASS_REGULAR_EXPRESSION "PASSED"
    FAIL_REGULAR_EXPRESSION "FAILED"
  )
endfunction()


#include the lumped headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lumped_include_headers)

ADD_SUBDIRECTORY(serial)
ADD_SUBDIRECTORY(mpi)
# ADD_SUBDIRECTORY(kokkos)
