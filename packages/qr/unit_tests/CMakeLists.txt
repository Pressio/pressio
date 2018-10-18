
#------------------------------------------------------------------------
# this function creates and adds an individual test executable
# this is needed for cases where we need to specify the number
# of MPI processes, for example to do testing for distributed objects
#
#  if more than 2 arguments, the third one HAS to be the num of ranks
#
#
function(add_unit_test test_name mainfile )

  SET(NUMMPI ${ARGV2})

  if(${NUMMPI} EQUAL 1)
    SET(COMM1 serial mpi)
  else(${NUMMPI} GREATER 1)
    SET(COMM1 mpi)
  endif()

  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    utest_${test_name}
    SOURCES ${test_name}.cc ${mainfile}.cc 
    # COMM just means that the test will be run if rompp is built with or without MPI
    # for MPI, the number of ranks to use is set by the
    # test itself when calling this function
    COMM ${COMM1}
    NUM_MPI_PROCS ${NUMMPI}
    PASS_REGULAR_EXPRESSION "PASSED"
    FAIL_REGULAR_EXPRESSION "FAILED"
  )
endfunction()

#------------------------------------------------------------------------
# to include the lumped headers to have access to core classes
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)

#------------------------------------------------------------------------
# include unit tests subfolders
#------------------------------------------------------------------------
ADD_SUBDIRECTORY(serial)

IF(HAVE_MPI)
  ADD_SUBDIRECTORY(mpi)
ENDIF()
