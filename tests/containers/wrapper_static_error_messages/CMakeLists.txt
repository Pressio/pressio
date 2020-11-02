
# the following macro creates and adds an individual SERIAL test that
# we use for testing that the build of an exe fails.
# The way this works is that the test itself is excluded from
# all the list of tests built regularly and we make
# ctest itself try to build the exe.
# When we run ctest, if the build fails, it means the test succeeds.

macro(add_serial_test_to_check_build_failure
    TESTNAME SRC GREPSTRING1 GREPSTRING2 PRIVATEDEF1)
  # set name of the executable
  set(testNameFinal containers_${TESTNAME})
  add_executable(${testNameFinal} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)

  # Avoid building these targets normally
  set_target_properties(${testNameFinal} PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    EXCLUDE_FROM_DEFAULT_BUILD TRUE)

  target_compile_definitions(${testNameFinal} PRIVATE ${PRIVATEDEF1})

  # Add the tests.  These invoke "cmake --build ..." which is a
  # cross-platform way of building the given target.
  # add_test(
  #   NAME ${testNameFinal}
  #   COMMAND ${CMAKE_COMMAND} --build . --target ${testNameFinal} --config $<CONFIGURATION>
  #   WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  # here we use runtest.cmake which will
  # - attempt to buold the executable for the test
  # - dump the log to file
  # - check that the file contains specific string
  add_test(
    NAME ${testNameFinal}
    COMMAND ${CMAKE_COMMAND}
    -DtargetName=${testNameFinal}
    -DtargetLogFile=${CMAKE_CURRENT_BINARY_DIR}/log_${testNameFinal}.out
    -DstringToGrepOne=${GREPSTRING1}
    -DstringToGrepTwo=${GREPSTRING2}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  # Expect these tests to fail (i.e. cmake --build should return a non-zero value)
  #set_tests_properties(${testNameFinal} PROPERTIES WILL_FAIL TRUE)

  # taken from: https://stackoverflow.com/questions/30155619/expected-build-failure-tests-in-cmake
endmacro()
#=====================================================================

##################
### for vector ###
##################

set(main main_for_vector)

set(errMgs 
  "You cannot wrap an Eigen dense matrix as a pressio::containers::Vector<>" 
  "You cannot wrap an Eigen sparse matrix as a pressio::containers::Vector<>" 
  "You cannot wrap an Epetra multivector as a pressio::containers::Vector<>"
  "You cannot wrap a Tpetra multivector as a pressio::containers::Vector<>"
  "You cannot wrap a Tpetra block multivector as a pressio::containers::Vector<>"
  "You cannot wrap a Kokkos 2d view as a pressio::containers::Vector<>"
  )

list(GET errMgs 0 err)
add_serial_test_to_check_build_failure(
  try_wrap_eig_dmat_with_vector ${main} ${err} " " DO_EIGEN_DENSE)

list(GET errMgs 1 err)
add_serial_test_to_check_build_failure(
  try_wrap_eig_spmat_with_vector ${main} ${err} " " DO_EIGEN_SPARSE)

if(PRESSIO_ENABLE_TPL_TRILINOS)
  list(GET errMgs 2 err)
  add_serial_test_to_check_build_failure(
    try_wrap_epetra_mv_with_vector ${main} ${err} " " DO_EPETRA_MV)

  list(GET errMgs 3 err)
  add_serial_test_to_check_build_failure(
    try_wrap_tpetra_mv_with_vector ${main} ${err} " " DO_TPETRA_MV)

  list(GET errMgs 4 err)
  add_serial_test_to_check_build_failure(
    try_wrap_tpetra_block_mv_with_vector ${main} ${err} " " DO_TPETRA_BLOCK_MV)
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  list(GET errMgs 5 err)
  add_serial_test_to_check_build_failure(
    try_wrap_Kokkos_dm_with_vector ${main} ${err} " " DO_KOKKOS_DM)
endif()


########################
### for multi vector ###
########################

set(main main_for_multivector)

set(errMgs 
  "You cannot wrap an Eigen vector as a pressio::containers::MultiVector<>" 
  "You cannot wrap an Eigen sparse matrix as a pressio::containers::MultiVector<>" 
  "You cannot wrap an Epetra vector as a pressio::containers::MultiVector<>"
  "You cannot wrap a Tpetra vector as a pressio::containers::MultiVector<>"
  "You cannot wrap a Tpetra block vector as a pressio::containers::MultiVector<>"
  "You cannot wrap a Kokkos 1d view as a pressio::containers::MultiVector<>"
  )

list(GET errMgs 0 err)
add_serial_test_to_check_build_failure(
  try_wrap_eig_vec_with_multivector ${main} ${err} " " DO_EIGEN_DENSE)

list(GET errMgs 1 err)
add_serial_test_to_check_build_failure(
  try_wrap_eig_spmat_with_multivector ${main} ${err} " " DO_EIGEN_SPARSE)

if(PRESSIO_ENABLE_TPL_TRILINOS)
  list(GET errMgs 2 err)
  add_serial_test_to_check_build_failure(
    try_wrap_epetra_vec_with_multivector ${main} ${err} " " DO_EPETRA_V)

  list(GET errMgs 3 err)
  add_serial_test_to_check_build_failure(
    try_wrap_tpetra_vec_with_multivector ${main} ${err} " " DO_TPETRA_V)

  list(GET errMgs 4 err)
  add_serial_test_to_check_build_failure(
    try_wrap_tpetra_block_vec_with_multivector ${main} ${err} " " DO_TPETRA_BLOCK_V)
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  list(GET errMgs 5 err)
  add_serial_test_to_check_build_failure(
    try_wrap_Kokkos_vec_with_multivector ${main} ${err} " " DO_KOKKOS_V)
endif()




########################
### for matrix ###
########################

set(main main_for_matrix)

set(errMgs 
  "You cannot wrap an Eigen vector as a pressio::containers::DenseMatrix<>" 
  "You cannot wrap an Eigen sparse matrix as a pressio::containers::DenseMatrix<>" 
  "You cannot wrap an Epetra vector as a pressio::containers::DenseMatrix<>"
  "You cannot wrap a Tpetra vector as a pressio::containers::DenseMatrix<>"
  "You cannot wrap a Tpetra block vector as a pressio::containers::DenseMatrix<>"
  "You cannot wrap a Kokkos 1d view as a pressio::containers::DenseMatrix<>"
  )

list(GET errMgs 0 err)
add_serial_test_to_check_build_failure(
  try_wrap_eig_vec_with_matrix ${main} ${err} " " DO_EIGEN_DENSE)

list(GET errMgs 1 err)
add_serial_test_to_check_build_failure(
  try_wrap_eig_spmat_with_matrix ${main} ${err} " " DO_EIGEN_SPARSE)

if(PRESSIO_ENABLE_TPL_TRILINOS)
  list(GET errMgs 2 err)
  add_serial_test_to_check_build_failure(
    try_wrap_epetra_vec_with_matrix ${main} ${err} " " DO_EPETRA_V)

  list(GET errMgs 3 err)
  add_serial_test_to_check_build_failure(
    try_wrap_tpetra_vec_with_matrix ${main} ${err} " " DO_TPETRA_V)

  list(GET errMgs 4 err)
  add_serial_test_to_check_build_failure(
    try_wrap_tpetra_block_vec_with_matrix ${main} ${err} " " DO_TPETRA_BLOCK_V)
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  list(GET errMgs 5 err)
  add_serial_test_to_check_build_failure(
    try_wrap_Kokkos_vec_with_matrix ${main} ${err} " " DO_KOKKOS_V)
endif()
