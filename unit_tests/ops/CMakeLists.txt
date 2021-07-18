
# to include google test fixtures
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../containers/fixtures)

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_native.cc)
  add_serial_utest(ops_eigen_native_objects "${SRC1}")

  #set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_expressions.cc)
  #add_serial_utest(ops_eigen_expressions "${SRC1}")
endif()

#if(PRESSIO_ENABLE_TPL_KOKKOS)
#add_subdirectory(kokkos)
#endif()
#if(PRESSIO_ENABLE_TPL_MPI AND PRESSIO_ENABLE_TPL_TRILINOS)
#add_subdirectory(trilinos)
#endif()
