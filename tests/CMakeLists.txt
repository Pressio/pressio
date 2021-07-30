
include_directories(${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/include)

if(PRESSIO_ENABLE_FUNCTIONAL_SMALL_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_small)
endif()

if(PRESSIO_ENABLE_FUNCTIONAL_MEDIUM_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_medium)
endif()
