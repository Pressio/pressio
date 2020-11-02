
configure_file(../../../basis/decoder_ncell100_romsize8_t30_dt_001_implicit_euler.pt basis.pt COPYONLY)

find_package(Torch REQUIRED)

set(testNameFinal rom_burgers1d_eigen_lspg_euler_gn_normal_eqns_discrete_time_api_libtorch_mapping)
add_executable(${testNameFinal} main.cc)
target_link_libraries(${testNameFinal} "${TORCH_LIBRARIES}")
set_property(TARGET ${testNameFinal} PROPERTY CXX_STANDARD 14)
target_include_directories(${testNameFinal} BEFORE PUBLIC "${TORCH_INCLUDE_DIRS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
  set_tests_properties(
    ${testNameFinal}
    PROPERTIES
    PASS_REGULAR_EXPRESSION "PASSED"
    FAIL_REGULAR_EXPRESSION "FAILED"
  )
