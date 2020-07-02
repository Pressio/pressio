
configure_file(../../basis/svd_basis_ncell20_t010_dt001_implicit_euler.txt basis.txt COPYONLY)
find_package(Torch REQUIRED)
#find_package(OpenMP REQUIRED)

set(testNameFinal rom_burgers1d_eigen_lspg_euler_gn_normal_eqns_residual_api_user_defined_mapping)
add_executable(${testNameFinal} main.cc)
target_link_libraries(${testNameFinal} "${TORCH_LIBRARIES}")
set_property(TARGET ${testNameFinal} PROPERTY CXX_STANDARD 14)
target_include_directories(${testNameFinal} SYSTEM BEFORE PUBLIC "/scratch/pressio_install/pressio_repos/libtorch/include")

add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
  set_tests_properties(
    ${testNameFinal}
    PROPERTIES
    PASS_REGULAR_EXPRESSION "PASSED"
    FAIL_REGULAR_EXPRESSION "FAILED"
  )

