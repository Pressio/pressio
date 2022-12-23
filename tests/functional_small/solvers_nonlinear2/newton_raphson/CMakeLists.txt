include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../problems)
if (PRESSIO_ENABLE_TPL_EIGEN)
add_serial_exe_and_test(newton_raphson_eigen ${ROOTNAME} eigen.cc "PASSED")
#add_serial_exe_and_test(eigen_ch_w_bad_jac_w_backtrack_strictly_decreasing_objective ${ROOTNAME} eigen_ch_w_bad_jac_w_backtrack_strictly_decreasing_objective.cc "PASSED")
endif()

#add_serial_exe_and_test(newton_raphson_custom_types ${ROOTNAME} custom_types.cc "PASSED")
