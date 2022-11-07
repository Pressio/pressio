
if (PRESSIO_ENABLE_TPL_EIGEN)
add_serial_exe_and_test(newton_raphson_eigen ${ROOTNAME} eigen.cc "PASSED")
add_serial_exe_and_test(newton_raphson_eigen_backtrack ${ROOTNAME} eigen_w_basic_backtrack.cc "PASSED")
endif()

#add_serial_exe_and_test(newton_raphson_custom_types ${ROOTNAME} custom_types.cc "PASSED")
