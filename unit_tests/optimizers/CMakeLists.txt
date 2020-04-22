
add_serial_utest(optimizers_param_list optimizers_parameter_list.cc)

if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_serial_utest(optimizers_param_list_to_rol optimizers_parameter_list_rol.cc)
endif()
