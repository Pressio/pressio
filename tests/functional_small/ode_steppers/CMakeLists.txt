
add_serial_utest(ode_steppers_concepts ${CMAKE_CURRENT_SOURCE_DIR}/concepts.cc)

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(FILENAME ode_forward_euler_simple_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_forward_euler_mass_matrix_manual_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_rk4_simple_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_ssprk3_simple_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_ab2_simple_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_all_schemes_mass_matrix_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})
endif()


set(FILENAME ode_forward_euler_simple_correctness_custom_types)
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

set(FILENAME ode_rk4_simple_correctness_custom_types)
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

set(FILENAME ode_ab2_simple_correctness_custom_types)
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})
