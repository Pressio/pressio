
add_serial_unit_test(vector_serial_eigen_dynamic_class)
add_serial_unit_test(vector_serial_eigen_dynamic_expr_templates)

IF(HAVE_TRILINOS)
	add_serial_unit_test(vector_serial_teuchos_dense_class)
	add_serial_unit_test(vector_serial_teuchos_dense_expr_templates)
ENDIF()

IF(HAVE_BLAZE)
	add_serial_unit_test(vector_serial_blaze_class)
ENDIF()

IF(HAVE_ARMADILLO)
	add_serial_unit_test(vector_serial_armadillo_class)
ENDIF()