
add_serial_unit_test(matrix_dense_serial_eigen_class)

IF(HAVE_TRILINOS)
	add_serial_unit_test(matrix_dense_serial_teuchos_class)
ENDIF()

