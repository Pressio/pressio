

add_serial_utest(containers_meta_basic .)
add_serial_utest(containers_meta_detect_typedefs .)
add_serial_utest(containers_config_defs_structs .)
add_serial_utest(containers_serial_printer .)
add_serial_utest(containers_vector_serial_eigen_dynamic_class .)
add_serial_utest(containers_vector_serial_eigen_dynamic_expr_templates .)
add_serial_utest(containers_matrix_vector_product_eigen .)
add_serial_utest(containers_matrix_matrix_product_eigen .)
add_serial_utest(containers_matrix_dense_serial_eigen_class .)
add_serial_utest(containers_multi_vector_eigen_product_eigen_vector .)
add_serial_utest(containers_sparse_matrix_product_multi_vector_eigen .)
add_serial_utest(containers_multi_vector_serial_eigen_dynamic_class .)
add_serial_utest(containers_multi_vector_eigen_dot_eigen_vector .)

# if(HAVE_TRILINOS)
# 	add_serial_utest(vector_serial_teuchos_dense_class . .)
# 	add_serial_utest(vector_serial_teuchos_dense_expr_templates . .)
# 	add_serial_utest(matrix_dense_serial_teuchos_class .)
# endif()

# if(HAVE_BLAZE)
# 	add_serial_utest(vector_serial_blaze_class .)
# endif()

# if(HAVE_ARMADILLO)
# 	add_serial_utest(vector_serial_armadillo_class .)
# endif()
