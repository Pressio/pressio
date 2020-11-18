
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(qr_eigen_dense_matrix_factorization_out_of_place
    qr_eigen_dense_matrix_factorization_out_of_place.cc )
  add_serial_utest(qr_eigen_dense_matrix_out_of_place_and_solve
    qr_eigen_dense_matrix_out_of_place_and_solve.cc)
endif()
