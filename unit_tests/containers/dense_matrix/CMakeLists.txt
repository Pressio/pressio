
set(exeRootName containers_dense_matrix)

#############
### EIGEN ###
#############
# column-major matrix dynamic
set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/eigen_dynamic.cc)
add_serial_utest(${exeRootName}_eigen_dynamic_colmajor ${SRC1})
target_compile_definitions(utest_${exeRootName}_eigen_dynamic_colmajor PRIVATE USE_COLMAJ)
# row-major matrix dynamic
set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/eigen_dynamic.cc)
add_serial_utest(${exeRootName}_eigen_dynamic_rowmajor ${SRC1})
target_compile_definitions(utest_${exeRootName}_eigen_dynamic_rowmajor PRIVATE USE_ROWMAJ)

# column-major matrix static 
set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/eigen_static.cc)
add_serial_utest(${exeRootName}_eigen_static_colmajor ${SRC1})
target_compile_definitions(utest_${exeRootName}_eigen_static_colmajor PRIVATE USE_COLMAJ)
# row-major matrix static
set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/eigen_static.cc)
add_serial_utest(${exeRootName}_eigen_static_rowmajor ${SRC1})
target_compile_definitions(utest_${exeRootName}_eigen_static_rowmajor PRIVATE USE_ROWMAJ)

