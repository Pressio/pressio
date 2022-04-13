
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bases_p_0.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bases_p_1.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bases_p_2.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bases_p_3.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bases_p_4.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bases_p_5.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bases_p_6.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bases_p_7.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/state_vec_rows_wrt_full_mesh_p_0.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/state_vec_rows_wrt_full_mesh_p_1.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/state_vec_rows_wrt_full_mesh_p_2.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/state_vec_rows_wrt_full_mesh_p_3.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/state_vec_rows_wrt_full_mesh_p_4.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/state_vec_rows_wrt_full_mesh_p_5.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/state_vec_rows_wrt_full_mesh_p_6.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/state_vec_rows_wrt_full_mesh_p_7.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tpetra_map_gids_0.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tpetra_map_gids_1.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tpetra_map_gids_2.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tpetra_map_gids_3.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tpetra_map_gids_4.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/check_proj.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/check_rec.py ${CMAKE_CURRENT_BINARY_DIR}  COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/rank_partition_mapping.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/odrom_tpetra.cc)
# add_utest_mpi(${TESTING_LEVEL}_odrom_trilinos_tpetra "${SRC1}" gTestMain_tpetra 5)

set(exename odrom_trilinos_tpetra)
add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/odrom_tpetra.cc)
target_link_libraries(${exename} ${MPI_CXX_LIBRARIES})

add_test(NAME ${TESTING_LEVEL}_odrom_trilinos_tpetra
COMMAND ${CMAKE_COMMAND}
-DEXENAME=$<TARGET_FILE:${exename}>
-P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
)

# set_tests_properties(${testname}
#   PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
#   FAIL_REGULAR_EXPRESSION "FAILED"
# )