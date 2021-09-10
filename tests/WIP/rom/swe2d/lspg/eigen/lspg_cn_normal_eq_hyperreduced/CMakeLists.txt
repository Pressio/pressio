
configure_file(../../../basis_hyperreduced/PhiSamplePlusStencil.txt PhiSamplePlusStencil.txt COPYONLY)
configure_file(../../../basis_hyperreduced/info_file.txt info_file.txt COPYONLY)
configure_file(../../../basis_hyperreduced/sample_mesh_gids.txt sample_mesh_gids.txt COPYONLY)
configure_file(../../../basis_hyperreduced/sample_mesh_plus_stencil_gids.txt sample_mesh_plus_stencil_gids.txt COPYONLY)
configure_file(../../../basis_hyperreduced/sample_mesh_relative_indices.txt sample_mesh_relative_indices.txt COPYONLY)
configure_file(../../../basis_hyperreduced/basis.txt basis.txt COPYONLY)

add_serial_exe_and_test(swe2d_eigen_lspg_crank_nicolson_gn_normal_eqns_hyper rom main.cc "PASSED")
