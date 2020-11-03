
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_subdirectory(arbitrary_data_structures)
  add_subdirectory(eigen)
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_subdirectory(epetra)
endif()
