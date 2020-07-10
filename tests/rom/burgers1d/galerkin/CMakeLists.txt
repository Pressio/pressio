
add_subdirectory(eigen)
add_subdirectory(arbitrary_data_structures)

if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_subdirectory(epetra)
endif()

