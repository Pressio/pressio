
TRIBITS_ADD_EXECUTABLE(
  generate_data
  SOURCES main.cc
  )

configure_file(./basis/basis_row189_col5.txt basis.txt COPYONLY)
