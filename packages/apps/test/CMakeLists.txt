
#include the lumped headers to access classes of solvers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lumped_include_headers)

ADD_SUBDIRECTORIES(
 burgers1d_epetra_wrapped
# burgers1d_eigen_wrapped
)
