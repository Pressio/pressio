
#include "pressio_containers.hpp"

int main(int argc, char *argv[])
{

#if defined DO_EIGEN_DENSE
  using T = Eigen::MatrixXd;
  using w_t = pressio::containers::Vector<T>;
  w_t(4);
#endif
#if defined DO_EIGEN_SPARSE
  using T = Eigen::SparseMatrix<double>;
  using w_t = pressio::containers::Vector<T>;
  w_t(4);
#endif

#if defined DO_EPETRA_MV
  using T = Epetra_MultiVector;
  using w_t = pressio::containers::Vector<T>;
  w_t(4);
#endif

#if defined DO_EPETRA_SM
  using T = Epetra_CrsMatrix;
  using w_t = pressio::containers::Vector<T>;
  w_t(4);
#endif


#if defined DO_KOKKOS_DM
  using T = Kokkos::View<double**>;
  using w_t = pressio::containers::Vector<T>;
  w_t(4);
#endif

#if defined DO_KOKKOS_SM
  using execution_space = Kokkos::DefaultExecutionSpace;
  using T = KokkosSparse::CrsMatrix<double, int, execution_space, void, int>;
  using w_t = pressio::containers::Vector<T>;
  w_t(4);
#endif

#if defined DO_TPETRA_MV
  using T = Tpetra::MultiVector<>;
  using w_t = pressio::containers::Vector<T>;
  w_t(4);
#endif

#if defined DO_TPETRA_SM
  using T = Tpetra::CrsMatrix<>;
  using w_t = pressio::containers::Vector<T>;
  w_t(4);
#endif


#if defined DO_TPETRA_BLOCK_MV
  using T = Tpetra::BlockMultiVector<>;
  using w_t = pressio::containers::Vector<T>;
  w_t(4);
#endif

  return 0;
}
