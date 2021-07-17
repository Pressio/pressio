
#include "pressio_containers.hpp"

int main(int argc, char *argv[])
{

#if defined DO_EIGEN_DENSE
  using T = Eigen::VectorXd;
  using w_t = pressio::containers::DenseMatrix<T>;
  w_t(4);
#endif
#if defined DO_EIGEN_SPARSE
  using T = Eigen::SparseMatrix<double>;
  using w_t = pressio::containers::DenseMatrix<T>;
  w_t(4);
#endif

#if defined DO_EPETRA_V
  using T = Epetra_Vector;
  using w_t = pressio::containers::DenseMatrix<T>;
  w_t(4);
#endif

#if defined DO_TPETRA_V
  using T = Tpetra::Vector<>;
  using w_t = pressio::containers::DenseMatrix<T>;
  w_t(4);
#endif

#if defined DO_TPETRA_BLOCK_V
  using T = Tpetra::BlockVector<>;
  using w_t = pressio::containers::DenseMatrix<T>;
  w_t(4);
#endif

#if defined DO_KOKKOS_V
  using T = Kokkos::View<double*>;
  using w_t = pressio::containers::DenseMatrix<T>;
  w_t(4);
#endif

  return 0;
}
