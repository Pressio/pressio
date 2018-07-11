
#include <gtest/gtest.h>
#include "experimental/solvers_linear_iterative_factory.hpp"
#include "matrix/concrete/core_matrix_dense_serial_eigen.hpp"
#include "matrix/concrete/core_matrix_sparse_serial_eigen.hpp"
#include "vector/concrete/core_vector_serial_eigen.hpp"


TEST(solvers_linear_iterative_factory, solversTestLinearIterativeFactoryEigen)
{
  // Typedefs and namespaces
  using namespace solvers;
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = core::matrix<matrix_n_t>;

  // Define linear system
  matrix_w_t A(3, 3);
  auto solver = LinearIterativeSolvers::createSolver<CG>(A);
}

