
#include <gtest/gtest.h>
#include "experimental/solvers_linear_iterative_factory.hpp"
#include "matrix/concrete/core_matrix_dense_serial_eigen.hpp"
#include "matrix/concrete/core_matrix_sparse_serial_eigen.hpp"
#include "vector/concrete/core_vector_serial_eigen.hpp"


TEST(solvers_linear_iterative_factory, solversTestLinearIterativeFactoryEigen)
{
  // Namespaces
  using namespace solvers;

  // Matrix typedefs
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = core::matrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = core::vector<vector_n_t>;

  // Define linear system
  vector_w_t b(3);
  matrix_w_t A(3, 3);

  // Solve linear system
  auto solver = LinearIterativeSolvers::createSolver<CG>(A);
  auto x = solver.solve(b);
}

