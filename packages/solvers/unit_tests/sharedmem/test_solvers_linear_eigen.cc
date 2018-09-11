
#include <gtest/gtest.h>

#include "experimental/solvers_linear_factory.hpp"
#include "CORE_MATRIX"
// #include "CORE_VECTOR"
// #include "matrix/concrete/core_matrix_dense_serial_eigen.hpp"
// #include "matrix/concrete/core_matrix_sparse_serial_eigen.hpp"
// #include "vector/concrete/core_vector_serial_eigen.hpp"

TEST(solvers_linear_iterative_eigen, solversTestLinearIterativeEigenCreateWithMatrixAndSolve)
{
  // Namespaces
  using namespace solvers;

  // Matrix typedefs
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = core::Matrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = core::Vector<vector_n_t>;

  // Define linear system
  vector_w_t b(2);
  (*b.data()) << 7.0, 1.0;

  matrix_w_t A(2, 2);
  A.data()->insert(0, 0) +=  3.0;
  A.data()->insert(0, 1) += -1.0;
  A.data()->insert(1, 0) +=  2.0;
  A.data()->insert(1, 1) +=  3.0;

  // Solve linear system
  auto solver = LinearSolvers::createIterativeSolver<linear::Bicgstab>(A);
  auto x = solver.solve(b);

  // Expectations
  EXPECT_NEAR( x[0],  2.0, 1e-14 );
  EXPECT_NEAR( x[1], -1.0, 1e-14 );
}


TEST(solvers_linear_iterative_eigen, solversTestLinearIterativeEigenCreateWithoutMatrixAndSolve)
{
  // Namespaces
  using namespace solvers;

  // Matrix typedefs
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = core::Matrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = core::Vector<vector_n_t>;

  // Define linear system
  vector_w_t b(2);
  (*b.data()) << 7.0, 1.0;

  matrix_w_t A(2, 2);
  A.data()->insert(0, 0) +=  3.0;
  A.data()->insert(0, 1) += -1.0;
  A.data()->insert(1, 0) +=  2.0;
  A.data()->insert(1, 1) +=  3.0;

  // Solve linear system
  auto solver = LinearSolvers::createIterativeSolver<linear::Bicgstab, matrix_w_t>();
  auto x = solver.solve(A, b);

  // Expectations
  EXPECT_NEAR( x[0],  2.0, 1e-14 );
  EXPECT_NEAR( x[1], -1.0, 1e-14 );
}
