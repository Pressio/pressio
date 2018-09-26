#include <gtest/gtest.h>

#include "experimental/solvers_linear_factory.hpp"
#include "CORE_MATRIX"


TEST(solvers_linear_dense_eigen, solversTestLinearDenseEigenCreateWithMatrixAndSolve)
{
  // Namespaces
  using namespace rompp;
  using namespace rompp::solvers;

  // Matrix typedefs
  using matrix_n_t = Eigen::MatrixXd;
  using matrix_w_t = core::Matrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = core::Vector<vector_n_t>;

  // Define linear system
  vector_w_t b(2);
  (*b.data()) << 7.0, 1.0;

  matrix_w_t A(2, 2);
  A(0, 0) +=  3.0;
  A(0, 1) += -1.0;
  A(1, 0) +=  2.0;
  A(1, 1) +=  3.0;

  // Solve linear system
  auto solver = LinearSolvers::createDirectSolver<linear::ColPivHouseholderQR>(A);
  auto x(solver.solve(b));

  // Expectations
  EXPECT_NEAR( x[0],  2.0, 1e-14 );
  EXPECT_NEAR( x[1], -1.0, 1e-14 );
}

TEST(solvers_linear_dense_eigen, solversTestLinearLeastSquareDenseEigenCreateWithMatrixAndSolve)
{
  // Namespaces
  using namespace rompp;
  using namespace rompp::solvers;

  // Matrix typedefs
  using matrix_n_t = Eigen::MatrixXd;
  using matrix_w_t = core::Matrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = core::Vector<vector_n_t>;

  // Define linear system
  vector_w_t b(3);
  (*b.data()) << 4.0, 1.0, 3.0;

  matrix_w_t A(3, 2);
  A(0, 0) += 1.0;
  A(0, 1) += 3.0;
  A(1, 0) += 2.0;
  A(1, 1) += 4.0;
  A(2, 0) += 1.0;
  A(2, 1) += 6.0;

  // Solve linear system
  auto solver = LinearSolvers::createDirectSolver<linear::ColPivHouseholderQR>(A);
  auto x(solver.solve(b));

  // Expectations
  EXPECT_NEAR( x[0], -0.377, 1e-3 );
  EXPECT_NEAR( x[1],  0.662, 1e-3 );
}
