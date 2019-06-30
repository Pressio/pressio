#include <gtest/gtest.h>

#include "experimental/solvers_linear_factory.hpp"
#include "CONTAINERS_MATRIX"


TEST(solvers_linear_least_square_factory_eigen, solversTestLinearLeastSquareFactoryEigen)
{
  // Namespaces
  using namespace rompp;
  using namespace solvers;

  // Matrix typedefs
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = containers::Matrix<matrix_n_t>;

  // Define matrix
  matrix_w_t A(3, 2);

  // Create solvers
  auto solver = LinearSolvers::createLeastSquareIterativeSolver(A);
}


TEST(solvers_linear_least_square_eigen, solversTestLinearLeastSquareEigenCreateWithMatrixAndSolve)
{
  // Namespaces
  using namespace rompp;
  using namespace solvers;

  // Matrix typedefs
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = containers::Matrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = containers::Vector<vector_n_t>;

  // Define linear system
  vector_w_t b(3);
  (*b.data()) << 4.0, 1.0, 3.0;

  matrix_w_t A(3, 2);
  A.data()->insert(0, 0) += 1.0;
  A.data()->insert(0, 1) += 3.0;
  A.data()->insert(1, 0) += 2.0;
  A.data()->insert(1, 1) += 4.0;
  A.data()->insert(2, 0) += 1.0;
  A.data()->insert(2, 1) += 6.0;

  // Solve linear system
  auto solver = LinearSolvers::createLeastSquareIterativeSolver(A);
  auto x = solver.solve(b);

  // Expectations
  EXPECT_NEAR( x[0], -0.377, 1e-3 );
  EXPECT_NEAR( x[1],  0.662, 1e-3 );
}
