
#include <gtest/gtest.h>
#include "CORE_MATRIX"
#include "SOLVERS_LINEAR"

TEST(solvers_linear_iterative, LSCGDenseEigen)
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

  vector_w_t y(2);

  // Solve linear system
  using solver_t = EigenIterative<linear::LSCG, matrix_w_t>;
  solver_t solver;
  solver.solve(A,b,y);

  // Expectations
  EXPECT_NEAR( y[0], -0.377, 1e-3 );
  EXPECT_NEAR( y[1],  0.662, 1e-3 );
}
