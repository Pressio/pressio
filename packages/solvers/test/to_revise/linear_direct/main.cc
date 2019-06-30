
#include <iostream>
#include "CONTAINERS_MATRIX"
#include "SOLVERS_LINEAR"

int main() {

  // Namespaces
  using namespace rompp;
  using namespace rompp::solvers;

  // Matrix typedefs
  using matrix_n_t = Eigen::MatrixXd;
  using matrix_w_t = containers::Matrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = containers::Vector<vector_n_t>;

  // Define linear system
  vector_w_t b(3);
  (*b.data()) << 7.5, 1.0, 0.0;

  matrix_w_t A(3, 3);
  A(0, 0) +=  1.0;
  A(0, 1) +=  1.0;
  A(0, 2) +=  2.0;
  A(1, 0) += -1.0;
  A(1, 1) +=  0.0;
  A(1, 2) +=  1.0;
  A(2, 0) +=  0.0;
  A(2, 1) +=  6.0;
  A(2, 2) +=  1.0;

  // Solve linear system
  auto solver = LinearSolvers::createDirectSolver<linear::ColPivHouseholderQR>(A);
  auto x = solver.solve(b);

  std::cout << "The solution of the linear system is: " << std::endl << *x.data() << std::endl;

  return 0;
}
