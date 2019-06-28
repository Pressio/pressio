
#include <iostream>
#include "ALGEBRA_MATRIX"
#include "SOLVERS_LINEAR"

int main() {

  // Namespaces
  using namespace rompp;
  using namespace solvers;

  // Matrix typedefs
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = algebra::Matrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = algebra::Vector<vector_n_t>;

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

  // Solve linear least squares system
  auto solver = LinearSolvers::createLeastSquareIterativeSolver(A);
  auto x = solver.solve(b);  

  std::cout << "The solution of the linear system is: " << std::endl << *x.data() << std::endl;

  return 0;
}
