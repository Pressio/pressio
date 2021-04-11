
#include "pressio_solvers.hpp"

int main()
{
  // Namespaces
  using namespace pressio;
  using namespace pressio::solvers;

  // Matrix typedefs
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = containers::SparseMatrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = containers::Vector<vector_n_t>;

  // Define linear system
  vector_w_t b(2);
  (*b.data()) << 7.0, 1.0;

  matrix_w_t A(2, 2);
  A.data()->insert(0, 0) +=  3.0;
  A.data()->insert(0, 1) += -1.0;
  A.data()->insert(1, 0) +=  2.0;
  A.data()->insert(1, 1) +=  3.0;

  vector_w_t y(2);

  // Solve linear system
  using solver_t = linear::Solver<linear::iterative::LSCG, matrix_w_t>;
  solver_t solver;
  solver.solve(A,b,y);

  std::string strOut = "PASSED";
  const auto e1 = std::abs(y(0) - (2.0)); 
  const auto e2 = std::abs(y(1) - (-1.0)); 
  if (e1>1e-13 or e2>1e-13) strOut = "FAILED";

  std::cout << *y.data() << std::endl;
  std::cout << strOut << std::endl;
  return 0;
}
