
#include "pressio_solvers_linear.hpp"

int main()
{
  using namespace pressio;
  using namespace pressio::linearsolvers;

  using matrix_t = Eigen::SparseMatrix<double>;
  using vector_t = Eigen::VectorXd;

  // Define linear system
  vector_t b(2);
  b << 7.0, 1.0;

  matrix_t A(2, 2);
  A.insert(0, 0) +=  3.0;
  A.insert(0, 1) += -1.0;
  A.insert(1, 0) +=  2.0;
  A.insert(1, 1) +=  3.0;

  vector_t y(2);
  using solver_t = Solver<iterative::LSCG, matrix_t>;
  solver_t solver;
  solver.solve(A,b,y);

  std::string strOut = "PASSED";
  const auto e1 = std::abs(y(0) - (2.0)); 
  const auto e2 = std::abs(y(1) - (-1.0)); 
  if (e1>1e-13 or e2>1e-13) strOut = "FAILED";
  std::cout << y << std::endl;
  std::cout << strOut << std::endl;

  return 0;
}
