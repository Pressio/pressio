
#include "pressio_solvers_linear.hpp"

int main()
{
  using namespace pressio;
  using namespace pressio::linearsolvers;

  using matrix_t = Eigen::MatrixXd;
  using vector_t = Eigen::VectorXd;

  vector_t b(3);
  b << 4.0, 1.0, 3.0;

  matrix_t A(3, 2);
  A(0, 0) = 1.0;
  A(0, 1) = 3.0;
  A(1, 0) = 2.0;
  A(1, 1) = 4.0;
  A(2, 0) = 1.0;
  A(2, 1) = 6.0;

  vector_t y(2);

#if defined(DO_LSCG)
  using tag = iterative::LSCG;
#elif defined(DO_ColPivHSQR)
  using tag = direct::ColPivHouseholderQR;
#elif defined(DO_HouseholderQR)
  using tag = direct::HouseholderQR;
#endif

  using solver_t = Solver<tag, matrix_t>;
  solver_t solver;
  solver.solve(A,b,y);

  std::string strOut = "PASSED";
  const auto e1 = std::abs(y(0) - (-0.377)); 
  const auto e2 = std::abs(y(1) - (0.662)); 
  if (e1>1e-3 or e2>1e-3) strOut = "FAILED";
  std::cout << y << std::endl;
  std::cout << strOut << std::endl;

  return 0;
}