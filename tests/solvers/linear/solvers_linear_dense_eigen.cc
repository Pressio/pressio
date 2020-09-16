
#include "pressio_solvers.hpp"

int main()
{
  // Namespaces
  using namespace pressio;
  using namespace pressio::solvers;

  // Matrix typedefs
  using matrix_n_t = Eigen::MatrixXd;
  using matrix_w_t = containers::Matrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = containers::Vector<vector_n_t>;

  // Define linear system
  vector_w_t b(3);
  (*b.data()) << 4.0, 1.0, 3.0;

  matrix_w_t A(3, 2);
  A(0, 0) = 1.0;
  A(0, 1) = 3.0;
  A(1, 0) = 2.0;
  A(1, 1) = 4.0;
  A(2, 0) = 1.0;
  A(2, 1) = 6.0;

  vector_w_t y(2);

#if defined(DO_LSCG)
  using tag = linear::iterative::LSCG;
#elif defined(DO_ColPivHSQR)
  using tag = linear::direct::ColPivHouseholderQR;
#elif defined(DO_HouseholderQR)
  using tag = linear::direct::HouseholderQR;
#endif

  using solver_t = linear::Solver<tag, matrix_w_t>;

  solver_t solver;
  solver.solve(A,b,y);

  std::string strOut = "PASSED";
  const auto e1 = std::abs(y[0] - (-0.377)); 
  const auto e2 = std::abs(y[1] - (0.662)); 
  if (e1>1e-3 or e2>1e-3) strOut = "FAILED";

  std::cout << *y.data() << std::endl;

  std::cout << strOut << std::endl;
}