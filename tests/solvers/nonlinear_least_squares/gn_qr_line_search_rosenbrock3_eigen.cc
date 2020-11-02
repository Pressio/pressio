
#include "pressio_solvers.hpp"
#include "./problems/eigen_rosenbrock_N3.hpp"


int main()
{
  std::string sentinel = "PASSED";

  using namespace pressio;
  using problem_t = solvers::test::EigenRosenbrock3;
  using state_w_t = typename problem_t::state_type;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;

  state_w_t x(3);
  x(0) = -1.5; x(1) = 1.1; x(2) = 1.2;

  using qr_solver_t = qr::QRSolver<mat_type, qr::Householder>;
  qr_solver_t qrSolver;

  auto GNSolver = pressio::solvers::nonlinear::createGaussNewtonQR(problem,x,qrSolver);
  GNSolver.setTolerance(1e-8);
  GNSolver.setMaxIterations(10);
  GNSolver.setUpdatingCriterion(pressio::solvers::nonlinear::update::armijo);

  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;
  if( std::abs(x(0) - 0.88235299629773) > 1e-6 ) sentinel = "FAILED";
  if( std::abs(x(1) - 0.70131382960973) > 1e-6 ) sentinel = "FAILED";
  if( std::abs(x(2) - 0.37498940601500) > 1e-6 ) sentinel = "FAILED";
  std::cout << sentinel << std::endl;
}
