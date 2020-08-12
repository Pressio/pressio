
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"
#include "./problems/solvers_utest_serial_rosenbrock_N3.hpp"

TEST(solvers_nonlinear_least_squares,
     gn_qr_line_search_armijo_rosenbrock3){

  using namespace pressio;
  using problem_t = solvers::test::Rosenbrock3;
  using state_w_t = typename problem_t::state_type;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;

  state_w_t x(3);
  x[0] = -1.5; x[1] = 1.1; x[2] = 1.2;

  using qr_solver_t = qr::QRSolver<mat_type, qr::Householder>;
  qr_solver_t qrSolver;
  // GaussNewton solver
  using solver = pressio::solvers::nonlinear::composeGaussNewtonQR_t<
    problem_t, pressio::solvers::nonlinear::armijoUpdate,
    qr_solver_t>;
  solver GNsolver(problem, x, qrSolver);
  GNsolver.setTolerance(1e-8);
  GNsolver.setMaxIterations(10);

  GNsolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 0.88235299629773, 1e-6 );
  EXPECT_NEAR( x(1), 0.70131382960973, 1e-6 );
  EXPECT_NEAR( x(2), 0.37498940601500, 1e-6 );
}
