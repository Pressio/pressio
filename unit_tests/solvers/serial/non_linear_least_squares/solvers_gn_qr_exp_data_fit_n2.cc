
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"
#include "./problems/solvers_utest_serial_expon_data_fit_n2.hpp"

TEST(solvers_nonlinear_least_squares, gn_qr_exp_data_fit_n2){
  using namespace pressio;

  using problem_t   = solvers::test::ExpDataFitN2;
  using state_w_t = typename problem_t::state_type;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;
  state_w_t x(2); x[0] = 2.0; x[1] = 0.25;

  using qr_solver_type = qr::QRSolver<mat_type, qr::Householder>;
  qr_solver_type qrSolver;

  // GaussNewton solver
  using solver = pressio::solvers::nonlinear::composeGaussNewtonQR_t<
    problem_t, pressio::solvers::nonlinear::DefaultUpdate,
    pressio::solvers::nonlinear::StopWhenCorrectionNormBelowTol,
    qr_solver_type>;
  solver GNsolver(problem, x, qrSolver);
  GNsolver.setTolerance(1e-8);
  GNsolver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-9 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-9 );

  // print summary from timers
  #ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  utils::TeuchosPerformanceMonitor::stackedTimersReportSerial();
  #endif
}

TEST(solvers_nonlinear_least_squares, gn_qr_only_2_steps_exp_data_fit_n2)
{
  using namespace pressio;

  using problem_t   = solvers::test::ExpDataFitN2;
  using state_w_t = typename problem_t::state_type;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;
  state_w_t x(2); x[0] = 2.0; x[1] = 0.25;

  using qr_solver_type = qr::QRSolver<mat_type, qr::Householder>;
  qr_solver_type qrSolver;

  // GaussNewton solver
  using solver = pressio::solvers::nonlinear::composeGaussNewtonQR_t<
    problem_t, pressio::solvers::nonlinear::DefaultUpdate,
    pressio::solvers::nonlinear::StopAfterMaxIters,
    qr_solver_type>;
  solver GNsolver(problem, x, qrSolver);

  // setting 2 max iters so the solver will exit afte 2 steps
  GNsolver.setMaxIterations(2);
  GNsolver.solve(problem, x);

  std::cout << std::setprecision(16) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.415361667771343 , 1e-8);
  EXPECT_NEAR( x(1), 0.2648293802571118 , 1e-8);
}
