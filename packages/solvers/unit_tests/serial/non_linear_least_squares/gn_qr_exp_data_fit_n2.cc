
#include <gtest/gtest.h>
#include "SOLVERS_NONLINEAR"
#include "QR_BASIC"
#include "./problems/solvers_utest_serial_expon_data_fit_n2.hpp"

TEST(solvers_nonlinear_least_squares, gn_qr_exp_data_fit_n2){
  using namespace rompp;

  using problem_t   = solvers::test::ExpDataFitN2;

  using state_w_t = typename problem_t::state_type;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;

  state_w_t x(2); x[0] = 2.0; x[1] = 0.25;

  // define type of QR and GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  solvers::iterative::GaussNewtonQR<problem_t, qr_type> solver(problem, x);
  solver.setTolerance(1e-8);
  solver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-9 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-9 );

  // print summary from timers
  #ifdef HAVE_TEUCHOS_TIMERS
  utils::TeuchosPerformanceMonitor::stackedTimersReportSerial();
  #endif
}


TEST(solvers_nonlinear_least_squares, gn_qr_only_2_steps_exp_data_fit_n2){
  using namespace rompp;

  using problem_t   = solvers::test::ExpDataFitN2;
  using state_w_t = typename problem_t::state_type;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;

  state_w_t x(2); x[0] = 2.0; x[1] = 0.25;

  // define type of QR and GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using converged_when_t
    = solvers::iterative::converged_when::completingNumMaxIters;
  solvers::iterative::GaussNewtonQR<problem_t, qr_type, converged_when_t> solver(problem, x);
  // setting 2 max iters so that in combination with the
  // above convergence method, the solver will exit afte 2 steps
  solver.setMaxIterations(2);
  solver.solve(problem, x);

  std::cout << std::setprecision(16) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.415361667771343 , 1e-8);
  EXPECT_NEAR( x(1), 0.2648293802571118 , 1e-8);
}


TEST(solvers_nonlinear_least_squares,
     gn_qr_only_pass_sys_type_exp_data_fit_n2){

  using namespace rompp;

  using problem_t   = solvers::test::ExpDataFitN2;
  using state_w_t = typename problem_t::state_type;
  using mat_t   = typename problem_t::jacobian_type;
  problem_t problem;

  state_w_t x(2); x[0] = 2.0; x[1] = 0.25;

  // define type of QR solver and GaussNewton solver
  using qr_algo		 = qr::Householder;
  using qr_type		 = qr::QRSolver<mat_t, qr_algo>;
  using converged_when_t = solvers::iterative::default_convergence;
  using gnsolver_t	 =
    solvers::iterative::GaussNewtonQR<qr_type, converged_when_t, problem_t>;
  gnsolver_t GNSolver(problem, x);
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-9 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-9 );
}
