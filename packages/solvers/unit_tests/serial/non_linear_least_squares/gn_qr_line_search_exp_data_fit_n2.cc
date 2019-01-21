
#include <gtest/gtest.h>
#include "SOLVERS_NONLINEAR"
#include "QR_BASIC"
#include "./problems/solvers_utest_serial_expon_data_fit_n2.hpp"

TEST(solvers_nonlinear_least_squares,
     gn_qr_line_search_armijo_exp_data_fit_n2){

  using namespace rompp;
  using problem_t   = solvers::test::ExpDataFitN2;
  using state_w_t = typename problem_t::state_type;
  using sc_t	  = double;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;

  // define type of QR and GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  solvers::iterative::GaussNewtonQRLineSearch<sc_t, qr_type,
					      lsearch_t> solver;
  solver.setTolerance(1e-8);

  state_w_t x(2); x[0] = 2.0; x[1] = 0.25;
  solver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-9 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-9 );
}


TEST(solvers_nonlinear_least_squares,
     gn_qr_line_search_armijo_only_2steps_exp_data_fit_n2){

  using namespace rompp;
  using problem_t   = solvers::test::ExpDataFitN2;
  using state_w_t = typename problem_t::state_type;
  using sc_t	  = double;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;

  // define type of QR and GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  using converged_when_t
    = solvers::iterative::converged_when::completingNumMaxIters;
  solvers::iterative::GaussNewtonQRLineSearch<sc_t, qr_type, lsearch_t,
					      converged_when_t> solver;
  // setting max iters so that in combination with the
  // above convergence method, the solver will exit after target steps
  solver.setMaxIterations(2);

  state_w_t x(2); x[0] = 2.0; x[1] = 0.25;
  solver.solve(problem, x);

  std::cout << std::setprecision(16) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.415361667771343 , 1e-8);
  EXPECT_NEAR( x(1), 0.2648293802571118 , 1e-8);
}


TEST(solvers_nonlinear_least_squares,
     gn_qr_line_search_armijo_pass_sys_type_exp_data_fit_n2){

  using namespace rompp;
  using problem_t   = solvers::test::ExpDataFitN2;
  using state_w_t = typename problem_t::state_type;
  using sc_t	  = double;
  using mat_t  = typename problem_t::jacobian_type;
  problem_t problem;

  state_w_t x(2); x[0] = 2.0; x[1] = 0.25;

  // define type of QR solver and GaussNewton solver
  using qr_algo		 = qr::Householder;
  using qr_type		 = qr::QRSolver<mat_t, qr_algo>;
  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  using converged_when_t = solvers::iterative::default_convergence;
  using gnsolver_t	 =
    solvers::iterative::GaussNewtonQRLineSearch<sc_t, qr_type,
						lsearch_t,
						converged_when_t,
						problem_t>;
  gnsolver_t GNSolver(problem, x);
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-9 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-9 );
}
