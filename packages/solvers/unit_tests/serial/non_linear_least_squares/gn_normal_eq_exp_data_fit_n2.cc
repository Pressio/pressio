
#include <gtest/gtest.h>
#include "SOLVERS_NONLINEAR"
#include "QR_BASIC"
#include "./problems/solvers_utest_serial_expon_data_fit_n2.hpp"

TEST(solvers_nonlinear_least_squares,
     gn_normal_eq_lscg_exp_data_fit_n2){

  using namespace rompp;

  using problem_t   = solvers::test::ExpDataFitN2;
  using state_t	    = typename problem_t::state_type;
  using sc_t	    = double;
  problem_t problem;
  state_t x(2); x[0] = 2.0; x[1] = 0.25;

  // linear solver type and GaussNewton solver
  using solver_tag = solvers::linear::iterative::LSCG;
  solvers::iterative::GaussNewton<sc_t, solver_tag,
				  solvers::EigenIterative> GNSolver;
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-7 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-7 );
}


TEST(solvers_nonlinear_least_squares,
     gn_normal_eq_pass_sys_type_lscg_exp_data_fit_n2){

  using namespace rompp;

  using problem_t   = solvers::test::ExpDataFitN2;
  using vec_t	  = typename problem_t::state_type;
  using mat_t	  = typename problem_t::jacobian_type;

  using state_t	  = vec_t;
  using hessian_t = mat_t;
  using sc_t	  = double;

  problem_t problem;
  state_t x(2); x[0] = 2.0; x[1] = 0.25;

  // define linear solver type and GaussNewton solver
  using solver_tag = solvers::linear::iterative::LSCG;
  using converged_when_t = solvers::iterative::default_convergence;
  using gn_t = solvers::iterative::GaussNewton<sc_t, solver_tag,
					       solvers::EigenIterative,
					       converged_when_t, problem_t,
					       hessian_t>;
  gn_t GNSolver(problem, x);
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-7 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-7 );
}
