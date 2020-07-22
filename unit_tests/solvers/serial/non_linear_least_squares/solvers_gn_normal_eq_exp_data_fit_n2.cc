
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"
#include "./problems/solvers_utest_serial_expon_data_fit_n2.hpp"

TEST(solvers_nonlinear_least_squares,
     gn_normal_eq_lscg_exp_data_fit_n2){

  using namespace pressio;

  using problem_t   = solvers::test::ExpDataFitN2;
  using state_t	    = typename problem_t::state_type;
  problem_t problem;
  state_t x(2); x[0] = 2.0; x[1] = 0.25;

  using hessian_t = containers::Matrix<Eigen::MatrixXd>;

  // linear solver type
  using solver_tag	= solvers::linear::iterative::LSCG;
  using linear_solver_t = solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  // GaussNewton solver
  using solver = pressio::solvers::nonlinear::composeGaussNewton_t<
    problem_t, pressio::solvers::nonlinear::DefaultUpdate, linear_solver_t>;
  solver GNSolver(problem, x, linSolver);
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-7 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-7 );
}


TEST(solvers_nonlinear_least_squares,
     gn_normal_eq_grad_norm_lscg_exp_data_fit_n2)
{

  using namespace pressio;

  using problem_t   = solvers::test::ExpDataFitN2;
  using vec_t	  = typename problem_t::state_type;
  using mat_t	  = typename problem_t::jacobian_type;
  using state_t	  = vec_t;
  using hessian_t = mat_t;

  problem_t problem;
  state_t x(2); x[0] = 2.0; x[1] = 0.25;

  // linear solver type
  using solver_tag	= solvers::linear::iterative::LSCG;
  using linear_solver_t = solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  // GaussNewton solver
  using solver = pressio::solvers::nonlinear::composeGaussNewton_t<
    problem_t, pressio::solvers::nonlinear::DefaultUpdate,
    linear_solver_t>;
  solver GNSolver(problem, x, linSolver);
  auto criterion = pressio::solvers::nonlinear::stop::whenGradientAbsoluteNormBelowTolerance;
  GNSolver.setStoppingCriterion(criterion);
  // 1e3 is chosen to test the convergence condition
  GNSolver.setTolerance(1e3);
  GNSolver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4153884777105201 , 1e-11);
  EXPECT_NEAR( x(1), 0.26879930127342189, 1e-11);
}


TEST(solvers_nonlinear_least_squares,
     gn_normal_eq_rel_grad_norm_lscg_exp_data_fit_n2){

  using namespace pressio;

  using problem_t   = solvers::test::ExpDataFitN2;
  using vec_t	  = typename problem_t::state_type;
  using mat_t	  = typename problem_t::jacobian_type;
  using state_t	  = vec_t;
  using hessian_t = mat_t;

  problem_t problem;
  state_t x(2); x[0] = 2.0; x[1] = 0.25;

  // linear solver type
  using solver_tag	= solvers::linear::iterative::LSCG;
  using linear_solver_t = solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  // GaussNewton solver
  using solver = pressio::solvers::nonlinear::composeGaussNewton_t<
    problem_t, pressio::solvers::nonlinear::DefaultUpdate,
    linear_solver_t>;
  solver GNSolver(problem, x, linSolver);
  GNSolver.setStoppingCriterion(solver::stop::whenGradientRelativeNormBelowTolerance);
  GNSolver.setTolerance(1e-5);
  GNSolver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.41728158794844, 1e-11);
  EXPECT_NEAR( x(1), 0.26465375115797, 1e-11);
}
