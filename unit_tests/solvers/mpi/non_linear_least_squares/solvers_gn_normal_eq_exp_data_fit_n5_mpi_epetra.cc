
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"
#include "./problems/solvers_utest_mpi_epetra_expon_data_fit_n5.hpp"

/*
 * test taken from:
 * http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
 * section 3.4
 * Data can be found at:
 * http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dedffj.f
 */


TEST(solvers_nonlin_lsq,
     gn_normal_eq_lscg_exp_data_fit_n5_mpi_epetra){
  using namespace pressio;

  using problem_t   = solvers::test::ExpDataFitN5;
  using state_t	    = typename problem_t::state_type;

  problem_t problem;
  state_t x(problem.numUn_);
  x[0] = 0.5;  x[1] = 1.5;
  x[2] = -1.0; x[3] = 0.01;
  x[4] = 0.02;

  using hessian_t = containers::Matrix<Eigen::MatrixXd>;

  // linear solver type
  using solver_tag  = solvers::linear::iterative::LSCG;
  using linear_solver_t = solvers::iterative::EigenIterative<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  using gn_t = solvers::iterative::GaussNewton<
  linear_solver_t, problem_t, hessian_t>;
  gn_t GNSolver(problem, x, linSolver);

  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<5; i++)
    EXPECT_NEAR( x(i), problem.trueS[i], 1e-6 );
}



TEST(solvers_nonlin_lsq,
     gn_normal_eq_lscg_line_search_armijo_pass_type_exp_data_fit_n5_mpi_epetra){
  using namespace pressio;

  using problem_t   = solvers::test::ExpDataFitN5;
  using state_t	    = typename problem_t::state_type;

  problem_t problem;
  state_t x(5);
  x[0] = 0.5;  x[1] = 1.5;
  x[2] = -1.0; x[3] = 0.01;
  x[4] = 0.02;

  using hessian_t = containers::Matrix<Eigen::MatrixXd>;

  // linear solver type
  using solver_tag  = solvers::linear::iterative::LSCG;
  using linear_solver_t = solvers::iterative::EigenIterative<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  using gn_t = solvers::iterative::GaussNewton<
  linear_solver_t, problem_t, hessian_t, lsearch_t>;
  gn_t GNSolver(problem, x, linSolver);
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<5; i++)
    EXPECT_NEAR( x(i), problem.trueS[i], 1e-6 );
}


TEST(solvers_nonlin_lsq,
     gn_normal_eq_pass_sys_type_lscg_exp_data_fit_n5_mpi_epetra){
  using namespace pressio;

  using problem_t   = solvers::test::ExpDataFitN5;
  using state_t	    = typename problem_t::state_type;
  using hessian_t   = containers::Matrix<Eigen::MatrixXd>;

  problem_t problem;
  state_t x(5);
  x[0] = 0.5;  x[1] = 1.5;
  x[2] = -1.0; x[3] = 0.01;
  x[4] = 0.02;

  using hessian_t = containers::Matrix<Eigen::MatrixXd>;

  // linear solver type
  using solver_tag  = solvers::linear::iterative::LSCG;
  using linear_solver_t = solvers::iterative::EigenIterative<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  using converged_when_t = solvers::iterative::default_convergence;
  using gn_t = solvers::iterative::GaussNewton<
  linear_solver_t, problem_t, hessian_t, converged_when_t>;
  gn_t GNSolver(problem, x, linSolver);

  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<5; i++)
    EXPECT_NEAR( x(i), problem.trueS[i], 1e-6 );
}
