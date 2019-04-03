
#include <gtest/gtest.h>
#include "SOLVERS_NONLINEAR"
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
  using namespace rompp;

  using problem_t   = solvers::test::ExpDataFitN5;
  using state_t	    = typename problem_t::state_type;
  using sc_t	    = double;

  problem_t problem;
  state_t x(problem.numUn_);
  x[0] = 0.5;  x[1] = 1.5;
  x[2] = -1.0; x[3] = 0.01;
  x[4] = 0.02;

  // linear solver type and GaussNewton solver
  using solver_tag = solvers::linear::iterative::LSCG;
  solvers::iterative::GaussNewton<sc_t, solver_tag,
  				  solvers::EigenIterative> GNSolver;
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<5; i++)
    EXPECT_NEAR( x(i), problem.trueS[i], 1e-6 );
}



TEST(solvers_nonlin_lsq,
     gn_normal_eq_lscg_line_search_armijo_pass_type_exp_data_fit_n5_mpi_epetra){
  using namespace rompp;

  using problem_t   = solvers::test::ExpDataFitN5;
  using state_t	    = typename problem_t::state_type;
  using sc_t	    = double;

  problem_t problem;
  state_t x(5);
  x[0] = 0.5;  x[1] = 1.5;
  x[2] = -1.0; x[3] = 0.01;
  x[4] = 0.02;

  // linear solver type and GaussNewton solver
  using solver_tag = solvers::linear::iterative::LSCG;
  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  using gn_t = solvers::iterative::GaussNewtonLineSearch<
    sc_t, solver_tag, solvers::EigenIterative, lsearch_t>;

  gn_t GNSolver;
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<5; i++)
    EXPECT_NEAR( x(i), problem.trueS[i], 1e-6 );
}


TEST(solvers_nonlin_lsq,
     gn_normal_eq_pass_sys_type_lscg_exp_data_fit_n5_mpi_epetra){
  using namespace rompp;

  using problem_t   = solvers::test::ExpDataFitN5;
  using state_t	    = typename problem_t::state_type;
  using sc_t	    = double;
  using hessian_t   = core::Matrix<Eigen::MatrixXd>;

  problem_t problem;
  state_t x(5);
  x[0] = 0.5;  x[1] = 1.5;
  x[2] = -1.0; x[3] = 0.01;
  x[4] = 0.02;

  // linear solver type and GaussNewton solver
  using solver_tag = solvers::linear::iterative::LSCG;
  using converged_when_t = solvers::iterative::default_convergence;
  using gn_t = solvers::iterative::GaussNewton<sc_t, solver_tag,
  				  solvers::EigenIterative,
				  converged_when_t,
					       problem_t, hessian_t>;
  gn_t GNSolver(problem, x);

  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<5; i++)
    EXPECT_NEAR( x(i), problem.trueS[i], 1e-6 );
}
