
#include <gtest/gtest.h>
#include "SOLVERS_ALL"
#include "./problems/solvers_utest_mpi_epetra_expon_data_fit_n11.hpp"

/*
 * test taken from:
 * http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
 * section 3.5
 * Data can be found at:
 * http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dgdffj.f
 */

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
TEST(solvers_nonlin_lsq,
     gn_qr_tsqr_exp_data_fit_n11_mpi_epetra){
  using namespace pressio;

  using problem_t   = solvers::test::ExpDataFitN11;
  using state_t	    = typename problem_t::state_type;
  using mat_type    = typename problem_t::jacobian_type;

  problem_t problem;
  state_t x(problem.numUn_);
  x(0) = 1.3;
  x(1) = 6.5e-1;
  x(2) = 6.5e-1;
  x(3) = 7.0e-1;
  x(4) = 6.0e-1;
  x(5) = 3.0;
  x(6) = 5.0;
  x(7) = 7.0;
  x(8) = 2.0;
  x(9) = 4.5;
  x(10) = 5.5;

  /* NOTE: this problem with GN only works with line search */
  // GaussNewton solver
  using qr_algo = qr::TSQR;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  using gn_t = solvers::iterative::GaussNewtonQR<
		    problem_t, qr_type, lsearch_t>;
  gn_t GNSolver(problem, x);

  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<problem.numUn_; i++)
    EXPECT_NEAR( x(i), problem.trueS[i], 1e-6 );
}
#endif


TEST(solvers_nonlin_lsq,
     gn_qr_householder_exp_data_fit_n5_mpi_epetra){
  using namespace pressio;

  using problem_t   = solvers::test::ExpDataFitN11;
  using state_t	    = typename problem_t::state_type;
  using mat_type    = typename problem_t::jacobian_type;

  problem_t problem;
  state_t x(problem.numUn_);
  x(0) = 1.3;
  x(1) = 6.5e-1;
  x(2) = 6.5e-1;
  x(3) = 7.0e-1;
  x(4) = 6.0e-1;
  x(5) = 3.0;
  x(6) = 5.0;
  x(7) = 7.0;
  x(8) = 2.0;
  x(9) = 4.5;
  x(10) = 5.5;

  // GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  using gn_t = solvers::iterative::GaussNewtonQR<
		    problem_t, qr_type, lsearch_t>;
  gn_t GNSolver(problem, x);

  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<problem.numUn_; i++)
    EXPECT_NEAR( x(i), problem.trueS[i], 1e-6 );
}
