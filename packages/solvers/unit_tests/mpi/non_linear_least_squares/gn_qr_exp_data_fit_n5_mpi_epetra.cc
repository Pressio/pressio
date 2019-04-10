
#include <gtest/gtest.h>
#include "SOLVERS_NONLINEAR"
#include "QR_BASIC"
#include "./problems/solvers_utest_mpi_epetra_expon_data_fit_n5.hpp"

/*
 * test taken from:
 * http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
 * section 3.4
 * Data can be found at:
 * http://ftp.mcs.anl.gov/pub/MINPACK-2/tprobs/dedffj.f
 */

#ifdef HAVE_TRILINOS
TEST(solvers_nonlin_lsq,
     gn_qr_tsqr_exp_data_fit_n5_mpi_epetra){
  using namespace rompp;

  using problem_t   = solvers::test::ExpDataFitN5;
  using state_t	    = typename problem_t::state_type;
  using mat_type    = typename problem_t::jacobian_type;

  problem_t problem;
  state_t x(5);
  x[0] = 0.5;  x[1] = 1.5;
  x[2] = -1.0; x[3] = 0.01;
  x[4] = 0.02;

  // GaussNewton solver
  // define type of QR and GaussNewton solver
  using qr_algo = qr::TSQR;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using gn_t    = solvers::iterative::GaussNewtonQR<problem_t, qr_type>;
  gn_t GNSolver(problem, x);

  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<5; i++)
    EXPECT_NEAR( x(i), problem.trueS[i], 1e-6 );
}
#endif


TEST(solvers_nonlin_lsq,
     gn_qr_householder_exp_data_fit_n5_mpi_epetra){
  using namespace rompp;

  using problem_t   = solvers::test::ExpDataFitN5;
  using state_t	    = typename problem_t::state_type;
  using mat_type    = typename problem_t::jacobian_type;

  problem_t problem;
  state_t x(5);
  x[0] = 0.5;  x[1] = 1.5;
  x[2] = -1.0; x[3] = 0.01;
  x[4] = 0.02;

  // GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using gn_t    = solvers::iterative::GaussNewtonQR<problem_t, qr_type>;
  gn_t GNSolver(problem, x);

  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  for(auto i=0; i<5; i++)
    EXPECT_NEAR( x(i), problem.trueS[i], 1e-6 );
}
