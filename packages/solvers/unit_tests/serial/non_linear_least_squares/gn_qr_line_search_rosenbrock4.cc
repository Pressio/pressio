
#include <gtest/gtest.h>
#include "SOLVERS_NONLINEAR"
#include "QR_BASIC"
#include "./problems/solvers_utest_serial_rosenbrock_N4.hpp"

TEST(solvers_nonlinear_least_squares,
     gn_qr_line_search_armijo_rosenbrock4){

  using namespace rompp;
  using problem_t = solvers::test::Rosenbrock4;
  using state_w_t = typename problem_t::state_type;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;

  state_w_t x(4);
  x[0] = -0.05; x[1] = 1.1; x[2] = 1.2; x[3] = 1.5;

  // define type of QR and GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  solvers::iterative::GaussNewtonQR<problem_t, qr_type,
					      lsearch_t> solver(problem, x);

  solver.setTolerance(1e-8);
  solver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 1., 1e-6 );
  EXPECT_NEAR( x(1), 1., 1e-6 );
  EXPECT_NEAR( x(2), 1., 1e-6 );
  EXPECT_NEAR( x(3), 1., 1e-6 );
}

TEST(solvers_nonlinear_least_squares,
     gn_qr_line_search_armijo_pass_sys_type_rosenbrock4){

  using namespace rompp;
  using problem_t = solvers::test::Rosenbrock4;
  using state_w_t = typename problem_t::state_type;
  using mat_type  = typename problem_t::jacobian_type;
  problem_t problem;

  state_w_t x(4);
  x[0] = -0.05; x[1] = 1.1; x[2] = 1.2; x[3] = 1.5;

  // define type of QR and GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  using converged_when_t = solvers::iterative::default_convergence;
  using gnsolver_t = solvers::iterative::GaussNewtonQR<
    qr_type, lsearch_t, converged_when_t, problem_t>;
  gnsolver_t solver(problem, x);

  solver.setTolerance(1e-8);
  solver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 1., 1e-6 );
  EXPECT_NEAR( x(1), 1., 1e-6 );
  EXPECT_NEAR( x(2), 1., 1e-6 );
  EXPECT_NEAR( x(3), 1., 1e-6 );
}
