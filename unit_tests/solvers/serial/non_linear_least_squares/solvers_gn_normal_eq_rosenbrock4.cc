
#include <gtest/gtest.h>
#include "SOLVERS_NONLINEAR"
#include "SOLVERS_EXPERIMENTAL"
#include "./problems/solvers_utest_serial_rosenbrock_N4.hpp"


TEST(solvers_nonlinear_least_squares, gn_hess_grad_api_rosenbrock4){
  using namespace pressio;

  Eigen::Vector4d stateDefApi;
  Eigen::Vector4d stateHessApi;

  // solve problem using regular GN
  {
    using problem_t = solvers::test::Rosenbrock4;
    using state_t   = typename problem_t::state_type;
    using hessian_t = containers::Matrix<Eigen::MatrixXd>;

    problem_t problem;
    state_t x(4);
    x[0] = -0.05; x[1] = 1.1; x[2] = 1.2; x[3] = 1.5;

    using solver_tag	= solvers::linear::direct::HouseholderQR;
    using linear_solver_t = solvers::direct::EigenDirect<solver_tag, hessian_t>;
    linear_solver_t linSolver;

    using gn_t = solvers::iterative::GaussNewton<linear_solver_t, problem_t, hessian_t>;
    gn_t GNSolver(problem, x, linSolver);
    GNSolver.setTolerance(1e-1);
    GNSolver.solve(problem, x);
    stateDefApi = *x.data();
    std::cout << std::setprecision(14) << stateDefApi << std::endl;
  }

  // solve problem using regular GN
  {
    using problem_t = solvers::test::Rosenbrock4HessGradApi;
    using state_t   = typename problem_t::state_type;
    using hessian_t = typename problem_t::hessian_type;

    problem_t problem;
    state_t x(4);
    x[0] = -0.05; x[1] = 1.1; x[2] = 1.2; x[3] = 1.5;

    // linear solver type
    using solver_tag	= solvers::linear::direct::HouseholderQR;
    using linear_solver_t = solvers::direct::EigenDirect<solver_tag, hessian_t>;
    linear_solver_t linSolver;

    // GaussNewton solver
    using gn_t = solvers::iterative::GaussNewton<linear_solver_t, problem_t>;
    gn_t GNSolver(problem, x, linSolver);
    GNSolver.setTolerance(1e-1);
    GNSolver.solve(problem, x);
    stateHessApi = *x.data();
    std::cout << std::setprecision(14) << stateHessApi << std::endl;
  }

  EXPECT_DOUBLE_EQ( stateDefApi(0), stateHessApi(0) );
  EXPECT_DOUBLE_EQ( stateDefApi(1), stateHessApi(1) );
  EXPECT_DOUBLE_EQ( stateDefApi(2), stateHessApi(2) );
  EXPECT_DOUBLE_EQ( stateDefApi(3), stateHessApi(3) );

  std::vector<double> gold = {1.00000001567414e+00, 9.99999999124769e-01,
			      9.99999996519930e-01, 9.99999988898883e-01};
  EXPECT_NEAR( stateDefApi(0), gold[0], 1e-6 );
  EXPECT_NEAR( stateDefApi(1), gold[1], 1e-6 );
  EXPECT_NEAR( stateDefApi(2), gold[2], 1e-6 );
  EXPECT_NEAR( stateDefApi(3), gold[3], 1e-6 );
}
