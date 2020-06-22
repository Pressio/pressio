
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"
#include "./problems/solvers_utest_serial_rosenbrock_N4.hpp"


TEST(solvers_nonlinear_least_squares, gn_res_jac_api_rosenbrock4){
  using namespace pressio;

  Eigen::Vector4d state;

  using problem_t = solvers::test::Rosenbrock4;
  using state_t   = typename problem_t::state_type;
  using hessian_t = containers::Matrix<Eigen::MatrixXd>;

  problem_t problem;
  state_t x(4);
  x[0] = -0.05; x[1] = 1.1; x[2] = 1.2; x[3] = 1.5;

  using lin_tag      = solvers::linear::direct::HouseholderQR;
  using lin_solver_t = solvers::linear::Solver<lin_tag, hessian_t>;
  lin_solver_t linSolver;

  using solver = pressio::solvers::nonlinear::composeGaussNewton_t<
    problem_t, pressio::solvers::nonlinear::DefaultUpdate,
    pressio::solvers::nonlinear::StopWhenCorrectionNormBelowTol,
    lin_solver_t>;
  solver GNSolver(problem, x, linSolver);


  GNSolver.setTolerance(1e-1);
  GNSolver.solve(problem, x);
  state = *x.data();
  std::cout << std::setprecision(14) << state << std::endl;

  std::vector<double> gold = {1.00000001567414e+00, 9.99999999124769e-01,
  			      9.99999996519930e-01, 9.99999988898883e-01};
  EXPECT_NEAR( state(0), gold[0], 1e-6 );
  EXPECT_NEAR( state(1), gold[1], 1e-6 );
  EXPECT_NEAR( state(2), gold[2], 1e-6 );
  EXPECT_NEAR( state(3), gold[3], 1e-6 );
}


TEST(solvers_nonlinear_least_squares, gn_hess_grad_api_rosenbrock4){
  using namespace pressio;

  Eigen::Vector4d state;

  using problem_t = solvers::test::Rosenbrock4HessGradApi;
  using state_t   = typename problem_t::state_type;
  using hessian_t = typename problem_t::hessian_type;

  problem_t problem;
  state_t x(4);
  x[0] = -0.05; x[1] = 1.1; x[2] = 1.2; x[3] = 1.5;

  // linear solver type
  using lin_tag      = solvers::linear::direct::HouseholderQR;
  using lin_solver_t = solvers::linear::Solver<lin_tag, hessian_t>;
  lin_solver_t linSolver;

  using solver = pressio::solvers::nonlinear::composeGaussNewton_t<
    problem_t, pressio::solvers::nonlinear::DefaultUpdate,
    pressio::solvers::nonlinear::StopWhenCorrectionNormBelowTol,
    lin_solver_t>;
  solver GNSolver(problem, x, linSolver);

  GNSolver.setTolerance(1e-1);
  GNSolver.solve(problem, x);
  state = *x.data();
  std::cout << std::setprecision(14) << state << std::endl;

  std::vector<double> gold = {1.00000001567414e+00, 9.99999999124769e-01,
  			      9.99999996519930e-01, 9.99999988898883e-01};
  EXPECT_NEAR( state(0), gold[0], 1e-6 );
  EXPECT_NEAR( state(1), gold[1], 1e-6 );
  EXPECT_NEAR( state(2), gold[2], 1e-6 );
  EXPECT_NEAR( state(3), gold[3], 1e-6 );
}
