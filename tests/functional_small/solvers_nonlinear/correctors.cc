
#include "mocks.hpp"
#include "pressio/solvers_nonlinear.hpp"

TEST(solvers_nonlinear, correctors_residual_jacobian)
{
  using ::testing::Return;
  using ::testing::_;

  using mock_system_t = MockSystemRJ<true>;
  mock_system_t sysObj;
  EXPECT_CALL(sysObj, createResidual()).Times(1);
  EXPECT_CALL(sysObj, createJacobian()).Times(1);
  EXPECT_CALL(sysObj, residual(_,_)).Times(2);
  EXPECT_CALL(sysObj, jacobian(_,_)).Times(2);

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using scalar_type = typename mock_system_t::scalar_type;
  using residual_type = typename mock_system_t::residual_type;
  using jacobian_type = typename mock_system_t::jacobian_type;

  using state_type = Eigen::VectorXd;
  state_type state(5);

  using mock_lin_solver_t = MockLinearSolver<jacobian_type, residual_type, state_type>;
  mock_lin_solver_t linSolver;
  EXPECT_CALL(linSolver, solve(_,_,_)).Times(2);

  using op_t = psimpl::ResidualJacobianOperators<residual_type, jacobian_type, scalar_type>;
  using corrector_t = psimpl::RJCorrector<op_t, state_type, mock_lin_solver_t &>;
  corrector_t corrector(sysObj, state, linSolver);

  scalar_type norm = {};
  corrector.computeOperators(sysObj, state, norm);
  const auto & J = corrector.jacobianCRef();
  const auto & r = corrector.residualCRef();
  linSolver.solve(J,r,state);

  corrector.computeOperators(sysObj, state, norm);
  linSolver.solve(J,r,state);
}

TEST(solvers_nonlinear, correctors_hessian_gradient)
{
  using ::testing::Return;
  using ::testing::_;

  using mock_system_t = MockSystemRJ<true>;
  mock_system_t sysObj;
  EXPECT_CALL(sysObj, createResidual()).Times(1);
  EXPECT_CALL(sysObj, createJacobian()).Times(1);
  EXPECT_CALL(sysObj, residual(_,_)).Times(2);
  EXPECT_CALL(sysObj, jacobian(_,_)).Times(2);

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using scalar_type = typename mock_system_t::scalar_type;
  using residual_type = typename mock_system_t::residual_type;
  using jacobian_type = typename mock_system_t::jacobian_type;

  using state_type = Eigen::VectorXd;
  using gradient_type = Eigen::VectorXd;
  using hessian_type = Eigen::MatrixXd;
  state_type state(5);

  using mock_lin_solver_t = MockLinearSolver<jacobian_type, residual_type, state_type>;
  mock_lin_solver_t linSolver;
  EXPECT_CALL(linSolver, solve(_,_,_)).Times(2);

  using op_t = psimpl::HessianGradientOperatorsRJApiNoWeighting<
    hessian_type, gradient_type, residual_type, jacobian_type, scalar_type>;
  using corrector_t = psimpl::HessianGradientCorrector<op_t, state_type, mock_lin_solver_t &>;
  corrector_t corrector(sysObj, state, linSolver);

  scalar_type norm = {};
  corrector.computeOperators(sysObj, state, norm);
  const auto & H = corrector.hessianCRef();
  const auto & g = corrector.gradientCRef();
  linSolver.solve(H,g,state);

  corrector.computeOperators(sysObj, state, norm);
  linSolver.solve(H,g,state);
}

