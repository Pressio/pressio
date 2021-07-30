
#include "mocks.hpp"
#include "pressio_solvers_nonlinear.hpp"

TEST(solvers_nonlinear, operators_residual_jacobian)
{
  using ::testing::Return;
  using ::testing::_;

  using mock_t = MockSystemRJ<true>;
  mock_t sysObj;
  EXPECT_CALL(sysObj, createResidual()).Times(1);
  EXPECT_CALL(sysObj, createJacobian()).Times(1);
  EXPECT_CALL(sysObj, residual(_,_)).Times(2);
  EXPECT_CALL(sysObj, jacobian(_,_)).Times(1);

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using scalar_type = typename mock_t::scalar_type;
  using residual_type = typename mock_t::residual_type;
  using jacobian_type = typename mock_t::jacobian_type;

  using op_t = psimpl::ResidualJacobianOperators<residual_type, jacobian_type, scalar_type>;
  Eigen::VectorXd state;
  op_t operators(sysObj, state);
  scalar_type norm = {};
  operators.computeOperators(sysObj, state, norm);

  const auto & r = operators.residualCRef();
  for (int i=0; i<r.size(); ++i) EXPECT_DOUBLE_EQ(r(i), 1.);

  scalar_type rNorm = {};
  operators.residualNorm(sysObj, state, rNorm);
  EXPECT_DOUBLE_EQ(rNorm, std::sqrt(5.));
}

TEST(solvers_nonlinear, operators_residual_jacobian_fused)
{
  using ::testing::Return;
  using ::testing::_;

  using mock_t = MockSystemRJFused<true>;
  mock_t sysObj;
  EXPECT_CALL(sysObj, createResidual()).Times(1);
  EXPECT_CALL(sysObj, createJacobian()).Times(1);
  EXPECT_CALL(sysObj, residualAndJacobian(_,_,_,_)).Times(2);

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using scalar_type = typename mock_t::scalar_type;
  using residual_type = typename mock_t::residual_type;
  using jacobian_type = typename mock_t::jacobian_type;

  using op_t = psimpl::ResidualJacobianOperators<residual_type, jacobian_type, scalar_type>;
  Eigen::VectorXd state;
  op_t operators(sysObj, state);
  scalar_type norm = {};
  operators.computeOperators(sysObj, state, norm);

  const auto & r = operators.residualCRef();
  for (int i=0; i<r.size(); ++i) EXPECT_DOUBLE_EQ(r(i), 1.);

  scalar_type rNorm = {};
  operators.residualNorm(sysObj, state, rNorm);
  EXPECT_DOUBLE_EQ(rNorm, std::sqrt(5.));
}

TEST(solvers_nonlinear, operators_hessian_gradient_hg_api_enable_jacobian_computation)
{
  using ::testing::Return;
  using ::testing::_;

  MockSystemHessGrad sysObj;
  EXPECT_CALL(sysObj, createGradient()).Times(1);
  EXPECT_CALL(sysObj, createHessian()).Times(1);
  EXPECT_CALL(sysObj, gradient(_,_,_,_,_)).Times(1);
  EXPECT_CALL(sysObj, hessian(_,_)).Times(1);
  EXPECT_CALL(sysObj, residualNorm(_,_,_)).Times(1);

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using scalar_type = typename MockSystemHessGrad::scalar_type;
  using gradient_type = typename MockSystemHessGrad::gradient_type;
  using hessian_type = typename MockSystemHessGrad::hessian_type;

  using op_t = psimpl::HessianGradientOperatorsHGApi<hessian_type, gradient_type, scalar_type>;
  Eigen::VectorXd state;
  op_t operators(sysObj, state);
  scalar_type norm = {};
  operators.computeOperators(sysObj, state, norm, true);

  // check gradient
  const auto & g = operators.gradientCRef();
  // the minus sign is because of the convention 
  for (int i=0; i<g.size(); ++i) EXPECT_DOUBLE_EQ(g(i), -3.4);

  // check hessian
  const auto & H = operators.hessianCRef();
  const auto goldH = sysObj.viewConcrete().goldHessian();
  EXPECT_TRUE(H.isApprox(goldH, 1e-10));

  // check norm
  scalar_type rNorm = {};
  operators.residualNorm(sysObj, state, rNorm);
  EXPECT_DOUBLE_EQ(rNorm, 8.);
}

TEST(solvers_nonlinear, operators_hessian_gradient_hg_api_disable_jacobian_computation)
{
  using ::testing::Return;
  using ::testing::_;

  MockSystemHessGrad sysObj;
  EXPECT_CALL(sysObj, createGradient()).Times(1);
  EXPECT_CALL(sysObj, createHessian()).Times(1);
  EXPECT_CALL(sysObj, gradient(_,_,_,_,_)).Times(1);
  EXPECT_CALL(sysObj, hessian(_,_)).Times(0);
  EXPECT_CALL(sysObj, residualNorm(_,_,_)).Times(1);

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using scalar_type = typename MockSystemHessGrad::scalar_type;
  using gradient_type = typename MockSystemHessGrad::gradient_type;
  using hessian_type = typename MockSystemHessGrad::hessian_type;

  using op_t = psimpl::HessianGradientOperatorsHGApi<hessian_type, gradient_type, scalar_type>;
  Eigen::VectorXd state;
  op_t operators(sysObj, state);
  scalar_type norm = {};
  operators.computeOperators(sysObj, state, norm, false);

  // check gradient
  const auto & g = operators.gradientCRef();
  // the minus sign is because of the convention 
  for (int i=0; i<g.size(); ++i) EXPECT_DOUBLE_EQ(g(i), -3.4);

  // check hessian
  const auto & H = operators.hessianCRef();
  for (int i=0; i<H.rows(); ++i){
    for (int j=0; j<H.cols(); ++j){
      EXPECT_DOUBLE_EQ(H(i,j), 0.);
    }
  }

  // check norm
  scalar_type rNorm = {};
  operators.residualNorm(sysObj, state, rNorm);
  EXPECT_DOUBLE_EQ(rNorm, 8.);
}


TEST(solvers_nonlinear, operators_hessian_gradient_hg_api_fused_enable_jacobian_computation)
{
  using ::testing::Return;
  using ::testing::_;

  MockSystemHessGradFused sysObj;
  EXPECT_CALL(sysObj, createGradient()).Times(1);
  EXPECT_CALL(sysObj, createHessian()).Times(1);
  EXPECT_CALL(sysObj, hessianAndGradient(_,_,_,_,_,_)).Times(1);
  EXPECT_CALL(sysObj, residualNorm(_,_,_)).Times(1);

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using scalar_type = typename MockSystemHessGradFused::scalar_type;
  using gradient_type = typename MockSystemHessGradFused::gradient_type;
  using hessian_type = typename MockSystemHessGradFused::hessian_type;

  using op_t = psimpl::HessianGradientOperatorsHGApi<hessian_type, gradient_type, scalar_type>;
  Eigen::VectorXd state;
  op_t operators(sysObj, state);
  scalar_type norm = {};
  operators.computeOperators(sysObj, state, norm, true);

  // check gradient
  const auto & g = operators.gradientCRef();
  // the minus sign is because of the convention 
  for (int i=0; i<g.size(); ++i) EXPECT_DOUBLE_EQ(g(i), -3.4);

  // check hessian
  const auto & H = operators.hessianCRef();
  const auto goldH = sysObj.viewConcrete().goldHessian();
  EXPECT_TRUE(H.isApprox(goldH, 1e-10));

  // check norm
  scalar_type rNorm = {};
  operators.residualNorm(sysObj, state, rNorm);
  EXPECT_DOUBLE_EQ(rNorm, 8.);
}


TEST(solvers_nonlinear, operators_hessian_gradient_hg_api_fused_disable_jacobian_computation)
{
  using ::testing::Return;
  using ::testing::_;

  MockSystemHessGradFused sysObj;
  EXPECT_CALL(sysObj, createGradient()).Times(1);
  EXPECT_CALL(sysObj, createHessian()).Times(1);
  EXPECT_CALL(sysObj, hessianAndGradient(_,_,_,_,_,_)).Times(1);
  EXPECT_CALL(sysObj, residualNorm(_,_,_)).Times(1);

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using scalar_type = typename MockSystemHessGradFused::scalar_type;
  using gradient_type = typename MockSystemHessGradFused::gradient_type;
  using hessian_type = typename MockSystemHessGradFused::hessian_type;

  using op_t = psimpl::HessianGradientOperatorsHGApi<hessian_type, gradient_type, scalar_type>;
  Eigen::VectorXd state;
  op_t operators(sysObj, state);
  scalar_type norm = {};
  operators.computeOperators(sysObj, state, norm, false);

  // check gradient
  const auto & g = operators.gradientCRef();
  // the minus sign is because of the convention 
  for (int i=0; i<g.size(); ++i) EXPECT_DOUBLE_EQ(g(i), -3.4);

  // check hessian
  const auto & H = operators.hessianCRef();
  for (int i=0; i<H.rows(); ++i){
    for (int j=0; j<H.cols(); ++j){
      EXPECT_DOUBLE_EQ(H(i,j), 0.);
    }
  }

  // check norm
  scalar_type rNorm = {};
  operators.residualNorm(sysObj, state, rNorm);
  EXPECT_DOUBLE_EQ(rNorm, 8.);
}



TEST(solvers_nonlinear, operators_hessian_gradient_rj_api_enable_jacobian_computation)
{
  using ::testing::Return;
  using ::testing::_;

  using mock_t = MockSystemRJ<true>; 
  mock_t sysObj;
  EXPECT_CALL(sysObj, createResidual()).Times(1);
  EXPECT_CALL(sysObj, createJacobian()).Times(1);
  EXPECT_CALL(sysObj, residual(_,_)).Times(2);
  EXPECT_CALL(sysObj, jacobian(_,_)).Times(1);

  using scalar_type = typename mock_t::scalar_type;
  using residual_type = typename mock_t::residual_type;
  using jacobian_type = typename mock_t::jacobian_type;
  using gradient_type = Eigen::VectorXd;
  using hessian_type = Eigen::MatrixXd;

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using op_t = psimpl::HessianGradientOperatorsRJApiNoWeighting<hessian_type, 
      gradient_type, residual_type, jacobian_type, scalar_type>;
  Eigen::VectorXd state(5);
  op_t operators(sysObj, state);

  scalar_type norm = {};
  operators.computeOperators(sysObj, state, norm, true);
  EXPECT_DOUBLE_EQ(norm, std::sqrt(5));

  // check gradient
  const auto & g = operators.gradientCRef();
  // the minus sign is because of the convention 
  EXPECT_DOUBLE_EQ(g(0), -8.);
  EXPECT_DOUBLE_EQ(g(1), -1.);
  EXPECT_DOUBLE_EQ(g(2), -7.);
  EXPECT_DOUBLE_EQ(g(3), -0.);
  EXPECT_DOUBLE_EQ(g(4), -13.);

  // check hessian
  const auto & H = operators.hessianCRef();
  const auto goldH = sysObj.viewConcrete().goldHessian();
  EXPECT_TRUE(H.isApprox(goldH, 1e-10));

  // check norm
  scalar_type rNorm = {};
  operators.residualNorm(sysObj, state, rNorm);
  EXPECT_DOUBLE_EQ(rNorm, std::sqrt(5.));
}


TEST(solvers_nonlinear, operators_hessian_gradient_rj_api_disable_jacobian_computation)
{
  using ::testing::Return;
  using ::testing::_;

  using mock_t = MockSystemRJ<true>; 
  mock_t sysObj;
  EXPECT_CALL(sysObj, createResidual()).Times(1);
  EXPECT_CALL(sysObj, createJacobian()).Times(1);
  EXPECT_CALL(sysObj, residual(_,_)).Times(2);
  EXPECT_CALL(sysObj, jacobian(_,_)).Times(0);

  using scalar_type = typename mock_t::scalar_type;
  using residual_type = typename mock_t::residual_type;
  using jacobian_type = typename mock_t::jacobian_type;
  using gradient_type = Eigen::VectorXd;
  using hessian_type = Eigen::MatrixXd;

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using op_t = psimpl::HessianGradientOperatorsRJApiNoWeighting<hessian_type, 
      gradient_type, residual_type, jacobian_type, scalar_type>;
  Eigen::VectorXd state(5);
  op_t operators(sysObj, state);

  scalar_type norm = {};
  operators.computeOperators(sysObj, state, norm, false);
  EXPECT_DOUBLE_EQ(norm, std::sqrt(5));

  // check hessian
  const auto & H = operators.hessianCRef();
  for (int i=0; i<H.rows(); ++i){
    for (int j=0; j<H.cols(); ++j){
      EXPECT_DOUBLE_EQ(H(i,j), 0.);
    }
  }

  // check gradient
  const auto & g = operators.gradientCRef();
  // the minus sign is because of the convention 
  EXPECT_DOUBLE_EQ(g(0), 0.);
  EXPECT_DOUBLE_EQ(g(1), 0.);
  EXPECT_DOUBLE_EQ(g(2), 0.);
  EXPECT_DOUBLE_EQ(g(3), 0.);
  EXPECT_DOUBLE_EQ(g(4), 0.);

  // check norm
  scalar_type rNorm = {};
  operators.residualNorm(sysObj, state, rNorm);
  EXPECT_DOUBLE_EQ(rNorm, std::sqrt(5.));
}


TEST(solvers_nonlinear, operators_hessian_gradient_weighted_rj_api_enable_jacobian_computation)
{
  using ::testing::Return;
  using ::testing::_;

  using mock_t = MockSystemRJ<true>; 
  mock_t sysObj;
  EXPECT_CALL(sysObj, createResidual()).Times(1);
  EXPECT_CALL(sysObj, createJacobian()).Times(1);
  EXPECT_CALL(sysObj, residual(_,_)).Times(2);
  EXPECT_CALL(sysObj, jacobian(_,_)).Times(1);

  using scalar_type = typename mock_t::scalar_type;
  using residual_type = typename mock_t::residual_type;
  using jacobian_type = typename mock_t::jacobian_type;
  using gradient_type = Eigen::VectorXd;
  using hessian_type = Eigen::MatrixXd;
  using weighting_type = MockWeightingOperator;

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using op_t = psimpl::WeightedHessianGradientOperatorsRJApi<hessian_type, 
      gradient_type, residual_type, jacobian_type, scalar_type, const weighting_type & >;
  Eigen::VectorXd state(5);
  weighting_type M;
  EXPECT_CALL(M, compute1(_,_)).Times(2);
  EXPECT_CALL(M, compute2(_,_)).Times(1);
  op_t operators(sysObj, state, M);

  scalar_type norm = {};
  operators.computeOperators(sysObj, state, norm, true);
  EXPECT_DOUBLE_EQ(norm, std::sqrt(15.));

  // check gradient
  const auto & g = operators.gradientCRef();
  // the minus sign is because of the convention 
  EXPECT_DOUBLE_EQ(g(0), -24.);
  EXPECT_DOUBLE_EQ(g(1), -3.);
  EXPECT_DOUBLE_EQ(g(2), -21.);
  EXPECT_DOUBLE_EQ(g(3), -0.);
  EXPECT_DOUBLE_EQ(g(4), -39.);

  // check hessian
  const auto & H = operators.hessianCRef();
  for (int j=0; j<H.cols(); j++){
    EXPECT_DOUBLE_EQ(H(0,j), 17.6);
    EXPECT_DOUBLE_EQ(H(1,j), 2.2);
    EXPECT_DOUBLE_EQ(H(2,j), 15.4);
    EXPECT_DOUBLE_EQ(H(3,j), 0.);
    EXPECT_DOUBLE_EQ(H(4,j), 28.6);
  }
  // check norm
  scalar_type rNorm = {};
  operators.residualNorm(sysObj, state, rNorm);
  EXPECT_DOUBLE_EQ(rNorm, std::sqrt(15.));
}

TEST(solvers_nonlinear, operators_hessian_gradient_weighted_rj_api_disable_jacobian_computation)
{
  using ::testing::Return;
  using ::testing::_;

  using mock_t = MockSystemRJ<true>; 
  mock_t sysObj;
  EXPECT_CALL(sysObj, createResidual()).Times(1);
  EXPECT_CALL(sysObj, createJacobian()).Times(1);
  EXPECT_CALL(sysObj, residual(_,_)).Times(2);
  EXPECT_CALL(sysObj, jacobian(_,_)).Times(0);

  using scalar_type = typename mock_t::scalar_type;
  using residual_type = typename mock_t::residual_type;
  using jacobian_type = typename mock_t::jacobian_type;
  using gradient_type = Eigen::VectorXd;
  using hessian_type = Eigen::MatrixXd;
  using weighting_type = MockWeightingOperator;

  namespace psimpl = pressio::nonlinearsolvers::impl; 
  using op_t = psimpl::WeightedHessianGradientOperatorsRJApi<hessian_type, 
      gradient_type, residual_type, jacobian_type, scalar_type, const weighting_type & >;
  Eigen::VectorXd state(5);
  weighting_type M;
  EXPECT_CALL(M, compute1(_,_)).Times(2);
  EXPECT_CALL(M, compute2(_,_)).Times(0);
  op_t operators(sysObj, state, M);

  scalar_type norm = {};
  operators.computeOperators(sysObj, state, norm, false);
  EXPECT_DOUBLE_EQ(norm, std::sqrt(15.));

  // check gradient
  const auto & g = operators.gradientCRef();
  // the minus sign is because of the convention 
  EXPECT_DOUBLE_EQ(g(0), 0.);
  EXPECT_DOUBLE_EQ(g(1), 0.);
  EXPECT_DOUBLE_EQ(g(2), 0.);
  EXPECT_DOUBLE_EQ(g(3), 0.);
  EXPECT_DOUBLE_EQ(g(4), 0.);

  // check hessian
  const auto & H = operators.hessianCRef();
  for (int i=0; i<H.rows(); i++){
   for (int j=0; j<H.cols(); j++){
    EXPECT_DOUBLE_EQ(H(i,j), 0.);
   }
  }
  // check norm
  scalar_type rNorm = {};
  operators.residualNorm(sysObj, state, rNorm);
  EXPECT_DOUBLE_EQ(rNorm, std::sqrt(15.));
}

