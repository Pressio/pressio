
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "reference_apps_for_testing.hpp"


TEST(ode_implicit_euler, traits){
  using namespace rompp;
  
  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  using nres_t = typename app_t::residual_type;
  using njac_t = typename app_t::jacobian_type;
  using state_t = core::Vector<nstate_t>;
  using res_t = core::Vector<nres_t>;
  using jac_t = core::Matrix<njac_t>;

  static_assert(
    ode::meta::is_legitimate_model_for_explicit_ode<app_t>::value, "");
  static_assert(
    ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
    
  using stepper_t = ode::ImplicitStepper<
    ode::ImplicitSteppersEnum::Euler, state_t, jac_t, app_t>;

  using impl_t = typename stepper_t::base_t;
  using traits = ode::details::traits<impl_t>;
  
  // static_assert(ode::meta::is_implicit_euler_residual_standard_policy<
  // 		traits::residual_policy_t>::value,
  // 		"");
  // static_assert(ode::meta::is_implicit_euler_jacobian_standard_policy<
  // 		typename traits::jacobian_policy_t>::value,
  // 		"");
  
  ::testing::StaticAssertTypeEq<typename
  				traits::state_t, state_t>();
  ::testing::StaticAssertTypeEq<typename
  				traits::residual_t,res_t>();
  ::testing::StaticAssertTypeEq<typename
  				traits::jacobian_t,jac_t>();
  ::testing::StaticAssertTypeEq<typename
  				traits::scalar_t,double>();
  ::testing::StaticAssertTypeEq<typename
  				traits::model_t,app_t>();
  static_assert( traits::order_value == 1, "");
}


TEST(ode_implicit_euler, numerics){
  using namespace rompp;
  
  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  using nresidual_t = typename app_t::residual_type;
  using njacobian_t = typename app_t::jacobian_type;
  app_t appObj;
  
  using state_t = core::Vector<nstate_t>;
  using res_t = core::Vector<nresidual_t>;
  using jac_t = core::Matrix<njacobian_t>;
  state_t y(3);//appObj.y0);
  y[0] = 1.; y[1] = 2.; y[2] = 3.;
  
  res_t r(3);
  appObj.residual(*y.data(), *r.data(), 0.0);
  //std::cout << std::setprecision(14) << *r.data();

  // define stepper
  using stepper_t = ode::ImplicitStepper<
    ode::ImplicitSteppersEnum::Euler, state_t, jac_t, app_t>;
  stepper_t stepperObj(appObj, y);

  // define solver
  using namespace rompp::solvers;
  auto solverO = NonLinearSolvers::createIterativeSolver<
    nonlinear::NewtonRaphson,linear::Bicgstab>();
  // solverO.setMaxIterations(500);
  // solverO.setMaxNonLinearIterations(500);
  // solverO.setTolerance(1e-6);
  // solverO.setNonLinearTolerance(1e-6);

  // integrate in time
  int nSteps = 2;
  double dt = 0.01;
  ode::integrateNSteps(stepperObj, y, 0.0, dt, nSteps, solverO);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  appObj.analyticAdvanceBackEulerNSteps(dt, nSteps);
  
  EXPECT_DOUBLE_EQ(y[0], appObj.y0[0]);
  EXPECT_DOUBLE_EQ(y[1], appObj.y0[1]);
  EXPECT_DOUBLE_EQ(y[2], appObj.y0[2]);  
}
