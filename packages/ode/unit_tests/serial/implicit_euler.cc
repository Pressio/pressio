
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

  // // the RESIDUAL standard policy is
  // using res_pol_t = ode::policy::ImplicitEulerResidualStandardPolicy<
  //   state_t, app_t, res_t>;
  // // the JACOBIAN standard policy is
  // using jac_pol_t = ode::policy::ImplicitEulerJacobianStandardPolicy<
  //   state_t, app_t, jac_t>;

  using stepper_t = ode::ImplicitStepper<
    ode::ImplicitEnum::Euler,
    state_t, res_t, jac_t, app_t, void>; /*aux stepper NOT needed for backEuler*/

  using traits = ode::details::traits<stepper_t>;
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


TEST(ode_implicit_euler, numericsStdPoliciesDefaultCreated){
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

  using stepper_t = ode::ImplicitStepper<
    ode::ImplicitEnum::Euler,
    state_t, res_t, jac_t, app_t>; /*aux stepper NOT needed for backEuler*/
  stepper_t stepperObj(appObj, y);

  // define solver
  using lin_solver_t = solvers::EigenIterative<solvers::linear::Bicgstab, jac_t>;
  solvers::NewtonRaphson<double, lin_solver_t> solverO;

  // integrate in time
  int nSteps = 2;
  double dt = 0.01;
  ode::integrateNSteps(stepperObj, y, 0.0, dt, nSteps, solverO);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  appObj.analyticAdvanceBackEulerNSteps(dt, nSteps);

  EXPECT_DOUBLE_EQ(y[0], appObj.y[0]);
  EXPECT_DOUBLE_EQ(y[1], appObj.y[1]);
  EXPECT_DOUBLE_EQ(y[2], appObj.y[2]);
}


TEST(ode_implicit_euler, guesserLambda){
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

  using stepper_t = ode::ImplicitStepper<
    ode::ImplicitEnum::Euler,
    state_t, res_t, jac_t, app_t>; /*aux stepper NOT needed for backEuler*/
  stepper_t stepperObj(appObj, y);

  // define solver
  using lin_algo_t = solvers::linear::Bicgstab;
  using lin_solver_t = solvers::EigenIterative<lin_algo_t, jac_t>;
  solvers::NewtonRaphson<double, lin_solver_t> solverO;

  // integrate in time

  auto testLambda = [](size_t step, double time, state_t & y){
  		      y[0] = -20.; y[1] = -20.; y[2] = -20.; };

  int nSteps = 1;
  double dt = 0.01;
  ode::integrateNSteps(stepperObj, y, 0.0, dt,
		       nSteps, solverO, testLambda);

  std::cout << std::setprecision(14) << *y.data() << "\n";
  // appObj.analyticAdvanceBackEulerNSteps(dt, nSteps);
  // EXPECT_DOUBLE_EQ(y[0], appObj.y[0]);
  // EXPECT_DOUBLE_EQ(y[1], appObj.y[1]);
  // EXPECT_DOUBLE_EQ(y[2], appObj.y[2]);
}


TEST(ode_implicit_euler, numericsStdResidualPolPassedByUser){
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

  //**********************
  // define policies and stepper
  //**********************
  using res_pol_t = ode::policy::ImplicitResidualStandardPolicy<
    state_t, app_t, res_t>;
  using jac_pol_t = ode::policy::ImplicitJacobianStandardPolicy<
    state_t, app_t, jac_t>;

  using stepper_t = ode::ImplicitStepper<
    ode::ImplicitEnum::Euler,
    state_t, res_t, jac_t, app_t, void, /*aux stepper NOT needed for backEuler*/
    res_pol_t, jac_pol_t>;
  stepper_t stepperObj(appObj, res_pol_t(), jac_pol_t(), y);

  //**********************
  // define solver
  //**********************
  using lin_solver_t = solvers::EigenIterative<solvers::linear::Bicgstab, jac_t>;
  solvers::NewtonRaphson<double, lin_solver_t> solverO;

  // integrate in time
  int nSteps = 2;
  double dt = 0.01;
  ode::integrateNSteps(stepperObj, y, 0.0, dt, nSteps, solverO);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  appObj.analyticAdvanceBackEulerNSteps(dt, nSteps);

  EXPECT_DOUBLE_EQ(y[0], appObj.y[0]);
  EXPECT_DOUBLE_EQ(y[1], appObj.y[1]);
  EXPECT_DOUBLE_EQ(y[2], appObj.y[2]);
}
