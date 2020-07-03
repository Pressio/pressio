
#include <gtest/gtest.h>
#include "pressio_ode.hpp"
#include "../reference_apps_for_testing.hpp"



TEST(ode_implicit_bdf2, numericsStdPoliciesDefaultCreated){
  using namespace pressio;

  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  using nveloc_t = typename app_t::velocity_type;
  using njacobian_t = typename app_t::jacobian_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  using res_t = containers::Vector<nveloc_t>;
  using jac_t = containers::Matrix<njacobian_t>;
  state_t y(3);
  *y.data() = appObj.getInitCond();

  // define auxiliary stepper
  using aux_stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::Euler, state_t, res_t, jac_t, app_t>;
  aux_stepper_t stepperAux(y, appObj);

  // bdf2 stepper
  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF2, state_t, res_t, jac_t, app_t, aux_stepper_t>;
  stepper_t stepperObj(y, appObj, stepperAux);

  // define solver
  using lin_solver_t = solvers::linear::Solver<
    solvers::linear::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;
  using nonlinear_solver_t = pressio::solvers::NewtonRaphson<stepper_t, lin_solver_t, double>;
  nonlinear_solver_t solverO(stepperObj, y, linSolverObj);

  // integrate in time
  ::pressio::ode::types::step_t nSteps = 4;
  double dt = 0.01;
  ode::integrateNSteps(stepperObj, y, 0.0, dt, nSteps, solverO);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  appObj.analyticAdvanceBackEulerNSteps(dt, 1);
  appObj.analyticAdvanceBDF2NSteps(dt, 3);
  std::cout << std::setprecision(14) << appObj.y << "\n";
  EXPECT_DOUBLE_EQ(y[0], appObj.y[0]);
  EXPECT_DOUBLE_EQ(y[1], appObj.y[1]);
  EXPECT_DOUBLE_EQ(y[2], appObj.y[2]);
}


TEST(ode_implicit_bdf2, numericsStdResidualPolPassedByUser){
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  using nveloc_t = typename app_t::velocity_type;
  using njacobian_t = typename app_t::jacobian_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  using res_t = containers::Vector<nveloc_t>;
  using jac_t = containers::Matrix<njacobian_t>;
  state_t y(3);
  *y.data() = appObj.getInitCond();

  // define auxiliary policies and stepper
  using res_pol_t
    = ode::implicitmethods::policy::ResidualStandardPolicy<state_t, app_t, res_t>;
  using jac_pol_t
    = ode::implicitmethods::policy::JacobianStandardPolicy<state_t, app_t, jac_t>;

  using aux_stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::Euler,
    state_t, res_t, jac_t, app_t, res_pol_t, jac_pol_t>;
  aux_stepper_t stepperAux(y, appObj, res_pol_t(), jac_pol_t());

  // stepper for BDF2
  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF2,
    state_t, res_t, jac_t, app_t, aux_stepper_t, res_pol_t, jac_pol_t>;
  stepper_t stepperObj(y, appObj, res_pol_t(), jac_pol_t(), stepperAux);

  // define solver
  using lin_solver_t = solvers::linear::Solver<
    solvers::linear::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;
  using nonlinear_solver_t = pressio::solvers::NewtonRaphson<stepper_t, lin_solver_t, double>;
  nonlinear_solver_t solverO(stepperObj, y, linSolverObj);

  // integrate in time
  ::pressio::ode::types::step_t nSteps = 4;
  double dt = 0.01;
  ode::integrateNSteps(stepperObj, y, 0.0, dt, nSteps, solverO);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  appObj.analyticAdvanceBackEulerNSteps(dt, 1);
  appObj.analyticAdvanceBDF2NSteps(dt, 3);
  std::cout << std::setprecision(14) << appObj.y << "\n";
  EXPECT_DOUBLE_EQ(y[0], appObj.y[0]);
  EXPECT_DOUBLE_EQ(y[1], appObj.y[1]);
  EXPECT_DOUBLE_EQ(y[2], appObj.y[2]);
}


// TEST(ode_implicit_bdf2, traits){
//   using namespace pressio;

//   using app_t = ode::testing::refAppForImpEigen;
//   using nstate_t = typename app_t::state_type;
//   using nveloc_t = typename app_t::velocity_type;
//   using njac_t = typename app_t::jacobian_type;
//   using state_t = containers::Vector<nstate_t>;
//   using res_t = containers::Vector<nveloc_t>;
//   using jac_t = containers::Matrix<njac_t>;

//   static_assert(
//     ode::meta::is_legitimate_model_for_explicit_ode<app_t>::value, "");
//   static_assert(
//     ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");

//   using aux_stepper_t = ode::ImplicitStepper<
//     ode::implicitmethods::Euler,
//     state_t, res_t, jac_t, app_t, void>; /*aux stepper NOT needed for backEuler*/

//   using stepper_t = ode::ImplicitStepper<
//     ode::implicitmethods::BDF2,
//     state_t, res_t, jac_t, app_t, aux_stepper_t>;

//   using traits = ode::details::traits<stepper_t>;

//   // static_assert(ode::meta::is_implicit_bdf2_residual_standard_policy<
//   //     typename traits::residual_policy_t>::value,
//   //     "");
//   // static_assert(ode::meta::is_implicit_bdf2_jacobian_standard_policy<
//   //     typename traits::jacobian_policy_t>::value,
//   //     "");

//   ::testing::StaticAssertTypeEq<typename
//          traits::aux_stepper_t, aux_stepper_t>();
//   ::testing::StaticAssertTypeEq<typename
//          traits::state_t, state_t>();
//   ::testing::StaticAssertTypeEq<typename
//          traits::residual_t,res_t>();
//   ::testing::StaticAssertTypeEq<typename
//          traits::jacobian_t,jac_t>();
//   ::testing::StaticAssertTypeEq<typename
//          traits::scalar_t,double>();
//   ::testing::StaticAssertTypeEq<typename
//          traits::system_t,app_t>();
//   static_assert( traits::order_value == 2, "");
// }


// TEST(ode_implicit_bdf2, numericsUserResidualDefaultJac){
//   using namespace pressio;
//   using app_t		= ode::testing::refAppForImpEigen;
//   using nstate_t	= typename app_t::state_type;
//   using nveloc_t	= typename app_t::velocity_type;
//   using njacobian_t	= typename app_t::jacobian_type;
//   app_t appObj;

//   using state_t = containers::Vector<nstate_t>;
//   using res_t = containers::Vector<nveloc_t>;
//   using jac_t = containers::Matrix<njacobian_t>;
//   state_t y(3);
//   *y.data() = appObj.getInitCond();

//   //**********************
//   // residual policy is user-defined (even if it is standard)
//   // jacobian_policy is standard and default-constructed

//   using res_pol_t = ode::implicitmethods::policy::ResidualStandardPolicy<state_t, app_t, res_t>;

//   using aux_stepper_t = ode::ImplicitStepper<
//     ode::implicitmethods::Euler,
//     state_t, res_t, jac_t, app_t,
//     void, /*aux stepper NOT needed for backEuler*/
//     res_pol_t>;
//   aux_stepper_t stepperAux(y, appObj, res_pol_t());

//   // stepper for BDF2
//   using stepper_t = ode::ImplicitStepper<
//     ode::implicitmethods::BDF2,
//     state_t, res_t, jac_t, app_t, aux_stepper_t, res_pol_t>;
//   stepper_t stepperObj(y, appObj, res_pol_t(), stepperAux);

//   // define solver
//   using lin_solver_t = solvers::linear::Solver<
//     solvers::linear::iterative::Bicgstab, jac_t>;
//   lin_solver_t linSolverObj;
//   using nonlinear_solver_t = pressio::solvers::NewtonRaphson<stepper_t, lin_solver_t, double>;
//   nonlinear_solver_t solverO(stepperObj, y, linSolverObj);

//   // integrate in time
//   ::pressio::ode::types::step_t nSteps = 4;
//   double dt = 0.01;
//   ode::integrateNSteps(stepperObj, y, 0.0, dt, nSteps, solverO);
//   std::cout << std::setprecision(14) << *y.data() << "\n";

//   appObj.analyticAdvanceBackEulerNSteps(dt, 1);
//   appObj.analyticAdvanceBDF2NSteps(dt, 3);
//   std::cout << std::setprecision(14) << appObj.y << "\n";
//   EXPECT_DOUBLE_EQ(y[0], appObj.y[0]);
//   EXPECT_DOUBLE_EQ(y[1], appObj.y[1]);
//   EXPECT_DOUBLE_EQ(y[2], appObj.y[2]);
// }
