
#include <gtest/gtest.h>
#include "pressio_ode_explicit.hpp"
#include "../reference_apps_for_testing.hpp"

TEST(ode_explicit_euler, numericsStdResidualPolDefaultCreated){
  using namespace pressio;
  using app_t = ode::testing::refAppEigen;
  using nstate_t = typename app_t::state_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  // using stepper_t = ode::ExplicitStepper<
  //   ode::explicitmethods::Euler, state_t, app_t>;
  // stepper_t stepperObj(y, appObj);
  auto stepperObj = ode::createForwardEulerStepper(y, appObj);

  // integrate in time
  double dt = 0.1;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 1);
  EXPECT_DOUBLE_EQ( y(0), 1.1);
  EXPECT_DOUBLE_EQ( y(1), 2.2);
  EXPECT_DOUBLE_EQ( y(2), 3.3);
  std::cout << std::setprecision(14) << *y.data();

  // integrate in time
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 1);
  EXPECT_DOUBLE_EQ( y(0), 1.21);
  EXPECT_DOUBLE_EQ( y(1), 2.42);
  EXPECT_DOUBLE_EQ( y(2), 3.63);
}

TEST(ode_explicit_euler, numericsStdResidualPolPassedByUser){
  using namespace pressio;
  using app_t = ode::testing::refAppEigen;
  using nstate_t = typename app_t::state_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  // the standard policy
  using res_std_pol_t = ode::explicitmethods::policy::VelocityStandardPolicy<state_t>;
  res_std_pol_t polObj;

  // using stepper_t = ode::ExplicitStepper<
  //   ode::explicitmethods::Euler, state_t, app_t, res_std_pol_t>;
  // stepper_t stepperObj(y, appObj, polObj);
  auto stepperObj = ode::createForwardEulerStepper(y, appObj, polObj);

  // integrate in time
  double dt = 0.1;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 1);
  EXPECT_DOUBLE_EQ( y(0), 1.1);
  EXPECT_DOUBLE_EQ( y(1), 2.2);
  EXPECT_DOUBLE_EQ( y(2), 3.3);
  std::cout << std::setprecision(14) << *y.data();

  // integrate in time
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 1);
  EXPECT_DOUBLE_EQ( y(0), 1.21);
  EXPECT_DOUBLE_EQ( y(1), 2.42);
  EXPECT_DOUBLE_EQ( y(2), 3.63);
}
