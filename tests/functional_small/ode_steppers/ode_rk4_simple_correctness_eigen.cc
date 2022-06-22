
#include <gtest/gtest.h>
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "testing_apps.hpp"

#define TEST_ODE_RK4_NUMERICS(y, stepperObj, appObj) \
  y(0) = 1.; y(1) = 2.; y(2) = 3.; \
  double dt = 0.1;				    \
  ode::advance_n_steps(stepperObj, y, 0.0, dt, pressio::ode::StepCount(1));  \
  appObj.analyticAdvanceRK4(dt);		    \
  EXPECT_DOUBLE_EQ(y(0), appObj.y(0));		    \
  EXPECT_DOUBLE_EQ(y(1), appObj.y(1));		    \
  EXPECT_DOUBLE_EQ(y(2), appObj.y(2));		    \

TEST(ode_explicit_steppers, rk4_system_reference)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  auto stepperObj = ode::create_runge_kutta4_stepper(appObj);
  TEST_ODE_RK4_NUMERICS(y, stepperObj, appObj);
}

TEST(ode_explicit_steppers, rk4_system_move)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  auto stepperObj = ode::create_runge_kutta4_stepper(app_t());
  TEST_ODE_RK4_NUMERICS(y, stepperObj, appObj);
}
