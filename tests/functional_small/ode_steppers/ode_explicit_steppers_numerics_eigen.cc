
#include <gtest/gtest.h>
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "testing_apps.hpp"

//
// EULER 
//
#define TEST_ODE_EULER_NUMERICS(y, stepperObj) \
  y(0) = 1.; y(1) = 2.; y(2) = 3.; \
  double dt = 0.1; \
  ode::advance_n_steps(stepperObj, y, 0.0, dt, 1);\
  EXPECT_DOUBLE_EQ( y(0), 1.1); \
  EXPECT_DOUBLE_EQ( y(1), 2.2); \
  EXPECT_DOUBLE_EQ( y(2), 3.3); \
  std::cout << std::setprecision(14) << y; \
  ode::advance_n_steps(stepperObj, y, 0.0, dt, 1); \
  EXPECT_DOUBLE_EQ( y(0), 1.21); \
  EXPECT_DOUBLE_EQ( y(1), 2.42); \
  EXPECT_DOUBLE_EQ( y(2), 3.63); \

TEST(ode, explicit_euler_system_reference)
{
  using namespace pressio;
  using app_t = ode::testing::refAppEigen;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  auto stepperObj = ode::create_forward_euler_stepper(y, appObj);
  TEST_ODE_EULER_NUMERICS(y, stepperObj);
}

TEST(ode, explicit_euler_system_move)
{
  using namespace pressio;
  using app_t = ode::testing::refAppEigen;
  using state_t = typename app_t::state_type;
  app_t appObj;

  state_t y(3);
  auto stepperObj = ode::create_forward_euler_stepper(y, std::move(appObj));
  TEST_ODE_EULER_NUMERICS(y, stepperObj);
}

//
// RK4
//
#define TEST_ODE_RK4_NUMERICS(y, stepperObj, appObj) \
  y(0) = 1.; y(1) = 2.; y(2) = 3.; \
  double dt = 0.1; \
  ode::advance_n_steps(stepperObj, y, 0.0, dt, 1); \
  std::cout << std::setprecision(14) << y; \
  appObj.analyticAdvanceRK4(dt); \
  EXPECT_DOUBLE_EQ(y(0), appObj.y(0)); \
  EXPECT_DOUBLE_EQ(y(1), appObj.y(1)); \
  EXPECT_DOUBLE_EQ(y(2), appObj.y(2)); \

TEST(ode, explicit_rk4_system_reference)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  auto stepperObj = ode::create_runge_kutta4_stepper(y,appObj);
  TEST_ODE_RK4_NUMERICS(y, stepperObj, appObj);
}

TEST(ode, explicit_rk4_system_move)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  auto stepperObj = ode::create_runge_kutta4_stepper(y,app_t());
  TEST_ODE_RK4_NUMERICS(y, stepperObj, appObj);
}

//
//SSPRK3
//
struct AppForSSPRK3
{
  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using velocity_type = state_type;

  void velocity(const state_type & y,
                const scalar_type time, 
                velocity_type & rhs) const
  {
    auto sz = y.size();
    for (decltype(sz) i=0; i<sz; i++){
      rhs[i] = y[i] + time;
    }
  };

  velocity_type createVelocity() const{
    return velocity_type(3);
  };
};
TEST(ode, explicit_ssprk3)
{
  using namespace pressio;
  using app_t = AppForSSPRK3;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.; 
  auto stepperObj = ode::create_ssp_runge_kutta3_stepper(y,appObj);
  double dt = 2.; 
  ode::advance_n_steps(stepperObj, y, 0.0, dt, 1);
  EXPECT_DOUBLE_EQ( y(0), 29./3.); 
  EXPECT_DOUBLE_EQ( y(1), 48./3.); 
  EXPECT_DOUBLE_EQ( y(2), 67./3.);
}


namespace {
struct AB2MyApp
{
  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using velocity_type = state_type;

  void velocity(const state_type & y,
    scalar_type t,
    velocity_type & f) const
  {
    f.setConstant(t);
  };

  velocity_type createVelocity() const
  {
    velocity_type R(3);
    return R;
  };
};

struct Collector
{
  void operator()(const ::pressio::ode::step_count_type & step,
      const double & time,
      const Eigen::VectorXd & y)
  {
    if (step==0){
      EXPECT_DOUBLE_EQ( y(0), 1.);
      EXPECT_DOUBLE_EQ( y(1), 2.);
      EXPECT_DOUBLE_EQ( y(2), 3.);
    }

    if (step==1){
      EXPECT_DOUBLE_EQ( y(0), 1.);
      EXPECT_DOUBLE_EQ( y(1), 2.);
      EXPECT_DOUBLE_EQ( y(2), 3.);
    }

    if (step==2){
      EXPECT_DOUBLE_EQ( y(0), 7.);
      EXPECT_DOUBLE_EQ( y(1), 8.);
      EXPECT_DOUBLE_EQ( y(2), 9.);
    }

    if (step==3){
      EXPECT_DOUBLE_EQ( y(0), 17.);
      EXPECT_DOUBLE_EQ( y(1), 18.);
      EXPECT_DOUBLE_EQ( y(2), 19.);
    }
  }
};
}

TEST(ode, explicit_ab2_system_reference)
{
  /*
    dy/dt = f
    where f returns [t t t]

    for AB2, step1 from t_0 -> t_1 is euler:
    y_1 = [1 2 3] + 2.*f(t0)
    f0 = [0 0 0]

    step2: from t_1 -> t_2
    y_2 = y_1 + dt*[ (3/2)*f(t_1) - (1/2)*f(t_0)]
        = [1 2 3] + 2*(3/2)*[2 2 2]
  = [7 8 9]

    step3:
    y_3 = y_2 + dt*[ (3/2)*f(t_2) - (1/2)*f(t_1)]
        = [7 8 9] + 2*[ (3/2)*[4 4 4] - (1/2)[2 2 2] ]
  = [17 18 19]
   */

  using namespace pressio;
  using app_t = AB2MyApp;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  auto stepperObj = ode::create_adams_bashforth2_stepper(y,appObj);

  double dt = 2.;
  Collector C;
  ode::advance_n_steps_and_observe(stepperObj, y, 0.0, dt, 3, C);
}

TEST(ode, explicit_ab2_custom_system_move)
{
  using namespace pressio;
  using app_t = AB2MyApp;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  auto stepperObj = ode::create_adams_bashforth2_stepper(y,std::move(appObj));
  double dt = 2.;
  Collector C;
  ode::advance_n_steps_and_observe(stepperObj, y, 0.0, dt, 3, C);
}