
#include <gtest/gtest.h>
#include "pressio_ode_explicit.hpp"
#include "testing_apps.hpp"

TEST(ode, explicit_euler_numericsStdResidualPolDefaultCreated){
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

TEST(ode, explicit_euler_numericsStdResidualPolPassedByUser){
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



TEST(ode, explicit_rk4_numericsStdResidualPolDefaultCreated)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  // using nveloc_t = typename app_t::velocity_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  // using stepper_t = ode::ExplicitStepper<
  //   ode::explicitmethods::RungeKutta4, state_t, app_t>;
  // stepper_t stepperObj(y, appObj);
  auto stepperObj = ode::createRungeKutta4Stepper(y, appObj);


  // // integrate in time
  double dt = 0.1;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 1);
  std::cout << std::setprecision(14) << *y.data();

  appObj.analyticAdvanceRK4(dt);

  EXPECT_DOUBLE_EQ(y(0), appObj.y(0));
  EXPECT_DOUBLE_EQ(y(1), appObj.y(1));
  EXPECT_DOUBLE_EQ(y(2), appObj.y(2));
}


TEST(ode, explicit_rk4_numericsStdResidualPolPassedByUser)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  // using nveloc_t = typename app_t::velocity_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  // using res_t = containers::Vector<nveloc_t>;

  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  // // the standard policy
  using res_std_pol_t = ::pressio::ode::explicitmethods::policy::VelocityStandardPolicy<state_t>;
  res_std_pol_t polObj;
  // using stepper_t = ode::ExplicitStepper<
  //   ode::explicitmethods::RungeKutta4, state_t,
  //   app_t, res_t, res_std_pol_t>;
  // stepper_t stepperObj(y, appObj, polObj);
  auto stepperObj = ode::createRungeKutta4Stepper(y, appObj, polObj);


  // integrate in time
  double dt = 0.1;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 1);
  std::cout << std::setprecision(14) << *y.data();

  appObj.analyticAdvanceRK4(dt);

  EXPECT_DOUBLE_EQ(y(0), appObj.y(0));
  EXPECT_DOUBLE_EQ(y(1), appObj.y(1));
  EXPECT_DOUBLE_EQ(y(2), appObj.y(2));
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
  void operator()(const ::pressio::ode::step_type & step,
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

TEST(ode, explicit_ab2_numericsStdResidualPolDefaultCreated)
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
  using nstate_t = typename app_t::state_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  auto stepperObj = ode::createAdamsBashforth2Stepper(y, appObj);

  double dt = 2.;
  Collector C;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 3, C);
}

TEST(ode, explicit_ab2numericsStdResidualPolPassedByUser)
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
  using nstate_t = typename app_t::state_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  using pol_t = ode::explicitmethods::policy::VelocityStandardPolicy<state_t>;
  pol_t polObj;
  auto stepperObj = ode::createAdamsBashforth2Stepper(y, appObj, polObj);

  double dt = 2.;
  Collector C;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 3, C);
}
