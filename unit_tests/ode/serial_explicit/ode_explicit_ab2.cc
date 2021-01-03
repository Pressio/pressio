
#include <gtest/gtest.h>
#include "pressio_ode.hpp"

struct MyApp
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
  void operator()(const ::pressio::ode::types::step_t & step,
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

TEST(ode_explicit_ab2, numericsStdResidualPolDefaultCreated)
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
  using app_t = MyApp;
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

TEST(ode_explicit_ab2, numericsStdResidualPolPassedByUser)
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
  using app_t = MyApp;
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
