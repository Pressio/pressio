
#include <gtest/gtest.h>
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "custom_independent_variable.hpp"

struct MyApp{
  using independent_variable_type = MyCustomTime;
  using state_type = Eigen::VectorXd;
  using rhs_type = state_type;

  state_type createState() const{
    state_type ret(3);
    ret.setZero();
    return ret;
  }

  rhs_type createRhs() const{
    rhs_type ret(3);
    ret.setZero();
    return ret;
  };

  void rhs(const state_type & y,
	   independent_variable_type /*unused*/,
	   rhs_type & f) const
  {
    f[0] = -10. * y[0];
    f[1] = -10. * y[1];
    f[2] = -10. * y[2];
  };
};

TEST(ode_explicit_steppers, euler_custom_ind_var_type)
{
  using namespace pressio;
  using app_t = MyApp;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  auto stepperObj = ode::create_forward_euler_stepper(appObj);

  y(0) = 1.; y(1) = 2.; y(2) = 3.;
  MyCustomTime t0{0.0};
  MyCustomTime dt{0.1};
  ode::advance_n_steps(stepperObj, y, t0,
		       dt, pressio::ode::StepCount(1));
  EXPECT_DOUBLE_EQ( y[0], 0.);
  EXPECT_DOUBLE_EQ( y[1], 0.);
  EXPECT_DOUBLE_EQ( y[2], 0.);
}
