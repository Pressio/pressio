#include <gtest/gtest.h>

#include "custom_types_and_ops_for_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressio/ode_steppers_explicit.hpp"

struct MyApp{
  using independent_variable_type   = ScalarType;
  using state_type    = VectorType;
  using right_hand_side_type = VectorType;

public:
  state_type createState() const{ return state_type(3); }

  void operator()(const state_type & y,
		  independent_variable_type /*unused*/,
		  right_hand_side_type & R) const
  {
    R[0] = -10. * y[0];
    R[1] = -10. * y[1];
    R[2] = -10. * y[2];
  };

  right_hand_side_type createRightHandSide() const{
    right_hand_side_type R(3);
    return R;
  };
};

TEST(ode, explicit_rk4_custom_types)
{
  using namespace pressio;
  using app_t     = MyApp;
  using state_t    = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  y[0] = 1.; y[1] = 2.; y[2] = 3.;

  auto stepperObj = ode::create_rk4_stepper(appObj);

  ScalarType dt = 0.1;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, ::pressio::ode::StepCount(1));
  EXPECT_DOUBLE_EQ( y[0], 0.375);
  EXPECT_DOUBLE_EQ( y[1], 0.75);
  EXPECT_DOUBLE_EQ( y[2], 1.125);
}
