
#include <gtest/gtest.h>
#include "pressio/ode_advancers.hpp"

using ScalarType = double;
using VectorType = std::vector<ScalarType>;

struct Stepper1
{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;

  void operator()(state_type & odeState,
		  pressio::ode::StepStartAt<independent_variable_type> currTime,
		  pressio::ode::StepCount /*unused*/,
		  pressio::ode::StepSize<independent_variable_type> /*unused*/)
  {
    auto timeval = currTime.get();
    for (std::size_t i=0; i<odeState.size(); i++){
      odeState[i] += timeval;
    }
  }
};

TEST(ode, advance_n_steps_const_dt_stepper)
{
  VectorType odeState(5);
  std::for_each(odeState.begin(), odeState.end(),
		[](ScalarType & val){val = 0.0; });

  Stepper1 stepper;
  ScalarType dt = 2.2;
  ScalarType t0 = 1.2;
  pressio::ode::advance_n_steps(stepper, odeState, t0, dt, pressio::ode::StepCount(4));
  std::for_each(odeState.begin(), odeState.end(),
		[](const ScalarType & val){ EXPECT_DOUBLE_EQ(val, 18.);});
}
