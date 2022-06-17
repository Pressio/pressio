
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
      pressio::ode::StepCount step,
      pressio::ode::StepSize<independent_variable_type> dt)
  {
    for (std::size_t i=0; i<odeState.size(); i++){
      odeState[i] += 1.;
    }
  }
};

class Observer1{
  public:
  using state_type = VectorType;
  template<class TimeType>
  void operator()(pressio::ode::StepCount step,
		  TimeType currtime,
		  const state_type & state) const
  {

    if (step.get()==0){
      EXPECT_DOUBLE_EQ(currtime, 0.);
      for (std::size_t i=0; i<state.size(); i++){
  EXPECT_DOUBLE_EQ(state[i], 0.);
      }
    }

    if (step.get()==1){
      EXPECT_DOUBLE_EQ(currtime, 2.);
      for (std::size_t i=0; i<state.size(); i++){
  EXPECT_DOUBLE_EQ(state[i], 1.);
      }
    }

    if (step.get()==2){
      EXPECT_DOUBLE_EQ(currtime, 4.);
      for (std::size_t i=0; i<state.size(); i++){
  EXPECT_DOUBLE_EQ(state[i], 2.);
      }
    }

    if (step.get()==3){
      EXPECT_DOUBLE_EQ(currtime, 6.);
      for (std::size_t i=0; i<state.size(); i++){
  EXPECT_DOUBLE_EQ(state[i], 3.);
      }
    }

  }
};

TEST(ode, advance_n_steps_const_dt_stepper_with_observer)
{
  VectorType odeState(5);
  std::for_each(odeState.begin(), odeState.end(), [](ScalarType & val){val = 0.0; });

  Observer1 obs;
  Stepper1 stepper;
  ScalarType dt = 2.;
  ScalarType t0 = 0.;
  pressio::ode::advance_n_steps(stepper, odeState, t0, dt, pressio::ode::StepCount(3), obs);
}
