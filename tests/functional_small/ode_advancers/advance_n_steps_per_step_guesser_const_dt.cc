
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
    auto timeval = currTime.get();
    for (std::size_t i=0; i<odeState.size(); i++){
      odeState[i] += timeval;
    }
  }
};

struct Guesser{
  void operator()(pressio::ode::StepCount currStep,
                  pressio::ode::StepStartAt<double> currTime,
                  VectorType & state) const
  {
    if(currStep.get()==1) {
      state = {4., 2., 10.};
    }
    if(currStep.get()==2) {
      state = {3., 0., 1.};
    }
    if(currStep.get()==3) {
      state = {5., -1., -2.};
    }

  }
};

class Observer{
  public:
  using state_type = VectorType;
  template<class TimeType>
  void operator()(pressio::ode::StepCount step,
		  TimeType currtime,
		  const state_type & state) const
  {

    if (step.get()==0){
      for (std::size_t i=0; i<state.size(); i++){
	EXPECT_DOUBLE_EQ(state[i], 0.);
      }
    }

    if (step.get()==1){
      EXPECT_DOUBLE_EQ(state[0], 5.2);
      EXPECT_DOUBLE_EQ(state[1], 3.2);
      EXPECT_DOUBLE_EQ(state[2], 11.2);
    }

    if (step.get()==2){
      EXPECT_DOUBLE_EQ(state[0], 6.2);
      EXPECT_DOUBLE_EQ(state[1], 3.2);
      EXPECT_DOUBLE_EQ(state[2], 4.2);
    }

    if (step.get()==3){
      EXPECT_DOUBLE_EQ(state[0], 10.2);
      EXPECT_DOUBLE_EQ(state[1], 4.2);
      EXPECT_DOUBLE_EQ(state[2], 3.2);
    }
  }
};

TEST(ode, advance_n_steps_guesser_const_dt_stepper)
{
  VectorType odeState(3);
  std::for_each(odeState.begin(), odeState.end(),
		[](ScalarType & val){val = 0.0; });

  Stepper1 stepper;
  ScalarType dt = 2.0;
  ScalarType t0 = 1.2;
  pressio::ode::advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt,
						      pressio::ode::StepCount(3),
						      Guesser(), Observer());

  EXPECT_DOUBLE_EQ(odeState[0], 10.2);
  EXPECT_DOUBLE_EQ(odeState[1], 4.2);
  EXPECT_DOUBLE_EQ(odeState[2], 3.2);
}
