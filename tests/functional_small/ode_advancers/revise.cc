
#include <gtest/gtest.h>
#include "pressio/ode_advancers.hpp"

using ScalarType = double;
using VectorType = std::vector<ScalarType>;

struct SimpleStepper
{
  using state_type = VectorType;
  using time_type = ScalarType;

  void operator()(state_type & odeState,
		  pressio::ode::TimeValue<time_type> currTime,
		  pressio::ode::StepCount step,
		  pressio::ode::TimeStepSize<time_type> dt)
  {
    auto timeval = currTime.get();
    for (std::size_t i=0; i<odeState.size(); i++){
      odeState[i] += timeval;
    }
  }
};

TEST(ode, explicit_advance_n_steps_fix_dt)
{
  VectorType odeState(5);
  std::for_each(odeState.begin(), odeState.end(),
		[](ScalarType & val){val = 0.0; });

  SimpleStepper stepper;
  ScalarType dt = 2.2;
  ScalarType t0 = 1.2;
  pressio::ode::advance_n_steps(stepper, odeState, t0, dt, pressio::ode::StepCount(4));
  std::for_each(odeState.begin(), odeState.end(),
		[](const ScalarType & val){ EXPECT_DOUBLE_EQ(val, 18.);});
}

////////////////////////////////////////////////////////////////

struct StepperForCheckingStepAndTime
{
  using state_type = VectorType;
  using time_type = ScalarType;

  void operator()(state_type & odeState,
		  pressio::ode::TimeValue<time_type> currTime,
		  pressio::ode::StepCount step,
		  pressio::ode::TimeStepSize<time_type> dt)
  {
    EXPECT_DOUBLE_EQ(dt.get(), 2.);
    if (step.get()==1){
      EXPECT_DOUBLE_EQ(currTime.get(), 0.);
    }
    if (step.get()==2){
      EXPECT_DOUBLE_EQ(currTime.get(), 2.);
    }
    if (step.get()==3){
      EXPECT_DOUBLE_EQ(currTime.get(), 4.);
    }
    if (step.get()==4){
      EXPECT_DOUBLE_EQ(currTime.get(), 6.);
    }
    if (step.get()==5){
      EXPECT_DOUBLE_EQ(currTime.get(), 8.);
    }
  }
};

TEST(ode, advance_n_steps_check_step_and_time)
{
  VectorType odeState(5);
  StepperForCheckingStepAndTime stepper;
  ScalarType dt = 2.;
  ScalarType t0 = 0.;
  pressio::ode::advance_n_steps(stepper, odeState, t0, dt, pressio::ode::StepCount(5));
}

////////////////////////////////////////////////////////////////

struct SimpleStepper2
{
  using state_type = VectorType;
  using time_type = ScalarType;

  void operator()(state_type & odeState,
		  pressio::ode::TimeValue<time_type> currTime,
		  pressio::ode::StepCount step,
		  pressio::ode::TimeStepSize<time_type> dt)
  {
    for (std::size_t i=0; i<odeState.size(); i++){
      odeState[i] += 1.;
    }
  }
};

class Observer2{
  public:
  using state_type = VectorType;
  template<class TimeType>
  void operator()(pressio::ode::StepCount step,
		  pressio::ode::TimeValue<TimeType> currtime,
		  const state_type & state) const
  {

    if (step.get()==0){
      EXPECT_DOUBLE_EQ(currtime.get(), 0.);
      for (std::size_t i=0; i<state.size(); i++){
	EXPECT_DOUBLE_EQ(state[i], 0.);
      }
    }

    if (step.get()==1){
      EXPECT_DOUBLE_EQ(currtime.get(), 2.);
      for (std::size_t i=0; i<state.size(); i++){
	EXPECT_DOUBLE_EQ(state[i], 1.);
      }
    }

    if (step.get()==2){
      EXPECT_DOUBLE_EQ(currtime.get(), 4.);
      for (std::size_t i=0; i<state.size(); i++){
	EXPECT_DOUBLE_EQ(state[i], 2.);
      }
    }

    if (step.get()==3){
      EXPECT_DOUBLE_EQ(currtime.get(), 6.);
      for (std::size_t i=0; i<state.size(); i++){
	EXPECT_DOUBLE_EQ(state[i], 3.);
      }
    }

  }
};

TEST(ode, explicit_advance_n_steps_fix_dt_check_collector_args)
{
  VectorType odeState(5);
  std::for_each(odeState.begin(), odeState.end(),
		[](ScalarType & val){val = 0.0; });

  Observer2 obs;
  SimpleStepper2 stepper;
  ScalarType dt = 2.;
  ScalarType t0 = 0.;
  pressio::ode::advance_n_steps(stepper, odeState, t0, dt, pressio::ode::StepCount(3), obs);
}

////////////////////////////////////////////////////////////////

// struct MyFakeStepper
// {
// 	using state_type = VectorType;
// 	using time_type = ScalarType;

//   template<typename solver_type>
//   void operator()(state_type & odeState,
// 		  const time_type & t,
// 		  const time_type & dt,
// 		  const typename pressio::ode::StepCount::value_type & step,
// 		  solver_type & solver)
//   {
//     for (int i=0; i<odeState.size(); i++){
//       odeState[i] += dt;
//     }
//   }
// };

// struct MyFakeSolver
// {
//   template<typename system_t, typename state_t>
//   void solve(const system_t & sys, state_t & state){}
// };

// TEST(ode, advance_n_steps_with_dt_setter_and_collector)
// {
//   /*
//     at step 0: [1,2,3]
//     at step 1: dt=2,   and leads to [3,4,5]
//     at step 2: dt=0.5, and leads to [3.5,4.5,5.5]

//    */

//   VectorType y(3);
//   y[0] = 1.0; y[1] = 2.0; y[2] = 3.0;

//   MyFakeStepper stepper;
//   MyFakeSolver solver;

//   std::string checkStr= "PASSED";

//   auto dtManager = [](const typename pressio::ode::StepCount::value_type & step,
// 		      const ScalarType & time,
// 		      ScalarType & dt)
// 		{
// 		  if(step==1) dt = 2.;
// 		  if(step==2) dt = 0.5;
// 		  if(step==3) dt = 1.;
// 		};

//   auto collector = [&checkStr](const typename pressio::ode::StepCount::value_type & step,
// 			       const ScalarType & time,
// 			       const VectorType & y)
// 		   {
// 		     using vec_t = std::vector<ScalarType>;
// 		     vec_t true0 = {1., 2., 3.};
// 		     vec_t true1 = {3., 4., 5.};
// 		     vec_t true2 = {3.5, 4.5, 5.5};
// 		     vec_t true3 = {4.5, 5.5, 6.5};

// 		     vec_t * trueY = nullptr;
// 		     if (step==0) trueY = &true0;
// 		     if (step==1) trueY = &true1;
// 		     if (step==2) trueY = &true2;
// 		     if (step==3) trueY = &true3;

// 		     for (auto i=0; i<3; ++i){
// 		       if( std::abs(y[i]-(*trueY)[i]) > 1e-11 )
// 			 checkStr = "FAILED";
// 		     }
// 		   };
//   pressio::ode::advance_n_steps_and_observe(stepper, y, 0., dtManager,
// 					    pressio::ode::StepCount(3),
// 					    collector, solver);
// }
