
#include <gtest/gtest.h>
#include "pressio/ode_advancers.hpp"

using ScalarType = double;
using VectorType = std::array<ScalarType,3>;

struct Stepper1{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    for (auto & it : state){ it += 1. + dt.get(); }
  }
};

struct AuxClass {
  ScalarType value() { return -8.0; }
  ScalarType value() const { return -9.0; }
};

struct AuxClass2{
  ScalarType value() { return -10.0; }
  ScalarType value()const { return -12.0; }
};

struct Stepper2a{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  AuxClass && aux)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    for (auto & it : state){ it += 2.1 + aux.value() + dt.get(); }
  }
};

struct Stepper2b{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  AuxClass & aux)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    for (auto & it : state){ it += 2.2 + aux.value() + dt.get(); }
  }
};

struct Stepper2c{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  const AuxClass & aux)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    for (auto & it : state){ it += 2.3 + aux.value() + dt.get(); }
  }
};

struct Stepper3a{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  AuxClass && aux,
		  AuxClass2 && aux2)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    for (auto & it : state){ it += 3.1 + aux.value() + aux2.value() + dt.get(); }
  }
};

struct Stepper3b{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  AuxClass && aux,
		  AuxClass2 & aux2)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    for (auto & it : state){ it += 3.2 + aux.value() + aux2.value() + dt.get(); }
  }
};

struct Stepper3c{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  AuxClass && aux,
		  const AuxClass2 & aux2)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    for (auto & it : state){ it += 3.3 + aux.value() + aux2.value() + dt.get(); }
  }
};

struct Guesser{
  void operator()(const pressio::ode::StepCount & /*unused*/,
		  pressio::ode::StepStartAt<double> tn,
		  VectorType & state) const
  { state.fill(10.0); }
};

struct Observer{
  void operator()(const pressio::ode::StepCount & /*unused*/,
		  const ScalarType & /*unused*/,
		  const VectorType & /*unused*/) const{}
};

struct DtSetter{
  void operator()(const pressio::ode::StepCount & /*unused*/,
                  const pressio::ode::StepStartAt<double> & /*unused*/,
                  pressio::ode::StepSize<double> & dt) const
  {
    dt = 4.0;
  }
};

void check_state_and_reset(VectorType & a, ScalarType value){
  EXPECT_DOUBLE_EQ(a[0], value);
  EXPECT_DOUBLE_EQ(a[1], value);
  EXPECT_DOUBLE_EQ(a[2], value);
  a.fill(0.);
}

TEST(ode, test)
{
  VectorType odeState; odeState.fill(0.);
  ScalarType dt = 2.0;
  ScalarType t0 = 1.0;
  const auto numSteps = pressio::ode::StepCount(1);
  Observer obs;
  DtSetter dtPol;
  Guesser guess;

  using namespace pressio::ode;

  {
    Stepper1 stepper;

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess);
    check_state_and_reset(odeState, 13.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess);
    check_state_and_reset(odeState, 15.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess, obs);
    check_state_and_reset(odeState, 13.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess, obs);
    check_state_and_reset(odeState, 15.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0,
					  [](const StepCount & /*unused*/,
					     const StepStartAt<double> & /*unused*/,
					     StepSize<double> & dt){ dt = 5.0; },
					  numSteps, guess);
    check_state_and_reset(odeState, 16.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess,
					  [](const StepCount & /*unused*/, ScalarType /*unused*/,
					     const VectorType & /*unused*/){});
    check_state_and_reset(odeState, 13.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0,
					  [](const StepCount & /*unused*/,
					     const StepStartAt<double> & /*unused*/,
					     StepSize<double> & dt){ dt=5.0; },
					  numSteps,
					  guess,
					  [](const StepCount & /*unused*/, ScalarType /*unused*/,
					     const VectorType & /*unused*/){});
    check_state_and_reset(odeState, 16.0);
  }

  {
    Stepper2a stepper;

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess, AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+2.0+10.);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess, AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+4.0+10.);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess, obs, AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+dt+10.);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess, obs, AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+4.0+10.);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0,
					  [](const StepCount & /*unused*/,
					     const StepStartAt<double> & /*unused*/,
					     StepSize<double> & dt){dt = 6.0;},
					  numSteps, guess,AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+6.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess,
		    [](const StepCount & /*unused*/, ScalarType /*unused*/,
		       const VectorType & /*unused*/){},
		    AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+2.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0,
					  [](const StepCount & /*unused*/,
					     const StepStartAt<double> & /*unused*/,
					     StepSize<double> & dt){dt=6.0;},
					  numSteps, guess,
					  [](const StepCount & /*unused*/, ScalarType /*unused*/,
					     const VectorType & /*unused*/){},
					  AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+6.0+10.0);
  }

  {
    Stepper2b stepper; AuxClass aux;

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess, aux);
    check_state_and_reset(odeState, -8.0+2.2+dt+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess, aux);
    check_state_and_reset(odeState, -8.0+2.2+4.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess,obs, aux);
    check_state_and_reset(odeState, -8.0+2.2+dt+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess,obs, aux);
    check_state_and_reset(odeState, -8.0+2.2+4.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0,
					  [](const StepCount & /*unused*/,
					     const StepStartAt<double> & /*unused*/,
					     StepSize<double> & dt){dt=5.0;},
					  numSteps, guess,aux);
    check_state_and_reset(odeState, -8.0+2.2+5.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess,
					  [](const StepCount & /*unused*/, ScalarType /*unused*/,
					     const VectorType & /*unused*/){},
					  aux);
    check_state_and_reset(odeState, -8.0+2.2+dt+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0,
					  [](const StepCount & /*unused*/,
					     const StepStartAt<double> & /*unused*/,
					     StepSize<double> & dt){dt=6.0;},
					  numSteps, guess,
					  [](const StepCount & /*unused*/, ScalarType /*unused*/,
					     const VectorType & /*unused*/){},
					  aux);
    check_state_and_reset(odeState, -8.0+2.2+6.0+10.0);
  }

  {
    Stepper2c stepper; const AuxClass aux;

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess,aux);
    check_state_and_reset(odeState, -9.0+2.3+dt+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess,aux);
    check_state_and_reset(odeState, -9.0+2.3+4.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess,obs, aux);
    check_state_and_reset(odeState, -9.0+2.3+dt+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess,obs, aux);
    check_state_and_reset(odeState, -9.0+2.3+4.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0,
					  [](const StepCount & /*unused*/,
					     const StepStartAt<double> & /*unused*/,
					     StepSize<double> & dt){dt=5.0;},
					  numSteps, guess,aux);
    check_state_and_reset(odeState, -9.0+2.3+5.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess,
					  [](const StepCount & /*unused*/,
					     ScalarType /*unused*/,
					     const VectorType & /*unused*/){},
					  aux);
    check_state_and_reset(odeState, -9.0+2.3+dt+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0,
					  [](const StepCount & /*unused*/,
					     const StepStartAt<double> & /*unused*/,
					     StepSize<double> & dt){dt=6.0;},
					  numSteps, guess,
					  [](const StepCount & /*unused*/,
					     ScalarType /*unused*/,
					     const VectorType & /*unused*/){},
					  aux);
    check_state_and_reset(odeState, -9.0+2.3+6.0+10.0);
  }

  {
    Stepper3a stepper; AuxClass aux; AuxClass2 aux2;

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess, AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+2.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess, AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+4.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess, obs, AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+2.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess, obs, AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+4.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess,
					  [](const StepCount & /*unused*/, ScalarType /*unused*/,
					     const VectorType & /*unused*/){},
					  AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+2.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0,
					  [](const StepCount & /*unused*/,
					     const StepStartAt<double> & /*unused*/,
					     StepSize<double> & dt){dt=5.0;},
					  numSteps, guess,
					  [](const StepCount & /*unused*/,
					     ScalarType /*unused*/,
					     const VectorType & /*unused*/){},
					  AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+5.0+10.0);
  }

  {
    Stepper3b stepper; AuxClass aux; AuxClass2 aux2;

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess, AuxClass{}, aux2);
    check_state_and_reset(odeState, -18.0+3.2+2.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess, AuxClass{}, aux2);
    check_state_and_reset(odeState, -18.0+3.2+4.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess, obs, AuxClass{}, aux2);
    check_state_and_reset(odeState, -18.0+3.2+2.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess, obs, AuxClass{}, aux2);
    check_state_and_reset(odeState, -18.0+3.2+4.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0,
		    [](const StepCount & /*unused*/, const StepStartAt<double> & /*unused*/,
		       StepSize<double> & dt){dt=0.0;},
					  numSteps, guess, AuxClass{}, aux2);
    check_state_and_reset(odeState, -18.0+3.2+10.0);
  }

  {
    Stepper3c stepper; AuxClass aux; const AuxClass2 aux2;

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess, AuxClass{}, aux2);
    check_state_and_reset(odeState, -20.0+3.3+2.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess, AuxClass{}, aux2);
    check_state_and_reset(odeState, -20.0+3.3+4.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dt, numSteps, guess, obs, AuxClass{}, aux2);
    check_state_and_reset(odeState, -20.0+3.3+2.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0, dtPol, numSteps, guess, obs, AuxClass{}, aux2);
    check_state_and_reset(odeState, -20.0+3.3+4.0+10.0);

    advance_n_steps_with_pre_step_guesser(stepper, odeState, t0,
		    [](const StepCount & /*unused*/, const StepStartAt<double> & /*unused*/,
		       StepSize<double> & dt){dt=0.0;},
					  numSteps, guess, AuxClass{}, aux2);
    check_state_and_reset(odeState, -20.0+3.3+10.0);
  }
}
