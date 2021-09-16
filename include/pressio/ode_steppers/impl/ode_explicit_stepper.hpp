/*
//@HEADER
// ************************************************************************
//
// ode_explicit_euler_stepper_impl.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef ODE_EXPLICIT_IMPL_ODE_EXPLICIT_STP_IMPL_HPP_
#define ODE_EXPLICIT_IMPL_ODE_EXPLICIT_STP_IMPL_HPP_

#include <array>

namespace pressio{ namespace ode{ namespace impl{

template<
  class ScalarType,
  class StateType,
  class SystemType,
  class VelocityType
  >
class ExplicitStepper
{

public:
  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;


private:
  StepScheme name_;
  const stepper_order_type order_;
  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
  std::vector<VelocityType> velocities_;
  StateType auxiliaryState_;

public:
  ExplicitStepper() = delete;
  ExplicitStepper(const ExplicitStepper &) = default;
  ExplicitStepper & operator=(const ExplicitStepper &) = delete;
  ExplicitStepper(ExplicitStepper &&) = default;
  ExplicitStepper & operator=(ExplicitStepper &&) = delete;
  ~ExplicitStepper() = default;

  ExplicitStepper(ode::ForwardEuler,
		  const StateType & state,
		  SystemType && systemObj)
    : name_(StepScheme::ForwardEuler),
      order_(1),
      systemObj_(std::forward<SystemType>(systemObj)),
      velocities_{systemObj.createVelocity()},
      auxiliaryState_{::pressio::ops::clone(state)}
  {}

  ExplicitStepper(ode::RungeKutta4,
		  const StateType & state,
		  SystemType && systemObj)
    : name_(StepScheme::RungeKutta4),
      order_(4),
      systemObj_(std::forward<SystemType>(systemObj)),
      velocities_{systemObj.createVelocity(),
		  systemObj.createVelocity(),
		  systemObj.createVelocity(),
		  systemObj.createVelocity()},
      auxiliaryState_{::pressio::ops::clone(state)}
  {}

  ExplicitStepper(ode::AdamsBashforth2,
		  const StateType & state,
		  SystemType && systemObj)
    : name_(StepScheme::AdamsBashforth2),
      order_(2),
      systemObj_(std::forward<SystemType>(systemObj)),
      velocities_{systemObj.createVelocity(), systemObj.createVelocity()},
      auxiliaryState_{::pressio::ops::clone(state)}
  {}

  ExplicitStepper(ode::SSPRungeKutta3,
		  const StateType & state,
		  SystemType && systemObj)
    : name_(StepScheme::SSPRungeKutta3),
      order_(3),
      systemObj_(std::forward<SystemType>(systemObj)),
      velocities_{systemObj.createVelocity()},
      auxiliaryState_{::pressio::ops::clone(state)}
  {}

public:
  stepper_order_type order() const
  {
    return order_;
  }

  template<class StepCountType>
  void operator()(StateType & odeState,
		  const ScalarType & time,
		  const ScalarType & dt,
		  const StepCountType & step)
  {
    if (name_ == ode::StepScheme::ForwardEuler){
      doStepImpl(ode::ForwardEuler(), odeState, time, dt, step);
    }
    else if (name_ == ode::StepScheme::RungeKutta4){
      doStepImpl(ode::RungeKutta4(), odeState, time, dt, step);
    }
    else if (name_ == ode::StepScheme::AdamsBashforth2){
      doStepImpl(ode::AdamsBashforth2(), odeState, time, dt, step);
    }
    else if (name_ == ode::StepScheme::SSPRungeKutta3){
      doStepImpl(ode::SSPRungeKutta3(), odeState, time, dt, step);
    }
  }

private:
  template<class StepCountType>
  void doStepImpl(ode::ForwardEuler,
		  StateType & odeState,
		  const ScalarType & time,
		  const ScalarType & dt,
		  const StepCountType & step)
  {
    PRESSIOLOG_DEBUG("euler forward stepper: do step");

    auto & rhs = velocities_[0];
    //eval RHS
    systemObj_.get().velocity(odeState, time, rhs);
    // y = y + dt * rhs
    constexpr auto one  = ::pressio::utils::Constants<ScalarType>::one();
    ::pressio::ops::update(odeState, one, rhs, dt);
  }

  template<class StepCountType>
  void doStepImpl(ode::AdamsBashforth2,
		  StateType & odeState,
		  const ScalarType & currentTime,
		  const ScalarType & dt,
		  const StepCountType & stepNumber)
  {
    PRESSIOLOG_DEBUG("adams-bashforth2 stepper: do step");

    // y_n+1 = y_n + dt*[ (3/2)*f(y_n, t_n) - (1/2)*f(y_n-1, t_n-1) ]

    const auto cfn   = ::pressio::utils::Constants<ScalarType>::threeOvTwo()*dt;
    const auto cfnm1 = ::pressio::utils::Constants<ScalarType>::negOneHalf()*dt;

    if (stepNumber==1){
      // use Euler forward or we could use something else here maybe RK4
      auto & rhs = velocities_[0];
      systemObj_.get().velocity(odeState, currentTime, rhs);

      constexpr auto one   = ::pressio::utils::Constants<ScalarType>::one();
      ::pressio::ops::update(odeState, one, rhs, dt);
    }
    else{
      auto & fn   = velocities_[0];
      auto & fnm1 = velocities_[1];
      // fn -> fnm1
      ::pressio::ops::deep_copy(fnm1, fn);

      systemObj_.get().velocity(odeState, currentTime, fn);
      constexpr auto one   = ::pressio::utils::Constants<ScalarType>::one();
      ::pressio::ops::update(odeState, one, fn, cfn, fnm1, cfnm1);
    }
  }


  template<class StepCountType>
  void doStepImpl(ode::SSPRungeKutta3,
		  StateType & odeState,
		  const ScalarType & time,
		  const ScalarType & dt,
		  const StepCountType & step)
  {
    PRESSIOLOG_DEBUG("ssprk3 stepper: do step");

    constexpr auto zero  = ::pressio::utils::Constants<ScalarType>::zero();
    constexpr auto one   = ::pressio::utils::Constants<ScalarType>::one();
    constexpr auto two   = ::pressio::utils::Constants<ScalarType>::two();
    constexpr auto three = ::pressio::utils::Constants<ScalarType>::three();
    constexpr auto four  = ::pressio::utils::Constants<ScalarType>::four();
    constexpr auto oneOvTwo = one/two;
    constexpr auto oneOvThree = one/three;
    constexpr auto twoOvThree = two/three;
    constexpr auto threeOvFour = three/four;
    constexpr auto fourInv = one/four;

    auto & rhs0 = velocities_[0];

    // see e.g. https://gkeyll.readthedocs.io/en/latest/dev/ssp-rk.html

    // rhs(u_n, t_n)
    systemObj_.get().velocity(odeState, time, rhs0);
    // u_1 = u_n + dt * rhs(u_n, t_n)
    ::pressio::ops::update(auxiliaryState_, zero,
                           odeState,     one,
                           rhs0,            dt);

    // rhs(u_1, t_n+dt)
    systemObj_.get().velocity(auxiliaryState_, time+dt, rhs0);
    // u_2 = 3/4*u_n + 1/4*u_1 + 1/4*dt*rhs(u_1, t_n+dt)
    ::pressio::ops::update(auxiliaryState_, fourInv,
		                       odeState,     threeOvFour,
                           rhs0,            fourInv*dt);

    // rhs(u_2, t_n + 0.5*dt)
    systemObj_.get().velocity(auxiliaryState_, time + oneOvTwo*dt, rhs0);
    // u_n+1 = 1/3*u_n + 2/3*u_2 + 2/3*dt*rhs(u_2, t_n+0.5*dt)
    ::pressio::ops::update(odeState,     oneOvThree,
		                       auxiliaryState_, twoOvThree,
                           rhs0,            twoOvThree*dt);
  }

  template<class StepCountType>
  void doStepImpl(ode::RungeKutta4,
		  StateType & odeState,
		  const ScalarType & t,
		  const ScalarType & dt,
		  const StepCountType & step)
  {
    PRESSIOLOG_DEBUG("rk4 stepper: do step");

    auto & rhs0 = velocities_[0];
    auto & rhs1 = velocities_[1];
    auto & rhs2 = velocities_[2];
    auto & rhs3 = velocities_[3];

    constexpr auto two  = ::pressio::utils::Constants<ScalarType>::two();
    constexpr auto three  = ::pressio::utils::Constants<ScalarType>::three();
    constexpr auto six  = two * three;

    const ScalarType dt_half = dt / two;
    const ScalarType t_phalf = t + dt_half;
    const ScalarType dt6 = dt / six;
    const ScalarType dt3 = dt / three;

    // stage 1: ytmp = y + rhs0*dt_half;
    systemObj_.get().velocity(odeState, t, rhs0);
    this->rk4_stage_update_impl(auxiliaryState_, odeState, rhs0, dt_half);

    // stage 2: ytmp = y + rhs1*dt_half;
    systemObj_.get().velocity(auxiliaryState_, t_phalf, rhs1);
    this->rk4_stage_update_impl(auxiliaryState_, odeState, rhs1, dt_half);

    // stage 3: ytmp = y + rhs2*dt;
    systemObj_.get().velocity(auxiliaryState_, t_phalf, rhs2);
    this->rk4_stage_update_impl(auxiliaryState_, odeState, rhs2, dt);

    // stage 4: y_n += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    systemObj_.get().velocity(auxiliaryState_,  t + dt, rhs3);
    this->rk4_stage_update_impl(odeState, rhs0, rhs1, rhs2, rhs3, dt6, dt3);
  }

  template<class rhs_t>
  void rk4_stage_update_impl(StateType & yIn,
		    const StateType & stateIn,
		    const rhs_t & rhsIn,
		    ScalarType dtValue)
  {
    constexpr auto zero  = ::pressio::utils::Constants<ScalarType>::zero();
    constexpr auto one  = ::pressio::utils::Constants<ScalarType>::one();
    ::pressio::ops::update(yIn, zero, stateIn, one, rhsIn, dtValue);
  }

  template<class rhs_t>
  void rk4_stage_update_impl(StateType & stateIn,
		    const rhs_t & rhsIn0,
        const rhs_t & rhsIn1,
		    const rhs_t & rhsIn2,
        const rhs_t & rhsIn3,
		    ScalarType dt6,
        ScalarType dt3)
  {
    constexpr auto one  = ::pressio::utils::Constants<ScalarType>::one();
    ::pressio::ops::update(stateIn, one,
			   rhsIn0, dt6,
			   rhsIn1, dt3,
			   rhsIn2, dt3,
			   rhsIn3, dt6);
  }

};

}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_EXPLICIT_IMPL_ODE_EXPLICIT_EULER_STEPPER_IMPL_HPP_
