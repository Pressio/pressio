/*
//@HEADER
// ************************************************************************
//
// ode_explicit_runge_kutta4_stepper_impl.hpp
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

#ifndef ODE_EXPLICIT_IMPL_ODE_EXPLICIT_SSP_RUNGE_KUTTA3_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_IMPL_ODE_EXPLICIT_SSP_RUNGE_KUTTA3_STEPPER_IMPL_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<
  class ScalarType,
  class StateType,
  class SystemType,
  class VelocityType
  >
class ExplicitSSPRungeKutta3Stepper
{
public:
  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  static constexpr stepper_order_type order_value = 3;

private:
  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
  std::array<VelocityType, 1> velocities_;
  StateType auxiliaryState_;

public:
  ExplicitSSPRungeKutta3Stepper() = delete;
  ExplicitSSPRungeKutta3Stepper(const ExplicitSSPRungeKutta3Stepper & other) = default;
  ExplicitSSPRungeKutta3Stepper & operator=(const ExplicitSSPRungeKutta3Stepper & other) = delete;
  ExplicitSSPRungeKutta3Stepper(ExplicitSSPRungeKutta3Stepper && other)  = default;
  ExplicitSSPRungeKutta3Stepper & operator=(ExplicitSSPRungeKutta3Stepper && other)  = delete;
  ~ExplicitSSPRungeKutta3Stepper() = default;

  ExplicitSSPRungeKutta3Stepper(const StateType & state,
                                 SystemType && systemObj)
    : systemObj_(std::forward<SystemType>(systemObj)),
      velocities_{systemObj.createVelocity()},
      auxiliaryState_{::pressio::ops::clone(state)}
  {}

public:
  stepper_order_type order() const
  {
    return order_value;
  }

  void doStep(StateType & odeSolution,
  	      const ScalarType & time,
  	      const ScalarType & dt,
  	      const step_count_type & step)
  {
    PRESSIOLOG_DEBUG("ssprk3 stepper: do step");

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
    systemObj_.get().velocity(odeSolution, time, rhs0);
    // u_1 = u_n + dt * rhs(u_n, t_n)
    ::pressio::ops::update(auxiliaryState_, odeSolution, one, rhs0, dt);

    // rhs(u_1, t_n+dt)
    systemObj_.get().velocity(auxiliaryState_, time+dt, rhs0);
    // u_2 = 3/4*u_n + 1/4*u_1 + 1/4*dt*rhs(u_1, t_n+dt)
    ::pressio::ops::update(auxiliaryState_, fourInv,
		   odeSolution, threeOvFour, rhs0, fourInv*dt);

    // rhs(u_2, t_n + 0.5*dt)
    systemObj_.get().velocity(auxiliaryState_, time + oneOvTwo*dt, rhs0);
    // u_n+1 = 1/3*u_n + 2/3*u_2 + 2/3*dt*rhs(u_2, t_n+0.5*dt)
    ::pressio::ops::update(odeSolution, oneOvThree,
		   auxiliaryState_, twoOvThree, rhs0, twoOvThree*dt);
  }
};

}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_EXPLICIT_IMPL_ODE_EXPLICIT_RUNGE_KUTTA4_STEPPER_IMPL_HPP_
