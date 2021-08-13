/*
//@HEADER
// ************************************************************************
//
// ode_explicit_adams_bashforth2_stepper_impl.hpp
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

#ifndef ODE_EXPLICIT_IMPL_ODE_EXPLICIT_ADAMS_BASHFORTH2_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_IMPL_ODE_EXPLICIT_ADAMS_BASHFORTH2_STEPPER_IMPL_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<
  class ScalarType,
  class StateType,
  class SystemType,
  class VelocityType
  >
class ExplicitAdamsBashforth2Stepper
{
public:
  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  static constexpr stepper_order_type order_value = 2;

private:
  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
  std::array<VelocityType, 2> velocities_;

public:
  ExplicitAdamsBashforth2Stepper() = delete;
  ExplicitAdamsBashforth2Stepper(const ExplicitAdamsBashforth2Stepper & other) = default;
  ExplicitAdamsBashforth2Stepper & operator=(const ExplicitAdamsBashforth2Stepper & other) = delete;
  ExplicitAdamsBashforth2Stepper(ExplicitAdamsBashforth2Stepper && other)  = default;
  ExplicitAdamsBashforth2Stepper & operator=(ExplicitAdamsBashforth2Stepper && other)  = delete;
  ~ExplicitAdamsBashforth2Stepper() = default;

  ExplicitAdamsBashforth2Stepper(const StateType & state,
                                 SystemType && systemObj)
    : systemObj_(std::forward<SystemType>(systemObj)),
      velocities_{systemObj.createVelocity(), systemObj.createVelocity()}
  {}

public:
  stepper_order_type order() const
  {
    return order_value;
  }

  template<class StepCountType>
  void operator()(StateType & odeSolution,
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
      systemObj_.get().velocity(odeSolution, currentTime, rhs);
      doUpdate1(odeSolution, rhs, dt);
    }
    else{
      auto & fn   = velocities_[0];
      auto & fnm1 = velocities_[1];
      // fn -> fnm1
      updateStoredVelocities(fnm1, fn);

      systemObj_.get().velocity(odeSolution, currentTime, fn);
      doUpdate2(odeSolution, fn, cfn, fnm1, cfnm1);
    }
  }

private:
  template<class T>
  void updateStoredVelocities(T & to, const T & from)
  {
    ::pressio::ops::deep_copy(to, from);
  }

  template<typename f_t>
  void doUpdate1(StateType & odeSolution,
	    const f_t & rhs,
	    const ScalarType & dt)
  {
    constexpr auto one   = ::pressio::utils::Constants<ScalarType>::one();
    ::pressio::ops::update(odeSolution, one, rhs, dt);
  }

  template<typename f_t>
  void doUpdate2(StateType & odeSolution,
	    const f_t & fn,   const ScalarType & cfn,
	    const f_t & fnm1, const ScalarType & cfnm1)
  {
    constexpr auto one   = ::pressio::utils::Constants<ScalarType>::one();
    ::pressio::ops::update(odeSolution, one, fn, cfn, fnm1, cfnm1);
  }
};

}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_EXPLICIT_IMPL_ODE_EXPLICIT_ADAMS_BASHFORTH2_STEPPER_IMPL_HPP_
