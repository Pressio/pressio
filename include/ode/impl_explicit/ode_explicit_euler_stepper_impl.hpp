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

#ifndef ODE_EXPLICIT_IMPL_ODE_EXPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_IMPL_ODE_EXPLICIT_EULER_STEPPER_IMPL_HPP_

namespace pressio{ namespace ode{ namespace explicitmethods{ namespace impl{

template<
  class scalar_type,
  class state_type,
  class system_type,
  class velocity_type,
  class velocity_policy_type,
  bool is_standard_policy
  >
class ExplicitEulerStepper
{

public:
  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  static constexpr stepper_order_type order_value = 1;

private:
  std::reference_wrapper<const system_type> systemObj_;
  ::pressio::utils::instance_or_reference_wrapper<velocity_policy_type> policy_;
  std::array<velocity_type, 1> velocities_;

public:
  ExplicitEulerStepper() = delete;
  ExplicitEulerStepper(const ExplicitEulerStepper &) = default;
  ExplicitEulerStepper & operator=(const ExplicitEulerStepper &) = delete;
  ExplicitEulerStepper(ExplicitEulerStepper &&) = default;
  ExplicitEulerStepper & operator=(ExplicitEulerStepper &&) = delete;
  ~ExplicitEulerStepper() = default;

  ExplicitEulerStepper(const state_type & state,
		       const system_type & systemObj,
		       velocity_policy_type && policy)
    : systemObj_(systemObj),
      policy_(std::forward<velocity_policy_type>(policy)),
      velocities_{policy_.get().create(systemObj)}
  {}

  // only enabled if policy standard 
  template <
    bool _is_standard_policy = is_standard_policy,
    mpl::enable_if_t<_is_standard_policy,int > = 0
    >
  ExplicitEulerStepper(const state_type & state,
		       const system_type & systemObj)
    : systemObj_(systemObj),
      policy_(),
      velocities_{policy_.get().create(systemObj)}
  {}

public:
  stepper_order_type order() const
  {
    return order_value;
  }

  void doStep(state_type & odeSolution,
	 const scalar_type & time,
	 const scalar_type & dt,
	 const step_count_type & step)
  {
    PRESSIOLOG_DEBUG("euler forward stepper: do step");

    auto & rhs = velocities_[0];
    //eval RHS
    policy_.get().compute(odeSolution, rhs, systemObj_.get(), time);
    // y = y + dt * rhs
    constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
    ::pressio::ops::update(odeSolution, one, rhs, dt);
  }
};

}}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_EXPLICIT_IMPL_ODE_EXPLICIT_EULER_STEPPER_IMPL_HPP_