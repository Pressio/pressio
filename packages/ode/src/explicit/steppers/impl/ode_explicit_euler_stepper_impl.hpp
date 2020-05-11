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

#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_EULER_STEPPER_IMPL_HPP_

namespace pressio{ namespace ode{ namespace explicitmethods{ namespace impl{

template<
  typename scalar_type,
  typename state_type,
  typename system_type,
  typename velocity_type,
  typename velocity_policy_type,
  typename ops_t
  >
class ExplicitEulerStepperImpl<scalar_type,
			       state_type,
			       system_type,
			       velocity_type,
			       velocity_policy_type,
			       ops_t>
{

  using this_t = ExplicitEulerStepperImpl< scalar_type,
					   state_type, system_type,
					   velocity_type,
					   velocity_policy_type,
					   ops_t>;

  using velocity_storage_t = VelocitiesContainer<velocity_type, 1>;
  using system_wrapper_t   = ::pressio::ode::impl::OdeSystemWrapper<system_type>;

  system_wrapper_t sys_;
  const velocity_policy_type & policy_;
  velocity_storage_t veloAuxStorage_;
  const ops_t * udOps_ = nullptr;

public:
  template <typename _ops_t = ops_t, mpl::enable_if_t< std::is_void<_ops_t>::value, int > = 0>
  ExplicitEulerStepperImpl(const system_type & model,
  			   const velocity_policy_type & policy_obj,
  			   const state_type & stateIn0,
  			   const velocity_type & f0)
    : sys_(model), policy_(policy_obj), veloAuxStorage_(f0)
  {}

  template <typename _ops_t = ops_t, mpl::enable_if_t< !std::is_void<_ops_t>::value, int > = 0>
  ExplicitEulerStepperImpl(const system_type & model,
  			   const velocity_policy_type & policy_obj,
  			   const state_type & stateIn0,
  			   const velocity_type & f0,
			   const _ops_t & udOps)
    : sys_(model), policy_(policy_obj), veloAuxStorage_(f0), udOps_(&udOps)
  {}

  ExplicitEulerStepperImpl() = delete;
  ~ExplicitEulerStepperImpl() = default;

  // copy cnstr
  ExplicitEulerStepperImpl(const ExplicitEulerStepperImpl & other)  = delete;
  // copy assignment
  ExplicitEulerStepperImpl & operator=(const ExplicitEulerStepperImpl & other)  = delete;
  // move cnstr
  ExplicitEulerStepperImpl(ExplicitEulerStepperImpl && other)  = delete;
  // move assign
  ExplicitEulerStepperImpl & operator=(ExplicitEulerStepperImpl && other)  = delete;

public:

  /* user does NOT provide custom ops, so we use ops */
  template<
    typename _ops_t = ops_t,
    typename _state_type = state_type
  >
  mpl::enable_if_t< std::is_void<_ops_t>::value >
  doStep(_state_type & stateInOut,
	      const scalar_type & time,
	      const scalar_type & dt,
	      const types::step_t & step)
  {
    auto & auxRhs0 = veloAuxStorage_(0);
    //eval RHS
    policy_(stateInOut, auxRhs0, sys_.get(), time);
    // y = y + dt * rhs
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();
    ::pressio::ops::do_update(stateInOut, one, auxRhs0, dt);
  }

  /* user provides custom ops */
  template<
    typename _ops_t = ops_t,
    typename _state_type = state_type
  >
  mpl::enable_if_t< 
    !std::is_void<_ops_t>::value and
    containers::meta::is_wrapper<_state_type>::value
  > 
  doStep(_state_type & stateInOut,
  	      const scalar_type & time,
  	      const scalar_type & dt,
  	      const types::step_t & step)
  {
    auto & auxRhs0 = veloAuxStorage_(0);
    policy_(stateInOut, auxRhs0, sys_.get(), time);
    // y = y + dt * rhs
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();
    udOps_->do_update(*stateInOut.data(), one, *auxRhs0.data(), dt);
  }
};

}}}}//end namespace pressio::ode::explicitmethods::impl
#endif
