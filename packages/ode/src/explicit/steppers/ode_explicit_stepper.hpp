/*
//@HEADER
// ************************************************************************
//
// ode_explicit_stepper.hpp
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

#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_

#include "../../ode_fwd.hpp"
#include "./impl/ode_explicit_euler_stepper_impl.hpp"
#include "./impl/ode_explicit_runge_kutta4_stepper_impl.hpp"

namespace pressio{ namespace ode{ namespace explicitmethods{

template<
  typename tag,
  typename state_type,
  typename ...Args
  >
class Stepper : public StepperBase<
  Stepper<tag, state_type, Args...>
  >
{

  using this_t		= Stepper<tag, state_type, Args...>;
  using base_t		= StepperBase<this_t>;
  // need to friend base to allow it to access the () operator below
  friend base_t;

  using mytraits	= ::pressio::ode::details::traits<this_t>;
  using scalar_type	= typename mytraits::scalar_t;
  using system_type	= typename mytraits::model_t;
  using velocity_type	= typename mytraits::velocity_t;
  using policy_t	= typename mytraits::velocity_policy_t;

  // this is the impl class type which holds all the implement details
  using impl_class_t	= typename mytraits::impl_t;
  impl_class_t myImpl_  = {};

public:
  Stepper()  = delete;
  ~Stepper() = default;

  // this is enabled all the time
  Stepper(state_type const	  & stateIn0,
	  const system_type	  & model,
	  const policy_t	  & policyObj)
    : myImpl_(model,
	      policyObj,
	      stateIn0,
	      policyObj(stateIn0, model,
			::pressio::utils::constants::zero<scalar_type>())
	      )
  {}

  // only enable if the policy is standard
  template <
    typename T = policy_t,
    ::pressio::mpl::enable_if_t<
      mpl::is_same<
  	T, policy_t
  	>::value
      > * = nullptr
    >
  Stepper(const	state_type & stateIn0,
	  const system_type & model)
    : myImpl_(model,
  	      T(),
  	      stateIn0,
  	      T()(stateIn0, model,
		  ::pressio::utils::constants::zero<scalar_type>())
  	      )
  {}

private:
  // the compute method is private because we want users to use
  // the () operator in the base, which in turn calls compute here
  template<typename ... Args2>
  void compute(Args2 && ... args){
    myImpl_.doStep( std::forward<Args2>(args)... );
  }

};//end class

}}} // end namespace pressio::ode::explicitmethods
#endif
