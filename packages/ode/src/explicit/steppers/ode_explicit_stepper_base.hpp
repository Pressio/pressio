/*
//@HEADER
// ************************************************************************
//
// ode_explicit_stepper_base.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_BASE_EXPLICIT_STEPPER_BASE_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_BASE_EXPLICIT_STEPPER_BASE_HPP_

#include "ode_explicit_stepper_traits.hpp"
#include "../policies/ode_is_legitimate_explicit_velocity_policy.hpp"
#include "../../ode_storage.hpp"
#include "../../ode_system_wrapper.hpp"

namespace pressio{ namespace ode{

/*
 * (1) constructors here should be private but we need
 * them public to enable interfacing with pybind11
 */

template<typename stepper_type>
class ExplicitStepperBase
{
private:
  using step_traits	= ode::details::traits<stepper_type>;
  using scalar_t	= typename step_traits::scalar_t;
  using state_t		= typename step_traits::state_t;
  using velocity_t	= typename step_traits::velocity_t;
  using model_t		= typename step_traits::model_t;
  using policy_t	= typename step_traits::velocity_policy_t;

  static_assert( meta::is_legitimate_explicit_state_type<state_t>::value,
  "OOPS: STATE_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");

  static_assert( meta::is_legitimate_explicit_velocity_type<velocity_t>::value,
  "OOPS: VELOCITY_TYPE IN SELECTED EXPLICIT STEPPER IS NOT VALID");

  static_assert( meta::is_legitimate_explicit_velocity_policy<
  		 policy_t>::value,
  "VELOCITY_POLICY NOT ADMISSIBLE: MAYBE NOT INHERITING FROM EXPLICIT POLICY BASE");

public:
  ExplicitStepperBase() = default;
  ~ExplicitStepperBase() = default;

public:
  typename step_traits::order_t order() const{
    return step_traits::order_value;
  }

  template<typename step_t>
  void operator()(state_t & yinout,
		  scalar_t t,
		  scalar_t dt,
		  step_t step){
    static_cast<stepper_type&>(*this).compute(yinout, t, dt, step);
  }

};//end class

}}//end namespace pressio::ode
#endif
