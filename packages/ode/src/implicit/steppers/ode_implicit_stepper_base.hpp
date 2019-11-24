/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_base.hpp
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

#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_

#include "../../ode_aux_states_container.hpp"
#include "../../ode_system_wrapper.hpp"

namespace pressio{ namespace ode{

/*
 * (1) constructors here should be private but we need
 * them public to enable interfacing with pybind11
 */
template<typename concrete_stepper_type>
class ImplicitStepperBase
{
  using traits			= typename details::traits<concrete_stepper_type>;
  using sc_t			= typename traits::scalar_t;
  using state_t			= typename traits::state_t;
  using residual_t		= typename traits::residual_t;
  using jacobian_t		= typename traits::jacobian_t;
  using standard_res_policy_t	= typename traits::standard_res_policy_t;
  using standard_jac_policy_t	= typename traits::standard_jac_policy_t;
  using residual_pol_t		= typename traits::residual_policy_t;
  using jacobian_pol_t		= typename traits::jacobian_policy_t;
  using system_t		= typename traits::system_t;
  using system_wrapper_t	= impl::OdeSystemWrapper<system_t>;

  //do checking here that things are as supposed
  static_assert( meta::is_legitimate_implicit_state_type<state_t>::value,
       "OOPS: STATE_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::is_legitimate_implicit_residual_type<residual_t>::value,
       "OOPS: RESIDUAL_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::is_legitimate_jacobian_type<jacobian_t>::value,
       "OOPS: JACOBIAN_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");

public:
  using aux_states_t = ::pressio::ode::AuxStatesContainer<false, state_t, traits::numAuxStates>;

  types::stepper_order_t order() const{
    return traits::order_value;
  }

  void residual(const state_t & odeState, residual_t & R) const{
    static_cast<const concrete_stepper_type &>(*this).residualImpl(odeState, R);
  }

  residual_t residual(const state_t & odeState) const{
    return static_cast<const concrete_stepper_type &>(*this).residualImpl(odeState);
  }

  void jacobian(const state_t & odeState, jacobian_t & J) const{
    static_cast<const concrete_stepper_type &>(*this).jacobianImpl(odeState, J);
  }

  jacobian_t jacobian(const state_t & odeState) const{
    return static_cast<const concrete_stepper_type &>(*this).jacobianImpl(odeState);
  }

protected:
  // procted because these are accessed only by children classes
  sc_t t_  = {};
  sc_t dt_ = {};
  types::step_t step_  = {};
  system_wrapper_t sys_;
  aux_states_t auxStates_;

  // conditionally set the type of the object knowing how to compute residual
  // if we have a standard policy, then it takes a copy
  // if we have a user-defined policy, we take a const & to it
  typename std::conditional<
    mpl::is_same<standard_res_policy_t, residual_pol_t>::value,
    const residual_pol_t,
    const residual_pol_t &
    >::type residual_obj_;

  // conditionally set the type of the object knowing how to compute jacobian
  // if we have a standard policy, then it takes a copy
  // if we have a user-defined policy, we take a const & to it
  typename std::conditional<
    mpl::is_same<standard_jac_policy_t, jacobian_pol_t>::value,
    const jacobian_pol_t,
    const jacobian_pol_t &
    >::type jacobian_obj_;

public:
  ImplicitStepperBase() = delete;
  ~ImplicitStepperBase() = default;

  ImplicitStepperBase(const state_t & stateIn0,
		      const system_t & model,
		      const residual_pol_t & resPolicyObj,
		      const jacobian_pol_t & jacPolicyObj)
    : sys_{model},
      auxStates_{stateIn0},
      residual_obj_{resPolicyObj},
      jacobian_obj_{jacPolicyObj}
  {}

  // cstr for standard residual and jacob policies
  template <
    typename T1 = standard_res_policy_t,
    typename T2 = standard_jac_policy_t,
    ::pressio::mpl::enable_if_t<
      mpl::is_same<T1, residual_pol_t>::value and
      mpl::is_same<T2, jacobian_pol_t>::value
      > * = nullptr
    >
  ImplicitStepperBase(const state_t & stateIn0,
  		      const system_t & model)
    : sys_{model},
      auxStates_{stateIn0},
      residual_obj_{},
      jacobian_obj_{}
  {}

  // cstr for standard jacob policies
  template <
    typename T2 = standard_jac_policy_t,
    ::pressio::mpl::enable_if_t<
      mpl::is_same<T2, jacobian_pol_t>::value
      > * = nullptr
    >
  ImplicitStepperBase(const state_t & stateIn0,
  		      const system_t & model,
  		      const residual_pol_t & resPolicyObj)
    : sys_{model},
      auxStates_{stateIn0},
      residual_obj_{resPolicyObj},
      jacobian_obj_{}
  {}

};//end class

}}//end namespace pressio::ode
#endif
