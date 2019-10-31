/*
//@HEADER
// ************************************************************************
//
// ode_implicit_residual_standard_policy_for_arbitrary_stepper.hpp
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

#ifndef ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_FOR_ARBITRARY_STEPPER_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_FOR_ARBITRARY_STEPPER_HPP_

#include "../../../ode_fwd.hpp"
#include "../base/ode_implicit_residual_policy_base.hpp"

namespace pressio{ namespace ode{ namespace policy{

template<
  typename state_type,
  typename system_type,
  typename residual_type
  >
class ImplicitResidualStandardPolicyForArbitraryStepper<
  state_type, system_type, residual_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::is_legitimate_implicit_residual_type<residual_type>::value and
    containers::meta::is_wrapper<state_type>::value and
    containers::meta::is_wrapper<residual_type>::value
    >
  >
  : public ImplicitResidualPolicyBase<
  ImplicitResidualStandardPolicyForArbitraryStepper<state_type, system_type, residual_type>>
{

  using this_t
  = ImplicitResidualStandardPolicyForArbitraryStepper<state_type, system_type, residual_type>;
  friend ImplicitResidualPolicyBase<this_t>;

public:
  ImplicitResidualStandardPolicyForArbitraryStepper() = default;
  ~ImplicitResidualStandardPolicyForArbitraryStepper() = default;

public:

  residual_type operator()(const state_type & odeCurrentState,
  			   const system_type & model) const{

    residual_type R(model.createTimeDiscreteResidualObject(*odeCurrentState.data()));
    return R;
  }

  //-------------------------------
  // specialize for n == 1
  //-------------------------------
  template <std::size_t n, typename scalar_type, mpl::enable_if_t<n==1> * = nullptr>
  void operator()(const state_type & odeCurrentState,
		  const ::pressio::ode::StatesContainer<state_type, n> & prevStates,
		  const system_type & model,
		  const scalar_type & t,
		  const scalar_type & dt,
		  const types::step_t & step,
		  residual_type & R) const{

    model.template timeDiscreteResidual(step, t, dt,
					*R.data(),
					*odeCurrentState.data(),
					*prevStates[0].data());
  }

  //-------------------------------
  // specialize for n == 2
  //-------------------------------
  template <std::size_t n, typename scalar_type, mpl::enable_if_t<n==2> * = nullptr>
  void operator()(const state_type & odeCurrentState,
		  const ::pressio::ode::StatesContainer<state_type, n> & prevStates,
		  const system_type & model,
		  const scalar_type & t,
		  const scalar_type & dt,
		  const types::step_t & step,
		  residual_type & R) const{

    model.template timeDiscreteResidual(step, t, dt,
					*R.data(),
					*odeCurrentState.data(),
					*prevStates[0].data(),
					*prevStates[1].data());
  }

  //-------------------------------
  // specialize for n == 3
  //-------------------------------
  template <std::size_t n, typename scalar_type, mpl::enable_if_t<n==3> * = nullptr>
  void operator()(const state_type & odeCurrentState,
		  const ::pressio::ode::StatesContainer<state_type, n> & prevStates,
		  const system_type & model,
		  const scalar_type & t,
		  const scalar_type & dt,
		  const types::step_t &  step,
		  residual_type & R) const{

    model.template timeDiscreteResidual(step, t, dt,
					*R.data(),
					*odeCurrentState.data(),
					*prevStates[0].data(),
					*prevStates[1].data(),
					*prevStates[2].data());
  }

};//end class

}}}//end namespace pressio::ode::policy
#endif










  // template <std::size_t n, typename scalar_type, mpl::enable_if_t<n==1> * = nullptr>
  // residual_type operator()(const state_type & odeCurrentState,
  // 			   const ::pressio::ode::StatesContainer<state_type, n> & prevStates,
  // 			   const system_type & model,
  // 			   const scalar_type & t,
  // 			   const scalar_type & dt,
  // 			   const types::step_t & step) const{

  //   residual_type R(model.template timeDiscreteResidual(step, t, dt,
  // 							*odeCurrentState.data(),
  // 							*prevStates[0].data() ));
  //   return R;
  // }

  // //-------------------------------
  // // specialize for n == 2
  // //-------------------------------
  // template <std::size_t n, typename scalar_type, mpl::enable_if_t<n==2> * = nullptr>
  // void operator()(const state_type & odeCurrentState,
  // 		  const ::pressio::ode::StatesContainer<state_type, n> & prevStates,
  // 		  const system_type & model,
  // 		  const scalar_type & t,
  // 		  const scalar_type & dt,
  // 		  const types::step_t & step,
  // 		  residual_type & R) const{

  //   model.template timeDiscreteResidual(step, t, dt,
  // 					*R.data(),
  // 					*odeCurrentState.data(),
  // 					*prevStates[0].data(),
  // 					*prevStates[1].data());
  // }

  // template <std::size_t n, typename scalar_type, mpl::enable_if_t<n==2> * = nullptr>
  // residual_type operator()(const state_type & odeCurrentState,
  // 			   const ::pressio::ode::StatesContainer<state_type, n> & prevStates,
  // 			   const system_type & model,
  // 			   const scalar_type & t,
  // 			   const scalar_type & dt,
  // 			   const types::step_t &  step) const{

  //   residual_type R( model.template timeDiscreteResidual(step, t, dt,
  // 							 *odeCurrentState.data(),
  // 							 *prevStates[0].data(),
  // 							 *prevStates[1].data() ));
  //   return R;
  // }

  // //-------------------------------
  // // specialize for n == 3
  // //-------------------------------
  // template <std::size_t n, typename scalar_type, mpl::enable_if_t<n==3> * = nullptr>
  // void operator()(const state_type & odeCurrentState,
  // 		  const ::pressio::ode::StatesContainer<state_type, n> & prevStates,
  // 		  const system_type & model,
  // 		  const scalar_type & t,
  // 		  const scalar_type & dt,
  // 		  const types::step_t &  step,
  // 		  residual_type & R) const{

  //   model.template timeDiscreteResidual(step, t, dt,
  // 					*R.data(),
  // 					*odeCurrentState.data(),
  // 					*prevStates[0].data(),
  // 					*prevStates[1].data(),
  // 					*prevStates[2].data());
  // }

  // template <std::size_t n, typename scalar_type, mpl::enable_if_t<n==3> * = nullptr>
  // residual_type operator()(const state_type & odeCurrentState,
  // 			   const ::pressio::ode::StatesContainer<state_type, n> & prevStates,
  // 			   const system_type & model,
  // 			   const scalar_type & t,
  // 			   const scalar_type & dt,
  // 			   const types::step_t &  step) const{

  //   residual_type R(model.template timeDiscreteResidual(step, t, dt,
  // 							*odeCurrentState.data(),
  // 							*prevStates[0].data(),
  // 							*prevStates[1].data(),
  // 							*prevStates[2].data() ));
  //   return R;
  // }
