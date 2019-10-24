/*
//@HEADER
// ************************************************************************
//
// ode_implicit_jacobian_standard_policy_for_arbitrary_stepper.hpp
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

#ifndef ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_FOR_ARBITRARY_STEPPER_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_FOR_ARBITRARY_STEPPER_HPP_

#include "../../../ode_fwd.hpp"
#include "../base/ode_jacobian_policy_base.hpp"

namespace pressio{ namespace ode{ namespace policy{

template<
  typename state_type,
  typename system_type,
  typename jacobian_type
  >
class ImplicitJacobianStandardPolicyForArbitraryStepper<
  state_type, system_type, jacobian_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::is_legitimate_jacobian_type<jacobian_type>::value and
    containers::meta::is_wrapper<state_type>::value and
    containers::meta::is_wrapper<jacobian_type>::value
    >
  > : public JacobianPolicyBase<ImplicitJacobianStandardPolicyForArbitraryStepper<
    state_type, system_type, jacobian_type> >
{

  using this_t
  = ImplicitJacobianStandardPolicyForArbitraryStepper<state_type, system_type, jacobian_type>;
  friend JacobianPolicyBase<this_t>;

public:
  ImplicitJacobianStandardPolicyForArbitraryStepper() = default;
  ~ImplicitJacobianStandardPolicyForArbitraryStepper() = default;

public:

  //-------------------------------
  // specialize for n == 1
  //-------------------------------
  template <int n, typename scalar_type, mpl::enable_if_t<n==1> * = nullptr>
  void operator()(const state_type & stateIn,
		  const std::array<state_type, 1> & oldStates,
		  const system_type & model,
		  const scalar_type & t,
		  const scalar_type & dt,
		  const types::step_t &  step,
		  jacobian_type & J) const
  {
    model.template timeDiscreteJacobian(step, t, dt,
					*J.data(),
					*stateIn.data(),
					*oldStates[0].data());
  }

  template <int n, typename scalar_type, mpl::enable_if_t<n==1> * = nullptr>
  jacobian_type operator()(const state_type & stateIn,
			   const std::array<state_type, 1> & oldStates,
  			   const system_type & model,
  			   const scalar_type & t,
  			   const scalar_type & dt,
			   const types::step_t & step) const
  {
    jacobian_type JJ(model.template timeDiscreteJacobian(step, t, dt, *stateIn.data(), *oldStates[0].data()));
    return JJ;
  }

  //-------------------------------
  // specialize for n == 2
  //-------------------------------
  template <int n, typename scalar_type, mpl::enable_if_t<n==2> * = nullptr>
  void operator()(const state_type & stateIn,
		  const std::array<state_type, 2> & oldStates,
		  const system_type & model,
		  const scalar_type & t,
		  const scalar_type & dt,
		  const types::step_t & step,
		  jacobian_type & J) const
  {
    model.template timeDiscreteJacobian(step, t, dt,
					*J.data(),
					*stateIn.data(),
					*oldStates[0].data(),
					*oldStates[1].data());
  }

  template <int n, typename scalar_type, mpl::enable_if_t<n==2> * = nullptr>
  jacobian_type operator()(const state_type & stateIn,
			   const std::array<state_type, 2> & oldStates,
  			   const system_type & model,
  			   const scalar_type &  t,
  			   const scalar_type &  dt,
			   const types::step_t & step) const
  {
    jacobian_type JJ(model.template timeDiscreteJacobian(step, t, dt,
							 *stateIn.data(),
							 *oldStates[0].data(),
							 *oldStates[1].data()));
    return JJ;
  }

  //-------------------------------
  // specialize for n == 3
  //-------------------------------
  template <int n, typename scalar_type, mpl::enable_if_t<n==3> * = nullptr>
  void operator()(const state_type & stateIn,
		  const std::array<state_type, 3> & oldStates,
		  const system_type & model,
		  const scalar_type & t,
		  const scalar_type & dt,
		  const types::step_t & step,
		  jacobian_type & J) const
  {
    model.template timeDiscreteJacobian(step, t, dt,
					*J.data(),
					*stateIn.data(),
					*oldStates[0].data(),
					*oldStates[1].data(),
					*oldStates[2].data());
  }

  template <int n, typename scalar_type, mpl::enable_if_t<n==3> * = nullptr>
  jacobian_type operator()(const state_type & stateIn,
			   const std::array<state_type, 3> & oldStates,
  			   const system_type & model,
  			   const scalar_type &  t,
  			   const scalar_type &  dt,
			   const types::step_t & step) const
  {
    jacobian_type JJ(model.timeDiscreteJacobian(step, t, dt,
						*stateIn.data(),
						*oldStates[0].data(),
						*oldStates[1].data(),
						*oldStates[2].data()));
    return JJ;
  }

};//end class

}}}//end namespace pressio::ode::policy
#endif
