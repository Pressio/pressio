/*
//@HEADER
// ************************************************************************
//
// ode_implicit_jacobian_standard_policy.hpp
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

#ifndef ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_HPP_

#include "../../../ode_fwd.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../../ode_time_discrete_jacobian.hpp"
#include "../../meta/ode_is_legitimate_implicit_state_type.hpp"
#include "../../meta/ode_is_legitimate_implicit_jacobian_type.hpp"

namespace pressio{ namespace ode{ namespace policy{

template<
  typename state_type,
  typename system_type,
  typename jacobian_type
  >
class ImplicitJacobianStandardPolicy<
  state_type, system_type, jacobian_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::is_legitimate_jacobian_type<jacobian_type>::value and
    containers::meta::is_wrapper<state_type>::value and
    containers::meta::is_wrapper<jacobian_type>::value
    >
  > : public JacobianPolicyBase<ImplicitJacobianStandardPolicy<
    state_type, system_type, jacobian_type> >
{

  using this_t = ImplicitJacobianStandardPolicy<state_type, system_type, jacobian_type>;
  friend JacobianPolicyBase<this_t>;

public:
  ImplicitJacobianStandardPolicy() = default;
  ~ImplicitJacobianStandardPolicy() = default;

public:

  template <
    ode::ImplicitEnum method, typename scalar_t
  >
  void operator()(const state_type & odeCurrentState,
		  jacobian_type & J,
		  const system_type & model,
		  const scalar_t & t,
		  const scalar_t & dt,
		  const types::step_t & step) const
  {
    model.jacobian( *odeCurrentState.data(), t, *J.data());
    ::pressio::ode::impl::time_discrete_jacobian<method>(J, dt);
  }

  template <
    ode::ImplicitEnum method,
    typename scalar_t
    >
  jacobian_type operator()(const state_type & odeCurrentState,
  			   const system_type & model,
  			   const scalar_t & t,
  			   const scalar_t & dt,
			   const types::step_t & step) const
  {
    jacobian_type JJ(model.jacobian(*odeCurrentState.data(), t));
    ::pressio::ode::impl::time_discrete_jacobian<method>(JJ, dt);
    return JJ;
  }

};//end class

}}}//end namespace pressio::ode::policy
#endif
