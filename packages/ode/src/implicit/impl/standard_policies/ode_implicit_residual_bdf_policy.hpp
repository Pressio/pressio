/*
//@HEADER
// ************************************************************************
//
// ode_implicit_residual_bdf_policy.hpp
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

#ifndef ODE_IMPLICIT_IMPL_STANDARD_POLICIES_ODE_IMPLICIT_RESIDUAL_BDF_POLICY_HPP_
#define ODE_IMPLICIT_IMPL_STANDARD_POLICIES_ODE_IMPLICIT_RESIDUAL_BDF_POLICY_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{ namespace policy{

template<typename state_type, typename residual_type>
class ResidualStandardPolicyBdf
{
  static_assert
  (::pressio::ode::constraints::implicit_state<state_type>::value,
   "Invalid state type for standard residual policy");

  static_assert
  (::pressio::ode::constraints::implicit_residual<residual_type>::value,
   "Invalid residual type for standard residual policy");

public:
  ResidualStandardPolicyBdf() = default;
  ResidualStandardPolicyBdf(const ResidualStandardPolicyBdf &) = default;
  ResidualStandardPolicyBdf & operator=(const ResidualStandardPolicyBdf &) = default;
  ResidualStandardPolicyBdf(ResidualStandardPolicyBdf &&) = default;
  ResidualStandardPolicyBdf & operator=(ResidualStandardPolicyBdf &&) = default;
  ~ResidualStandardPolicyBdf() = default;

public:
  template <typename system_type>
  residual_type create(const system_type & system) const
  {
    static_assert
      (::pressio::ode::constraints::continuous_time_system_with_user_provided_jacobian
       <system_type>::value,
       "system type must meet the continuous time api");

    residual_type R(system.createVelocity());
    return R;
  }

  template <
    class ode_tag,
    class stencil_states_type,
    class system_type,
    class scalar_type
    >
  mpl::enable_if_t<
    std::is_same<ode_tag, ::pressio::ode::implicitmethods::BDF1>::value or
    std::is_same<ode_tag, ::pressio::ode::implicitmethods::BDF2>::value
    >
  compute(const state_type & predictedState,
	  const stencil_states_type & stencilStatesManager,
	  const system_type & system,
	  const scalar_type & rhsEvaluationTime,
	  const scalar_type & dt,
	  const types::step_t & step,
	  residual_type & R) const
  {
    static_assert
      (::pressio::ode::constraints::continuous_time_system_with_user_provided_jacobian
       <system_type>::value,
       "system type must meet the continuous time api");

    try{
      system.velocity(*predictedState.data(), rhsEvaluationTime, *R.data());
      ::pressio::ode::impl::discrete_time_residual(predictedState,
						   R, stencilStatesManager,
						   dt, ode_tag());
    }
    catch (::pressio::eh::velocity_failure_unrecoverable const & e){
      throw ::pressio::eh::residual_evaluation_failure_unrecoverable();
    }
  }
};

}}}}//end namespace pressio::ode::implicitmethods::policy
#endif  // ODE_IMPLICIT_IMPL_STANDARD_POLICIES_ODE_IMPLICIT_RESIDUAL_BDF_POLICY_HPP_
