/*
//@HEADER
// ************************************************************************
//
// ode_implicit_residual_crank_nicolson_policy.hpp
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

#ifndef ODE_IMPLICIT_IMPL_STANDARD_POLICIES_ODE_IMPLICIT_RESIDUAL_CRANK_NICOLSON_POLICY_HPP_
#define ODE_IMPLICIT_IMPL_STANDARD_POLICIES_ODE_IMPLICIT_RESIDUAL_CRANK_NICOLSON_POLICY_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{ namespace policy{

template<typename state_type, typename residual_type>
class ResidualStandardPolicyCrankNicolson
{
  static_assert
  (::pressio::ode::constraints::implicit_state<state_type>::value,
   "Invalid state type for standard residual policy");

  static_assert
  (::pressio::ode::constraints::implicit_residual<residual_type>::value,
   "Invalid residual type for standard residual policy");

  mutable ::pressio::ode::types::step_t stepTracker_ = -1;

public:
  ResidualStandardPolicyCrankNicolson() = default;
  ResidualStandardPolicyCrankNicolson(const ResidualStandardPolicyCrankNicolson &) = default;
  ResidualStandardPolicyCrankNicolson & operator=(const ResidualStandardPolicyCrankNicolson &) = default;
  ResidualStandardPolicyCrankNicolson(ResidualStandardPolicyCrankNicolson &&) = default;
  ResidualStandardPolicyCrankNicolson & operator=(ResidualStandardPolicyCrankNicolson &&) = default;
  ~ResidualStandardPolicyCrankNicolson() = default;

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
    class stencil_velocities_type,
    class system_type,
    class scalar_type
    >
  mpl::enable_if_t<
    std::is_same<ode_tag, ::pressio::ode::implicitmethods::CrankNicolson>::value
    >
  compute(const state_type & predictedState,
	  const stencil_states_type & stencilStates,
	  const system_type & system,
	  const scalar_type & t_np1,
	  const scalar_type & dt,
	  const types::step_t & step,
	  stencil_velocities_type & stencilVelocities,
	  residual_type & R) const
  {
    static_assert
    (constraints::continuous_time_system_with_user_provided_jacobian
     <system_type>::value,
     "system type must meet the continuous time api");

    if (stepTracker_ != step){
      auto & f_n = stencilVelocities(::pressio::ode::n());
      auto & state_n = stencilStates(::pressio::ode::n());
      const auto tn = t_np1-dt;
      system.velocity(*state_n.data(), tn, *f_n.data());
    }

    auto & f_np1 = stencilVelocities(::pressio::ode::nPlusOne());
    system.velocity(*predictedState.data(), t_np1, *f_np1.data());
    ::pressio::ode::impl::discrete_time_residual
    	(predictedState, R, stencilStates, stencilVelocities,
	 dt, ode_tag());


    // if (stepTracker_ != step){
    //   auto & state_n = stencilDataManager.stateAt(::pressio::ode::n());
    //   const auto tn = t_np1-dt;
    //   system.velocity(*state_n.data(), tn, *f_n_.data());
    // }
    // system.velocity(*predictedState.data(), t_np1, *R.data());
    // ::pressio::ode::impl::discrete_time_residual
    // 	(predictedState, R, stencilDataManager, f_n_, dt, ode_tag());

    // if (step==1 and stepTracker_ == -1)
    // {
    //   auto & state_n = stencilDataManager(::pressio::ode::n());
    //   const auto tn = t_np1-dt;
    //   system.velocity(*state_n.data(), tn, *f_n_.data());
    // }
    // if (stepTracker_ != step and step!=1){
    //   ::pressio::ops::deep_copy(f_n_, f_np1_);
    // }
    // system.velocity(*predictedState.data(), t_np1, *f_np1_.data());
    // ::pressio::ode::impl::discrete_time_residual
    // 	(predictedState, R, stencilDataManager, f_n_, f_np1_, dt, ode_tag());

    stepTracker_ = step;
  }
};

}}}}//end namespace pressio::ode::implicitmethods::policy
#endif  // ODE_IMPLICIT_IMPL_STANDARD_POLICIES_ODE_IMPLICIT_RESIDUAL_CRANK_NICOLSON_POLICY_HPP_
