/*
//@HEADER
// ************************************************************************
//
// rom_masked_jacobian_policy.hpp
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

#ifndef ROM_DECORATORS_ROM_MASKED_JACOBIAN_POLICY_HPP_
#define ROM_DECORATORS_ROM_MASKED_JACOBIAN_POLICY_HPP_

namespace pressio{ namespace rom{ namespace decorator{


template <typename maskable_policy>
class MaskedJacobianPolicy : public maskable_policy
{

public:
  using typename maskable_policy::apply_jac_return_t;
  using maskable_policy::JJ_;
  using maskable_policy::fomStatesMngr_;

public:
  MaskedJacobianPolicy() = delete;
  MaskedJacobianPolicy(const maskable_policy & obj) : maskable_policy(obj)
  {}

  template <typename ... Args>
  MaskedJacobianPolicy(Args && ... args) : maskable_policy(std::forward<Args>(args)...)
  {}

public:
  template <typename app_t>
  apply_jac_return_t create(const app_t & app) const
  {
    JJ_ = maskable_policy::create(app);
    return apply_jac_return_t(app.createApplyMaskResult(*JJ_.data()));
  }

  // unsteady case
  template <
    typename stepper_tag,
    typename state_t,
    typename prev_states_mgr,
    typename app_t,
    typename time_type,
    typename jac_t
  >
  void compute(const state_t & state,
      const prev_states_mgr & prevStatesMgr,
		  const app_t & app,
		  const time_type & time,
		  const time_type & dt,
		  const ::pressio::ode::types::step_t & step,
		  jac_t & odeJJ) const
  {
    maskable_policy::template compute<stepper_tag>(state, prevStatesMgr, app, time, dt, step, JJ_);
    app.applyMask(*JJ_.data(), time, *odeJJ.data());
  }

  // steady case
  template <typename state_t, typename jac_t, typename app_t>
  void compute(const state_t & state, jac_t & odeJJ, const app_t & app) const
  {
    maskable_policy::compute(state, JJ_, app);
    app.applyMask(*JJ_.data(), *odeJJ.data());    
  }

};//end class

}}} //end namespace pressio::rom::decorator
#endif  // ROM_DECORATORS_ROM_MASKED_JACOBIAN_POLICY_HPP_
