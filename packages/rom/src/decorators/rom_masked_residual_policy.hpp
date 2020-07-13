/*
//@HEADER
// ************************************************************************
//
// rom_masked_residual_policy.hpp
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

#ifndef ROM_DECORATORS_ROM_MASKED_RESIDUAL_POLICY_HPP_
#define ROM_DECORATORS_ROM_MASKED_RESIDUAL_POLICY_HPP_

namespace pressio{ namespace rom{ namespace decorator{


template <typename maskable_policy>
class MaskedResidualPolicy : public maskable_policy
{
  using typename maskable_policy::residual_t;
  using typename maskable_policy::ud_ops_t;
  using maskable_policy::R_;
  using maskable_policy::udOps_;
  using maskable_policy::fomStatesMngr_;

public:
  MaskedResidualPolicy() = delete;
  MaskedResidualPolicy(const maskable_policy & obj) : maskable_policy(obj)
  {}

  template <typename ... Args>
  MaskedResidualPolicy(Args && ... args) : maskable_policy(std::forward<Args>(args)...)
  {}

public:
  template <typename app_t>
  residual_t create(const app_t	& app) const
  {
    R_ = maskable_policy::create(app);
    return residual_t(app.createApplyMaskResult(*R_.data()));
  }

  //-------------------------------
  // unsteady case
  //-------------------------------
  template <
    typename stepper_tag,
    typename state_t,
    typename prev_states_t,
    typename app_t,
    typename time_type,
    typename norm_value_type
  >
  void compute(const state_t & state,
		  const prev_states_t & prevStates,
  		const app_t	& app,
		  const time_type & time,
		  const time_type & dt,
		  const ::pressio::ode::types::step_t & step,
		  residual_t & R,
		  ::pressio::Norm normKind,
		  norm_value_type & normValue) const
  {
    maskable_policy::template compute<stepper_tag>(state, prevStates, app, time, 
        dt, step, R_, normKind, normValue);

    app.applyMask(*R_.data(), time, *R.data());

    // need to compute norm since after the mask the object is a different 
    comouteNormOfMaskedResidual(R, normKind, normValue);
  }

  //-------------------------------
  // steady case
  //-------------------------------
  template <typename state_t, typename app_t, typename norm_value_type>
  void compute(const state_t & state,
      residual_t & R,
      const app_t & app,
      ::pressio::Norm normKind,
      norm_value_type & normValue) const
  {
    maskable_policy::compute(state, R_, app, normKind, normValue);

    app.applyMask(*R_.data(), *R.data());

    // need to compute norm since after the mask the object is a different 
    comouteNormOfMaskedResidual(R, normKind, normValue);
  }

private:
  template <typename norm_value_type, typename _ud_ops_t = ud_ops_t>
  ::pressio::mpl::enable_if_t< !std::is_void<_ud_ops_t>::value >
  comouteNormOfMaskedResidual(const residual_t & R, 
                              ::pressio::Norm normKind, 
                              norm_value_type & normValue) const
  {
    if (normKind == ::pressio::Norm::L2)
      normValue = udOps_->norm2(*R_.data());
    else if (normKind == ::pressio::Norm::L1)
      normValue = udOps_->norm1(*R_.data());
    else
      throw std::runtime_error("Invalid norm kind for lspg unsteady residual policy");
  }

  template <typename norm_value_type, typename _ud_ops_t = ud_ops_t>
  ::pressio::mpl::enable_if_t< std::is_void<_ud_ops_t>::value >
  comouteNormOfMaskedResidual(const residual_t & R, 
                              ::pressio::Norm normKind, 
                              norm_value_type & normValue) const
  {
    if (normKind == ::pressio::Norm::L2)
      normValue = ::pressio::ops::norm2(R_);
    else if (normKind == ::pressio::Norm::L1)
      normValue = ::pressio::ops::norm1(R_);
    else
      throw std::runtime_error("Invalid norm kind for lspg unsteady residual policy");
  }

};//end class

}}} //end namespace pressio::rom::decorator
#endif  // ROM_DECORATORS_ROM_MASKED_RESIDUAL_POLICY_HPP_
