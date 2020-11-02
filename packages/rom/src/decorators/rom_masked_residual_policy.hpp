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

template <typename masker_t, typename maskable_policy>
class MaskedResidualPolicy : public maskable_policy
{
  using typename maskable_policy::residual_t;
  mutable residual_t R_;
  std::reference_wrapper<const masker_t> maskerObj_;

public:
  MaskedResidualPolicy() = delete;
  MaskedResidualPolicy(const MaskedResidualPolicy &) = default;
  MaskedResidualPolicy & operator=(const MaskedResidualPolicy &) = default;
  MaskedResidualPolicy(MaskedResidualPolicy &&) = default;
  MaskedResidualPolicy & operator=(MaskedResidualPolicy &&) = default;
  ~MaskedResidualPolicy() = default;

  template <typename fom_system_t, typename ... Args>
  MaskedResidualPolicy(const masker_t & maskerObj,
		       const fom_system_t & fomObj,
		       Args && ... args)
    : maskable_policy(std::forward<Args>(args)...),
      R_(maskable_policy::create(fomObj)),
      maskerObj_(maskerObj)
  {}

public:
  template <typename fom_system_t>
  residual_t create(const fom_system_t & systemObj) const
  {
    return residual_t(maskerObj_.get().createApplyMaskResult(*R_.data()));
  }

  //-------------------------------
  // unsteady case
  //-------------------------------
  template <
    typename stepper_tag,
    typename state_t,
    typename prev_states_t,
    typename fom_system_t,
    typename time_type
    >
  void compute(const state_t & state,
	       const prev_states_t & prevStates,
	       const fom_system_t & systemObj,
	       const time_type & time,
	       const time_type & dt,
	       const ::pressio::ode::types::step_t & step,
	       residual_t & R) const
  {
    maskable_policy::template compute<stepper_tag>
      (state, prevStates, systemObj, time, dt, step, R_);

    maskerObj_.get().applyMask(*R_.data(), time, *R.data());
  }

  //-------------------------------
  // steady case
  //-------------------------------
  template <
    typename state_t,
    typename fom_system_t
    >
  void compute(const state_t & state,
  	       residual_t & R,
  	       const fom_system_t & systemObj) const
  {
    maskable_policy::compute(state, R_, systemObj);
    maskerObj_.get().applyMask(*R_.data(), *R.data());
  }

};

}}} //end namespace pressio::rom::decorator
#endif  // ROM_DECORATORS_ROM_MASKED_RESIDUAL_POLICY_HPP_
