/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_problem_generator.hpp
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

#ifndef PRESSIO_ROM_LSPG_STEADY_PROBLEM_GENERATOR_HPP_
#define PRESSIO_ROM_LSPG_STEADY_PROBLEM_GENERATOR_HPP_

#include "rom_lspg_steady_type_generator_default.hpp"

namespace pressio{ namespace rom{

template <typename lspg_problem>
class LSPGSteadyProblemGenerator<
  lspg_problem > : lspg_problem
{

public:
  using typename lspg_problem::fom_t;
  using typename lspg_problem::scalar_t;
  using typename lspg_problem::fom_native_state_t;
  using typename lspg_problem::fom_state_t;
  using typename lspg_problem::fom_velocity_t;

  using typename lspg_problem::lspg_state_t;
  using typename lspg_problem::decoder_t;
  using typename lspg_problem::fom_state_reconstr_t;
  using typename lspg_problem::fom_states_data;

  using typename lspg_problem::lspg_matrix_t;
  using typename lspg_problem::fom_eval_rhs_policy_t;
  using typename lspg_problem::fom_apply_jac_policy_t;
  using typename lspg_problem::lspg_residual_policy_t;
  using typename lspg_problem::lspg_jacobian_policy_t;
  using typename lspg_problem::lspg_system_t;

private:
  fom_eval_rhs_policy_t		rhsQuerier_;
  fom_apply_jac_policy_t	applyJacobQuerier_;
  fom_state_t			fomStateReference_;
  fom_state_reconstr_t		fomStateReconstructor_;
  fom_velocity_t		fomVelocityRef_;
  fom_states_data		fomStates_;
  lspg_matrix_t			jPhiMatrix_;
  lspg_residual_policy_t	residualPolicy_;
  lspg_jacobian_policy_t	jacobianPolicy_;
  lspg_system_t			systemObj_;

public:
  lspg_system_t & getSystemRef(){
    return systemObj_;
  }

  const fom_state_reconstr_t & getFomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

public:
  LSPGSteadyProblemGenerator(const fom_t	& appObj,
			     const fom_native_state_t & yFomRefNative,
			     const decoder_t	& decoder,
			     lspg_state_t	& yROM)
    : rhsQuerier_{},
      applyJacobQuerier_{},
      fomStateReference_(yFomRefNative),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomVelocityRef_( rhsQuerier_.evaluate(appObj, fomStateReference_) ),
      fomStates_(fomStateReference_, fomStateReconstructor_),
      jPhiMatrix_(applyJacobQuerier_.evaluate(appObj, fomStateReference_,
					  decoder.getReferenceToJacobian())),
      residualPolicy_(fomVelocityRef_, fomStates_, rhsQuerier_),
      jacobianPolicy_(fomStates_, applyJacobQuerier_, jPhiMatrix_, decoder),
      systemObj_(appObj, residualPolicy_, jacobianPolicy_)
  {}

};

}}//end namespace pressio::rom
#endif
