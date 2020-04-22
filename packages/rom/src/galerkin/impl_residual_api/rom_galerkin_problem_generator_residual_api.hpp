/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_problem_generator_residual_api.hpp
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

#ifndef PRESSIO_ROM_EXP_GALERKIN_PROBLEM_GENERATOR_RESIDUAL_API_IMPL_HPP_
#define PRESSIO_ROM_EXP_GALERKIN_PROBLEM_GENERATOR_RESIDUAL_API_IMPL_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <
  template <class ...> class galerkin_type,
  typename fom_type,
  typename stepper_tag,
  typename rom_state_type,
  typename rom_jacobian_type,
  typename ...Args
  >
class ProblemGeneratorResidualApi
{

public:
  using problem_t	= galerkin_type<fom_type, stepper_tag, rom_state_type, rom_jacobian_type, Args...>;

  using fom_t		= typename problem_t::fom_t;
  using scalar_t	= typename problem_t::scalar_t;
  using fom_nat_state_t = typename problem_t::fom_native_state_t;
  using fom_state_t	= typename problem_t::fom_state_t;

  using decoder_t	= typename problem_t::decoder_t;
  using fom_reconstr_t	= typename problem_t::fom_state_reconstr_t;
  using fom_states_data	= typename problem_t::fom_states_data;
  using ud_ops_t	= typename problem_t::ud_ops_t;

  using rom_state_t	= typename problem_t::rom_state_t;

  using fom_r_querier_policy_t = typename problem_t::fom_residual_querier_policy_t;
  using fom_apply_jac_policy_t = typename problem_t::fom_apply_jac_policy_t;

  using residual_policy_t = typename problem_t::residual_policy_t;
  using jacobian_policy_t = typename problem_t::jacobian_policy_t;

  using stepper_t = typename problem_t::stepper_t;

private:
  ::pressio::ode::types::step_t step0_;
  scalar_t t0_;
  scalar_t dt0_;
  fom_r_querier_policy_t residualQuerier_;
  fom_apply_jac_policy_t applyJacobQuerier_;
  fom_state_t fomStateReference_;
  fom_reconstr_t fomStateReconstructor_;
  fom_states_data fomStates_;

  residual_policy_t residualPolicy_;
  jacobian_policy_t jacobianPolicy_;
  stepper_t stepperObj_;

public:
  stepper_t & getStepperRef(){
    return stepperObj_;
  }

  const fom_reconstr_t & getFomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

public:
  ProblemGeneratorResidualApi() = delete;

  ProblemGeneratorResidualApi(const fom_t & appObj,
			      const fom_nat_state_t & fomStateReferenceNative,
			      decoder_t	 & decoder,
			      rom_state_t & yROM,
			      scalar_t	t0)
    : step0_{},
      t0_{t0},
      dt0_{},
      residualQuerier_{},
      applyJacobQuerier_{},
      fomStateReference_(fomStateReferenceNative),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomStates_(fomStateReconstructor_, fomStateReference_),
      // policies
      residualPolicy_(fomStates_, residualQuerier_, decoder, appObj),
      jacobianPolicy_(fomStates_, applyJacobQuerier_, decoder, appObj),
      // stepper
      stepperObj_(yROM, appObj, residualPolicy_, jacobianPolicy_)
  {}

};

}}}}//end namespace
#endif
