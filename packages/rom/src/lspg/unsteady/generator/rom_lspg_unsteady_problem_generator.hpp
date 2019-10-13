/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_problem_generator.hpp
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

#ifndef PRESSIO_ROM_LSPG_UNSTEADY_PROBLEM_GENERATOR_HPP_
#define PRESSIO_ROM_LSPG_UNSTEADY_PROBLEM_GENERATOR_HPP_

#include "rom_lspg_unsteady_problem_type_generator_default.hpp"
#include "rom_lspg_unsteady_problem_type_generator_masked.hpp"
#include "rom_lspg_unsteady_problem_type_generator_preconditioned.hpp"

namespace pressio{ namespace rom{

template <typename lspg_problem>
class LSPGUnsteadyProblemGenerator<lspg_problem>
  : public lspg_problem
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
  using typename lspg_problem::fom_velocity_data;
  using typename lspg_problem::ud_ops_t;

  using typename lspg_problem::lspg_matrix_t;
  using typename lspg_problem::fom_eval_velocity_policy_t;
  using typename lspg_problem::fom_apply_jac_policy_t;
  using typename lspg_problem::lspg_residual_policy_t;
  using typename lspg_problem::lspg_jacobian_policy_t;

  using typename lspg_problem::aux_stepper_t;
  using typename lspg_problem::lspg_stepper_t;

private:
  fom_eval_velocity_policy_t	rhsEv_;
  fom_apply_jac_policy_t	ajacEv_;
  fom_state_t			yFomRef_;
  fom_state_reconstr_t		fomStateReconstructor_;
  fom_velocity_t		rFomRef_;
  fom_states_data		fomStates_;
  fom_velocity_data		fomRhs_;
  lspg_matrix_t			romMat_;
  lspg_residual_policy_t	resObj_;
  lspg_jacobian_policy_t	jacObj_;

  /* here we use conditional type if auxiliary stepper is non-void,
   * otherwise we set it to a dummy type and we dont construct it */
  typename std::conditional<
    std::is_void<aux_stepper_t>::value,
    utils::impl::empty, aux_stepper_t>::type auxStepperObj_;

  // actual stepper object
  lspg_stepper_t		stepperObj_;

public:
  lspg_stepper_t & getStepperRef(){
    return stepperObj_;
  }

  const fom_state_reconstr_t & getFomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

public:

  /* - aux stepper NOT needed
   * - ud_ops_t != void
  */
  template <
    typename T = aux_stepper_t,
    typename T2 = ud_ops_t,
    typename ::pressio::mpl::enable_if_t<
      std::is_void<T>::value and
      !std::is_void<T2>::value
      > * = nullptr
  >
  LSPGUnsteadyProblemGenerator(const fom_t	 & appObj,
			       const fom_native_state_t & yFomRefNative,
			       decoder_t	 & decoder,
			       lspg_state_t	 & yROM,
			       scalar_t		 t0,
			       const T2		 & udOps)
    : rhsEv_{},
      ajacEv_{},
      yFomRef_(yFomRefNative),
      fomStateReconstructor_(yFomRef_, decoder),
      rFomRef_( rhsEv_.evaluate(appObj, yFomRef_, t0) ),
      fomStates_(yFomRef_, fomStateReconstructor_),
      fomRhs_(rFomRef_),
      romMat_(ajacEv_.evaluate(appObj, yFomRef_,
      			       decoder.getReferenceToJacobian(), t0)),
      resObj_(fomStates_, fomRhs_, rhsEv_, udOps),
      jacObj_(fomStates_, ajacEv_, romMat_, decoder, udOps),
      auxStepperObj_{},
      stepperObj_(yROM, appObj, resObj_, jacObj_)
  {
#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    std::cout << std::endl;
    std::cout << "LSPGProbGen" << std::endl;
    std::cout << "yFomRef_ " << yFomRef_.data() << std::endl;
    std::cout << "rFomRef_ " << rFomRef_.data() << std::endl;
    std::cout << "romMat_ " << romMat_.data() << std::endl;
    std::cout << std::endl;
#endif
  }


  /*
   * - aux stepper NOT needed
   * - ud_ops_t = void
  */
  template <
    typename T = aux_stepper_t,
    typename T2 = ud_ops_t,
  typename ::pressio::mpl::enable_if_t<
      std::is_void<T>::value and
      std::is_void<T2>::value
      > * = nullptr
  >
  LSPGUnsteadyProblemGenerator(const fom_t	 & appObj,
  			       const fom_native_state_t & yFomRefNative,
  			       decoder_t	 & decoder,
  			       lspg_state_t	 & yROM,
  			       scalar_t		 t0)
    : rhsEv_{},
      ajacEv_{},
      yFomRef_(yFomRefNative),
      fomStateReconstructor_(yFomRef_, decoder),
      rFomRef_( rhsEv_.evaluate(appObj, yFomRef_, t0) ),
      fomStates_(yFomRef_, fomStateReconstructor_),
      fomRhs_(rFomRef_),
      romMat_(ajacEv_.evaluate(appObj, yFomRef_, decoder.getReferenceToJacobian(), t0)),
      resObj_(fomStates_, fomRhs_, rhsEv_),
      jacObj_(fomStates_, ajacEv_, romMat_, decoder),
      auxStepperObj_{},
      stepperObj_(yROM, appObj, resObj_, jacObj_)
  {}


  /* sfinae here for when we need aux stepper */
  template <
    typename T = aux_stepper_t,
    typename ::pressio::mpl::enable_if_t<
      std::is_void<T>::value == false
      > * = nullptr
    >
  LSPGUnsteadyProblemGenerator(const fom_t	 & appObj,
  			       const fom_native_state_t & yFomRefNative,
  			       const decoder_t	 & decoder,
  			       lspg_state_t	 & yROM,
  			       scalar_t		 t0)
    : rhsEv_{},
      ajacEv_{},
      yFomRef_(yFomRefNative),
      fomStateReconstructor_(yFomRef_, decoder),
      rFomRef_( rhsEv_.evaluate(appObj, yFomRef_, t0) ),
      fomStates_(yFomRef_, fomStateReconstructor_),
      fomRhs_(rFomRef_),
      romMat_(ajacEv_.evaluate(appObj, yFomRef_,
  			       decoder.getReferenceToJacobian(), t0)),
      resObj_(fomStates_, fomRhs_, rhsEv_),
      jacObj_(fomStates_, ajacEv_, romMat_, decoder),
      auxStepperObj_(yROM, appObj, resObj_, jacObj_),
      stepperObj_(yROM, appObj, resObj_, jacObj_, auxStepperObj_)
  {}

};

}}//end namespace pressio::rom
#endif
