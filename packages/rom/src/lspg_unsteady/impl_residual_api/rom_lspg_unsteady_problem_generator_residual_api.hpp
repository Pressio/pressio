/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_problem_generator_residual_api.hpp
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

#ifndef PRESSIO_ROM_LSPG_UNSTEADY_PROBLEM_GENERATOR_RESIDUAL_API_HPP_
#define PRESSIO_ROM_LSPG_UNSTEADY_PROBLEM_GENERATOR_RESIDUAL_API_HPP_

#include "rom_lspg_unsteady_problem_type_generator_default_residual_api.hpp"

namespace pressio{ namespace rom{ namespace impl{

template <
  template <::pressio::ode::ImplicitEnum, class, class, class ...> class lspg_t,
  ::pressio::ode::ImplicitEnum odeName,
  typename fom_t,
  typename lspg_state_t,
  typename ...Args
  >
class LSPGUnsteadyProblemGeneratorResidualAPI
{
  // ::pressio::rom::meta::model_meets_residual_api_for_unsteady_lspg<
  //     typename lspg_problem_t::fom_t
  //     >::value

// public:
//   using fom_t			= typename lspg_problem_t::fom_t;
//   using scalar_t		= typename lspg_problem_t::scalar_t;
  // using fom_native_state_t	= typename lspg_problem_t::fom_native_state_t;
  // using fom_state_t		= typename lspg_problem_t::fom_state_t;
  // using fom_velocity_t		= typename lspg_problem_t::fom_velocity_t;
  // using lspg_state_t		= typename lspg_problem_t::lspg_state_t;
  // using decoder_t		= typename lspg_problem_t::decoder_t;
  // using fom_state_reconstr_t	= typename lspg_problem_t::fom_state_reconstr_t;
  // using fom_states_data		= typename lspg_problem_t::fom_states_data;
  // using ud_ops_t		= typename lspg_problem_t::ud_ops_t;
  // using lspg_matrix_t		= typename lspg_problem_t::lspg_matrix_t;
  // using fom_eval_velo_policy_t	= typename lspg_problem_t::fom_eval_velocity_policy_t;
  // using fom_apply_jac_policy_t	= typename lspg_problem_t::fom_apply_jac_policy_t;
  // using lspg_residual_policy_t	= typename lspg_problem_t::lspg_residual_policy_t;
  // using lspg_jacobian_policy_t	= typename lspg_problem_t::lspg_jacobian_policy_t;
  // using aux_stepper_t		= typename lspg_problem_t::aux_stepper_t;
  // using lspg_stepper_t		= typename lspg_problem_t::lspg_stepper_t;

// private:
//   fom_eval_velo_policy_t	veloQuerier_;
//   fom_apply_jac_policy_t	applyJacobQuerier_;
//   fom_state_t			fomStateReference_;
//   fom_state_reconstr_t		fomStateReconstructor_;
//   fom_velocity_t		fomVelocityRef_;
//   fom_states_data		fomStates_;
//   lspg_matrix_t			jPhiMatrix_;
//   lspg_residual_policy_t	residualPolicy_;
//   lspg_jacobian_policy_t	jacobianPolicy_;

//   /* here we use conditional type if auxiliary stepper is non-void,
//    * otherwise we set it to a dummy type and we dont construct it */
//   typename std::conditional<
//     std::is_void<aux_stepper_t>::value,
//     utils::impl::empty, aux_stepper_t
//     >::type auxStepperObj_;

//   // actual stepper object
//   lspg_stepper_t  stepperObj_;


// public:
//   lspg_stepper_t & getStepperRef(){
//     return stepperObj_;
//   }

//   const fom_state_reconstr_t & getFomStateReconstructorCRef() const{
//     return fomStateReconstructor_;
//   }

// public:

//   /*- aux stepper = void
//    * - ud_ops_t = void */
//   template <
//     typename _aux_stepper_t = aux_stepper_t,
//     typename _ud_ops_t = ud_ops_t,
//     typename ::pressio::mpl::enable_if_t<
//       std::is_void<_aux_stepper_t>::value and
//       std::is_void<_ud_ops_t>::value
//       > * = nullptr
//   >
//   LSPGUnsteadyProblemGenerator(const fom_t	 & appObj,
//   			       const fom_native_state_t & fomStateReferenceNative,
//   			       decoder_t	 & decoder,
//   			       lspg_state_t	 & yROM,
//   			       scalar_t		 t0)
//     : veloQuerier_{},
//       applyJacobQuerier_{},
//       fomStateReference_(fomStateReferenceNative),
//       fomStateReconstructor_(fomStateReference_, decoder),
//       fomVelocityRef_( veloQuerier_.evaluate(appObj, fomStateReference_, t0) ),
//       fomStates_(fomStateReference_, fomStateReconstructor_),
//       jPhiMatrix_(applyJacobQuerier_.evaluate(appObj, fomStateReference_, decoder.getReferenceToJacobian(), t0)),
//       // here we pass a fom velocity object to the residual policy to use it to initialize the residual data
//       // since the lspg residual is of same type and size of the fom velocity (this is true w and w/o hyperreduction)
//       residualPolicy_(fomVelocityRef_, fomStates_, veloQuerier_),
//       jacobianPolicy_(fomStates_, applyJacobQuerier_, jPhiMatrix_, decoder),
//       auxStepperObj_{},
//       stepperObj_(yROM, appObj, residualPolicy_, jacobianPolicy_)
//   {}

};

}}}//end namespace pressio::rom::impl
#endif
