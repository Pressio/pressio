/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_masked_problem_continuous_time_api.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_ROM_LSPG_UNSTEADY_MASKED_PROBLEM_CONTINUOUS_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_ROM_LSPG_UNSTEADY_MASKED_PROBLEM_CONTINUOUS_TIME_API_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template<typename ...Args>
class MaskedProblemContinuousTimeApi
{
public:
  using this_t = MaskedProblemContinuousTimeApi<Args...>;
  using traits = ::pressio::rom::details::traits<this_t>;

  using fom_system_t		= typename traits::fom_system_t;
  using scalar_t		= typename traits::scalar_t;
  using fom_native_state_t	= typename traits::fom_native_state_t;
  using fom_state_t		= typename traits::fom_state_t;
  using fom_velocity_t		= typename traits::fom_velocity_t;
  using lspg_state_t		= typename traits::lspg_state_t;
  using decoder_t		= typename traits::decoder_t;
  using fom_state_reconstr_t	= typename traits::fom_state_reconstr_t;
  using fom_states_manager_t	= typename traits::fom_states_manager_t;
  using masker_t		= typename traits::masker_t;
  using ud_ops_t		= typename traits::ud_ops_t;
  using lspg_matrix_t		= typename traits::lspg_matrix_t;
  using residual_policy_t	= typename traits::residual_policy_t;
  using jacobian_policy_t	= typename traits::jacobian_policy_t;
  using aux_stepper_t		= typename traits::aux_stepper_t;
  using stepper_t		= typename traits::stepper_t;

private:
  const fom_state_t		fomStateReference_;
  const fom_velocity_t		fomVelocityRef_;
  const fom_state_reconstr_t	fomStateReconstructor_;
  fom_states_manager_t		fomStatesMngr_;
  lspg_matrix_t			jPhiMatrix_;
  residual_policy_t	residualPolicy_;
  jacobian_policy_t	jacobianPolicy_;

  /* here we use conditional type if auxiliary stepper is non-void,
   * otherwise we set it to a dummy type and we dont construct it */
  typename std::conditional<
    std::is_void<aux_stepper_t>::value,
    ::pressio::utils::impl::empty, aux_stepper_t
    >::type auxStepperObj_ = {};

  // actual stepper object
  stepper_t  stepperObj_;

public:
  stepper_t & stepperRef(){
    return stepperObj_;
  }

  const fom_native_state_t & currentFomState() const{
    return *fomStatesMngr_.currentFomStateCRef().data();
  }

  const fom_state_reconstr_t & fomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

public:
  MaskedProblemContinuousTimeApi() = delete;
  MaskedProblemContinuousTimeApi(const MaskedProblemContinuousTimeApi &) = default;
  MaskedProblemContinuousTimeApi & operator=(const MaskedProblemContinuousTimeApi &) = default;
  MaskedProblemContinuousTimeApi(MaskedProblemContinuousTimeApi &&) = default;
  MaskedProblemContinuousTimeApi & operator=(MaskedProblemContinuousTimeApi &&) = default;
  ~MaskedProblemContinuousTimeApi() = default;

  /* specialize for:
   * - the fom_system_t is regular c++
   * - aux stepper is NOT needed (e.g. for BDF1)
   * - ud_ops_t == void
   */
  template <
    typename _fom_system_t = fom_system_t,
    typename _aux_stepper_t = aux_stepper_t,
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      !::pressio::ops::predicates::is_object_pybind<_fom_system_t>::value and
      std::is_void<_aux_stepper_t>::value and
      std::is_void<_ud_ops_t>::value,
      int > = 0
    >
  MaskedProblemContinuousTimeApi(const _fom_system_t	& fomSystemObj,
				 const decoder_t & decoder,
				 const lspg_state_t & romStateIn,
				 const fom_native_state_t & fomStateReferenceNative,
				 const masker_t & maskerObj)
    : fomStateReference_(fomStateReferenceNative),
      fomVelocityRef_(fomSystemObj.createVelocity()),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      //
      jPhiMatrix_(fomSystemObj.createApplyJacobianResult
		  ( *decoder.jacobianCRef().data() )),
      //
      // here we pass a fom velocity object to the residual policy to
      // use it to initialize the residual data
      // since the lspg residual is of same type and size of the fom velocity
      // (this is true w and w/o hyperreduction)
      residualPolicy_(maskerObj, fomVelocityRef_, fomStatesMngr_),
      jacobianPolicy_(maskerObj, fomStatesMngr_, jPhiMatrix_, decoder),
      auxStepperObj_{},
      stepperObj_(romStateIn, fomSystemObj, residualPolicy_, jacobianPolicy_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }

  /* specialize for:
   * - the fom_system_t is regular c++
   * - aux stepper is needed (e.g. for BDF2)
   * - ud_ops_t == void
   */
  template <
    typename _fom_system_t = fom_system_t,
    typename _aux_stepper_t = aux_stepper_t,
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      !::pressio::ops::predicates::is_object_pybind<_fom_system_t>::value and
      !std::is_void<_aux_stepper_t>::value and
      std::is_void<_ud_ops_t>::value,
      int > = 0
    >
  MaskedProblemContinuousTimeApi(const _fom_system_t & fomSystemObj,
				 const decoder_t & decoder,
				 const lspg_state_t & romStateIn,
				 const fom_native_state_t & fomStateReferenceNative,
				 const masker_t & maskerObj)
    : fomStateReference_(fomStateReferenceNative),
      fomVelocityRef_(fomSystemObj.createVelocity()),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomStateReference_),
      //
      jPhiMatrix_(fomSystemObj.createApplyJacobianResult
		  (*decoder.jacobianCRef().data())),
      //
      // here we pass a fom velocity object to the residual policy to
      // use it to initialize the residual data
      // since the lspg residual is of same type and size of the fom velocity
      // (this is true w and w/o hyperreduction)
      residualPolicy_(maskerObj, fomVelocityRef_, fomStatesMngr_),
      jacobianPolicy_(maskerObj, fomStatesMngr_, jPhiMatrix_, decoder),
      auxStepperObj_(romStateIn, fomSystemObj, residualPolicy_, jacobianPolicy_),
      stepperObj_(romStateIn, fomSystemObj, residualPolicy_,
		  jacobianPolicy_, auxStepperObj_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }

};

}}}}}
#endif  // ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_ROM_LSPG_UNSTEADY_MASKED_PROBLEM_CONTINUOUS_TIME_API_HPP_
