/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_policy_residual.hpp
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

#ifndef ROM_IMPL_TO_REVISE_LSPG_UNSTEADY_RJ_POLICY_MASKED_HPP_
#define ROM_IMPL_TO_REVISE_LSPG_UNSTEADY_RJ_POLICY_MASKED_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class IndVarType,
  class ReducedStateType,
  class LspgResidualType,
  class LspgJacobianType,
  class TrialSpaceType,
  class FomSystemType,
  class RhsMaskerType,
  class JacobianActionMaskerType,
  class HypRedUpdaterType
  >
class LspgUnsteadyResidualJacobianPolicyMasked
{
  using basis_type = typename TrialSpaceType::basis_type;

  // deduce the unmasked types
  using unmasked_fom_rhs_type = typename FomSystemType::right_hand_side_type;
  using unmasked_fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_type const &>()));

  // deduce the masked types
  using masked_fom_rhs_type = typename RhsMaskerType::result_type;
  using masked_fom_jac_action_result_type = typename JacobianActionMaskerType::result_type;

public:
  // required
  using independent_variable_type = IndVarType;
  using state_type    = ReducedStateType;
  using residual_type = LspgResidualType;
  using jacobian_type = LspgJacobianType;

public:
  LspgUnsteadyResidualJacobianPolicyMasked() = delete;
  LspgUnsteadyResidualJacobianPolicyMasked(const TrialSpaceType & trialSpace,
					   const FomSystemType & fomSystem,
					   LspgFomStatesManager<TrialSpaceType> & fomStatesManager,
					   const RhsMaskerType & rhsMasker,
					   const JacobianActionMaskerType & jaMasker,
					   const HypRedUpdaterType & hrUpdater)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomStatesManager_(fomStatesManager),
      rhsMasker_(rhsMasker),
      jaMasker_(jaMasker),
      hypRedUpdater_(hrUpdater),
      fomStateHelperInstance_(trialSpace.createFullState()),
      unMaskedFomRhs_(fomSystem.createRightHandSide()),
      unMaskedFomJacAction_(fomSystem.createApplyJacobianResult(trialSpace_.get().viewBasis()))
  {}

public:
  state_type createState() const{
    // this needs to create an instance of the reduced state
    return trialSpace_.get().createReducedState();
  }

  residual_type createResidual() const{
    // for lspg, a residual instance can be contructed from the masked rhs
    auto tmp = fomSystem_.get().createResidual();
    ::pressio::ops::set_zero(tmp);
    return rhsMasker_.get().createApplyMaskResult(tmp);
  }

  jacobian_type createJacobian() const{
    const auto phi = trialSpace_.get().viewBasis();
    auto tmp = fomSystem_.get().createApplyJacobianResult(phi);
    ::pressio::ops::set_zero(J);
    return jaMasker_.get().createApplyMaskResult(tmp);
  }

  template <class StencilStatesContainerType, class StencilRhsContainerType>
  void operator()(::pressio::ode::StepScheme odeSchemeName,
		  const state_type & predictedReducedState,
		  const StencilStatesContainerType & reducedStatesStencilManager,
		  StencilRhsContainerType & fomRhsStencilManger,
		  const ::pressio::ode::StepEndAt<IndVarType> & rhsEvaluationTime,
		  ::pressio::ode::StepCount step,
		  const ::pressio::ode::StepSize<IndVarType> & dt,
		  residual_type & R,
		  jacobian_type & J,
		  bool computeJacobian) const
  {

    if (odeSchemeName == ::pressio::ode::StepScheme::BDF1){
      (*this).template compute_impl_bdf
	(odeSchemeName, predictedReducedState, reducedStatesStencilManager,
	 fomRhsStencilManger, rhsEvaluationTime.get(),
	 dt.get(), step.get(), R, J, computeJacobian);
    }

    else if (odeSchemeName == ::pressio::ode::StepScheme::BDF2){

    }

    else{
      throw std::runtime_error("Only BDF1 currently impl for default unstedy LSPG");
    }
  }

private:
  template <class StencilStatesContainerType, class StencilRhsContainerType>
  void compute_impl_bdf(::pressio::ode::StepScheme odeSchemeName,
			const state_type & predictedReducedState,
			const StencilStatesContainerType & reducedStatesStencilManager,
			StencilRhsContainerType & /*unused*/,
			const IndVarType & rhsEvaluationTime,
			const IndVarType & dt,
			const typename ::pressio::ode::StepCount::value_type & step,
			residual_type & R,
			jacobian_type & J,
			bool computeJacobian) const
  {

    /* the FOM state for the prediction has to be always recomputed
       regardless of whether the currentStepNumber changes since
       we might be inside a subiteration of the non-linear solve
       where the time step does not change but the predicted state does
     */
    fomStatesManager_.get().reconstructAtWithoutStencilUpdate(predictedReducedState,
							      ::pressio::ode::nPlusOne());
    const auto & fomStateAt_np1 = fomStatesManager_(::pressio::ode::nPlusOne());

    /* previous FOM states should only be recomputed when the time step changes.
       The method below does not recompute all previous states, but only
       recomputes the n-th state and updates/shifts back all the other
       FOM states stored. */
    if (stepTracker_ != step){
      const auto & lspgStateAt_n = reducedStatesStencilManager(::pressio::ode::n());
      fomStatesManager_.get().reconstructAtWithStencilUpdate(lspgStateAt_n,
							     ::pressio::ode::n());
      stepTracker_ = step;
    }

    /*
      BDF1 residual :

         R(y_n+1) = cnp1*y_n+1 + cn*y_n + cf*f(t_n+1, y_n+1)

      which we do follows:
      1. R = masked(f(t_n+1, y_n+1))
      2. fomStateHelpInstance_ = cnp1*y_np1 + cn*y_n
      3. call the hypRedUpdater to handle the rest

      for BDF2, we have:

         R(y_n+1) = cnp1*y_n+1 + cn*y_n + cnm1*y_n-1 + cf*f(t_n+1, y_n+1)

      so only difference is step 2
      2. fomStateHelpInstance_ = cnp1*y_np1 + cn*y_n + cnm1*y_n-1
    */

    // step 1 (same for both)
    fomSystem_.get().rightHandSide(fomStateAt_np1, rhsEvaluationTime, unMaskedFomRhs_);
    rhsMasker_(unMaskedFomRhs_, R);

    // step 2 is different
    const auto & fomStateAt_n = fomStatesManager_(::pressio::ode::n());
    using fom_state_type = typename FomSystemType::state_type;
    using sc_t = typename ::pressio::Traits<fom_state_type>::scalar_type;
    constexpr auto zero = ::pressio::utils::Constants<sc_t>::zero();
    constexpr auto one = ::pressio::utils::Constants<sc_t>::one();

    sc_t cf = {};
    if (odeSchemeName == ::pressio::ode::StepScheme::BDF1)
    {
      constexpr auto cnp1 = ::pressio::ode::constants::bdf1<sc_t>::c_np1_;
      constexpr auto cn   = ::pressio::ode::constants::bdf1<sc_t>::c_n_;
      cf = ::pressio::ode::constants::bdf1<sc_t>::c_f_ * dt;

      ::pressio::ops::update(fomStateHelperInstance_, zero,
			     fomStateAt_np1, cnp1,
			     fomStateAt_n, cn);
    }
    else if (odeSchemeName == ::pressio::ode::StepScheme::BDF2)
    {
      constexpr auto cnp1 = ::pressio::ode::constants::bdf2<sc_t>::c_np1_;
      constexpr auto cn   = ::pressio::ode::constants::bdf2<sc_t>::c_n_;
      constexpr auto cnm1 = ::pressio::ode::constants::bdf2<sc_t>::c_nm1_;
      cf = ::pressio::ode::constants::bdf2<sc_t>::c_f_ * dt;

      const auto & fomStateAt_nm1 = fomStatesManager_(::pressio::ode::nMinusOne());
      ::pressio::ops::update(fomStateHelperInstance_, zero,
			     fomStateAt_np1, cnp1,
			     fomStateAt_n, cn,
			     fomStateAt_nm1, cnm1);
    }

    // step 3
    hypRedUpdater_.get().updateSampleMeshOperandWithStencilMeshOne
      (R, cf, fomStateHelperInstance_, one);


    if (computeJacobian)
    {
      /* for BDF1 and BDF2 this means:

	 lspgJac = decoderJac + dt*coeff*masked(J*decoderJac)

	 where J is the d(fomrhs)/dy and coeff depends on the scheme.
      */

      // first, store J*phi into J
      const auto & phi = trialSpace_.get().viewBasis();
      fomSystem_.get().applyJacobian(fomStateAt_np1, phi, rhsEvaluationTime, unMaskedFomJacAction_);
      jaMasker_(unMaskedFomJacAction_, J);

      using sc_t = typename ::pressio::Traits<typename TrialSpaceType::basis_type>::scalar_type;
      const auto one = ::pressio::utils::Constants<sc_t>::one();
      IndVarType cf = dt;
      if (odeSchemeName == ::pressio::ode::StepScheme::BDF1){
	cf *= ::pressio::ode::constants::bdf1<sc_t>::c_f_;
      }
      else if (odeSchemeName == ::pressio::ode::StepScheme::BDF2){
	cf *= ::pressio::ode::constants::bdf2<sc_t>::c_f_;
      }

      hypRedUpdater_.get().updateSampleMeshOperandWithStencilMeshOne(J, cf, phi, one);
    }
  }

private:
  using raw_step_type = typename ::pressio::ode::StepCount::value_type;
  static_assert(std::is_signed<raw_step_type>::value, "");
  // storedStep is used to keep track of which step we are at.
  // used to decide if we need to update/recompute the previous FOM states or not.
  // To avoid recomputing previous FOM states if we are not in a new time step.
  mutable raw_step_type stepTracker_ = -1;

  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<LspgFomStatesManager<TrialSpaceType>> fomStatesManager_;
  mutable typename FomSystemType::state_type fomStateHelperInstance_;
  std::reference_wrapper<const HypRedUpdaterType> hypRedUpdater_;

  // masker
  std::reference_wrapper<const RhsMaskerType> rhsMasker_;
  std::reference_wrapper<const JacobianActionMaskerType> jaMasker_;
  // UNMASKED objects
  mutable unmasked_fom_rhs_type unMaskedFomRhs_;
  mutable unmasked_fom_jac_action_result_type unMaskedFomJacAction_;
};

}}}
#endif  // ROM_IMPL_TO_REVISE_LSPG_UNSTEADY_RJ_POLICY_MASKED_HPP_
