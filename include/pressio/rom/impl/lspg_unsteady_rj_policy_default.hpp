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

#ifndef ROM_IMPL_ROM_LSPG_UNSTEADY_DEFAULT_RJ_POLICY_HPP_
#define ROM_IMPL_ROM_LSPG_UNSTEADY_DEFAULT_RJ_POLICY_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class IndVarType,
  class ReducedStateType,
  class LspgResidualType,
  class LspgJacobianType,
  class TrialSubspaceType,
  class FomSystemType
  >
class LspgUnsteadyResidualJacobianPolicy
{
public:
  // required
  using independent_variable_type = IndVarType;
  using state_type    = ReducedStateType;
  using residual_type = LspgResidualType;
  using jacobian_type = LspgJacobianType;

public:
  LspgUnsteadyResidualJacobianPolicy() = delete;
  LspgUnsteadyResidualJacobianPolicy(const TrialSubspaceType & trialSubspace,
				     const FomSystemType & fomSystem,
				     LspgFomStatesManager<TrialSubspaceType> & fomStatesManager)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomStatesManager_(fomStatesManager)
  {}

public:
  state_type createState() const{
    // this needs to create an instance of the reduced state
    return trialSubspace_.get().createReducedState();
  }

  residual_type createResidual() const{
    residual_type R(fomSystem_.get().createRightHandSide());
    return R;
  }

  jacobian_type createJacobian() const{
    const auto phi = trialSubspace_.get().basisOfTranslatedSpace();
    return fomSystem_.get().createApplyJacobianResult(phi);
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
	<ode::BDF1>(predictedReducedState, reducedStatesStencilManager,
		    fomRhsStencilManger, rhsEvaluationTime.get(),
		    dt.get(), step.get(), R, J, computeJacobian);
    }

    else if (odeSchemeName == ::pressio::ode::StepScheme::BDF2){
      (*this).template compute_impl_bdf
	<ode::BDF2>(predictedReducedState, reducedStatesStencilManager,
		    fomRhsStencilManger, rhsEvaluationTime.get(),
		    dt.get(), step.get(), R, J, computeJacobian);
    }

    else{
      throw std::runtime_error("Only BDF1 currently impl for default unstedy LSPG");
    }
  }

private:
  template <class OdeTag, class StencilStatesContainerType, class StencilRhsContainerType>
  void compute_impl_bdf(const state_type & predictedReducedState,
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

    //
    // always compute residual
    //
    fomSystem_.get().rightHandSide(fomStateAt_np1, rhsEvaluationTime, R);
    // default lspg does not do anything special, so I can use the
    // ode functions for computing the discrete residual
    ::pressio::ode::impl::discrete_residual(OdeTag(), fomStateAt_np1,
					    R, fomStatesManager_.get(), dt);

    //
    // deal with jacobian if needed
    //
    if (computeJacobian){
      // lspg Jac looks something like: lspgJac = decoderJac + dt*coeff*d(fomrhs)/dy*phi
      // where J is the d(fomrhs)/dy and coeff depends on the scheme.

      // first, store J*phi into J
      const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
      fomSystem_.get().applyJacobian(fomStateAt_np1, phi, rhsEvaluationTime, J);

      // second, we just need to update J properly
      using basis_sc_t = typename ::pressio::Traits<
	typename TrialSubspaceType::basis_matrix_type>::scalar_type;
      const auto one = ::pressio::utils::Constants<basis_sc_t>::one();
      IndVarType factor = {};
      if (std::is_same<OdeTag, ode::BDF1>::value){
	factor = dt*::pressio::ode::constants::bdf1<IndVarType>::c_f_;
      }
      else if (std::is_same<OdeTag, ode::BDF2>::value){
	factor = dt*::pressio::ode::constants::bdf2<IndVarType>::c_f_;
      }

      ::pressio::ops::update(J, factor, phi, one);
    }

  }

private:
  using raw_step_type = typename ::pressio::ode::StepCount::value_type;
  static_assert(std::is_signed<raw_step_type>::value, "");
  // storedStep is used to keep track of which step we are at.
  // used to decide if we need to update/recompute the previous FOM states or not.
  // To avoid recomputing previous FOM states if we are not in a new time step.
  mutable raw_step_type stepTracker_ = -1;

  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<LspgFomStatesManager<TrialSubspaceType>> fomStatesManager_;
};

}}}
#endif  // ROM_IMPL_ROM_LSPG_UNSTEADY_POLICY_RESIDUAL_HPP_


  // template <
  //   class LspgStateType,
  //   class LspgStencilStatesContainerType,
  //   class LspgStencilVelocitiesContainerType,
  //   class ScalarType>
  // void operator()(::pressio::ode::StepScheme name,
  // 		  const LspgStateType & lspgState,
  // 		  const LspgStencilStatesContainerType & lspgStencilStates,
  // 		  LspgStencilVelocitiesContainerType & lspgStencilVelocities,
  // 		  const ScalarType & t_np1,
  // 		  const ScalarType & dt,
  // 		  const ::pressio::ode::step_count_type & currentStepNumber,
  // 		  residual_type & lspgResidual) const
  // {

  //   if (name == ::pressio::ode::StepScheme::BDF1){
  //     (*this).template compute_impl_bdf<ode::BDF1>
  // 	(lspgState, lspgStencilStates, lspgStencilVelocities,
  // 	 t_np1, dt, currentStepNumber, lspgResidual);
  //   }
  //   else if (name == ::pressio::ode::StepScheme::BDF2){
  //     (*this).template compute_impl_bdf<ode::BDF2>
  // 	(lspgState, lspgStencilStates, lspgStencilVelocities,
  // 	 t_np1, dt, currentStepNumber, lspgResidual);
  //   }
  //   else if (name == ::pressio::ode::StepScheme::CrankNicolson){
  //     (*this).compute_impl_cn(lspgState, lspgStencilStates,lspgStencilVelocities,
  // 			      t_np1, dt, currentStepNumber, lspgResidual);
  //   }
  // }

//   template <
//     class LspgStateType,
//     class LspgStencilStatesContainerType,
//     class LspgStencilVelocitiesContainerType,
//     class ScalarType
//     >
//   void compute_impl_cn(const LspgStateType & lspgState,
// 		       const LspgStencilStatesContainerType & lspgStencilStates,
// 		       // for CN, stencilVelocities holds f_n+1 and f_n
// 		       LspgStencilVelocitiesContainerType & lspgStencilVelocities,
// 		       const ScalarType & t_np1,
// 		       const ScalarType & dt,
// 		       const ::pressio::ode::step_count_type & currentStepNumber,
// 		       residual_type & lspgResidual) const
//   {
//     PRESSIOLOG_DEBUG("residual policy with compute_cn_impl");

//     fomStatesMngr_.get().reconstructAt(lspgState, ::pressio::ode::nPlusOne());

//     if (storedStep_ != currentStepNumber){
//       const auto & lspgStateAt_n = lspgStencilStates(::pressio::ode::n());
//       fomStatesMngr_.get().reconstructAtAndUpdatePrevious(lspgStateAt_n,
// 							  ::pressio::ode::n());
//       storedStep_ = currentStepNumber;

//       // if the step changed, I need to compute f(y_n, t_n)
//       const auto tn = t_np1-dt;
//       auto & f_n = lspgStencilVelocities(::pressio::ode::n());
//       const auto & fomState_n = fomStatesMngr_(::pressio::ode::n());
//       fomSystem_.get().velocity(fomState_n, tn, f_n);
//     }

//     // always compute f(y_n+1, t_n+1)
//     const auto & fomState_np1 = fomStatesMngr_(::pressio::ode::nPlusOne());
//     auto & f_np1 = lspgStencilVelocities(::pressio::ode::nPlusOne());
//     fomSystem_.get().velocity(fomState_np1, t_np1, f_np1);

//     ::pressio::ode::impl::discrete_time_residual
// 	(fomState_np1, lspgResidual, fomStatesMngr_.get(),
// 	 lspgStencilVelocities, dt, ::pressio::ode::CrankNicolson());
//   }
