/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_residual_policy_continuous_time_api.hpp
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

#ifndef ROM_LSPG_UNSTEADY_HYPRED_POLICY_RESIDUAL_HPP_
#define ROM_LSPG_UNSTEADY_HYPRED_POLICY_RESIDUAL_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <
  class ResidualType,
  class FomStatesManagerType,
  class FomSystemType,
  class HypRedOperatorUpdater
  >
class UnsteadyHypRedResidualPolicy
{
public:
  // required
  using residual_type = ResidualType;

public:
  UnsteadyHypRedResidualPolicy() = delete;
  UnsteadyHypRedResidualPolicy(const UnsteadyHypRedResidualPolicy &) = default;
  UnsteadyHypRedResidualPolicy & operator=(const UnsteadyHypRedResidualPolicy &) = delete;
  UnsteadyHypRedResidualPolicy(UnsteadyHypRedResidualPolicy &&) = default;
  UnsteadyHypRedResidualPolicy & operator=(UnsteadyHypRedResidualPolicy &&) = delete;
  ~UnsteadyHypRedResidualPolicy() = default;

  UnsteadyHypRedResidualPolicy(const FomSystemType & fomSystem,
			       FomStatesManagerType & fomStatesMngr,
			       const HypRedOperatorUpdater & hrUpdater)
    : fomStatesMngr_(fomStatesMngr),
      fomSystem_(fomSystem),
      hypredOperatorUpdater_(hrUpdater),
      fomStateHelperInstance_(::pressio::ops::clone(fomStatesMngr(::pressio::ode::nPlusOne())))
  {}

public:
  residual_type create() const{
    residual_type R(fomSystem_.get().createVelocity());
    ::pressio::ops::set_zero(R);
    return R;
  }

  template <
    class LspgStateType,
    class LspgStencilStatesContainerType,
    class LspgStencilVelocitiesContainerType,
    class ScalarType>
  void operator()(::pressio::ode::StepScheme name,
		  const LspgStateType & lspgState,
		  const LspgStencilStatesContainerType & lspgStencilStates,
		  LspgStencilVelocitiesContainerType & lspgStencilVelocities,
		  const ScalarType & t_np1,
		  const ScalarType & dt,
		  const ::pressio::ode::step_count_type & currentStepNumber,
		  residual_type & lspgResidual) const
  {

    /* the currrent FOM has to be recomputed every time regardless of
     * whether the currentStepNumber changes since we might be inside a non-linear solve
     * where the time step does not change but this residual method
     * is called multiple times.
     */
    fomStatesMngr_.get().reconstructAt(lspgState, ::pressio::ode::nPlusOne());

    if (name == ::pressio::ode::StepScheme::BDF1)
    {
      (*this).template compute_impl_bdf<ode::BDF1>
	(name, lspgState, lspgStencilStates, lspgStencilVelocities,
	 t_np1, dt, currentStepNumber, lspgResidual);
    }
    else if (name == ::pressio::ode::StepScheme::BDF2)
    {
      (*this).template compute_impl_bdf<ode::BDF2>
	(name, lspgState, lspgStencilStates, lspgStencilVelocities,
	 t_np1, dt, currentStepNumber, lspgResidual);
    }
    else if (name == ::pressio::ode::StepScheme::CrankNicolson){
      (*this).compute_impl_cn(lspgState, lspgStencilStates,lspgStencilVelocities,
			      t_np1, dt, currentStepNumber, lspgResidual);
    }
  }

private:
  template <
    class StepperTag,
    class LspgStateType,
    class LspgStencilStatesContainerType,
    class LspgStencilVelocitiesContainerType,
    class ScalarType
    >
  void compute_impl_bdf(::pressio::ode::StepScheme name,
			const LspgStateType & lspgState,
			const LspgStencilStatesContainerType & lspgStencilStates,
			LspgStencilVelocitiesContainerType & lspgStencilVelocities,
			const ScalarType & t_np1,
			const ScalarType & dt,
			const ::pressio::ode::step_count_type & currentStepNumber,
			residual_type & lspgResidual) const
  {

    /* previous FOM states should only be recomputed when the time step changes.
     * The method below does not recompute all previous states, but only
     * recomputes the n-th state and updates/shifts back all the other
     * FOM states stored. */
    if (storedStep_ != currentStepNumber){
      const auto & lspgStateAt_n = lspgStencilStates(::pressio::ode::n());
      fomStatesMngr_.get().reconstructAtAndUpdatePrevious(lspgStateAt_n,
							  ::pressio::ode::n());
      storedStep_ = currentStepNumber;
    }

    const auto & fomStateAt_np1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    fomSystem_.get().velocity(fomStateAt_np1, t_np1, lspgResidual);

    if (name == ::pressio::ode::StepScheme::BDF1)
    {
      /*R(y_n+1) = cnp1*y_n+1 + cn*y_n + cf*f(t_n+1, y_n+1)

	so we do:
	1. lspgResidual = f(t_n+1, y_n+1) (see call to fom velocity above)

	2. fomStateHelpInstance_ = cnp1*y_np1 + cn*y_n

	3. we call the combiner to handle the rest
      */
      const auto & fomStateAt_n = fomStatesMngr_(::pressio::ode::n());

      constexpr auto zero = ::pressio::utils::Constants<ScalarType>::zero();
      constexpr auto cnp1 = ::pressio::ode::constants::bdf1<ScalarType>::c_np1_;
      constexpr auto cn   = ::pressio::ode::constants::bdf1<ScalarType>::c_n_;
      ::pressio::ops::update(fomStateHelperInstance_, zero,
			     fomStateAt_np1, cnp1,
			     fomStateAt_n, cn);

      constexpr auto one = ::pressio::utils::Constants<ScalarType>::one();
      const auto cf = ::pressio::ode::constants::bdf1<ScalarType>::c_f_ * dt;
      hypredOperatorUpdater_.get().updateSampleMeshOperandWithStencilMeshOne(lspgResidual, cf,
								       fomStateHelperInstance_, one);
    }

    else if (name == ::pressio::ode::StepScheme::BDF2)
    {
      /*R(y_n+1) = cnp1*y_n+1 + cn*y_n + cnm1*y_n-1 + cf*f(t_n+1, y_n+1)

	so we do:
	1. lspgResidual = f(t_n+1, y_n+1) (see call to fom velocity above)

	2. fomStateHelpInstance_ = cnp1*y_np1 + cn*y_n + cnm1*y_n-1

	3. we call the combiner to handle the rest
      */

      constexpr auto cnp1 = ::pressio::ode::constants::bdf2<ScalarType>::c_np1_;
      constexpr auto cn   = ::pressio::ode::constants::bdf2<ScalarType>::c_n_;
      constexpr auto cnm1 = ::pressio::ode::constants::bdf2<ScalarType>::c_nm1_;
      const auto cf	  = ::pressio::ode::constants::bdf2<ScalarType>::c_f_ * dt;

      const auto & fomStateAt_n   = fomStatesMngr_(::pressio::ode::n());
      const auto & fomStateAt_nm1 = fomStatesMngr_(::pressio::ode::nMinusOne());
      constexpr auto zero = ::pressio::utils::Constants<ScalarType>::zero();
      ::pressio::ops::update(fomStateHelperInstance_, zero,
			     fomStateAt_np1, cnp1,
			     fomStateAt_n, cn,
			     fomStateAt_nm1, cnm1);

      constexpr auto one = ::pressio::utils::Constants<ScalarType>::one();
      hypredOperatorUpdater_.get().updateSampleMeshOperandWithStencilMeshOne(lspgResidual, cf,
								       fomStateHelperInstance_, one);
    }
  }

  template <
    class LspgStateType,
    class LspgStencilStatesContainerType,
    class LspgStencilVelocitiesContainerType,
    class ScalarType
    >
  void compute_impl_cn(const LspgStateType & lspgState,
		       const LspgStencilStatesContainerType & lspgStencilStates,
		       // for CN, stencilVelocities holds f_n+1 and f_n
		       LspgStencilVelocitiesContainerType & lspgStencilVelocities,
		       const ScalarType & t_np1,
		       const ScalarType & dt,
		       const ::pressio::ode::step_count_type & currentStepNumber,
		       residual_type & lspgResidual) const
  {

    throw std::runtime_error("Lspg hyp-red residual policy for CN missing");

    // PRESSIOLOG_DEBUG("residual policy with compute_cn_impl");
    // if (storedStep_ != currentStepNumber){
    //   const auto & lspgStateAt_n = lspgStencilStates(::pressio::ode::n());
    //   fomStatesMngr_.get().reconstructAtAndUpdatePrevious(lspgStateAt_n,
    // 							  ::pressio::ode::n());
    //   storedStep_ = currentStepNumber;

    //   // if the step changed, I need to compute f(y_n, t_n)
    //   const auto tn = t_np1-dt;
    //   auto & f_n = lspgStencilVelocities(::pressio::ode::n());
    //   const auto & fomState_n = fomStatesMngr_(::pressio::ode::n());
    //   fomSystem_.get().velocity(fomState_n, tn, f_n);
    // }

    // // always compute f(y_n+1, t_n+1)
    // const auto & fomState_np1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    // auto & f_np1 = lspgStencilVelocities(::pressio::ode::nPlusOne());
    // fomSystem_.get().velocity(fomState_np1, t_np1, f_np1);

    // ::pressio::ode::impl::discrete_time_residual
    // 	(fomState_np1, lspgResidual, fomStatesMngr_.get(),
    // 	 lspgStencilVelocities, dt, ::pressio::ode::CrankNicolson());
  }

protected:
  // storedStep is used to keep track of which step we are at.
  // used to decide if we need to update/recompute the previous FOM states or not.
  // To avoid recomputing previous FOM states if we are not in a new time step.
  mutable int32_t storedStep_ = -1;

  std::reference_wrapper<FomStatesManagerType> fomStatesMngr_;
  std::reference_wrapper<const FomSystemType> fomSystem_;

  std::reference_wrapper<const HypRedOperatorUpdater> hypredOperatorUpdater_;
  mutable typename FomStatesManagerType::value_type fomStateHelperInstance_;
};

}}}}
#endif  // ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_CONTINUOUS_TIME_API_HPP_
