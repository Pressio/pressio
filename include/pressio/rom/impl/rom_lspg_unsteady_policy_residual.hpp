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

#ifndef ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_CONTINUOUS_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_CONTINUOUS_TIME_API_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <bool is_cont_time, class ResidualType, class FomStatesManagerType, class FomSystemType>
class UnsteadyResidualPolicy
{
public:
  // required
  using residual_type = ResidualType;

public:
  UnsteadyResidualPolicy() = delete;
  UnsteadyResidualPolicy(const UnsteadyResidualPolicy &) = default;
  UnsteadyResidualPolicy & operator=(const UnsteadyResidualPolicy &) = delete;
  UnsteadyResidualPolicy(UnsteadyResidualPolicy &&) = default;
  UnsteadyResidualPolicy & operator=(UnsteadyResidualPolicy &&) = delete;
  ~UnsteadyResidualPolicy() = default;

  UnsteadyResidualPolicy(const FomSystemType & fomSystem,
			 FomStatesManagerType & fomStatesMngr)
    : fomStatesMngr_(fomStatesMngr),
      fomSystem_(fomSystem)
  {}

public:
  template<bool _is_cont_time = is_cont_time>
  mpl::enable_if_t<_is_cont_time, residual_type>
  create() const{
    residual_type R(fomSystem_.get().createVelocity());
    ::pressio::ops::set_zero(R);
    return R;
  }

  template <
    class StepperTag,
    class LspgStateType,
    class LspgStencilStatesContainerType,
    class LspgStencilVelocitiesContainerType,
    class ScalarType
    >
  mpl::enable_if_t<
     std::is_same<StepperTag, ::pressio::ode::BDF1>::value or
     std::is_same<StepperTag, ::pressio::ode::BDF2>::value
    >
  compute(const LspgStateType & lspgState,
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

    /* previous FOM states should only be recomputed when the time step changes.
     * The method below does not recompute all previous states, but only
     * recomputes the n-th state and updates/shifts back all the other
     * FOM states stored. */
    if (storedStep_ != currentStepNumber){
      fomStatesMngr_.get().reconstructWithStencilUpdate(lspgStencilStates(::pressio::ode::n()));
      storedStep_ = currentStepNumber;
    }

    const auto & fomState = fomStatesMngr_(::pressio::ode::nPlusOne());
    fomSystem_.get().velocity(fomState, t_np1, lspgResidual);
    ::pressio::ode::impl::discrete_time_residual(fomState, lspgResidual,
						 fomStatesMngr_.get(), dt,
						 StepperTag());
  }

  template <
    class StepperTag,
    class LspgStateType,
    class LspgStencilStatesContainerType,
    class LspgStencilVelocitiesContainerType,
    class ScalarType
    >
  mpl::enable_if_t<
    std::is_same<StepperTag, ::pressio::ode::CrankNicolson>::value
    >
  compute(const LspgStateType & lspgState,
	  const LspgStencilStatesContainerType & lspgStencilStates,
	  LspgStencilVelocitiesContainerType & lspgStencilVelocities,
	  const ScalarType & t_np1,
	  const ScalarType & dt,
	  const ::pressio::ode::step_count_type & currentStepNumber,
	  residual_type & lspgResidual) const
  {
    // this->compute_cn_impl<StepperTag>(lspgState, lspgResidual,
    // 				       lspgStencilStates, lspgStencilVelocities,
    // 				       t_np1, dt, currentStepNumber);
  }

private:
  // template <
  //   class StepperTag,
  //   class LspgStateType,
  //   class LspgStencilStatesContainerType,
  //   class LspgStencilVelocitiesContainerType,
  //   class ScalarType
  //   >
  // void compute_impl(const LspgStateType & lspgState,
  // 	       residual_type & lspgResidual,
  // 	       const LspgStencilStatesContainerType & lspgStencilStates,
  // 	       LspgStencilVelocitiesContainerType & lspgStencilVelocities,
  // 	       const ScalarType & t_np1,
  // 	       const ScalarType & dt,
  // 	       const ::pressio::ode::step_count_type & currentStepNumber) const
  // {

  // }

//   template <
//     class StepperTag,
//     class lspg_state_t,
//     class stencil_states_t,
//     class fom_system_t,
//     class scalar_t,
//     class stencil_velocities_t
//   >
//   void compute_cn_impl(const lspg_state_t & lspgState,
// 		       residual_type & lspgResidual,
// 		       const stencil_states_t & lspgStencilStates,
// 		       const fom_system_t & fomSystemObj,
// 		       const scalar_t & t_np1,
// 		       const scalar_t & dt,
// 		       const ::pressio::ode::step_count_type & currentStepNumber,
// 		       // for CN, stencilVelocities holds f_n+1 and f_n
// 		       stencil_velocities_t & stencilVelocities) const
//   {
//     PRESSIOLOG_DEBUG("residual policy with compute_cn_impl");

//     fomStatesMngr_.get().reconstructAt(lspgState, ::pressio::ode::nPlusOne());

//     if (storedStep_ != currentStepNumber){
//       fomStatesMngr_.get().reconstructWithStencilUpdate(lspgStencilStates(ode::n()));
//       storedStep_ = currentStepNumber;

//       // if the step changed, I need to compute f(y_n, t_n)
//       const auto tn = t_np1-dt;
//       auto & f_n = stencilVelocities(::pressio::ode::n());
//       const auto & fomState_n = fomStatesMngr_(::pressio::ode::n());
//       fomSystemObj.velocity(*fomState_n.data(), tn, *f_n.data());
//     }

//     // always compute f(y_n+1, t_n+1)
//     auto & f_np1 = stencilVelocities(::pressio::ode::nPlusOne());
//     const auto & fomState_np1 = fomStatesMngr_(::pressio::ode::nPlusOne());
//     fomSystemObj.velocity(*fomState_np1.data(), t_np1, *f_np1.data());

//     ::pressio::rom::lspg::impl::unsteady::time_discrete_residual
// 	<StepperTag>(fomStatesMngr_.get(), stencilVelocities, lspgResidual, dt);
//   }

protected:
  // storedStep is used to keep track of which step we are at.
  // used to decide if we need to update/recompute the previous FOM states or not.
  // To avoid recomputing previous FOM states if we are not in a new time step.
  mutable int32_t storedStep_ = {};

  std::reference_wrapper<FomStatesManagerType> fomStatesMngr_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
};

}}}}
#endif  // ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_CONTINUOUS_TIME_API_HPP_
