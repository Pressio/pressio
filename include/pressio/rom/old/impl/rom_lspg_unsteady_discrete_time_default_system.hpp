/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_discrete_time_default_system.hpp
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

#ifndef ROM_IMPL_ROM_LSPG_UNSTEADY_DISCRETE_TIME_DEFAULT_SYSTEM_HPP_
#define ROM_IMPL_ROM_LSPG_UNSTEADY_DISCRETE_TIME_DEFAULT_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <
  class ScalarType,
  class LspgStateType,
  class LspgResidualType,
  class LspgJacobianType,
  class FomSystemType,
  class DecoderType,
  class FomStatesManagerType
  >
struct DiscreteTimeDefaultSystem
{
  using scalar_type = ScalarType;
  using state_type  = LspgStateType;
  using discrete_time_residual_type = LspgResidualType;
  using discrete_time_jacobian_type = LspgJacobianType;

  DiscreteTimeDefaultSystem() = delete;
  DiscreteTimeDefaultSystem(const DiscreteTimeDefaultSystem &) = default;
  DiscreteTimeDefaultSystem & operator=(const DiscreteTimeDefaultSystem &) = delete;
  DiscreteTimeDefaultSystem(DiscreteTimeDefaultSystem &&) = default;
  DiscreteTimeDefaultSystem & operator=(DiscreteTimeDefaultSystem &&) = delete;
  ~DiscreteTimeDefaultSystem() = default;

  DiscreteTimeDefaultSystem(const FomSystemType & fomSystem,
			    DecoderType & decoder,
			    FomStatesManagerType & fomStatesMngr)
    : fomSystem_(fomSystem),
      fomStatesMngr_(fomStatesMngr),
      decoderObj_(decoder),
      phi_(decoder.jacobianCRef())
  {}

  discrete_time_residual_type createDiscreteTimeResidual() const{
    discrete_time_residual_type R(fomSystem_.get().createDiscreteTimeResidual());
    return R;
  }

  discrete_time_jacobian_type createDiscreteTimeJacobian() const{
    const auto & phi = decoderObj_.get().jacobianCRef();
    discrete_time_jacobian_type J(fomSystem_.get().createApplyDiscreteTimeJacobianResult(phi));
    return J;
  }

  template<class StepCountType>
  void discreteTimeResidual(StepCountType currentStepNumber,
                            scalar_type time_np1,
                            scalar_type dt,
                            discrete_time_residual_type & R,
                            const state_type & lspg_state_np1,
			    const state_type & lspg_state_n) const
  {

    doFomStatesReconstruction(currentStepNumber, lspg_state_np1, lspg_state_n);

    const auto & ynp1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesMngr_(::pressio::ode::n());
    try{
      fomSystem_.get().discreteTimeResidual(currentStepNumber, time_np1, dt,
					    R, ynp1, yn);
    }
    catch (::pressio::eh::DiscreteTimeResidualFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

  template<class StepCountType>
  void discreteTimeResidual(StepCountType currentStepNumber,
			    scalar_type time_np1,
			    scalar_type dt,
			    discrete_time_residual_type & R,
			    const state_type & lspg_state_np1,
			    const state_type & lspg_state_n,
			    const state_type & lspg_state_nm1) const
  {
    doFomStatesReconstruction(currentStepNumber,
			      lspg_state_np1, lspg_state_n, lspg_state_nm1);

    const auto & ynp1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesMngr_(::pressio::ode::n());
    const auto & ynm1 = fomStatesMngr_(::pressio::ode::nMinusOne());
    try{
      fomSystem_.get().discreteTimeResidual(currentStepNumber, time_np1, dt,
					    R, ynp1, yn, ynm1);
    }
    catch (::pressio::eh::DiscreteTimeResidualFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

  template<class StepCountType>
  void discreteTimeJacobian(StepCountType currentStepNumber,
			    scalar_type time_np1,
			    scalar_type dt,
			    discrete_time_jacobian_type & J,
			    const state_type & lspg_state_np1,
			    const state_type & lspg_state_n) const
  {
    // here we assume that the current state has already been reconstructd
    // by the residual policy. So we do not recompute the FOM state.
    // Maybe we should find a way to ensure this is the case.
    //fomStatesMngr_.get().template reconstructCurrentFomState(romState);

    // update Jacobian of decoder
    decoderObj_.get().updateJacobian(lspg_state_np1);

    const auto & ynp1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesMngr_(::pressio::ode::n());

    fomSystem_.get().applyDiscreteTimeJacobian(currentStepNumber, time_np1,
					       dt, phi_.get(), J, ynp1, yn);
  }

  template<class StepCountType>
  void discreteTimeJacobian(StepCountType currentStepNumber,
			    scalar_type time_np1,
			    scalar_type dt,
			    discrete_time_jacobian_type & J,
			    const state_type & lspg_state_np1,
			    const state_type & lspg_state_n,
			    const state_type & lspg_state_nm1) const
  {
    // here we assume that the current state has already been reconstructd
    // by the residual policy. So we do not recompute the FOM state.
    // Maybe we should find a way to ensure this is the case.
    //fomStatesMngr_.get().template reconstructCurrentFomState(romState);

    // update Jacobian of decoder
    decoderObj_.get().updateJacobian(lspg_state_np1);

    const auto & ynp1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesMngr_(::pressio::ode::n());
    const auto & ynm1 = fomStatesMngr_(::pressio::ode::nMinusOne());

    fomSystem_.get().applyDiscreteTimeJacobian(currentStepNumber, time_np1,
					       dt, phi_.get(), J, ynp1, yn, ynm1);
  }

private:
  void doFomStatesReconstruction(const int32_t & step_number,
				 const state_type & lspg_state_np1) const
  {
    fomStatesMngr_.get().reconstructAt(lspg_state_np1, ::pressio::ode::nPlusOne());
  }

  void doFomStatesReconstruction(const int32_t & step_number,
				 const state_type & lspg_state_np1,
				 const state_type & lspg_state_n) const
  {
    /* the FOM state corresponding to the new predicted state has to be
     * recomputed every time regardless of the time step chaning or not,
     *  since we might be inside a non-linear solve
     * where the time step does not change but this residual method
     * is called multiple times. */
    fomStatesMngr_.get().reconstructAt(lspg_state_np1, ::pressio::ode::nPlusOne());

    /* previous FOM states should only be recomputed when the time step changes.
     * The method below does not recompute all previous states, but only
     * recomputes the n-th state and updates/shifts back all the other
     * FOM states stored. */
    if (storedStep_ != step_number){
      fomStatesMngr_.get().reconstructAtAndUpdatePrevious(lspg_state_n,
							  ::pressio::ode::n());
      storedStep_ = step_number;
    }
  }

  void doFomStatesReconstruction(const int32_t & step_number,
				 const state_type & lspg_state_np1,
				 const state_type & lspg_state_n,
				 const state_type & lspg_state_nm1) const
  {
    (void)lspg_state_nm1;
    doFomStatesReconstruction(step_number, lspg_state_np1, lspg_state_n);
  }

protected:
  // storedStep is used to keep track of which step we are doing.
  // This is used to decide whether we need to update/recompute the previous
  // FOM states or not. Since it does not make sense to recompute previous
  // FOM states if we are not in a new time step.
  mutable int32_t storedStep_ = {};

  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<FomStatesManagerType> fomStatesMngr_;
  std::reference_wrapper<DecoderType> decoderObj_;
  std::reference_wrapper<const typename DecoderType::jacobian_type> phi_;
};

}}}}//end namespace pressio::rom::lspg::impl
#endif  // ROM_IMPL_ROM_LSPG_UNSTEADY_DISCRETE_TIME_DEFAULT_SYSTEM_HPP_
