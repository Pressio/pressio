/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_hyper_reduced_jacobian_policy_continuous_time_api.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_HYPER_REDUCED_JACOBIAN_POLICY_CONTINUOUS_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_HYPER_REDUCED_JACOBIAN_POLICY_CONTINUOUS_TIME_API_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template<
  typename fom_states_manager_t,
  typename apply_jac_return_type,
  typename decoder_type,
  typename sample_to_stencil_t
  >
class HypRedJacobianPolicyContinuousTimeApi
{

public:
  using data_type = apply_jac_return_type;

public:
  HypRedJacobianPolicyContinuousTimeApi() = delete;
  HypRedJacobianPolicyContinuousTimeApi(const HypRedJacobianPolicyContinuousTimeApi &) = default;
  HypRedJacobianPolicyContinuousTimeApi & operator=(const HypRedJacobianPolicyContinuousTimeApi &) = delete;
  HypRedJacobianPolicyContinuousTimeApi(HypRedJacobianPolicyContinuousTimeApi &&) = default;
  HypRedJacobianPolicyContinuousTimeApi & operator=(HypRedJacobianPolicyContinuousTimeApi &&) = delete;
  ~HypRedJacobianPolicyContinuousTimeApi() = default;

  HypRedJacobianPolicyContinuousTimeApi(fom_states_manager_t & fomStatesMngr,
					decoder_type & decoder,
					const sample_to_stencil_t & sTosInfo)
    : fomStatesMngr_(fomStatesMngr),
      decoderObj_(decoder),
      decoderJacobian_(decoder.jacobianCRef()),
      sTosInfo_(sTosInfo)
  {}

public:
  template <typename fom_system_t>
  apply_jac_return_type create(const fom_system_t & fomObj) const
  {
    const auto & decJac = decoderJacobian_.get();
    return apply_jac_return_type(fomObj.createApplyJacobianResult(*decJac.data()));
  }

  template <
    typename stepper_tag,
    typename stencil_states_t,
    typename lspg_state_t,
    typename lspg_jac_t,
    typename fom_system_t,
    typename scalar_t
    >
  void compute(const lspg_state_t & romState,
	       const stencil_states_t & stencilStates,
	       const fom_system_t & fomSystemObj,
	       const scalar_t & timeAtNextStep,
	       const scalar_t & dt,
	       const ::pressio::ode::types::step_t & currentStepNumber,
	       lspg_jac_t & romJac) const
  {
    // since this is for hyp-red, I need to make sure the sTosInfo
    // is consistent with romJacobian
    assert(sTosInfo_.get().extent(0) == romJac.extent(0));

    this->compute_impl<stepper_tag>(romState, romJac, fomSystemObj,
				    timeAtNextStep, dt, currentStepNumber);
  }

private:
  template <
    typename stepper_tag,
    typename lspg_state_t,
    typename lspg_jac_t,
    typename fom_system_t,
    typename scalar_t
    >
  void compute_impl(const lspg_state_t & romState,
		    lspg_jac_t & romJac,
		    const fom_system_t & fomSystemObj,
		    const scalar_t   & timeAtNextStep,
		    const scalar_t   & dt,
		    const ::pressio::ode::types::step_t & currentStepNumber) const
  {
    // here we assume that the current state has already been reconstructd
    // by the residual policy. So we do not recompute the FOM state.
    // Maybe we should find a way to ensure this is the case.
    // fomStatesMngr_.get().reconstructCurrentFomState(romState);

    // update Jacobian of decoder
    decoderObj_.get().updateJacobian(romState);

    const auto & fomState = fomStatesMngr_(::pressio::ode::nPlusOne());
    const auto & basis = decoderObj_.get().jacobianCRef();
    fomSystemObj.applyJacobian(*fomState.data(), *basis.data(),
			       timeAtNextStep, *romJac.data());

    ::pressio::rom::lspg::impl::unsteady::time_discrete_jacobian<
      stepper_tag>(romJac, dt, decoderJacobian_.get(), sTosInfo_.get());
  }

protected:
  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
  std::reference_wrapper<decoder_type> decoderObj_ = {};
  std::reference_wrapper<const typename decoder_type::jacobian_type> decoderJacobian_ = {};
  std::reference_wrapper<const sample_to_stencil_t> sTosInfo_;
};

}}}}}
#endif  // ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_HYPER_REDUCED_JACOBIAN_POLICY_CONTINUOUS_TIME_API_HPP_
