/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_jacobian_policy_discrete_time_api.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_DISCRETE_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_DISCRETE_TIME_API_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template<
  typename fom_states_manager_t,
  typename apply_jac_return_type,
  typename decoder_type
  >
class JacobianPolicyDiscreteTimeApi
{

public:
  using apply_jac_return_t = apply_jac_return_type;

public:
  JacobianPolicyDiscreteTimeApi() = delete;
  JacobianPolicyDiscreteTimeApi(const JacobianPolicyDiscreteTimeApi &) = default;
  JacobianPolicyDiscreteTimeApi & operator=(const JacobianPolicyDiscreteTimeApi &) = default;
  JacobianPolicyDiscreteTimeApi(JacobianPolicyDiscreteTimeApi &&) = default;
  JacobianPolicyDiscreteTimeApi & operator=(JacobianPolicyDiscreteTimeApi &&) = default;
  ~JacobianPolicyDiscreteTimeApi() = default;

  JacobianPolicyDiscreteTimeApi(fom_states_manager_t & fomStatesMngr,
				const decoder_type & decoder)
    : decoderObj_(decoder), fomStatesMngr_(fomStatesMngr){}

public:
  template <typename fom_system_t>
  apply_jac_return_t create(const fom_system_t & fomSystemObj) const
  {
    // // this is only called once
    const auto & phi = decoderObj_.get().jacobianCRef();
    apply_jac_return_t romJac
      (fomSystemObj.createApplyDiscreteTimeJacobianResult(*phi.data()));
    return romJac;
  }

  template <
    typename ode_tag,
    typename lspg_state_t,
    typename lspg_prev_states_t,
    typename lspg_jac_t,
    typename fom_system_t,
    typename scalar_t
    >
  void compute(const lspg_state_t & romState,
	       const lspg_prev_states_t	& romPrevStates,
	       const fom_system_t & fomSystemObj,
	       const scalar_t & time,
	       const scalar_t & dt,
	       const ::pressio::ode::types::step_t & step,
	       lspg_jac_t & romJac) const
  {
    this->compute_impl(romState, romPrevStates,
		       fomSystemObj, time, dt, step, romJac);
  }

private:
  // we have here n = 1 prev rom states
  template<
    typename lspg_state_t, typename lspg_prev_states_t, typename fom_system_t,
    typename scalar_t, typename lspg_jac_t
    >
  mpl::enable_if_t< lspg_prev_states_t::size()==1 >
  compute_impl(const lspg_state_t			& romState,
		    const lspg_prev_states_t		& romPrevStates,
  		    const fom_system_t		        & fomSystemObj,
  		    const scalar_t			& time,
  		    const scalar_t			& dt,
		    const ::pressio::ode::types::step_t	& step,
		    lspg_jac_t				& romJac) const
  {
    // here we assume that the current state has already been reconstructd
    // by the residual policy. So we do not recompute the FOM state.
    // Maybe we should find a way to ensure this is the case.
    //fomStatesMngr_.get().template reconstructCurrentFomState(romState);

    // update Jacobian of decoder
    decoderObj_.get().updateJacobian(romState);

    const auto & phi = decoderObj_.get().jacobianCRef();
    const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
    const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();
    fomSystemObj.applyDiscreteTimeJacobian(step, time, dt,
					   *phi.data(),
					   *romJac.data(),
					   *yn.data(),
					   *ynm1.data());
  }

  // we have here n = 2 prev rom states
  template<
    typename lspg_state_t, typename lspg_prev_states_t, typename fom_system_t,
    typename scalar_t, typename lspg_jac_t
  >
  mpl::enable_if_t< lspg_prev_states_t::size()==2 >
  compute_impl(const lspg_state_t			& romState,
		    const lspg_prev_states_t		& romPrevStates,
  		    const fom_system_t		        & fomSystemObj,
  		    const scalar_t			& time,
  		    const scalar_t			& dt,
		    const ::pressio::ode::types::step_t	& step,
		    lspg_jac_t				& romJac) const
  {
    // here we assume that the current state has already been reconstructd
    // by the residual policy. So we do not recompute the FOM state.
    // Maybe we should find a way to ensure this is the case.
    //fomStatesMngr_.get().template reconstructCurrentFomState(romState);

    // update Jacobian of decoder
    decoderObj_.get().updateJacobian(romState);

    const auto & phi = decoderObj_.get().jacobianCRef();
    const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
    const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();
    const auto & ynm2 = fomStatesMngr_.get().fomStatePrevStepCRef();

    fomSystemObj.applyDiscreteTimeJacobian(step, time, dt,
					   *phi.data(),
					   *romJac.data(),
					   *yn.data(),
					   *ynm1.data(),
					   *ynm2.data());
  }

protected:
  std::reference_wrapper<const decoder_type> decoderObj_;
  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
};

}}}}}//end namespace pressio::rom::lspg::unsteady::impl
#endif  // ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_DISCRETE_TIME_API_HPP_
