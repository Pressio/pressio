/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_jacobian_policy_discrete_time_api.hpp
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

#ifndef ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_JACOBIAN_POLICY_DISCRETE_TIME_API_HPP_
#define ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_JACOBIAN_POLICY_DISCRETE_TIME_API_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template<
  typename rom_jacobian_type,
  typename fom_apply_jacobian_ret_type,
  typename decoder_type,
  typename fom_states_manager_t
  >
class JacobianPolicyDiscreteTimeApi
{
public:
  using scalar_t =
    typename ::pressio::containers::details::traits<rom_jacobian_type>::scalar_t;
  using rom_jacobian_t = rom_jacobian_type;
  using time_step_type = ::pressio::ode::types::step_t;

public:
  JacobianPolicyDiscreteTimeApi() = delete;
  JacobianPolicyDiscreteTimeApi(const JacobianPolicyDiscreteTimeApi &) = default;
  JacobianPolicyDiscreteTimeApi & operator=(const JacobianPolicyDiscreteTimeApi &) = default;
  JacobianPolicyDiscreteTimeApi(JacobianPolicyDiscreteTimeApi &&) = default;
  JacobianPolicyDiscreteTimeApi & operator=(JacobianPolicyDiscreteTimeApi &&) = default;
  ~JacobianPolicyDiscreteTimeApi() = default;

  template< typename app_t>
  JacobianPolicyDiscreteTimeApi(fom_states_manager_t & fomStatesMngr,
				const decoder_type & decoder,
				const app_t & appObj)
    : fomStatesMngr_(fomStatesMngr),
      phi_(decoder.getReferenceToJacobian()),
      fomApplyJac_(appObj.createApplyDiscreteTimeJacobianResult(*phi_.get().data()))
  {}

public:
  template <typename fom_system_t>
  rom_jacobian_t create(const fom_system_t & fomSystemObj) const
  {
    const auto nRows = phi_.get().extent(1);
    const auto nCols = phi_.get().extent(1);
    rom_jacobian_t romJac(nRows, nCols);
    return romJac;
  }

  template <
    typename ode_tag,
    typename rom_state_t,
    typename rom_prev_states_t,
    typename fom_system_t
    >
  void compute(const rom_state_t	& romState,
	       const rom_prev_states_t	& romPrevStates,
	       const fom_system_t	& fomSystemObj,
	       const scalar_t		& time,
	       const scalar_t		& dt,
	       const time_step_type	& step,
	       rom_jacobian_t		& romJac) const
  {
    this->compute_impl(romState, romPrevStates, fomSystemObj,
		       time, dt, step, romJac);
  }

private:
  // we have here n = 1 prev rom states
  template<
    typename rom_state_t, typename rom_prev_states_t, typename fom_system_t,
    typename scalar_t, typename rom_jac_t
  >
  mpl::enable_if_t< rom_prev_states_t::size()==1 >
  compute_impl(const rom_state_t	& romState,
	       const rom_prev_states_t	& romPrevStates,
	       const fom_system_t	& fomSystemObj,
	       const scalar_t		& time,
	       const scalar_t		& dt,
	       const time_step_type	& step,
	       rom_jac_t		& romJac) const
  {
    // here we assume that the residual policy already:
    // - reconstucted the previous states
    // - called the decoder to update the jacobian

    const auto & yn   = fomStatesMngr_.get().getCRefToCurrentFomState();
    const auto & ynm1 = fomStatesMngr_.get().getCRefToFomStatePrevStep();
    ::pressio::rom::queryFomApplyDiscreteTimeJacobian(yn, ynm1, fomSystemObj,
						      time, dt, step, phi_.get(),
						      fomApplyJac_);

    constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
                            one, phi_.get(), fomApplyJac_, zero, romJac);
  }

  // we have here n = 2 prev rom states
  template<
    typename rom_state_t, typename rom_prev_states_t, typename fom_system_t,
    typename scalar_t, typename rom_jac_t
    >
  mpl::enable_if_t< rom_prev_states_t::size()==2 >
  compute_impl(const rom_state_t	& romState,
	       const rom_prev_states_t	& romPrevStates,
	       const fom_system_t       & fomSystemObj,
	       const scalar_t		& time,
	       const scalar_t		& dt,
	       const time_step_type	& step,
	       rom_jac_t		& romJac) const
  {
    // here we assume that the residual policy already:
    // - reconstucted the previous states
    // - called the decoder to update the jacobian

    const auto & yn   = fomStatesMngr_.get().getCRefToCurrentFomState();
    const auto & ynm1 = fomStatesMngr_.get().getCRefToFomStatePrevStep();
    const auto & ynm2 = fomStatesMngr_.get().getCRefToFomStatePrevStep();
    ::pressio::rom::queryFomApplyDiscreteTimeJacobian(yn, ynm1, ynm2,
						      fomSystemObj, time,
						      dt, step, phi_.get(),
						      fomApplyJac_);

    constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
                            one, phi_.get(), fomApplyJac_, zero, romJac);
  }

private:
  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
  std::reference_wrapper<const typename decoder_type::jacobian_type> phi_;
  mutable fom_apply_jacobian_ret_type fomApplyJac_;
};

}}}}//end namespace
#endif  // ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_JACOBIAN_POLICY_DISCRETE_TIME_API_HPP_
