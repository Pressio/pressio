/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_residual_policy_residual_api.hpp
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

#ifndef ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_RESIDUAL_api_HPP_
#define ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_RESIDUAL_api_HPP_

#include "../../rom_fwd.hpp"
#include "../../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "../../rom_data_fom_states.hpp"

namespace pressio{ namespace rom{ namespace impl{

template <
  typename residual_type,
  typename fom_states_data_type,
  typename fom_querier_policy
  >
class LSPGUnsteadyResidualPolicyResidualApi
  : public ode::policy::ImplicitResidualPolicyBase<
      LSPGUnsteadyResidualPolicyResidualApi<residual_type,
					    fom_states_data_type,
					    fom_querier_policy>>,
    protected fom_querier_policy
{

public:
  using this_t = LSPGUnsteadyResidualPolicyResidualApi<residual_type,
						       fom_states_data_type,
						       fom_querier_policy>;
  friend ode::policy::ImplicitResidualPolicyBase<this_t>;

  static constexpr bool isResidualPolicy_ = true;
  using residual_t = residual_type;

public:
  LSPGUnsteadyResidualPolicyResidualApi() = delete;
  ~LSPGUnsteadyResidualPolicyResidualApi() = default;

  LSPGUnsteadyResidualPolicyResidualApi(const residual_t & RIn,
					fom_states_data_type & fomStatesIn,
					const fom_querier_policy & fomQuerierFunctor)
    : R_{RIn},
      fomStates_(fomStatesIn),
      fom_querier_policy(fomQuerierFunctor){}

public:
  template <
    int n,
    typename lspg_state_t,
    typename fom_t,
    typename scalar_t
  >
  void operator()(const lspg_state_t		   & romY,
		  residual_t			   & romR,
  		  const std::array<lspg_state_t,n> & romOldYs,
  		  const fom_t			   & app,
		  scalar_t			   t,
		  scalar_t			   dt,
		  ::pressio::ode::types::step_t step) const
  {
    std::cout << " empty LSPGUnsteadyResidualPolicyResidualApi" << std::endl;
    //this->compute_impl<n>(romY, romR, romOldYs, app, t, dt, step);
  }

  template <
    int n,
    typename lspg_state_t,
    typename fom_t,
    typename scalar_t
    >
  residual_t operator()(const lspg_state_t		& romY,
			const std::array<lspg_state_t,n>  & romOldYs,
			const fom_t			& app,
			scalar_t			t,
			scalar_t			dt,
			::pressio::ode::types::step_t   step) const
  {
    std::cout << " empty LSPGUnsteadyResidualPolicyResidualApi" << std::endl;
    //this->compute_impl<n>(romY, R_, romOldYs, app, t, dt, step);
    return R_;
  }

// private:
//   template <
//     int n,
//     typename lspg_state_t,
//     typename fom_t,
//     typename scalar_t
//   >
//   void compute_impl(const lspg_state_t		     & romY,
// 		    residual_t			     & romR,
// 		    const std::array<lspg_state_t,n> & romOldYs,
// 		    const fom_t			     & app,
// 		    scalar_t			     t,
// 		    scalar_t			     dt,
// 		    ::pressio::ode::types::step_t   step) const
//   {
//     fomStates_.template reconstructCurrentFomState(romY);
//     fomStates_.template reconstructFomOldStates<n>(romOldYs);

//     // fom_querier_policy::evaluate<n>(app,
//     // 				    fomStates_.getCRefToFomState(),
//     // 				    fomStates_.getCRefToFomOldStates(), romR, t);
//   }

protected:
  mutable residual_t R_ = {};
  fom_states_data_type & fomStates_;

};//end class

}}}//end namespace pressio::rom::impl
#endif
