/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_jacobian_policy.hpp
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

#ifndef ROM_LSPG_STEADY_JACOBIAN_POLICY_HPP_
#define ROM_LSPG_STEADY_JACOBIAN_POLICY_HPP_

#include "../../rom_fwd.hpp"
#include "../../rom_data_fom_states.hpp"

namespace pressio{ namespace rom{

template<
  typename fom_states_data,
  typename apply_jac_return_type,
  typename fom_apply_jac_policy,
  typename decoder_type
  >
class LSPGSteadyJacobianPolicy
  : protected fom_apply_jac_policy{

protected:
  using this_t = LSPGSteadyJacobianPolicy<
  fom_states_data, apply_jac_return_type,
  fom_apply_jac_policy, decoder_type>;

public:
  static constexpr bool isResidualPolicy_ = false;
  using apply_jac_return_t = apply_jac_return_type;

public:
  LSPGSteadyJacobianPolicy() = delete;
  ~LSPGSteadyJacobianPolicy() = default;

  LSPGSteadyJacobianPolicy(fom_states_data	& fomStates,
			   const fom_apply_jac_policy	& applyJacFunctor,
			   const apply_jac_return_type	& applyJacObj,
			   const decoder_type		& decoder)
    : fomStates_(fomStates),
      fom_apply_jac_policy(applyJacFunctor),
      JJ_(applyJacObj),
      decoderObj_(decoder){}

public:

  template <
    typename lspg_state_t,
    typename lspg_jac_t,
    typename app_t
  >
  void operator()(const lspg_state_t & romY,
		  lspg_jac_t	     & romJJ,
  		  const app_t	     & app) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg apply jac");
#endif

    // todo: this is not needed if jacobian is called after resiudal
    // because residual takes care of reconstructing the fom state
    //    timer->start("reconstruct fom state");
    fomStates_.template reconstructCurrentFomState(romY);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom apply jac");
#endif

    const auto & basis = decoderObj_.getReferenceToJacobian();
    fom_apply_jac_policy::evaluate(app, fomStates_.getCRefToFomState(), basis, romJJ);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom apply jac");
    timer->stop("lspg apply jac");
#endif
  }

  template <typename lspg_state_t, typename app_t>
  apply_jac_return_t operator()(const lspg_state_t & romY,
				const app_t	   & app) const
  {
    (*this).template operator()(romY, JJ_, app);
    return JJ_;
  }

protected:
  mutable apply_jac_return_t JJ_   = {};
  const decoder_type & decoderObj_ = {};
  fom_states_data & fomStates_;
};

}}//end namespace pressio::rom
#endif
