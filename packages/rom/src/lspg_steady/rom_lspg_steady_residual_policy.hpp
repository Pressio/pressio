/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_residual_policy.hpp
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

#ifndef ROM_LSPG_STEADY_RESIDUAL_POLICY_HPP_
#define ROM_LSPG_STEADY_RESIDUAL_POLICY_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace steady{

template <
  typename residual_type,
  typename fom_states_data,
  typename fom_eval_rhs_policy
  >
class ResidualPolicy : protected fom_eval_rhs_policy
{

public:
  using this_t = ResidualPolicy<residual_type, fom_states_data,	fom_eval_rhs_policy>;

  static constexpr bool isResidualPolicy_ = true;
  using residual_t = residual_type;

public:
  ResidualPolicy() = delete;
  ~ResidualPolicy() = default;

  ResidualPolicy(const residual_t & RIn,
		 fom_states_data & fomStatesIn,
		 const fom_eval_rhs_policy & fomEvalRhsFunctor)
    : fom_eval_rhs_policy(fomEvalRhsFunctor),
      R_{RIn},
      fomStates_(fomStatesIn){}

public:
  template <
    typename lspg_state_t,
    typename lspg_residual_t,
    typename fom_t>
  void operator()(const lspg_state_t	& romState,
		  lspg_residual_t	& romR,
  		  const fom_t		& app) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg residual");
#endif

    fomStates_.template reconstructCurrentFomState(romState);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif

    fom_eval_rhs_policy::evaluate(app, fomStates_.getCRefToCurrentFomState(), romR);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
    timer->stop("lspg residual");
#endif
  }

  template <typename lspg_state_t, typename fom_t>
  residual_t operator()(const lspg_state_t & romState,
		       const fom_t	  & app) const
  {
    (*this).template operator()(romState, R_, app);
    return R_;
  }

protected:
  mutable residual_t R_ = {};
  fom_states_data & fomStates_;

};//end class

}}}}//end namespace pressio::rom::lspg::steady
#endif
