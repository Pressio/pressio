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

#include "../../rom_fwd.hpp"
#include "../../rom_data_fom_rhs.hpp"
#include "../../rom_data_fom_states.hpp"

namespace pressio{ namespace rom{

template <
  typename fom_states_data,
  typename fom_rhs_data,
  typename fom_eval_rhs_policy
  >
class LSPGSteadyResidualPolicy
  : protected fom_states_data,
    protected fom_rhs_data,
    protected fom_eval_rhs_policy{

protected:
  // protected because we might have decorators of this class
  using this_t = LSPGSteadyResidualPolicy<fom_states_data,
					  fom_rhs_data,
					  fom_eval_rhs_policy>;

  using fom_states_data::yFom_;
  using fom_rhs_data::fomRhs_;

public:
  static constexpr bool isResidualPolicy_ = true;
  using typename fom_rhs_data::fom_rhs_t;

public:
  LSPGSteadyResidualPolicy() = delete;
  ~LSPGSteadyResidualPolicy() = default;

  LSPGSteadyResidualPolicy(const fom_states_data     & fomStates,
			   const fom_rhs_data	     & fomResids,
			   const fom_eval_rhs_policy & fomEvalRhsFunctor)
    : fom_states_data(fomStates),
      fom_rhs_data(fomResids),
      fom_eval_rhs_policy(fomEvalRhsFunctor){}

public:
  template <typename lspg_state_t,
	    typename lspg_residual_t,
	    typename fom_t>
  void operator()(const lspg_state_t	& romY,
		  lspg_residual_t	& romR,
  		  const fom_t		& app) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg residual");
#endif

    fom_states_data::template reconstructCurrentFomState(romY);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
    fom_eval_rhs_policy::evaluate(app, yFom_, romR);
    timer->stop("fom eval rhs");
#else
    fom_eval_rhs_policy::evaluate(app, yFom_, romR);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg residual");
#endif
  }

  template <typename lspg_state_t, typename fom_t>
  fom_rhs_t operator()(const lspg_state_t & romY,
			 const fom_t	    & app) const
  {
    (*this).template operator()(romY, fomRhs_, app);
    return fomRhs_;
  }

};//end class

}}//end namespace pressio::rom
#endif
