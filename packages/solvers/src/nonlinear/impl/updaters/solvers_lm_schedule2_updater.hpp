/*
//@HEADER
// ************************************************************************
//
// solvers_lm_schedule2_updater.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_LM_SCHEDULE2_UPDATER_HPP_
#define SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_LM_SCHEDULE2_UPDATER_HPP_

#include "solvers_lm_gain_factor.hpp"

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template<typename state_t>
class LMSchedule2Updater : public BaseUpdater
{
  using scalar_t = typename ::pressio::containers::details::traits<state_t>::scalar_t;

  LMGainFactor<state_t> gainFactorEval_;

  using cnst		   = pressio::utils::constants<scalar_t>;
  const scalar_t rho1_	   = static_cast<scalar_t>(0.2);
  const scalar_t rho2_     = static_cast<scalar_t>(0.8);
  const scalar_t beta_	   = cnst::two();
  const scalar_t gammaInv_ = cnst::one()/cnst::three();
  const scalar_t tau_	   = cnst::one();

public:
  LMSchedule2Updater() = delete;
  LMSchedule2Updater(LMSchedule2Updater const &) = default;
  LMSchedule2Updater & operator=(LMSchedule2Updater const &) = default;
  LMSchedule2Updater(LMSchedule2Updater &&) = default;
  LMSchedule2Updater & operator=(LMSchedule2Updater &&) = default;
  ~LMSchedule2Updater() = default;

  LMSchedule2Updater(const state_t & state)
    : gainFactorEval_(state){}

public:
  void resetForNewCall() final
  {
    gainFactorEval_.resetForNewCall();
  }

  template<typename system_t, typename solver_mixin_t>
  void updateState(const system_t & sys,
		   state_t & state,
		   solver_mixin_t & solver)
  {
    PRESSIOLOG_DEBUG("nonlinsolver: lm2 update");
    constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
    constexpr auto ten  = static_cast<scalar_t>(10);
    constexpr auto seven  = static_cast<scalar_t>(7);
    constexpr auto negSeven  = ::pressio::utils::constants<scalar_t>::negOne()*seven;
    const auto tenToSev  = std::pow(ten, seven);
    const auto tenToNegSev  = std::pow(ten, negSeven);

    const auto mu	    = solver.lmDampParam();
    const auto & correction = solver.correctionCRef();

    // *** compute gain factor (rho) ***
    const auto rho = gainFactorEval_.compute(sys, state, mu, solver);

    // *** update mu ***
    if (rho < rho1_){
      solver.setLMDampParam(std::min(mu*beta_, tenToSev) );
    }

    if (rho > rho2_){
      solver.setLMDampParam(std::max( tenToNegSev, mu*gammaInv_) );
    }

    if (rho > 0){
      ::pressio::ops::update(state, one, correction, one);
    }
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_LM_SCHEDULE2_UPDATER_HPP_
