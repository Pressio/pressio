/*
//@HEADER
// ************************************************************************
//
// solvers_lm_schedule1.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_LM_SCHEDULE1_HPP_
#define SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_LM_SCHEDULE1_HPP_

#include "solvers_lm_gain_factor.hpp"

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template<typename StateType>
class LMSchedule1Updater
{
  using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
  LMGainFactor<StateType> gainFactorEval_;

  using cnst = pressio::utils::Constants<scalar_type>;
  const scalar_type beta_     = cnst::two();
  const scalar_type gammaInv_ = cnst::one()/cnst::three();
  const scalar_type p_	   = cnst::three();
  const scalar_type tau_	   = cnst::one();
  scalar_type nu_		   = cnst::two();

public:
  LMSchedule1Updater() = delete;
  LMSchedule1Updater(LMSchedule1Updater const &) = default;
  LMSchedule1Updater & operator=(LMSchedule1Updater const &) = default;
  LMSchedule1Updater(LMSchedule1Updater &&) = default;
  LMSchedule1Updater & operator=(LMSchedule1Updater &&) = default;
  ~LMSchedule1Updater() = default;

  LMSchedule1Updater(const StateType & state)
    : gainFactorEval_(state){}

public:
  void reset(){
    PRESSIOLOG_DEBUG("lm1 update reset");
    nu_ = cnst::two();
    gainFactorEval_.resetForNewCall();
  }

  template<typename SystemType, typename SolverMixinType>
  void operator()(const SystemType & sys,
		   StateType & state,
		   SolverMixinType & solver)
  {
    PRESSIOLOG_DEBUG("nonlinsolver: lm1 update");
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    constexpr auto two  = ::pressio::utils::Constants<scalar_type>::two();

    const auto mu	    = solver.lmDampParam();
    const auto & correction = solver.correctionCRef();

    // *** compute gain factor (rho) ***
    const auto rho = gainFactorEval_.compute(sys, state, mu, solver);

    // *** update mu ***
    if (rho > 0){
      PRESSIOLOG_DEBUG("lm1 update: rho>0");
      ::pressio::ops::update(state, one, correction, one);
      scalar_type mu_rat = one - (beta_ - one)*std::pow(two*rho - one, p_);
      mu_rat = std::max(mu_rat, gammaInv_);
      solver.setLMDampParam(mu*mu_rat);
      nu_ = beta_;
      PRESSIOLOG_DEBUG("lm1 update: rho>0, {}, {}", nu_, mu*mu_rat);
    }
    else{
      solver.setLMDampParam(mu*nu_);
      nu_ = two*nu_;
      PRESSIOLOG_DEBUG("lm1 update: rho<=0, {}, {}", nu_, mu*nu_);
    }
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_LM_SCHEDULE1_HPP_
