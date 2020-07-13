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

#ifndef SOLVERS_NONLINEAR_IMPL_UPDATE_MIXINS_SOLVERS_LM_SCHEDULE2_UPDATER_HPP_
#define SOLVERS_NONLINEAR_IMPL_UPDATE_MIXINS_SOLVERS_LM_SCHEDULE2_UPDATER_HPP_

#include "solvers_lm_gain_factor.hpp"

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template<typename scalar_t, typename state_t, typename T>
class LMSchedule2Updater : public T
{
private:
  LMGainFactor<scalar_t, state_t> gainFactorEval_;

  using cnst		   = pressio::utils::constants<scalar_t>;
  const scalar_t rho1_	   = static_cast<scalar_t>(0.2);
  const scalar_t rho2_     = static_cast<scalar_t>(0.8);
  const scalar_t beta_	   = cnst::two();
  const scalar_t gammaInv_ = cnst::one()/cnst::three();
  const scalar_t tau_	   = cnst::one();

public:
  template <typename system_t, typename ...Args>
  LMSchedule2Updater(const system_t & sys, const state_t & state, Args &&... args)
    : gainFactorEval_(state), T(sys, state, std::forward<Args>(args)...){}

  template<typename system_t>
  void updateState(const system_t & sys, state_t & state)
  {
    constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
    constexpr auto two  = ::pressio::utils::constants<scalar_t>::two();
    constexpr auto ten  = static_cast<scalar_t>(10);
    constexpr auto seven  = static_cast<scalar_t>(7);
    constexpr auto negSeven  = -seven;
    constexpr auto tenToSev  = std::pow(ten, seven);
    constexpr auto tenToNegSev  = std::pow(ten, negSeven);

    const auto mu	    = T::getLMDampParam();
    const auto & correction = T::viewCorrection();

    // *** compute gain factor (rho) ***
    const auto rho = gainFactorEval_.compute(sys, state, mu, *this);

    // *** update mu ***
    if (rho < rho1_){
      T::setLMDampParam( std::min(mu*beta_, tenToSev) );
    }

    if (rho > rho2_){
      T::setLMDampParam( std::max( tenToNegSev, mu*gammaInv_) );
    }

    if (rho > 0){
      ::pressio::ops::do_update(state, one, correction, one);
    }
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_UPDATE_MIXINS_SOLVERS_LM_SCHEDULE2_UPDATER_HPP_
