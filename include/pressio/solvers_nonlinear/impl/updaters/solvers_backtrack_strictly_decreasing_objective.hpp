/*
//@HEADER
// ************************************************************************
//
// solvers_backtrack_strictly_decreasing_objective.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_BACKTRACK_STRICLTY_DECREASING_OBJECTIVE_HPP_
#define SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_BACKTRACK_STRICLTY_DECREASING_OBJECTIVE_HPP_

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template <typename StateType>
class BacktrackStrictlyDecreasingObjectiveUpdater
{
  StateType trialState_;

public:
  BacktrackStrictlyDecreasingObjectiveUpdater() = delete;
  BacktrackStrictlyDecreasingObjectiveUpdater(BacktrackStrictlyDecreasingObjectiveUpdater const &) = default;
  BacktrackStrictlyDecreasingObjectiveUpdater & operator=(BacktrackStrictlyDecreasingObjectiveUpdater const &) = default;
  BacktrackStrictlyDecreasingObjectiveUpdater(BacktrackStrictlyDecreasingObjectiveUpdater &&) = default;
  BacktrackStrictlyDecreasingObjectiveUpdater & operator=(BacktrackStrictlyDecreasingObjectiveUpdater &&) = default;
  ~BacktrackStrictlyDecreasingObjectiveUpdater() = default;

  BacktrackStrictlyDecreasingObjectiveUpdater(const StateType & state)
    : trialState_(::pressio::ops::clone(state))
  {
    setTrialStateToZero();
  }

public:
  void reset(){ /* no op */ }

  template<typename system_t, typename solver_mixin_t>
  void operator()(const system_t & system,
		  StateType & state,
		  solver_mixin_t & solver)
  {
    PRESSIOLOG_DEBUG("BacktrackStrictlyDecreasingObjective update");
    setTrialStateToZero();
    const auto & p_k   = solver.correctionCRef();
    auto fx_k    = solver.residualNormCurrentCorrectionStep();

    using scalar_type = decltype(fx_k);
    constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
    auto alpha = static_cast<scalar_type>(1);
    fx_k = std::pow(fx_k, ::pressio::utils::Constants<scalar_type>::two());

    PRESSIOLOG_DEBUG("starting backtracking");
    scalar_type ftrial = {};
    bool done = false;
    while (not done)
    {
      if (std::abs(alpha) <= 0.001){
        /*
        Presently set an exit alpha drops below 0.001; anything smaller 
        than this is probably unreasonable. Note that this quantity doesn't depend on 
        the dimension or magnitude of the state
        */
	PRESSIOLOG_DEBUG("alpha = {:6e}, too small, exiting line search", alpha);
        alpha = 0.;
	done = true;
      }

      // update : trialState = x_k + alpha*p_k
      constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
      ::pressio::ops::update(trialState_, zero, state, one, p_k, alpha);

      solver.residualNorm(system, trialState_, ftrial);
      ftrial = std::pow(ftrial, ::pressio::utils::Constants<scalar_type>::two());

      if (ftrial <= fx_k){
	PRESSIOLOG_DEBUG("backtrack; condition satisfied with alpha= {:6e}", alpha);
	done = true;
        ::pressio::ops::update(state, one, p_k, alpha);
      }
      /* convectional way to backtrack:alpha_l+1 = 0.5 * alpha_l */
      if (!done) alpha *= 0.5;
    }//while
  }

private:
  void setTrialStateToZero(){
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
    ::pressio::ops::fill(trialState_, zero);
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_BACKTRACK_STRICLTY_DECREASING_OBJECTIVE_HPP_
