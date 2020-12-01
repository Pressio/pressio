/*
//@HEADER
// ************************************************************************
//
// solvers_armijo_updater.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_ARMIJO_UPDATER_HPP_
#define SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_ARMIJO_UPDATER_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template <typename state_t>
class ArmijoUpdater : public BaseUpdater
{
  using sc_t = typename ::pressio::containers::details::traits<state_t>::scalar_t;

  state_t trialState_;
  const sc_t beta_  = 0.0001;

public:
  ArmijoUpdater() = delete;
  ArmijoUpdater(ArmijoUpdater const &) = default;
  ArmijoUpdater & operator=(ArmijoUpdater const &) = default;
  ArmijoUpdater(ArmijoUpdater &&) = default;
  ArmijoUpdater & operator=(ArmijoUpdater &&) = default;
  ~ArmijoUpdater() = default;

  ArmijoUpdater(const state_t & state)
    : trialState_(state)
  {
    constexpr auto zero = ::pressio::utils::constants<sc_t>::zero();
    ::pressio::ops::fill(trialState_, zero);
  }

public:
  void resetForNewCall() final{ /* no op */ }

  template<typename system_t, typename solver_mixin_t>
  void updateState(const system_t & system,
		   state_t & state,
		   solver_mixin_t & solver)
  {
    PRESSIOLOG_DEBUG("armijo update");

    constexpr auto one = ::pressio::utils::constants<sc_t>::one();
    auto alpha = static_cast<sc_t>(1);

    ::pressio::ops::fill(trialState_,
			 ::pressio::utils::constants<sc_t>::zero());

    // https://people.maths.ox.ac.uk/hauser/hauser_lecture2.pdf

    // Armijo rule says to backtrack alpha until:
    // f(x_k + alpha_l * p_k) - f(x_k) <= alpha_l * beta * dot(g_k, p_k)
    //
    // where
    // k = the GN step
    // l = indexes the armijo backtracking stages
    // p_k is the correction at GN k-th step
    // g_k is the gradient at GN k-th step

    const auto & p_k   = solver.correctionCRef();
    const auto & g_k   = solver.gradientCRef();
    auto fx_k    = solver.residualNormCurrentCorrectionStep();
    fx_k = std::pow(fx_k, ::pressio::utils::constants<sc_t>::two());
    const auto gkDotpk = ::pressio::ops::dot(g_k, p_k);

    PRESSIOLOG_DEBUG("starting backtracking"); //: alpha = {}", alpha);
    sc_t ftrial = {};
    bool done = false;
    while (not done)
    {
      // update : trialState = x_k + alpha*p_k
      ::pressio::ops::update(trialState_, state, one, p_k, alpha);

      // compute rhs_l = alpha_l * beta * dot(g_k, p_k)
      const auto rhs = alpha * beta_ * gkDotpk;

      // eval f(x_k + alpha_l * p_k)
      solver.residualNorm(system, trialState_, ftrial);
      ftrial = std::pow(ftrial, ::pressio::utils::constants<sc_t>::two());

      // lhs = f(x_k + alpha_l * p_k) - f(x_k)
      const auto lhs = ftrial - fx_k;

      //::pressio::utils::io::print_stdout(" f(x_k+alpha*p_k) =", ftrial, "\n");
      PRESSIOLOG_DEBUG("alpha = {:5f}: (f_trial-f) = {:6e}, rhs = {:6e}", alpha, lhs, rhs);

      if (lhs <= rhs){
	PRESSIOLOG_DEBUG("condition satisfied: f_trial-f <= rhs, exiting");
	done = true;
      }

      // exit when abs(fytrail-fy) < eps, leave eps = 1e-14 for now
      // change later with some machine epsilon
      if (std::abs(lhs) <= 1e-14){
	PRESSIOLOG_DEBUG("change of obj. function too small, terminating");
	done = true;
      }

      /* convectional way to backtrack:alpha_l+1 = 0.5 * alpha_l */
      if (!done) alpha *= 0.5;
    }//while

    // solution update: state = state + alpha*p_k
    ::pressio::ops::update(state, one, p_k, alpha);
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_ARMIJO_UPDATER_HPP_
