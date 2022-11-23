/*
//@HEADER
// ************************************************************************
//
// solvers_armijo.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_ARMIJO_HPP_
#define SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_ARMIJO_HPP_

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template <typename StateType>
class ArmijoUpdater
{
  StateType trialState_;

public:
  ArmijoUpdater() = delete;
  ArmijoUpdater(ArmijoUpdater const &) = default;
  ArmijoUpdater & operator=(ArmijoUpdater const &) = default;
  ArmijoUpdater(ArmijoUpdater &&) = default;
  ArmijoUpdater & operator=(ArmijoUpdater &&) = default;
  ~ArmijoUpdater() = default;

  ArmijoUpdater(const StateType & state)
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
    PRESSIOLOG_DEBUG("Armijo update");

    setTrialStateToZero();

    // https://people.maths.ox.ac.uk/hauser/hauser_lecture2.pdf

    // Armijo rule says to backtrack alpha until:
    // f(x_k + alpha_l * p_k) - f(x_k) <= alpha_l * beta * dot(g_k, p_k)
    //
    // where
    // k = the GN step
    // l = indexes the Armijo backtracking stages
    // p_k is the correction at GN k-th step
    // g_k is the gradient at GN k-th step

    const auto & p_k   = solver.correctionCRef();
    const auto & g_k   = solver.gradientCRef();
    auto fx_k    = solver.residualNormCurrentCorrectionStep();

    using scalar_type = decltype(fx_k);
    constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
    auto alpha = static_cast<scalar_type>(1);
    const scalar_type beta_ = static_cast<scalar_type>(0.0001);
    fx_k = std::pow(fx_k, ::pressio::utils::Constants<scalar_type>::two());
    const auto gkDotpk = ::pressio::ops::dot(g_k, p_k);

    PRESSIOLOG_DEBUG("starting backtracking");
    scalar_type ftrial = {};
    while (true)
    {
      if (std::abs(alpha) <= static_cast<scalar_type>(0.001)){
	PRESSIOLOG_DEBUG("alpha = {:6e}, too small, exiting line search", alpha);
	break;
      }

      // update : trialState = x_k + alpha*p_k
      constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
      ::pressio::ops::update(trialState_, zero, state, one, p_k, alpha);

      // compute rhs_l = alpha_l * beta * dot(g_k, p_k)
      const auto rhs = alpha * beta_ * gkDotpk;

      // eval f(x_k + alpha_l * p_k)
      solver.residualNorm(system, trialState_, ftrial);
      ftrial = std::pow(ftrial, ::pressio::utils::Constants<scalar_type>::two());

      // lhs = f(x_k + alpha_l * p_k) - f(x_k)
      const auto lhs = ftrial - fx_k;
      PRESSIOLOG_DEBUG("alpha = {:5f}: (f_trial-f) = {:6e}, rhs = {:6e}", alpha, lhs, rhs);

      if (lhs <= rhs){
	PRESSIOLOG_DEBUG("condition satisfied: f_trial-f <= rhs, exiting");
	// solution update: state = state + alpha*p_k
	::pressio::ops::update(state, one, p_k, alpha);
	break;
      }

      // exit when abs(fytrail-fy) < eps, leave eps = 1e-14 for now
      // change later with some machine epsilon
      if (std::abs(lhs) <= 1e-14){
	PRESSIOLOG_DEBUG("obj. function change too small, terminating");
	break;
      }

      /* convectional way to backtrack:alpha_l+1 = 0.5 * alpha_l */
      alpha *= 0.5;

    }//while not done
  }

private:
  void setTrialStateToZero(){
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
    ::pressio::ops::fill(trialState_, zero);
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_UPDATERS_SOLVERS_ARMIJO_HPP_
