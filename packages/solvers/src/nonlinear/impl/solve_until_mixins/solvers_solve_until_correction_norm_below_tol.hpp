/*
//@HEADER
// ************************************************************************
//
// solvers_solve_until_correction_norm_below_tol.hpp
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


#ifndef SOLVERS_NONLINEAR_IMPL_SOLVE_UNTIL_MIXINS_SOLVERS_SOLVE_UNTIL_CORRECTION_NORM_BELOW_TOL_HPP_
#define SOLVERS_NONLINEAR_IMPL_SOLVE_UNTIL_MIXINS_SOLVERS_SOLVE_UNTIL_CORRECTION_NORM_BELOW_TOL_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template<typename sc_t, typename T>
class SolveUntilCorrectionNormBelowTol
  : public T,
    public IterativeBase< SolveUntilCorrectionNormBelowTol<sc_t, T>, sc_t >
{
  using this_t = SolveUntilCorrectionNormBelowTol<sc_t, T>;
  using iterative_base_t = IterativeBase<this_t, sc_t>;
  // friend iterative_base_t so it can access my private methods
  friend iterative_base_t;
  using typename iterative_base_t::iteration_t;

  iteration_t iStep_ = {};
  #ifdef PRESSIO_ENABLE_DEBUG_PRINT
  ::pressio::solvers::nonlinear::impl::NonlinearLeastSquaresDefaultMetricsPrinter<sc_t> solverStatusPrinter = {};
  #endif
public:
  SolveUntilCorrectionNormBelowTol() = delete;

  template <typename ...Args>
  SolveUntilCorrectionNormBelowTol(Args &&... args)
    : T(std::forward<Args>(args)...){}

  template<typename system_t, typename state_t>
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  mpl::enable_if_t< !::pressio::containers::predicates::is_array_pybind<state_t>::value>
#else
  void
#endif
  solve(const system_t & sys, state_t & state)
  {
    this->solveImpl(sys, state);
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  template<typename system_t, typename state_t>
  mpl::enable_if_t<::pressio::containers::predicates::is_array_pybind<state_t>::value>
  solve(const system_t & sys, state_t & state)
  {
    // here we want to view the state since we want to modify its data,
    // which is numpy array owned by the user inside their Python code.
    // upon exit of this function, the original state is changed.
    ::pressio::containers::Vector<state_t> stateView(state, ::pressio::view());
    this->solveImpl(sys, stateView);
  }
#endif

private:
  template<typename system_t, typename state_t>
  void solveImpl(const system_t & sys, state_t & state)
  {
    sc_t initialNorm = {};
    sc_t relativeNorm = {};
    iStep_ = 0;

    // Order is as follows:
    //   1.) Compute operators at step n and obtain correction
    //   2.) Compute statistics pertaining to step n and print
    //   3.) check convergence criteria
    //   4.) Update state if needed

    while (++iStep_ <= iterative_base_t::maxIters_)
    {
      // 1.
      T::computeCorrection(sys, state);

      // 2.
      const auto correctionNorm = T::correctionNormCurrentCorrectionStep();
      if (iStep_==1) initialNorm = correctionNorm;
      relativeNorm = correctionNorm/initialNorm;

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
      solverStatusPrinter.givenCorrectionNormsPrintRest(*this, iStep_, correctionNorm, relativeNorm);
#endif

      // 3.
      if (correctionNorm < iterative_base_t::tolerance_)
	break;

      // 4.
      T::updateState(sys, state);
    }
  }

  iteration_t getNumIterationsExecutedImpl() const {
    return iStep_;
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVE_UNTIL_MIXINS_SOLVERS_SOLVE_UNTIL_CORRECTION_NORM_BELOW_TOL_HPP_