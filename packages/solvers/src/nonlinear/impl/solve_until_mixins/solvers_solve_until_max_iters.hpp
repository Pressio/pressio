/*
//@HEADER
// ************************************************************************
//
// solvers_solve_until_max_iters.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_SOLVE_UNTIL_MIXINS_SOLVERS_SOLVE_UNTIL_MAX_ITERS_HPP_
#define SOLVERS_NONLINEAR_IMPL_SOLVE_UNTIL_MIXINS_SOLVERS_SOLVE_UNTIL_MAX_ITERS_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template<typename sc_t, typename T>
class SolveUntilMaxIters
  : public T,
    public IterativeBase< SolveUntilMaxIters<sc_t, T>, sc_t >
{
  using this_t = SolveUntilMaxIters<sc_t, T>;
  using iterative_base_t = IterativeBase<this_t, sc_t>;
  // friend iterative_base_t so it can access my private methods
  friend iterative_base_t;
  using typename iterative_base_t::iteration_t;
  #ifdef PRESSIO_ENABLE_DEBUG_PRINT
  ::pressio::solvers::nonlinear::impl::NonlinearLeastSquaresDefaultMetricsPrinter<sc_t> solverStatusPrinter = {};
  #endif 

public:
  SolveUntilMaxIters() = delete;

  template <typename ...Args>
  SolveUntilMaxIters(Args &&... args)
    : T(std::forward<Args>(args)...){}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    iteration_t iStep = 0;
    T::computeCorrection(sys, state);
    while (++iStep <= iterative_base_t::maxIters_)
    {
  #ifdef PRESSIO_ENABLE_DEBUG_PRINT
      solverStatusPrinter.print(*this, iStep);
  #endif 
      T::updateState(sys, state);
      T::computeCorrection(sys, state);

    }
  }

private:
  iteration_t getNumIterationsExecutedImpl() const {
    return iterative_base_t::maxIters_;
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVE_UNTIL_MIXINS_SOLVERS_SOLVE_UNTIL_MAX_ITERS_HPP_
