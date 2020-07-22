/*
//@HEADER
// ************************************************************************
//
// solver.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_
#define SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template<typename T, typename sc_t>
class Solver
  : public T,
    public IterativeBase< Solver<T, sc_t>, sc_t >
{

  using this_t = Solver<T, sc_t>;
  using iterative_base_t = IterativeBase<this_t, sc_t>;
  friend iterative_base_t;
  using typename iterative_base_t::iteration_t;
  using printer_t = NonlinearLeastSquaresDefaultMetricsPrinter<sc_t>;

public:
  enum class stop
    {
     whenCorrectionAbsoluteNormBelowTolerance, // this is the default
     whenCorrectionRelativeNormBelowTolerance,
     whenResidualAbsoluteNormBelowTolerance,
     whenResidualRelativeNormBelowTolerance,
     whenGradientAbsoluteNormBelowTolerance,
     whenGradientRelativeNormBelowTolerance,
     afterMaxIters
    };

private:
  stop stopping_ = this_t::stop::whenCorrectionAbsoluteNormBelowTolerance;
  iteration_t iStep_ = {};
  std::array<sc_t, 6> norms_;
#ifdef PRESSIO_ENABLE_DEBUG_PRINT
   printer_t solverStatusPrinter = {};
#endif

public:
  Solver() = delete;

  template <typename ...Args>
  Solver(stop stopping, Args &&... args)
    : T(std::forward<Args>(args)...),
      stopping_(stopping_){}

  template <typename ...Args>
  Solver(Args &&... args)
    : T(std::forward<Args>(args)...){}

  void setStoppingCriterion(stop value){
    stopping_ = value;
  }

  stop getStoppingCriterion() const{
    return stopping_;
  }

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
    sc_t residualNorm0 = {};
    sc_t correctionNorm0 = {};
    sc_t gradientNorm0 = {};

    iStep_ = 0;
    while (++iStep_ <= iterative_base_t::maxIters_)
    {
      // 1.
      T::computeCorrection(sys, state);

      // 2.
      const auto correctionNorm = T::correctionNormCurrentCorrectionStep();
      const auto residualNorm	= T::residualNormCurrentCorrectionStep();
      const auto gradientNorm	= T::gradientNormCurrentCorrectionStep();
      if (iStep_==1) {
	residualNorm0   = residualNorm;
	correctionNorm0 = correctionNorm;
	gradientNorm0   = gradientNorm;
      }

      norms_[0] = correctionNorm;
      norms_[1] = correctionNorm/correctionNorm0;
      norms_[2] = residualNorm;
      norms_[3] = residualNorm/residualNorm0;
      norms_[4] = gradientNorm;
      norms_[5] = gradientNorm/gradientNorm0;

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
      solverStatusPrinter.print(*this, iStep_,
				norms_[0], norms_[1],
				norms_[2], norms_[3],
				norms_[4], norms_[5]);
#endif

      // 3.
      if (stop(iStep_)) break;

      // 4.
      T::updateState(sys, state);
    }
  }

  bool stop(const iteration_t & iStep) const
  {
    switch (stopping_)
      {
      case stop::afterMaxIters:
    	return iStep == iterative_base_t::maxIters_;

      case stop::whenCorrectionAbsoluteNormBelowTolerance:
    	return norms_[0] < iterative_base_t::tolerance_;
      case stop::whenCorrectionRelativeNormBelowTolerance:
    	return norms_[1] < iterative_base_t::tolerance_;

      case stop::whenResidualAbsoluteNormBelowTolerance:
    	return norms_[2] < iterative_base_t::tolerance_;
      case stop::whenResidualRelativeNormBelowTolerance:
    	return norms_[3] < iterative_base_t::tolerance_;

      case stop::whenGradientAbsoluteNormBelowTolerance:
    	return norms_[4] < iterative_base_t::tolerance_;
      case stop::whenGradientRelativeNormBelowTolerance:
    	return norms_[5] < iterative_base_t::tolerance_;

      default:
	return false;
      };
  }

  iteration_t getNumIterationsExecutedImpl() const {
    return iStep_;
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_
