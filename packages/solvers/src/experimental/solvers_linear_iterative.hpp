/*
//@HEADER
// ************************************************************************
//
// solvers_linear_iterative.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_HPP

#include "../solvers_ConfigDefs.hpp"
#include "solvers_linear_base.hpp"


namespace pressio {
namespace solvers {

// Forward declarations
struct LinearSolvers;


/**
 * @brief Class that implements a linear iterative solver
 */
template<
  typename SolverT,
  typename MatrixT
>
class LinearIterativeSolver
  : public LinearSolverBase<
      SolverT,
      MatrixT,
      LinearIterativeSolver<
        SolverT,
        MatrixT
      >
    >
{

  private:

    friend LinearSolvers;
    typedef LinearSolverBase<SolverT, MatrixT, LinearIterativeSolver<SolverT, MatrixT>> base_type;


  public:

    LinearIterativeSolver(LinearIterativeSolver&& other) :
      base_type(std::move(other)), maxIters_(other.maxIters_), tolerance_(other.tolerance_) {}


    template <typename T>
    auto _solve(const T& b)
     -> decltype(this->getSolver()->solve(b)) {
      auto solver = this->getSolver();
      solver->setMaxIterations(this->getMaxIterations());
      solver->setTolerance(this->getTolerance());
      return solver->solve(b);
    }


    inline containers::default_types::uint getMaxIterations() {
      return maxIters_;
    }


    void setMaxIterations(containers::default_types::uint maxIters) {
      maxIters_ = maxIters;
    }


    inline double getTolerance() {
      return tolerance_;
    }


    void setTolerance(double tolerance) {
      tolerance_ = tolerance;
    }


  protected:

    LinearIterativeSolver() :
      base_type(), maxIters_(100), tolerance_(1.0e-6) {};


    LinearIterativeSolver(std::shared_ptr<SolverT> solver) :
      base_type(solver), maxIters_(100), tolerance_(1.0e-6) {};


  private:

    containers::default_types::uint maxIters_;
    double tolerance_;
};

} //end namespace solvers
} //end namespace pressio

#endif
