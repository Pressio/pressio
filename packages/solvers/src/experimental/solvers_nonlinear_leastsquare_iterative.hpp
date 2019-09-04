/*
//@HEADER
// ************************************************************************
//
// solvers_nonlinear_leastsquare_iterative.hpp
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
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_LEASTSQUARE_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_LEASTSQUARE_ITERATIVE_HPP

#include <type_traits>
#include "solvers_nonlinear_base.hpp"
#include "solvers_nonlinear_iterative_helper.hpp"
#include "../solvers_ConfigDefs.hpp"


namespace pressio {
namespace solvers {

class NonLinearSolvers; // Fwd declaration

template <
  typename PolicyT,
  typename LSolverT
>
class NonLinearLeastSquareIterativeSolver
  : public NonLinearIterativeSolverHelper,
    public NonLinearSolverBase<
      NonLinearLeastSquareIterativeSolver<
        PolicyT,
        LSolverT
      >
    >
{

  private:

    friend NonLinearSolvers;
    typedef NonLinearSolverBase<NonLinearLeastSquareIterativeSolver<PolicyT, LSolverT>> base_type;


  public:


    /**
     * Implements the solve method for a non linear solver.
     */
    template <
      typename PrecT,
      typename NormT,
      typename SystemT,
      typename VectorT
    >
    auto solve_(const SystemT& sys, const VectorT& b) 
    -> decltype(
        PolicyT::template solve<LSolverT, PrecT, NormT>(sys, b, 1, 1, 1e-2, 1e-2, 1.)
      ){

      double tolerance = this->getTolerance();
      double nonLinearTolerance = this->getNonLinearTolerance();
      containers::default_types::uint maxIterations = this->getMaxIterations();
      containers::default_types::uint maxNonLinearIterations = this->getMaxNonLinearIterations();

      return PolicyT::template solve<LSolverT, PrecT, NormT>(sys, b, maxIterations, 
                maxNonLinearIterations, tolerance, nonLinearTolerance, lambda_);
    }


    /**
     * Get the value of lambda used in nonlinear least-square algorithms.
     */
    double getLambda() {
      return lambda_;
    }


  protected:

    NonLinearLeastSquareIterativeSolver() : NonLinearIterativeSolverHelper(), base_type(), lambda_(1.0) {}


  private:

    double lambda_;
};

} // end namespace solvers
} // end namespace pressio

#endif
