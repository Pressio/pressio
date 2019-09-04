/*
//@HEADER
// ************************************************************************
//
// solvers_linear_direct.hpp
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

#ifndef SOLVERS_EXPERIMENTAL_LINEAR_DIRECT_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_DIRECT_HPP

#include "solvers_linear_base.hpp"


namespace pressio{
namespace solvers{

// Forward declarations
struct LinearSolvers;


/**
 * Class that implements a linear direct solver.
 *
 * @tparam SolverT linear dense solver.
 * @tparam MatrixT matrix defining the linear system.
 * @tparam PolicyT policy that implements the solution algorithm.
 */
template<
  typename SolverT,
  typename MatrixT
//  typename PolicyT
>
class LinearDirectSolver
  : public LinearSolverBase<
      SolverT,
      MatrixT,
      LinearDirectSolver<
        SolverT,
        MatrixT
      >
    >
{

  private:

    friend LinearSolvers;
    typedef LinearSolverBase<SolverT, MatrixT, LinearDirectSolver<SolverT, MatrixT>> base_type;


  public:

    LinearDirectSolver(LinearDirectSolver&& other) : base_type(std::move(other)) {}


    template <typename T>
    auto _solve(const T& b) 
      -> decltype(this->getSolver()->solve(b)) {
      return this->getSolver()->solve(b);
    }


  public:

    LinearDirectSolver() : base_type() {};

    LinearDirectSolver(std::shared_ptr<SolverT> solver) : base_type(solver) {};

};

} //end namespace solvers

}//end namespace pressio
#endif
