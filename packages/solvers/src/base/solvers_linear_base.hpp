/*
//@HEADER
// ************************************************************************
//
// solvers_linear_base.hpp
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

#ifndef SOLVERS_LINEAR_BASE_HPP
#define SOLVERS_LINEAR_BASE_HPP

namespace pressio{ namespace solvers{

/**
 * Base class for linear solver
 *
 * @section DESCRIPTION
 *
 * This class defines the public interface for a linear solver class.
 */
template<typename MatrixT, typename Derived>
struct LinearBase {

  LinearBase() = default;
  LinearBase(const LinearBase&) = delete;
  ~LinearBase() = default;

  template <typename CompatibleMatrixT>
  void resetLinearSystem(const CompatibleMatrixT& A)
  {
    static_assert(solvers::meta::are_matrix_compatible<
		  MatrixT, CompatibleMatrixT>::value,
		  "Matrix types not compatible");

    static_cast<Derived&>(*this).resetLinearSystemImpl(A);
  }

  template <typename VectorT>
  void solve(const VectorT & b, VectorT& x) {
    static_cast<Derived&>(*this).solveImpl(b, x);
  }

  template <typename VectorT>
  void solve(const MatrixT & A, const VectorT & b, VectorT& x) {
    static_cast<Derived&>(*this).solveImpl(A, b, x);
  }

  // this method allows solver to overwrite matrix
  template <typename VectorT>
  void solveAllowMatOverwrite(MatrixT & A, const VectorT & b, VectorT& x) {
    static_cast<Derived&>(*this).solveAllowMatOverwriteImpl(A, b, x);
  }

};

}}//end namespace pressio::solvers
#endif
