/*
//@HEADER
// ************************************************************************
//
// solvers_linear_wrapper_armadillo.hpp
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
#ifndef SOLVERS_LINEAR_WRAPPER_ARMADILLO_HPP
#define SOLVERS_LINEAR_WRAPPER_ARMADILLO_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../../../containers/src/meta/containers_meta_detection_idiom.hpp"


namespace pressio {
namespace solvers {

template <typename MatrixT>
class SolversLinearDirectWrapperArmadillo {

  public:

    /**
     * Update the matrix that describes the linear system.
     */
    void resetLinearSystem(const MatrixT& other) {
      A = other;
    }


    /**
     * Solve the linear system
     */
    template <typename VectorT>
    VectorT solve(const VectorT& b) {
      return this->template _solveImpl<MatrixT>(b);
    }


  private:

    /**
     * Implement the solve method for a dense matrix.
     */
    template <
      typename DMatrixT,
      typename VectorT,
      typename ::pressio::mpl::enable_if_t<
        containers::details::traits<DMatrixT>::is_dense
      >* = nullptr
    >
    VectorT _solveImpl(const VectorT& b) const {
      return b;//VectorT(arma::solve(*A.data(), *b.data()));
    }


    /**
     * Implement the solve method for a sparse matrix.
     */
    template <
      typename DMatrixT,
      typename VectorT,
      typename ::pressio::mpl::enable_if_t<
        containers::details::traits<DMatrixT>::is_sparse
      >* = nullptr
    >
    VectorT _solveImpl(const VectorT& b) const {
      return b; // VectorT(arma::spsolve(*A.data(), *b.data()));
    }


  private:
    MatrixT A;
};


} // end namespace solvers
} // end namespace pressio

#endif
