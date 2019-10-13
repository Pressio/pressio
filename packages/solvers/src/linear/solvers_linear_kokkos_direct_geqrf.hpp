/*
//@HEADER
// ************************************************************************
//
// solvers_linear_kokkos_direct_geqrf.hpp
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

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#ifndef SOLVERS_LINEAR_KOKKOS_DIRECT_GEQRF_HPP
#define SOLVERS_LINEAR_KOKKOS_DIRECT_GEQRF_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../base/solvers_linear_base.hpp"
#include "solvers_linear_traits.hpp"
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_BLAS.hpp>
#endif

#if defined PRESSIO_ENABLE_TPL_KOKKOS and defined KOKKOS_ENABLE_CUDA
#include <cuda_runtime.h>
#include <cusolverDn.h>
#endif

namespace pressio { namespace solvers { namespace direct{

template<typename MatrixT>
class KokkosDirect<::pressio::solvers::linear::direct::geqrf, MatrixT>
  : public LinearBase<
  ::pressio::solvers::linear::direct::geqrf,
  MatrixT,
  KokkosDirect<::pressio::solvers::linear::direct::geqrf, MatrixT>
  >
{
public:
  static_assert( ::pressio::containers::meta::is_dense_matrix_wrapper_kokkos<MatrixT>::value or
  		 ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<MatrixT>::value,
  		 "Kokkos direct dense solver expects either (a) dense matrix wrapper or a (b) multi-vector wrapper, both wrapping a rank=2 Kokkos View");

  using solver_t	= ::pressio::solvers::linear::direct::geqrf;
  using this_t          = KokkosDirect<solver_t, MatrixT>;
  using matrix_type	= MatrixT;
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using scalar_t        = typename containers::details::traits<MatrixT>::scalar_t;
  using exe_space       = typename containers::details::traits<MatrixT>::execution_space;

  using base_t  = LinearBase<solver_t, MatrixT, this_t>;
  friend base_t;
  using solver_traits   = linear::details::traits<solver_t>;

  static_assert( solver_traits::kokkos_enabled == true,
  		 "the native solver must suppport kokkos to use in KokkosDirect");
  static_assert( solver_traits::direct == true,
  		 "the native solver must be direct to use in KokkosDirect");

public:

  KokkosDirect(){
#ifdef KOKKOS_ENABLE_CUDA
    auto cusolverStatus = cusolverDnCreate(&cuDnHandle_);
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
#endif
  }

  KokkosDirect(const KokkosDirect &) = delete;

  ~KokkosDirect(){
#ifdef KOKKOS_ENABLE_CUDA
    auto cusolverStatus = cusolverDnDestroy(cuDnHandle_);
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
#endif
  }


private:

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  /*
   * enable if:
   * the matrix has layout left (i.e. column major)
   * T is a kokkos vector wrapper
   * T has host execution space
   * T and MatrixT have same execution space
   */
  template <
    typename _MatrixT = MatrixT,
    typename T,
    mpl::enable_if_t<
      mpl::is_same<
	typename ::pressio::containers::details::traits<_MatrixT>::layout,
	Kokkos::LayoutLeft
	>::value
      and
      ::pressio::containers::meta::is_vector_wrapper_kokkos<T>::value
      and
      ::pressio::containers::details::traits<T>::has_host_execution_space
      and
      mpl::is_same<
	typename containers::details::traits<T>::execution_space,
	typename containers::details::traits<_MatrixT>::execution_space
	>::value
      > * = nullptr
  >
  void solveAllowMatOverwriteImpl(_MatrixT & A, const T& b, T & y)
  {
    assert(A.rows() == b.size() );
    assert(A.cols() == y.size() );
    // gerts is for square matrices
    assert(A.rows() == A.cols() );

    // only one rhs because this is only enabled if T is a vector wrapper
    constexpr int nRhs = 1;

    // just use n, since rows == cols
    const auto n = A.rows();

    // to store the return code of the function
    int info = 0;

    if (tau_.size() != n)
      tau_.resize(n);

    if (lwork_ == -1){
      // this means we need to query what lwork should be
      // and resize the work array "work_" properly
      lpk_.GEQRF(n, n, A.data()->data(), n, tau_.data(), work_.data(), lwork_, &info);
      work_.resize(work_[0]);
      lwork_ = work_[0];
    }
    // do QR
    lpk_.GEQRF(n, n, A.data()->data(), n, tau_.data(), work_.data(), lwork_, &info);
    assert(info == 0);
    // std::cout << " info " << info
    // 	      << " lwork " << lwork_
    // 	      << " workSz " << work_.size() << std::endl;

    // Note that A is overwritten starting from here!

    // we need to deep copy b into y and pass y to ormqr
    // because it is overwritten with Q^T b
    Kokkos::deep_copy(*y.data(), *b.data());

    // compute Q^T b
    const char side = 'L';
    const char trans = 'T';
    lpk_.ORMQR(side, trans, n, nRhs, n,
    	       A.data()->data(),
    	       n,
    	       tau_.data(),
    	       y.data()->data(),
    	       n,
    	       work_.data(),
    	       lwork_,
    	       &info);
    //std::cout << " info-ormqr " << info << std::endl;

    // solver R y = Q^T b
    const char uplo = 'U';
    const char diag = 'N';
    constexpr scalar_t alpha = ::pressio::utils::constants::one<scalar_t>();
    blas_.TRSM(Teuchos::ESide::LEFT_SIDE,
    	       Teuchos::EUplo::UPPER_TRI,
    	       Teuchos::ETransp::NO_TRANS,
    	       Teuchos::EDiag::NON_UNIT_DIAG,
    	       n, nRhs, alpha,
    	       A.data()->data(), n, y.data()->data(), n);
  }
#endif


private:

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  Teuchos::LAPACK<int, scalar_t> lpk_;
  Teuchos::BLAS<int, scalar_t> blas_;

  // lwork is used for the geqrf, but if you read the doc,
  // if lwork == -1, then geqrf does a query of what is needed.
  // more details are shown in the code above on how we use lwork
  int lwork_ = -1;

  std::vector<scalar_t> work_ = {0};
  std::vector<scalar_t> tau_ = {};
#endif

#if defined PRESSIO_ENABLE_TPL_KOKKOS and defined KOKKOS_ENABLE_CUDA
  cusolverDnHandle_t cuDnHandle_;
#endif
};

}}} // end namespace pressio::solvers::direct
#endif
#endif
