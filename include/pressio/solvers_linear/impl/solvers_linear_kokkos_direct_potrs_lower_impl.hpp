/*
//@HEADER
// ************************************************************************
//
// solvers_linear_kokkos_direct_potrs_lower_impl.hpp
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

#ifndef SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_KOKKOS_DIRECT_POTRS_LOWER_IMPL_HPP_
#define SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_KOKKOS_DIRECT_POTRS_LOWER_IMPL_HPP_

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#endif

#if defined PRESSIO_ENABLE_TPL_KOKKOS and defined KOKKOS_ENABLE_CUDA
#include <cuda_runtime.h>
#include <cusolverDn.h>
#endif

namespace pressio { namespace linearsolvers{ namespace impl{

template<typename MatrixType>
class KokkosDirect<::pressio::linearsolvers::direct::potrsL, MatrixType>
{
public:

  using solver_tag	    = ::pressio::linearsolvers::direct::potrsL;
  using solver_traits   = ::pressio::linearsolvers::Traits<solver_tag>;
  using this_type          = KokkosDirect<solver_tag, MatrixType>;
  using matrix_type	    = MatrixType;
  using scalar_type        = typename MatrixType::value_type;
  using exe_space       = typename MatrixType::traits::execution_space;

  static_assert( solver_traits::kokkos_enabled == true,
  		 "the native solver must suppport kokkos to use in KokkosDirect");
  static_assert( solver_traits::direct == true,
  		 "the native solver must be direct to use in KokkosDirect");

public:
  KokkosDirect() = default;
  KokkosDirect(const KokkosDirect &) = delete;
  ~KokkosDirect() = default;


// because this uses teuchos lapack wrapper
#ifdef PRESSIO_ENABLE_TPL_TRILINOS

  /*
   * enable if:
   * the matrix has layout left (i.e. column major)
   * T is a kokkos vector 
   * has host execution space
   * T and MatrixType have same execution space
   */
  template < typename _MatrixType = MatrixType, typename T>
  mpl::enable_if_t<
    mpl::is_same<typename _MatrixType::traits::array_layout, Kokkos::LayoutLeft>::value 
    and ::pressio::is_vector_kokkos<T>::value 
    /*and ::pressio::containers::details::traits<T>::has_host_execution_space and*/
    and mpl::is_same<typename T::traits::execution_space, typename _MatrixType::traits::execution_space>::value
  >
  solve(const _MatrixType & A, const T& b, T & y)
  {
    const auto Aext0 = A.extent(0);
    const auto Aext1 = A.extent(1);
    if (Aext0 != auxMat_.extent(0) or Aext1 != auxMat_.extent(1))
    {
      Kokkos::resize(auxMat_, Aext0, Aext1);
	  }

    Kokkos::deep_copy(auxMat_, A);
    this->solveAllowMatOverwrite(auxMat_, b, y);
  }


  /*
   * enable if:
   * the matrix has layout left (i.e. column major)
   * T is a kokkos vector 
   * has host execution space
   * T and MatrixType have same execution space
   */
  template < typename _MatrixType = MatrixType, typename T>
  mpl::enable_if_t<
    mpl::is_same<typename _MatrixType::traits::array_layout, Kokkos::LayoutLeft>::value 
    and ::pressio::is_vector_kokkos<T>::value 
    /*and ::pressio::containers::details::traits<T>::has_host_execution_space and*/
    and mpl::is_same<typename T::traits::execution_space, typename _MatrixType::traits::execution_space>::value
  >
  solveAllowMatOverwrite(_MatrixType & A, const T& b, T & y)
  {
    assert(A.extent(0) == b.extent(0) );
    assert(A.extent(1) == y.extent(0) );
    // potrs is for symmetric pos def
    assert(A.extent(0) == A.extent(1) );

    // only one rhs
    constexpr int nRhs = 1;
    // just use n, since rows == cols
    const auto n = A.extent(0);

    // Cholesky factorization
    int info = 0;
    lpk_.POTRF(uplo_, n, A.data(), n, &info);
    assert(info == 0);

    // we need to deep copy b into y and pass y
    // because we overwrite the RHS in place with the solution
    Kokkos::deep_copy(y, b);

    lpk_.POTRS(uplo_, n, nRhs, A.data(), n, y.data(), y.extent(0), &info);
    assert(info == 0);
  }
#endif

  const char uplo_ = 'L';

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  Teuchos::LAPACK<int, scalar_type> lpk_;
  MatrixType auxMat_ = {};
#endif
};

}}} // end namespace pressio::linearsolvers::impl
#endif  // SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_KOKKOS_DIRECT_POTRS_LOWER_IMPL_HPP_
