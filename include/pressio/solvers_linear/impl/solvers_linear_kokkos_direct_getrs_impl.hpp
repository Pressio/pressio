/*
//@HEADER
// ************************************************************************
//
// solvers_linear_kokkos_direct_getrs_impl.hpp
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

#ifndef SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_KOKKOS_DIRECT_GETRS_IMPL_HPP_
#define SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_KOKKOS_DIRECT_GETRS_IMPL_HPP_

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
class KokkosDirect<::pressio::linearsolvers::direct::getrs, MatrixType>
{
public:

  using solver_tag	    = ::pressio::linearsolvers::direct::getrs;
  using this_type          = KokkosDirect<solver_tag, MatrixType>;
  using matrix_type	    = MatrixType;
  using scalar_type        = typename MatrixType::value_type;
  using exe_space       = typename MatrixType::traits::execution_space;
  using solver_traits   = ::pressio::linearsolvers::Traits<solver_tag>;

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
    /*::pressio::containers::details::traits<T>::has_host_execution_space and*/
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
    /*::pressio::containers::details::traits<T>::has_host_execution_space and*/
    and mpl::is_same<typename T::traits::execution_space, typename _MatrixType::traits::execution_space>::value
  >
  solveAllowMatOverwrite(_MatrixType & A, const T& b, T & y)
  {
    assert(A.extent(0) == b.extent(0) );
    assert(A.extent(1) == y.extent(0) );
    // gerts is for square matrices
    assert(A.extent(0) == A.extent(1) );

    // only one rhs
    constexpr int nRhs = 1;
    // just use n, since rows == cols
    const auto n = A.extent(0);

    int info = 0;
    const int ipivSz = n;
    // this needs to be moved out, does not make sense to construct it every time
    std::vector<int> ipiv(ipivSz);

    // LU factorize using GETRF
    lpk_.GETRF(n, n, A.data(), n, ipiv.data(), &info);
    assert(info == 0);

    // we need to deep copy b into y and pass y
    // because getrs overwrite the RHS in place with the solution
    Kokkos::deep_copy(y, b);

    const char ct = 'N';
    lpk_.GETRS(ct, n, nRhs, A.data(), n, ipiv.data(), y.data(), y.extent(0), &info);
    assert(info == 0);
  }
#endif


#if defined PRESSIO_ENABLE_TPL_KOKKOS and defined KOKKOS_ENABLE_CUDA
  /*
   * enable if:
   * the matrix has layout left (i.e. column major)
   * T is a kokkos vector 
   * has CUDA execution space
   * T and MatrixType have same execution space
   */
  template < typename _MatrixType = MatrixType, typename T>
  mpl::enable_if_t<
    mpl::is_same<typename wrapped_type::traits::array_layout, Kokkos::LayoutLeft>::value 
    and ::pressio::is_vector_kokkos<T>::value 
    /*and ::pressio::containers::details::traits<T>::has_host_execution_space and*/
    and mpl::is_same<typename T::traits::execution_space, typename _MatrixType::traits::execution_space>::value
    >
  solveAllowMatOverwrite(_MatrixType & A, const T& b, T & y)
  {
    assert(A.extent(0) == b.extent(0) );
    assert(A.extent(1) == y.extent(0) );
    // gerts is for square matrices
    assert(A.extent(0) == A.extent(1) );

    // only one rhs
    constexpr int nRhs = 1;

    // n = nRows = nCols: because it is square matrix
    const auto n = A.extent(0);

    cusolverStatus_t cusolverStatus;
    //cusolverDnHandle_t handle;
    cudaError cudaStatus;
    int Lwork = 0;

    // cuDnHandle already created in constructor

    // compute buffer size and prep.memory
    cusolverStatus = cusolverDnDgetrf_bufferSize(cuDnHandle_, n, n, A.data(), n, &Lwork);
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);

    // for now, working buffers are stored as Kokkos arrays but
    // maybe later we can use directly cuda allocations
    using k1d_d = Kokkos::View<scalar_type*, Kokkos::LayoutLeft, exe_space>;
    using k1di_d = Kokkos::View<int*, Kokkos::LayoutLeft, exe_space>;
    k1d_d work_d("d_work", Lwork);
    k1di_d pivot_d("d_pivot", n);
    k1di_d info_d("d_info", 1);

    cusolverStatus = cusolverDnDgetrf(cuDnHandle_, n, n, A.data(), n,
                                      work_d.data(),
				                              pivot_d.data(),
    				                          info_d.data());
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);

    // we need to deep copy b into y and pass y
    // because getrs overwrites the RHS in place with the solution
    Kokkos::deep_copy(y, b);

    cusolverStatus = cusolverDnDgetrs(cuDnHandle_, CUBLAS_OP_N,
    				      n, nRhs,
    				      A.data(), n,
    				      pivot_d.data(),
    				      y.data(), n,
    				      info_d.data());
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);

    // make sure the solver kernel is done before exiting
    cudaStatus = cudaDeviceSynchronize();
  }
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  Teuchos::LAPACK<int, scalar_type> lpk_;

  MatrixType auxMat_ = {};
#endif

#if defined PRESSIO_ENABLE_TPL_KOKKOS and defined KOKKOS_ENABLE_CUDA
  cusolverDnHandle_t cuDnHandle_;
#endif
};

}}} // end namespace pressio::linearsolvers::impl
#endif  // SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_KOKKOS_DIRECT_GETRS_IMPL_HPP_
