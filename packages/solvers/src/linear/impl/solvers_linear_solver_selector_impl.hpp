/*
//@HEADER
// ************************************************************************
//
// solvers_linear_solver_selector_impl.hpp
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

#ifndef SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_SOLVER_SELECTOR_IMPL_HPP_
#define SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_SOLVER_SELECTOR_IMPL_HPP_

#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "solvers_linear_eigen_direct_impl.hpp"
#include "solvers_linear_eigen_iterative_impl.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "solvers_linear_kokkos_direct_geqrf_impl.hpp"
#include "solvers_linear_kokkos_direct_getrs_impl.hpp"
#include "solvers_linear_kokkos_direct_potrs_lower_impl.hpp"
#include "solvers_linear_kokkos_direct_potrs_upper_impl.hpp"
#endif

namespace pressio{ namespace solvers{ namespace linear { namespace impl{

template<typename tag, typename MatrixT, typename enable = void>
struct LinearSolverSelector{
  using type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<typename tag, typename MatrixT>
struct LinearSolverSelector<
  tag, MatrixT,
  mpl::enable_if_t<
    ::pressio::solvers::linear::details::traits<tag>::iterative and
    (::pressio::containers::predicates::is_dense_matrix_wrapper_eigen<MatrixT>::value or
     ::pressio::containers::predicates::is_sparse_matrix_wrapper_eigen<MatrixT>::value or
     ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<MatrixT>::value)
    >
  >
{
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using solver_traits   = linear::details::traits<tag>;
  using type = ::pressio::solvers::linear::impl::EigenIterative<tag, MatrixT>;
};

template<typename tag, typename MatrixT>
struct LinearSolverSelector<
  tag, MatrixT,
  mpl::enable_if_t<
    ::pressio::solvers::linear::details::traits<tag>::direct and
    (::pressio::containers::predicates::is_dense_matrix_wrapper_eigen<MatrixT>::value or
     ::pressio::containers::predicates::is_sparse_matrix_wrapper_eigen<MatrixT>::value or
     ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<MatrixT>::value)    >
  >
{
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using solver_traits   = linear::details::traits<tag>;
  using type = ::pressio::solvers::linear::impl::EigenDirect<tag, MatrixT>;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template<typename tag, typename MatrixT>
struct LinearSolverSelector<
  tag, MatrixT,
  mpl::enable_if_t<
    ::pressio::solvers::linear::details::traits<tag>::direct and
    (::pressio::containers::predicates::is_dense_matrix_wrapper_kokkos<MatrixT>::value or
     ::pressio::containers::predicates::is_multi_vector_wrapper_kokkos<MatrixT>::value)
    >
  >
{
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using solver_traits   = linear::details::traits<tag>;
  using type = ::pressio::solvers::linear::impl::KokkosDirect<tag, MatrixT>;
};
#endif

}}}}// end namespace pressio::solvers::linear::impl
#endif  // SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_SOLVER_SELECTOR_IMPL_HPP_
