/*
//@HEADER
// ************************************************************************
//
// solvers_linear_solver.hpp
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

#ifndef SOLVERS_LINEAR_SOLVER_HPP_
#define SOLVERS_LINEAR_SOLVER_HPP_

#include "solvers/src/linear/impl/solvers_linear_eigen_direct.hpp"
#include "solvers/src/linear/impl/solvers_linear_eigen_iterative.hpp"
#include "solvers/src/linear/impl/solvers_linear_kokkos_direct_geqrf.hpp"
#include "solvers/src/linear/impl/solvers_linear_kokkos_direct_getrs.hpp"
#include "solvers/src/linear/impl/solvers_linear_kokkos_direct_potrs_lower.hpp"
#include "solvers/src/linear/impl/solvers_linear_kokkos_direct_potrs_upper.hpp"

namespace pressio{ namespace solvers{ namespace linear {

namespace impl{

template<typename tag, typename MatrixT, typename enable = void>
struct _LinearSolverSelector{
  using type = void;
};

template<typename tag, typename MatrixT>
struct _LinearSolverSelector<
  tag, MatrixT,
  mpl::enable_if_t<
    ::pressio::solvers::linear::details::traits<tag>::iterative and
    (::pressio::containers::meta::is_matrix_wrapper_eigen<MatrixT>::value or
     ::pressio::containers::meta::is_multi_vector_wrapper_eigen<MatrixT>::value)
    >
  >
{
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using solver_traits   = linear::details::traits<tag>;
  using type = ::pressio::solvers::linear::impl::EigenIterative<tag, MatrixT>;
};

template<typename tag, typename MatrixT>
struct _LinearSolverSelector<
  tag, MatrixT,
  mpl::enable_if_t<
    ::pressio::solvers::linear::details::traits<tag>::direct and
    (::pressio::containers::meta::is_matrix_wrapper_eigen<MatrixT>::value or
     ::pressio::containers::meta::is_multi_vector_wrapper_eigen<MatrixT>::value)
    >
  >
{
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using solver_traits   = linear::details::traits<tag>;
  using type = ::pressio::solvers::linear::impl::EigenDirect<tag, MatrixT>;
};

template<typename tag, typename MatrixT>
struct _LinearSolverSelector<
  tag, MatrixT,
  mpl::enable_if_t<
    ::pressio::solvers::linear::details::traits<tag>::direct and
    (::pressio::containers::meta::is_dense_matrix_wrapper_kokkos<MatrixT>::value or
     ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<MatrixT>::value)
    >
  >
{
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using solver_traits   = linear::details::traits<tag>;
  using type = ::pressio::solvers::linear::impl::KokkosDirect<tag, MatrixT>;
};

}// end namespace pressio::solvers::linear::impl


template <typename ... Args>
using Solver = typename impl::_LinearSolverSelector<Args...>::type;


}}}//end namespace pressio::solvers::linear
#endif
