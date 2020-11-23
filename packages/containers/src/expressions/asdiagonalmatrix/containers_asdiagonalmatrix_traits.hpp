/*
//@HEADER
// ************************************************************************
//
// containers_expressions_traits.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_ASDIAGONALMATRIX_TRAITS_HPP_
#define CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_ASDIAGONALMATRIX_TRAITS_HPP_

namespace pressio{ namespace containers{ namespace details{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename v_type>
struct traits<
  ::pressio::containers::expressions::AsDiagonalMatrixExpr<v_type>,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_dynamic_vector_wrapper_eigen<v_type>::value
    >
  >
  : public containers_shared_traits<
  ::pressio::containers::expressions::AsDiagonalMatrixExpr<v_type>,
  typename details::traits<v_type>::wrapped_t,
  //isVector=false, isMatrix=true, isMv=false
  false, true, false,
  WrappedPackageIdentifier::Eigen,
  true
  >,
    public matrix_shared_traits<false>
{
  static constexpr auto wrapped_matrix_identifier=WrappedMatrixIdentifier::DenseEigen;
  static constexpr bool is_static = true;
  static constexpr bool is_dynamic = false;
  using wrapped_t = typename traits<v_type>::wrapped_t;
  using scalar_t  = typename traits<v_type>::scalar_t;
  using ordinal_t = typename traits<v_type>::ordinal_t;
  using size_t    = ordinal_t;

  // conditiona ref type because native expression returns by value when object is const
  using reference_t = typename std::conditional<
    std::is_const<v_type>::value, scalar_t, scalar_t &
  >::type;

  using const_reference_t = typename std::conditional<
    std::is_const<v_type>::value, scalar_t, scalar_t const &
    >::type;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template <typename v_type>
struct traits<
  ::pressio::containers::expressions::AsDiagonalMatrixExpr<v_type>,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_vector_wrapper_tpetra<v_type>::value
    >
  >
  : public containers_shared_traits<
  ::pressio::containers::expressions::AsDiagonalMatrixExpr<v_type>,
  typename details::traits<v_type>::wrapped_t,
  //isVector=false, isMatrix=true, isMv=false
  false, true, false,
  WrappedPackageIdentifier::Trilinos,
  true
  >,
    public matrix_shared_traits<false>
{
  static constexpr auto wrapped_matrix_identifier=WrappedMatrixIdentifier::Tpetra;
  static constexpr bool is_static = true;
  static constexpr bool is_dynamic = false;
  using wrapped_t = typename traits<v_type>::wrapped_t;
  using scalar_t  = typename traits<v_type>::scalar_t;
  using local_ordinal_t  = typename traits<v_type>::local_ordinal_t;
  using global_ordinal_t = typename traits<v_type>::global_ordinal_t;
  using size_t = global_ordinal_t;
};
#endif

}}}//end namespace pressio::containers::details
#endif  // CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_ASDIAGONALMATRIX_TRAITS_HPP_
