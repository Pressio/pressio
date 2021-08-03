/*
//@HEADER
// ************************************************************************
//
// containers_span_traits.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_SPAN_CONTAINERS_SPAN_TRAITS_HPP_
#define CONTAINERS_EXPRESSIONS_SPAN_CONTAINERS_SPAN_TRAITS_HPP_

namespace pressio{ namespace expressions{ namespace impl{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename VectorType>
struct span_traits<
  SpanExpr<VectorType>,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dynamic_vector_eigen<VectorType>::value
    >
  >
  : public ContainersSharedTraits<PackageIdentifier::Eigen, true, 1>
{
  static constexpr bool is_static = true;
  static constexpr bool is_dynamic = !is_static;

  using vec_remove_cv_t = typename std::remove_cv<VectorType>::type;
  using scalar_type  = typename ::pressio::Traits<vec_remove_cv_t>::scalar_type;
  using ordinal_type = typename ::pressio::Traits<vec_remove_cv_t>::ordinal_type;
  using size_type    = typename ::pressio::Traits<vec_remove_cv_t>::size_type;

  // type of the native expression
  using _native_expr_type =
    decltype(
     std::declval<VectorType>().segment(ordinal_type{}, ordinal_type{})
     );

  using _const_native_expr_type =
    decltype(
     std::declval<const VectorType>().segment(ordinal_type{}, ordinal_type{})
     );

  using native_expr_type = typename std::conditional<
    std::is_const<VectorType>::value,
    _const_native_expr_type,
    _native_expr_type
    >::type;

  using reference_type = decltype( std::declval<_native_expr_type>()(0) );
  using const_reference_type = decltype( std::declval<_const_native_expr_type>()(0) );
};
#endif


#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename VectorType>
struct span_traits<
  SpanExpr<VectorType>,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_vector_kokkos<VectorType>::value
    >
  >
  : public ContainersSharedTraits<PackageIdentifier::Kokkos, true, 1>
{
  static constexpr bool is_static = true;
  static constexpr bool is_dynamic = !is_static;

  using vec_remove_cv_t = typename std::remove_cv<VectorType>::type;
  using scalar_type	    = typename ::pressio::Traits<vec_remove_cv_t>::scalar_type;
  using execution_space = typename ::pressio::Traits<vec_remove_cv_t>::execution_space;
  using memory_space	  = typename ::pressio::Traits<vec_remove_cv_t>::memory_space;
  using device_type	    = typename ::pressio::Traits<vec_remove_cv_t>::device_type;
  using ordinal_type	  = typename ::pressio::Traits<vec_remove_cv_t>::ordinal_type;
  using reference_type  = typename ::pressio::Traits<vec_remove_cv_t>::reference_type;
  using size_type	      = typename ::pressio::Traits<vec_remove_cv_t>::size_type;
  using pair_type       = std::pair<size_type, size_type>;

  using native_expr_type =
    decltype(
      Kokkos::subview(std::declval<VectorType>(), std::declval<pair_type>())
     );

  // using _const_native_expr_type =
  //   decltype(
  //    Kokkos::subview(std::declval<const VectorType>(), std::declval<pair_type>())
  //    );
  // using native_expr_type = typename std::conditional<
  //   std::is_const<VectorType>::value,
  //   _const_native_expr_type,
  //   _native_expr_type
  //   >::type;
  // using const_data_return_type = native_expr_type const *;
  // using data_return_type = native_expr_type *;
};
#endif

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template <typename T>
// struct traits<
//   ::pressio::expressions::SpanExpr<T>,
//   ::pressio::mpl::enable_if_t<
//     ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<T>::value
//     >
//   >
//   : public ContainersSharedTraits<
//   typename details::traits<T>::wrapped_t, PackageIdentifier::Pybind,
//   true,
//   1
//   >
// {
//   static constexpr bool is_static = true;
//   static constexpr bool is_dynamic  = !is_static;
//   using wrapped_t = typename traits<T>::wrapped_t;
//   using scalar_t  = typename traits<T>::scalar_t;
//   using ordinal_t = typename traits<T>::ordinal_t;
//   using size_t    = ordinal_t;
//   using reference_t =  scalar_t &;
//   using const_reference_t = scalar_t const &;
// };
// #endif

}}}
#endif  // CONTAINERS_EXPRESSIONS_SPAN_CONTAINERS_SPAN_TRAITS_HPP_
