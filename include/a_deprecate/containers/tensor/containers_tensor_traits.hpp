/*
//@HEADER
// ************************************************************************
//
// containers_tensor_traits.hpp
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

#ifndef CONTAINERS_TENSOR_CONTAINERS_TENSOR_TRAITS_HPP_
#define CONTAINERS_TENSOR_CONTAINERS_TENSOR_TRAITS_HPP_

namespace pressio{ namespace containers{ namespace details{

template<typename wrapped_type>
struct traits<
  Tensor<1, wrapped_type>,
    mpl::enable_if_t<
      containers::predicates::is_array_pybind<wrapped_type>::value
    >
  >
{
  static constexpr WrappedTensorIdentifier wrapped_tensor_identifier =
    WrappedTensorIdentifier::Pybind;
  static constexpr WrappedPackageIdentifier wrapped_package_identifier =
    WrappedPackageIdentifier::Pybind;

  using wrapped_t = wrapped_type;

  static constexpr int rank = 1;
  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;
  static constexpr bool is_vector	= false;
  static constexpr bool is_matrix	= false;
  static constexpr bool is_multi_vector = false;
  static constexpr bool is_shared_mem	= true;
  static constexpr bool is_distributed	= false;

  using scalar_t  = typename wrapped_type::value_type;
  using ordinal_t = pybind11::ssize_t;
  using size_t	  = ordinal_t;
  using mut_proxy_t = decltype( std::declval<wrapped_type &>().mutable_unchecked() );
  using proxy_t	    = decltype( std::declval<const wrapped_type &>().unchecked() );
  using const_data_return_t = wrapped_type const *;
  using data_return_t	    = wrapped_type *;
  using reference_t	    = scalar_t &;
  using const_reference_t   = scalar_t const &;

  using span_ret_t   = expressions::SpanExpr<Tensor<1,wrapped_type>>;
  using span_const_ret_t = expressions::SpanExpr< const Tensor<1,wrapped_type>>;
  using asdiagonalmatrix_ret_t = expressions::AsDiagonalMatrixExpr<Tensor<1,wrapped_type>>;
  using asdiagonalmatrix_const_ret_t = expressions::AsDiagonalMatrixExpr<const Tensor<1,wrapped_type>>;
};


template<typename wrapped_type>
struct traits<
  Tensor<2, wrapped_type>,
    mpl::enable_if_t<
      containers::predicates::is_fstyle_array_pybind<wrapped_type>::value
    >
  >
{
  static constexpr WrappedTensorIdentifier wrapped_tensor_identifier =
    WrappedTensorIdentifier::Pybind;
  static constexpr WrappedPackageIdentifier wrapped_package_identifier =
    WrappedPackageIdentifier::Pybind;

  using wrapped_t = wrapped_type;

  static constexpr int rank = 2;
  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;
  static constexpr bool is_vector	= false;
  static constexpr bool is_matrix	= false;
  static constexpr bool is_multi_vector = false;
  static constexpr bool is_shared_mem	= true;
  static constexpr bool is_distributed	= false;

  using scalar_t  = typename wrapped_type::value_type;
  using ordinal_t = pybind11::ssize_t;
  using size_t	  = ordinal_t;
  using mut_proxy_t = decltype( std::declval<wrapped_type &>().mutable_unchecked() );
  using proxy_t	    = decltype( std::declval<const wrapped_type &>().unchecked() );
  using const_data_return_t = wrapped_type const *;
  using data_return_t	    = wrapped_type *;
  using reference_t	    = scalar_t &;
  using const_reference_t   = scalar_t const &;

  using diag_ret_t = expressions::DiagExpr<Tensor<2,wrapped_type>>;
  using diag_const_ret_t = expressions::DiagExpr<const Tensor<2,wrapped_type>>;
  using subspan_ret_t   = expressions::SubspanExpr<Tensor<2,wrapped_type>>;
  using subspan_const_ret_t = expressions::SubspanExpr< const Tensor<2,wrapped_type>>;

  using unsuited_layout_wrapped_t = pybind11::array_t<scalar_t, pybind11::array::c_style>;
};

template<typename wrapped_type>
struct traits<
  Tensor<3, wrapped_type>,
    mpl::enable_if_t<
      containers::predicates::is_fstyle_array_pybind<wrapped_type>::value
    >
  >
{
  static constexpr WrappedTensorIdentifier wrapped_tensor_identifier =
    WrappedTensorIdentifier::Pybind;
  static constexpr WrappedPackageIdentifier wrapped_package_identifier =
    WrappedPackageIdentifier::Pybind;

  using wrapped_t = wrapped_type;

  static constexpr int rank = 3;
  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;
  static constexpr bool is_vector	= false;
  static constexpr bool is_matrix	= false;
  static constexpr bool is_multi_vector = false;
  static constexpr bool is_shared_mem	= true;
  static constexpr bool is_distributed	= false;

  using scalar_t  = typename wrapped_type::value_type;
  using ordinal_t = pybind11::ssize_t;
  using size_t	  = ordinal_t;
  using mut_proxy_t = decltype( std::declval<wrapped_type &>().mutable_unchecked() );
  using proxy_t	    = decltype( std::declval<const wrapped_type &>().unchecked() );
  using const_data_return_t = wrapped_type const *;
  using data_return_t	    = wrapped_type *;
  using reference_t	    = scalar_t &;
  using const_reference_t   = scalar_t const &;

  using unsuited_layout_wrapped_t = pybind11::array_t<scalar_t, pybind11::array::c_style>;
};

}}} //end namespace pressio::containers::details
#endif  // CONTAINERS_TENSOR_CONTAINERS_TENSOR_TRAITS_HPP_
