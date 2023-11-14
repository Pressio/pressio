/*
//@HEADER
// ************************************************************************
//
// span_traits.hpp
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

#ifndef EXPRESSIONS_IMPL_SPAN_TRAITS_HPP_
#define EXPRESSIONS_IMPL_SPAN_TRAITS_HPP_

namespace pressio{ namespace expressions{ namespace impl{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename VectorType>
class SpanTraits<
  SpanExpr<VectorType>,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dynamic_vector_eigen<VectorType>::value
    >
  > : public ::pressio::Traits<VectorType>
{
private:
  using _ordinal_type = typename VectorType::StorageIndex;

  using _native_expr_type =
    decltype(
     std::declval<VectorType>().segment(_ordinal_type{}, _ordinal_type{})
     );

  using _const_native_expr_type =
    decltype(
     std::declval<std::add_const_t<VectorType>>().segment(_ordinal_type{}, _ordinal_type{})
     );

public:
  using ordinal_type = _ordinal_type;

  using native_expr_type = std::conditional_t<
    std::is_const_v<VectorType>,
    _const_native_expr_type,
    _native_expr_type
    >;

  using reference_type = std::conditional_t<
    std::is_const_v<VectorType>,
    const typename VectorType::Scalar &,
    typename VectorType::Scalar &
    >;
};
#endif


#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename VectorType>
class SpanTraits<
  SpanExpr<VectorType>,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_vector_kokkos<VectorType>::value
    >
  > : public ::pressio::Traits<VectorType>
{
private:
  using _ordinal_type = typename VectorType::traits::size_type;
  using _pair_type = std::pair<_ordinal_type, _ordinal_type>;

public:
  using ordinal_type = _ordinal_type;
  using native_expr_type =
    decltype(
      Kokkos::subview(std::declval<VectorType>(), std::declval<_pair_type>())
     );

  using reference_type = typename VectorType::reference_type;
};
#endif

}}

}
#endif  // EXPRESSIONS_IMPL_SPAN_TRAITS_HPP_
