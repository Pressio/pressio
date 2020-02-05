/*
//@HEADER
// ************************************************************************
//
// containers_vector_span_expression.hpp
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

#ifndef CONTAINERS_VECTOR_SPAN_EXPRESSION_HPP_
#define CONTAINERS_VECTOR_SPAN_EXPRESSION_HPP_

#include "../vector/containers_vector_meta.hpp"
#include "containers_expression_base.hpp"

namespace pressio{ namespace containers{ namespace expressions{

template <typename vector_t, typename scalar_type>
struct SpanExpr<
  vector_t, scalar_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_dynamic_vector_wrapper_eigen<vector_t>::value
    >
  >
  : public BaseExpr< SpanExpr<vector_t, scalar_type> >
{
  using native_t   = typename ::pressio::containers::details::traits<vector_t>::wrapped_t;

private:
  vector_t & vecObj_;
  std::size_t startIndex_;
  std::size_t extent_ = {};

public:
  SpanExpr() = delete;
  ~SpanExpr() = default;
  SpanExpr(const SpanExpr & other) = default;
  SpanExpr(SpanExpr && other) = default;
  SpanExpr & operator=(const SpanExpr & other) = default;
  SpanExpr & operator=(SpanExpr && other) = default;

  SpanExpr(vector_t & objIn,
	   const std::size_t startIndexIn,
	   const std::size_t extentIn)
    : vecObj_(objIn), startIndex_(startIndexIn), extent_(extentIn)
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.size() );
    assert( extent_ <= objIn.size() );
  }

  SpanExpr(vector_t & objIn,
	   std::pair<std::size_t, std::size_t> indexRange)
    : vecObj_(objIn),
      startIndex_(std::get<0>(indexRange)),
      extent_(std::get<1>(indexRange)-startIndex_)
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.size() );
    assert( extent_ <= objIn.size() );
  }

  std::size_t const & size() const{ return extent_; }

  scalar_type & operator[](std::size_t i)
  {
    assert(i < extent_);
    return vecObj_(startIndex_+i);
  }

  scalar_type const & operator[](std::size_t i) const
  {
    assert(i < extent_);
    return vecObj_(startIndex_+i);
  }

  auto operator()()
    -> decltype(vecObj_.data()->segment(startIndex_, extent_))
  {
    return vecObj_.data()->segment(startIndex_, extent_);
  }

  auto operator()() const
    -> decltype( std::declval<const native_t>().segment(startIndex_, extent_))
  {
    return static_cast<native_t const *>(vecObj_.data())->segment(startIndex_, extent_);
  }
};




#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename vector_t, typename scalar_type>
struct SpanExpr<
  vector_t, scalar_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<vector_t>::value
    >
  >
  : public BaseExpr< SpanExpr<vector_t, scalar_type> >
{
  using native_t   = typename ::pressio::containers::details::traits<vector_t>::wrapped_t;

private:
  vector_t & vecObj_;
  std::size_t startIndex_;
  std::size_t extent_ = {};
  using pair_t = std::pair<std::size_t, std::size_t>;

public:
  SpanExpr() = delete;
  ~SpanExpr() = default;
  SpanExpr(const SpanExpr & other) = default;
  SpanExpr(SpanExpr && other) = default;
  SpanExpr & operator=(const SpanExpr & other) = default;
  SpanExpr & operator=(SpanExpr && other) = default;

  SpanExpr(vector_t & objIn,
	   const std::size_t startIndexIn,
	   const std::size_t extentIn)
    : vecObj_(objIn), startIndex_(startIndexIn), extent_(extentIn)
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.size() );
    assert( extent_ <= objIn.size() );
  }

  SpanExpr(vector_t & objIn,
	   std::pair<std::size_t, std::size_t> indexRange)
    : vecObj_(objIn),
      startIndex_(std::get<0>(indexRange)),
      extent_(std::get<1>(indexRange)-startIndex_)
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.size() );
    assert( extent_ <= objIn.size() );
  }

  std::size_t const & size() const{ return extent_; }

  // TODO: enable only on host
  scalar_type & operator[](std::size_t i)
  {
    assert(i < extent_);
    return vecObj_(startIndex_+i);
  }

  // TODO: enable only on host
  scalar_type const & operator[](std::size_t i) const
  {
    assert(i < extent_);
    return vecObj_(startIndex_+i);
  }

  auto operator()()
    -> decltype
    (
     Kokkos::subview(*vecObj_.data(), std::declval<pair_t>())
     )
  {
    return Kokkos::subview(*vecObj_.data(),
			   std::make_pair(startIndex_, startIndex_+extent_));
  }

  auto operator()() const
    -> decltype
    (
     Kokkos::subview(*vecObj_.data(), std::declval<pair_t>())
     )
  {
    return Kokkos::subview(*vecObj_.data(),
			   std::make_pair(startIndex_, startIndex_+extent_));
  }

  // auto operator()() const
  //   -> decltype( std::declval<const native_t>().segment(startIndex_, extent_))
  // {
  //   return static_cast<native_t const *>(vecObj_.data())->segment(startIndex_, extent_);
  // }
};
#endif


}}} //end namespace pressio::containers::expressions
#endif
