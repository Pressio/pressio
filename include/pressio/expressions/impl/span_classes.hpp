/*
//@HEADER
// ************************************************************************
//
// span_classes.hpp
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

#ifndef EXPRESSIONS_IMPL_SPAN_CLASSES_HPP_
#define EXPRESSIONS_IMPL_SPAN_CLASSES_HPP_

namespace pressio{ namespace expressions{ namespace impl{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename VectorType>
class SpanExpr<
  VectorType,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dynamic_vector_eigen<VectorType>::value
    >
  >
{

  using mytraits      = SpanTraits<SpanExpr<VectorType>>;
  using reference_t   = typename mytraits::reference_type;
  using native_expr_t = typename mytraits::native_expr_type;

  VectorType * operand_ = nullptr;
  std::size_t startIndex_ = {};
  std::size_t extent_ = {};
  native_expr_t nativeExprObj_;

public:
  SpanExpr(VectorType & operand,
	   std::size_t startIndex,
	   std::size_t endIndex)
    : operand_(&operand),
      startIndex_(startIndex),
      extent_(endIndex-startIndex_),
      nativeExprObj_(operand_->segment(startIndex_, extent_))
  {
    assert( startIndex >= 0 && startIndex <= endIndex );
    assert( endIndex <= std::size_t(operand.size()) );
  }

  std::size_t extent(std::size_t i) const{
    return (i == 0) ? extent_ : std::size_t(1);
  }

  native_expr_t const & native() const{ return nativeExprObj_; }
  native_expr_t & native(){ return nativeExprObj_; }

  reference_t operator()(std::size_t i) const{
    assert(i < extent_);
    return (*operand_)(startIndex_ + i);
  }

  reference_t operator[](std::size_t i) const{
    assert(i < extent_);
    return (*operand_)(startIndex_ + i);
  }

  auto data() const { return operand_; }
};
#endif


#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename VectorType>
class SpanExpr<
  VectorType,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_vector_kokkos<VectorType>::value
    >
  >
{
  static_assert(Kokkos::SpaceAccessibility<
		typename VectorType::execution_space,
		Kokkos::HostSpace>::accessible,
		"span is currently only valid for a host-accessible Kokkos View");

  using traits = SpanTraits<SpanExpr<VectorType>>;
  using native_expr_t = typename traits::native_expr_type;
  using ref_t = typename traits::reference_type; //decltype( std::declval<native_expr_t>()(0) );

  VectorType * operand_;
  std::size_t startIndex_ = {};
  std::size_t extent_ = {};
  native_expr_t nativeExprObj_;

public:
  SpanExpr(VectorType & operand,
	   std::size_t startIndex,
	   std::size_t endIndex)
    : operand_(&operand),
      startIndex_(startIndex),
      extent_(endIndex-startIndex_),
      nativeExprObj_(Kokkos::subview(*operand_, std::make_pair(startIndex_, startIndex_+extent_)))
  {
    assert( startIndex >= 0 && startIndex <= endIndex );
    assert( endIndex <= operand.extent(0) );
  }

  auto data() const { return operand_; }

  std::size_t extent(std::size_t i) const{
    return (i == 0) ? extent_ : std::size_t(1);
  }

  native_expr_t const & native() const{ return nativeExprObj_; }
  native_expr_t & native(){ return nativeExprObj_; }

  ref_t operator()(std::size_t i) const{
    assert(i < extent_);
    return nativeExprObj_(i);
  }
  ref_t operator[](std::size_t i) const{
    assert(i < extent_);
    return nativeExprObj_(i);
  }
};
#endif

}}}
#endif  // EXPRESSIONS_IMPL_SPAN_CLASSES_HPP_
