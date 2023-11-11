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
struct SpanExpr<
  VectorType,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dynamic_vector_eigen<VectorType>::value
    >
  >
{

  using this_t = SpanExpr<VectorType>;
  using mytraits = SpanTraits<this_t>;
  using ord_t = typename mytraits::ordinal_type;
  using size_t = typename mytraits::ordinal_type;
  using ref_t = typename mytraits::reference_type;
  using const_ref_t = typename mytraits::const_reference_type;
  using native_expr_t = typename mytraits::native_expr_type;

private:
  std::reference_wrapper<VectorType> operand_;
  std::size_t startIndex_;
  std::size_t extent_ = {};
  native_expr_t nativeExprObj_;

public:
  SpanExpr() = delete;
  SpanExpr(const SpanExpr & other) = default;
  SpanExpr & operator=(const SpanExpr & other) = delete;
  SpanExpr(SpanExpr && other) = default;
  SpanExpr & operator=(SpanExpr && other) = delete;
  ~SpanExpr() = default;

  SpanExpr(VectorType & operand,
	   std::size_t startIndex,
	   std::size_t endIndex)
    : operand_(operand),
      startIndex_(startIndex),
      extent_(endIndex-startIndex_),
      nativeExprObj_(operand_.get().segment(startIndex_, extent_))
  {
    assert( startIndex >= 0 && startIndex <= endIndex );
    assert( endIndex <= std::size_t(operand.size()) );
  }

public:
  size_t extent(size_t i) const{
    return (i == 0) ? std::size_t(extent_) : std::size_t(1);
  }

  native_expr_t const & native() const{
    return nativeExprObj_;
  }

  native_expr_t & native(){
    return nativeExprObj_;
  }

  ref_t operator()(std::size_t i)
  {
    assert(i < (std::size_t)extent_);
    return nativeExprObj_(i);
  }

  const_ref_t operator()(std::size_t i) const
  {
    assert(i < (std::size_t)extent_);
    return nativeExprObj_(i);
  }

};
#endif


#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename VectorType>
struct SpanExpr<
  VectorType,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_vector_kokkos<VectorType>::value
    >
  >
{
  using this_t = SpanExpr<VectorType>;
  using traits = SpanTraits<this_t>;
  using size_t = typename VectorType::traits::size_type;
  using pair_t = typename traits::pair_type;
  using native_expr_t = typename traits::native_expr_type;
  using ref_t = decltype( std::declval<native_expr_t>()(0) );

private:
  std::reference_wrapper<VectorType> operand_;
  std::size_t startIndex_;
  std::size_t extent_ = {};
  native_expr_t nativeExprObj_;

public:
  SpanExpr() = delete;
  SpanExpr(const SpanExpr & other) = default;
  SpanExpr & operator=(const SpanExpr & other) = delete;
  SpanExpr(SpanExpr && other) = default;
  SpanExpr & operator=(SpanExpr && other) = delete;
  ~SpanExpr() = default;

  SpanExpr(VectorType & operand,
	   std::size_t startIndex,
	   std::size_t endIndex)
    : operand_(operand),
      startIndex_(startIndex),
      extent_(endIndex-startIndex_),
      nativeExprObj_(Kokkos::subview(operand_.get(), std::make_pair(startIndex_, startIndex_+extent_)))
  {
    assert( startIndex >= 0 && startIndex <= endIndex );
    assert( endIndex <= operand.extent(0) );
  }

public:
  size_t extent(size_t i) const{
    return (i == 0) ? std::size_t(extent_) : std::size_t(1);
  }

  native_expr_t const & native() const{
    return nativeExprObj_;
  }

  native_expr_t & native(){
    return nativeExprObj_;
  }

  // I need to be careful with subscripting
  // because for kokkos the following would be legal:

  // using kv_t	      = Kokkos::View<double *>;
  // using pressio_v_t = pressio::containers::Vector<kv_t>;
  // kv_t a(..);
  // const pressio_v_t aw(a);
  // auto s = pressio::containers::span(aw,...);
  // s(0) = 1.1;

  template<typename _VectorType = VectorType>
  mpl::enable_if_t<
    std::is_same<typename _VectorType::memory_space, Kokkos::HostSpace>::value,
    ref_t
    >
  operator()(size_t i) const
  {
    assert(i < extent_);
    return nativeExprObj_(i);
  }
};
#endif

}}}
#endif  // EXPRESSIONS_IMPL_SPAN_CLASSES_HPP_
