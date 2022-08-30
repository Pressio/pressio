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
  using size_t = typename mytraits::size_type;
  using ref_t = typename mytraits::reference_type;
  using const_ref_t = typename mytraits::const_reference_type;
  using native_expr_t = typename mytraits::native_expr_type;

private:
  std::reference_wrapper<VectorType> vecObj_;
  ord_t startIndex_;
  ord_t extent_ = {};
  native_expr_t nativeExprObj_;

public:
  SpanExpr() = delete;
  SpanExpr(const SpanExpr & other) = default;
  SpanExpr & operator=(const SpanExpr & other) = delete;
  SpanExpr(SpanExpr && other) = default;
  SpanExpr & operator=(SpanExpr && other) = delete;
  ~SpanExpr() = default;

  SpanExpr(VectorType & objIn,
	   const ord_t startIndexIn,
	   const ord_t extentIn)
    : vecObj_(objIn),
      startIndex_(startIndexIn),
      extent_(extentIn),
      nativeExprObj_(vecObj_.get().segment(startIndex_, extent_))
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.size() );
    assert( extent_ <= objIn.size() );
  }

  SpanExpr(VectorType & objIn,
	   std::pair<ord_t, ord_t> indexRange)
    : vecObj_(objIn),
      startIndex_(std::get<0>(indexRange)),
      extent_(std::get<1>(indexRange)-startIndex_),
      nativeExprObj_(vecObj_.get().segment(startIndex_, extent_))
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.size() );
    assert( extent_ <= objIn.size() );
  }

public:
  size_t extent(size_t i) const{
    assert(i==0);
    (void) i;
    return extent_;
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
  using ord_t = typename traits::ordinal_type;
  using size_t = typename traits::size_type;
  using pair_t = typename traits::pair_type;
  using ref_t = typename traits::reference_type;
  using native_expr_t = typename traits::native_expr_type;

private:
  std::reference_wrapper<VectorType> vecObj_;
  size_t startIndex_;
  size_t extent_ = {};
  native_expr_t nativeExprObj_;

public:
  SpanExpr() = delete;
  SpanExpr(const SpanExpr & other) = default;
  SpanExpr & operator=(const SpanExpr & other) = delete;
  SpanExpr(SpanExpr && other) = default;
  SpanExpr & operator=(SpanExpr && other) = delete;
  ~SpanExpr() = default;

  SpanExpr(VectorType & objIn,
	   const size_t startIndexIn,
	   const size_t extentIn)
    : vecObj_(objIn),
      startIndex_(startIndexIn),
      extent_(extentIn),
      nativeExprObj_(Kokkos::subview(vecObj_.get(),std::make_pair(startIndex_, startIndex_+extent_)))
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.extent(0) );
    assert( extent_ <= objIn.extent(0) );
  }

  SpanExpr(VectorType & objIn,
	   pair_t indexRange)
    : vecObj_(objIn),
      startIndex_(std::get<0>(indexRange)),
      extent_(std::get<1>(indexRange)-startIndex_),
      nativeExprObj_(Kokkos::subview(vecObj_.get(), std::make_pair(startIndex_, startIndex_+extent_)))
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.extent(0) );
    assert( extent_ <= objIn.extent(0) );
  }

public:
  size_t extent(size_t i) const{
    assert(i==0);
    (void) i;
    return extent_;
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
    std::is_same<typename traits::memory_space, Kokkos::HostSpace>::value,
    ref_t
    >
  operator()(size_t i) const
  {
    assert(i < extent_);
    return nativeExprObj_(i);
  }
};
#endif


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename VectorType>
struct SpanExpr<
  VectorType,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_array_pybind<VectorType>::value
    >
  >
{
  using this_t = SpanExpr<VectorType>;
  using traits = SpanTraits<this_t>;
  using size_t = typename traits::size_type;
  using ref_t = typename traits::reference_type;
  using const_ref_t = typename traits::const_reference_type;
  using pair_t = std::pair<std::size_t, std::size_t>;

private:
  std::reference_wrapper<VectorType> vecObj_;
  size_t startIndex_;
  size_t extent_ = {};

public:
  SpanExpr() = delete;
  SpanExpr(const SpanExpr & other) = default;
  SpanExpr & operator=(const SpanExpr & other) = delete;
  SpanExpr(SpanExpr && other) = default;
  SpanExpr & operator=(SpanExpr && other) = delete;
  ~SpanExpr() = default;

  SpanExpr(VectorType & objIn,
	   const size_t startIndexIn,
	   const size_t extentIn)
    : vecObj_(objIn),
      startIndex_(startIndexIn),
      extent_(extentIn)
  {
    assert(objIn.ndim()==1);
    assert(startIndex_ >= 0 and startIndex_ < objIn.shape(0));
    assert(extent_ <= objIn.shape(0));
  }

  SpanExpr(VectorType & objIn,
	   pair_t indexRange)
    : vecObj_(objIn),
      startIndex_(std::get<0>(indexRange)),
      extent_(std::get<1>(indexRange)-startIndex_)
  {
    assert(objIn.ndim()==1);
    assert(startIndex_ >= 0 and startIndex_ < objIn.shape(0));
    assert(extent_ <= objIn.shape(0));
  }

public:
  int ndim() const{
    return 1;
  }

  size_t extent(size_t i) const{
    assert(i==0);
    return extent_;
  }

  // non-const subscripting
  ref_t operator()(size_t i)
  {
    assert(i < (size_t)extent_);
    return vecObj_(startIndex_+i);
  }

  // const subscripting
  const_ref_t operator()(size_t i) const
  {
    assert(i < (size_t)extent_);
    return vecObj_(startIndex_+i);
  }
};
#endif

}}}
#endif  // EXPRESSIONS_IMPL_SPAN_CLASSES_HPP_
