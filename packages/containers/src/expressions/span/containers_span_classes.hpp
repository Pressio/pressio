/*
//@HEADER
// ************************************************************************
//
// containers_span_classes.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_SPAN_CONTAINERS_SPAN_CLASSES_HPP_
#define CONTAINERS_EXPRESSIONS_SPAN_CONTAINERS_SPAN_CLASSES_HPP_

namespace pressio{ namespace containers{ namespace expressions{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename vector_t>
struct SpanExpr<
  vector_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_dynamic_vector_wrapper_eigen<vector_t>::value
    >
  >
{
  using this_t = SpanExpr<vector_t>;
  using traits = typename details::traits<this_t>;
  using sc_t = typename traits::scalar_t;
  using ord_t = typename traits::ordinal_t;
  using size_t = typename traits::size_t;

  using ref_t = typename traits::reference_t;
  using const_ref_t = typename traits::const_reference_t;

  using native_expr_t = typename traits::native_expr_t;
  using data_return_t = typename traits::data_return_t;
  using const_data_return_t = typename traits::const_data_return_t;

private:
  std::reference_wrapper<vector_t> vecObj_;
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

  SpanExpr(vector_t & objIn,
	   const ord_t startIndexIn,
	   const ord_t extentIn)
    : vecObj_(objIn),
      startIndex_(startIndexIn),
      extent_(extentIn),
      nativeExprObj_(vecObj_.get().data()->segment(startIndex_, extent_))
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.extent(0) );
    assert( extent_ <= objIn.extent(0) );
  }

  SpanExpr(vector_t & objIn,
	   std::pair<ord_t, ord_t> indexRange)
    : vecObj_(objIn),
      startIndex_(std::get<0>(indexRange)),
      extent_(std::get<1>(indexRange)-startIndex_),
      nativeExprObj_(vecObj_.get().data()->segment(startIndex_, extent_))
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.extent(0) );
    assert( extent_ <= objIn.extent(0) );
  }

public:
  size_t extent(size_t i) const{
    assert(i==0);
    return extent_;
  }

  const_data_return_t data() const{
    return &nativeExprObj_;
  }

  data_return_t data(){
    return &nativeExprObj_;
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
template <typename vector_t>
struct SpanExpr<
  vector_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_vector_wrapper_kokkos<vector_t>::value
    >
  >
{
  using this_t = SpanExpr<vector_t>;
  using traits = typename details::traits<this_t>;
  using sc_t = typename traits::scalar_t;
  using ord_t = typename traits::ordinal_t;
  using size_t = typename traits::size_t;
  using pair_t = typename traits::pair_t;

  using ref_t = typename traits::reference_t;
  using const_ref_t = typename traits::const_reference_t;

  using native_expr_t = typename traits::native_expr_t;
  using data_return_t = typename traits::data_return_t;
  using const_data_return_t = typename traits::const_data_return_t;

private:
  std::reference_wrapper<vector_t> vecObj_;
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

  SpanExpr(vector_t & objIn,
	   const size_t startIndexIn,
	   const size_t extentIn)
    : vecObj_(objIn),
      startIndex_(startIndexIn),
      extent_(extentIn),
      nativeExprObj_
      (Kokkos::subview(*vecObj_.get().data(),
		       std::make_pair(startIndex_, startIndex_+extent_)))
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.extent(0) );
    assert( extent_ <= objIn.extent(0) );
  }

  SpanExpr(vector_t & objIn,
	   pair_t indexRange)
    : vecObj_(objIn),
      startIndex_(std::get<0>(indexRange)),
      extent_(std::get<1>(indexRange)-startIndex_),
      nativeExprObj_
      (Kokkos::subview(*vecObj_.get().data(),
		       std::make_pair(startIndex_, startIndex_+extent_)))
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.extent(0) );
    assert( extent_ <= objIn.extent(0) );
  }

public:
  size_t extent(size_t i) const{
    assert(i==0);
    return extent_;
  }

  const_data_return_t data() const{
    return &nativeExprObj_;
  }

  data_return_t data(){
    return &nativeExprObj_;
  }

  // non-const subscripting
  /*
    I need to be careful with non-const subscripting
    because for kokkos the following would be legal:

    using kv_t	      = Kokkos::View<double *>;
    using pressio_v_t = pressio::containers::Vector<kv_t>;
    kv_t a(..);
    const pressio_v_t aw(a);
    auto s = pressio::containers::span(aw,...);
    s(0) = 1.1;

    which works because for kokkos we can assign a const view.
    but we do NOT wwant this since aw is const.
   */
  template<typename _vector_t = vector_t>
  mpl::enable_if_t<
    !std::is_const<typename std::remove_reference<_vector_t>::type>::value and
    std::is_same<typename traits::memory_space, Kokkos::HostSpace>::value,
    ref_t
    >
  operator()(size_t i)
  {
    assert(i < extent_);
    return nativeExprObj_(i);
  }

  // const subscripting
  template<typename _vector_t = vector_t>
  mpl::enable_if_t<
    std::is_same<typename traits::memory_space, Kokkos::HostSpace>::value,
    const_ref_t
    >
  operator()(size_t i) const
  {
    assert(i < extent_);
    return nativeExprObj_(i);
  }

  // template<typename _vector_t = vector_t>
  // mpl::enable_if_t<
  //   !std::is_const<typename std::remove_reference<_vector_t>::type>::value and
  //   std::is_same<typename traits::memory_space, Kokkos::HostSpace>::value,
  //   ref_t
  //   >
  // operator()(size_t i)
  // {
  //   return (*this)[i];
  // }

  // template<typename _vector_t = vector_t>
  // mpl::enable_if_t<
  //   std::is_same<typename traits::memory_space, Kokkos::HostSpace>::value,
  //   const_ref_t
  //   >
  // operator()(size_t i) const
  // {
  //   return (*this)[i];
  // }
};
#endif


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename vector_t>
struct SpanExpr<
  vector_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<vector_t>::value
    >
  >
{
  using this_t = SpanExpr<vector_t>;
  using traits = typename details::traits<this_t>;
  using sc_t = typename traits::scalar_t;
  using ord_t = typename traits::ordinal_t;
  using size_t = typename traits::size_t;
  using ref_t = typename traits::reference_t;
  using const_ref_t = typename traits::const_reference_t;
  using pair_t = std::pair<std::size_t, std::size_t>;

private:
  std::reference_wrapper<vector_t> vecObj_;
  size_t startIndex_;
  size_t extent_ = {};

public:
  SpanExpr() = delete;
  SpanExpr(const SpanExpr & other) = default;
  SpanExpr & operator=(const SpanExpr & other) = delete;
  SpanExpr(SpanExpr && other) = default;
  SpanExpr & operator=(SpanExpr && other) = delete;
  ~SpanExpr() = default;

  SpanExpr(vector_t & objIn,
	   const size_t startIndexIn,
	   const size_t extentIn)
    : vecObj_(objIn),
      startIndex_(startIndexIn),
      extent_(extentIn)
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.extent(0) );
    assert( extent_ <= objIn.extent(0) );
  }

  SpanExpr(vector_t & objIn,
	   pair_t indexRange)
    : vecObj_(objIn),
      startIndex_(std::get<0>(indexRange)),
      extent_(std::get<1>(indexRange)-startIndex_)
  {
    assert( startIndex_ >= 0 and startIndex_ < objIn.extent(0) );
    assert( extent_ <= objIn.extent(0) );
  }

public:
  size_t extent(size_t i) const{
    assert(i==0);
    return extent_;
  }

  // non-const subscripting
  template<typename _vector_t = vector_t>
  mpl::enable_if_t<
    !std::is_const<typename std::remove_reference<_vector_t>::type>::value,
    ref_t
    >
  operator()(size_t i)
  {
    assert(i < (size_t)extent_);
    return vecObj_.get()(startIndex_+i);
  }

  // const subscripting
  const_ref_t operator()(size_t i) const
  {
    assert(i < (size_t)extent_);
    return vecObj_.get()(startIndex_+i);
  }
};
#endif

}}} //end namespace pressio::containers::expressions
#endif  // CONTAINERS_EXPRESSIONS_SPAN_CONTAINERS_SPAN_CLASSES_HPP_
