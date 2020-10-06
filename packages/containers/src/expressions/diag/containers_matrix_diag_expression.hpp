/*
//@HEADER
// ************************************************************************
//
// containers_matrix_diag_expression.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_DIAG_CONTAINERS_MATRIX_DIAG_EXPRESSION_HPP_
#define CONTAINERS_EXPRESSIONS_DIAG_CONTAINERS_MATRIX_DIAG_EXPRESSION_HPP_

namespace pressio{ namespace containers{ namespace expressions{

template <typename matrix_t>
struct DiagExpr<
  matrix_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_dense_matrix_wrapper_eigen<matrix_t>::value
    >
  >
{
  using this_t = DiagExpr<matrix_t>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using size_t = typename mytraits::size_t;

  using ref_t = typename mytraits::reference_t;
  using const_ref_t = typename mytraits::const_reference_t;

  using native_expr_t = typename mytraits::native_expr_t;
  using data_return_t = typename mytraits::data_return_t;
  using const_data_return_t = typename mytraits::const_data_return_t;

  using pair_t = std::pair<std::size_t, std::size_t>;

private:
  std::reference_wrapper<matrix_t> matObj_;
  native_expr_t nativeExprObj_;
  size_t numRows_ = {};
  size_t numCols_ = {};
  size_t extent_ = {};

public:
  DiagExpr() = delete;

  DiagExpr(const DiagExpr & other) = default;
  DiagExpr & operator=(const DiagExpr & other) = delete;

  DiagExpr(DiagExpr && other) = default;
  DiagExpr & operator=(DiagExpr && other) = delete;
  ~DiagExpr() = default;

  DiagExpr(matrix_t & matObjIn)
    : matObj_(matObjIn),
    nativeExprObj_(matObj_.get().data()->diagonal()),
    numRows_(matObj_.get().data()->rows()),
    numCols_(matObj_.get().data()->cols()),
    extent_(matObj_.get().data()->rows())
  {
    assert(numRows_ == numCols_);
  }

  size_t extent() const{
    return extent_;
  }

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

  ref_t operator()(size_t i)
  {
    assert(i < (size_t)extent_);
    return nativeExprObj_(i);
  }

  const_ref_t operator()(size_t i) const
  {
    assert(i < (size_t)extent_);
    return nativeExprObj_(i);
  }

  ref_t operator[](size_t i)
  {
    assert(i < (size_t)extent_);
    return nativeExprObj_(i);
  }

  const_ref_t operator[](size_t i) const
  {
    assert(i < (size_t)extent_);
    return nativeExprObj_(i);
  }
};


#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename matrix_t>
struct DiagExpr<
  matrix_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_dense_matrix_wrapper_kokkos<matrix_t>::value
    >
  >
{
  using this_t = DiagExpr<matrix_t>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using size_t = typename mytraits::size_t;

  using ref_t = typename mytraits::reference_t;
  using const_ref_t = typename mytraits::const_reference_t;

  //using native_expr_t = typename mytraits::native_expr_t;
  ///using data_return_t = typename mytraits::data_return_t;
  //using const_data_return_t = typename mytraits::const_data_return_t;

  using pair_t = std::pair<std::size_t, std::size_t>;

private:
  std::reference_wrapper<matrix_t> matObj_;
  //native_expr_t nativeExprObj_;
  size_t numRows_ = {};
  size_t numCols_ = {};
  size_t extent_ = {};

public:
  DiagExpr() = delete;

  DiagExpr(const DiagExpr & other) = default;
  DiagExpr & operator=(const DiagExpr & other) = delete;

  DiagExpr(DiagExpr && other) = default;
  DiagExpr & operator=(DiagExpr && other) = delete;

  ~DiagExpr() = default;

  DiagExpr(matrix_t & matObjIn)
    : matObj_(matObjIn),
      /*nativeExprObj_(),*/
      numRows_(matObj_.get().extent(0)),
      numCols_(matObj_.get().extent(1)),
      extent_(matObj_.get().extent(0))
  {
    assert(numRows_ == numCols_);
  }

  size_t extent() const{
    return extent_;
  }

  size_t extent(size_t i) const{
    assert(i==0);
    return extent_;
  }

  matrix_t & getUnderlyingObject(){
    return matObj_;
  }

  const matrix_t & getUnderlyingObject() const{
    return matObj_;
  }

  // const_data_return_t data() const{
  //   return &nativeExprObj_;
  // }
  // data_return_t data(){
  //   return &nativeExprObj_;
  // }
  ref_t operator[](size_t i)
  {
    assert(i < (size_t)extent_);
    return matObj_(i,i);
  }

  const_ref_t operator[](size_t i) const
  {
    assert(i < (size_t)extent_);
    return matObj_(i,i);
  }

  ref_t operator()(size_t i)
  {
    assert(i < (size_t)extent_);
    return matObj_(i,i);
  }

  const_ref_t operator()(size_t i) const
  {
    assert(i < (size_t)extent_);
    return matObj_(i,i);
  }
};
#endif

}}} //end namespace pressio::containers::expressions
#endif  // CONTAINERS_EXPRESSIONS_DIAG_CONTAINERS_MATRIX_DIAG_EXPRESSION_HPP_
