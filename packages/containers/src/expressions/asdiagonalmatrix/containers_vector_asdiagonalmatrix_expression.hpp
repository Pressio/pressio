/*
//@HEADER
// ************************************************************************
//
// containers_vector_asdiagonalmatrix_expression.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_VECTOR_ASDIAGONALMATRIX_EXPRESSION_HPP_
#define CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_VECTOR_ASDIAGONALMATRIX_EXPRESSION_HPP_

namespace pressio{ namespace containers{ namespace expressions{

template <typename T>
struct AsDiagonalMatrixExpr<
  T,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_sharedmem_vector_wrapper<T>::value
    >
  >
{
  using this_t = AsDiagonalMatrixExpr<T>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using size_t = typename mytraits::size_t;
  using ref_t = typename mytraits::reference_t;
  using const_ref_t = typename mytraits::const_reference_t;

private:
  std::reference_wrapper<T> vecObj_;
  size_t extent_ = {};

public:
  AsDiagonalMatrixExpr() = delete;
  AsDiagonalMatrixExpr(const AsDiagonalMatrixExpr & other) = default;
  AsDiagonalMatrixExpr & operator=(const AsDiagonalMatrixExpr & other) = delete;
  AsDiagonalMatrixExpr(AsDiagonalMatrixExpr && other) = default;
  AsDiagonalMatrixExpr & operator=(AsDiagonalMatrixExpr && other) = delete;
  ~AsDiagonalMatrixExpr() = default;

  AsDiagonalMatrixExpr(T & objIn)
    : vecObj_(objIn),
      extent_(objIn.extent(0))
  {}

public:
  size_t extent(size_t i) const{
    assert(i==0 or i==1);
    return extent_;
  }

  const T * pressioObj() const{
    return &vecObj_.get();
  }
  T * pressioObj(){
    return &vecObj_.get();
  }

  template<class _T = T>
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_sharedmem_host_accessible_vector_wrapper<_T>::value, ref_t
    >
  operator()(size_t i, size_t j)
  {
    assert(i==j);
    assert(j < extent_);
    return vecObj_(i);
  }

  template<class _T = T>
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_sharedmem_host_accessible_vector_wrapper<_T>::value, const_ref_t
    >
  operator()(size_t i, size_t j) const
  {
    assert(i==j);
    assert(j < extent_);
    return vecObj_(i);
  }
};

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template <typename T>
struct AsDiagonalMatrixExpr<
  T,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_vector_wrapper_tpetra<T>::value or
    ::pressio::containers::predicates::is_vector_wrapper_epetra<T>::value or
    ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<T>::value
    >
  >
{
  using this_t = AsDiagonalMatrixExpr<T>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using lo_t = typename mytraits::local_ordinal_t;
  using go_t = typename mytraits::global_ordinal_t;

private:
  std::reference_wrapper<T> vecObj_;
  go_t extent_ = {};
  lo_t extentLocal_ = {};

public:
  AsDiagonalMatrixExpr() = delete;
  AsDiagonalMatrixExpr(const AsDiagonalMatrixExpr & other) = default;
  AsDiagonalMatrixExpr & operator=(const AsDiagonalMatrixExpr & other) = delete;
  AsDiagonalMatrixExpr(AsDiagonalMatrixExpr && other) = default;
  AsDiagonalMatrixExpr & operator=(AsDiagonalMatrixExpr && other) = delete;
  ~AsDiagonalMatrixExpr() = default;

  AsDiagonalMatrixExpr(T & objIn)
    : vecObj_(objIn),
      extent_(objIn.extent(0)),
      extentLocal_(objIn.extentLocal(0))
  {}

public:
  go_t extent(size_t i) const{
    assert(i==0 or i==1);
    return extent_;
  }

  lo_t extentLocal(size_t i) const{
    assert(i==0 or i==1);
    return extentLocal_;
  }

  const T * pressioObj() const{
    return &vecObj_.get();
  }
  T * pressioObj(){
    return &vecObj_.get();
  }
};
#endif

}}} //end namespace pressio::containers::expressions
#endif  // CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_VECTOR_ASDIAGONALMATRIX_EXPRESSION_HPP_