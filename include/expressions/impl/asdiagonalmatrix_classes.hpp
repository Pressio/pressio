/*
//@HEADER
// ************************************************************************
//
// containers_asdiagonalmatrix_classes.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_ASDIAGONALMATRIX_CLASSES_HPP_
#define CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_ASDIAGONALMATRIX_CLASSES_HPP_

namespace pressio{ namespace expressions{ namespace impl{

#if defined PRESSIO_ENABLE_TPL_EIGEN or defined PRESSIO_ENABLE_TPL_PYBIND11
template <typename VectorType>
struct AsDiagonalMatrixExpr<
  VectorType,
  ::pressio::mpl::enable_if_t<
#ifdef PRESSIO_ENABLE_TPL_EIGEN
    ::pressio::is_dynamic_vector_eigen<
      typename std::remove_cv<VectorType>::type
     >::value
#endif
#if defined PRESSIO_ENABLE_TPL_EIGEN and defined PRESSIO_ENABLE_TPL_PYBIND11
    or
#endif
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    ::pressio::is_rank1_tensor_pybind<
     typename std::remove_cv<VectorType>::type
     >::value
#endif
    >
  >
{
  using this_t = AsDiagonalMatrixExpr<VectorType>;
  using mytraits = asdiagmatrix_traits<this_t>;
  using sc_t = typename mytraits::scalar_type;
  using size_t = typename mytraits::size_type;
  using ref_t = typename mytraits::reference_type;
  using const_ref_t = typename mytraits::const_reference_type;

private:
  std::reference_wrapper<VectorType> vecObj_;
  size_t extent_ = {};

public:
  AsDiagonalMatrixExpr() = delete;
  AsDiagonalMatrixExpr(const AsDiagonalMatrixExpr & other) = default;
  AsDiagonalMatrixExpr & operator=(const AsDiagonalMatrixExpr & other) = delete;
  AsDiagonalMatrixExpr(AsDiagonalMatrixExpr && other) = default;
  AsDiagonalMatrixExpr & operator=(AsDiagonalMatrixExpr && other) = delete;
  ~AsDiagonalMatrixExpr() = default;

  AsDiagonalMatrixExpr(VectorType & objIn)
    : vecObj_(objIn),
      extent_(objIn.size())
  {}

public:
  size_t extent(size_t i) const{
    assert(i==0 or i==1);
    return extent_;
  }

  // const VectorType * pressioObj() const{
  //   return &vecObj_.get();
  // }
  // VectorType * pressioObj(){
  //   return &vecObj_.get();
  // }

  ref_t operator()(size_t i, size_t j)
  {
    assert(i==j and j < extent_);
    return vecObj_(i);
  }

  const_ref_t operator()(size_t i, size_t j) const
  {
    assert(i==j and j < extent_);
    return vecObj_(i);
  }
};
#endif

// #ifdef PRESSIO_ENABLE_TPL_TRILINOS
// template <typename T>
// struct AsDiagonalMatrixExpr<
//   T,
//   ::pressio::mpl::enable_if_t<
//     ::pressio::is_vector_tpetra<typename std::remove_cv<VectorType>::type>::value or
//     ::pressio::is_vector_epetra<typename std::remove_cv<VectorType>::type>::value or
//     ::pressio::is_vector_tpetra_block<typename std::remove_cv<VectorType>::type>::value
//     >
//   >
// {
//   using this_t = AsDiagonalMatrixExpr<VectorType>;
//   using traits = typename details::traits<this_t>;
//   using sc_t = typename traits::scalar_t;
//   using lo_t = typename traits::local_ordinal_t;
//   using go_t = typename traits::global_ordinal_t;

// private:
//   std::reference_wrapper<VectorType> vecObj_;
//   go_t extent_ = {};
//   lo_t extentLocal_ = {};

// public:
//   AsDiagonalMatrixExpr() = delete;
//   AsDiagonalMatrixExpr(const AsDiagonalMatrixExpr & other) = default;
//   AsDiagonalMatrixExpr & operator=(const AsDiagonalMatrixExpr & other) = delete;
//   AsDiagonalMatrixExpr(AsDiagonalMatrixExpr && other) = default;
//   AsDiagonalMatrixExpr & operator=(AsDiagonalMatrixExpr && other) = delete;
//   ~AsDiagonalMatrixExpr() = default;

//   template <typename _T = T, 
//     typename mpl::enable_if_t<
//       ::pressio::is_vector_epetra<_T>::value
//     > * = nullptr
//     >
//   AsDiagonalMatrixExpr(_T & objIn)
//     : vecObj_(objIn),
//       extent_(objIn.GlobalLength()),
//       extentLocal_(objIn.MyLength())
//   {}

//   template <typename _T = T, 
//     typename mpl::enable_if_t<
//       ::pressio::is_vector_tpetra<_T>::value
//     > * = nullptr
//     >
//   AsDiagonalMatrixExpr(_T & objIn)
//     : vecObj_(objIn),
//       extent_(objIn.getGlobalLength()),
//       extentLocal_(objIn.getLocalLength())
//   {}

//   template <typename _T = T, 
//     typename mpl::enable_if_t<
//       ::pressio::is_vector_tpetra_block<_T>::value
//     > * = nullptr
//     >
//   AsDiagonalMatrixExpr(_T & objIn)
//     : vecObj_(objIn),
//       extent_(objIn.getMap()->getGlobalNumElements()),
//       extentLocal_(objIn.getMap()->getNodeNumElements())
//   {}

// public:
//   go_t extent(size_t i) const{
//     assert(i==0 or i==1);
//     return extent_;
//   }

//   lo_t extentLocal(size_t i) const{
//     assert(i==0 or i==1);
//     return extentLocal_;
//   }

//   const VectorType * pressioObj() const{
//     return &vecObj_.get();
//   }
//   VectorType * pressioObj(){
//     return &vecObj_.get();
//   }
// };
// #endif

}}}
#endif  // CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_ASDIAGONALMATRIX_CLASSES_HPP_
