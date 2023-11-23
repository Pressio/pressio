/*
//@HEADER
// ************************************************************************
//
// diag_classes.hpp
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

#ifndef EXPRESSIONS_IMPL_DIAGONAL_CLASSES_HPP_
#define EXPRESSIONS_IMPL_DIAGONAL_CLASSES_HPP_

namespace pressio{ namespace expressions{ namespace impl{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename MatrixType>
class DiagonalExpr<
  MatrixType,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dense_matrix_eigen<MatrixType>::value
    >
  >
{
  using traits = DiagonalTraits<DiagonalExpr<MatrixType>>;
  using reference_t = typename traits::reference_type;
  using native_expr_t = typename traits::native_expr_type;

private:
  MatrixType * operand_;
  native_expr_t nativeExprObj_;
  std::size_t numRows_ = {};
  std::size_t numCols_ = {};
  std::size_t extent_ = {};

public:
  explicit DiagonalExpr(MatrixType & matObjIn)
    : operand_(&matObjIn),
      nativeExprObj_(operand_->diagonal()),
      numRows_(operand_->rows()),
      numCols_(operand_->cols()),
      extent_(operand_->rows())
  {
    assert(numRows_ == numCols_);
  }

  std::size_t extent(std::size_t i) const{
    return (i < 1) ? extent_ : std::size_t(1);
  }
  native_expr_t const & native() const{ return nativeExprObj_; }
  native_expr_t & native(){ return nativeExprObj_; }

  reference_t operator()(std::size_t i) const{
    assert(i < extent_);
    return (*operand_)(i,i);
  }

  reference_t operator[](std::size_t i) const{
    assert(i < extent_);
    return (*operand_)(i,i);
  }

  auto data() const { return operand_; }
  // ref_t operator()(size_t i){
  //   assert(i < (size_t)extent_);
  //   return nativeExprObj_(i);
  // }
  // const_ref_t operator()(size_t i) const {
  //   assert(i < (size_t)extent_);
  //   return nativeExprObj_(i);
  // }
};
#endif


#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename MatrixType>
class DiagonalExpr<
  MatrixType,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dense_matrix_kokkos<MatrixType>::value
    >
  >
{
  static_assert(Kokkos::SpaceAccessibility<
		typename MatrixType::execution_space,
		Kokkos::HostSpace>::accessible,
		"diagonal is currently only valid for a host-accessible Kokkos View");

  using traits	= DiagonalTraits<DiagonalExpr<MatrixType>>;
  using native_expr_t	= typename traits::native_expr_type;
  using ref_t = typename traits::reference_type;

private:
  MatrixType * matObj_;
  native_expr_t nativeExprObj_;
  std::size_t extent_ = {};

  using natexpr_layout = typename native_expr_t::traits::array_layout;
  // for now leave this assert, then remove later
  static_assert
  (std::is_same<natexpr_layout, Kokkos::LayoutStride>::value,
   "The layout for the native type of the diagonal kokkos expression does not match the strided layout expected");

public:
  explicit DiagonalExpr(MatrixType & M)
    : matObj_(&M),
      nativeExprObj_(M.data(), natexpr_layout(M.extent(0), M.stride(0)+M.stride(1))),
      extent_(M.extent(0))
  {
    // make sure the diagonal is taken on a square matrix
    assert(M.extent(0) == M.extent(1));
  }

public:
  std::size_t extent(std::size_t i) const{
    return (i < 1) ? extent_ : std::size_t(1);
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
#endif  // EXPRESSIONS_IMPL_DIAGONAL_CLASSES_HPP_
