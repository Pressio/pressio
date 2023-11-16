/*
//@HEADER
// ************************************************************************
//
// subspan_classes.hpp
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

#ifndef EXPRESSIONS_IMPL_COLUMN_CLASSES_HPP_
#define EXPRESSIONS_IMPL_COLUMN_CLASSES_HPP_

namespace pressio{ namespace expressions{ namespace impl{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename MatrixType>
class ColumnExpr<
  MatrixType,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dense_matrix_eigen<MatrixType>::value
    >
  >
{

  using traits = ColumnTraits<ColumnExpr<MatrixType>>;
  using native_expr_t = typename traits::native_expr_type;
  using reference_t   = typename traits::reference_type;

private:
  MatrixType * operand_;
  std::size_t colIndex_;
  std::size_t numRows_;
  native_expr_t nativeExprObj_;

public:
  ColumnExpr(MatrixType & matObjIn, std::size_t colIndex)
    : operand_(&matObjIn),
      colIndex_(colIndex),
      numRows_(matObjIn.rows()),
      nativeExprObj_(operand_->col(colIndex))
  {
    assert( colIndex_ >= 0 and colIndex_ < std::size_t(matObjIn.cols()) );
  }

public:
  auto data() const { return operand_; }

  std::size_t extent(std::size_t i) const{
    if (i == 0) {
      return numRows_;
    } else {
      return std::size_t(1);
    }
  }

  native_expr_t const & native() const{ return nativeExprObj_; }
  native_expr_t & native(){ return nativeExprObj_; }

  reference_t operator()(std::size_t i) const{
    assert(i < numRows_);
    return (*operand_)(i, colIndex_);
  }

  reference_t operator[](std::size_t i) const{
    assert(i < numRows_);
    return (*operand_)(i, colIndex_);
  }
};
#endif

}}}
#endif  // EXPRESSIONS_IMPL_COLUMN_CLASSES_HPP_
