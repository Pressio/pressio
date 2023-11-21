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
template <typename T>
class ColumnExpr<
  T,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dense_matrix_eigen<T>::value
    >
  >
{

  using traits = ColumnTraits<ColumnExpr<T>>;
  using native_expr_t = typename traits::native_expr_type;
  using reference_t   = typename traits::reference_type;

  T * operand_;
  std::size_t colIndex_;
  std::size_t numRows_;
  native_expr_t nativeExprObj_;

public:
  ColumnExpr(T & matObjIn, std::size_t colIndex)
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


#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template <typename T>
class ColumnExpr<
  T, ::pressio::mpl::enable_if_t<
       ::pressio::is_multi_vector_tpetra<T>::value > >
{

  using traits = ColumnTraits<ColumnExpr<T>>;
  using sc_t   = typename T::scalar_type;
  using go_t   = typename T::global_ordinal_type;
  using lo_t   = typename T::local_ordinal_type;
  using node_t = typename T::node_type;
  using vec_t  = Tpetra::Vector<sc_t, lo_t, go_t, node_t>;

  T * operand_;
  std::size_t colIndex_;
  std::size_t numRowsLocal_;
  std::size_t numRowsGlobal_;
  vec_t colVec_;

public:
  explicit ColumnExpr(T & matObjIn, std::size_t colIndex)
    : operand_(&matObjIn)
    , colIndex_(colIndex)
    , numRowsLocal_(matObjIn.getLocalLength())
    , numRowsGlobal_(matObjIn.getGlobalLength())
    , colVec_(*matObjIn.getVectorNonConst(colIndex))
  {
    assert( colIndex_ >= 0 &&
	    colIndex_ < std::size_t(matObjIn.getNumVectors()) );
  }

public:
  auto data() const { return operand_; }

  std::size_t extentLocal(std::size_t i) const{
    if (i == 0) { return numRowsLocal_; } else { return std::size_t(1); }
  }

  std::size_t extentGlobal(std::size_t i) const{
    if (i == 0) { return numRowsGlobal_; } else { return std::size_t(1); }
  }

  vec_t native() const{ return colVec_; }
  vec_t native(){ return colVec_; }
};

template <typename T>
class ColumnExpr<
  T, ::pressio::mpl::enable_if_t<
       ::pressio::is_multi_vector_tpetra_block<T>::value > >
{
  using traits = ColumnTraits<ColumnExpr<T>>;

  T * operand_;
  std::size_t colIndex_;
  std::size_t numRowsGlobal_;
  typename T::mv_type tpetraMv_;
  typename traits::native_expr_type colVec_;

public:
  explicit ColumnExpr(T & matObjIn, std::size_t colIndex)
    : operand_(&matObjIn)
    , colIndex_(colIndex)
    , numRowsGlobal_(matObjIn.getMap()->getGlobalNumElements())
    , tpetraMv_(*matObjIn.getMultiVectorView().getVectorNonConst(colIndex))
    , colVec_(tpetraMv_, *(matObjIn.getMap()), matObjIn.getBlockSize())
  {
    assert( colIndex_ >= 0 &&
	    colIndex_ < std::size_t(matObjIn.getNumVectors()) );
  }

public:
  auto data() const { return operand_; }

  std::size_t extentGlobal(std::size_t i) const{
    if (i == 0) { return numRowsGlobal_; } else { return std::size_t(1); }
  }

  typename traits::native_expr_type native() const{ return colVec_; }
  typename traits::native_expr_type native(){ return colVec_; }
};

#endif

}}}
#endif  // EXPRESSIONS_IMPL_COLUMN_CLASSES_HPP_
