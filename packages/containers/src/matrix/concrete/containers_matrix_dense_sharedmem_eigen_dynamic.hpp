/*
//@HEADER
// ************************************************************************
//
// containers_matrix_dense_sharedmem_eigen_dynamic.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef CONTAINERS_MATRIX_CONCRETE_MATRIX_DENSE_SHAREDMEM_EIGEN_DYNAMIC_HPP_
#define CONTAINERS_MATRIX_CONCRETE_MATRIX_DENSE_SHAREDMEM_EIGEN_DYNAMIC_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_resizable_base.hpp"
#include "../../shared_base/containers_container_subscriptable_base.hpp"

#include "../base/containers_matrix_base.hpp"
#include "../base/containers_matrix_sharedmem_base.hpp"
#include "../base/containers_matrix_dense_sharedmem_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::meta::is_dense_dynamic_matrix_eigen<
		 wrapped_type>::value>
	     >
  : public ContainerBase< Matrix<wrapped_type>, wrapped_type >,
    public MatrixBase< Matrix<wrapped_type> >,
    public MatrixSharedMemBase< Matrix<wrapped_type> >,
    public MatrixDenseSharedMemBase< Matrix<wrapped_type> >,
    public ContainerSubscriptable2DBase<
     Matrix<wrapped_type>,
     typename details::traits<Matrix<wrapped_type>>::scalar_t,
     typename details::traits<Matrix<wrapped_type>>::ordinal_t>,
    public ContainerResizableBase<Matrix<wrapped_type>, 2>
{

  using derived_t = Matrix<wrapped_type>;
  using mytraits = typename details::traits<derived_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;
  using der_t = typename mytraits::derived_t;

public:
  Matrix() = default;

  Matrix(ord_t nrows, ord_t ncols) {
    this->resize(nrows,ncols);
  }

  explicit Matrix(const wrap_t & other)
    : data_(other){}

  explicit Matrix(const sc_t * datain)
    : data_(datain){}

  ~Matrix() = default;

public:
  sc_t & operator() (ord_t row, ord_t col){
    // check if we are withinbound
    return data_(row,col);
  }

  sc_t const & operator() (ord_t row, ord_t col) const{
    // check if we are withinbound
    return data_(row,col);
  }

  derived_t & operator+=(const derived_t & other) {
    assert( other.rows() == this->rows() );
    assert( other.cols() == this->cols() );
    this->data_ += *other.data();
    return *this;
  }

  derived_t & operator-=(const derived_t & other) {
    assert( other.rows() == this->rows() );
    assert( other.cols() == this->cols() );
    this->data_ -= *other.data();
    return *this;
  }

private:
  wrap_t * dataImpl(){
    return &data_;
  };

  wrap_t const * dataImpl() const{
    return &data_;
  };

  void addToDiagonalImpl(sc_t value) {
    // check matrix is diagonal
    assert(this->rows()==this->cols());
    for (ord_t ir=0; ir<this->rows(); ir++)
      data_(ir,ir) += value;
  };

  ord_t rowsImpl() const{
    return data_.rows();
  }

  ord_t colsImpl() const{
    return data_.cols();
  }

  void setZeroImpl() {
    data_.setConstant(static_cast<sc_t>(0));
  }

  void resizeImpl(ord_t nrows, ord_t ncols){
    data_.resize(nrows, ncols);
    this->setZero();
  }

  void scaleImpl(sc_t value) {
    for (auto i=0; i<this->rows(); i++)
      for (auto j=0; j<this->cols(); j++)
	data_(i,j) *= value;
  }

private:
  friend ContainerBase< derived_t, wrapped_type >;
  friend MatrixBase< derived_t >;
  friend MatrixSharedMemBase< derived_t >;
  friend MatrixDenseSharedMemBase< derived_t >;
  friend ContainerSubscriptable2DBase< derived_t, sc_t, ord_t>;
  friend ContainerResizableBase<derived_t, 2>;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif
