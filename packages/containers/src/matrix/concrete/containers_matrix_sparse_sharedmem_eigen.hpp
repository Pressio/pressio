/*
//@HEADER
// ************************************************************************
//
// containers_matrix_sparse_sharedmem_eigen.hpp
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

#ifndef CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_SPARSE_SHAREDMEM_EIGEN_HPP_
#define CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_SPARSE_SHAREDMEM_EIGEN_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::predicates::is_sparse_matrix_eigen<
		 wrapped_type >::value>
	     >
  : public MatrixSharedMemBase< Matrix<wrapped_type>>
{

  using derived_t = Matrix<wrapped_type>;
  using mytraits = details::traits<derived_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;
  using der_t = typename mytraits::derived_t;

public:
  Matrix() = delete;

  explicit Matrix(ord_t nrows, ord_t ncols) {
    data_.resize(nrows, ncols);
    data_.makeCompressed();
  }

  // row-major matrix constructor
  template <
    typename U = ord_t,
    mpl::enable_if_t<
     !std::is_void<U>::value and 
     mytraits::is_row_major==1, int> = 0
    >
  explicit Matrix(U nrows, U ncols, U nonZerosPerRow) {
    data_.resize(nrows, ncols);
    if( nonZerosPerRow > ncols )
      throw std::runtime_error(
    "SPARSE MATRIX CNTR: estimated nonzeros larger then size of cols");
    data_.reserve(Eigen::VectorXi::Constant(nrows,nonZerosPerRow));
    data_.makeCompressed();
  }

  // col-major matrix constructor
  template <
    typename U = ord_t,
    mpl::enable_if_t<
     !std::is_void<U>::value and
      mytraits::is_row_major==0, int> = 0
    >
  explicit Matrix(U nrows, U ncols, U nonZerosPerCol) {
    data_.resize(nrows, ncols);
    if( nonZerosPerCol > nrows )
      throw std::runtime_error(
    "SPARSE MATRIX CNTR: estimated nonzeros larger then size of rows");
    data_.reserve(Eigen::VectorXi::Constant(ncols,nonZerosPerCol));
    data_.makeCompressed();
  }

  explicit Matrix(const wrap_t & other) : data_(other){
    data_.makeCompressed();
  }

  ~Matrix() = default;

  // note here that we return by copy
  sc_t operator() (ord_t row, ord_t col) const{
    // eigen returns 0 if the item is zero
    return data_.coeff(row,col);
  }

  derived_t & operator+=(const derived_t & other) {
    assert(haveCompatibleDimensions(*this, other) );
    this->data_ += *other.data();
    return *this;
  }

  derived_t & operator-=(const derived_t & other) {
    assert(haveCompatibleDimensions(*this, other) );
    this->data_ -= *other.data();
    return *this;
  }

  wrap_t const * data() const{
    return &data_;
  };

  wrap_t * data(){
    return &data_;
  };

private:
  friend MatrixSharedMemBase<derived_t>;
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif  // CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_SPARSE_SHAREDMEM_EIGEN_HPP_
