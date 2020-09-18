/*
//@HEADER
// ************************************************************************
//
// containers_matrix_dense_sharedmem_eigen_dynamic.hpp
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

#ifndef CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_DENSE_SHAREDMEM_EIGEN_DYNAMIC_HPP_
#define CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_DENSE_SHAREDMEM_EIGEN_DYNAMIC_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class DenseMatrix<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::predicates::is_dense_dynamic_matrix_eigen<
		 wrapped_type>::value>
	     >
  : public DenseMatrixSharedMemBase< DenseMatrix<wrapped_type> >
{

  using derived_t = DenseMatrix<wrapped_type>;
  using mytraits = typename details::traits<derived_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;
  using der_t = typename mytraits::derived_t;

public:
  DenseMatrix() = default;

  DenseMatrix(ord_t nrows, ord_t ncols)
    : data_(nrows, ncols){
    data_.setConstant(static_cast<sc_t>(0));
  }

  explicit DenseMatrix(const wrap_t & other)
    : data_(other){}

  explicit DenseMatrix(const sc_t * datain)
    : data_(datain){}

  ~DenseMatrix() = default;

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
    assert( other.extent(0) == this->extent(0) );
    assert( other.extent(1) == this->extent(1) );
    this->data_ += *other.data();
    return *this;
  }

  derived_t & operator-=(const derived_t & other) {
    assert( other.extent(0) == this->extent(0) );
    assert( other.extent(1) == this->extent(1) );
    this->data_ -= *other.data();
    return *this;
  }

public:
  wrap_t * data(){
    return &data_;
  };

  wrap_t const * data() const{
    return &data_;
  };

  ord_t extent(ord_t i) const {
    assert(i==0 or i==1);
    return (i==0) ? data_.rows() : data_.cols();
  }

private:
  friend DenseMatrixSharedMemBase< derived_t >;
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif  // CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_DENSE_SHAREDMEM_EIGEN_DYNAMIC_HPP_
