/*
//@HEADER
// ************************************************************************
//
// containers_matrix_sharedmem_pybind11.hpp
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

#ifndef CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_SHAREDMEM_PYBIND11_HPP_
#define CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_SHAREDMEM_PYBIND11_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class DenseMatrix<
  wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_array_pybind<wrapped_type>::value
    >
  >
  : public DenseMatrixSharedMemBase< DenseMatrix<wrapped_type> >
{

  using this_t	    = DenseMatrix<wrapped_type>;
  using mytraits    = typename details::traits<this_t>;
  using sc_t	    = typename mytraits::scalar_t;
  using ord_t	    = typename mytraits::ordinal_t;
  using wrap_t	    = typename mytraits::wrapped_t;
  // using ref_t	    = typename mytraits::reference_t;
  // using const_ref_t = typename mytraits::const_reference_t;
  // using mut_proxy_t = typename mytraits::mut_proxy_t;
  // using proxy_t	    = typename mytraits::proxy_t;

public:
  DenseMatrix() = delete;

  explicit DenseMatrix(std::size_t ext1, std::size_t ext2)
    : data_({ext1, ext2}){}

  explicit DenseMatrix(const wrap_t & src)
    : data_{ wrap_t(const_cast<wrap_t &>(src).request()) }
  {
    // src must be a matrix to be wraped into a matrix
    assert( data_.ndim() == 2 );

    // copy data from src to this
    auto proxy = data_.mutable_unchecked();
    const auto srcPx = src.unchecked();
    for (ord_t i=0; i<src.shape(0); ++i){
      for (ord_t j=0; j<src.shape(1); ++j){
	proxy(i,j) = srcPx(i,j);
      }
    }
  }

  // use only if you know what you are doing
  // it is currently used only in specific places
  DenseMatrix(wrap_t src, ::pressio::view)
    : data_{src}{
    assert( data_.ndim() == 2 );
  }

  // copy cnstr
  DenseMatrix(DenseMatrix const & other)
    : data_({ other.extent(0), other.extent(1) })
  {
    assert( other.data_.ndim() == 2 );
    // copy data from src to this
    auto proxy = data_.mutable_unchecked();
    const auto srcPx = other.data_.unchecked();
    for (ord_t i=0; i<other.extent(0); ++i){
      for (ord_t j=0; j<other.extent(1); ++j){
	proxy(i,j) = srcPx(i,j);
      }
    }
  }

  // copy assignment
  DenseMatrix & operator=(const DenseMatrix & other){
    if (&other != this){
      assert(this->extent(0) == other.extent(0));
      assert(this->extent(1) == other.extent(1));

      // copy data from src to this
      auto proxy = data_.mutable_unchecked();
      const auto srcPx = other.data_.unchecked();
      for (ord_t i=0; i<other.extent(0); ++i){
  	for (ord_t j=0; j<other.extent(1); ++j){
  	  proxy(i,j) = srcPx(i,j);
  	}
      }
    }
    return *this;
  }

  // // move cnstr and assign
  // DenseMatrix(DenseMatrix && o);
  // DenseMatrix & operator=(DenseMatrix && o) = delete;

  // destructor
  ~DenseMatrix(){};

public:
  ord_t extent(ord_t i) const {
    assert( i <= 1 );
    return data_.shape(i);
  }

  // ref_t operator()(ord_t i, ord_t j){
  //   return mutProxy_(i,j);
  // };

  // const_ref_t operator()(ord_t i, ord_t j) const{
  //   return proxy_(i,j);
  // };

  wrap_t const * data() const{
    return &data_;
  }

  wrap_t * data(){
    return &data_;
  }

private:
  friend DenseMatrixSharedMemBase< this_t >;
  wrap_t data_ = {};
};//end class

}}//end namespace pressio::containers
#endif  // CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_SHAREDMEM_PYBIND11_HPP_