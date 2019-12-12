/*
//@HEADER
// ************************************************************************
//
// containers_multi_vector_sharedmem_eigen_dynamic.hpp
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

#ifndef CONTAINERS_MULTIVECTOR_CONCRETE_MULTIVECTOR_SHAREDMEM_EIGEN_DYNAMIC_HPP_
#define CONTAINERS_MULTIVECTOR_CONCRETE_MULTIVECTOR_SHAREDMEM_EIGEN_DYNAMIC_HPP_

#include "../base/containers_multi_vector_sharedmem_base.hpp"
// #include "../../meta/containers_native_multi_vector_meta.hpp"
#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_resizable_base.hpp"
#include "../../shared_base/containers_container_subscriptable_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class MultiVector<wrapped_type,
		  ::pressio::mpl::enable_if_t<
		    meta::is_dynamic_multi_vector_eigen<
		      wrapped_type>::value>
		  >
  : public ContainerBase< MultiVector<wrapped_type>, wrapped_type >,
    public MultiVectorSharedMemBase< MultiVector<wrapped_type> >,
    public ContainerSubscriptable2DBase<
     MultiVector<wrapped_type>,
     typename details::traits<MultiVector<wrapped_type>>::scalar_t,
     typename details::traits<MultiVector<wrapped_type>>::ordinal_t>,
    public ContainerResizableBase<MultiVector<wrapped_type>, 2>
{

private:
  using this_t = MultiVector<wrapped_type>;
  using sc_t = typename details::traits<this_t>::scalar_t;
  using ord_t = typename details::traits<this_t>::ordinal_t;
  using wrap_t = typename details::traits<this_t>::wrapped_t;

public:
  MultiVector() = delete;

  MultiVector(ord_t length, ord_t numVectors)
    : data_(length, numVectors){
    this->setZero();
  }

  explicit MultiVector(const wrap_t & other)
    : data_(other){}

  // copy cnstr
  MultiVector(MultiVector const & other) = default;
  // copy assignment
  MultiVector & operator=(const MultiVector & other) = default;
  // move cnstr
  MultiVector(MultiVector && o) = default;
  // move assignment
  MultiVector & operator=(MultiVector && other) = default;
  // destructor
  ~MultiVector() = default;


public:
  sc_t & operator()(ord_t irow, ord_t iVec){
    assert(iVec < this->numVectors() );
    assert(irow < this->length() );
    return data_(irow, iVec);
  }

  sc_t const & operator()(ord_t irow, ord_t iVec)const{
    assert(iVec < this->numVectors() );
    assert(irow < this->length() );
    return data_(irow, iVec);
  }

private:

  ord_t numVectorsImpl() const{
    return data_.cols();
  }

  ord_t lengthImpl() const {
    return data_.rows();
  };

  wrap_t * dataImpl(){
    return &data_;
  };

  wrap_t const * dataImpl() const{
    return &data_;
  };

  wrap_t dataCpImpl() const{
    return data_;
  };

  bool emptyImpl() const {
    return this->length()==0 ? true : false;
  }

  void scaleImpl(sc_t & factor){
    data_.coeffs() *= factor;
  };

  void setZeroImpl() {
    data_.setConstant(static_cast<sc_t>(0));
  }

  void resizeImpl(ord_t newlength, ord_t nVec){
    data_.resize(newlength, nVec);
    this->setZero();
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend MultiVectorSharedMemBase< this_t >;
  friend ContainerSubscriptable2DBase< this_t, sc_t, ord_t>;
  friend ContainerResizableBase<this_t, 2>;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif
