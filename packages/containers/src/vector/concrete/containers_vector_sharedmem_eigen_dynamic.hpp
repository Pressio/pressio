/*
//@HEADER
// ************************************************************************
//
// containers_vector_sharedmem_eigen_dynamic.hpp
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

#ifndef CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_EIGEN_DYNAMIC_HPP_
#define CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_EIGEN_DYNAMIC_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_resizable_base.hpp"
#include "../../shared_base/containers_container_subscriptable_base.hpp"
#include "../../shared_base/containers_container_printable_base.hpp"
#include "../base/containers_vector_sharedmem_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::meta::is_dynamic_vector_eigen<wrapped_type>::value
	       >
	     >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public VectorSharedMemBase< Vector<wrapped_type> >,
    public ContainerResizableBase<Vector<wrapped_type>, 1>,
    public ContainerSubscriptable1DBase<
	       Vector<wrapped_type>,
	       typename details::traits<Vector<wrapped_type>>::scalar_t,
	       typename details::traits<Vector<wrapped_type>>::ordinal_t>,
    public ContainerPrintable1DBase<
	       Vector<wrapped_type>,
	       typename details::traits<Vector<wrapped_type>>::ordinal_t>{

  using this_t = Vector<wrapped_type>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename  mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;

public:
  Vector() = default;

  ~Vector() = default;

  explicit Vector(ord_t insize){
    this->resize(insize);
    this->setZero();
  }

  explicit Vector(const wrap_t & src)
    : data_(src){}

  Vector(this_t const & other)
    : data_(*other.data()){
  }


public:

  // constructor from any expression, force evaluation
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  explicit Vector(const T & expr){
    this->resize(expr.size());
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] = expr(i);
  }

  // assignment from any expression, force evaluation
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  this_t & operator=(const T & expr){
    if(this->size() != expr.size())
      this->resize(expr.size());
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] = expr(i);
    return *this;
  }

  // assignment with other vector of same type
  this_t & operator=(const this_t & other){
    if(this->size() != other.size())
       this->resize( other.size() );
    data_ = *other.data();
    return *this;
  }

  // assignment with value
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
	      std::is_same<T, sc_t>::value> * = nullptr>
  this_t & operator=(const T value){
    for (ord_t i = 0; i != data_.size(); ++i)
      data_[i] = value;
    return *this;
  }


  // compound assignment from expression template
  // this += expr
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  this_t & operator+=(const T & expr) {
    assert( expr.size() == this->size() );
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] += expr(i);
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this += b
  this_t & operator+=(const this_t & other) {
    assert( other.size() == this->size() );
    this->data_ += *other.data();
    return *this;
  }


  // compound assignment from expression template
  // this -= expr
  template <typename T,
  	    ::pressio::mpl::enable_if_t<
  	      T::is_vector_expression> * = nullptr>
  this_t & operator-=(const T & expr) {
    assert( expr.size() == this->size() );
    for (ord_t i = 0; i != expr.size(); ++i)
      data_[i] -= expr(i);
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this -= b
  this_t & operator-=(const this_t & other) {
    assert( other.size() == this->size() );
    this->data_ -= *other.data();
    return *this;
  }

public:

  sc_t & operator [] (ord_t i){
    return data_(i);
  };

  sc_t const & operator [] (ord_t i) const{
    return data_(i);
  };

  sc_t & operator()(ord_t i){
    return data_[i];
  };
  sc_t const & operator()(ord_t i) const{
    return data_[i];
  };

private:

  template <typename stream_t>
  void printImpl(stream_t & os, char c, ord_t nIn) const{
    assert(nIn <= this->size());
    auto nToPrint = (nIn==-1) ? this->size() : nIn;
    ::pressio::utils::impl::setStreamPrecision<stream_t, sc_t>(os);

    if (c=='d')
      this->printVertically(os, nToPrint);
    if (c=='f')
      this->printFlatten(os, nToPrint);
  }

  template <typename stream_t>
  void printVertically(stream_t & os, ord_t nToPrint) const{
    for (auto i=0; i<nToPrint; i++)
      os << data_[i] << "\n";
    os << std::endl;
  }

  template <typename stream_t>
  void printFlatten(stream_t & os, ord_t nToPrint) const{
    for (auto i=0; i<nToPrint; i++)
      os << data_[i] << " ";
    os << std::endl;
  }

  void matchLayoutWithImpl(const this_t & other){
    this->resize( other.size() );
  }

  wrap_t const * dataImpl() const{
    return &data_;
  }

  void scaleImpl(sc_t value) {
    data_ *= value;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  void putScalarImpl(sc_t value) {
    data_.setConstant(value);
  }

  void setZeroImpl() {
    this->putScalarImpl( static_cast<sc_t>(0) );
  }

  bool emptyImpl() const{
    return this->size()==0 ? true : false;
  }

  ord_t sizeImpl() const {
    return (data_.rows()==1) ? data_.cols() : data_.rows();
  }

  void resizeImpl(ord_t val){
    data_.resize(val);
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorSharedMemBase< this_t >;
  friend ContainerResizableBase<this_t, 1>;
  friend ContainerSubscriptable1DBase<this_t, sc_t, ord_t>;
  friend ContainerPrintable1DBase<this_t, ord_t>;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif
