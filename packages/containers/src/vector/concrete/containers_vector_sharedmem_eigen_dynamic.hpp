/*
//@HEADER
// ************************************************************************
//
// containers_vector_sharedmem_eigen_dynamic.hpp
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

#ifndef CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_EIGEN_DYNAMIC_HPP_
#define CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_EIGEN_DYNAMIC_HPP_

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "ROL_Vector.hpp"
#endif

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::predicates::is_dynamic_vector_eigen<wrapped_type>::value
	       >
	     >
  : public VectorSharedMemBase< Vector<wrapped_type> >
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  , public ROL::Vector< typename wrapped_type::Scalar>
#endif
{

  using this_t = Vector<wrapped_type>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename  mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;
  using ref_t = typename mytraits::reference_t;
  using const_ref_t = typename mytraits::const_reference_t;

public:
  Vector() = default;

  explicit Vector(ord_t insize){
    data_.resize(insize);
    data_.setConstant( ::pressio::utils::constants<sc_t>::zero() );
  }
  explicit Vector(const wrap_t & src) : data_(src){}

  // copy cnstr
  Vector(Vector const & other) = default;
  // copy assignment
  Vector & operator=(const Vector & other) = default;
  // move cnstr
  Vector(Vector && o) = default;
  // move assignment
  Vector & operator=(Vector && other) = default;
  // destructor
  ~Vector() = default;

public:
  // assignment with value
  this_t & operator=(const sc_t & value){
    for (ord_t i = 0; i != data_.size(); ++i)
      data_[i] = value;
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this += b
  this_t & operator+=(const this_t & other) {
    assert( other.extent(0) == this->extent(0) );
    this->data_ += *other.data();
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this -= b
  this_t & operator-=(const this_t & other) {
    assert( other.extent(0) == this->extent(0) );
    this->data_ -= *other.data();
    return *this;
  }

public:

  ref_t operator [] (ord_t i){
    return data_(i);
  };
  const_ref_t operator [] (ord_t i) const{
    return data_(i);
  };

  ref_t operator()(ord_t i){
    return data_[i];
  };
  const_ref_t operator()(ord_t i) const{
    return data_[i];
  };

  wrap_t const * data() const{
    return &data_;
  }

  wrap_t * data(){
    return &data_;
  }

  wrap_t dataCp(){
    return data_;
  }

  bool empty() const{
    return this->extent(0)==0 ? true : false;
  }

  ord_t extent(ord_t i) const {
    assert( i == 0 );
    return (data_.rows()==1) ? data_.cols() : data_.rows();
  }

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  /* the following methods are needed to enable using this class
     for the Rol optimization */
  int dimension() const final{
    return data_.size();
  }
  void plus( const ROL::Vector<sc_t> &x ) final{
    const this_t & ex = static_cast<const this_t&>(x);
    data_ += *ex.data();
  }
  void axpy( const sc_t alpha, const ROL::Vector<sc_t> &x ) final{
    const this_t & ex = static_cast<const this_t&>(x);
    data_ += alpha * (*ex.data());
  }
  void scale( const sc_t alpha ) final{
    data_ *= alpha;
  }
  void zero() final{
    data_.setConstant(static_cast<sc_t>(0));
  }
  sc_t dot( const ROL::Vector<sc_t> &x ) const final{
    const this_t & ex = static_cast<const this_t&>(x);
    return data_.dot( *ex.data() );
  }
  sc_t norm() const final{
    sc_t result = 0.0;
    for (decltype(this->extent(0)) i=0; i<this->extent(0); i++)
      result += data_[i]*data_[i];
    return std::sqrt(result);
  }
  void set( const ROL::Vector<sc_t> &x ) final{
    const this_t & ex = static_cast<const this_t&>(x);
    data_ = *ex.data();
  }
  ROL::Ptr<ROL::Vector<sc_t> > clone() const final{
    return ROL::makePtr<this_t>(data_);
  }
  ROL::Ptr<ROL::Vector<sc_t>> basis( const int i ) const final{
    auto  b_ptr = clone();
    auto& b_ref = static_cast<this_t&>(*b_ptr);
    b_ref.zero();
    b_ref[i] = sc_t(1);
    return b_ptr;
  }
#endif

private:
  friend VectorSharedMemBase< this_t >;
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif  // CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_EIGEN_DYNAMIC_HPP_
