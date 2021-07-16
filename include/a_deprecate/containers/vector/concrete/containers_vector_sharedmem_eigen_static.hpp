/*
//@HEADER
// ************************************************************************
//
// containers_vector_sharedmem_eigen_static.hpp
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

#ifndef CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_EIGEN_STATIC_HPP_
#define CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_EIGEN_STATIC_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<
  wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_static_vector_eigen<wrapped_type>::value
    >
  >
{
public:
  using this_t = Vector<wrapped_type>;
  using traits = details::traits<this_t>;
  using sc_t = typename traits::scalar_t;
  using ord_t = typename  traits::ordinal_t;
  using wrap_t = typename traits::wrapped_t;
  using ref_t = typename traits::reference_t;
  using const_ref_t = typename traits::const_reference_t;

public:
  Vector(){
    data_.setConstant(0);
  }

  explicit Vector(const wrap_t & src) : data_(src){}

  Vector(wrap_t && src) : data_(std::move(src)){}

  // copy cnstr
  Vector(Vector const & other) = default;
  // delete copy assign to force usage of ops::deep_copy
  Vector & operator=(Vector const & other) = delete;

  /* move semantics for static vectors do not invalidate
     the moved-from object, they just copy the data, see:
     https://gitlab.com/libeigen/eigen/-/issues/2000
  */
  // move cnstr
  Vector(Vector && other) = default;
  // move assignment
  Vector & operator=(Vector && other) = default;

  // destructor
  ~Vector() = default;

public:
  wrap_t const * data() const{
    return &data_;
  }

  wrap_t * data(){
    return &data_;
  }

  ref_t operator()(ord_t i){
    assert(i < this->extent(0));
    return data_(i);
  };
  const_ref_t operator()(ord_t i) const{
    assert(i < this->extent(0));
    return data_(i);
  };

  [[deprecated("Use operator() instead.")]]
  ref_t operator[](ord_t i){
    return (*this)(i);
  };
  [[deprecated("Use operator() instead.")]]
  const_ref_t operator[](ord_t i) const{
    return (*this)(i);
  };

  bool empty() const{
    return this->size()==0 ? true : false;
  }

  ord_t extent(ord_t i) const {
    assert( i == 0 );
    return (data_.rows()==1) ? data_.cols() : data_.rows();
  }

private:
  wrap_t data_ = {};

};
}}//end namespace pressio::containers
#endif  // CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_EIGEN_STATIC_HPP_
