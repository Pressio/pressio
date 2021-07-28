/*
//@HEADER
// ************************************************************************
//
// containers_vector_sharedmem_pybind11.hpp
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

#ifndef CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_PYBIND11_HPP_
#define CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_PYBIND11_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<
  wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_array_pybind<wrapped_type>::value
    >
  >
{

  using this_t	    = Vector<wrapped_type>;
  using traits      = details::traits<this_t>;
  using sc_t	    = typename traits::scalar_t;
  using ord_t	    = typename traits::ordinal_t;
  using wrap_t	    = typename traits::wrapped_t;
  using data_r_t    = typename traits::data_return_t;
  using const_data_r_t  = typename traits::const_data_return_t;
  using ref_t      = typename traits::reference_t;
  using const_ref_t = typename traits::const_reference_t;
  using mut_proxy_t = typename traits::mut_proxy_t;
  using proxy_t	    = typename traits::proxy_t;

public:
  Vector() = delete;

  template <
    typename int_t,
    mpl::enable_if_t<
     std::is_integral<int_t>::value, std::size_t
     > = 0
    >
  explicit Vector(int_t insize)
    : data_(insize)
  {
    auto proxy = data_.mutable_unchecked();
    for (ord_t i=0; i<insize; ++i)
      proxy(i) = static_cast<sc_t>(0);
  }

  explicit Vector(const wrap_t & src)
    : data_{ wrap_t(const_cast<wrap_t &>(src).request()) }
  {
    // src must be a vector to be wraped into a vector
    assert( data_.ndim() == 1 );

    // copy data from src to this
    auto proxy = data_.mutable_unchecked();
    const auto srcPx = src.unchecked();
    for (ord_t i=0; i<src.size(); ++i)
      proxy(i) = srcPx(i);
  }

  // use only if you know what you are doing
  // it is currently used only in specific places
  Vector(wrap_t src, ::pressio::view)
    : data_{src}
  {
    assert( data_.ndim() == 1 );
  }

  // copy cnstr
  Vector(Vector const & other)
    : data_(other.extent(0))
  {
    assert( other.data_.ndim() == 1 );
    // copy data from src to this
    auto proxy = data_.mutable_unchecked();
    const auto srcPx = other.data_.unchecked();
    for (ord_t i=0; i<other.extent(0); ++i)
      proxy(i) = srcPx(i);
  }

  // delete copy assign to force usage of ops::deep_copy
  Vector & operator=(const Vector & other) = delete;
  //   if (&other != this){
  //     assert( other.ndim() == 1 );
  //     assert(this->extent(0) == other.extent(0));

  //     // copy data from src to this
  //     auto proxy = data_.mutable_unchecked();
  //     const auto srcPx = other.data_.unchecked();
  //     for (ord_t i=0; i<other.extent(0); ++i)
  // 	proxy_(i) = srcPx(i);
  //   }
  //   return *this;
  // }

  // move cnstr and assign
  Vector(Vector && other) = default;
  Vector & operator=(Vector && o) = delete;

  // destructor
  ~Vector(){};

public:
  const_data_r_t data() const{
    return &data_;
  }

  data_r_t data(){
    return &data_;
  }

  bool empty() const{
    return this->extent(0)==0 ? true : false;
  }

  proxy_t proxy() const{
    return data_.unchecked();
  }

  mut_proxy_t proxy(){
    return data_.mutable_unchecked();
  }

  ord_t extent(ord_t i) const {
    assert( i == 0 );
    return data_.shape(0);
  }

  ref_t operator()(ord_t i){
    assert(i < this->extent(0));
    return data_(i);
  };

  const_ref_t operator()(ord_t i) const{
    assert(i < this->extent(0));
    return data_(i);
  };

  // [[deprecated("Use operator() instead.")]]
  // ref_t operator[](ord_t i){
  //   return (*this)(i);
  // };
  // [[deprecated("Use operator() instead.")]]
  // const_ref_t operator[](ord_t i) const{
  //   return (*this)(i);
  // };

private:
  wrap_t data_ = {};
};

}}//end namespace pressio::containers
#endif  // CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_SHAREDMEM_PYBIND11_HPP_
