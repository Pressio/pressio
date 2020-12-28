/*
//@HEADER
// ************************************************************************
//
// containers_vector_distributed_epetra.hpp
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

#ifndef CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_DISTRIBUTED_EPETRA_HPP_
#define CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_DISTRIBUTED_EPETRA_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<
  wrapped_type,
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_vector_epetra<
      wrapped_type>::value
    >
  >
{

public:
  using this_t = Vector<wrapped_type>;
  using traits = details::traits<this_t>;
  using sc_t = typename details::traits<this_t>::scalar_t;
  using LO_t = typename details::traits<this_t>::local_ordinal_t;
  using GO_t = typename details::traits<this_t>::global_ordinal_t;
  using map_t = typename details::traits<this_t>::data_map_t;

public:
  // default cnstr
  Vector() = delete;

  // cnstrs
  explicit Vector(const map_t & mapobj) : data_(mapobj){}
  explicit Vector(const wrapped_type & vecobj) : data_(vecobj){}

  // copy cnstr
  Vector(Vector const & other) = default;

  // delete copy assign to force usage of ops::deep_copy
  Vector & operator=(Vector const & other) = delete;

  // move cnstr
  Vector(Vector && other) = default;
  // move assignment
  Vector & operator=(Vector && other) = default;
  // destructor
  ~Vector() = default;

public:
  wrapped_type const * data() const{
    return &data_;
  }

  wrapped_type * data(){
    return &data_;
  }

  sc_t & operator()(LO_t i){
    assert(i < this->extentLocal(0));
    return data_[i];
  };
  sc_t const & operator()(LO_t i) const{
    assert(i < this->extentLocal(0));
    return data_[i];
  };

  [[deprecated("Use operator() instead.")]]
  sc_t & operator[](LO_t i){
    assert(i < this->extentLocal(0));
    return data_[i];
  };
  [[deprecated("Use operator() instead.")]]
  sc_t const & operator[](LO_t i) const{
    assert(i < this->extentLocal(0));
    return data_[i];
  };

  GO_t extent(std::size_t i) const{
    assert(i==0);
    return data_.GlobalLength();
  }

  LO_t extentLocal(std::size_t i) const{
    assert(i==0);
    return data_.MyLength();
  }

private:
  wrapped_type data_ = {};

};

}}//end namespace pressio::containers
#endif  // CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_DISTRIBUTED_EPETRA_HPP_
