/*
//@HEADER
// ************************************************************************
//
// containers_vector_distributed_tpetra.hpp
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

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#ifndef CONTAINERS_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_TPETRA_HPP_
#define CONTAINERS_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_TPETRA_HPP_

#include <MatrixMarket_Tpetra.hpp>

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<wrapped_type,
	     typename
	     std::enable_if<
	       meta::is_vector_tpetra<
		 wrapped_type>::value
	       >::type
	     >
  : public VectorDistributedBase< Vector<wrapped_type> >
{

  using this_t = Vector<wrapped_type>;
  using sc_t = typename details::traits<this_t>::scalar_t;
  using LO_t = typename details::traits<this_t>::local_ordinal_t;
  using GO_t = typename details::traits<this_t>::global_ordinal_t;
  using device_t = typename details::traits<this_t>::device_t;
  using der_t = this_t;
  using wrap_t = typename details::traits<this_t>::wrapped_t;
  using map_t = typename details::traits<this_t>::data_map_t;
  using mpicomm_t = typename details::traits<this_t>::communicator_t;

public:
  // default cnstr
  Vector() = delete;

  // cnstrs
  explicit Vector(const wrap_t & vecobj)
    // use the deep_copy cnstr
    : data_(vecobj, Teuchos::Copy){}

  explicit Vector(Teuchos::RCP<const map_t> mapO)
    : data_(mapO){}

  // here we do not default the copy and move because if we did that,
  // it would use the tpetra copy/move which have view semantics
  // which is not what we want here (for the time being)

  // copy cnstr
  Vector(Vector const & other) : data_(*other.data(), Teuchos::Copy){}
  // copy assignment
  Vector & operator=(const Vector & other){
    if (&other != this){
      assert(this->extentLocal(0) == other.extentLocal(0));
      data_.assign( *other.data() );
    }
    return *this;
  }

  // move cnstr
  Vector(Vector && other) : data_(*other.data(), Teuchos::Copy){}

  // move assignment
  Vector & operator=(Vector && other){
    assert(this->extentLocal(0) == other.extentLocal(0));
    data_.assign( *other.data() );
    return *this;
  }

  // destructor
  ~Vector() = default;

public:
  // compound assignment when type(b) = type(this)
  // this += b
  this_t & operator+=(const this_t & other) {
    this->data_.update(1.0, *other.data(), 1.0 );
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this -= b
  this_t & operator-=(const this_t & other) {
    this->data_.update(-1.0, *other.data(), 1.0 );
    return *this;
  }

  void print(std::string tag) const{
    Tpetra::MatrixMarket::Writer<wrap_t>::writeDense
      (std::cout << std::setprecision(15), data_, tag, tag);
  }

  wrap_t const * data() const{
    return &data_;
  }

  wrap_t * data(){
    return &data_;
  }

  wrap_t dataCp(){
    return data_;
  }

  GO_t extent(std::size_t i) const{
    assert(i==0);
    return data_.getGlobalLength();
  }

  LO_t extentLocal(std::size_t i) const{
    assert(i==0);
    return data_.getLocalLength();
  }

private:
  friend VectorDistributedBase< this_t >;
  wrap_t data_ = {};
};//end class

}}//end namespace pressio::containers
#endif
#endif
