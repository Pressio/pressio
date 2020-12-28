/*
//@HEADER
// ************************************************************************
//
// containers_vector_distributed_tpetra_block.hpp
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

#ifndef CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_DISTRIBUTED_TPETRA_BLOCK_HPP_
#define CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_DISTRIBUTED_TPETRA_BLOCK_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<
  wrapped_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_vector_tpetra_block<
      wrapped_type
      >::value
    >
  >
{
public:
  using this_t	= Vector<wrapped_type>;
  using traits	= details::traits<this_t>;
  using sc_t	= typename details::traits<this_t>::scalar_t;
  using LO_t	= typename details::traits<this_t>::local_ordinal_t;
  using GO_t	= typename details::traits<this_t>::global_ordinal_t;
  using device_t= typename details::traits<this_t>::device_t;
  using map_t	= typename details::traits<this_t>::data_map_t;

public:
  Vector() = delete;

  // something seems wrong with
  //    data_(other.data_, Teuchos::Copy){}
  // for tpetra block, so don't use that.
  // just construct the vecrtor using map and block size.

  explicit Vector(const wrapped_type & vecobj)
    : data_( *vecobj.getMap(), vecobj.getBlockSize())
  {
    // just a trick to copy data
    data_.update(::pressio::utils::constants<sc_t>::one(),
     vecobj,
     ::pressio::utils::constants<sc_t>::zero());
  }

  explicit Vector(wrapped_type && vecobj)
    : Vector(vecobj){}

  explicit Vector(const map_t & mapO, LO_t blockSize)
    : data_(mapO, blockSize){}


  // copy cnstr delegating (for now) to the one above
  Vector(Vector const & other) : Vector(*other.data()){}

  // delete copy assign to force usage of ops::deep_copy
  Vector & operator=(const Vector & other) = delete;


  /*something is wrong with the move for blockVector,
    I get all sorts of strange behaviors.

    for example this:
    using vec_t = tpetra::BlockVector<>;
    vec_t a(*map, 3);
    a.putScalar(1.);
    ASSERT_EQ( a.getVectorView().norm1(), 45. );
    vec_t b(std::move(a));
    ASSERT_EQ( b.getVectorView().norm1(), 45. );

    for now let's not implement move for tpetra block wrapper.
  */
  // move cnstr
  Vector(Vector && other) : Vector(*other.data()){}
  //   : data_(other.data_, Teuchos::Copy){}

  // move assignment
  Vector & operator=(Vector && other)
  {
    this->data_.update(::pressio::utils::constants<sc_t>::one(),
           *other.data(),
           ::pressio::utils::constants<sc_t>::zero() );
    return *this;
  }

  ~Vector() = default;

public:
  wrapped_type const * data() const{
    return &data_;
  }

  wrapped_type * data(){
    return &data_;
  }

  wrapped_type dataCp(){
    return data_;
  }

  GO_t extent(std::size_t i) const{
    assert(i==0);
    return data_.getMap()->getGlobalNumElements();
  }

  LO_t extentLocal(std::size_t i) const{
    assert(i==0);
    return data_.getMap()->getNodeNumElements();
  }

private:
  wrapped_type data_ = {};

};

}}//end namespace pressio::containers
#endif  // CONTAINERS_VECTOR_CONCRETE_CONTAINERS_VECTOR_DISTRIBUTED_TPETRA_BLOCK_HPP_
