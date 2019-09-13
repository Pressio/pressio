/*
//@HEADER
// ************************************************************************
//
// containers_multi_vector_arbitrary.hpp
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

#ifndef CONTAINERS_MULTIVECTOR_CONCRETE_MULTIVECTOR_ARBITRARYHPP_
#define CONTAINERS_MULTIVECTOR_CONCRETE_MULTIVECTOR_ARBITRARYHPP_

#include "../../shared_base/containers_container_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class MultiVector<
  wrapped_type,
  mpl::enable_if_t<
    details::traits<MultiVector<wrapped_type>>::wrapped_multi_vector_identifier
    == details::WrappedMultiVectorIdentifier::Arbitrary
    >
  >
  : public ContainerBase< MultiVector<wrapped_type>, wrapped_type >
{

  using this_t = MultiVector<wrapped_type>;

public:
  MultiVector() = delete;
  ~MultiVector() = delete;

  template <typename ...Args>
  MultiVector(Args && ... args)
    : data_( std::forward<Args>(args)... ){}

  explicit MultiVector(const wrap_t & vecobj)
    : data_(vecobj){}

  MultiVector(this_t const & other)
    : data_(*other.data()){}

private:
  wrapped_type const * dataImpl() const{
    return &data_;
  }

  wrapped_type * dataImpl(){
    return &data_;
  }

private:
  friend ContainerBase< this_t, wrapped_type >;

  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers

#endif
