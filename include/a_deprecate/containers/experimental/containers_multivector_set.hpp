/*
//@HEADER
// ************************************************************************
//
// containers_multivector_set.hpp
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

#ifndef CONTAINERS_EXPERIMENTAL_CONTAINERS_MULTIVECTOR_SET_HPP_
#define CONTAINERS_EXPERIMENTAL_CONTAINERS_MULTIVECTOR_SET_HPP_

namespace pressio{ namespace containers{

namespace experimental{
template<typename T>
class MultiVectorSet
{
public:
  using traits = ::pressio::containers::details::traits<MultiVectorSet<T>>;
  using mv_type	= ::pressio::containers::MultiVector<T>;
  using multivector_type = mv_type;
  using data_type = std::vector<mv_type>;

private:
  data_type data_;

public:
  std::size_t size() const{ return data_.size(); }

  mv_type & operator()(std::size_t i){
    return data_[i];
  }

  mv_type const & operator()(std::size_t i) const{
    return data_[i];
  }

public:
  template <
  typename _mv_type = mv_type,
  mpl::enable_if_t<std::is_default_constructible<_mv_type>::value, int> = 0
  >
  MultiVectorSet(){};

  template<typename ...Args>
  MultiVectorSet(std::size_t n, Args && ...args)
    : data_(n, mv_type(std::forward<Args>(args)...))
  {}

  MultiVectorSet(MultiVectorSet const & other) = default;
  MultiVectorSet & operator=(MultiVectorSet const & other) = default;
  MultiVectorSet(MultiVectorSet && other) = default;
  MultiVectorSet & operator=(MultiVectorSet && other) = default;
  ~MultiVectorSet() = default;
};

}// end namespace experimental

namespace details{

template<typename T>
struct traits<
  ::pressio::containers::experimental::MultiVectorSet<T>
>
{
  using mv_t = typename ::pressio::containers::experimental::MultiVectorSet<T>::mv_type;
  using scalar_t = typename traits<mv_t>::scalar_t;
};
}// end namespace details

}} //end namespace pressio::containers
#endif  // CONTAINERS_EXPERIMENTAL_CONTAINERS_MULTIVECTOR_SET_HPP_
