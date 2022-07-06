/*
//@HEADER
// ************************************************************************
//
// ops_clone.hpp
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

#ifndef OPS_PYBIND11_OPS_CLONE_HPP_
#define OPS_PYBIND11_OPS_CLONE_HPP_

namespace pressio{ namespace ops{

template <class T>
::pressio::mpl::enable_if_t<::pressio::is_array_pybind<T>::value, T>
clone(const T & src)
{
  if (src.ndim()==1)
  {

    const pybind11::ssize_t ext0 = ::pressio::ops::extent(src, 0);
    T result({ext0});
    using ord_t = typename ::pressio::Traits<T>::ordinal_type;
    for (ord_t i=0; i<::pressio::ops::extent(src, 0); ++i){
      result(i) = src(i);
    }
    return result;
  }
  else if (src.ndim()==2)
  {

    const pybind11::ssize_t ext0 = ::pressio::ops::extent(src, 0);
    const pybind11::ssize_t ext1 = ::pressio::ops::extent(src, 1);

    T result({ext0,  ext1});
    using ord_t = typename ::pressio::Traits<T>::ordinal_type;
    for (ord_t i=0; i<::pressio::ops::extent(src, 0); ++i){
      for (ord_t j=0; j<::pressio::ops::extent(src, 1); ++j){
	result(i,j) = src(i,j);
      }
    }
    return result;
  }
  else{

    throw std::runtime_error("clone: rank-3 not impl yet");
  }

  return T{};
}

}} //end namespace
#endif  // OPS_PYBIND11_OPS_CLONE_HPP_
