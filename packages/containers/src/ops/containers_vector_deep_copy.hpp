/*
//@HEADER
// ************************************************************************
//
// containers_vector_deep_copy.hpp
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

#ifndef CONTAINERS_CONTAINER_OPS_VECTOR_DEEP_COPY_HPP_
#define CONTAINERS_CONTAINER_OPS_VECTOR_DEEP_COPY_HPP_

#include "../vector/containers_vector_meta.hpp"
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

namespace pressio{ namespace containers{ namespace ops{

//--------------------------------------------------------------------------
// for wrappers, because we overload the = operator
// and for now we do NOT have view semantics
//--------------------------------------------------------------------------
template<
  typename T,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper<T>::value
    > * = nullptr
  >
void deep_copy(const T & src, T & dest){
  dest = src;
}

template<
  typename T1, typename T2,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper<T1>::value and
    ::pressio::containers::meta::is_expression<T2>::value
    > * = nullptr
  >
void deep_copy(const T1 & src, T2 & dest){
  for (auto i=0; i<src.size(); ++i)
    dest[i] = src[i];
}


//--------------------------------------------------------------------------
// enable for pybind11::array_t
//--------------------------------------------------------------------------
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template<
  typename T,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind11<T>::value
    > * = nullptr
  >
void deep_copy(const T & src, T & dest){

  if (src.ndim() > 1){
    throw std::runtime_error("containers::ops::deep_copy: v.ndims()!=1. this operation currently supported for vectors only");
  }

  const auto vsz = src.size();
  assert(vsz == dest.size());
  auto dest_proxy = dest.mutable_unchecked();
  auto src_proxy  = src.unchecked();
  for (std::size_t i=0; i<(std::size_t)vsz; ++i){
    dest_proxy(i) = src_proxy(i);
  }
}
#endif


}}}//end namespace pressio::containers::ops
#endif
