/*
//@HEADER
// ************************************************************************
//
// ops_abs.hpp
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

#ifndef OPS_PYBIND11_OPS_ABS_HPP_
#define OPS_PYBIND11_OPS_ABS_HPP_

namespace pressio{ namespace ops{

// to = abs(from)
template <class T1, class T2>
::pressio::mpl::enable_if_t<
  (::pressio::Traits<T1>::package_identifier == PackageIdentifier::Pybind and
   ::pressio::Traits<T2>::package_identifier == PackageIdentifier::Pybind)
  >
abs(T1 & to, const T2 & from)
{
  assert(to.ndim() == from.ndim());

  using size_t = typename ::pressio::Traits<T1>::ordinal_type;
  assert(extent(to,0)==extent(from,0));
  if (to.ndim() == 1)
  {
    for (size_t i=0; i<extent(to, 0); ++i){
      to(i) = std::abs(from(i));
    }
  }
  else{
    throw std::runtime_error("abs: rank-2, rank-3 not impl yet");
  }
}

}}//end namespace pressio::ops
#endif  // OPS_PYBIND11_OPS_ABS_HPP_
