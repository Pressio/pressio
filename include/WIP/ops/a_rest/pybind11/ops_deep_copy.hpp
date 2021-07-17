/*
//@HEADER
// ************************************************************************
//
// ops_deep_copy.hpp
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

#ifndef OPS_PYBIND11_OPS_DEEP_COPY_HPP_
#define OPS_PYBIND11_OPS_DEEP_COPY_HPP_

namespace pressio{ namespace ops{

//**************
//*** RANK-1 ***
//**************
template<typename T1, typename T2>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<T1>::value and
  ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<T2>::value
  >
deep_copy(T1 & dest, const T2 & src)
{
  assert( dest.extent(0) == src.extent(0) );
  for (std::size_t i=0; i<dest.extent(0); ++i){
    dest(i) = src(i);
  }
}

template<typename T1, typename T2>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<T1>::value
  >
deep_copy(T1 & dest,
	  const ::pressio::containers::expressions::SpanExpr<T2> & src)
{
  assert( dest.extent(0) == src.extent(0) );
  for (std::size_t i=0; i<dest.extent(0); ++i){
    dest(i) = src(i);
  }
}

template<typename T1, typename T2>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<T2>::value
  >
deep_copy(::pressio::containers::expressions::SpanExpr<T1> & dest,
	  const T2 & src)
{
  assert( dest.extent(0) == src.extent(0) );
  for (std::size_t i=0; i<dest.extent(0); ++i){
    dest(i) = src(i);
  }
}

template<typename T1, typename T2>
::pressio::mpl::enable_if_t<
  T1::traits::wrapped_package_identifier == ::pressio::containers::details::WrappedPackageIdentifier::Pybind and
  T2::traits::wrapped_package_identifier == ::pressio::containers::details::WrappedPackageIdentifier::Pybind
  >
deep_copy(::pressio::containers::expressions::SpanExpr<T1> & dest,
	  const ::pressio::containers::expressions::SpanExpr<T2> & src)
{
  assert( dest.extent(0) == src.extent(0) );
  for (std::size_t i=0; i<dest.extent(0); ++i){
    dest(i) = src(i);
  }
}


//**************
//*** RANK-2 ***
//**************
template<typename T1, typename T2>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T1>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T2>::value
  >
deep_copy(T1 & dest, const T2 & src)
{
  assert( dest.extent(0) == src.extent(0) );
  assert( dest.extent(1) == src.extent(1) );
  for (std::size_t j=0; j<dest.extent(1); ++j){
    for (std::size_t i=0; i<dest.extent(0); ++i){
      dest(i,j) = src(i,j);
    }
  }
}

template<typename T1, typename T2>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_cstyle_rank2_tensor_wrapper_pybind<T1>::value and
  ::pressio::containers::predicates::is_cstyle_rank2_tensor_wrapper_pybind<T2>::value
  >
deep_copy(T1 & dest, const T2 & src)
{
  assert( dest.extent(0) == src.extent(0) );
  assert( dest.extent(1) == src.extent(1) );
  for (std::size_t i=0; i<dest.extent(0); ++i){
    for (std::size_t j=0; j<dest.extent(1); ++j){
      dest(i,j) = src(i,j);
    }
  }
}

//**************
//*** RANK-3 ***
//**************
template<typename T1, typename T2>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank3_tensor_wrapper_pybind<T1>::value and
  ::pressio::containers::predicates::is_fstyle_rank3_tensor_wrapper_pybind<T2>::value
  >
deep_copy(T1 & dest, const T2 & src)
{
  assert( dest.extent(0) == src.extent(0) );
  assert( dest.extent(1) == src.extent(1) );
  assert( dest.extent(2) == src.extent(2) );
  for (std::size_t j=0; j<dest.extent(1); ++j){
    for (std::size_t i=0; i<dest.extent(0); ++i){
      dest(i,j) = src(i,j);
    }
  }
}

}}//end namespace pressio::ops
#endif  // OPS_PYBIND11_OPS_DEEP_COPY_HPP_
