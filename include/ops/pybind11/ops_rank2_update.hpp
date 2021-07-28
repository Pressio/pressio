/*
//@HEADER
// ************************************************************************
//
// ops_rank2_update.hpp
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

#ifndef OPS_PYBIND11_OPS_RANK2_UPDATE_HPP_
#define OPS_PYBIND11_OPS_RANK2_UPDATE_HPP_

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// computing:  M = a * M + b * M1
//----------------------------------------------------------------------
template<typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T1>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T2>::value
  >
update(T1 & M, scalar_t a,
       const T2 & M1, scalar_t b)
{
  using size_type = typename T1::traits::size_t;
  assert( impl::_matching_extents(M,M1) );

  for (size_type j=0; j<M.extent(1); ++j){
    for (size_type i=0; i<M.extent(0); ++i){
      M(i,j) = a*M(i,j) + b*M1(i,j);
    }
  }
}

template<typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T1>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T2>::value
  >
update(T1 & M, const T2 & M1, const scalar_t b)
{
  using size_type = typename T1::traits::size_t;
  assert( impl::_matching_extents(M,M1) );

  for (size_type j=0; j<M.extent(1); ++j){
    for (size_type i=0; i<M.extent(0); ++i){
      M(i,j) = b*M1(i,j);
    }
  }
}

//----------------------------------------------------------------------
//  overloads for computing:
//  M = a * M + b * M1 + c * M2
//----------------------------------------------------------------------
template<typename T, typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T1>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T2>::value
  >
update(T & M, scalar_t a,
       const T1 & M1, scalar_t b,
       const T2 & M2, scalar_t c)
{
  using size_type = typename T1::traits::size_t;
  assert( impl::_matching_extents(M,M1,M2) );

  for (size_type j=0; j<M.extent(1); ++j){
    for (size_type i=0; i<M.extent(0); ++i){
      M(i,j) = a*M(i,j) + b*M1(i,j) + c*M2(i,j);
    }
  }
}

//----------------------------------------------------------------------
//  overloads for computing:
//  M = b * M1 + c * M2
//----------------------------------------------------------------------
template<typename T, typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T1>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T2>::value
  >
update(T & M,
       const T1 & M1, const scalar_t &b,
       const T2 & M2, const scalar_t &c)
{
  using size_type = typename T1::traits::size_t;
  assert( impl::_matching_extents(M,M1,M2) );

  for (size_type j=0; j<M.extent(1); ++j){
    for (size_type i=0; i<M.extent(0); ++i){
      M(i,j) = b*M1(i,j) + c*M2(i,j);
    }
  }
}

//----------------------------------------------------------------------
//  overloads for computing:
//	M = a * M + b * M1 + c * M2 + d * M3
//----------------------------------------------------------------------
template<typename T, typename T1, typename T2, typename T3, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T1>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T2>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T3>::value
  >
update(T & M, const scalar_t &a,
       const T1 & M1, const scalar_t &b,
       const T2 & M2, const scalar_t &c,
       const T3 & M3, const scalar_t &d)
{
  using size_type = typename T1::traits::size_t;
  assert( impl::_matching_extents(M,M1,M2,M3) );

  for (size_type j=0; j<M.extent(1); ++j){
    for (size_type i=0; i<M.extent(0); ++i){
      M(i,j) = a*M(i,j) + b*M1(i,j) + c*M2(i,j) + d*M3(i,j);
    }
  }
}

//----------------------------------------------------------------------
//  overloads for computing:
//  M = a * M + b * M1 + c * M2 + d * M3 + e * M4
//----------------------------------------------------------------------
template<typename T, typename T1, typename T2, typename T3, typename T4, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T1>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T2>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T3>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<T4>::value
  >
update(T & M, const scalar_t &a,
       const T1 & M1, const scalar_t &b,
       const T2 & M2, const scalar_t &c,
       const T3 & M3, const scalar_t &d,
       const T4 & M4, const scalar_t &e)
{
  using size_type = typename T1::traits::size_t;
  assert( impl::_matching_extents(M,M1,M2,M3,M4) );

  for (size_type j=0; j<M.extent(1); ++j){
    for (size_type i=0; i<M.extent(0); ++i){
      M(i,j) = a*M(i,j) + b*M1(i,j) + c*M2(i,j) + d*M3(i,j) + e*M4(i,j);
    }
  }
}

}}//end namespace pressio::ops
#endif  // OPS_PYBIND11_OPS_RANK2_UPDATE_HPP_
