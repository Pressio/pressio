/*
//@HEADER
// ************************************************************************
//
// ops_multi_vector_update.hpp
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

#ifndef OPS_EIGEN_OPS_MULTI_VECTOR_UPDATE_HPP_
#define OPS_EIGEN_OPS_MULTI_VECTOR_UPDATE_HPP_

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
//  overloads for computing: MV = a * MV + b * MV1
// where MV is an eigen multivector wrapper
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  containers::predicates::is_multi_vector_wrapper_eigen<T>::value
>
update(T & mv, const scalar_t &a,
       const T & mv1, const scalar_t &b)
{
  assert( mv.extent(0) == mv1.extent(0) );
  assert( mv.extent(1) == mv1.extent(1) );
  *mv.data() = a * (*mv.data()) + b * (*mv1.data());
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T>::value
>
update(T & mv, const T & mv1, const scalar_t & b)
{
  assert( mv.extent(0) == mv1.extent(0) );
  assert( mv.extent(1) == mv1.extent(1) );
  *mv.data() = b * (*mv1.data());
}

//----------------------------------------------------------------------
//  overload for: MV = a * MV + b * MV1 + c * MV2
//----------------------------------------------------------------------
template<typename T, typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T1>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T2>::value
  >
update(T & mv, const scalar_t &a,
       const T1 & mv1, const scalar_t &b,
       const T2 & mv2, const scalar_t &c)
{
  static_assert
    (::pressio::containers::predicates::are_scalar_compatible<T,T1,T2>::value,
     "types T,T1,T2 are not scalar compatible");
  (*mv.data()) = a*(*mv.data()) + b*(*mv1.data()) + c*(*mv2.data());
}

template<typename T, typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T1>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T2>::value
  >
update(T & mv,
       const T1 & mv1, const scalar_t &b,
       const T2 & mv2, const scalar_t &c)
{
  static_assert
    (::pressio::containers::predicates::are_scalar_compatible<T,T1,T2>::value,
     "types T,T1,T2 are not scalar compatible");
  (*mv.data()) = b*(*mv1.data()) + c*(*mv2.data());
}

//----------------------------------------------------------------------
//  overload for: MV = a * MV + b * MV1 + c * MV2 + d * MV3
//----------------------------------------------------------------------
template<typename T, typename T1, typename T2, typename T3, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T1>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T2>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T3>::value
  >
update(T & mv, const scalar_t &a,
       const T1 & mv1, const scalar_t &b,
       const T2 & mv2, const scalar_t &c,
       const T3 & mv3, const scalar_t &d)
{
  static_assert
    (::pressio::containers::predicates::are_scalar_compatible<T,T1,T2,T3>::value,
     "types are not scalar compatible");
  (*mv.data()) = a*(*mv.data()) + b*(*mv1.data()) + c*(*mv2.data()) + d*(*mv3.data());
}

template<typename T, typename T1, typename T2, typename T3, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T1>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T2>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T3>::value
  >
update(T & mv,
       const T1 & mv1, const scalar_t &b,
       const T2 & mv2, const scalar_t &c,
       const T3 & mv3, const scalar_t &d)
{
  static_assert
    (::pressio::containers::predicates::are_scalar_compatible<T,T1,T2>::value,
     "types are not scalar compatible");
  (*mv.data()) = b*(*mv1.data()) + c*(*mv2.data()) + d*(*mv3.data());
}

//----------------------------------------------------------------------
//  overload for: MV = a * MV + b * MV1 + c * MV2 + d * MV3 + e*MV4
//----------------------------------------------------------------------
template<typename T, typename T1, typename T2, typename T3, typename T4, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T1>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T2>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T3>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T4>::value
  >
update(T & mv, const scalar_t &a,
       const T1 & mv1, const scalar_t &b,
       const T2 & mv2, const scalar_t &c,
       const T3 & mv3, const scalar_t &d,
       const T4 & mv4, const scalar_t &e)
{
  static_assert
    (::pressio::containers::predicates::are_scalar_compatible<T,T1,T2,T3,T4>::value,
     "types are not scalar compatible");
  (*mv.data()) = a*(*mv.data()) + b*(*mv1.data()) +
    c*(*mv2.data()) + d*(*mv3.data()) + e*(*mv4.data());
}

template<typename T, typename T1, typename T2, typename T3, typename T4, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T1>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T2>::value and
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T3>::value and
    ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T4>::value
  >
update(T & mv,
       const T1 & mv1, const scalar_t &b,
       const T2 & mv2, const scalar_t &c,
       const T3 & mv3, const scalar_t &d,
       const T4 & mv4, const scalar_t &e)
{
  static_assert
    (::pressio::containers::predicates::are_scalar_compatible<T,T1,T2,T4>::value,
     "types are not scalar compatible");
  (*mv.data()) = b*(*mv1.data()) + c*(*mv2.data()) +
    d*(*mv3.data()) + e*(*mv4.data());
}

}}//end namespace pressio::ops
#endif  // OPS_EIGEN_OPS_MULTI_VECTOR_UPDATE_HPP_
