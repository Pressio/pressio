/*
//@HEADER
// ************************************************************************
//
// ops_vector_update.hpp
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

#ifndef OPS_EIGEN_OPS_VECTOR_UPDATE_HPP_
#define OPS_EIGEN_OPS_VECTOR_UPDATE_HPP_

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<typename T, typename T1, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::traits<T>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T1>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T>::rank == 1 and
  ::pressio::traits<T1>::rank == 1 
  >
update(T & v, const scalar_t a, const T1 & v1, const scalar_t b)
{
  static_assert
    (::pressio::are_scalar_compatible<T,T1>::value,
     "vector types T and T1 in ops/src/eigen/ops_vector_update.hpp are not scalar compatible");
  using ord_t = typename ::pressio::traits<T1>::ordinal_type;
  for (ord_t i=0; i< ::pressio::ops::extent(v, 0); ++i){
    v(i) = a*v(i) + b*v1(i);
  }
}

template<typename T, typename T1, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::traits<T>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T1>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T>::rank == 1 and
  ::pressio::traits<T1>::rank == 1
  >
update(T & v, const T1 & v1, const scalar_t  b)
{
  static_assert
    (::pressio::are_scalar_compatible<T,T1>::value,
     "vector types T and T1 in ops/src/eigen/ops_vector_update.hpp are not scalar compatible");

  using ord_t = typename ::pressio::traits<T1>::ordinal_type;
  for (ord_t i=0; i< ::pressio::ops::extent(v, 0); ++i){
    v(i) = b*v1(i);
  }
}

//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<typename T, typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::traits<T>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T1>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T2>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T>::rank == 1 and
  ::pressio::traits<T1>::rank == 1 and 
  ::pressio::traits<T2>::rank == 1 
  >
update(T & v, const scalar_t &a,
       const T1 & v1, const scalar_t &b,
       const T2 & v2, const scalar_t &c)
{
  static_assert
    (::pressio::are_scalar_compatible<T,T1,T2>::value,
     "vector types T,T1,T2 in ops/src/eigen/ops_vector_update.hpp are not scalar compatible");
  using ord_t = typename ::pressio::traits<T1>::ordinal_type;

  for (ord_t i=0; i< ::pressio::ops::extent(v, 0); ++i){
    v(i) = a*v(i) + b*v1(i) + c*v2(i);
  }
}

template<typename T, typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::traits<T>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T1>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T2>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T>::rank == 1 and
  ::pressio::traits<T1>::rank == 1 and 
  ::pressio::traits<T2>::rank == 1 
  >
update(T & v,
       const T1 & v1, const scalar_t &b,
       const T2 & v2, const scalar_t &c)
{
  static_assert
    (::pressio::are_scalar_compatible<T,T1,T2>::value,
     "vector types T,T1,T2 in ops/src/eigen/ops_vector_update.hpp are not scalar compatible");
  using ord_t = typename ::pressio::traits<T1>::ordinal_type;

  for (ord_t i=0; i< ::pressio::ops::extent(v, 0); ++i){
    v(i) = b*v1(i) + c*v2(i);
  }
}


//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<typename T,
         typename T1,
         typename T2,
         typename T3,
         typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::traits<T>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T1>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T2>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T3>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T>::rank == 1 and
  ::pressio::traits<T1>::rank == 1 and 
  ::pressio::traits<T2>::rank == 1 and 
  ::pressio::traits<T3>::rank == 1 
  >
update(T  & v, const scalar_t &a,
       const T1 & v1, const scalar_t &b,
       const T2 & v2, const scalar_t &c,
       const T3 & v3, const scalar_t &d)
{
  static_assert
    (::pressio::are_scalar_compatible<T,T1,T2,T3>::value,
     "vector types T,T1,T2,T3 in ops/src/eigen/ops_vector_update.hpp are not scalar compatible");
  using ord_t = typename ::pressio::traits<T1>::ordinal_type;

  for (ord_t i=0; i< ::pressio::ops::extent(v, 0); ++i){
    v(i) = a*v(i) + b*v1(i) + c*v2(i) + d*v3(i);
  }
}

template<typename T,
         typename T1,
         typename T2,
         typename T3,
         typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::traits<T>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T1>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T2>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T3>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T>::rank == 1 and
  ::pressio::traits<T1>::rank == 1 and 
  ::pressio::traits<T2>::rank == 1 and 
  ::pressio::traits<T3>::rank == 1 
  >
update(T & v,
       const T1 & v1, const scalar_t &b,
       const T2 & v2, const scalar_t &c,
       const T3 & v3, const scalar_t &d)
{
  static_assert
    (::pressio::are_scalar_compatible<T,T1,T2,T3>::value,
     "vector types T,T1,T2,T3 in ops/src/eigen/ops_vector_update.hpp are not scalar compatible");
  using ord_t = typename ::pressio::traits<T1>::ordinal_type;

  for (ord_t i=0; i< ::pressio::ops::extent(v, 0); ++i){
    v(i) = b*v1(i) + c*v2(i) + d*v3(i);
  }
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template< typename T,
          typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::traits<T>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T1>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T2>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T3>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T4>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T>::rank == 1 and
  ::pressio::traits<T1>::rank == 1 and 
  ::pressio::traits<T2>::rank == 1 and 
  ::pressio::traits<T3>::rank == 1 and 
  ::pressio::traits<T4>::rank == 1 
  >
update(T & v, const scalar_t &a,
       const T1 & v1, const scalar_t &b,
       const T2 & v2, const scalar_t &c,
       const T3 & v3, const scalar_t &d,
       const T4 & v4, const scalar_t &e)
{
  static_assert
    (::pressio::are_scalar_compatible<T,T1,T2,T3,T4>::value,
     "vector types T,T1,T2,T3,T4 in ops/src/eigen/ops_vector_update.hpp are not scalar compatible");
  using ord_t = typename ::pressio::traits<T1>::ordinal_type;

  for (ord_t i=0; i< ::pressio::ops::extent(v, 0); ++i){
    v(i) = a*v(i) + b*v1(i) + c*v2(i) + d*v3(i) + e*v4(i);
  }
}

template<typename T,
         typename T1,
         typename T2,
         typename T3,
         typename T4,
         typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::traits<T>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T1>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T2>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T3>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T4>::package_identifier == PackageIdentifier::Eigen and
  ::pressio::traits<T>::rank == 1 and
  ::pressio::traits<T1>::rank == 1 and 
  ::pressio::traits<T2>::rank == 1 and 
  ::pressio::traits<T3>::rank == 1 and 
  ::pressio::traits<T4>::rank == 1 
  >
update(T & v,
       const T1 & v1, const scalar_t &b,
       const T2 & v2, const scalar_t &c,
       const T3 & v3, const scalar_t &d,
       const T4 & v4, const scalar_t &e)
{
  static_assert
    (::pressio::are_scalar_compatible<T,T1,T2,T3,T4>::value,
     "vector types T,T1,T2,T3,T4 in ops/src/eigen/ops_vector_update.hpp are not scalar compatible");
  using ord_t = typename ::pressio::traits<T1>::ordinal_type;

  for (ord_t i=0; i< ::pressio::ops::extent(v, 0); ++i){
    v(i) = b*v1(i) + c*v2(i) + d*v3(i) + e*v4(i);
  }
}


}}//end namespace pressio::ops
#endif  // OPS_EIGEN_OPS_VECTOR_UPDATE_HPP_
