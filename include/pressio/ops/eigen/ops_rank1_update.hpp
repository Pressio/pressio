/*
//@HEADER
// ************************************************************************
//
// ops_rank1_update.hpp
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

#ifndef OPS_EIGEN_OPS_RANK1_UPDATE_HPP_
#define OPS_EIGEN_OPS_RANK1_UPDATE_HPP_

namespace pressio{ namespace ops{

/*
   below we constrain in all cases via is_convertible
   because the implementations are using Eigen native operations
   which are based on expressions and require
   coefficients to be convertible to scalar types of the vectors operands
 */

//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<class T, class T1, class a_Type, class b_Type>
::pressio::mpl::enable_if_t<
     ::pressio::all_have_traits_and_same_scalar<T, T1>::value
  && ::pressio::package_identifier<T>::value == PackageIdentifier::Eigen
  && ::pressio::package_identifier<T1>::value == PackageIdentifier::Eigen
  && ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && std::is_convertible<a_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v,         const a_Type & a,
       const T1 & v1, const b_Type & b)
{

  auto & v_n = impl::get_native(v);
  const auto & v_n1 = impl::get_native(v1);
  v_n = a*v_n + b*v_n1;
}

//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<
  class T, class T1, class T2,
  class a_Type, class b_Type, class c_Type
  >
::pressio::mpl::enable_if_t<
     ::pressio::all_have_traits_and_same_scalar<T, T1, T2>::value
  && ::pressio::package_identifier<T>::value == PackageIdentifier::Eigen
  && ::pressio::package_identifier<T1>::value == PackageIdentifier::Eigen
  && ::pressio::package_identifier<T2>::value == PackageIdentifier::Eigen
  && ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  && std::is_convertible<a_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<c_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v,         const a_Type &a,
       const T1 & v1, const b_Type &b,
       const T2 & v2, const c_Type &c)
{

  auto & v_n = impl::get_native(v);
  const auto & v_n1 = impl::get_native(v1);
  const auto & v_n2 = impl::get_native(v2);
  v_n = a*v_n + b*v_n1 + c*v_n2;
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<
  class T, class T1, class T2, class T3,
  class a_Type, class b_Type, class c_Type, class d_Type
  >
::pressio::mpl::enable_if_t<
     ::pressio::all_have_traits_and_same_scalar<T, T1, T2, T3>::value
  && ::pressio::package_identifier<T>::value == PackageIdentifier::Eigen
  && ::pressio::package_identifier<T1>::value == PackageIdentifier::Eigen
  && ::pressio::package_identifier<T2>::value == PackageIdentifier::Eigen
  && ::pressio::package_identifier<T3>::value == PackageIdentifier::Eigen
  && ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  && ::pressio::Traits<T3>::rank == 1
  && std::is_convertible<a_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<c_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<d_Type, typename ::pressio::Traits<T>::scalar_type>::value

  >
update(T & v,         const a_Type &a,
       const T1 & v1, const b_Type &b,
       const T2 & v2, const c_Type &c,
       const T3 & v3, const d_Type &d)
{

  auto & v_n = impl::get_native(v);
  const auto & v_n1 = impl::get_native(v1);
  const auto & v_n2 = impl::get_native(v2);
  const auto & v_n3 = impl::get_native(v3);
  v_n = a*v_n + b*v_n1 + c*v_n2 + d*v_n3;
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<
  class T, class T1, class T2, class T3, class T4,
  class a_Type, class b_Type, class c_Type, class d_Type, class e_Type
  >
::pressio::mpl::enable_if_t<
     ::pressio::all_have_traits_and_same_scalar<T, T1, T2, T3, T4>::value
  && ::pressio::package_identifier<T>::value == PackageIdentifier::Eigen
  && ::pressio::package_identifier<T1>::value == PackageIdentifier::Eigen
  && ::pressio::package_identifier<T2>::value == PackageIdentifier::Eigen
  && ::pressio::package_identifier<T3>::value == PackageIdentifier::Eigen
  && ::pressio::package_identifier<T4>::value == PackageIdentifier::Eigen
  && ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  && ::pressio::Traits<T3>::rank == 1
  && ::pressio::Traits<T4>::rank == 1
  && std::is_convertible<a_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<c_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<d_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<e_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v,         const a_Type &a,
       const T1 & v1, const b_Type &b,
       const T2 & v2, const c_Type &c,
       const T3 & v3, const d_Type &d,
       const T4 & v4, const e_Type &e)
{

  auto & v_n = impl::get_native(v);
  const auto & v_n1 = impl::get_native(v1);
  const auto & v_n2 = impl::get_native(v2);
  const auto & v_n3 = impl::get_native(v3);
  const auto & v_n4 = impl::get_native(v4);
  v_n = a*v_n + b*v_n1 + c*v_n2 + d*v_n3 + e*v_n4;
}

}}//end namespace pressio::ops
#endif  // OPS_EIGEN_OPS_RANK1_UPDATE_HPP_
