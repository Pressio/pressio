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

#ifndef OPS_EPETRA_OPS_RANK1_UPDATE_HPP_
#define OPS_EPETRA_OPS_RANK1_UPDATE_HPP_

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<typename T, typename T1,
         typename a_Type, typename b_Type>
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  // TPL/container specific
  && ::pressio::is_vector_epetra<T>::value
  && ::pressio::is_vector_epetra<T1>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<a_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v,         const a_Type a,
       const T1 & v1, const b_Type b)
{
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  const auto zero = static_cast<scalar_t>(0);
  const scalar_t a_(a);
  const scalar_t b_(b);

  if (b_ == zero) {
    ::pressio::ops::scale(v, a_);
    return;
  }

  for (int i=0; i<v.MyLength(); ++i) {
    if (a_ == zero) {
      v[i] = b_*v1[i];
    } else {
      v[i] = a_*v[i] + b_*v1[i];
    }
  }
}

//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<
  class T, class T1, class T2,
  class a_Type, class b_Type, class c_Type
  >
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  // TPL/container specific
  && ::pressio::is_vector_epetra<T>::value
  && ::pressio::is_vector_epetra<T1>::value
  && ::pressio::is_vector_epetra<T2>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1, T2>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<a_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<c_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v,         const a_Type &a,
       const T1 & v1, const b_Type &b,
       const T2 & v2, const c_Type &c)
{
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v2, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  const auto zero = static_cast<scalar_t>(0);
  const scalar_t a_(a);
  const scalar_t b_(b);
  const scalar_t c_(c);

  if (b_ == zero) {
    ::pressio::ops::update(v, a, v2, c);
  } else if (c_ == zero) {
    ::pressio::ops::update(v, a, v1, b);
  } else {
    for (int i=0; i<v.MyLength(); ++i) {
      if (a_ == zero) {
        v[i] = b_*v1[i] + c_*v2[i];
      } else {
        v[i] = a_*v[i] + b_*v1[i] + c_*v2[i];
      }
    }
  }
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<
  class T, class T1, class T2, class T3,
  class a_Type, class b_Type, class c_Type, class d_Type
  >
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  && ::pressio::Traits<T3>::rank == 1
  // TPL/container specific
  && ::pressio::is_vector_epetra<T>::value
  && ::pressio::is_vector_epetra<T1>::value
  && ::pressio::is_vector_epetra<T2>::value
  && ::pressio::is_vector_epetra<T3>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1, T2, T3>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
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
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v2, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v3, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  const auto zero = static_cast<scalar_t>(0);
  const scalar_t a_(a);
  const scalar_t b_(b);
  const scalar_t c_(c);
  const scalar_t d_(d);

  if (b_ == zero) {
    ::pressio::ops::update(v, a, v2, c, v3, d);
  } else if (c_ == zero) {
    ::pressio::ops::update(v, a, v1, b, v3, d);
  } else if (d_ == zero) {
    ::pressio::ops::update(v, a, v1, b, v2, c);
  } else {
    for (int i=0; i<v.MyLength(); ++i) {
      if (a_ == zero) {
        v[i] = b_*v1[i] + c_*v2[i] + d_*v3[i];
      } else {
        v[i] = a_*v[i] + b_*v1[i] + c_*v2[i] + d_*v3[i];
      }
    }
  }
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<
  class T, class T1, class T2, class T3, class T4,
  class a_Type, class b_Type, class c_Type, class d_Type, class e_Type
  >
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  && ::pressio::Traits<T3>::rank == 1
  && ::pressio::Traits<T4>::rank == 1
  // TPL/container specific
  && ::pressio::is_vector_epetra<T>::value
  && ::pressio::is_vector_epetra<T1>::value
  && ::pressio::is_vector_epetra<T2>::value
  && ::pressio::is_vector_epetra<T3>::value
  && ::pressio::is_vector_epetra<T4>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1, T2, T3, T4>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
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
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v2, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v3, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v4, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  const auto zero = static_cast<scalar_t>(0);
  const scalar_t a_(a);
  const scalar_t b_(b);
  const scalar_t c_(c);
  const scalar_t d_(d);
  const scalar_t e_(e);

  if (b_ == zero) {
    ::pressio::ops::update(v, a, v2, c, v3, d, v4, e);
  } else if (c_ == zero) {
    ::pressio::ops::update(v, a, v1, b, v3, d, v4, e);
  } else if (d_ == zero) {
    ::pressio::ops::update(v, a, v1, b, v2, c, v4, e);
  } else if (e_ == zero) {
    ::pressio::ops::update(v, a, v1, b, v2, c, v3, d);
  } else {
    for (int i = 0; i < v.MyLength(); ++i) {
      if (a_ == zero) {
        v[i] = b_*v1[i] + c_*v2[i] + d_*v3[i] + e_*v4[i];
      } else {
        v[i] = a_*v[i] + b_*v1[i] + c_*v2[i] + d_*v3[i] + e_*v4[i];
      }
    }
  }
}

}}//end namespace pressio::ops
#endif  // OPS_EPETRA_OPS_RANK1_UPDATE_HPP_
