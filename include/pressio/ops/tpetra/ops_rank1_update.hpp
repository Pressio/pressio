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

#ifndef OPS_TPETRA_OPS_RANK1_UPDATE_HPP_
#define OPS_TPETRA_OPS_RANK1_UPDATE_HPP_

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<class T, class T1, class a_Type, class b_Type>
::pressio::mpl::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  // TPL/container specific
  && ::pressio::is_vector_tpetra<T>::value
  && ::pressio::is_vector_tpetra<T1>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<a_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v,         const a_Type & a,
       const T1 & v1, const b_Type & b)
{
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  sc_t a_{a};
  sc_t b_{b};
  v.update(b_, v1, a_);
}

//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<
  class T, class T1, class T2,
  class a_Type, class b_Type, class c_Type
  >
::pressio::mpl::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  // TPL/container specific
  && ::pressio::is_vector_tpetra<T>::value
  && ::pressio::is_vector_tpetra<T1>::value
  && ::pressio::is_vector_tpetra<T2>::value
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

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  constexpr auto one  = ::pressio::utils::Constants<sc_t>::one();
  sc_t a_{a};
  sc_t b_{b};
  sc_t c_{c};
  v.update(b_, v1, a_); // v = v*a + b * v1
  v.update(c_, v2, one); // add c*v2
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
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  && ::pressio::Traits<T3>::rank == 1
  // TPL/container specific
  && ::pressio::is_vector_tpetra<T>::value
  && ::pressio::is_vector_tpetra<T1>::value
  && ::pressio::is_vector_tpetra<T2>::value
  && ::pressio::is_vector_tpetra<T3>::value
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

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  constexpr auto one  = ::pressio::utils::Constants<sc_t>::one();
  sc_t a_{a};
  sc_t b_{b};
  sc_t c_{c};
  sc_t d_{d};
  v.update(b_, v1, a_); // v = a*v + b*v1
  v.update(c_, v2, one); // add c*v2
  v.update(d_, v3, one); // add d*v3
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
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  && ::pressio::Traits<T3>::rank == 1
  && ::pressio::Traits<T4>::rank == 1
  // TPL/container specific
  && ::pressio::is_vector_tpetra<T>::value
  && ::pressio::is_vector_tpetra<T1>::value
  && ::pressio::is_vector_tpetra<T2>::value
  && ::pressio::is_vector_tpetra<T3>::value
  && ::pressio::is_vector_tpetra<T4>::value
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

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  constexpr auto one  = ::pressio::utils::Constants<sc_t>::one();
  sc_t a_{a};
  sc_t b_{b};
  sc_t c_{c};
  sc_t d_{d};
  sc_t e_{e};
  v.update(b_, v1, a_); // v = a*v + b*v1
  v.update(c_, v2, one); // add c*v2
  v.update(d_, v3, one); // add d*v3
  v.update(e_, v4, one); // add e*v4
}

}}//end namespace pressio::ops
#endif  // OPS_TPETRA_OPS_RANK1_UPDATE_HPP_
