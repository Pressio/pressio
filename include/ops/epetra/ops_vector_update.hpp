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

#ifndef OPS_EPETRA_OPS_VECTOR_UPDATE_HPP_
#define OPS_EPETRA_OPS_VECTOR_UPDATE_HPP_

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_epetra<T>::value
  >
update(T & v, const scalar_t a, const T & v1, const scalar_t b)
{
  using int_t = typename ::pressio::containers::details::traits<T>::local_ordinal_t;
  for (int_t i=0; i<v.extentLocal(0); ++i)
    v(i) = a*v(i) + b*v1(i);
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_epetra<T>::value
  >
update(T & v, const T & v1, const scalar_t  b)
{
  using int_t = typename ::pressio::containers::details::traits<T>::local_ordinal_t;
  for (int_t i=0; i<v.extentLocal(0); ++i)
    v(i) = b*v1(i);
}


//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_epetra<T>::value
  >
update(T & v, const scalar_t &a,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c)
{
  using int_t = typename ::pressio::containers::details::traits<T>::local_ordinal_t;
  for (int_t i=0; i<v.extentLocal(0); ++i)
    v(i) = a*v(i) + b*v1(i) + c*v2(i);
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_epetra<T>::value
  >
update(T & v,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c)
{
  using int_t = typename ::pressio::containers::details::traits<T>::local_ordinal_t;
  for (int_t i=0; i<v.extentLocal(0); ++i)
    v(i) = b*v1(i) + c*v2(i);
}


//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_epetra<T>::value
  >
update(T & v, const scalar_t &a,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d)
{
  using int_t = typename ::pressio::containers::details::traits<T>::local_ordinal_t;
  for (int_t i=0; i<v.extentLocal(0); ++i)
    v(i) = a*v(i) + b*v1(i) + c*v2(i) + d*v3(i);
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_epetra<T>::value
  >
update(T & v,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d)
{
  using int_t = typename ::pressio::containers::details::traits<T>::local_ordinal_t;
  for (int_t i=0; i<v.extentLocal(0); ++i)
    v(i) = b*v1(i) + c*v2(i) + d*v3(i);
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_epetra<T>::value
  >
update(T & v, const scalar_t &a,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d,
	  const T & v4, const scalar_t &e)
{
  using int_t = typename ::pressio::containers::details::traits<T>::local_ordinal_t;
  for (int_t i=0; i<v.extentLocal(0); ++i)
    v(i) = a*v(i) + b*v1(i) + c*v2(i) + d*v3(i) + e*v4(i);
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_epetra<T>::value
  >
update(T & v,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d,
	  const T & v4, const scalar_t &e)
{
  using int_t = typename ::pressio::containers::details::traits<T>::local_ordinal_t;
  for (int_t i=0; i<v.extentLocal(0); ++i)
    v(i) = b*v1(i) + c*v2(i) + d*v3(i) + e*v4(i);
}


}}//end namespace pressio::ops
#endif  // OPS_EPETRA_OPS_VECTOR_UPDATE_HPP_
