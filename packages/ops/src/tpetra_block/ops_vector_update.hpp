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

#ifndef OPS_TPETRA_BLOCK_OPS_VECTOR_UPDATE_HPP_
#define OPS_TPETRA_BLOCK_OPS_VECTOR_UPDATE_HPP_

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<T>::value
  >
update(T & v, const scalar_t a,
	  const T & v1, const scalar_t b)
{
  v.data()->update(b, *v1.data(), a); // v = a*v + b * v1
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<T>::value
  >
update(T & v, const T & v1, const scalar_t b)
{
  constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
  v.data()->update(b, *v1.data(), zero); // v = b * v1
}


//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<T>::value
  >
update(T & v, const scalar_t &a,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c)
{
  constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
  v.data()->update(b, *v1.data(), a); // v = v + b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<T>::value
  >
update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c)
{
  constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
  constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();

  v.data()->update(b, *v1.data(), zero); // v = b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
}


//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<T>::value
  >
update(T & v, const scalar_t &a,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d)
{
  constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();

  v.data()->update(b, *v1.data(), a); // v = a*v + b*v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<T>::value
  >
update(T & v,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d)
{
  constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
  constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();

  v.data()->update(b, *v1.data(), zero); // v = b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
}


//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<T>::value
  >
update(T & v, const scalar_t &a,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d,
	  const T & v4, const scalar_t &e)
{
  constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();

  v.data()->update(b, *v1.data(), a); // v = a*v + b*v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
  v.data()->update(e, *v4.data(), one); // add e*v4
}

template<typename T, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<T>::value
  >
update(T & v,
	  const T & v1, const scalar_t &b,
	  const T & v2, const scalar_t &c,
	  const T & v3, const scalar_t &d,
	  const T & v4, const scalar_t &e)
{
  constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
  constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();

  v.data()->update(b, *v1.data(), zero); // v = b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
  v.data()->update(e, *v4.data(), one); // add e*v4
}


}}//end namespace pressio::ops
#endif  // OPS_TPETRA_BLOCK_OPS_VECTOR_UPDATE_HPP_
