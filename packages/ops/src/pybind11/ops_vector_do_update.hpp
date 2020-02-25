/*
//@HEADER
// ************************************************************************
//
// ops_vector_do_update.hpp
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

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#ifndef OPS_SRC_OPS_PYBIND11_VECTOR_DO_UPDATE_HPP_
#define OPS_SRC_OPS_PYBIND11_VECTOR_DO_UPDATE_HPP_

namespace pressio{ namespace ops{


//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v, scalar_t a, const T & v1, scalar_t b){
  using int_t = typename ::pressio::containers::details::traits<T>::ordinal_t;
  for (int_t i=0; i<v.extent(0); ++i)
    v[i] = a*v[i] + b*v1[i];
}

template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v, scalar_t a, const T & v1, scalar_t b){
  // make sure this is a vector
  assert(v.ndim()==1 && v1.ndim()==1);
  const auto vsz = v.size();
  auto v_proxy = v.mutable_unchecked();
  const auto v1_proxy = v1.unchecked();
  for (std::size_t i=0; i<(std::size_t)vsz; ++i){
    v_proxy(i) = a*v_proxy(i) + b*v1_proxy(i);
  }
}


template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v, const T & v1, const scalar_t b){
  using int_t = typename ::pressio::containers::details::traits<T>::ordinat_t;
  for (int_t i=0; i<v.extent(0); ++i)
    v[i] = b*v1[i];
}

template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v, const T & v1, const scalar_t b){
  // make sure this is a vector
  assert(v.ndim()==1 && v1.ndim()==1);
  assert(v.size()==v1.size());
  const auto vsz = v.size();
  auto v_proxy = v.mutable_unchecked();
  auto v1_proxy = v1.unchecked();
  for (std::size_t i=0; i<(std::size_t)vsz; ++i){
    v_proxy(i) = b*v1_proxy(i);
  }
}

//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){

  using int_t = typename ::pressio::containers::details::traits<T>::ordinal_t;
  for (int_t i=0; i<v.extent(0); ++i)
    v[i] = a*v[i] + b*v1[i] + c*v2[i];
}

template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){

  // make sure this is a vector
  assert(v.ndim()==1 && v1.ndim()==1 && v2.ndim()==1);
  assert(v.size()==v1.size() && v1.size()==v2.size());
  const auto vsz = v.size();
  auto v_proxy   = v.mutable_unchecked();
  auto v1_proxy	 = v1.unchecked();
  auto v2_proxy  = v2.unchecked();
  for (std::size_t i=0; i<(std::size_t)vsz; ++i){
    v_proxy(i) = a*v_proxy(i) + b*v1_proxy(i) + c*v2_proxy(i);
  }
}


template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){

  using int_t = typename ::pressio::containers::details::traits<T>::ordinal_t;
  for (int_t i=0; i<v.extent(0); ++i)
    v[i] = b*v1[i] + c*v2[i];
}

template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){

  // make sure this is a vector
  assert(v.ndim()==1 && v1.ndim()==1 && v2.ndim()==1);
  assert(v.size()==v1.size() && v1.size()==v2.size());
  const auto vsz = v.size();
  auto v_proxy   = v.mutable_unchecked();
  auto v1_proxy	 = v1.unchecked();
  auto v2_proxy  = v2.unchecked();
  for (std::size_t i=0; i<(std::size_t)vsz; ++i){
    v_proxy(i) = b*v1_proxy(i) + c*v2_proxy(i);
  }
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d){

  using int_t = typename ::pressio::containers::details::traits<T>::ordinal_t;
  for (int_t i=0; i<v.extent(0); ++i)
    v[i] = a*v[i] + b*v1[i] + c*v2[i] + d*v3[i];
}

template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d){

  // make sure this is a vector
  assert(v.ndim()==1 && v1.ndim()==1 && v2.ndim()==1 && v3.ndim()==1);
  assert(v.size()==v1.size() && v1.size()==v2.size());
  assert(v2.size()==v3.size());
  const auto vsz = v.size();
  auto v_proxy   = v.mutable_unchecked();
  auto v1_proxy	 = v1.unchecked();
  auto v2_proxy  = v2.unchecked();
  auto v3_proxy	 = v3.unchecked();
  for (std::size_t i=0; i<(std::size_t)vsz; ++i){
    v_proxy(i) = a*v_proxy(i) + b*v1_proxy(i) + c*v2_proxy(i) + d*v3_proxy(i);
  }
}

template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d){

  using int_t = typename ::pressio::containers::details::traits<T>::ordinal_t;
  for (int_t i=0; i<v.extent(0); ++i)
    v[i] = b*v1[i] + c*v2[i] + d*v3[i];
}

template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d){

  // make sure this is a vector
  assert(v.ndim()==1 && v1.ndim()==1 && v2.ndim()==1 && v3.ndim()==1);
  assert(v.size()==v1.size() && v1.size()==v2.size());
  assert(v2.size()==v3.size());
  const auto vsz = v.size();
  auto v_proxy   = v.mutable_unchecked();
  auto v1_proxy	 = v1.unchecked();
  auto v2_proxy  = v2.unchecked();
  auto v3_proxy	 = v3.unchecked();

  for (std::size_t i=0; i<(std::size_t)vsz; ++i){
    v_proxy(i) = b*v1_proxy(i) + c*v2_proxy(i) + d*v3_proxy(i);
  }
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){

  using int_t = typename ::pressio::containers::details::traits<T>::ordinal_t;
  for (int_t i=0; i<v.extent(0); ++i)
    v[i] = a*v[i] + b*v1[i] + c*v2[i] + d*v3[i] + e*v4[i];
}

template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){

  // make sure this is a vector
  assert(v.ndim()==1 && v1.ndim()==1 && v2.ndim()==1 && v3.ndim()==1 && v4.ndim()==1);
  assert(v.size()==v1.size() && v1.size()==v2.size());
  assert(v2.size()==v3.size() && v3.size()==v4.size());
  const auto vsz = v.size();
  auto v_proxy   = v.mutable_unchecked();
  auto v1_proxy	 = v1.unchecked();
  auto v2_proxy  = v2.unchecked();
  auto v3_proxy	 = v3.unchecked();
  auto v4_proxy  = v4.unchecked();
  for (std::size_t i=0; i<(std::size_t)vsz; ++i){
    v_proxy(i) = a*v_proxy(i) + b*v1_proxy(i) + c*v2_proxy(i) + d*v3_proxy(i) + e*v4_proxy(i);
  }
}

template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){

  using int_t = typename ::pressio::containers::details::traits<T>::ordinal_t;
  for (int_t i=0; i<v.extent(0); ++i)
    v[i] = b*v1[i] + c*v2[i] + d*v3[i] + e*v4[i];
}

template<
  typename T, typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){

  // make sure this is a vector
  assert(v.ndim()==1 && v1.ndim()==1 && v2.ndim()==1 && v3.ndim()==1 && v4.ndim()==1);
  assert(v.size()==v1.size() && v1.size()==v2.size());
  assert(v2.size()==v3.size() && v3.size()==v4.size());
  const auto vsz = v.size();
  auto v_proxy   = v.mutable_unchecked();
  auto v1_proxy	 = v1.unchecked();
  auto v2_proxy  = v2.unchecked();
  auto v3_proxy	 = v3.unchecked();
  auto v4_proxy  = v4.unchecked();
  for (std::size_t i=0; i<(std::size_t)vsz; ++i){
    v_proxy(i) = b*v1_proxy(i) + c*v2_proxy(i) + d*v3_proxy(i) + e*v4_proxy(i);
  }
}

}}//end namespace pressio::ops
#endif
#endif
