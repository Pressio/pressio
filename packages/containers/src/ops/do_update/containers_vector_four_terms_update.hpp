/*
//@HEADER
// ************************************************************************
//
// containers_vector_four_terms_update.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef CONTAINERS_CONTAINER_OPS_VECTOR_FOUR_TERMS_UPDATE_HPP_
#define CONTAINERS_CONTAINER_OPS_VECTOR_FOUR_TERMS_UPDATE_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"

#ifdef HAVE_KOKKOS
#include "containers_vector_do_update_kokkos_functors.hpp"
#endif

#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------

namespace pressio{ namespace containers{ namespace ops{


//---------------------------------------------------------------
// enable for vectors supporting expression templates
//---------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){
  v = a*v + b*v1 + c*v2 + d*v3 + e*v4;
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){
  v = b*v1 + c*v2 + d*v3 + e*v4;
}



//--------------------------------------------------------------------------
// enable for pybind11::array_t
//--------------------------------------------------------------------------
#ifdef HAVE_PYBIND11
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size()
      and vsz != v3.size() and vsz != v4.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = a*v.at(i) + b*v1.at(i) + c*v2.at(i) + d*v3.at(i) + e*v4.at(i);
  }
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size()
      and vsz != v3.size() and vsz != v4.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = b*v1.at(i) + c*v2.at(i) + d*v3.at(i) + e*v4.at(i);
  }
}
#endif



//-----------------------------------------------------------------------------
// enable for tpetra and tpetra block vectors
//-----------------------------------------------------------------------------
#ifdef HAVE_TRILINOS
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value or
    ::pressio::containers::meta::is_vector_wrapper_tpetra_block<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e)
{
  constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();

  v.data()->update(b, *v1.data(), a); // v = a*v + b*v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
  v.data()->update(e, *v4.data(), one); // add e*v4
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value or
    ::pressio::containers::meta::is_vector_wrapper_tpetra_block<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e)
{
  constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  v.data()->update(b, *v1.data(), zero); // v = b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
  v.data()->update(e, *v4.data(), one); // add e*v4
}
#endif


//--------------------------------------------------------------------------
// enable for Kokkos wrappers
//--------------------------------------------------------------------------
#ifdef HAVE_KOKKOS
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::containers::ops::impl::DoUpdateFourTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), *v4.data(), a, b, c, d, e);
  Kokkos::parallel_for(v.size(), F);
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::containers::ops::impl::DoUpdateFourTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), *v4.data(), b, c, d, e);
  Kokkos::parallel_for(v.size(), F);
}
#endif

}}}//end namespace pressio::containers::ops
#endif
