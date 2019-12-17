/*
//@HEADER
// ************************************************************************
//
// containers_subspan.hpp
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

#ifndef CONTAINERS_SUBSPAN_HPP_
#define CONTAINERS_SUBSPAN_HPP_

#include "containers_expression_base.hpp"
#include "../matrix/containers_matrix_meta.hpp"
#include "containers_matrix_subspan_expression.hpp"

namespace pressio{ namespace containers{

template <typename T>
mpl::enable_if_t<
  meta::is_matrix_wrapper<T>::value,
  typename details::traits<T>::subspan_const_ret_t
  >
subspan(const T & obj,
	const typename details::traits<T>::subspan_const_ret_t::interval_t & rowRangeIn,
	const typename details::traits<T>::subspan_const_ret_t::interval_t & colRangeIn)
{
  using return_t = typename details::traits<T>::subspan_const_ret_t;
  return return_t(obj, rowRangeIn, colRangeIn);
}

template <typename T>
mpl::enable_if_t<
  meta::is_matrix_wrapper<T>::value,
  typename details::traits<T>::subspan_ret_t
  >
subspan(T & obj,
	const typename details::traits<T>::subspan_ret_t::interval_t & rowRangeIn,
	const typename details::traits<T>::subspan_ret_t::interval_t & colRangeIn)
{
  using return_t = typename details::traits<T>::subspan_ret_t;
  return return_t(obj, rowRangeIn, colRangeIn);
}


}} //end namespace pressio::containers

#endif
