/*
//@HEADER
// ************************************************************************
//
// ops_pow.hpp
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

#ifndef OPS_TPETRA_BLOCK_OPS_POW_HPP_
#define OPS_TPETRA_BLOCK_OPS_POW_HPP_

namespace pressio{ namespace ops{

// y = |x|^exponent, expo>0
template <typename T1, typename T2>
::pressio::mpl::enable_if_t<
  (   ::pressio::is_vector_tpetra_block<T1>::value
   || ::pressio::is_expression_column_acting_on_tpetra_block<T1>::value)
  && (::pressio::is_vector_tpetra_block<T2>::value
   || ::pressio::is_expression_column_acting_on_tpetra_block<T2>::value)
  >
abs_pow(T1 & y,
	const T2 & x,
	const typename ::pressio::Traits<T1>::scalar_type & exponent)
{

  using sc_t = typename ::pressio::Traits<T1>::scalar_type;

  assert(extent(x,0) == extent(y,0));
  assert(exponent > ::pressio::utils::Constants<sc_t>::zero());
  if (exponent < ::pressio::utils::Constants<sc_t>::zero()){
    throw std::runtime_error("this overload is only for exponent > 0");
  }

  auto x_tp = impl::get_underlying_tpetra_object(x);
  auto y_tp = impl::get_underlying_tpetra_object(y);
  ::pressio::ops::abs_pow(y_tp, x_tp, exponent);
}

// y = |x|^exponent, expo<0
template <typename T1, typename T2>
::pressio::mpl::enable_if_t<
  (::pressio::is_vector_tpetra_block<T1>::value
   || ::pressio::is_expression_column_acting_on_tpetra_block<T1>::value)
  && (::pressio::is_vector_tpetra_block<T2>::value
   || ::pressio::is_expression_column_acting_on_tpetra_block<T2>::value)
  >
abs_pow(T1 & y,
	const T2 & x,
	const typename ::pressio::Traits<T1>::scalar_type & exponent,
	const typename ::pressio::Traits<T1>::scalar_type & eps)
{

  using sc_t = typename ::pressio::Traits<T1>::scalar_type;

  assert(extent(x,0) == extent(y,0));
  assert(exponent < ::pressio::utils::Constants<sc_t>::zero());
  if (exponent > ::pressio::utils::Constants<sc_t>::zero()){
    throw std::runtime_error("this overload is only for exponent < 0");
  }

  auto x_tp = impl::get_underlying_tpetra_object(x);
  auto y_tp = impl::get_underlying_tpetra_object(y);
  ::pressio::ops::abs_pow(y_tp, x_tp, exponent, eps);
}

template <typename T>
::pressio::mpl::enable_if_t<
  ::pressio::is_vector_tpetra_block<T>::value
  || ::pressio::is_expression_column_acting_on_tpetra_block<T>::value
  >
pow(T & x,
    const typename ::pressio::Traits<T>::scalar_type & exponent)
{
  auto x_tp = impl::get_underlying_tpetra_object(x);
  ::pressio::ops::pow(x_tp, exponent);
}

}}//end namespace pressio::ops
#endif  // OPS_TPETRA_BLOCK_OPS_POW_HPP_
