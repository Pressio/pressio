/*
//@HEADER
// ************************************************************************
//
// containers_vector_subtract_operator_distributed.hpp
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

#ifndef CONTAINERS_VECTOR_VECTOR_SUBTRACT_OPERATOR_DISTRIBUTED_HPP_
#define CONTAINERS_VECTOR_VECTOR_SUBTRACT_OPERATOR_DISTRIBUTED_HPP_

#include "../containers_vector_traits.hpp"
#include "../containers_vector_meta.hpp"
#include "containers_vector_distributed_binary_expression_templates.hpp"
#include "../../containers_expression_templates_operators.hpp"

namespace pressio{ namespace containers{


// T1: expre, T2: vector:
// example: a*3 - b
template <typename T1,
	  typename T2,
	  ::pressio::mpl::enable_if_t<
  meta::is_admissible_vec_for_dist_expression<T2>::value
	    > * = nullptr>
auto operator-(const T1 & u, const T2 & v)
  -> decltype(
      containers::exprtemplates::DistributedVectorBinaryExp<
      containers::exprtemplates::subtract_,
      T1, T2,
      typename containers::details::traits<T2>::scalar_t,
      typename containers::details::traits<T2>::local_ordinal_t>(u, v)
    )
{
  using sc_t = typename containers::details::traits<T2>::scalar_t;
  using LO_t = typename containers::details::traits<T2>::local_ordinal_t;

  return containers::exprtemplates::DistributedVectorBinaryExp<
    containers::exprtemplates::subtract_,T1, T2, sc_t, LO_t>(u, v);
}


//-----------------------------------------------------
// T1: expre, T2: expr:
// example: a*3 - b*21
template <typename T1,
	  typename T2,
	  ::pressio::mpl::enable_if_t<
  exprtemplates::is_distributed_vector_expression<T1>::value &&
  exprtemplates::is_distributed_vector_expression<T2>::value
	    > * = nullptr>
auto operator-(const T1 & u, const T2 & v)
  -> decltype(
      containers::exprtemplates::DistributedVectorBinaryExp<
      containers::exprtemplates::subtract_,
      T1, T2,
      typename T2::sc_type,
      typename T2::LO_type>(u, v)
    )
{
  using sc_t = typename T2::sc_type;
  using LO_t = typename T2::LO_type;

  return containers::exprtemplates::DistributedVectorBinaryExp<
    containers::exprtemplates::subtract_,T1, T2, sc_t, LO_t>(u, v);
}


//-----------------------------------------------------
// T1: vector, T2: expr:
// example: a - b*21
template <typename T1,
	  typename T2,
	  ::pressio::mpl::enable_if_t<
  meta::is_admissible_vec_for_dist_expression<T1>::value &&
  exprtemplates::is_distributed_vector_expression<T2>::value
	    > * = nullptr>
auto operator-(const T1 & u, const T2 & v)
  -> decltype(
      containers::exprtemplates::DistributedVectorBinaryExp<
      containers::exprtemplates::subtract_,
      T1, T2,
      typename T2::sc_type,
      typename T2::LO_type>(u, v)
    )
{
  using sc_t = typename T2::sc_type;
  using LO_t = typename T2::LO_type;

  return containers::exprtemplates::DistributedVectorBinaryExp<
    containers::exprtemplates::subtract_,T1, T2, sc_t, LO_t>(u, v);
}


}}//end namespace pressio::containers
#endif
