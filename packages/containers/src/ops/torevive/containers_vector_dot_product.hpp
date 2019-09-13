/*
//@HEADER
// ************************************************************************
//
// containers_vector_dot_product.hpp
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

#ifndef CONTAINERS_VECTOR_OPERATIONS_VECTOR_DOT_PRODUCT_HPP_
#define CONTAINERS_VECTOR_OPERATIONS_VECTOR_DOT_PRODUCT_HPP_

#include "../../meta/containers_vector_meta.hpp"
#include "../concrete/containers_vector_sharedmem_eigen.hpp"

namespace pressio{
namespace containers{
namespace vec_ops{

template <typename vec_a_type,
	  typename vec_b_type,
	  typename res_t,
	  ::pressio::mpl::enable_if_t<
	   details::traits<vec_a_type>::is_vector &&
	   details::traits<vec_a_type>::isEigen && 
	   details::traits<vec_b_type>::is_vector &&
	   details::traits<vec_b_type>::isEigen &&
	   std::is_same<
	     typename details::traits<vec_a_type>::scalar_t,
	     typename details::traits<vec_b_type>::scalar_t
	     >::value &&
	   std::is_same<
	     typename details::traits<vec_a_type>::scalar_t,
	     res_t>::value
	   > * = nullptr
	  >
void dot(const vec_a_type & vecA,
	 const vec_b_type & vecB,
	 res_t & result)
{
  assert(vecA.size() == vecB.size());
  result = vecA.data()->dot(*vecB.data());
}


} // end namespace vec_ops
} // end namespace containers
}//end namespace pressio
#endif
