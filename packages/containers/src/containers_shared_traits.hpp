/*
//@HEADER
// ************************************************************************
//
// containers_shared_traits.hpp
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

#ifndef CONTAINERS_SHARED_TRAITS_HPP_
#define CONTAINERS_SHARED_TRAITS_HPP_

#include "containers_wrapped_types_enum.hpp"

namespace pressio{ namespace containers{ namespace details {

template<
  typename container_T,
  typename wrapped_T,
  bool is_vector_t,
  bool is_matrix_t,
  bool is_multi_vector_t,
  WrappedPackageIdentifier wpid,
  bool is_shared_mem_t,
  bool is_static_t
  >
struct containers_shared_traits{

  using wrapped_t = wrapped_T;
  using derived_t = container_T;

  static constexpr WrappedPackageIdentifier wrapped_package_identifier = wpid;

  static constexpr bool is_vector = is_vector_t;
  static constexpr bool is_matrix = is_matrix_t;
  static constexpr bool is_multi_vector = is_multi_vector_t;
  static constexpr bool is_shared_mem = is_shared_mem_t;
  static constexpr bool is_static = is_static_t;
  static constexpr bool is_dynamic = !is_static_t;

  // by default, any container is not admissible to expr templates
  // the ones that are, will overwrite this
  static constexpr bool is_admissible_for_expression_templates = false;
};

/// common traits of matrices
template<bool is_sparse_t>
struct matrix_shared_traits{
  static constexpr bool is_sparse = is_sparse_t;
  static constexpr bool is_dense = !is_sparse_t;
};

}}} // end namespace pressio::containers::details
#endif
