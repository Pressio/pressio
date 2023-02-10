/*
//@HEADER
// ************************************************************************
//
// ops.hpp
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

#ifndef PRESSIO_OPS_HPP_
#define PRESSIO_OPS_HPP_

#include "./mpl.hpp"
#include "./utils.hpp"
#include "./type_traits.hpp"
#include "./expressions.hpp"

namespace pressio{
struct transpose{};
struct nontranspose{};
enum class Norm{Undefined, L1, L2};

namespace ops{
// matching_extents needs extent but we need to fwd
// because it is used by, e.g., deep_copy
template<class ...> struct matching_extents;
}//end namespace ops
}//end namespace pressio

#include "ops/ops_get_native.hpp"
#include "ops/ops_known_data_type.hpp"
#include "ops/ops_ordinal_type.hpp"

// Eigen
#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "ops/eigen/ops_clone.hpp"
#include "ops/eigen/ops_extent.hpp"
#include "ops/eigen/ops_deep_copy.hpp"
#include "ops/eigen/ops_abs.hpp"
#include "ops/eigen/ops_scale.hpp"
#include "ops/eigen/ops_set_zero.hpp"
#include "ops/eigen/ops_fill.hpp"
#include "ops/eigen/ops_resize.hpp"
#include "ops/eigen/ops_add_to_diagonal.hpp"
#include "ops/eigen/ops_min_max.hpp"
#include "ops/eigen/ops_norms.hpp"
#include "ops/eigen/ops_dot.hpp"
#include "ops/eigen/ops_pow.hpp"
#include "ops/eigen/ops_rank1_update.hpp"
#include "ops/eigen/ops_rank2_update.hpp"
#include "ops/eigen/ops_elementwise_multiply.hpp"
#include "ops/eigen/ops_level2.hpp"
#include "ops/eigen/ops_level3.hpp"
#endif

// Kokkos
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "ops/kokkos/ops_clone.hpp"
#include "ops/kokkos/ops_extent.hpp"
#include "ops/kokkos/ops_deep_copy.hpp"
#include "ops/kokkos/ops_set_zero.hpp"
#include "ops/kokkos/ops_scale.hpp"
#include "ops/kokkos/ops_fill.hpp"
#include "ops/kokkos/ops_resize.hpp"
#include "ops/kokkos/ops_abs.hpp"
#include "ops/kokkos/ops_norms_vector.hpp"
#include "ops/kokkos/ops_dot.hpp"
#include "ops/kokkos/ops_pow.hpp"
#include "ops/kokkos/ops_vector_update.hpp"
#include "ops/kokkos/ops_rank2_update.hpp"
#include "ops/kokkos/ops_elementwise_multiply.hpp"
#include "ops/kokkos/ops_level2.hpp"
#include "ops/kokkos/ops_level3.hpp"
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
// Teuchos
#include "ops/teuchos/ops_extent.hpp"
#include "ops/teuchos/ops_set_zero.hpp"
#include "ops/teuchos/ops_scale.hpp"
#include "ops/teuchos/ops_fill.hpp"
#include "ops/teuchos/ops_level2.hpp"

// Tpetra
#include "ops/tpetra/ops_clone.hpp"
#include "ops/tpetra/ops_extent.hpp"
#include "ops/tpetra/ops_deep_copy.hpp"
#include "ops/tpetra/ops_set_zero.hpp"
#include "ops/tpetra/ops_fill.hpp"
#include "ops/tpetra/ops_abs.hpp"
#include "ops/tpetra/ops_dot.hpp"
#include "ops/tpetra/ops_norms.hpp"
#include "ops/tpetra/ops_pow.hpp"
#include "ops/tpetra/ops_rank1_update.hpp"
#include "ops/tpetra/ops_elementwise_multiply.hpp"
#include "ops/tpetra/ops_multi_vector_update.hpp"
#include "ops/tpetra/ops_level2.hpp"
#include "ops/tpetra/ops_level3.hpp"

// Tpetra block
#include "ops/tpetra_block/ops_clone.hpp"
#include "ops/tpetra_block/ops_extent.hpp"
#include "ops/tpetra_block/ops_deep_copy.hpp"
#include "ops/tpetra_block/ops_set_zero.hpp"
#include "ops/tpetra_block/ops_fill.hpp"
#include "ops/tpetra_block/ops_abs.hpp"
#include "ops/tpetra_block/ops_dot.hpp"
#include "ops/tpetra_block/ops_norms.hpp"
#include "ops/tpetra_block/ops_pow.hpp"
#include "ops/tpetra_block/ops_rank1_update.hpp"
#include "ops/tpetra_block/ops_elementwise_multiply.hpp"
#include "ops/tpetra_block/ops_multi_vector_update.hpp"
#include "ops/tpetra_block/ops_level2.hpp"
#include "ops/tpetra_block/ops_level3.hpp"

// Epetra
#include "ops/epetra/ops_clone.hpp"
#include "ops/epetra/ops_extent.hpp"
#include "ops/epetra/ops_deep_copy.hpp"
#include "ops/epetra/ops_scale.hpp"
#include "ops/epetra/ops_set_zero.hpp"
#include "ops/epetra/ops_fill.hpp"
#include "ops/epetra/ops_abs.hpp"
#include "ops/epetra/ops_dot.hpp"
#include "ops/epetra/ops_min_max.hpp"
#include "ops/epetra/ops_norms.hpp"
#include "ops/epetra/ops_pow.hpp"
#include "ops/epetra/ops_rank1_update.hpp"
#include "ops/epetra/ops_elementwise_multiply.hpp"
#include "ops/epetra/ops_multi_vector_update.hpp"
#include "ops/epetra/ops_level2.hpp"
#include "ops/epetra/ops_level3.hpp"
#endif //PRESSIO_ENABLE_TPL_TRILINOS

// keep this last
#include "ops/ops_matching_extents.hpp"

#endif
