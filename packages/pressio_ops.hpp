/*
//@HEADER
// ************************************************************************
//
// pressio_ops.hpp
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

#include "pressio_mpl.hpp"
#include "pressio_utils.hpp"
#include "pressio_containers.hpp"

#include "ops/src/ops_fwd.hpp"

// predicates
#include "ops/src/predicates/ops_has_method_deep_copy.hpp"
#include "ops/src/predicates/ops_has_method_set_zero.hpp"
#include "ops/src/predicates/ops_has_method_axpy.hpp"
#include "ops/src/predicates/ops_has_method_norm1.hpp"
#include "ops/src/predicates/ops_has_method_norm2.hpp"
#include "ops/src/predicates/ops_has_void_method_product_mat_mat.hpp"
#include "ops/src/predicates/ops_has_nonvoid_method_product_mat_mat.hpp"
#include "ops/src/predicates/ops_has_void_method_product_mat_vec.hpp"
#include "ops/src/predicates/ops_has_method_add_to_diagonal.hpp"
#include "ops/src/predicates/ops_has_method_scale.hpp"
#include "ops/src/predicates/ops_has_method_update_one_term.hpp"
#include "ops/src/predicates/ops_has_method_update_two_terms.hpp"
#include "ops/src/predicates/ops_has_method_update_three_terms.hpp"
#include "ops/src/predicates/ops_has_method_update_four_terms.hpp"

// ops_is_object_pybind: not within preproc direc because we need to use it
// even when pybind is disabled, in which case it will always be false
#include "ops/src/predicates/ops_is_object_pybind.hpp"

#include "ops/src/constraints/ops_sharedmem_host_subscriptable_rank1_container.hpp"
#include "ops/src/constraints/ops_sharedmem_host_subscriptable_rank2_container.hpp"

#include "ops/src/ops_matching_extents_impl.hpp"

// Eigen
#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "ops/src/constraints/ops_container_eigen_with_native_data_access.hpp"
#include "ops/src/constraints/ops_rank1_container_eigen_with_native_data_access.hpp"
#include "ops/src/constraints/ops_rank2_container_eigen_with_native_data_access.hpp"
#include "ops/src/eigen/ops_abs.hpp"
#include "ops/src/eigen/ops_set_zero.hpp"
#include "ops/src/eigen/ops_scale.hpp"
#include "ops/src/eigen/ops_fill.hpp"
#include "ops/src/eigen/ops_resize.hpp"
#include "ops/src/eigen/ops_deep_copy.hpp"
#include "ops/src/eigen/ops_add_to_diagonal.hpp"
#include "ops/src/eigen/ops_min_max.hpp"
#include "ops/src/eigen/ops_level2.hpp"
#include "ops/src/eigen/ops_level3.hpp"
#include "ops/src/eigen/ops_multi_vector_update.hpp"
#include "ops/src/eigen/ops_norms_vector.hpp"
#include "ops/src/eigen/ops_dot.hpp"
#include "ops/src/eigen/ops_vector_update.hpp"
#include "ops/src/eigen/ops_pow.hpp"
#include "ops/src/eigen/ops_elementwise_multiply.hpp"
#endif

// Kokkos
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "ops/src/constraints/ops_container_kokkos_with_native_data_access.hpp"
#include "ops/src/constraints/ops_rank1_container_kokkos_with_native_data_access.hpp"
#include "ops/src/constraints/ops_rank2_container_kokkos_with_native_data_access.hpp"
#include "ops/src/kokkos/ops_abs.hpp"
#include "ops/src/kokkos/ops_set_zero.hpp"
#include "ops/src/kokkos/ops_scale.hpp"
#include "ops/src/kokkos/ops_fill.hpp"
#include "ops/src/kokkos/ops_deep_copy.hpp"
#include "ops/src/kokkos/ops_level2.hpp"
#include "ops/src/kokkos/ops_level3.hpp"
#include "ops/src/kokkos/ops_norms_vector.hpp"
#include "ops/src/kokkos/ops_vector_update_kokkos_functors.hpp"
#include "ops/src/kokkos/ops_vector_update.hpp"
#include "ops/src/kokkos/ops_multi_vector_update.hpp"
#include "ops/src/kokkos/ops_dot.hpp"
#include "ops/src/kokkos/ops_elementwise_multiply.hpp"
#include "ops/src/kokkos/ops_pow.hpp"
#endif

// Epetra
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "ops/src/epetra/ops_abs.hpp"
#include "ops/src/epetra/ops_set_zero.hpp"
#include "ops/src/epetra/ops_fill.hpp"
#include "ops/src/epetra/ops_deep_copy.hpp"
#include "ops/src/epetra/ops_min_max_vector.hpp"
#include "ops/src/epetra/ops_level2.hpp"
#include "ops/src/epetra/ops_level3.hpp"
#include "ops/src/epetra/ops_norms_vector.hpp"
#include "ops/src/epetra/ops_vector_update.hpp"
#include "ops/src/epetra/ops_dot.hpp"
#include "ops/src/epetra/ops_pow.hpp"
#include "ops/src/epetra/ops_elementwise_multiply.hpp"

// teuchos
#include "ops/src/teuchos/ops_set_zero.hpp"
#include "ops/src/teuchos/ops_fill.hpp"
#include "ops/src/teuchos/ops_deep_copy.hpp"
#include "ops/src/teuchos/ops_norms_vector.hpp"
#include "ops/src/teuchos/ops_vector_update.hpp"
#include "ops/src/teuchos/ops_level2.hpp"

// Tpetra
#include "ops/src/tpetra/ops_abs.hpp"
#include "ops/src/tpetra/ops_set_zero.hpp"
#include "ops/src/tpetra/ops_fill.hpp"
#include "ops/src/tpetra/ops_deep_copy.hpp"
#include "ops/src/tpetra/ops_level2.hpp"
#include "ops/src/tpetra/ops_level3.hpp"
#include "ops/src/tpetra/ops_norms_vector.hpp"
#include "ops/src/tpetra/ops_vector_update.hpp"
#include "ops/src/tpetra/ops_multi_vector_update.hpp"
#include "ops/src/tpetra/ops_dot.hpp"
#include "ops/src/tpetra/ops_pow.hpp"
#include "ops/src/tpetra/ops_elementwise_multiply.hpp"

// Tpetra block
#include "ops/src/tpetra_block/ops_abs.hpp"
#include "ops/src/tpetra_block/ops_set_zero.hpp"
#include "ops/src/tpetra_block/ops_fill.hpp"
#include "ops/src/tpetra_block/ops_deep_copy.hpp"
#include "ops/src/tpetra_block/ops_level2.hpp"
#include "ops/src/tpetra_block/ops_level3.hpp"
#include "ops/src/tpetra_block/ops_norms_vector.hpp"
#include "ops/src/tpetra_block/ops_vector_update.hpp"
#include "ops/src/tpetra_block/ops_multi_vector_update.hpp"
#include "ops/src/tpetra_block/ops_dot.hpp"
#include "ops/src/tpetra_block/ops_pow.hpp"
#include "ops/src/tpetra_block/ops_elementwise_multiply.hpp"
#endif

// pybind11
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include "ops/src/constraints/ops_rank1_container_pybind.hpp"
#include "ops/src/pybind11/ops_abs.hpp"
#include "ops/src/pybind11/ops_fill.hpp"
#include "ops/src/pybind11/ops_set_zero.hpp"
#include "ops/src/pybind11/ops_deep_copy.hpp"
#include "ops/src/pybind11/ops_rank1_update.hpp"
#include "ops/src/pybind11/ops_rank2_update.hpp"
#include "ops/src/pybind11/ops_rank3_update.hpp"
#include "ops/src/pybind11/ops_level2.hpp"
#include "ops/src/pybind11/ops_level3.hpp"
#include "ops/src/pybind11/ops_scale.hpp"
#include "ops/src/pybind11/ops_dot.hpp"
#include "ops/src/pybind11/ops_pow.hpp"
#include "ops/src/pybind11/ops_norms_vector.hpp"
#include "ops/src/pybind11/ops_elementwise_multiply.hpp"
#endif

#endif
