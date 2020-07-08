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

#include "ops/src/ops_ConfigDefs.hpp"
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
#include "ops/src/predicates/ops_has_method_do_update_one_term.hpp"
#include "ops/src/predicates/ops_has_method_do_update_two_terms.hpp"
#include "ops/src/predicates/ops_has_method_do_update_three_terms.hpp"
#include "ops/src/predicates/ops_has_method_do_update_four_terms.hpp"

// Eigen
#include "ops/src/eigen/ops_set_zero.hpp"
#include "ops/src/eigen/ops_scale.hpp"
#include "ops/src/eigen/ops_fill.hpp"
#include "ops/src/eigen/ops_resize.hpp"
#include "ops/src/eigen/ops_deep_copy.hpp"
#include "ops/src/eigen/ops_add_to_diagonal.hpp"
#include "ops/src/eigen/ops_min_max.hpp"
#include "ops/src/eigen/ops_mat_prod_vec.hpp"
#include "ops/src/eigen/ops_mat_prod_mat.hpp"
#include "ops/src/eigen/ops_multi_vector_do_update.hpp"
#include "ops/src/eigen/ops_norms_vector.hpp"
#include "ops/src/eigen/ops_vec_dot_vec.hpp"
#include "ops/src/eigen/ops_vector_do_update.hpp"
#include "ops/src/eigen/ops_elementwise_multiply.hpp"

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
// Kokkos
#include "ops/src/kokkos/ops_set_zero.hpp"
#include "ops/src/kokkos/ops_scale.hpp"
#include "ops/src/kokkos/ops_fill.hpp"
#include "ops/src/kokkos/ops_deep_copy.hpp"
#include "ops/src/kokkos/ops_mvec_prod_mvec.hpp"
#include "ops/src/kokkos/ops_mvec_prod_vec.hpp"
#include "ops/src/kokkos/ops_norms_vector.hpp"
#include "ops/src/kokkos/ops_vector_do_update_kokkos_functors.hpp"
#include "ops/src/kokkos/ops_vector_do_update.hpp"
#include "ops/src/kokkos/ops_multi_vector_do_update.hpp"
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
// Epetra
#include "ops/src/epetra/ops_set_zero.hpp"
#include "ops/src/epetra/ops_fill.hpp"
#include "ops/src/epetra/ops_deep_copy.hpp"
#include "ops/src/epetra/ops_min_max_vector.hpp"
#include "ops/src/epetra/ops_mvec_prod_mvec.hpp"
#include "ops/src/epetra/ops_mvec_prod_vec.hpp"
#include "ops/src/epetra/ops_norms_vector.hpp"
#include "ops/src/epetra/ops_vector_do_update.hpp"
// teuchos
#include "ops/src/teuchos/ops_set_zero.hpp"
#include "ops/src/teuchos/ops_fill.hpp"
#include "ops/src/teuchos/ops_deep_copy.hpp"
#include "ops/src/teuchos/ops_norms_vector.hpp"
#include "ops/src/teuchos/ops_vector_do_update.hpp"
#include "ops/src/teuchos/ops_mat_prod_vec.hpp"
// Tpetra
#include "ops/src/tpetra/ops_set_zero.hpp"
#include "ops/src/tpetra/ops_fill.hpp"
#include "ops/src/tpetra/ops_deep_copy.hpp"
#include "ops/src/tpetra/ops_mvec_prod_mvec.hpp"
#include "ops/src/tpetra/ops_mvec_prod_vec.hpp"
#include "ops/src/tpetra/ops_norms_vector.hpp"
#include "ops/src/tpetra/ops_vector_do_update.hpp"
#include "ops/src/tpetra/ops_multi_vector_do_update.hpp"
// Tpetra block
#include "ops/src/tpetra_block/ops_set_zero.hpp"
#include "ops/src/tpetra_block/ops_fill.hpp"
#include "ops/src/tpetra_block/ops_deep_copy.hpp"
#include "ops/src/tpetra_block/ops_mvec_prod_mvec.hpp"
#include "ops/src/tpetra_block/ops_mvec_prod_vec.hpp"
#include "ops/src/tpetra_block/ops_norms_vector.hpp"
#include "ops/src/tpetra_block/ops_vector_do_update.hpp"
#include "ops/src/tpetra_block/ops_multi_vector_do_update.hpp"
#endif

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
// pybind11
#include "ops/src/pybind11/ops_set_zero.hpp"
#include "ops/src/pybind11/ops_deep_copy.hpp"
#include "ops/src/pybind11/ops_norms_vector.hpp"
#include "ops/src/pybind11/ops_vector_do_update.hpp"
#include "ops/src/pybind11/ops_mat_prod_vec.hpp"
#endif

#endif
