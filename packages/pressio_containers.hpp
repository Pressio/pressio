/*
//@HEADER
// ************************************************************************
//
// pressio_containers.hpp
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

#ifndef PRESSIO_CONTAINERS_HPP_
#define PRESSIO_CONTAINERS_HPP_

/*************** IMPORTANT **************
The order below matters because this is how we make sure
headers are included (and classes found) in the proper order.
****************************************/

// pressio dependencies
#include "pressio_mpl.hpp"
#include "pressio_utils.hpp"

// generic
#include "containers/src/containers_ConfigDefs.hpp"
#include "containers/src/containers_fwd.hpp"
#include "containers/src/containers_wrapped_types_enum.hpp"
#include "containers/src/containers_shared_traits.hpp"

// common predicates
#include "containers/src/predicates/typedefs/containers_has_communicator_typedef.hpp"
#include "containers/src/predicates/typedefs/containers_has_data_map_typedef.hpp"
#include "containers/src/predicates/typedefs/containers_has_global_ordinal_typedef.hpp"
#include "containers/src/predicates/typedefs/containers_has_local_ordinal_typedef.hpp"
#include "containers/src/predicates/typedefs/containers_has_ordinal_typedef.hpp"
#include "containers/src/predicates/typedefs/containers_has_scalar_typedef.hpp"
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "containers/src/predicates/containers_is_teuchos_rcp.hpp"
#endif

//-----------------------------------------
// predicates for native types detection
//-----------------------------------------
// for pybind we don't have a distinction for vec/mat/mv, because pybind arrays
// can have dimensions. so there is not yet a way to detect if it is a vector at compile time.
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include "containers/src/predicates/native_types_detection/containers_native_pybind_array.hpp"
#endif

//*** vector ****
#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "containers/src/predicates/native_types_detection/containers_native_eigen_vector.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "containers/src/predicates/native_types_detection/containers_native_kokkos_vector.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "containers/src/predicates/native_types_detection/containers_native_epetra_vector.hpp"
#include "containers/src/predicates/native_types_detection/containers_native_teuchos_vector.hpp"
#include "containers/src/predicates/native_types_detection/containers_native_tpetra_block_vector.hpp"
#include "containers/src/predicates/native_types_detection/containers_native_tpetra_vector.hpp"
#endif

//*** matrix ****
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "containers/src/predicates/native_types_detection/containers_native_kokkos_dense_matrix.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "containers/src/predicates/native_types_detection/containers_native_epetra_dense_matrix.hpp"
#include "containers/src/predicates/native_types_detection/containers_native_teuchos_dense_matrix.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "containers/src/predicates/native_types_detection/containers_native_eigen_dense_matrix.hpp"
#include "containers/src/predicates/native_types_detection/containers_native_eigen_sparse_matrix.hpp"
#endif

//*** multi vector ****
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "containers/src/predicates/native_types_detection/containers_native_kokkos_multi_vector.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "containers/src/predicates/native_types_detection/containers_native_epetra_multi_vector.hpp"
#include "containers/src/predicates/native_types_detection/containers_native_tpetra_block_multi_vector.hpp"
#include "containers/src/predicates/native_types_detection/containers_native_tpetra_multi_vector.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "containers/src/predicates/native_types_detection/containers_native_eigen_multi_vector.hpp"
#endif

// arbitrary must be at end because they depend on the above
#include "containers/src/predicates/native_types_detection/containers_native_arbitrary_vector.hpp"
#include "containers/src/predicates/native_types_detection/containers_native_arbitrary_dense_matrix.hpp"
#include "containers/src/predicates/native_types_detection/containers_native_arbitrary_multi_vector.hpp"

//-------------------------------------------
// include actual Vector/Matrix/MV classes
//-------------------------------------------
#include "containers/src/vector/pressio_containers_vector_include.hpp"
#include "containers/src/dense_matrix/pressio_containers_dense_matrix_include.hpp"
#include "containers/src/sparse_matrix/pressio_containers_sparse_matrix_include.hpp"
#include "containers/src/multi_vector/pressio_containers_multi_vector_include.hpp"

//-------------------------------------------
// expressions
//-------------------------------------------
#include "containers/src/expressions/pressio_containers_expressions_include.hpp"

//-------------------------------------------
// others
//-------------------------------------------
#include "containers/src/predicates/containers_is_wrapper.hpp"
#include "containers/src/predicates/containers_are_wrappers.hpp"
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "containers/src/predicates/containers_have_matching_exe_space.hpp"
#endif
#include "containers/src/predicates/containers_are_scalar_compatible.hpp"
#include "containers/src/collection/containers_static_collection.hpp"

#endif
