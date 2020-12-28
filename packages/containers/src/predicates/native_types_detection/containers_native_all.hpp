/*
//@HEADER
// ************************************************************************
//
// containers_native_all.hpp
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

#ifndef CONTAINERS_PREDICATES_NATIVE_TYPES_DETECTION_CONTAINERS_NATIVE_ALL_HPP_
#define CONTAINERS_PREDICATES_NATIVE_TYPES_DETECTION_CONTAINERS_NATIVE_ALL_HPP_

//-----------------------------------------
// predicates for native types detection
//-----------------------------------------
// for pybind arrays there is no way to detect if a native type
// it is a vector or something else at compile time. We only have an "array"
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include "containers_native_pybind_array.hpp"
#endif

//*** vector ****
#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "containers_native_eigen_vector.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "containers_native_kokkos_vector.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "containers_native_epetra_vector.hpp"
#include "containers_native_teuchos_vector.hpp"
#include "containers_native_tpetra_block_vector.hpp"
#include "containers_native_tpetra_vector.hpp"
#endif

//*** matrix ****
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "containers_native_kokkos_dense_matrix.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "containers_native_epetra_dense_matrix.hpp"
#include "containers_native_teuchos_dense_matrix.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "containers_native_eigen_dense_matrix.hpp"
#include "containers_native_eigen_sparse_matrix.hpp"
#endif

//*** multi vector ****
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "containers_native_kokkos_multi_vector.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "containers_native_epetra_multi_vector.hpp"
#include "containers_native_tpetra_block_multi_vector.hpp"
#include "containers_native_tpetra_multi_vector.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "containers_native_eigen_multi_vector.hpp"
#endif

// arbitrary must be at end because they depend on the above
#include "containers_native_arbitrary_vector.hpp"
#include "containers_native_arbitrary_dense_matrix.hpp"
#include "containers_native_arbitrary_multi_vector.hpp"


#endif  // CONTAINERS_PREDICATES_NATIVE_TYPES_DETECTION_CONTAINERS_NATIVE_ALL_HPP_
