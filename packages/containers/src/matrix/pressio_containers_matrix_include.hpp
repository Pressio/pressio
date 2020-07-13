/*
//@HEADER
// ************************************************************************
//
// pressio_containers_matrix_include.hpp
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

#ifndef CONTAINERS_MATRIX_PRESSIO_CONTAINERS_MATRIX_INCLUDE_HPP_
#define CONTAINERS_MATRIX_PRESSIO_CONTAINERS_MATRIX_INCLUDE_HPP_

// predicates native
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./predicates/containers_native_kokkos_matrix_meta.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./predicates/containers_native_trilinos_matrix_meta.hpp"
#endif
#include "./predicates/containers_native_eigen_matrix_meta.hpp"
#include "./predicates/containers_native_arbitrary_matrix_meta.hpp"

// predicates wrappers
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./predicates/containers_is_dense_matrix_wrapper_epetra.hpp"
#include "./predicates/containers_is_dense_matrix_wrapper_teuchos.hpp"
#include "./predicates/containers_is_sparse_matrix_wrapper_epetra.hpp"
#include "./predicates/containers_is_sparse_matrix_wrapper_tpetra.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./predicates/containers_is_dense_matrix_wrapper_kokkos.hpp"
#include "./predicates/containers_is_sparse_matrix_wrapper_kokkos.hpp"
#include "./predicates/containers_is_matrix_wrapper_kokkos.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include "./predicates/containers_is_matrix_wrapper_pybind11.hpp"
#endif
#include "./predicates/containers_is_sparse_matrix_wrapper_eigen.hpp"
#include "./predicates/containers_is_dense_matrix_wrapper_eigen.hpp"
#include "./predicates/containers_is_matrix_wrapper_eigen.hpp"
#include "./predicates/containers_is_matrix_wrapper_arbitrary.hpp"
#include "./predicates/containers_is_matrix_wrapper.hpp"

// traits
#include "./containers_matrix_traits.hpp"

// concrete
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include "./concrete/containers_matrix_sharedmem_pybind11.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./concrete/containers_matrix_dense_distributed_epetra.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_teuchos_serial.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./concrete/containers_matrix_dense_sharedmem_kokkos.hpp"
#endif
#include "./concrete/containers_matrix_arbitrary.hpp"
#include "./concrete/containers_matrix_sparse_sharedmem_eigen.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_eigen_dynamic.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_eigen_static.hpp"

#endif  // CONTAINERS_MATRIX_PRESSIO_CONTAINERS_MATRIX_INCLUDE_HPP_
