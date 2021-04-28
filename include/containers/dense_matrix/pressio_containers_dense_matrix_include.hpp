/*
//@HEADER
// ************************************************************************
//
// pressio_containers_dense_matrix_include.hpp
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

#ifndef CONTAINERS_DENSE_MATRIX_PRESSIO_CONTAINERS_DENSE_MATRIX_INCLUDE_HPP_
#define CONTAINERS_DENSE_MATRIX_PRESSIO_CONTAINERS_DENSE_MATRIX_INCLUDE_HPP_

/* WARNING: the inclusion order below matters: */

//--------
// traits
//--------
#include "./containers_dense_matrix_traits.hpp"

//---------------
// concrete types
//---------------
#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "./concrete/containers_matrix_dense_sharedmem_eigen_dynamic.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_eigen_static.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./concrete/containers_matrix_dense_distributed_epetra.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_teuchos_serial.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./concrete/containers_matrix_dense_sharedmem_kokkos.hpp"
#endif
// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// #include "./concrete/containers_matrix_dense_sharedmem_pybind11.hpp"
// #endif
#include "./concrete/containers_matrix_dense_arbitrary.hpp"

//----------------------
// wrapper predicates
//----------------------
#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "./wrapper_predicates/containers_is_dense_matrix_wrapper_eigen.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./wrapper_predicates/containers_is_dense_matrix_wrapper_epetra.hpp"
#include "./wrapper_predicates/containers_is_dense_matrix_wrapper_teuchos.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./wrapper_predicates/containers_is_dense_matrix_wrapper_kokkos.hpp"
#endif
// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// #include "./wrapper_predicates/containers_is_dense_matrix_wrapper_pybind11.hpp"
// #endif
#include "./wrapper_predicates/containers_is_dense_matrix_wrapper_arbitrary.hpp"
#include "./wrapper_predicates/containers_is_dense_matrix_wrapper.hpp"

#endif  // CONTAINERS_DENSE_MATRIX_PRESSIO_CONTAINERS_DENSE_MATRIX_INCLUDE_HPP_
