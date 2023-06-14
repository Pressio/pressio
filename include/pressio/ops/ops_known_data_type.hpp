/*
//@HEADER
// ************************************************************************
//
// ops_get_native.hpp
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

#ifndef OPS_OPS_KNOWN_DATA_TYPE_HPP_
#define OPS_OPS_KNOWN_DATA_TYPE_HPP_

namespace pressio{ namespace ops{


template<class T>
struct is_known_data_type
{
  static constexpr bool v1 =
#ifdef PRESSIO_ENABLE_TPL_EIGEN
    is_vector_eigen<T>::value
    || is_dense_matrix_eigen<T>::value
    || is_sparse_matrix_eigen<T>::value;
#else
  false;
#endif

  static constexpr bool v2 =
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
    is_vector_kokkos<T>::value
    || is_dense_matrix_kokkos<T>::value;
#else
  false;
#endif

  static constexpr bool v3 =
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
    is_vector_tpetra<T>::value
    || is_multi_vector_tpetra<T>::value
    || is_vector_tpetra_block<T>::value
    || is_multi_vector_tpetra_block<T>::value
    || is_vector_epetra<T>::value
    || is_multi_vector_epetra<T>::value;
#else
    false;
#endif

  static constexpr bool value = v1 || v2 || v3;
};

template<class ... Args>
struct are_known_data_types
{
  static constexpr bool value = ::pressio::mpl::all_of<
    is_known_data_type, Args...>::value;
};

}}
#endif  // OPS_OPS_KNOWN_DATA_TYPE_HPP_
