/*
//@HEADER
// ************************************************************************
//
// solvers_get_matrix_size_helper.hpp
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

#ifndef SOLVERS_GET_MATRIX_SIZE_HELPER_HPP
#define SOLVERS_GET_MATRIX_SIZE_HELPER_HPP

namespace pressio{ namespace solvers{ namespace impl{

template <typename T, typename enable = void>
struct MatrixGetSizeHelper;

template <typename T>
struct MatrixGetSizeHelper<
  T,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper<T>::value and
    containers::details::traits<T>::is_shared_mem == false
    >
  >{
  static auto globalRows(const T & A) -> decltype(A.extent(0)){
    return A.extent(0);
  }
  static auto globalCols(const T & A) -> decltype(A.numVectors()){
    return A.numVectors();
  }
};


template <typename T>
struct MatrixGetSizeHelper<
  T,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_matrix_wrapper<T>::value and
    containers::details::traits<T>::is_shared_mem == false
    >
  >{
  static auto globalRows(const T & A) -> decltype(A.extent(0)){
    return A.extent(0);
  }
  static auto globalCols(const T & A) -> decltype(A.extent(1)){
    return A.extent(1);
  }
};

template <typename T>
struct MatrixGetSizeHelper<
  T,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_matrix_wrapper<T>::value and
    containers::details::traits<T>::is_shared_mem == true
    >
  >{
  static auto globalRows(const T & A) -> decltype(A.extent(0)){
    return A.extent(0);
  }
  static auto globalCols(const T & A) -> decltype(A.extent(1)){
    return A.extent(1);
  }
};

template <typename T>
struct MatrixGetSizeHelper<
  T,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper<T>::value and
    containers::details::traits<T>::is_shared_mem == true
    >
  >{
  static auto globalRows(const T & A) -> decltype(A.extent()){
    return A.extent();
  }
  static auto globalCols(const T & A) -> decltype(A.numVectors()){
    return A.numVectors();
  }
};

}}} //end namespace pressio::solvers::impl
#endif
