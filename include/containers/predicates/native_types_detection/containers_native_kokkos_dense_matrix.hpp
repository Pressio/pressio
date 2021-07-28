/*
//@HEADER
// ************************************************************************
//
// containers_native_kokkos_dense_matrix.hpp
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

#ifndef CONTAINERS_PREDICATES_NATIVE_TYPES_DETECTION_CONTAINERS_NATIVE_KOKKOS_DENSE_MATRIX_HPP_
#define CONTAINERS_PREDICATES_NATIVE_TYPES_DETECTION_CONTAINERS_NATIVE_KOKKOS_DENSE_MATRIX_HPP_

#include <KokkosSparse_CrsMatrix.hpp>

namespace pressio{ namespace containers{ namespace predicates {

// T is a static kokkos view if T:
// - is a view with rank==2
// - the number of runtime determined dimensions == 0
template <typename T, typename enable = void>
struct is_static_dense_matrix_kokkos : std::false_type {};

template <typename T>
struct is_static_dense_matrix_kokkos<
  T,
  ::pressio::mpl::enable_if_t<
    Kokkos::is_view<T>::value &&
    T::traits::rank==2 &&
    T::traits::rank_dynamic==0
    >
  > : std::true_type{};

// -------------------------------------------------
template <typename T, typename enable = void>
struct is_dynamic_dense_matrix_kokkos : std::false_type {};

template <typename T>
struct is_dynamic_dense_matrix_kokkos<
  T,
  ::pressio::mpl::enable_if_t<
    Kokkos::is_view<T>::value &&
    T::traits::rank==2 &&
    T::traits::rank_dynamic!=0
    >
  > : std::true_type{};

// -------------------------------------------------
template <typename T, typename enable = void>
struct is_dense_matrix_kokkos : std::false_type {};

template <typename T>
struct is_dense_matrix_kokkos<
  T,
  ::pressio::mpl::enable_if_t<
    is_static_dense_matrix_kokkos<T>::value or
    is_dynamic_dense_matrix_kokkos<T>::value
    >
  > : std::true_type{};

}}}//end namespace pressio::containers::predicates
#endif  // CONTAINERS_PREDICATES_NATIVE_TYPES_DETECTION_CONTAINERS_NATIVE_KOKKOS_DENSE_MATRIX_HPP_
