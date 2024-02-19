/*
//@HEADER
// ************************************************************************
//
// type_traits.hpp
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

#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "traits_shared.hpp"

namespace pressio { namespace traits { namespace test {

template <
  typename T,
  typename traits = pressio::Traits<T>
>
void test_kokkos_container_traits()
{
  // traits and shared predicates
  // constexpr bool is_dynamic = T::traits::rank_dynamic != 0;
  test_container_traits<
    T,
    T::traits::rank,
    typename T::traits::value_type
  >();

  // negative checks (cross-package)
  test_is_not_eigen_container<T>();
  test_is_not_teuchos_container<T>();
  test_is_not_tpetra_container<T>();
  test_is_not_tpetra_block_container<T>();
}

//*******************************
// Kokkos vector
//*******************************

/*
  Verify values of Kokkos matrix traits and relatad predicates
*/
template <
  typename T,
  typename traits = pressio::Traits<T>
>
void test_kokkos_vector_type_traits()
{
  // traits and shared predicates
  test_kokkos_container_traits<T>();

  // vector predicates
  constexpr bool is_dynamic = T::traits::rank_dynamic != 0;
  static_assert(pressio::is_dynamic_vector_kokkos<T>::value == is_dynamic,"");
  static_assert(pressio::is_static_vector_kokkos<T>::value == !is_dynamic,"");
  static_assert(pressio::is_vector_kokkos<T>::value,"");

  // negative checks (within Kokkos)
  static_assert(pressio::is_dense_matrix_kokkos<T>::value == false,"");
  static_assert(pressio::is_static_dense_matrix_kokkos<T>::value == false,"");
  static_assert(pressio::is_dynamic_dense_matrix_kokkos<T>::value == false,"");
}

TEST(type_traits, kokkos_vector) {
  test_kokkos_vector_type_traits<
    Kokkos::View<double*>
  >();
  test_kokkos_vector_type_traits<
    Kokkos::View<float[32]>
  >();
}

//*******************************
// Kokkos dense matrix
//*******************************

/*
  Verify values of Kokkos matrix traits and relatad predicates
*/
template <
  typename T,
  typename traits = pressio::Traits<T>
>
void test_kokkos_matrix_type_traits()
{
  // traits and shared predicates
  test_kokkos_container_traits<T>();

  // static_assert(std::is_same<
  //       typename T::traits::array_layout,
  //       Kokkos::LayoutLeft
  //     >::value, "");

  // dense matrix predicates
  constexpr bool is_dynamic = T::traits::rank_dynamic != 0;
  static_assert(pressio::is_static_dense_matrix_kokkos<T>::value == !is_dynamic,"");
  static_assert(pressio::is_dynamic_dense_matrix_kokkos<T>::value == is_dynamic,"");
  static_assert(pressio::is_dense_matrix_kokkos<T>::value,"");

  // negative checks (within Kokkos)
  static_assert(pressio::is_dynamic_vector_kokkos<T>::value == false,"");
  static_assert(pressio::is_static_vector_kokkos<T>::value == false,"");
  static_assert(pressio::is_vector_kokkos<T>::value == false,"");
}

TEST(type_traits, kokkos_matrix) {
  test_kokkos_matrix_type_traits<
    Kokkos::View<double**>
  >();
  test_kokkos_matrix_type_traits<
    Kokkos::View<float[32][32]>
  >();
}

}}} // pressio::traits::test
