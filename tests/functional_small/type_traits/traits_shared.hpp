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
namespace pressio { namespace traits { namespace test {
// -------------------------------------------------

/*
  Verify package predicates
*/

template <typename T>
void test_is_not_eigen_container() {
#ifdef PRESSIO_ENABLE_TPL_EIGEN
  static_assert(!pressio::is_native_container_eigen<T>::value, "");
  static_assert(!pressio::is_vector_eigen<T>::value, "");
  static_assert(!pressio::is_dynamic_vector_eigen<T>::value, "");
  static_assert(!pressio::is_static_vector_eigen<T>::value, "");
  static_assert(!pressio::is_dynamic_row_vector_eigen<T>::value, "");
  static_assert(!pressio::is_static_row_vector_eigen<T>::value, "");
  static_assert(!pressio::is_dynamic_column_vector_eigen<T>::value, "");
  static_assert(!pressio::is_static_column_vector_eigen<T>::value, "");
  static_assert(!pressio::is_dense_matrix_eigen<T>::value, "");
  static_assert(!pressio::is_sparse_matrix_eigen<T>::value, "");
  static_assert(!pressio::is_static_dense_matrix_eigen<T>::value, "");
  static_assert(!pressio::is_dynamic_dense_matrix_eigen<T>::value, "");
  static_assert(!pressio::is_dense_row_major_matrix_eigen<T>::value, "");
  static_assert(!pressio::is_dense_col_major_matrix_eigen<T>::value, "");
#endif
}

template <typename T>
void test_is_not_teuchos_container() {
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  static_assert(!pressio::is_dense_vector_teuchos<T>::value,"");
  static_assert(!pressio::is_dense_matrix_teuchos<T>::value,"");
#endif
}

template <typename T>
void test_is_not_epetra_container() {
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  static_assert(!pressio::is_vector_epetra<T>::value, "");
  static_assert(!pressio::is_multi_vector_epetra<T>::value, "");
#endif
}

template <typename T>
void test_is_not_tpetra_container() {
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  static_assert(!pressio::is_vector_tpetra<T>::value, "");
  static_assert(!pressio::is_multi_vector_tpetra<T>::value, "");
  static_assert(!pressio::is_vector_tpetra_block<T>::value, "");
  static_assert(!pressio::is_multi_vector_tpetra_block<T>::value, "");
#endif
}

template <typename T>
void test_is_not_kokkos_container() {
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
  static_assert(!pressio::is_native_container_kokkos<T>::value,"");
  static_assert(!pressio::is_dynamic_vector_kokkos<T>::value,"");
  static_assert(!pressio::is_static_vector_kokkos<T>::value,"");
  static_assert(!pressio::is_vector_kokkos<T>::value,"");
  static_assert(!pressio::is_static_dense_matrix_kokkos<T>::value,"");
  static_assert(!pressio::is_dynamic_dense_matrix_kokkos<T>::value,"");
  static_assert(!pressio::is_dense_matrix_kokkos<T>::value,"");
#endif
}

// -------------------------------------------------

/*
    Verifies traits common for all containers
*/
template <
  typename T,
  int rank,
  typename Scalar,
  typename traits = pressio::Traits<T>
>
void test_container_traits()
{
  static_assert(traits::rank == rank, "rank is different than expected");
  testing::StaticAssertTypeEq<typename traits::scalar_type, Scalar>();
}

// -------------------------------------------------
}}} // pressio::traits::test