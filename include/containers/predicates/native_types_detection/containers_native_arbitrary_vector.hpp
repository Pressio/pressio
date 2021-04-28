/*
//@HEADER
// ************************************************************************
//
// containers_native_arbitrary_vector.hpp
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

#ifndef CONTAINERS_PREDICATES_NATIVE_TYPES_DETECTION_CONTAINERS_NATIVE_ARBITRARY_VECTOR_HPP_
#define CONTAINERS_PREDICATES_NATIVE_TYPES_DETECTION_CONTAINERS_NATIVE_ARBITRARY_VECTOR_HPP_

namespace pressio{ namespace containers{ namespace predicates {

/*
  T is admissible to be wrapped as an arbitrary vector iff it is:
  - not one the supported vector types
  - not one of the supported dense matrices
  - not one of the supported multivectors
  - not already a wrapper

  NOTE that these checks are necessary to guard against cases like this:
    using T = Eigen::MatrixXd;
    using v_t = pressio::containers::Vector<T>;

    If we were to just check that T is not an eigen vector,
    then this case would still compile. But this is wrong.
    This is why we need to make sure T is not a matrix or multivector.
*/

template <typename T, typename enable = void>
struct is_admissible_as_vector_arbitrary : std::false_type
{
  static_assert
  (!containers::predicates::is_wrapper<T>::value,
   "You cannot wrap a pressio container as a pressio::containers::Vector<>.");

#ifdef PRESSIO_ENABLE_TPL_EIGEN
  static_assert
  (!containers::predicates::is_dense_matrix_eigen<T>::value and
   !containers::predicates::is_admissible_as_dynamic_multi_vector_eigen<T>::value,
   "You cannot wrap an Eigen dense matrix as a pressio::containers::Vector<>.");

  static_assert
  (!containers::predicates::is_sparse_matrix_eigen<T>::value,
   "You cannot wrap an Eigen sparse matrix as a pressio::containers::Vector<>.");
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
  static_assert
  (!containers::predicates::is_dense_matrix_kokkos<T>::value and
   !containers::predicates::is_admissible_as_multi_vector_kokkos<T>::value,
   "You cannot wrap a Kokkos 2d view as a pressio::containers::Vector<>.");
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  static_assert
  (!containers::predicates::is_multi_vector_epetra<T>::value and
   !containers::predicates::is_admissible_as_dense_matrix_epetra<T>::value,
   "You cannot wrap an Epetra multivector as a pressio::containers::Vector<>.");

  static_assert
  (!containers::predicates::is_dense_matrix_teuchos<T>::value
   or !containers::predicates::is_dense_matrix_teuchos_rcp<T>::value,
   "You cannot wrap a Teuchos dense matrix as a pressio::containers::Vector<>.");

  static_assert
  (!containers::predicates::is_multi_vector_tpetra_block<T>::value,
   "You cannot wrap a Tpetra block multivector as a pressio::containers::Vector<>.");

  static_assert
  (!containers::predicates::is_multi_vector_tpetra<T>::value,
   "You cannot wrap a Tpetra multivector as a pressio::containers::Vector<>.");
#endif
};

template <typename T>
struct is_admissible_as_vector_arbitrary<
  T,
  ::pressio::mpl::enable_if_t<
    !std::is_void<T>::value
    and !containers::predicates::is_wrapper<T>::value
#ifdef PRESSIO_ENABLE_TPL_EIGEN
    and !containers::predicates::is_vector_eigen<T>::value
    and !containers::predicates::is_dense_matrix_eigen<T>::value
    and !containers::predicates::is_sparse_matrix_eigen<T>::value
    and !containers::predicates::is_admissible_as_multi_vector_eigen<T>::value
#endif
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    and !containers::predicates::is_array_pybind<T>::value
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
    and !containers::predicates::is_vector_kokkos<T>::value
    and !containers::predicates::is_dense_matrix_kokkos<T>::value
    and !containers::predicates::is_admissible_as_multi_vector_kokkos<T>::value
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
    and !containers::predicates::is_vector_epetra<T>::value
    and !containers::predicates::is_multi_vector_epetra<T>::value
    and !containers::predicates::is_admissible_as_dense_matrix_epetra<T>::value
    //
    and !containers::predicates::is_dense_vector_teuchos<T>::value
    and !containers::predicates::is_dense_matrix_teuchos<T>::value
    and !containers::predicates::is_dense_matrix_teuchos_rcp<T>::value
    //
    and !containers::predicates::is_vector_tpetra_block<T>::value
    and !containers::predicates::is_multi_vector_tpetra_block<T>::value
    //
    and !containers::predicates::is_vector_tpetra<T>::value
    and !containers::predicates::is_multi_vector_tpetra<T>::value
#endif
    >
  > : std::true_type{};


}}}//end namespace pressio::containers::predicates
#endif  // CONTAINERS_PREDICATES_NATIVE_TYPES_DETECTION_CONTAINERS_NATIVE_ARBITRARY_VECTOR_HPP_
