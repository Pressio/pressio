/*
//@HEADER
// ************************************************************************
//
// traits_tpl.hpp
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

#ifndef TYPE_TRAITS_TRAITS_TPL_HPP_
#define TYPE_TRAITS_TRAITS_TPL_HPP_

namespace pressio {

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename T>
struct is_native_container_eigen {
  static constexpr auto value = ::pressio::is_vector_eigen<T>::value
    || ::pressio::is_dense_matrix_eigen<T>::value
    || ::pressio::is_sparse_matrix_eigen<T>::value;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename T>
struct is_native_container_kokkos {
  static constexpr auto value = ::pressio::is_vector_kokkos<T>::value
    || ::pressio::is_dense_matrix_kokkos<T>::value;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename T>
struct Traits< T,
  std::enable_if_t< is_vector_eigen<T>::value > >
{
  static constexpr int rank = 1;
  using scalar_type = typename T::Scalar;
};
#endif //PRESSIO_ENABLE_TPL_EIGEN

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename T>
struct Traits< T,
  std::enable_if_t< is_vector_kokkos<T>::value > >
{
  static constexpr int rank = 1;
  using scalar_type = typename T::traits::value_type;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename T>
struct Traits< T,
   std::enable_if_t<
    is_vector_tpetra<T>::value || is_vector_tpetra_block<T>::value > >
{
  static constexpr int rank = 1;
  using scalar_type = typename T::impl_scalar_type;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename T>
struct Traits< T,
  std::enable_if_t<
    is_dense_vector_teuchos<T>::value > >
{
  static constexpr int rank = 1;
  using scalar_type = typename T::scalarType;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename T>
struct Traits<
  T, std::enable_if_t<
    is_dense_matrix_eigen<T>::value >
  >
{
  static constexpr int rank = 2;
  using scalar_type = typename T::Scalar;
};

#endif //PRESSIO_ENABLE_TPL_EIGEN

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename T>
struct Traits<
  T, std::enable_if_t<
      is_dense_matrix_kokkos<T>::value >
  >
{
  static constexpr int rank = 2;
  using scalar_type = typename T::traits::value_type;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename T>
struct Traits<
  T, std::enable_if_t<
    is_dense_matrix_teuchos<T>::value >
  >
{
  static constexpr int rank = 2;
  using scalar_type = typename T::scalarType;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename T>
struct Traits<
  T, std::enable_if_t<
    is_multi_vector_tpetra<T>::value ||
    is_multi_vector_tpetra_block<T>::value >
  >
{
  static constexpr int rank = 2;
  using scalar_type = typename T::impl_scalar_type;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename T>
struct Traits<
  T, std::enable_if_t<
    is_sparse_matrix_eigen<T>::value >
  >
{
  static constexpr int rank = 2;
  using scalar_type = typename T::Scalar;

private:
  using ordinal_type_ = typename T::StorageIndex;

  static_assert(
    std::is_integral<ordinal_type_>::value &&
    std::is_signed<ordinal_type_>::value,
    "ordinal type for indexing eigen sparse matrix has to be signed"
  );
};
#endif //PRESSIO_ENABLE_TPL_EIGEN


} // namespace pressio
#endif  // TYPE_TRAITS_TRAITS_TPL_HPP_
