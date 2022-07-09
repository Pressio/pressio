/*
//@HEADER
// ************************************************************************
//
// containers_multi_vector_traits.hpp
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

namespace pressio { namespace impl {

#ifdef PRESSIO_ENABLE_TPL_TRILINOS

//*******************************
// Common traits
//*******************************

template<
    int Rank,
    typename Scalar,
    typename LocalOrdinal,
    typename GlobalOrdinal,
    typename DataMap,
    typename Comm
>
struct TrilinosTraits
  : public ::pressio::impl::ContainersSharedTraits<
      PackageIdentifier::Trilinos,
      false,
      Rank
    >,
    public ::pressio::impl::ScalarTrait<Scalar>,
    public ::pressio::impl::OrdinalTrait<LocalOrdinal>
{
};

//*******************************
// Epetra traits
//*******************************
template<int Rank>
using EpetraTraits = TrilinosTraits<
  Rank,
  double, int, int,
  Epetra_BlockMap,
  Epetra_Comm
>;

//*******************************
// Tpetra traits
//*******************************
template<typename T>
struct TrilinosCommType
{
  using type = decltype(
    std::declval<
      typename T::map_type
    >().getComm()
  );
};

template<typename T, int Rank>
struct TpetraTraits
  : public ::pressio::impl::TrilinosTraits<
      Rank,
      typename T::impl_scalar_type,
      typename T::local_ordinal_type,
      typename T::global_ordinal_type,
      typename T::map_type,
      typename TrilinosCommType<T>::type
    >
{};

//*******************************
// Teuchos traits
//*******************************
template<typename T, int Rank>
struct TeuchosTraits
  : public ::pressio::impl::ContainersSharedTraits<
      PackageIdentifier::Trilinos,
      true,
      Rank
    >,
    public ::pressio::impl::ScalarTrait<typename T::scalarType>,
    public ::pressio::impl::OrdinalTrait<typename T::ordinalType>
{
};

#endif // PRESSIO_ENABLE_TPL_TRILINOS

//*******************************
// Kokkos traits
//*******************************
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <
  typename T,
  int Rank,
  bool is_static = T::traits::rank_dynamic == 0 // static has no runtime dimensions
>
struct KokkosTraits
  : public ::pressio::impl::ContainersSharedTraits<
      PackageIdentifier::Kokkos,
      true,
      Rank
    >,
    public ::pressio::impl::OrdinalTrait<
      typename T::traits::size_type
    >,
    public ::pressio::impl::ScalarTrait<
      typename T::traits::value_type
    >
{};
#endif // PRESSIO_ENABLE_TPL_KOKKOS

//*******************************
// Eigen vector helpers
//*******************************
#ifdef PRESSIO_ENABLE_TPL_EIGEN

template <typename T, int Rank>
struct EigenTraits
  : public ::pressio::impl::ContainersSharedTraits<
      PackageIdentifier::Eigen,
      true,
      Rank
    >,
    public ::pressio::impl::ScalarTrait<typename T::Scalar>,
    public ::pressio::impl::OrdinalTrait<typename T::StorageIndex>
{};

#endif // PRESSIO_ENABLE_TPL_EIGEN

// device_type is just an (execution space, memory space) pair.
// defined as: Kokkos::Device<execution_space, memory_space>
// so from the device we can get the device execution and memory space
template <typename T, typename Enabled=void> struct device_type;

template <typename T, typename Enabled=void> struct execution_space;

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename T>
struct device_type<
  T,
  ::pressio::mpl::enable_if_t<
    is_vector_kokkos<T>::value
    || ::pressio::is_dense_matrix_kokkos<T>::value
    >
  >
{
  using type = typename T::traits::device_type;
};

template <typename T>
struct execution_space<
  T,
  ::pressio::mpl::enable_if_t<
    is_vector_kokkos<T>::value
    || ::pressio::is_dense_matrix_kokkos<T>::value
    >
  >
{
  using type = typename T::traits::execution_space;
};
#endif


#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template <typename T>
struct device_type<
  T,
  ::pressio::mpl::enable_if_t<
    is_multi_vector_tpetra<T>::value
    || is_multi_vector_tpetra_block<T>::value
    || is_vector_tpetra<T>::value
    || is_vector_tpetra_block<T>::value
    >
  >
{
  using type = typename T::device_type;
};
#endif

}} // namespace pressio::impl
#endif  // TYPE_TRAITS_TRAITS_TPL_HPP_
