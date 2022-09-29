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
    public ::pressio::impl::OrdinalTrait<LocalOrdinal>,
    public ::pressio::impl::DynamicAllocTrait
{
  using local_ordinal_type = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using data_map_type = DataMap;
  using size_type  = GlobalOrdinal;
  using communicator_type = Comm;
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

template<typename T, typename enabled = void> struct TpetraExtraTraits {};

template<typename T>
struct TpetraExtraTraits<
  T,
  mpl::enable_if_t< // NOT Tpetra::BlockVector
    is_vector_tpetra<T>::value or
    is_multi_vector_tpetra<T>::value
  >
>
{
  using dual_view_type = typename T::dual_view_type;
  using dot_type = typename T::dot_type;
  using mag_type = typename T::mag_type;
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
    >,
    public ::pressio::impl::TpetraExtraTraits<T>
{
  // using const_data_return_t = T const *;
  // using data_return_t = T *;

  /* node is a Tpetra concept, defined as:
   * node_type = ::Kokkos::Compat::KokkosDeviceWrapperNode<execution_space>;
   * where memory space is taken from the execution_space
   */
  using node_type = typename T::node_type;

  // device_type is just an (execution space, memory space) pair.
  // defined as: Kokkos::Device<execution_space, memory_space>
  // so from the device we can get the device execution and memory space
  using device_type = typename T::device_type;
  using device_mem_space_type = typename device_type::memory_space;
  using device_exec_space_type = typename device_type::execution_space;

  // store types for host
  using host_mem_space_type = typename Kokkos::HostSpace::memory_space;
  using host_exec_space_type = typename Kokkos::HostSpace::execution_space;
};

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
    public ::pressio::impl::OrdinalTrait<typename T::ordinalType>,
    public ::pressio::impl::DynamicAllocTrait
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
      typename T::traits::value_type,
      typename T::reference_type
    >,
    public ::pressio::impl::AllocTrait<is_static>
{
  using layout_type       = typename T::traits::array_layout;
  using execution_space   = typename T::traits::execution_space;
  using memory_space      = typename T::traits::memory_space;
  using device_type       = typename T::traits::device_type;
  using memory_traits     = typename T::traits::memory_traits;
  using host_mirror_space = typename T::traits::host_mirror_space;
  using host_mirror_t     = typename T::host_mirror_type;

  static constexpr bool has_host_execution_space =
    (false
     #ifdef KOKKOS_ENABLE_SERIAL
     || std::is_same<execution_space, Kokkos::Serial>::value
     #endif
     #ifdef KOKKOS_ENABLE_OPENMP
     || std::is_same<execution_space, Kokkos::OpenMP>::value
     #endif
     );
};
#endif // PRESSIO_ENABLE_TPL_KOKKOS

//*******************************
// Eigen vector helpers
//*******************************
#ifdef PRESSIO_ENABLE_TPL_EIGEN

template <typename T, typename enabled = void> struct EigenVectorIdentifier {};

#define _IMPL_EIGEN_VECTOR_VECTOR_IDENTIFIER(is_vector, id) \
template <typename T> \
struct EigenVectorIdentifier< \
  T, \
  mpl::enable_if_t<::pressio::is_vector<T>::value> \
> { \
  static constexpr VectorIdentifier vector_identifier = VectorIdentifier::id; \
};

_IMPL_EIGEN_VECTOR_VECTOR_IDENTIFIER(is_static_row_vector_eigen,     EigenRowStatic)
_IMPL_EIGEN_VECTOR_VECTOR_IDENTIFIER(is_static_column_vector_eigen,  EigenColStatic)
_IMPL_EIGEN_VECTOR_VECTOR_IDENTIFIER(is_dynamic_row_vector_eigen,    EigenRowDynamic)
_IMPL_EIGEN_VECTOR_VECTOR_IDENTIFIER(is_dynamic_column_vector_eigen, EigenColDynamic)

template <typename T>
using EigenVectorAllocTrait = ::pressio::impl::AllocTrait<
  is_static_vector_eigen<T>::value
>;

template <typename T>
using EigenMatrixAllocTrait = ::pressio::impl::AllocTrait<
  T::RowsAtCompileTime != Eigen::Dynamic &&
  T::ColsAtCompileTime != Eigen::Dynamic
>;

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


}} // namespace pressio::impl
#endif  // TYPE_TRAITS_TRAITS_TPL_HPP_
