/*
//@HEADER
// ************************************************************************
//
// traits_vector.hpp
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

#ifndef TYPE_TRAITS_TRAITS_VECTOR_HPP_
#define TYPE_TRAITS_TRAITS_VECTOR_HPP_

namespace pressio{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
//*******************************
// Eigen STATIC ROW vector
//*******************************
template <typename T>
struct Traits<
  T,
  mpl::enable_if_t<
    is_static_row_vector_eigen<T>::value
    >
  >
  : public ContainersSharedTraits<PackageIdentifier::Eigen, true, 1>
{

  static constexpr VectorIdentifier vector_identifier = VectorIdentifier::EigenRowStatic;
  static constexpr bool is_static = true;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_type   = typename T::Scalar;
  using ordinal_type  = typename T::StorageIndex;
  using size_type     = ordinal_type;
  using reference_type = scalar_type &;
  using const_reference_type = scalar_type const &;
};

//*******************************
// Eigen STATIC COLUMN vector
//*******************************
template <typename T>
struct Traits<
  T,
  mpl::enable_if_t<
    is_static_column_vector_eigen<T>::value
    >
  >
  : public ContainersSharedTraits<PackageIdentifier::Eigen, true, 1>
{

  static constexpr VectorIdentifier vector_identifier = VectorIdentifier::EigenColStatic;
  static constexpr bool is_static = true;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_type   = typename T::Scalar;
  using ordinal_type  = typename T::StorageIndex;
  using size_type     = ordinal_type;
  using reference_type = scalar_type &;
  using const_reference_type = scalar_type const &;
};

//*******************************
// Eigen DYNAMIC ROW vector
//*******************************
template <typename T>
struct Traits<
  T,
  mpl::enable_if_t<
    is_dynamic_row_vector_eigen<T>::value
    >
  >
  : public ContainersSharedTraits<PackageIdentifier::Eigen, true, 1>
{

  static constexpr VectorIdentifier vector_identifier = VectorIdentifier::EigenRowDynamic;
  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_type   = typename T::Scalar;
  using ordinal_type  = typename T::StorageIndex;
  using size_type     = ordinal_type;
  using reference_type = scalar_type &;
  using const_reference_type = scalar_type const &;
};

//*******************************
// Eigen DYNAMIC COLUMN vector
//*******************************
template <typename T>
struct Traits<
  T,
  ::pressio::mpl::enable_if_t<
    is_dynamic_column_vector_eigen<T>::value
    >
  >
  : public ContainersSharedTraits<PackageIdentifier::Eigen, true, 1>
{

  static constexpr VectorIdentifier vector_identifier = VectorIdentifier::EigenColDynamic;
  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_type   = typename T::Scalar;
  using ordinal_type  = typename T::StorageIndex;
  using size_type     = ordinal_type;
  using reference_type = scalar_type &;
  using const_reference_type = scalar_type const &;
};
#endif //PRESSIO_ENABLE_TPL_EIGEN


//*******************************
// Kokkos view, rank==1
//*******************************
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename T>
struct Traits<
  T,
  ::pressio::mpl::enable_if_t<
    is_vector_kokkos<T>::value
    >
  >
  : public ContainersSharedTraits<PackageIdentifier::Kokkos, true, 1>
{

  // static view if the number of runtime determined dimensions == 0
  static constexpr bool is_static = T::traits::rank_dynamic==0;
  static constexpr bool is_dynamic  = !is_static;

  static constexpr VectorIdentifier
  vector_identifier =
    is_static ? VectorIdentifier::KokkosStatic :
    VectorIdentifier::KokkosDynamic;

  using scalar_type      = typename T::traits::value_type;
  using layout_type      = typename T::traits::array_layout;
  using ordinal_type     = typename T::traits::size_type;
  using size_type        = typename T::traits::size_type;
  using execution_space  = typename T::traits::execution_space;
  using memory_space     = typename T::traits::memory_space;
  using device_type      = typename T::traits::device_type;
  using memory_traits    = typename T::traits::memory_traits;
  using host_mirror_space= typename T::traits::host_mirror_space;
  using host_mirror_type = typename T::host_mirror_type;

  using reference_type = typename T::reference_type;

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
#endif

//*******************************
// for tpetra vector
//*******************************
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename T>
struct Traits<
  T,
  mpl::enable_if_t<
    is_vector_tpetra<T>::value
    >
  >
  : public ContainersSharedTraits<PackageIdentifier::Trilinos, false, 1>
{

  static constexpr VectorIdentifier vector_identifier = VectorIdentifier::Tpetra;
  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_type   = typename T::impl_scalar_type;
  using local_ordinal_type  = typename T::local_ordinal_type;
  using global_ordinal_type = typename T::global_ordinal_type;
  using data_map_type   = typename T::map_type;
  using size_type     = global_ordinal_type;
  // using const_data_return_t = T const *;
  // using data_return_t = T *;

  // node is a Tpetra concept, defined as:
  // node_type = ::Kokkos::Compat::KokkosDeviceWrapperNode<execution_space>;
  // where memory space is taken from the execution_space
  ///
  using node_type = typename T::node_type;
  using dual_view_type = typename T::dual_view_type;
  // device_type is just an (execution space, memory space) pair.
  // defined as: Kokkos::Device<execution_space, memory_space>
  // so from the device we can get the device execution and memory space
  using device_type = typename T::device_type;
  using device_mem_space_type = typename device_type::memory_space;
  using device_exec_space_type = typename device_type::execution_space;
  using host_mem_space_type = typename Kokkos::HostSpace::memory_space;
  using host_exec_space_type = typename Kokkos::HostSpace::execution_space;

  using dot_type = typename T::dot_type;
  using mag_type = typename T::mag_type;
  using communicator_type = decltype(std::declval<data_map_type>().getComm());
};
#endif

// *******************************
// for epetra vector
// *******************************
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename T>
struct Traits<
  T,
  mpl::enable_if_t<
    is_vector_epetra<T>::value
    >
  >
  : public ContainersSharedTraits<PackageIdentifier::Trilinos, false, 1>
{

  static constexpr VectorIdentifier vector_identifier = VectorIdentifier::Epetra;
  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_type  = double;
  using local_ordinal_type  = int;
  using global_ordinal_type = int;
  using size_type    = global_ordinal_type;
  using data_map_type  = Epetra_BlockMap;
  using communicator_type   = Epetra_Comm;
};
#endif

//*******************************
// for block tpetra vector
//*******************************
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename T>
struct Traits<
  T,
  mpl::enable_if_t<
    is_vector_tpetra_block<T>::value
    >
  >
  : public ContainersSharedTraits<PackageIdentifier::Trilinos, false, 1>
{

  static constexpr VectorIdentifier vector_identifier = VectorIdentifier::TpetraBlock;
  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_type = typename T::impl_scalar_type;
  using local_ordinal_type = typename T::local_ordinal_type;
  using global_ordinal_type = typename T::global_ordinal_type;
  using data_map_type = typename T::map_type;
  using size_type    = global_ordinal_type;
  using mag_type = scalar_type;

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
  using communicator_type = decltype(std::declval<data_map_type>().getComm());
};
#endif

//*******************************
// for teuchos serial dense vec
//*******************************
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename T>
struct Traits<
  T,
  mpl::enable_if_t<
    is_dense_vector_teuchos<T>::value
    >
  >
  : public ContainersSharedTraits<PackageIdentifier::Trilinos, true, 1>
{

  static constexpr VectorIdentifier vector_identifier = VectorIdentifier::TeuchosSerialDense;
  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_type  = typename T::scalarType;
  using ordinal_type = typename T::ordinalType;
  using size_type    = ordinal_type;
  using reference_type = scalar_type &;
  using const_reference_type = scalar_type const &;
};
#endif

}
#endif  // TYPE_TRAITS_TRAITS_VECTOR_HPP_
