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

#ifndef CONTAINERS_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_
#define CONTAINERS_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_

namespace pressio{ namespace containers{ namespace details{

/********************************
an arbitrary multi vector is one
for which a user must provide ops
*******************************/
template <typename wrapped_type>
struct traits<
  MultiVector<
    wrapped_type,
    mpl::enable_if_t<
      !containers::meta::is_dynamic_multi_vector_eigen<wrapped_type>::value
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
      and
      !containers::meta::is_multi_vector_epetra<wrapped_type>::value and
      !containers::meta::is_multi_vector_tpetra_block<wrapped_type>::value and
      !containers::meta::is_multi_vector_tpetra<wrapped_type>::value
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
      and
      !containers::meta::is_multi_vector_kokkos<wrapped_type>::value
#endif
      >
    >
  >
  : public containers_shared_traits<MultiVector<wrapped_type>,
				    wrapped_type,
				    false, false, true,
				    WrappedPackageIdentifier::Arbitrary,
				    false>
{

  using wrapped_t = wrapped_type;
  using derived_t = MultiVector<wrapped_t>;
  using scalar_t  = typename wrapped_type::value_type;
  using value_t   = typename wrapped_type::value_type;
  using size_t   = typename wrapped_type::size_type;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;
  using data_cp_return_t = wrapped_type;

  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Arbitrary;

  static constexpr WrappedPackageIdentifier
  wrapped_package_identifier = WrappedPackageIdentifier::Arbitrary;

  static constexpr bool is_vector = false;
  static constexpr bool is_matrix = false;
  static constexpr bool is_multi_vector = true;
};


#ifdef PRESSIO_ENABLE_TPL_TRILINOS
//*******************************
// for epetra multivector
//*******************************
template<typename wrapped_type>
struct traits<
  MultiVector<
    wrapped_type,
    ::pressio::mpl::enable_if_t<
      meta::is_multi_vector_epetra<
	wrapped_type
	>::value
      >
    >
  >
  : public containers_shared_traits<MultiVector<wrapped_type>,
				    wrapped_type,
				    false, false, true,
				    WrappedPackageIdentifier::Trilinos,
				    false>
{
  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Epetra;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;
  using data_cp_return_t = wrapped_type;

  using scalar_t = double;
  using local_ordinal_t = int;
  using global_ordinal_t = int;
  using size_t    = global_ordinal_t;

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using data_map_t = Epetra_BlockMap;
  using communicator_t = Epetra_Comm;
};


//*******************************
// for tpetra multivector
//*******************************
template<typename wrapped_type>
struct traits<
  MultiVector<
    wrapped_type,
    ::pressio::mpl::enable_if_t<
      meta::is_multi_vector_tpetra<
	wrapped_type
	>::value
      >
    >
  >
  : public containers_shared_traits<MultiVector<wrapped_type>,
				    wrapped_type,
				    false, false, true,
				    WrappedPackageIdentifier::Trilinos,
				    false>
{
  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Tpetra;

  using scalar_t = typename wrapped_type::impl_scalar_type;
  using local_ordinal_t = typename wrapped_type::local_ordinal_type;
  using global_ordinal_t = typename wrapped_type::global_ordinal_type;
  using data_map_t = typename wrapped_type::map_type;
  using size_t    = global_ordinal_t;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;
  using data_cp_return_t = wrapped_type;

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  /* node is a Tpetra concept, defined as:
   * node_type = ::Kokkos::Compat::KokkosDeviceWrapperNode<execution_space>;
   * where memory space is taken from the execution_space
   */
  using node_t = typename wrapped_type::node_type;
  using dual_view_t = typename wrapped_type::dual_view_type;
  // device_type is just an (execution space, memory space) pair.
  // defined as: Kokkos::Device<execution_space, memory_space>
  // so from the device we can get the device execution and memory space
  using device_t = typename wrapped_type::device_type;
  using device_mem_space_t = typename device_t::memory_space;
  using device_exec_space_t = typename device_t::execution_space;
  // store types for host
  using host_mem_space_t = typename Kokkos::HostSpace::memory_space;
  using host_exec_space_t = typename Kokkos::HostSpace::execution_space;

  using dot_t = typename wrapped_type::dot_type;
  using mag_t = typename wrapped_type::mag_type;
  using communicator_t = decltype(std::declval<data_map_t>().getComm());
};


//*******************************
// for block tpetra multivector
//*******************************
template<typename wrapped_type>
struct traits<
  MultiVector<
    wrapped_type,
    ::pressio::mpl::enable_if_t<
      meta::is_multi_vector_tpetra_block<
	wrapped_type
	>::value
      >
    >
  >
  : public containers_shared_traits<MultiVector<wrapped_type>,
				    wrapped_type,
				    false, false, true,
				    WrappedPackageIdentifier::Trilinos,
				    false>
{
  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::TpetraBlock;

  using scalar_t = typename wrapped_type::impl_scalar_type;
  using local_ordinal_t = typename wrapped_type::local_ordinal_type;
  using global_ordinal_t = typename wrapped_type::global_ordinal_type;
  using data_map_t = typename wrapped_type::map_type;
  using size_t    = global_ordinal_t;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;
  using data_cp_return_t = wrapped_type;

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  /* node is a Tpetra concept, defined as:
   * node_type = ::Kokkos::Compat::KokkosDeviceWrapperNode<execution_space>;
   * where memory space is taken from the execution_space
   */
  using node_t = typename wrapped_type::node_type;

  // device_type is just an (execution space, memory space) pair.
  // defined as: Kokkos::Device<execution_space, memory_space>
  // so from the device we can get the device execution and memory space
  using device_t = typename wrapped_type::device_type;
  using device_mem_space_t = typename device_t::memory_space;
  using device_exec_space_t = typename device_t::execution_space;
  // store types for host
  using host_mem_space_t = typename Kokkos::HostSpace::memory_space;
  using host_exec_space_t = typename Kokkos::HostSpace::execution_space;

  using communicator_t = decltype(std::declval<data_map_t>().getComm());
};
#endif


//*******************************
// for eigen dynamic multivector
//*******************************
template<typename wrapped_type>
struct traits<
  MultiVector<
    wrapped_type,
    ::pressio::mpl::enable_if_t<
      meta::is_dynamic_multi_vector_eigen<wrapped_type>::value
      >
    >
  >
  : public containers_shared_traits<
  MultiVector<wrapped_type>,
  wrapped_type, false, false, true,
  WrappedPackageIdentifier::Eigen, true
  >
{
  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Eigen;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;
  using data_cp_return_t = wrapped_type;

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_t  = typename wrapped_type::Scalar;
  using ordinal_t = typename wrapped_type::StorageIndex;
  using size_t    = ordinal_t;

  static constexpr bool is_admissible_for_expression_templates = true;
};


//*******************************
// Kokkos multi vector
//*******************************
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename wrapped_type>
struct traits<
  MultiVector<
    wrapped_type,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_multi_vector_kokkos<
	wrapped_type
	>::value
      >
    >
  >
  : public containers_shared_traits<
  MultiVector<wrapped_type>,
  wrapped_type,
  false, false, true,
  WrappedPackageIdentifier::Kokkos,
  true //true because kokkos is for shared mem
  >
{

  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Kokkos;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;
  using data_cp_return_t = wrapped_type;

  // static view if the number of runtime determined dimensions == 0
  static constexpr bool is_static = wrapped_type::traits::rank_dynamic==0;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_t	  = typename wrapped_type::traits::value_type;
  using layout		  = typename wrapped_type::traits::array_layout;
  using ordinal_t	  = typename wrapped_type::traits::size_type;
  using size_t		  = ordinal_t;
  using execution_space	  = typename wrapped_type::traits::execution_space;
  using memory_space	  = typename wrapped_type::traits::memory_space;
  using device_type	  = typename wrapped_type::traits::device_type;
  using memory_traits	  = typename wrapped_type::traits::memory_traits;
  using host_mirror_space = typename wrapped_type::traits::host_mirror_space;
  using host_mirror_t     = typename wrapped_type::host_mirror_type;
};
#endif

}}}//end namespace pressio::containers::details
#endif
