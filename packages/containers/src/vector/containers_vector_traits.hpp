/*
//@HEADER
// ************************************************************************
//
// containers_vector_traits.hpp
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

#ifndef CONTAINERS_VECTOR_VECTOR_TRAITS_HPP_
#define CONTAINERS_VECTOR_VECTOR_TRAITS_HPP_

namespace pressio{ namespace containers{ namespace details{

/********************************
an arbitrary vector is one
for which a user must provide ops
*******************************/
template <typename wrapped_type>
struct traits<
  Vector<
    wrapped_type,
    mpl::enable_if_t<
      !containers::meta::is_vector_eigen<wrapped_type>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      and !containers::meta::is_array_pybind<wrapped_type>::value
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
      and !containers::meta::is_vector_kokkos<wrapped_type>::value
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
      and !containers::meta::is_vector_epetra<wrapped_type>::value
      and !containers::meta::is_dense_vector_teuchos<wrapped_type>::value
      and !containers::meta::is_vector_tpetra_block<wrapped_type>::value
      and !containers::meta::is_vector_tpetra<wrapped_type>::value
#endif
      >
    >
  >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Undefined,
            false>
{

  using wrapped_t = wrapped_type;
  using derived_t = Vector<wrapped_t>;
  using scalar_t  = typename wrapped_type::value_type;
  using value_t   = typename wrapped_type::value_type;
  using size_t    = typename wrapped_type::size_type;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::Arbitrary;

  static constexpr WrappedPackageIdentifier
  wrapped_package_identifier = WrappedPackageIdentifier::Arbitrary;

  static constexpr bool is_vector = true;
  static constexpr bool is_matrix = false;
  static constexpr bool is_multi_vector = false;
};


//*******************************
// Eigen STATIC ROW vector
//*******************************
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     typename
		     std::enable_if<
		       containers::meta::is_static_row_vector_eigen<
			 wrapped_type
			 >::value
		       >::type
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Eigen,
				    true>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::EigenRowStatic;

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;
  using size_t    = ordinal_t;
  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

  static constexpr bool is_static = true;
  static constexpr bool is_dynamic  = !is_static;

  using reference_t = scalar_t &;
  using const_reference_t = scalar_t const &;
};


//*******************************
// Eigen STATIC COLUMN vector
//*******************************
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     typename
		     std::enable_if<
		       containers::meta::is_static_column_vector_eigen<
			 wrapped_type
			 >::value
		       >::type
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Eigen,
				    true>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::EigenColStatic;

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;
  using size_t    = ordinal_t;
  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

  static constexpr bool is_static = true;
  static constexpr bool is_dynamic  = !is_static;

  using reference_t = scalar_t &;
  using const_reference_t = scalar_t const &;
};


//*******************************
// Eigen DYNAMIC ROW vector
//*******************************
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     typename
		     std::enable_if<
		       containers::meta::is_dynamic_row_vector_eigen<
			 wrapped_type
			 >::value
		       >::type
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Eigen,
				    true>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::EigenRowDynamic;

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;
  using size_t    = ordinal_t;
  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

  using reference_t = scalar_t &;
  using const_reference_t = scalar_t const &;

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using span_ret_t	 = expressions::SpanExpr<Vector<wrapped_type>>;
  using span_const_ret_t = expressions::SpanExpr< const Vector<wrapped_type>>;
};


//*******************************
// Eigen DYNAMIC COLUMN vector
//*******************************
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     typename
		     std::enable_if<
		       containers::meta::is_dynamic_column_vector_eigen<
			 wrapped_type
			 >::value
		       >::type
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Eigen,
				    true>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::EigenColDynamic;

  using scalar_t	 = typename wrapped_type::Scalar;
  using ordinal_t	 = int;
  using size_t    = ordinal_t;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;
  using reference_t = scalar_t &;
  using const_reference_t = scalar_t const &;

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  using span_ret_t	 = expressions::SpanExpr<Vector<wrapped_type>>;
  using span_const_ret_t = expressions::SpanExpr< const Vector<wrapped_type>>;
};


//*******************************
// for teuchos serial dense vec
//*******************************
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename wrapped_type>
struct traits<Vector<wrapped_type,
	  typename
	  std::enable_if<
	    containers::meta::is_dense_vector_teuchos<
	      wrapped_type>::value
	    >::type
	  >
	>
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
			       WrappedPackageIdentifier::Trilinos,
				    true>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::TeuchosSerialDense;

  using scalar_t = typename wrapped_type::scalarType;
  using ordinal_t = typename wrapped_type::ordinalType;
  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;
  using size_t    = int;
  using reference_t = scalar_t &;
  using const_reference_t = scalar_t const &;
  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;
};
#endif


//*******************************
// for epetra vector
//*******************************
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename wrapped_type>
struct traits<Vector<wrapped_type,
	  typename
	  std::enable_if<
	    containers::meta::is_vector_epetra<
	      wrapped_type>::value
	    >::type
	  >
	>
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
			       WrappedPackageIdentifier::Trilinos,
				    false>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::Epetra;

  using scalar_t = double;
  using local_ordinal_t = int;
  using global_ordinal_t = int;
  using size_t    = global_ordinal_t;
  using data_map_t = Epetra_BlockMap;
  using communicator_t = Epetra_Comm;
  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;
  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

};
#endif


//*******************************
// for tpetra vector
//*******************************
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename wrapped_type>
struct traits<Vector<wrapped_type,
	  typename
	  std::enable_if<
	    containers::meta::is_vector_tpetra<
	      wrapped_type>::value
	    >::type
	  >
	>
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
			       WrappedPackageIdentifier::Trilinos,
				    false>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::Tpetra;

  using scalar_t = typename wrapped_type::impl_scalar_type;
  using local_ordinal_t = typename wrapped_type::local_ordinal_type;
  using global_ordinal_t = typename wrapped_type::global_ordinal_type;
  using data_map_t = typename wrapped_type::map_type;
  using size_t    = global_ordinal_t;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

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

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;
};
#endif


//*******************************
// Kokkos vector
//*******************************
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
	  ::pressio::mpl::enable_if_t<
	    containers::meta::is_vector_kokkos<
	      wrapped_type>::value
	    >
	  >
	>
  : public containers_shared_traits<
  Vector<wrapped_type>,
  wrapped_type,
  true, false, false,
  WrappedPackageIdentifier::Kokkos,
  true //true because kokkos is for shared mem
  >
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::Kokkos;

  using scalar_t	   = typename wrapped_type::traits::value_type;
  using layout		   = typename wrapped_type::traits::array_layout;
  using ordinal_t	   = typename wrapped_type::traits::size_type;
  using size_t    = ordinal_t;
  using execution_space    = typename wrapped_type::traits::execution_space;
  using memory_space	   = typename wrapped_type::traits::memory_space;
  using device_type	   = typename wrapped_type::traits::device_type;
  using device_t	   = typename wrapped_type::traits::device_type;
  using memory_traits	   = typename wrapped_type::traits::memory_traits;
  using host_mirror_space  = typename wrapped_type::traits::host_mirror_space;
  using host_mirror_t      = typename wrapped_type::host_mirror_type;

  using span_ret_t	   = expressions::SpanExpr<Vector<wrapped_type>>;
  using span_const_ret_t   = expressions::SpanExpr<const Vector<wrapped_type>>;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

  using reference_t = scalar_t &;
  using const_reference_t = scalar_t const &;

  // static view if the number of runtime determined dimensions == 0
  static constexpr bool is_static = wrapped_type::traits::rank_dynamic==0;
  static constexpr bool is_dynamic  = !is_static;

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
// for block tpetra vector
//*******************************
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename wrapped_type>
struct traits<
  Vector<wrapped_type,
	 typename
	 std::enable_if<
	   containers::meta::is_vector_tpetra_block<
	     wrapped_type
	     >::value
	   >::type
	 >
  >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
			       WrappedPackageIdentifier::Trilinos,
				    false>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::TpetraBlock;

  using scalar_t = typename wrapped_type::impl_scalar_type;
  using local_ordinal_t = typename wrapped_type::local_ordinal_type;
  using global_ordinal_t = typename wrapped_type::global_ordinal_type;
  using data_map_t = typename wrapped_type::map_type;
  using size_t    = global_ordinal_t;

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;
  using mag_t = scalar_t;

  /* node is a Tpetra concept, defined as:
   * node_type = ::Kokkos::Compat::KokkosDeviceWrapperNode<execution_space>;
   * where memory space is taken from the execution_space
   */
  using node_t = typename wrapped_type::node_type;

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

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
// Pybind array
//*******************************
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename wrapped_type>
struct traits<
  Vector<
    wrapped_type,
    mpl::enable_if_t<
      containers::meta::is_array_pybind<wrapped_type>::value
      >
    >
  >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Pybind,
				    true>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::Pybind;

  using scalar_t	 = typename wrapped_type::value_type;
  using ordinal_t	 = std::size_t;
  using size_t		 = ordinal_t;

  using mut_proxy_t = decltype( std::declval<wrapped_type &>().mutable_unchecked() );
  using proxy_t	    = decltype( std::declval<const wrapped_type &>().unchecked() );

  using const_data_return_t = wrapped_type const *;
  using data_return_t = wrapped_type *;

  using reference_t = scalar_t &;
  using const_reference_t = scalar_t const &;

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;

  // using span_ret_t	 = expressions::SpanExpr<Vector<wrapped_type>>;
  // using span_const_ret_t = expressions::SpanExpr< const Vector<wrapped_type>>;
};
#endif

}}}//end namespace pressio::containers::details
#endif
