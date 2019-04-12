
#ifndef CORE_VECTOR_VECTOR_TRAITS_HPP_
#define CORE_VECTOR_VECTOR_TRAITS_HPP_

#include "../core_forward_declarations.hpp"
#include "../meta/core_native_vector_meta.hpp"
// #include "../meta/core_meta_detect_typedefs.hpp"
// #include "../meta/core_meta_detect_operators.hpp"
#include "../core_shared_traits.hpp"

namespace rompp{ namespace core{ namespace details{

//*******************************
// Eigen STATIC ROW vector
//*******************************
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     typename
		     std::enable_if<
		       core::meta::is_static_row_vector_eigen<
			 wrapped_type
			 >::value
		       >::type
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Eigen,
				    true, true>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::EigenRowStatic;

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;
};


//*******************************
// Eigen STATIC COLUMN vector
//*******************************
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     typename
		     std::enable_if<
		       core::meta::is_static_column_vector_eigen<
			 wrapped_type
			 >::value
		       >::type
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Eigen,
				    true, true>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::EigenColStatic;

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;
};


//*******************************
// Eigen DYNAMIC ROW vector
//*******************************
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     typename
		     std::enable_if<
		       core::meta::is_dynamic_row_vector_eigen<
			 wrapped_type
			 >::value
		       >::type
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Eigen,
				    true, false>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::EigenRowDynamic;

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;
};


//*******************************
// Eigen DYNAMIC COLUMN vector
//*******************************
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     typename
		     std::enable_if<
		       core::meta::is_dynamic_column_vector_eigen<
			 wrapped_type
			 >::value
		       >::type
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Eigen,
				    true, false>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::EigenColDynamic;

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;
};


//*******************************
// Armadillo column vector
//*******************************
#ifdef HAVE_ARMADILLO
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     ::rompp::mpl::enable_if_t<
		       core::meta::is_column_vector_armadillo<
			 wrapped_type>::value
		       >
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Armadillo,
				    true, false>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier =
    WrappedVectorIdentifier::ArmadilloCol;

  using scalar_t = typename wrapped_type::elem_type;
  using ordinal_t = unsigned long;
};
#endif


//*******************************
// Armadillo row vector
//*******************************
#ifdef HAVE_ARMADILLO
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     ::rompp::mpl::enable_if_t<
		       core::meta::is_row_vector_armadillo<
			 wrapped_type>::value
		       >
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Armadillo,
				    true, false>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier =
    WrappedVectorIdentifier::ArmadilloRow;

  using scalar_t = typename wrapped_type::elem_type;
  using ordinal_t = unsigned long;
};
#endif


//*******************************
// Blaze dynamic vector
//*******************************
#ifdef HAVE_BLAZE
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     ::rompp::mpl::enable_if_t<
		       core::meta::is_dynamic_vector_blaze<
			 wrapped_type>::value
		       >
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Blaze,
				    true, false>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::BlazeDynamic;

  using scalar_t = typename wrapped_type::ElementType;
  using ordinal_t = unsigned long;
};
#endif


//*******************************
// for teuchos serial dense vec
//*******************************
#ifdef HAVE_TRILINOS
template<typename wrapped_type>
struct traits<Vector<wrapped_type,
	  typename
	  std::enable_if<
	    core::meta::is_dense_vector_teuchos<
	      wrapped_type>::value
	    >::type
	  >
	>
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
			       WrappedPackageIdentifier::Trilinos,
				    true, false>{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::TeuchosSerialDense;

  using scalar_t = typename wrapped_type::scalarType;
  using ordinal_t = typename wrapped_type::ordinalType;
};
#endif


//*******************************
// for epetra vector
//*******************************
#ifdef HAVE_TRILINOS
template<typename wrapped_type>
struct traits<Vector<wrapped_type,
	  typename
	  std::enable_if<
	    core::meta::is_vector_epetra<
	      wrapped_type>::value
	    >::type
	  >
	>
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
			       WrappedPackageIdentifier::Trilinos,
				    false, false>{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::Epetra;

  using scalar_t = default_types::epetra_scalar_t;
  using local_ordinal_t = core::default_types::epetra_lo_t;
  using global_ordinal_t = core::default_types::epetra_go_t1;
  using data_map_t = Epetra_BlockMap;
  using communicator_t = Epetra_Comm;
};
#endif


//*******************************
// for tpetra vector
//*******************************
#ifdef HAVE_TRILINOS
template<typename wrapped_type>
struct traits<Vector<wrapped_type,
	  typename
	  std::enable_if<
	    core::meta::is_vector_tpetra<
	      wrapped_type>::value
	    >::type
	  >
	>
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
			       WrappedPackageIdentifier::Trilinos,
				    false, false>{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::Tpetra;

  using scalar_t = typename wrapped_type::impl_scalar_type;
  using local_ordinal_t = typename wrapped_type::local_ordinal_type;
  using global_ordinal_t = typename wrapped_type::global_ordinal_type;
  using data_map_t = typename wrapped_type::map_type;

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
#endif


//*******************************
// Kokkos vector
//*******************************
#ifdef HAVE_TRILINOS
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
	  ::rompp::mpl::enable_if_t<
	    core::meta::is_vector_kokkos<
	      wrapped_type>::value
	    >
	  >
	>
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Kokkos,
				    true, wrapped_type::traits::rank_dynamic==0>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::Kokkos;

  using scalar_t = typename wrapped_type::traits::value_type;
  using ordinal_t = int;

  using  execution_space = typename wrapped_type::traits::execution_space;
  using  memory_space = typename wrapped_type::traits::memory_space;
  using  device_type = typename wrapped_type::traits::device_type;
  using  memory_traits = typename wrapped_type::traits::memory_traits;
  using  host_mirror_space = typename wrapped_type::traits::host_mirror_space;
};
#endif


//*******************************
// for a std vector
//*******************************
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
    typename
    std::enable_if<
      core::meta::is_vector_stdlib<
	wrapped_type>::value
      >::type
    >
  >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
			       WrappedPackageIdentifier::CppStdLib,
				    true, false>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::CppStdLib;

  using scalar_t = typename wrapped_type::value_type;
  using ordinal_t = core::default_types::local_ordinal_t;
};


}}}//end namespace rompp::core::details
#endif
