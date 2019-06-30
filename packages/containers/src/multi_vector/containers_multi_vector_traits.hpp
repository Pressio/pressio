
#ifndef CONTAINERS_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_
#define CONTAINERS_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_

#include "../containers_fwd.hpp"
#include "../containers_shared_traits.hpp"
#include "./meta/containers_native_eigen_multi_vector_meta.hpp"
#include "./meta/containers_native_epetra_multi_vector_meta.hpp"
#include "./meta/containers_native_tpetra_multi_vector_meta.hpp"
#include "./meta/containers_native_tpetra_block_multi_vector_meta.hpp"

namespace rompp{ namespace containers{ namespace details{

#ifdef HAVE_TRILINOS

/********************************
an arbitrary multi vector is one
for which a user must provide ops
*******************************/
template <typename wrapped_type>
struct traits<
  MultiVector<
    wrapped_type,
    mpl::enable_if_t<
      !containers::meta::is_dynamic_multi_vector_eigen<wrapped_type>::value and
      !containers::meta::is_multi_vector_epetra<wrapped_type>::value and
      !containers::meta::is_multi_vector_tpetra_block<wrapped_type>::value and
      !containers::meta::is_multi_vector_tpetra<wrapped_type>::value
      >
    >
  > {

  using wrapped_t = wrapped_type;
  using derived_t = MultiVector<wrapped_t>;

  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Arbitrary;

  static constexpr WrappedPackageIdentifier
  wrapped_package_identifier = WrappedPackageIdentifier::Arbitrary;

  static constexpr bool is_vector = false;
  static constexpr bool is_matrix = false;
  static constexpr bool is_multi_vector = true;

  // by default, any container is not admissible to expr templates
  // the ones that are, will overwrite this
  static constexpr bool is_admissible_for_expression_templates = false;
};

//*******************************
// for epetra multivector
//*******************************
template<typename wrapped_type>
struct traits<MultiVector<wrapped_type,
      typename std::enable_if<
       meta::is_multi_vector_epetra<wrapped_type
      >::value>::type>
     >
  : public containers_shared_traits<MultiVector<wrapped_type>,
				    wrapped_type,
				    false, false, true,
				    WrappedPackageIdentifier::Trilinos,
				    false, false>
{
  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Epetra;

  using scalar_t = default_types::epetra_scalar_t;
  using local_ordinal_t = containers::default_types::epetra_lo_t;
  using global_ordinal_t = containers::default_types::epetra_go_t1;
  using data_map_t = Epetra_BlockMap;
  using communicator_t = Epetra_Comm;
};
#endif


#ifdef HAVE_TRILINOS
//*******************************
// for tpetra multivector
//*******************************
template<typename wrapped_type>
struct traits<MultiVector<wrapped_type,
      typename std::enable_if<
       meta::is_multi_vector_tpetra<wrapped_type
      >::value>::type>
     >
  : public containers_shared_traits<MultiVector<wrapped_type>,
				    wrapped_type,
				    false, false, true,
				    WrappedPackageIdentifier::Trilinos,
				    false, false>
{
  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Tpetra;

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


#ifdef HAVE_TRILINOS
//*******************************
// for block tpetra multivector
//*******************************
template<typename wrapped_type>
struct traits<
  MultiVector<
    wrapped_type,
    ::rompp::mpl::enable_if_t<
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
				    false, false>
{
  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::TpetraBlock;

  using scalar_t = typename wrapped_type::impl_scalar_type;
  using local_ordinal_t = typename wrapped_type::local_ordinal_type;
  using global_ordinal_t = typename wrapped_type::global_ordinal_type;
  using data_map_t = typename wrapped_type::map_type;

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
// for eigen multivector
//*******************************
template<typename wrapped_type>
struct traits<MultiVector<wrapped_type,
      typename std::enable_if<
       meta::is_dynamic_multi_vector_eigen<wrapped_type
      >::value>::type>
     >
  : public containers_shared_traits<MultiVector<wrapped_type>,
            wrapped_type,
            false, false, true,
            WrappedPackageIdentifier::Eigen, true,
	   ( wrapped_type::RowsAtCompileTime != Eigen::Dynamic &&
	     wrapped_type::ColsAtCompileTime != Eigen::Dynamic )>
{
  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Eigen;

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;

  // static constexpr bool is_static =
  //   ( wrapped_type::RowsAtCompileTime != Eigen::Dynamic &&
  //     wrapped_type::ColsAtCompileTime != Eigen::Dynamic );
  // static constexpr bool is_dynamic = !is_static;
};

}}}//end namespace rompp::containers::details
#endif
