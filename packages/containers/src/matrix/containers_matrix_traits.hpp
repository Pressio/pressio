
#ifndef CONTAINERS_MATRIX_MATRIX_TRAITS_HPP_
#define CONTAINERS_MATRIX_MATRIX_TRAITS_HPP_

#include "../containers_fwd.hpp"
#include "../containers_shared_traits.hpp"
#include "./meta/containers_native_eigen_matrix_meta.hpp"
#ifdef HAVE_TRILINOS
#include "./meta/containers_native_trilinos_matrix_meta.hpp"
#include "./meta/containers_native_kokkos_matrix_meta.hpp"
#endif

namespace pressio{ namespace containers{ namespace details{

/********************************
an arbitrary matrix is one
for which a user must provide ops
*******************************/
template <typename wrapped_type>
struct traits<
  Matrix<
    wrapped_type,
    mpl::enable_if_t<
      !containers::meta::is_dense_matrix_eigen<wrapped_type>::value and
      !containers::meta::is_sparse_matrix_eigen<wrapped_type>::value
#ifdef HAVE_TRILINOS
      and
      !containers::meta::is_sparse_matrix_epetra<wrapped_type>::value and
      !containers::meta::is_dense_matrix_epetra<wrapped_type>::value and
      !containers::meta::is_dense_matrix_teuchos<wrapped_type>::value and
      !containers::meta::is_dense_matrix_teuchos_rcp<wrapped_type>::value and
      !containers::meta::is_sparse_matrix_tpetra<wrapped_type>::value and
      !containers::meta::is_sparse_matrix_kokkos<wrapped_type>::value and
      !containers::meta::is_dense_matrix_kokkos<wrapped_type>::value
#endif
      >
    >
  > {

  using wrapped_t = wrapped_type;
  using derived_t = Matrix<wrapped_t>;

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::Arbitrary;

  static constexpr WrappedPackageIdentifier
  wrapped_package_identifier = WrappedPackageIdentifier::Arbitrary;

  static constexpr bool is_vector = false;
  static constexpr bool is_matrix = true;
  static constexpr bool is_multi_vector = false;

  // by default, any container is not admissible to expr templates
  // the ones that are, will overwrite this
  static constexpr bool is_admissible_for_expression_templates = false;
};


//***********************************
// eigen dense matrix
//***********************************
template <typename wrapped_type>
struct traits< Matrix<
    wrapped_type,
    typename
    std::enable_if<
      containers::meta::is_dense_matrix_eigen<
	wrapped_type>::value
      >::type
    >
  >
  : public containers_shared_traits<Matrix<wrapped_type>,
				    wrapped_type, false, true, false,
				    WrappedPackageIdentifier::Eigen,true,
				    ( wrapped_type::RowsAtCompileTime != Eigen::Dynamic &&
				      wrapped_type::ColsAtCompileTime != Eigen::Dynamic )>,
    public matrix_shared_traits<false>
{

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::DenseEigen;

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;

  // static constexpr bool is_static =
  //   ( wrapped_type::RowsAtCompileTime != Eigen::Dynamic &&
  //     wrapped_type::ColsAtCompileTime != Eigen::Dynamic );
};


//***********************************
// eigen sparse matrix
//***********************************
template <typename wrapped_type>
struct traits< Matrix<
    wrapped_type,
    typename
    std::enable_if<
      containers::meta::is_sparse_matrix_eigen<
	wrapped_type
	>::value
      >::type
    >
  >
  : public containers_shared_traits<Matrix<wrapped_type>,
				    wrapped_type, false, true, false,
				    WrappedPackageIdentifier::Eigen, true, false>,
    public matrix_shared_traits<true>
{
  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::SparseEigen;

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = typename wrapped_type::StorageIndex;
  //  ordinal has to be integral and signed
  static_assert( std::is_integral<ordinal_t>::value &&
  		 std::is_signed<ordinal_t>::value,
  "ordinal type for indexing eigen sparse matrix has to be signed");

  static constexpr bool is_row_major = wrapped_type::IsRowMajor;
  static constexpr bool is_col_major = !is_row_major;
};


//***********************************
// epetra sparse distributed matrix
//***********************************
#ifdef HAVE_TRILINOS
template <typename wrapped_type>
struct traits<Matrix
   <wrapped_type,
    typename
    std::enable_if<
      containers::meta::is_sparse_matrix_epetra<
	wrapped_type
	>::value
      >::type
    >
  >
  : public containers_shared_traits<Matrix<wrapped_type>,
				    wrapped_type, false, true, false,
				    WrappedPackageIdentifier::Trilinos,
				    false, false>,
    public matrix_shared_traits<true>
{

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::CrsEpetra;

  using scalar_t = default_types::epetra_scalar_t;
  using local_ordinal_t = containers::default_types::epetra_lo_t;
  using global_ordinal_t = containers::default_types::epetra_go_t1;
  using communicator_t = Epetra_Comm;
  using row_map_t = Epetra_Map;
  using col_map_t = Epetra_Map;
  using range_map_t = Epetra_Map;
  using domain_map_t = Epetra_Map;
  using crs_graph_t = Epetra_CrsGraph;
};
#endif


//**********************************
// for teuchos serial dense matrix
//**********************************
#ifdef HAVE_TRILINOS
template<typename wrapped_type>
struct traits<Matrix<wrapped_type,
	  typename
	  std::enable_if<
	    containers::meta::is_dense_matrix_teuchos<
	      wrapped_type>::value
	    >::type
	  >
	>
  : public containers_shared_traits<Matrix<wrapped_type>,
				    wrapped_type,
				    false, true, false,
				    WrappedPackageIdentifier::Trilinos,
				    true, false>,
    public matrix_shared_traits<false>
{

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::TeuchosSerialDense;

  using scalar_t = typename wrapped_type::scalarType;
  using ordinal_t = typename wrapped_type::ordinalType;
};
#endif


//***********************************
// epetra dense distributed matrix
//***********************************
#ifdef HAVE_TRILINOS
template <typename wrapped_type>
struct traits<Matrix
   <wrapped_type,
    typename
    std::enable_if<
      containers::meta::is_dense_matrix_epetra<
	wrapped_type
	>::value
      >::type
    >
  >
  : public containers_shared_traits<Matrix<wrapped_type>,
				    wrapped_type, false, true, false,
				    WrappedPackageIdentifier::Trilinos,
				    false, false>,
    public matrix_shared_traits<false>
{
  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::DenseEpetra;

  using scalar_t = default_types::epetra_scalar_t;
  using local_ordinal_t = containers::default_types::epetra_lo_t;
  using global_ordinal_t = containers::default_types::epetra_go_t1;
  using communicator_t = Epetra_Comm;
  using row_map_t = Epetra_BlockMap;

  // // a multivector can also be seen as literally a multi-vector.
  // // But here is intended as distributed dense matrix.
  // // So we set actsAsMultiVector = false, actingAsDenseMatrix = true
  // static constexpr int actingAsMultiVector = 0;
  // static constexpr int actingAsDenseMatrix = 1;
};
#endif


#ifdef HAVE_TRILINOS
//*******************************
// for tpetra crs matrix
//*******************************
template<typename wrapped_type>
struct traits<Matrix<wrapped_type,
      typename std::enable_if<
       meta::is_sparse_matrix_tpetra<wrapped_type
      >::value>::type>
     >
  : public containers_shared_traits<Matrix<wrapped_type>,
            wrapped_type, false, true, false,
	    WrappedPackageIdentifier::Trilinos,false, false>,
    public matrix_shared_traits<true>
{
  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::SparseTpetra;

  static constexpr int is_static = 0;

  using scalar_t = typename wrapped_type::impl_scalar_type;
  using local_ordinal_t = typename wrapped_type::local_ordinal_type;
  using global_ordinal_t = typename wrapped_type::global_ordinal_type;
  using row_map_t = typename wrapped_type::map_type;
  using col_map_t = typename wrapped_type::map_type;
  using range_map_t = typename wrapped_type::map_type;
  using domain_map_t = typename wrapped_type::map_type;

  /* node is a Tpetra concept, defined as:
   * node_type = ::Kokkos::Compat::KokkosDeviceWrapperNode<execution_space>;
   * where memory space is taken from the execution_space
   */
  using node_t = typename wrapped_type::node_type;
  //using dual_view_t = typename wrapped_type::dual_view_type;
  // device_type is just an (execution space, memory space) pair.
  // defined as: Kokkos::Device<execution_space, memory_space>
  // so from the device we can get the device execution and memory space
  using device_t = typename wrapped_type::device_type;
  using device_mem_space_t = typename device_t::memory_space;
  using device_exec_space_t = typename device_t::execution_space;
  // store types for host
  using host_mem_space_t = typename Kokkos::HostSpace::memory_space;
  using host_exec_space_t = typename Kokkos::HostSpace::execution_space;

  using mag_t = typename wrapped_type::mag_type;
  using communicator_t = decltype(std::declval<wrapped_type>().getComm());
};
#endif


//*******************************
// Kokkos crs matrix
//*******************************
#ifdef HAVE_TRILINOS
template <typename wrapped_type>
struct traits<
  Matrix<
    wrapped_type,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_sparse_matrix_kokkos<
	wrapped_type
	>::value
      >
    >
  >
  : public containers_shared_traits<
  Matrix<wrapped_type>,
  wrapped_type,
  false, true, false,
  WrappedPackageIdentifier::Kokkos,
  true, //true because kokkos is for shared mem
  // the values of the crs matrix are stored in a 1d dynamic view,
  // so set crs kokkos matrix to be dynamic, but maybe we should decide in
  // a different way
  false
  >
{

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::CrsKokkos;

  using scalar_t	= typename wrapped_type::value_type;
  using ordinal_t	= typename wrapped_type::ordinal_type;
  using execution_space = typename wrapped_type::execution_space;
  using device_type	= typename wrapped_type::device_type;
  using memory_traits	= typename wrapped_type::memory_traits;
  using size_type	= typename wrapped_type::size_type;
};
#endif




//*******************************
// Kokkos dense matrix
//*******************************
#ifdef HAVE_TRILINOS
template <typename wrapped_type>
struct traits<
  Matrix<
    wrapped_type,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_dense_matrix_kokkos<
	wrapped_type
	>::value
      >
    >
  >
  : public containers_shared_traits<
  Matrix<wrapped_type>,
  wrapped_type,
  false, true, false,
  WrappedPackageIdentifier::Kokkos,
  true, //true because kokkos is for shared mem
  // static view if the number of runtime determined dimensions == 0
  wrapped_type::traits::rank_dynamic==0
  >
{

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::DenseKokkos;

  using scalar_t	  = typename wrapped_type::traits::value_type;
  using layout		  = typename wrapped_type::traits::array_layout;
  using ordinal_t	  = typename wrapped_type::traits::size_type;
  using execution_space   = typename wrapped_type::traits::execution_space;
  using memory_space	  = typename wrapped_type::traits::memory_space;
  using device_type	  = typename wrapped_type::traits::device_type;
  using memory_traits	  = typename wrapped_type::traits::memory_traits;
  using host_mirror_space = typename wrapped_type::traits::host_mirror_space;
  using host_mirror_t     = typename wrapped_type::host_mirror_type;
};
#endif


}}}//end namespace pressio::containers::details
#endif
