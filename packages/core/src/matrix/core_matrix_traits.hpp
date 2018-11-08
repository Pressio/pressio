
#ifndef CORE_MATRIX_MATRIX_TRAITS_HPP_
#define CORE_MATRIX_MATRIX_TRAITS_HPP_

#include "../core_forward_declarations.hpp"
#include "../meta/core_native_matrix_meta.hpp"
#include "../meta/core_meta_detect_typedefs.hpp"
#include "../meta/core_meta_detect_operators.hpp"
#include "../core_shared_traits.hpp"

namespace rompp{
namespace core{
namespace details{

  
//***********************************
// eigen dense matrix 
//***********************************
template <typename wrapped_type>
struct traits< Matrix<
    wrapped_type,
    typename
    std::enable_if<
      core::meta::is_matrix_dense_sharedmem_eigen<
	wrapped_type>::value
      >::type
    >
  >
  : public containers_shared_traits<Matrix<wrapped_type>,
				    wrapped_type, false, true, false,
				    WrappedPackageIdentifier::Eigen,true>,
    public matrix_shared_traits<false>
{

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::DenseEigen;
  
  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;

  static constexpr bool is_static =
    ( wrapped_type::RowsAtCompileTime != Eigen::Dynamic &&
      wrapped_type::ColsAtCompileTime != Eigen::Dynamic );
};


//***********************************
// eigen sparse matrix 
//***********************************
template <typename wrapped_type>
struct traits< Matrix<
    wrapped_type,
    typename
    std::enable_if<
      core::meta::is_matrix_sparse_sharedmem_eigen<
	wrapped_type
	>::value
      >::type
    >
  >
  : public containers_shared_traits<Matrix<wrapped_type>,
				    wrapped_type, false, true, false,
				    WrappedPackageIdentifier::Eigen, true>,
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
  static constexpr bool is_static = false;
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
      core::meta::is_matrix_sparse_distributed_epetra<
	wrapped_type
	>::value
      >::type
    >
  >
  : public containers_shared_traits<Matrix<wrapped_type>,
				    wrapped_type, false, true, false,
				    WrappedPackageIdentifier::Trilinos,false>,
    public matrix_shared_traits<true>
{

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::CrsEpetra;

  using scalar_t = default_types::epetra_scalar_t;
  using local_ordinal_t = core::default_types::epetra_lo_t;
  using global_ordinal_t = core::default_types::epetra_go_t1;
  using communicator_t = Epetra_Comm;
  using row_map_t = Epetra_Map;
  using col_map_t = Epetra_Map;
  using range_map_t = Epetra_Map;
  using domain_map_t = Epetra_Map;
  using crs_graph_t = Epetra_CrsGraph;

  static constexpr int is_static = 0;
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
      core::meta::is_matrix_dense_distributed_epetra<
	wrapped_type
	>::value
      >::type
    >
  >
  : public containers_shared_traits<Matrix<wrapped_type>,
				    wrapped_type, false, true, false,
				    WrappedPackageIdentifier::Trilinos,false>,
    public matrix_shared_traits<false>
{
  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::DenseEpetra;

  using scalar_t = default_types::epetra_scalar_t;
  using local_ordinal_t = core::default_types::epetra_lo_t;
  using global_ordinal_t = core::default_types::epetra_go_t1;
  using communicator_t = Epetra_Comm;
  using row_map_t = Epetra_BlockMap;

  // // a multivector can also be seen as literally a multi-vector.
  // // But here is intended as distributed dense matrix.
  // // So we set actsAsMultiVector = false, actingAsDenseMatrix = true
  // static constexpr int actingAsMultiVector = 0;
  // static constexpr int actingAsDenseMatrix = 1;

  static constexpr bool is_static = false;
};
#endif
  

  
//***********************************
// based on std::vector<std::vector<>>
//***********************************
template <typename wrapped_type>
struct traits< Matrix<
    wrapped_type,
    typename
    std::enable_if<
      core::meta::is_matrix_dense_sharedmem_stdlib<
	wrapped_type
	>::value
      >::type
    >
  >
  : public containers_shared_traits<Matrix<wrapped_type>,
				    wrapped_type, false, true, false,
				    WrappedPackageIdentifier::CppStdLib, true>,
    public matrix_shared_traits<false>
{

  static constexpr WrappedMatrixIdentifier
  wrapped_matrix_identifier = WrappedMatrixIdentifier::CppStdLib;
 
  using scalar_t = typename wrapped_type::value_type::value_type;
  using ordinal_t = int;
  static constexpr bool is_static = false;
};

  
}//end namespace details
}//end namespace core

}//end namespace rompp
#endif
