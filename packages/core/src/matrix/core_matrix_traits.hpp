
#ifndef CORE_MATRIX_MATRIX_TRAITS_HPP_
#define CORE_MATRIX_MATRIX_TRAITS_HPP_

#include "../core_forward_declarations.hpp"
#include "../meta/core_matrix_meta.hpp"
#include "../meta/core_meta_detect_typedefs.hpp"
#include "../meta/core_meta_detect_operators.hpp"

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
{

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;
  using wrapped_t = wrapped_type;
  using derived_t = Matrix<wrapped_t>;

  static constexpr int isMatrix = 1;
  static constexpr int isEigen = 1;
  static constexpr int isEpetra = 0;
  static constexpr int isDense = 1;
  static constexpr int isSparse = !isDense;
  static constexpr int isSharedMem = 1;
  static constexpr int isVector = !isMatrix;
  static constexpr int isDistributed = !isSharedMem;
  static constexpr int isStdlib = 0;
  static constexpr int isStatic =
    ( wrapped_t::RowsAtCompileTime != Eigen::Dynamic &&
      wrapped_t::ColsAtCompileTime != Eigen::Dynamic );

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
{

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = typename wrapped_type::StorageIndex;
  //  ordinal has to be integral and signed
  static_assert( std::is_integral<ordinal_t>::value &&
  		 std::is_signed<ordinal_t>::value,
  "ordinal type for indexing eigen sparse matrix has to be signed"
  		 );
  
  using wrapped_t = wrapped_type;
  using derived_t = Matrix<wrapped_t>;

  static constexpr int isMatrix = 1;
  static constexpr int isRowMajor = wrapped_type::IsRowMajor;
  static constexpr int isColMajor = !isRowMajor;
  static constexpr int isEigen = 1;
  static constexpr int isEpetra = 0;
  static constexpr int isDense = 0;
  static constexpr int isSparse = 1;
  static constexpr int isSharedMem = 1;
  static constexpr int isVector = !isMatrix;
  static constexpr int isDistributed = !isSharedMem;
  static constexpr int isStdlib = 0;
  static constexpr int isStatic = 0;
};

  

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
{

  using scalar_t = typename wrapped_type::value_type::value_type;
  using ordinal_t = int;
  using wrapped_t = wrapped_type;
  using derived_t = Matrix<wrapped_t>;

  static constexpr int isMatrix = 1;
  static constexpr int isStdlib = 1;
  static constexpr int isEigen = 0;
  static constexpr int isEpetra = 0;
  static constexpr int isDense = 1;
  static constexpr int isSparse = !isDense;
  static constexpr int isSharedMem = 1;
  static constexpr int isVector = !isMatrix;
  static constexpr int isDistributed = !isSharedMem;
  static constexpr int isStatic = 0;
};


//***********************************
// epetra sparse distributed matrix 
//***********************************
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
{

  using wrapped_t = wrapped_type;
  using derived_t = Matrix<wrapped_t>;
  using scalar_t = defaultTypes::epetra_scalar_t;
  using local_ordinal_t = core::defaultTypes::epetra_lo_t;
  using global_ordinal_t = core::defaultTypes::epetra_go_t1;
  using communicator_t = Epetra_Comm;
  using row_map_t = Epetra_Map;
  using col_map_t = Epetra_Map;
  using range_map_t = Epetra_Map;
  using domain_map_t = Epetra_Map;
  using crs_graph_t = Epetra_CrsGraph;

  static constexpr int isMatrix = 1;
  static constexpr int isEpetra = 1;
  static constexpr int isSparse = 1;
  static constexpr int isDistributed = 1;
  static constexpr int isEigen = 0;
  static constexpr int isDense = 0;
  static constexpr int isSharedMem = 0;
  static constexpr int isVector = 0;
  static constexpr int isStdlib = 0;
  static constexpr int isStatic = 0;
};


//***********************************
// epetra dense distributed matrix 
//***********************************
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
{

  using wrapped_t = wrapped_type;
  using derived_t = Matrix<wrapped_t>;
  using scalar_t = defaultTypes::epetra_scalar_t;
  using local_ordinal_t = core::defaultTypes::epetra_lo_t;
  using global_ordinal_t = core::defaultTypes::epetra_go_t1;
  using communicator_t = Epetra_Comm;
  using row_map_t = Epetra_BlockMap;

  static constexpr int isMatrix = 1;
  static constexpr int isEpetra = 1;
  static constexpr int isSparse = 0;
  static constexpr int isDistributed = 1;

  // a multivector can also be seen as literally a multi-vector.
  // But here is intended as distributed dense matrix.
  // So we set actsAsMultiVector = false, actingAsDenseMatrix = true
  static constexpr int actingAsMultiVector = 0;
  static constexpr int actingAsDenseMatrix = 1;

  static constexpr int isEigen = 0;
  static constexpr int isDense = 1;
  static constexpr int isSharedMem = 0;
  static constexpr int isVector = 0;
  static constexpr int isStdlib = 0;
  static constexpr int isStatic = 0;
};
  
    
////////////////////////////
}//end namespace details
////////////////////////////

  
namespace meta{

template <typename T, typename enable = void>
struct is_core_matrix_wrapper : std::false_type {};

template <typename T>
struct is_core_matrix_wrapper< T,
		       typename
		       std::enable_if<
			 core::details::traits<T>::isMatrix==1
			 >::type
		       > : std::true_type{};

#define STATIC_ASSERT_IS_CORE_MATRIX_WRAPPER(TYPE) \
  static_assert( core::meta::is_core_matrix_wrapper<TYPE>::value, \
		 "THIS_IS_NOT_A_CORE_MATRIX_WRAPPER")
  
}//end meta

  
}//end namespace core
#endif
