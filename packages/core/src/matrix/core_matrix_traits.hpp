
#ifndef CORE_MATRIX_TRAITS_HPP_
#define CORE_MATRIX_TRAITS_HPP_

#include "../core_forward_declarations.hpp"
#include "./meta/core_matrix_meta.hpp"
#include <vector>
#include <Eigen/Core>

namespace core{
namespace details{

//***********************************
// eigen dense matrix 
//***********************************
template <typename wrapped_type>
struct traits<matrix<wrapped_type,
		     typename
		     std::enable_if<
		       core::meta::is_matrixDenseSerialEigen<
			                  wrapped_type>::value
				    >::type
		     >
	     >
{
  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;
  using wrapped_t = wrapped_type;
  using derived_t = matrix<wrapped_t>;

  static constexpr int isMatrix = 1;
  static constexpr int isEigen = 1;
  static constexpr int isDense = 1;
  static constexpr int isSparse = !isDense;
  static constexpr int isSerial = 1;
  static constexpr int isVector = !isMatrix;
  static constexpr int isDistributed = !isSerial;
  static constexpr int isStdlib = 0;
  static constexpr int isStatic =
    ( wrapped_t::RowsAtCompileTime != Eigen::Dynamic &&
      wrapped_t::ColsAtCompileTime != Eigen::Dynamic );
};


//***********************************
// eigen sparse matrix 
//***********************************
template <typename wrapped_type>
struct traits<matrix<wrapped_type,
		     typename
		     std::enable_if<
		       core::meta::is_matrixSparseSerialEigen<
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
  using derived_t = matrix<wrapped_t>;

  static constexpr int isMatrix = 1;
  static constexpr int isRowMajor = wrapped_type::IsRowMajor;
  static constexpr int isColMajor = !isRowMajor;
  static constexpr int isEigen = 1;
  static constexpr int isDense = 0;
  static constexpr int isSparse = 1;
  static constexpr int isSerial = 1;
  static constexpr int isVector = !isMatrix;
  static constexpr int isDistributed = !isSerial;
  static constexpr int isStdlib = 0;
  static constexpr int isStatic = 0;
};

  

//***********************************
// based on std::vector<std::vector<>>
//***********************************
template <typename wrapped_type>
struct traits<matrix<wrapped_type,
		     typename
		     std::enable_if<
		       core::meta::is_matrixDenseSerialStdlib<
			                          wrapped_type
						 >::value
				    >::type
		     >
	     >
{
  using scalar_t = typename wrapped_type::value_type::value_type;
  using ordinal_t = int;
  using wrapped_t = wrapped_type;
  using derived_t = matrix<wrapped_t>;

  static constexpr int isMatrix = 1;
  static constexpr int isStdlib = 1;
  static constexpr int isEigen = 0;
  static constexpr int isDense = 1;
  static constexpr int isSparse = !isDense;
  static constexpr int isSerial = 1;
  static constexpr int isVector = !isMatrix;
  static constexpr int isDistributed = !isSerial;
  static constexpr int isStatic = 0;
};

  
}//end namespace details


namespace meta{

template <typename T, typename enable = void>
struct is_coreMatrixWrapper : std::false_type {};

template <typename T>
struct is_coreMatrixWrapper< T,
		       typename
		       std::enable_if<
			 core::details::traits<T>::isMatrix==1
			 >::type
		       > : std::true_type{};

#define STATIC_ASSERT_IS_CORE_MATRIX_WRAPPER(TYPE) \
  static_assert( core::meta::is_coreMatrixWrapper<TYPE>::value, \
		 "THIS_IS_NOT_A_CORE_MATRIX_WRAPPER")
  
}//end meta

  
}//end namespace core
#endif
