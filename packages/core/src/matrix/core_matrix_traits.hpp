
#ifndef CORE_MATRIX_TRAITS_HPP_
#define CORE_MATRIX_TRAITS_HPP_

#include "core_forward_declarations.hpp"
#include "core_matrix_meta.hpp"
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
  enum {
    isMatrix = 1,
    isEigen = 1,
    isDense = 1,
    isSparse = !isDense,
    isSerial = 1,
    isVector = !isMatrix,
    isDistributed = !isSerial,
    isStdlib = 0,
    isStatic =  ( wrapped_t::RowsAtCompileTime != Eigen::Dynamic && 
		  wrapped_t::ColsAtCompileTime != Eigen::Dynamic )
  };
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
  enum {
    isMatrix = 1,
    isEigen = 1,
    isDense = 0,
    isSparse = 1,
    isSerial = 1,
    isVector = !isMatrix,
    isDistributed = !isSerial,
    isStdlib = 0,
    isStatic = 0
  };
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
  enum {
    isMatrix = 1,
    isStdlib = 1,
    isEigen = 0,
    isDense = 1,
    isSparse = !isDense,
    isSerial = 1,
    isVector = !isMatrix,
    isDistributed = !isSerial,
    isStatic = 0
  };
};

  
  
}//end namespace details
}//end namespace core

#endif
