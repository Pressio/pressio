
#ifndef CORE_MATRIX_TRAITS_HPP_
#define CORE_MATRIX_TRAITS_HPP_

#include "core_forward_declarations.hpp"
#include "core_matrix_meta.hpp"
#include <vector>
#include <Eigen/Core>


namespace core{
namespace details{


//***********************************
// eigen dense matrix specialization
//***********************************
template <typename wrapped_type>
struct traits<matrix<wrapped_type,
		     typename
		     std::enable_if<core::meta::is_denseMatrixEigen<wrapped_type>::value
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
    isSerial = 1,
    isVector = !isMatrix,
    isDistributed = !isSerial,
    isStdlib = 0,
    isStatic =  ( wrapped_t::RowsAtCompileTime != Eigen::Dynamic && 
		  wrapped_t::ColsAtCompileTime != Eigen::Dynamic )
  };
};


//***********************************
// matrix specialization for matrix 
// based on std::vector<std::vector<>>
//***********************************
template <typename wrapped_type>
struct traits<matrix<wrapped_type,
		     typename
		     std::enable_if<core::meta::is_stdlibMatrix<wrapped_type>::value
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
    isStdlib = 1,
    isEigen = 0,
    isDense = 1,
    isSerial = 1,
    isVector = !isMatrix,
    isDistributed = !isSerial,
    isStatic = 0
  };
};

  
  
}//end namespace details
}//end namespace core

#endif
