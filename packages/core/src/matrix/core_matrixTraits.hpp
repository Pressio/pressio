
#ifndef CORE_MATRIXTRAITS_HPP
#define CORE_MATRIXTRAITS_HPP

#include "forwardDeclarations.hpp"
#include "meta.hpp"
#include <vector>
#include <Eigen/Core>


namespace core{
namespace details{

  // //*******************************
  // // std::vector-based matrix wrapper 
  // //******************************* 
  // template <typename scalar_type,
  // 	    typename ordinal_type>
  // struct traits<matrix<std::vector<std::vector<scalar_type>>,
  // 		       scalar_type,
  // 		       ordinal_type> >
  // {
  //   using scalar_t = scalar_type;
  //   using ordinal_t = ordinal_type;
  //   using wrapped_t = std::vector<std::vector<scalar_t>>;
  //   using derived_t = vector<wrapped_t,scalar_t,ordinal_t>;
  //   enum {
  //     isMatrix=1,
  //     isVector = !isMatrix,
  //     isSerial = 1,
  //     isDistributed = !isSerial
  //   };
  // };

  //*******************************
  // eigen matrix wrapper 
  //******************************* 
  template <typename wrapped_type>
  struct traits<matrix<wrapped_type,
		       typename
		       std::enable_if< core::meta::is_matrixEigen<wrapped_type>::value >::type
		       >
		>     		       
  {
    using scalar_t = typename wrapped_type::Scalar;
    using ordinal_t = int;
    using wrapped_t = wrapped_type;
    using derived_t = matrix<wrapped_t>;
    enum {
      isMatrix=1,
      isVector = !isMatrix,
      isSerial = 1,
      isDistributed = !isSerial,
      isEigen = 1,
      isStatic =  ( wrapped_t::RowsAtCompileTime != Eigen::Dynamic && 
		    wrapped_t::ColsAtCompileTime != Eigen::Dynamic )
    };
  };
  
  
}//end namespace details
}//end namespace core

#endif
