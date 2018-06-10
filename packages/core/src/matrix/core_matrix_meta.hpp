
#ifndef CORE_MATRIX_META_HPP_
#define CORE_MATRIX_META_HPP_

#include "meta/core_meta_basic.hpp"
//#include "meta/core_meta_detect_operators.hpp"
//#include "meta/core_meta_detect_typedefs.hpp"
#include "vector/core_vector_meta.hpp"


namespace core{
namespace meta {


  template<typename T>
  struct is_matrixEigen : std::integral_constant<bool,!is_vectorEigen<T>::value> {};


  template <typename T, typename enable = void>
  struct is_denseMatrixEigen : std::false_type {};

  
  /* a type is a dense eigen matrix if
     (a) it is not a eigen vector, 
     (b) matches a eigen matrix sized at compile time
     (c) its rows and cols are fully dynamic
  */
  template<typename T>
  struct is_denseMatrixEigen<T, typename
  			        std::enable_if<is_matrixEigen<T>::value &&
  					       (std::is_same<T,
  						             Eigen::Matrix<typename T::Scalar,
  						                           T::RowsAtCompileTime,
  						                           T::ColsAtCompileTime
  						                           >
  						            >::value ||
  						std::is_same<T,
  						             Eigen::Matrix<typename T::Scalar,
  						                           Eigen::Dynamic,
  						                           Eigen::Dynamic>
  						            >::value
						)
  					       >::type
  			     > : std::true_type{};


  template <typename T, typename enable = void>
  struct is_stdlibMatrix : std::false_type {};

  template <typename T>
  struct  is_stdlibMatrix<T,
  			  typename
  			  std::enable_if<std::is_same<T,std::vector<std::vector<typename
  										T::value_type::value_type
  										>
  								    >
  						      >::value
  					 >::type
  			 > : std::true_type{};
  
    
} // namespace meta
} // namespace core
#endif
