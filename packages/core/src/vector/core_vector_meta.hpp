
#ifndef CORE_VECTOR_META_HPP_
#define CORE_VECTOR_META_HPP_

#include "../meta/core_meta_basic.hpp"
#include "../meta/core_meta_detect_operators.hpp"
#include "../meta/core_meta_detect_typedefs.hpp"
#include <Eigen/Dense>
#include <vector>


namespace core{
namespace meta {

  // check if type is an Eigen vector
  // this is true when at least one dimension is 1
  
  template <typename T, typename enable = void>
  struct is_vectorEigen : std::false_type {};

  template <typename T>
  struct is_vectorEigen< T,
			 typename
			 std::enable_if<
			   std::is_same<T,
					Eigen::Matrix<typename T::Scalar,
						      1,
						      T::ColsAtCompileTime
						      >
					>::value
			   >::type
  			 > : std::true_type{};

  template <typename T>
  struct is_vectorEigen< T,
			 typename
			 std::enable_if<
			   std::is_same<T,
					Eigen::Matrix<typename T::Scalar,
						      T::RowsAtCompileTime,
						      1>
					>::value
			   >::type
  			 > : std::true_type{};

  template <typename T, typename enable = void>
  struct is_stdlibVector : std::false_type {};

  template <typename T>
  struct is_stdlibVector<T,
  			 typename
  			 std::enable_if<
			   std::is_same<T,
					std::vector<typename T::value_type>
					>::value &&
			   // we do not want to have vector<vector<...>>
			   // so we need to check that the T::value_type is a
			   // scalar type or integral type or complex
			   (std::is_floating_point<typename T::value_type>::value ||
			    std::is_integral<typename T::value_type>::value ||
			    is_stdComplex<typename T::value_type>::value
			    )
			   >::type
  			 > : std::true_type{};

   
} // namespace meta
} // namespace core
#endif
