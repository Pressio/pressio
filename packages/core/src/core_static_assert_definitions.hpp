
#ifndef CORE_STATIC_ASSERT_DEFINITIONS_HPP_
#define CORE_STATIC_ASSERT_DEFINITIONS_HPP_

#include "vector/core_vector_meta.hpp"


namespace core{

// // static assertion failing if the type is not a vector
// #define STATIC_ASSERT_IS_VECTOR(TYPE) \
//   static_assert( TYPE::isVector, \
// 		 "NOT_A_VECTOR")
// #define STATIC_ASSERT_IS_NOT_VECTOR(TYPE) \
//   static_assert( TYPE::isVector, \
// 		 "IS_A_VECTOR")
// static assertion failing if the type is not a matrix
// #define STATIC_ASSERT_IS_MATRIX(TYPE) \
//   static_assert( TYPE::isMatrix, \
// 		 "NOT_A_MATRIX")

  
#define STATIC_ASSERT_IS_VECTOR_EIGEN(TYPE) \
  static_assert( core::meta::is_vectorEigen<TYPE>::value, \
		 "NOT_A_VECTOR_FROM_EIGEN")
#define STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(TYPE) \
  static_assert( !core::meta::is_vectorEigen<TYPE>::value, \
		 "IS_A_VECTOR_FROM_EIGEN")

#define STATIC_ASSERT_IS_STDLIB_VECTOR(TYPE) \
  static_assert( core::meta::is_stdlibVector<TYPE>::value, \
		 "NOT_A_STDLIB_VECTOR")
#define STATIC_ASSERT_IS_NOT_STDLIB_VECTOR(TYPE) \
  static_assert( !core::meta::is_stdlibVector<TYPE>::value, \
		 "IS_A_STDLIB_VECTOR")

  

  
} // end namespace core
#endif
