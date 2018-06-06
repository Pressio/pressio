
#ifndef CORE_STATIC_ASSERT_HPP_
#define CORE_STATIC_ASSERT_HPP_

#include <type_traits>


namespace core{

// static assertion failing if the type is not a vector
#define STATIC_ASSERT_IS_VECTOR(TYPE) \
  std::static_assert( TYPE::isVector, \
		      THIS_IS_NOT_A_VECTOR)

// static assertion failing if the type is not a matrix
#define STATIC_ASSERT_IS_MATRIX(TYPE) \
  std::static_assert( TYPE::isMatrix, \
		      THIS_IS_NOT_A_MATRIX)

  
} // end namespace core
#endif
