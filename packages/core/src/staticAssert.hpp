
#ifndef CORE_STATIC_ASSERT_HPP
#define CORE_STATIC_ASSERT_HPP

#include <type_traits>



/* 
 * - in ROMPP_STATIC_ASSERT(CONDITION,MSG) the parameter CONDITION must be a compile time boolean
 *    expression, and MSG an enum listed in struct internal::static_assertion<true>
 *
 *  - currently ROMPP_STATIC_ASSERT can only be used in function scope
 */

namespace core{

// if native static_assert is enabled, let's use it
#define ROMPP_STATIC_ASSERT(X,MSG) static_assert(X,#MSG);


namespace details{

  template<bool condition>
  struct static_assertion {};

  template<>
  struct static_assertion<true>{
    enum {
      YOU_MIXED_VECTORS_OF_DIFFERENT_SIZES,
      OUT_OF_RANGE_ACCESS
    };
  };

} // end namespace details


/*
 // static assertion failing if the two vector expression types are not compatible (same fixed-size or dynamic size)
 #define ROMPP_STATIC_ASSERT_SAME_VECTOR_SIZE(TYPE0,TYPE1) \
   ROMPP_STATIC_ASSERT( int(TYPE0::SizeAtCompileTime)==int(TYPE1::SizeAtCompileTime),\
     YOU_MIXED_VECTORS_OF_DIFFERENT_SIZES)
*/

  
} // end namespace core


#endif // ROMPP_STATIC_ASSERT_H
