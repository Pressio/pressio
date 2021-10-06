#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"

namespace pressio { namespace testing {

TEST(type_traits, method_detect)
{
  struct A {
    int size()          const { return 0; };
    int size(int dim)   const { return 0; };
    int extent(int dim) const { return 0; };
  };

  static_assert(pressio::has_method_extent<A>::value);
  static_assert(pressio::has_method_size<A>::value);
  static_assert(pressio::has_method_size_with_arg<A>::value);
}

}} // pressio::testing