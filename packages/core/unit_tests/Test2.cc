
#include <gtest/gtest.h>
#include "meta/core_meta_basic.hpp"
#include "meta/core_meta_detect_operators.hpp"
#include "meta/core_meta_detect_typedefs.hpp"

TEST(core_meta, isDefaultConstructible)
{
  class A{
  public:
    A() = default;
  };
  class B{
  public:
    B(){}
  };
  class C{
  public:
    C() = delete;
  };

  EXPECT_EQ( core::meta::is_default_constructible<A>::value, true);
  EXPECT_EQ( core::meta::is_default_constructible<B>::value, true);
  EXPECT_EQ( core::meta::is_default_constructible<C>::value, false);
}
