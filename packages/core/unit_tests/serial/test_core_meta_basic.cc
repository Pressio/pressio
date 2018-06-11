
#include <gtest/gtest.h>
#include "meta/core_meta_basic.hpp"
#include "meta/core_meta_detect_operators.hpp"

TEST(core_meta_basic, isDefaultConstructible)
{
  class A{
  public:
    A() = default;
  };
  EXPECT_EQ( core::meta::is_default_constructible<A>::value, true);

  class B{
  public:
    B(){}
  };
  EXPECT_EQ( core::meta::is_default_constructible<B>::value, true);

  class C{
  public:
    C() = delete;
  };
  EXPECT_EQ( core::meta::is_default_constructible<C>::value, false);
}



TEST(core_meta_basic, isComplexNumber)
{
  using t1 = std::complex<float>;
  static_assert( core::meta::is_stdComplex<t1>::value, "should be complex" );
  using t2 = std::complex<double>;
  static_assert( core::meta::is_stdComplex<t2>::value, "should be complex" );
  using t3 = std::complex<long double>;
  static_assert( core::meta::is_stdComplex<t3>::value, "should be complex" );
  using t4 = double;
  static_assert( !core::meta::is_stdComplex<t4>::value, "should not be complex" );
}

