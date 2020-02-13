
#include <gtest/gtest.h>
#include "pressio_mpl.hpp"

TEST(containers_meta_basic, isDefaultConstructible){
  using namespace pressio;

  class A{
  public:
    A() = default;
  };
  EXPECT_EQ( mpl::is_default_constructible<A>::value, true);

  class B{
  public:
    B(){}
  };
  EXPECT_EQ( mpl::is_default_constructible<B>::value, true);

  class C{
  public:
    C() = delete;
  };
  EXPECT_EQ( mpl::is_default_constructible<C>::value, false);
}

TEST(containers_meta_basic, isComplexNumber){
  using namespace pressio;

  using t1 = std::complex<float>;
  static_assert( mpl::is_std_complex<t1>::value, "should be complex" );
  using t2 = std::complex<double>;
  static_assert( mpl::is_std_complex<t2>::value, "should be complex" );
  using t3 = std::complex<long double>;
  static_assert( mpl::is_std_complex<t3>::value, "should be complex" );
  using t4 = double;
  static_assert( !mpl::is_std_complex<t4>::value, "should not be complex" );
}

