
#include <gtest/gtest.h>
#include "CORE_BASIC"


TEST(core_meta_detect_operators, operatorAddDetecting)
{
  using namespace rompp;

  class A{
   public:
    A operator+(const A & other) { return A(); }
  };
  EXPECT_EQ( core::meta::has_add_op<A>::value, true);

  class B{};
  EXPECT_EQ( core::meta::has_add_op<B>::value, false);

  class C{
  public:
    C operator+(const C & other) const { return C(); }
  };
  EXPECT_EQ( core::meta::has_add_op<C>::value, true);

  class D{
  public:
    D operator+(D & other) { return D(); }
  };
  EXPECT_EQ( core::meta::has_add_op<D>::value, false);

  class E{
  public:
    E operator+(E other) { return E(); }
  };  
  EXPECT_EQ( core::meta::has_add_op<E>::value, true);
}


TEST(core_meta_detect_operators, operatorSubtractDetecting)
{
  using namespace rompp;

  class A{
   public:
    A operator-(const A & other) { return A(); }
  };
  EXPECT_EQ( core::meta::has_subtract_op<A>::value, true);

  class B{};
  EXPECT_EQ( core::meta::has_subtract_op<B>::value, false);

  class C{
  public:
    C operator-(const C & other) const { return C(); }
  };
  EXPECT_EQ( core::meta::has_subtract_op<C>::value, true);

  class D{
  public:
    D operator-(D & other) { return D(); }
  };
  EXPECT_EQ( core::meta::has_subtract_op<D>::value, false);

  class E{
  public:
    E operator-(E other) { return E(); }
  };  
  EXPECT_EQ( core::meta::has_subtract_op<E>::value, true);
}


TEST(core_meta_detect_operators, operatorStarDetecting)
{
  using namespace rompp;

  class A{
   public:
    A operator*(const A & other) { return A(); }
  };
  EXPECT_EQ( core::meta::has_star_op<A>::value, true);

  class B{};
  EXPECT_EQ( core::meta::has_star_op<B>::value, false);

  class C{
  public:
    C operator*(const C & other) const { return C(); }
  };
  EXPECT_EQ( core::meta::has_star_op<C>::value, true);

  class D{
  public:
    D operator*(D & other) { return D(); }
  };
  EXPECT_EQ( core::meta::has_star_op<D>::value, false);

  class E{
  public:
    E operator*(E other) { return E(); }
  };
  EXPECT_EQ( core::meta::has_star_op<E>::value, true);
}



TEST(core_meta_detect_operators, operatorCompAssignDetecting)
{
  using namespace rompp;

  class A{
   public:
    A & operator+=(const A & other) { return *this; }
  };
  EXPECT_EQ( core::meta::has_comp_assign_plus_op<A>::value, true);

  class B{};
  EXPECT_EQ( core::meta::has_comp_assign_plus_op<B>::value, false);

  class D{
  public:
    D & operator+=(D & other) { return *this; }
  };
  EXPECT_EQ( core::meta::has_comp_assign_plus_op<D>::value, false);

  class E{
  public:
    E & operator+=(E other) { return *this; }
  };
  EXPECT_EQ( core::meta::has_comp_assign_plus_op<E>::value, true);
}



TEST(core_meta_detect_operators, operatorCompAssignMinusDetecting)
{
  using namespace rompp;

  class A{
   public:
    A & operator-=(const A & other) { return *this; }
  };
  EXPECT_EQ( core::meta::has_comp_assign_minus_op<A>::value, true);

  class B{};
  EXPECT_EQ( core::meta::has_comp_assign_minus_op<B>::value, false);

  class D{
  public:
    D & operator-=(D & other) { return *this; }
  };
  EXPECT_EQ( core::meta::has_comp_assign_minus_op<D>::value, false);

  class E{
  public:
    E & operator-=(E other) { return *this; }
  };
  EXPECT_EQ( core::meta::has_comp_assign_minus_op<E>::value, true);
}
