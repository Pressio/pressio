
#include <gtest/gtest.h>
#include "MPL_ALL"

TEST(containers_meta_detect_operators, operatorAddDetecting)
{
  using namespace pressio;

  class A{
   public:
    A operator+(const A & other) { return A(); }
  };
  EXPECT_EQ( mpl::has_addition_op<A>::value, true);

  class B{};
  EXPECT_EQ( mpl::has_addition_op<B>::value, false);

  class C{
  public:
    C operator+(const C & other) const { return C(); }
  };
  EXPECT_EQ( mpl::has_addition_op<C>::value, true);

  class D{
  public:
    D operator+(D & other) { return D(); }
  };
  EXPECT_EQ( mpl::has_addition_op<D>::value, false);

  class E{
  public:
    E operator+(E other) { return E(); }
  };  
  EXPECT_EQ( mpl::has_addition_op<E>::value, true);
}


TEST(containers_meta_detect_operators, operatorSubtractDetecting)
{
  using namespace pressio;

  class A{
   public:
    A operator-(const A & other) { return A(); }
  };
  EXPECT_EQ( mpl::has_subtraction_op<A>::value, true);

  class B{};
  EXPECT_EQ( mpl::has_subtraction_op<B>::value, false);

  class C{
  public:
    C operator-(const C & other) const { return C(); }
  };
  EXPECT_EQ( mpl::has_subtraction_op<C>::value, true);

  class D{
  public:
    D operator-(D & other) { return D(); }
  };
  EXPECT_EQ( mpl::has_subtraction_op<D>::value, false);

  class E{
  public:
    E operator-(E other) { return E(); }
  };  
  EXPECT_EQ( mpl::has_subtraction_op<E>::value, true);
}


TEST(containers_meta_detect_operators, operatorStarDetecting)
{
  using namespace pressio;

  class A{
   public:
    A operator*(const A & other) { return A(); }
  };
  EXPECT_EQ( mpl::has_multiplication_op<A>::value, true);

  class B{};
  EXPECT_EQ( mpl::has_multiplication_op<B>::value, false);

  class C{
  public:
    C operator*(const C & other) const { return C(); }
  };
  EXPECT_EQ( mpl::has_multiplication_op<C>::value, true);

  class D{
  public:
    D operator*(D & other) { return D(); }
  };
  EXPECT_EQ( mpl::has_multiplication_op<D>::value, false);

  class E{
  public:
    E operator*(E other) { return E(); }
  };
  EXPECT_EQ( mpl::has_multiplication_op<E>::value, true);
}



TEST(containers_meta_detect_operators, operatorCompAssignDetecting)
{
  using namespace pressio;

  class A{
   public:
    A & operator+=(const A & other) { return *this; }
  };
  EXPECT_EQ( mpl::has_addition_assign_op<A>::value, true);

  class B{};
  EXPECT_EQ( mpl::has_addition_assign_op<B>::value, false);

  class D{
  public:
    D & operator+=(D & other) { return *this; }
  };
  EXPECT_EQ( mpl::has_addition_assign_op<D>::value, false);

  class E{
  public:
    E & operator+=(E other) { return *this; }
  };
  EXPECT_EQ( mpl::has_addition_assign_op<E>::value, true);
}



TEST(containers_meta_detect_operators, operatorCompAssignMinusDetecting)
{
  using namespace pressio;

  class A{
   public:
    A & operator-=(const A & other) { return *this; }
  };
  EXPECT_EQ( mpl::has_subtraction_assign_op<A>::value, true);

  class B{};
  EXPECT_EQ( mpl::has_subtraction_assign_op<B>::value, false);

  class D{
  public:
    D & operator-=(D & other) { return *this; }
  };
  EXPECT_EQ( mpl::has_subtraction_assign_op<D>::value, false);

  class E{
  public:
    E & operator-=(E other) { return *this; }
  };
  EXPECT_EQ( mpl::has_subtraction_assign_op<E>::value, true);
}
