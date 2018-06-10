
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


TEST(core_meta, operatorPlusDetecting)
{
  class A{
   public:
    A operator+(const A & other) { return A(); }
  };

  class B{};

  class C{
  public:
    C operator+(const C & other) const { return C(); }
  };

  class D{
  public:
    D operator+(D & other) { return D(); }
  };

  class E{
  public:
    E operator+(E other) { return E(); }
  };
  
  
  EXPECT_EQ( core::meta::has_add_op<A>::value, true);
  EXPECT_EQ( core::meta::has_add_op<B>::value, false);
  EXPECT_EQ( core::meta::has_add_op<C>::value, true);
  EXPECT_EQ( core::meta::has_add_op<D>::value, false);
  EXPECT_EQ( core::meta::has_add_op<E>::value, true);
}


TEST(core_meta, operatorMinusDetecting)
{
  class A{
   public:
    A operator-(const A & other) { return A(); }
  };

  class B{};

  class C{
  public:
    C operator-(const C & other) const { return C(); }
  };

  class D{
  public:
    D operator-(D & other) { return D(); }
  };

  class E{
  public:
    E operator-(E other) { return E(); }
  };
  
  
  EXPECT_EQ( core::meta::has_diff_op<A>::value, true);
  EXPECT_EQ( core::meta::has_diff_op<B>::value, false);
  EXPECT_EQ( core::meta::has_diff_op<C>::value, true);
  EXPECT_EQ( core::meta::has_diff_op<D>::value, false);
  EXPECT_EQ( core::meta::has_diff_op<E>::value, true);
}


TEST(core_meta, operatorStarDetecting)
{
  class A{
   public:
    A operator*(const A & other) { return A(); }
  };

  class B{};

  class C{
  public:
    C operator*(const C & other) const { return C(); }
  };

  class D{
  public:
    D operator*(D & other) { return D(); }
  };

  class E{
  public:
    E operator*(E other) { return E(); }
  };
  
  
  EXPECT_EQ( core::meta::has_star_op<A>::value, true);
  EXPECT_EQ( core::meta::has_star_op<B>::value, false);
  EXPECT_EQ( core::meta::has_star_op<C>::value, true);
  EXPECT_EQ( core::meta::has_star_op<D>::value, false);
  EXPECT_EQ( core::meta::has_star_op<E>::value, true);
}



TEST(core_meta, operatorCompAssignDetecting)
{
  class A{
   public:
    A & operator+=(const A & other) { return *this; }
  };

  class B{};

  class D{
  public:
    D & operator+=(D & other) { return *this; }
  };

  class E{
  public:
    E & operator+=(E other) { return *this; }
  };
  
  
  EXPECT_EQ( core::meta::has_comp_assign_plus_op<A>::value, true);
  EXPECT_EQ( core::meta::has_comp_assign_plus_op<B>::value, false);
  EXPECT_EQ( core::meta::has_comp_assign_plus_op<D>::value, false);
  EXPECT_EQ( core::meta::has_comp_assign_plus_op<E>::value, true);
}



TEST(core_meta, operatorCompAssignMinusDetecting)
{
  class A{
   public:
    A & operator-=(const A & other) { return *this; }
  };

  class B{};

  class D{
  public:
    D & operator-=(D & other) { return *this; }
  };

  class E{
  public:
    E & operator-=(E other) { return *this; }
  };
  
  
  EXPECT_EQ( core::meta::has_comp_assign_minus_op<A>::value, true);
  EXPECT_EQ( core::meta::has_comp_assign_minus_op<B>::value, false);
  EXPECT_EQ( core::meta::has_comp_assign_minus_op<D>::value, false);
  EXPECT_EQ( core::meta::has_comp_assign_minus_op<E>::value, true);
}



// TEST(core_meta, operatorsPlusDetecting)
// {
//   class A{
//    public:
//     A() = default;

//     A operator+(const A & other) {
//       A res(other.size());
//       *res.data() = this->data_ + *other.data();
//       return res;
//     }

//     A operator-(const A & other) {
//       A res(other.size());
//       *res.data() = this->data_ - *other.data();
//       return res;
//     }

//     A operator*(const A & other) {
//       A res(other.size());
//       *res.data() = this->data_ * (*other.data());
//       return res;
//     }

//     A & operator+=(const A & other) {
//       this->data_ += *other.data();
//       return *this;
//     }

//     A & operator-=(const A & other) {
//       this->data_ -= *other.data();
//       return *this;
//     }
//   };

//   // EXPECT_EQ( core::meta::is_default_constructible<A>::value, true);
//   // EXPECT_EQ( core::meta::is_default_constructible<B>::value, true);
//   // EXPECT_EQ( core::meta::is_default_constructible<C>::value, false);
// }
