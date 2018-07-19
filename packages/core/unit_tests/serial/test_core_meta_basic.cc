#include <gtest/gtest.h>
#include "CORE_VECTOR"
#include "CORE_MATRIX"


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
  static_assert( core::meta::is_std_complex<t1>::value, "should be complex" );
  using t2 = std::complex<double>;
  static_assert( core::meta::is_std_complex<t2>::value, "should be complex" );
  using t3 = std::complex<long double>;
  static_assert( core::meta::is_std_complex<t3>::value, "should be complex" );
  using t4 = double;
  static_assert( !core::meta::is_std_complex<t4>::value, "should not be complex" );
}


TEST(core_meta_basic, inheritanceVector)
{
  using eigV_t = core::Vector<Eigen::Matrix<double,4,1>>;
  using base_t1 = core::VectorSerialBase<eigV_t>;

  using stdV_t = core::Vector<std::vector<double>>; 
  using base_t2 = core::VectorSerialBase<stdV_t>;

  using epeV_t = core::Vector<Epetra_Vector>;
  using base_t3 = core::VectorDistributedBase<epeV_t>;
  
  static_assert(core::meta::publicly_inherits_from<eigV_t,base_t1>::value, "");
  static_assert(core::meta::publicly_inherits_from<eigV_t,base_t2>::value==false, "");
  static_assert(core::meta::publicly_inherits_from<eigV_t,base_t3>::value==false, "");

  static_assert(core::meta::publicly_inherits_from<stdV_t,base_t1>::value==false, "");
  static_assert(core::meta::publicly_inherits_from<stdV_t,base_t2>::value==true, "");
  static_assert(core::meta::publicly_inherits_from<stdV_t,base_t3>::value==false, "");

  static_assert(core::meta::publicly_inherits_from<epeV_t,base_t1>::value==false, "");
  static_assert(core::meta::publicly_inherits_from<epeV_t,base_t2>::value==false, "");
  static_assert(core::meta::publicly_inherits_from<epeV_t,base_t3>::value==true, "");
}



TEST(core_meta_basic, hasSizeMethod)
{
  struct foo{
    int size() const{
      return 3;
    };
  };
  static_assert(core::meta::has_size_method<foo>::value==true,"");

  struct foo2{
    int size(){
      return 3;
    };
  };
  static_assert(core::meta::has_size_method<foo2>::value==true,"");

  struct foo3{
    void size(){
    };
  };
  static_assert(core::meta::has_size_method<foo3>::value==false,"");
  
}
