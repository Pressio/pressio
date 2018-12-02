#include <gtest/gtest.h>
#include "CORE_VECTOR"
#include "CORE_MATRIX"


TEST(core_meta_basic, isDefaultConstructible){
  using namespace rompp;

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

TEST(core_meta_basic, isComplexNumber){
  using namespace rompp;

  
  using t1 = std::complex<float>;
  static_assert( core::meta::is_std_complex<t1>::value, "should be complex" );
  using t2 = std::complex<double>;
  static_assert( core::meta::is_std_complex<t2>::value, "should be complex" );
  using t3 = std::complex<long double>;
  static_assert( core::meta::is_std_complex<t3>::value, "should be complex" );
  using t4 = double;
  static_assert( !core::meta::is_std_complex<t4>::value, "should not be complex" );
}


TEST(core_meta_basic, inheritanceVector){
  
  using namespace rompp;

  using eigV_t = core::Vector<Eigen::Matrix<double,4,1>>;
  using base_t1 = core::VectorSharedMemBase<eigV_t>;

  using stdV_t = core::Vector<std::vector<double>>; 
  using base_t2 = core::VectorSharedMemBase<stdV_t>;

#ifdef HAVE_TRILINOS
  using epeV_t = core::Vector<Epetra_Vector>;
  using base_t3 = core::VectorDistributedBase<epeV_t>;
#endif

  static_assert(core::meta::publicly_inherits_from<eigV_t,base_t1>::value, "");
  static_assert(core::meta::publicly_inherits_from<eigV_t,base_t2>::value==false, "");
#ifdef HAVE_TRILINOS
  static_assert(core::meta::publicly_inherits_from<eigV_t,base_t3>::value==false, "");
#endif
  static_assert(core::meta::publicly_inherits_from<stdV_t,base_t1>::value==false, "");
  static_assert(core::meta::publicly_inherits_from<stdV_t,base_t2>::value==true, "");
#ifdef HAVE_TRILINOS
  static_assert(core::meta::publicly_inherits_from<stdV_t,base_t3>::value==false, "");
#endif
#ifdef HAVE_TRILINOS
  static_assert(core::meta::publicly_inherits_from<epeV_t,base_t1>::value==false, "");
  static_assert(core::meta::publicly_inherits_from<epeV_t,base_t2>::value==false, "");
  static_assert(core::meta::publicly_inherits_from<epeV_t,base_t3>::value==true, "");
#endif
}


#ifdef HAVE_TRILINOS
TEST(core_meta_basic, isTeuchosRCP){
  using namespace rompp;

  class foo{
    int a_ = 0;
    public:
      foo(int a) : a_(a) {};
  };

  using foo_t1 = foo;
  using foo_t2 = foo *;
  using foo_t3 = std::shared_ptr<foo>;
  using foo_t4 = Teuchos::RCP<foo>;
  using foo_t5 = Teuchos::RCP<const foo>;

  EXPECT_EQ( core::meta::is_teuchos_rcp_ptr<foo_t1>::value, false);
  EXPECT_EQ( core::meta::is_teuchos_rcp_ptr<foo_t2>::value, false);
  EXPECT_EQ( core::meta::is_teuchos_rcp_ptr<foo_t3>::value, false);
  EXPECT_EQ( core::meta::is_teuchos_rcp_ptr<foo_t4>::value, true); 
  EXPECT_EQ( core::meta::is_teuchos_rcp_ptr<foo_t5>::value, true); 
}
#endif


TEST(core_meta_basic, hasSizeMethod){
  
  using namespace rompp;

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
