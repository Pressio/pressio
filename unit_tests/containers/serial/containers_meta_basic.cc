#include <gtest/gtest.h>
#include "CONTAINERS_VECTOR"
#include "CONTAINERS_MATRIX"

TEST(containers_meta_basic, inheritanceVector){
  
  using namespace pressio;

  using eigV_t = containers::Vector<Eigen::Matrix<double,4,1>>;
  using base_t1 = containers::VectorSharedMemBase<eigV_t>;

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  using epeV_t = containers::Vector<Epetra_Vector>;
  using base_t3 = containers::VectorDistributedBase<epeV_t>;
#endif

  static_assert( mpl::publicly_inherits_from<eigV_t,base_t1>::value, "");
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  static_assert( mpl::publicly_inherits_from<eigV_t,base_t3>::value==false, "");
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  static_assert( mpl::publicly_inherits_from<epeV_t,base_t1>::value==false, "");
  static_assert( mpl::publicly_inherits_from<epeV_t,base_t3>::value==true, "");
#endif
}


#ifdef PRESSIO_ENABLE_TPL_TRILINOS
TEST(containers_meta_basic, isTeuchosRCP){
  using namespace pressio;

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

  EXPECT_EQ( containers::meta::is_teuchos_rcp<foo_t1>::value, false);
  EXPECT_EQ( containers::meta::is_teuchos_rcp<foo_t2>::value, false);
  EXPECT_EQ( containers::meta::is_teuchos_rcp<foo_t3>::value, false);
  EXPECT_EQ( containers::meta::is_teuchos_rcp<foo_t4>::value, true); 
  EXPECT_EQ( containers::meta::is_teuchos_rcp<foo_t5>::value, true); 
}
#endif


TEST(containers_meta_basic, hasSizeMethod){
  
  using namespace pressio;

  struct foo{
    int size() const{
      return 3;
    };
  };
  static_assert(containers::meta::has_size_method<foo>::value==true,"");

  struct foo2{
    int size(){
      return 3;
    };
  };
  static_assert(containers::meta::has_size_method<foo2>::value==true,"");

  struct foo3{
    void size(){
    };
  };
  static_assert(containers::meta::has_size_method<foo3>::value==false,"");
  
}
