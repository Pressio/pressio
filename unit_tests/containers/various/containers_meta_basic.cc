#include <gtest/gtest.h>
#include "pressio_containers.hpp"

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

  EXPECT_EQ( containers::predicates::is_teuchos_rcp<foo_t1>::value, false);
  EXPECT_EQ( containers::predicates::is_teuchos_rcp<foo_t2>::value, false);
  EXPECT_EQ( containers::predicates::is_teuchos_rcp<foo_t3>::value, false);
  EXPECT_EQ( containers::predicates::is_teuchos_rcp<foo_t4>::value, true);
  EXPECT_EQ( containers::predicates::is_teuchos_rcp<foo_t5>::value, true);
}
