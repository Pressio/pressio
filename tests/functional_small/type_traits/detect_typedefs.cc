
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"

TEST(type_traits, scalarTypedefDetect)
{
  using namespace pressio;

  class A{
  public:
    using scalar_type = double;
    scalar_type x;
  };
  static_assert( pressio::has_scalar_typedef<A>::value, "");

  class B{
  public:
    using scalar_t = double;
    scalar_t x;
  };
  static_assert( pressio::has_scalar_typedef<B>::value == false, "");

  struct C{
    using scalar_type = double;
    scalar_type x;
  };
  static_assert( pressio::has_scalar_typedef<C>::value, "");
}


TEST(type_traits, ordinalTypedefDetect)
{
  using namespace pressio;
  class A{
  public:
    using ordinal_type = int;
    ordinal_type x;
  };
  EXPECT_EQ( pressio::has_ordinal_typedef<A>::value, true);

  class B{
  public:
    using ordinal_t = int;
    ordinal_t x;
  };
  EXPECT_EQ( pressio::has_ordinal_typedef<B>::value, false);

  struct C{
    using ordinal_type = int;
    ordinal_type x;
  };
  EXPECT_EQ( pressio::has_ordinal_typedef<C>::value, true);
}


TEST(type_traits, localglobalOrdinalTypedefDetect)
{
  using namespace pressio;

  class A{
  public:
    using local_ordinal_type = int;
    local_ordinal_type x;
    using global_ordinal_type = int;
    global_ordinal_type x2;
  };
  EXPECT_EQ( pressio::has_local_ordinal_typedef<A>::value, true);
  EXPECT_EQ( pressio::has_global_ordinal_typedef<A>::value, true);

  class B{
  public:
    using local_ordinal_t = int;
    local_ordinal_t x;
    using global_ordinal_t = int;
    global_ordinal_t x2;
  };
  EXPECT_EQ( pressio::has_local_ordinal_typedef<B>::value, false);
  EXPECT_EQ( pressio::has_global_ordinal_typedef<B>::value, false);

  struct C{
    using local_ordinal_type = int;
    local_ordinal_type x;
    using global_ordinal_type = int;
    global_ordinal_type x2;
  };
  EXPECT_EQ( pressio::has_local_ordinal_typedef<C>::value, true);
  EXPECT_EQ( pressio::has_global_ordinal_typedef<C>::value, true);
}


TEST(type_traits, mapCommTypedefDetect)
{
  using namespace pressio;

  class A{
  public:
    using data_map_type = int;
    data_map_type x;
    using communicator_type = int;
    communicator_type x2;
  };
  EXPECT_EQ( pressio::has_data_map_typedef<A>::value, true);
  EXPECT_EQ( pressio::has_communicator_typedef<A>::value, true);

  class B{
  public:
    using data_map_t = int;
    data_map_t x;
    using communicator_t = int;
    communicator_t x2;
  };
  EXPECT_EQ( pressio::has_data_map_typedef<B>::value, false);
  EXPECT_EQ( pressio::has_communicator_typedef<B>::value, false);

  struct C{
    using data_map_type = int;
    data_map_type x;
    using communicator_type = int;
    communicator_type x2;
  };
  EXPECT_EQ( pressio::has_data_map_typedef<C>::value, true);
  EXPECT_EQ( pressio::has_communicator_typedef<C>::value, true);
}
