
#include <gtest/gtest.h>
#include "CORE_BASIC"

TEST(core_meta_detect_typedefs, scalarTypedefDetect)
{
  class A{
  public:
    using scalar_type = double;
    scalar_type x;
  };
  static_assert( core::meta::is_detected<core::meta::has_scalar_typedef, A>::value, "");

  class B{
  public:
    using scalar_t = double;
    scalar_t x;
  };
  static_assert( core::meta::is_detected<core::meta::has_scalar_typedef, B>::value == false, "");

  struct C{
    using scalar_type = double;
    scalar_type x;
  };
  static_assert( core::meta::is_detected<core::meta::has_scalar_typedef, C>::value, "");
}


TEST(core_meta_detect_typedefs, ordinalTypedefDetect)
{
  class A{
  public:
    using ordinal_type = int;
    ordinal_type x;
  };
  EXPECT_EQ( core::meta::has_ordinal_typedef<A>::value, true);

  class B{
  public:
    using ordinal_t = int;
    ordinal_t x;
  };
  EXPECT_EQ( core::meta::has_ordinal_typedef<B>::value, false);

  struct C{
    using ordinal_type = int;
    ordinal_type x;
  };
  EXPECT_EQ( core::meta::has_ordinal_typedef<C>::value, true);  
}


TEST(core_meta_detect_typedefs, localglobalOrdinalTypedefDetect)
{
  class A{
  public:
    using local_ordinal_type = int;
    local_ordinal_type x;
    using global_ordinal_type = int;
    global_ordinal_type x2;
  };
  EXPECT_EQ( core::meta::has_local_ordinal_typedef<A>::value, true);
  EXPECT_EQ( core::meta::has_global_ordinal_typedef<A>::value, true);

  class B{
  public:
    using local_ordinal_t = int;
    local_ordinal_t x;
    using global_ordinal_t = int;
    global_ordinal_t x2;
  };
  EXPECT_EQ( core::meta::has_local_ordinal_typedef<B>::value, false);
  EXPECT_EQ( core::meta::has_global_ordinal_typedef<B>::value, false);

  struct C{
    using local_ordinal_type = int;
    local_ordinal_type x;
    using global_ordinal_type = int;
    global_ordinal_type x2;
  };
  EXPECT_EQ( core::meta::has_local_ordinal_typedef<C>::value, true);  
  EXPECT_EQ( core::meta::has_global_ordinal_typedef<C>::value, true);
}


TEST(core_meta_detect_typedefs, mapCommTypedefDetect)
{
  class A{
  public:
    using data_map_type = int;
    data_map_type x;
    using communicator_type = int;
    communicator_type x2;
  };
  EXPECT_EQ( core::meta::has_data_map_typedef<A>::value, true);
  EXPECT_EQ( core::meta::has_mpi_comm_typedef<A>::value, true);

  class B{
  public:
    using data_map_t = int;
    data_map_t x;
    using communicator_t = int;
    communicator_t x2;
  };
  EXPECT_EQ( core::meta::has_data_map_typedef<B>::value, false);
  EXPECT_EQ( core::meta::has_mpi_comm_typedef<B>::value, false);

  struct C{
    using data_map_type = int;
    data_map_type x;
    using communicator_type = int;
    communicator_type x2;
  };
  EXPECT_EQ( core::meta::has_data_map_typedef<C>::value, true);  
  EXPECT_EQ( core::meta::has_mpi_comm_typedef<C>::value, true);
}

