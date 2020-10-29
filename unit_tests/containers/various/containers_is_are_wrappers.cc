
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using eigvec_d_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;


TEST(containers_meta, is_wrapper)
{
  using myv1_t = pressio::containers::Vector<eigvec_d_t>;
  static_assert( pressio::containers::predicates::is_wrapper<myv1_t>::value, "");
}

TEST(containers_meta, are_wrappers)
{
  using myv1_t = pressio::containers::Vector<eigvec_d_t>;
  using myv2_t = myv1_t;
  using myv3_t = myv1_t;
  static_assert( pressio::containers::predicates::are_wrappers<myv1_t, myv2_t, myv3_t>::value, "");
}

TEST(containers_meta, are_wrappers2)
{
  using myv1_t = pressio::containers::Vector<eigvec_d_t>;
  using myv2_t = myv1_t;
  static_assert( !pressio::containers::predicates::are_wrappers<myv1_t, myv2_t, eigvec_d_t>::value, "");
}

TEST(containers_meta, are_wrappers3)
{
  using myv1_t = pressio::containers::Vector<eigvec_d_t>;
  static_assert( pressio::containers::predicates::are_wrappers<myv1_t>::value, "");
}
