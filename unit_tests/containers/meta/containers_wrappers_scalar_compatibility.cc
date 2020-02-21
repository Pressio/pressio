#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using eigvec_d_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using eigvec_i_t = Eigen::Matrix<int, Eigen::Dynamic, 1>;
using eigvec_f_t = Eigen::Matrix<float, Eigen::Dynamic, 1>;

using eigm_d_t = Eigen::Matrix<double, -1, -1>;
using eigm_i_t = Eigen::Matrix<int, -1, -1>;
using eigm_f_t = Eigen::Matrix<float, -1, -1>;


TEST(containers_meta, two_vector_scalar_compatible){
  using myv1_t = pressio::containers::Vector<eigvec_d_t>;
  using myv2_t = myv1_t;
  using myv3_t = pressio::containers::Vector<eigvec_i_t>;
  static_assert( pressio::containers::meta::are_scalar_compatible<myv1_t, myv2_t>::value, "");
  static_assert( !pressio::containers::meta::are_scalar_compatible<myv1_t, myv3_t>::value, "");
}

TEST(containers_meta, three_vector_scalar_compatible){
  using myv1_t = pressio::containers::Vector<eigvec_d_t>;
  using myv2_t = myv1_t;
  using myv3_t = myv2_t;
  static_assert( pressio::containers::meta::are_scalar_compatible<myv1_t, myv2_t, myv3_t>::value, "");

  using myv4_t = pressio::containers::Vector<eigvec_i_t>;
  static_assert( !pressio::containers::meta::are_scalar_compatible<myv1_t, myv2_t, myv4_t>::value, "");

  using myv5_t = pressio::containers::Vector<eigvec_f_t>;
  static_assert( !pressio::containers::meta::are_scalar_compatible<myv1_t, myv2_t, myv5_t>::value, "");
}

TEST(containers_meta, four_vector_scalar_compatible){
  using myv1_t = pressio::containers::Vector<eigvec_d_t>;
  using myv2_t = myv1_t;
  using myv3_t = myv2_t;
  using myv4_t = pressio::containers::Vector<eigvec_i_t>;

  static_assert( pressio::containers::meta::are_scalar_compatible<myv1_t, myv2_t,
                myv3_t, myv2_t>::value, "");
  static_assert( !pressio::containers::meta::are_scalar_compatible<myv1_t, myv2_t,
                myv3_t, myv4_t>::value, "");
}
