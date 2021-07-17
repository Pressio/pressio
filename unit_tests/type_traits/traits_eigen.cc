#include <gtest/gtest.h>
#include "pressio_type_traits.hpp"

TEST(type_traits, eigen_dynamic_vector)
{
  using T = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using traits = pressio::traits<T>;

  static_assert(std::is_same<typename traits::scalar_type, double>::value, "");
}

TEST(type_traits, two_vector_scalar_compatible)
{
	using eigvec_d_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
	using eigvec_i_t = Eigen::Matrix<int, Eigen::Dynamic, 1>;

  using myv1_t = eigvec_d_t;
  using myv2_t = myv1_t;
  using myv3_t = eigvec_i_t;
  static_assert( pressio::are_scalar_compatible<myv1_t, myv2_t>::value, "");
  static_assert( !pressio::are_scalar_compatible<myv1_t, myv3_t>::value, "");
}

TEST(type_traits, three_vector_scalar_compatible)
{
using eigvec_d_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using eigvec_i_t = Eigen::Matrix<int, Eigen::Dynamic, 1>;
using eigvec_f_t = Eigen::Matrix<float, Eigen::Dynamic, 1>;

  using myv1_t = eigvec_d_t;
  using myv2_t = myv1_t;
  using myv3_t = myv2_t;
  static_assert( pressio::are_scalar_compatible<myv1_t, myv2_t, myv3_t>::value, "");

  using myv4_t = eigvec_i_t;
  static_assert( !pressio::are_scalar_compatible<myv1_t, myv2_t, myv4_t>::value, "");

  using myv5_t = eigvec_f_t;
  static_assert( !pressio::are_scalar_compatible<myv1_t, myv2_t, myv5_t>::value, "");
}

TEST(type_traits, four_vector_scalar_compatible)
{
	using eigvec_d_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
	using eigvec_i_t = Eigen::Matrix<int, Eigen::Dynamic, 1>;

  using myv1_t = eigvec_d_t;
  using myv2_t = myv1_t;
  using myv3_t = myv2_t;
  using myv4_t = eigvec_i_t;

  static_assert( pressio::are_scalar_compatible<myv1_t, myv2_t,
                myv3_t, myv2_t>::value, "");
  static_assert( !pressio::are_scalar_compatible<myv1_t, myv2_t,
                myv3_t, myv4_t>::value, "");
}




// TEST(containers_meta_basic, isTeuchosRCP)
// {
//   using namespace pressio;

//   class foo{
//     int a_ = 0;
//     public:
//       foo(int a) : a_(a) {};
//   };

//   using foo_t1 = foo;
//   using foo_t2 = foo *;
//   using foo_t3 = std::shared_ptr<foo>;
//   using foo_t4 = Teuchos::RCP<foo>;
//   using foo_t5 = Teuchos::RCP<const foo>;

//   EXPECT_EQ( containers::predicates::is_teuchos_rcp<foo_t1>::value, false);
//   EXPECT_EQ( containers::predicates::is_teuchos_rcp<foo_t2>::value, false);
//   EXPECT_EQ( containers::predicates::is_teuchos_rcp<foo_t3>::value, false);
//   EXPECT_EQ( containers::predicates::is_teuchos_rcp<foo_t4>::value, true);
//   EXPECT_EQ( containers::predicates::is_teuchos_rcp<foo_t5>::value, true);
// }
