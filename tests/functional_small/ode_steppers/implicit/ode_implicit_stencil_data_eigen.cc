
#include <gtest/gtest.h>
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"

template<class T>
bool all_equal_to(const T a, double value)
{
  for (int i=0; i<a.size(); ++i){
    if (a(i)!=value)
      return false;
  }
  return true;
}

TEST(ode, implicit_stencil_velocities_constructor)
{
  using T = Eigen::VectorXd;
  T a(5);

  // data should NOT reference a
  pressio::ode::ImplicitStencilVelocitiesContainerStatic<T, 1> data1(a);
  const auto & v11 = data1(pressio::ode::nPlusOne());
  EXPECT_EQ(data1.size(), 1);
  EXPECT_TRUE(v11.data() != a.data());

  pressio::ode::ImplicitStencilVelocitiesContainerStatic<T, 2> data2(a);
  const auto & v21 = data2(pressio::ode::nPlusOne());
  const auto & v22 = data2(pressio::ode::n());
  EXPECT_EQ(data2.size(), 2);
  EXPECT_TRUE(v21.data() != a.data());
  EXPECT_TRUE(v22.data() != a.data());

  pressio::ode::ImplicitStencilVelocitiesContainerStatic<T, 3> data3(a);
  const auto & v31 = data3(pressio::ode::nPlusOne());
  const auto & v32 = data3(pressio::ode::n());
  const auto & v33 = data3(pressio::ode::nMinusOne());
  EXPECT_EQ(data3.size(), 3);
  EXPECT_TRUE(v31.data() != a.data());
  EXPECT_TRUE(v32.data() != a.data());
  EXPECT_TRUE(v33.data() != a.data());

  pressio::ode::ImplicitStencilVelocitiesContainerStatic<T, 4> data4(a);
  const auto & v41 = data4(pressio::ode::nPlusOne());
  const auto & v42 = data4(pressio::ode::n());
  const auto & v43 = data4(pressio::ode::nMinusOne());
  const auto & v44 = data4(pressio::ode::nMinusTwo());
  EXPECT_EQ(data4.size(), 4);
  EXPECT_TRUE(v41.data() != a.data());
  EXPECT_TRUE(v42.data() != a.data());
  EXPECT_TRUE(v43.data() != a.data());
  EXPECT_TRUE(v44.data() != a.data());
}


TEST(ode, implicit_stencil_velocities_data)
{
  using T = Eigen::VectorXd;
  T a(5);

  pressio::ode::ImplicitStencilVelocitiesContainerStatic<T, 1> data1(a);
  auto & v11 = data1(pressio::ode::nPlusOne());
  v11.setConstant(1.);
  const auto & v11r = data1(pressio::ode::nPlusOne());
  EXPECT_TRUE(all_equal_to(v11r, 1.));

  pressio::ode::ImplicitStencilVelocitiesContainerStatic<T, 2> data2(a);
  auto & v21 = data2(pressio::ode::nPlusOne());
  auto & v22 = data2(pressio::ode::n());
  v21.setConstant(2.);
  v22.setConstant(3.);
  const auto & v21r = data2(pressio::ode::nPlusOne());
  const auto & v22r = data2(pressio::ode::n());
  EXPECT_TRUE(all_equal_to(v21r, 2.));
  EXPECT_TRUE(all_equal_to(v22r, 3.));

  pressio::ode::ImplicitStencilVelocitiesContainerStatic<T, 3> data3(a);
  auto & v31 = data3(pressio::ode::nPlusOne());
  auto & v32 = data3(pressio::ode::n());
  auto & v33 = data3(pressio::ode::nMinusOne());
  v31.setConstant(2.);
  v32.setConstant(3.);
  v33.setConstant(4.);
  const auto & v31r = data3(pressio::ode::nPlusOne());
  const auto & v32r = data3(pressio::ode::n());
  const auto & v33r = data3(pressio::ode::nMinusOne());
  EXPECT_TRUE(all_equal_to(v31r, 2.));
  EXPECT_TRUE(all_equal_to(v32r, 3.));
  EXPECT_TRUE(all_equal_to(v33r, 4.));

  pressio::ode::ImplicitStencilVelocitiesContainerStatic<T, 4> data4(a);
  auto & v41 = data4(pressio::ode::nPlusOne());
  auto & v42 = data4(pressio::ode::n());
  auto & v43 = data4(pressio::ode::nMinusOne());
  auto & v44 = data4(pressio::ode::nMinusTwo());
  v41.setConstant(2.);
  v42.setConstant(3.);
  v43.setConstant(4.);
  v44.setConstant(5.);
  const auto & v41r = data4(pressio::ode::nPlusOne());
  const auto & v42r = data4(pressio::ode::n());
  const auto & v43r = data4(pressio::ode::nMinusOne());
  const auto & v44r = data4(pressio::ode::nMinusTwo());
  EXPECT_TRUE(all_equal_to(v41r, 2.));
  EXPECT_TRUE(all_equal_to(v42r, 3.));
  EXPECT_TRUE(all_equal_to(v43r, 4.));
  EXPECT_TRUE(all_equal_to(v44r, 5.));
}


TEST(ode, implicit_stencil_states_constructor)
{
  using T = Eigen::VectorXd;
  T a(5);

  // data should NOT reference a
  pressio::ode::ImplicitStencilStatesContainerStatic<T, 1> data1(a);
  const auto & v11 = data1(pressio::ode::n());
  EXPECT_EQ(data1.size(), 1);
  EXPECT_TRUE(v11.data() != a.data());

  pressio::ode::ImplicitStencilStatesContainerStatic<T, 2> data2(a);
  const auto & v21 = data2(pressio::ode::n());
  const auto & v22 = data2(pressio::ode::nMinusOne());
  EXPECT_EQ(data2.size(), 2);
  EXPECT_TRUE(v21.data() != a.data());
  EXPECT_TRUE(v22.data() != a.data());

  pressio::ode::ImplicitStencilStatesContainerStatic<T, 3> data3(a);
  const auto & v31 = data3(pressio::ode::n());
  const auto & v32 = data3(pressio::ode::nMinusOne());
  const auto & v33 = data3(pressio::ode::nMinusTwo());
  EXPECT_EQ(data3.size(), 3);
  EXPECT_TRUE(v31.data() != a.data());
  EXPECT_TRUE(v32.data() != a.data());
  EXPECT_TRUE(v33.data() != a.data());

  pressio::ode::ImplicitStencilStatesContainerStatic<T, 4> data4(a);
  const auto & v41 = data4(pressio::ode::n());
  const auto & v42 = data4(pressio::ode::nMinusOne());
  const auto & v43 = data4(pressio::ode::nMinusTwo());
  const auto & v44 = data4(pressio::ode::nMinusThree());
  EXPECT_EQ(data4.size(), 4);
  EXPECT_TRUE(v41.data() != a.data());
  EXPECT_TRUE(v42.data() != a.data());
  EXPECT_TRUE(v43.data() != a.data());
  EXPECT_TRUE(v44.data() != a.data());
}


TEST(ode, implicit_stencil_states_data)
{
  using T = Eigen::VectorXd;
  T a(5);

  pressio::ode::ImplicitStencilStatesContainerStatic<T, 1> data1(a);
  auto & v11 = data1(pressio::ode::n());
  v11.setConstant(1.);
  const auto & v11r = data1(pressio::ode::n());
  EXPECT_TRUE(all_equal_to(v11r, 1.));

  pressio::ode::ImplicitStencilStatesContainerStatic<T, 2> data2(a);
  auto & v21 = data2(pressio::ode::n());
  auto & v22 = data2(pressio::ode::nMinusOne());
  v21.setConstant(2.);
  v22.setConstant(3.);
  const auto & v21r = data2(pressio::ode::n());
  const auto & v22r = data2(pressio::ode::nMinusOne());
  EXPECT_TRUE(all_equal_to(v21r, 2.));
  EXPECT_TRUE(all_equal_to(v22r, 3.));

  pressio::ode::ImplicitStencilStatesContainerStatic<T, 3> data3(a);
  auto & v31 = data3(pressio::ode::n());
  auto & v32 = data3(pressio::ode::nMinusOne());
  auto & v33 = data3(pressio::ode::nMinusTwo());
  v31.setConstant(2.);
  v32.setConstant(3.);
  v33.setConstant(4.);
  const auto & v31r = data3(pressio::ode::n());
  const auto & v32r = data3(pressio::ode::nMinusOne());
  const auto & v33r = data3(pressio::ode::nMinusTwo());
  EXPECT_TRUE(all_equal_to(v31r, 2.));
  EXPECT_TRUE(all_equal_to(v32r, 3.));
  EXPECT_TRUE(all_equal_to(v33r, 4.));

  pressio::ode::ImplicitStencilStatesContainerStatic<T, 4> data4(a);
  auto & v41 = data4(pressio::ode::n());
  auto & v42 = data4(pressio::ode::nMinusOne());
  auto & v43 = data4(pressio::ode::nMinusTwo());
  auto & v44 = data4(pressio::ode::nMinusThree());
  v41.setConstant(2.);
  v42.setConstant(3.);
  v43.setConstant(4.);
  v44.setConstant(5.);
  const auto & v41r = data4(pressio::ode::n());
  const auto & v42r = data4(pressio::ode::nMinusOne());
  const auto & v43r = data4(pressio::ode::nMinusTwo());
  const auto & v44r = data4(pressio::ode::nMinusThree());
  EXPECT_TRUE(all_equal_to(v41r, 2.));
  EXPECT_TRUE(all_equal_to(v42r, 3.));
  EXPECT_TRUE(all_equal_to(v43r, 4.));
  EXPECT_TRUE(all_equal_to(v44r, 5.));
}