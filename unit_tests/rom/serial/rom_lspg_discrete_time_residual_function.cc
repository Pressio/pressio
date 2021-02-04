
#include <gtest/gtest.h>
#include "pressio_rom.hpp"

template<typename T>
struct StatesManagerHelp
{
  T ynp1_;
  T yn_;
  T ynm1_;

  StatesManagerHelp() = delete;
  StatesManagerHelp(int N)
    : ynp1_(N), yn_(N), ynm1_(N)
  {
    for (auto i=0; i<N; ++i)
    {
      ynp1_(i) = (double) i;
      yn_(i) = 2.*ynp1_(i);
      ynm1_(i) = 3.*ynp1_(i);
    }
  }

  const T & fomStateAt(::pressio::ode::nPlusOne ) const{ return ynp1_; }
  const T & fomStateAt(::pressio::ode::n ) const{ return yn_; }
  const T & fomStateAt(::pressio::ode::nMinusOne ) const{ return ynm1_; }
};

TEST(rom_lspg_time_discrete_residual_eigen, bdf1)
{
  using eig_vec = Eigen::Matrix<double, -1, 1>;
  using vec_t  = pressio::containers::Vector<eig_vec>;

  int N = 8;
  StatesManagerHelp<vec_t> statesMgr(N);

  vec_t R(8);
  R.data()->setConstant(2.);

  using tag = ::pressio::ode::implicitmethods::Euler;
  pressio::rom::lspg::impl::unsteady::time_discrete_residual<tag>(statesMgr, R, 1.1);

  EXPECT_DOUBLE_EQ(R(0), 0.-0.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(1), 1.-2.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(2), 2.-4.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(3), 3.-6.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(4), 4.-8.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(5), 5.-10.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(6), 6.-12.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(7), 7.-14.-1.1*2.);
}

TEST(rom_lspg_time_discrete_residual_eigen, bdf1_with_sample_indices)
{
  using eig_vec = Eigen::Matrix<double, -1, 1>;
  using vec_t  = pressio::containers::Vector<eig_vec>;

  int N = 8;
  StatesManagerHelp<vec_t> statesMgr(N);

  vec_t R(4);
  R.data()->setConstant(2.);

  vec_t indices(4);
  indices(0) = 1;
  indices(1) = 4;
  indices(2) = 5;
  indices(3) = 7;

  using tag = ::pressio::ode::implicitmethods::Euler;
  pressio::rom::lspg::impl::unsteady::time_discrete_residual<tag>(statesMgr, R, 1.1, indices);

  EXPECT_DOUBLE_EQ(R(0), 1.-2.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(1), 4.-8.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(2), 5.-10.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(3), 7.-14.-1.1*2.);
}
