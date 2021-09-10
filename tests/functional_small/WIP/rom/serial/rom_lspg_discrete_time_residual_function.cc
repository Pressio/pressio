
#include <gtest/gtest.h>
#include "pressio_rom.hpp"


template<typename T>
struct VelocitiesManagerHelp
{
  T ynp1_;
  T yn_;
  T ynm1_;

  VelocitiesManagerHelp() = delete;
  VelocitiesManagerHelp(int N, double s0=1.)
    : ynp1_(N), yn_(N), ynm1_(N)
  {
    for (auto i=0; i<N; ++i)
    {
      ynp1_(i) = (double) i*s0;
      yn_(i) = 2.*ynp1_(i);
      ynm1_(i) = 3.*ynp1_(i);
    }
  }

  const T & operator()(::pressio::ode::nPlusOne ) const{ return ynp1_; }
  const T & operator()(::pressio::ode::n ) const{ return yn_; }
  const T & operator()(::pressio::ode::nMinusOne ) const{ return ynm1_; }
};

template<typename T>
struct StatesManagerHelp
{
  T ynp1_;
  T yn_;
  T ynm1_;

  StatesManagerHelp() = delete;
  StatesManagerHelp(int N, double s0=1.)
    : ynp1_(N), yn_(N), ynm1_(N)
  {
    for (auto i=0; i<N; ++i)
    {
      ynp1_(i) = (double) i*s0;
      yn_(i) = 2.*ynp1_(i);
      ynm1_(i) = 3.*ynp1_(i);
    }
  }

  const T & fomStateAt(::pressio::ode::nPlusOne ) const{ return ynp1_; }
  const T & fomStateAt(::pressio::ode::n ) const{ return yn_; }
  const T & fomStateAt(::pressio::ode::nMinusOne ) const{ return ynm1_; }
  const T & operator()(::pressio::ode::nPlusOne ) const{ return ynp1_; }
  const T & operator()(::pressio::ode::n ) const{ return yn_; }
  const T & operator()(::pressio::ode::nMinusOne ) const{ return ynm1_; }

};

TEST(rom_lspg_time_discrete_residual_eigen, bdf1)
{
  using eig_vec = Eigen::Matrix<double, -1, 1>;
  using vec_t  = pressio::containers::Vector<eig_vec>;

  int N = 8;
  StatesManagerHelp<vec_t> statesMgr(N);

  vec_t R(8);
  R.data()->setConstant(2.);

  using tag = ::pressio::ode::BDF1;
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

  using tag = ::pressio::ode::BDF1;
  pressio::rom::lspg::impl::unsteady::time_discrete_residual<tag>(statesMgr, R, 1.1, indices);

  EXPECT_DOUBLE_EQ(R(0), 1.-2.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(1), 4.-8.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(2), 5.-10.-1.1*2.);
  EXPECT_DOUBLE_EQ(R(3), 7.-14.-1.1*2.);
}

TEST(rom_lspg_time_discrete_residual_eigen, bdf2_with_sample_indices)
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

  using tag = ::pressio::ode::ode::BDF2;
  pressio::rom::lspg::impl::unsteady::time_discrete_residual<tag>(statesMgr, R, 1.1, indices);

  constexpr double cnp1 = 1;
  constexpr double cn   = -4./3.;
  constexpr double cnm1 = 1./3.;
  constexpr double cf = 2./3.;

  EXPECT_DOUBLE_EQ(R(0), cnp1*1.+cn*2. +cnm1*3. -1.1*cf*2.);
  EXPECT_DOUBLE_EQ(R(1), cnp1*4.+cn*8. +cnm1*12.-1.1*cf*2.);
  EXPECT_DOUBLE_EQ(R(2), cnp1*5.+cn*10.+cnm1*15.-1.1*cf*2.);
  EXPECT_DOUBLE_EQ(R(3), cnp1*7.+cn*14.+cnm1*21.-1.1*cf*2.);
}

TEST(rom_lspg_time_discrete_residual_eigen, cranknicolson_with_sample_indices)
{
  using eig_vec = Eigen::Matrix<double, -1, 1>;
  using vec_t  = pressio::containers::Vector<eig_vec>;

  int N = 8;
  StatesManagerHelp<vec_t> statesMgr(N);
  VelocitiesManagerHelp<vec_t> velocitiesMgr(4,5.);

  vec_t R(4);
  R.data()->setConstant(0.);

  vec_t indices(4);
  indices(0) = 1;
  indices(1) = 4;
  indices(2) = 5;
  indices(3) = 7;

  using tag = ::pressio::ode::CrankNicolson;
  pressio::rom::lspg::impl::unsteady::time_discrete_residual<tag>(statesMgr, velocitiesMgr,R, 1.1, indices);

  constexpr double cnp1 = 1;
  constexpr double cn   = -1;
  constexpr double cf = 0.5;
  EXPECT_DOUBLE_EQ(R(0), cnp1*1.+cn*2. -1.1*cf*(5*0 + 5*2*0));
  EXPECT_DOUBLE_EQ(R(1), cnp1*4.+cn*8. -1.1*cf*(5*1 + 5*2*1));
  EXPECT_DOUBLE_EQ(R(2), cnp1*5.+cn*10.-1.1*cf*(5*2 + 5*2*2));
  EXPECT_DOUBLE_EQ(R(3), cnp1*7.+cn*14.-1.1*cf*(5*3 + 5*2*3));
}
