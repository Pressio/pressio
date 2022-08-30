
#include <gtest/gtest.h>
#include "pressio_optimizers.hpp"

TEST(optimizers_param_list, getDefault)
{
  using scalar_t = double;
  using opt_param_t = pressio::optimizers::Parameters<scalar_t>;
  opt_param_t MyPars;
  EXPECT_EQ(MyPars.maxIterations(), 100); //100 is default
  EXPECT_DOUBLE_EQ( MyPars.getGradientNormOptimalityTolerance(), 1e-10);
  EXPECT_DOUBLE_EQ( MyPars.getStepNormOptimalityTolerance(), 1e-12);
}

TEST(optimizers_param_list, setGet)
{
  using scalar_t = double;
  using opt_param_t = pressio::optimizers::Parameters<scalar_t>;
  opt_param_t MyPars;
  MyPars.setMaxIterations(33);
  MyPars.setGradientNormOptimalityTolerance(0.12345);
  MyPars.setStepNormOptimalityTolerance(1e-2);

  EXPECT_EQ(MyPars.maxIterations(), 33);
  EXPECT_DOUBLE_EQ( MyPars.getGradientNormOptimalityTolerance(), 0.12345);
  EXPECT_DOUBLE_EQ( MyPars.getStepNormOptimalityTolerance(), 0.01);
}
