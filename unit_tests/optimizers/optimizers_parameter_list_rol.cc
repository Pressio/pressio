
#include <gtest/gtest.h>
#include "pressio_optimizers.hpp"

TEST(optimizers_param_list, convertToRol)
{
  using scalar_t = double;
  using opt_param_t = pressio::optimizers::Parameters<scalar_t>;
  opt_param_t MyPars;
  MyPars.setGradientNormOptimalityTolerance(0.12345);
  MyPars.setStepNormOptimalityTolerance(1e-2);

  ROL::ParameterList rolParList;
  ::pressio::optimizers::convertToRolParameterList(MyPars, rolParList);

  {
    const auto entry = rolParList.sublist("Status Test").getEntry("Iteration Limit");
    const auto value = Teuchos::getValue<int>(entry);
    EXPECT_EQ(value, MyPars.maxIterations());
  }

  {
    const auto entry = rolParList.sublist("Status Test").getEntry("Gradient Tolerance");
    const auto value = Teuchos::getValue<scalar_t>(entry);
    EXPECT_EQ(value, 0.12345);
  }

  {
    const auto entry = rolParList.sublist("Status Test").getEntry("Step Tolerance");
    const auto value = Teuchos::getValue<scalar_t>(entry);
    EXPECT_EQ(value, 0.01);
  }
}
