#include <gtest/gtest.h>

#include <type_traits>

#include "experimental/solvers_nonlinear_factory.hpp"
#include "experimental/solvers_nonlinear_traits.hpp"
#include "experimental/solvers_policy_nonlinear_iterative.hpp"


TEST(solvers_factory_nonlinear, solversFactoryTestValidSolver)
{
  // Namespaces
  using namespace pressio;
  using namespace solvers;

  // Create linear solver using a valid solver type
  auto solver = NonLinearSolvers::createIterativeSolver<nonlinear::NewtonRaphson, linear::Bicgstab>();

  // Expectations
  auto value = std::is_same<decltype(solver), NonLinearIterativeSolver<SolversNonLinearIterativeNewtonRaphsonPolicy, linear::Bicgstab>>::value;
  EXPECT_EQ( value, true);
}
