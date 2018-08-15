#include <gtest/gtest.h>

#include <type_traits>

#include "experimental/solvers_nonlinear_factory.hpp"
#include "experimental/solvers_nonlinear_traits.hpp"


TEST(solvers_factory_nonlinear, solversFactoryTestValidSolver)
{
  // Namespaces
  using namespace solvers;

  // Create linear solver using a valid solver type
  auto solver = NonlinearSolvers::createSolver<nonlinear::NewtonRaphson>();
   
  // Expectations
  auto value = std::is_same<decltype(solver), NonLinearSolverBase>::value;
  EXPECT_EQ( value, true);
}


TEST(solvers_factory_nonlinear, solversFactoryTestInvalidSolver)
{
  // Namespace
  using namespace solvers;

  // Create linear solver using an invalid solver type
  ASSERT_DEATH(NonlinearSolvers::createSolver<void>(), "Error: the nonlinear solver selected is not available or its name was mispelt");
}