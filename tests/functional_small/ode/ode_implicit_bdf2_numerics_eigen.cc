
#include <gtest/gtest.h>
#include "pressio/solvers.hpp"
#include "pressio/ode_implicit.hpp"
#include "testing_apps.hpp"

TEST(ode, implicit_bdf2_policy_default_created)
{
  using namespace pressio;
  using problem_t = ode::testing::refAppForImpEigen;
  using state_t = typename problem_t::state_type;
  problem_t problemObj;
  state_t y(3);
  y = problemObj.getInitCond();

  auto stepperObj = ode::create_bdf2_stepper(problemObj, y);

  using jac_t = typename problem_t::jacobian_type;
  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver = nonlinearsolvers::create_newton_raphson(stepperObj,y,linSolverObj);

  // integrate in time
  ode::step_count_type nSteps = 4;
  double dt = 0.01;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, nSteps, NonLinSolver);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  problemObj.analyticAdvanceBackEulerNSteps(dt, 1);
  problemObj.analyticAdvanceBDF2NSteps(dt, 3);
  std::cout << std::setprecision(14) << problemObj.y << "\n";
  EXPECT_DOUBLE_EQ(y(0), problemObj.y(0));
  EXPECT_DOUBLE_EQ(y(1), problemObj.y(1));
  EXPECT_DOUBLE_EQ(y(2), problemObj.y(2));
}

TEST(ode, implicit_bdf2_policy_default_created_specify_types)
{
  using namespace pressio;
  using problem_t = ode::testing::refAppForImpEigen;
  using state_t = typename problem_t::state_type;
  problem_t problemObj;
  state_t y = problemObj.getInitCond();

  using res_t = typename problem_t::velocity_type;
  using jac_t = typename problem_t::jacobian_type;
  auto stepperObj = ode::create_bdf2_stepper_partial_deduction<res_t, jac_t>(problemObj, y);

  using jac_t = typename problem_t::jacobian_type;
  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver = nonlinearsolvers::create_newton_raphson(stepperObj,y,linSolverObj);

  ode::step_count_type nSteps = 4;
  double dt = 0.01;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, nSteps, NonLinSolver);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  problemObj.analyticAdvanceBackEulerNSteps(dt, 1);
  problemObj.analyticAdvanceBDF2NSteps(dt, 3);
  std::cout << std::setprecision(14) << problemObj.y << "\n";
  EXPECT_DOUBLE_EQ(y(0), problemObj.y(0));
  EXPECT_DOUBLE_EQ(y(1), problemObj.y(1));
  EXPECT_DOUBLE_EQ(y(2), problemObj.y(2));
}

TEST(ode, implicit_bdf2_custom_policy)
{
  using namespace pressio;
  using problem_t = ode::testing::refAppForImpEigen;
  using state_t = typename problem_t::state_type;
  problem_t problemObj;
  state_t y = problemObj.getInitCond();

  using res_t = typename problem_t::velocity_type;
  using jac_t = typename problem_t::jacobian_type;
  using res_pol_t = ode::impl::ResidualStandardPolicyBdf<state_t, res_t>;
  using jac_pol_t = ode::impl::JacobianStandardPolicyBdf<state_t, jac_t>;
  auto stepperObj = ode::create_bdf2_stepper(problemObj, y, res_pol_t(), jac_pol_t());

  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver = nonlinearsolvers::create_newton_raphson(stepperObj,y,linSolverObj);

  // integrate in time
  ode::step_count_type nSteps = 4;
  double dt = 0.01;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, nSteps, NonLinSolver);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  problemObj.analyticAdvanceBackEulerNSteps(dt, 1);
  problemObj.analyticAdvanceBDF2NSteps(dt, 3);
  std::cout << std::setprecision(14) << problemObj.y << "\n";
  EXPECT_DOUBLE_EQ(y(0), problemObj.y(0));
  EXPECT_DOUBLE_EQ(y(1), problemObj.y(1));
  EXPECT_DOUBLE_EQ(y(2), problemObj.y(2));
}

TEST(ode, implicit_bdf2_custom_policy_specify_types)
{
  using namespace pressio;
  using problem_t = ode::testing::refAppForImpEigen;
  using state_t = typename problem_t::state_type;
  using res_t = typename problem_t::velocity_type;
  using jac_t = typename problem_t::jacobian_type;
  problem_t problemObj;
  state_t y = problemObj.getInitCond();

  using res_pol_t = ode::impl::ResidualStandardPolicyBdf<state_t, res_t>;
  using jac_pol_t = ode::impl::JacobianStandardPolicyBdf<state_t, jac_t>;
  auto stepperObj = ode::create_bdf2_stepper_partial_deduction<res_t, jac_t>(problemObj, y, res_pol_t(), jac_pol_t());

  // define solver
  using lin_solver_t = linearsolvers::Solver<
    linearsolvers::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver = pressio::nonlinearsolvers::create_newton_raphson(stepperObj,y,linSolverObj);

  // integrate in time
  ::pressio::ode::step_count_type nSteps = 4;
  double dt = 0.01;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, nSteps, NonLinSolver);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  problemObj.analyticAdvanceBackEulerNSteps(dt, 1);
  problemObj.analyticAdvanceBDF2NSteps(dt, 3);
  std::cout << std::setprecision(14) << problemObj.y << "\n";
  EXPECT_DOUBLE_EQ(y(0), problemObj.y(0));
  EXPECT_DOUBLE_EQ(y(1), problemObj.y(1));
  EXPECT_DOUBLE_EQ(y(2), problemObj.y(2));
}
