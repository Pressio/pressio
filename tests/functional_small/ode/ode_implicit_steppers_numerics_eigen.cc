
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"
#include "pressio_ode_implicit.hpp"
#include "testing_apps.hpp"

TEST(ode, implicit_euler_policy_default_created)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  using res_t = typename app_t::velocity_type;
  using jac_t = typename app_t::jacobian_type;
  app_t appObj;

  state_t y(3);
  y = appObj.getInitCond();

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF1, state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  using lin_solver_t = linearsolvers::Solver<
    linearsolvers::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;

  auto NonLinSolver = pressio::nonlinearsolvers::create_newton_raphson(
    stepperObj,y,linSolverObj);

  ::pressio::ode::step_count_type nSteps = 2;
  double dt = 0.01;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, nSteps, NonLinSolver);
  std::cout << std::setprecision(14) << y << "\n";

  appObj.analyticAdvanceBackEulerNSteps(dt, nSteps);
  EXPECT_DOUBLE_EQ(y(0), appObj.y(0));
  EXPECT_DOUBLE_EQ(y(1), appObj.y(1));
  EXPECT_DOUBLE_EQ(y(2), appObj.y(2));
}

TEST(ode, implicit_euler_guesserLambda)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  using res_t = typename app_t::velocity_type;
  using jac_t = typename app_t::jacobian_type;
  app_t appObj;

  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF1, state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  // define solver
  using lin_algo_t = linearsolvers::iterative::Bicgstab;
  using lin_solver_t = linearsolvers::Solver<lin_algo_t, jac_t>;
  lin_solver_t linSolverObj;

  auto NonLinSolver = pressio::nonlinearsolvers::create_newton_raphson(stepperObj,y,linSolverObj);
  NonLinSolver.setMaxIterations(0);

  // integrate in time
  const auto testLambda = [](const ode::step_count_type & step,
             const double & time,
             state_t & yIn)
          {
            yIn(0) = -22.; yIn(1) = -26.; yIn(2) = -28.;
          };

  double dt = 0.01;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, 1, testLambda, NonLinSolver);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  EXPECT_DOUBLE_EQ(y(0), -22.0);
  EXPECT_DOUBLE_EQ(y(1), -26.0);
  EXPECT_DOUBLE_EQ(y(2), -28.0);
}

TEST(ode, implicit_euler_custom_policy)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  using res_t = typename app_t::velocity_type;
  using jac_t = typename app_t::jacobian_type;
  app_t appObj;
  state_t y(3);
  y = appObj.getInitCond();

  // define policies and stepper
  using res_pol_t = ode::impl::ResidualStandardPolicyBdf<state_t, res_t>;
  using jac_pol_t = ode::impl::JacobianStandardPolicyBdf<state_t, jac_t>;

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF1, state_t, res_t, jac_t, app_t, res_pol_t, jac_pol_t>;
  stepper_t stepperObj(y, appObj, res_pol_t(), jac_pol_t());

  //**********************
  // define solver
  //**********************
  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver = pressio::nonlinearsolvers::create_newton_raphson(
      stepperObj,y,linSolverObj);

  // integrate in time
  ::pressio::ode::step_count_type nSteps = 2;
  double dt = 0.01;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, nSteps, NonLinSolver);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  appObj.analyticAdvanceBackEulerNSteps(dt, nSteps);

  EXPECT_DOUBLE_EQ(y(0), appObj.y(0));
  EXPECT_DOUBLE_EQ(y(1), appObj.y(1));
  EXPECT_DOUBLE_EQ(y(2), appObj.y(2));
}



namespace
{
struct CustomUpdate
{
  void reset(){}

  template<typename system_t, typename state_t, typename solver_t>
  void operator()(const system_t & sys,
      state_t & state,
      solver_t & solver)
  {
    PRESSIOLOG_DEBUG("custom update");
    const auto & correction = solver.correctionCRef();
    std::cout << *state.data() << std::endl;
    ::pressio::ops::update(state, 1., correction, 0.);
    std::cout << *state.data() << std::endl;
  }
};
}

TEST(ode, implicit_euler_policy_default_created_custom_update)
{
  // this test is trivial, just aimed at testing
  // we can pass the custom update via the ode integrate
  // to the solver
  // the custom update sets correction to zero so that solution should never change

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::trace});

  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  using res_t = typename app_t::velocity_type;
  using jac_t = typename app_t::jacobian_type;
  app_t appObj;

  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 4.;

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF1, state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  // define solver
  using lin_solver_t = linearsolvers::Solver<
      linearsolvers::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;

  auto NonLinSolver = pressio::nonlinearsolvers::create_newton_raphson(stepperObj,y,linSolverObj);
  NonLinSolver.setMaxIterations(1);

  // integrate in time
  double dt = 0.01;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, 2, NonLinSolver, CustomUpdate{});
  std::cout << std::setprecision(14) << *y.data() << "\n";

  EXPECT_DOUBLE_EQ(y(0), 1.);
  EXPECT_DOUBLE_EQ(y(1), 2.);
  EXPECT_DOUBLE_EQ(y(2), 4.);
  pressio::log::finalize();
}


TEST(ode, implicit_euler_guesserLambdaCustomUpdate)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  using res_t = typename app_t::velocity_type;
  using jac_t = typename app_t::jacobian_type;
  app_t appObj;

  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF1, state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  using lin_algo_t = linearsolvers::iterative::Bicgstab;
  using lin_solver_t = linearsolvers::Solver<lin_algo_t, jac_t>;
  lin_solver_t linSolverObj;

  auto NonLinSolver = pressio::nonlinearsolvers::create_newton_raphson(stepperObj,y,linSolverObj);
  NonLinSolver.setMaxIterations(0);

  // integrate in time
  const auto testLambda = [](const ode::step_count_type & step,
             const double & time, 
             state_t & yIn)
          {
            yIn(0) = -22.; yIn(1) = -26.; yIn(2) = -28.;
          };

  double dt = 0.01;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, 5, testLambda, NonLinSolver, CustomUpdate{});
  std::cout << std::setprecision(14) << *y.data() << "\n";

  EXPECT_DOUBLE_EQ(y(0), -22.0);
  EXPECT_DOUBLE_EQ(y(1), -26.0);
  EXPECT_DOUBLE_EQ(y(2), -28.0);
}


TEST(ode, implicit_bdf2_policy_default_created)
{
  using namespace pressio;

  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  using res_t = typename app_t::velocity_type;
  using jac_t = typename app_t::jacobian_type;
  app_t appObj;
  state_t y(3);
  y = appObj.getInitCond();

  // define auxiliary stepper
  using aux_stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF1, state_t, res_t, jac_t, app_t>;
  aux_stepper_t stepperAux(y, appObj);

  // bdf2 stepper
  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF2, state_t, res_t, jac_t, app_t, aux_stepper_t>;
  stepper_t stepperObj(y, appObj, stepperAux);

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

  appObj.analyticAdvanceBackEulerNSteps(dt, 1);
  appObj.analyticAdvanceBDF2NSteps(dt, 3);
  std::cout << std::setprecision(14) << appObj.y << "\n";
  EXPECT_DOUBLE_EQ(y(0), appObj.y(0));
  EXPECT_DOUBLE_EQ(y(1), appObj.y(1));
  EXPECT_DOUBLE_EQ(y(2), appObj.y(2));
}


TEST(ode, implicit_bdf2_custom_policy)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  using res_t = typename app_t::velocity_type;
  using jac_t = typename app_t::jacobian_type;
  app_t appObj;
  state_t y(3);
  y = appObj.getInitCond();

  // define auxiliary policies and stepper
  using res_pol_t = ode::impl::ResidualStandardPolicyBdf<state_t, res_t>;
  using jac_pol_t = ode::impl::JacobianStandardPolicyBdf<state_t, jac_t>;

  using aux_stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF1,
    state_t, res_t, jac_t, app_t, res_pol_t, jac_pol_t>;
  aux_stepper_t stepperAux(y, appObj, res_pol_t(), jac_pol_t());

  // stepper for BDF2
  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF2,
    state_t, res_t, jac_t, app_t, aux_stepper_t, res_pol_t, jac_pol_t>;
  stepper_t stepperObj(y, appObj, res_pol_t(), jac_pol_t(), stepperAux);

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

  appObj.analyticAdvanceBackEulerNSteps(dt, 1);
  appObj.analyticAdvanceBDF2NSteps(dt, 3);
  std::cout << std::setprecision(14) << appObj.y << "\n";
  EXPECT_DOUBLE_EQ(y(0), appObj.y(0));
  EXPECT_DOUBLE_EQ(y(1), appObj.y(1));
  EXPECT_DOUBLE_EQ(y(2), appObj.y(2));
}


//==================================================
//==================================================
// CrankNicolson
//==================================================
//==================================================
namespace
{
struct CNTestApp
{
  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

  void velocity(const state_type & y,
    scalar_type t,
    velocity_type & R) const
  {
    auto sz = y.size();
    for (decltype(sz) i=0; i<sz; i++)
      R[i] = y[i];
  };

  void jacobian(const state_type & yIn,
      scalar_type t,
    jacobian_type & JJ) const
  {

  }

  velocity_type createVelocity() const
  {
    velocity_type R(3);
    return R;
  };
  jacobian_type createJacobian() const
  {
    jacobian_type JJ(3,3);
    return JJ;
  };
};
}

TEST(ode, implicit_crank_nicolson_constructor)
{
  using namespace pressio;
  using app_t = CNTestApp;
  using state_t = typename app_t::state_type;
  using res_t = typename app_t::velocity_type;
  using jac_t = typename app_t::jacobian_type;
  app_t appObj;
  state_t y(3);

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::CrankNicolson, state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  EXPECT_EQ(stepperObj.order(), 2);
}

