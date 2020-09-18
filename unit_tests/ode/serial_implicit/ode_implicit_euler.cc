
#include <gtest/gtest.h>
#include "pressio_ode.hpp"
#include "../reference_apps_for_testing.hpp"


TEST(ode_implicit_euler, numericsStdPoliciesDefaultCreated){
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  using nveloc_t = typename app_t::velocity_type;
  using njacobian_t = typename app_t::jacobian_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  using res_t = containers::Vector<nveloc_t>;
  using jac_t = containers::SparseMatrix<njacobian_t>;
  state_t y(3);
  *y.data() = appObj.getInitCond();

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::Euler,
    state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  // define solver
  using lin_solver_t = solvers::linear::Solver<
      solvers::linear::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;

  using nl_solver_t = pressio::solvers::nonlinear::composeNewtonRaphson_t<
    stepper_t, pressio::solvers::nonlinear::DefaultUpdate,
    lin_solver_t>;
  nl_solver_t NonLinSolver(stepperObj, y, linSolverObj);

  // integrate in time
  ::pressio::ode::types::step_t nSteps = 2;
  double dt = 0.01;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 2, NonLinSolver);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  appObj.analyticAdvanceBackEulerNSteps(dt, nSteps);

  EXPECT_DOUBLE_EQ(y[0], appObj.y[0]);
  EXPECT_DOUBLE_EQ(y[1], appObj.y[1]);
  EXPECT_DOUBLE_EQ(y[2], appObj.y[2]);
}


TEST(ode_implicit_euler, guesserLambda){
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  using nveloc_t = typename app_t::velocity_type;
  using njacobian_t = typename app_t::jacobian_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  using res_t = containers::Vector<nveloc_t>;
  using jac_t = containers::SparseMatrix<njacobian_t>;
  state_t y(3);//appObj.y0);
  y[0] = 1.; y[1] = 2.; y[2] = 3.;

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::Euler, state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  // define solver
  using lin_algo_t = solvers::linear::iterative::Bicgstab;
  using lin_solver_t = solvers::linear::Solver<lin_algo_t, jac_t>;
  lin_solver_t linSolverObj;

  using nl_solver_t = pressio::solvers::nonlinear::composeNewtonRaphson_t<
    stepper_t, pressio::solvers::nonlinear::DefaultUpdate,
    lin_solver_t>;
  nl_solver_t NonLinSolver(stepperObj, y, linSolverObj);
  NonLinSolver.setMaxIterations(0);
  // using nonlinear_solver_t = pressio::solvers::NewtonRaphson<stepper_t, lin_solver_t, double>;
  // nonlinear_solver_t solverO(stepperObj, y, linSolverObj);

  // integrate in time
  const auto testLambda = [](const ode::types::step_t & step,
  			     const double & time,
  			     state_t & yIn)
  			  {
  			    yIn[0] = -22.; yIn[1] = -26.; yIn[2] = -28.;
  			  };

  double dt = 0.01;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 1, NonLinSolver, testLambda);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  EXPECT_DOUBLE_EQ(y[0], -22.0);
  EXPECT_DOUBLE_EQ(y[1], -26.0);
  EXPECT_DOUBLE_EQ(y[2], -28.0);
}


TEST(ode_implicit_euler, numericsStdResidualPolPassedByUser){
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  using nveloc_t = typename app_t::velocity_type;
  using njacobian_t = typename app_t::jacobian_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  using res_t = containers::Vector<nveloc_t>;
  using jac_t = containers::SparseMatrix<njacobian_t>;
  state_t y(3);
  *y.data() = appObj.getInitCond();

  //**********************
  // define policies and stepper
  //**********************
  using res_pol_t = ode::implicitmethods::policy::ResidualStandardPolicy<state_t, app_t, res_t>;
  using jac_pol_t = ode::implicitmethods::policy::JacobianStandardPolicy<state_t, app_t, jac_t>;

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::Euler, state_t, res_t, jac_t, app_t, res_pol_t, jac_pol_t>;
  stepper_t stepperObj(y, appObj, res_pol_t(), jac_pol_t());

  //**********************
  // define solver
  //**********************
  using lin_solver_t = solvers::linear::Solver<
      solvers::linear::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;

  using nl_solver_t = pressio::solvers::nonlinear::composeNewtonRaphson_t<
    stepper_t, pressio::solvers::nonlinear::DefaultUpdate,
    lin_solver_t>;
  nl_solver_t NonLinSolver(stepperObj, y, linSolverObj);

  // integrate in time
  ::pressio::ode::types::step_t nSteps = 2;
  double dt = 0.01;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, nSteps, NonLinSolver);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  appObj.analyticAdvanceBackEulerNSteps(dt, nSteps);

  EXPECT_DOUBLE_EQ(y[0], appObj.y[0]);
  EXPECT_DOUBLE_EQ(y[1], appObj.y[1]);
  EXPECT_DOUBLE_EQ(y[2], appObj.y[2]);
}

