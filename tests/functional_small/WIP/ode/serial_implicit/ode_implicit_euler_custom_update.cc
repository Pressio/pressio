
#include <gtest/gtest.h>
#include "pressio_ode_implicit.hpp"
#include "../reference_apps_for_testing.hpp"

namespace{
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

TEST(ode_implicit_euler, numericsStdPoliciesDefaultCreatedCustomUpdate)
{
  // this test is trivial, just aimed at testing
  // we can pass the custom update via the ode integrate
  // to the solver
  // the custom update sets correction to zero so that solution should never change

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::trace});

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
  y(0) = 1.; y(1) = 2.; y(2) = 4.;

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::Euler,
    state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  // define solver
  using lin_solver_t = solvers::linear::Solver<
      solvers::linear::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;

  auto NonLinSolver = pressio::solvers::nonlinear::createNewtonRaphson(stepperObj,y,linSolverObj);
  NonLinSolver.setMaxIterations(1);

  // integrate in time
  double dt = 0.01;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 2, NonLinSolver, CustomUpdate{});
  std::cout << std::setprecision(14) << *y.data() << "\n";

  EXPECT_DOUBLE_EQ(y(0), 1.);
  EXPECT_DOUBLE_EQ(y(1), 2.);
  EXPECT_DOUBLE_EQ(y(2), 4.);
  pressio::log::finalize();
}


TEST(ode_implicit_euler, guesserLambdaCustomUpdate)
{
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
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::Euler, state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  // define solver
  using lin_algo_t = solvers::linear::iterative::Bicgstab;
  using lin_solver_t = solvers::linear::Solver<lin_algo_t, jac_t>;
  lin_solver_t linSolverObj;

  auto NonLinSolver = pressio::solvers::nonlinear::createNewtonRaphson(stepperObj,y,linSolverObj);
  NonLinSolver.setMaxIterations(0);

  // integrate in time
  const auto testLambda = [](const ode::types::step_t & step,
  			     const double & time,
  			     state_t & yIn)
  			  {
  			    yIn(0) = -22.; yIn(1) = -26.; yIn(2) = -28.;
  			  };

  double dt = 0.01;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 5, testLambda, NonLinSolver, CustomUpdate{});
  std::cout << std::setprecision(14) << *y.data() << "\n";

  EXPECT_DOUBLE_EQ(y(0), -22.0);
  EXPECT_DOUBLE_EQ(y(1), -26.0);
  EXPECT_DOUBLE_EQ(y(2), -28.0);
}
