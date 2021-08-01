
#include <gtest/gtest.h>
#include "pressio_ode_implicit.hpp"

namespace{
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

struct MySolver
{

  template <class T, class state_t>
  void solve(T &, state_t & state)
  {

  }
};
}

TEST(ode_implicit_crank_nicolson, constructor)
{
  using namespace pressio;
  using app_t = CNTestApp;
  using nstate_t = typename app_t::state_type;
  using nveloc_t = typename app_t::velocity_type;
  using njacobian_t = typename app_t::jacobian_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  using res_t = containers::Vector<nveloc_t>;
  using jac_t = containers::SparseMatrix<njacobian_t>;
  state_t y(3);

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::CrankNicolson,
    state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  EXPECT_EQ(stepperObj.order(), 2);
  // // define solver
  // using lin_solver_t = solvers::linear::Solver<
  //     solvers::linear::iterative::Bicgstab, jac_t>;
  // lin_solver_t linSolverObj;

  // auto NonLinSolver = pressio::solvers::nonlinear::createNewtonRaphson(stepperObj,y,linSolverObj);

  // // integrate in time
  // ::pressio::ode::step_type nSteps = 2;
  // double dt = 0.01;
  // ode::advanceNSteps(stepperObj, y, 0.0, dt, 2, NonLinSolver);
  // std::cout << std::setprecision(14) << *y.data() << "\n";

  // appObj.analyticAdvanceBackEulerNSteps(dt, nSteps);

  // EXPECT_DOUBLE_EQ(y(0), appObj.y(0));
  // EXPECT_DOUBLE_EQ(y(1), appObj.y(1));
  // EXPECT_DOUBLE_EQ(y(2), appObj.y(2));
}
