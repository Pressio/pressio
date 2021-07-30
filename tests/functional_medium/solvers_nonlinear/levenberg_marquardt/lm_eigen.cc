#include <iostream>
#include "pressio_solvers.hpp"

struct NonLinearLeastSquareSystem
{
  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using residual_type = Eigen::VectorXd;
  using jacobian_type = Eigen::Matrix<double, -1, -1>;

  void residual(const state_type& x, residual_type & res) const
  {
    res(0) = x(0) - x(1)*(2. - x(1)*(5. - x(1)) ) - 13.;
    res(1) = x(0) - x(1)*(14. - x(1)*(1. + x(1)) ) - 29.;
  }

  void jacobian(const state_type& x, jacobian_type & jac) const {
    jac.coeffRef(0,0) = 1.;
    jac.coeffRef(0,1) = -x(1)*(2.*x(1) - 5.) + (5. - x(1))*x(1) - 2.;
    jac.coeffRef(1,0) = 1.;
    jac.coeffRef(1,1) = x(1)*(x(1) + 1.) - (-2.*x(1) - 1.)*x(1) - 14.;
  }

  residual_type createResidual() const {
    return residual_type(2);
  }

  jacobian_type createJacobian() const {
    return jacobian_type(2, 2);
  }
};

int main()
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::trace});

  std::string checkStr = "PASSED";
  using vector_t = Eigen::VectorXd;
  using hessian_t  = Eigen::Matrix<double, -1, -1>;

  using solver_tag   = pressio::linearsolvers::iterative::LSCG;
  using linear_solver_t = pressio::linearsolvers::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  NonLinearLeastSquareSystem sys;

  // // LM with default update
  vector_t x0(2);
  {
    x0(0) = 0.5;
    x0(1) = -2.;
    auto solver = pressio::nonlinearsolvers::createLevenbergMarquardt(sys, x0, linSolverObj);
    solver.setUpdatingCriterion(pressio::nonlinearsolvers::update::LMSchedule1);
    solver.setTolerance(1e-15);
    solver.solve(sys, x0);
  }

  // LM with a different update
  vector_t x1(2);
  {
    x1(0) = 0.5;
    x1(1) = -2.;
    auto solver = pressio::nonlinearsolvers::createLevenbergMarquardt(sys, x0, linSolverObj);
    solver.setTolerance(1e-15);
    solver.setUpdatingCriterion(pressio::nonlinearsolvers::update::LMSchedule2);
    solver.solve(sys, x1);
  }

  // LM with default update, multiple solves
  // if we call solve twice in a row starting from same initial condition
  // we should get the same answer. if we don't it means the solver does
  // not reset its parameters correctly.
  vector_t x2a(2);
  vector_t x2b(2);
  {
    auto solver = pressio::nonlinearsolvers::createLevenbergMarquardt(sys, x0, linSolverObj);
    solver.setMaxIterations(4);
    solver.setUpdatingCriterion(pressio::nonlinearsolvers::update::LMSchedule1);

    x2a(0) = 0.5; x2a(1) = -2.;
    solver.solve(sys, x2a);

    x2b(0) = 0.5; x2b(1) = -2.;
    solver.solve(sys, x2b);

    const bool b1 = abs(x2a(0) - x2b(0)) < 1e-13;
    const bool b2 = abs(x2a(1) - x2b(1)) < 1e-13;
    if (!b1 or !b2){
      checkStr = "FAILED";
      std::cout << "LM sequential solves failed" << std::endl;
    }
  }

  // check solution
  vector_t xstar(2);
  xstar(0) = 11.412779;
  xstar(1) = -0.896805;

  for (int i=0; i< 2; i++){
    if (abs((x0)(i) - (xstar)(i)) > 1e-6){
      checkStr = "FAILED";
      std::cout << "Default policy failed" << std::endl;
    }
    if (abs((x1)(i) - (xstar)(i)) > 1e-6){
      checkStr = "FAILED";
      std::cout << "Policy 2 failed" << std::endl;
    }
  }

  std::cout << checkStr << std::endl;
  pressio::log::finalize();
  return 0;
}
