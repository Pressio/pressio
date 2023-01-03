
#include "pressio/solvers.hpp"

struct ValidSystem
{
  using state_type = Eigen::VectorXd;
  using residual_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

  state_type createState() const {
    return state_type(2);
  }

  residual_type createResidual() const {
    return residual_type(2);
  }

  jacobian_type createJacobian() const {
    return jacobian_type(2, 2);
  }

  void residual(const state_type& x,
                residual_type& res) const
  {
    res(0) =  x(0)*x(0)*x(0) + x(1) - 1.0;
    res(1) = -x(0) + x(1)*x(1)*x(1) + 1.0;
  }

  void jacobian(const state_type& x, jacobian_type& jac) const {
    jac.coeffRef(0, 0) = 3.0*x(0)*x(0);
    jac.coeffRef(0, 1) =  1.0;
    jac.coeffRef(1, 0) = -1.0;
    jac.coeffRef(1, 1) = 3.0*x(1)*x(1);
  }
};

int main()
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  using namespace pressio;
  using problem_t  = ValidSystem;
  using state_t    = problem_t::state_type;
  using jacobian_t = problem_t::jacobian_type;

  problem_t sys;
  // my solution vector
  state_t y(2);
  // initial condition
  y(0) = 0.001; y(1) = 0.0001;

  // linear system
  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::LSCG, jacobian_t>;
  lin_solver_t linearSolverObj;

  auto nonLinSolver = nonlinearsolvers::create_newton_raphson(sys, linearSolverObj);
  nonLinSolver.solve(sys, y);

  std::string strOut = "PASSED";
  const auto e1 = std::abs(y(0) - (1.));
  const auto e2 = std::abs(y(1) - (0.));
  if (e1>1e-8 or e2>1e-8) strOut = "FAILED";
  std::cout <<  strOut << std::endl;
  std::cout << y << std::endl;

  pressio::log::finalize();
  return 0;
}
