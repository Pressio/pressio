//
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
    res(0) =  x(0)*x(0) + x(1)*x(1) - 4.0;
    res(1) = x(0)*x(1)  - 1.0;
  }

  void jacobian(const state_type& x, jacobian_type& jac) const {
    jac.coeffRef(0, 0) = 2*x(0);
    jac.coeffRef(0, 1) = 2*x(1);
    // Have incorrect entries to mimic requirement for line search
    jac.coeffRef(1, 0) = 0.1*x(1);
    jac.coeffRef(1, 1) = 0.1*x(0);
  }
};

int main()
{
  pressio::log::initialize(pressio::logto::fileAndTerminal, "log.txt");
  pressio::log::setVerbosity({pressio::log::level::trace, pressio::log::level::trace});

  using namespace pressio;
  using problem_t  = ValidSystem;
  using state_t    = problem_t::state_type;
  using jacobian_t = problem_t::jacobian_type;

  problem_t sys;
  // my solution vector
  state_t y(2);
  // initial condition
  y(0) = 0.3; y(1) = 0.4;

  // linear system
  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::LSCG, jacobian_t>;
  lin_solver_t linearSolverObj;

  auto NonLinSolver = nonlinearsolvers::create_newton_raphson(sys, linearSolverObj);
  NonLinSolver.setUpdatingCriterion(pressio::nonlinearsolvers::Update::BasicBacktrack);
  NonLinSolver.solve(sys, y);

  std::string strOut = "PASSED";
  const auto e1 = std::abs(y(0) - (1.93185165));
  const auto e2 = std::abs(y(1) - (0.51763809));
  if (e1>1e-7 or e2>1e-7) strOut = "FAILED";

  std::cout <<  strOut << std::endl;
  std::cout << y << std::endl;
  pressio::log::finalize();
  return 0;
}
