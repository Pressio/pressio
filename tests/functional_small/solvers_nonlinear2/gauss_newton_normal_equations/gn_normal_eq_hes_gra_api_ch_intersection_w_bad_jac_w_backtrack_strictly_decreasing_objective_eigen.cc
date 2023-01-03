//
#include "pressio/solvers.hpp"
#include "eigen_ch_intersection_with_bad_jacobian.hpp"

int main()
{
  pressio::log::initialize(pressio::logto::fileAndTerminal, "log.txt");
  pressio::log::setVerbosity({pressio::log::level::trace, pressio::log::level::trace});

  using namespace pressio;
  using problem_t  = solvers::test::CircleHyperbolaIntersectionWithBadJacobianSystem4HessGradApi;
  using state_t    = problem_t::state_type;
  using jacobian_t = problem_t::hessian_type;

  problem_t sys;
  // my solution vector
  state_t y(2);
  // initial condition
  y(0) = 0.3; y(1) = 0.4;

  // linear system
  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::LSCG, jacobian_t>;
  lin_solver_t linearSolverObj;

  auto NonLinSolver = pressio::nonlinearsolvers::create_gauss_newton(sys, linearSolverObj);
  NonLinSolver.setUpdatingCriterion(pressio::nonlinearsolvers::Update::BacktrackStrictlyDecreasingObjective);
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
