
#include <gtest/gtest.h>
#include "SOLVERS_NONLINEAR"
#include "./problems/solvers_utest_serial_rosenbrock_N4.hpp"

TEST(solvers_nonlinear_least_squares, gn_conv_conditions){

  using namespace pressio;

  // does not matter which problem we use here
  using problem_t   = solvers::test::Rosenbrock4;
  using vec_t	  = typename problem_t::state_type;
  using mat_t	  = typename problem_t::jacobian_type;
  using state_t	  = vec_t;
  using hessian_t = mat_t;

  problem_t problem;
  state_t x(4);
  x[0] = -0.05; x[1] = 1.1; x[2] = 1.2; x[3] = 1.5;

  // linear solver type
  using solver_tag	= solvers::linear::iterative::LSCG;
  using linear_solver_t = solvers::iterative::EigenIterative<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  {
    // GaussNewton solver
    using converged_when_t_a = solvers::iterative::converged_when::completingNumMaxIters;

    using gn_t = solvers::iterative::GaussNewton
      <linear_solver_t, converged_when_t_a, problem_t, hessian_t>;
    gn_t GNSolver(problem, x, linSolver);
    GNSolver.setTolerance(1e-8);
    GNSolver.setMaxIterations(1);
    GNSolver.solve(problem, x);

    const auto solverConvStr = GNSolver.getConvergenceConditionDescription();
    EXPECT_EQ(solverConvStr, "complete max iters");
  }

  using norm_t = solvers::L1Norm;

  {
    // GaussNewton solver
    using converged_when_t_a = solvers::iterative::converged_when::absoluteNormCorrectionBelowTol<norm_t>;

    using gn_t = solvers::iterative::GaussNewton
      <linear_solver_t, converged_when_t_a, problem_t, hessian_t>;
    gn_t GNSolver(problem, x, linSolver);
    GNSolver.setTolerance(1e-8);
    GNSolver.setMaxIterations(1);
    GNSolver.solve(problem, x);

    const auto solverConvStr = GNSolver.getConvergenceConditionDescription();
    EXPECT_EQ(solverConvStr, "||dy|| < tol");
  }

  {
    // GaussNewton solver
    using converged_when_t_a = solvers::iterative::converged_when::absoluteNormResidualBelowTol<norm_t>;

    using gn_t = solvers::iterative::GaussNewton
      <linear_solver_t, converged_when_t_a, problem_t, hessian_t>;
    gn_t GNSolver(problem, x, linSolver);
    GNSolver.setTolerance(1e-8);
    GNSolver.setMaxIterations(1);
    GNSolver.solve(problem, x);

    const auto solverConvStr = GNSolver.getConvergenceConditionDescription();
    EXPECT_EQ(solverConvStr, "||R|| < tol");
  }


  {
    // GaussNewton solver
    using converged_when_t_a = solvers::iterative::converged_when::relativeNormResidualBelowTol<norm_t>;

    using gn_t = solvers::iterative::GaussNewton
      <linear_solver_t, converged_when_t_a, problem_t, hessian_t>;
    gn_t GNSolver(problem, x, linSolver);
    GNSolver.setTolerance(1e-8);
    GNSolver.setMaxIterations(1);
    GNSolver.solve(problem, x);

    const auto solverConvStr = GNSolver.getConvergenceConditionDescription();
    EXPECT_EQ(solverConvStr, "||R||(r) < tol");
  }

  {
    // GaussNewton solver
    using converged_when_t_a = solvers::iterative::converged_when::absoluteNormGradientBelowTol<norm_t>;

    using gn_t = solvers::iterative::GaussNewton
      <linear_solver_t, converged_when_t_a, problem_t, hessian_t>;
    gn_t GNSolver(problem, x, linSolver);
    GNSolver.setTolerance(1e-8);
    GNSolver.setMaxIterations(1);
    GNSolver.solve(problem, x);

    const auto solverConvStr = GNSolver.getConvergenceConditionDescription();
    EXPECT_EQ(solverConvStr, "||J^T R|| < tol");
  }

  {
    // GaussNewton solver
    using converged_when_t_a = solvers::iterative::converged_when::relativeNormGradientBelowTol<norm_t>;

    using gn_t = solvers::iterative::GaussNewton
      <linear_solver_t, converged_when_t_a, problem_t, hessian_t>;
    gn_t GNSolver(problem, x, linSolver);
    GNSolver.setTolerance(1e-8);
    GNSolver.setMaxIterations(1);
    GNSolver.solve(problem, x);

    const auto solverConvStr = GNSolver.getConvergenceConditionDescription();
    EXPECT_EQ(solverConvStr, "||J^T R||(r) < tol");
  }

}
