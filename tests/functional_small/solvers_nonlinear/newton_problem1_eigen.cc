
#include <gtest/gtest.h>
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_newton.hpp"
#include "./problems/problem1.hpp"

template<class SystemType, bool logOn = false>
void run_impl(int reps, bool callSolveWithJustState = true)
{
  if constexpr(logOn){
    pressio::log::initialize(pressio::logto::terminal);
    pressio::log::setVerbosity({pressio::log::level::debug});
  }

  using namespace pressio;
  using problem_t  = SystemType;
  using state_t    = typename problem_t::state_type;
  using jacobian_t = typename problem_t::jacobian_type;

  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::LSCG, jacobian_t>;
  lin_solver_t linearSolverObj;

  problem_t sys;
  state_t y(2);
  auto nonLinSolver = create_newton_solver(sys, linearSolverObj);

  if (reps){
    y(0) = 0.001; y(1) = 0.0001;

    if (callSolveWithJustState){
      nonLinSolver.solve(y);
    }else{
      nonLinSolver.solve(sys, y);
    }
    const auto e1 = std::abs(y(0) - (1.));
    const auto e2 = std::abs(y(1) - (0.));
    std::cout << y << std::endl;
    ASSERT_TRUE(e1<1e-8);
    ASSERT_TRUE(e2<1e-8);
  }
  else{

    for (int i=0; i<reps; ++i){
      y(0) = 0.001; y(1) = 0.0001;
      if (callSolveWithJustState){
	nonLinSolver.solve(y);
      }else{
	nonLinSolver.solve(sys, y);
      }

      const auto e1 = std::abs(y(0) - (1.));
      const auto e2 = std::abs(y(1) - (0.));
      ASSERT_TRUE(e1<1e-8);
      ASSERT_TRUE(e2<1e-8);
    }
  }

  if constexpr(logOn){
    pressio::log::finalize();
  }
}

TEST(solvers_nonlinear, problem1){
  run_impl<pressio::solvers::test::Problem1, true>(1, false);
}

TEST(solvers_nonlinear, problem1_repeated_solve){
  run_impl<pressio::solvers::test::Problem1>(100, false);
}

TEST(solvers_nonlinear, problem1_call_solve_with_only_state){
  run_impl<pressio::solvers::test::Problem1, true>(1, true);
}

TEST(solvers_nonlinear, problem1_repeated_solve_call_solve_with_only_state){
  run_impl<pressio::solvers::test::Problem1>(100, true);
}

