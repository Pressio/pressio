
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear.hpp"

namespace{
struct MySystem
{
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;;
  using residual_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

  residual_type createResidual() const { return residual_type(2); }
  jacobian_type createJacobian() const { return jacobian_type(2, 2); }

  void residual(const state_type& x, residual_type& res) const{}
  void jacobian(const state_type& x, jacobian_type& jac) const{}
};
}

TEST(solvers_nonlinear, tolerances)
{
  MySystem sysObj;
  typename MySystem::state_type y(2);

  using lin_solver_t = pressio::linearsolvers::Solver<
    pressio::linearsolvers::iterative::LSCG, typename MySystem::jacobian_type>;
  lin_solver_t linearSolverObj;
  auto solver = pressio::nonlinearsolvers::create_newton_raphson(sysObj, y, linearSolverObj);

  // if we don't specify anything, the tolerance should be defaulted to 0.000001
  const double defaultTol = 0.000001;
  ASSERT_EQ( solver.correctionAbsoluteTolerance(), defaultTol );
  ASSERT_EQ( solver.correctionRelativeTolerance(), defaultTol );
  ASSERT_EQ( solver.residualAbsoluteTolerance(),   defaultTol );
  ASSERT_EQ( solver.residualRelativeTolerance(),   defaultTol );
  ASSERT_EQ( solver.gradientAbsoluteTolerance(),   defaultTol );
  ASSERT_EQ( solver.gradientRelativeTolerance(),   defaultTol );

  // set each and check
  solver.setCorrectionAbsoluteTolerance(0.01);
  solver.setCorrectionRelativeTolerance(0.02);
  solver.setResidualAbsoluteTolerance(0.03);
  solver.setResidualRelativeTolerance(0.04);
  solver.setGradientAbsoluteTolerance(0.05);
  solver.setGradientRelativeTolerance(0.06);
  ASSERT_EQ( solver.correctionAbsoluteTolerance(), 0.01 );
  ASSERT_EQ( solver.correctionRelativeTolerance(), 0.02 );
  ASSERT_EQ( solver.residualAbsoluteTolerance(),   0.03 );
  ASSERT_EQ( solver.residualRelativeTolerance(),   0.04 );
  ASSERT_EQ( solver.gradientAbsoluteTolerance(),   0.05 );
  ASSERT_EQ( solver.gradientRelativeTolerance(),   0.06 );
}
