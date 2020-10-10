
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"

TEST(solvers, nonlinear_tolerances)
{
  struct system
  {
    using matrix_n_t = Eigen::SparseMatrix<double>;
    using matrix_w_t = pressio::containers::SparseMatrix<matrix_n_t>;
    using vector_n_t = Eigen::VectorXd;
    using vector_w_t = pressio::containers::Vector<vector_n_t>;

    using scalar_type     = double;
    using state_type	= vector_w_t;
    using residual_type	= state_type;
    using jacobian_type	= matrix_w_t;

    residual_type createResidual() const { return residual_type(2); }
    jacobian_type createJacobian() const { return jacobian_type(2, 2); }

    void residual(const state_type& x, residual_type& res) const{}
    void jacobian(const state_type& x, jacobian_type& jac) const{}
  };

  system sysObj;
  typename system::state_type y(2);

  // it does not matter here what types we use, 
  using lin_solver_t = pressio::solvers::linear::Solver<
    pressio::solvers::linear::iterative::LSCG, typename system::jacobian_type>;
  lin_solver_t linearSolverObj;
  auto solver = pressio::solvers::nonlinear::createNewtonRaphson(sysObj, y, linearSolverObj);

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
