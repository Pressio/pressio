#include <iostream>
#include "pressio_solvers.hpp"

struct NonLinearLeastSquareSystem 
{
  using matrix_n_t  = Eigen::Matrix<double, -1, -1>;
  using matrix_w_t = pressio::containers::Matrix<matrix_n_t>;
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = pressio::containers::Vector<vector_n_t>;
  typedef vector_w_t vector_type;
  typedef matrix_w_t matrix_type;

  using scalar_type = double;
  using state_type = vector_w_t;
  using residual_type = vector_w_t;
  using jacobian_type = matrix_w_t;

  void residual(const state_type& x, 
    residual_type & res,
		::pressio::Norm normKind, 
    scalar_type & normResidual) const 
  {
    res[0] = x[0] - x[1]*(2. - x[1]*(5. - x[1]) ) - 13.;
    res[1] = x[0] - x[1]*(14. - x[1]*(1. + x[1]) ) - 29.;

    if (normKind == pressio::Norm::L2)
      normResidual = res.data()->norm();
  }

  void jacobian(const state_type& x, jacobian_type & jac) const {
    jac.data()->coeffRef(0,0) = 1.;
    jac.data()->coeffRef(0,1) = -x[1]*(2.*x[1] - 5.) + (5. - x[1])*x[1] - 2.;
    jac.data()->coeffRef(1,0) = 1.;
    jac.data()->coeffRef(1,1) = x[1]*(x[1] + 1.) - (-2.*x[1] - 1.)*x[1] - 14.;
  }

  residual_type createResidual() const {
    return residual_type(2);
  }

  jacobian_type createJacobian() const {
    return jacobian_type(2, 2);
  }
};

int main() {
  std::string checkStr = "PASSED";
  // Namespaces
  using namespace pressio;
  using namespace pressio::solvers;

  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = containers::Vector<vector_n_t>;
  using eig_dyn_mat  = Eigen::Matrix<double, -1, -1>;
  using hessian_t  = pressio::containers::Matrix<eig_dyn_mat>;

  using solver_tag   = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t = pressio::solvers::linear::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  using system_t = NonLinearLeastSquareSystem;
  NonLinearLeastSquareSystem sys;

  // LM with default update
  vector_w_t x0(2);
  {
    x0[0] = 0.5;
    x0[1] = -2.;

    using lmsolver = pressio::solvers::nonlinear::composeLevenbergMarquardt_t<
      system_t, pressio::solvers::nonlinear::LMDefaultUpdate,
      linear_solver_t>;
    lmsolver solver1(sys, x0, linSolverObj);

    solver1.setTolerance(1e-15);
    solver1.solve(sys, x0);
  }

  // LM with a different update
  vector_w_t x1(2);
  {
    x1[0] = 0.5;
    x1[1] = -2.;

    using lmsolver = pressio::solvers::nonlinear::composeLevenbergMarquardt<
      system_t, pressio::solvers::nonlinear::LMUpdateSchedule2,
      linear_solver_t>::type;
    lmsolver solver2(sys, x0, linSolverObj);
    solver2.setTolerance(1e-15);
    solver2.solve(sys, x1);
  }

  // check solution
  vector_w_t xstar(2);
  xstar[0] = 11.412779;
  xstar[1] = -0.896805;

  for (int i=0; i< 2; i++){
    if (abs((*x0.data())(i) - (*xstar.data())(i)) > 1e-6){
      checkStr = "FAILED";
      std::cout << "Default policy failed" << std::endl;
    }
    if (abs((*x1.data())(i) - (*xstar.data())(i)) > 1e-6){
      checkStr = "FAILED";
      std::cout << "Policy 2 failed" << std::endl;
    }
  }

  std::cout << checkStr << std::endl;
  return 0;
}
