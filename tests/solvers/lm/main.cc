#include <iostream>
#include "pressio_solvers.hpp"

struct NonLinearLeastSquareSystem {
    // Matrix typedefs
    //using matrix_n_t = Eigen::SparseMatrix<double>;
    using matrix_n_t  = Eigen::Matrix<double, -1, -1>;
    using matrix_w_t = pressio::containers::Matrix<matrix_n_t>;

    // Vector typedefs
    using vector_n_t = Eigen::VectorXd;
  
    using vector_w_t = pressio::containers::Vector<vector_n_t>;

    typedef vector_w_t vector_type;
    typedef matrix_w_t matrix_type;


    using scalar_type = double;
    using state_type = vector_w_t;
    using residual_type = vector_w_t;
    using jacobian_type = matrix_w_t; 

    void residual(const vector_w_t& x, vector_w_t& res) const {
      res[0] = x[0] - x[1]*(2. - x[1]*(5. - x[1]) ) - 13.;
      res[1] = x[0] - x[1]*(14. - x[1]*(1. + x[1]) ) - 29.;

    }


    vector_w_t residual(const vector_w_t& x) const {
      vector_w_t res(2);
      this->residual(x, res);
      return res;
    }


    void jacobian(const vector_w_t& x, matrix_w_t& jac) const {
        jac.data()->coeffRef(0,0) = 1;
        jac.data()->coeffRef(0,1) = -x[1]*(2.*x[1] - 5.) + (5. - x[1])*x[1] - 2.;
        jac.data()->coeffRef(1,0) = 1.;
        jac.data()->coeffRef(1,1) = x[1]*(x[1] + 1.) - (-2.*x[1] - 1.)*x[1] - 14.;
    }


    matrix_w_t jacobian(const vector_w_t& x) const {
      matrix_w_t jac(2, 2);
      this->jacobian(x, jac);
      return jac;
    }

    const double coeff[5] = {1.0, 2.0, 4.0, 5.0, 8.0};
    const double vansw[5] = {3.2939, 4.2699, 7.1749, 9.3008, 20.259};
};

int main() {

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

  using lm_schedule_policy_tag = pressio::solvers::iterative::lm::SchedulePolicy2;
  using system_t = NonLinearLeastSquareSystem;
  using lmsolver_t   = pressio::solvers::nonlinear::LM<system_t, linear_solver_t,lm_schedule_policy_tag>;
  vector_w_t x0(2);
  x0[0] = 0.5;
  x0[1] = -2.;
  NonLinearLeastSquareSystem sys;
  lmsolver_t solver(sys, x0, linSolverObj);
  solver.solve(sys, x0);

  std::cout << "The solution of the nonlinear system is: " << std::endl << *x0.data() << std::endl;
  return 0;
}
