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
        jac.data()->coeffRef(0,0) = 1.;
        jac.data()->coeffRef(0,1) = -x[1]*(2.*x[1] - 5.) + (5. - x[1])*x[1] - 2.;
        jac.data()->coeffRef(1,0) = 1.;
        jac.data()->coeffRef(1,1) = x[1]*(x[1] + 1.) - (-2.*x[1] - 1.)*x[1] - 14.;
    }


    matrix_w_t jacobian(const vector_w_t& x) const {
      matrix_w_t jac(2, 2);
      this->jacobian(x, jac);
      return jac;
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

  using lm_schedule_policy_tag1 = pressio::solvers::iterative::lm::SchedulePolicyDefault;
  using lm_schedule_policy_tag2 = pressio::solvers::iterative::lm::SchedulePolicy2;
  // We can also build a policy and pass it to the solver. For example, let's say we want to run 
  // with lm_schedule_policy_tag_1 but with custom coefficients, we can do:
  //    "pressio::solvers::iterative::impl::LMSchedule<lm_schedule_policy_tag1,double> LMSchedule(2.,3.,0.2,0.8,1.);"
  // and pass to the solver below
  //
  using system_t = NonLinearLeastSquareSystem;
  using lmsolver_t1   = pressio::solvers::nonlinear::LM<system_t, linear_solver_t,lm_schedule_policy_tag1>;
  using lmsolver_t2   = pressio::solvers::nonlinear::LM<system_t, linear_solver_t,lm_schedule_policy_tag2>;

  vector_w_t x0(2);
  x0[0] = 0.5;
  x0[1] = -2.;
  NonLinearLeastSquareSystem sys;
  lmsolver_t1 solver1(sys, x0, linSolverObj);
  // if we wanted to pass a schedule policy, solver1(sys, x0, linSolverObj,LMSchedule);
  solver1.setTolerance(1e-15);
  solver1.solve(sys, x0);


  vector_w_t x1(2);
  x1[0] = 0.5;
  x1[1] = -2.;
  lmsolver_t2 solver2(sys, x0, linSolverObj);
  solver2.setTolerance(1e-15);
  solver2.solve(sys, x1);

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
