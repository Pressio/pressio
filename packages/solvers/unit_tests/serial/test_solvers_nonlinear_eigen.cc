
#include <gtest/gtest.h>

#include "SOLVERS_EXP"
#include "CORE_VECTOR"
#include "CORE_MATRIX"

TEST(solvers_nonlinear, simpleTest)
{
  using native_state_t = Eigen::VectorXd;
  using native_jac_t = Eigen::Matrix<double,-1,-1>;
  using state_t = core::vector<native_state_t>;
  using jac_t = core::matrix<native_jac_t>;
  using scalar_t = double;

  struct app{
    void residual(const state_t & x, state_t & res){
      res[0] = x[0]*x[0]*x[0] + x[1] - 1.;
      res[1] = -x[0] + x[1]*x[1]*x[1] + 1.;
    }
    void jacobian(const state_t & x, jac_t & jac){
      jac(0,0) = 3.0*x[0]*x[0];
      jac(0,1) = 1.;
      jac(1,0) = -1;
      jac(1,1) = 3.*x[1]*x[1];
    }  
  };

  // linear solver
  using lin_solve_t =
    solvers::experimental::linearSolver<jac_t,state_t,state_t>;
  lin_solve_t ls;
    
  // nonlinear solver
  using nonlin_solve_t =
    solvers::experimental::newtonRaphson<state_t,state_t,jac_t,lin_solve_t>;
  nonlin_solve_t nonls(ls);

  state_t y(2);
  y[0]=0.0015; y[1]=1.5;
  app obj;
  nonls.solve(y,obj);
  std::cout << y[0] << " " << y[1] << std::endl;
}
