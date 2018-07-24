
#include <gtest/gtest.h>

#include <iostream>

#include "experimental/solvers_norms.hpp"
#include "matrix/concrete/core_matrix_sparse_serial_eigen.hpp"
#include "vector/concrete/core_vector_serial_eigen.hpp"
#include "experimental/solvers_nonlinear_iterative_factory.hpp"

struct NonLinearSystem {

    // Matrix typedefs

    using matrix_n_t = Eigen::SparseMatrix<double>;
    using matrix_w_t = core::Matrix<matrix_n_t>;

    // Vector typedefs
    using vector_n_t = Eigen::VectorXd;
    using vector_w_t = core::Vector<vector_n_t>;

    typedef vector_w_t state_type;
    typedef matrix_w_t matrix_type;


    void residual(const vector_w_t& x, vector_w_t& res) const {
      res[0] =  x[0]*x[0]*x[0] + x[1] - 1.0;
      res[1] = -x[0] + x[1]*x[1]*x[1] + 1.0;
    }


    auto residual(const vector_w_t& x) const {
      vector_w_t res(2);
      res[0] =  x[0]*x[0]*x[0] + x[1] - 1.0;
      res[1] = -x[0] + x[1]*x[1]*x[1] + 1.0;
      return res;
    }

    
    void jacobian(const vector_w_t& x, matrix_w_t& jac) const {
      jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
      jac.data()->coeffRef(0, 1) =  1.0;
      jac.data()->coeffRef(1, 0) = -1.0;
      jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1]; 
    }


    auto jacobian(const vector_w_t& x) const {
      matrix_w_t jac(2, 2);
      jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
      jac.data()->coeffRef(0, 1) =  1.0;
      jac.data()->coeffRef(1, 0) = -1.0;
      jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1]; 
      return jac;
    }
};


TEST(solvers_nonlinear_iterative_eigen, solversTestNonLinearIterativeEigenRungeKutta)
{
  // Namespaces
  using namespace solvers;

  // Define a system to solve
  NonLinearSystem system;
  core::Vector<Eigen::VectorXd> b(2);

  // Initialize b
  b[0] = 0.15;
  b[1] = 0.5; 

  // Solve nonlinear system using 
  auto solver = NonlinearIterativeSolvers::createSolver<nonlinear::NewtonRaphson>();
  auto y = solver.solve<linear::Bicgstab>(system, b);

  // Expectations
  EXPECT_NEAR( y[0],  1.0, 1e-8 );
  EXPECT_NEAR( y[1],  0.0, 1e-8 );
 
}


#if 0
TEST(solvers_nonlinear, simpleTest)
{
  using native_state_t = Eigen::VectorXd;
  using native_jac_t = Eigen::Matrix<double,-1,-1>;
  using state_t = core::Vector<native_state_t>;
  using jac_t = core::Matrix<native_jac_t>;
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
#endif


