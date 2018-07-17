
#include <gtest/gtest.h>

#include "matrix/concrete/core_matrix_sparse_serial_eigen.hpp"
#include "vector/concrete/core_vector_serial_eigen.hpp"
#include "experimental/solvers_nonlinear_iterative_factory.hpp"


TEST(solvers_nonlinear_iterative_eigen, solversTestNonLinearIterativeEigenRungeKutta)
{
  // Namespaces
  using namespace solvers;

  struct System {

    // Matrix typedefs
    using matrix_n_t = Eigen::SparseMatrix<double>;
    using matrix_w_t = core::matrix<matrix_n_t>;

    // Vector typedefs
    using vector_n_t = Eigen::VectorXd;
    using vector_w_t = core::vector<vector_n_t>;

    typedef vector_w_t state_type;
    typedef matrix_w_t matrix_type;


    void residual(const vector_w_t& x, vector_w_t& res) {
      res[0] =  x[0]*x[0]*x[0] + x[1] - 1.0;
      res[1] = -x[0] + x[1]*x[1]*x[1] + 1.0;
    }


    auto residual(const vector_w_t& x) {
      vector_w_t res(2);
      res[0] =  x[0]*x[0]*x[0] + x[1] - 1.0;
      res[1] = -x[0] + x[1]*x[1]*x[1] + 1.0;
      return res;
    }

    
    void jacobian(const vector_w_t& x, matrix_w_t& jac) {
      jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
      jac.data()->coeffRef(0, 1) =  1.0;
      jac.data()->coeffRef(1, 0) = -1.0;
      jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1]; 
    }


    auto jacobian(const vector_w_t& x) {
      matrix_w_t jac(2, 2);
      jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
      jac.data()->coeffRef(0, 1) =  1.0;
      jac.data()->coeffRef(1, 0) = -1.0;
      jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1]; 
      return jac;
    }
  };

  // Define a system to solve
  System system; 

  // Solve nonlinear system using 
  auto solver = NonlinearIterativeSolvers::createSolver<nonlinear::NewtonRaphson>(system);
 // auto x = solver.solve<linear::CG>(b);
  
  // Expectations
  // EXPECT_NEAR( x[0],  1.0, 1e-14 );
//  EXPECT_NEAR( x[1], -1.0, 1e-14 );
}

#if 0
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

#endif
