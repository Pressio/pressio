
#include <gtest/gtest.h>

#include <iostream>

#include "CONTAINERS_ALL"
// #include "matrix/concrete/containers_matrix_sparse_serial_eigen.hpp"
// #include "vector/concrete/containers_vector_serial_eigen.hpp"

#include <iostream>
#include "experimental/solvers_l2_vector_norm.hpp"
#include "experimental/solvers_nonlinear_base.hpp"
#include "experimental/solvers_nonlinear_traits.hpp"
#include "experimental/solvers_linear_factory.hpp"
#include "experimental/solvers_nonlinear_factory.hpp"
#include "experimental/solvers_policy_nonlinear_iterative.hpp"
#include "experimental/solvers_policy_nonlinear_leastsquare_iterative.hpp"


struct ValidSystemLeastSquares {

    // Matrix typedefs
    using matrix_n_t = Eigen::SparseMatrix<double>;
    using matrix_w_t = pressio::containers::Matrix<matrix_n_t>;

    // Vector typedefs
    using vector_n_t = Eigen::VectorXd;
    using vector_w_t = pressio::containers::Vector<vector_n_t>;

    typedef vector_w_t vector_type;
    typedef matrix_w_t matrix_type;


    void residual(const vector_w_t& x, vector_w_t& res) const {
      for (int i = 0; i < 5; i++) {
        res[i] = x[0] * exp(x[1] * coeff[i]) - vansw[i];
      }
    }


    vector_w_t residual(const vector_w_t& x) const {
      vector_w_t res(5);
      this->residual(x, res);
      return res;
    }


    void jacobian(const vector_w_t& x, matrix_w_t& jac) const {
      for (int i = 0; i < 5; i++) {
        double expval = exp(x[1] * coeff[i]);
        jac.data()->coeffRef(i, 0) = expval;
        jac.data()->coeffRef(i, 1) = x[0]*coeff[i]*expval;
      }
    }


    matrix_w_t jacobian(const vector_w_t& x) const {
      matrix_w_t jac(5, 2);
      this->jacobian(x, jac);
      return jac;
    }

    const double coeff[5] = {1.0, 2.0, 4.0, 5.0, 8.0};
    const double vansw[5] = {3.2939, 4.2699, 7.1749, 9.3008, 20.259};
};


TEST(solvers_non_linear_least_square_base, solversBaseSolveTest)
{
  using namespace pressio;
  using namespace pressio::solvers;

  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = containers::Vector<vector_n_t>;

  auto solver = NonLinearSolvers::createNonLinearIterativeLeastSquareSolver<nonlinearleastsquare::LevenbergMarquardt, linear::Bicgstab>();

  vector_w_t x0(2);
  x0[0] = 2.50;
  x0[1] = 0.25;

  ValidSystemLeastSquares sys;
  auto x = solver.solve(sys, x0);

  EXPECT_NEAR( x[0],  2.541, 1e-3 );
  EXPECT_NEAR( x[1],  0.259, 1e-3 );
}

/*
TEST(solvers_non_linear_base, solversBaseBadSolveTest)
{
  using namespace pressio;
  using namespace pressio::solvers;

  auto solver = NonLinearSolvers::createIterativeSolver<nonlinear::NewtonRaphson, linear::Bicgstab>();

  double left; int right;

  ASSERT_DEATH(solver.solve(left, right), "Error: either the nonlinear system or the solution hint is invalid.");
}


TEST(solvers_non_linear_base, solversNewtonRaphsonSolve_Test)
{
  using namespace pressio;
  using namespace pressio::solvers;

  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = containers::Vector<vector_n_t>;

  auto solver = NonLinearSolvers::createIterativeSolver<nonlinear::NewtonRaphson, linear::Bicgstab>();

  vector_w_t b(2);
  b[0] = 0.4;
  b[1] = 0.5;

  ValidSystem sys;
  auto y = solver.solve<linear::DefaultPreconditioner, L2Norm>(sys, b);

  EXPECT_NEAR( y[0],  1.0, 1e-8 );
  EXPECT_NEAR( y[1],  0.0, 1e-8 );
}
*/














/*
struct NonLinearSystem {

    // Matrix typedefs

    using matrix_n_t = Eigen::SparseMatrix<double>;
    using matrix_w_t = containers::Matrix<matrix_n_t>;

    // Vector typedefs
    using vector_n_t = Eigen::VectorXd;
    using vector_w_t = containers::Vector<vector_n_t>;

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
  containers::Vector<Eigen::VectorXd> b(2);

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
  using state_t = containers::Vector<native_state_t>;
  using jac_t = containers::Matrix<native_jac_t>;
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
*/
