
#include <gtest/gtest.h>

#include <iostream>

#include "CONTAINERS_ALL"
// #include "matrix/concrete/containers_matrix_sparse_serial_eigen.hpp"
// #include "vector/concrete/containers_vector_serial_eigen.hpp"

#include "experimental/solvers_l2_vector_norm.hpp"
#include "experimental/solvers_nonlinear_base.hpp"
#include "experimental/solvers_nonlinear_traits.hpp"
#include "experimental/solvers_linear_factory.hpp"
#include "experimental/solvers_nonlinear_factory.hpp"
#include "experimental/solvers_policy_nonlinear_iterative.hpp"


struct ValidSystem {

    // Matrix typedefs
    using matrix_n_t = Eigen::SparseMatrix<double>;
    using matrix_w_t = rompp::containers::Matrix<matrix_n_t>;

    // Vector typedefs
    using vector_n_t = Eigen::VectorXd;
    using vector_w_t = rompp::containers::Vector<vector_n_t>;

    typedef vector_w_t vector_type;
    typedef matrix_w_t matrix_type;


    void residual(const vector_w_t& x, vector_w_t& res) const {
      res[0] =  x[0]*x[0]*x[0] + x[1] - 1.0;
      res[1] = -x[0] + x[1]*x[1]*x[1] + 1.0;
    }


    vector_w_t residual(const vector_w_t& x) const {
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


    matrix_w_t jacobian(const vector_w_t& x) const {
      matrix_w_t jac(2, 2);
      jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
      jac.data()->coeffRef(0, 1) =  1.0;
      jac.data()->coeffRef(1, 0) = -1.0;
      jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1];
      return jac;
    }
};


TEST(solvers_nonlinear_base, solversBaseGettersTest)
{
  using namespace rompp;
  using namespace rompp::solvers;

  auto solver = NonLinearSolvers::createIterativeSolver<nonlinear::NewtonRaphson, linear::Bicgstab>();

  auto x = solver.getMaxIterations();
  auto xNL = solver.getMaxNonLinearIterations();

  auto tol = solver.getTolerance();
  auto tolNL = solver.getNonLinearTolerance();

  EXPECT_EQ(x, 100);
  EXPECT_EQ(xNL, 100);
  EXPECT_NEAR(tol, 1.0e-5, 1.0e-8);
  EXPECT_NEAR(tolNL, 1.0e-5, 1.0e-8);
}


TEST(solvers_non_linear_base, solversBaseSettersTest)
{
  using namespace rompp;
  using namespace rompp::solvers;

  auto solver = NonLinearSolvers::createIterativeSolver<nonlinear::NewtonRaphson, linear::Bicgstab>();

  solver.setMaxIterations(200);
  solver.setMaxNonLinearIterations(222);

  solver.setTolerance(-2.0e-5);
  solver.setNonLinearTolerance(-2.0e-5);

  auto x = solver.getMaxIterations();
  auto xNL = solver.getMaxNonLinearIterations();

  auto tol = solver.getTolerance();
  auto tolNL = solver.getNonLinearTolerance();

  EXPECT_EQ(x, 200);
  EXPECT_EQ(xNL, 222);
  EXPECT_NEAR(tol, 2.0e-5, 1.0e-8);
  EXPECT_NEAR(tolNL, 2.0e-5, 1.0e-8);
}


TEST(solvers_non_linear_base, solversBaseSolveTest)
{
  using namespace rompp;
  using namespace rompp::solvers;

  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = containers::Vector<vector_n_t>;

  auto solver = NonLinearSolvers::createIterativeSolver<nonlinear::NewtonRaphson, linear::Bicgstab>();

  vector_w_t b(2);
  b[0] = 0.4;
  b[1] = 0.5;

  ValidSystem sys;
  auto y = solver.solve(sys, b);

  EXPECT_NEAR( y[0],  1.0, 1e-8 );
  EXPECT_NEAR( y[1],  0.0, 1e-8 );
}


// TEST(solvers_non_linear_base, solversBaseBadSolveTest)
// {
//   using namespace rompp;
//   using namespace rompp::solvers;

//   auto solver = NonLinearSolvers::createIterativeSolver<nonlinear::NewtonRaphson, linear::Bicgstab>();

//   double left; int right;

//   ASSERT_DEATH(solver.solve(left, right), "Error: either the nonlinear system or the solution hint is invalid.");
// }


TEST(solvers_non_linear_base, solversNewtonRaphsonSolve_Test)
{
  using namespace rompp;
  using namespace rompp::solvers;

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
