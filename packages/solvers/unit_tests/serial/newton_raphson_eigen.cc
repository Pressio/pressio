
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "SOLVERS_NONLINEAR"

struct ValidSystem {
  // Matrix typedefs
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = rompp::core::Matrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = rompp::core::Vector<vector_n_t>;

  using state_type = vector_w_t;
  using residual_type = state_type;
  using jacobian_type = matrix_w_t;
  using scalar_type = double;

  void residual(const state_type& x, residual_type& res) const {
    res[0] =  x[0]*x[0]*x[0] + x[1] - 1.0;
    res[1] = -x[0] + x[1]*x[1]*x[1] + 1.0;
  }

  residual_type residual(const state_type& x) const {
    residual_type res(2);
    res[0] =  x[0]*x[0]*x[0] + x[1] - 1.0;
    res[1] = -x[0] + x[1]*x[1]*x[1] + 1.0;
    return res;
  }

  void jacobian(const state_type& x, jacobian_type& jac) const {
    jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
    jac.data()->coeffRef(0, 1) =  1.0;
    jac.data()->coeffRef(1, 0) = -1.0;
    jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1];
  }

  jacobian_type jacobian(const state_type& x) const {
    jacobian_type jac(2, 2);
    jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
    jac.data()->coeffRef(0, 1) =  1.0;
    jac.data()->coeffRef(1, 0) = -1.0;
    jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1];
    return jac;
  }
};


TEST(solvers_nonlinear, NewtonRaphsonEigen)
{
  using namespace rompp;
  using namespace rompp::solvers;

  using problem_t  = ValidSystem;
  using state_t	   = problem_t::state_type;
  using jacobian_t = problem_t::jacobian_type;
  using scalar_t   = double;

  problem_t sys;

  // linear system
  using lin_solver_t = iterative::EigenIterative<linear::iterative::LSCG, jacobian_t>;
  // nonlinear system
  NewtonRaphson<scalar_t, lin_solver_t> solver;

  // my solution vector
  state_t y(2);
  solver.solve(sys, y);

  EXPECT_NEAR( y[0],  1.0, 1e-8 );
  EXPECT_NEAR( y[1],  0.0, 1e-8 );
}






// TEST(solvers_nonlinear_base, solversBaseGettersTest)
// {
//   using namespace rompp;
//   using namespace rompp::solvers;

//   auto solver = NonLinearSolvers::createIterativeSolver<nonlinear::NewtonRaphson, linear::Bicgstab>();

//   auto x = solver.getMaxIterations();
//   auto xNL = solver.getMaxNonLinearIterations();

//   auto tol = solver.getTolerance();
//   auto tolNL = solver.getNonLinearTolerance();

//   EXPECT_EQ(x, 100);
//   EXPECT_EQ(xNL, 100);
//   EXPECT_NEAR(tol, 1.0e-5, 1.0e-8);
//   EXPECT_NEAR(tolNL, 1.0e-5, 1.0e-8);
// }


// TEST(solvers_non_linear_base, solversBaseSettersTest)
// {
//   using namespace rompp;
//   using namespace rompp::solvers;

//   auto solver = NonLinearSolvers::createIterativeSolver<nonlinear::NewtonRaphson, linear::Bicgstab>();

//   solver.setMaxIterations(200);
//   solver.setMaxNonLinearIterations(222);

//   solver.setTolerance(-2.0e-5);
//   solver.setNonLinearTolerance(-2.0e-5);

//   auto x = solver.getMaxIterations();
//   auto xNL = solver.getMaxNonLinearIterations();

//   auto tol = solver.getTolerance();
//   auto tolNL = solver.getNonLinearTolerance();

//   EXPECT_EQ(x, 200);
//   EXPECT_EQ(xNL, 222);
//   EXPECT_NEAR(tol, 2.0e-5, 1.0e-8);
//   EXPECT_NEAR(tolNL, 2.0e-5, 1.0e-8);
// }


// TEST(solvers_non_linear_base, solversBaseSolveTest)
// {
//   using namespace rompp;
//   using namespace rompp::solvers;

//   using vector_n_t = Eigen::VectorXd;
//   using vector_w_t = core::Vector<vector_n_t>;

//   auto solver = NonLinearSolvers::createIterativeSolver<nonlinear::NewtonRaphson, linear::Bicgstab>();

//   vector_w_t b(2);
//   b[0] = 0.4;
//   b[1] = 0.5;

//   ValidSystem sys;
//   auto y = solver.solve(sys, b);

//   EXPECT_NEAR( y[0],  1.0, 1e-8 );
//   EXPECT_NEAR( y[1],  0.0, 1e-8 );
// }


// TEST(solvers_non_linear_base, solversBaseBadSolveTest)
// {
//   using namespace rompp;
//   using namespace rompp::solvers;

//   auto solver = NonLinearSolvers::createIterativeSolver<nonlinear::NewtonRaphson, linear::Bicgstab>();

//   double left; int right;

//   ASSERT_DEATH(solver.solve(left, right), "Error: either the nonlinear system or the solution hint is invalid.");
// }
