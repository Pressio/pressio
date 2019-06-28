
#include <gtest/gtest.h>
#include "ALGEBRA_ALL"
#include "SOLVERS_NONLINEAR"
#include "QR_BASIC"

struct NonLinearLeastSquareSystem {

  using jacobian_w_t = rompp::algebra::Matrix<Eigen::MatrixXd>;
  using state_w_t = rompp::algebra::Vector<Eigen::VectorXd>;
  using state_type = state_w_t;
  using residual_type = state_type;
  using jacobian_type =  jacobian_w_t;
  static constexpr int n = 2;

  void residual(const state_type& x, residual_type & res) const {
    auto x1 = x[0];
    auto x2 = x[1];

    res[0] = 10.*(x2 - x1*x1);
    res[1] = (1.-x1);
  }

  residual_type residual(const state_type& x) const {
    residual_type res(n);
    this->residual(x, res);
    return res;
  }

  void jacobian(const state_type & x, jacobian_type & jac) const {
    auto x1 = x[0];

    (*jac.data())(0,0) = -20.*x1;
    (*jac.data())(0,1) = 10.;

    (*jac.data())(1,0) = -1.;
    (*jac.data())(1,1) = 0.0;
  }

  jacobian_type jacobian(const state_type& x) const {
    jacobian_type jac(n, 2);
    this->jacobian(x, jac);
    return jac;
  }
};


TEST(solvers_nonlinear_least_squares, rosenbrockGNQRLineSearch){
  using namespace rompp;
  using state_w_t = algebra::Vector<Eigen::VectorXd>;
  using sc_t	  = double;
  using mat_type  = typename NonLinearLeastSquareSystem::jacobian_w_t;
  NonLinearLeastSquareSystem problem;

  // define type of QR and GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  solvers::iterative::GaussNewtonQRLineSearch<sc_t, qr_type,
					      lsearch_t> solver;
  solver.setTolerance(1e-8);

  state_w_t x(2); x[0] = -1.9; x[1] = 2.;
  solver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;
  // EXPECT_NEAR( x(0), 2.4173449278229, 1e-9 );
  // EXPECT_NEAR( x(1), 0.26464986197941, 1e-9 );
}


// TEST(solvers_nonlinear_least_squares, gaussNewtonQRLineSearchDoOnly2Steps){
//   using namespace rompp;
//   using state_w_t = algebra::Vector<Eigen::VectorXd>;
//   using sc_t	  = double;
//   using mat_type  = typename NonLinearLeastSquareSystem::jacobian_w_t;
//   NonLinearLeastSquareSystem problem;

//   // define type of QR and GaussNewton solver
//   using qr_algo = qr::Householder;
//   using qr_type = qr::QRSolver<mat_type, qr_algo>;
//   using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
//   using converged_when_t
//     = solvers::iterative::converged_when::completingNumMaxIters;
//   solvers::iterative::GaussNewtonQRLineSearch<sc_t, qr_type, lsearch_t,
// 					      converged_when_t> solver;
//   // setting max iters so that in combination with the
//   // above convergence method, the solver will exit after target steps
//   solver.setMaxIterations(2);

//   state_w_t x(2); x[0] = 2.0; x[1] = 0.25;
//   solver.solve(problem, x);

//   std::cout << std::setprecision(16) << *x.data() << std::endl;
//   EXPECT_NEAR( x(0), 2.415361667771343 , 1e-8);
//   EXPECT_NEAR( x(1), 0.2648293802571118 , 1e-8);
// }


// TEST(solvers_nonlinear_least_squares, gaussNewtonQRPassTypes){
//   using namespace rompp;
//   using problem_t = NonLinearLeastSquareSystem;
//   using state_w_t = typename problem_t::state_w_t;
//   using sc_t	  = double;
//   using mat_t	  = typename problem_t::jacobian_w_t;

//   problem_t problem;
//   state_w_t x(2); x[0] = 2.0; x[1] = 0.25;

//   // define type of QR solver and GaussNewton solver
//   using qr_algo		 = qr::Householder;
//   using qr_type		 = qr::QRSolver<mat_t, qr_algo>;
//   using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
//   using converged_when_t = solvers::iterative::default_convergence;
//   using gnsolver_t	 =
//     solvers::iterative::GaussNewtonQRLineSearch<sc_t, qr_type,
// 						lsearch_t,
// 						converged_when_t,
// 						problem_t>;
//   gnsolver_t GNSolver(problem, x);
//   GNSolver.setTolerance(1e-8);
//   GNSolver.solve(problem, x);
//   std::cout << std::setprecision(14) << *x.data() << std::endl;
//   EXPECT_NEAR( x(0), 2.4173449278229, 1e-9 );
//   EXPECT_NEAR( x(1), 0.26464986197941, 1e-9 );
// }
