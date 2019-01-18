
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "SOLVERS_NONLINEAR"
#include "QR_BASIC"

struct NonLinearLeastSquareSystem {
  using state_type = rompp::core::Vector<Eigen::VectorXd>;
  using residual_type = state_type;
  using jacobian_type = rompp::core::Matrix<Eigen::MatrixXd>;
  static constexpr int n = 8;
  const double times_[n] = {1.,2.,3.,4.,
			   5.,6.,7.,8};
  const double y_[n] = {3.29, 4.27, 5.3, 7.1,
		       10.1, 9.8, 16.1, 20.2};

  void residual(const state_type& x, residual_type & res) const {
    for (auto i = 0; i < n; i++) {
      res[i] = x[0] * exp(x[1]*times_[i]) - y_[i];
    }
  }

  residual_type residual(const state_type& x) const {
    residual_type res(n);
    this->residual(x, res);
    return res;
  }

  void jacobian(const state_type & x, jacobian_type & jac) const {
    for (int i = 0; i < n; i++) {
      double expval = exp(x[1] * times_[i]);
      (*jac.data())(i,0) = expval;
      (*jac.data())(i,1) = x[0]*times_[i]*expval;
    }
  }

  jacobian_type jacobian(const state_type& x) const {
    jacobian_type jac(n, 2);
    this->jacobian(x, jac);
    return jac;
  }
};


TEST(solvers_nonlinear_least_squares, gaussNewtonNormalEq){
  using namespace rompp;

  using problem_t   = NonLinearLeastSquareSystem;
  using state_t	    = typename problem_t::state_type;
  using sc_t	    = double;

  problem_t problem;
  state_t x(2); x[0] = 2.0; x[1] = 0.25;

  // define linear solver type and GaussNewton solver
  using solver_tag = solvers::linear::LSCG;
  solvers::iterative::GaussNewton<sc_t, solver_tag,
				  solvers::EigenIterative> GNSolver;
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-7 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-7 );
}


TEST(solvers_nonlinear_least_squares, gaussNewtonNormalEqPassSystemType){
  using namespace rompp;

  using problem_t = NonLinearLeastSquareSystem;
  using vec_t	  = typename problem_t::state_type;
  using mat_t	  = typename problem_t::jacobian_type;

  using state_t	  = vec_t;
  using hessian_t = mat_t;
  using sc_t	  = double;

  problem_t problem;
  state_t x(2); x[0] = 2.0; x[1] = 0.25;

  // define linear solver type and GaussNewton solver
  using solver_tag = solvers::linear::LSCG;
  using converged_when_t = solvers::iterative::default_convergence;
  using gn_t = solvers::iterative::GaussNewton<sc_t, solver_tag,
					       solvers::EigenIterative,
					       converged_when_t, problem_t,
					       hessian_t>;
  gn_t GNSolver(problem, x);
  GNSolver.setTolerance(1e-8);
  GNSolver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-7 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-7 );
}
