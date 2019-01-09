
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "SOLVERS_NONLINEAR"

struct NonLinearLeastSquareSystem {
  //using jacobian_w_t = rompp::core::Matrix<Eigen::SparseMatrix<double>>;
  using jacobian_w_t = rompp::core::Matrix<Eigen::MatrixXd>;
  using state_w_t = rompp::core::Vector<Eigen::VectorXd>;

  using state_type = state_w_t;
  using residual_type = state_type;
  using jacobian_type =  jacobian_w_t;

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
      (*jac.data())(i,0) = expval; //->coeffRef(i, 0) = expval;
      (*jac.data())(i,1) = x[0]*times_[i]*expval;

      //jac.data()->coeffRef(i,0) = expval;
      //jac.data()->coeffRef(i, 1) = x[0]*times_[i]*expval;
    }
  }

  jacobian_type jacobian(const state_type& x) const {
    jacobian_type jac(n, 2);
    this->jacobian(x, jac);
    return jac;
  }
};


TEST(solvers_nonlinear_least_squares, gaussNewtonQR)
{
  using namespace rompp;

  using state_w_t = core::Vector<Eigen::VectorXd>;
  using sc_t	  = double;
  using mat_type  = typename NonLinearLeastSquareSystem::jacobian_w_t;

  // my system object
  NonLinearLeastSquareSystem problem;

  // define type of QR and GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  rompp::solvers::GaussNewtonQR<sc_t, qr_type> solver;
  solver.setTolerance(1e-8);

  // initial state
  state_w_t x(2);
  x[0] = 2.0;
  x[1] = 0.25;

  // solve
  solver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-9 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-9 );
}



TEST(solvers_nonlinear_least_squares, gaussNewtonQRPassTypes)
{
  using namespace rompp;

  using problem_t = NonLinearLeastSquareSystem;
  using state_w_t = typename problem_t::state_w_t;
  using sc_t	  = double;
  using mat_t	  = typename problem_t::jacobian_w_t;

  // my system object
  problem_t problem;

  // initial state
  state_w_t x(2);
  x[0] = 2.0;
  x[1] = 0.25;

  // define type of QR solver and GaussNewton solver
  using qr_algo	   = qr::Householder;
  using qr_type	   = qr::QRSolver<mat_t, qr_algo>;
  using gnsolver_t = solvers::GaussNewtonQR<sc_t, qr_type, problem_t>;
  gnsolver_t GNSolver(problem, x);
  GNSolver.setTolerance(1e-8);
  // solve
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 2.4173449278229, 1e-9 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-9 );
}
