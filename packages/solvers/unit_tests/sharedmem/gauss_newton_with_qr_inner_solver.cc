
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "SOLVERS_NONLINEAR"

struct NonLinearLeastSquareSystem {
  
  using jacobian_w_t = rompp::core::Matrix<Eigen::SparseMatrix<double>>;
  using state_w_t = rompp::core::Vector<Eigen::VectorXd>;
  using state_type = state_w_t;
  using jacobian_type =  jacobian_w_t;
  using vector_type =  state_type;
  using matrix_type =  jacobian_type;
  
  static constexpr int n = 8;
  const double times_[n] = {1.,2.,3.,4.,
			   5.,6.,7.,8};
  const double y_[n] = {3.29, 4.27, 5.3, 7.1,
		       10.1, 9.8, 16.1, 20.2};

  void residual(const state_w_t& x, state_w_t& res) const {
    for (auto i = 0; i < n; i++) {
      res[i] = x[0] * exp(x[1]*times_[i]) - y_[i];
    }
  }

  state_w_t residual(const state_w_t& x) const {
    state_w_t res(n);
    this->residual(x, res);
    return res;
  }

  void jacobian(const state_w_t& x, jacobian_w_t& jac) const {
    for (int i = 0; i < n; i++) {
      double expval = exp(x[1] * times_[i]);
      jac.data()->coeffRef(i, 0) = expval;
      jac.data()->coeffRef(i, 1) = x[0]*times_[i]*expval;
    }
  }

  jacobian_w_t jacobian(const state_w_t& x) const {
    jacobian_w_t jac(n, 2);
    this->jacobian(x, jac);
    return jac;
  }
};

TEST(solvers_nonlinear_least_squares, gaussNewtonQR)
{
  using namespace rompp;
  using namespace rompp::solvers;

  using state_w_t = core::Vector<Eigen::VectorXd>;
  auto solver = NonLinearSolvers::createNonLinearIterativeLeastSquareQRBasedSolver<
    nonlinearleastsquare::GaussNewtonQR>();

  state_w_t x0(2);
  x0[0] = 2.0;
  x0[1] = 0.25;
 
  NonLinearLeastSquareSystem sys;
  static_assert( details::system_traits<NonLinearLeastSquareSystem>::is_system, "");
  auto x = solver.solve(sys, x0);
  std::cout << std::setprecision(14) << *x.data() << std::endl;

  EXPECT_NEAR( x(0), 2.4173449278229, 1e-9 );
  EXPECT_NEAR( x(1), 0.26464986197941, 1e-9 );
}
