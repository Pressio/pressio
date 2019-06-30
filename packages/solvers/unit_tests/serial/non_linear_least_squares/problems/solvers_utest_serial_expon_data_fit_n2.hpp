
#if not defined SOLVERS_UTEST_SERIAL_EXPONENTIAL_DATA_FIT_N2_HPP_
#define SOLVERS_UTEST_SERIAL_EXPONENTIAL_DATA_FIT_N2_HPP_

#include "CONTAINERS_VECTOR"
#include "CONTAINERS_MATRIX"

namespace rompp{ namespace solvers{ namespace test{

struct ExpDataFitN2 {
  using scalar_type = double;
  using state_type	= rompp::containers::Vector<Eigen::VectorXd>;
  using residual_type	= state_type;
  using jacobian_type	= rompp::containers::Matrix<Eigen::MatrixXd>;

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

}}} //end namespace rompp::solvers::test
#endif
