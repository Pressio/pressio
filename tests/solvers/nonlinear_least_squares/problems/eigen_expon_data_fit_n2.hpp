
#ifndef SOLVERS_TESTS_EIGEN_EXPONENTIAL_DATA_FIT_N2_HPP_
#define SOLVERS_TESTS_EIGEN_EXPONENTIAL_DATA_FIT_N2_HPP_

#include "pressio_solvers.hpp"

namespace pressio{ namespace solvers{ namespace test{

struct EigenExpDataFitN2 
{
  using scalar_type	= double;
  using state_type	= pressio::containers::Vector<Eigen::VectorXd>;
  using residual_type	= state_type;
  using jacobian_type	= pressio::containers::DenseMatrix<Eigen::MatrixXd>;

  static constexpr int n = 8;
  const double times_[n] = {1.,2.,3.,4., 5.,6.,7.,8};
  const double y_[n] = {3.29, 4.27, 5.3, 7.1, 10.1, 9.8, 16.1, 20.2};

  residual_type createResidual() const{return residual_type(n);}
  jacobian_type createJacobian() const{return jacobian_type(n,2);}

  void residual(const state_type& x, residual_type & res) const
  {
    for (auto i = 0; i < n; i++) {
      res[i] = x[0] * exp(x[1]*times_[i]) - y_[i];
    }
    // if (normKind == pressio::Norm::L2) normResidual = res.data()->norm();
    // if (normKind == pressio::Norm::L1) normResidual = res.data()->lpNorm<1>();
  }

  void jacobian(const state_type & x, jacobian_type & jac) const {
    for (int i = 0; i < n; i++) {
      double expval = exp(x[1] * times_[i]);
      (*jac.data())(i,0) = expval;
      (*jac.data())(i,1) = x[0]*times_[i]*expval;
    }
  }

};

}}} //end namespace pressio::solvers::test
#endif
