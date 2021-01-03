// Note: This problem is mgh10 from 
// https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
#ifndef SOLVERS_TESTS_EIGEN_NIST_MGH10_HPP_
#define SOLVERS_TESTS_EIGEN_NIST_MGH10_HPP_

#include "pressio_solvers.hpp"

namespace pressio{ namespace solvers{ namespace test{

struct EigenNISTmgh10 
{
  using scalar_type	= double;
  using state_type	= pressio::containers::Vector<Eigen::VectorXd>;
  using residual_type	= state_type;
  using jacobian_type	= pressio::containers::DenseMatrix<Eigen::MatrixXd>;

  static constexpr int n = 16;
  const double times_[n] = 
    {
      5.00E+01,
      5.50E+01,
      6.00E+01,
      6.50E+01,
      7.00E+01,
      7.50E+01,
      8.00E+01,
      8.50E+01,
      9.00E+01,
      9.50E+01,
      1.00E+02,
      1.05E+02,
      1.10E+02,
      1.15E+02,
      1.20E+02,
      1.25E+02
    };
  const double y_[n] = 
    {
      3.4780E+04,
      2.8610E+04,
      2.3650E+04,
      1.9630E+04,
      1.6370E+04,
      1.3720E+04,
      1.1540E+04,
      9.7440E+03,
      8.2610E+03,
      7.0300E+03,
      6.0050E+03,
      5.1470E+03,
      4.4270E+03,
      3.8200E+03,
      3.3070E+03,
      2.8720E+03
    };

  residual_type createResidual() const{return residual_type(n);}
  jacobian_type createJacobian() const{return jacobian_type(n,3);}

  void residual(const state_type& x, residual_type & res) const
  {
    for (auto i = 0; i < n; i++) {
      res(i) = x(0) * exp(x(1)/(times_[i] + x(2))) - y_[i];
    }
    // if (normKind == pressio::Norm::L2) normResidual = res.data()->norm();
    // if (normKind == pressio::Norm::L1) normResidual = res.data()->lpNorm<1>();
  }

  void jacobian(const state_type & x, jacobian_type & jac) const {
    for (int i = 0; i < n; i++) {
      double expval = exp(x(1)/(times_[i] + x(2)));
      double tplusx2 = times_[i] + x(2);
      (*jac.data())(i,0) = expval;
      (*jac.data())(i,1) = x(0)*expval/tplusx2;
      (*jac.data())(i,2) = -1*x(0)*expval*x(1)/(tplusx2*tplusx2);
    }
  }

};

}}} //end namespace pressio::solvers::test
#endif
