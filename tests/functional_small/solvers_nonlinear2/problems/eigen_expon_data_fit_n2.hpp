
#ifndef SOLVERS_TESTS_EIGEN_EXPONENTIAL_DATA_FIT_N2_HPP_
#define SOLVERS_TESTS_EIGEN_EXPONENTIAL_DATA_FIT_N2_HPP_

#include "pressio/solvers.hpp"

namespace pressio{ namespace solvers{ namespace test{

template<class scalar_t = double>
struct EigenExpDataFitN2 
{
  using scalar_type	= scalar_t;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::MatrixXd;

  const double times_[8] = {1.,2.,3.,4., 5.,6.,7.,8};
  const double y_[8] = {3.29, 4.27, 5.3, 7.1, 10.1, 9.8, 16.1, 20.2};

  state_type createState() const{return state_type(2);}
  residual_type createResidual() const{return residual_type(8);}
  jacobian_type createJacobian() const{return jacobian_type(8,2);}

  void residual(const state_type& x, residual_type & res) const
  {
    for (auto i = 0; i < res.size(); i++) {
      res(i) = x(0) * exp(x(1)*times_[i]) - y_[i];
    }
  }

  void jacobian(const state_type & x, jacobian_type & jac) const {
    for (int i = 0; i < 8; i++) {
      double expval = exp(x(1) * times_[i]);
      jac(i,0) = expval;
      jac(i,1) = x(0)*times_[i]*expval;
    }
  }

};

}}} //end namespace pressio::solvers::test
#endif
