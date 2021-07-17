
#ifndef SOLVERS_TESTS_EIGEN_ROSENBROCK_N3_HPP_
#define SOLVERS_TESTS_EIGEN_ROSENBROCK_N3_HPP_

#include "pressio_solvers.hpp"

namespace pressio{ namespace solvers{ namespace test{

struct EigenRosenbrock3 
{

  using scalar_type = double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::MatrixXd;

  residual_type createResidual() const{return residual_type(4);}
  jacobian_type createJacobian() const{return jacobian_type(4, 3);}

  void residual(const state_type& x, residual_type & res) const
  {
    auto x1 = x(0);
    auto x2 = x(1);
    auto x3 = x(2);

    res(0) = 10.*(x3 - x2*x2);
    res(1) = 10.*(x2 - x1*x1);
    res(2) = (1.-x1);
    res(3) = (1.-x2);
  }

  void jacobian(const state_type & x, jacobian_type & JJ) const {
    auto x1 = x(0);
    auto x2 = x(1);
    JJ.setZero();

    JJ(0,1) = -20.*x2;
    JJ(0,2) = 10.;
    JJ(1,0) = -20.*x1;
    JJ(1,1) = 10.;
    JJ(2,0) = -1.;
    JJ(3,1) = -1.;
  }

};

}}} //end namespace pressio::solvers::test
#endif
