
#ifndef SOLVERS_TESTS_EIGEN_ROSENBROCK_N3_HPP_
#define SOLVERS_TESTS_EIGEN_ROSENBROCK_N3_HPP_

#include "pressio_solvers.hpp"

namespace pressio{ namespace solvers{ namespace test{

struct EigenRosenbrock3 
{

  using eig_dyn_mat	= Eigen::MatrixXd;
  using eig_dyn_vec	= Eigen::VectorXd;
  using jacobian_w_t	= pressio::containers::DenseMatrix<eig_dyn_mat>;
  using state_w_t	= pressio::containers::Vector<eig_dyn_vec>;

  using scalar_type = double;
  using state_type	= state_w_t;
  using residual_type	= state_type;
  using jacobian_type	= jacobian_w_t;
  static constexpr int nf = 4; // num functions
  static constexpr int nv = 3; // num variables

  residual_type createResidual() const {
    return residual_type(nf);
  }

  jacobian_type createJacobian() const {
    return jacobian_type(nf, nv);
  }

  void residual(const state_type& x,
		residual_type & res,
		::pressio::Norm normKind,
		scalar_type & normResidual) const
  {
    auto x1 = x[0];
    auto x2 = x[1];
    auto x3 = x[2];

    res[0] = 10.*(x3 - x2*x2);
    res[1] = 10.*(x2 - x1*x1);
    res[2] = (1.-x1);
    res[3] = (1.-x2);

    if (normKind == pressio::Norm::L2) normResidual = res.data()->norm();
    if (normKind == pressio::Norm::L1) normResidual = res.data()->lpNorm<1>();
  }

  void jacobian(const state_type & x, jacobian_type & jac) const {
    auto x1 = x[0];
    auto x2 = x[1];
    auto & JJ = *jac.data();
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
