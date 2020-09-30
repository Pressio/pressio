
#ifndef SOLVERS_TESTS_EIGEN_ROSENBROCK_N4_HPP_
#define SOLVERS_TESTS_EIGEN_ROSENBROCK_N4_HPP_

#include "pressio_solvers.hpp"

namespace pressio{ namespace solvers{ namespace test{

struct EigenRosenbrock4Impl
{
  using eig_dyn_mat	= Eigen::MatrixXd;
  using eig_dyn_vec	= Eigen::VectorXd;
  using jacobian_w_t	= pressio::containers::DenseMatrix<eig_dyn_mat>;
  using state_w_t	= pressio::containers::Vector<eig_dyn_vec>;

  using scalar_type = double;
  using state_type	= state_w_t;
  using residual_type	= state_type;
  using jacobian_type	= jacobian_w_t;

  static constexpr int nf = 6; // num functions
  static constexpr int nv = 4; // num variables

  residual_type createResidual() const{return residual_type(nf);}
  jacobian_type createJacobian() const{return jacobian_type(nf, nv);}

  // void residualNorm(const state_type & state,
  // 		    pressio::Norm normKind,
  // 		    scalar_type & resNorm) const
  // {
  //   // here I can create one R every time, because performance does not matter
  //   // but it would be better to create a R only once
  //   auto R = createResidualObject(state);
  //   residual(state, R, normKind, resNorm);
  // }

  void residual(const state_type& x, residual_type & res) const
  {
    auto x1 = x[0];
    auto x2 = x[1];
    auto x3 = x[2];
    auto x4 = x[3];
    res[0] = 10.*(x4 - x3*x3);
    res[1] = 10.*(x3 - x2*x2);
    res[2] = 10.*(x2 - x1*x1);
    res[3] = (1.-x1);
    res[4] = (1.-x2);
    res[5] = (1.-x3);
    // if (normKind == pressio::Norm::L2) normResidual = res.data()->norm();
    // if (normKind == pressio::Norm::L1) normResidual = res.data()->lpNorm<1>();
  }

  void jacobian(const state_type & x, jacobian_type & jac) const 
  {
    auto x1 = x[0];
    auto x2 = x[1];
    auto x3 = x[2];
    auto & JJ = *jac.data();
    JJ.setZero();

    JJ(0,2) = -20.*x3;
    JJ(0,3) = 10.;
    JJ(1,1) = -20.*x2;
    JJ(1,2) = 10.;
    JJ(2,0) = -20.*x1;
    JJ(2,1) = 10.;
    JJ(3,0) = -1.;
    JJ(4,1) = -1.;
    JJ(5,2) = -1.;
  }

};

using EigenRosenbrock4 = EigenRosenbrock4Impl;

struct EigenRosenbrock4HessGradApi
{
  using eig_dyn_mat	= Eigen::MatrixXd;
  using eig_dyn_vec	= Eigen::VectorXd;
  using jacobian_w_t	= pressio::containers::DenseMatrix<eig_dyn_mat>;
  using state_w_t	= pressio::containers::Vector<eig_dyn_vec>;

  using scalar_type	= double;
  using state_type	= state_w_t;
  using hessian_type	= pressio::containers::DenseMatrix<eig_dyn_mat>;
  using gradient_type	= state_type;

  static constexpr int nf = 6; // num functions
  static constexpr int nv = 4; // num variables

  EigenRosenbrock4Impl rosImpl;

public:
  hessian_type createHessian() const{
    // this only constructs empty object
    return hessian_type(nv, nv);
  }

  gradient_type createGradient() const{
    // this only constructs empty object
    return gradient_type(nv);
  }

  void residualNorm(const state_type & state,
		    pressio::Norm normKind,
		    scalar_type & normResidual) const
  {
    auto R = rosImpl.createResidual();
    rosImpl.residual(state, R);//, normKind, normResidual);
    if (normKind == pressio::Norm::L2) normResidual = R.data()->norm();
    if (normKind == pressio::Norm::L1) normResidual = R.data()->lpNorm<1>();    
  }

  void hessianAndGradient(const state_type & x,
			  hessian_type & hess,
			  gradient_type & grad,
			  pressio::Norm normType,
			  scalar_type & residualNorm,
        bool recomputeJacobian) const
  {
    auto J = rosImpl.createJacobian();
    rosImpl.jacobian(x, J);
    *hess.data() = J.data()->transpose() * (*J.data());

    auto R = rosImpl.createResidual();
    rosImpl.residual(x, R);//, normType, residualNorm);

    *grad.data() = J.data()->transpose() * (*R.data());

    if (normType == ::pressio::Norm::L2) residualNorm = R.data()->norm();
    if (normType == ::pressio::Norm::L1) residualNorm = R.data()->lpNorm<1>();
  }
};

}}} //end namespace pressio::solvers::test
#endif
