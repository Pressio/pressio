
#if not defined SOLVERS_UTEST_SERIAL_ROSENBROCK_N4_HPP_
#define SOLVERS_UTEST_SERIAL_ROSENBROCK_N4_HPP_

#include "pressio_solvers.hpp"

namespace pressio{ namespace solvers{ namespace test{

struct Rosenbrock4Impl{
  using eig_dyn_mat	= Eigen::MatrixXd;
  using eig_dyn_vec	= Eigen::VectorXd;
  using jacobian_w_t	= pressio::containers::Matrix<eig_dyn_mat>;
  using state_w_t	= pressio::containers::Vector<eig_dyn_vec>;

  using scalar_type = double;
  using state_type	= state_w_t;
  using residual_type	= state_type;
  using jacobian_type	= jacobian_w_t;

  static constexpr int nf = 6; // num functions
  static constexpr int nv = 4; // num variables

  void residual(const state_type& x, residual_type & res) const {
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
  }

  residual_type residual(const state_type& x) const {
    residual_type res(nf);
    this->residual(x, res);
    return res;
  }

  void jacobian(const state_type & x, jacobian_type & jac) const {
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

  jacobian_type jacobian(const state_type& x) const {
    jacobian_type jac(nf, nv);
    this->jacobian(x, jac);
    return jac;
  }
};


using Rosenbrock4 = Rosenbrock4Impl;



struct Rosenbrock4HessGradApi{
  using eig_dyn_mat	= Eigen::MatrixXd;
  using eig_dyn_vec	= Eigen::VectorXd;
  using jacobian_w_t	= pressio::containers::Matrix<eig_dyn_mat>;
  using state_w_t	= pressio::containers::Vector<eig_dyn_vec>;

  using scalar_type = double;
  using state_type	= state_w_t;
  using hessian_type	= pressio::containers::Matrix<eig_dyn_mat>;
  using gradient_type	= state_type;

  static constexpr int nf = 6; // num functions
  static constexpr int nv = 4; // num variables

  Rosenbrock4Impl rosImpl;

  void computeHessianAndGradient(const state_type & x,
  				 hessian_type & hess,
  				 gradient_type & grad,
  				 const pressio::solvers::Norm & normType,
  				 scalar_type & residualNorm) const{
    auto J = rosImpl.jacobian(x);
    *hess.data() = J.data()->transpose() * (*J.data());
    const auto R = rosImpl.residual(x);
    *grad.data() = J.data()->transpose() * (*R.data());
    if (normType == ::pressio::solvers::Norm::L2)
      residualNorm = R.data()->norm();
    else
      throw std::runtime_error("RosenbrockN4 only supports L2 norm");
  }

  hessian_type createHessianObject(const state_type & x) const{
    // this only construct empty objects
    hessian_type hess(nv, nv);
    return hess;
  }

  gradient_type createGradientObject(const state_type & x) const{
    // this only construct empty objects
    gradient_type grad(nv);
    return grad;
  }
};

}}} //end namespace pressio::solvers::test
#endif
