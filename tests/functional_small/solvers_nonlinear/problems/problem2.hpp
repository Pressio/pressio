
#ifndef NONLINEAR_SOLVERS_TESTS_PROBLEM2_HPP_
#define NONLINEAR_SOLVERS_TESTS_PROBLEM2_HPP_

#include "pressio/solvers.hpp"

namespace pressio{ namespace solvers{ namespace test{

struct Problem2
{
  using state_type = Eigen::VectorXd;
  using residual_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

  state_type createState() const {
    state_type a(2);
    a.setZero();
    return a;
  }

  residual_type createResidual() const {
    residual_type a(2);
    a.setZero();
    return a;
  }

  jacobian_type createJacobian() const {
    jacobian_type a(2, 2);
    a.setZero();
    return a;
  }

  void residual(const state_type& x,
                residual_type& res) const
  {
    res(0) = x(0)*x(0) + x(1)*x(1) - 4.0;
    res(1) = x(0)*x(1) - 1.0;
  }

  void jacobian(const state_type& x, jacobian_type& jac) const {
    jac.coeffRef(0, 0) = 2*x(0);
    jac.coeffRef(0, 1) = 2*x(1);
    // Have incorrect entries so that line search is required
    jac.coeffRef(1, 0) = 0.1*x(1);
    jac.coeffRef(1, 1) = 0.1*x(0);
  }
};

}}} //end namespace pressio::solvers::test
#endif


// struct CircleHyperbolaIntersectionWithBadJacobianSystem4HessGradApi
// {
//   using scalar_type	= double;
//   using state_type	= Eigen::VectorXd;
//   using hessian_type	= Eigen::MatrixXd;
//   using gradient_type	= state_type;
//   using residual_norm_type = double;

//   CircleHyperbolaIntersectionWithBadJacobianSystem mySystem;

// public:
//   state_type createState() const{return state_type(2);}

//   hessian_type createHessian() const{
//     return hessian_type(2, 2);
//   }

//   gradient_type createGradient() const{
//     return gradient_type(2);
//   }

//   void residualNorm(const state_type & state,
// 		    pressio::Norm normKind,
// 		    residual_norm_type & normResidual) const
//   {
//     auto R = mySystem.createResidual();
//     mySystem.residual(state, R);//, normKind, normResidual);
//     if (normKind == pressio::Norm::L2) normResidual = R.norm();
//     if (normKind == pressio::Norm::L1) normResidual = R.lpNorm<1>();
//     std::cout << R[0] << " " << R[1] << " : " << normResidual << "\n";
//   }

//   void hessianAndGradient(const state_type & x,
// 			  hessian_type & hess,
// 			  gradient_type & grad,
// 			  pressio::Norm normType,
// 			  residual_norm_type & residualNorm,
//         bool /*recomputeJacobian*/) const
//   {
//     auto J = mySystem.createJacobian();
//     mySystem.jacobian(x, J);
//     hess = J.transpose() * J;

//     auto R = mySystem.createResidual();
//     mySystem.residual(x, R);

//     grad = J.transpose() * R;

//     if (normType == ::pressio::Norm::L2) residualNorm = R.norm();
//     if (normType == ::pressio::Norm::L1) residualNorm = R.lpNorm<1>();
//   }
// };
