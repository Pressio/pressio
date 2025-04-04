
#ifndef NONLINEAR_SOLVERS_TESTS_PROBLEM5_HPP_
#define NONLINEAR_SOLVERS_TESTS_PROBLEM5_HPP_

namespace pressio{ namespace solvers{ namespace test{

template<class scalar_t = double>
struct Problem5a
{
  using scalar_type = double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::MatrixXd;

  state_type createState() const {
    state_type a(4);
    a.setZero();
    return a;
  }

  residual_type createResidual() const {
    residual_type a(6);
    a.setZero();
    return a;
  }

  jacobian_type createJacobian() const {
    jacobian_type a(6, 4);
    a.setZero();
    return a;
  }

  void residualAndJacobian(const state_type& x,
			   residual_type& res,
			   std::optional<jacobian_type*> Jin) const
  {
    auto x1 = x(0);
    auto x2 = x(1);
    auto x3 = x(2);
    auto x4 = x(3);

    res(0) = 10.*(x4 - x3*x3);
    res(1) = 10.*(x3 - x2*x2);
    res(2) = 10.*(x2 - x1*x1);
    res(3) = (1.-x1);
    res(4) = (1.-x2);
    res(5) = (1.-x3);

    if (Jin){
      auto & J = *Jin.value();
      J.setZero();
      J(0,2) = -20.*x3;
      J(0,3) = 10.;
      J(1,1) = -20.*x2;
      J(1,2) = 10.;
      J(2,0) = -20.*x1;
      J(2,1) = 10.;
      J(3,0) = -1.;
      J(4,1) = -1.;
      J(5,2) = -1.;
    }
  }
};

}}} //end namespace pressio::solvers::test
#endif


// struct EigenRosenbrock4HessGradApi
// {
//   using scalar_type	= double;
//   using state_type	= Eigen::VectorXd;
//   using hessian_type	= Eigen::MatrixXd;
//   using gradient_type	= state_type;
//   using residual_norm_type = double;

//   EigenRosenbrock4Impl rosImpl;

// public:
//   state_type createState() const{return state_type(4);}

//   hessian_type createHessian() const{
//     return hessian_type(4, 4);
//   }

//   gradient_type createGradient() const{
//     return gradient_type(4);
//   }

//   void residualNorm(const state_type & state,
// 		    pressio::Norm normKind,
// 		    residual_norm_type & normResidual) const
//   {
//     auto R = rosImpl.createResidual();
//     rosImpl.residual(state, R);//, normKind, normResidual);
//     if (normKind == pressio::Norm::L2) normResidual = R.norm();
//     if (normKind == pressio::Norm::L1) normResidual = R.lpNorm<1>();
//   }

//   void gradient(const state_type & x,
// 		gradient_type & grad,
// 		pressio::Norm normType,
// 		residual_norm_type & residualNorm,
// 		bool recomputeJacobian) const
//   {
//     assert(recomputeJacobian);
//     auto J = rosImpl.createJacobian();
//     rosImpl.jacobian(x, J);

//     auto R = rosImpl.createResidual();
//     rosImpl.residual(x, R);
//     grad = J.transpose() * R;

//     if (normType == ::pressio::Norm::L2) residualNorm = R.norm();
//     if (normType == ::pressio::Norm::L1) residualNorm = R.lpNorm<1>();
//   }

//   void hessian(const state_type & x,
// 	       hessian_type & hess) const
//   {
//     auto J = rosImpl.createJacobian();
//     rosImpl.jacobian(x, J);
//     hess = J.transpose() * J;
//   }
// };

// struct EigenRosenbrock4HessGradApiFused
// {
//   using scalar_type	= double;
//   using state_type	= Eigen::VectorXd;
//   using hessian_type	= Eigen::MatrixXd;
//   using gradient_type	= state_type;
//   using residual_norm_type = double;

//   EigenRosenbrock4Impl rosImpl;

// public:
//   state_type createState() const{return state_type(4);}

//   hessian_type createHessian() const{
//     return hessian_type(4, 4);
//   }

//   gradient_type createGradient() const{
//     return gradient_type(4);
//   }

//   void residualNorm(const state_type & state,
// 		    pressio::Norm normKind,
// 		    residual_norm_type & normResidual) const
//   {
//     auto R = rosImpl.createResidual();
//     rosImpl.residual(state, R);//, normKind, normResidual);
//     if (normKind == pressio::Norm::L2) normResidual = R.norm();
//     if (normKind == pressio::Norm::L1) normResidual = R.lpNorm<1>();
//   }

//   void hessianAndGradient(const state_type & x,
// 			  hessian_type & hess,
// 			  gradient_type & grad,
// 			  pressio::Norm normType,
// 			  residual_norm_type & residualNorm,
// 			  bool /*recomputeJacobian*/) const
//   {
//     auto J = rosImpl.createJacobian();
//     rosImpl.jacobian(x, J);
//     hess = J.transpose() * J;

//     auto R = rosImpl.createResidual();
//     rosImpl.residual(x, R);

//     grad = J.transpose() * R;

//     if (normType == ::pressio::Norm::L2) residualNorm = R.norm();
//     if (normType == ::pressio::Norm::L1) residualNorm = R.lpNorm<1>();
//   }
// };
