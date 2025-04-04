
#ifndef NONLINEAR_SOLVERS_TESTS_PROBLEM3_HPP_
#define NONLINEAR_SOLVERS_TESTS_PROBLEM3_HPP_

namespace pressio{ namespace solvers{ namespace test{

template<class scalar_t = double>
struct Problem3
{
  using scalar_type	= scalar_t;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::MatrixXd;

  const double times_[8] = {1.,2.,3.,4., 5.,6.,7.,8};
  const double y_[8] = {3.29, 4.27, 5.3, 7.1, 10.1, 9.8, 16.1, 20.2};


  state_type createState() const {
    state_type a(2);
    a.setZero();
    return a;
  }

  residual_type createResidual() const {
    residual_type a(8);
    a.setZero();
    return a;
  }

  jacobian_type createJacobian() const {
    jacobian_type a(8, 2);
    a.setZero();
    return a;
  }

  void residualAndJacobian(const state_type& x,
			   residual_type& res,
			   std::optional<jacobian_type*> Jin) const
  {
    auto * jac = Jin.value_or(nullptr);

    for (int i = 0; i < 8; i++){
      const double expval = exp(x(1) * times_[i]);
      res(i) = x(0) * expval - y_[i];

      if (jac){
	(*jac)(i,0) = expval;
	(*jac)(i,1) = x(0)*times_[i]*expval;
      }
    }
  }
};

}}} //end namespace pressio::solvers::test
#endif
