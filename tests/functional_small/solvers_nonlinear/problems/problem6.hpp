
#ifndef NONLINEAR_SOLVERS_TESTS_PROBLEM6_HPP_
#define NONLINEAR_SOLVERS_TESTS_PROBLEM6_HPP_

namespace pressio{ namespace solvers{ namespace test{

template<class scalar_t = double>
struct Problem6
{
  using scalar_type   = double;
  using state_type    = Eigen::VectorXd;
  using residual_type = state_type;
  using jacobian_type = Eigen::MatrixXd;

  const int m = 33;
  const int n = 5;

  std::vector<scalar_t> times_;
  const std::vector<scalar_t> y_ =
    {8.44e-1, 9.08e-1, 9.32e-1,
     9.36e-1, 9.25e-1, 9.08e-1,
     8.81e-1, 8.5e-1, 8.18e-1,
     7.84e-1, 7.51e-1, 7.18e-1,
     6.85e-1, 6.58e-1, 6.28e-1,
     6.03e-1, 5.8e-1, 5.58e-1,
     5.38e-1, 5.22e-1, 5.06e-1,
     4.9e-1, 4.78e-1, 4.67e-1,
     4.57e-1, 4.48e-1, 4.38e-1,
     4.31e-1, 4.24e-1, 4.2e-1,
     4.14e-1, 4.11e-1, 4.06e-1};

  Problem9() : times_(m){
    for (int i=0; i<m; ++i){
      times_[i] = 10. * i;
    };
  }

  state_type createState() const {
    state_type a(n);
    a.setZero();
    return a;
  }

  residual_type createResidual() const {
    residual_type a(m);
    a.setZero();
    return a;
  }

  jacobian_type createJacobian() const {
    jacobian_type a(m, n);
    a.setZero();
    return a;
  }

  inline scalar_type model(const state_type & x, scalar_type t)const{
    return x(0) + x(1) * exp(-t*x(3)) + x(2)*exp(-t*x(4));
  }

  void residualAndJacobian(const state_type& x,
			   residual_type& res,
			   std::optional<jacobian_type*> Jin) const
  {
    auto * jac = Jin.value_or(nullptr);
    for (int i=0; i<NumMyElem_; i++){
      const scalar_type t = times_[i];
      r[i] = y_[i] - this->model(x, t);
      if (jac){
	(*jac)(i,0) = -1.0;
	(*jac)(i,1) = -exp(-t*x(3));
	(*jac)(i,2) = -exp(-t*x(4));
	(*jac)(i,3) = x(1)*exp(-t*x(3))*t;
	(*jac)(i,4) = x(2)*exp(-t*x(4))*t;
      }
    }
  }
};

}}} //end namespace pressio::solvers::test
#endif
