
#ifndef NONLINEAR_SOLVERS_TESTS_PROBLEM9_HPP_
#define NONLINEAR_SOLVERS_TESTS_PROBLEM9_HPP_

namespace pressio{ namespace solvers{ namespace test{

template<class scalar_t = double>
struct Problem9
{
  using scalar_type = double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::MatrixXd;

  const int m = 65;
  const int n = 11;
  std::vector<scalar_t> times_;
  const std::vector<scalar_t> y_ =
    {1.366, 1.191, 1.112, 1.013,
     9.91e-1, 8.85e-1, 8.31e-1, 8.47e-1,
     7.86e-1, 7.25e-1, 7.46e-1, 6.79e-1,
     6.08e-1, 6.55e-1, 6.16e-1, 6.06e-1,
     6.02e-1, 6.26e-1, 6.51e-1, 7.24e-1,
     6.49e-1, 6.49e-1, 6.94e-1, 6.44e-1,
     6.24e-1, 6.61e-1, 6.12e-1, 5.58e-1,
     5.33e-1, 4.95e-1, 5.0e-1,  4.23e-1,
     3.95e-1, 3.75e-1, 3.72e-1, 3.91e-1,
     3.96e-1, 4.05e-1, 4.28e-1, 4.29e-1,
     5.23e-1, 5.62e-1, 6.07e-1, 6.53e-1,
     6.72e-1, 7.08e-1, 6.33e-1, 6.68e-1,
     6.45e-1, 6.32e-1, 5.91e-1, 5.59e-1,
     5.97e-1, 6.25e-1, 7.39e-1, 7.1e-1,
     7.29e-1, 7.2e-1, 6.36e-1, 5.81e-1,
     4.28e-1, 2.92e-1, 1.62e-1, 9.8e-2,
     5.4e-2};

  Problem9() : times_(m){
    for (int i=0; i<m; ++i){
      times_[i] = i/10.;
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

  scalar_type model(const state_type & x, scalar_type t)const{
    auto temp1 = exp(-x(4)*t);
    auto temp2 = exp(-x(5) * (t-x(8)) * (t-x(8)) );
    auto temp3 = exp(-x(6) * (t-x(9)) * (t-x(9)) );
    auto temp4 = exp(-x(7) * (t-x(10)) * (t-x(10)) );
    return x(0)*temp1 + x(1)*temp2 + x(2)*temp3 + x(3)*temp4;
  }

  void residualAndJacobian(const state_type& x,
			   residual_type& r,
			   std::optional<jacobian_type*> Jin) const
  {
    auto * jac = Jin.value_or(nullptr);

    for (int i=0; i<m; i++){
      const scalar_type t = times_[i];
      const auto temp1 = exp(-x(4)*t);
      const auto temp2 = exp(-x(5) * (t-x(8)) * (t-x(8)) );
      const auto temp3 = exp(-x(6) * (t-x(9)) * (t-x(9)) );
      const auto temp4 = exp(-x(7) * (t-x(10)) * (t-x(10)) );

      r[i] = y_[i] - this->model(x, t);

      if (jac){
	(*jac)(i,0)  = -temp1;
	(*jac)(i,1)  = -temp2;
	(*jac)(i,2)  = -temp3;
	(*jac)(i,3)  = -temp4;
	(*jac)(i,4)  = x(0)*temp1*t;
	(*jac)(i,5)  = x(1) * (t-x(8))  * (t-x(8))  * temp2;
	(*jac)(i,6)  = x(2) * (t-x(9))  * (t-x(9))  * temp3;
	(*jac)(i,7)  = x(3) * (t-x(10)) * (t-x(10)) * temp4;
	(*jac)(i,8)  = -2.*x(1)*x(5)*(t-x(8))*temp2;
	(*jac)(i,9)  = -2.*x(2)*x(6)*(t-x(9))*temp3;
	(*jac)(i,10) = -2.*x(3)*x(7)*(t-x(10))*temp4;
      }
    }
  }
};

}}} //end namespace pressio::solvers::test
#endif
