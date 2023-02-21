
#ifndef NONLINEAR_SOLVERS_TESTS_PROBLEM8_HPP_
#define NONLINEAR_SOLVERS_TESTS_PROBLEM8_HPP_

namespace pressio{ namespace solvers{ namespace test{

template<class scalar_t = double>
struct Problem8
{
  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::MatrixXd;

  const double times_[16] =
    {
      5.00E+01,
      5.50E+01,
      6.00E+01,
      6.50E+01,
      7.00E+01,
      7.50E+01,
      8.00E+01,
      8.50E+01,
      9.00E+01,
      9.50E+01,
      1.00E+02,
      1.05E+02,
      1.10E+02,
      1.15E+02,
      1.20E+02,
      1.25E+02
    };
  const double y_[16] =
    {
      3.4780E+04,
      2.8610E+04,
      2.3650E+04,
      1.9630E+04,
      1.6370E+04,
      1.3720E+04,
      1.1540E+04,
      9.7440E+03,
      8.2610E+03,
      7.0300E+03,
      6.0050E+03,
      5.1470E+03,
      4.4270E+03,
      3.8200E+03,
      3.3070E+03,
      2.8720E+03
    };

  state_type createState() const{return state_type(3);}
  residual_type createResidual() const{return residual_type(16);}
  jacobian_type createJacobian() const{return jacobian_type(16,3);}

  void residual(const state_type& x, residual_type & res) const
  {
    for (auto i = 0; i < res.size(); i++) {
      res(i) = x(0) * exp(x(1)/(times_[i] + x(2))) - y_[i];
    }
  }

  void jacobian(const state_type & x, jacobian_type & jac) const {
    for (int i = 0; i < 16; i++) {
      double expval = exp(x(1)/(times_[i] + x(2)));
      double tplusx2 = times_[i] + x(2);
      jac(i,0) = expval;
      jac(i,1) = x(0)*expval/tplusx2;
      jac(i,2) = -1*x(0)*expval*x(1)/(tplusx2*tplusx2);
    }
  }

};

}}} //end namespace pressio::solvers::test
#endif
