
#if not defined SOLVERS_UTEST_SERIAL_ROSENBROCK_N4_HPP_
#define SOLVERS_UTEST_SERIAL_ROSENBROCK_N4_HPP_

#include "CORE_VECTOR"
#include "CORE_MATRIX"

namespace rompp{ namespace solvers{ namespace test{

struct Rosenbrock4 {

  using eig_dyn_mat	= Eigen::MatrixXd;
  using eig_dyn_vec	= Eigen::VectorXd;
  using jacobian_w_t	= rompp::core::Matrix<eig_dyn_mat>;
  using state_w_t	= rompp::core::Vector<eig_dyn_vec>;

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

}}} //end namespace rompp::solvers::test
#endif
