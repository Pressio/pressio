
#if not defined SOLVERS_UTEST_SERIAL_ROSENBROCK_N3_HPP_
#define SOLVERS_UTEST_SERIAL_ROSENBROCK_N3_HPP_

#include "CONTAINERS_VECTOR"
#include "CONTAINERS_MATRIX"

namespace pressio{ namespace solvers{ namespace test{

struct Rosenbrock3 {

  using eig_dyn_mat	= Eigen::MatrixXd;
  using eig_dyn_vec	= Eigen::VectorXd;
  using jacobian_w_t	= pressio::containers::Matrix<eig_dyn_mat>;
  using state_w_t	= pressio::containers::Vector<eig_dyn_vec>;

  using scalar_type = double;
  using state_type	= state_w_t;
  using residual_type	= state_type;
  using jacobian_type	= jacobian_w_t;
  static constexpr int nf = 4; // num functions
  static constexpr int nv = 3; // num variables

  void residual(const state_type& x, residual_type & res) const {
    auto x1 = x[0];
    auto x2 = x[1];
    auto x3 = x[2];

    res[0] = 10.*(x3 - x2*x2);
    res[1] = 10.*(x2 - x1*x1);
    res[2] = (1.-x1);
    res[3] = (1.-x2);
  }

  residual_type residual(const state_type& x) const {
    residual_type res(nf);
    this->residual(x, res);
    return res;
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

  jacobian_type jacobian(const state_type& x) const {
    jacobian_type jac(nf, nv);
    this->jacobian(x, jac);
    return jac;
  }
};

}}} //end namespace pressio::solvers::test
#endif
