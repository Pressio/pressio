
#ifndef NONLINEAR_SOLVERS_TESTS_PROBLEM1_HPP_
#define NONLINEAR_SOLVERS_TESTS_PROBLEM1_HPP_

namespace pressio{ namespace solvers{ namespace test{

struct Problem1A
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
    jacobian_type a(2,2);
    a.setZero();
    return a;
  }

  void residual(const state_type& x,
                residual_type& res) const {
    res(0) =  x(0)*x(0)*x(0) + x(1) - 1.0;
    res(1) = -x(0) + x(1)*x(1)*x(1) + 1.0;
  }

  void jacobian(const state_type& x, jacobian_type& jac) const {
    jac.coeffRef(0, 0) = 3.0*x(0)*x(0);
    jac.coeffRef(0, 1) =  1.0;
    jac.coeffRef(1, 0) = -1.0;
    jac.coeffRef(1, 1) = 3.0*x(1)*x(1);
  }
};

struct Problem1B
{
  using state_type = Eigen::VectorXd;
  using residual_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

  state_type createState() const {
    return state_type(2);
  }

  residual_type createResidual() const {
    return residual_type(2);
  }

  jacobian_type createJacobian() const {
    return jacobian_type(2, 2);
  }

  void residual(const state_type& x, residual_type& res) const{
    res(0) =  x(0)*x(0)*x(0) + x(1) - 1.0;
    res(1) = -x(0) + x(1)*x(1)*x(1) + 1.0;
  }

  void residualAndJacobian(const state_type& x,
         residual_type& res,
         jacobian_type& jac) const
  {
    residual(x, res);
    jac.coeffRef(0, 0) = 3.0*x(0)*x(0);
    jac.coeffRef(0, 1) =  1.0;
    jac.coeffRef(1, 0) = -1.0;
    jac.coeffRef(1, 1) = 3.0*x(1)*x(1);
  }
};

}}} //end namespace pressio::solvers::test
#endif
