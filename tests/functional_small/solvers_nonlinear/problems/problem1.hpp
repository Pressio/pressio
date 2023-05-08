
#ifndef NONLINEAR_SOLVERS_TESTS_PROBLEM1_HPP_
#define NONLINEAR_SOLVERS_TESTS_PROBLEM1_HPP_

namespace pressio{ namespace solvers{ namespace test{

struct Problem1
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

  void residualAndJacobian(const state_type& x,
			   residual_type& res,
#ifdef PRESSIO_ENABLE_CXX17
			   std::optional<jacobian_type*> Jin) const
#else
                           jacobian_type* Jin) const
#endif
  {
    res(0) =  x(0)*x(0)*x(0) + x(1) - 1.0;
    res(1) = -x(0) + x(1)*x(1)*x(1) + 1.0;

    if (Jin){
#ifdef PRESSIO_ENABLE_CXX17
      auto & jac = *Jin.value();
#else
      auto & jac = *Jin;
#endif
      jac.coeffRef(0, 0) = 3.0*x(0)*x(0);
      jac.coeffRef(0, 1) =  1.0;
      jac.coeffRef(1, 0) = -1.0;
      jac.coeffRef(1, 1) = 3.0*x(1)*x(1);
    }
  }
};


}}} //end namespace pressio::solvers::test
#endif
