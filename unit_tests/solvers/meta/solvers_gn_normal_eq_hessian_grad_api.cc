
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "SOLVERS_EXPERIMENTAL"

struct System {
  using matrix_n_t = Eigen::MatrixXd;
  using matrix_w_t = pressio::containers::Matrix<matrix_n_t>;
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = pressio::containers::Vector<vector_n_t>;

  using scalar_type     = double;
  using state_type	= vector_w_t;
  using hessian_type	= matrix_w_t;
  using gradient_type	= vector_w_t;

  void computeHessianAndGradient(const state_type & x,
				 hessian_type & hess,
				 gradient_type & grad,
				 const pressio::solvers::Norm & normType,
				 scalar_type & residualNorm) const;
  hessian_type createHessianObject(const state_type & x) const;
  gradient_type createGradientObject(const state_type & x) const;
};

TEST(solvers_meta, gn_normeq_hess_grad_api){
  using namespace pressio::solvers::meta::experimental;
  using sys_t   = System;

  static_assert(system_meets_gn_hessian_gradient_api<sys_t>::value, "");
}
