
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"

struct ValidSystemA {
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = pressio::containers::SparseMatrix<matrix_n_t>;
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = pressio::containers::Vector<vector_n_t>;

  using scalar_type     = double;
  using state_type	= vector_w_t;
  using residual_type	= state_type;
  using jacobian_type	= matrix_w_t;

  residual_type createResidual() const;
  jacobian_type createJacobian() const;

  void residual(const state_type& x, residual_type & res) const;
  void jacobian(const state_type& x, jacobian_type & jac) const;
};

struct NonValidSystemA1 {
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = pressio::containers::SparseMatrix<matrix_n_t>;
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = pressio::containers::Vector<vector_n_t>;

  using scalar_type     = double;
  using state_type	= vector_w_t;
  using residual_type	= state_type;
  using jacobian_type	= matrix_w_t;

  //residual_type createResidual() const; // missing on purpose to fail assert
  jacobian_type createJacobian() const;
  void residual(const state_type& x, residual_type & res) const;
  void jacobian(const state_type& x, jacobian_type & jac) const;
};

struct ValidSystemB {
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = pressio::containers::SparseMatrix<matrix_n_t>;
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = pressio::containers::Vector<vector_n_t>;

  using scalar_type     = double;
  using state_type	= vector_w_t;
  using residual_type	= state_type;
  using jacobian_type	= matrix_w_t;

  residual_type createResidual() const;
  jacobian_type createJacobian() const;

  void residualAndJacobian(const state_type& x,
			   residual_type & res, 
         jacobian_type & J,
			   bool recomputeJacobian) const;
};


TEST(solvers_meta, system_admissible_res_jac_api){
  using namespace pressio;
  using system_t   = ValidSystemA;
  static_assert(solvers::constraints::system_residual_jacobian<system_t>::value, "");
  static_assert(!solvers::constraints::system_fused_residual_jacobian<system_t>::value, "");
}

TEST(solvers_meta, system_non_admissible_res_jac_api){
  using namespace pressio;
  using system_t   = NonValidSystemA1;
  static_assert(!solvers::constraints::system_residual_jacobian<system_t>::value, "");
  static_assert(!solvers::constraints::system_fused_residual_jacobian<system_t>::value, "");
}

TEST(solvers_meta, system_admissible_fused_res_jac_api){
  using namespace pressio;
  using system_t   = ValidSystemB;
  static_assert(!solvers::constraints::system_residual_jacobian<system_t>::value, "");
  static_assert(solvers::constraints::system_fused_residual_jacobian<system_t>::value, "");
}

