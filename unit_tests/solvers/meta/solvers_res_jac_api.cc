
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"

struct ValidSystemA {
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = pressio::containers::Matrix<matrix_n_t>;
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = pressio::containers::Vector<vector_n_t>;

  using scalar_type     = double;
  using state_type	= vector_w_t;
  using residual_type	= state_type;
  using jacobian_type	= matrix_w_t;

  residual_type createResidualObject(const state_type& x) const;
  jacobian_type createJacobianObject(const state_type& x) const;

  void residual(const state_type& x, residual_type & res,
		::pressio::solvers::Norm normKind, scalar_type & normResidual) const;

  void jacobian(const state_type& x, jacobian_type & jac) const;
};

struct NonValidSystemA1 {
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = pressio::containers::Matrix<matrix_n_t>;
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = pressio::containers::Vector<vector_n_t>;

  using scalar_type     = double;
  using state_type	= vector_w_t;
  using residual_type	= state_type;
  using jacobian_type	= matrix_w_t;

  //residual_type createResidualObject(const state_type& x) const;
  jacobian_type createJacobianObject(const state_type& x) const;

  void residual(const state_type& x, residual_type & res,
		::pressio::solvers::Norm normKind, scalar_type & normResidual) const;

  void jacobian(const state_type& x, jacobian_type & jac) const;
};

struct ValidSystemB {
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = pressio::containers::Matrix<matrix_n_t>;
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = pressio::containers::Vector<vector_n_t>;

  using scalar_type     = double;
  using state_type	= vector_w_t;
  using residual_type	= state_type;
  using jacobian_type	= matrix_w_t;

  residual_type createResidualObject(const state_type& x) const;
  jacobian_type createJacobianObject(const state_type& x) const;

  void residualNorm(const state_type & state,
		    pressio::solvers::Norm normKind,
		    scalar_type & resNorm) const;

  void residualAndJacobian(const state_type& x,
			   residual_type & res, jacobian_type & J,
			   ::pressio::solvers::Norm normKind,
			   scalar_type & normResidual) const;
};


TEST(solvers_meta, system_admissible_res_jac_api){
  using namespace pressio;
  using system_t   = ValidSystemA;
  static_assert(solvers::meta::system_meets_residual_jacobian_api<system_t>::value, "");
  static_assert(!solvers::meta::system_meets_fused_residual_jacobian_api<system_t>::value, "");
}

TEST(solvers_meta, system_non_admissible_res_jac_api){
  using namespace pressio;
  using system_t   = NonValidSystemA1;
  static_assert(!solvers::meta::system_meets_residual_jacobian_api<system_t>::value, "");
  static_assert(!solvers::meta::system_meets_fused_residual_jacobian_api<system_t>::value, "");
}

TEST(solvers_meta, system_admissible_fused_res_jac_api){
  using namespace pressio;
  using system_t   = ValidSystemB;
  static_assert(!solvers::meta::system_meets_residual_jacobian_api<system_t>::value, "");
  static_assert(solvers::meta::system_meets_fused_residual_jacobian_api<system_t>::value, "");
}


// TEST(solvers_meta, detect_residual_methods){
//   using namespace pressio;
//   using system_t   = ValidSystemA;
//   static_assert(solvers::meta::system_has_needed_residual_methods
// 		<system_t,
// 		typename system_t::state_type,
// 		typename system_t::residual_type
// 		>::value, "");

//   static_assert(solvers::meta::has_residual_method_callable_with_one_arg<
// 			system_t, typename system_t::state_type>::value,
// 		"system does not have residual with one arg");

//   static_assert(solvers::meta::has_residual_method_callable_with_one_arg<
// 			system_t, double>::value == false,
// 		"system has residual method with one arg with wrong type");

//   static_assert(solvers::meta::has_residual_method_callable_with_two_args<
// 		system_t,
// 		typename system_t::state_type,
// 		typename system_t::residual_type
// 		>::value,
// 		"system does not have residual with two args");

//   static_assert(solvers::meta::has_residual_method_callable_with_two_args<
// 		system_t,
// 		typename system_t::state_type,
// 		double
// 		>::value == false,
// 		"system has residual with two args wrong type");
// }

// TEST(solvers_meta, detect_jacobian_methods){
//   using namespace pressio;
//   using system_t   = ValidSystemA;

//   static_assert(solvers::meta::has_jacobian_method_callable_with_one_arg<
// 		system_t,
// 		typename system_t::state_type
// 		>::value,
// 		"system does not have jacobian with one arg");

//   static_assert(solvers::meta::has_jacobian_method_callable_with_two_args<
// 		system_t,
// 		typename system_t::state_type,
// 		typename system_t::jacobian_type
// 		>::value,
// 		"system does not have jacobian with two args");

//   static_assert(solvers::meta::system_has_needed_jacobian_methods
// 		<system_t,
// 		typename system_t::state_type,
// 		typename system_t::jacobian_type
// 		>::value, "");
// }
