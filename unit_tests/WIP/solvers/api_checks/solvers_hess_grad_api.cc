
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"

struct ValidSystemA {
  using scalar_type     = double;
  using state_type	= Eigen::VectorXd;
  using hessian_type	= Eigen::SparseMatrix<double>;
  using gradient_type	= Eigen::VectorXd;

  hessian_type createHessian() const;
  gradient_type createGradient() const;

  void residualNorm(const state_type & state,
  		    pressio::Norm normKind,
  		    scalar_type & resNorm) const;

  void gradient(const state_type &,
		gradient_type &,
		::pressio::Norm normKind,
		scalar_type & normResidual,
		bool recomputeJacobian) const;

  void hessian(const state_type &, hessian_type &) const;
};

struct ValidSystemB {
  using scalar_type     = double;
  using state_type	= Eigen::VectorXd;
  using hessian_type	= Eigen::SparseMatrix<double>;
  using gradient_type	= Eigen::VectorXd;

  hessian_type createHessian() const;
  gradient_type createGradient() const;

  void residualNorm(const state_type & state,
  		    pressio::Norm normKind,
  		    scalar_type & resNorm) const;

  void hessianAndGradient(const state_type &, hessian_type &,
			  gradient_type &,
			  ::pressio::Norm normKind,
			  scalar_type & normResidual,
			  bool recomputeJacobian) const;
};

TEST(solvers_meta, system_admissible_hes_gra_api){
  using namespace pressio;
  using system_t   = ValidSystemA;
  static_assert(nonlinearsolvers::constraints::system_hessian_gradient<system_t>::value, "");
  static_assert(!nonlinearsolvers::constraints::system_fused_hessian_gradient<system_t>::value, "");
}


TEST(solvers_meta, system_admissible_fused_hes_gra_api){
  using namespace pressio;
  using system_t   = ValidSystemB;
  static_assert(!nonlinearsolvers::constraints::system_hessian_gradient<system_t>::value, "");
  static_assert(nonlinearsolvers::constraints::system_fused_hessian_gradient<system_t>::value, "");
}
