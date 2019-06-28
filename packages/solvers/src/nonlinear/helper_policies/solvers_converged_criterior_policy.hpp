
#ifndef SOLVERS_IMPL_CONVERGE_CRITERION_POLICY_HPP
#define SOLVERS_IMPL_CONVERGE_CRITERION_POLICY_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../ALGEBRA_OPS"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <typename conv_tag>
struct IsConvergedHelper;

template <>
struct IsConvergedHelper<converged_when::completingNumMaxIters>{
  static constexpr char const * description_ = "complete max iters";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y, const state_t & dy,
		       scalar_t norm_dy, scalar_t norm_r,
		       scalar_t norm_r0, step_t step,
		       step_t maxIters, scalar_t tol) {
    return step==maxIters;
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::absoluteNormCorrectionBelowTol<norm_t>>{

  static constexpr char const * description_ = "||dy|| < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y, const state_t & dy,
		       scalar_t norm_dy, scalar_t norm_r,
		       scalar_t norm_r0, step_t step,
		       step_t maxIters, scalar_t tol) {
    return (norm_dy<tol);
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::absoluteNormResidualBelowTol<norm_t>>{

  static constexpr char const * description_ = "||R|| < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y, const state_t & dy,
		       scalar_t norm_dy, scalar_t norm_r,
		       scalar_t norm_r0, step_t step,
		       step_t maxIters, scalar_t tol) {
    return (norm_r<tol);
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::relativeNormResidualBelowTol<norm_t>>{

  static constexpr char const * description_ = "||R||(r) < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y, const state_t & dy,
		       scalar_t norm_dy, scalar_t norm_r,
		       scalar_t norm_r0, step_t step,
		       step_t maxIters, scalar_t tol) {
    return (norm_r/norm_r0<tol);
  }
};


}}}} //end namespace rompp::solvers::iterative::impl
#endif
