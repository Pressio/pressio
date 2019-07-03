
#ifndef SOLVERS_IMPL_CONVERGE_CRITERION_POLICY_HPP
#define SOLVERS_IMPL_CONVERGE_CRITERION_POLICY_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../CONTAINERS_OPS"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <typename conv_tag>
struct IsConvergedHelper;

template <>
struct IsConvergedHelper<converged_when::completingNumMaxIters>{
  static constexpr char const * description_ = "complete max iters";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return step==maxIters;
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::absoluteNormCorrectionBelowTol<norm_t>>{

  static constexpr char const * description_ = "||dy|| < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return (norm_dy<tol);
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::absoluteNormResidualBelowTol<norm_t>>{

  static constexpr char const * description_ = "||R|| < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return (norm_r<tol);
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::relativeNormResidualBelowTol<norm_t>>{

  static constexpr char const * description_ = "||R||(r) < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return (norm_r/norm_r0<tol);
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::absoluteNormProjectedResidualBelowTol<norm_t>>{

  static constexpr char const * description_ = "||P^T R|| < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return (norm_proj_r<tol);
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::relativeNormProjectedResidualBelowTol<norm_t>>{

  static constexpr char const * description_ = "||P^T R||(r) < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return (norm_proj_r/norm_proj_r0<tol);
  }
};



}}}} //end namespace pressio::solvers::iterative::impl
#endif
