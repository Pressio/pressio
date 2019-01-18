
#ifndef SOLVERS_IMPL_CONVERGE_CRITERION_MIXIN_HPP
#define SOLVERS_IMPL_CONVERGE_CRITERION_MIXIN_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../CORE_OPS"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <typename conv_tag>
struct IsConvergedHelper;

template <>
struct IsConvergedHelper<converged_when::completingNumMaxIters>{
  static constexpr char const * description_ = "complete max iters";

  template <typename state_t, typename step_t, typename scalar_t>
  bool operator()(const state_t & y, const state_t & dy,
		  scalar_t norm_dy, step_t step,
		  step_t maxIters, scalar_t tol) const{
    return step==maxIters;
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::absoluteNormCorrectionBelowTol<norm_t>>{

  static constexpr char const * description_ = "norm(dx) < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  bool operator()(const state_t & y, const state_t & dy,
		  scalar_t norm_dy, step_t step,
		  step_t maxIters, scalar_t tol) const{
    return (norm_dy<tol);
  }
};


}}}} //end namespace rompp::solvers::iterative::impl
#endif
