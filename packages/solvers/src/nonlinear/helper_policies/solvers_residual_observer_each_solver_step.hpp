
#ifndef SOLVERS_IMPL_RESIDUAL_OBSERVER_EACH_SOLVER_STEP_HPP
#define SOLVERS_IMPL_RESIDUAL_OBSERVER_EACH_SOLVER_STEP_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../CORE_OPS"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <typename observer_t, typename residual_t, typename = void>
struct has_method_observe_residual_each_step
  : std::false_type{};

template <typename observer_t, typename residual_t>
struct has_method_observe_residual_each_step<
  observer_t, residual_t,
  core::meta::void_t<
    decltype
    (
     std::declval<observer_t>().observeResidualEachStep
     (
      std::declval<const residual_t &>(),
      std::declval<int>()
      )
     )
    >
  > : std::true_type{};
// ---------------------------------------------------------------


template<typename observer_t,
	 typename residual_t,
	 typename = void>
struct ResidualObserverEachSolverStep{

  void operator()(const observer_t * obs,
		  const residual_t & resid,
		  int step) const{
    // no op
  }
};


template<typename observer_t,
	 typename residual_t>
struct ResidualObserverEachSolverStep<
  observer_t, residual_t,
  core::meta::enable_if_t<
    has_method_observe_residual_each_step<
      observer_t, residual_t
     >::value
    >
  >{

  void operator()(const observer_t * obs,
		  const residual_t & resid,
		  int gn_step) const{
    obs->observeResidualEachStep(gn_step, resid);
  }
};


}}}} //end namespace rompp::solvers::iterative::impl
#endif
