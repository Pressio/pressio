
#ifndef SOLVERS_IMPL_RESIDUAL_OBSERVER_WHEN_SOLVER_CONVERGED_HPP
#define SOLVERS_IMPL_RESIDUAL_OBSERVER_WHEN_SOLVER_CONVERGED_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../CORE_OPS"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <typename observer_t, typename residual_t, typename = void>
struct has_method_observe_residual_when_solver_converged
  : std::false_type{};

template <typename observer_t, typename residual_t>
struct has_method_observe_residual_when_solver_converged<
  observer_t, residual_t,
  core::meta::void_t<
    decltype
    (
     std::declval<observer_t>().observeResidualWhenSolverConverged
     (
      std::declval<const residual_t &>()
      )
     )
    >
  > : std::true_type{};
// ---------------------------------------------------------------


template<typename observer_t,
	 typename residual_t,
	 typename = void>
struct ResidualObserverWhenSolverConverged{

  void operator()(const observer_t * obs,
		  const residual_t & resid) const{
    // no op
  }
};


template<typename observer_t,
	 typename residual_t>
struct ResidualObserverWhenSolverConverged<
  observer_t, residual_t,
  core::meta::enable_if_t<
    has_method_observe_residual_when_solver_converged<
      observer_t, residual_t
     >::value
    >
  >{

  void operator()(const observer_t * obs,
		  const residual_t & resid) const{
    obs->observeResidualWhenSolverConverged(resid);
  }
};


}}}} //end namespace rompp::solvers::iterative::impl
#endif
