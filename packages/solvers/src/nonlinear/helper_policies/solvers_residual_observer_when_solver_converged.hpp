
#ifndef SOLVERS_IMPL_RESIDUAL_OBSERVER_WHEN_SOLVER_CONVERGED_HPP
#define SOLVERS_IMPL_RESIDUAL_OBSERVER_WHEN_SOLVER_CONVERGED_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../meta/solvers_is_legitimate_residual_observer_when_converged.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template<typename observer_t,
	 typename residual_t,
	 typename = void>
struct ResidualObserverWhenSolverConverged{

  static void evaluate(const observer_t * obs,
		       const residual_t & resid) {
    // no op
  }
};


template<typename observer_t,
	 typename residual_t>
struct ResidualObserverWhenSolverConverged<
  observer_t, residual_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::solvers::meta::is_legitimate_residual_observer_when_solver_converged<
      observer_t, residual_t
     >::value
    >
  >{

  static void evaluate(const observer_t * obs,
		       const residual_t & resid) {
    obs->observeResidualWhenSolverConverged(resid);
  }
};


}}}} //end namespace pressio::solvers::iterative::impl
#endif
